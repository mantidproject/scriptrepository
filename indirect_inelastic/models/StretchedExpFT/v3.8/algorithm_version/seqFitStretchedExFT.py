from __future__ import print_function
import sys
import string
import re
import numpy as np
import mantid.api as mapi
import mantid.simpleapi as msapi
import mantid.kernel as mkrnl

class FitWorker(object):
    """
    Object in charge to do a sequential fit
    """

    def __init__(self, alg):
        """
        Object in charge of doing the minimization
        :param alg: Algorithm calling the object
        """
        self._alg = alg
        self._fittedQ = list()

    def _GetQsInQENSData(self):
        """
        Attempt to extract or compute the Q values from a workspace
        :param workspace: matrix workspace
        :return: python list of Q values
        :except: unable to extract the list of Q values
        """
        # Check if the vertical axis has units of Momentum transfer
        workspace = self._alg.getProperty("DataWorkspace").value
        number_spectra = workspace.getNumberHistograms()
        try:
            axis = workspace.getAxis(1)
        except:
            raise RuntimeError("Emtpy vertical axis")
        if axis.getUnit().unitID() == "MomentumTransfer":
            try:
                qvalues = axis.extractValues()
            except:
                raise RuntimeError("Unable to extract the Q values")
            # convert to point values if qvalues is histogram data
            if len(qvalues) == number_spectra+1:
                qvalues = (qvalues[1:]+qvalues[:-1])/2
            qvalues = qvalues.tolist()
        else:
            # Attempt to compute the Q values
            qvalues=list()
            try:
                for index in range(number_spectra):
                    detector = workspace.getDetector(index)
                    efixed = workspace.getEFixed(detector.getID())  # in meV
                    wavelength=9.044567/np.sqrt(efixed)  # in Angstroms
                    usignTheta = 0.5 * workspace.detectorTwoTheta(detector)
                    qvalue = (4*np.pi/wavelength)*np.sin(usignTheta)  # in inverse Angstroms
                    qvalues.append(qvalue)
            except:
                raise RuntimeError("Unable to compute the Q values")
        return qvalues

    def _prefix(self):
        """
        Return the root name for all output workspaces
        :return: prefix string
        """
        return "seq"+self._alg.getProperty("OutputName").value

    def groupWorkspaces(self, output_prefixes):
        """
        Gather the workspaces output from Fit algorithm into three WorkspaceGroups
        :param inputWs: list of workspace root names to be grouped
        """
        for suffix in ("_NormalisedCovarianceMatrix", "_Parameters", "_Workspace"):
            inputWs = [prefix+suffix for prefix in output_prefixes]
            msapi.GroupWorkspaces(InputWorkspaces=inputWs, OutputWorkspace=self._prefix()+suffix)

    def cleanPreviousFitResults(self):
        """
        Remove previous existing group workspace.
        """
        for suffix in ("_NormalisedCovarianceMatrix", "_Parameters", "_Workspace"):
            outputWorkspace = self._prefix()+suffix
            if mapi.AnalysisDataService.doesExist(outputWorkspace):
                workspaceGroup = mapi.mtd[outputWorkspace]
                for workspace in workspaceGroup:
                    print(workspace.name())
                    mapi.AnalysisDataService.remove(workspace.name())
                mapi.AnalysisDataService.remove(outputWorkspace)
        for suffix in ("_Qdependence", "_QQdependence"):
            outputWorkspace = self._prefix()+suffix
            if mapi.AnalysisDataService.doesExist(outputWorkspace):
                mapi.AnalysisDataService.remove(outputWorkspace)
        self._fittedQ = list()

    def doFit(self):
        """
        Carry out sequential fit
        :return:
        """
        self.cleanPreviousFitResults()
        qvals = self._GetQsInQENSData()
        mm = self._alg._extractMaxIterations(self._alg.getProperty("seqMinimizer").value)
        selectedwi = self._alg.getProperty("SelectedWI").value
        output_prefixes = list()
        for wi in selectedwi:
            wi = int(wi)  # Avoid numpy.int64 to int pitfall
            output_prefix = self._prefix()+str(wi)
            # fitreport is a tuple with the following contents, by index
            # 0: success info-string
            # 1: Chi2
            # 2: reference to the table workspace containing correlations between pair of parameters
            # 3: reference to the table workspace containing optimal parameter values
            # 4: reference to MatrixWorkspace containing the data, fit, and difference curves
            print(self._alg._applyGuess(wi=wi))
            try:
                fitreport = msapi.Fit(self._alg._applyGuess(wi=wi),
                                      InputWorkspace=self._alg.getPropertyValue("DataWorkspace"),
                                      WorkspaceIndex=wi, CreateOutput=1,
                                      startX=self._alg.getProperty("EnergyRange").value[0],
                                      endX=self._alg.getProperty("EnergyRange").value[1],
                                      Minimizer=mm["Minimizer"], MaxIterations=mm["MaxIterations"],
                                      OutputCompositeMembers=True, ConvolveMembers=True,
                                      Output=output_prefix)
                self._alg._updateGuess(fitreport)
                output_prefixes.append(output_prefix)
                self._fittedQ.append(qvals[wi])
            except:
                break
        output_prefixes.reverse()  # recall we fitted from highest to lowest workspace index
        self._fittedQ.reverse()
        self.groupWorkspaces(output_prefixes)

    def extractOptimalParameter(self, parname, compare_method="=="):
        """
        Extract value and indeterminancy for a parameter from the parameter WorkspaceGroup
        :param parname: name of the parameter of interest
        :return: list of values and indeterminacies of the parameter
        """
        values = list()
        errors = list()
        wsg = mapi.mtd[self._prefix()+"_Parameters"]  # WorkspaceGroup containing TableWorkspace parameters
        # Iterate over the workspaces
        for parmtable in wsg:
            for row in parmtable:
                comparison = False
                if compare_method == "==":
                    comparison = (parname == row['Name'])
                elif compare_method == "in":
                    comparison = (parname in row['Name'])
                if comparison:
                    values.append(row['Value'])
                    errors.append(row['Error'])
                    break
        return values, errors
    
    def dependQ(self):
        """
        Collect Q-dependence on the optimal parameters of interest
        :return: table workspaces for Q and Q^2 dependence
        """
        names = list()
        Y = list()
        E = list()
        for prop, key in self._alg._propvalue2key.items():
            names.append(prop)
            values, errors = self.extractOptimalParameter(key.replace('_','.'))
            Y += values
            E += errors
        values, errors = self.extractOptimalParameter("Cost", compare_method="in")
        names.append("Cost"); Y += values; E += errors
        nhistograms = len(self._alg._propvalue2key) + 1
        Q = np.tile(self._fittedQ, nhistograms)  # only selected ones
        msapi.CreateWorkspace(DataX=Q, dataY=Y, dataE=E,
                              NSpec=nhistograms, UnitX="MomentumTransfer",
                              WorkspaceTitle="Q dependence of optimal parameters",
                              VerticalAxisUnit="Text", VerticalAxisValues=names,
                              OutputWorkspace=self._prefix()+"_Qdependence")
        msapi.CreateWorkspace(DataX=Q*Q, dataY=Y, dataE=E,
                              NSpec=nhistograms, UnitX="QSquared",
                              WorkspaceTitle="Q-squared dependence of optimal parameters",
                              VerticalAxisUnit="Text", VerticalAxisValues=names,
                              OutputWorkspace=self._prefix()+"_QQdependence")
        # Sort by ascending Q. This is neccessary if SelectedWI not in ascending order
        msapi.SortXAxis(InputWorkspace=self._prefix()+"_Qdependence", OutputWorkspace=self._prefix()+"_Qdependence")
        msapi.SortXAxis(InputWorkspace=self._prefix()+"_QQdependence", OutputWorkspace=self._prefix()+"_QQdependence")
