'''
  For problems with this algorithm, please either write a:
  - post in the Mantid forum <http://forum.mantidproject.org>
  - message in the "Contact us" page <http://www.mantidproject.org/Contact>
  - email to <mantid-help@mantidproject.org>

  Algorithm to global fit of QENS data to the jump-diffusion model by Teixeira
  Fitting model. In this case:
        Convolution( A*Resolution, x*Delta + (1-x)*TeixeiraWaterSQE ) + LinearBackground
    with 0<x<1 is the fraction of the elastic intensity
  Parameters DiffCoeff and Tau (residence time) of fit function TeixeiraWaterSQE are
  the same for all spectra. All other fitting  parameters are different for each spectrum
'''
from mantid.simpleapi import *
from mantid.api import DataProcessorAlgorithm, AlgorithmFactory, MatrixWorkspaceProperty, \
                       WorkspaceGroupProperty, Progress
from mantid.kernel import Direction
from mantid import logger
import re
from copy import copy
import sys
import numpy as np
from os.path import join as pjoin

class TeixeiraWaterFit(DataProcessorAlgorithm):

    def category(self):
        return "Workflow\\MIDAS"


    def summary (self):
        return "Calculates the nth moment of y(q,w)"

    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty("DataWorkspace", "", Direction.Input),
                             doc="Data Workspace to use.")
        self.declareProperty(MatrixWorkspaceProperty("ResolutionWorkspace", "", Direction.Input),
                             doc="Resolution Workspace to use.")
        self.declareProperty(name='EnergyMin', defaultValue=-0.5,
                             doc='Minimum energy in meV for fit')
        self.declareProperty(name='EnergyMax', defaultValue=0.5,
                             doc='Maximum energy in meV for fit')
        self.declareProperty(name='OutputName', defaultValue='Output',
                             doc='Name for output workspace')


    def PyExec(self):

        workflow_prog = Progress(self, start=0.0, end=1.0, nreports=20)

        data = self.getPropertyValue('DataWorkspace')
        logger.information('Data workspace %s' % data)
        data_ws = mtd[data]
        resolution = self.getPropertyValue('ResolutionWorkspace')
        logger.information('Resolution workspace %s' % resolution)
        resolution_ws = mtd[resolution]
        output = self.getProperty('OutputName').value

        # Find out the Q-values from the loaded data
        workflow_prog.report('Getting Q values')
        numHist = data_ws.getNumberHistograms()
        Qs=list()  # will store the Q-values
        for idx in range(numHist):
            detector = data_ws.getDetector(idx) 
            efixed = data_ws.getEFixed(detector.getID())  # in meV
            wavelength=9.044567/np.sqrt(efixed)  # in Angstroms
            usignTheta = 0.5 * data_ws.detectorTwoTheta(detector)
            Q = (4*np.pi/wavelength)*np.sin(usignTheta)  # in inverse Angstroms
            Qs.append(Q)
        nQ=len(Qs)
        logger.information('Sample %i Q values from %f to %f' % (nQ, Qs[0], Qs[nQ-1]))

        # Our initial guesses are DiffCoeff=1.0 and Tau=1.0
        workflow_prog.report('Getting initial guesses')
        single_model_template="""(composite=Convolution,NumDeriv=true;
        name=TabulatedFunction,Workspace=%s,WorkspaceIndex=0,Scaling=1,Shift=0,XScaling=1;
        (name=DeltaFunction,Height=1.5,Centre=0;
        name=TeixeiraWaterSQE,Q=_Q_,Height=1.0,Tau=1.0,DiffCoeff=1.0,Centre=0;
        ties=(f1.Centre=f0.Centre)));
        name=LinearBackground,A0=0,A1=0""" % resolution_ws
        single_model_template = re.sub('[\s+]', '', single_model_template)  # remove whitespaces and such

        # Do the fit using only these workspaces
        #selected_wi = [ 0, 3, 4, 7] # select a few workspace indexes
        selected_wi = range(len(Qs))  # use all spectra for the fit
        nWi=len(selected_wi)  # number of selected spectra for fitting

        # Energy range over which we do the fitting.
        minE = self.getProperty('EnergyMin').value
        maxE = self.getProperty('EnergyMax').value
        logger.information('Energy range from %f to %f' % (minE, maxE))

        # Create the string representation of the global model (for selected spectra):
        global_model="composite=MultiDomainFunction,NumDeriv=true;"
        for wi in selected_wi:
            Q=Qs[wi]
            single_model = single_model_template.replace("_Q_", str(Q))  # insert Q-value
            global_model += "(composite=CompositeFunction,NumDeriv=true,$domains=i;{0});\n".format(single_model)

        # Tie DiffCoeff and Tau since they are global parameters
        ties=['='.join(["f{0}.f0.f1.f1.DiffCoeff".format(di) for di in reversed(range(nWi))]),
        '='.join(["f{0}.f0.f1.f1.Tau".format(wi) for wi in reversed(range(nWi))]) ]
        global_model += "ties=("+','.join(ties)+')'  # tie Radius

        # Now relate each domain(i.e. spectrum) to each single-spectrum model
        domain_model = dict()
        domain_index = 0
        for wi in selected_wi:
            if domain_index == 0:
                domain_model.update({"InputWorkspace": data_ws.name(), "WorkspaceIndex": str(wi),
                             "StartX": str(minE), "EndX": str(maxE)})
            else:
                di = str(domain_index)
                domain_model.update({"InputWorkspace_"+di: data_ws.name(), "WorkspaceIndex_"+di: str(wi),
                             "StartX_"+di: str(minE), "EndX_"+di: str(maxE)})
            domain_index += 1

        # Invoke the Fit algorithm using global_model and domain_model:
        workflow_prog.report('Fitting')
        output_workspace = output + '_' + data_ws.name()
        Fit(Function=global_model,
            Output=output_workspace,
            CreateOutput=True,
            MaxIterations=500,
            **domain_model)

        # Save Q-dependencies of the optimized parameters
        parameters_workspace = mtd[output_workspace + "_Parameters"]
        aliases = {"f0.f1.f0.Height": "Ha", "f0.f1.f1.Height": "Hb", "f0.f1.f1.Tau": "tau", "f0.f1.f1.DiffCoeff": "diffCoeff"}
        optimals = dict()
        for alias in aliases.values():
            optimals[alias] = list()
        names = aliases.keys()
        for row in parameters_workspace:
        # name of the fitting parameter for this particular row
            matches = re.search("^f(\d+)\.(.*)", row['Name']) # for instance, f0.f0.f1.f0.Height
            if matches:
                iq, name = matches.groups()
                if name in names:
                    optimals[aliases[name]].append(row["Value"])
        # Calculate EISF as the fraction of the total intensity multiplying the Delta-Dirac
        # (The data is water at ambient conditions, thus EISF turns out very small)
        optimals["EISF"]=  [optimals["Ha"][i]/(optimals["Ha"][i]+optimals["Hb"][i]) for i in range(nWi)]
        optimals["qvalues"]=[Qs[wi] for wi in selected_wi]  # list the Q-values for the selected spectra
        # Save parameters Q-dependence to workspace "glofit_Qdependencies"
        dataY = optimals["EISF"] + optimals["diffCoeff"] + optimals["tau"]
        dataX = optimals["qvalues"] * 3
        CreateWorkspace(OutputWorkspace=output + '_Qdependencies',
	                    DataX=dataX,
                        UnitX="MomentumTransfer",
                        DataY=dataY,
                        NSpec=3,
                        WorkspaceTitle="Q-dependence of parameters",
                        VerticalAxisUnit="Text",
                        VerticalAxisValues=["EISF", "diffCoeff", "tau"],
                        EnableLogging=False)
        # Save parameters Q-dependence to workspace "glofit_Qdependencies"
        dataX = np.array(dataX) ** 2
        CreateWorkspace(OutputWorkspace=output + '_Q2dependencies',
	                    DataX=dataX,
                        UnitX="QSquared",
                        DataY=dataY, NSpec=3,
                        WorkspaceTitle="Q squared-dependence of parameters",
                        VerticalAxisUnit="Text",
                        VerticalAxisValues=["EISF", "diffCoeff", "tau"],
                        EnableLogging=False)

# Register algorithm with Mantid
AlgorithmFactory.subscribe(TeixeiraWaterFit)
