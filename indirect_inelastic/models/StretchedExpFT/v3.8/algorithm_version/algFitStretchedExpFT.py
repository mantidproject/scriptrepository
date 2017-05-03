'''
  For problems with this algorithm, please either write a:
  - post in the Mantid forum <http://forum.mantidproject.org>
  - message in the "Contact us" page <http://www.mantidproject.org/Contact>
  - email to <mantid-help@mantidproject.org>
'''

import mantid.api as mapi
import mantid.simpleapi as msapi
import mantid.kernel as mkrnl
import string
import os
import sys
import re
import numpy as np
from collections import OrderedDict
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import seqFitStretchedExFT
reload(seqFitStretchedExFT)
#import gloFitStretchedExFTBnew
#reload(gloFitStretchedExFTBnew)

class FitStretchedExpFT(mapi.PythonAlgorithm):

    _scheme2mod = {"Sequential": seqFitStretchedExFT}  # , "Global": gloFitStretchedExFTnew}
    #  String representation of the fit model, with templated parameters and ties

    _fitTpl = """
(composite=Convolution,FixResolution=false,NumDeriv=true;
   name=TabulatedFunction,Workspace=${resolution},WorkspaceIndex=$wi,
        Scaling=${f0_f0_Scaling},Shift=0,XScaling=1,
        ties=(XScaling=1);
  (composite=CompositeFunction,NumDeriv=true;
     name=DeltaFunction,Height=${f0_f1_f0_Height},Centre=0,
          constraints=(0<Height<1),ties=(Centre=0${fixf0_f1_f0_Height});
     name=StretchedExpFT,Tau=${f0_f1_f1_Tau},Beta=${f0_f1_f1_Beta},Centre=0,
          constraints=(0<Tau,0<Beta),
          ties=(Centre=0${fixf0_f1_f1_Beta});
     ties=(f1.Height=1-f0.Height)
   )
);
name=LinearBackground,A0=${f1_A0},A1=${f1_A1}
"""

    # Relation between algorithm property fix and
    # (property parameter, global parameter name, local parameter name)
    _propfix2key = {"FixEI": ("EISF", "f0_f1_f0_Height", "Height"),
                    "FixBeta": ("Beta", "f0_f1_f1_Beta", "Beta")
                   }

    # Relation between algorithm property value and global parameter name
    _propvalue2key = OrderedDict([("Intensity", "f0_f0_Scaling"),
                                  ("EISF", "f0_f1_f0_Height"),
                                  ("Tau", "f0_f1_f1_Tau"),
                                  ("Beta", "f0_f1_f1_Beta"),
                                  ("A0", "f1_A0"),
                                  ("A1", "f1_A1")
                                  ])

    def category(self):
        return "Optimization\\FitAlgorithms"

    def summary (self):
        return "Fit QENS data resulting from a streched exponential decay in the time domain"

    def PyInit(self):

        self.declareProperty(MatrixWorkspaceProperty("DataWorkspace", "",
                                                     direction=mkrnl.Direction.Input),
                             doc="Data Workspace to use.")
        self.declareProperty(MatrixWorkspaceProperty("ResolutionWorkspace", "",
                                                     direction=mkrnl.Direction.Input),
                             doc="Resolution Workspace to use")

        # Options
        fitEngineOptions = "Fit Engine Options"
        fitSchemes = ["Sequential", "Global"]
        self.declareProperty("FitScheme", "Sequential", mkrnl.StringListValidator(fitSchemes),
                             direction=mkrnl.Direction.Input, doc="Type of fit scheme")
        self.setPropertyGroup("FitScheme", fitEngineOptions)
        self.declareProperty("seqMinimizer", "FABADA,ConvergenceCriteria=0.02,MaxIterations=2000", direction=mkrnl.Direction.Input,
                             doc="Minimizer for sequential fit")  # Levenberg-Marquardt,MaxIterations=2000
        self.setPropertyGroup("seqMinimizer", fitEngineOptions)
        fitSchemeCondition = mkrnl.EnabledWhenProperty("FitScheme",mkrnl.PropertyCriterion.IsEqualTo,"Global")
        self.declareProperty("gloMinimizer", "FABADA,ConvergenceCriteria=0.02,MaxIterations=1000",
                             direction=mkrnl.Direction.Input,
                             doc="Minimizer for global fit")
        self.setPropertyGroup("gloMinimizer", fitEngineOptions)
        self.setPropertySettings("gloMinimizer", fitSchemeCondition)
        self.declareProperty(mkrnl.IntArrayProperty("SelectedWI", values=[],
                                                      direction=mkrnl.Direction.Input),
                             doc="Selected workspace indexes to fit, in descending order, leave empty to select all")
        self.setPropertyGroup("SelectedWI", fitEngineOptions)
        er_validator = mkrnl.FloatArrayLengthValidator(2)
        self.declareProperty(mkrnl.FloatArrayProperty("EnergyRange", values=[-0.1, 0.1],
                                                      validator=er_validator, direction=mkrnl.Direction.Input),
                             doc="Fitting energy range (in meV)")
        self.setPropertyGroup("EnergyRange", fitEngineOptions)

        # Initial Guess
        initialGuessTitle = "Initial guess for spectrum with highest Q"
        self.declareProperty(name="Intensity", defaultValue=1.0, direction=Direction.Input,
                             doc="Overall intensity, area under the model")
        self.setPropertyGroup("Intensity", initialGuessTitle)
        self.declareProperty(name="EISF", defaultValue=0.5, direction=Direction.Input,
                             doc="List of elastic incoherent structure factors, between 0 and 1.")
        self.setPropertyGroup("EISF", initialGuessTitle)
        self.declareProperty(name="FixEI", defaultValue=False, direction=Direction.Input,
                             doc="Fix intensity of the elastic line for all spectra?")
        self.setPropertyGroup("FixEI", initialGuessTitle)
        self.declareProperty(name="Tau", defaultValue=50.0,
                             doc="relaxation time (ps) for the highest selected workspace index")
        self.setPropertyGroup("Tau", initialGuessTitle)
        self.declareProperty(name="Beta", defaultValue=1.0,
                             doc="Beta exponent for the highest selected workspace index")
        self.setPropertyGroup("Beta", initialGuessTitle)
        self.declareProperty(name="FixBeta", defaultValue=False, direction=Direction.Input,
                             doc="Fix value of Beta for all spectra?")
        self.setPropertyGroup("FixBeta", initialGuessTitle)
        self.declareProperty(name="A0", defaultValue=0.0,
                             doc="Flat component of the linear background")
        self.setPropertyGroup("A0", initialGuessTitle)
        self.declareProperty(name="A1", defaultValue=0.0,
                             doc="Slope of the linear background")
        self.setPropertyGroup("A1", initialGuessTitle)
        self.declareProperty(name="OutputName", defaultValue="fit", direction=Direction.Input,
                             doc="prefix for output workspaces")
        self.declareProperty(name="DryRun", defaultValue=False, direction=Direction.Input,
                             doc="Output the string representation of the fit function using the initial guess but do not carry out any fitting")

    def PyExec(self):
        workflow_prog = Progress(self, start=0.0, end=1.0, nreports=20)

        # Substitute the resolution workspace name
        template = string.Template(self._fitTpl)
        self._fitTpl = template.safe_substitute({"resolution":self.getPropertyValue('ResolutionWorkspace')})

        # Find out selected workspaces and sort from lowest to highest
        selectedwi = self.getProperty("SelectedWI").value
        if not selectedwi.any():
            selectedwi = np.array(range(mapi.mtd[self.getPropertyValue("DataWorkspace")].getNumberHistograms()))
        selectedwi.sort()
        selectedwi = selectedwi[::-1] # from Highest to lowest
        self.setProperty("SelectedWI", selectedwi)
        self._applyFixes()
        if self.getProperty("DryRun").value:
            self._fitTpl = self._applyGuess(wi=self.getProperty("SelectedWI").value[0])  # LAST of the selected workspaces
            mkrnl.logger.notice("FUNCTION STRING FOR THE LAST SELECTED SPECTRUM\n########\n"+self._fitTpl+"\n########")
        else:
            fit_scheme_mod = FitStretchedExpFT._scheme2mod[self.getProperty("FitScheme").value]  # seq or glob module

            print dir(fit_scheme_mod)
            print os.path.abspath(fit_scheme_mod.__file__)
            fitter = fit_scheme_mod.FitWorker(self)  # call sequential or global fit
            fitter.doFit()
            fitter.dependQ()  # create workspaces with Q-dependence

    def _extractMaxIterations(self, minimizer_info, maxiterations=1000):
        """
        Split the minimizer info string into the maxiterations and the rest
        :param minimizer_info: string containing all info about the minimizer, plus MaxIterations
        :return: two entry dictionary with MaxIterations and rest of the minimizer info string
        :except: no maxiterations found, set to the default
        """
        try:
            found = re.search('(,\s*MaxIterations\s*=\s*\d+)', minimizer_info).group(1)
            mi = int(re.search('(\d+)', found).group(1))
            minimizer = minimizer_info.replace(found, "")
            return {"MaxIterations": mi, "Minimizer": minimizer}
        except AttributeError:
            mkrnl.logger.error(
                "No MaxIterations in minimizer infor string. Setting MaxIterations=" + str(maxiterations))
            return {"MaxIterations": maxiterations, "Minimizer": minimizer_info}

    def _applyFixes(self):
        """
        Introduce ties for fixed parameters and substitute with initial guess values
        """
        template=string.Template(self._fitTpl)
        # Introduce the ties
        kwargs={}
        for propf in ("FixEI", "FixBeta"):
            propv, key, local = self._propfix2key[propf]
            if self.getProperty(propf).value:
                value = str(self.getProperty(propv).value)
                kwargs[key] = value  # substitute global parameter name with value
                kwargs["fix"+key] = ",{0}={1}".format(local,value)  # insert tie
            else:
                kwargs["fix"+key] = ""  # remove tie template
        self._fitTpl = template.safe_substitute(**kwargs)

    def _applyGuess(self, wi=-1):
        """
        Substitute the initial guess values
        :return: function string with initial guess values
        """
        template = string.Template(self._fitTpl)
        # Find initial guess values
        kwargs = {key:str(self.getProperty(prop).value) for prop,key in self._propvalue2key.items()}
        if wi >= 0:
            kwargs.update({"wi":str(wi)})
        return template.safe_substitute(**kwargs)  # this self._fitTpl is owned by object, not class

    def _updateGuess(self, fitreport):
        """
        Scan a table workspace, and update the algorithm properties of the initial guess.
        Note: special treatment for parameter tau because of its exponential Q-dependence
        :param fitreport: output of algorithm Fit
        :return:
        """
        # Find out the table containing the parameters
        parmtable = None
        for item in fitreport:
            if type(item)==type(msapi.CreateEmptyTableWorkspace()):
                if item.columnCount()==3:
                    parmtable = item
                    break
        # Iterate over the table contents
        for row in parmtable:
            name = row['Name'].replace('.', '_')  # name of the fitting parameter
            for propvalue, key in self._propvalue2key.items():
                if name == key:
                    factor = 1.0
                    if name == "f0_f1_f1_Tau":
                        factor = 2.0  # Q increases fast when going towards smaller Q values
                    self.setProperty(propvalue, factor*row['Value'])

# Register algorithm with Mantid
AlgorithmFactory.subscribe(FitStretchedExpFT)
