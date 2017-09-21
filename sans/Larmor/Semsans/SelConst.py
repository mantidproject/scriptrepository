# -*- coding: utf-8 -*-
import numpy as np
from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *


class SelConst(PythonAlgorithm):
    """A Mantid implementation of the SelConst algorithm"""
    def category(self):
        """What categories to file SelConst under in the GUI menu."""
        return "SESANS"

    def PyInit(self):
        """The initial setup of the algorithm"""
        # Tube1=wksp=Input Tube2=wksp=Input Tube3=wksp=Input dist=4.0 thickness=5e-3
        self.declareProperty(
            WorkspaceProperty(name="Tube1",
                defaultValue="",
                direction=Direction.Input))
        self.declareProperty(
            WorkspaceProperty(name="Tube2",
                defaultValue="",
                direction=Direction.Input))
        self.declareProperty(
            WorkspaceProperty(name="Tube3",
                defaultValue="",
                direction=Direction.Input))
        self.declareProperty("dist", defaultValue=4.0)
        self.declareProperty("thickness", defaultValue=0.005)
        self.declareProperty("show_fits", defaultValue=False)
        self.declareProperty("show_quality", defaultValue=False)
        self.declareProperty("result", defaultValue=0.0, direction=Direction.Output)

    def PyExec(self):
        """The actual running of the algorithm"""
        Tube1 = self.getProperty("Tube1").value
        Tube2 = self.getProperty("Tube2").value
        Tube3 = self.getProperty("Tube3").value
        dist = self.getProperty("dist").value
        thickness = self.getProperty("thickness").value
        show_fits = self.getProperty("show_fits").value
        show_quality = self.getProperty("show_quality").value

        freqs = []
        errs = []
        for run in [Tube1, Tube2, Tube3]:
            wsTemp = FFT(run)
            x = wsTemp.readX(2)
            y = wsTemp.readY(2)
            max_arg = np.nanargmax(y)
            amp = y[max_arg]
            self.log().notice(str(x[max_arg]))
            self.log().notice(str(amp))

            model = "name=UserFunction,Formula=e*cos(x*f),f={},e={}"
            Fit(Function=model.format(x[max_arg]*2*np.pi, amp),
                InputWorkspace=run, StartX=3, EndX=7.5, CreateOutput=True)

            result = mtd[run.getName()+"_Parameters"]

            freqs.append(np.abs(result.column(1)[1]))
            errs.append(np.abs(result.column(1)[1]))
            DeleteWorkspace(wsTemp)
            if not show_fits:
                DeleteWorkspace(run.getName()+"_Parameters")
                DeleteWorkspace(run.getName()+"_NormalisedCovarianceMatrix")
                DeleteWorkspace(run.getName()+"_Workspace")

        xs = np.arange(3)*thickness
        wlin = CreateWorkspace(xs, freqs, errs)
        Fit(Function='name=LinearBackground',
            InputWorkspace='wlin', CreateOutput=True)
        result = 0.1 * mtd["wlin_Parameters"].column(1)[1] * dist / 2 / np.pi
        if not show_quality:
            DeleteWorkspace("wlin_Parameters")
            DeleteWorkspace("wlin_NormalisedCovarianceMatrix")
            DeleteWorkspace("wlin_Workspace")
            DeleteWorkspace("wlin")
        self.setProperty("result", result)



AlgorithmFactory.subscribe(SelConst)
