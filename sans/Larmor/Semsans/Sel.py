# -*- coding: utf-8 -*-
from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *


class Patterson(PythonAlgorithm):
    """A Mantid implementation of the Patterson algorithm"""
    def category(self):
        """What categories to file Patterson under in the GUI menu."""
        return "SESANS"

    def PyInit(self):
        """The initial setup of the algorithm"""
        # pol=wksp=Input const=10 OutputWorkspace
        self.declareProperty(
            WorkspaceProperty(name="pol",
                              defaultValue="",
                              direction=Direction.Input))
        self.declareProperty("const", defaultValue=10)
        self.declareProperty(
            WorkspaceProperty(name="OutputWorkspace",
                              defaultValue="",
                              direction=Direction.Output))

    def PyExec(self):
        """The actual running of the algorithm"""
        pol = self.getProperty("pol").value
        const = self.getProperty("const").value

        wksp = CloneWorkspace(pol)
        x = wksp.extractX()
        x = (x[:, 1:]+x[:, :-1])/2
        wksp = Logarithm(wksp)
        for i in range(x.shape[0]):
            wksp.setY(i, wksp.readY(i)/x[i]**2)
        wksp = ConvertUnits(wksp, "SpinEchoLength", EFixed=const)

        self.setProperty("OutputWorkspace", wksp)
        DeleteWorkspace(wksp)


AlgorithmFactory.subscribe(Patterson)
