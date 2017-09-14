# -*- coding: utf-8 -*-
from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *


class SumToShim(PythonAlgorithm):
    """A Mantid implementation of the SumToShim algorithm"""
    def category(self):
        """What categories to file SumToShim under in the GUI menu."""
        return "SESANS"

    def PyInit(self):
        """The initial setup of the algorithm"""
        # rnum=20401 OutputWorkspace
        self.declareProperty("rnum", defaultValue=20401)
        self.declareProperty(
            WorkspaceProperty(name="OutputWorkspace",
                defaultValue="",
                direction=Direction.Output))

    def PyExec(self):
        """The actual running of the algorithm"""
        rnum = self.getProperty("rnum").value
        try:
            wtemp = Load(str(rnum), LoadMonitors=True)
            RebinToWorkspace('wtemp_1', 'wtemp_monitors_1',
                             PreserveEvents=False,
                             OutputWorkspace='wtemp_1')
            RebinToWorkspace('wtemp_2', 'wtemp_monitors_1',
                             PreserveEvents=False,
                             OutputWorkspace='wtemp_2')
            wtemp_1 = ConjoinWorkspaces('wtemp_monitors_1', 'wtemp_1')
            wtemp_2 = ConjoinWorkspaces('wtemp_monitors_2', 'wtemp_2')
        except:
            wtemp_monitors = Load(rnum)

        wtempShim = mtd['wtemp_monitors_1'] + mtd['wtemp_monitors_2']
        self.setProperty("OutputWorkspace", wtempShim)
        DeleteWorkspaces([wtempShim, "wtemp_monitors"])

AlgorithmFactory.subscribe(SumToShim)
