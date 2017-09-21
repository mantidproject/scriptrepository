# -*- coding: utf-8 -*-
from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *


class Polarisation(PythonAlgorithm):
    def category(self):
        return "SESANS"

    def PyInit(self):
        self.declareProperty("RunNumber", defaultValue=0)
        self.declareProperty("Binning", defaultValue="1.8,0.1,8.0")
        self.declareProperty(
            WorkspaceProperty(
                name="MaskFile",
                direction=Direction.Input,
                defaultValue=r'MaskWorkspace'))
        self.declareProperty(
            WorkspaceProperty(name="OutputWorkspace",
                              defaultValue="",
                              direction=Direction.Output))

    def PyExec(self):
        w1 = Load(str(self.getProperty("RunNumber").value), LoadMonitors=True)
        w1mon = ExtractSingleSpectrum("w1_monitors", 0)
        w1 = ConvertUnits("w1", "Wavelength", AlignBins=1)
        w1mon = ConvertUnits(w1mon, "Wavelength")
        
        binning = self.getProperty("Binning").value
        w1 = Rebin(w1, binning, PreserveEvents=False)
        w1mon = Rebin(w1mon, binning)
        w1 = w1 / w1mon

        Mask_Tube = self.getProperty("MaskFile").value
        MaskDetectors(w1, MaskedWorkspace=Mask_Tube)
        w1 = SumSpectra(w1)

        up = mtd["w1_1"]
        dn = mtd["w1_2"]
        pol = (up-dn)/(up+dn)
        self.setProperty("OutputWorkspace", pol)
        DeleteWorkspaces([pol, "w1_1", "w1_2", "w1_monitors_1",
                          "w1_monitors_2", "w1mon_1", "w1mon_2"])

AlgorithmFactory.subscribe(Polarisation)
