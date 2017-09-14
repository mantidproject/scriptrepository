# -*- coding: utf-8 -*-
from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *


class Polarisation(PythonAlgorithm):
    def category(self):
        return "SESANS"

    def PyInit(self):
        self.declareProperty("RunNumber", defaultValue=0)
        self.declareProperty(
            FileProperty(
                name="MaskFile", action=FileAction.Load,
                defaultValue=r'\\isis\inst$\NDXLARMOR\User\Users\Edler\May_2017\Mask_Tube39.xml'))
        self.declareProperty(
            WorkspaceProperty(name="OutputWorkspace",
                              defaultValue="",
                              direction=Direction.Output))

    def PyExec(self):
        w1 = Load(str(self.getProperty("RunNumber").value), LoadMonitors=True)
        w1mon = ExtractSingleSpectrum("w1_monitors", 0)
        w1 = ConvertUnits("w1", "Wavelength", AlignBins=1)
        w1mon = ConvertUnits(w1mon, "Wavelength")
        w1 = Rebin(w1, "1.8,0.1,8", PreserveEvents=False)
        w1mon = Rebin(w1mon, "1.8,0.1,8")
        w1 = w1 / w1mon

        mask = self.getProperty("MaskFile").value
        Mask_Tube = LoadMask("LARMOR", mask)
        MaskDetectors(w1, MaskedWorkspace=Mask_Tube)
        w1 = SumSpectra(w1)

        up = mtd["w1_1"]
        dn = mtd["w1_2"]
        pol = (up-dn)/(up+dn)
        self.setProperty("OutputWorkspace", pol)
        DeleteWorkspaces([pol, "Mask_Tube", "w1_1", "w1_2", "w1_monitors_1",
                          "w1_monitors_2", "w1mon_1", "w1mon_2"])

AlgorithmFactory.subscribe(Polarisation)
