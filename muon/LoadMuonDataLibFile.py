# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
from mantid.api import *
from mantid.kernel import *
import matplotlib.pyplot as plt
import numpy as np
import h5py

# Calibrate muon dead time from MuonDataLib histogram file

class LoadMuonDataLibFile(PythonAlgorithm):
    def category(self):
        return "Muon"

    def PyInit(self):
        self.declareProperty(FileProperty(name="InputFile",defaultValue="",action=FileAction.Load,extensions=[".nxs"]))
        self.declareProperty("Instrument","",StringListValidator(["","HIFI","MUSR"])) # ,help="Instrument definition to add"
        self.declareProperty(WorkspaceProperty("OutputWorkspace","",Direction.Output))

    def PyExec(self):
        fh=h5py.File(self.getProperty("InputFile").value,"r")
        dat=fh["/raw_data_1/instrument/detector_1/counts"]
        print (dat.shape)
        
        nper=dat.shape[0]
        if nper>1:
            raise NotImplementedError("Periods in use")
        nvec=dat.shape[1]
        ntbins=dat.shape[2]
        tbins=fh["/raw_data_1/instrument/detector_1/raw_time"]
        goodfrm=fh["/raw_data_1/good_frames"][0]
        
        instr=self.getProperty("Instrument").value
        if instr:
            w=CreateSimulationWorkspace(Instrument=instr,BinParams=(0,1,ntbins))
        else:
            w=WorkspaceFactory.create("Workspace2D",XLength=len(tbins),YLength=ntbins,NVectors=nvec)
        for i in range(nvec):
            w.dataX(i)[:]=tbins
            w.dataY(i)[:]=dat[0,i,:]
            w.dataE(i)[:]=np.sqrt(dat[0,i,:])
        w.getRun().addProperty("goodfrm",IntPropertyWithValue("goodfrm",int(goodfrm)),replace=True)
        
        self.setProperty("OutputWorkspace",w)
        
AlgorithmFactory.subscribe(LoadMuonDataLibFile)
