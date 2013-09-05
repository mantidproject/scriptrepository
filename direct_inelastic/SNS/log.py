"""
Quicly check logs for multiple runs
"""
from numpy import *
from string import *
runs=arange(22579,23388,1)
prefix="/SNS/SEQ/IPTS-6179/data/SEQ_"
suffix="_event.nxs"
print "  Run   Temp       Angle    PC 	HIST	Energy"
for r in runs:
	CreateSingleValuedWorkspace(OutputWorkspace="w",DataValue="0")
	LoadNexusLogs(Workspace="w",Filename=prefix+str(r)+suffix,OverwriteLogs=True)
	w=mtd["w"]
	runinfo=w.getRun()
	temp=array(runinfo["SampleTemp"].value).mean()
	angle=array(runinfo["CCR13VRot"].value).mean()
	pc=array(runinfo["proton_charge"].value).sum()
	energy=array(runinfo["EnergyRequest"].value).mean()
	LoadNexusMonitors(OutputWorkspace='w1',Filename=prefix+str(r)+suffix)
	nhistmon=mtd['w1'].getNumberHistograms()
	print r,temp,angle,pc,nhistmon,energy
	DeleteWorkspace('w')
	DeleteWorkspace('w1')
