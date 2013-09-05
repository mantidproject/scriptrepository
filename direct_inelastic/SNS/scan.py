"""
Rocking scan (+ and - directions)
"""
import numpy as np

ang_name='CCR12Rot'
ang_deriv='CCR12RotDeriv'
step=0.3
filename="/home/andrei/Desktop/ARCS_12924_event.nxs"
groupFile='/home/andrei/Desktop/group1.xml'

Load(Filename=filename,OutputWorkspace="w")
FilterBadPulses(InputWorkspace='w',OutputWorkspace='w')
AddLogDerivative(InputWorkspace='w',LogName=ang_name,NewLogName=ang_deriv)
w=mtd["w"]
psi=np.array(w.getRun()[ang_name].value)
minpsi=min(psi)
maxpsi=max(psi)
ang_list_plus=[]
ang_list_minus=[]
int_plus=[]
int_minus=[]
e_plus=[]
e_minus=[]


for i in np.arange(minpsi,maxpsi-step,step):
	print i
	FilterByLogValue(InputWorkspace="w",OutputWorkspace="wi",LogName=ang_name,MinimumValue=i,MaximumValue=i+step)
	if (mtd["wi"].getNumberEvents()>1):
		FilterByLogValue(InputWorkspace="wi",OutputWorkspace="wip",LogName=ang_deriv,MinimumValue=0.,MaximumValue=1e100)
		FilterByLogValue(InputWorkspace="wi",OutputWorkspace="wim",LogName=ang_deriv,MinimumValue=-1e100,MaximumValue=0)
		GroupDetectors(InputWorkspace="wip",OutputWorkspace="wip",MapFile=groupFile)
		GroupDetectors(InputWorkspace="wim",OutputWorkspace="wim",MapFile=groupFile)
		wip=mtd["wip"]
		psi=np.array(wip.getRun()[ang_name].value)
		pc=wip.getRun().getProtonCharge()
		if (pc>0):
			ang_list_plus.append(psi.mean())
			int_plus.append((wip.readY(0)[0])/pc)
			e_plus.append((wip.readE(0)[0])/pc)
		wim=mtd["wim"]
		psi=np.array(wim.getRun()[ang_name].value)
		pc=wim.getRun().getProtonCharge()
		if (pc>0):
			ang_list_minus.append(psi.mean())
			int_minus.append((wim.readY(0)[0])/pc)
			e_minus.append((wim.readE(0)[0])/pc)


CreateWorkspace(OutputWorkspace="Pscan",DataX=np.array(ang_list_plus),DataY=np.array(int_plus),DataE=np.array(e_plus))
CreateWorkspace(OutputWorkspace="Mscan",DataX=np.array(ang_list_minus),DataY=np.array(int_minus),DataE=np.array(e_minus))

all_scans=plotSpectrum(("Pscan","Mscan"),0,True)
