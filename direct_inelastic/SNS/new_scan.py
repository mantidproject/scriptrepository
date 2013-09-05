"""
Updated version of rocking curve scan.py, using FilterEvents
"""
import numpy as np
ang_name='CCR12Rot'
ang_deriv='CCR12RotDeriv'

t_off=0.0
step=0.2
ang_bins=np.arange(160-step/2,169+step/1.5,step)

runs=range(12956,12966+1)			#list of runs that have the same condition
datadir="/SNS/ARCS/IPTS-4041/0/"		#Data directory	
run_num = str(runs[0])
filename=datadir+run_num+"/NeXus/ARCS_"+run_num+"_event.nxs"
LoadEventNexus(filename,OutputWorkspace="w_all")

for run in runs[1:]:
	run_num = str(run)
	filename=datadir+run_num+"/NeXus/ARCS_"+run_num+"_event.nxs"
	LoadEventNexus(filename,OutputWorkspace="temp_ws")
	Plus("w_all","temp_ws",OutputWorkspace="w_all")	
FilterBadPulses(InputWorkspace='w_all',OutputWorkspace='w')
ChangeLogTime(InputWorkspace='w',OutputWorkspace='w',LogName=ang_name,TimeOffset=t_off)
AddLogDerivative(InputWorkspace='w',LogName=ang_name,NewLogName=ang_deriv)
# Get mask if needed
LoadNexus(Filename=r'/SNS/users/vua/Mantid/Angle_binning_reduction/mask5_exclude.nxs',OutputWorkspace='MaskWorkspace')
MaskDetectors("w",MaskedWorkspace="MaskWorkspace")
w=mtd["w"]
rate = np.absolute(np.array(w.getRun()[ang_deriv].value)).mean()
print "Rate(deg/sec)):",rate
ang_list_plus=[]
ang_list_minus=[]
int_plus=[]
int_minus=[]
e_plus=[]
e_minus=[]
pc_plus=[]
pc_minus=[]


GenerateEventsFilter(InputWorkspace='w',OutputWorkspace='splp',SplittersInformationWorkspace='info',LogName=ang_name,MinimumLogValue='159.9',MaximumLogValue='170',LogValueInterval='0.2',FilterLogValueByChangingDirection='Increase')
FilterEvents(InputWorkspace='w',OutputWorkspaceBaseName='wsplit_pos',SplittersInformationWorkspace='info',InputSplittersWorkspace='splp',FilterByPulseTime='1',GroupWorkspaces='1')
DeleteWorkspace('wsplit_pos_-1')
SumSpectra(InputWorkspace='wsplit_pos',OutputWorkspace='desum',IncludeMonitors='0')


w=mtd['desum']
for i in range(w.getNumberOfEntries()):
	wip=w.getItem(i)
	pc=wip.getRun().getProtonCharge()
	if (pc>0):
		psi=np.array(wip.getRun()[ang_name].value)
		ang_list_plus.append(psi.mean())
		cts = wip.getNumberEvents()
		int_plus.append(cts/pc)
		e_plus.append(np.sqrt(cts)/pc)
		pc_plus.append(pc)

GenerateEventsFilter(InputWorkspace='w',OutputWorkspace='splm',SplittersInformationWorkspace='info',LogName=ang_name,MinimumLogValue='159.9',MaximumLogValue='170',LogValueInterval='0.2',FilterLogValueByChangingDirection='Decrease')
FilterEvents(InputWorkspace='w',OutputWorkspaceBaseName='wsplit_neg',SplittersInformationWorkspace='info',InputSplittersWorkspace='splm',FilterByPulseTime='1',GroupWorkspaces='1')
DeleteWorkspace('wsplit_neg_-1')
SumSpectra(InputWorkspace='wsplit_neg',OutputWorkspace='desum',IncludeMonitors='0')


w=mtd['desum']
for i in range(w.getNumberOfEntries()):
	wip=w.getItem(i)
	pc=wip.getRun().getProtonCharge()
	if (pc>0):
		psi=np.array(wip.getRun()[ang_name].value)
		ang_list_minus.append(psi.mean())
		cts = wip.getNumberEvents()
		int_minus.append(cts/pc)
		e_minus.append(np.sqrt(cts)/pc)
		pc_minus.append(pc)
		
CreateWorkspace(OutputWorkspace="Pscan",DataX=np.array(ang_list_plus),DataY=np.array(int_plus),DataE=np.array(e_plus))
CreateWorkspace(OutputWorkspace="Mscan",DataX=np.array(ang_list_minus),DataY=np.array(int_minus),DataE=np.array(e_minus))
CreateWorkspace(OutputWorkspace="Ppc",DataX=np.array(ang_list_minus),DataY=np.array(pc_plus))
CreateWorkspace(OutputWorkspace="Mpc",DataX=np.array(ang_list_minus),DataY=np.array(pc_minus))

all_scans=plotSpectrum(("Pscan","Mscan"),0,True)
pc_scans=plotSpectrum(("Ppc","Mpc"),0)

xp=np.array(ang_list_plus)
xm=np.array(ang_list_minus)
yp=np.array(int_plus)
ym=np.array(int_minus)
cp=(xp*yp).sum()/yp.sum()
cm=(xm*ym).sum()/ym.sum()
print "Center-of-Mass - Plus:",cp,"Minus:",cm,"Difference:",cp-cm
print "Estimated time offset:", t_off + (cp-cm)/rate/2
