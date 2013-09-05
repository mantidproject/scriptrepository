"""
Simplified version of reduction using DGSReduction algorithm
"""
#!/usr/bin/env python

import os
import sys
from numpy import *

mantid_root = "/opt/mantidnightly"
mantid_bin = sys.path.append(os.path.join(mantid_root, "bin"))

from mantid.simpleapi import *


from string import *

def ShiftTime(WName,lg_name):
	"""
	shift the time in a given log to match the time in the proton charge log"
	"""
	H_IN = mtd[WName]
	PC =  H_IN.getRun()['proton_charge'].firstTime()
	#print "P="+str(PC)+"\n"
	P =  H_IN.getRun()[lg_name].firstTime()
	#print "P="+str(P)+"\n"
	Tdiff = PC-P
	Tdiff_num = Tdiff.total_milliseconds()*1E-3
	#print "Tdiff="+str(Tdiff_num)+"\n"
	ChangeLogTime(InputWorkspace=WName, OutputWorkspace = WName, LogName = lg_name, TimeOffset = Tdiff_num)




runs=arange(35385,35520)		#First last+1
#data directory
IPTS_directory='/SNS/SEQ/IPTS-8138/'
for irun in runs:
	filename=IPTS_directory+'data/SEQ_'+str(irun)+'_event.nxs'
	w=Load(filename)
	
	wm=LoadNexusMonitors(filename)	
	alg=GetEi(wm)
	Ei=alg[0]
	erange=str(-0.2*Ei)+','+str(0.005*Ei)+','+str(0.95*Ei)
	#change log times
	for x in w.getRun().keys():
		if x not in ['duration','proton_charge','start_time','run_title','run_start','run_number','gd_prtn_chrg','end_time']:
			try:
				ShiftTime('w',x)
			except:
				pass
	#T0 chopper problem
	valC3=w.getRun()['Phase3'].getStatistics().median
	FilterByLogValue(InputWorkspace=w,OutputWorkspace=w,LogName='Phase3',MinimumValue=valC3-0.15,MaximumValue=valC3+0.15)

	# adjust time for pack C17 wired backward
	ChangeBinOffset(InputWorkspace=w,OutputWorkspace=w,Offset=500,IndexMin=54272,IndexMax=55295) 

	# FILTERS OUT LOW-POWER EVENTS BELOW 10% LEVEL
	FilterBadPulses(InputWorkspace=w,OutputWorkspace=w,LowerCutoff=10)				

#GroupingFile='/SNS/SEQ/IPTS-6225/shared/SEQ_4x1.xml'

	DgsReduction(SampleInputWorkspace='w',IncidentBeamNormalisation='ByCurrent',OutputWorkspace='reduced',
		DetectorVanadiumInputFile=IPTS_directory+'shared/autoreduce/van.nx5',UseProcessedDetVan=1,EnergyTransferRange=erange )

	r=mtd['reduced']
	psi=r.getRun()['CCR16Rot'].getStatistics().mean
	outfilename=IPTS_directory+'shared/Rotreduced/SEQ_'+str(irun)+'.nxspe'
	print irun,psi
	irun=irun+1
	SaveNXSPE(r,outfilename,r.getRun()['Ei'].value,psi,'1')

