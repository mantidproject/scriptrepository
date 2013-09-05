"""
Load a white beam data, convert to MD, find peaks, remove those belonging to Al, and index the rest
"""
from numpy import *
from string import *

def removeAlpeaks(wsname,sigma=0.1,dmin=0,dmax=1e10):
	#good values for dmin,dmax are 1 and 5
	AlD=4.04/array([sqrt(3.),2.,sqrt(8.)])
	peakstoremove=[]
	w=mtd[wsname]
	for i in range(w.getNumberPeaks()):
		p=w.getPeak(i)
		if (p.getDetectorID()<0) or (p.getDSpacing()<dmin) or (p.getDSpacing()>dmax) :
			peakstoremove.append(i)
		else:
			d=abs(AlD-p.getDSpacing()).min()
			if d<sigma:
				peakstoremove.append(i)
	DeleteTableRows(TableWorkspace=wsname,Rows=str(peakstoremove).replace(']','').replace('[',''))
	
	
Load(Filename='/SNS/SEQ/IPTS-6047/data/SEQ_17932_event.nxs',OutputWorkspace='SEQ_17932_event')
SetGoniometer(Workspace='SEQ_17932_event',Axis0='CCR13Vrot,0,1,0,1')
ConvertToDiffractionMDWorkspace(InputWorkspace='SEQ_17932_event',OutputWorkspace='d_Md',OutputDimensions='Q (sample frame)')
FindPeaksMD(InputWorkspace='d_Md',MaxPeaks='50',DensityThresholdFactor='1000000',OutputWorkspace='peaks')
removeAlpeaks('peaks',dmin=1,dmax=5)
FindUBUsingLatticeParameters('peaks',8.455,8.455,8.455,90,90,90,Tolerance=0.2)
IndexPeaks('peaks')
