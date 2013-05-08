"""
Temperature dependence of scattering intensity between Qmin,Qmax. Detectors are groupped 16x4, the groups belonging to data/background are selected below
"""

from numpy import *
from string import *

IPTSpath='/SNS/SEQ/IPTS-6225'
runsup=range(32298,32357)
runsdown=arange(32357,32411)
runsdown=append(runsdown,arange(32412,32513))

Qmin=0.7
Qmax=1.5
#backgroundGroups=[921,922,925,926,929,930,933,934]
backgroundGroups=[921,926,929,934,937,942,945,950,953,958]
dataGroups=[923,924,931,932]

alltemps=[]
allsignal=[]
allerror=[]
allbkg=[]
allbkgerr=[]

for r in runsup:
	filename=IPTSpath+'/data/SEQ_'+str(r)+'_event.nxs'
	w=Load(filename)																	# Mantid algorithm to load data
	w=ChangeBinOffset(w,500,IndexMin=54272,IndexMax=55295) 							# Mantid algorithm to fix timing in miswired det. bank C17 (add 500 microsec) 
	w=NormaliseByCurrent(w)															# Mantid algorithm to normalize data by current
	w=ConvertUnits(w,'MomentumTransfer','Elastic')									# Mantid algorithm to convert time-of-flight to Q assuming scattering is in elastic channel
	w=Rebin(w,str(Qmin)+','+str(Qmax-Qmin)+','+str(Qmax),PreserveEvents=0)			# Mantid algorithm to rebin into a single bin of length 
	wg=GroupDetectors(w,'/SNS/SEQ/IPTS-6225/shared/SEQ_16x4.xml',Behaviour='Sum')		# Mantid algorithm to group detector pixels using SEQ_16x4.xml mask, Behaviour='Sum' overrides the default ('Average')
	bkg=0.
	signal=0.
	for i in backgroundGroups:
		bkg+=wg.readY(i)[0]
	for i in dataGroups:
		signal+=wg.readY(i)[0]
	errbkg=sqrt(bkg)/len(backgroundGroups)
	bkg=bkg/len(backgroundGroups)
	err=sqrt(signal)/len(dataGroups)
	signal=signal/len(dataGroups)
	temp=(wg.getRun()['SampleTemp']).getStatistics().mean
	print temp,signal,err,bkg,errbkg
	alltemps.append(temp)
	allsignal.append(signal)
	allerror.append(err)
	allbkg.append(bkg)
	allbkgerr.append(errbkg)
	
z=zip(alltemps,allsignal,allerror,allbkg,allbkgerr)
z.sort()
alltemps,allsignal,allerror,allbkg,allbkgerr=zip(*z)

Intup=CreateWorkspace(alltemps,allsignal,allerror,1,YUnitLabel='Counts/muAh')
BGup=CreateWorkspace(alltemps,allbkg,allbkgerr,1,YUnitLabel='Counts/muAh')
Inetup=Intup-BGup

alltemps=[]
allsignal=[]
allerror=[]
allbkg=[]
allbkgerr=[]

for r in runsdown:
	filename=IPTSpath+'/data/SEQ_'+str(r)+'_event.nxs'
	w=Load(filename)																	# Mantid algorithm to load data
	w=ChangeBinOffset(w,500,IndexMin=54272,IndexMax=55295) 							# Mantid algorithm to fix timing in miswired det. bank C17 (add 500 microsec) 
	w=NormaliseByCurrent(w)															# Mantid algorithm to normalize data by current
	w=ConvertUnits(w,'MomentumTransfer','Elastic')									# Mantid algorithm to convert time-of-flight to Q assuming scattering is in elastic channel
	w=Rebin(w,str(Qmin)+','+str(Qmax-Qmin)+','+str(Qmax),PreserveEvents=0)			# Mantid algorithm to rebin into a single bin of length 
	wg=GroupDetectors(w,'/SNS/SEQ/IPTS-6225/shared/SEQ_16x4.xml',Behaviour='Sum')		# Mantid algorithm to group detector pixels using SEQ_16x4.xml mask, Behaviour='Sum' overrides the default ('Average')
	bkg=0.
	signal=0.
	for i in backgroundGroups:
		bkg+=wg.readY(i)[0]
	for i in dataGroups:
		signal+=wg.readY(i)[0]
	errbkg=sqrt(bkg)/len(backgroundGroups)
	bkg=bkg/len(backgroundGroups)
	err=sqrt(signal)/len(dataGroups)
	signal=signal/len(dataGroups)
	temp=(wg.getRun()['SampleTemp']).getStatistics().mean
	print temp,signal,err,bkg,errbkg
	alltemps.append(temp)
	allsignal.append(signal)
	allerror.append(err)
	allbkg.append(bkg)
	allbkgerr.append(errbkg)
	
z=zip(alltemps,allsignal,allerror,allbkg,allbkgerr)
z.sort()
alltemps,allsignal,allerror,allbkg,allbkgerr=zip(*z)

Intdown=CreateWorkspace(alltemps,allsignal,allerror,1,YUnitLabel='Counts/muAh')
BGdown=CreateWorkspace(alltemps,allbkg,allbkgerr,1,YUnitLabel='Counts/muAh')
Inetdown=Intdown-BGdown


pl=plot((Intup,Intdown,BGup,BGdown),0)
pld=plot((Inetup,Inetdown),0)

