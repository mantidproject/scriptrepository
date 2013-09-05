"""
Plot temperature dependence for data in a set of runs
"""
from numpy import *
from string import *

IPTSpath='/SNS/SEQ/IPTS-6225'
runsTup=range(32752,32846)

#filename=IPTSpath+'/data/SEQ_32846_event.nxs'
#wm=LoadNexusMonitors(filename)
#alg=GetEi(wm)
#Efixed=alg[0]
#t0=-alg[3]
Efixed=49.83
t0=-29.64

roi1=LoadMask('/SNS/SEQ/IPTS-6225/shared/Tdep50meV/roi_1p5_0p0.xml','SEQUOIA')
roi2=LoadMask('/SNS/SEQ/IPTS-6225/shared/Tdep50meV/roi_1p5_0p5.xml','SEQUOIA')

Emin=-1
Emax=1

allTemps=[]
allIntegratedSignal1=[]
allPeakSignal1=[]
allIntegratedSignal2=[]
allPeakSignal2=[]

for r in runsTup:
	filename=IPTSpath+'/data/SEQ_'+str(r)+'_event.nxs'
	w1=Load(filename)								# Mantid algorithm to load data

	w1=ChangeBinOffset(w1,t0)
	w1=ChangeBinOffset(w1,500,IndexMin=54272,IndexMax=55295)			# Mantid algorithm to fix timing in miswired det. bank C17 (add 500 microsec) 
	w1=NormaliseByCurrent(w1)							# Mantid algorithm to normalize data by current
	w1=ConvertUnits(w1,'DeltaE','Direct',Efixed)
	w1=Rebin(w1,str(Emin)+','+str(Emax-Emin)+','+str(Emax),PreserveEvents=0)	# Mantid algorithm to rebin into a single bin of length 
	w2=CloneWorkspace(w1)

	MaskDetectors(w1,MaskedWorkspace=roi1)
	s1=SumSpectra(w1)
	Imax1=w1.extractY().max()
	signal1=s1.readY(0)[0]
	MaskDetectors(w2,MaskedWorkspace=roi2)
	s2=SumSpectra(w2)
	Imax2=w2.extractY().max()
	signal2=s2.readY(0)[0]

	temp=(w1.getRun()['SampleTemp']).getStatistics().mean
	print temp,signal1,Imax1,signal2,Imax2
	allTemps.append(temp)
	allIntegratedSignal1.append(signal1)
	allPeakSignal1.append(Imax1)
	allIntegratedSignal2.append(signal2)
	allPeakSignal2.append(Imax2)
	
z=zip(allTemps,allPeakSignal1,allIntegratedSignal1,allPeakSignal2,allIntegratedSignal2)
z.sort()
allTemps,allPeakSignal1,allIntegratedSignal1,allPeakSignal2,allIntegratedSignal2=zip(*z)

Int_1p5_0p0=CreateWorkspace(allTemps,allIntegratedSignal1,sqrt(allIntegratedSignal1),1,YUnitLabel='Counts/muAh')
Peak_1p5_0p0=CreateWorkspace(allTemps,allPeakSignal1,sqrt(allPeakSignal1),1,YUnitLabel='Counts/muAh')

Int_1p5_0p5=CreateWorkspace(allTemps,allIntegratedSignal2,sqrt(allIntegratedSignal2),1,YUnitLabel='Counts/muAh')
Peak_1p5_0p5=CreateWorkspace(allTemps,allPeakSignal2,sqrt(allPeakSignal2),1,YUnitLabel='Counts/muAh')

pl=plot((Int_1p5_0p0,Int_1p5_0p5),0)
pld=plot((Peak_1p5_0p0,Peak_1p5_0p5),0)
