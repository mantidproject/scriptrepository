# takes wkspin as string of run number 
# returns ei, as value and rebin parameners as string
import numpy as np
from peakdet import *
from mantid.simpleapi import *


def autoEi(WkspIn,BkgdGen=None,monspecin=None,BkgdLevel=None):
	#monitor spectrum in spectrum number
	#pars for MARI
	
	monspec=monspecin
	
	
	if mtd.doesExist(WkspIn)==True:
		#data is passed as an existing workspace
		tmp_Monitors=ExtractSingleSpectrum(WkspIn,monspec-1)
	else:	
		#load single spectrum for monitor
		Load(Filename=WkspIn,OutputWorkspace='tmp_Monitors',Cache=r'Never',LoadLogFiles='0',SpectrumMin=monspec,SpectrumMax=monspec)
		#get some distances for later
		tmp_Monitors=mtd['tmp_Monitors']
	
	
	instrumentObj=tmp_Monitors.getInstrument()
	sampleObj=instrumentObj.getSample()

	sourceObj=instrumentObj.getSource()

	detObj=instrumentObj.getDetector(4101)
	L1=sourceObj.getDistance(sampleObj)
	L2=sampleObj.getDistance(detObj)
	
	ConvertToDistribution(Workspace='tmp_Monitors')
	#convert to energy & rebin to make life simple later the rebin is required for mari as monitor 2 always has a sit load of gammas at low tofs
	#these need to be ignored rebin can be removed for monitor 3 but the peak find is less reliable as there is always more noise in the m3
	ConvertUnits(InputWorkspace='tmp_Monitors',OutputWorkspace='tmp_Monitors',Target='Energy')
	Rebin(InputWorkspace='tmp_Monitors',OutputWorkspace='tmp_Monitors',Params='1,1,2000',PreserveEvents='0')
	
	#extract x and y data from specrum & delete mantid wksp goto point data in x
	dat=mtd['tmp_Monitors']
	x=dat.extractX()		
	y=dat.extractY()
	DeleteWorkspace('tmp_Monitors')
	y=y[0]
	x=x[0]

	xx=(x[1:len(x)]+x[0:len(x)-1])/2

	[maximal,miniaml]=peakdet(y, y[0]/2 ,xx)
	#find the biggest peak in the returned list of max data
	max1=maximal[:,1].argmax()
	#get the correspoiding energy
	ei=maximal[max1,0]
	print ei
	#generate some rebin parameters
	#deltaE is set at 10 points per resolution function of 2%Ei this is to make sofqw3 work. 
	#There is an issue with the overlapping polygon rebin that results in bad output if the energy transfer vector is 
	# is not over sampled
	deltaEi=(ei*.02)/10.0
	rebin=str(-ei*.5)+','+str(deltaEi)+','+str(ei*.95)
	#print rebin
	#print type(rebin)
	
	if BkgdGen == True:
		if BkgdLevel > 0:
			emin=BkgdLevel*ei   #minimum energy is with 70% energy loss
		else:
			emin=0.2*ei
		
		lam=(81.81/ei)**0.5
		lam_max=(81.81/emin)**0.5
		tsam=252.82*lam*L1   #time at sample
		#tmon2=252.82*lam*23.5 #time to monitor 6 on LET
		tmax=tsam+(252.82*lam_max*L2) #maximum time to measure inelastic signal to
		BkgdRange=[tmax,19500]
		return ei,rebin,BkgdRange
	else:
		return ei,rebin



#ei,rebin=autoEi('18322')
