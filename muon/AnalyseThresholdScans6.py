from __future__ import print_function
from mantid.kernel import *
from mantid.api import *
import math
 
class AnalyseThresholdScans(PythonAlgorithm):
 
    def PyInit(self):
        # Declare properties
	# first and last runs to be loaded as for PlotAsymmetryByLogValue
	self.declareProperty(FileProperty("FirstFile","",FileAction.Load, extensions = ["nxs"]),doc="First of a continuous sequence of runs measured at different threshold settings")
	self.declareProperty("NumberOfRuns",1000000,doc="Number of runs, default continues until sequence stops. All must have threshold in Comment line")
	# use as many files as entries in thresholds list
	#self.declareFileProperty("LastFile","",FileAction.Load, Exts = ["nxs"])
	# assume thresholds are set in Comment line "thresholds  123mV" and in increasing order
	# self.declareListProperty("Thresholds",[25.0,37.0,50.0,62.0,75.0,87.0,100.0,125.0,150.0],Description="Comma separated list of the thresholds used")
	# if TF set to 0 then assume ZF / high LF data
	self.declareProperty("TransverseFieldGuess",20.0,doc="Hint for fit, or 0 for ZF/LF non-relaxing data, -1 for integral counts")
	self.declareProperty("FitDataFrom",20.0,doc="Start time for fit, may need to avoid high count rate distortions. Always fits to end of histograms.")
	self.declareProperty("FixMuonLifetime",True,doc="Could allow muon lifetime to vary with uncalibrated TDC and high stats runs")
	self.declareProperty("FixBackgroundToZero",True,doc="Could fit background rate (random, light-leak) with high stats only")
	self.declareProperty("Grouped",False,doc="Group detectors before fitting")

    def category(self):
        return 'Muon'

    def PyExec(self):
        # Run the algorithm
        firstfn = self.getProperty('FirstFile').value
	i=firstfn.rindex('.')
	j=i-1
	while firstfn[j-1].isdigit():
		j=j-1
	firstnum=int(firstfn[j:i])
	k=j-1
	while firstfn[k-1].isalpha():
		k=k-1
	#print "old run number was ",n0
	#voltlist=self.getProperty('Thresholds')
		
#path="\\\\hifi\\data\\cycle_10_3\\"
#path2="\\\\hifi\\data\\"
#inst="hifi000"
#extn=".nxs"
	bguess=self.getProperty('TransverseFieldGuess').value
	nmax=self.getProperty('NumberOfRuns').value
	t1=self.getProperty('FitDataFrom').value
	fml=self.getProperty('FixMuonLifetime').value
	fbz=self.getProperty('FixBackgroundToZero').value
	agroup=self.getProperty('Grouped').value
	
# table workspace would be nice but there's no way to create it?!
# also no way to extract Comment line from a run?
# tb=mtd.getTableWorkspace("threshold_runs")

#runlist= [15487, 15488, 15489, 15490, 15491, 15492, 15493, 15494, 15495, 30261]
#voltlist=[25.0,  37.0,  50.0,  62.0,  75.0,  87.0,  100.0, 125.0, 150.0, 50.0]

# columns are Run No, threshold(mV)
# first: group all data and fit to get B? take guess provided by user
# then: group across runs and get phase (and freq) for each spectrum
# then: individual fits for bg, asym, deadtime, countrate only.

# dead time correction and background need high stats to be meaningful

# allow "muon lifetime" to be variable too if TDC timebase is slightly off

# load the workspaces and add them keeping individuals too
	voltlist=[]
	runhandles=[]
	runnames=[]
	vi=0
	oldvolts=-1
	while(vi<nmax and oldvolts<999999):
	  #runlist.append(vi+firstnum)
	  thisrunname=firstfn[k:j]+str(vi+firstnum).zfill(i-j)
	  thisrunpath=firstfn[:j]+str(vi+firstnum).zfill(i-j)+firstfn[i:]
	  
	  #self.progress(0.25*vi/numruns)
	  try:
	    print("calling LoadMuonNexus with ",thisrunpath," and ",thisrunname)
	    thisoneX=LoadMuonNexus(filename=thisrunpath,OutputWorkspace=thisrunname,AutoGroup=agroup)
	  except:
	    print("no more runs to load")
	    break
	  thisone=thisoneX[0]
	  print("Comment field: [",thisone.getComment(),"]")
	  comnt=thisone.getComment()
	  mv=comnt.find("mV")
	  try:
		volts=int(comnt[11:mv])
	  except:
		volts=-1.0
	  if(mv>0 and (volts>oldvolts or nmax<1000000) and comnt[0:11]=="thresholds "):
	    voltlist.append(volts)
	    runnames.append(thisrunname)
	    runhandles.append(thisone)
	    print("threshold incrementing from ",oldvolts," to ",volts)
	    oldvolts=volts
	  else:
	    oldvolts=999999
	    break
	  if(vi==0):
	    sumd=CloneWorkspace(InputWorkspace=runnames[vi],OutputWorkspace="Summed")
	    proc=ConvertToPointData(InputWorkspace=runnames[vi],OutputWorkspace="rates")
	    #proc=mtd.getMatrixWorkspace("rates")
	    procA=ConvertToPointData(InputWorkspace=runnames[vi],OutputWorkspace="Asyms")
	    #procA=mtd.getMatrixWorkspace("Asyms")
	    procB=ConvertToPointData(InputWorkspace=runnames[vi],OutputWorkspace="BGs")
	    #procB=mtd.getMatrixWorkspace("BGs")
	    procD=ConvertToPointData(InputWorkspace=runnames[vi],OutputWorkspace="DTs")
	    #procD=mtd.getMatrixWorkspace("DTs")
	    #sumd=mtd.getMatrixWorkspace("summed")
	  else:
	    sumd += thisone
	  vi=vi+1

	numruns=len(voltlist)
	print("found sequence of ",numruns," runs with thresholds between ",voltlist[0]," and ",voltlist[-1]," mV")

# fitting..

# loop through spectra
	nhtp=sumd.getNumberHistograms()
# nhtp=3 # testing
	for i in range(0,nhtp):
	  #self.progress(0.25+0.7*i/nhtp)
# fit the summed runs first (if TF)
	  if(bguess>0):
		  if(i>=sumd.getNumberHistograms()/2):
		    phguess=3.14
		  else:
		    phguess=6.28
		  nguess=sumd.dataY(i)[64]
		  #ff=mantid.createAlgorithm("Fit") # had mtd.createAlgorithm
		  #ff.setPropertyValue("Function","name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*f+p))*exp(-x/tau)+b)+d), n="+str(nguess)+", a=0.2, f="+str(bguess*0.01355*3.14159*2.0)+", p="+str(phguess)+", b=0.01, d=0.001, tau=2.19703")
		  #ff.setPropertyValue("InputWorkspace","summed")
		  #ff.setPropertyValue("StartX",str(t1))
		  #ff.setPropertyValue("EndX","32")
		  #ff.setPropertyValue("Output","fitres")
		  #print "Predicting f=",str(bguess*0.01355*3.14159*2.0)," and p=",str(phguess)
		  #ff.setPropertyValue("WorkspaceIndex",str(i))
		  if(fml):
			if(fbz):
				#ff.setPropertyValue("Ties","tau=2.19703,b=0.0")
				ties="tau=2.19703,b=0.0"
			else:
				#ff.setPropertyValue("Ties","tau=2.19703")
				ties="tau=2.19703"
		  else:
			if(fbz):
				#ff.setPropertyValue("Ties","b=0.0")
				ties="b=0.0"
			else:
				ties=""
		  #ff.execute()
		  ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*f+p))*exp(-x/tau)+b)+d), n="+str(nguess)+", a=0.2, f="+str(bguess*0.01355*3.14159*2.0)+", p="+str(phguess)+", b=0.01, d=0.001, tau=2.19703",InputWorkspace="Summed",StartX=str(t1),EndX=32,Output="fitres",Ties=ties,WorkspaceIndex=str(i))
		  fr=mtd["fitres_Parameters"]
		  afitted=fr.cell("Value",1)
		  bfitted=fr.cell("Value",2)
		  phfitted=fr.cell("Value",3)
		  taufitted=fr.cell("Value",4)
		  if(fr.cell("Value",1)<0):
		   phfitted=phfitted+3.14159265 # correct for phase jumping negative
		   afitted=-afitted
		  print("hist ",i," freq=",bfitted," ampl=",afitted," phase ",phfitted," tau=", taufitted)
		  if(fml):
			  taufitted=2.19703
		  for vi in range(0,len(voltlist)):
	    # now fit the individual runs' spectra
		    #nfst0=runhandles[vi].getSampleDetails().getLogData("ICPEVENT").value
		    #gfi=nfst0.rfind("GF")
		    #gfie=nfst0.rfind("RF")
		    #nframes=int(nfst0[gfi+2:gfie])
		    nframes=runhandles[vi].getSampleDetails().getLogData("goodfrm").value
		    tres=runhandles[vi].readX(i)[2]-runhandles[vi].readX(i)[1]  # microseconds
		    #ff=mtd.createAlgorithm("Fit")
		    nguess=runhandles[vi].dataY(i)[64]
		    #ff.setPropertyValue("Function","name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+"+str(phfitted)+"))*exp(-x/"+str(taufitted)+")+b)+d), n="+str(nguess)+", a=0.2, b=0.01, d=0.001")
		    #ff.setPropertyValue("InputWorkspace",runnames[vi])
		    #ff.setPropertyValue("WorkspaceIndex",str(i))
		    #ff.setPropertyValue("StartX",str(t1))
		    #ff.setPropertyValue("EndX","32")
		    #ff.setPropertyValue("Output","fitres")
		    if(fbz):
			#ff.setPropertyValue("Ties","b=0.0")
			ties="b=0.0"
		    else:
			    ties=""
		    #ff.execute()
		    ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+"+str(phfitted)+"))*exp(-x/"+str(taufitted)+")+b)+d), n="+str(nguess)+", a=0.2, b=0.01, d=0.001",InputWorkspace=runnames[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,Output="fitres",Ties=ties)
		    fr=mtd["fitres_Parameters"]
		    Nfitted=fr.cell("Value",0)
		    NfitE=fr.cell("Error",0)
		    Afitted=fr.cell("Value",1)
		    AfitE=fr.cell("Error",1)
		    Bfitted=fr.cell("Value",2)
		    BfitE=fr.cell("Error",2)
		    Dfitted=fr.cell("Value",3)
		    DfitE=fr.cell("Error",3)
	    #print "fitted individual run ",vi," to give asym=",Afitted," backgd=",Bfitted," and deadtime=",Dfitted
		    proc.dataX(i)[vi]=voltlist[vi]
		    proc.dataY(i)[vi]=Nfitted/nframes/tres # counts/us/frame
		    proc.dataE(i)[vi]=NfitE/nframes/tres
		    procA.dataX(i)[vi]=voltlist[vi]
		    procA.dataY(i)[vi]=Afitted
		    procA.dataE(i)[vi]=AfitE
		    procB.dataX(i)[vi]=voltlist[vi]
		    procB.dataY(i)[vi]=Bfitted/nframes/tres*1.0E6 # units of counts per second 
		    procB.dataE(i)[vi]=BfitE/nframes/tres*1.0E6
		    procD.dataX(i)[vi]=voltlist[vi]
		    procD.dataY(i)[vi]=Dfitted*nframes*tres # dead time in microseconds
		    procD.dataE(i)[vi]=DfitE*nframes*tres
	  #print "finished fitting spectrum ",i
	  elif (bguess==0):
		  for vi in range(0,len(voltlist)):
	    # now the individual runs' spectra
		    #nfst0=runhandles[vi].getSampleDetails().getLogData("ICPEVENT").value
		    #gfi=nfst0.rfind("GF")
		    #gfie=nfst0.rfind("RF")
		    #nframes=int(nfst0[gfi+2:gfie])
		    nframes=runhandles[vi].getSampleDetails().getLogData("goodfrm").value
		    tres=runhandles[vi].readX(i)[2]-runhandles[vi].readX(i)[1]  # microseconds
		    #ff=mtd.createAlgorithm("Fit")
		    nguess=runhandles[vi].dataY(i)[64]
		    #ff.setPropertyValue("Function","name=UserFunction, Formula=1/(1/(n*exp(-x/tau)+b)+d), n="+str(nguess)+", b=0.01, d=0.001, tau=2.19703")
		    #ff.setPropertyValue("InputWorkspace",runnames[vi])
		    #ff.setPropertyValue("WorkspaceIndex",str(i))
		    #ff.setPropertyValue("StartX",str(t1))
		    #ff.setPropertyValue("EndX","32")
		    #ff.setPropertyValue("Output","fitres")
		    if(fml):
			if(fbz):
				#ff.setPropertyValue("Ties","tau=2.19703,b=0.0")
				ties="tau=2.19703,b=0.0"
			else:
				#ff.setPropertyValue("Ties","tau=2.19703")
				ties="tau=2.19703"
		    else:
			if(fbz):
				#ff.setPropertyValue("Ties","b=0.0")
				ties="b=0.0"
			else:
				ties=""
		    #ff.execute()
		    ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*exp(-x/tau)+b)+d), n="+str(nguess)+", b=0.01, d=0.001, tau=2.19703",InputWorkspace=runnames[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,Output="fitres",Ties=ties)
		    fr=mtd["fitres_Parameters"]
		    Nfitted=fr.cell("Value",0)
		    NfitE=fr.cell("Error",0)
		    Bfitted=fr.cell("Value",2)
		    BfitE=fr.cell("Error",2)
		    Dfitted=fr.cell("Value",3)
		    DfitE=fr.cell("Error",3)
	    #print "fitted individual run ",vi," to give asym=",Afitted," backgd=",Bfitted," and deadtime=",Dfitted
		    proc.dataX(i)[vi]=voltlist[vi]
		    proc.dataY(i)[vi]=Nfitted/nframes/tres # counts/us/frame
		    proc.dataE(i)[vi]=NfitE/nframes/tres
		    procA.dataX(i)[vi]=voltlist[vi]
		    procA.dataY(i)[vi]=0.0 # no asym possible for LF/ZF runs
		    procA.dataE(i)[vi]=0.0
		    procB.dataX(i)[vi]=voltlist[vi]
		    procB.dataY(i)[vi]=Bfitted/nframes/tres*1.0E6 # units of counts per second 
		    procB.dataE(i)[vi]=BfitE/nframes/tres*1.0E6
		    procD.dataX(i)[vi]=voltlist[vi]
		    procD.dataY(i)[vi]=Dfitted*nframes*tres # dead time in microseconds
		    procD.dataE(i)[vi]=DfitE*nframes*tres
	  else:
		  # integral counts only
		  for vi in range(0,len(voltlist)):
		    #nfst0=runhandles[vi].getSampleDetails().getLogData("ICPEVENT").value
		    #gfi=nfst0.rfind("GF")
		    #gfie=nfst0.rfind("RF")
		    #nframes=int(nfst0[gfi+2:gfie])
		    nframes=runhandles[vi].getSampleDetails().getLogData("goodfrm").value
		    tres=runhandles[vi].readX(i)[2]-runhandles[vi].readX(i)[1]  # microseconds. Not needed for integral?
		    sum=0
		    hx=runhandles[vi].readX(i)
		    hy=runhandles[vi].readY(i)
		    for tt in range(len(hy)):
			    if(hx[tt+1]>t1):
				    sum=sum+hy[tt]
		    # integral from t1 only
		    # true total counts per microsecond at t=0 is sum*exp(t1/tau)/tau
		    proc.dataX(i)[vi]=voltlist[vi]
		    proc.dataY(i)[vi]=sum*math.exp(t1/2.19703)/2.19703/nframes # counts/us/frame
		    proc.dataE(i)[vi]=math.sqrt(sum)*math.exp(t1/2.19703)/2.19703/nframes
		    procA.dataX(i)[vi]=voltlist[vi]
		    procA.dataY(i)[vi]=0.0 # no asym possible for LF/ZF runs
		    procA.dataE(i)[vi]=0.0
		    procB.dataX(i)[vi]=voltlist[vi]
		    procB.dataY(i)[vi]=0.0 # BG not fitted 
		    procB.dataE(i)[vi]=0.0
		    procD.dataX(i)[vi]=voltlist[vi]
		    procD.dataY(i)[vi]=0.0 # deadtime not fitted
		    procD.dataE(i)[vi]=0.0
  #sleep(0.5)
	  for vi in range(len(voltlist),len(proc.dataX(i))):
	    proc.dataX(i)[vi]=1001.0+vi
	    procA.dataX(i)[vi]=1001.0+vi
	    procB.dataX(i)[vi]=1001.0+vi
	    procD.dataX(i)[vi]=1001.0+vi


# tidy up workspace "proc"
	CropWorkspace(InputWorkspace="rates",OutputWorkspace="rates",XMin=0.0,XMax=998.0)
	CropWorkspace(InputWorkspace="Asyms",OutputWorkspace="Asyms",XMin=0.0,XMax=998.0)
	CropWorkspace(InputWorkspace="BGs",OutputWorkspace="BGs",XMin=0.0,XMax=998.0)
	CropWorkspace(InputWorkspace="DTs",OutputWorkspace="DTs",XMin=0.0,XMax=998.0)
# workspace "results" has 3 "time" bins, rate, mean pulse height and width
	CloneWorkspace(InputWorkspace="rates",OutputWorkspace="results")
	resw=mtd["results"]
# fit N versus volts to erf() or something similar
# omit <25mV point2 from fit as it may have noise as well as positrons
	for i in range(0,nhtp):
	  #self.progress(0.95+0.05*i/proc.getNumberHistograms())
	  #ff=mtd.createAlgorithm("Fit")
	  #ff.setPropertyValue("Function","name=UserFunction, Formula=rate/(1+exp(x/lam - 2.75)), rate=5.0, lam=20.0")
	  #ff.setPropertyValue("InputWorkspace","rates")
	  #ff.setPropertyValue("StartX","24.0")
	  #ff.setPropertyValue("EndX","999.0")
	  #ff.setPropertyValue("Output","fitres")
  #ff.setPropertyValue("Ties","PeakCentre=0.0")
	  #ff.setPropertyValue("WorkspaceIndex",str(i))
	  #ff.execute()
	  ff=Fit(Function="name=UserFunction, Formula=rate/(1+exp(x/lam - 2.75)), rate=5.0, lam=20.0",InputWorkspace="rates",StartX="24.0",EndX="999.0",Output="fitres",WorkspaceIndex=str(i))
	  fr=mtd["fitres_Parameters"]
	  for j in range(0,3):
	    resw.dataX(i)[j]=j*1.0
	    resw.dataY(i)[j]=fr.cell("Value",j)
	    try:
	      resw.dataE(i)[j]=fr.cell("Error",j)
	    except:
	      print("missing Error column from fit")
	      resw.dataE(i)[j]=0.0
	  resw.dataX(i)[3]=5.0
	  #print "fitted pulse height spectrum ",i," with mean=",fr.getDouble("Value",1)
  #time.sleep(0.5)

	CropWorkspace(InputWorkspace="results",OutputWorkspace="results",XMin=0.0,XMax=3.0)
	#self.progress(1.0)
#print "last hist had nframes=",nframes," and time res=",tres
	#print "Finished!"
# open Instrument View?

AlgorithmFactory.subscribe(AnalyseThresholdScans)
