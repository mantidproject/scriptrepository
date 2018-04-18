from mantid.kernel import *
from mantid.api import *
import math
import numpy
import re

class AnalyseThresholdScans(PythonAlgorithm):

	def PyInit(self):
		# Declare properties
		# first and last runs to be loaded as for PlotAsymmetryByLogValue
		self.declareProperty(FileProperty("FirstFile","",FileAction.Load, extensions = ["nxs"]),doc="First of a continuous sequence of runs measured at different threshold settings")
		self.declareProperty("NumberOfRuns",1000000,doc="Number of runs, default continues until sequence stops. All must have threshold in Comment line")
		self.declareProperty(FloatArrayProperty("ThresholdList",direction=Direction.Input),doc="Threshold values if not given in run comments, overrides NumberOfRuns")
		# use as many files as entries in thresholds list
		#self.declareFileProperty("LastFile","",FileAction.Load, Exts = ["nxs"])
		# assume thresholds are set in Comment line "thresholds  123mV" and in increasing order
		# self.declareListProperty("Thresholds",[25.0,37.0,50.0,62.0,75.0,87.0,100.0,125.0,150.0],Description="Comma separated list of the thresholds used")
		# if TF set to 0 then assume ZF / high LF data
		self.declareProperty("TransverseFieldGuess",20.0,doc="Hint for fit, or 0 for ZF/LF non-relaxing data, -1 for integral counts")
		self.declareProperty("FitDataFrom",0.2,doc="Start time for fit, may need to avoid high count rate distortions. Always fits to end of histograms.")
		self.declareProperty("FixMuonLifetime",True,doc="Could allow muon lifetime to vary with uncalibrated TDC and high stats runs")
		self.declareProperty("FixBackgroundToZero",True,doc="Could fit background rate (random, light-leak) with high stats only")
		self.declareProperty("IndividualPhase",False,doc="Allow phase to vary with threshold as well as between detectors")
		self.declareProperty("Grouped",False,doc="Group detectors before fitting")
		self.declareProperty(WorkspaceProperty("Rates","",Direction.Output),doc="Rates as function of threshold")
		self.declareProperty(ITableWorkspaceProperty("ResultTable","",Direction.Output),doc="Positron rates and pulse heights")
		self.declareProperty(WorkspaceProperty("Asymmetries","",Direction.Output,PropertyMode.Optional),doc="Asymmetry as function of threshold (TF only)")
		self.declareProperty(WorkspaceProperty("Phases","",Direction.Output,PropertyMode.Optional),doc="Phase as function of threshold (TF only)")
		self.declareProperty(WorkspaceProperty("DeadTimes","",Direction.Output,PropertyMode.Optional),doc="Apparent dead time as function of threshold")
		self.declareProperty(WorkspaceProperty("Backgrounds","",Direction.Output,PropertyMode.Optional),doc="Background counts as function of threshold")
		self.declareProperty(FileProperty("IDF","",FileAction.OptionalLoad, extensions = ["xml"]),doc="IDF for the result files, for Instrument View")


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

		bguess=self.getProperty('TransverseFieldGuess').value
		nmax=self.getProperty('NumberOfRuns').value
		t1=self.getProperty('FitDataFrom').value
		fml=self.getProperty('FixMuonLifetime').value
		fbz=self.getProperty('FixBackgroundToZero').value
		agroup=self.getProperty('Grouped').value
		indPhases=self.getProperty("IndividualPhase").value

		thresholdsGiven=self.getProperty("ThresholdList").value
		if(len(thresholdsGiven)>0):
			nmax=len(thresholdsGiven)
	
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
		thUnit="mV" # default, can override
		while(vi<nmax and oldvolts<999999):
			#runlist.append(vi+firstnum)
			thisrunname=firstfn[k:j]+str(vi+firstnum).zfill(i-j)
			thisrunpath=firstfn[:j]+str(vi+firstnum).zfill(i-j)+firstfn[i:]
	  
			try:
				self.log().information("calling LoadMuonNexus with "+thisrunpath+" and "+thisrunname)
				thisoneX=LoadMuonNexus(filename=thisrunpath,OutputWorkspace=thisrunname,AutoGroup=agroup)
			except:
				self.log().information("no more runs to load")
				break
			thisone=thisoneX[0]
			self.log().information("Comment field: ["+thisone.getComment()+"]")
			if(len(thresholdsGiven)>0):
				mo=True
				volts=thresholdsGiven[vi]
			else:
				comnt=thisone.getComment()
				pat=r"^(?P<name>.*?)(?P<num>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s*(?P<unit>\S*)"
				mo=re.search(pat,comnt)
				if(mo):
					if(mo.group("name")!="" and "thr" not in mo.group("name").lower()):
						mo=None
				if(mo):
					volts=float(mo.group("num"))
					if(mo.group("unit") and vi==0):
						thUnit=mo.group("unit") # select on 1st only
				else:
					volts=-1.0
			if(mo and (volts>oldvolts or nmax<1000000) ):
				voltlist.append(volts)
				runnames.append(thisrunname)
				runhandles.append(thisone)
				self.log().information("threshold incrementing from "+str(oldvolts)+" to "+str(volts))
				oldvolts=volts
			else:
				oldvolts=999999
				break
			if(vi==0):
				sumd=CloneWorkspace(InputWorkspace=runnames[vi],OutputWorkspace="Summed")
			else:
				sumd += thisone
			vi=vi+1

		numruns=len(voltlist)
		self.log().notice("found sequence of "+str(numruns)+" runs with thresholds between "+str(voltlist[0])+" and "+str(voltlist[-1])+thUnit)

		# fitting..

		# loop through spectra
		nhtp=sumd.getNumberHistograms()
		proc=WorkspaceFactory.create("Workspace2D",NVectors=nhtp,XLength=numruns,YLength=numruns) # rates
		procA=WorkspaceFactory.create("Workspace2D",NVectors=nhtp,XLength=numruns,YLength=numruns) # Asyms
		procB=WorkspaceFactory.create("Workspace2D",NVectors=nhtp,XLength=numruns,YLength=numruns) # BGs
		procD=WorkspaceFactory.create("Workspace2D",NVectors=nhtp,XLength=numruns,YLength=numruns) # DTs
		procP=WorkspaceFactory.create("Workspace2D",NVectors=nhtp,XLength=numruns,YLength=numruns) # phases, if used
		# output tables and workspaces - give them names so that MaskDetectors, etc will work
		self.setProperty("Rates",proc)
		if(not self.getProperty("Asymmetries").isDefault):
			self.setProperty("Asymmetries",procA)
		if(not self.getProperty("Backgrounds").isDefault):
			self.setProperty("Backgrounds",procB)
		if(not self.getProperty("DeadTimes").isDefault):
			self.setProperty("DeadTimes",procD)
		if(not self.getProperty("Phases").isDefault):
			self.setProperty("Phases",procP)
		# X axes are thresholds in mv
		for w in (proc,procA,procB,procD,procP):
			na=NumericAxis.create(numruns)
			for i,v in enumerate(voltlist):
				na.setValue(i,v)
				lbl=na.setUnit("Label")
				lbl.setLabel("Threshold",thUnit)
			w.replaceAxis(0,na)
		proc.setYUnitLabel("Counts per frame")
		procA.setYUnitLabel("Asymmetry")
		procB.setYUnitLabel("Background counts/frame")
		procD.setYUnitLabel("Dead time (us)")
		procP.setYUnitLabel("Phase (rad)")
		if(not self.getProperty("IDF").isDefault):
			# CHRONUS, etc: load alternative IDF. Do now since running this afterwards clears the mask flags!
			for pp in [proc,procA,procB,procD,procP]:
				iloader=AlgorithmManager.createUnmanaged("LoadInstrument")
				iloader.setChild(True) # doesn't need a name
				iloader.initialize()
				iloader.setProperty("Workspace", pp)
				iloader.setProperty("RewriteSpectraMap", True)
				iloader.setProperty("Filename", self.getProperty("IDF").value)
				iloader.execute()
			# LoadInstrument(Workspace=proc,RewriteSpectraMap=True,Filename=self.getProperty("IDF").value)

		prog=Progress(self,start=0.0,end=1.0,nreports=nhtp)

		for i in range(0,nhtp):
			# fit the summed runs first (if TF)
			# check for dead detector and mask it
			if(numpy.sum(sumd.dataY(i))<1):
				for pp in [proc,procA,procB,procD,procP]:
					masker = AlgorithmManager.createUnmanaged("MaskDetectors")
					masker.setChild(True) # doesn't need a name
					masker.initialize()
					masker.setProperty("Workspace", pp) 
					masker.setProperty("WorkspaceIndexList", [i]) # or the appropriate index
					masker.execute()

			else:
				if(bguess>0):
					if(i>=sumd.getNumberHistograms()/2):
						phguess=3.14
					else:
						phguess=6.28
					nguess=numpy.amax(sumd.dataY(i))
					#ff=mantid.createAlgorithm("Fit") # had mtd.createAlgorithm
					#ff.setPropertyValue("Function","name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*f+p))*exp(-x/tau)+b)+d), n="+str(nguess)+", a=0.2, f="+str(bguess*0.01355*3.14159*2.0)+", p="+str(phguess)+", b=0.01, d=0.001, tau=2.19703")
					#ff.setPropertyValue("InputWorkspace","summed")
					#ff.setPropertyValue("StartX",str(t1))
					#ff.setPropertyValue("EndX","32")
					#ff.setPropertyValue("Output","fitres")
					#print "Predicting f=",str(bguess*0.01355*3.14159*2.0)," and p=",str(phguess)
					#ff.setPropertyValue("WorkspaceIndex",str(i))
					ties=[]
					if(fml):
						ties.append("tau=2.19703")
					if(fbz):
						ties.append("b=0.0")
					#ff.execute()
					ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*f+p))*exp(-x/tau)+b)+d), n="+str(nguess)+", a=0.2, f="+str(bguess*0.01355*3.14159*2.0)+", p="+str(phguess)+", b=0.01, d=0.001, tau=2.19703",
						InputWorkspace="Summed",StartX=str(t1),EndX=32,Output="fitres",Ties=",".join(ties),WorkspaceIndex=str(i))
					fr=ff[3] # mtd["fitres_Parameters"]
					afitted=fr.cell("Value",1)
					bfitted=fr.cell("Value",2)
					phfitted=fr.cell("Value",3)
					taufitted=fr.cell("Value",4)
					if(fr.cell("Value",1)<0):
						phfitted=phfitted+3.14159265 # correct for phase jumping negative
						afitted=-afitted
					self.log().information("hist "+str(i)+" freq="+str(bfitted)+" ampl="+str(afitted)+" phase "+str(phfitted)+" tau="+str(taufitted))
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
						nguess=numpy.amax(runhandles[vi].dataY(i))
						#ff.setPropertyValue("Function","name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+"+str(phfitted)+"))*exp(-x/"+str(taufitted)+")+b)+d), n="+str(nguess)+", a=0.2, b=0.01, d=0.001")
						#ff.setPropertyValue("InputWorkspace",runnames[vi])
						#ff.setPropertyValue("WorkspaceIndex",str(i))
						#ff.setPropertyValue("StartX",str(t1))
						#ff.setPropertyValue("EndX","32")
						#ff.setPropertyValue("Output","fitres")
						ties=[]
						if(fbz):
							ties.append("b=0.0")
						if(not indPhases):
							ties.append("p="+str(phfitted))
						#ff.execute()
						ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+p))*exp(-x/"+str(taufitted)+")+b)+d), n="+str(nguess)+", a=0.2, b=0.01, d=0.001, p="+str(phfitted),
							InputWorkspace=runnames[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,Output="fitres",Ties=",".join(ties))
						fr=ff[3] # mtd["fitres_Parameters"]
						Nfitted=fr.cell("Value",0)
						NfitE=fr.cell("Error",0)
						Afitted=fr.cell("Value",1)
						AfitE=fr.cell("Error",1)
						Pfitted=fr.cell("Value",2)
						PfitE=fr.cell("Error",2)
						Bfitted=fr.cell("Value",3)
						BfitE=fr.cell("Error",3)
						Dfitted=fr.cell("Value",4)
						DfitE=fr.cell("Error",4)
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
						procP.dataX(i)[vi]=voltlist[vi]
						procP.dataY(i)[vi]=Pfitted
						procP.dataE(i)[vi]=PfitE
						#print "finished fitting spectrum ",i
				elif (bguess==0):
					for vi in range(0,len(voltlist)):
						nframes=runhandles[vi].getSampleDetails().getLogData("goodfrm").value
						tres=runhandles[vi].readX(i)[2]-runhandles[vi].readX(i)[1]  # microseconds
						#ff=mtd.createAlgorithm("Fit")
						nguess=numpy.amax(runhandles[vi].dataY(i))
						#ff.setPropertyValue("Function","name=UserFunction, Formula=1/(1/(n*exp(-x/tau)+b)+d), n="+str(nguess)+", b=0.01, d=0.001, tau=2.19703")
						#ff.setPropertyValue("InputWorkspace",runnames[vi])
						#ff.setPropertyValue("WorkspaceIndex",str(i))
						#ff.setPropertyValue("StartX",str(t1))
						#ff.setPropertyValue("EndX","32")
						#ff.setPropertyValue("Output","fitres")
						ties=[]
						if(fml):
							ties.append("tau=2.19703")
						if(fbz):
							ties.append("b=0.0")
						#ff.execute()
						ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*exp(-x/tau)+b)+d), n="+str(nguess)+", b=0.01, d=0.001, tau=2.19703",
							InputWorkspace=runnames[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,Output="fitres",Ties=",".join(ties))
						fr=ff[3] # mtd["fitres_Parameters"]
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
			prog.report()

		resw=WorkspaceFactory.createTable()
		resw.addColumn("int","Spectrum",1)
		resw.addColumn("double","CountsPerFrame",2)
		resw.addColumn("double","ErrCPF",5)
		resw.addColumn("double","PulseHeight",2)
		resw.addColumn("double","ErrPH",5)
		resw.addColumn("double","ChiSq",2)
		self.setProperty("ResultTable",resw)
		# fit N versus volts to erf() or something similar
		# omit <25mV point2 from fit as it may have noise as well as positrons
		# though if voltlist range is all <25 then assume relative (CHRONUS) and fit everything
		if(voltlist[-1]>50):
			x0=24.0
		else:
			x0=0.0
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
			ff=Fit(Function="name=UserFunction, Formula=rate/(1+exp(x/lam - 2.75)), rate=5.0, lam=20.0",InputWorkspace=proc,StartX=x0,EndX="999.0",Output="fitres",WorkspaceIndex=str(i))
			fr=ff[3] # mtd["fitres_Parameters"]
			resw.addRow((i,fr.cell("Value",0),fr.cell("Error",0),fr.cell("Value",1),fr.cell("Error",1),fr.cell("Value",2)))

		# output workspaces already set...

AlgorithmFactory.subscribe(AnalyseThresholdScans)
