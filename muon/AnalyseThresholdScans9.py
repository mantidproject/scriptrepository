from mantid.kernel import *
from mantid.api import *
from mantid.simpleapi import *
import math
import numpy
import re
import time

class AnalyseThresholdScans(PythonAlgorithm):

	timeSetup=0.0
	callsSetup=0
	timeExec=0.0
	callsExec=0
	timeSetPar=0.0
	callsSetPar=0

	def doChildAlg(self,AlgName,**InputProps):
		if isinstance(AlgName,str):
			#t1=time.clock()
			myalg=AlgorithmManager.create(AlgName) # brand new disposable one
			#t2=time.clock()
			if(AlgName=="Fit"):
			#	self.timeSetup += (t2-t1)
				self.callsSetup+=1
		else:
			myalg=AlgName # given recyclable algorithm
		myalg.setChild(True)
		priorities=myalg.orderedProperties()
		#t1=time.clock()
		for pri in priorities:
			if(pri in InputProps):
				myalg.setProperty(pri,InputProps[pri])
				del InputProps[pri]
		myalg.setProperties(InputProps) # rest in any order
		#t2=time.clock()
		if(AlgName=="Fit"):
			#self.timeSetPar += (t2-t1)
			self.callsSetPar+=1
		#t1=time.clock()
		myalg.execute()
		#t2=time.clock()
		if(AlgName=="Fit"):
			#self.timeExec += (t2-t1)
			self.callsExec+=1
		myres={}
		for op in myalg.outputProperties():
			myres[op]=myalg.getProperty(op).value
		return myres

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
		self.declareProperty("TransverseFieldGuess",20.0,doc="Hint for fit, or 0 for ZF/LF non-relaxing data, -1 for integral counts, -2 for rate and chisq")
		self.declareProperty("LambdaGuess",0.0,doc="positive to hint at TF lambda, negative to set fixed value")
		self.declareProperty("FitDataFrom",0.2,doc="Start time for fit, may need to avoid high count rate distortions. Always fits to end of histograms.")
		self.declareProperty("FixMuonLifetime",True,doc="Could allow muon lifetime to vary with uncalibrated TDC and high stats runs")
		self.declareProperty("FixBackgroundToZero",True,doc="Could fit background rate (random, light-leak) with high stats only")
		self.declareProperty("IndividualPhase",False,doc="Allow phase to vary with threshold as well as between detectors")
		self.declareProperty("Grouped",False,doc="Group detectors before fitting")
		self.declareProperty(ITableWorkspaceProperty("GroupingTable","",Direction.Input,PropertyMode.Optional),doc="Alternative grouping table e.g. by ring")
		self.declareProperty("BackgroundBefore",999.0,"Time before anything arrives, for background measurement on CHRONUS")
		self.declareProperty("PositronPeakTime",0.0,"Time at which a positron peak may be expected")
		self.declareProperty("MinThresholdForSum",0.0,"Exclude thresholds below this for initial summing and fitting TF")
		self.declareProperty(WorkspaceProperty("Rates","",Direction.Output),doc="Rates as function of threshold")
		self.declareProperty(ITableWorkspaceProperty("ResultTable","",Direction.Output),doc="Positron rates and pulse heights")
		self.declareProperty(WorkspaceProperty("Asymmetries","",Direction.Output,PropertyMode.Optional),doc="Asymmetry as function of threshold (TF only)")
		self.declareProperty(WorkspaceProperty("Phases","",Direction.Output,PropertyMode.Optional),doc="Phase as function of threshold (TF only)")
		self.declareProperty(WorkspaceProperty("DeadTimes","",Direction.Output,PropertyMode.Optional),doc="Apparent dead time as function of threshold")
		self.declareProperty(WorkspaceProperty("ChiSquares","",Direction.Output,PropertyMode.Optional),doc="Chi-squared of fit as function of threshold")
		self.declareProperty(WorkspaceProperty("Backgrounds","",Direction.Output,PropertyMode.Optional),doc="Background counts as function of threshold")
		self.declareProperty(WorkspaceProperty("Positrons","",Direction.Output,PropertyMode.Optional),doc="Positron peak area as function of threshold")
		self.declareProperty(FileProperty("IDF","",FileAction.OptionalLoad, extensions = ["xml"]),doc="IDF for the result files, for Instrument View")


	def category(self):
		return 'Muon'

	def PyExec(self):
		# Run the algorithm
		self.timeSetup=0.0
		self.callsSetup=0
		self.timeExec=0.0
		self.callsExec=0
		self.timeSetPar=0.0
		self.callsSetPar=0


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
		lamguess=self.getProperty('LambdaGuess').value
		nmax=self.getProperty('NumberOfRuns').value
		t1=self.getProperty('FitDataFrom').value
		fml=self.getProperty('FixMuonLifetime').value
		fbz=self.getProperty('FixBackgroundToZero').value
		agroup=self.getProperty('Grouped').value
		indPhases=self.getProperty("IndividualPhase").value
		minsum=self.getProperty("MinThresholdForSum").value
		bbefore=self.getProperty("BackgroundBefore").value
		pktime=self.getProperty("PositronPeakTime").value

		thresholdsGiven=self.getProperty("ThresholdList").value
		if(len(thresholdsGiven)>0):
			nmax=len(thresholdsGiven)
	
		# first: group all data and fit to get B? take guess provided by user
		# then: group across runs and get phase (and freq) for each spectrum
		# then: individual fits for bg, asym, deadtime, countrate only.

		# dead time correction and background need high stats to be meaningful

		# allow "muon lifetime" to be variable too if TDC timebase is slightly off

		# progress. Loading phase is 0.1, main fitting 0.85, final analysis 0.05
		# Progress in Load assumes 100 runs unless number given
		# awaiting correct way to use this...
		#progSoFar=
		
		# prepare recyclable child algorithms for speed
		loadAlg="LoadMuonNexus" # AlgorithmManager.create("LoadMuonNexus")
		loadInstAlg=AlgorithmManager.create("LoadInstrument")
		groupDetAlg=AlgorithmManager.create("GroupDetectors")
		asymAlg=AlgorithmManager.create("AsymmetryCalc")
		mgrpAlg=AlgorithmManager.create("MuonGroupDetectors")
		# Fit is fussy!
		fitAlg="Fit" # AlgorithmManager.create("Fit")
		maskAlg=AlgorithmManager.create("MaskDetectors")
		plusAlg=AlgorithmManager.create("Plus")
		rebunchAlg=AlgorithmManager.create("Rebunch")

		# load the workspaces and add them keeping individuals too
		voltlist=[]
		runhandles=[]
		runnames=[]
		vi=0
		oldvolts=-1
		thUnit="mV" # default, can override
		sumd=None
		while(vi<nmax and oldvolts<999999):
			#runlist.append(vi+firstnum)
			thisrunname=firstfn[k:j]+str(vi+firstnum).zfill(i-j)
			thisrunpath=firstfn[:j]+str(vi+firstnum).zfill(i-j)+firstfn[i:]
	  
			try:
				self.log().notice("calling LoadMuonNexus with "+thisrunpath+" and "+thisrunname)
				if(self.getProperty("GroupingTable").isDefault):
					thisoneX=self.doChildAlg(loadAlg,Filename=thisrunpath,OutputWorkspace=thisrunname,AutoGroup=agroup,DetectorGroupingTable="dgt")
					dgt=thisoneX["DetectorGroupingTable"]
				else:
					thisoneX1=self.doChildAlg(loadAlg,Filename=thisrunpath,OutputWorkspace=thisrunname,AutoGroup=False,DetectorGroupingTable="dgt0")
					if(not self.getProperty("IDF").isDefault):
						thisoneY=self.doChildAlg(loadInstAlg,Workspace=thisoneX1["OutputWorkspace"],RewriteSpectraMap=False,Filename=self.getProperty("IDF").value)
					dgt=self.getProperty("GroupingTable").value
					glist=[dgt.cell(iii,0) for iii in range(dgt.rowCount())]
					gpat=",".join(["+".join([str(y-1) for y in x]) for x in glist])
					thisoneX=self.doChildAlg(groupDetAlg,InputWorkspace=thisoneX1["OutputWorkspace"],GroupingPattern=gpat,OutputWorkspace=thisrunname)
					#thisoneX=self.doChildAlg("MuonGroupDetectors",InputWorkspace=thisoneX1["OutputWorkspace"],DetectorGroupingTable=dgt,OutputWorkspace=thisrunname)
				#thisoneX=LoadMuonNexus(filename=thisrunpath,OutputWorkspace=thisrunname,AutoGroup=agroup,DetectorGroupingTable="dgt")
			except IOError as e:				
				self.log().notice("no more runs to load; "+str(e))
				break
			thisone=thisoneX["OutputWorkspace"]

			# rebunch if bins are too small
			if thisone.dataX(0)[1]-thisone.dataX(0)[0]<0.015:
				self.log().notice("rebunching run")
				rb=self.doChildAlg(rebunchAlg,InputWorkspace=thisone,NBunch=16,OutputWorkspace=thisrunname)
				thisone=rb["OutputWorkspace"]
				
			self.log().notice("Comment field: ["+thisone.getComment()+"]")
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
				self.log().notice("threshold incrementing from "+str(oldvolts)+" to "+str(volts))
				oldvolts=volts
			else:
				oldvolts=999999
				break
			# exclude low threshold runs which might have background
			if(volts>=minsum):
				if(sumd is None):
					sumd=thisone
				else:
					sumd=self.doChildAlg(plusAlg,LHSWorkspace=sumd,RHSWorkspace=thisone,OutputWorkspace="sumd")["OutputWorkspace"]
					#sumd = sumd + thisone # note creates new, not adding into original ws!
			vi=vi+1

		numruns=len(voltlist)
		self.log().notice("found sequence of "+str(numruns)+" runs with thresholds between "+str(voltlist[0])+" and "+str(voltlist[-1])+thUnit)

		# sort runs in order of increasing voltage leading to neater plots later
		ordering=numpy.argsort(voltlist)
		voltlist=[voltlist[ordering[i]] for i in range(len(ordering))]
		runnames=[runnames[ordering[i]] for i in range(len(ordering))]
		runhandles=[runhandles[ordering[i]] for i in range(len(ordering))]
		
		# fitting..

		# loop through spectra
		nhtp=sumd.getNumberHistograms()
		#nhtp=10 # subset for testing
		proc=WorkspaceFactory.create(sumd,NVectors=nhtp,XLength=numruns,YLength=numruns) # rates
		procA=WorkspaceFactory.create(sumd,NVectors=nhtp,XLength=numruns,YLength=numruns) # Asyms
		procB=WorkspaceFactory.create(sumd,NVectors=nhtp,XLength=numruns,YLength=numruns) # BGs
		procD=WorkspaceFactory.create(sumd,NVectors=nhtp,XLength=numruns,YLength=numruns) # DTs
		procP=WorkspaceFactory.create(sumd,NVectors=nhtp,XLength=numruns,YLength=numruns) # phases, if used
		procX=WorkspaceFactory.create(sumd,NVectors=nhtp,XLength=numruns,YLength=numruns) # chi squared, if used
		if(not self.getProperty("Positrons").isDefault):
			doPositrons=True
			procPos=WorkspaceFactory.create(sumd,NVectors=nhtp,XLength=numruns,YLength=numruns) # phases, if used
			self.setProperty("Positrons",procPos)
		else:
			procPos=None
			doPositrons=False

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
		if(not self.getProperty("ChiSquares").isDefault):
			self.setProperty("ChiSquares",procX)
		# X axes are thresholds in mv
		for w in (proc,procA,procB,procD,procP,procPos,procX):
			na=NumericAxis.create(numruns)
			for i,v in enumerate(voltlist):
				na.setValue(i,v)
				lbl=na.setUnit("Label")
				lbl.setLabel("Threshold",thUnit)
			if(w):
				w.replaceAxis(0,na)
		proc.setYUnitLabel("Counts per frame")
		procA.setYUnitLabel("Asymmetry")
		procB.setYUnitLabel("Background counts/frame")
		procD.setYUnitLabel("Dead time (us)")
		procP.setYUnitLabel("Phase (rad)")
		procX.setYUnitLabel("Chi squared")
		if(doPositrons):
			procPos.setYUnitLabel("Positron peak counts/frame")
		# IDF already reloaded into raw data and workspaces based on it
		#if(not self.getProperty("IDF").isDefault):
			# CHRONUS, etc: load alternative IDF. Do now since running this afterwards clears the mask flags!
			#for pp in [proc,procA,procB,procD,procP,procPos,procX]:
				#if(pp):
					#self.doChildAlg("LoadInstrument",Workspace=pp,RewriteSpectraMap=False,Filename=self.getProperty("IDF").value)
					#iloader=AlgorithmManager.createUnmanaged("LoadInstrument")
					#iloader.setChild(True) # doesn't need a name
					#iloader.initialize()
					#iloader.setProperty("Workspace", pp)
					#iloader.setProperty("RewriteSpectraMap", True)
					#iloader.setProperty("Filename", self.getProperty("IDF").value)
					#iloader.execute()
			# LoadInstrument(Workspace=proc,RewriteSpectraMap=True,Filename=self.getProperty("IDF").value)

		prog=Progress(self,start=0.0,end=1.0,nreports=nhtp)

		# initial fit to summed data to get actual frequency and lambda
		if(bguess>0):
			if(agroup):
				#bulkasym=AsymmetryCalc(sumd,Alpha=1.0)
				asyc=self.doChildAlg(asymAlg,InputWorkspace=sumd,Alpha=1.0,ForwardSpectra=[0],BackwardSpectra=[nhtp-1],OutputWorkspace="asyc")
				#asyc=AlgorithmManager.create("AsymmetryCalc")
				#asyc.setChild(True)
				#asyc.setProperty("InputWorkspace",sumd)
				#asyc.setProperty("Alpha",1.0)
				#asyc.execute
				bulkasym=asyc["OutputWorkspace"]
			else:
				#bulkgrp=MuonGroupDetectors(sumd,dgt)
				mgd=self.doChildAlg(mgrpAlg,InputWorkspace=sumd,DetectorGroupingTable=dgt,OutputWorkspace="mgd")
				#mgd=AlgorithmManager.create("MuonGroupDetectors")
				#mgd.setChild(True)
				#mgd.setProperty("InputWorkspace",sumd)
				#mgd.setProperty("DetectorGroupingTable",dgt)
				#mgd.execute
				bulkgrp=mgd["OutputWorkspace"]
				#bulkasym=AsymmetryCalc(bulkgrp,Alpha=1.0)
				asyc=self.doChildAlg(asymAlg,InputWorkspace=bulkgrp,Alpha=1.0,OutputWorkspace="asyc")
				#asyc=AlgorithmManager.create("AsymmetryCalc")
				#asyc.setChild(True)
				#asyc.setProperty("InputWorkspace",bulkgrp)
				#asyc.setProperty("Alpha",1.0)
				#asyc.execute
				bulkasym=asyc["OutputWorkspace"]
			if(lamguess<0):
				edo=ExpDecayOsc(A=0.2,Frequency=bguess*0.01355,Phi=0.1,Lambda=-lamguess)
			else:
				edo=ExpDecayOsc(A=0.2,Frequency=bguess*0.01355,Phi=0.1,Lambda=lamguess)
				edo.fix("Lambda")
			bulkfunc=edo + FlatBackground(A0=0.0)
			#fitalg=AlgorithmManager.create("Fit")
			#fitalg.setChild(True)
			#fitalg.setProperty("InputWorkspace",bulkasym)
			#fitalg.setProperty("Function",bulkfunc)
			#fitalg.setProperty("")
			fitalg=self.doChildAlg(fitAlg,InputWorkspace=bulkasym,Function=str(bulkfunc),StartX=t1,EndX=t1+12.0,CreateOutput=True,OutputParametersOnly=True,Output="bulk")
			#ff=Fit(InputWorkspace=bulkasym,Function=bulkfunc,StartX=t1,EndX=t1+12.0,CreateOutput=True,Output="bulk")
			fr=dict([(r["Name"],r["Value"]) for r in fitalg["OutputParameters"]])
			bfitted=fr["f0.Frequency"]*2*math.pi
			lamfitted=fr["f0.Lambda"]
			phizero=fr["f0.Phi"]
			self.log().notice("Fitted frequency "+str(bfitted)+" Mrad/s, lambda "+str(lamfitted)+" us-1")

		for i in range(0,nhtp):
			# fit the summed runs first (if TF)
			# check for dead detector and mask it
			if(numpy.sum(sumd.dataY(i))<1):
				for pp in [proc,procA,procB,procD,procP,procPos,procX]:
					if(pp):
						pp.dataX(i)[:]=voltlist # do need X values!
						self.doChildAlg(maskAlg,Workspace=pp,WorkspaceIndexList=[i])
						#masker = AlgorithmManager.createUnmanaged("MaskDetectors")
						#masker.setChild(True) # doesn't need a name
						#masker.initialize()
						#masker.setProperty("Workspace", pp) 
						#masker.setProperty("WorkspaceIndexList", [i]) # or the appropriate index
						#masker.execute()

			else:
				if(bguess>0):
					if(not fml) or (not indPhases): # need to fit this histogram, summed, to get the phase and/or tau
						if(i+1) in dgt.cell(0,0): # forward detector list
							phguess=phizero
						else:
							phguess=phizero+math.pi
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
						#print "fit function will be '"+"name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+p)*exp(-x*"+str(lamfitted)+"))*exp(-x/tau)+b)+d), n="+str(nguess)+", a=0.2, p="+str(phguess)+", b=0.01, d=0.001, tau=2.19703"+"'"
						ff=self.doChildAlg(fitAlg,Function="name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+p)*exp(-x*"+str(lamfitted)+"))*exp(-x/tau)+b)+d), n="+str(nguess)+", a=0.2, p="+str(phguess)+", b=0.01, d=0.001, tau=2.19703",
							InputWorkspace=sumd,StartX=t1,EndX=t1+12,CreateOutput=True,OutputParametersOnly=True,Ties=",".join(ties),WorkspaceIndex=i)
						#ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+p)*exp(-x/"+str(lamfitted)+"))*exp(-x/tau)+b)+d), n="+str(nguess)+", a=0.2, p="+str(phguess)+", b=0.01, d=0.001, tau=2.19703",
						#	InputWorkspace="Summed",StartX=t1,EndX=t1+12,Output="fitresSpec",Ties=",".join(ties),WorkspaceIndex=i)
						#for rrr in ff["OutputParameters"]:
						#	print rrr
						fr=dict([(r["Name"],r["Value"]) for r in ff["OutputParameters"]])
						afitted=fr["a"]
						#bfitted=fr["b"]
						phfitted=fr["p"]
						taufitted=fr["tau"]
						if(afitted<0):
							phfitted=phfitted+math.pi # correct for phase jumping negative
							afitted=-afitted
						phfitted=phfitted-2*math.pi*math.floor((phfitted-phizero)/2/math.pi+0.25) # wrap to range -pi/2 .. +3pi/2 with respect to phizero
						self.log().notice("hist "+str(i)+" freq="+str(bfitted)+" ampl="+str(afitted)+" phase "+str(phfitted)+"(from "+str(phguess)+") tau="+str(taufitted))
					else:
						if(fml):
							taufitted=2.19703
						if(i+1) in dgt.cell(0,0): # forward detector list
							phguess=phizero
						else:
							phguess=phizero+math.pi
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
						# maybe fit background from before t0 and fix it
						elif(bbefore<32):
							ff=self.doChildAlg(fitAlg,Function=str(FlatBackground()),InputWorkspace=runhandles[vi],WorkspaceIndex=i,StartX=-999.0,EndX=bbefore,CreateOutput=True,OutputParametersOnly=True)
							#ff=Fit(Function=FlatBackground(),InputWorkspace=runnames[vi],WorkspaceIndex=i,StartX=-999.0,EndX=bbefore)
							frbg=dict([(r["Name"],r["Value"]) for r in ff["OutputParameters"]])
							ties.append("b="+str(frbg["A0"]))
						if(not indPhases):
							ties.append("p="+str(phfitted))
						#ff.execute()
						ff=self.doChildAlg(fitAlg,Function="name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+p)*exp(-x*"+str(lamfitted)+"))*exp(-x/"+str(taufitted)+")+b)+d), n="+str(nguess)+", a=0.2, b=0.01, d=0.001, p="+str(phguess),
							InputWorkspace=runhandles[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,CreateOutput=True,OutputParametersOnly=True,Ties=",".join(ties))
						#ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*(1+a*cos(x*"+str(bfitted)+"+p)*exp(-x*"+str(lamfitted)+"))*exp(-x/"+str(taufitted)+")+b)+d), n="+str(nguess)+", a=0.2, b=0.01, d=0.001, p="+str(phguess),
						#	InputWorkspace=runnames[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,Output="fitres",Ties=",".join(ties))
						fr=dict([(r["Name"],r["Value"]) for r in ff["OutputParameters"]])
						fe=dict([(r["Name"],r["Error"]) for r in ff["OutputParameters"]])
						Nfitted=fr["n"]
						NfitE=fe["n"]
						Afitted=fr["a"]
						AfitE=fe["a"]
						Pfitted=fr["p"]
						PfitE=fe["p"]
						Bfitted=fr["b"]
						BfitE=fe["b"]
						Dfitted=fr["d"]
						DfitE=fe["d"]
						# check phase
						if (Afitted<0):
							Afitted=-Afitted
							Pfitted=Pfitted+math.pi
						Pfitted=Pfitted-2*math.pi*math.floor((Pfitted-phizero)/2/math.pi+0.25) # wrap to range -pi/2 .. +3pi/2
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
						procX.dataX(i)[vi]=voltlist[vi]
						procX.dataY(i)[vi]=fr["Cost function value"]
						procX.dataE(i)[vi]=0
						# and also fit the positron peak (fixed position/width Gaussian)
						if(doPositrons):
							ff=self.doChildAlg(fitAlg,Function=str(FlatBackground()),InputWorkspace=runhandles[vi],WorkspaceIndex=i,StartX=pktime-0.05,EndX=pktime+0.05,CreateOutput=True,OutputParametersOnly=True)
							#ff=Fit(Function=FlatBackground(),InputWorkspace=runnames[vi],WorkspaceIndex=i,StartX=pktime-0.05,EndX=pktime+0.05,CreateOutput=True,Output="posfit")
							pfr=dict([(r["Name"],r["Value"]) for r in ff["OutputParameters"]])
							pfe=dict([(r["Name"],r["Error"]) for r in ff["OutputParameters"]])
							procPos.dataX(i)[vi]=voltlist[vi]
							procPos.dataY(i)[vi]=(pfr["A0"]-fr["b"])/nframes # counts/bin/frame
							procPos.dataE(i)[vi]=pfe["A0"]/nframes # assume BG error is small by comparison
													
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
						ff=self.doChildAlg(fitAlg,Function="name=UserFunction, Formula=1/(1/(n*exp(-x/tau)+b)+d), n="+str(nguess)+", b=0.01, d=0.001, tau=2.19703",
							InputWorkspace=runhandles[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,CreateOutput=True,OutputParametersOnly=True,Ties=",".join(ties))
						#ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*exp(-x/tau)+b)+d), n="+str(nguess)+", b=0.01, d=0.001, tau=2.19703",
						#	InputWorkspace=runnames[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,Output="fitres",Ties=",".join(ties))
						fr=ff["OutputParameters"] # mtd["fitres_Parameters"]
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
						procX.dataX(i)[vi]=voltlist[vi]
						procX.dataY(i)[vi]=fr.cell("Value",4)
						procX.dataE(i)[vi]=0
				elif (bguess==-2):
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
						#ff.execute()
						ff=self.doChildAlg(fitAlg,Function="name=Chebyshev,n=10,StartX="+str(t1)+",EndX="+str(t1+4),
							InputWorkspace=runhandles[vi],WorkspaceIndex=str(i),StartX=t1,EndX=t1+4,CreateOutput=True,OutputParametersOnly=True)
						#ff=Fit(Function="name=UserFunction, Formula=1/(1/(n*exp(-x/tau)+b)+d), n="+str(nguess)+", b=0.01, d=0.001, tau=2.19703",
						#	InputWorkspace=runnames[vi],WorkspaceIndex=str(i),StartX=str(t1),EndX=32,Output="fitres",Ties=",".join(ties))
						fr=dict([(r["Name"],r["Value"]) for r in ff["OutputParameters"]])
						fe=dict([(r["Name"],r["Error"]) for r in ff["OutputParameters"]])
						proc.dataX(i)[vi]=voltlist[vi]
						proc.dataY(i)[vi]=fr["A0"]/nframes/tres # counts/us/frame, mean including any BG
						proc.dataE(i)[vi]=fe["A0"]/nframes/tres
						procA.dataX(i)[vi]=voltlist[vi]
						procA.dataY(i)[vi]=0.0 # no asym possible for LF/ZF runs
						procA.dataE(i)[vi]=0.0
						procB.dataX(i)[vi]=voltlist[vi]
						procB.dataY(i)[vi]=0.0 
						procB.dataE(i)[vi]=0.0
						procD.dataX(i)[vi]=voltlist[vi]
						procD.dataY(i)[vi]=0.0
						procD.dataE(i)[vi]=0.0
						procX.dataX(i)[vi]=voltlist[vi]
						procX.dataY(i)[vi]=fr["Cost function value"]
						procX.dataE(i)[vi]=0
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
					
			#print "voltlist is",voltlist
			#print "runnames in order are ",runnames
			#print "Start X for fit is ",x0
			#print "x values of data set are",proc.dataX(i)
			if (not numpy.all(numpy.isfinite(proc.dataX(i))) ) or (not numpy.all(numpy.isfinite(proc.dataY(i))) ) or (not numpy.all(numpy.isfinite(proc.dataE(i))) ):
				self.log().warning("about to fit model to spectrum "+str(i))
				self.log().warning("proc.dataX="+str(proc.dataX(i)))
				self.log().warning("proc.dataY="+str(proc.dataY(i)))
				self.log().warning("proc.dataE="+str(proc.dataE(i)))
				print ("Bad data in histogram",i,"Will crash Fit...")
			else:
				ff=self.doChildAlg(fitAlg,Function="name=UserFunction, Formula=rate/(1+exp(x/lam - 2.75)), rate=5.0, lam=20.0",InputWorkspace=proc,StartX=x0,EndX="999.0",CreateOutput=True,OutputParametersOnly=True,WorkspaceIndex=str(i))
				#ff=Fit(Function="name=UserFunction, Formula=rate/(1+exp(x/lam - 2.75)), rate=5.0, lam=20.0",InputWorkspace=proc,StartX=x0,EndX="999.0",Output="fitModel",WorkspaceIndex=str(i))
				fr=ff["OutputParameters"] # mtd["fitres_Parameters"]
				resw.addRow((i,fr.cell("Value",0),fr.cell("Error",0),fr.cell("Value",1),fr.cell("Error",1),fr.cell("Value",2)))

		# output workspaces already set...
		if(self.callsSetup>0):
			self.log().notice("created "+str(self.callsSetup)+" new Fits, "+str(self.timeSetup/self.callsSetup)+" seconds each")
		else:
			self.log().notice("Don't seem to have created any Fit algorithms")
		if(self.callsSetPar>0):
			self.log().notice("set Fit properties "+str(self.callsSetPar)+" times, "+str(self.timeSetPar/self.callsSetPar)+" seconds each batch")
		if(self.callsExec>0):
			self.log().notice("executed Fit "+str(self.callsExec)+" times, "+str(self.timeExec/self.callsExec)+" seconds each")

AlgorithmFactory.subscribe(AnalyseThresholdScans)
