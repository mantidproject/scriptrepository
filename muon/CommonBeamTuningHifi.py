"""
Scan a magnet, take data on three machines, analyse it, plot results
"""
from __future__ import print_function

# HIFI only version! (for now)

# deal with Limits on Blocks. Does cset(block, value=too_large) raise an exception?

# return to default/original at end of scan
# auxiliary tables named xx_Reference with cols "name" "value"
# columns for all magnets
# "last setpoint" table named "SE_Block_Cache", same format
CACHENAME="SE_Block_Cache"

# select test mode
# DAQ=True (simulated), "abort" (real DAE but don't save) or False (real working system)
TESTMODE_DAQ=False
TESTMODE_BLOCKS=False

# script top level (maybe not Algorithm unless Workflow is debugged)
# stage 0 get input
# multi section scan pars permitted, as Rebin() - x1 dx12 x2 dx23 x3

# stage 1 run ScanAndMeasureMagnet
# output an intermediate table? has headings "scanvar" "inst1" "inst2" and rows <x> <run1> <run2>
# can edit this table to insert old data and go on to run AnalyseTuningScan

# stage 2 run AnalyseTuningScan

# table "comment" fields gives starting currents/slits and scan variable?

# if old table given:
# calculate extra points so that no gap exceeds specified step. Generally have equal sized gaps in between old points. Don't repeat old points.
# re-run all analysis (old and new points).

# table column names
# left column type x, named for Magnet (or other tunable)
# for instruments, always prefix with inst name, then...
# Run Rate RateErr Asym AsymErr Alpha AlphaErr Prompt PromptErr DoubleP DoublePErr
# for cameras, suffixes
# File Xpos XposErr Ypos YposErr Width WidthErr Height HeightErr Ellip EllipErr Intens IntensErr Backgd BackgdErr

# auxiliary table for scan conditions, named <table>_Conditions
# three columns
# col 1 for scanned name, magnet and slit names, instrument and camera names (run number ranges)
# col 2 for starting value
# col 3 for end value (same if that magnet not scanned)

from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty
from mantid.kernel import *
from mantid.simpleapi import *
from logger import Logger
#import inspect
#from genie_python.genie_startup import *
import os
import re
import time
import numpy
import math

####################################
# SETUP PARAMETERS
####################################

# camera control subroutine(s)
# params are fileformat, lastRunNum, flag -1=initialise 0=just started data taking, 1=about to end
# return tuple (run num,flag) for this time, valid on 2nd call. flag=True if good data taken
def getRecentFile(fmt,lastRunNum,flag):
	# continuous dumping of files, use latest one
	if(flag==-1):
		dir=os.path.dirname(fmt)
		matcher=re.compile("^"+os.path.basename(fmt).replace("%d","([0-9]+)")+"$",re.IGNORECASE)
		num=0
		allFiles=os.listdir(dir)
		for f in allFiles:
			m=matcher.match(f)
			if(m):
				num=max(num,int(m.group(1)))
		return (num,False) # highest numbered matching file
	else:
		N=lastRunNum
		while (os.access(fmt % (N+1), os.F_OK)):
			N=N+1
		if(flag==1):
			if(N-lastRunNum>=2):
				return (N,True)
			else:
				return (N,False)
		else:
			return (N,False)

def processHifiCam(filename):
	img=LoadBeamCameraFile(SourceFile=filename, FilterNoiseLevel='50', BinSize='1', OutputWorkspace='IMG',EnableLogging=False)
	#img=LoadFITS(Filename=filename, LoadAsRectImg=True, FilterNoiseLevel='50', BinSize='1', Scale='80', OutputWorkspace='IMG')
	squash=MaskAndFlattenCameraImage(img,120,590,80,440,EnableLogging=False)
	
	try:
		fplist=Fit(Function="name=TwoDimGaussianPlusBG,Intensity=40000",
			InputWorkspace='squash',
			Output='squash',
			OutputCompositeMembers=True,
			StartX=0,
			EndX=100000000,
			EnableLogging=False)
	
		#print fplist
		(YesNo, CostFn, CovarMat, Params,Workspace,Function,CostFn2)=fplist	
		(X0,Y0,XSig,YSig,Skew,Background,Intens,CostF2)=Params.column(1)
		(eX0,eY0,eXSig,eYSig,eSkew,eBackground,eIntens,eCostF2)=Params.column(2)
	except Exception as e:
		print("error fitting or processing? ",e)
		(X0,Y0,XSig,YSig,Skew,Background,Intens,CostF2)=[0.0]*8
		(eX0,eY0,eXSig,eYSig,eSkew,eBackground,eIntens,eCostF2)=[0.0]*8
		
	return (X0,eX0,Y0,eY0,XSig,eXSig,YSig,eYSig,Skew,eSkew,Intens,eIntens,Background,eBackground)
	
WAITFOR_TIMEOUT_SECONDS = 60
MAGNETS=["SEPTUM_A","UQ13A","UQ14A","UQ15A"]
#MAGNETS=["UQ1","UQ2","UQ3","UQ4","UQ5","UQ6","UQ7","UQ8","UQ9","UQ10","UQ11","UQ12","UB1","UB2","SEPARATOR","SEPTUM_A","SEPTUM_C","UQ13A","UQ14A","UQ15A","UQ13B","UQ14B","UQ13C","UQ14C","UQ15C"] # example
MAGNETSETS=["H131415A","V131415A","M131415A"]
#MAGNETSETS=["H123","V123","M123","H456","V456","M456","H789","V789","M789","H101112","V101112","M101112","H131415A","V131415A","M131415A","H1314B","V1314B","H131415C","V131415C","M131415C","Momentum"]
OTHERS=[] # ["Mslits"] not yet connected this way

BLOCKNAMES={"UQ1":"UQ1_CURR","UQ2":"UQ2_CURR","UQ3":"UQ3_CURR",
	"UQ4":"UQ4_CURR","UQ5":"UQ5_CURR","UQ6":"UQ6_CURR","UQ7":"UQ7_CURR",
	"UQ8":"UQ8_CURR","UQ9":"UQ9_CURR","UQ10":"UQ10_CURR","UQ11":"UQ11_CURR",
	"UQ12":"UQ12_CURR","UB1":"UB1_CURR","UB2":"UB2_CURR","SEPARATOR":"SEPARATOR_CURR",
	"SEPTUM_A":"SEPTUM_A_CURR","SEPTUM_C":"SEPTUM_C_CURR","UQ13A":"UQ13A_CURR",
	"UQ14A":"UQ14A_CURR","UQ15A":"UQ15A_CURR","UQ13B":"UQ13B_CURR",
	"UQ14B":"UQ14B_CURR","UQ13C":"UQ13C_CURR","UQ14C":"UQ14C_CURR","UQ15C":"UQ15C_CURR"}

# *********** override for tests
#MAGNETS=["SEPTUM_A_CURR","SEPTUM_C_CURR"]
#OTHERS=[]

# full list for convenience
TUNABLES=MAGNETS+MAGNETSETS+OTHERS

DPARS={} # how to analyse data
#DPARS["MUSR"]={"machine":"NDXMUSR","data":"//musr/data/musr%08d.nxs","tgbegin":0.5,"tgend":10.0,"p1begin":-0.23,"p1end":-0.1,"promptbegin":-0.4,"promptend":-0.3,"pulse":2}
#DPARS["EMU"]={"machine":"NDXEMU","data":"//emu/data/emu%08d.nxs","tgbegin":0.5,"tgend":10.0,"p1begin":0.1,"p1end":0.23,"promptbegin":-0.2,"promptend":0.0,"pulse":1}
DPARS["HIFI"]={"machine":"NDXHIFI","data":"//hifi/data/hifi%08d.nxs","tgbegin":0.5,"tgend":10.0,"p1begin":0.1,"p1end":0.23,"promptbegin":-0.2,"promptend":0.0,"pulse":1}
MACHINES0 = list(DPARS.keys())

CPARS={} # how to get camera data
CPARS["HIFI"]={"data":"//ndlt626/Users/Public/Documents/sxvh9usb/HIFI Sept2018/IMG%d.fit","finder":getRecentFile,"processor":processHifiCam}
#CPARS["MUSR"]={finder=TakeMusrPic,processor=processMusrCam} # filenames specified by user
CAMERAS0=list(CPARS.keys())


if(TESTMODE_DAQ is True):
	DPARS["MUSR"]["data"]="//musr/data/musr%08d.nxs"
	DPARS["EMU"]["data"]="//emu/data/Cycle 14_3/emu%08d.nxs"
	DPARS["HIFI"]["data"]="//hifi/data/cycle_14_3/hifi%08d.nxs"
	class MultipleDaeController:
		def __init__(self):
			self.startrunnums={"NDXMUSR":51909,"NDXEMU":50690,"NDXHIFI":77061} # B1B2 tune starting 20 Mar 2015 13:52 approx
			self.runnums={}
			self.states={}
			self.frames={}
		def add_machine(self,mac,timeout):
			self.states[mac]="SETUP"
			self.frames[mac]=0
			self.runnums[mac]=self.startrunnums[mac]
			print("registered "+str(mac))
		def run_on_all(self,code,**kwds):
			if(code=="get_run_state"):
				return self.states
			if(code=="begin"):
				print("doing BEGIN, "+str(kwds))
				for m in list(self.states.keys()):
					self.states[m]="RUNNING"
					self.frames[m]=0
				if("waitfor_frames" in kwds):
					print("waiting for frames="+str(kwds["waitfor_frames"]))
					time.sleep(1)
					for m in list(self.states.keys()):
						self.frames[m]=kwds["waitfor_frames"]+1
			if(code=="get_run_number"):
				return self.runnums
			if(code=="end"):
				print("ending runs")
				for m in list(self.states.keys()):
					self.states[m]="SETUP"
					self.runnums[m]+=1
			if(code=="abort"):
				print("aborting runs")
				for m in list(self.states.keys()):
					self.states[m]="SETUP"
else:
	from multiple_dae_controller import MultipleDaeController
	import inspect

if(TESTMODE_BLOCKS):
	blockvals={}
	for b in MAGNETS+OTHERS:
		blockvals[BLOCKNAMES[b]]=50.0
	def cset(block=None,value=None,wait=False,lowlimit=-99999.9,highlimit=99999.9,**args):
		if(len(args)==1):
			block=list(args.keys())[0]
			value=list(args.values())[0]
		print("setting "+str(block)+" to "+str(value))
		blockvals[block]=value
		if(wait):
			print("waiting for value to stabilise")
			time.sleep(1)
	def cshow(block):
		print("CSHOW: "+str(blockvals[block]))
		return None
	def cget(block):
		if block in blockvals:
			return {"name":"TEST:"+block,"value":blockvals[block],"lowlimit":0.0,"highlimit":1.0}
		else:
			return None
	def waitfor_block(block,lowlimit=-99999.9,highlimit=99999.9):
		print("waiting for "+block+" to stabilise")
		time.sleep(1)
else:
	from genie_python.genie_startup import *
	set_instrument("IN:MUONFE:")
	
def cached_cset(block,val,wait=False,lowlimit=-99999.9,highlimit=99999.9):
	# translate magnet name to block name
	# set the block
	# also update the cache
	ca=mtd[CACHENAME]
	names=ca.column(0)
	args={BLOCKNAMES[block]:val} # ,"wait":wait
	cset(**args)
	if block in names:
		ca.setCell(names.index(block),1,val) # only on successful cset
	
def cached_get_block(block):
	# translate magnet name to block name
	# get block val from cache
	# get actual value
	# get value from table tt (_Reference) if supplied
	# return cached value if all close
	# or raise exception
	print("doing cget("+block+")")
	rbk=cget(BLOCKNAMES[block])["value"]
	ca=mtd[CACHENAME]
	cnames=ca.column(0)
	if(block in cnames):
		cval=ca.column(1)[cnames.index(block)]
		if(abs(rbk-cval)>1.0):
			raise Exception("Readback of "+block+"("+str(rbk)+") is not close to cached setpoint("+str(cval)+")")
		return cval
	else:
		return rbk # uncached

def raw_cget(block):
	# translate magnet name to block name
	return cget(BLOCKNAMES[block])["value"]


def setTune(startPars,magName,tuneVal):
	toSet={}
	if( (magName in MAGNETS) or (magName in OTHERS) ): # direct setting
		toSet[magName]=tuneVal
	elif(magName=="H123"): # delta amps (largest change) for all triplet options
		toSet["UQ1"]=startPars["UQ1"]+tuneVal*0.25
		toSet["UQ2"]=startPars["UQ2"]+tuneVal*1.0
		toSet["UQ3"]=startPars["UQ3"]+tuneVal*0.25
	elif(magName=="V123"):
		toSet["UQ1"]=startPars["UQ1"]+tuneVal*1.0
		toSet["UQ2"]=startPars["UQ2"]+tuneVal*0.4
		toSet["UQ3"]=startPars["UQ3"]+tuneVal*0.6
	elif(magName=="M123"):
		toSet["UQ1"]=startPars["UQ1"]+tuneVal*0.5
		toSet["UQ2"]=startPars["UQ2"]+tuneVal*-0.5
		toSet["UQ3"]=startPars["UQ3"]+tuneVal*-1.0
	elif(magName=="H456"):
		toSet["UQ4"]=startPars["UQ4"]+tuneVal*1.0
		toSet["UQ5"]=startPars["UQ5"]+tuneVal*0.67
		toSet["UQ6"]=startPars["UQ6"]+tuneVal*1.0
	elif(magName=="V456"):
		toSet["UQ4"]=startPars["UQ4"]+tuneVal*0.25
		toSet["UQ5"]=startPars["UQ5"]+tuneVal*1.0
		toSet["UQ6"]=startPars["UQ6"]+tuneVal*0.25
	elif(magName=="M456"):
		toSet["UQ4"]=startPars["UQ4"]+tuneVal*-1.0
		toSet["UQ6"]=startPars["UQ6"]+tuneVal*1.0
	elif(magName=="H789"):
		toSet["UQ7"]=startPars["UQ7"]+tuneVal*0.25
		toSet["UQ8"]=startPars["UQ8"]+tuneVal*1.0
		toSet["UQ9"]=startPars["UQ9"]+tuneVal*0.25
	elif(magName=="V789"):
		toSet["UQ7"]=startPars["UQ7"]+tuneVal*1.0
		toSet["UQ8"]=startPars["UQ8"]+tuneVal*0.67
		toSet["UQ9"]=startPars["UQ9"]+tuneVal*1.0
	elif(magName=="M789"):
		toSet["UQ7"]=startPars["UQ7"]+tuneVal*1.0
		toSet["UQ9"]=startPars["UQ9"]+tuneVal*-1.0
	elif(magName=="H101112"):
		toSet["UQ10"]=startPars["UQ10"]+tuneVal*1.0
		toSet["UQ11"]=startPars["UQ11"]+tuneVal*0.67
		toSet["UQ12"]=startPars["UQ12"]+tuneVal*1.0
	elif(magName=="V101112"):
		toSet["UQ10"]=startPars["UQ10"]+tuneVal*0.25
		toSet["UQ11"]=startPars["UQ11"]+tuneVal*1.0
		toSet["UQ12"]=startPars["UQ12"]+tuneVal*0.25
	elif(magName=="M101112"):
		toSet["UQ10"]=startPars["UQ10"]+tuneVal*-1.0
		toSet["UQ12"]=startPars["UQ12"]+tuneVal*1.0
	elif(magName=="H131415A"):
		toSet["UQ13A"]=startPars["UQ13A"]+tuneVal*0.2
		toSet["UQ14A"]=startPars["UQ14A"]+tuneVal*1.0
		toSet["UQ15A"]=startPars["UQ15A"]+tuneVal*0.4
	elif(magName=="V131415A"):
		toSet["UQ13A"]=startPars["UQ13A"]+tuneVal*0.6
		toSet["UQ14A"]=startPars["UQ14A"]+tuneVal*0.4
		toSet["UQ15A"]=startPars["UQ15A"]+tuneVal*1.0
	elif(magName=="M131415A"):
		toSet["UQ13A"]=startPars["UQ13A"]+tuneVal*-1.0
		toSet["UQ15A"]=startPars["UQ15A"]+tuneVal*1.0
	elif(magName=="H131415C"):
		toSet["UQ13C"]=startPars["UQ13C"]+tuneVal*0.25
		toSet["UQ14C"]=startPars["UQ14C"]+tuneVal*1.0
		toSet["UQ15C"]=startPars["UQ15C"]+tuneVal*0.25
	elif(magName=="V131415C"):
		toSet["UQ13C"]=startPars["UQ13C"]+tuneVal*0.6
		toSet["UQ14C"]=startPars["UQ14C"]+tuneVal*0.4
		toSet["UQ15C"]=startPars["UQ15C"]+tuneVal*1.0
	elif(magName=="M131415C"):
		toSet["UQ13B"]=startPars["UQ13B"]+tuneVal*-1.0
		toSet["UQ14C"]=startPars["UQ14C"]+tuneVal*-0.5
		toSet["UQ15C"]=startPars["UQ15C"]+tuneVal*0.5
	elif(magName=="H1314B"):
		toSet["UQ13B"]=startPars["UQ13B"]+tuneVal*1.0
		toSet["UQ14B"]=startPars["UQ14B"]+tuneVal*0.33
	elif(magName=="V1314B"):
		toSet["UQ13B"]=startPars["UQ13B"]+tuneVal*0.16
		toSet["UQ14B"]=startPars["UQ14B"]+tuneVal*1.0
	elif(magName=="Momentum"): # units of "B1"
		for mag in MAGNETS:
			if(mag=="SEPARATOR"):
				toSet["SEPARATOR"]=startPars["SEPARATOR"]*startPars["UB1"]/tuneVal
			else:
				toSet[mag]=startPars[mag]/startPars["UB1"]*tuneVal
	return toSet

class Beamline_cset(PythonAlgorithm):
	def PyInit(self):
		self.declareProperty("Operation","",StringListValidator(["Set to value","Refresh","Refresh All","Set All from Table"]))
		self.declareProperty("Magnet","",StringListValidator(MAGNETS+OTHERS))
		self.declareProperty("Value",0.0)
		self.declareProperty(ITableWorkspaceProperty("Table","",direction=Direction.Input,optional=PropertyMode.Optional))
	
	def category(self):
		return "Muon\\Tuning"

	def PyExec(self):
		if(CACHENAME not in mtd):
			# create and initialise?
			tt2=WorkspaceFactory.createTable()
			tt2.addColumn("str","MagnetName",6)
			tt2.addColumn("double","Reference",2)
			for mag1 in MAGNETS+OTHERS:
				tt2.addRow([mag1,-1.0]) # initially unknown
			AnalysisDataService.addOrReplace(CACHENAME,tt2)
			self.log().notice("Creating block cache")
		op=self.getProperty("Operation").value
		block=self.getProperty("Magnet").value
		if(op=="Set to value"):
			val=self.getProperty("Value").value
			cached_cset(block,val,wait=True)
		elif(op=="Refresh"):
			val=raw_cget(block) # bypass cache, get readback
			cached_cset(block,val)
		elif(op=="Refresh All"):
			for block in MAGNETS+OTHERS:
				val=raw_cget(block) # bypass cache, get readback
				cached_cset(block,val)
		elif(op=="Set All from Table"):
			# fairly free options in table. Meant for x_Reference tables, but can use Parameters table from a fit for example. Not all values need be set!
			t=self.getProperty("Table").value
			tblks=t.column(0)
			tvals=t.column(1)
			for block in tblks:
				if(block in (MAGNETS+OTHERS)):
					val=tvals[tblks.index(block)] # table entry
					cached_cset(block,val,wait=True)
				else:
					self.log().notice("Ignoring table row "+block)
		else:
			raise Exception("Unknown operation")
AlgorithmFactory.subscribe(Beamline_cset)

class TuneMagnet(DataProcessorAlgorithm):
	""" workflow algorithm, scan a magnet, analyse results, plot curve """
	def PyInit(self):
		#self.declareProperty(StringArrayProperty("Machines",MACHINES0))
		for mac in MACHINES0:
			self.declareProperty(mac,True,doc="Measure data from this instrument")
		for cam in CAMERAS0:
			self.declareProperty(cam+"cam",True,doc="Measure spot on this camera") # eg //laptop/directory/picture%d.fit"
		self.declareProperty("Magnet","",StringListValidator(TUNABLES))
		self.declareProperty(FloatArrayProperty("ScanRange",[0.0,1.0,10.0]),doc="first,step,[midpt,step2]*n,end, or blank just to re-analyse old results")
		self.declareProperty("Frames",5000)
		self.declareProperty(ITableWorkspaceProperty("OldResultTable","",direction=Direction.Input,optional=PropertyMode.Optional),doc="Already measured points to append to") # will copy from this table and append new points
		self.declareProperty(ITableWorkspaceProperty("ResultTable","",direction=Direction.Output)) # out; will analyse all new and old runs; string only here
		self.declareProperty("DoDeadTimeCorr",True)
		self.declareProperty("MeasureTF20asym",False)
		self.declareProperty("MeasurePromptPeaks",False)
		self.declareProperty("MeasureDoublePulse",False)
		self.declareProperty(FileProperty("AutoSaveDir","",action=FileAction.OptionalDirectory),doc="Somewhere to save the results")
	def category(self):
		return "Muon\\Tuning"

	def CalcScanValues(self,ScanRange,ExistingPoints,tolerance=0.3):
		#self.declareProperty(FloatArrayProperty("ScanRange",[0.0,1.0,10.0]),doc="first,step,[midpt,step2]*n,end")
		#self.declareProperty(FloatArrayProperty("ExistingPoints",[1.5,2.5]),doc="already measured, not to duplicate")
		#self.declareProperty(FloatArrayProperty("ScanValues",[0.0]),direction=Direction.Output,doc="chosen values to scan")
		# fill ranges
		# end and break points done exactly unless duplicates
		# don't add new points if within tolerance*stepsize of an existing point
		ScanR2=list(ScanRange)
		newPts=[]
		tol2=0.001 # default for a single point
		while(len(ScanR2)>2):
			ExPts=sorted(ExistingPoints+newPts) # existing, and new already added points
			tol2=tolerance*ScanR2[1]
			morePts=numpy.arange(ScanR2[0],ScanR2[2],ScanR2[1])
			for morePt in morePts:
				i=numpy.searchsorted(ExPts,morePt) # i points to old point just above new one
				if((i==0 or abs(ExPts[i-1]-morePt)>tol2) and(i==len(ExPts) or abs(ExPts[i]-morePt)>tol2)):
					newPts.append(morePt)
			del ScanR2[0:2]
		if(len(ScanR2)==1):
			ExPts=sorted(ExistingPoints+newPts) # existing, and new break points
			i=numpy.searchsorted(ExPts,ScanR2[0]) # i points to old point just above new one
			if((i==0 or abs(ExPts[i-1]-ScanR2[0])>tol2) and(i==len(ExPts) or abs(ExPts[i]-ScanR2[0])>tol2)):
				newPts.append(ScanR2[0])
		return newPts

	def PyExec(self):
		# get old scan values
		tunable=self.getProperty("Magnet").value

		if not self.getProperty("OldResultTable").isDefault:
			oldTab=self.getProperty("OldResultTable").value
			oldScanVals=oldTab.column(0)
			if(oldTab.getColumnNames()[0] != tunable):
				raise Exception("The old table was scanning a different magnet, can't use it")
			# check _References table matches current settings (allow one magnet being tuned to differ)
			tt2=mtd[self.getProperty("OldResultTable").valueAsStr+"_Reference"]
			for r in tt2:
				if(r["MagnetName"] != tunable and abs(r["Reference"]-cached_get_block(r["MagnetName"]))>0.01):
					raise Exception("Magnet "+r["MagnetName"]+" is not at the value used last time")
			
		else:
			oldTab=None
			oldScanVals=[]
		newScanPars=self.getProperty("ScanRange").value
		#print "newScanPars=", newScanPars
		newScanPts=self.CalcScanValues(newScanPars,oldScanVals)
		if(len(newScanPts)==0):
			self.log().notice("No values to scan, or already measured")
			runNumList=oldTab
			oldTab=None
		else:
			MACHINES=[]
			for mac in MACHINES0:
				if(self.getProperty(mac).value):
					MACHINES.append(DPARS[mac]["machine"])
			CAMERAS=[]
			for cam in CAMERAS0:
				if(self.getProperty(cam+"cam").value):
					CAMERAS.append(cam)
			#runNumList=ScanMagnet(Machines=MACHINES,Cameras=CAMERAS,Magnet=Tunable,ScanValues=newScanPts,Frames=self.getProperty("Frames").value,RunNumList="__NewRuns")
			scmag=self.createChildAlgorithm("ScanMagnet",startProgress=0.0,endProgress=0.9)
			scmag.setProperty("Machines",MACHINES)
			scmag.setProperty("Cameras",CAMERAS)
			scmag.setProperty("Magnet",tunable)
			scmag.setProperty("ScanValues",newScanPts)
			scmag.setProperty("Frames",self.getProperty("Frames").value)
			scmag.setPropertyValue("RunNumList", "__NewRuns")
			#scmag.setAlwaysStoreInADS(True)
			scmag.execute()
			runNumList = scmag.getProperty("RunNumList").value
			AnalysisDataService.addOrReplace("__NewRuns", runNumList)

		ats=self.createChildAlgorithm("AnalyseTuningScan",startProgress=0.9,endProgress=1.0)
		if(oldTab):
			ats.setProperty("OldTable",oldTab)
		ats.setProperty("NewRuns",runNumList)
		ats.setProperty("OutputTable",self.getProperty("ResultTable").valueAsStr)
		ats.setProperty("DoDeadTimeCorr",self.getProperty("DoDeadTimeCorr").value)
		ats.setProperty("MeasureTF20asym",self.getProperty("MeasureTF20asym").value)
		ats.setProperty("MeasurePromptPeaks",self.getProperty("MeasurePromptPeaks").value)
		ats.setProperty("MeasureDoublePulse",self.getProperty("MeasureDoublePulse").value)
		#ats.setAlwaysStoreInADS(True)
		ats.execute()
		optab=ats.getProperty("OutputTable").value
		self.setProperty("ResultTable",optab)
		
		if(not self.getProperty("AutoSaveDir").isDefault):
			asdir=self.getProperty("AutoSaveDir").value
			AnalysisDataService.addOrReplace(self.getProperty("ResultTable").valueAsStr, optab)
			SaveNexus(InputWorkspace=self.getProperty("ResultTable").valueAsStr,Filename=os.path.join(asdir,self.getProperty("ResultTable").valueAsStr+".nxs"))
			SaveNexus(InputWorkspace=self.getProperty("ResultTable").valueAsStr+"_Reference",Filename=os.path.join(asdir,self.getProperty("ResultTable").valueAsStr)+"_Reference.nxs")
			
AlgorithmFactory.subscribe(TuneMagnet)

class ChooseBestValue(PythonAlgorithm):
	""" given scan table and best value of X, calculate and set magnet values """
	def PyInit(self):
		self.declareProperty(ITableWorkspaceProperty("ResultTable","",direction=Direction.Input))
		self.declareProperty("BestValue",-1000.0,direction=Direction.InOut,doc="Specify best X coord, or default to let algorithm guess") # default is to guess based on rates
		self.declareProperty("Parameter","",StringListValidator(["Rate","Spot"]))
	def category(self):
		return "Muon\\Tuning"
		
	def PyExec(self):
		x=self.getProperty("BestValue").value
		tab=self.getProperty("ResultTable").value
		what=self.getProperty("Parameter").value
		mag=tab.getColumnNames()[0]
		if(x==-1000):
			if(what=="Rate"):
				# guess best
				rmax=0.0
				for row in tab:
					rates=[]
					for (key,val) in list(row.items()):
						if(key[-4:]=="Rate"):
							rates.append(val)
					rthis=sum(rates)
					if(rthis>rmax):
						x=row[mag]
						rmax=rthis
			elif(what=="Spot"):
				# smallest spot, 1st set of params only
				rmin=99999.9
				for row in tab:
					width=-1
					height=-1
					for (key,val) in list(row.items()):
						if(key[-5:]=="Width"):
							width=val
						if(key[-6:]=="Height"):
							height=val
					rthis=math.sqrt(width**2+height**2)
					if(width>0 and height>0 and rthis<rmin):
						x=row[mag]
						rmin=rthis
			self.log().notice("Found best value "+mag+"="+str(x))
		refvals={}
		for mag0 in MAGNETS+OTHERS:
			refvals[mag0]=cached_get_block(mag0)
		toSet=setTune(refvals,mag,x)
		try:
			for tsb,tsv in list(toSet.items()):
				#print "doing cached_cset(",tsb,",",tsv,",wait=True)"
				cached_cset(tsb,tsv,wait=True) # all at once with one "Wait"
				#print "done the cached cset."
			self.setProperty("BestValue",x)
		except:
			self.log().warning("Couldn't set one of the magnets, resetting all other values")
			for tsb in list(toSet.keys()):
				tsv=refvals[tsb]
				cached_cset(tsb,tsv,wait=True) # all at once with one "Wait"
			self.setProperty("BestValue",-1000.0)
			toSet={}
		for mag0 in MAGNETS+OTHERS:
			if(mag0 in toSet):
				self.log().notice(mag0+" changed from "+str(refvals[mag0])+" to "+str(toSet[mag0]))
			else:
				self.log().notice(mag0+" unchanged, still "+str(refvals[mag0]))
AlgorithmFactory.subscribe(ChooseBestValue)

class ScanMagnet(PythonAlgorithm):
	def PyInit(self):
		self.declareProperty(StringArrayProperty("Machines",[])) # NDX.. names
		self.declareProperty(StringArrayProperty("Cameras",[])) # names
		self.declareProperty("Magnet","",StringListValidator(TUNABLES))
		self.declareProperty(FloatArrayProperty("ScanValues",[0.0]),doc="list of values to scan")
		self.declareProperty("Frames",5000) # only if New Runs
		self.declareProperty(ITableWorkspaceProperty("RunNumList","",direction=Direction.Output)) # output table, new x vals and run numbers only
	def category(self):
		return "Muon\\Tuning"
		
	def PyExec(self):
		# catch errors in cset and omit the out of range points from the result table
		# sanity check that the algorithm isn't already running!
		if (len(AlgorithmManager.runningInstancesOf(self.name()))>1):
			raise Exception("TuneMagnet is already running")
		
		MACHINES=self.getProperty("Machines").value
		CAMERAS=self.getProperty("Cameras").value
		
		Mag=self.getProperty("Magnet").value
		scanvals=self.getProperty("ScanValues").value
		frames=self.getProperty("Frames").value
		tt=WorkspaceFactory.createTable()
		tt.addColumn("double",Mag,1)
		for mac in MACHINES:
			mac0=mac.replace("NDX","",1)
			tt.addColumn("int",mac0+"Run",6)
		for cam in CAMERAS:
			tt.addColumn("int",cam+"File",6)

		Prog=Progress(self,start=0.0,end=1.0,nreports=2*len(scanvals))
		frames=self.getProperty("Frames").value
		controller = MultipleDaeController()
		for m in MACHINES:
			controller.add_machine(m, WAITFOR_TIMEOUT_SECONDS)
		controller.run_on_all("set_waitfor_sleep_interval",5)
		statuses=controller.run_on_all("get_run_state")
		errstr=""
		for m,status in list(statuses.items()):
			if(status != "SETUP"):
				errstr=errstr+" "+m+"="+status
		if(errstr != ""):
			raise Exception("Instruments are not all in SETUP state:"+errstr)
		# initialise cameras
		cnums=[0]*len(CAMERAS)
		for i,cam in enumerate(CAMERAS):
			cs=CPARS[cam]["finder"](CPARS[cam]["data"],0,-1)
			cnums[i]=cs[0]
			self.log().notice("camera "+str(i)+" starting file="+str(cnums[i]))

		refvals={}
		for mag0 in MAGNETS+OTHERS:
			refvals[mag0]=cached_get_block(mag0)

		for x in scanvals:
			try:
				toSet=setTune(refvals,Mag,x)
				#for (mag,value) in toSet.items():
				#	cset(mag,value,wait=True)
				for tsb,tsv in list(toSet.items()):
					self.log().notice("doing cached_cset("+str(tsb)+","+str(tsv)+")")
					cached_cset(tsb,tsv,wait=True) # all at once with one "Wait"
					#print "done cached cset."
			except:
				self.log().warning("Couldn't set magnets for x="+str(x))
			else: # successful cset, take data, report errors here
				# camera call 1
				for i,cam in enumerate(CAMERAS):
					cs=CPARS[cam]["finder"](CPARS[cam]["data"],cnums[i],0)
					cnums[i]=cs[0]
					#print "camera",i,"file now=",cnums[i]
				Prog.report("Data collection")
				controller.run_on_all("begin",waitfor_frames=frames) # assumes waits for this many frames to be reached on all before continuing
				runs0=controller.run_on_all("get_run_number") # new method which will return a list of the run numbers in progress?
				runs={}
				for m,r in list(runs0.items()):
					runs[m]=r[0] # real get_run_number returns tuple of number and an empty string!
				self.log().notice("runs in progress are"+str(runs))
				# ensure runs really have finished
				finished=False
				abandon=False
				while(not finished):
					finished=True
					statuses=controller.run_on_all("get_run_state")
					errstr=""
					for m,status in list(statuses.items()):
						if(status != "RUNNING"):
							errstr=errstr+" "+m+"="+status
					if(errstr != ""):
						self.log().warning("Instruments are not all in RUNNING state:"+errstr+"; abandoning scan and saving previous completed runs")
						abandon=True # to get out of for() scan loop
						break # out of while loop
					allframes=controller.run_on_all("get_frames")
					for m,ff in list(allframes.items()):
						self.log().notice(str(m)+" is at "+str(ff)+"frames out of "+str(frames))
						if(ff<frames):
							finished=False
					if(not finished):
						self.log().notice("Waiting again. Hit 'Abort' or 'End' on a DAE to abandon scan")
						waitagain=controller.run_on_all("waitfor_frames",frames)
				if(abandon):
					break # out of for loop
				# camera call 2
				for i,cam in enumerate(CAMERAS):
					cs=[0,False]
					ci=100
					while((not cs[1]) and (ci>0)):
						cs=CPARS[cam]["finder"](CPARS[cam]["data"],cnums[i],1)
						if(not cs[1]):
							self.log().notice("camera "+cam+" has not taken a good picture yet, waiting for it")
							time.sleep(5)
							ci=ci-1
					self.log().notice("chosen picture "+str(cs[0])+" from camera "+cam)
					cnums[i]=cs[0]
				if(TESTMODE_DAQ=="abort"):
					controller.run_on_all("abort") # end but don't fill the archive disks.
					for mac in list(runs.keys()):
						runs[mac]-=1 # current run not saved, use the previous run instead!
				else:
					controller.run_on_all("end") # end and save them all.
				self.log().notice("measured at x="+str(x)+" into runs "+str(runs))
				runs2=[]
				for mac in MACHINES:
					runs2.append(runs[mac])
				Prog.report("Data collection")
				tt.addRow([x]+runs2+cnums)

		for tsb in list(toSet.keys()):
			tsv=refvals[tsb]
			self.log().notice("doing cached_cset("+str(tsb)+","+str(tsv)+")")
			cached_cset(tsb,tsv,wait=True) # all at once with one "Wait"
		self.setProperty("RunNumList",tt)
		# create reference file
		tt2=WorkspaceFactory.createTable()
		tt2.addColumn("str","MagnetName",6)
		tt2.addColumn("double","Reference",2)
		for mag1 in MAGNETS+OTHERS:
			tt2.addRow([mag1,refvals[mag1]])
		AnalysisDataService.addOrReplace(self.getProperty("RunNumList").valueAsStr+"_Reference",tt2)
AlgorithmFactory.subscribe(ScanMagnet)

class AnalyseTuningScan(PythonAlgorithm):
	def PyInit(self):
		#self.declareProperty(StringArrayProperty("Machines",MACHINES0))
		#self.declareProperty("Magnet","",StringListValidator(TUNABLES))
		#self.declareProperty("LowerLimit",0.0)
		#self.declareProperty("UpperLimit",100.0)
		#self.declareProperty("Steps",10)
		#self.declareProperty("Frames",5000) # only if New Runs
		self.declareProperty(ITableWorkspaceProperty("OldTable","",direction=Direction.Input,optional=PropertyMode.Optional),doc="If already measured, run numbers on each machine") # table with old run numbers vs X to re-analyse; old analysis will be redone if present
		self.declareProperty(ITableWorkspaceProperty("NewRuns","",direction=Direction.Input),doc="run numbers on each machine") # table with newly measured run numbers to analyse
		self.declareProperty(ITableWorkspaceProperty("OutputTable","",direction=Direction.Output),doc="Full analysis as requested")
		self.declareProperty("DoDeadTimeCorr",True)
		self.declareProperty("MeasureTF20asym",False)
		self.declareProperty("MeasurePromptPeaks",False)
		self.declareProperty("MeasureDoublePulse",False)
	def category(self):
		return "Muon\\Tuning"

	def PyExec(self):
		self.setAlgStartupLogging=False # attempt to reduce unnecessary log messages
		doDead=self.getProperty("DoDeadTimeCorr").value
		doTF=self.getProperty("MeasureTF20asym").value
		doPrompt=self.getProperty("MeasurePromptPeaks").value
		do2nd=self.getProperty("MeasureDoublePulse").value

		# determine machine/analysis options and column names and build new table header
		tt=WorkspaceFactory.createTable()
		colNames=self.getProperty("NewRuns").value.getColumnNames()
		macNames=[]
		camNames=[]
		tt.addColumn("double",colNames[0],1) # scan variable
		for cc in colNames:
			if(cc.endswith("Run")):
				mac=cc[:-3]
				macNames.append(mac)
				tt.addColumn("int",mac+"Run",6)
				tt.addColumn("double",mac+"Rate",2)
				tt.addColumn("double",mac+"RateErr",5)
				if(doTF):
					tt.addColumn("double",mac+"Alpha",2)
					tt.addColumn("double",mac+"AlphaErr",5)
					tt.addColumn("double",mac+"Asym",2)
					tt.addColumn("double",mac+"AsymErr",5)
				if(doPrompt):
					tt.addColumn("double",mac+"Prompt",2)
					tt.addColumn("double",mac+"PromptErr",5)
				if(do2nd):
					tt.addColumn("double",mac+"DoubleP",2)
					tt.addColumn("double",mac+"DoublePErr",5)
			if(cc.endswith("File")):
				cam=cc[:-4]
				camNames.append(cam)
				tt.addColumn("int",cam+"File",6)
				tt.addColumn("double",cam+"Xpos",2)
				tt.addColumn("double",cam+"XposErr",5)
				tt.addColumn("double",cam+"Ypos",2)
				tt.addColumn("double",cam+"YposErr",5)
				tt.addColumn("double",cam+"Width",2)
				tt.addColumn("double",cam+"WidthErr",5)
				tt.addColumn("double",cam+"Height",2)
				tt.addColumn("double",cam+"HeightErr",5)
				tt.addColumn("double",cam+"Ellip",2)
				tt.addColumn("double",cam+"EllipErr",5)
				tt.addColumn("double",cam+"Intens",2)
				tt.addColumn("double",cam+"IntensErr",5)
				tt.addColumn("double",cam+"Backgd",2)
				tt.addColumn("double",cam+"BackgdErr",5)
		
		#print "grand list of columns:",tt.getColumnNames()
		rowsToDo=[]
		if(not self.getProperty("OldTable").isDefault):
			for r in self.getProperty("OldTable").value:
				rowsToDo.append(r)
		for r in self.getProperty("NewRuns").value:
			rowsToDo.append(r)
		# determine MACHINES from column headings *Run and *File
		#MACHINES=self.getProperty("Machines").value
		
		# determine magnet scanned by column 1 in all tables, actual name just to be copied to output
		#Mag=self.getProperty("Magnet").value
		Mag=colNames[0]

		# clone reference file for output (parent alg already checked consistency)
		tt2=WorkspaceFactory.createTable()
		tt2.addColumn("str","MagnetName",6)
		tt2.addColumn("double","Reference",2)
		tt3=mtd[self.getProperty("NewRuns").valueAsStr+"_Reference"]
		for r in tt3:
			tt2.addRow([r["MagnetName"],r["Reference"]])
		AnalysisDataService.addOrReplace(self.getProperty("OutputTable").valueAsStr+"_Reference",tt2)

		Prog=Progress(self,start=0.0,end=1.0,nreports=len(rowsToDo))
		for r in rowsToDo:
			tr=[r[Mag]]
			for mac in macNames:
				f=DPARS[mac]["data"] % r[mac+"Run"]
				toload=10
				while toload>0:
					try:
						wss=LoadMuonNexus(Filename=f,DeadTimeTable="dttab",DetectorGroupingTable="groups",EnableLogging=False)
						toload=0
					except:
						self.log().notice("file"+str(f)+"isn't in the archive yet, waiting 10 seconds")
						time.sleep(10)
						toload = toload-1
				ws=wss[0]
				frames=wss[0].getRun().getProperty("goodfrm").value
				if(doDead):
					# correct for dead time
					ws=ApplyDeadTimeCorr(InputWorkspace=ws,DeadTimeTable=wss[4],EnableLogging=False)
				# rate by default
				grouped=MuonGroupDetectors(InputWorkspace=ws,DetectorGroupingTable=wss[5],EnableLogging=False)
				if(doTF):
					alpha=AlphaCalc(grouped,FirstGoodValue=0.5,LastGoodValue=10.0,EnableLogging=False)
					self.log().notice("First guess of alpha="+str(alpha))
					asy=AsymmetryCalc(grouped,Alpha=alpha,EnableLogging=False)
					fguess=ws.getRun().getProperty("sample_magn_field").value*0.01355
					fitdat=Fit(Function="name=ExpDecayOsc,A=0.2,Lambda=0.01,Frequency="+str(fguess)+",Phi=0.1;name=FlatBackground,A0=0.0",InputWorkspace=asy,StartX=DPARS[mac]["tgbegin"],EndX=DPARS[mac]["tgend"],Output="Fitted",EnableLogging=False)
					alpha=alpha+(2.0*alpha*fitdat[3].cell("Value",4)) # refined value, redo fit
					self.log().notice("Second guess of alpha="+str(alpha)+" using A0="+str(fitdat[3].cell("Value",4)))
					asy=AsymmetryCalc(grouped,Alpha=alpha,EnableLogging=False)
					fitdat=Fit(Function="name=ExpDecayOsc,A=0.2,Lambda=0.01,Frequency="+str(fguess)+",Phi=0.1;name=FlatBackground,A0=0.0",InputWorkspace=asy,StartX=DPARS[mac]["tgbegin"],EndX=DPARS[mac]["tgend"],Output="Fitted",EnableLogging=False)
					alpha=alpha+(2.0*alpha*fitdat[3].cell("Value",4)) # refined value, redo fit
					self.log().notice("Third guess of alpha="+str(alpha)+" using A0="+str(fitdat[3].cell("Value",4)))
					asy=AsymmetryCalc(grouped,Alpha=alpha,EnableLogging=False)
					fitdat=Fit(Function="name=ExpDecayOsc,A=0.2,Lambda=0.01,Frequency="+str(fguess)+",Phi=0.1;name=FlatBackground,A0=0.0",InputWorkspace=asy,StartX=DPARS[mac]["tgbegin"],EndX=DPARS[mac]["tgend"],Output="Fitted",EnableLogging=False)
					alpha=alpha+(2.0*alpha*fitdat[3].cell("Value",4))
					self.log().notice("Fourth guess of alpha="+str(alpha)+" using A0="+str(fitdat[3].cell("Value",4)))
				else:
					alpha=1.0
				sf=ExtractSingleSpectrum(InputWorkspace=grouped,WorkspaceIndex=0,EnableLogging=False)
				sb=ExtractSingleSpectrum(InputWorkspace=grouped,WorkspaceIndex=1,EnableLogging=False)
				sba=Scale(InputWorkspace=sb,Factor=alpha,Operation="Multiply",EnableLogging=False)
				s=Plus(LHSWorkspace=sf,RHSWorkspace=sba,EnableLogging=False)
				integ=Integration(InputWorkspace=s,RangeLower=DPARS[mac]["tgbegin"],RangeUpper=DPARS[mac]["tgend"],EnableLogging=False)
				rate=integ.dataY(0)[0]/frames
				rateErr=integ.dataE(0)[0]/frames
				tr.extend([r[mac+"Run"],rate,rateErr])
				if(doTF):
					# check that Phi is close to 0, Asym might be negative if muons are polarised the wrong way. Or save phi too?
					mAsym=fitdat[3].cell("Value",0) # from xx_Parameters table
					eAsym=fitdat[3].cell("Error",0)
					eAlpha=fitdat[3].cell("Error",4)*2*alpha
					if(math.cos(fitdat[3].cell("Value",3)) < 0):
						mAsym=-mAsym
					tr.extend([alpha,eAlpha,mAsym,eAsym])
				if(doPrompt):
					prws=Integration(InputWorkspace=ws,RangeLower=DPARS[mac]["promptbegin"],RangeUpper=DPARS[mac]["promptend"],EnableLogging=False)
					tr.extend([prws.dataY(0)[0]/frames,prws.dataE(0)[0]/frames])
					# integrate a suggested region or regions
					# give absolute rate as counts/frame
					pass
				if(do2nd):
					# contamination factor, 0.0=pure correct single pulse, 1.0=equal double pulse, >1 = wrong pulse!
					ecs=RemoveExpDecay(s,EnableLogging=False)
					p1ws=Integration(InputWorkspace=ecs,RangeLower=DPARS[mac]["p1begin"],RangeUpper=DPARS[mac]["p1end"],EnableLogging=False)
					p12ws=Integration(InputWorkspace=ecs,RangeLower=DPARS[mac]["tgbegin"],RangeUpper=DPARS[mac]["tgend"],EnableLogging=False)
					r1=p1ws.dataY(0)[0]/(p1ws.dataX(0)[1]-p1ws.dataX(0)[0])
					r12=p12ws.dataY(0)[0]/(p12ws.dataX(0)[1]-p12ws.dataX(0)[0])
					if(DPARS[mac]["pulse"]==1):
						tr.extend([(r12-r1)/r1,0.0]) # do error analysis!
					elif(DPARS[mac]["pulse"]==2):
						tr.extend([r1/(r12-r1),0.0])
					# integrate region between pulses, avoiding prompt pulses = P1
					# also get final count rate extrapolated back = P12
					# for MuSR: fraction = P1/(P12-P1)
					# for EMU/HIFI: fraction=(P12-P1)/P1
					# if TF data, having got alpha, calc "sum" = F+alpha*B rather than F+B.
			for cam in camNames:
				# load camera file and analyse it, append to tr[]
				parset=CPARS[cam]["processor"](CPARS[cam]["data"] % r[mac+"File"])
				tr.append(r[mac+"File"])
				tr.extend(parset)
			#print "row ready to add:",tr
			tt.addRow(tr)
			Prog.report("Data analysis")
		self.setProperty("OutputTable",tt)
AlgorithmFactory.subscribe(AnalyseTuningScan)
