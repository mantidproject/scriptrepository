from mantid.kernel import *
from mantid.simpleapi import *
#import math
 
class Create2DALCMap(PythonAlgorithm):
 
	def PyInit(self):
		# Declare properties
		# first and last runs to be loaded as for PlotAsymmetryByLogValue
		self.declareProperty(FileProperty("FirstFile","",FileAction.Load, extensions = ["nxs"]),doc="First of a sequence of runs measured while scanning field or something")
		self.declareProperty("NumberOfRuns",1000000,doc="Number of runs, default continues until a log value is out of sequence.")
		# otherwise use this many runs and sort them into increasing log value order so will merge interleaved scans
		self.declareProperty("StartTime",20.0,doc="Start time for result map")
		self.declareProperty("EndTime",20.0,doc="End time for result map")
		self.declareProperty("NBunch",10,doc="Bunching factor to improve stats")
		self.declareProperty("LogSelection","Field_Main_Target",doc="Which log value to use for Y axis")
		self.declareProperty(WorkspaceProperty("OutputWorkspace","",Direction.Output),doc="Results!")

	def category(self):
		return 'Muon'

	def PyExec(self):
		# Run the algorithm
	
		wsOutput = self.getPropertyValue("OutputWorkspace")
		
		firstfn = self.getProperty('FirstFile').value
		nmax=self.getProperty("NumberOfRuns").value
		StartTime=self.getProperty("StartTime").value
		EndTime=self.getProperty("EndTime").value
		NBunch=self.getProperty("NBunch").value
		logname=self.getProperty("LogSelection").value
		
		progRep=Progress(self,start=0.0,end=1.0,nreports=nmax)
		
		i=firstfn.rindex('.')
		j=i-1
		while firstfn[j-1].isdigit():
			j=j-1
		firstnum=int(firstfn[j:i])
		k=j-1
		while firstfn[k-1].isalpha():
			k=k-1

		filelist=[]
		vi=0
		periods=False
		oldlogval=-1000000.0
		while(vi<nmax and oldlogval<999999):
		  #runlist.append(vi+firstnum)
		  thisrunname=firstfn[k:j]+str(vi+firstnum).zfill(i-j)
		  thisrunpath=firstfn[:j]+str(vi+firstnum).zfill(i-j)+firstfn[i:]
		  try:
			  ws0=LoadMuonNexus(Filename=thisrunpath,OutputWorkspace='__ALCmap_raw')
			  ws=ws0[0]
			  try:
				  try:
					ylogval=ws[0].getRun().getProperty(logname).value[-1]
				  except:
					ylogval=float(ws[0].getRun().getProperty(logname).value)
				  if(len(filelist)==0):
					  periods=len(ws)
			  except:
				  try:
					ylogval=ws.getRun().getProperty(logname).value[-1]
				  except:
					ylogval=float(ws.getRun().getProperty(logname).value)			  
			  if(ylogval<=oldlogval and nmax>=1000000):
				 print "Out of sequence run found, ending"
				 oldlogval=10000001.0 # end now
			  else:
				  self.log().information(logname+"="+str(ylogval))
				  wsname1='__ALCmap'+str(vi)
				  Rebunch(InputWorkspace='__ALCmap_raw',OutputWorkspace='__ALCmap_fine',NBunch=NBunch)
				  AsymmetryCalc(InputWorkspace='__ALCmap_fine',OutputWorkspace='__ALCmap_asym',ForwardSpectra='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32',BackwardSpectra='33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64',Alpha='1.0')
				  ws=CropWorkspace(InputWorkspace='__ALCmap_asym',OutputWorkspace=wsname1,XMin=StartTime,XMax=EndTime)
				  print logname+"="+str(ylogval)+" increment="+str(ylogval-oldlogval)
				  oldlogval=ylogval
				  filelist.append((ylogval,wsname1))
				  vi=vi+1
				  DeleteWorkspace("__ALCmap_raw")
				  DeleteWorkspace("__ALCmap_fine")
				  DeleteWorkspace("__ALCmap_asym")
		  except Exception as e:
			  print e
			  print e.__repr__()
			  print "run out of data"
			  oldlogval=10000001.0 # end now, run maybe not found?
		  progRep.report()
		# assemble conjoined run
		filelist.sort()
		
		# the axis
		if(periods):
			na=[None]*periods
			for i in range(periods):
				na[i]=NumericAxis.create(len(filelist))
		else:
			na = NumericAxis.create(len(filelist))

		for loopIndex in range(len(filelist)):
			(yval,wsn)=filelist[loopIndex]
			if(periods):
				for i in range(periods):
					na[i].setValue(loopIndex,yval)
			else:
				na.setValue(loopIndex,yval)
			if (loopIndex>0):
				ConjoinWorkspaces(InputWorkspace1=wsOutput,InputWorkspace2=wsn,CheckOverlapping=False)
			else:
				RenameWorkspace(InputWorkspace=wsn,OutputWorkspace=wsOutput)
			if wsn in mtd:
				DeleteWorkspace(Workspace=wsn)

		# na.setUnit("TOF")

		wsOut = mtd[wsOutput]
		#replace the spectrun axis
		if(periods):
			for i in range(periods):
				wsOut[i].replaceAxis(1,na[i])
		else:
			wsOut.replaceAxis(1,na)


		self.setProperty("OutputWorkspace",wsOut)

AlgorithmFactory.subscribe(Create2DALCMap)
