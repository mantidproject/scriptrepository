""" SplitLaserPowerLog
Author: James Lord Aug 2017
Take a continuous time series log in separate text file and extract the values corresponding to a range of runs
Build a table of run number, value, rms variation.
"""
import datetime
import numpy
# read Laser Power log and slice up to match runs
class SplitLaserPowerLog(PythonAlgorithm):

	def category(self):
		return 'Muon;DataHandling'

	def PyInit(self):
		self.declareProperty(FileProperty(name="FirstMuonFile",defaultValue="",action=FileAction.Load,extensions = ["nxs"]))
		self.declareProperty(FileProperty(name="LastMuonFile",defaultValue="",action=FileAction.Load,extensions = ["nxs"]))
		self.declareProperty(FileProperty(name="LaserPowerLog",defaultValue="",action=FileAction.Load,extensions = ["txt"]))
		self.declareProperty(ITableWorkspaceProperty('Results','',Direction.Output),"Mean power versus run number")

	def PyExec(self):
		# create output table
		oTab=WorkspaceFactory.createTable()
		oTab.addColumn("int","Run",1)
		oTab.addColumn("float","Power",2)
		oTab.addColumn("float","Error",5)
		
		# open and process power log
		pl=self.getProperty("LaserPowerLog").value
		doneHeader=False
		refTime=None
		timestamps=[]
		powers=[]
		with open(pl,"r") as file:
			for line in file:
				if (doneHeader):
					(t,p)=list(map(float,line.split()))
					timestamps.append(t)
					powers.append(p)
				else:
					if(";Logged:" in line):
						#print "extracting time from <"+line+">"
						refTime=datetime.datetime.strptime(line.strip(),";Logged:%d/%m/%Y at %H:%M:%S")
					if("Channel A" in line) and not (";" in line):
						doneHeader=True
		N=len(timestamps)
		timestamps=numpy.array(timestamps)
		powers=numpy.array(powers)
		
		#print "timestamps go from ",timestamps[0],"to",timestamps[-1]
		
		# load Muon files and get dates. Add lines to table
		file1=self.getProperty("FirstMuonFile").value
		file9=self.getProperty("LastMuonFile").value
		if(file1 != file9):
			i1=file1.rindex('.')
			j1=i1-1
			while file1[j1-1].isdigit():
				j1=j1-1
			firstnum=int(file1[j1:i1])
			i9=file9.rindex('.')
			j9=i9-1
			while file9[j9-1].isdigit():
				j9=j9-1
			lastnum=int(file9[j9:i9])
			if(file1[:j9] != file9[:j9]):
				raise Exception("Files from different directories or instruments")
			if(file1[i1:] != file9[i9:]):
				raise Exception("Files of different types")
			if(i1-j1 != i9-j9):
				raise Exception("File numbering error")
			if(lastnum < firstnum):
				if(firstnum != 0):
					raise Exception("Run numbers must increase")
				else:
					# run zero: load temporary file instead
					lastnum=10^(i1-j1)-1 # highest possible number as place holder
			#		runpaths.append(firstfn[:j]+str(n+firstnum).zfill(i-j)+firstfn[i:])
		else: # load any single file even if it's not numbered, eg Temporary File.
			firstnum=0
			lastnum=0
		# arrays
		# loop and load files. Absolute numbers for now.
		gotToCurrentRun=False
		prog_reporter = Progress(self,start=0.0,end=1.0,nreports=lastnum+1-firstnum)

		for ff in range(firstnum,lastnum+1):
			if(file1 != file9):
				thispath=file1[:j1]+str(ff).zfill(i1-j1)+file1[i1:]
			else:
				thispath=file1
			# cope with Periods - have Workspace Group and extra return values from LoadMuonNexus
			#print "loading from ",thispath
			try:
				wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__CopyLogsTmp")
			except:
				thispath=file1[:j1]+"auto_A.tmp" # should really try all of them and pick newest!
				wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__CopyLogsTmp")
				gotToCurrentRun=True
			#print len(wstuple)
			ws=wstuple[0] # plain run
			if isinstance(ws,mantid.api.WorkspaceGroup):
				ws=ws[0] # first period, should overlap with same log data in each
			prog_reporter.report("Loading")
			#print "run has times ",ws.getRun().getProperty("run_start").value," - ",ws.getRun().getProperty("run_end").value
			thisstart=datetime.datetime.strptime(ws.getRun().getProperty("run_start").value,"%Y-%m-%dT%H:%M:%S")
			thisend=datetime.datetime.strptime(ws.getRun().getProperty("run_end").value,"%Y-%m-%dT%H:%M:%S")
			st1=(thisstart-refTime).total_seconds()
			st2=(thisend-refTime).total_seconds()
			# look up times
			it1,it2=numpy.searchsorted(timestamps,(st1,st2))
			#print "run:",ff," slice ",it1,"to",it2," out of ",N
			if(it2==0 or it1==N):
				raise IndexError("Run "+str(ff)+" is outside the power log time span")
			if(it1==it2):
				raise KeyError("Run "+str(ff)+" was too quick to get a log point")
			if(it1==0 or it2==N):
				self.logger().warn("Run "+str(ff)+" is at the start or end, are points missing?")
			meanPower=numpy.mean(powers[it1:it2])
			if(it1+1==it2):
				self.logger().warn("Run "+str(ff)+" has only one log point, no statistics possible")
				rmsPower=0.0
			else:
				rmsPower=numpy.std(powers[it1:it2])
			oTab.addRow((ff,meanPower,rmsPower))
			
		self.setProperty("Results",oTab)
		
AlgorithmFactory.subscribe(SplitLaserPowerLog)