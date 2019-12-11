from __future__ import print_function
import time
import datetime
import math
import numpy
from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty
from mantid.kernel import *
from mantid.simpleapi import *

class CopyLogsToHist(PythonAlgorithm):
	# hard coded for HiFI data
	# copy the log values "DetectorTemp1...DetectorTemp63" to the corresponding histograms, in place of the data
	# time values are now hours relative to midnight on the day the first run started.
	# Due to SECI the times are not necessarily equal from one spectrum to another
	# logs truncated to the length of the shortest one (so MatrixWorkspace is not too ragged)
	# option to interpolate the log to equal bin sizes
	# batch processing of a block of consecutive runs
	# last run can be numbered 0, or beyond the end of the saved runs: will go to end of saved runs and also load the latest temporary file and then stop there
	# use first run and last run the same to load any one file of any name (even not numbered)
	# non-logged detectors are interpolated between their neighbours (times too, if not interpolating time)
	# for example use Instrument View or 2D Colour Plot on the output workspace
	def PyInit(self):
		#self.declareProperty(WorkspaceProperty("InputWS","",Direction.Input),"Source run")
		self.declareProperty(FileProperty(name="FirstFile",defaultValue="",action=FileAction.Load,extensions = ["nxs"]))
		self.declareProperty(FileProperty(name="LastFile",defaultValue="",action=FileAction.OptionalLoad,extensions = ["nxs"]))
		self.declareProperty(WorkspaceProperty("OutputWS","",Direction.Output),"With histograms filled with logs")
		self.declareProperty("Relative",0,doc="+1 for first point as reference, -1 for last, 2 for mean, 3 to remove average time dependence,4=take out scatter from 1st")
		self.declareProperty("InterpolatePts",-1,doc="Number of time bins to interpolate into")
		self.declareProperty("InterpolateDelta",-1.0,doc="Width (s) of time bins to interpolate into")
		self.declareProperty("OmitUnmonitoredSpectra",False,doc="Remove spectra without their own thermometers instead of interpolating from neighbours")
		self.declareProperty("RunningOnly",False,doc="Remove log points before run starts")
		self.declareProperty("AutoConvKtoC",False,doc="Auto conversion K to C")
		self.declareProperty(StringArrayProperty(name="ExtraLogs",direction=Direction.Input,values=""))

	def category(self):
		return "Muon"
		
	def PyExec(self):
		
		# "Det_Water","Temp_Air1","Temp_Air2"
		LogNames=["DetectorTemp1","DetectorTemp3","DetectorTemp5","DetectorTemp7","DetectorTemp9","DetectorTemp11","DetectorTemp13","DetectorTemp15","DetectorTemp17","DetectorTemp19","DetectorTemp21","DetectorTemp23","DetectorTemp25","DetectorTemp27","DetectorTemp29","DetectorTemp31","DetectorTemp33","DetectorTemp35","DetectorTemp37","DetectorTemp39","DetectorTemp41","DetectorTemp43","DetectorTemp45","DetectorTemp47","DetectorTemp49","DetectorTemp51","DetectorTemp53","DetectorTemp55","DetectorTemp57","DetectorTemp59","DetectorTemp61","DetectorTemp63"]
		ExtraLogs=self.getProperty("ExtraLogs").value
		RunningOnly=self.getProperty("RunningOnly").value
		LogNames.extend(ExtraLogs)
		print("logs to be got: ",LogNames)
		#ws=self.getProperty("InputWS").value
		file1=self.getProperty("FirstFile").value
		file9=self.getProperty("LastFile").value
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
		times=[0]*len(LogNames)
		temps=[0]*len(LogNames)
		for i in range(len(LogNames)):
			times[i]=numpy.zeros(0)
			temps[i]=numpy.zeros(0)
		# loop and load files. Absolute numbers for now.
		gotToCurrentRun=False
		prog_reporter = Progress(self,start=0.0,end=1.0,nreports=lastnum+1-firstnum)
		for ff in range(firstnum,lastnum+1):
			if(file1 != file9):
				thispath=file1[:j1]+str(ff).zfill(i1-j1)+file1[i1:]
			else:
				thispath=file1
			# cope with Periods - have Workspace Group and extra return values from LoadMuonNexus
			try:
				wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__CopyLogsTmp")
			except:
				thispath=file1[:j1]+"auto_A.tmp" # should really try all of them and pick newest!
				wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__CopyLogsTmp")
				gotToCurrentRun=True
			try:
				ws=wstuple[0][0] # first period if there are some
			except:
				ws=wstuple[0] # plain run
			prog_reporter.report("Loading")
			if(gotToCurrentRun and int(ws.getRun().getProperty("run_number").value) != ff):
				print("temporary file is for the previous run, probably no new one started")
				break
			if(ff==firstnum):
				#begin=datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value,"%Y-%m-%dT%H:%M:%S")[0:6]))
				begin=datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value.strip(),"%Y-%m-%dT%H:%M:%S")[0:3])) # start of day
			for ss in range(len(LogNames)):
				try:
					#log=ws.getRun().getLogData("DetectorTemp"+str(2*ss+1))
					log=ws.getRun().getLogData(LogNames[ss])
					nnp=log.size()
					tconv=numpy.zeros(nnp)
					for tt in range(nnp):
						tconv[tt]=(datetime.datetime(*(time.strptime(str(log.times[tt]).strip(),"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds()
					if(RunningOnly):
						thisbegin=(datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value.strip(),"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds()
						ri=numpy.searchsorted(tconv,thisbegin,side='right')+1 # throw away any point at t==0.0 and one more
					else:
						ri=0
					times[ss]=numpy.append(times[ss],tconv[ri:])
					temps[ss]=numpy.append(temps[ss],log.value[ri:])
				except Exception as ee:
					print(ee)
					print("no log in this run for ",LogNames[ss]," inserting gap")
					# one NaN point at the start of this run
					tconv=numpy.array([(datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value.strip(),"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds()])
					times[ss]=numpy.append(times[ss],tconv)
					temps[ss]=numpy.append(temps[ss],numpy.nan)
			print("finished file ",thispath)
			if(gotToCurrentRun):
				break
		autoKC=self.getProperty("AutoConvKtoC").value
		if(autoKC):
			# K to C conversion (approximate): values >200 are K. Outputs all in C.
			for ss in range(len(LogNames)):
				if(temps[ss][-1]>200):
					temps[ss]=temps[ss]-273.15
		relflag=self.getProperty("Relative").value
		unmon=self.getProperty("OmitUnmonitoredSpectra").value
		if(unmon):
			skip=1
		else:
			skip=2
		print("time zero for result is ", begin)
		n=1
		t0m= 100000000.0
		t1m=-100000000.0
		for d in range(len(LogNames)):
			n=max(n,len(times[d]))
			t0m=min(t0m,times[d][0])
			t1m=max(t1m,times[d][-1])
			print("detector ",2*d+1," has ",len(times[d])," time stamps and ",len(temps[d])," readings")
			if(len(times[d]) != len(temps[d]) ):
				raise Exception("Something went wrong reading logs!")
		print("longest log has ",n," points")
		print("logs span ",(t1m-t0m)/3600.0," hours run time")
		interpflag = False
		ipts=self.getProperty("InterpolatePts").value
		if(ipts>1):
			interpflag=True
			interpdt=numpy.linspace(t0m/3600.0,t1m/3600.0,num=ipts)
			n=len(interpdt)
			print("interpolating into bins of width ",(interpdt[1]-interpdt[0])*3600.0," seconds")
		else:
			idelta=self.getProperty("InterpolateDelta").value
			if(idelta>0):
				interpflag=True
				t0m=math.ceil(t0m/idelta)*idelta # points are on round minute/hour boundaries where appropriate
				interpdt=numpy.arange(t0m/3600.0,t1m/3600.0,idelta/3600.0)		
				n=len(interpdt)
				print("interpolating into ",n," bins")
		ows=WorkspaceFactory.create(ws,len(LogNames)*skip,n,n)
		ows.setYUnitLabel("Temperature / K")
		# ows.getAxis(0).setUnit("hours") would be nice...
		for d in range(0,len(LogNames)):
			#log=ws.getRun().getLogData("DetectorTemp"+str(d+1))
			if(interpflag):
				nthis=len(times[d])
			else:
				nthis=n
			newX=numpy.zeros(nthis)
			newY=numpy.zeros(nthis)
			newE=numpy.zeros(n)
			reference=0
			if(relflag == 1):
				reference=temps[d][0]
			if(relflag == -1):
				reference=temps[d][-1]
			if(relflag == 2 or relflag == 3):
				reference=numpy.mean(temps[d][0:nthis])
			if(relflag==4):
				reference=temps[d][0]-numpy.mean(temps[:][0])
			for i in range(len(times[d])):
				newX[i]=times[d][i]/3600.0 # hours
				try:
					newY[i]=float(temps[d][i])-reference
				except:
					if (temps[d][i]=="ON"): # logical values
						newY[i]=1.0
					elif (temps[d][i]=="OFF"):
						newY[i]=0.0
					else:
						newY[i]=0-01.0
			for i in range(len(times[d]),nthis):
				newX[i]=times[d][len(times[d])-1]/3600.0 # hours
				newY[i]=newY[len(times[d])-1]
			dxp=ows.dataX(d*skip)
			dyp=ows.dataY(d*skip)
			dep=ows.dataE(d*skip)
			if(interpflag):
				dxp[...]=interpdt
				dyp[...]=numpy.interp(interpdt,newX,newY)
				dep[...]=newE
			else:
				dxp[...]=newX
				dyp[...]=newY
				dep[...]=newE
		if(relflag==3 and interpflag):
			refcurve=numpy.zeros(n)
			for d in range(0,len(LogNames)):
				dyp=ows.dataY(d*skip)
				refcurve+=dyp
			refcurve /= (len(LogNames)*1.0)
			#refcurve -= numpy.mean(refcurve)
			for d in range(0,len(LogNames)):
				dyp=ows.dataY(d*skip)
				dyp[...]=dyp[...]-refcurve
			#print "hist ",d," has ",len(log.value)," points with ",dxp[n-2]-dxp[1]," span"
		if(unmon):
			# spare spectra never existed but renumber the rest
			for d in range(0,len(LogNames)):
				#ows.getAxis(1).setValue(d,d*2+1) # obsolete non working method
				ows.getSpectrum(d).setSpectrumNo(d*2+1)
		else:
			for d in range(1,len(LogNames)*2,2):
				d1=d-1
				d2=d+1
				if(d==31):
					d2=0
				if(d==63):
					d2=32
				if(d>63):
					d1=d-1
					d2=d-1 # duplicate "extras", don't interpolate
				ows.dataX(d)[...]=(ows.dataX(d1)+ows.dataX(d2))/2
				ows.dataY(d)[...]=(ows.dataY(d1)+ows.dataY(d2))/2
				ows.dataE(d)[...]=(ows.dataE(d1)+ows.dataE(d2))/2
		
		self.setProperty("OutputWS",ows)
		DeleteWorkspace(ws)

AlgorithmFactory.subscribe(CopyLogsToHist)
