"""*WIKI* 

This algorithm calculates integrated counts for each spectrum (normalised per frame) from a range of files to be loaded and places into a single workspace. The "time" axis is now in clock time (hours from 00:00 on the day the first of the runs began) with the time between runs considered as belonging to the run just ended, apart from the last one.

*WIKI*"""
from __future__ import (absolute_import, division, print_function, unicode_literals)

import time
import datetime
import math
import numpy
from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty
from mantid.kernel import *
from mantid.simpleapi import *

class DetectorRateTimeSeries(PythonAlgorithm):
	# integrate each spectrum and calcuate counts per frame
	# assemble into complete time series using run_start times (and run_end of the last run)
	# time values are now hours relative to midnight on the day the first run started.
	# batch processing of a block of consecutive runs
	# last run can be numbered 0, or beyond the end of the saved runs: will go to end of saved runs and also load the latest temporary file and then stop there
	# for example use Instrument View or 2D Colour Plot on the output workspace
	# integrate whole time range for dead or unstable channels
	# integrate late times only to find noisy channels
	# uses LoadMuonNexus but in principle will work for any type of data
	def PyInit(self):
		#self.declareProperty(WorkspaceProperty("InputWS","",Direction.Input),"Source run")
		self.declareProperty(FileProperty(name="FirstFile",defaultValue="",action=FileAction.Load,extensions = ["nxs"]))
		self.declareProperty(FileProperty(name="LastFile",defaultValue="",action=FileAction.OptionalLoad,extensions = ["nxs"]))
		self.declareProperty(WorkspaceProperty("OutputWS","",Direction.Output),"Summary rates versus time")
		self.declareProperty("IntegrateTimeFrom",0.0)
		self.declareProperty("IntegrateTimeTo",32.0)		
		self.declareProperty("NormaliseByCurrent",False)

	def category(self):
		return "Muon"
		
	def PyExec(self):
		
		#ws=self.getProperty("InputWS").value
		file1=self.getProperty("FirstFile").value
		file9=self.getProperty("LastFile").value
		tfrom=self.getProperty("IntegrateTimeFrom").value
		tto=self.getProperty("IntegrateTimeTo").value
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
			try:
				wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__CopyLogsTmp")
			except:
				thispath=file1[:j1]+"auto_A.tmp" # should really try all of them and pick newest!
				wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__CopyLogsTmp")
				gotToCurrentRun=True
			# print( len(wstuple))
			ws=wstuple[0] # plain run
			if isinstance(ws,mantid.api.WorkspaceGroup):
				ws=ws[0] # first period, should overlap with same log data in each
			prog_reporter.report("Loading")
			if(gotToCurrentRun and int(ws.getRun().getProperty("run_number").value) != ff):
				print("temporary file is for the previous run, probably no new one started")
				break
			if(ff==firstnum):
				times=numpy.zeros(1)
				temps=[0]*ws.getNumberHistograms()
				for i in range(ws.getNumberHistograms()):
					temps[i]=numpy.zeros(0)
				#begin=datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value,"%Y-%m-%dT%H:%M:%S")[0:6]))
				begin=datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value,"%Y-%m-%dT%H:%M:%S")[0:3])) # start of day
				nbc=self.getProperty("NormaliseByCurrent").value
				if(nbc):
					nbclog=None
					for nbclog0 in ["Beam_Current","Beamlog_TS1_Current","beamlog_current"]: # consistency would be nice!
						if(ws.getRun().hasProperty(nbclog0)):
							nbclog=nbclog0
					if(nbclog is None):
						raise Exception("Beam current doesn't seem to be logged on this instrument under any of the names I know")
			ws2=Integration(ws,tfrom,tto)
			if(nbc):
				bclog=ws.getRun().getProperty(nbclog).value
				bcon=numpy.extract(bclog>100.0,bclog) # exclude beam-off and base rate
				bc=numpy.mean(bcon) # accidentally includes beam current before run
			else:
				bc=1.0
			for ss in range(ws2.getNumberHistograms()):
				temps[ss]=numpy.append(temps[ss],ws2.readY(ss)[0]/ws.getRun().getProperty("goodfrm").value/bc)
			times[-1]=((datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_start").value,"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds())/3600
			times=numpy.append(times,((datetime.datetime(*(time.strptime(ws.getRun().getProperty("run_end").value,"%Y-%m-%dT%H:%M:%S")[0:6]))-begin).total_seconds())/3600)
			print("finished file ",thispath)
			DeleteWorkspace(ws2)
			if(gotToCurrentRun):
				break
		print("time zero for result is ", begin)
		print("size might be X=",len(times)," and Y=",len(temps[0]))
		ows=WorkspaceFactory.create(ws,NVectors=ws.getNumberHistograms(),XLength=len(times),YLength=len(temps[0]))
		ows.setYUnitLabel("Counts per Frame")
		# ows.getAxis(0).setUnit("hours") would be nice...
		newE=numpy.zeros(len(temps[0]))
		for d in range(0,ws.getNumberHistograms()):
			ows.dataX(d)[...]=times
			ows.dataY(d)[...]=temps[d]
			ows.dataE(d)[...]=newE
		
		self.setProperty("OutputWS",ows)
		DeleteWorkspace(ws)

AlgorithmFactory.subscribe(DetectorRateTimeSeries)
