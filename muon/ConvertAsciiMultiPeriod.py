import time
import datetime
import math
import numpy
from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty
from mantid.kernel import *
from mantid.simpleapi import *

class ConvertAsciiMultiPeriod(PythonAlgorithm):
	# muon data
	# create multi period ascii file
	# col 1=time
	# col 2=F(P1) 3=B(P1) 4=F(P2) 5=B(P2) etc
	# no headers
	# optional dead time
	def PyInit(self):
		self.declareProperty(FileProperty(name="FirstFile",defaultValue="",action=FileAction.Load,extensions = ["nxs"]))
		self.declareProperty(FileProperty(name="LastFile",defaultValue="",action=FileAction.OptionalLoad,extensions = ["nxs"]))
		self.declareProperty("CorrectDeadTime",True,doc="Apply dead time corrections")
		self.declareProperty("Grouped",True,doc="Group data")
		self.declareProperty(FileProperty(name="OutputDirectory",defaultValue="",action=FileAction.Directory),doc="Files put here, named to match source files with extension .dat")

	def category(self):
		return "Muon"
		
	def PyExec(self):
		
		file1=self.getProperty("FirstFile").value
		file9=self.getProperty("LastFile").value
		doDeadTimes=self.getProperty("CorrectDeadTime").value
		group=self.getProperty("Grouped").value
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
		# analyse output file
		fileout1=self.getProperty("OutputDirectory").value
		prog_reporter = Progress(self,start=0.0,end=1.0,nreports=lastnum+1-firstnum)
		for ff in range(firstnum,lastnum+1):
			if(file1 != file9):
				thispath=file1[:j1]+str(ff).zfill(i1-j1)+file1[i1:]
			else:
				thispath=file1
			# cope with Periods - have Workspace Group and extra return values from LoadMuonNexus
			wstuple=LoadMuonNexus(filename=thispath,OutputWorkspace="__ConvAsciiTmp",DeadTimeTable="__ConvAsciiDead",DetectorGroupingTable="__ConvAsciiGroup")
			try:
				ws=wstuple[0] # first period if there are some
				Deads=wstuple[4]
				Groups=wstuple[5]
				L=len(ws)
			except:
				ws=(wstuple[0],) # plain run
				Deads=(wstuple[4],)
				Groups=(wstuple[5],)
			prog_reporter.report("Processing")
			nper=len(ws)
			if (group):
				fb=MuonGroupDetectors(InputWorkspace=ws[0],DetectorGroupingTable=Groups[i])
				nh=fb.getNumberHistograms()
			else:
				nh=ws[0].getNumberHistograms()
			outArray=numpy.zeros([len(ws[0].dataY(0)),nper*nh*2+1])
			outArray[:,0]=(ws[0].dataX(0)[:-1]+ws[0].dataX(0)[1:])/2
			for i in range(nper):
				if(doDeadTimes):
					dtcdat=ApplyDeadTimeCorr(ws[i],Deads[i])
				else:
					dtcdat=ws[i]
				if group:
					fb=MuonGroupDetectors(InputWorkspace=dtcdat,DetectorGroupingTable=Groups[i])
				else:
					fb=dtcdat
				for j in range(nh):
					outArray[:,i*nh*2+j*2+1]=fb.dataY(j)
					outArray[:,i*nh*2+j*2+2]=fb.dataE(j)
				#outArray[:,i*nh*2+3]=fb.dataY(1)
				#outArray[:,i*nh*2+4]=fb.dataE(1)
			outpath=os.path.join(fileout1,os.path.splitext(os.path.split(thispath)[1])[0]+".dat")
			fmtlist=['%9.1f']*(nper*2*nh)
			fmtlist.insert(0,'%6.3f')
			numpy.savetxt(outpath,outArray,fmt=fmtlist)
			for w in ws:
				DeleteWorkspace(w)
			for g in Groups:
				DeleteWorkspace(g)
			for d in Deads:
				DeleteWorkspace(d)
			# print "finished file ",thispath

AlgorithmFactory.subscribe(ConvertAsciiMultiPeriod)
