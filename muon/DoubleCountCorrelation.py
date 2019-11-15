""" For Dark Counts and other sparse runs 
Generate a NSpec*NSpec workspace with 1 added to every bin where there is a non-zero number of counts in both of those spectra in the same time bin
now just multiples bin contents, same for really sparse workspaces but meaningful if more counts present
Filter on time range if requested (time applies to S1 only)
Optional offset to bin in S2 for delayed correlation
Requires time bins to match across spectra (but doesn't check)
Upgrade to Event Mode would be ideal here!
if given Workspace Group then process all workspaces individually and sum results? or use MergeRuns on output WorkspaceGroup"""
from __future__ import print_function
import numpy

from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty

class DoubleCountCorrelation(PythonAlgorithm):

	def PyInit(self):
		self.declareProperty(WorkspaceProperty("InputWorkspace","",Direction.Input))
		self.declareProperty("StartTimeWindow",-1000000.0,Direction.Input)
		self.declareProperty("EndTimeWindow",1000000.0,Direction.Input)
		self.declareProperty("BinOffset",0,Direction.Input)
		self.declareProperty("SuppressSelfCorrelation",True,Direction.Input)
		self.declareProperty(WorkspaceProperty("OutputWorkspace","",Direction.Output))

	def category(self):
		return "Diagnostics;Muon"

	def PyExec(self):
		inWS = self.getProperty("InputWorkspace").value
		t1 = self.getProperty("StartTimeWindow").value
		t2 = self.getProperty("EndTimeWindow").value
		dt = self.getProperty("BinOffset").value
		ssc = self.getProperty("SuppressSelfCorrelation").value
		Nspec=inWS.getNumberHistograms()
		Xlen=len(inWS.readX(0))
		Ylen=len(inWS.readY(0))
		print("lengths ",Xlen,Ylen)
		
		ows=WorkspaceFactory.create("Workspace2D",Nspec,Nspec+1,Nspec)
		for j in range(Nspec):
			ows.dataX(j)[:]=numpy.linspace(0.5,0.5+Nspec,Nspec+1)
		#ows.replaceAxis(0,inWS.getAxis(1))
		#ows.replaceAxis(1,inWS.getAxis(1))
		ows.setYUnitLabel("Double Counts")
		print(inWS.dataX(0))
		i1=numpy.searchsorted(inWS.dataX(0),t1,side='left')
		i2=numpy.searchsorted(inWS.dataX(0),t2,side='right')-1
		print("requested bins [",i1,":",i2,"]")
		if(i1+dt < 0):
			i1=-dt
		if(i2+dt > Ylen):
			i2=Ylen-dt
		print("processing bins [",i1,":",i2,"]")
		if(i1<0 or i2<=i1 or i2>Ylen):
			raise Exception("Oops, bin calculation mistake")
		for s1 in range(Nspec):
			for s2 in range(Nspec):
				if(s1 == s2 and dt==0 and ssc):
					ows.dataY(s1)[s2]=0 # self correlation not wanted
				else:
					ows.dataY(s1)[s2]=numpy.sum(inWS.dataY(s1)[i1:i2] * inWS.dataY(s2)[i1+dt:i2+dt])

		self.setProperty("OutputWorkspace",ows)

AlgorithmFactory.subscribe(DoubleCountCorrelation)
