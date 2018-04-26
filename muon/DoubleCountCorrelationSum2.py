import numpy
""" For Dark Counts and other sparse runs 
Correlates spectra looking for counts in S1 at time T1 and S2 at T2+dt
now just multiples bin contents, same for really sparse workspaces but meaningful if more counts present
Filter on time range if requested (time applies to S1 only)
Scan dT and generate auto correlation time plot
Spectra are: 0=self, 1=nnL, 2=nnR, 3=nnnL, 4=nnnR, 5=elsewhere same bank, 6=opposite bank F->B, 7=opposite B->F, 8=total other
Requires time bins to match across spectra (but doesn't check)
Upgrade to Event Mode would be ideal here! """

from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty

class DoubleCountCorrelationSum(PythonAlgorithm):

	def PyInit(self):
		self.declareProperty(WorkspaceProperty("InputWorkspace","",Direction.Input))
		self.declareProperty("StartTimeWindow",-1000000.0,Direction.Input)
		self.declareProperty("EndTimeWindow",1000000.0,Direction.Input)
		self.declareProperty("MaxBinOffset",10,IntBoundedValidator(lower=1,upper=16384))
		self.declareProperty(WorkspaceProperty("OutputWorkspace","",Direction.Output))

	def category(self):
		return "Diagnostics;Muon"

	def PyExec(self):
		inWS = self.getProperty("InputWorkspace").value
		t1 = self.getProperty("StartTimeWindow").value
		t2 = self.getProperty("EndTimeWindow").value
		dtm = self.getProperty("MaxBinOffset").value
		Nspec=inWS.getNumberHistograms()
		Xlen=len(inWS.readX(0))
		Ylen=len(inWS.readY(0))
		
		ows=WorkspaceFactory.create("Workspace2D",9,2*dtm+1,2*dtm+1)
		axis0 = ows.getAxis(0)
		unit = axis0.setUnit("Label")
		unit.setLabel("Time","microsecond")

		for j in range(9):
			ows.dataX(j)[:]=numpy.linspace(-dtm*(inWS.dataX(0)[1]-inWS.dataX(0)[0]),dtm*(inWS.dataX(0)[1]-inWS.dataX(0)[0]),2*dtm+1) # also assumes uniform bin sizes!
			ows.dataY(j)[:]=numpy.zeros(2*dtm+1)
			ows.dataE(j)[:]=numpy.zeros(2*dtm+1)

		ows.setYUnitLabel("Double Counts")
		axis1=TextAxis.create(9)
		axis1.setLabel(0,'self')
		axis1.setLabel(1,'neighbour (cw)')
		axis1.setLabel(2,'neighbour (ac)')
		axis1.setLabel(3,'next nearest (cw)')
		axis1.setLabel(4,'next nearest (ac)')
		axis1.setLabel(5,'same bank')
		axis1.setLabel(6,'other bank F->B')
		axis1.setLabel(7,'other bank F->B')
		axis1.setLabel(8,'total others')
		ows.replaceAxis(1,axis1)

		i1=numpy.searchsorted(inWS.dataX(0),t1,side='left')
		i2=numpy.searchsorted(inWS.dataX(0),t2,side='right')-1
		if(i1-dtm < 0):
			i1=-dtm
		if(i2+dtm > Ylen):
			i2=Ylen-dtm
		print "processing bins [",i1,":",i2,"]"
		if(i1<0 or i2<=i1 or i2>Ylen):
			raise Exception("Oops, bin calculation mistake")
		half=Nspec/2 # assume 2 banks
		for j,dt in enumerate(range(-dtm,dtm+1)):
			for s1 in range(Nspec):
				for s2 in range(Nspec):
					if(s1/half == s2/half):
						if(s1==s2):
							ss=0
						elif((s1%half - s2%half == 1 and s1<half) or (s1%half - s2%half == -1 and s1>=half)):
							ss=1 # nn, cw?
						elif((s1%half - s2%half == 1 and s1>=half) or (s1%half - s2%half == -1 and s1<half)):
							ss=2 # nn, acw?
						elif((s1%half - s2%half == 2 and s1<half) or (s1%half - s2%half == -2 and s1>=half)):
							ss=3 # nnn, cw
						elif((s1%half - s2%half == 2 and s1>=half) or (s1%half - s2%half == -2 and s1<half)):
							ss=4 # nnn, acw
						else:
							ss=5 # elsewhere same bank
					else:
						if(s1<s2):
							ss=6 # other bank, F then B
						else:
							ss=7 # other bank, B then F
					cc=numpy.sum(inWS.dataY(s1)[i1:i2] * inWS.dataY(s2)[i1+dt:i2+dt])
					ows.dataY(ss)[j]=ows.dataY(ss)[j]+cc
					if(ss>0):
						ows.dataY(8)[j]=ows.dataY(8)[j]+cc

		self.setProperty("OutputWorkspace",ows)

AlgorithmFactory.subscribe(DoubleCountCorrelationSum)
