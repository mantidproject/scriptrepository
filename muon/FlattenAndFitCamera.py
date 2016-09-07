from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *
import numpy
import math

class MaskAndFlattenCameraImage(PythonAlgorithm):
	def category(self):
		return "Muon"

	def PyInit(self):
		self.declareProperty(WorkspaceProperty("InputWorkspace","",Direction.Input))
		self.declareProperty(WorkspaceProperty("OutputWorkspace","",Direction.Output))
		self.declareProperty("Xmin",0,Direction.Input)
		self.declareProperty("Xmax",-1,Direction.Input)
		self.declareProperty("Ymin",0,Direction.Input)
		self.declareProperty("Ymax",-1,Direction.Input) # half open interval like Python arrays

	def PyExec(self):
		iws=self.getProperty("InputWorkspace").value
		x1=self.getProperty("Xmin").value
		x2=self.getProperty("Xmax").value
		y1=self.getProperty("Ymin").value
		y2=self.getProperty("Ymax").value
		# defaults
		if(self.getProperty("Xmin").isDefault):
			x1=0
		if(self.getProperty("Xmax").isDefault):
			x2=len(iws.dataY(0))
		if(self.getProperty("Ymin").isDefault):
			y1=0
		if(self.getProperty("Ymax").isDefault):
			y2=iws.getNumberHistograms()

		if(x1<0 or x2<=x1 or x2>len(iws.dataY(0))):
			print "X1=",x1," X2=",x2," Xlimit=",len(iws.dataY(0))
			raise Exception("X range error")
		if(y1<0 or y2<=y1 or y2>iws.getNumberHistograms()):
			print "Y1=",y1," Y2=",y2," Ylimit=",iws.getNumberHistograms()
			raise Exception("Y range error")

		ows=WorkspaceFactory.create("Workspace2D",NVectors=1,XLength=(x2-x1)*(y2-y1),YLength=(x2-x1)*(y2-y1))
		n=0
		for j in range(y1,y2):
			for i in range(x1,x2):
				xx=i+j*10000
				ows.dataX(0)[n]=xx
				ows.dataY(0)[n]=iws.dataY(j)[i]
				ows.dataE(0)[n]=iws.dataE(j)[i]
				n=n+1
		self.setProperty("OutputWorkspace",ows)
		
AlgorithmFactory.subscribe(MaskAndFlattenCameraImage)

class UnflattenCameraImage(PythonAlgorithm):
	def category(self):
		return "Muon"

	def PyInit(self):
		self.declareProperty(WorkspaceProperty("InputWorkspace","",Direction.Input))
		self.declareProperty("SpectrumNumber",0,Direction.Input) # 0 for plain flat file, 1 or 2 for Calc or Diff from fit results
		self.declareProperty(WorkspaceProperty("OutputWorkspace","",Direction.Output))
		self.declareProperty("Scale",1.0/80.0,Direction.Input)
	
	def PyExec(self):
		iws=self.getProperty("InputWorkspace").value
		spec=self.getProperty("SpectrumNumber").value
		xy1=int(iws.dataX(spec)[0])
		x1=xy1 % 10000
		y1=(xy1-x1)/10000
		xy2=int(iws.dataX(spec)[-1])
		x2=xy2 % 10000
		y2=((xy2-x2)/10000 ) + 1
		x2=x2+1
		if((x2-x1)*(y2-y1) != len(iws.dataX(spec))):
			print "guessed x1=",x1," x2=",x2," y1=",y1," y2=",y2," length=", len(iws.dataX(spec))
			raise Exception("Can't unflatten axes")
		scale=self.getProperty("Scale").value

		ows=WorkspaceFactory.create("Workspace2D",NVectors=(y2-y1),XLength=(x2-x1),YLength=(x2-x1))
		n=0
		for j in range(y1,y2):
			for i in range(x1,x2):
				ows.dataX(j-y1)[i-x1]=i*scale
				ows.dataY(j-y1)[i-x1]=iws.dataY(spec)[n]
				ows.dataE(j-y1)[i-x1]=iws.dataE(spec)[n]
				n=n+1
				
		na = NumericAxis.create(y2-y1)
		for i in range(y2-y1):
			na.setValue(i,(i+y1)*scale)
		lbl=na.setUnit("Label")
		lbl.setLabel("height","cm")
		ows.replaceAxis(1,na)

		lbl=ows.getAxis(0).setUnit("Label")
		lbl.setLabel("width","cm")
		ows.setYUnitLabel("brightness")
		
		self.setProperty("OutputWorkspace",ows)

AlgorithmFactory.subscribe(UnflattenCameraImage)

class TwoDimGaussianPlusBG(IFunction1D):
	def init(self):
		self.declareParameter("X0",4.0) # all sizes in cm
		self.declareParameter("Y0",3.0)
		self.declareParameter("XSig",1.0)
		self.declareParameter("YSig",1.0)
		self.declareParameter("Skew",0.1)
		self.declareParameter("Background",1500.0) # per pixel
		self.declareParameter("Intensity",600.0) # integrated

	def function1D(self,xvals):
		X0=self.getParameterValue("X0")
		Y0=self.getParameterValue("Y0")
		XSig=self.getParameterValue("XSig")
		YSig=self.getParameterValue("YSig")
		Skew=self.getParameterValue("Skew")
		BG=self.getParameterValue("Background")
		Intens=self.getParameterValue("Intensity")
		scale=1.0/80.0 # hard coded! Could be Attribute

		realXi=numpy.mod(xvals,10000)
		realY=(xvals-realXi)/10000*scale
		realX=realXi*scale
		realXn=(realX-X0)/XSig
		realYn=(realY-Y0)/YSig
		
		Ampl=Intens/XSig/YSig/math.sqrt(1+Skew**2/4) # peak
		y=BG+Ampl*numpy.exp(-(realXn**2)-(realYn**2)-realXn*realYn*Skew)
		return y

FunctionFactory.subscribe(TwoDimGaussianPlusBG)

class DrawGaussianEllipse(PythonAlgorithm):
	def category(self):
		return "Muon"

	def PyInit(self):
		self.declareProperty(ITableWorkspaceProperty("FitTable","",Direction.Input))
		self.declareProperty(WorkspaceProperty("OutputCurves","",Direction.Output))
		self.declareProperty("MajorAxis",0.0,Direction.Output)
		self.declareProperty("MinorAxis",0.0,Direction.Output)
		self.declareProperty("Rotation",0.0,Direction.Output)
	
	def PyExec(self):
		tab=self.getProperty("FitTable").value
		for i,lab in enumerate(tab.column(0)):
			if(lab=="X0"):
				X0=tab.column(1)[i]
			if(lab=="Y0"):
				Y0=tab.column(1)[i]
			if(lab=="XSig"):
				XSig=tab.column(1)[i]
			if(lab=="YSig"):
				YSig=tab.column(1)[i]
			if(lab=="Skew"):
				Skew=tab.column(1)[i]
			if(lab=="Background"):
				BG=tab.column(1)[i]
			if(lab=="Intensity"):
				Intens=tab.column(1)[i]

		# ellipse at 1 sigma
		# intersects X axis at X-XSig and X+XSig
		# etc
		# or:
		# x=X0+XM*cos(phi)+XL*sin(phi); y=Y0+YM*cos(phi)+YL*sin(phi)
		# for phi=0..2pi, where
		# XM=Major*cos(rot)
		# YM=Major*sin(rot)
		# XL=-Minor*sin(rot)
		# YL=Minor*cos(rot)
		
		# from mathworld...
		a=1/(XSig**2)
		b=Skew/(XSig*YSig)/2.0
		c=1/(YSig**2)
		# d=f=0
		# g=-1
		
		K1=2*(-b**2+a*c)
		K2a=(b**2-a*c)*(math.sqrt((a-c)**2+4*b**2)-(a+c))
		K2b=(b**2-a*c)*(-math.sqrt((a-c)**2+4*b**2)-(a+c))
		
		ax1=math.sqrt(K1/K2a)
		ax2=math.sqrt(K1/K2b)
		if(ax1>ax2):
			major=ax1
			minor=ax2
		else:
			major=ax2
			minor=ax1
		if(b==0 and a<c):
			phi=0.0
		elif(b==0):
			phi=math.pi/2.0
		# elif(a<c):
		#	phi=0.5*math.atan2(2*b,a-c) # math.cot**-1((a-c)/2/b)
		else:
			phi=math.pi/2.0+0.5*math.atan2(2*b,a-c) # math.cot**-1((a-c)/2/b)
		
		self.setProperty("MinorAxis",minor)
		print "major axis length: ",major
		self.setProperty("MajorAxis",major)
		print "minor axis length: ",minor
		self.setProperty("Rotation",phi)
		print "rotation (degrees anticlockwise from x): ",phi*180.0/math.pi
		XM=major*math.cos(phi)
		XL=-minor*math.sin(phi)
		YM=major*math.sin(phi)
		YL=minor*math.cos(phi)
	
		phidat=numpy.linspace(0.0,2*math.pi,37)
		cosdat=numpy.cos(phidat)
		sindat=numpy.sin(phidat)
		xelip=X0+XM*cosdat+XL*sindat
		yelip=Y0+YM*cosdat+YL*sindat
		
		# axes
		lindat=numpy.linspace(-1.0,1.0,37)
		minorX=X0+XL*lindat
		minorY=Y0+YL*lindat
		majorX=X0+XM*lindat
		majorY=Y0+YM*lindat

		ows=WorkspaceFactory.create("Workspace2D",NVectors=3,XLength=37,YLength=37)
		ows.dataX(0)[:]=xelip
		ows.dataY(0)[:]=yelip
		ows.dataX(1)[:]=majorX
		ows.dataY(1)[:]=majorY
		ows.dataX(2)[:]=minorX
		ows.dataY(2)[:]=minorY
		na = TextAxis.create(3)
		na.setLabel(0,"Ellipse")
		na.setLabel(1,"Major axis")
		na.setLabel(2,"Minor axis")
		ows.replaceAxis(1,na)

		self.setProperty("OutputCurves",ows)
		
AlgorithmFactory.subscribe(DrawGaussianEllipse)
