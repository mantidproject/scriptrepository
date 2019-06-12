"""*WIKI* 
Run this script on the Mantid Script window
Analyses a single image from Trius Monochrome beamline camera
Uses classes from FlattenAndFitCamera.py and LoadBeamCameraFile.py
Written by Koji Yokoyama on 23 Jan. 2019
Hint: you will need to set: a scaling factor, initial values for fitting
*WIKI*"""

import numpy
import math
import sys

image_filename = "C:/your_path/vsm-1.5.fit"
par = -1.5
par_label = "vsm" # "par" is only for the result table
result_table = "result_table"
image_scale = 198.6 # Pixels per cm.  198.6 for EMU, and 80.0 for HiFi images (1/2 pixels in HiFi)
mask = [0,1391,0,1039] # 120,590,80,440 for HiFi 0,1391,0,1039 for Emu
clear_windows = "yes" # Close existing windows? "yes" or ""

# clears windows
if clear_windows == "yes":
	for w in windows():
	  if w.inherits('MultiLayer'):
	    w.close()
else:
	pass

# defines fit function and initial parameters for fitting
class TwoDimGaussianPlusBG(IFunction1D):
	def init(self):
		self.declareParameter("X0",3.6) # all sizes in cm
		self.declareParameter("Y0",2.6)
		self.declareParameter("XSig",1.0)
		self.declareParameter("YSig",1.0)
		self.declareParameter("Skew",0.1)
		self.declareParameter("Background",1400.0) # counts per pixel
		self.declareParameter("Intensity",200.0) # integrated

	def function1D(self,xvals):
		X0=self.getParameterValue("X0")
		Y0=self.getParameterValue("Y0")
		XSig=self.getParameterValue("XSig")
		YSig=self.getParameterValue("YSig")
		Skew=self.getParameterValue("Skew")
		BG=self.getParameterValue("Background")
		Intens=self.getParameterValue("Intensity")
		scale=1.0/image_scale # hard coded! Could be Attribute

		realXi=numpy.mod(xvals,10000)
		realY=(xvals-realXi)/10000*scale
		realX=realXi*scale
		realXn=(realX-X0)/XSig
		realYn=(realY-Y0)/YSig
		
		Ampl=Intens/XSig/YSig/math.sqrt(1+Skew**2/4) # peak
		y=BG+Ampl*numpy.exp(-(realXn**2)-(realYn**2)-realXn*realYn*Skew)
		return y

FunctionFactory.subscribe(TwoDimGaussianPlusBG)

img=LoadBeamCameraFile(SourceFile=image_filename, FilterNoiseLevel='50', BinSize='1', OutputWorkspace='IMG', Scale=str(image_scale))

starttime=time.mktime(time.strptime(str(img.getRun().startTime()).strip(),"%Y-%m-%dT%H:%M:%S"))
endtime=time.mktime(time.strptime(str(img.getRun().endTime()).strip(),"%Y-%m-%dT%H:%M:%S"))
print "capture duration:", endtime-starttime, "sec \n Fit parameters:"

squash=MaskAndFlattenCameraImage(img,*mask)
# Do fitting
fplist=Fit(Function="name=TwoDimGaussianPlusBG",
	InputWorkspace='squash',
	Output='squash',
	OutputCompositeMembers=True,
	StartX=0,
	EndX=100000000,
	OutputNormalisedCovarianceMatrix='squash_NormalisedCovarianceMatrix',
	OutputParameters='squash_Parameters',
	OutputWorkspace='squash_Workspace',
	MaxIterations=50)

# print fplist
Params = fplist.OutputParameters
CostFn = fplist.OutputChi2overDoF
fitWorkspace = fplist.OutputWorkspace
(EllipseCurves,Major,Minor,Phi)=DrawGaussianEllipse(Params)

(X0,Y0,XSig,YSig,Skew,Background,Intens,CostF2)=Params.column(1)
(eX0,eY0,eXSig,eYSig,eSkew,eBackground,eIntens,eCostF2)=Params.column(2)
SFac=math.sqrt(1-Skew**2/4)
PeakInt=Intens/XSig/YSig*SFac
XOverall=XSig/SFac
YOverall=YSig/SFac
Amplitude = Intens/XSig/YSig/math.sqrt(1+Skew**2/4) # peak

print "X centre = ", X0
print "Y centre = ", Y0
print "Skew = ", Skew
print "Background = ", Background
print "Intensity = ", Intens
print "Amplitude = ", Amplitude
print "area1=",XOverall*YSig, " area2=",Major*Minor

# Generates the result table with the parameter "par" and calculated fit params
try:
	tt = mtd["result_table"]
except:
	print "Workspace result_table doesn't exit. Creating one ..."
	tt = WorkspaceFactory.createTable()
	columns = [par_label,"X0","eX0","Y0","eY0","XSig","eXSig","YSig","eYSig","Skew","eSkew","Background","eBackground","Intens","eIntens","Major","Minor","Phi"]
	for co in columns:
		if(co==par_label):
			tt.addColumn("double",co,1)
		elif(co[0]=="e"):
			tt.addColumn("double",co,5)
		else:
			tt.addColumn("double",co,2)
row = [par,X0,eX0,Y0,eY0,XSig,eXSig,YSig,eYSig,Skew,eSkew,Background,eBackground,Intens,eIntens,Major,Minor,Phi*180/math.pi]
tt.addRow(row)
AnalysisDataService.addOrReplace(result_table,tt)

# plot data with the Gaussian ellipse curve
RightScaleMin = Background
RightScaleMax = Background + Amplitude
LeftScaleMin = -0.2
LeftScaleMax = 6.8
BottomScaleMin = -0.2
BottomScaleMax = 8.8

g=plot2D(img)
g2=plotSpectrum(EllipseCurves,[0,1,2])
mergePlots(g,g2)
lay=g.activeLayer()
lay.setAxisScale(Layer.Left, LeftScaleMin, LeftScaleMax)
lay.setAxisScale(Layer.Bottom, BottomScaleMin,BottomScaleMax)
lay.setAxisScale(Layer.Right, RightScaleMin, RightScaleMax)
# plot 2D fit function with the Gaussian ellipse curve
fit_img = UnflattenCameraImage(InputWorkspace=fitWorkspace, OutputWorkspace='FITIMG', Scale=1/image_scale, SpectrumNumber=1)
g3=plot2D(fit_img)
g4=plotSpectrum(EllipseCurves,[0,1,2])
mergePlots(g3,g4)
lay3=g3.activeLayer()
lay3.setAxisScale(Layer.Left, LeftScaleMin, LeftScaleMax)
lay3.setAxisScale(Layer.Bottom, BottomScaleMin,BottomScaleMax)
lay3.setAxisScale(Layer.Right, RightScaleMin, RightScaleMax)
