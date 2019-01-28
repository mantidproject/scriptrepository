"""*WIKI* 
Run this script on the Mantid Script window
Analyses a series of images from Trius Monochrome beamline camera
The images and nxs run files are taken while changing a beamline parameter e.g. Vert_Steer_Mag
The camera and SECI work independently.
The image and nxs files are compared using their time stamps.
So make sure your camera is connected to network during the measurement.
Uses classes from FlattenAndFitCamera.py and LoadBeamCameraFile.py
Written by Koji Yokoyama on 28 Jan. 2019
Hint: you will need to set: a scaling factor, initial values for fitting
*WIKI*"""

import time
import datetime
import numpy as np
import math

muonfile = "//hifi/data/cycle_16_3/HIFI%08d.nxs"
camfile = "C:/your_path/IMG%d.FIT"
export_image = "yes" # converts graphs into JPG files 
export_img_file_directory = "C:/your_path/"
piclist = range(2815,2850 +1)
runlist = range(108048,108052 +1) # run numbers from the measurement
beamline_par1 = "Steer_VSM1" # Refer to the sample logs to find the parameter name
resultname = beamline_par1

# ***change these parameters depending on your local time***
format_paramlog = "%Y-%m-%dT%H:%M:%S.000000000" # in normal time "%Y-%m-%dT%H:%M:%S.000000000", in summer time: "%Y-%m-%dT%H:%M:%S.000000000+0100"
adjust = 0 # this is 0 in normal time and 3600 in summer time.

# cam parameters for analysis
image_scale = 198.6 # Pixels per cm.  198.6 for EMU, and 80.0 for HiFi images (1/2 pixels in HiFi)
mask = [0,1391,0,1039] # 120,590,80,440 for HiFi 0,1391,0,1039 for Emu

clear_windows = "" # Close existing windows? "yes" or ""


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

		realXi=np.mod(xvals,10000)
		realY=(xvals-realXi)/10000*scale
		realX=realXi*scale
		realXn=(realX-X0)/XSig
		realYn=(realY-Y0)/YSig
		
		Ampl=Intens/XSig/YSig/math.sqrt(1+Skew**2/4) # peak
		y=BG+Ampl*np.exp(-(realXn**2)-(realYn**2)-realXn*realYn*Skew)
		return y
FunctionFactory.subscribe(TwoDimGaussianPlusBG)

# clears windows
if clear_windows == "yes":
	for w in windows():
	  if w.inherits('MultiLayer'):
	    w.close()
else:
	pass

# reads runs and makes a dictionary of runs and their sample logs
results={}
for run in runlist:
	print "loading ", run
	wstuple=LoadMuonNexus(Filename=muonfile % run)
	try:
		ws = wstuple[0][0] # first period if there are some
		nperiod = 2
	except:
		ws = wstuple[0] # for plain run
		nperiod = 1

	thisres={}
	if(run==runlist[0]):
		# absolute start time in sec when measurement started
		time0 = time.mktime(time.strptime(ws.run().getLogData("run_start").value,"%Y-%m-%dT%H:%M:%S")) 
	
	run_start = time.mktime(time.strptime(ws.run().getLogData("run_start").value,"%Y-%m-%dT%H:%M:%S"))
	run_end = time.mktime(time.strptime(ws.run().getLogData("run_end").value,"%Y-%m-%dT%H:%M:%S"))
	thisres["begin"] = run_start - time0
	thisres["end"] = run_end - time0
	
	events=0
	for i in range(ws.getNumberHistograms()):
		events+=np.sum(ws.dataY(i))
	thisres["rate"]=events/ws.getRun().getProperty("goodfrm").value
	events=0
	for i in range(ws.getNumberHistograms()/2):
		events+=np.sum(ws.dataY(i))
	thisres["rateF"]=events/ws.getRun().getProperty("goodfrm").value
	events=0
	for i in range(ws.getNumberHistograms()/2,ws.getNumberHistograms()):
		events+=np.sum(ws.dataY(i))
	thisres["rateB"]=events/ws.getRun().getProperty("goodfrm").value
	
	# calculates (for weighted average) the beamline parameter that you are changing
	date_time_array = map(str,np.array(ws.run().getLogData(beamline_par1).times))
	relative_times = []
	for date_time in date_time_array:
		relative_times.append(time.mktime(time.strptime(date_time,format_paramlog)) - run_start - adjust)
	value_array = map(float,np.array(ws.run().getLogData(beamline_par1).value))
	relative_times_reduced = []
	value_array_reduced = []
	for idx, ielem in enumerate(relative_times):
		if ielem >= 0:
			relative_times_reduced.append(ielem)
			value_array_reduced.append(value_array[idx])
	w = np.diff(relative_times_reduced, n=1)	
	beamline_par1_ave = np.average(value_array_reduced[:-1], weights = w) # takes a weighted average
	print beamline_par1, " relative_times_reduced =", relative_times_reduced
	print beamline_par1, " value_array_reduced =", value_array_reduced
	
	thisres["Run"] = run
	thisres[beamline_par1] = beamline_par1_ave
	thisres["Slits"] = ws.getRun().getProperty("Slits").value[-1]
	results[run] = thisres
	
print results

# Creates the result table
tt=WorkspaceFactory.createTable()
columns=["Run","Image",beamline_par1,"Slits","rate","X0","eX0","Y0","eY0","XSig","eXSig","YSig","eYSig","Skew","eSkew","Background","eBackground","Intens","eIntens","Major","Minor"]
for co in columns:
	if(co=="Run" or co=="Image"):
		tt.addColumn("int",co,6)
	elif(co[0]=="e"):
		tt.addColumn("double",co,5)
	else:
		tt.addColumn("double",co,2)

# Based on the runs analysed above, now analyses image files
for i in piclist:
	print "processing", (camfile % i)
	img=LoadBeamCameraFile(SourceFile=(camfile % i), FilterNoiseLevel='50', BinSize='1', OutputWorkspace='IMG', Scale=str(image_scale))
	# Find times when the image was taken, relative to time0
	starttime = time.mktime(time.strptime(str(img.getRun().startTime()).strip(),"%Y-%m-%dT%H:%M:%S")) - time0
	endtime = time.mktime(time.strptime(str(img.getRun().endTime()).strip(),"%Y-%m-%dT%H:%M:%S")) - time0
	print "capture duration:", endtime-starttime, "sec \n "
	
	for rn,r in results.items():
		if(starttime>r["begin"] and endtime<r["end"]): # find out which beamline parameter this image corresponding to
			print "image",i,"taken during run",rn
		
			squash=MaskAndFlattenCameraImage(img,*mask)
			
			try:
				fplist=Fit(Function="name=TwoDimGaussianPlusBG",
					InputWorkspace='squash',
					Output='squash',
					OutputCompositeMembers=True,
					StartX=0,
					EndX=100000000,
					OutputNormalisedCovarianceMatrix='squash_NormalisedCovarianceMatrix',
					OutputParameters='squash_Parameters',
					OutputWorkspace='squash_Workspace',
					MaxIterations=20)
				
				#print fplist				
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
				# print "Skew = ", Skew
				# print "Background = ", Background
				# print "Intensity = ", Intens
				# print "Amplitude = ", Amplitude
				# print "area1=",XOverall*YSig, " area2=",Major*Minor

				r["X0"]=X0
				r["Y0"]=Y0
				r["XSig"]=XSig
				r["YSig"]=YSig
				r["Skew"]=Skew
				r["Background"]=Background
				r["Intens"]=Intens
				r["eX0"]=eX0
				r["eY0"]=eY0
				r["eXSig"]=eXSig
				r["eYSig"]=eYSig
				r["eSkew"]=eSkew
				r["eBackground"]=eBackground
				r["eIntens"]=eIntens
				r["Major"]=Major
				r["Minor"]=Minor
				r["Image"]=i
				r["PeakInt"]=PeakInt
				r["XOverall"]=XSig/SFac
				r["YOverall"]=YSig/SFac
				
				if(export_image == "yes"):
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
					lay.setAxisScale(Layer.Bottom, BottomScaleMin, BottomScaleMax)
					lay.setAxisScale(Layer.Right, RightScaleMin, RightScaleMax)
					titl = "Data " + beamline_par1 + " = " + str("{:.2f}".format(r[beamline_par1]))
					lay.setTitle(titl)
					lay.exportImage(export_img_file_directory + "IMG" + str(i) + ".JPG")
					g.close()
					# plot fit image with the Gaussian ellipse curve
					fit_img = UnflattenCameraImage(InputWorkspace=fitWorkspace, OutputWorkspace='FITIMG', Scale=1/image_scale, SpectrumNumber=1)
					g3=plot2D(fit_img)
					g4=plotSpectrum(EllipseCurves,[0,1,2])
					mergePlots(g3,g4)
					lay3=g3.activeLayer()
					lay3.setAxisScale(Layer.Left, LeftScaleMin, LeftScaleMax)
					lay3.setAxisScale(Layer.Bottom, BottomScaleMin, BottomScaleMax)
					lay3.setAxisScale(Layer.Right, RightScaleMin, RightScaleMax)
					titl = "Fit " + beamline_par1 + " = " + str("{:.2f}".format(r[beamline_par1]))
					lay3.setTitle(titl)
					lay3.exportImage(export_img_file_directory + "IMG" + str(i) + "_FIT.JPG")
					g3.close()
	
			except:
				print "fit failed for image :( ",i # don't add to r, another image might be ok	

# put the calculated values in the result table
for r in results.values():
	row=map(lambda x:r[x],columns)
	tt.addRow(row)

AnalysisDataService.addOrReplace(resultname,tt)
