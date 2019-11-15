# load FITS file from HIFI Starlight camera into Workspace
# set start and end times from header
# optional filtering out of noise points (over specified delta from all 4 neighbours)
# optional binning
# set scale of image
from __future__ import print_function
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *
import numpy

class LoadBeamCameraFile(PythonAlgorithm):
	def category(self):
		return 'Muon'

	def PyInit(self):
		self.declareProperty(FileProperty(name="SourceFile",defaultValue="",action=FileAction.Load,extensions=["FIT"]))
		self.declareProperty("FilterNoiseLevel",0,doc="Threshold to remove noisy pixels, try 50")
		self.declareProperty("BinSize",1,doc="Rebunch n*n on both axes")
		self.declareProperty("Scale",80.0,doc="Pixels per cm")
		self.declareProperty(WorkspaceProperty('OutputWorkspace','',Direction.Output))
		
	def PyExec(self):
		fnam=self.getProperty("SourceFile").value
		fh=open(fnam,'rb')
		# read FITS header looking for lines
		width=0
		height=0
		nbuf=0
		datestr=""
		timestr=""
		exposurelen=0
		for nlr in range(299):
			if(nbuf<=0):
				lines = fh.read(80*36)
				nbuf=36
			line=lines[80*(36-nbuf):80*(36-nbuf)+80]
			nbuf=nbuf-1
			if(line[0:3]=="END"):
				print("found the end of the header")
				break
			if(line[0:6]=="NAXIS1"):
				print("found axis1 line: ",line)
				width=int(line[9:30])
			if(line[0:6]=="NAXIS2"):
				print("found axis2 line: ",line)
				height=int(line[9:30])
			if("DATE-OBS" in line):
				datestr=line.split("'")[1]
			if("TIME-OBS" in line):
				timestr=line.split("'")[1]
			if("EXPTIME" in line):
				exposurelen=float(line[9:])
			if('END' in line):
				print("nearly got the end in [",line,"]")

		print("read ",nlr," potential header lines")
		fh.seek(0,1) # back "here", flush read ahead buffers
		if(height==0 or width==0):
			raise Exception("Failed to identify image size")

		dt=numpy.dtype('>i2')
		flatimage=numpy.fromfile(fh,dtype=dt,count=width*height)
		image=numpy.reshape(flatimage,[width,height],order="F")
		image=image+32768

		scale=1.0/self.getProperty("Scale").value

		# filter?
		filt=self.getProperty("FilterNoiseLevel").value
		if(filt>0):
			image1=numpy.zeros([width+2,height+2])
			image1[0,:]=-65536 # bad pixels at edge will differ from this
			image1[:,0]=-65536
			image1[-1,:]=-65536
			image1[:,-1]=-65536
			image1[1:-1,1:-1]=image
			image2=numpy.array(image1)
			good=numpy.ones([width+2,height+2])
			good[0,:]=0 # don't copy from border
			good[:,0]=0
			good[-1,:]=0
			good[:,-1]=0
			image2[0,:]=0 # no stray weight from border either
			image2[:,0]=0
			image2[-1,:]=0
			image2[:,-1]=0
			for i in range(1,width+1):
				for j in range(1,height+1):
					if(abs(image1[i,j]-image1[i-1,j])>filt and abs(image1[i,j]-image1[i+1,j])>filt and abs(image1[i,j]-image1[i,j-1])>filt and abs(image1[i,j]-image1[i,j+1])>filt):
						image2[i,j]=0
						good[i,j]=0
			image3=numpy.array(image1)
			for i in range(1,width+1):
				for j in range(1,height+1):
					if(good[i,j]==0):
						ngn=good[i-1,j]+good[i+1,j]+good[i,j-1]+good[i,j+1]
						if(ngn>0): # at least one good pixel to copy from, otherwise leave alone
							image3[i,j]=(image2[i-1,j]+image2[i+1,j]+image2[i,j-1]+image2[i,j+1])/ngn
			image=image3[1:-1,1:-1]

		# resize?
		bin=self.getProperty("BinSize").value
		ys=1.0
		if(bin>1):
			newWid=int(width/bin)
			newHei=int(height/bin)
			# pad is not (yet) in the old numpy bundled with Mantid
			#imageP=numpy.pad(image,((1,0),(1,0)),mode='constant')
			imageP=numpy.insert(image,0,0,0)
			imCS=numpy.cumsum(imageP,axis=0)
			imCSB=imCS[0:(newWid+1)*bin:bin,:] # 0:(newHei+1)*bin:bin
			imCSUX=imCSB[1:,:]-imCSB[:-1,:]
			imageP=numpy.insert(imCSUX,0,0,1)
			imCS=numpy.cumsum(imageP,axis=1)
			imCSB=imCS[:,0:(newHei+1)*bin:bin] # 
			image=imCSB[:,1:]-imCSB[:,:-1]
			width=newWid
			height=newHei
			scale=scale*bin
			ys=1.0/(bin*bin)
				
		ws=WorkspaceFactory.create("Workspace2D",NVectors=height,XLength=width+1,YLength=width)
		ws.setDistribution(True) # always polarisation or similar normalised value, not raw or simulated counts
		xaxdat=numpy.linspace(0.0,width*scale,width+1)
		yaxdat=numpy.linspace(0.0,(height-1)*scale,height)
		
		for i in range(height):
			ws.dataX(i)[:]=xaxdat
			ws.dataY(i)[:]=image[:,height-1-i]*ys
		na = NumericAxis.create(height)
		for i in range(height):
			na.setValue(i,yaxdat[i])
		lbl=na.setUnit("Label")
		lbl.setLabel("height","cm")
		ws.replaceAxis(1,na)

		lbl=ws.getAxis(0).setUnit("Label")
		lbl.setLabel("width","cm")
		ws.setYUnitLabel("brightness")

		if(timestr != "" and datestr != ""):
			starttime=DateAndTime(datestr+"T"+timestr)
			endtime=DateAndTime(int(starttime.total_nanoseconds()+exposurelen*1000000000))
			ws.getRun().setStartAndEndTime(starttime,endtime)

		self.setProperty('OutputWorkspace',ws)

AlgorithmFactory.subscribe(LoadBeamCameraFile)
