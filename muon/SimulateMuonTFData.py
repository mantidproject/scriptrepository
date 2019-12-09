"""*WIKI* 

Algorithm SetupSimInstrument defines which instrument is to be modelled, and chooses sample offsets to define alpha, etc.
It outputs a setup table for use by SimulateMuonTFData.

Algorithm SimulateMuonTFData generates a simulated muon run workspace, based on a model or a user specified function, field value and counting time (statistics).
It includes the frequency response, dead time effects, and random variations of detector position, alpha, etc as defined in the setup table.
Use the same input table each time to simulate successive runs in the same experiment.

*WIKI*"""

import time
import datetime
import math
import numpy
import scipy
from mantid.api import * # PythonAlgorithm, registerAlgorithm, WorkspaceProperty
from mantid.kernel import *
from mantid.simpleapi import *

class SetupSimInstrument(PythonAlgorithm):
	def PyInit(self):
		# instrument name (for detector positions and num of histograms) - pass to CreateSimulationWorkspace()
		self.declareProperty("Instrument","",validator=StringListValidator(["","MUSR","EMU","HIFI","DEVA","ARGUS","CHRONUS"]),doc="to find IDF")
		self.declareProperty("Hint","",StringListValidator(["","Standard LF","Rotated for TF","T20 Coils"]),doc="Suggested field and instrument setup")
		self.declareProperty("FieldAxis","",doc="Override guess of magnetic field direction with +-xyz or 3-vector") # +-x,y,z or 3-vector
		self.declareProperty("MuonSpin","",doc="Override guess of initial muon spin with +-xyz or 3-vector")
		# num of bins. Assume start at t=0
		self.declareProperty("NBins",2048,"Number of data bins in a histogram")
		# bin width (ns)
		self.declareProperty("BinWidth",0.016,"Bin width in microseconds")
		# time zero
		self.declareProperty("T0",0.3654,"Muon arrival time relative to histogram start")
		self.declareProperty("T0Scatter",0.005,"T0 RMS variation from detector to detector")
		self.declareProperty("DetPosScatter",2.0,"p-p position variation (mm) in 3D from detector to detector, affects phases")
		# pulse parameters
		self.declareProperty("ProtonPulseWidth",0.1,"FWHM of muon pulse, for count build up around t0 and frequency response")
		self.declareProperty("PionLifetime",0.026,"Pion lifetime, for frequency response")
		# full asymmetry (z)
		self.declareProperty("FullAsym",0.30,"Full asymmetry with polarisation pointing at a detector")
		self.declareProperty("FullAsymScatter",0.02,"Random p-p variation detector to detector")
		# sample offset, defines alpha
		self.declareProperty("SampleX",1.0,"Sample/beam x position in mm, adjusts alpha (x)")
		self.declareProperty("SampleY",1.0,"Sample/beam y position in mm, adjusts alpha (y)")
		self.declareProperty("SampleZ",1.0,"Sample z position in m, adjusts alpha")
		self.declareProperty("CountRate",50.0,"Total in Mevents/hour, for good frame and dead time estimation, 40Hz beam assumed")
		self.declareProperty("CountScatter",1.3,"count rate variation (ratio P-P) detector to detector")
		self.declareProperty("DeadTime",0.01,"Mean dead time (us)")
		self.declareProperty("DeadTimeScatter",0.005,"Dead time variation p-p (us)")
		self.declareProperty("DeadDets",0,"How many to mark as completely dead")
		# output a definition table
		self.declareProperty(ITableWorkspaceProperty('InstrumentSetupTable','',Direction.Output))

	def category(self):
		return "Muon"
		
	def PyExec(self):
		# table has columns detnum ax ay az t0 countrate
		# title of column detnum is instrument name!
		# extra rows column 0 <= 0:
		# -1: field axis (ax,ay,az)
		# -2: muon spin (ax,ay,az)
		# -3: time bins (-t0, dt, nbins,pulsewidth)
		# -4: sample location (ax,ay,az)? Not needed
		inst=self.getProperty("Instrument").value
		cdws=self.createChildAlgorithm("CreateSimulationWorkspace")
		cdws.setProperty("Instrument",inst)
		cdws.setProperty("BinParams","0,1,32")
		cdws.execute()
		dws=cdws.getProperty("OutputWorkspace").value
		#tidy_dws=True
		self.log().notice("Looking up default IDF for "+inst)
		
		# axis-guessing code from CreatePhaseQuadTable2
		axvec=None
		idvec=None
		main_field_direction=""
		# standard coil directions/polarities
		# coords: x to left as seen by beam, i.e. towards EPB for south side instruments, y upwards, z downstream (except MUSR in TF, when z towards HIFI, x downstream, y still up)
		MuSRmain=V3D(0,0,-1) # ok
		MuSRbeamLF=V3D(0,0.1,-1) # ok
		MuSRbeamTF=V3D(-1,0.1,0) # ok
		MuSRT20coil=V3D(1,0,0) # ok
		EMULF=V3D(0,0,-1) # ok
		EMUTF=V3D(0,-1,0) # ok
		EMUbeam=V3D(0.05,0.1,-1) # ok
		HIFImain=V3D(0,0,1) # ok
		HIFITFx=V3D(1,0,0) # ok
		HIFITFy=V3D(0,1,0) # ok
		HIFILFz=V3D(0,0,-1) # ok
		HIFIbeam=V3D(-0.05,0.1,-1) # fairly close
		
		hint=self.getProperty("Hint").value
		#self.log().notice("The field hint is "+hint)
		if(inst=="MUSR"):
			if(hint=="Standard LF"):
				idvec=MuSRbeamLF
				axvec=MuSRmain
				self.log().notice("Setting up MUSR in LF, main coils")
				main_field_direction="_L"
			elif(hint=="T20 Coils"):
				idvec=MuSRbeamLF
				axvec=MuSRT20coil
				self.log().notice("Setting up MUSR in LF, T20 coils")
				main_field_direction="_L"
			elif(hint=="Rotated for TF"):
				idvec=MuSRbeamTF
				axvec=MuSRmain
				self.log().notice("Setting up MUSR in TF, main coils")
				main_field_direction="_T"
			else: # guess based on magnets, will fail and not set anything if no sample run available
				try:
					if(dws.getRun().getProperty("a_selected_magnet").value[-1] == "T20 Coils"):
						idvec=MuSRbeamLF
						axvec=MuSRT20coil
						self.log().notice("Setting up MUSR in LF, T20 coils")
						main_field_direction="_L"
					elif(dws.getRun().getProperty("a_selected_magnet").value[-1] == "Danfysik"): # assume TF orientation
						idvec=MuSRbeamTF
						axvec=MuSRmain
						self.log().notice("Setting up MUSR, assume TF orientation, main coils")
						main_field_direction="_T"
					else:
						self.log().notice("Don't know about MUSR magnet "+dws.getRun().getProperty("a_selected_magnet").value[-1])
						main_field_direction="_T" # guess
				except:
					self.log().warning("Can't get information about field direction, try Hint or axes")
					main_field_direction="_T" # guess
		elif(inst=="EMU"):
			if(hint=="Standard LF"):
				idvec=EMUbeam
				axvec=EMULF
				self.log().notice("Setting up EMU, main coils")
			elif(hint=="T20 Coils"):
				idvec=EMUbeam
				axvec=EMUTF
				self.log().notice("Setting up EMU, TF coils")
			elif(hint=="Rotated for TF"):
				raise Exception("EMU can't be rotated")
			else: # guess based on magnets, will fail and not set anything if no sample run available
				try:
					if(dws.getRun().getProperty("a_selected_magnet").value[-1] == "T20 Coils"):
						idvec=EMUbeam
						axvec=EMUTF
						self.log().notice("Setting up EMU, TF coils")
					elif(dws.getRun().getProperty("a_selected_magnet").value[-1] == "Danfysik"):
						idvec=EMUbeam
						axvec=EMULF
						self.log().notice("Setting up EMU, LF coils")
					else:
						self.log().notice("Don't know about EMU magnet "+dws.getRun().getProperty("a_selected_magnet").value[-1])
				except:
					self.log().warning("Can't get information about field direction, try Hint or axes")
		elif(inst=="HIFI"):
			self.log().notice("It's HIFI...")
			if(hint=="Standard LF"):
				idvec=HIFIbeam
				axvec=HIFImain
				self.log().notice("Setting up HIFI, main LF coils")
			elif(hint=="T20 Coils"):
				idvec=HIFIbeam
				axvec=HIFITFy
				self.log().notice("Setting up HIFI, TF y coils")
			elif(hint=="Rotated for TF"):
				raise Exception("HIFI can't be rotated")
			else: # guess based on magnets, will fail and not set anything if no sample run available
				try:
					if(dws.getRun().getProperty("a_selected_magnet").value[-1] == "X Transverse"):
						idvec=HIFIbeam
						axvec=HIFITFx
						self.log().notice("Setting up HIFI, TF x coils")
					elif(dws.getRun().getProperty("a_selected_magnet").value[-1] == "Y Transverse"):
						idvec=HIFIbeam
						axvec=HIFITFy
						self.log().notice("Setting up HIFI, TF y coils")
					elif(dws.getRun().getProperty("a_selected_magnet").value[-1] == "Z Sweep"):
						idvec=HIFIbeam
						axvec=HIFILFz
						self.log().notice("Setting up HIFI, low LF z coils")
					elif(dws.getRun().getProperty("a_selected_magnet").value[-1] == "Main Longitudinal"):
						idvec=HIFIbeam
						axvec=HIFImain
						self.log().notice("Setting up HIFI, main LF coils")
					else:
						self.log().notice("Don't know about HIFI magnet "+dws.getRun().getProperty("a_selected_magnet").value[-1])
				except:
					self.log().warning("Can't get information about field direction, try Hint or axes")
		elif(inst=="ARGUS"):
			self.log().warning("Guess not set up for ARGUS, specify axes")
		elif(inst=="CHRONUS"):
			self.log().warning("Guess not set up for CHRONUS, specify axes")
		else:
			self.log().warning("No idea about instrument "+inst)

		if not self.getProperty("FieldAxis").isDefault:
			ax=self.getProperty("FieldAxis").value
			axnames={"+x":"1,0,0","+y":"0,1,0","+z":"0,0,1","-x":"-1,0,0","-y":"0,-1,0","-z":"0,0,-1"}
			if(ax in axnames):
				ax=axnames[ax]
			axvec=V3D(*list(map(float,ax.split(","))))
		axvec=axvec*(1.0/axvec.norm())

		if not self.getProperty("MuonSpin").isDefault:
			id=self.getProperty("MuonSpin").value
			if(id in axnames):
				id=axnames[id]
			idvec=V3D(*list(map(float,id.split(","))))
		idvec=idvec*(1.0/idvec.norm())
		
		nr=dws.getNumberHistograms()

		tab=WorkspaceFactory.createTable()
		tab.addColumn("int",inst+main_field_direction)
		tab.addColumn("double","ax")
		tab.addColumn("double","ay")
		tab.addColumn("double","az")
		tab.addColumn("double","t0")
		tab.addColumn("double","countrate")
		tab.addColumn("double","deadtime")
		tab.addRow([-3,self.getProperty("T0").value,self.getProperty("BinWidth").value,self.getProperty("NBins").value,self.getProperty("ProtonPulseWidth").value,self.getProperty("PionLifetime").value,0.0])
		tab.addRow([-2,idvec.X(),idvec.Y(),idvec.Z(),0.0,0.0,0.0])
		tab.addRow([-1,axvec.X(),axvec.Y(),axvec.Z(),0.0,0.0,0.0])

		deads=numpy.random.choice(nr, self.getProperty("DeadDets").value,replace=False)

		DetPosScatter=self.getProperty("DetPosScatter").value * 0.001
		SamplePosOffset=V3D(
			self.getProperty("SampleX").value * 0.001,
			self.getProperty("SampleY").value * 0.001,
			self.getProperty("SampleZ").value * 0.001 )

		amplitudes=self.getProperty("FullAsym").value + numpy.random.uniform(-1.0,1.0,size=nr)*self.getProperty("FullAsymScatter").value
		timezeros=self.getProperty("T0Scatter").value * numpy.random.uniform(-1.0,1.0,size=nr)
		countrates= numpy.random.uniform(1.0,self.getProperty("CountScatter").value,size=nr)
		for d in deads:
			countrates[d]=0.0
		totalcount=self.getProperty("CountRate").value * 1.E6/3600.0/40.0 # counts/frame
		countrates=countrates*totalcount/numpy.sum(countrates)

		dtmean=self.getProperty("DeadTime").value
		dtvar=self.getProperty("DeadTimeScatter").value
		deadtime=numpy.random.uniform(dtmean-dtvar/2,dtmean+dtvar/2,size=nr)

		for d in range(nr):
			if(d in deads):
				tab.addRow([d,0.0,0.0,0.0,0.0,0.0])
			else:
				det=V3D(*numpy.random.uniform(low=-DetPosScatter,high=DetPosScatter,size=3))+dws.getDetector(d).getPos() -dws.getInstrument().getSample().getPos() - SamplePosOffset
				det=det*(amplitudes[d]/det.norm())
				tab.addRow([d,det.X(),det.Y(),det.Z(),timezeros[d],countrates[d],deadtime[d]])

		del dws
		self.setProperty('InstrumentSetupTable',tab)
				
AlgorithmFactory.subscribe(SetupSimInstrument)

class SimulateMuonTFData(PythonAlgorithm):
	def PyInit(self):
		# instrument definition table
		self.declareProperty(ITableWorkspaceProperty('InstrumentSetupTable','',Direction.Input))
		self.declareProperty("SampleModel","Calibration",StringListValidator(["Calibration","ExpDecay","GaussDecay","Shallow Donor","Flux Line lattice","Antiferromagnet","UserDefined"]))
		self.declareProperty("Field",20.0,"Magnetic field in Gauss")
		self.declareProperty("Width",10.0,"Lineshape parameter 1: sigma/lambda (us-1), A (MHz), lambda(nm), B_int(G)")
		self.declareProperty("Width2",1.0,"Lineshape parameter 2: D (MHz), correl length (nm), modulation of B_int(G)")
		self.declareProperty("Shift",1.0,"Knight shift parameter in Gauss, for predefined models")
		self.declareProperty("LambdaSig",0.05,"Relaxation rate applied to signal, for predefined models")
		self.declareProperty("Background",0.2,"Fraction of diamagnetic background, for predefined models")
		self.declareProperty("LambdaBG",0.01,"Relaxation rate applied to background, for predefined models")
		# user relaxation (polarisation) function (Z)
		self.declareProperty("Gz","math.exp(-0.1*t)","Longitudinal polarisation function, t=time in microseconds")
		# optional user X and Y polarisation for TF (needs detector positions)
		self.declareProperty("Gx","math.cos(5.0*2*math.pi*t)","Transverse (X) polarisation function, t=time in microseconds")
		self.declareProperty("Gy","math.sin(5.0*2*math.pi*t)","Transverse (Y) polarisation function, t=time in microseconds")
		# total Mevents
		self.declareProperty("Mevents",10.0,"Total Mevents in run (nominal)")
		self.declareProperty("CountRate",0.0,"Rate in Mevents/hour, overriding SetupSimInstrument value")
		# output workspace
		self.declareProperty(WorkspaceProperty("OutputWorkspace","",Direction.Output),"The Simulated Data File")

	def category(self):
		return "Muon"

	def PyExec(self):
		# get stuff
		tab=self.getProperty('InstrumentSetupTable').value
		inst=list(tab.keys())[0]
		orient=None
		if "_" in inst:
			(inst1,orient)=inst.split("_")
		else:
			inst1=inst
		counters=0.0 # total counting capacity
		for r in tab:
			if(r[inst]==-3):
				T0=r["ax"]
				dT=r["ay"]
				NBins=r["az"]
				pwidth=r["t0"]/2.0
				pilife=r["countrate"]
			elif(r[inst]==-2):
				idvec=V3D(r["ax"],r["ay"],r["az"])
			elif(r[inst]==-1):
				axvec=V3D(r["ax"],r["ay"],r["az"])
			else:
				counters += r["countrate"]
		frames=self.getProperty("Mevents").value * 1.E6/counters
		if(self.getProperty("CountRate").value>0):
			frames=frames*(counters*40.0*3600/1.E6)/self.getProperty("CountRate").value
		frames=max(int(frames),1) # non-zero integer!
		# precalculate pulse shape function (kernel)
		subsample=int(dT/min(pilife,pwidth)*20)+1 # at least 20 points per pion life
		self.log().notice("Sub sampling "+str(subsample))
		dT2=dT/subsample
		pslen=int((4*pilife+1.1*pwidth)/dT2+1) # half length of kernel, "sufficient"
		ktime=numpy.linspace(start=-dT2*pslen,stop=dT2*pslen,num=2*pslen+1)
		proton=numpy.maximum(1.0-(ktime/pwidth)**2,0.0)
		kedges=numpy.linspace(start=-dT2*(pslen+0.5),stop=dT2*(pslen+0.5),num=2*pslen+2)
		pionInteg=numpy.minimum(numpy.exp(-kedges/pilife-1.0),1.0)
		pion=pionInteg[1:]-pionInteg[:-1]
		pulseshape=scipy.convolve(proton,pion,mode="same")
		pulseshape=pulseshape/numpy.sum(pulseshape) # normalised!
		
		cdws=self.createChildAlgorithm("CreateSimulationWorkspace")
		cdws.setProperty("Instrument",inst1)
		cdws.setProperty("BinParams",[-T0,dT,NBins*dT-T0])
		cdws.setProperty("UnitX","TOF")
		cdws.execute()
		ws=cdws.getProperty("OutputWorkspace").value
		lbl=ws.getAxis(0).setUnit("Label")
		lbl.setLabel("Time","microsecond")
		#binCentres=(ws.dataX(0)[1:]+ws.dataX(0)[:-1])/2.0
		binCentres=numpy.linspace(-T0+dT/2,-T0+(NBins*subsample-1)*dT2+dT/2,NBins*subsample) # with sub sampling, 1st bin is the centre of the 1st full size one
		N0=self.getProperty("Mevents").value * 1.E6 /counters * (dT/2.19703) # counts per bin (per unit counts/frame factor) close to T0
		
		model=self.getProperty("SampleModel").value
		modelset={"Calibration":None,"ExpDecay":self.sampleExpOsc,"GaussDecay":self.sampleGaussOsc,"Shallow Donor":self.sampleShallowDonor,"Flux Line lattice":self.sampleFluxLineLattice,"Antiferromagnet":self.sampleAntiferromagnet,"UserDefined":None}
		modfunc=modelset[model]
		
		field=self.getProperty("Field").value
		ws.getRun().addProperty("sample_magn_field",field,replace=True)
		ws.getRun().addProperty("goodfrm",frames,replace=True)
		if(orient=="L"):
			ws.getRun().addProperty("main_field_direction","Longitudinal",replace=True)
		if(orient=="T"):
			ws.getRun().addProperty("main_field_direction","Transverse",replace=True)
		width=self.getProperty("Width").value
		width2=self.getProperty("Width2").value
		shift=self.getProperty("Shift").value
		LambdaSig=self.getProperty("LambdaSig").value
		LambdaBG=self.getProperty("LambdaBG").value
		bgfrac=self.getProperty("Background").value
		Gz=self.getProperty("Gz").value
		Gx=self.getProperty("Gx").value
		Gy=self.getProperty("Gy").value

		idperp=(idvec-axvec*idvec.scalar_prod(axvec)*(1.0/axvec.norm()**2)) # component of muon spin perp to field
		idparal=idvec-idperp # component of muon spin along field
		ipd=idperp.cross_prod(axvec) # component of muon spin perp to field, rotated round field by 90 deg
		#self.log().notice("idparal is "+str(idparal))
		#self.log().notice("idperp is "+str(idperp))
		#self.log().notice("idrot is "+str(ipd))

		# precalculate the more difficult distributions if required, for speed
		if (modfunc):
			precalc=modfunc(None,field+shift,width,width2,None)

		for r in tab:
			i=r[inst]
			if(i>=0):
				thisX=binCentres+r["t0"]
				thisAmpl=V3D(r["ax"],r["ay"],r["az"])
				amplPar=thisAmpl.scalar_prod(idparal)
				amplPerp=thisAmpl.scalar_prod(idperp)
				amplRot=thisAmpl.scalar_prod(ipd)
				#self.log().notice("detector {0}, amplitude factors z={1},x={2},y={3}".format(i,amplPar,amplPerp,amplRot) )
				gx=numpy.cos(thisX*field*0.01355*2*math.pi)*numpy.exp(-thisX*LambdaBG)
				gy=numpy.sin(thisX*field*0.01355*2*math.pi)*numpy.exp(-thisX*LambdaBG)
				gz=numpy.exp(-thisX*LambdaBG)
				if(model!="Calibration"):
					if(model=="UserDefined"):
						gx1=numpy.zeros([len(thisX)])
						gy1=numpy.zeros([len(thisX)])
						gz1=numpy.zeros([len(thisX)])
						for ix in range(NBins):
							t=thisX[ix]
							if(Gx):
								gx1[i]=eval(Gx)
							else:
								gx1[i]=1.0
							if(Gy):
								gy1[i]=eval(Gy)
							else:
								gy1[i]=0.0
							if(Gz):
								gz1[i]=eval(Gz)
							else:
								gz1[i]=1.0
					else:
						gx1,gy1,gz1=modfunc(thisX,field+shift,width,width2,precalc)
					gx=gx1*(1-bgfrac)*numpy.exp(-thisX*LambdaSig)+gx*bgfrac
					gy=gy1*(1-bgfrac)*numpy.exp(-thisX*LambdaSig)+gy*bgfrac
					gz=gz1*(1-bgfrac)*numpy.exp(-thisX*LambdaSig)+gz*bgfrac
				# now gz is the component to preserve along field axis
				# gx is the component of initial spin perp to field to preserve along that axis
				# gy is the component that appears along the third axis
				#self.log().notice("gx span {0} to {1}".format(numpy.amin(gx),numpy.amax(gx)))
				#self.log().notice("gy span {0} to {1}".format(numpy.amin(gy),numpy.amax(gy)))
				#self.log().notice("gz span {0} to {1}".format(numpy.amin(gz),numpy.amax(gz)))
				ampl=1.0+amplPar*gz+amplPerp*gx+amplRot*gy
				muflux=N0*numpy.exp(-thisX/2.19703)*r["countrate"] # total counts per bin for this detector, without asymmetry modulation
				muflux=numpy.where(thisX>0,muflux,0)
				#for ix in range(len(muflux)):
				#	if(thisX[ix]<-pwidth):
				#		muflux[ix]=1.E-10
				#	elif(thisX[ix]<pwidth):
				#		muflux[ix]=muflux[ix]*(0.5+0.75*(thisX[ix]/pwidth)-0.25*(thisX[ix]/pwidth)**3)
				#self.log().notice("muflux span {0} to {1}".format(numpy.amin(muflux),numpy.amax(muflux)))
				posflux=muflux*ampl
				# could add delta function for prompt signal?
				# frequency response and pulse shape
				posflux2=scipy.convolve(posflux,pulseshape,mode="same") # counts per bin for this detector, with asymmetry and convolution
				# dead time
				posflux3=posflux2[::subsample]/(1+posflux2[::subsample]*r["deadtime"]/frames/dT)
				#self.log().notice("posflux span {0} to {1}".format(numpy.amin(posflux),numpy.amax(posflux)))
				ws.dataY(i)[:]=numpy.random.poisson(posflux3)
				ws.dataE(i)[:]=numpy.sqrt(ws.dataY(i))
			
		self.setProperty("OutputWorkspace",ws)

	# functions for example data
	# input: x array and parameters
	# output: 0 and 90 deg spectra
	# if called with x=None, calculate and return a field distribution in some form, to be passed back later
	def sampleGaussOsc(self,x,field,width,width2,precalc):
		if (x is None):
			return None
		Omega=field*0.01355*2*math.pi
		envel=numpy.exp(-(x*width)**2)
		return (numpy.cos(x*Omega)*envel, numpy.sin(x*Omega)*envel,1.0)

	def sampleExpOsc(self,x,field,width,width2,precalc):
		if (x is None):
			return None
		Omega=field*0.01355*2*math.pi
		envel=numpy.exp(-(x*width))
		return (numpy.cos(x*Omega)*envel, numpy.sin(x*Omega)*envel,1.0)
	
	def sampleShallowDonor(self,x,field, width,width2,precalc):
		if (x is None):
			return None
		if(width2!=0):
			# anisotropic, width2=D
			# integral over 2D is uniform in cos(theta)
			# need cos^2(theta)
			yr=numpy.zeros([len(x)])
			yi=numpy.zeros([len(x)])
			for i in range(200):
				b1=(field*0.01355+width/2.0+width2*(((i+0.5)/200.0)**2-0.5))*2*math.pi
				b2=(field*0.01355-width/2.0-width2*(((i+0.5)/200.0)**2-0.5))*2*math.pi
				yr=yr+numpy.cos(b1*x)+numpy.cos(b2*x)
				yi=yi+numpy.sin(b1*x)+numpy.sin(b2*x)
			return (yr/400.0, yi/400.0, 1.0)
		else:
			b1=(field+width)*0.01355*2*math.pi
			b2=(field-width)*0.01355*2*math.pi
			yr=numpy.cos(b1*x)+numpy.cos(b2*x)
			yi=numpy.sin(b1*x)+numpy.sin(b2*x)
			return (yr/2.0, yi/2.0,1.0)

	def sampleFluxLineLattice(self,x,field, width, width2,precalc):
		if (x is None):
			self.log().notice("Precalculating the flux line lattice")
			# useful to pre-calculate the distribution
			precalc=numpy.zeros([50*50])
			# numpy array of frequencies (all amplitudes equal so not stored)
			# Fourier components b_G = exp(-xi^2*|G|^2/4)/(1+lambda^2*|G|^2)
			# lattice spacing a=sqrt( (sqrt(3)*Phi0 / (2*B_mean) ) (triangular)
			# B(r) = B_mean * sum_lattice (b_G exp(-i G.r))
			# width -> lambda (nm)
			# width2 -> correlation length (nm)
			# spacing
			Phi0_SI=2.067833E-15 # Wb (T.m2)
			Phi0 = Phi0_SI * 1.E18 * 1.E4 # into G. nm^2
			a=math.sqrt( (math.sqrt(3.0)*Phi0)/(2.0*field) )
			astar=2*math.pi/a
			xi_sq_4=width2**2/4.0
			lamsq=width**2
			# sum over 1/4 of a primitive cell (muon locations):
			for ix in range(50):
				for iy in range(50):
					rx=(ix+0.5)/100.0 # multiple of a
					ry=(iy+0.5)/100.0
					b=0
					# sum over enough Fourier components
					for ki in range(-20,21):
						for kj in range(-20,21):
							dot=(rx*ki+ry*kj+0.5*(rx*kj+ry*ki))*2*math.pi
							Gsqr=(ki**2+ki*kj+kj**2)*astar**2
							b_G=math.exp(-xi_sq_4*Gsqr)/(1.0+lamsq*Gsqr)
							b=b+b_G*math.cos(dot)
					f=b*field*0.01355*2*math.pi
					precalc[ix+iy*50]=f
			return precalc
		else:
			yr=numpy.zeros([len(x)])
			yi=numpy.zeros([len(x)])
			for i in range(len(precalc)):
				yr=yr+numpy.cos(x*precalc[i])
				yi=yi+numpy.sin(x*precalc[i])
		return (yr/len(precalc),yi/len(precalc),1.0)
		
	def sampleAntiferromagnet(self,x,field,width,width2,precalc):
		if (x is None):
			return None
		# AF powder sample, zero field
		# ignore supplied "field" parameter
		# internal field = width
		# incommensurate field modulation to +-width2
		if (width2>0):
			y=numpy.zeros([len(x)])
			for i in range(180):
				y=y+numpy.cos(x*((width+width2*math.cos((i+0.5)*math.pi/180.0))*0.01355*2*math.pi))
			y=y/180.0
		else:
			y=numpy.cos(x*(width*0.01355*2*math.pi))
		return (y, 0.0, y )
				
AlgorithmFactory.subscribe(SimulateMuonTFData)
