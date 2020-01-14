from __future__ import print_function, unicode_literals

## Quantum - a program for solving spin evolution of the muon
## Author: James Lord
## Version 1.04, August 2018
import numpy
import math
#import time
from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import *
from mantid.geometry import *
#from collections import Counter
import re
import numbers
#import decimal
#import sys
import traceback
numpy.set_printoptions(linewidth=100)

from quantumtools import *
from quantumtabletools import *

# model strings for description
# parts separated by ;
# spins=mu,e,H uses look up table and populates I and gamma
# spins=1/2,1/2,3 sets I only
# gamma(0)=nnn
# axis(i)=1,1,1 not needed
# A(i)=nnn or A(i,j)=nnn if more than one e (or none). A(mu)=nnn also allowed, or a(e1,e2)
# D(i)=nnn
# E(i)=nnn
# or A(i,j)=aaa,ddd,[x,y,z],eee,[x2,y2,z2]
# etc (as implemented)

# orientation definition string
# "LFUniform=n"
# "LFRandom=n"
# "TFUniform=n"
# "TFRandom=n"
# "LF=1,0,0"
# "TF="1,0,0, 0,1,0"
# "Bmag=10234.0"

def ParseAndApplyModelString(string,B,beam,det):
	# return (spmat,Ham,rho0,scint)
	# assumes B known
	pass

def ParseAndApplyModelString_notZ(string):
	# return (spmat,Ham)
	# omits Zeeman component
	pass

def ParseAndApplyModelString_addZ(string,B,beam,det,Ham0):
	# return (spmat,Ham,rho0,scint) adding in Zeeman to Ham0
	pass

def ParseOrientationString(string):
	# iterator?
	pass


# def ParseStringToDict(st,di0={}): now in quantumtabletools

def ParseTableWorkspaceToDict(tw,di0={}):
	# column 0 is param names (text, inc index)
	# column 1 is values (ideally text, free)
	# implied "=" between columns
	# omit empty rows
	# no semicolon or double "=" lines possible now
	textlists=("measure","spins","detectspin","initialspin","tzero","recycle","brf","morespin","crystal","goniometer")
	freetext=("fitfunction","loop0par","loop1par","fit0par","fit1par","fit2par","fit3par","fit4par","fit5par","fit6par","fit7par","fit8par","fit9par") # do not parse at all (are on one line)
	di=dict(di0) # copy
	text = list(zip(list(map(str,tw.column(0))),list(map(str,tw.column(1)))))
	for (c1,c2) in text:
		if(c1 != "" and c1[0] != "#"):
			if(c1 in freetext):
				di[c1]=c2
			else:
				val=c2.split(",")
				kvl0s=c1.split("(")[0]
				if(kvl0s not in textlists):
					val=list(map(float,val))
				di[c1]=val
	return di

# def calc_eqm(rates): now in quantumtabletools

# def EnumerateSites(p,d,keystr): now in quantumtabletools
# def ParseAndIterateModel(pars): now in quantumtabletools

# def SelectOrientIterator(iterator,which): now in quantumtabletools

# def processor_TimeSpectra(pars,ybins,ebins,dest,asym): now in quantumtabletools
#def processor_TimeSpectraPQ(pars,ybins,ebins,dest,asym):
#def processor_IntegralAsym(pars,ybins,ebins,dest,asym):
#def processor_FittedCurve(pars,ybins,ebins,dest,asym):
#def tidyup_FittedCurve(pars,x0,x1):
#def processor_moments(pars,ybins,ebins,dest,bigomega,bigccos,bigcsin,numave):
#def processor_freqspec(pars,ybins,ebins,dest,bigomega,bigccos,bigcsin,numave):
#def processor_breitrabi(pars,ybins,ebins,dest,bigomega,bigccos,bigcsin,numave):

# def PreParseLoop(pars,hadaxis0,prog=None): now in quantumtabletools
# def ParseMeasureType(pars,prog=None): now in quantumtabletools
# def ParseMeasureType(pars,prog=None): now in quantumtabletools
# def ParseAndIterateLoop(pars): now in quantumtabletools

# def RunModelledSystem(pars0,prog=None): now in quantumtabletools
# def GetUserAxisName(axname): now in quantumtabletools
# def GetUserYAxisName(axname): now in quantumtabletools

class QuantumTableDrivenSimulation(PythonAlgorithm):
	def name(self):
		return 'SimulateBasedOnTable'
	def category(self):
		return 'Muon\Quantum'

	def PyInit(self):
		self.declareProperty(ITableWorkspaceProperty('ModelTable','',Direction.Input))
		self.declareProperty(WorkspaceProperty('Results','',Direction.Output))

	def PyExec(self):
		prog=Progress(self,0.0,1.0,100)
		tw=self.getProperty('ModelTable').value
		pars=ParseTableWorkspaceToDict(tw)
		axes,ybins,ebins = RunModelledSystem(pars,prog)
		if(len(axes)==1):
			# 1D result, create dummy Spectrum axis
			axes.append(["Simulation"])
		ns=len(axes[1])
		ws=WorkspaceFactory.create("Workspace2D",NVectors=ns,XLength=len(axes[0]),YLength=ybins.shape[1])
		ws.setDistribution(True) # always polarisation or similar normalised value, not raw or simulated counts
		for s in range(ns):
			ws.dataX(s)[:]=axes[0]
			#print "reduction???",ws.dataY(s)[:].shape,ybins[s,:].shape
			ws.dataY(s)[:]=ybins[s,:]
			ws.dataE(s)[:]=ebins[s,:]
		if(isinstance(axes[1][0], numbers.Number)):
			na = NumericAxis.create(ns)
			for i in range(ns):
				na.setValue(i,axes[1][i])
			#na.setUnit("TOF")
			#for i in range(ns):
			#	if(axes[1][i] != 0):
			#		e10=int(math.log10(abs(axes[1][i])))-3
			#		digits=decimal.Decimal((0,(1,),e10))
			#	else:
			#		e10=0
			#		digits=decimal.Decimal('0.0')
			#	print "rounding ",axes[1][i]," with ",e10,digits," to ",float(decimal.Decimal(axes[1][i]).quantize(digits))
			#	na.setValue(i,float(decimal.Decimal(axes[1][i]).quantize(digits))) # round to nearest decimal value to prevent formatting .9999999999 or .00000000001
			#y2=numpy.around(axes[1],decimals=6)
			#for i in range(ns):
			#	na.setValue(i,y2[i])
			#e1=1.1-sys.float_info.epsilon*10
			#for i in range(ns):
			#	na.setValue(i,e1)
			#	e1=e1+sys.float_info.epsilon
			yvar,yunit=GetUserAxisName(pars.get("axis1name"))
			if(yvar==""):
				na.setUnit(yunit)
			else:
				lbl=na.setUnit("Label")
				lbl.setLabel(yvar,yunit)
		else:
			na = TextAxis.create(ns)
			for i in range(ns):
				na.setLabel(i,str(axes[1][i]))
		ws.replaceAxis(1,na)
		xvar,xunit=GetUserAxisName(pars.get("axis0name"))
		if(xvar==""):
			ws.getAxis(0).setUnit(xunit)
		else:
			lbl=ws.getAxis(0).setUnit("Label")
			lbl.setLabel(xvar,xunit)
		ws.setYUnitLabel(GetUserYAxisName(pars.get("measure")[0]))
		self.setProperty('Results',ws)

AlgorithmFactory.subscribe(QuantumTableDrivenSimulation)

def InsertFitPars(pars,fps):
	# set fit parameters to fps = [fp0,fp1,...]
	for i in range(len(fps)):
		# find what it means
		fpnl=pars["fit"+str(i)+"par"].split(";")
		for fpn0 in fpnl:
			fpn=fpn0.strip()
			fp=fps[i]
			if(fpn[0]=="-"):
				fp=-fps[i]
				fpn=fpn[1:]
				#print "setting neg ",fpn," to ",fp
			#else:
				#print "setting pos ",fpn," to ",fp
			if(fpn[0]=="^"):
				fp=10.0**fps[i]
				fpn=fpn[1:]
			fpnp=fpn.split("[")
			if(len(fpnp)==1):
				# simple par=value, always overwrite
				#print "setting simple par ",lpnp[0],"=",x
				pars[fpnp[0]]=(fp,)
			else:
				# par[..]=value, need to loop up previous one and adjust
				fpi=fpnp[1].strip("[]")
				if fpnp[0] in pars:
					d2=pars[fpnp[0]][:]
				else:
					raise Exception("Scanning one component of "+fpnp[0]+" but rest is undefined")
				fpi2=fpi.split(",")
				if(len(fpi2)==1):
					d2[int(fpi2[0])]=fp
				elif(len(fpi2)==2): # 2 of 3 components allowing scanning of angle. Use for A axes, LF axis, etc
					d2[int(fpi2[0])]=math.cos(fp*math.pi/180.0)
					d2[int(fpi2[1])]=math.sin(fp*math.pi/180.0)
				else:
					raise Exception("Loop varying more than 2 components")
				#print "setting complex par ",lpnp[0],"=",d2
				pars[fpnp[0]]=d2

def RunModelledSystemCached(pars,key,cache={},counter=[0]):
	# keep cache of old "key" values and their results
	# in principle could key on most or all of "pars" but that may be too slow to check for a match, and some irrelevant stuff in there
	#content of cache is [ybins,lastTimeUsed]
	# don't store Axes and error values as they're not used in the fit anyway (but return if available)
	# pass None as progress indicator, as Fits do anyway
	counter[0]+=1
	if key in cache:
		pt=cache[key]
		pt[1]=counter[0]
		return None,pt[0],None
	else:
		if(len(cache)>200):
			# clear some cache entries, least recently used
			stale=sorted(list(cache.items()),key=lambda x: x[1][1])
			for obsolete in stale[0:50]:
				del cache[obsolete[0]]
		ax,yb,eb=RunModelledSystem(pars,None)
		cache[key]=[yb,counter[0]]
		return ax,yb,eb

class QuantumTableDrivenFunction1(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))
		
FunctionFactory.subscribe(QuantumTableDrivenFunction1)

class QuantumTableDrivenFunction2(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1",2.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			P2=self.getParameterValue("P1")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,P2))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,P2),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))
		
FunctionFactory.subscribe(QuantumTableDrivenFunction2)

class QuantumTableDrivenFunction3(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1",2.0)
		self.declareParameter("P2",3.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			P2=self.getParameterValue("P1")
			P3=self.getParameterValue("P2")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,P2,P3))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,P2,P3),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction3)

class QuantumTableDrivenFunction3SD(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1plusP2",2.0)
		self.declareParameter("P1minusP2",3.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P0=self.getParameterValue("P0")
			P1p2=self.getParameterValue("P1plusP2")
			P1m2=self.getParameterValue("P1minusP2")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P0,(P1p2+P1m2)/2,(P1p2-P1m2)/2))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P0,(P1p2+P1m2)/2,(P1p2-P1m2)/2),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction3SD)

class QuantumTableDrivenFunction4(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1",2.0)
		self.declareParameter("P2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	def setAttributeValue(self,name,value):
		if name == "Table":
			self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			P2=self.getParameterValue("P1")
			P3=self.getParameterValue("P2")
			P4=self.getParameterValue("P3")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,P2,P3,P4))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,P2,P3,P4),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction4)

class QuantumTableDrivenFunction4SD(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1plusP2",2.0)
		self.declareParameter("P1minusP2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P0=self.getParameterValue("P0")
			P1p2=self.getParameterValue("P1plusP2")
			P1m2=self.getParameterValue("P1minusP2")
			P3=self.getParameterValue("P3")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P0,(P1p2+P1m2)/2,(P1p2-P1m2)/2,P3))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])
			
			# cache control
			key=((P0,(P1p2+P1m2)/2,(P1p2-P1m2)/2,P3),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction4SD)

class QuantumTableDrivenFunction5(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1",2.0)
		self.declareParameter("P2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("P4",5.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			P2=self.getParameterValue("P1")
			P3=self.getParameterValue("P2")
			P4=self.getParameterValue("P3")
			P5=self.getParameterValue("P4")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,P2,P3,P4,P5))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,P2,P3,P4,P5),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction5)

class QuantumTableDrivenFunction5SD(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1plusP2",2.0)
		self.declareParameter("P1minusP2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("P4",4.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P0=self.getParameterValue("P0")
			P1p2=self.getParameterValue("P1plusP2")
			P1m2=self.getParameterValue("P1minusP2")
			P3=self.getParameterValue("P3")
			P4=self.getParameterValue("P4")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P0,(P1p2+P1m2)/2,(P1p2-P1m2)/2,P3,P4))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])
			
			# cache control
			key=((P0,(P1p2+P1m2)/2,(P1p2-P1m2)/2,P3,P4),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction5SD)

class QuantumTableDrivenFunction6(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1",2.0)
		self.declareParameter("P2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("P4",5.0)
		self.declareParameter("P5",6.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			P2=self.getParameterValue("P1")
			P3=self.getParameterValue("P2")
			P4=self.getParameterValue("P3")
			P5=self.getParameterValue("P4")
			P6=self.getParameterValue("P5")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,P2,P3,P4,P5,P6))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,P2,P3,P4,P5,P6),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction6)

class QuantumTableDrivenFunction6SD(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1plusP2",2.0)
		self.declareParameter("P1minusP2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("P4",4.0)
		self.declareParameter("P5",4.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P0=self.getParameterValue("P0")
			P1p2=self.getParameterValue("P1plusP2")
			P1m2=self.getParameterValue("P1minusP2")
			P3=self.getParameterValue("P3")
			P4=self.getParameterValue("P4")
			P5=self.getParameterValue("P5")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P0,(P1p2+P1m2)/2,(P1p2-P1m2)/2,P3,P4,P5))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])
			
			# cache control
			key=((P0,(P1p2+P1m2)/2,(P1p2-P1m2)/2,P3,P4,P5),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction6SD)

class QuantumTableDrivenFunction7(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1",2.0)
		self.declareParameter("P2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("P4",5.0)
		self.declareParameter("P5",6.0)
		self.declareParameter("P6",7.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			P2=self.getParameterValue("P1")
			P3=self.getParameterValue("P2")
			P4=self.getParameterValue("P3")
			P5=self.getParameterValue("P4")
			P6=self.getParameterValue("P5")
			P7=self.getParameterValue("P6")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,P2,P3,P4,P5,P6,P7))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,P2,P3,P4,P5,P6,P7),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction7)

class QuantumTableDrivenFunction8(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1",2.0)
		self.declareParameter("P2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("P4",5.0)
		self.declareParameter("P5",6.0)
		self.declareParameter("P6",7.0)
		self.declareParameter("P7",8.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			P2=self.getParameterValue("P1")
			P3=self.getParameterValue("P2")
			P4=self.getParameterValue("P3")
			P5=self.getParameterValue("P4")
			P6=self.getParameterValue("P5")
			P7=self.getParameterValue("P6")
			P8=self.getParameterValue("P7")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,P2,P3,P4,P5,P6,P7,P8))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,P2,P3,P4,P5,P6,P7,P8),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction8)

class QuantumTableDrivenFunction9(IFunction1D):
	def category(self):
		return 'Muon'

	def init(self):
		#self.declareAttribute("Table","Tab")
		self.declareParameter("P0",1.0)
		self.declareParameter("P1",2.0)
		self.declareParameter("P2",3.0)
		self.declareParameter("P3",4.0)
		self.declareParameter("P4",5.0)
		self.declareParameter("P5",6.0)
		self.declareParameter("P6",7.0)
		self.declareParameter("P7",8.0)
		self.declareParameter("P8",9.0)
		self.declareParameter("Scale",1.0) # always have variable asymmetry and background, could tie if not needed
		self.declareParameter("Baseline",0.0)

	#def setAttributeValue(self,name,value):
	#	if name == "Table":
	#		self._table = mtd[value]
			
	def function1D(self,xvals):
		try:
			P1=self.getParameterValue("P0")
			P2=self.getParameterValue("P1")
			P3=self.getParameterValue("P2")
			P4=self.getParameterValue("P3")
			P5=self.getParameterValue("P4")
			P6=self.getParameterValue("P5")
			P7=self.getParameterValue("P6")
			P8=self.getParameterValue("P7")
			P9=self.getParameterValue("P8")
			scal=self.getParameterValue("Scale")
			base=self.getParameterValue("Baseline")
			pars=ParseTableWorkspaceToDict(mtd["Tab"]) # self._table
			InsertFitPars(pars,(P1,P2,P3,P4,P5,P6,P7,P8,P9))
			# case 1: mtype="timespectra", X axis of data is time, copy X bins of data into time for simulation
			# case 2: mtype="integral", X axis of data means scan some parameter such as field. Copy into loop0.
			# either case, pre-define "axis0" with the data's X. 


			if(pars["measure"][0]=="timespectra"):
				xbins=numpy.zeros(len(xvals)+1)
				xbins[1:-1]=(xvals[1:]+xvals[:-1])/2.0 # inner bin boundaries
				xbins[0]=xbins[1]*2-xbins[2]
				xbins[-1]=xbins[-2]*2-xbins[-3] # bin boundaries, end bins set to same width as next ones in		
				pars["axis0"]=xbins
				pars["axis0extra"]=(1,)
				
			elif(pars["measure"][0]=="integral"):
				pars["axis0"]=numpy.array(xvals) # point data type for field, etc
				pars["axis0extra"]=(0,)

			else:
				raise Exception("Can't fit with measure type"+pars["measure"][0])

			# cache control
			key=((P1,P2,P3,P4,P5,P6,P7,P8,P9),len(xvals),xvals[0])
			axes,ybins,ebins = RunModelledSystemCached(pars,key)
			
			return ybins[0,:]*scal+base
		except Exception as e:
			traceback.print_exc()
			return numpy.array([1.E30]*len(xvals))

FunctionFactory.subscribe(QuantumTableDrivenFunction9)

class ReplaceFitParsInTable(PythonAlgorithm):
	def category(self):
		return 'Muon\Quantum'

	def PyInit(self):
		self.declareProperty(ITableWorkspaceProperty('OrigModelTable','',Direction.Input))
		self.declareProperty(ITableWorkspaceProperty('FitParamTable','',Direction.Input))
		self.declareProperty("Component","",doc='Which component of a Composite Function to use, if ambiguous')
		self.declareProperty(ITableWorkspaceProperty('UpdatedModelTable','',Direction.Output))

	def PyExec(self):
		Fittab=self.getProperty('OrigModelTable').value
		pars=ParseTableWorkspaceToDict(Fittab)
		FuncTab=self.getProperty('FitParamTable').value
		if(FuncTab.columnCount() != 3):
			raise Exception("Doesn't look like a fit parameter table")
		names=FuncTab.column(0)
		values=FuncTab.column(1)
		# look for a set of entries with common prefix and suffix P0,P1,..Pn,Scale,Baseline
		# group them first
		groupedNames={}
		for (i,name) in enumerate(names):
			splts=name.rsplit(".",1)
			if(len(splts)==1):
				prefix=''
				pname=splts[0]
			else:
				prefix,pname=splts
			if(prefix not in groupedNames):
				groupedNames[prefix]={}
			groupedNames[prefix][pname]=values[i]
		#print "processed parameter set is ",groupedNames
		comp=self.getProperty('Component').value
		if(comp):
			try:
				selGroup=groupedNames[comp]
			except:
				raise Exception("Component "+comp+" is not in this table")
		else:
			# search for the first one with right set of keys
			selGroup=None
			for g in list(groupedNames.values()):
				print("examining param group ",g)
				if('P0' in g):
					if(selGroup):
						raise Exception("More than one possible Quantum function found - please specify")
					else:
						selGroup=g
			if(not(selGroup)):
				raise Exception("Didn't find a Quantum fit function")
		pardict={}
		for (pnam,val) in list(selGroup.items()):
			if(pnam[0]=="P" and pnam[1].isdigit()):
				if(pnam[2:7]=="plusP"):
					pnum1=int(pnam[1])
					pnum2=int(pnam[7])
					if(pnum1 in pardict and pnum2 in pardict):
						pardict[pnum1]=pardict[pnum1]+float(val)/2.0
						pardict[pnum2]=pardict[pnum2]+float(val)/2.0
					else:
						pardict[pnum1]=float(val)/2.0
						pardict[pnum2]=float(val)/2.0
				elif(pnam[2:8]=="minusP"):
					pnum1=int(pnam[1])
					pnum2=int(pnam[8])
					if(pnum1 in pardict and pnum2 in pardict):
						pardict[pnum1]=pardict[pnum1]+float(val)/2.0
						pardict[pnum2]=pardict[pnum2]-float(val)/2.0
					else:
						pardict[pnum1]=float(val)/2.0
						pardict[pnum2]=-float(val)/2.0
				else:
					pnum=int(pnam[1:])
					pardict[pnum]=float(val)
		if(min(pardict.keys())==0 and max(pardict.keys())==len(pardict)-1):
			# good, conv to list
			fps=[float("NaN")]*len(pardict)
			for (i,v) in list(pardict.items()):
				fps[i]=v
		else:
			print(pardict)
			print(min(pardict.keys()),max(pardict.keys()),len(pardict))
			raise Exception("Some parameters seem to be missing")
		if(numpy.any(numpy.isnan(fps))):
			raise Exception("A parameter is NaN - not updating the table")
		if(numpy.any(numpy.isinf(fps))):
			raise Exception("A parameter is Inf - not updating the table") 
		InsertFitPars(pars,fps)
		# rewrite the table, based on the original one to keep the ordering
		# include any blank or comment lines
		nTab=WorkspaceFactory.createTable() # CreateEmptyTableWorkspace()
		nTab.addColumn("str","First")
		nTab.addColumn("str","Second")
		for i in range(Fittab.rowCount()):
			c1=Fittab.column(0)[i]
			c2=Fittab.column(1)[i]
			try:
				c2=pars[c1.strip()]
				if(isinstance(c2,list) or isinstance(c2,tuple)):
					c2=",".join(map(str,c2))
			except:
				pass
			nTab.addRow([c1,c2])
		self.setProperty("UpdatedModelTable",nTab)

AlgorithmFactory.subscribe(ReplaceFitParsInTable)
