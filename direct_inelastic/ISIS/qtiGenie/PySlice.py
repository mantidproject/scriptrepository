#from utils import *
from mantid.simpleapi import *
from mantid.kernel import funcreturns
from mantidplot import *

try:
    from DirectEnergyConversion import *
    import dgreduce    
except ImportError:
    from Direct.DirectEnergyConversion import *
    import Direct.dgreduce as dgreduce

import time as time

import inspect
import numpy



class data2D():
	global active_layer, graph
	def __init__(self):
		self.data=0
		self.init_params()
	
	def init_params(self):
		self.linecol=2
		

	def rebinproj(self,wksp_in,qbin):
		data=sqw(wksp_in,qbin)
		self.data=data
		
	def colour(self,col):
		active_layer.setCurveLineColor(1,col)
	
	def cutE(self,Qmin,Qmax):
		cut(self.data,Qmin,Qmax,along='e')
	
	def cutQ(self,Emin,Emax):
		cut(self.data,min,max,along='q')
	
	def plotOverE(self,Qmin,Qmax):
		cut(self.data,Qmin,Qmax,along='e',over=True)
	
	def cutQ(self,Emin,Emax):
		cut(self.data,min,max,along='q',over=True)
	
	def setX(self,xscale_min,xscale_max):
		setscale(xscale_min,xscale_max,'x')
	
	def setY(self,yscale_min,yscale_max):
		setscale(yscale_min,yscale_max,'y')
	
	def display(self,*args):
			
		dat=self.data.importMatrixWorkspace().plotGraph2D()
		#set z scale 
		if len(args)==1:
			ll=dat.activeLayer()
			ll.setAxisScale(Layer.Right,0,args[1])
		if len(args)==2:
			ll=g2d.activeLayer()
			ll.setAxisScale(Layer.Right,args[1],args[2])
			
def sqw(wksp_in,qbin):
	"""
	convert to SmodQw assume direct geom requires string of rebin parameters
	sqw(w1,'0,.1,12')
	"""
	try:
		n,r=funcreturns.lhs_info('both')
		wksp_out=r[0]
		ei= (wksp_in.getRun().getLogData("Ei").value)

		SofQW2(wksp_in, OutputWorkspace=wksp_out,QAxisBinning=qbin,EMode="Direct",EFixed=str(ei))
		Transpose(InputWorkspace=wksp_out,OutputWorkspace=wksp_out)
		return mtd[wksp_out]
	except:
		print 'no output workpsace defined'

	
def cut(wksp,min,max,**kwargs):
	'''
	Take a cut out of a wksp created by sqw
	assumes |Q| is X and Energy is Y
	keywords
	along=q
	along=e
	'''
	global active_layer, graph
	
	if kwargs.has_key('along') and kwargs.get('along')=='e' or kwargs.has_key('along') and kwargs.get('along')=='E':
	
		# get the x vector which in principle is Q
		xvec=wksp.readX(0)
		#return elements in the region of interest
		#there must be a more elegant way than this!
		minvec=numpy.nonzero(xvec>=min)
		maxvec=numpy.nonzero(xvec<=max)
		minel=minvec[0][0]
		maxel=maxvec[0][numpy.size(maxvec)-1]
		
		numpix=maxel-minel
		out=numpy.zeros(wksp.getNumberHistograms())
		outerr=numpy.zeros(wksp.getNumberHistograms())
		print numpix
		for i in range(wksp.getNumberHistograms()):
			out[i]=numpy.sum(wksp.readY(i)[minel:maxel])/numpix
			
			
			outerr[i]=(numpy.sqrt(numpy.sum(wksp.readE(i)[minel:maxel])**2))/numpix
			
		tmp=transpose(wksp)
		xvec=tmp.readX(0)
		#put cut data into mantid wksp
		CreateWorkspace(OutputWorkspace='Cut_from_'+str(wksp),DataX=list(xvec),DataY=list(out),DataE=list(outerr),VerticalAxisValues="Energy Transfer",WorkspaceTitle='Cut from '+str(wksp)+' between '+str(min)+' and '+str(max)+' invA')
		#plotSpectrum('Cut_from_'+str(wksp), 1,'true')
		label='Cut-from-'+str(wksp)
		table=fillqtitable(wksp,xvec,out,outerr,label)
		labely="Intensity arb. units"
		labelx="Energy Transfer meV"
		title='Cut from '+str(wksp)+' between '+str(min)+' and '+str(max)+' invA'
		if kwargs.has_key('over') and kwargs.get('over')==True:
			
			plotOverFromQtiTable(table,labelx,labely,title)
		else:
			plotFromQtiTable(table,labelx,labely,title)
	
	if kwargs.has_key('along') and kwargs.get('along')=='q' or kwargs.has_key('along') and kwargs.get('along')=='Q':
		wksp=transpose(wksp)
		# get the x vector which in principle is Q
		xvec=wksp.readX(0)
		#return elements in the region of interest
		#there must be a more elegant way than this!
		minvec=numpy.nonzero(xvec>=min)
		maxvec=numpy.nonzero(xvec<=max)
		minel=minvec[0][0]
		maxel=maxvec[0][numpy.size(maxvec)-1]
		
		numpix=maxel-minel
		out=numpy.zeros(wksp.getNumberHistograms())
		outerr=numpy.zeros(wksp.getNumberHistograms())
		print numpix
		for i in range(wksp.getNumberHistograms()):
			out[i]=numpy.sum(wksp.readY(i)[minel:maxel])/numpix
			
			
			outerr[i]=(numpy.sqrt(numpy.sum(wksp.readE(i)[minel:maxel])**2))/numpix
			
		tmp=transpose(wksp)
		xvec=tmp.readX(0)
		
		#put cut data into mantid wksp
		CreateWorkspace(OutputWorkspace='Cut_from_'+str(wksp),DataX=list(xvec),DataY=list(out),DataE=list(outerr),VerticalAxisValues="Momentum Transfer",WorkspaceTitle='Cut from '+str(wksp)+' between '+str(min)+' and '+str(max)+' meV')
		#plotSpectrum('Cut_from_'+str(wksp), 1,'true')
		label='Cut-from-'+str(wksp)
		table=fillqtitable(wksp,xvec,out,outerr,label)
		labely="Intensity arb. units"
		labelx="Momentum Transfer |Q|"
		title='Cut from '+str(wksp)+' between '+str(min)+' and '+str(max)+' meV'
		if kwargs.has_key('over') and kwargs.get('over')==True:
			
			plotOverFromQtiTable(table,labelx,labely,title)
		else:
			plotFromQtiTable(table,labelx,labely,title)

def plotFromQtiTable(tb_in,labelx,labely,title):
	global active_layer, graph
	graph=plot(tb_in,(1,2,3),2)
	active_layer = graph.activeLayer()
	#l.setLineColor(1, 2)
	active_layer.setTitle(title)
	active_layer.setAxisTitle(0, labelx)
	active_layer.setAxisTitle(2, labely)
	#l.setScale(0,0,50)# set y scale
	#l.setScale(2,0,50)#set xscale

def plotOverFromQtiTable(tb_in,labelx,labely,title):
	global active_layer, graph
	active_layer.addCurves(tb_in,(1,2,3),2)
	#l.setLineColor(1, 2)
	#l.setScale(0,0,50)# set y scale
	#l.setScale(2,0,50)#set xscale

#def makeCurrent():
#	global active_layer, graph
#	gr=currentGraph()

def fillqtitable(wksp_in,xvec,datY,daterr,label):
	outdat=createqtiTable(label,wksp_in.getNumberHistograms()+1)
	jj=1
	for i in range(wksp_in.getNumberHistograms()-1):
		outdat.setCell(1,jj,xvec[jj])
		outdat.setCell(2,jj,datY[jj])
		outdat.setCell(3,jj,daterr[jj])
		jj=jj+1
	return outdat
def setscale(min,max,dir):
	if dir == 'x':
		active_layer.setScale(2,min,max)#set xscale
	if dir == 'y':
		active_layer.setScale(0,min,max)#set xscale

def transpose(wksp_in):
	"""
	transpose workspace
	"""
	n,r=funcreturns.lhs_info('both')
	wksp_out=r[0]
	Transpose(InputWorkspace=wksp_in,OutputWorkspace=wksp_out)
	return mtd[wksp_out]

def createqtiTable(*args):
#create a qti table of length arg1 with name arg0
	if len(args)==0:
		out=newTable()
	if len(args)==2:
		out=newTable(args[0],args[1],3)
		out.setColumnRole(3, 5)
	return out 
