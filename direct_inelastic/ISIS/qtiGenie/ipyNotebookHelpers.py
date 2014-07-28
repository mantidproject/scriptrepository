try:
  from mantidplot import *
except ImportError:
  pass
from mantid import *
from mantid.simpleapi import *
from mantid.kernel import funcreturns

# (NEW!) To communicate with Plotly's server, sign in with credentials file
import plotly.plotly as py  

# (NEW!) Useful Python/Plotly tools
import plotly.tools as tls   

# (NEW!) Graph objects to piece together plots
from plotly.graph_objs import Data, Layout, Figure

import numpy  # numpy for math functions and arrays

#import matplotlib
import matplotlib.pyplot as plt   # side-stepping mpl's backend
#%matplotlib inline

import mpld3


class mplFig:
#a package for data for matplotlib
	def __init__(self,data):
		
		if type(data) != list:
			if mtd.doesExist(data)==True:# extract data from workspace for plotting. else assume its a x,y,e list
				cutDataWorkSpace=mtd[data]
			
				if cutDataWorkSpace.getNumberHistograms() == 1: # data is 1D

					evals=cutDataWorkSpace.extractE()
					xvals=cutDataWorkSpace.extractX()
					yvals=cutDataWorkSpace.extractY()

					NewXvals=self.centerBins(xvals)
					yvals=yvals[0]
					evals=evals[0]
		
					self.data=[NewXvals,yvals,evals]
			
				if cutDataWorkSpace.getNumberHistograms() > 1: # data is 2D
					cutDataWorkSpace=ReplaceSpecialValues(cutDataWorkSpace,NaNValue='0',InfinityValue='0')
					y=cutDataWorkSpace.getAxis(1) #need a fix for this hard coded value
					
					ydata=y.extractValues()
					xdata=cutDataWorkSpace.extractX()
					zdata=cutDataWorkSpace.extractY()
					
					zshape=zdata.shape

					yshape=ydata.shape

					xshape=xdata.shape
					
					if zshape[1] != xshape[0]:
						xdata=self.centerBins(xdata[0])
						
					if zshape[0] != yshape[0]:
						ydata=self.centerBins(ydata[0])
						
					self.data=[xdata,ydata,zdata]
		
		else:
			#input data is already a list of xyz or xye
			print 'data is a list '
			self.data=data#x,y,e or ,x,y,zsignal
		
		self.init_params()

	def init_params(self):
		self.linecol=2
		self.figure_dict={}
		self.fignum=0
		self.antialiasing=True #antialiasing
		self.marker='o' #marker type
		self.linestyle='-' # line style
		self.color='red'#line color
		self.linewidth='1' #line width
		self.markerfacecolor='red' #marker face color
		self.markeredgecolor='black' #marker edge color
		self.markersize=3 # marker size
		self.markeredgesize=1 #marker edge size
		self.xlim='' #xlimits
		self.ylim='' #ylimits
		self.xlabel='X Axis' #xlabel
		self.ylabel='Y Axis' #ylabel
		self.grid=True # grid
		self.vmin=''
		self.vmax=''
		self.numBins=15




	def mplPlotLine(self):

		mpld3.enable_notebook()
		plt.errorbar(self.data[0],self.data[1],self.data[2],
		aa=self.antialiasing,
		marker=self.marker,
		ls=self.linestyle,
		c=self.color,
		lw=self.linewidth,
		mfc=self.markerfacecolor,
		mec=self.markeredgecolor,
		ms=self.markersize,
		mew=self.markeredgesize)

		plt.xlim(self.xlim[0],self.xlim[1])
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)
		plt.grid(self.grid)

	def centerBins(dataIn):

		newVals=( dataIn + numpy.roll(dataIn,-1) )/2 # calculate the bin center 
		newVals = numpy.delete(newVals,-1) # remove the last element which is junk

		return newVals
	
	def centerBins(self,dataIn):

		newVals=( dataIn + numpy.roll(dataIn,-1) )/2 # calculate the bin center 
		newVals = numpy.delete(newVals,-1) # remove the last element which is junk

		return newVals

	def mplPlot2D(self):	
	#newdata=Transpose(ww.data)
		mpld3.disable_notebook()
		from matplotlib.ticker import MaxNLocator
		from matplotlib.colors import Normalize
		z=self.data[2]
		#normLev=Normalize(self.vmin,self.vmax)
		#levels = MaxNLocator(nbins=self.numBins).tick_values(z.min(), z.max())
		z[z>self.vmax]=numpy.nan
		z[z<self.vmin]=numpy.nan
		#plt.pcolor(self.data[0],self.data[1],self.data[2],vmin=self.vmin, vmax=self.vmax)
		plt.pcolor(self.data[0],self.data[1],self.data[2],
		vmin=self.vmin, 
		vmax=self.vmax,
		)
		#plt.colorbar()



