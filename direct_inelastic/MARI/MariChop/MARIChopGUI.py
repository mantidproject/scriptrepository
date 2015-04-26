import sys
from MARIChopUI import Ui_MainWindow
from PyQt4 import QtCore, uic,QtGui
import time as time
#from mantidplotpy import *
#import dgreduce
import inspect
import numpy
from mantidplot import *
from mantid import *
from mantid.simpleapi import *
import math
import MariChop as MarChop
#class MainWindow(QtGui.QMainWindow, Ui_MainWindow):

class MainWindow(QtGui.QMainWindow):

	def __init__(self, parent=None):
		QtGui.QMainWindow.__init__(self,parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		QtCore.QObject.connect(self.ui.calc, QtCore.SIGNAL("clicked()"), self.calc )
		QtCore.QObject.connect(self.ui.ei, QtCore.SIGNAL("returnPressed()"), self.ei)
		QtCore.QObject.connect(self.ui.actionMari, QtCore.SIGNAL("triggered()"), self.mari )
		QtCore.QObject.connect(self.ui.actionMaps, QtCore.SIGNAL("triggered()"), self.maps )
		QtCore.QObject.connect(self.ui.actionMerlin, QtCore.SIGNAL("triggered()"), self.merlin )
		QtCore.QObject.connect(self.ui.actionGadolinium, QtCore.SIGNAL("triggered()"), self.gd )
		QtCore.QObject.connect(self.ui.actionSloppy, QtCore.SIGNAL("triggered()"), self.sloppy )
		QtCore.QObject.connect(self.ui.actionA, QtCore.SIGNAL("triggered()"), self.a )
		QtCore.QObject.connect(self.ui.actionRType, QtCore.SIGNAL("triggered()"), self.r )
		
		
		
		self.inst=''
		self.chop=''
		self.ei=0.0
		
	def maps(self):
		self.inst='maps'
	def merlin(self):
		self.inst='merlin'
	def mari(self):
		self.inst='mari'
		print 'MARI selected'
	
	def gd(self):
		self.chop='g'
		print 'Gd chopper selected'
	def sloppy(self):
		self.chop='s'
		print 'Sloppy chopper selected'
	def a(self):
		self.chop='a'
		print 'A chopper selected'
	def r(self):
		self.chop='r'
		print 'Relaxed chopper selected'
		
	def ei(self):
		self.ei=float(self.ui.ei.text())
		
	def calc(self):
		self.ei=float(self.ui.ei.text())
		MarChop.setchoptype(self.inst,self.chop)
		speed=[]
		incidentFlux=[]
		percent=[]
		WkspNames=[]
		string='Calculation for '+self.inst+' '+self.chop+' chopper at '+self.ui.ei.text()+'meV'
		self.ui.disp.addItem(string)
		string_title = 'Frequency [Hz]    Flux [n/s/cm^2]  Resolution [meV]  Resolution[ % ei]'                    
		self.ui.disp.addItem(string_title)        
		for freq in range(50,650,50):
			van_el,van,flux=MarChop.calculate(self.ei,freq,all=True)
			string1 = '{0}\t\t{1:.2f}\t{2:.2f}\t{3:.2f}'.format(freq,float(flux),float(van_el),((van_el/self.ei)*100))
			self.ui.disp.addItem(string1)
			fac=.99
			en_lo=1
			eps_min=en_lo
			eps_max=fac*self.ei +(1-fac)*en_lo
			eeps=range(int(eps_min),int(eps_max),1)
			dat=list(van)
			dat.reverse()
			CreateWorkspace(OutputWorkspace='Energy transfer resolution at'+str(freq)+'Hz',DataX=eeps,DataY=dat,DataE=list(van*0),WorkspaceTitle='Resolution at'+str(self.ei)+'meV',UnitX='Energy Transfer [meV]',YUnitLabel='resolution [meV]')
			WkspNames.append('Energy transfer resolution at'+str(freq)+'Hz')
			speed.append(freq)
			incidentFlux.append(flux.item())
			percent.append((van_el/self.ei)*100)
			
		GroupWorkspaces(OutputWorkspace='EnergyTransfer',InputWorkspaces=WkspNames)
		CreateWorkspace(OutputWorkspace="Resolution",DataX=speed,DataY=percent,DataE=(percent*0),WorkspaceTitle=self.chop+' chopper at '+str(self.ei)+'meV Resolution in % of'+str(self.ei)+'meV',UnitX='Chopper Frequency',YUnitLabel='Elastic line resolution % of Ei')
		CreateWorkspace(OutputWorkspace="Flux",DataX=speed,DataY=incidentFlux,DataE=(percent*0),WorkspaceTitle='Flux for '+self.chop+' chopper at '+str(self.ei)+'meV',UnitX='Chopper Frequency',YUnitLabel='Incident Flux n/s/cm^-2')

def qapp():
	if QtGui.QApplication.instance():
		app = QtGui.QApplication.instance()
	else:
		app = QtGui.QApplication(sys.argv)
	return app
 
app = qapp()
reducer = MainWindow()
reducer.show()
app.exec_()
