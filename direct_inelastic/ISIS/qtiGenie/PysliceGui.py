
import sys
from PySliceUI2 import Ui_MainWindow
from PyQt4 import QtCore, uic,QtGui
#import MantidFramework 
#MantidFramework.mtd.initialise()
from mantid.simpleapi import *
#from DirectEnergyConversion import *
import time as time
#from mantidplotpy import *
#import dgreduce
import inspect
import numpy
from mantidplot import *
from mantid import *
from mantid.simpleapi import *
from PySlice2 import *
from qtiGenie import *
from autoEi import *
#class MainWindow(QtGui.QMainWindow, Ui_MainWindow):

class MainWindow(QtGui.QMainWindow):

	def __init__(self, parent=None):
		QtGui.QMainWindow.__init__(self,parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		QtCore.QObject.connect(self.ui.plot, QtCore.SIGNAL("clicked()"), self.plot )
		QtCore.QObject.connect(self.ui.Reduce, QtCore.SIGNAL("clicked()"), self.Reduce )
		QtCore.QObject.connect(self.ui.plotOver, QtCore.SIGNAL("clicked()"), self.plotOver )
		QtCore.QObject.connect(self.ui.calcProj, QtCore.SIGNAL("clicked()"), self.calcProj )
		QtCore.QObject.connect(self.ui.display, QtCore.SIGNAL("clicked()"), self.display )
		QtCore.QObject.connect(self.ui.refresh, QtCore.SIGNAL("clicked()"), self.refreshlist )
		QtCore.QObject.connect(self.ui.cutMin, QtCore.SIGNAL("returnPressed()"), self.cutMin)
		QtCore.QObject.connect(self.ui.cutMax, QtCore.SIGNAL("returnPressed()"), self.cutMax)
		QtCore.QObject.connect(self.ui.delta, QtCore.SIGNAL("returnPressed()"), self.delta)
		QtCore.QObject.connect(self.ui.intMin, QtCore.SIGNAL("returnPressed()"), self.intMin)
		QtCore.QObject.connect(self.ui.intMax, QtCore.SIGNAL("returnPressed()"), self.intMax)
		QtCore.QObject.connect(self.ui.WkspIn, QtCore.SIGNAL("activated()"), self.setwksp)
		QtCore.QObject.connect(self.ui.actionDensityOfstates, QtCore.SIGNAL("triggered()"), self.pdos )
		QtCore.QObject.connect(self.ui.actionBoseFactor, QtCore.SIGNAL("triggered()"), self.bose )
		QtCore.QObject.connect(self.ui.removeWksp, QtCore.SIGNAL("clicked()"), self.removeWksp )
		QtCore.QObject.connect(self.ui.removePlot, QtCore.SIGNAL("clicked()"), self.removePlot )
		
		#reduction UI commands
		
		QtCore.QObject.connect(self.ui.InstName, QtCore.SIGNAL("activated(int)"), self.setInst )
		QtCore.QObject.connect(self.ui.NormMethod, QtCore.SIGNAL("activated(int)"), self.setNormMethod )
		
		QtCore.QObject.connect(self.ui.BkgSwitch,  QtCore.SIGNAL("stateChanged(int)"), self.BkgSwitchOn)
		
		QtCore.QObject.connect(self.ui.AutoEiChbox,  QtCore.SIGNAL("stateChanged(int)"), self.AutoEiOn)
		QtCore.QObject.connect(self.ui.FixEi,  QtCore.SIGNAL("stateChanged(int)"), self.FixEiOn)
		QtCore.QObject.connect(self.ui.FixMonitorSpectrum,  QtCore.SIGNAL("stateChanged(int)"), self.MonitorSpec)
		
		QtCore.QObject.connect(self.ui.AbsNormSwitch,  QtCore.SIGNAL("stateChanged(int)"), self.AbsNormOn)
		QtCore.QObject.connect(self.ui.sumRuns,  QtCore.SIGNAL("stateChanged(int)"), self.sumRunsOn)
		
		
		
		self.smooth=0
		self.cutMin=0
		self.cutMax=0
		self.delta=0
		self.intMin=0
		self.intMax=0
		self.scalemin=''
		self.scalemax=''
		self.selectedWksp=''
		self.selectedFig=''
		self.data_dict={}
		self.masterFigureDict={}
		self.fignum=0
		self.currentWkspNames=[]
		
		self.ui.axis.insertItem(1,'|Q|')
		self.ui.axis.insertItem(2,'E')
		
		#stuff for the reduction UI
		self.inst=''
		self.BkgSwitch=False
		self.AutoEi=False
		self.FixEi=False
		self.EiVal=[]
		self.MonitorSpec=False
		self.MonitorSpectrumNumber=[]
		self.Normalisation='Current'
		self.CurrentMapFile=''
		
		self.DetVanNum=''
		self.runNum=''
		self.WkspOutName=''
		
		self.AbsNormStatus=False
		
		self.MonoDetVanNum=''
		self.MonorunNum=''
		self.sampleMass=[]
		self.RMM=[]
		self.sumRunsCheck=False
		self.BkgdRange=[]
		self.ui.mapfile.setText('mari_res2013')
		self.ui.MonSpecNumber.setText('2')
		self.shortname=''
		
		
		self.currentWkspNames=mtd.getObjectNames()
		iter=0
		for item in self.currentWkspNames: 
			self.ui.WkspIn.insertItem(iter,item)
			iter=iter+1
		
	
	def setInst(self):
	
		self.inst=str(self.ui.InstName.currentText())
		print self.inst 
		
		if self.inst=='Mari':
			self.shortname='MAR'
			self.exten='.raw'
			self.ui.mapfile.setText('mari_res2013')
	
	def setNormMethod(self):
	
		self.Normalisation=str(self.ui.NormMethod.currentText())
		print self.Normalisation
		
		
	def BkgSwitchOn(self):
	
		if self.BkgSwitch==False:
			self.BkgSwitch =True
			self.BkgdRange=[]
			print 'Background subtraction on'
		else:
			self.BkgSwitch =False
			self.BkgdRange=[]
			print 'Background subtraction off'
	
	def sumRunsOn(self):
	
		if self.sumRunsCheck==False:
			self.sumRunsCheck =True
			print 'Reduction will sum runs'
		else:
			self.sumRunsCheck =False
			print 'Reduction will not sum runs'
			
	def AutoEiOn(self):
	
		if self.AutoEi==False:
			self.AutoEi =True
			print 'AutoEi on'
		else:
			self.AutoEi =False
			self.EiVal=double(self.ui.EiGuess.text())
			print 'AutoEi off: Ei guess is', str(self.EiVal),' meV'
			
	
	def FixEiOn(self):
	
		if self.FixEi==False:
			self.FixEi =True
			self.EiVal=double(self.ui.EiGuess.text())
			print 'FixEi on set to :', str(self.EiVal),' meV'
		else:
			self.FixEi =False
			print 'FixEi off'
	
	
	def MonitorSpec(self):
	
		if self.MonitorSpec==False:
			self.MonitorSpec =True
			self.MonitorSpectrumNumber=int(self.ui.MonSpecNumber.text())
			print 'Monitor spectrum set to :', str(self.MonitorSpectrumNumber)
		else:
			self.MonitorSpec =False
			print 'Monitor spectrum set to default'
	
	def AbsNormOn(self):
	
		if self.AbsNormStatus==False:
			self.AbsNormStatus =True
			print 'Absolute Normalisation on'
		else:
			self.AbsNormStatus =False
			print 'Absolute Normalisation off'
	def getRuns(self):
		WB=str(self.ui.DetVanNum.text())
 		Run=str(self.ui.RunNum.text())
 		
 		if len(Run.split(','))>1:
 			#case for a comma separated list
 			#generate a list of strings
 			Runs=Run.split(',')
 		
 		if len(Run.split(':'))>1:
 			#case for a colon separated list
 			#generate a list of strings from range
 			tmp=Run.split(':')
 			Runs=range(int(tmp[0]),int(tmp[1])+1)
 		
 		if len(Run.split(','))==1 & len(Run.split(':'))==1:
			if Run=='current':
				Runs='current'
			else:

 				#for a single run return a 1 element list
 				tmp=int(Run)
 				Runs=[1]
 				Runs[0]=tmp

 		print Runs, type(Run)
 		return WB,Runs
 	
 	def Reduce(self):
 	
 		if self.AbsNormStatus==True:
 			self.ReduceArbsolute()
 		else:
 			self.ReduceArbitary()
 		
 		
 	def ReduceArbsolute(self):
 		#absolute normalisation will be preformed
 		#try:
		WB,Runs=self.getRuns()
		print Runs, type(Runs)
		
		iliad_setup(self.shortname)
		mapfile=str(self.ui.mapfile.text())
		cal_file=self.shortname+WB+self.exten
		
		monoWB=str(self.ui.DetVanNum_2.text())
 		MonoRun=str(self.ui.MonoRunNum.text())
		
		sampleMass=double(self.ui.sampleMass.text())
		RMMmass=double(self.ui.RMMmass.text())
		
		monovan_mapfile=mapfile
		
		if self.sumRunsCheck==True:
		#sum runs
			Load(Filename=Runs[0],OutputWorkspace='SummedWksp',Cache=r'Never',LoadLogFiles='0')
			SummedWksp=mtd['SummedWksp']
			for i in range(1,len(Runs)):
				Load(Filename=Runs[i],OutputWorkspace='tmp',Cache=r'Never',LoadLogFiles='0')
				tmp=mtd['tmp']
				SummedWksp=SummedWksp+tmp
			DeleteWorkspace('tmp')
			#overwrite the runs list with 1 entry corresponding to the summed data
			Runs=['SummedWksp']
		
		for Run in Runs:
		
			if self.AutoEi==True:
				ei,rebin_params,self.BkgdRange=autoEi(str(Run),BkgdGen=True,monspecin=int(self.ui.MonSpecNumber.text()))
				#set the absolute van integration range
				tmp=rebin_params.split(',')
				monovanreb=[tmp[0],tmp[2]]
				
			else:
				ei=self.EiVal
				deltaEi=(ei*.01)/20.0
				rebin_params=str(-ei*.5)+','+str(deltaEi)+','+str(ei*.95)
				tmp=rebin_params.split(',')
				monovanreb=[tmp[0],tmp[2]]
			
			
			if self.FixEi == True:
				if self.BkgSwitch==True:
					w1=iliad_abs(str(WB),str(Run),str(MonoRun),str(monoWB),RMMmass,sampleMass,ei,rebin_params,mapfile,monovan_mapfile,det_cal_file=cal_file,fixei=True,norm_method=self.Normalisation,save_format='',abs_units_van_range=monovanreb,background=True,bkgd_range=self.BkgdRange)
					RenameWorkspace(InputWorkspace='w1',OutputWorkspace=str(Run)+'reduced')
				else:
					w1=iliad_abs(WB,str(Run),MonoRun,monoWB,RMMmass,sampleMass,ei,rebin_params,mapfile,monovan_mapfile,det_cal_file=cal_file,background=False,fixei=True,norm_method=self.Normalisation,save_format='',abs_units_van_range=monovanreb)
					RenameWorkspace(InputWorkspace='w1',OutputWorkspace=str(Run)+'reduced')
			else:
				if self.BkgSwitch==True:
					w1=iliad_abs(str(WB),str(Run),str(MonoRun),str(monoWB),RMMmass,sampleMass,ei,rebin_params,mapfile,monovan_mapfile,det_cal_file=cal_file,norm_method=self.Normalisation,save_format='',abs_units_van_range=monovanreb,background=True,bkgd_range=self.BkgdRange)
					RenameWorkspace(InputWorkspace='w1',OutputWorkspace=str(Run)+'reduced')
				else:
					w1=iliad_abs(WB,str(Run),MonoRun,monoWB,RMMmass,sampleMass,ei,rebin_params,mapfile,monovan_mapfile,det_cal_file=cal_file,background=False,norm_method=self.Normalisation,save_format='',abs_units_van_range=monovanreb)
					RenameWorkspace(InputWorkspace='w1',OutputWorkspace=str(Run)+'reduced')
		
		
		WkspOut=str(self.ui.outputWksp.text())
		
		if len(WkspOut)>0:
			RenameWorkspace(InputWorkspace=str(Run)+'reduced',OutputWorkspace=str(self.ui.outputWksp.text()))
		else:
			RenameWorkspace(InputWorkspace=str(Run)+'reduced',OutputWorkspace='mar'+str(Run)+'reduced')
			
		String=':) Run '+str(Run)+' reduced with absolute normalisation successfully ei guess ='+str(ei)+'meV'+' data rebinned with Emin,deltaE,Emax of '+rebin_params
		self.ui.ReduceOutput.addItem(String)
		#except:
		#	String=':( Something is wrong, terribly wrong ' 
		#	self.ui.ReduceOutput.addItem(String)
 	
 		
 	def ReduceArbitary(self):
 	
 		#try:
		WB,Runs=self.getRuns()
		print Runs, type(Runs)
		
		if Runs=='current':
			#get live data
			StartLiveData(UpdateEvery='0',Instrument='MARI',AccumulationMethod="Replace", OutputWorkspace='CurrentRun')
			Runs=['CurrentRun']

		iliad_setup(self.shortname)
		mapfile=str(self.ui.mapfile.text())
		cal_file=self.shortname+WB+self.exten
		
		
		if self.sumRunsCheck==True:
		#sum runs
			Load(Filename=Runs[0],OutputWorkspace='SummedWksp',Cache=r'Never',LoadLogFiles='0')
			SummedWksp=mtd['SummedWksp']
			for i in range(1,len(Runs)):
				Load(Filename=Runs[i],OutputWorkspace='tmp',Cache=r'Never',LoadLogFiles='0')
				tmp=mtd['tmp']
				SummedWksp=SummedWksp+tmp
			DeleteWorkspace('tmp')
			#overwrite the runs list with 1 entry corresponding to the summed data
			Runs=['SummedWksp']
		
		for Run in Runs:
		
			if self.AutoEi==True:
				ei,rebin_params,self.BkgdRange=autoEi(str(Run),BkgdGen=True,monspecin=int(self.ui.MonSpecNumber.text())) 
			else:
				ei=self.EiVal
				deltaEi=(ei*.01)/20.0
				rebin_params=str(-ei*.5)+','+str(deltaEi)+','+str(ei*.95)	
		
			#w1=iliad(WB,Run,ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method='current')
			if self.FixEi == True:
				if self.BkgSwitch==True:
					w1=iliad(WB,str(Run),ei,rebin_params,mapfile,det_cal_file=cal_file,fixei=True,norm_method=self.Normalisation,save_format='',background=self.BkgSwitch,bkgd_range=self.BkgdRange)
					RenameWorkspace(InputWorkspace='w1',OutputWorkspace=str(Run)+'reduced')
				else:
					w1=iliad(WB,str(Run),ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method=self.Normalisation,background=False,fixei=True,save_format='')
					RenameWorkspace(InputWorkspace='w1',OutputWorkspace=str(Run)+'reduced')
			else:
				if self.BkgSwitch==True:
					w1=iliad(WB,str(Run),ei,rebin_params,mapfile,det_cal_file=cal_file,norm_method=self.Normalisation,save_format='',background=self.BkgSwitch,bkgd_range=self.BkgdRange)
					RenameWorkspace(InputWorkspace='w1',OutputWorkspace=str(Run)+'reduced')
				else:
					w1=iliad(WB,str(Run),ei,rebin_params,mapfile,det_cal_file=cal_file,background=False,norm_method=self.Normalisation,save_format='')
					RenameWorkspace(InputWorkspace='w1',OutputWorkspace=str(Run)+'reduced')
		
		WkspOut=str(self.ui.outputWksp.text())
		
		if len(WkspOut)>0:
			RenameWorkspace(InputWorkspace=str(Run)+'reduced',OutputWorkspace=str(self.ui.outputWksp.text()))
		else:
			RenameWorkspace(InputWorkspace=str(Run)+'reduced',OutputWorkspace='mar'+str(Run)+'reduced')
			
		String=':) Run '+str(Run)+' reduced successfully ei guess ='+str(ei)+'meV'+' data rebinned with Emin,deltaE,Emax of '+rebin_params
		self.ui.ReduceOutput.addItem(String)
		#except:
		#	String=':( Something is wrong, terribly wrong ' 
		#	self.ui.ReduceOutput.addItem(String)
			
 	def setwksp(self):
 		print 'high'
 	
 	def removeWksp(self):
 		tmp=self.ui.wkspList.selectedItems()[0].text()
 		row= self.ui.wkspList.currentRow()
 		self.ui.wkspList.takeItem(row)
 		Wksp = tmp.split(':')[0]
 		
 		tmpdata=self.data_dict.get(str(Wksp))
 		DeleteWorkspace(tmpdata.data)
 		DeleteWorkspace(tmpdata.data_disp)
 		del self.data_dict[str(Wksp)]
 	
 	def removePlot(self):
 		tmp=self.ui.FigureList.selectedItems()[0].text()
 		row= self.ui.FigureList.currentRow()
 		self.ui.FigureList.takeItem(row)
 		fig = tmp.split(':')[0]
 		del self.masterFigureDict[int(fig)]
 		
 	def bose(self):
 		temp=float(self.ui.Temperature.text())
 		
 		if temp >=0.0:
 			print 'Correct for Bose factor'
 			tmp=self.ui.wkspList.selectedItems()[0].text()
 			row= self.ui.wkspList.currentRow()
 			self.ui.wkspList.takeItem(row)
 			
			self.selectedWksp = tmp.split(':')[0]
			string=tmp+'Bose Factor Applied with T='+" %.2f" % temp+'K'
			self.ui.wkspList.insertItem(row,string)
			data=self.data_dict.get(str(self.selectedWksp))
			data.boseFac(temp)
 		else:
 			print 'Enter temperature in kelvin for run'
 			
 		
 	def pdos(self):
 		temp=float(self.ui.Temperature.text())
 		dbwfac=float(self.ui.DBWFac.text())
 		
 		if temp >=0.0:
 			print 'Calculate Density of States'
 			tmp=self.ui.wkspList.selectedItems()[0].text()
 			row= self.ui.wkspList.currentRow()
 			self.ui.wkspList.takeItem(row)
			self.selectedWksp = tmp.split(':')[0]
			string=tmp+'Density of states calculated with T='+" %.2f" % temp+'K'
			self.ui.wkspList.insertItem(row,string)
			data=self.data_dict.get(str(self.selectedWksp))
			data.gofe(temp,dbwfac)	
		else:
 			print 'Enter temperature in kelvin for run'
 			
 	def refreshlist(self): 
		self.ui.WkspIn.clear() 
		
 		self.currentWkspNames=mtd.getObjectNames()
 		
		iter=0
		#for item in self.currentWkspNames: 
		#	self.ui.WkspIn.insertItem(iter,item)
		#	iter=iter+1
		#this finds only workspaces with the appropriate units for calc projections to mod q.
		# the hack is to discount existing workspaces that are in e trans but have SQW in the name. 
		#if the data2D class is changed such that the names for transformed workspaces are different this bit of the GUI will break
		
		allNames=mtd.getObjectNames()
		for i in range(len(allNames)):
			tmp=mtd[allNames[i]]
			try:
				if tmp.getAxis(0).getUnit().caption() == 'Energy transfer' and tmp.getName().find('SQW')==-1:
					print allNames[i] +' has units of energy transfer'
					if tmp.getNumberHistograms() > 1:
						self.ui.WkspIn.insertItem(iter,allNames[i])
					iter=iter+1
			except:
				print allNames[i]  +' is not a matrix workspace'

		
		
		
	def cutMin(self):
		self.cutMin=self.ui.cutMin.text()
		#print self.cutMin
	def cutMax(self):
		self.cutMax=self.ui.cutMax.text()
		#print self.cutMax
	def delta(self):
		self.delta=self.ui.delta.text()
		#print self.delta
	def intMin(self):
		self.intMin=self.ui.intMin.text()
		#print self.intMin
	def intMax(self):
		self.intMax=self.ui.intMax.text()
		#print self.intMax
		
	def plot(self):
		
		#try:
		self.intMax=float(self.ui.intMax.text())
		self.intMin=float(self.ui.intMin.text())
		self.delta=float(self.ui.delta.text())
		self.cutMin=float(self.ui.cutMin.text())
		self.cutMax=float(self.ui.cutMax.text())
		self.scalemin=(self.ui.scalemin.text())
		self.scalemax=(self.ui.scalemax.text())
		
		#print self.cutMin
		#print self.cutMax
		#print self.delta
		#print self.intMin
		#print self.intMax
		#print self.selectedWksp
		try:
			dir=self.ui.axis.currentText()
			tmp=self.ui.wkspList.selectedItems()[0].text()
			self.selectedWksp = tmp.split(':')[0]
			data=self.data_dict.get(str(self.selectedWksp))
		except:
			print 'select workspace to plot from'
		if dir=='E':
			print 'Cut along E'
			self.fignum=self.fignum+1
			data.ECut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,shoelace=True)
			entry=data.figure_dict.get(data.fignum)
			
			fighandle=entry[1]
			if len(self.scalemax) and len(self.scalemin)>0:
				self.scaleY(fighandle)
			
			self.masterFigureDict.setdefault(self.fignum,entry)
			#entry[0]is title entry[1]is handle
			string=str(self.fignum)+':'+entry[0]
			self.ui.FigureList.addItem(string)
			
		if dir =='|Q|':
			print 'Cut along |Q|'
			self.fignum=self.fignum+1
			data.QCut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,shoelace=True)
			entry=data.figure_dict.get(data.fignum)
			fighandle=entry[1]
			if len(self.scalemax) and len(self.scalemin)>0:
				self.scaleY(fighandle)
			self.masterFigureDict.setdefault(self.fignum,entry)
			#entry[0]is title entry[1]is handle
			string=str(self.fignum)+':'+entry[0]
			self.ui.FigureList.addItem(string)
			
		#except:
		 #return
		 
	def plotOver(self):
		
		
		self.intMax=float(self.ui.intMax.text())
		self.intMin=float(self.ui.intMin.text())
		self.delta=float(self.ui.delta.text())
		self.cutMin=float(self.ui.cutMin.text())
		self.cutMax=float(self.ui.cutMax.text())
		self.scalemin=(self.ui.scalemin.text())
		self.scalemax=(self.ui.scalemax.text())
		#print self.cutMin
		#print self.cutMax
		#print self.delta
		#print self.intMin
		#print self.intMax
		#print self.selectedWksp
		dir=self.ui.axis.currentText()
		
		try:
			tmp=self.ui.wkspList.selectedItems()[0].text()
			self.selectedWksp = tmp.split(':')[0]
			data=self.data_dict.get(str(self.selectedWksp))
			tmp=self.ui.FigureList.selectedItems()[0].text()
			self.selectedFig = tmp.split(':')[0]
			#print self.selectedFig
			figHandle=self.masterFigureDict.get(int(self.selectedFig))[1]
			#print self.selectedFig
			#print self.masterFigureDict
			#print figHandle
		except:
			print 'Select both workspace to plot from and figure to plot onto'
		if dir=='E':
			print 'Plot over along E'
			
			data.ECut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,over=True,Handle=figHandle,shoelace=True)
			
			entry=data.figure_dict.get(data.fignum)
			fighandle=entry[1]
			if len(self.scalemax) and len(self.scalemin)>0:
				self.scaleY(fighandle)
			
			del self.masterFigureDict[int(self.selectedFig)]
			self.masterFigureDict.setdefault(int(self.selectedFig),entry)
			
		if dir =='|Q|':
			print 'Plot Over along |Q|'
			data.QCut(self.intMin,self.intMax,self.cutMin,self.delta,self.cutMax,over=True,Handle=figHandle,shoelace=True)
			
			entry=data.figure_dict.get(data.fignum)
			fighandle=entry[1]
			if len(self.scalemax) and len(self.scalemin)>0:
				self.scaleY(fighandle)
			del self.masterFigureDict[int(self.selectedFig)]
			self.masterFigureDict.setdefault(int(self.selectedFig),entry)
			
		
	def scaleY(self,handle):
		min=(self.ui.scalemin.text())
		max=(self.ui.scalemax.text())
		active_layer = handle.activeLayer()
		active_layer.setScale(0,float(min),float(max))
			 
	def display(self):
		smooth= self.ui.smooth.text()
		IMin= self.ui.DispMin.text()
		IMax= self.ui.DispMax.text()
		
		if len(IMax)>0 and len(IMin)==0:
			IMin='0'
		
		try:
			tmp=self.ui.wkspList.selectedItems()[0].text()
			self.selectedWksp = tmp.split(':')[0]
			#print self.selectedWksp
			data=self.data_dict.get(str(self.selectedWksp))
			
			
			if len(smooth)>0 and len(IMin)>0 and len(IMax)>0:
				data.smooth(int(smooth),int(IMin),int(IMax))
				
			if len(smooth)>0 and len(IMin)==0 and len(IMax)==0:
				data.smooth(int(smooth))
				
			if len(smooth)==0 and len(IMin)>0 and len(IMax)>0:	
				data.display(int(IMin),int(IMax))
				
			if len(smooth)==0 and len(IMin)==0 and len(IMax)==0:	
				data.display()
				
		except:
			print 'Select data to display'
				
	def calcProj(self):
	
		wkspName=str(self.ui.WkspIn.currentText())
		
		if self.data_dict.get(str(wkspName)) != None:
			print 'S(qw) data from ',wkspName,'exists, deleting existing instance'
			tmpdata=self.data_dict.get(str(wkspName))
 			DeleteWorkspace(tmpdata.data)
 			DeleteWorkspace(tmpdata.data_disp)
 			del self.data_dict[str(wkspName)]
 			#now deal with the displayed text
 			numLines=self.ui.wkspList.count()
 			for i in range(numLines):
 				tmp= self.ui.wkspList.item(i).text()
				wkspExists=tmp.split(':')[0]
				
				if wkspExists==wkspName:
					self.ui.wkspList.takeItem(i)
					
				
		data=mtd[str(self.ui.WkspIn.currentText())]
		numHist=data.getNumberHistograms()
		#find the max and min angles
		maxTheta=0.1
		for i in range(numHist):
			tmp=(data.detectorTwoTheta(data.getDetector(i)))
			if tmp > maxTheta:
				maxTheta=tmp
		
		minTheta=360
		for i in range(numHist):
			tmp=(data.detectorTwoTheta(data.getDetector(i)))
			if tmp < minTheta:
				minTheta=tmp
				
		#calculate extents in mod q for the data and generate rebin parameters
		
		ei= data.getRun().getLogData('Ei').value
		ki=sqrt(ei/2.07)
		kf=sqrt((ei-0.0)/2.07)
		Qx=ki*1-cos(minTheta)*kf
		
		Qy=-(sin(minTheta)*cos(0))*ki
		Qz=-(sin(minTheta)*sin(0))*ki
		Qmin=sqrt(Qx**2+Qy**2+Qz**2)
		#print Qx,Qy,Qz
		
		Qx=ki*1-cos(maxTheta)*kf
		Qy=-(sin(maxTheta)*cos(0))*ki
		Qz=-(sin(maxTheta)*sin(0))*ki
		Qmax=sqrt(Qx**2+Qy**2+Qz**2)
		
		#deltaQ=(Qmax-Qmin)/(numHist-50)
		#set the deltaQ to be the elastic line intrinsic resolution
		
		#approximate deltaQ given the number of histograms
		if not self.ui.DeltaQ.text():
			deltaQ=(Qmax-Qmin)/(numHist)
		else:
			deltaQ=float(self.ui.DeltaQ.text())
			
			
		text=str(self.ui.WkspIn.currentText())
		string1=text+': ei= '+" %.2f" % ei+'meV qmin= '+" %.2f" % Qmin+' qmax = '+" %.2f" % Qmax+' deltaQ = '+" %.2f" % deltaQ
		
		
		
		rebin=str(Qmin)+','+str(deltaQ)+','+str(Qmax)
		self.ui.wkspList.addItem(string1)
		#print self.ui.WkspIn.currentText()
		instring=str(self.ui.WkspIn.currentText())
		data=data2D(instring)
		#newVarName=instring+data
		#new_var_name='some_new_variable_name'
		#exec '%s = tmp'%new_var_name
		data.rebinproj(rebin)
		self.data_dict.setdefault(instring,data)
		
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
