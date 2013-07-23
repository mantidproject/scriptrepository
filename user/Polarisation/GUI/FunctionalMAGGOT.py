from PyQt4 import QtCore, QtGui
import os
from mantidsimple import *
from numpy import *

up='up'
down='down'
yes='y'
no='n'
#'hello \n')#

Hours=1
Minutes=60
Seconds=3600
first='1'

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MAGGOTWindow(object):
    def setupUi(self, MAGGOTWindow):
        MAGGOTWindow.setObjectName(_fromUtf8("MAGGOT"))
        MAGGOTWindow.resize(1089, 1001)
        self.centralwidget = QtGui.QWidget(MAGGOTWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
	
        self.MainTab = QtGui.QTabWidget(self.centralwidget)
        self.MainTab.setGeometry(QtCore.QRect(20, 40, 1041, 851))
        self.MainTab.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.MainTab.setObjectName(_fromUtf8("MainTab"))
	
        self.NMR = QtGui.QWidget()
        self.NMR.setObjectName(_fromUtf8("NMR"))
	
        self.File_Info = QtGui.QGroupBox(self.NMR)
        self.File_Info.setGeometry(QtCore.QRect(20, 10, 311, 111))
        self.File_Info.setObjectName(_fromUtf8("File_Info"))
	
        self.gridLayoutWidget_6 = QtGui.QWidget(self.File_Info)
        self.gridLayoutWidget_6.setGeometry(QtCore.QRect(10, 20, 291, 81))
        self.gridLayoutWidget_6.setObjectName(_fromUtf8("gridLayoutWidget_6"))
	
        self.File_Info_Grid = QtGui.QGridLayout(self.gridLayoutWidget_6)
        self.File_Info_Grid.setMargin(0)
        self.File_Info_Grid.setObjectName(_fromUtf8("File_Info_Grid"))
	
        self.NMRFileLocation = QtGui.QLabel(self.gridLayoutWidget_6)
        self.NMRFileLocation.setObjectName(_fromUtf8("NMRFileLocation"))
	
        self.File_Info_Grid.addWidget(self.NMRFileLocation, 0, 0, 1, 1)
        self.File_Info_Grid_2 = QtGui.QGridLayout()
        self.File_Info_Grid_2.setObjectName(_fromUtf8("File_Info_Grid_2"))
	
        self.browsebutton = QtGui.QPushButton(self.gridLayoutWidget_6)
        self.browsebutton.setObjectName(_fromUtf8("browsebutton"))
	#
	QtCore.QObject.connect(self.browsebutton,QtCore.SIGNAL("clicked()"),self.selectFile)
	#
        self.File_Info_Grid_2.addWidget(self.browsebutton, 0, 1, 1, 1)
        self.FilePath = QtGui.QLineEdit(self.gridLayoutWidget_6)
        self.FilePath.setObjectName(_fromUtf8("FilePath"))
        self.File_Info_Grid_2.addWidget(self.FilePath, 0, 0, 1, 1)
        self.File_Name = QtGui.QLabel(self.gridLayoutWidget_6)
        self.File_Name.setObjectName(_fromUtf8("File_Name"))
        self.File_Info_Grid_2.addWidget(self.File_Name, 1, 0, 1, 1)
        self.File_Name_2 = QtGui.QLineEdit(self.gridLayoutWidget_6)
        self.File_Name_2.setObjectName(_fromUtf8("File_Name_2"))
        self.File_Info_Grid_2.addWidget(self.File_Name_2, 1, 1, 1, 1)
        self.File_Info_Grid.addLayout(self.File_Info_Grid_2, 1, 0, 1, 1)
	###################################################
        self.SetUp_FIDMax = QtGui.QTabWidget(self.NMR)
        self.SetUp_FIDMax.setGeometry(QtCore.QRect(20, 130, 311, 111))
        self.SetUp_FIDMax.setObjectName(_fromUtf8("SetUp_FIDMax"))
        self.SetUp_FID = QtGui.QWidget()
        self.SetUp_FID.setObjectName(_fromUtf8("SetUp_FID"))
        self.gridLayoutWidget_5 = QtGui.QWidget(self.SetUp_FID)
        self.gridLayoutWidget_5.setGeometry(QtCore.QRect(10, 10, 281, 31))
        self.gridLayoutWidget_5.setObjectName(_fromUtf8("gridLayoutWidget_5"))
        self.SetUp_FID_Grid = QtGui.QGridLayout(self.gridLayoutWidget_5)
        self.SetUp_FID_Grid.setContentsMargins(0, 0, 0, 10)
        self.SetUp_FID_Grid.setObjectName(_fromUtf8("SetUp_FID_Grid"))
        self.SetUp_Freq_Val = QtGui.QSpinBox(self.gridLayoutWidget_5)
        self.SetUp_Freq_Val.setMaximum(9999)
        self.SetUp_Freq_Val.setSingleStep(1)
        self.SetUp_Freq_Val.setProperty(_fromUtf8("value"), 200)
        self.SetUp_Freq_Val.setObjectName(_fromUtf8("SetUp_Freq_Val"))
        self.SetUp_FID_Grid.addWidget(self.SetUp_Freq_Val, 0, 2, 1, 1)
        self.SetUp_Freq = QtGui.QLabel(self.gridLayoutWidget_5)
        self.SetUp_Freq.setObjectName(_fromUtf8("SetUp_Freq"))
	#
	#QtCore.QObject.connect(self.SetUp_Freq_Val,QtCore.SIGNAL("valueChanged(int)"),self.Frequency_F)
	#
	######################################################
        self.SetUp_FID_Grid.addWidget(self.SetUp_Freq, 0, 1, 1, 1)
        self.SetUp_FIDMax.addTab(self.SetUp_FID, _fromUtf8(""))
        self.SetUp_Max = QtGui.QWidget()
        self.SetUp_Max.setObjectName(_fromUtf8("SetUp_Max"))
        self.gridLayoutWidget_2 = QtGui.QWidget(self.SetUp_Max)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(20, 10, 261, 71))
        self.gridLayoutWidget_2.setObjectName(_fromUtf8("gridLayoutWidget_2"))
        self.SetUp_Max_Grid = QtGui.QGridLayout(self.gridLayoutWidget_2)
        self.SetUp_Max_Grid.setMargin(0)
        self.SetUp_Max_Grid.setObjectName(_fromUtf8("SetUp_Max_Grid"))
        self.Smoothing = QtGui.QCheckBox(self.gridLayoutWidget_2)
        self.Smoothing.setObjectName(_fromUtf8("Smoothing"))
	#
	#QtCore.QObject.connect( self.Smoothing, QtCore.SIGNAL("toggled(bool)"), self.Smoothing_F)
	#
        self.SetUp_Max_Grid.addWidget(self.Smoothing, 0, 0, 1, 2)
        self.Step_Tol = QtGui.QLabel(self.gridLayoutWidget_2)
        self.Step_Tol.setObjectName(_fromUtf8("Step_Tol"))
        self.SetUp_Max_Grid.addWidget(self.Step_Tol, 1, 0, 1, 1)
        self.A_Tol = QtGui.QLabel(self.gridLayoutWidget_2)
        self.A_Tol.setObjectName(_fromUtf8("A_Tol"))
        self.SetUp_Max_Grid.addWidget(self.A_Tol, 2, 0, 1, 1)
        self.Step_Tol_Val = QtGui.QDoubleSpinBox(self.gridLayoutWidget_2)
        self.Step_Tol_Val.setDecimals(5)
        self.Step_Tol_Val.setMaximum(1.0)
        self.Step_Tol_Val.setSingleStep(1e-05)
        self.Step_Tol_Val.setProperty(_fromUtf8("value"), 0.0001)
        self.Step_Tol_Val.setObjectName(_fromUtf8("Step_Tol_Val"))
	#
	#QtCore.QObject.connect(self.Step_Tol_Val,QtCore.SIGNAL("valueChanged(int)"),self.Step_Tol_F)
	#
        self.SetUp_Max_Grid.addWidget(self.Step_Tol_Val, 1, 1, 1, 1)
        self.A_Tol_Val = QtGui.QDoubleSpinBox(self.gridLayoutWidget_2)
        self.A_Tol_Val.setDecimals(4)
        self.A_Tol_Val.setMaximum(1.0)
        self.A_Tol_Val.setSingleStep(0.0001)
        self.A_Tol_Val.setProperty(_fromUtf8("value"), 0.003)
        self.A_Tol_Val.setObjectName(_fromUtf8("A_Tol_Val"))
	#
	#QtCore.QObject.connect(self.A_Tol_Val,QtCore.SIGNAL("valueChanged(int)"),self.A_Tol_F)
	#
        self.SetUp_Max_Grid.addWidget(self.A_Tol_Val, 2, 1, 1, 1)
        self.SetUp_FIDMax.addTab(self.SetUp_Max, _fromUtf8(""))
	######################################################
	#Single, Multi or LPF
	######################################################
        self.SingleorMulti = QtGui.QTabWidget(self.NMR)
        self.SingleorMulti.setGeometry(QtCore.QRect(20, 260, 311, 391))
        self.SingleorMulti.setObjectName(_fromUtf8("SingleorMulti"))
        self.SingleFit = QtGui.QWidget()
        self.SingleFit.setObjectName(_fromUtf8("SingleFit"))
        self.verticalLayoutWidget = QtGui.QWidget(self.SingleFit)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(10, 10, 281, 51))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.SingleGrid = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.SingleGrid.setMargin(0)
        self.SingleGrid.setObjectName(_fromUtf8("SingleGrid"))
        self.FileNo = QtGui.QLabel(self.verticalLayoutWidget)
        self.FileNo.setObjectName(_fromUtf8("FileNo"))
        self.SingleGrid.addWidget(self.FileNo)
        self.FileNo_2 = QtGui.QSpinBox(self.verticalLayoutWidget)
        self.FileNo_2.setMaximum(9999)
        self.FileNo_2.setSingleStep(1)
        self.FileNo_2.setProperty(_fromUtf8("value"), 1)
        self.FileNo_2.setObjectName(_fromUtf8("FileNo_2"))
	#
	#QtCore.QObject.connect(self.FileNo_2,QtCore.SIGNAL("valueChanged(int)"),self.FileNo_2_F)
	#
        self.SingleGrid.addWidget(self.FileNo_2)
        self.SF_FIDorMax = QtGui.QTabWidget(self.SingleFit)
        self.SF_FIDorMax.setGeometry(QtCore.QRect(10, 70, 281, 80))
        self.SF_FIDorMax.setObjectName(_fromUtf8("SF_FIDorMax"))
        self.SF_FID = QtGui.QWidget()
        self.SF_FID.setObjectName(_fromUtf8("SF_FID"))
        self.FitSF_FID = QtGui.QPushButton(self.SF_FID)
        self.FitSF_FID.setGeometry(QtCore.QRect(10, 10, 221, 23))
        self.FitSF_FID.setObjectName(_fromUtf8("FitSF_FID"))
        self.SF_FIDorMax.addTab(self.SF_FID, _fromUtf8(""))
	#
	QtCore.QObject.connect(self.FitSF_FID,QtCore.SIGNAL("clicked()"),self.MAGGOT_SingleFit)
	#
        self.SF_Max = QtGui.QWidget()
        self.SF_Max.setObjectName(_fromUtf8("SF_Max"))
        self.RunSF_Max = QtGui.QPushButton(self.SF_Max)
        self.RunSF_Max.setGeometry(QtCore.QRect(10, 10, 211, 23))
        self.RunSF_Max.setObjectName(_fromUtf8("RunSF_Max"))
	#
	#QtCore.QObject.connect(self.RunSF_Max,QtCore.SIGNAL("clicked()"),self.MAX_SingleFit)
	#
        self.SF_FIDorMax.addTab(self.SF_Max, _fromUtf8(""))
        self.SingleFitResults = QtGui.QGroupBox(self.SingleFit)
        self.SingleFitResults.setGeometry(QtCore.QRect(10, 160, 281, 131))
        self.SingleFitResults.setObjectName(_fromUtf8("SingleFitResults"))
        self.gridLayoutWidget_8 = QtGui.QWidget(self.SingleFitResults)
        self.gridLayoutWidget_8.setGeometry(QtCore.QRect(10, 20, 261, 100))
        self.gridLayoutWidget_8.setObjectName(_fromUtf8("gridLayoutWidget_8"))
        self.SingleFitResults_2 = QtGui.QGridLayout(self.gridLayoutWidget_8)
        self.SingleFitResults_2.setMargin(0)
        self.SingleFitResults_2.setObjectName(_fromUtf8("SingleFitResults_2"))
        self.SFResFreq = QtGui.QLabel(self.gridLayoutWidget_8)
        self.SFResFreq.setObjectName(_fromUtf8("SFResFreq"))
        self.SingleFitResults_2.addWidget(self.SFResFreq, 1, 0, 1, 1)
        self.T2 = QtGui.QLabel(self.gridLayoutWidget_8)
        self.T2.setObjectName(_fromUtf8("T2"))
        self.SingleFitResults_2.addWidget(self.T2, 2, 0, 1, 1)
        self.SFResFreq_2 = QtGui.QLineEdit(self.gridLayoutWidget_8)
	
        self.SFResFreq_2.setObjectName(_fromUtf8("SFResFreq_2"))
        self.SingleFitResults_2.addWidget(self.SFResFreq_2, 1, 2, 1, 1)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.SingleFitResults_2.addItem(spacerItem, 0, 1, 1, 1)
        self.SFAmplitude_2 = QtGui.QLineEdit(self.gridLayoutWidget_8)
	
        self.SFAmplitude_2.setObjectName(_fromUtf8("SFAmplitude_2"))
        self.SingleFitResults_2.addWidget(self.SFAmplitude_2, 0, 2, 1, 1)
        self.T2_2 = QtGui.QLineEdit(self.gridLayoutWidget_8)
	
        self.T2_2.setObjectName(_fromUtf8("T2_2"))
        self.SingleFitResults_2.addWidget(self.T2_2, 2, 2, 1, 1)
        self.SFAmplitude = QtGui.QLabel(self.gridLayoutWidget_8)
        self.SFAmplitude.setObjectName(_fromUtf8("SFAmplitude"))
        self.SingleFitResults_2.addWidget(self.SFAmplitude, 0, 0, 1, 1)
        self.Phase = QtGui.QLabel(self.gridLayoutWidget_8)
        self.Phase.setObjectName(_fromUtf8("Phase"))
        self.SingleFitResults_2.addWidget(self.Phase, 3, 0, 1, 1)
        self.Phase_2 = QtGui.QLineEdit(self.gridLayoutWidget_8)
	
        self.Phase_2.setObjectName(_fromUtf8("Phase_2"))
        self.SingleFitResults_2.addWidget(self.Phase_2, 3, 2, 1, 1)
	
	######################################################
        self.SingleorMulti.addTab(self.SingleFit, _fromUtf8(""))
        self.MultiFit = QtGui.QWidget()
        self.MultiFit.setObjectName(_fromUtf8("MultiFit"))
        self.gridLayoutWidget_3 = QtGui.QWidget(self.MultiFit)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(10, 10, 281, 151))
        self.gridLayoutWidget_3.setObjectName(_fromUtf8("gridLayoutWidget_3"))
        self.MultiGrid = QtGui.QGridLayout(self.gridLayoutWidget_3)
        self.MultiGrid.setMargin(0)
        self.MultiGrid.setObjectName(_fromUtf8("MultiGrid"))
        self.TimeDiff = QtGui.QLabel(self.gridLayoutWidget_3)
        self.TimeDiff.setObjectName(_fromUtf8("TimeDiff"))
        self.MultiGrid.addWidget(self.TimeDiff, 2, 1, 1, 1)
        self.DiagMode = QtGui.QCheckBox(self.gridLayoutWidget_3) #Diagnostic Mode Check Box
		
        self.DiagMode.setObjectName(_fromUtf8("DiagMode"))
        self.MultiGrid.addWidget(self.DiagMode, 6, 1, 1, 2)
        self.Units = QtGui.QComboBox(self.gridLayoutWidget_3) #Units selection
	
        self.Units.setObjectName(_fromUtf8("Units"))
        self.Units.addItem(_fromUtf8(""))
        self.Units.addItem(_fromUtf8(""))
        self.Units.addItem(_fromUtf8(""))
        self.MultiGrid.addWidget(self.Units, 2, 2, 1, 1)
        self.UP_Down_LPP = QtGui.QComboBox(self.gridLayoutWidget_3) #Fitting options
	
        self.UP_Down_LPP.setObjectName(_fromUtf8("UP_Down_LPP"))
        self.UP_Down_LPP.addItem(_fromUtf8(""))
        self.UP_Down_LPP.addItem(_fromUtf8(""))
        self.UP_Down_LPP.addItem(_fromUtf8(""))
        self.UP_Down_LPP.addItem(_fromUtf8(""))
        self.MultiGrid.addWidget(self.UP_Down_LPP, 1, 2, 1, 1)
        self.TimeBetweenNMR = QtGui.QDoubleSpinBox(self.gridLayoutWidget_3) #Time between pulses
	
        self.TimeBetweenNMR.setSingleStep(0.01)
        self.TimeBetweenNMR.setProperty(_fromUtf8("value"), 1.0)
        self.TimeBetweenNMR.setObjectName(_fromUtf8("TimeBetweenNMR"))
        self.MultiGrid.addWidget(self.TimeBetweenNMR, 3, 1, 1, 2)
        self.MultiGrid_3 = QtGui.QGridLayout()
        self.MultiGrid_3.setObjectName(_fromUtf8("MultiGrid_3"))
        self.InitialNMR_2 = QtGui.QLabel(self.gridLayoutWidget_3)
        self.InitialNMR_2.setObjectName(_fromUtf8("InitialNMR_2"))
        self.MultiGrid_3.addWidget(self.InitialNMR_2, 0, 0, 1, 1)
        self.FinalNMR_2 = QtGui.QLabel(self.gridLayoutWidget_3)
        self.FinalNMR_2.setToolTip(_fromUtf8(""))
        self.FinalNMR_2.setWhatsThis(_fromUtf8(""))
        self.FinalNMR_2.setObjectName(_fromUtf8("FinalNMR_2"))
        self.MultiGrid_3.addWidget(self.FinalNMR_2, 1, 0, 1, 1)
        self.InitialNMR = QtGui.QSpinBox(self.gridLayoutWidget_3) #first file to fit
	
        self.InitialNMR.setMaximum(99999)
        self.InitialNMR.setProperty(_fromUtf8("value"), 1)
        self.InitialNMR.setObjectName(_fromUtf8("InitialNMR"))
        self.MultiGrid_3.addWidget(self.InitialNMR, 0, 1, 1, 1)
        self.FinalNMR = QtGui.QSpinBox(self.gridLayoutWidget_3) #last file to fit
	
        self.FinalNMR.setMaximum(99999)
        self.FinalNMR.setProperty(_fromUtf8("value"), 10000)
        self.FinalNMR.setObjectName(_fromUtf8("FinalNMR"))
        self.MultiGrid_3.addWidget(self.FinalNMR, 1, 1, 1, 1)
        self.MultiGrid.addLayout(self.MultiGrid_3, 4, 1, 1, 2)
        self.Fit_to = QtGui.QLabel(self.gridLayoutWidget_3)
        self.Fit_to.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.Fit_to.setAlignment(QtCore.Qt.AlignCenter)
        self.Fit_to.setObjectName(_fromUtf8("Fit_to"))
        self.MultiGrid.addWidget(self.Fit_to, 1, 1, 1, 1)
        self.MF_FIDorMax = QtGui.QTabWidget(self.MultiFit)
        self.MF_FIDorMax.setGeometry(QtCore.QRect(10, 170, 281, 61))
        self.MF_FIDorMax.setObjectName(_fromUtf8("MF_FIDorMax"))
        self.MF_FID = QtGui.QWidget()
        self.MF_FID.setObjectName(_fromUtf8("MF_FID"))
        self.Fit_MF_FID = QtGui.QPushButton(self.MF_FID) #fit multi button
	#
	QtCore.QObject.connect(self.Fit_MF_FID,QtCore.SIGNAL("clicked()"),self.MAGGOT)
	#
        self.Fit_MF_FID.setGeometry(QtCore.QRect(10, 10, 211, 23))
        self.Fit_MF_FID.setObjectName(_fromUtf8("Fit_MF_FID"))
        self.MF_FIDorMax.addTab(self.MF_FID, _fromUtf8(""))
        self.MF_Max = QtGui.QWidget()
        self.MF_Max.setObjectName(_fromUtf8("MF_Max"))
        self.Run_MF_Max = QtGui.QPushButton(self.MF_Max) #run multi max
	
        self.Run_MF_Max.setGeometry(QtCore.QRect(10, 10, 211, 23))
        self.Run_MF_Max.setObjectName(_fromUtf8("Run_MF_Max"))
        self.MF_FIDorMax.addTab(self.MF_Max, _fromUtf8(""))
        self.MF_Results = QtGui.QGroupBox(self.MultiFit)
        self.MF_Results.setGeometry(QtCore.QRect(10, 240, 281, 111))
        self.MF_Results.setObjectName(_fromUtf8("MF_Results"))
        self.gridLayoutWidget = QtGui.QWidget(self.MF_Results)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 20, 261, 80))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.MF_Results_2 = QtGui.QGridLayout(self.gridLayoutWidget)
        self.MF_Results_2.setMargin(0)
        self.MF_Results_2.setObjectName(_fromUtf8("MF_Results_2"))
	#Results
        self.LPP_2 = QtGui.QLineEdit(self.gridLayoutWidget)
	
        self.LPP_2.setObjectName(_fromUtf8("LPP_2"))
        self.MF_Results_2.addWidget(self.LPP_2, 2, 2, 1, 1)
        self.MF_Amplitude_2 = QtGui.QLineEdit(self.gridLayoutWidget)
	
        self.MF_Amplitude_2.setObjectName(_fromUtf8("MF_Amplitude_2"))
        self.MF_Results_2.addWidget(self.MF_Amplitude_2, 0, 2, 1, 1)
        self.LPP = QtGui.QLabel(self.gridLayoutWidget)
        self.LPP.setObjectName(_fromUtf8("LPP"))
        self.MF_Results_2.addWidget(self.LPP, 2, 0, 1, 1)
        self.MF_Amplitude = QtGui.QLabel(self.gridLayoutWidget)
        self.MF_Amplitude.setObjectName(_fromUtf8("MF_Amplitude"))
        self.MF_Results_2.addWidget(self.MF_Amplitude, 0, 0, 1, 1)
        self.T1orTup_2 = QtGui.QLineEdit(self.gridLayoutWidget)
	
        self.T1orTup_2.setObjectName(_fromUtf8("T1orTup_2"))
        self.MF_Results_2.addWidget(self.T1orTup_2, 1, 2, 1, 1)
        self.T1orTup = QtGui.QLabel(self.gridLayoutWidget)
        self.T1orTup.setObjectName(_fromUtf8("T1orTup"))
        self.MF_Results_2.addWidget(self.T1orTup, 1, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.MF_Results_2.addItem(spacerItem1, 0, 1, 1, 1)
	
	############################################################
        self.SingleorMulti.addTab(self.MultiFit, _fromUtf8(""))
        self.LPF = QtGui.QWidget()
        self.LPF.setObjectName(_fromUtf8("LPF"))
        self.layoutWidget = QtGui.QWidget(self.LPF)
        self.layoutWidget.setGeometry(QtCore.QRect(20, 150, 271, 51))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.LPF_Grid = QtGui.QGridLayout(self.layoutWidget)
        self.LPF_Grid.setMargin(0)
        self.LPF_Grid.setObjectName(_fromUtf8("LPF_Grid"))
        self.FLips = QtGui.QLabel(self.layoutWidget)
        self.FLips.setObjectName(_fromUtf8("FLips"))
        self.LPF_Grid.addWidget(self.FLips, 0, 0, 1, 1)
        self.LPF_FinalNMR = QtGui.QSpinBox(self.layoutWidget) #number of flips
	
        self.LPF_FinalNMR.setMaximum(99999)
        self.LPF_FinalNMR.setProperty(_fromUtf8("value"), 200)
        self.LPF_FinalNMR.setObjectName(_fromUtf8("LPF_FinalNMR"))
        self.LPF_Grid.addWidget(self.LPF_FinalNMR, 0, 1, 1, 1)
        self.LPP_Value = QtGui.QDoubleSpinBox(self.layoutWidget)# loss per pulse
	
        self.LPP_Value.setDecimals(6)
        self.LPP_Value.setMinimum(0.0)
        self.LPP_Value.setMaximum(1.0)
        self.LPP_Value.setSingleStep(1e-06)
        self.LPP_Value.setProperty(_fromUtf8("value"), 0.0005)
        self.LPP_Value.setObjectName(_fromUtf8("LPP_Value"))
        self.LPF_Grid.addWidget(self.LPP_Value, 1, 1, 1, 1)
        self.Use_LPP = QtGui.QCheckBox(self.layoutWidget)#use lpp?
	
        self.Use_LPP.setObjectName(_fromUtf8("Use_LPP"))
        self.LPF_Grid.addWidget(self.Use_LPP, 1, 0, 1, 1)
        self.Flip_File_Info = QtGui.QGroupBox(self.LPF)
        self.Flip_File_Info.setGeometry(QtCore.QRect(10, 10, 291, 131))
        self.Flip_File_Info.setObjectName(_fromUtf8("Flip_File_Info"))
        self.gridLayoutWidget_7 = QtGui.QWidget(self.Flip_File_Info)
        self.gridLayoutWidget_7.setGeometry(QtCore.QRect(10, 20, 271, 98))
        self.gridLayoutWidget_7.setObjectName(_fromUtf8("gridLayoutWidget_7"))
        self.Flip_File_Info_2 = QtGui.QGridLayout(self.gridLayoutWidget_7)
        self.Flip_File_Info_2.setMargin(0)
        self.Flip_File_Info_2.setObjectName(_fromUtf8("Flip_File_Info_2"))
        self.NMRFileLocation_2 = QtGui.QLabel(self.gridLayoutWidget_7)
        self.NMRFileLocation_2.setObjectName(_fromUtf8("NMRFileLocation_2"))
        self.Flip_File_Info_2.addWidget(self.NMRFileLocation_2, 0, 0, 1, 1)
        self.Flip_File_Info_3 = QtGui.QGridLayout()
        self.Flip_File_Info_3.setObjectName(_fromUtf8("Flip_File_Info_3"))
        self.browsebutton_2 = QtGui.QPushButton(self.gridLayoutWidget_7)
	#
	QtCore.QObject.connect(self.browsebutton_2,QtCore.SIGNAL("clicked()"),self.LPFselectFile)
	#
        self.browsebutton_2.setObjectName(_fromUtf8("browsebutton_2"))
        self.Flip_File_Info_3.addWidget(self.browsebutton_2, 0, 1, 1, 1)
        self.PreFlipName_2 = QtGui.QLineEdit(self.gridLayoutWidget_7)
	
        self.PreFlipName_2.setObjectName(_fromUtf8("PreFlipName_2"))
        self.Flip_File_Info_3.addWidget(self.PreFlipName_2, 1, 1, 1, 1)
        self.PreFlipName = QtGui.QLabel(self.gridLayoutWidget_7)
        self.PreFlipName.setObjectName(_fromUtf8("PreFlipName"))
        self.Flip_File_Info_3.addWidget(self.PreFlipName, 1, 0, 1, 1)
        self.FilePath_2 = QtGui.QLineEdit(self.gridLayoutWidget_7)
	
        self.FilePath_2.setObjectName(_fromUtf8("FilePath_2"))
        self.Flip_File_Info_3.addWidget(self.FilePath_2, 0, 0, 1, 1)
        self.PostFlipName = QtGui.QLabel(self.gridLayoutWidget_7)
        self.PostFlipName.setObjectName(_fromUtf8("PostFlipName"))
        self.Flip_File_Info_3.addWidget(self.PostFlipName, 2, 0, 1, 1)
        self.PostFlipName_2 = QtGui.QLineEdit(self.gridLayoutWidget_7)
	
        self.PostFlipName_2.setObjectName(_fromUtf8("PostFlipName_2"))
        self.Flip_File_Info_3.addWidget(self.PostFlipName_2, 2, 1, 1, 1)
        self.Flip_File_Info_2.addLayout(self.Flip_File_Info_3, 1, 0, 1, 1)
        self.LPF_FIDorMax = QtGui.QTabWidget(self.LPF)
        self.LPF_FIDorMax.setGeometry(QtCore.QRect(10, 210, 281, 80))
        self.LPF_FIDorMax.setObjectName(_fromUtf8("LPF_FIDorMax"))
        self.LPF_FID = QtGui.QWidget()
        self.LPF_FID.setObjectName(_fromUtf8("LPF_FID"))
        self.Fit_LPF_FID = QtGui.QPushButton(self.LPF_FID)#fit button
	#
	QtCore.QObject.connect(self.Fit_LPF_FID,QtCore.SIGNAL("clicked()"),self.LPF_f)
	#
        self.Fit_LPF_FID.setGeometry(QtCore.QRect(10, 10, 211, 23))
        self.Fit_LPF_FID.setObjectName(_fromUtf8("Fit_LPF_FID"))
        self.LPF_FIDorMax.addTab(self.LPF_FID, _fromUtf8(""))
        self.LPF_Max = QtGui.QWidget()
        self.LPF_Max.setObjectName(_fromUtf8("LPF_Max"))
        self.Run_LPF_Max = QtGui.QPushButton(self.LPF_Max)#run max
	
        self.Run_LPF_Max.setGeometry(QtCore.QRect(10, 10, 211, 23))
        self.Run_LPF_Max.setObjectName(_fromUtf8("Run_LPF_Max"))
        self.LPF_FIDorMax.addTab(self.LPF_Max, _fromUtf8(""))
        self.LPF_Results = QtGui.QGroupBox(self.LPF)
        self.LPF_Results.setGeometry(QtCore.QRect(10, 290, 281, 61))
        self.LPF_Results.setObjectName(_fromUtf8("LPF_Results"))
        self.gridLayoutWidget_9 = QtGui.QWidget(self.LPF_Results)
        self.gridLayoutWidget_9.setGeometry(QtCore.QRect(10, 20, 261, 31))
        self.gridLayoutWidget_9.setObjectName(_fromUtf8("gridLayoutWidget_9"))
        self.LPF_Results_2 = QtGui.QGridLayout(self.gridLayoutWidget_9)
        self.LPF_Results_2.setMargin(0)
        self.LPF_Results_2.setObjectName(_fromUtf8("LPF_Results_2"))
        self.LPF_Result = QtGui.QLineEdit(self.gridLayoutWidget_9)#loss per flip result
	
        self.LPF_Result.setObjectName(_fromUtf8("LPF_Result"))
        self.LPF_Results_2.addWidget(self.LPF_Result, 0, 2, 1, 1)
        self.LPF_Res_Label = QtGui.QLabel(self.gridLayoutWidget_9)
        self.LPF_Res_Label.setObjectName(_fromUtf8("LPF_Res_Label"))
        self.LPF_Results_2.addWidget(self.LPF_Res_Label, 0, 0, 1, 1)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.LPF_Results_2.addItem(spacerItem2, 0, 1, 1, 1)
	
	########################################################################
        self.SingleorMulti.addTab(self.LPF, _fromUtf8(""))
        self.MainTab.addTab(self.NMR, _fromUtf8(""))
        self.NeutronData = QtGui.QWidget()
        self.NeutronData.setObjectName(_fromUtf8("NeutronData"))
        self.Neutron = QtGui.QTabWidget(self.NeutronData)
        self.Neutron.setGeometry(QtCore.QRect(0, 0, 1031, 301))
        self.Neutron.setObjectName(_fromUtf8("Neutron"))
        self.Pol = QtGui.QWidget()
        self.Pol.setObjectName(_fromUtf8("Pol"))
        self.gridLayoutWidget_4 = QtGui.QWidget(self.Pol)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 10, 178, 171))
        self.gridLayoutWidget_4.setObjectName(_fromUtf8("gridLayoutWidget_4"))
        self.Pol_Grid = QtGui.QGridLayout(self.gridLayoutWidget_4)
        self.Pol_Grid.setMargin(0)
        self.Pol_Grid.setObjectName(_fromUtf8("Pol_Grid"))
        self.WaveL = QtGui.QLabel(self.gridLayoutWidget_4)
        self.WaveL.setObjectName(_fromUtf8("WaveL"))
        self.Pol_Grid.addWidget(self.WaveL, 3, 0, 1, 1)
        self.PL_Val = QtGui.QDoubleSpinBox(self.gridLayoutWidget_4)
	
        self.PL_Val.setDecimals(3)
        self.PL_Val.setSingleStep(0.001)
        self.PL_Val.setProperty(_fromUtf8("value"), 10.0)
        self.PL_Val.setObjectName(_fromUtf8("PL_Val"))
        self.Pol_Grid.addWidget(self.PL_Val, 1, 0, 1, 2)
        self.CalibrationNeut = QtGui.QLabel(self.gridLayoutWidget_4)
        self.CalibrationNeut.setObjectName(_fromUtf8("CalibrationNeut"))
        self.Pol_Grid.addWidget(self.CalibrationNeut, 4, 0, 1, 1)
        self.InitialNeut = QtGui.QLabel(self.gridLayoutWidget_4)
        self.InitialNeut.setObjectName(_fromUtf8("InitialNeut"))
        self.Pol_Grid.addWidget(self.InitialNeut, 5, 0, 1, 1)
        self.CalibrationNeut_No = QtGui.QSpinBox(self.gridLayoutWidget_4)
        self.CalibrationNeut_No.setMaximum(99999)
        self.CalibrationNeut_No.setObjectName(_fromUtf8("CalibrationNeut_No"))
        self.Pol_Grid.addWidget(self.CalibrationNeut_No, 4, 1, 1, 1)
        self.PL = QtGui.QLabel(self.gridLayoutWidget_4)
        self.PL.setObjectName(_fromUtf8("PL"))
        self.Pol_Grid.addWidget(self.PL, 0, 0, 1, 2)
        self.InitialNeut_No = QtGui.QSpinBox(self.gridLayoutWidget_4)
        self.InitialNeut_No.setMaximum(99999)
        self.InitialNeut_No.setObjectName(_fromUtf8("InitialNeut_No"))
        self.Pol_Grid.addWidget(self.InitialNeut_No, 5, 1, 1, 1)
        self.WaveL_Val = QtGui.QDoubleSpinBox(self.gridLayoutWidget_4)
        self.WaveL_Val.setMaximum(30.0)
        self.WaveL_Val.setSingleStep(0.01)
        self.WaveL_Val.setProperty(_fromUtf8("value"), 5.0)
        self.WaveL_Val.setObjectName(_fromUtf8("WaveL_Val"))
        self.Pol_Grid.addWidget(self.WaveL_Val, 3, 1, 1, 1)
        self.FinalNeut_No = QtGui.QSpinBox(self.gridLayoutWidget_4)
        self.FinalNeut_No.setMaximum(99999)
        self.FinalNeut_No.setObjectName(_fromUtf8("FinalNeut_No"))
        self.Pol_Grid.addWidget(self.FinalNeut_No, 6, 1, 1, 1)
        self.FinalNeut = QtGui.QLabel(self.gridLayoutWidget_4)
        self.FinalNeut.setObjectName(_fromUtf8("FinalNeut"))
        self.Pol_Grid.addWidget(self.FinalNeut, 6, 0, 1, 1)
        self.NeutTUnits = QtGui.QLabel(self.gridLayoutWidget_4)
        self.NeutTUnits.setObjectName(_fromUtf8("NeutTUnits"))
        self.Pol_Grid.addWidget(self.NeutTUnits, 7, 0, 1, 1)
        self.UnitSelect = QtGui.QComboBox(self.gridLayoutWidget_4)
        self.UnitSelect.setModelColumn(0)
        self.UnitSelect.setObjectName(_fromUtf8("UnitSelect"))
        self.UnitSelect.addItem(_fromUtf8(""))
        self.UnitSelect.addItem(_fromUtf8(""))
        self.UnitSelect.addItem(_fromUtf8(""))
        self.UnitSelect.addItem(_fromUtf8(""))
        self.Pol_Grid.addWidget(self.UnitSelect, 7, 1, 1, 1)
        self.groupBox_2 = QtGui.QGroupBox(self.Pol)
        self.groupBox_2.setGeometry(QtCore.QRect(200, 10, 181, 211))
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.checkBox_3 = QtGui.QCheckBox(self.groupBox_2)
        self.checkBox_3.setGeometry(QtCore.QRect(10, 20, 131, 17))
        self.checkBox_3.setObjectName(_fromUtf8("checkBox_3"))
        self.Neutron.addTab(self.Pol, _fromUtf8(""))
        self.tab_12 = QtGui.QWidget()
        self.tab_12.setObjectName(_fromUtf8("tab_12"))
        self.Neutron.addTab(self.tab_12, _fromUtf8(""))
        self.tab_11 = QtGui.QWidget()
        self.tab_11.setObjectName(_fromUtf8("tab_11"))
        self.Neutron.addTab(self.tab_11, _fromUtf8(""))
        self.MainTab.addTab(self.NeutronData, _fromUtf8(""))
        MAGGOTWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MAGGOTWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1089, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MAGGOTWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MAGGOTWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MAGGOTWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MAGGOTWindow)
        self.MainTab.setCurrentIndex(0)
        self.SetUp_FIDMax.setCurrentIndex(1)
        self.SingleorMulti.setCurrentIndex(2)
        self.SF_FIDorMax.setCurrentIndex(1)
        self.MF_FIDorMax.setCurrentIndex(1)
        self.LPF_FIDorMax.setCurrentIndex(1)
        self.Neutron.setCurrentIndex(0)
        self.UnitSelect.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MAGGOTWindow)

    def retranslateUi(self, MAGGOTWindow):
        MAGGOTWindow.setWindowTitle(QtGui.QApplication.translate("MAGGOTWindow", "MAGGOTWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.File_Info.setTitle(QtGui.QApplication.translate("MAGGOTWindow", "File Information", None, QtGui.QApplication.UnicodeUTF8))
        self.NMRFileLocation.setText(QtGui.QApplication.translate("MAGGOTWindow", "NMR File Location", None, QtGui.QApplication.UnicodeUTF8))
        self.browsebutton.setText(QtGui.QApplication.translate("MAGGOTWindow", "Browse Files", None, QtGui.QApplication.UnicodeUTF8))
        self.File_Name.setText(QtGui.QApplication.translate("MAGGOTWindow", "File Name", None, QtGui.QApplication.UnicodeUTF8))
        self.SetUp_Freq.setText(QtGui.QApplication.translate("MAGGOTWindow", "Frequency", None, QtGui.QApplication.UnicodeUTF8))
        self.SetUp_FIDMax.setTabText(self.SetUp_FIDMax.indexOf(self.SetUp_FID), QtGui.QApplication.translate("MAGGOTWindow", "Setup FID NMR", None, QtGui.QApplication.UnicodeUTF8))
        self.Smoothing.setText(QtGui.QApplication.translate("MAGGOTWindow", "Use Smoothing?", None, QtGui.QApplication.UnicodeUTF8))
        self.Step_Tol.setText(QtGui.QApplication.translate("MAGGOTWindow", "Step Tolerance", None, QtGui.QApplication.UnicodeUTF8))
        self.A_Tol.setText(QtGui.QApplication.translate("MAGGOTWindow", "Amplitude Tolerance", None, QtGui.QApplication.UnicodeUTF8))
        self.SetUp_FIDMax.setTabText(self.SetUp_FIDMax.indexOf(self.SetUp_Max), QtGui.QApplication.translate("MAGGOTWindow", "Setup Max Value NMR", None, QtGui.QApplication.UnicodeUTF8))
        self.FileNo.setText(QtGui.QApplication.translate("MAGGOTWindow", "File Number", None, QtGui.QApplication.UnicodeUTF8))
        self.FitSF_FID.setText(QtGui.QApplication.translate("MAGGOTWindow", "Fit", None, QtGui.QApplication.UnicodeUTF8))
        self.SF_FIDorMax.setTabText(self.SF_FIDorMax.indexOf(self.SF_FID), QtGui.QApplication.translate("MAGGOTWindow", "Fit FID NMR", None, QtGui.QApplication.UnicodeUTF8))
        self.RunSF_Max.setText(QtGui.QApplication.translate("MAGGOTWindow", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.SF_FIDorMax.setTabText(self.SF_FIDorMax.indexOf(self.SF_Max), QtGui.QApplication.translate("MAGGOTWindow", "Find NMR Max", None, QtGui.QApplication.UnicodeUTF8))
        self.SingleFitResults.setTitle(QtGui.QApplication.translate("MAGGOTWindow", "Results", None, QtGui.QApplication.UnicodeUTF8))
        self.SFResFreq.setText(QtGui.QApplication.translate("MAGGOTWindow", "Frequency(Hz)", None, QtGui.QApplication.UnicodeUTF8))
        self.T2.setText(QtGui.QApplication.translate("MAGGOTWindow", "T2(s)", None, QtGui.QApplication.UnicodeUTF8))
        self.SFAmplitude.setText(QtGui.QApplication.translate("MAGGOTWindow", "Amplitude(V)", None, QtGui.QApplication.UnicodeUTF8))
        self.Phase.setText(QtGui.QApplication.translate("MAGGOTWindow", "Phase", None, QtGui.QApplication.UnicodeUTF8))
        self.SingleorMulti.setTabText(self.SingleorMulti.indexOf(self.SingleFit), QtGui.QApplication.translate("MAGGOTWindow", "Single-Fit", None, QtGui.QApplication.UnicodeUTF8))
        self.TimeDiff.setText(QtGui.QApplication.translate("MAGGOTWindow", "Measurement Gap", None, QtGui.QApplication.UnicodeUTF8))
        self.DiagMode.setText(QtGui.QApplication.translate("MAGGOTWindow", "Diagnostic Mode", None, QtGui.QApplication.UnicodeUTF8))
        self.Units.setItemText(0, QtGui.QApplication.translate("MAGGOTWindow", "Hours", None, QtGui.QApplication.UnicodeUTF8))
        self.Units.setItemText(1, QtGui.QApplication.translate("MAGGOTWindow", "Minutes", None, QtGui.QApplication.UnicodeUTF8))
        self.Units.setItemText(2, QtGui.QApplication.translate("MAGGOTWindow", "Seconds", None, QtGui.QApplication.UnicodeUTF8))
        self.UP_Down_LPP.setItemText(0, QtGui.QApplication.translate("MAGGOTWindow", "Spin Up", None, QtGui.QApplication.UnicodeUTF8))
        self.UP_Down_LPP.setItemText(1, QtGui.QApplication.translate("MAGGOTWindow", "T1 Decay", None, QtGui.QApplication.UnicodeUTF8))
        self.UP_Down_LPP.setItemText(2, QtGui.QApplication.translate("MAGGOTWindow", "Loss Per Pulse", None, QtGui.QApplication.UnicodeUTF8))
        self.UP_Down_LPP.setItemText(3, QtGui.QApplication.translate("MAGGOTWindow", "Nothing", None, QtGui.QApplication.UnicodeUTF8))
        self.InitialNMR_2.setText(QtGui.QApplication.translate("MAGGOTWindow", "Initial File", None, QtGui.QApplication.UnicodeUTF8))
        self.FinalNMR_2.setText(QtGui.QApplication.translate("MAGGOTWindow", "Final File", None, QtGui.QApplication.UnicodeUTF8))
        self.Fit_to.setText(QtGui.QApplication.translate("MAGGOTWindow", "Fit to", None, QtGui.QApplication.UnicodeUTF8))
        self.Fit_MF_FID.setText(QtGui.QApplication.translate("MAGGOTWindow", "Fit", None, QtGui.QApplication.UnicodeUTF8))
        self.MF_FIDorMax.setTabText(self.MF_FIDorMax.indexOf(self.MF_FID), QtGui.QApplication.translate("MAGGOTWindow", "Fit FID NMR", None, QtGui.QApplication.UnicodeUTF8))
        self.Run_MF_Max.setText(QtGui.QApplication.translate("MAGGOTWindow", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.MF_FIDorMax.setTabText(self.MF_FIDorMax.indexOf(self.MF_Max), QtGui.QApplication.translate("MAGGOTWindow", "Find NMR Max", None, QtGui.QApplication.UnicodeUTF8))
       
        self.MF_Results.setTitle(QtGui.QApplication.translate("MAGGOTWindow", "Results", None, QtGui.QApplication.UnicodeUTF8))
        self.LPP.setText(QtGui.QApplication.translate("MAGGOTWindow", "LPP(%)", None, QtGui.QApplication.UnicodeUTF8))
        self.MF_Amplitude.setText(QtGui.QApplication.translate("MAGGOTWindow", "Amplitude(V)", None, QtGui.QApplication.UnicodeUTF8))
        self.T1orTup.setText(QtGui.QApplication.translate("MAGGOTWindow", "T1/Tup(Time)", None, QtGui.QApplication.UnicodeUTF8))
        self.SingleorMulti.setTabText(self.SingleorMulti.indexOf(self.MultiFit), QtGui.QApplication.translate("MAGGOTWindow", "Multi-Fit", None, QtGui.QApplication.UnicodeUTF8))
        self.FLips.setText(QtGui.QApplication.translate("MAGGOTWindow", "Flips", None, QtGui.QApplication.UnicodeUTF8))
        self.Use_LPP.setText(QtGui.QApplication.translate("MAGGOTWindow", "Correct for LPP?", None, QtGui.QApplication.UnicodeUTF8))
        self.Flip_File_Info.setTitle(QtGui.QApplication.translate("MAGGOTWindow", "Flipping File Information", None, QtGui.QApplication.UnicodeUTF8))
        self.NMRFileLocation_2.setText(QtGui.QApplication.translate("MAGGOTWindow", "File Location", None, QtGui.QApplication.UnicodeUTF8))
        self.browsebutton_2.setText(QtGui.QApplication.translate("MAGGOTWindow", "Browse Files", None, QtGui.QApplication.UnicodeUTF8))
        self.PreFlipName.setText(QtGui.QApplication.translate("MAGGOTWindow", "Pre-Flip File Name", None, QtGui.QApplication.UnicodeUTF8))
        self.PostFlipName.setText(QtGui.QApplication.translate("MAGGOTWindow", "Post-Flip File Name", None, QtGui.QApplication.UnicodeUTF8))
        self.Fit_LPF_FID.setText(QtGui.QApplication.translate("MAGGOTWindow", "Fit", None, QtGui.QApplication.UnicodeUTF8))
        self.LPF_FIDorMax.setTabText(self.LPF_FIDorMax.indexOf(self.LPF_FID), QtGui.QApplication.translate("MAGGOTWindow", "Fit FID NMR", None, QtGui.QApplication.UnicodeUTF8))
        self.Run_LPF_Max.setText(QtGui.QApplication.translate("MAGGOTWindow", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.LPF_FIDorMax.setTabText(self.LPF_FIDorMax.indexOf(self.LPF_Max), QtGui.QApplication.translate("MAGGOTWindow", "Find NMR Max", None, QtGui.QApplication.UnicodeUTF8))
        self.LPF_Results.setTitle(QtGui.QApplication.translate("MAGGOTWindow", "Results", None, QtGui.QApplication.UnicodeUTF8))
        self.LPF_Res_Label.setText(QtGui.QApplication.translate("MAGGOTWindow", "LPF(%)", None, QtGui.QApplication.UnicodeUTF8))
        self.SingleorMulti.setTabText(self.SingleorMulti.indexOf(self.LPF), QtGui.QApplication.translate("MAGGOTWindow", "Loss Per Flip", None, QtGui.QApplication.UnicodeUTF8))
        self.MainTab.setTabText(self.MainTab.indexOf(self.NMR), QtGui.QApplication.translate("MAGGOTWindow", "NMR Data", None, QtGui.QApplication.UnicodeUTF8))
        self.WaveL.setText(QtGui.QApplication.translate("MAGGOTWindow", "Wavelength (Ang)", None, QtGui.QApplication.UnicodeUTF8))
        self.CalibrationNeut.setText(QtGui.QApplication.translate("MAGGOTWindow", "Calibration File", None, QtGui.QApplication.UnicodeUTF8))
        self.InitialNeut.setText(QtGui.QApplication.translate("MAGGOTWindow", "Initial File", None, QtGui.QApplication.UnicodeUTF8))
        self.PL.setText(QtGui.QApplication.translate("MAGGOTWindow", "Polariser Pressure-Length (bar.cm)", None, QtGui.QApplication.UnicodeUTF8))
        self.FinalNeut.setText(QtGui.QApplication.translate("MAGGOTWindow", "Final File", None, QtGui.QApplication.UnicodeUTF8))
        self.NeutTUnits.setText(QtGui.QApplication.translate("MAGGOTWindow", "Units", None, QtGui.QApplication.UnicodeUTF8))
        self.UnitSelect.setItemText(0, QtGui.QApplication.translate("MAGGOTWindow", "Hours", None, QtGui.QApplication.UnicodeUTF8))
        self.UnitSelect.setItemText(1, QtGui.QApplication.translate("MAGGOTWindow", "Minutes", None, QtGui.QApplication.UnicodeUTF8))
        self.UnitSelect.setItemText(2, QtGui.QApplication.translate("MAGGOTWindow", "Seconds", None, QtGui.QApplication.UnicodeUTF8))
        self.UnitSelect.setItemText(3, QtGui.QApplication.translate("MAGGOTWindow", "Days", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_2.setTitle(QtGui.QApplication.translate("MAGGOTWindow", " Analyser Cell", None, QtGui.QApplication.UnicodeUTF8))
        self.checkBox_3.setText(QtGui.QApplication.translate("MAGGOTWindow", "Include Analyser Cell?", None, QtGui.QApplication.UnicodeUTF8))
        self.Neutron.setTabText(self.Neutron.indexOf(self.Pol), QtGui.QApplication.translate("MAGGOTWindow", "Polarisation Measurements", None, QtGui.QApplication.UnicodeUTF8))
        self.Neutron.setTabText(self.Neutron.indexOf(self.tab_12), QtGui.QApplication.translate("MAGGOTWindow", "Cell Measurements", None, QtGui.QApplication.UnicodeUTF8))
        self.Neutron.setTabText(self.Neutron.indexOf(self.tab_11), QtGui.QApplication.translate("MAGGOTWindow", "Guide Efficiency", None, QtGui.QApplication.UnicodeUTF8))
        self.MainTab.setTabText(self.MainTab.indexOf(self.NeutronData), QtGui.QApplication.translate("MAGGOTWindow", "Neutron Data", None, QtGui.QApplication.UnicodeUTF8))

    def selectFile(self):
	
	#global filepath , filename
	dialog=QtGui.QFileDialog()
	try:
		PathandName = dialog.getOpenFileName(caption='Select NMR Files',directory=r'\\Britannic\3He\Cells')
	except:
		PathandName = dialog.getOpenFileName(caption='Select NMR Files',directory=r'C:/')
	#Extracts File Path
	for i in range(1,len(PathandName)):
		if PathandName[len(PathandName)-i] <> '/' :
			PathandName=PathandName
		else:
			break
	
	filepath=PathandName[0:(len(PathandName)-i+1)]
	#Extracts File Name
	for j in range(1,len(PathandName)):
		if PathandName[len(PathandName)-j] <> '_' :
			PathandName=PathandName
		else:
			break
	filename=PathandName[(len(PathandName)+1-i):(len(PathandName)+1-j)]
	
        self.FilePath.setText(filepath)
	self.File_Name_2.setText(filename)
	self.FileNo_2.setValue(int(PathandName[(len(PathandName)+1-j):len(PathandName)]))
	
    def LPFselectFile(self):
	
	#global filepath , filename
	dialog=QtGui.QFileDialog()
	try:
		PathandName = dialog.getOpenFileName(caption='Select NMR Files',directory=r'\\Britannic\3He\Cells')
	except:
		PathandName = dialog.getOpenFileName(caption='Select NMR Files',directory=r'C:/')
	#Extracts File Path
	for i in range(1,len(PathandName)):
		if PathandName[len(PathandName)-i] <> '/' :
			PathandName=PathandName
		else:
			break
	
	filepath=PathandName[0:(len(PathandName)-i+1)]
	#Extracts File Name
	for j in range(1,len(PathandName)):
		if PathandName[len(PathandName)-j] <> '_' :
			PathandName=PathandName
		else:
			break
	filename=PathandName[(len(PathandName)+1-i):(len(PathandName)+1-j)]
	
        self.FilePath_2.setText(filepath)
	self.PreFlipName_2.setText(filename)
	self.PostFlipName_2.setText(filename)
	
	
	#Extracts frequency
    #def Frequency_F(self):
	#global Freq
	#Freq = self.SetUp_Freq_Val.value()
	#self.File_Name_2.setText(str(Freq))
	
	#def test(self):
	#self.FilePath.setText(str(self.SetUp_Freq_Val.value())) #Can use results of other functions as arguments for other functions! supposedly...
	
	#Reads whether to use smoothing
  #  def Smoothing_F(self):
	    #if self.Smoothing.isChecked() == True :
		    #self.File_Name_2.setText(self.Smoothing.isChecked())
		  #  Smoothing = yes
	  #  else:
		    #self.File_Name_2.setText(self.Smoothing.isChecked())
		   # Smoothing = no
	#Reads Step Tolerance	    
   # def Step_Tol_F(self):
	#    Step_Tol = self.Step_Tol_Val.value()
	#Reads Amplitude Tolerance
   # def A_Tol_F(self):
	#   A_Tol = self.A_Tol_Val.value()
###############################################################
    #def FileNo_2_F(self):
	    
	  #  fileno = self.FileNo_2.value()

    def MAGGOT_SingleFit(self): #Fits a single NMR pulse file
    
	filepath=self.FilePath.text()
	filename=self.File_Name_2.text()
	fileno = self.FileNo_2.value()
	Freq = self.SetUp_Freq_Val.value()
	
	LoadAscii(Filename=filepath+filename+str(fileno),OutputWorkspace=filename+str(fileno),Unit='Time')
	Fit(Function='name=UserFunction,Formula=A*exp(-(x/T2))*cos(2*'+str(pi)+'*(f*x+p)),A=0.0004,T2=0.1,f='+str(Freq)+'.0,p=0.000000',InputWorkspace=filename+str(fileno) ,Output=filename+str(fileno)+'res',StartX='0.002',EndX='0.19999')
	
	table = mtd[str(filename)+str(fileno)+'res_Parameters']
	self.SFAmplitude_2.setText(str(table.cell(0,1)))
	self.T2_2.setText(str(table.cell(1,1)))
	self.SFResFreq_2.setText(str(table.cell(2,1)))
	self.Phase_2.setText(str(table.cell(3,1)))
	

    def MAGGOT(self):	#Basic fitting algorithm for polarising a cell
##Open files in appropriate place to write fit parameters to	
	try:
		A_NMR= open(r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\A_NMR', 'w')
		NMRDIAG= open(r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\NMR Diagnostics.csv', 'w')
		ISIS='yes'
	except:
		A_NMR= open(r'C:\MantidInstall\logs\\A_NMR', 'w')
		NMRDIAG= open(r'C:\MantidInstall\logs\\NMR Diagnostics.csv', 'w')
		ISIS='no'
##read filepath,name,frequency and file range from GUI
	filepath=self.FilePath.text()
	filename=self.File_Name_2.text()
	Freq = self.SetUp_Freq_Val.value()
	#print 'Freq'+ str(Freq)
	n0=self.InitialNMR.value()
	n=self.FinalNMR.value()
##load and fit range of files
	try:
		for i in range(n0,n+1):
			LoadAscii(Filename=str(filepath)+str(filename)+str(i),OutputWorkspace=str(filename),Unit='Time')
##fit first file with input params
			if i == self.InitialNMR.value():
				Fit(Function='name=UserFunction,Formula=A*exp(-(x/T2))*cos(2*'+str(pi)+'*(f*x+p)),A=0.0004,T2=0.1,f='+str(Freq)+'.0,p=0.000000',InputWorkspace=str(filename)  ,Output=str(filename)+'res',StartX='0.002',EndX='0.19999')
##fit subsequent with params from previous fit				
			else:
				Fit(Function='name=UserFunction,Formula=A*exp(-(x/T2))*cos(2*'+str(pi)+'*(f*x+p)),A=0.0004,T2=0.1,f='+ str(abs(table1.cell(2,1)))+',p='+ str(abs(table1.cell(3,1))),InputWorkspace=str(filename)  ,Output=str(filename)+'res',StartX='0.002',EndX='0.19999')
				
			table1 = mtd[str(filename) +'res_Parameters']
##write params to file(daig mode writes all and normal mode only writes Amplitude and amplitude error)			
			if self.DiagMode.isChecked() == True : 
				if i==n0:
					NMRDIAG.write('T'+','+'A'+','+'Aer' + ',' + 'T2'+',' + 'T2er' + ',' + 'f' +','+ 'fer'+',' + 'p' + ','+'per'+ '\n')
				NMRDIAG.write(str(i*self.TimeBetweenNMR.value()) +','+ str(abs(table1.cell(0,1)))+','+ str(table1.cell(0,2))+','+ str(table1.cell(1,1))+','+ str(table1.cell(1,2))+','+ str(table1.cell(2,1))+','+ str(table1.cell(2,2))+','+ str(table1.cell(3,1))+',' + str(table1.cell(3,2))+'\n')
				
			else:
				A_NMR.write(str(i*self.TimeBetweenNMR.value())+'	'+ str(abs(table1.cell(0,1)))+'	' + str(table1.cell(0,2))+'\n')
							
			self.MF_Results.setTitle(QtGui.QApplication.translate("MAGGOTWindow", "Fit complete for file "+str(i)+".", None, QtGui.QApplication.UnicodeUTF8))
##error if no files fit or all files fitted
	except:
		if i == n0:
			print 'Error: Initial filepath incorrect.'
		elif n0 < i <= n+1 :
			print 'All NMR files fitted.'
##close params file
	NMRDIAG.close()
	A_NMR.close()		
	self.MF_Results.setTitle(QtGui.QApplication.translate("MAGGOTWindow", "Results", None, QtGui.QApplication.UnicodeUTF8))
##choose function for fitting polarisation rate behaviour
	if self.UP_Down_LPP.currentText() == 'Spin Up':
		Function='name=UserFunction,Formula=A*(1-exp( -(x+T0)/T1)),T0=0.011,A=0.000286216,T1=10.00,constraints=(0<=T1)'
	
	elif self.UP_Down_LPP.currentText() == 'T1 Decay' or 'Loss Per Pulse':
		Function='name=UserFunction,Formula=A*(exp( -(x+T0)/T1)),T0=0.011,A=0.000106216,T1=50.00,constraints=(0<=T1)'
	
	#elif self.UP_Down_LPP.currentText() == 'Loss Per Pulse':
		#Function='name=UserFunction,Formula=A*(exp( -(x+T0)/T1)),T0=0.011,A=0.000286216,T1=50.00,constraints=(0<=T1)'
	
	else:
		fun = 0.000017239
	#if Smoothing == yes:
	#	LoadAscii(Filename=r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\ANMRsmooth',OutputWorkspace='T1_Data')
	#	Fit(Function,InputWorkspace='T1_Data',Output=File+'res',StartX='0',EndX=str(i*T))
	#	table = mtd[File +'res_Parameters']
		
	#elif Smoothing == no :
##load parameter files from fitting
	if self.DiagMode.isChecked() == True :
		if ISIS == 'no':
			LoadAscii(Filename=r'C:\MantidInstall\logs\\NMR Diagnostics.csv',OutputWorkspace='Rate_Data',Unit='Time')
		else:
			LoadAscii(Filename=r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\NMR Diagnostics.csv',OutputWorkspace='Rate_Data',Unit='Time')
	else:
		if ISIS == 'no':
			LoadAscii(Filename=r'C:\MantidInstall\logs\\A_NMR',OutputWorkspace='Rate_Data',Unit='Time')
		else:
			LoadAscii(Filename=r'\\Britannic\3he\NMR\1 Current NMR Data\1Extracted Fit Data\\A_NMR',OutputWorkspace='Rate_Data',Unit='Time')
##fit using already chosen function
	Fit(Function,InputWorkspace='Rate_Data',Output=str(filename)+'res',StartX='0',EndX=str(i*self.TimeBetweenNMR.value()))
	table = mtd[str(filename) +'res_Parameters']
##write results of relevant fit into GUI	
	if self.UP_Down_LPP.currentText() == 'Spin Up':
		self.MF_Amplitude_2.setText(str(table.cell(0,1)))
		self.T1orTup.setText(QtGui.QApplication.translate("MAGGOTWindow", "Tup(Time)", None, QtGui.QApplication.UnicodeUTF8))
		self.T1orTup_2.setText(str(table.cell(2,1)))
		
	elif self.UP_Down_LPP.currentText() == 'T1 Decay':
		self.MF_Amplitude_2.setText(str(table.cell(0,1)))
		self.T1orTup.setText(QtGui.QApplication.translate("MAGGOTWindow", "T1(Time)", None, QtGui.QApplication.UnicodeUTF8))
		self.T1orTup_2.setText(str(table.cell(2,1)))
		
	elif self.UP_Down_LPP.currentText() == 'Loss Per Pulse':
		self.MF_Amplitude_2.setText(str(table.cell(0,1)))
		self.LPP_2.setText(str(exp(-1/table.cell(2,1)) ) )
		self.LPP_Value.setValue(exp(-1/table.cell(2,1)) )
	
	
    def LPF_f(self):

	p=0 #number of pulses with amplitude between preaverage and postaverage
	PreSum=0
	PostSum=0
	PreAverage=0
	PostAverage=0
	Frequency= self.SetUp_Freq_Val.value()
	Filepath = (self.FilePath_2.text())
	PreFile = (self.PreFlipName_2.text())
	PostFile = (self.PostFlipName_2.text())
	PathandName= Filepath + PreFile
	#print PathandName
	
	try:	
		for i in range (1, (1000)):
			LoadAscii(Filename=Filepath+PreFile+str(i),OutputWorkspace=PreFile)
			Fit(Function='name=UserFunction,Formula=A*exp(-(x/T2))*cos(2*'+str(pi)+'*(f*x+p)),A=0.0004,T2=0.1,f='+str(Frequency)+'.0,p=0.000000',InputWorkspace=PreFile  ,Output=PreFile+'res',StartX='0.002',EndX='0.19999',CalcErrors=True)
			table = mtd[str(PreFile)+'res_Parameters']	
			print 'Amplitude is ' + str(abs(table.cell(0,1))) + ' mV.'	
			PreSum=PreSum+abs(table.cell(0,1))
		
	except: 
		if i==1:
			print 'Error: only 1 PreFile attempted'
		else: 
			#print i
			print PreSum
			PreAverage = float(PreSum)/(float(i-1))
			print PreAverage

	try:
		for i in range (1, (1000)):	
			LoadAscii(Filename=Filepath+PostFile+str(i),OutputWorkspace=PostFile)
			Fit(Function='name=UserFunction,Formula=A*exp(-(x/T2))*cos(2*'+str(pi)+'*(f*x+p)),A=0.0004,T2=0.1,f='+str(Frequency)+'.0,p=0.000000',InputWorkspace=PostFile  ,Output=PostFile+'res',StartX='0.002',EndX='0.19999',CalcErrors=True)
			table = mtd[str(PostFile)+'res_Parameters']
			print 'Amplitude is ' + str(abs(table.cell(0,1))) + ' mV.'
			PostSum=PostSum+abs(table.cell(0,1))
	except:
		if i==1:
			print 'Error: only 1 PostFile attempted'
		else:
			#print i
			print PostSum
			PostAverage = float(PostSum)/(float(i-1))
			print PostAverage
	
	if self.Use_LPP.isChecked() == True : 
		try: 
			for i in range (1, 1000):
				LoadAscii(Filename=Filepath+PreFile+str(i),OutputWorkspace=PreFile)
				Fit(Function='name=UserFunction,Formula=A*exp(-(x/T2))*cos(2*'+str(pi)+'*(f*x+p)),A=0.0004,T2=0.1,f='+str(Frequency)+'.0,p=0.000000',InputWorkspace=PreFile  ,Output=PreFile+'res',StartX='0.002',EndX='0.19999',CalcErrors=True)
				table = mtd[str(PreFile)+'res_Parameters']
				if (abs(table.cell(0,1))) > PostAverage and (abs(table.cell(0,1))) < PreAverage: 
					p=p+1
				else:
					p=p
		except: 
			if i==1:
				print 'Error: only 1 PreFile attempted'
		print 'p from first is '+str(p)
		try:
			for i in range (1, 1000):
				LoadAscii(Filename=Filepath+PostFile+str(i),OutputWorkspace=PostFile)
				Fit(Function='name=UserFunction,Formula=A*exp(-(x/T2))*cos(2*'+str(pi)+'*(f*x+p)),A=0.0004,T2=0.1,f='+str(Frequency)+'.0,p=0.000000',InputWorkspace=PostFile  ,Output=PostFile+'res',StartX='0.002',EndX='0.19999',CalcErrors=True)
				table = mtd[str(PostFile)+'res_Parameters']
				if (abs(table.cell(0,1))) > PostAverage and (abs(table.cell(0,1))) < PreAverage :
					p=p+1
				else:
					p=p
		except: 
			if i==1:
				print 'Error: only 1 PostFile attempted'
		print 'p from second is '+str(p)
	else: 
		p=0
		
	
	F=((1/(float(1-self.LPP_Value.value())**float(p)))*((float(PostSum))/(float(PreSum))))**(1/float(self.LPF_FinalNMR.value()))
	print 'Loss per flip is ' +str(1-F)
	self.LPF_Result.setText(str(1-F))