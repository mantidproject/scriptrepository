# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'PysliceUi_2.ui'
#
# Created: Thu Nov 21 12:01:18 2013
#      by: PyQt4 UI code generator 4.10
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(831, 512)
        MainWindow.setToolTip(_fromUtf8(""))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.DetVanNum = QtGui.QLineEdit(self.tab)
        self.DetVanNum.setGeometry(QtCore.QRect(140, 10, 113, 22))
        self.DetVanNum.setObjectName(_fromUtf8("DetVanNum"))
        self.RunNum = QtGui.QLineEdit(self.tab)
        self.RunNum.setGeometry(QtCore.QRect(140, 50, 113, 22))
        self.RunNum.setToolTip(_fromUtf8(""))
        self.RunNum.setObjectName(_fromUtf8("RunNum"))
        self.label_16 = QtGui.QLabel(self.tab)
        self.label_16.setGeometry(QtCore.QRect(10, 10, 131, 16))
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.label_17 = QtGui.QLabel(self.tab)
        self.label_17.setGeometry(QtCore.QRect(10, 50, 131, 16))
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.Reduce = QtGui.QPushButton(self.tab)
        self.Reduce.setGeometry(QtCore.QRect(0, 130, 114, 71))
        self.Reduce.setObjectName(_fromUtf8("Reduce"))
        self.ReduceOutput = QtGui.QListWidget(self.tab)
        self.ReduceOutput.setGeometry(QtCore.QRect(10, 210, 781, 221))
        self.ReduceOutput.setObjectName(_fromUtf8("ReduceOutput"))
        self.sumRuns = QtGui.QCheckBox(self.tab)
        self.sumRuns.setGeometry(QtCore.QRect(270, 50, 141, 20))
        self.sumRuns.setObjectName(_fromUtf8("sumRuns"))
        self.label_18 = QtGui.QLabel(self.tab)
        self.label_18.setGeometry(QtCore.QRect(10, 80, 131, 41))
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.outputWksp = QtGui.QLineEdit(self.tab)
        self.outputWksp.setGeometry(QtCore.QRect(140, 90, 113, 22))
        self.outputWksp.setToolTip(_fromUtf8(""))
        self.outputWksp.setText(_fromUtf8(""))
        self.outputWksp.setObjectName(_fromUtf8("outputWksp"))
        self.tabWidget_2 = QtGui.QTabWidget(self.tab)
        self.tabWidget_2.setGeometry(QtCore.QRect(410, 0, 381, 201))
        self.tabWidget_2.setObjectName(_fromUtf8("tabWidget_2"))
        self.tab_3 = QtGui.QWidget()
        self.tab_3.setObjectName(_fromUtf8("tab_3"))
        self.BkgSwitch = QtGui.QCheckBox(self.tab_3)
        self.BkgSwitch.setGeometry(QtCore.QRect(10, 30, 181, 20))
        self.BkgSwitch.setObjectName(_fromUtf8("BkgSwitch"))
        self.AutoEiChbox = QtGui.QCheckBox(self.tab_3)
        self.AutoEiChbox.setGeometry(QtCore.QRect(10, 50, 71, 20))
        self.AutoEiChbox.setObjectName(_fromUtf8("AutoEiChbox"))
        self.FixEi = QtGui.QCheckBox(self.tab_3)
        self.FixEi.setGeometry(QtCore.QRect(10, 70, 61, 20))
        self.FixEi.setObjectName(_fromUtf8("FixEi"))
        self.EiGuess = QtGui.QLineEdit(self.tab_3)
        self.EiGuess.setGeometry(QtCore.QRect(80, 70, 41, 21))
        self.EiGuess.setText(_fromUtf8(""))
        self.EiGuess.setObjectName(_fromUtf8("EiGuess"))
        self.MonSpecNumber = QtGui.QLineEdit(self.tab_3)
        self.MonSpecNumber.setGeometry(QtCore.QRect(80, 95, 41, 20))
        self.MonSpecNumber.setText(_fromUtf8(""))
        self.MonSpecNumber.setObjectName(_fromUtf8("MonSpecNumber"))
        self.FixMonitorSpectrum = QtGui.QCheckBox(self.tab_3)
        self.FixMonitorSpectrum.setGeometry(QtCore.QRect(10, 90, 71, 31))
        self.FixMonitorSpectrum.setObjectName(_fromUtf8("FixMonitorSpectrum"))
        self.NormMethod = QtGui.QComboBox(self.tab_3)
        self.NormMethod.setGeometry(QtCore.QRect(10, 120, 81, 26))
        self.NormMethod.setObjectName(_fromUtf8("NormMethod"))
        self.NormMethod.addItem(_fromUtf8(""))
        self.NormMethod.addItem(_fromUtf8(""))
        self.label_23 = QtGui.QLabel(self.tab_3)
        self.label_23.setGeometry(QtCore.QRect(90, 125, 91, 16))
        self.label_23.setObjectName(_fromUtf8("label_23"))
        self.mapfile = QtGui.QLineEdit(self.tab_3)
        self.mapfile.setGeometry(QtCore.QRect(10, 150, 131, 21))
        self.mapfile.setText(_fromUtf8(""))
        self.mapfile.setObjectName(_fromUtf8("mapfile"))
        self.label_24 = QtGui.QLabel(self.tab_3)
        self.label_24.setGeometry(QtCore.QRect(150, 150, 60, 16))
        self.label_24.setObjectName(_fromUtf8("label_24"))
        self.InstName = QtGui.QComboBox(self.tab_3)
        self.InstName.setGeometry(QtCore.QRect(10, 0, 111, 26))
        self.InstName.setObjectName(_fromUtf8("InstName"))
        self.InstName.addItem(_fromUtf8(""))
        self.InstName.addItem(_fromUtf8(""))
        self.InstName.addItem(_fromUtf8(""))
        self.InstName.addItem(_fromUtf8(""))
        self.InstName.addItem(_fromUtf8(""))
        self.label_25 = QtGui.QLabel(self.tab_3)
        self.label_25.setGeometry(QtCore.QRect(130, 73, 60, 16))
        self.label_25.setObjectName(_fromUtf8("label_25"))
        self.tabWidget_2.addTab(self.tab_3, _fromUtf8(""))
        self.tab_4 = QtGui.QWidget()
        self.tab_4.setObjectName(_fromUtf8("tab_4"))
        self.AbsNormSwitch = QtGui.QCheckBox(self.tab_4)
        self.AbsNormSwitch.setGeometry(QtCore.QRect(10, 10, 221, 20))
        self.AbsNormSwitch.setObjectName(_fromUtf8("AbsNormSwitch"))
        self.DetVanNum_2 = QtGui.QLineEdit(self.tab_4)
        self.DetVanNum_2.setGeometry(QtCore.QRect(140, 50, 113, 22))
        self.DetVanNum_2.setObjectName(_fromUtf8("DetVanNum_2"))
        self.label_19 = QtGui.QLabel(self.tab_4)
        self.label_19.setGeometry(QtCore.QRect(10, 80, 131, 16))
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.label_20 = QtGui.QLabel(self.tab_4)
        self.label_20.setGeometry(QtCore.QRect(10, 50, 131, 31))
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.MonoRunNum = QtGui.QLineEdit(self.tab_4)
        self.MonoRunNum.setGeometry(QtCore.QRect(140, 80, 113, 22))
        self.MonoRunNum.setToolTip(_fromUtf8(""))
        self.MonoRunNum.setObjectName(_fromUtf8("MonoRunNum"))
        self.sampleMass = QtGui.QLineEdit(self.tab_4)
        self.sampleMass.setGeometry(QtCore.QRect(140, 110, 113, 22))
        self.sampleMass.setText(_fromUtf8(""))
        self.sampleMass.setObjectName(_fromUtf8("sampleMass"))
        self.label_21 = QtGui.QLabel(self.tab_4)
        self.label_21.setGeometry(QtCore.QRect(10, 140, 131, 31))
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.label_22 = QtGui.QLabel(self.tab_4)
        self.label_22.setGeometry(QtCore.QRect(10, 110, 131, 16))
        self.label_22.setObjectName(_fromUtf8("label_22"))
        self.RMMmass = QtGui.QLineEdit(self.tab_4)
        self.RMMmass.setGeometry(QtCore.QRect(140, 140, 113, 22))
        self.RMMmass.setToolTip(_fromUtf8(""))
        self.RMMmass.setText(_fromUtf8(""))
        self.RMMmass.setObjectName(_fromUtf8("RMMmass"))
        self.tabWidget_2.addTab(self.tab_4, _fromUtf8(""))
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.label_6 = QtGui.QLabel(self.tab_2)
        self.label_6.setGeometry(QtCore.QRect(534, 61, 62, 16))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.label = QtGui.QLabel(self.tab_2)
        self.label.setGeometry(QtCore.QRect(410, 92, 31, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.plot = QtGui.QPushButton(self.tab_2)
        self.plot.setGeometry(QtCore.QRect(554, 142, 91, 32))
        self.plot.setObjectName(_fromUtf8("plot"))
        self.removePlot = QtGui.QPushButton(self.tab_2)
        self.removePlot.setGeometry(QtCore.QRect(294, 380, 150, 32))
        self.removePlot.setObjectName(_fromUtf8("removePlot"))
        self.label_7 = QtGui.QLabel(self.tab_2)
        self.label_7.setGeometry(QtCore.QRect(409, 163, 31, 16))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.scalemin = QtGui.QLineEdit(self.tab_2)
        self.scalemin.setGeometry(QtCore.QRect(383, 178, 81, 22))
        self.scalemin.setObjectName(_fromUtf8("scalemin"))
        self.label_14 = QtGui.QLabel(self.tab_2)
        self.label_14.setGeometry(QtCore.QRect(299, 220, 101, 16))
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.cutMax = QtGui.QLineEdit(self.tab_2)
        self.cutMax.setGeometry(QtCore.QRect(469, 112, 81, 22))
        self.cutMax.setObjectName(_fromUtf8("cutMax"))
        self.label_10 = QtGui.QLabel(self.tab_2)
        self.label_10.setGeometry(QtCore.QRect(669, 330, 81, 21))
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.frame = QtGui.QFrame(self.tab_2)
        self.frame.setGeometry(QtCore.QRect(319, 3, 301, 80))
        self.frame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtGui.QFrame.Raised)
        self.frame.setObjectName(_fromUtf8("frame"))
        self.label_3 = QtGui.QLabel(self.tab_2)
        self.label_3.setGeometry(QtCore.QRect(584, 92, 30, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.removeWksp = QtGui.QPushButton(self.tab_2)
        self.removeWksp.setGeometry(QtCore.QRect(3, 410, 150, 32))
        self.removeWksp.setObjectName(_fromUtf8("removeWksp"))
        self.label_4 = QtGui.QLabel(self.tab_2)
        self.label_4.setGeometry(QtCore.QRect(322, 142, 71, 20))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.intMin = QtGui.QLineEdit(self.tab_2)
        self.intMin.setGeometry(QtCore.QRect(383, 142, 81, 22))
        self.intMin.setObjectName(_fromUtf8("intMin"))
        self.DispMin = QtGui.QLineEdit(self.tab_2)
        self.DispMin.setGeometry(QtCore.QRect(362, 40, 71, 20))
        self.DispMin.setObjectName(_fromUtf8("DispMin"))
        self.scalemax = QtGui.QLineEdit(self.tab_2)
        self.scalemax.setGeometry(QtCore.QRect(469, 177, 81, 22))
        self.scalemax.setObjectName(_fromUtf8("scalemax"))
        self.label_2 = QtGui.QLabel(self.tab_2)
        self.label_2.setGeometry(QtCore.QRect(490, 92, 31, 16))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_13 = QtGui.QLabel(self.tab_2)
        self.label_13.setGeometry(QtCore.QRect(456, 60, 41, 20))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.label_9 = QtGui.QLabel(self.tab_2)
        self.label_9.setGeometry(QtCore.QRect(674, 284, 62, 16))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.label_12 = QtGui.QLabel(self.tab_2)
        self.label_12.setGeometry(QtCore.QRect(9, 90, 101, 16))
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.FigureList = QtGui.QListWidget(self.tab_2)
        self.FigureList.setGeometry(QtCore.QRect(299, 240, 361, 141))
        self.FigureList.setObjectName(_fromUtf8("FigureList"))
        self.Temperature = QtGui.QLineEdit(self.tab_2)
        self.Temperature.setGeometry(QtCore.QRect(668, 260, 71, 22))
        self.Temperature.setObjectName(_fromUtf8("Temperature"))
        self.intMax = QtGui.QLineEdit(self.tab_2)
        self.intMax.setGeometry(QtCore.QRect(469, 142, 81, 22))
        self.intMax.setObjectName(_fromUtf8("intMax"))
        self.WkspIn = QtGui.QComboBox(self.tab_2)
        self.WkspIn.setGeometry(QtCore.QRect(13, 4, 171, 26))
        self.WkspIn.setStatusTip(_fromUtf8(""))
        self.WkspIn.setObjectName(_fromUtf8("WkspIn"))
        self.DispMax = QtGui.QLineEdit(self.tab_2)
        self.DispMax.setGeometry(QtCore.QRect(441, 40, 71, 20))
        self.DispMax.setObjectName(_fromUtf8("DispMax"))
        self.refresh = QtGui.QPushButton(self.tab_2)
        self.refresh.setGeometry(QtCore.QRect(9, 30, 114, 32))
        self.refresh.setObjectName(_fromUtf8("refresh"))
        self.axis = QtGui.QComboBox(self.tab_2)
        self.axis.setGeometry(QtCore.QRect(303, 109, 71, 26))
        self.axis.setObjectName(_fromUtf8("axis"))
        self.label_11 = QtGui.QLabel(self.tab_2)
        self.label_11.setGeometry(QtCore.QRect(664, 350, 91, 21))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.label_15 = QtGui.QLabel(self.tab_2)
        self.label_15.setGeometry(QtCore.QRect(379, 60, 41, 20))
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.genScript = QtGui.QPushButton(self.tab_2)
        self.genScript.setGeometry(QtCore.QRect(662, 370, 91, 40))
        self.genScript.setObjectName(_fromUtf8("genScript"))
        self.calcProj = QtGui.QPushButton(self.tab_2)
        self.calcProj.setGeometry(QtCore.QRect(200, 0, 101, 61))
        self.calcProj.setObjectName(_fromUtf8("calcProj"))
        self.plotOver = QtGui.QPushButton(self.tab_2)
        self.plotOver.setGeometry(QtCore.QRect(554, 172, 91, 32))
        self.plotOver.setObjectName(_fromUtf8("plotOver"))
        self.delta = QtGui.QLineEdit(self.tab_2)
        self.delta.setGeometry(QtCore.QRect(559, 112, 81, 22))
        self.delta.setObjectName(_fromUtf8("delta"))
        self.display = QtGui.QPushButton(self.tab_2)
        self.display.setGeometry(QtCore.QRect(357, 10, 241, 21))
        self.display.setObjectName(_fromUtf8("display"))
        self.label_5 = QtGui.QLabel(self.tab_2)
        self.label_5.setGeometry(QtCore.QRect(305, 177, 81, 20))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.cutMin = QtGui.QLineEdit(self.tab_2)
        self.cutMin.setGeometry(QtCore.QRect(383, 112, 81, 22))
        self.cutMin.setObjectName(_fromUtf8("cutMin"))
        self.wkspList = QtGui.QListWidget(self.tab_2)
        self.wkspList.setGeometry(QtCore.QRect(9, 110, 271, 301))
        self.wkspList.setObjectName(_fromUtf8("wkspList"))
        self.smooth = QtGui.QLineEdit(self.tab_2)
        self.smooth.setGeometry(QtCore.QRect(522, 39, 71, 22))
        self.smooth.setObjectName(_fromUtf8("smooth"))
        self.label_8 = QtGui.QLabel(self.tab_2)
        self.label_8.setGeometry(QtCore.QRect(499, 162, 31, 16))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.DBWFac = QtGui.QLineEdit(self.tab_2)
        self.DBWFac.setGeometry(QtCore.QRect(670, 310, 71, 22))
        self.DBWFac.setObjectName(_fromUtf8("DBWFac"))
        self.DeltaQ = QtGui.QLineEdit(self.tab_2)
        self.DeltaQ.setGeometry(QtCore.QRect(210, 80, 81, 22))
        self.DeltaQ.setObjectName(_fromUtf8("DeltaQ"))
        self.label_26 = QtGui.QLabel(self.tab_2)
        self.label_26.setGeometry(QtCore.QRect(200, 60, 111, 20))
        self.label_26.setObjectName(_fromUtf8("label_26"))
        self.tabWidget.addTab(self.tab_2, _fromUtf8(""))
        self.horizontalLayout.addWidget(self.tabWidget)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menuBar = QtGui.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 831, 22))
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuCalculations = QtGui.QMenu(self.menuBar)
        self.menuCalculations.setObjectName(_fromUtf8("menuCalculations"))
        MainWindow.setMenuBar(self.menuBar)
        self.actionDensityOfstates = QtGui.QAction(MainWindow)
        self.actionDensityOfstates.setObjectName(_fromUtf8("actionDensityOfstates"))
        self.actionBoseFactor = QtGui.QAction(MainWindow)
        self.actionBoseFactor.setObjectName(_fromUtf8("actionBoseFactor"))
        self.menuCalculations.addAction(self.actionDensityOfstates)
        self.menuCalculations.addAction(self.actionBoseFactor)
        self.menuBar.addAction(self.menuCalculations.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(1)
        self.tabWidget_2.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        MainWindow.setTabOrder(self.cutMin, self.cutMax)
        MainWindow.setTabOrder(self.cutMax, self.delta)
        MainWindow.setTabOrder(self.delta, self.intMin)
        MainWindow.setTabOrder(self.intMin, self.intMax)
        MainWindow.setTabOrder(self.intMax, self.scalemin)
        MainWindow.setTabOrder(self.scalemin, self.scalemax)
        MainWindow.setTabOrder(self.scalemax, self.plot)
        MainWindow.setTabOrder(self.plot, self.DetVanNum)
        MainWindow.setTabOrder(self.DetVanNum, self.removeWksp)
        MainWindow.setTabOrder(self.removeWksp, self.Reduce)
        MainWindow.setTabOrder(self.Reduce, self.DispMin)
        MainWindow.setTabOrder(self.DispMin, self.removePlot)
        MainWindow.setTabOrder(self.removePlot, self.FigureList)
        MainWindow.setTabOrder(self.FigureList, self.Temperature)
        MainWindow.setTabOrder(self.Temperature, self.ReduceOutput)
        MainWindow.setTabOrder(self.ReduceOutput, self.WkspIn)
        MainWindow.setTabOrder(self.WkspIn, self.DispMax)
        MainWindow.setTabOrder(self.DispMax, self.refresh)
        MainWindow.setTabOrder(self.refresh, self.axis)
        MainWindow.setTabOrder(self.axis, self.genScript)
        MainWindow.setTabOrder(self.genScript, self.calcProj)
        MainWindow.setTabOrder(self.calcProj, self.plotOver)
        MainWindow.setTabOrder(self.plotOver, self.RunNum)
        MainWindow.setTabOrder(self.RunNum, self.display)
        MainWindow.setTabOrder(self.display, self.wkspList)
        MainWindow.setTabOrder(self.wkspList, self.smooth)
        MainWindow.setTabOrder(self.smooth, self.DBWFac)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.label_16.setText(_translate("MainWindow", "Detector Vanadium", None))
        self.label_17.setText(_translate("MainWindow", "Run numbers", None))
        self.Reduce.setText(_translate("MainWindow", "Reduce", None))
        self.sumRuns.setText(_translate("MainWindow", "sum these runs", None))
        self.label_18.setText(_translate("MainWindow", "Output workpspace\n"
"name", None))
        self.BkgSwitch.setText(_translate("MainWindow", "BackGround Subtraction", None))
        self.AutoEiChbox.setText(_translate("MainWindow", "AutoEi", None))
        self.FixEi.setText(_translate("MainWindow", "FixEi", None))
        self.FixMonitorSpectrum.setText(_translate("MainWindow", "Monitor\n"
"Spec", None))
        self.NormMethod.setItemText(0, _translate("MainWindow", "Current", None))
        self.NormMethod.setItemText(1, _translate("MainWindow", "Monitor-1", None))
        self.label_23.setText(_translate("MainWindow", "Normalisation", None))
        self.label_24.setText(_translate("MainWindow", "MapFile", None))
        self.InstName.setItemText(0, _translate("MainWindow", "Select Inst", None))
        self.InstName.setItemText(1, _translate("MainWindow", "Mari", None))
        self.InstName.setItemText(2, _translate("MainWindow", "Merlin", None))
        self.InstName.setItemText(3, _translate("MainWindow", "Maps", None))
        self.InstName.setItemText(4, _translate("MainWindow", "LET", None))
        self.label_25.setText(_translate("MainWindow", "EiGuess", None))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_3), _translate("MainWindow", "Reduction options", None))
        self.AbsNormSwitch.setText(_translate("MainWindow", "Perform absolute normalisation", None))
        self.label_19.setText(_translate("MainWindow", "Mono Van Run", None))
        self.label_20.setText(_translate("MainWindow", "Detector Vanadium\n"
"for mono run", None))
        self.label_21.setText(_translate("MainWindow", "Mass of a\n"
"formula unit", None))
        self.label_22.setText(_translate("MainWindow", "Sample Mass", None))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_4), _translate("MainWindow", "Absolute Normalisation", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Reduction", None))
        self.label_6.setText(_translate("MainWindow", "Smooth", None))
        self.label.setText(_translate("MainWindow", "Min", None))
        self.plot.setToolTip(_translate("MainWindow", "Plot cut", None))
        self.plot.setText(_translate("MainWindow", "Plot", None))
        self.removePlot.setText(_translate("MainWindow", "Remove figure", None))
        self.label_7.setText(_translate("MainWindow", "Min", None))
        self.scalemin.setToolTip(_translate("MainWindow", "Minimum valaue for cut integration", None))
        self.label_14.setText(_translate("MainWindow", "Figure List ", None))
        self.cutMax.setToolTip(_translate("MainWindow", "Maximum valaue for cut axis", None))
        self.label_10.setText(_translate("MainWindow", "MeanSquare", None))
        self.frame.setToolTip(_translate("MainWindow", "Controls for displaying Q, W map", None))
        self.label_3.setText(_translate("MainWindow", "Step", None))
        self.removeWksp.setText(_translate("MainWindow", "Remove workspace", None))
        self.label_4.setText(_translate("MainWindow", "Integrate ", None))
        self.intMin.setToolTip(_translate("MainWindow", "Minimum valaue for cut integration", None))
        self.scalemax.setToolTip(_translate("MainWindow", "Maximum valaue for cut integration", None))
        self.label_2.setText(_translate("MainWindow", "Max", None))
        self.label_13.setText(_translate("MainWindow", "I Max", None))
        self.label_9.setText(_translate("MainWindow", "Temp[K]", None))
        self.label_12.setText(_translate("MainWindow", "Workspace List ", None))
        self.FigureList.setToolTip(_translate("MainWindow", "List of Plots", None))
        self.Temperature.setToolTip(_translate("MainWindow", "Temperature of Run Units K", None))
        self.intMax.setToolTip(_translate("MainWindow", "Maximum valaue for cut integration", None))
        self.WkspIn.setToolTip(_translate("MainWindow", "Avalailable Workspaces", None))
        self.refresh.setToolTip(_translate("MainWindow", "Refresh available workspaces", None))
        self.refresh.setText(_translate("MainWindow", "Refresh List", None))
        self.axis.setToolTip(_translate("MainWindow", "Axis to cut along", None))
        self.label_11.setText(_translate("MainWindow", "Displacement", None))
        self.label_15.setText(_translate("MainWindow", "I Min", None))
        self.genScript.setText(_translate("MainWindow", "GenScript", None))
        self.calcProj.setText(_translate("MainWindow", "Calculate \n"
" Projections", None))
        self.plotOver.setToolTip(_translate("MainWindow", "Plot cut over existing plot (Must seleect a plot below)", None))
        self.plotOver.setText(_translate("MainWindow", "Plot Over", None))
        self.delta.setToolTip(_translate("MainWindow", "Step valaue for cut axis rebin", None))
        self.display.setText(_translate("MainWindow", "display", None))
        self.label_5.setText(_translate("MainWindow", "Y Axis scale", None))
        self.cutMin.setToolTip(_translate("MainWindow", "Minimum valaue for cut axis", None))
        self.wkspList.setToolTip(_translate("MainWindow", "List of transformed workspaces", None))
        self.label_8.setText(_translate("MainWindow", "Max", None))
        self.DBWFac.setToolTip(_translate("MainWindow", "Mean square displacement [Angstrom]", None))
        self.DeltaQ.setToolTip(_translate("MainWindow", "Minimum valaue for cut axis", None))
        self.label_26.setText(_translate("MainWindow", "DeltaQ (optional)", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Visulisation", None))
        self.menuCalculations.setTitle(_translate("MainWindow", "Calculations", None))
        self.actionDensityOfstates.setText(_translate("MainWindow", "DensityOfstates", None))
        self.actionBoseFactor.setText(_translate("MainWindow", "BoseFactor", None))

