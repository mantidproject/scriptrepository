# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mslice\plotting\plot_window\plot_window.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
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
        MainWindow.resize(800, 591)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuInformation = QtGui.QMenu(self.menubar)
        self.menuInformation.setObjectName(_fromUtf8("menuInformation"))
        self.menuIntensity = QtGui.QMenu(self.menubar)
        self.menuIntensity.setObjectName(_fromUtf8("menuIntensity"))
        self.menuMake_Script = QtGui.QMenu(self.menubar)
        self.menuMake_Script.setObjectName(_fromUtf8("menuMake_Script"))
        self.menuHelp = QtGui.QMenu(self.menubar)
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.toolBar = QtGui.QToolBar(MainWindow)
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.action_save_image = QtGui.QAction(MainWindow)
        self.action_save_image.setObjectName(_fromUtf8("action_save_image"))
        self.actionSave_To_File = QtGui.QAction(MainWindow)
        self.actionSave_To_File.setObjectName(_fromUtf8("actionSave_To_File"))
        self.actionDump_To_Console = QtGui.QAction(MainWindow)
        self.actionDump_To_Console.setObjectName(_fromUtf8("actionDump_To_Console"))
        self.actionZoom_In = QtGui.QAction(MainWindow)
        self.actionZoom_In.setObjectName(_fromUtf8("actionZoom_In"))
        self.actionZoom_Out = QtGui.QAction(MainWindow)
        self.actionZoom_Out.setObjectName(_fromUtf8("actionZoom_Out"))
        self.actionDataCursor = QtGui.QAction(MainWindow)
        self.actionDataCursor.setCheckable(True)
        self.actionDataCursor.setObjectName(_fromUtf8("actionDataCursor"))
        self.actionKeep = QtGui.QAction(MainWindow)
        self.actionKeep.setCheckable(True)
        self.actionKeep.setObjectName(_fromUtf8("actionKeep"))
        self.actionMakeCurrent = QtGui.QAction(MainWindow)
        self.actionMakeCurrent.setCheckable(True)
        self.actionMakeCurrent.setChecked(True)
        self.actionMakeCurrent.setObjectName(_fromUtf8("actionMakeCurrent"))
        self.actionPlotOptions = QtGui.QAction(MainWindow)
        self.actionPlotOptions.setObjectName(_fromUtf8("actionPlotOptions"))
        self.actionToggleLegends = QtGui.QAction(MainWindow)
        self.actionToggleLegends.setObjectName(_fromUtf8("actionToggleLegends"))
        self.menuMake_Script.addAction(self.actionSave_To_File)
        self.menuMake_Script.addAction(self.actionDump_To_Console)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuInformation.menuAction())
        self.menubar.addAction(self.menuIntensity.menuAction())
        self.menubar.addAction(self.menuMake_Script.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.toolBar.addAction(self.action_save_image)
        self.toolBar.addAction(self.actionZoom_In)
        self.toolBar.addAction(self.actionZoom_Out)
        self.toolBar.addAction(self.actionDataCursor)
        self.toolBar.addAction(self.actionPlotOptions)
        self.toolBar.addAction(self.actionToggleLegends)
        self.toolBar.addAction(self.actionKeep)
        self.toolBar.addAction(self.actionMakeCurrent)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.menuInformation.setTitle(_translate("MainWindow", "Information", None))
        self.menuIntensity.setTitle(_translate("MainWindow", "Intensity", None))
        self.menuMake_Script.setTitle(_translate("MainWindow", "Make Script", None))
        self.menuHelp.setTitle(_translate("MainWindow", "Help", None))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar", None))
        self.action_save_image.setText(_translate("MainWindow", "Save Plot", None))
        self.action_save_image.setToolTip(_translate("MainWindow", "Save Plot As Image", None))
        self.actionSave_To_File.setText(_translate("MainWindow", "Save To File", None))
        self.actionDump_To_Console.setText(_translate("MainWindow", "Dump To Console", None))
        self.actionDump_To_Console.setShortcut(_translate("MainWindow", "Ctrl+D", None))
        self.actionZoom_In.setText(_translate("MainWindow", "Zoom In", None))
        self.actionZoom_Out.setText(_translate("MainWindow", "Zoom Out", None))
        self.actionDataCursor.setText(_translate("MainWindow", "Data Cursor", None))
        self.actionKeep.setText(_translate("MainWindow", "Keep", None))
        self.actionKeep.setToolTip(_translate("MainWindow", "Keep", None))
        self.actionMakeCurrent.setText(_translate("MainWindow", "Make Current", None))
        self.actionPlotOptions.setText(_translate("MainWindow", "Plot Options", None))
        self.actionToggleLegends.setText(_translate("MainWindow", "Toggle Legends", None))

