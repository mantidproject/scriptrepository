# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mslice\app\mainwindow.ui'
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
        MainWindow.resize(740, 689)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.tabWidget3 = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget3.setObjectName(_fromUtf8("tabWidget3"))
        self.tab_5 = QtGui.QWidget()
        self.tab_5.setObjectName(_fromUtf8("tab_5"))
        self.gridLayout_6 = QtGui.QGridLayout(self.tab_5)
        self.gridLayout_6.setObjectName(_fromUtf8("gridLayout_6"))
        self.wgtSlice = SliceWidget(self.tab_5)
        self.wgtSlice.setObjectName(_fromUtf8("wgtSlice"))
        self.gridLayout_6.addWidget(self.wgtSlice, 0, 0, 1, 1)
        self.tabWidget3.addTab(self.tab_5, _fromUtf8(""))
        self.tab_6 = QtGui.QWidget()
        self.tab_6.setObjectName(_fromUtf8("tab_6"))
        self.gridLayout_7 = QtGui.QGridLayout(self.tab_6)
        self.gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
        self.wgtCut = CutWidget(self.tab_6)
        self.wgtCut.setObjectName(_fromUtf8("wgtCut"))
        self.gridLayout_7.addWidget(self.wgtCut, 0, 0, 1, 1)
        self.tabWidget3.addTab(self.tab_6, _fromUtf8(""))
        self.gridLayout.addWidget(self.tabWidget3, 2, 0, 1, 1)
        self.tabWidget2 = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget2.setObjectName(_fromUtf8("tabWidget2"))
        self.tab_3 = QtGui.QWidget()
        self.tab_3.setObjectName(_fromUtf8("tab_3"))
        self.gridLayout_4 = QtGui.QGridLayout(self.tab_3)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        spacerItem = QtGui.QSpacerItem(666, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_4.addItem(spacerItem, 0, 1, 1, 1)
        self.wgtPowder = PowderWidget(self.tab_3)
        self.wgtPowder.setObjectName(_fromUtf8("wgtPowder"))
        self.gridLayout_4.addWidget(self.wgtPowder, 0, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 137, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_4.addItem(spacerItem1, 1, 0, 1, 1)
        self.tabWidget2.addTab(self.tab_3, _fromUtf8(""))
        self.gridLayout.addWidget(self.tabWidget2, 1, 0, 1, 1)
        self.tabWidget1 = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget1.setObjectName(_fromUtf8("tabWidget1"))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.gridLayout_2 = QtGui.QGridLayout(self.tab)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.wgtWorkspacemanager = WorkspaceManagerWidget(self.tab)
        self.wgtWorkspacemanager.setObjectName(_fromUtf8("wgtWorkspacemanager"))
        self.gridLayout_2.addWidget(self.wgtWorkspacemanager, 0, 0, 1, 1)
        self.tabWidget1.addTab(self.tab, _fromUtf8(""))
        self.gridLayout.addWidget(self.tabWidget1, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget3.setCurrentIndex(0)
        self.tabWidget2.setCurrentIndex(0)
        self.tabWidget1.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "Mantid Mslice", None))
        self.tabWidget3.setTabText(self.tabWidget3.indexOf(self.tab_5), _translate("MainWindow", "Slice", None))
        self.tabWidget3.setTabText(self.tabWidget3.indexOf(self.tab_6), _translate("MainWindow", "Cut", None))
        self.tabWidget2.setTabText(self.tabWidget2.indexOf(self.tab_3), _translate("MainWindow", "Powder", None))
        self.tabWidget1.setTabText(self.tabWidget1.indexOf(self.tab), _translate("MainWindow", "Workspace Manager", None))

from mslice.widgets.cut.cut import CutWidget
from mslice.widgets.projection.powder.powder import PowderWidget
from mslice.widgets.slice.slice import SliceWidget
from mslice.widgets.workspacemanager.workspacemanager import WorkspaceManagerWidget
