# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mslice\widgets\status\status.ui'
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

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(661, 233)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_14 = QtGui.QGridLayout()
        self.gridLayout_14.setObjectName(_fromUtf8("gridLayout_14"))
        self.btnStatusLog = QtGui.QPushButton(Form)
        self.btnStatusLog.setObjectName(_fromUtf8("btnStatusLog"))
        self.gridLayout_14.addWidget(self.btnStatusLog, 0, 2, 1, 1)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_14.addItem(spacerItem, 0, 1, 1, 1)
        self.labelStatus_2 = QtGui.QLabel(Form)
        self.labelStatus_2.setObjectName(_fromUtf8("labelStatus_2"))
        self.gridLayout_14.addWidget(self.labelStatus_2, 0, 0, 1, 1)
        self.txtStatusLog = QtGui.QTextEdit(Form)
        self.txtStatusLog.setReadOnly(True)
        self.txtStatusLog.setObjectName(_fromUtf8("txtStatusLog"))
        self.gridLayout_14.addWidget(self.txtStatusLog, 1, 0, 1, 3)
        self.gridLayout.addLayout(self.gridLayout_14, 0, 0, 2, 1)
        self.groupBox_11 = QtGui.QGroupBox(Form)
        font = QtGui.QFont()
        font.setBold(False)
        font.setUnderline(False)
        font.setWeight(50)
        self.groupBox_11.setFont(font)
        self.groupBox_11.setObjectName(_fromUtf8("groupBox_11"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.groupBox_11)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.gridLayout_17 = QtGui.QGridLayout()
        self.gridLayout_17.setObjectName(_fromUtf8("gridLayout_17"))
        self.label_132 = QtGui.QLabel(self.groupBox_11)
        font = QtGui.QFont()
        font.setUnderline(False)
        self.label_132.setFont(font)
        self.label_132.setObjectName(_fromUtf8("label_132"))
        self.gridLayout_17.addWidget(self.label_132, 0, 0, 1, 2)
        self.prgStatusMemory = QtGui.QProgressBar(self.groupBox_11)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.prgStatusMemory.sizePolicy().hasHeightForWidth())
        self.prgStatusMemory.setSizePolicy(sizePolicy)
        self.prgStatusMemory.setProperty("value", 24)
        self.prgStatusMemory.setObjectName(_fromUtf8("prgStatusMemory"))
        self.gridLayout_17.addWidget(self.prgStatusMemory, 0, 2, 1, 1)
        self.labelMaxMem_2 = QtGui.QLabel(self.groupBox_11)
        self.labelMaxMem_2.setObjectName(_fromUtf8("labelMaxMem_2"))
        self.gridLayout_17.addWidget(self.labelMaxMem_2, 1, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(77, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_17.addItem(spacerItem1, 1, 1, 1, 2)
        self.verticalLayout_3.addLayout(self.gridLayout_17)
        self.gridLayout_16 = QtGui.QGridLayout()
        self.gridLayout_16.setObjectName(_fromUtf8("gridLayout_16"))
        self.label_164 = QtGui.QLabel(self.groupBox_11)
        font = QtGui.QFont()
        font.setUnderline(False)
        self.label_164.setFont(font)
        self.label_164.setObjectName(_fromUtf8("label_164"))
        self.gridLayout_16.addWidget(self.label_164, 0, 0, 1, 1)
        self.prgStatusCPUUsage = QtGui.QProgressBar(self.groupBox_11)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.prgStatusCPUUsage.sizePolicy().hasHeightForWidth())
        self.prgStatusCPUUsage.setSizePolicy(sizePolicy)
        self.prgStatusCPUUsage.setProperty("value", 24)
        self.prgStatusCPUUsage.setObjectName(_fromUtf8("prgStatusCPUUsage"))
        self.gridLayout_16.addWidget(self.prgStatusCPUUsage, 0, 1, 2, 1)
        self.labelCPUCount_2 = QtGui.QLabel(self.groupBox_11)
        self.labelCPUCount_2.setObjectName(_fromUtf8("labelCPUCount_2"))
        self.gridLayout_16.addWidget(self.labelCPUCount_2, 1, 0, 2, 1)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_16.addItem(spacerItem2, 2, 1, 1, 1)
        self.verticalLayout_3.addLayout(self.gridLayout_16)
        self.gridLayout.addWidget(self.groupBox_11, 0, 1, 1, 1)
        self.groupBox_7 = QtGui.QGroupBox(Form)
        font = QtGui.QFont()
        font.setUnderline(False)
        self.groupBox_7.setFont(font)
        self.groupBox_7.setObjectName(_fromUtf8("groupBox_7"))
        self.gridLayout_15 = QtGui.QGridLayout(self.groupBox_7)
        self.gridLayout_15.setObjectName(_fromUtf8("gridLayout_15"))
        self.labelElapsedTime_2 = QtGui.QLabel(self.groupBox_7)
        self.labelElapsedTime_2.setObjectName(_fromUtf8("labelElapsedTime_2"))
        self.gridLayout_15.addWidget(self.labelElapsedTime_2, 1, 0, 1, 2)
        self.label_108 = QtGui.QLabel(self.groupBox_7)
        self.label_108.setObjectName(_fromUtf8("label_108"))
        self.gridLayout_15.addWidget(self.label_108, 0, 0, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(277, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_15.addItem(spacerItem3, 1, 2, 1, 1)
        self.prgStatusProgress = QtGui.QProgressBar(self.groupBox_7)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.prgStatusProgress.sizePolicy().hasHeightForWidth())
        self.prgStatusProgress.setSizePolicy(sizePolicy)
        self.prgStatusProgress.setProperty("value", 24)
        self.prgStatusProgress.setObjectName(_fromUtf8("prgStatusProgress"))
        self.gridLayout_15.addWidget(self.prgStatusProgress, 0, 1, 1, 3)
        self.gridLayout.addWidget(self.groupBox_7, 1, 1, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.btnStatusLog.setText(_translate("Form", "Save Log", None))
        self.labelStatus_2.setText(_translate("Form", "Status Area", None))
        self.groupBox_11.setTitle(_translate("Form", "System Status", None))
        self.label_132.setText(_translate("Form", "Memory Usage", None))
        self.labelMaxMem_2.setText(_translate("Form", "Max Mem:", None))
        self.label_164.setText(_translate("Form", "CPU Usage", None))
        self.labelCPUCount_2.setText(_translate("Form", "CPU count:", None))
        self.groupBox_7.setTitle(_translate("Form", "Application Status", None))
        self.labelElapsedTime_2.setText(_translate("Form", "Elapsed Time: 0 H  0 M", None))
        self.label_108.setText(_translate("Form", "Progress", None))

