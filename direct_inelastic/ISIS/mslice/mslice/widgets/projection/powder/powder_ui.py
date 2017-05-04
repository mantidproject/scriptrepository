# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mslice\widgets\projection\powder\powder.ui'
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
        Form.resize(400, 131)
        self.gridLayout_2 = QtGui.QGridLayout(Form)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.groupBox = QtGui.QGroupBox(Form)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.verticalLayout_4 = QtGui.QVBoxLayout(self.groupBox)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setSizeConstraint(QtGui.QLayout.SetNoConstraint)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_208 = QtGui.QLabel(self.groupBox)
        self.label_208.setObjectName(_fromUtf8("label_208"))
        self.gridLayout.addWidget(self.label_208, 0, 0, 1, 1)
        self.cmbPowderU1 = QtGui.QComboBox(self.groupBox)
        self.cmbPowderU1.setObjectName(_fromUtf8("cmbPowderU1"))
        self.gridLayout.addWidget(self.cmbPowderU1, 0, 1, 1, 1)
        self.label_210 = QtGui.QLabel(self.groupBox)
        self.label_210.setObjectName(_fromUtf8("label_210"))
        self.gridLayout.addWidget(self.label_210, 1, 0, 1, 1)
        self.cmbPowderU2 = QtGui.QComboBox(self.groupBox)
        self.cmbPowderU2.setObjectName(_fromUtf8("cmbPowderU2"))
        self.gridLayout.addWidget(self.cmbPowderU2, 1, 1, 1, 1)
        self.label = QtGui.QLabel(self.groupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 1, 2, 1, 1)
        self.cmbPowderUnits = QtGui.QComboBox(self.groupBox)
        self.cmbPowderUnits.setObjectName(_fromUtf8("cmbPowderUnits"))
        self.gridLayout.addWidget(self.cmbPowderUnits, 1, 3, 1, 1)
        self.btnPowderCalculateProjection = QtGui.QPushButton(self.groupBox)
        self.btnPowderCalculateProjection.setObjectName(_fromUtf8("btnPowderCalculateProjection"))
        self.gridLayout.addWidget(self.btnPowderCalculateProjection, 2, 1, 1, 1)
        self.verticalLayout_4.addLayout(self.gridLayout)
        self.gridLayout_2.addWidget(self.groupBox, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.groupBox.setTitle(_translate("Form", "Projection", None))
        self.label_208.setText(_translate("Form", "u1", None))
        self.label_210.setText(_translate("Form", "u2", None))
        self.label.setText(_translate("Form", "Units", None))
        self.btnPowderCalculateProjection.setText(_translate("Form", "Calculate Projection", None))

