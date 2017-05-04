# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mslice\widgets\slice\slice.ui'
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
        Form.resize(400, 147)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_10 = QtGui.QGridLayout()
        self.gridLayout_10.setObjectName(_fromUtf8("gridLayout_10"))
        self.label_10 = QtGui.QLabel(Form)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.gridLayout_10.addWidget(self.label_10, 0, 0, 1, 1)
        self.cmbSliceXAxis = QtGui.QComboBox(Form)
        self.cmbSliceXAxis.setObjectName(_fromUtf8("cmbSliceXAxis"))
        self.gridLayout_10.addWidget(self.cmbSliceXAxis, 0, 1, 1, 1)
        self.label_12 = QtGui.QLabel(Form)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.gridLayout_10.addWidget(self.label_12, 0, 2, 1, 1)
        self.lneSliceXStart = QtGui.QLineEdit(Form)
        self.lneSliceXStart.setObjectName(_fromUtf8("lneSliceXStart"))
        self.gridLayout_10.addWidget(self.lneSliceXStart, 0, 3, 1, 1)
        self.label_267 = QtGui.QLabel(Form)
        self.label_267.setObjectName(_fromUtf8("label_267"))
        self.gridLayout_10.addWidget(self.label_267, 0, 4, 1, 1)
        self.lneSliceXEnd = QtGui.QLineEdit(Form)
        self.lneSliceXEnd.setObjectName(_fromUtf8("lneSliceXEnd"))
        self.gridLayout_10.addWidget(self.lneSliceXEnd, 0, 5, 1, 2)
        self.label_268 = QtGui.QLabel(Form)
        self.label_268.setObjectName(_fromUtf8("label_268"))
        self.gridLayout_10.addWidget(self.label_268, 0, 7, 1, 1)
        self.lneSliceXStep = QtGui.QLineEdit(Form)
        self.lneSliceXStep.setText(_fromUtf8(""))
        self.lneSliceXStep.setObjectName(_fromUtf8("lneSliceXStep"))
        self.gridLayout_10.addWidget(self.lneSliceXStep, 0, 8, 1, 1)
        self.label_11 = QtGui.QLabel(Form)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridLayout_10.addWidget(self.label_11, 1, 0, 1, 1)
        self.cmbSliceYAxis = QtGui.QComboBox(Form)
        self.cmbSliceYAxis.setObjectName(_fromUtf8("cmbSliceYAxis"))
        self.gridLayout_10.addWidget(self.cmbSliceYAxis, 1, 1, 1, 1)
        self.label_271 = QtGui.QLabel(Form)
        self.label_271.setObjectName(_fromUtf8("label_271"))
        self.gridLayout_10.addWidget(self.label_271, 1, 2, 1, 1)
        self.lneSliceYStart = QtGui.QLineEdit(Form)
        self.lneSliceYStart.setObjectName(_fromUtf8("lneSliceYStart"))
        self.gridLayout_10.addWidget(self.lneSliceYStart, 1, 3, 1, 1)
        self.label_270 = QtGui.QLabel(Form)
        self.label_270.setObjectName(_fromUtf8("label_270"))
        self.gridLayout_10.addWidget(self.label_270, 1, 4, 1, 1)
        self.lneSliceYEnd = QtGui.QLineEdit(Form)
        self.lneSliceYEnd.setObjectName(_fromUtf8("lneSliceYEnd"))
        self.gridLayout_10.addWidget(self.lneSliceYEnd, 1, 5, 1, 2)
        self.label_255 = QtGui.QLabel(Form)
        self.label_255.setObjectName(_fromUtf8("label_255"))
        self.gridLayout_10.addWidget(self.label_255, 1, 7, 1, 1)
        self.lneSliceYStep = QtGui.QLineEdit(Form)
        self.lneSliceYStep.setText(_fromUtf8(""))
        self.lneSliceYStep.setObjectName(_fromUtf8("lneSliceYStep"))
        self.gridLayout_10.addWidget(self.lneSliceYStep, 1, 8, 1, 1)
        self.label_13 = QtGui.QLabel(Form)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.gridLayout_10.addWidget(self.label_13, 2, 0, 1, 2)
        self.label_269 = QtGui.QLabel(Form)
        self.label_269.setObjectName(_fromUtf8("label_269"))
        self.gridLayout_10.addWidget(self.label_269, 2, 2, 1, 1)
        self.lneSliceIntensityStart = QtGui.QLineEdit(Form)
        self.lneSliceIntensityStart.setObjectName(_fromUtf8("lneSliceIntensityStart"))
        self.gridLayout_10.addWidget(self.lneSliceIntensityStart, 2, 3, 1, 1)
        self.label_266 = QtGui.QLabel(Form)
        self.label_266.setObjectName(_fromUtf8("label_266"))
        self.gridLayout_10.addWidget(self.label_266, 2, 4, 1, 1)
        self.lneSliceIntensityEnd = QtGui.QLineEdit(Form)
        self.lneSliceIntensityEnd.setObjectName(_fromUtf8("lneSliceIntensityEnd"))
        self.gridLayout_10.addWidget(self.lneSliceIntensityEnd, 2, 5, 1, 2)
        self.label_14 = QtGui.QLabel(Form)
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.gridLayout_10.addWidget(self.label_14, 3, 0, 1, 2)
        self.lneSliceSmoothing = QtGui.QLineEdit(Form)
        self.lneSliceSmoothing.setEnabled(False)
        self.lneSliceSmoothing.setObjectName(_fromUtf8("lneSliceSmoothing"))
        self.gridLayout_10.addWidget(self.lneSliceSmoothing, 3, 2, 1, 2)
        self.cmbSliceColormap = QtGui.QComboBox(Form)
        self.cmbSliceColormap.setObjectName(_fromUtf8("cmbSliceColormap"))
        self.gridLayout_10.addWidget(self.cmbSliceColormap, 3, 8, 1, 1)
        self.label_15 = QtGui.QLabel(Form)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.gridLayout_10.addWidget(self.label_15, 3, 7, 1, 1)
        self.btnSliceDisplay = QtGui.QPushButton(Form)
        self.btnSliceDisplay.setObjectName(_fromUtf8("btnSliceDisplay"))
        self.gridLayout_10.addWidget(self.btnSliceDisplay, 4, 8, 1, 1)
        self.rdoSliceNormToOne = QtGui.QCheckBox(Form)
        self.rdoSliceNormToOne.setObjectName(_fromUtf8("rdoSliceNormToOne"))
        self.gridLayout_10.addWidget(self.rdoSliceNormToOne, 2, 7, 1, 2)
        self.gridLayout.addLayout(self.gridLayout_10, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.label_10.setText(_translate("Form", "x", None))
        self.label_12.setText(_translate("Form", "from", None))
        self.label_267.setText(_translate("Form", "to", None))
        self.label_268.setText(_translate("Form", "step", None))
        self.label_11.setText(_translate("Form", "y", None))
        self.label_271.setText(_translate("Form", "from", None))
        self.label_270.setText(_translate("Form", "to", None))
        self.label_255.setText(_translate("Form", "step", None))
        self.label_13.setText(_translate("Form", "Intensity", None))
        self.label_269.setText(_translate("Form", "from", None))
        self.label_266.setText(_translate("Form", "to", None))
        self.label_14.setText(_translate("Form", "Smoothing", None))
        self.label_15.setText(_translate("Form", "Colormap", None))
        self.btnSliceDisplay.setText(_translate("Form", "Display", None))
        self.rdoSliceNormToOne.setText(_translate("Form", "Norm to 1", None))

