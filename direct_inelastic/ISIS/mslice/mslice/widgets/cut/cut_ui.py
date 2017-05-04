# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mslice\widgets\cut\cut.ui'
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
        Form.resize(402, 150)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_12 = QtGui.QGridLayout()
        self.gridLayout_12.setObjectName(_fromUtf8("gridLayout_12"))
        self.label_251 = QtGui.QLabel(Form)
        self.label_251.setObjectName(_fromUtf8("label_251"))
        self.gridLayout_12.addWidget(self.label_251, 0, 0, 1, 1)
        self.cmbCutAxis = QtGui.QComboBox(Form)
        self.cmbCutAxis.setObjectName(_fromUtf8("cmbCutAxis"))
        self.gridLayout_12.addWidget(self.cmbCutAxis, 0, 1, 1, 1)
        self.label_4 = QtGui.QLabel(Form)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_12.addWidget(self.label_4, 0, 2, 1, 1)
        self.lneCutStart = QtGui.QLineEdit(Form)
        self.lneCutStart.setObjectName(_fromUtf8("lneCutStart"))
        self.gridLayout_12.addWidget(self.lneCutStart, 0, 3, 1, 1)
        self.label_257 = QtGui.QLabel(Form)
        self.label_257.setObjectName(_fromUtf8("label_257"))
        self.gridLayout_12.addWidget(self.label_257, 0, 4, 1, 1)
        self.lneCutEnd = QtGui.QLineEdit(Form)
        self.lneCutEnd.setObjectName(_fromUtf8("lneCutEnd"))
        self.gridLayout_12.addWidget(self.lneCutEnd, 0, 5, 1, 1)
        self.label_256 = QtGui.QLabel(Form)
        self.label_256.setObjectName(_fromUtf8("label_256"))
        self.gridLayout_12.addWidget(self.label_256, 0, 6, 1, 1)
        self.lneCutStep = QtGui.QLineEdit(Form)
        self.lneCutStep.setText(_fromUtf8(""))
        self.lneCutStep.setObjectName(_fromUtf8("lneCutStep"))
        self.gridLayout_12.addWidget(self.lneCutStep, 0, 7, 1, 1)
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_12.addWidget(self.label_2, 1, 1, 1, 1)
        self.label_258 = QtGui.QLabel(Form)
        self.label_258.setObjectName(_fromUtf8("label_258"))
        self.gridLayout_12.addWidget(self.label_258, 1, 2, 1, 1)
        self.lneCutIntegrationStart = QtGui.QLineEdit(Form)
        self.lneCutIntegrationStart.setObjectName(_fromUtf8("lneCutIntegrationStart"))
        self.gridLayout_12.addWidget(self.lneCutIntegrationStart, 1, 3, 1, 1)
        self.label_259 = QtGui.QLabel(Form)
        self.label_259.setObjectName(_fromUtf8("label_259"))
        self.gridLayout_12.addWidget(self.label_259, 1, 4, 1, 1)
        self.lneCutIntegrationEnd = QtGui.QLineEdit(Form)
        self.lneCutIntegrationEnd.setObjectName(_fromUtf8("lneCutIntegrationEnd"))
        self.gridLayout_12.addWidget(self.lneCutIntegrationEnd, 1, 5, 1, 1)
        self.label_253 = QtGui.QLabel(Form)
        self.label_253.setObjectName(_fromUtf8("label_253"))
        self.gridLayout_12.addWidget(self.label_253, 1, 6, 1, 1)
        self.lneCutIntegrationWidth = QtGui.QLineEdit(Form)
        self.lneCutIntegrationWidth.setEnabled(True)
        self.lneCutIntegrationWidth.setText(_fromUtf8(""))
        self.lneCutIntegrationWidth.setObjectName(_fromUtf8("lneCutIntegrationWidth"))
        self.gridLayout_12.addWidget(self.lneCutIntegrationWidth, 1, 7, 1, 1)
        self.label_3 = QtGui.QLabel(Form)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_12.addWidget(self.label_3, 2, 1, 1, 1)
        self.label_264 = QtGui.QLabel(Form)
        self.label_264.setObjectName(_fromUtf8("label_264"))
        self.gridLayout_12.addWidget(self.label_264, 2, 2, 1, 1)
        self.lneEditCutIntensityStart = QtGui.QLineEdit(Form)
        self.lneEditCutIntensityStart.setObjectName(_fromUtf8("lneEditCutIntensityStart"))
        self.gridLayout_12.addWidget(self.lneEditCutIntensityStart, 2, 3, 1, 1)
        self.label_265 = QtGui.QLabel(Form)
        self.label_265.setObjectName(_fromUtf8("label_265"))
        self.gridLayout_12.addWidget(self.label_265, 2, 4, 1, 1)
        self.lneCutIntensityEnd = QtGui.QLineEdit(Form)
        self.lneCutIntensityEnd.setObjectName(_fromUtf8("lneCutIntensityEnd"))
        self.gridLayout_12.addWidget(self.lneCutIntensityEnd, 2, 5, 1, 1)
        self.label_5 = QtGui.QLabel(Form)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_12.addWidget(self.label_5, 3, 1, 1, 1)
        self.lneCutSmoothing = QtGui.QLineEdit(Form)
        self.lneCutSmoothing.setEnabled(False)
        self.lneCutSmoothing.setObjectName(_fromUtf8("lneCutSmoothing"))
        self.gridLayout_12.addWidget(self.lneCutSmoothing, 3, 2, 1, 2)
        self.btnCutPlot = QtGui.QPushButton(Form)
        self.btnCutPlot.setObjectName(_fromUtf8("btnCutPlot"))
        self.gridLayout_12.addWidget(self.btnCutPlot, 3, 6, 1, 2)
        self.btnCutSaveAscii = QtGui.QPushButton(Form)
        self.btnCutSaveAscii.setEnabled(False)
        self.btnCutSaveAscii.setObjectName(_fromUtf8("btnCutSaveAscii"))
        self.gridLayout_12.addWidget(self.btnCutSaveAscii, 4, 0, 1, 2)
        self.btnCutSaveToWorkspace = QtGui.QPushButton(Form)
        self.btnCutSaveToWorkspace.setEnabled(True)
        self.btnCutSaveToWorkspace.setObjectName(_fromUtf8("btnCutSaveToWorkspace"))
        self.gridLayout_12.addWidget(self.btnCutSaveToWorkspace, 4, 3, 1, 3)
        self.btnCutPlotOver = QtGui.QPushButton(Form)
        self.btnCutPlotOver.setObjectName(_fromUtf8("btnCutPlotOver"))
        self.gridLayout_12.addWidget(self.btnCutPlotOver, 4, 6, 1, 2)
        self.rdoCutNormToOne = QtGui.QCheckBox(Form)
        self.rdoCutNormToOne.setObjectName(_fromUtf8("rdoCutNormToOne"))
        self.gridLayout_12.addWidget(self.rdoCutNormToOne, 2, 6, 1, 2)
        self.gridLayout.addLayout(self.gridLayout_12, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.label_251.setText(_translate("Form", "along", None))
        self.label_4.setText(_translate("Form", "from", None))
        self.label_257.setText(_translate("Form", "to", None))
        self.label_256.setText(_translate("Form", "step", None))
        self.label_2.setText(_translate("Form", "Integrating", None))
        self.label_258.setText(_translate("Form", "from", None))
        self.label_259.setText(_translate("Form", "to", None))
        self.label_253.setText(_translate("Form", "width", None))
        self.label_3.setText(_translate("Form", "Intensity", None))
        self.label_264.setText(_translate("Form", "from", None))
        self.label_265.setText(_translate("Form", "to", None))
        self.label_5.setText(_translate("Form", "Smoothing", None))
        self.btnCutPlot.setText(_translate("Form", "Plot", None))
        self.btnCutSaveAscii.setText(_translate("Form", "Save Ascii", None))
        self.btnCutSaveToWorkspace.setText(_translate("Form", "Save to Workspace", None))
        self.btnCutPlotOver.setText(_translate("Form", "Plot Over", None))
        self.rdoCutNormToOne.setText(_translate("Form", "Norm to 1", None))

