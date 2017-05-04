# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mslice\widgets\workspacemanager\workspacemanager.ui'
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
        Form.resize(529, 225)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.btnWorkspaceSave = QtGui.QPushButton(Form)
        self.btnWorkspaceSave.setObjectName(_fromUtf8("btnWorkspaceSave"))
        self.gridLayout.addWidget(self.btnWorkspaceSave, 1, 0, 1, 1)
        self.btnWorkspaceRemove = QtGui.QPushButton(Form)
        self.btnWorkspaceRemove.setObjectName(_fromUtf8("btnWorkspaceRemove"))
        self.gridLayout.addWidget(self.btnWorkspaceRemove, 1, 1, 1, 1)
        spacerItem = QtGui.QSpacerItem(186, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 1, 2, 1, 1)
        self.groupBox_6 = QtGui.QGroupBox(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_6.sizePolicy().hasHeightForWidth())
        self.groupBox_6.setSizePolicy(sizePolicy)
        self.groupBox_6.setMinimumSize(QtCore.QSize(120, 0))
        self.groupBox_6.setFlat(True)
        self.groupBox_6.setObjectName(_fromUtf8("groupBox_6"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox_6)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.btnLoad = QtGui.QPushButton(self.groupBox_6)
        self.btnLoad.setObjectName(_fromUtf8("btnLoad"))
        self.gridLayout_2.addWidget(self.btnLoad, 0, 0, 1, 1)
        self.btnWorkspaceCompose = QtGui.QPushButton(self.groupBox_6)
        self.btnWorkspaceCompose.setEnabled(False)
        self.btnWorkspaceCompose.setObjectName(_fromUtf8("btnWorkspaceCompose"))
        self.gridLayout_2.addWidget(self.btnWorkspaceCompose, 1, 0, 1, 1)
        self.btnRename = QtGui.QPushButton(self.groupBox_6)
        self.btnRename.setObjectName(_fromUtf8("btnRename"))
        self.gridLayout_2.addWidget(self.btnRename, 2, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem1, 3, 0, 1, 1)
        self.gridLayout.addWidget(self.groupBox_6, 0, 3, 1, 1)
        self.lstWorkspaces = QtGui.QListWidget(Form)
        self.lstWorkspaces.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.lstWorkspaces.setObjectName(_fromUtf8("lstWorkspaces"))
        self.gridLayout.addWidget(self.lstWorkspaces, 0, 0, 1, 3)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Form", None))
        self.btnWorkspaceSave.setText(_translate("Form", "Save", None))
        self.btnWorkspaceRemove.setText(_translate("Form", "Remove", None))
        self.groupBox_6.setTitle(_translate("Form", "Workspace Options", None))
        self.btnLoad.setText(_translate("Form", "Load", None))
        self.btnWorkspaceCompose.setText(_translate("Form", "Compose", None))
        self.btnRename.setText(_translate("Form", "Rename", None))

