# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'option_Dlg.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_OptionDialog(object):
    def setupUi(self, OptionDialog):
        OptionDialog.setObjectName("OptionDialog")
        OptionDialog.resize(517, 374)
        self.buttonBox = QtWidgets.QDialogButtonBox(OptionDialog)
        self.buttonBox.setGeometry(QtCore.QRect(170, 330, 341, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.ed_save_cvt = QtWidgets.QCheckBox(OptionDialog)
        self.ed_save_cvt.setGeometry(QtCore.QRect(40, 40, 211, 23))
        self.ed_save_cvt.setObjectName("ed_save_cvt")
        self.ed_save_rdf = QtWidgets.QCheckBox(OptionDialog)
        self.ed_save_rdf.setGeometry(QtCore.QRect(40, 80, 211, 23))
        self.ed_save_rdf.setObjectName("ed_save_rdf")
        self.ed_save_converg = QtWidgets.QCheckBox(OptionDialog)
        self.ed_save_converg.setGeometry(QtCore.QRect(40, 120, 211, 23))
        self.ed_save_converg.setObjectName("ed_save_converg")

        self.retranslateUi(OptionDialog)
        self.buttonBox.accepted.connect(OptionDialog.accept)
        self.buttonBox.rejected.connect(OptionDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(OptionDialog)

    def retranslateUi(self, OptionDialog):
        _translate = QtCore.QCoreApplication.translate
        OptionDialog.setWindowTitle(_translate("OptionDialog", "Save options"))
        self.ed_save_cvt.setText(_translate("OptionDialog", "save CVT on each iteration"))
        self.ed_save_rdf.setText(_translate("OptionDialog", "save RDF on each iteration"))
        self.ed_save_converg.setText(_translate("OptionDialog", "save convergence criteria"))
