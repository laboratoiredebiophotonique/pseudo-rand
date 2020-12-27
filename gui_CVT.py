# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui_CVT.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.setEnabled(True)
        MainWindow.resize(930, 650)
        font = QtGui.QFont()
        font.setFamily("Comic Sans MS")
        font.setBold(True)
        font.setWeight(75)
        MainWindow.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../../../Images/favicon.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setStyleSheet("")
        self.centralWidget = QtWidgets.QWidget(MainWindow)
        self.centralWidget.setObjectName("centralWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralWidget)
        self.gridLayout.setContentsMargins(10, 10, 10, 10)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setObjectName("gridLayout")
        self.frame = QtWidgets.QFrame(self.centralWidget)
        self.frame.setEnabled(True)
        self.frame.setStyleSheet("")
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.frame)
        self.gridLayout_2.setContentsMargins(10, 10, 10, 10)
        self.gridLayout_2.setSpacing(10)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.ed_boundary = QtWidgets.QComboBox(self.frame)
        self.ed_boundary.setEditable(False)
        self.ed_boundary.setObjectName("ed_boundary")
        self.ed_boundary.addItem("")
        self.ed_boundary.addItem("")
        self.ed_boundary.addItem("")
        self.gridLayout_2.addWidget(self.ed_boundary, 0, 3, 1, 1)
        self.lab_delta = QtWidgets.QLabel(self.frame)
        self.lab_delta.setEnabled(True)
        self.lab_delta.setMinimumSize(QtCore.QSize(61, 0))
        self.lab_delta.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lab_delta.setObjectName("lab_delta")
        self.gridLayout_2.addWidget(self.lab_delta, 0, 4, 1, 1)
        self.lab_mainfield = QtWidgets.QLabel(self.frame)
        self.lab_mainfield.setMinimumSize(QtCore.QSize(0, 30))
        self.lab_mainfield.setMaximumSize(QtCore.QSize(200, 16777215))
        self.lab_mainfield.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lab_mainfield.setObjectName("lab_mainfield")
        self.gridLayout_2.addWidget(self.lab_mainfield, 0, 0, 1, 1)
        self.lab_working_dir = QtWidgets.QLabel(self.frame)
        self.lab_working_dir.setEnabled(True)
        self.lab_working_dir.setMinimumSize(QtCore.QSize(61, 0))
        self.lab_working_dir.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lab_working_dir.setObjectName("lab_working_dir")
        self.gridLayout_2.addWidget(self.lab_working_dir, 1, 4, 1, 1)
        self.ed_geometry = QtWidgets.QComboBox(self.frame)
        self.ed_geometry.setEditable(False)
        self.ed_geometry.setObjectName("ed_geometry")
        self.ed_geometry.addItem("")
        self.ed_geometry.addItem("")
        self.ed_geometry.addItem("")
        self.ed_geometry.addItem("")
        self.gridLayout_2.addWidget(self.ed_geometry, 1, 3, 1, 1)
        self.bt_working_dir = QtWidgets.QPushButton(self.frame)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/icons/folder_2"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_working_dir.setIcon(icon1)
        self.bt_working_dir.setObjectName("bt_working_dir")
        self.gridLayout_2.addWidget(self.bt_working_dir, 1, 5, 1, 1)
        self.ed_mainfield = QtWidgets.QDoubleSpinBox(self.frame)
        self.ed_mainfield.setMinimumSize(QtCore.QSize(0, 22))
        self.ed_mainfield.setDecimals(0)
        self.ed_mainfield.setMinimum(1.0)
        self.ed_mainfield.setMaximum(500.0)
        self.ed_mainfield.setSingleStep(10.0)
        self.ed_mainfield.setProperty("value", 5.0)
        self.ed_mainfield.setObjectName("ed_mainfield")
        self.gridLayout_2.addWidget(self.ed_mainfield, 0, 1, 1, 1)
        self.ed_delta = QtWidgets.QSlider(self.frame)
        self.ed_delta.setMaximum(100)
        self.ed_delta.setSingleStep(8)
        self.ed_delta.setProperty("value", 10)
        self.ed_delta.setOrientation(QtCore.Qt.Horizontal)
        self.ed_delta.setTickPosition(QtWidgets.QSlider.NoTicks)
        self.ed_delta.setObjectName("ed_delta")
        self.gridLayout_2.addWidget(self.ed_delta, 0, 5, 1, 1)
        self.lab_boundary = QtWidgets.QLabel(self.frame)
        self.lab_boundary.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lab_boundary.setObjectName("lab_boundary")
        self.gridLayout_2.addWidget(self.lab_boundary, 0, 2, 1, 1)
        self.lab_geometry = QtWidgets.QLabel(self.frame)
        self.lab_geometry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lab_geometry.setObjectName("lab_geometry")
        self.gridLayout_2.addWidget(self.lab_geometry, 1, 2, 1, 1)
        self.lab_init_file = QtWidgets.QLabel(self.frame)
        self.lab_init_file.setEnabled(True)
        self.lab_init_file.setMinimumSize(QtCore.QSize(30, 30))
        self.lab_init_file.setMaximumSize(QtCore.QSize(230, 30))
        self.lab_init_file.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.lab_init_file.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.lab_init_file.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lab_init_file.setObjectName("lab_init_file")
        self.gridLayout_2.addWidget(self.lab_init_file, 2, 0, 1, 2)
        self.lab_periode = QtWidgets.QLabel(self.frame)
        self.lab_periode.setMinimumSize(QtCore.QSize(0, 21))
        self.lab_periode.setMaximumSize(QtCore.QSize(200, 16777215))
        self.lab_periode.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lab_periode.setObjectName("lab_periode")
        self.gridLayout_2.addWidget(self.lab_periode, 1, 0, 1, 1)
        self.bt_option = QtWidgets.QPushButton(self.frame)
        self.bt_option.setMaximumSize(QtCore.QSize(16777215, 30))
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/icons/options"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_option.setIcon(icon2)
        self.bt_option.setObjectName("bt_option")
        self.gridLayout_2.addWidget(self.bt_option, 2, 5, 1, 1)
        self.ed_periode = QtWidgets.QDoubleSpinBox(self.frame)
        self.ed_periode.setMinimumSize(QtCore.QSize(0, 21))
        self.ed_periode.setDecimals(0)
        self.ed_periode.setMinimum(100.0)
        self.ed_periode.setMaximum(10000000.0)
        self.ed_periode.setSingleStep(10.0)
        self.ed_periode.setProperty("value", 500.0)
        self.ed_periode.setObjectName("ed_periode")
        self.gridLayout_2.addWidget(self.ed_periode, 1, 1, 1, 1)
        self.ed_init_file = QtWidgets.QLabel(self.frame)
        self.ed_init_file.setMinimumSize(QtCore.QSize(0, 30))
        self.ed_init_file.setMaximumSize(QtCore.QSize(500, 100))
        self.ed_init_file.setFrameShape(QtWidgets.QFrame.Box)
        self.ed_init_file.setObjectName("ed_init_file")
        self.gridLayout_2.addWidget(self.ed_init_file, 2, 2, 1, 2)
        self.gridLayout.addWidget(self.frame, 0, 0, 1, 2)
        self.frame_2 = QtWidgets.QFrame(self.centralWidget)
        self.frame_2.setMinimumSize(QtCore.QSize(0, 50))
        self.frame_2.setMaximumSize(QtCore.QSize(16777215, 50))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.frame_2)
        self.gridLayout_3.setContentsMargins(10, 10, 10, 10)
        self.gridLayout_3.setSpacing(6)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.bt_stop = QtWidgets.QPushButton(self.frame_2)
        self.bt_stop.setEnabled(False)
        self.bt_stop.setMaximumSize(QtCore.QSize(16777215, 100))
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/icons/stop_1"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_stop.setIcon(icon3)
        self.bt_stop.setObjectName("bt_stop")
        self.gridLayout_3.addWidget(self.bt_stop, 0, 2, 1, 1)
        self.bt_pause = QtWidgets.QPushButton(self.frame_2)
        self.bt_pause.setEnabled(False)
        self.bt_pause.setMaximumSize(QtCore.QSize(16777215, 100))
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/icons/pause_1"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_pause.setIcon(icon4)
        self.bt_pause.setObjectName("bt_pause")
        self.gridLayout_3.addWidget(self.bt_pause, 0, 1, 1, 1)
        self.bt_start = QtWidgets.QPushButton(self.frame_2)
        self.bt_start.setMaximumSize(QtCore.QSize(16777215, 30))
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(":/icons/start_1"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.bt_start.setIcon(icon5)
        self.bt_start.setObjectName("bt_start")
        self.gridLayout_3.addWidget(self.bt_start, 0, 0, 1, 1)
        self.gridLayout.addWidget(self.frame_2, 1, 0, 1, 2)
        self.mlp_cvt = MplWidget(self.centralWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mlp_cvt.sizePolicy().hasHeightForWidth())
        self.mlp_cvt.setSizePolicy(sizePolicy)
        self.mlp_cvt.setStyleSheet("background-color: rgb(252, 233, 79);")
        self.mlp_cvt.setObjectName("mlp_cvt")
        self.gridLayout.addWidget(self.mlp_cvt, 2, 0, 2, 1)
        self.mlp_converg = MplWidget(self.centralWidget)
        self.mlp_converg.setStyleSheet("background-color: rgb(252, 233, 79);")
        self.mlp_converg.setObjectName("mlp_converg")
        self.gridLayout.addWidget(self.mlp_converg, 2, 1, 1, 1)
        self.mlp_rdf = MplWidget(self.centralWidget)
        self.mlp_rdf.setStyleSheet("background-color: rgb(252, 233, 79);")
        self.mlp_rdf.setObjectName("mlp_rdf")
        self.gridLayout.addWidget(self.mlp_rdf, 3, 1, 1, 1)
        self.mlp_cvt.raise_()
        self.mlp_converg.raise_()
        self.mlp_rdf.raise_()
        self.frame.raise_()
        self.frame_2.raise_()
        MainWindow.setCentralWidget(self.centralWidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "RandomLight"))
        self.ed_boundary.setItemText(0, _translate("MainWindow", "Periodic"))
        self.ed_boundary.setItemText(1, _translate("MainWindow", "Free"))
        self.ed_boundary.setItemText(2, _translate("MainWindow", "Mirror"))
        self.lab_delta.setText(_translate("MainWindow", "Delta :"))
        self.lab_mainfield.setText(_translate("MainWindow", "MainField length:"))
        self.lab_working_dir.setText(_translate("MainWindow", "Work directory:"))
        self.ed_geometry.setCurrentText(_translate("MainWindow", "Random"))
        self.ed_geometry.setItemText(0, _translate("MainWindow", "Random"))
        self.ed_geometry.setItemText(1, _translate("MainWindow", "Square"))
        self.ed_geometry.setItemText(2, _translate("MainWindow", "Hexagonal_Compact"))
        self.ed_geometry.setItemText(3, _translate("MainWindow", "From_File"))
        self.bt_working_dir.setText(_translate("MainWindow", "Working directory"))
        self.ed_mainfield.setSuffix(_translate("MainWindow", " µm"))
        self.ed_delta.setWhatsThis(_translate("MainWindow", "<html><head/><body><p>Bonjour</p></body></html>"))
        self.lab_boundary.setText(_translate("MainWindow", "Boundary cond:"))
        self.lab_geometry.setText(_translate("MainWindow", "Geometry:"))
        self.lab_init_file.setWhatsThis(_translate("MainWindow", "<html><head/><body><p><br/></p></body></html>"))
        self.lab_init_file.setText(_translate("MainWindow", "Inital points from file:"))
        self.lab_periode.setText(_translate("MainWindow", "Mean period:"))
        self.bt_option.setText(_translate("MainWindow", "Save options..."))
        self.ed_periode.setSuffix(_translate("MainWindow", " nm"))
        self.ed_init_file.setText(_translate("MainWindow", "None"))
        self.bt_stop.setText(_translate("MainWindow", "Stop"))
        self.bt_pause.setText(_translate("MainWindow", "Pause"))
        self.bt_start.setText(_translate("MainWindow", "Start from scratch"))
from mplwidget import MplWidget
import ressources
