#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pyuic5 gui_CVT.ui -o gui_CVT.py

import os

import matplotlib.patches as patches
import numpy as np
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog
# from PyQt5.QtWidgets import QMessageBox
from scipy.spatial import Voronoi, voronoi_plot_2d

import gui_CVT
import option_dlg
from CVT import PseudoRand


class MyOptionDlg(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = option_dlg.Ui_OptionDialog()
        self.ui.setupUi(self)


class MyWindow(QtWidgets.QMainWindow, gui_CVT.Ui_MainWindow):

    def __init__(self, parent=None):
        super(MyWindow, self).__init__(parent)
        self.setupUi(self)

        self.L = self.ed_mainfield.value() * 1e3
        self.d = self.ed_periode.value()
        self.nbpt = int(np.ceil(2 / np.sqrt(3) * self.L ** 2 / self.d ** 2))
        self.iteration = 0
        self.delta = 0.1
        self.init_filename = None
        self.working_dir = os.getcwd()
        self.from_scratch = True

        self.ed_mainfield.valueChanged.connect(self.set_nbpt)
        self.ed_periode.valueChanged.connect(self.set_nbpt)
        self.ed_boundary.activated[str].connect(self.set_boundary)
        self.ed_delta.valueChanged.connect(self.set_delta)
        self.ed_geometry.activated[str].connect(self.set_geometry)
        self.bt_working_dir.clicked.connect(self.set_working_dir)
        self.bt_option.clicked.connect(self.set_option)
        self.bt_start.clicked.connect(self.start)
        self.bt_pause.clicked.connect(self.pause)
        self.bt_stop.clicked.connect(self.stop)
        self.ed_geometry.setCurrentIndex(0)
        self.ed_boundary.setCurrentIndex(0)
        self.lab_init_file.setHidden(True)
        self.ed_init_file.setHidden(True)
        self.bt_working_dir.setText(os.path.basename(self.working_dir))

        self.boundary = self.ed_boundary.currentText()
        self.geometry = self.ed_geometry.currentText()

        self.myCVT = PseudoRand(nbpt=self.nbpt, datafile=self.init_filename, initial=self.geometry,
                                boundary=self.boundary, delta=self.delta)
        self.set_nbpt()
        self.save_cvt = False
        self.save_rdf = False
        self.save_converg = False
        self.timer = QtCore.QTimer()
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.update_cvt)

    def set_win_title(self):
        title = "RandomLight_" + self.boundary + "_" + self.geometry + "_" + str(self.nbpt)
        self.setWindowTitle(title)

    def set_working_dir(self):
        value = QFileDialog.getExistingDirectory(self, "Open Directory", self.working_dir, QFileDialog.ShowDirsOnly)
        if value:
            self.working_dir = value
        self.bt_working_dir.setText(os.path.basename(self.working_dir))

    def set_nbpt(self):
        self.from_scratch = True
        self.L = self.ed_mainfield.value() * 1e3
        if self.geometry == 'Square':
            self.d = self.ed_periode.value()
            nb = int(np.rint(self.L/self.d))
            self.nbpt = nb*nb

        elif self.geometry == 'From_File':
            self.nbpt = self.myCVT.nbpt
        else:
            self.d = self.ed_periode.value()
            self.nbpt = int(np.ceil(2 / np.sqrt(3) * self.L ** 2 / self.d ** 2))

        self.set_win_title()
        self.myCVT.d = self.d/self.L
        self.myCVT.set_nbpt(self.nbpt)
        self.myCVT.set_initial(self.geometry)
        self.plot_cvt()

    def set_delta(self):
        self.delta = self.ed_delta.value() / 100.0
        self.myCVT.set_delta(self.delta)
        self.plot_cvt()

    def set_boundary(self, value):
        if not self.from_scratch:
            self.set_nbpt()
        self.boundary = value
        self.myCVT.set_boundary(value)
        self.set_win_title()
        if self.boundary == 'Free':
            self.ed_delta.setEnabled(False)
        else:
            self.ed_delta.setEnabled(True)
        self.plot_cvt()

    def set_geometry(self, value):
        if value == 'Random':
            self.ed_periode.setEnabled(True)
            self.lab_periode.setEnabled(True)
            self.lab_init_file.setHidden(True)
            self.ed_init_file.setHidden(True)
            self.geometry = value
        elif value == 'Square':
            self.ed_periode.setEnabled(True)
            self.lab_periode.setEnabled(True)
            self.lab_init_file.setHidden(True)
            self.ed_init_file.setHidden(True)
            self.geometry = value
        elif value == "Hexagonal_Compact":
            self.ed_periode.setEnabled(False)
            self.lab_periode.setEnabled(False)
            self.lab_init_file.setHidden(True)
            self.ed_init_file.setHidden(True)
            self.ed_init_file.setHidden(True)
            self.geometry = value
        elif value == "From_File":
            self.init_filename, _ = QFileDialog.getOpenFileName(self, "Data File", "", "Data Files(*.csv)")
            if self.init_filename:
                self.myCVT.set_datafile(self.init_filename)
                self.ed_periode.setEnabled(False)
                self.lab_periode.setEnabled(False)
                self.lab_init_file.setHidden(False)
                self.ed_init_file.setHidden(False)
                self.ed_init_file.setText(os.path.basename(self.init_filename))
                self.geometry = value
            else:
                self.ed_geometry.setCurrentText(self.geometry)
        self.from_scratch = True
        self.myCVT.set_initial(self.geometry)
        self.set_nbpt()
        self.set_win_title()
        self.plot_cvt()

    def set_option(self):
        dlg = MyOptionDlg(self)
        dlg.ui.ed_save_cvt.setChecked(self.save_cvt)
        dlg.ui.ed_save_rdf.setChecked(self.save_rdf)
        dlg.ui.ed_save_converg.setChecked(self.save_converg)
        if dlg.exec():
            self.save_cvt = dlg.ui.ed_save_cvt.isChecked()
            self.save_rdf = dlg.ui.ed_save_rdf.isChecked()
            self.save_converg = dlg.ui.ed_save_converg.isChecked()

    def start(self):
        if not self.from_scratch:
            self.myCVT = PseudoRand(nbpt=self.nbpt, datafile=self.init_filename, initial=self.geometry,
                                    boundary=self.boundary, delta=self.delta)
        self.from_scratch = False
        self.iteration = 0
        self.bt_start.setEnabled(False)
        self.bt_pause.setEnabled(True)
        self.bt_stop.setEnabled(True)
        self.ed_mainfield.setEnabled(False)
        self.ed_periode.setEnabled(False)
        self.ed_boundary.setEnabled(False)
        self.ed_geometry.setEnabled(False)
        self.bt_working_dir.setEnabled(False)
        self.ed_delta.setEnabled(False)
        self.bt_option.setEnabled(False)

        self.mlp_cvt.canvas.ax.cla()
        self.mlp_converg.canvas.ax.cla()
        self.timer.start()

    def pause(self):
        self.bt_start.setEnabled(False)
        self.bt_stop.setEnabled(True)
        self.ed_mainfield.setEnabled(False)
        self.ed_periode.setEnabled(False)
        self.ed_boundary.setEnabled(False)
        self.ed_geometry.setEnabled(False)
        self.bt_working_dir.setEnabled(False)
        if self.bt_pause.text() == 'Pause':
            self.bt_pause.setText('Continue')
            self.timer.stop()
        else:
            self.bt_pause.setText('Pause')
            self.timer.start()

    def stop(self):
        self.bt_start.setEnabled(True)
        self.bt_pause.setEnabled(False)
        self.bt_pause.setText('Pause')
        self.bt_stop.setEnabled(False)
        self.ed_mainfield.setEnabled(True)
        self.ed_periode.setEnabled(True)
        self.ed_boundary.setEnabled(True)
        self.ed_geometry.setEnabled(True)
        self.bt_working_dir.setEnabled(True)
        self.ed_delta.setEnabled(True)
        self.bt_option.setEnabled(True)
        self.timer.stop()

    def update_cvt(self):
        # Drop off the first y element, append a new one.
        self.iteration += 1
        self.myCVT.iteration()
        self.plot_cvt()
        self.plot_converg()
        # Trigger the canvas to update and redraw.
        self.mlp_cvt.canvas.draw()
        self.mlp_converg.canvas.draw()
        self.myCVT.update()

    def plot_cvt(self):
        data = self.myCVT

        mybox = patches.Rectangle((-0.5, -0.5), 1.0, 1.0, linewidth=1, edgecolor='y', facecolor='y', alpha=0.5)
        self.mlp_cvt.canvas.ax.cla()
        self.mlp_cvt.canvas.ax.plot(data.site[:, 0] - 0.5, data.site[:, 1] - 0.5, 'ro', label='Total')
        voronoi_plot_2d(Voronoi(data.pts - [0.5, 0.5]), show_vertices=False, line_alpha=0.2, ax=self.mlp_cvt.canvas.ax)
        self.mlp_cvt.canvas.ax.add_patch(mybox)
        self.mlp_cvt.canvas.ax.set_aspect('equal', 'box')
        self.mlp_cvt.canvas.ax.set_xlim(-0.5 - self.delta, 0.5 + self.delta)
        self.mlp_cvt.canvas.ax.set_ylim(-0.5 - self.delta, 0.5 + self.delta)
        self.mlp_cvt.canvas.ax.set_title('CVT iteration %d, nb_pts = %d' % (self.iteration, len(data.pts)))
        self.mlp_cvt.canvas.draw()

    def plot_converg(self):
        data = self.myCVT
        y = data.dist_mean * self.L
        y_err = data.dist_sigma * self.L
        # dy = np.abs(np.append(np.nan, np.diff(y)))
        iteration = np.arange(0, len(y), 1)
        self.mlp_converg.canvas.ax.errorbar(iteration, y, yerr=y_err, fmt='bo-', ecolor='c', capsize=5)
        self.mlp_converg.canvas.ax.grid(True)
        self.mlp_converg.canvas.ax.set_title(r'distance : %d $\pm$ %d' % (y[-1], y_err[-1]))


if __name__ == '__main__':
    import sys

    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
