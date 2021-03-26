#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os

import matplotlib.patches as patches
from PIL import Image
import numpy as np
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog, QWidget, QVBoxLayout

# from PyQt5.QtWidgets import QMessageBox
from scipy.spatial import Voronoi, voronoi_plot_2d

import gui_CVT
import gui_dlg_option
import gui_dlg_post_option
from mplwidget import MplWidget
from CVT import PseudoRand
from RDF import rdf2d
from utils import create_png, fft2d, create_gds


class MyOptionDlg(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = gui_dlg_option.Ui_OptionDialog()
        self.ui.setupUi(self)
        self.working_dir = os.getcwd()
        self.ui.bt_working_dir.clicked.connect(self.set_working_dir)

    def set_working_dir(self):
        value = QFileDialog.getExistingDirectory(
            self, "Open Directory", self.working_dir, QFileDialog.ShowDirsOnly
        )
        if value:
            self.working_dir = value
        self.ui.bt_working_dir.setText(os.path.basename(self.working_dir))


class MyPostOptionDlg(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = gui_dlg_post_option.Ui_Dialog()
        self.ui.setupUi(self)


class GraphWindow(QWidget):
    def __init__(self, title="what is my name?"):
        super().__init__()
        self.setWindowTitle(title)
        layout = QVBoxLayout()
        self.mpl_graph = MplWidget()
        layout.addWidget(self.mpl_graph)
        self.setLayout(layout)


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
        self.ax = 50
        self.ay = 50
        self.pitch = 10

        self.ed_mainfield.valueChanged.connect(self.set_nbpt)
        self.ed_periode.valueChanged.connect(self.set_nbpt)
        self.ed_boundary.activated[str].connect(self.set_boundary)
        self.ed_delta.valueChanged.connect(self.set_delta)
        self.ed_geometry.activated[str].connect(self.set_geometry)
        self.bt_working_dir.clicked.connect(self.set_working_dir)
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

        self.myCVT = PseudoRand(
            nbpt=self.nbpt,
            datafile=self.init_filename,
            initial=self.geometry,
            boundary=self.boundary,
            delta=self.delta,
        )
        self.set_nbpt()
        self.rdf_distance = self.rdf_value = np.array([])
        self.save_cvt_data = False
        self.save_cvt_each_iteration = 1
        self.save_rdf_data = False
        self.save_converg = False
        self.save_png = False
        self.save_gds = False
        self.save_fft_img = False
        self.save_fft_data = False
        self.save_fht_img = False
        self.save_fht_data = False

        self.win_png = GraphWindow("PNG image")
        self.win_fft = GraphWindow("FFT image")
        self.win_fht = GraphWindow("FHT image")

        self.timer = QtCore.QTimer()
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.update_cvt)

    def set_win_title(self):
        title = (
            "RandomLight_" + self.boundary + "_" + self.geometry + "_" + str(self.nbpt)
        )
        self.setWindowTitle(title)

    def save_data(self, data, data_name, iteration=None, extension="csv"):
        filename = (
            data_name + "_" + self.boundary + "_" + self.geometry + "_" + str(self.nbpt)
        )
        if iteration:
            filename += "_" + iteration
        full_filename = os.path.join(self.working_dir, filename + "." + extension)
        np.savetxt(full_filename, data, delimiter=",")

    def set_working_dir(self):
        value = QFileDialog.getExistingDirectory(
            self, "Open Directory", self.working_dir, QFileDialog.ShowDirsOnly
        )
        if value:
            self.working_dir = value
        self.bt_working_dir.setText(os.path.basename(self.working_dir))

    def set_nbpt(self):
        self.L = self.ed_mainfield.value() * 1e3
        if self.geometry == "Square":
            self.d = self.ed_periode.value()
            nb = int(np.rint(self.L / self.d))
            self.nbpt = nb * nb

        elif self.geometry == "From_File":
            self.nbpt = self.myCVT.nbpt
        else:
            self.d = self.ed_periode.value()
            self.nbpt = int(np.ceil(2 / np.sqrt(3) * self.L ** 2 / self.d ** 2))

        self.set_win_title()
        self.myCVT.d = self.d / self.L
        self.myCVT.set_nbpt(self.nbpt)
        self.myCVT.set_initial(self.geometry)
        self.plot_cvt()
        self.plot_converg()
        self.plot_rdf()

    def set_delta(self):
        self.delta = self.ed_delta.value() / 100.0
        self.myCVT.set_delta(self.delta)
        self.plot_cvt()
        self.plot_rdf()

    def set_boundary(self, value):
        self.boundary = value
        self.myCVT.set_boundary(value)
        self.set_win_title()
        if self.boundary == "Free":
            self.ed_delta.setEnabled(False)
        else:
            self.ed_delta.setEnabled(True)
        self.plot_cvt()
        self.plot_rdf()

    def set_geometry(self, value):
        if value == "Random":
            self.ed_periode.setEnabled(True)
            self.lab_periode.setEnabled(True)
            self.lab_init_file.setHidden(True)
            self.ed_init_file.setHidden(True)
            self.geometry = value
        elif value == "Square":
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
            self.init_filename, _ = QFileDialog.getOpenFileName(
                self, "Data File", "", "Data Files(*.csv)"
            )
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
        self.myCVT.set_initial(self.geometry)
        self.set_nbpt()
        self.set_win_title()
        self.plot_cvt()
        self.plot_converg()
        self.plot_rdf()

    def set_option(self):
        dlg = MyOptionDlg(self)
        dlg.ui.ed_save_cvt_data.setChecked(self.save_cvt_data)
        dlg.ui.ed_each_cvt.setValue(self.save_cvt_each_iteration)
        dlg.ui.ed_save_rdf_data.setChecked(self.save_rdf_data)
        dlg.ui.ed_save_converg.setChecked(self.save_converg)
        dlg.working_dir = self.working_dir
        dlg.ui.bt_working_dir.setText(os.path.basename(dlg.working_dir))
        if dlg.exec():
            self.save_cvt_data = dlg.ui.ed_save_cvt_data.isChecked()
            self.save_cvt_each_iteration = dlg.ui.ed_each_cvt.value()
            self.save_rdf_data = dlg.ui.ed_save_rdf_data.isChecked()
            self.save_converg = dlg.ui.ed_save_converg.isChecked()
            self.working_dir = dlg.working_dir

    def set_post_option(self):
        dlg = MyPostOptionDlg(self)
        dlg.ui.ed_PNG.setChecked(self.save_png)
        dlg.ui.ed_GDS.setChecked(self.save_gds)
        dlg.ui.ed_FFT_img.setChecked(self.save_fft_img)
        dlg.ui.ed_FFT_data.setChecked(self.save_fft_data)
        dlg.ui.ed_FHT_img.setChecked(self.save_fht_img)
        dlg.ui.ed_FHT_data.setChecked(self.save_fht_data)
        dlg.ui.ed_mainfield.setText(self.L.__str__())
        dlg.ui.ed_distance.setText(self.d.__str__())
        dlg.ui.ed_ax.setValue(self.ax)
        dlg.ui.ed_ay.setValue(self.ay)
        dlg.ui.ed_pitch.setValue(self.pitch)

        if dlg.exec():
            self.save_png = dlg.ui.ed_PNG.isChecked()
            self.save_gds = dlg.ui.ed_GDS.isChecked()
            self.save_fft_img = dlg.ui.ed_FFT_img.isChecked()
            self.save_fft_data = dlg.ui.ed_FFT_data.isChecked()
            self.save_fht_img = dlg.ui.ed_FHT_img.isChecked()
            self.save_fht_data = dlg.ui.ed_FHT_data.isChecked()
            self.ax = dlg.ui.ed_ax.value()
            self.ay = dlg.ui.ed_ay.value()
            self.pitch = dlg.ui.ed_pitch.value()

    def start(self):
        self.iteration = 0
        self.myCVT.reset_convergence()
        self.bt_start.setEnabled(False)
        self.bt_pause.setEnabled(True)
        self.bt_stop.setEnabled(True)
        self.ed_mainfield.setEnabled(False)
        self.ed_periode.setEnabled(False)
        self.ed_boundary.setEnabled(False)
        self.ed_geometry.setEnabled(False)
        self.bt_working_dir.setEnabled(False)
        self.ed_delta.setEnabled(False)
        self.set_option()
        self.save_data(self.myCVT.pts, "CVT", "start")
        self.timer.start()

    def pause(self):
        self.bt_start.setEnabled(False)
        self.bt_stop.setEnabled(True)
        self.ed_mainfield.setEnabled(False)
        self.ed_periode.setEnabled(False)
        self.ed_boundary.setEnabled(False)
        self.ed_geometry.setEnabled(False)
        self.bt_working_dir.setEnabled(False)

        if self.bt_pause.text() == "Pause":
            self.bt_pause.setText("Continue")
            self.timer.stop()
        else:
            self.bt_pause.setText("Pause")
            self.timer.start()

    def stop(self):
        self.bt_start.setEnabled(True)
        self.bt_pause.setEnabled(False)
        self.bt_pause.setText("Pause")
        self.bt_stop.setEnabled(False)
        self.ed_mainfield.setEnabled(True)
        self.ed_periode.setEnabled(True)
        self.ed_boundary.setEnabled(True)
        self.ed_geometry.setEnabled(True)
        self.bt_working_dir.setEnabled(True)
        self.ed_delta.setEnabled(True)
        self.timer.stop()
        self.save_data(self.myCVT.pts, "CVT", "end")
        if self.save_converg:
            iteration = np.arange(0, len(self.myCVT.dist_mean), 1, dtype=int)
            data = np.vstack((iteration, self.myCVT.dist_mean, self.myCVT.dist_sigma)).T
            self.save_data(data, "MeanDistance", "end")
        if self.save_rdf_data:
            data = np.vstack((self.rdf_distance, self.rdf_value)).T
            self.save_data(data, "RDF", "end")
        data_png = create_png(self.myCVT.pts, self.L, self.ax, self.ay, self.pitch)
        self.set_post_option()
        if self.save_png:
            im = Image.fromarray(data_png)
            filename = (
                "Image_" + self.boundary + "_" + self.geometry + "_" + str(self.nbpt)
            )
            full_filename = os.path.join(self.working_dir, filename + ".png")
            im.save(full_filename)
            self.win_png.mpl_graph.canvas.ax.set_title("PNG image")
            self.win_png.mpl_graph.canvas.ax.imshow(data_png)
            self.win_png.show()

        if self.save_gds:
            filename = (
                "Layout_" + self.boundary + "_" + self.geometry + "_" + str(self.nbpt)
            )
            full_filename = os.path.join(self.working_dir, filename + ".gds")
            create_gds(self.myCVT.pts, self.L, self.ax, self.ay, full_filename)

        if self.save_fft_img or self.save_fft_data:
            data_freq, data_fft = fft2d(data_png)
            if self.save_fft_img:
                data_fft_img = np.log(np.abs(data_fft) ** 2)
                self.win_fft.mpl_graph.canvas.ax.set_title("FFT image")
                self.win_fft.mpl_graph.canvas.ax.pcolor(data_fft_img)
                self.win_fft.show()
                im = Image.fromarray(data_fft_img).convert("L")
                filename = (
                    "FFT_" + self.boundary + "_" + self.geometry + "_" + str(self.nbpt)
                )
                full_filename = os.path.join(self.working_dir, filename + ".png")
                print(type(data_fft_img))
                im.save(full_filename)
            if self.save_fft_data:
                self.save_data(np.real(data_fft), "FFT_Real")
                self.save_data(np.imag(data_fft), "FFT_Imag")

        if self.save_fht_img or self.save_fht_data:
            if self.save_fft_img:
                print("Save FHT img")
                self.win_fht.show()
            if self.save_fft_data:
                print("Save FHT data")

    def update_cvt(self):
        self.iteration += 1
        self.myCVT.iteration()
        self.plot_cvt()
        self.plot_converg()
        self.plot_rdf()
        self.myCVT.update()
        iteration = np.remainder(self.iteration, self.save_cvt_each_iteration)
        if self.save_cvt_data and iteration == 0:
            self.save_data(self.myCVT.pts, "CVT", str(self.iteration))

    def plot_cvt(self):
        data = self.myCVT
        mybox = patches.Rectangle(
            (0.0, 0.0),
            self.L,
            self.L,
            linewidth=1,
            edgecolor="y",
            facecolor="y",
            alpha=0.5,
        )
        self.mlp_cvt.canvas.ax.cla()
        self.mlp_cvt.canvas.ax.plot(
            data.site[:, 0] * self.L, data.site[:, 1] * self.L, "ro", label="Total"
        )
        voronoi_plot_2d(
            Voronoi(data.pts * self.L),
            show_vertices=False,
            line_alpha=0.2,
            ax=self.mlp_cvt.canvas.ax,
            point_size=13,
        )
        self.mlp_cvt.canvas.ax.add_patch(mybox)
        self.mlp_cvt.canvas.ax.set_aspect("equal", "box")
        self.mlp_cvt.canvas.ax.set_xlim(
            -self.delta * self.L, (1.0 + self.delta) * self.L
        )
        self.mlp_cvt.canvas.ax.set_ylim(
            -self.delta * self.L, (1.0 + self.delta) * self.L
        )
        self.mlp_rdf.canvas.ax.set_xlabel("Main field (nm)")
        self.mlp_rdf.canvas.ax.set_ylabel("Main field (nm)")
        self.mlp_cvt.canvas.ax.set_title(
            "CVT iteration %d, nb_pts = %d" % (self.iteration, len(data.pts))
        )
        self.mlp_cvt.canvas.draw()

    def plot_rdf(self):
        data = self.myCVT.pts
        self.rdf_distance, self.rdf_value = rdf2d(data, dr=0.005)
        self.mlp_rdf.canvas.ax.cla()
        self.mlp_rdf.canvas.ax.plot(self.rdf_distance * self.L, self.rdf_value, "bo-")
        self.mlp_rdf.canvas.ax.set_title("Radial Distribution Function (RDF)")
        self.mlp_rdf.canvas.ax.set_xlabel("Distance (nm)")
        self.mlp_rdf.canvas.ax.set_ylabel("RDF (a.u.)")
        self.mlp_rdf.canvas.ax.grid(True)
        self.mlp_rdf.canvas.draw()

    def plot_converg(self):
        data = self.myCVT
        y = data.dist_mean * self.L
        y_err = data.dist_sigma * self.L
        # dy = np.abs(np.append(np.nan, np.diff(y)))
        iteration = np.arange(0, len(y), 1)
        self.mlp_converg.canvas.ax.cla()
        self.mlp_converg.canvas.ax.errorbar(
            iteration, y, yerr=y_err, fmt="bo-", ecolor="c", capsize=5
        )
        self.mlp_rdf.canvas.ax.set_xlabel("Iteration")
        self.mlp_rdf.canvas.ax.set_ylabel("Mean distance (nm)")
        self.mlp_converg.canvas.ax.grid(True)
        self.mlp_converg.canvas.ax.set_title(
            r"Mean distance : %d $\pm$ %d" % (y[-1], y_err[-1])
        )
        self.mlp_converg.canvas.draw()


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
