# -*- coding: utf-8 -*-(
import numpy as np
import gdspy
from PIL import Image, ImageDraw


def create_png(data, scale, ax, ay, pitch):
    scale /= pitch
    pts = np.ceil(data * scale).astype('int')
    w = h = np.ceil(scale).astype('int')
    ax = np.ceil(ax / pitch).astype('int')
    ay = np.ceil(ay / pitch).astype('int')

    img_png = Image.new("L", (w, h), "white")
    png_draw = ImageDraw.Draw(img_png)

    for pt in pts:
        x_upper_left = pt[0] - ax
        y_upper_left = pt[1] - ay
        x_bottom_right = pt[0] + ax
        y_bottom_right = pt[1] + ay

        png_draw.ellipse([x_upper_left, y_upper_left, x_bottom_right, y_bottom_right], fill="black")
    return np.asarray(img_png)


def create_gds(data, scale, ax, ay, fullname):
    pts = np.ceil(data * scale).astype('int')
    my_gds = gdspy.GdsLibrary(unit=1e-09)
    cell = my_gds.new_cell('NP_layer')
    for pt in pts:
        motif = gdspy.Round((pt[0], pt[1]), (ax, ay), layer=1, number_of_points=24)
        cell.add(motif)
    my_gds.write_gds(fullname)
    return fullname


def fft2d(data):
    h, w = data.shape
    tmp_h = 2 * 2 ** (np.ceil(np.log2(h)).astype('int'))
    tmp_w = 2 * 2 ** (np.ceil(np.log2(w)).astype('int'))
    length = max(tmp_h, tmp_w)
    data_fft = data[0:length, 0:length]

    fft_2d = np.fft.fft2(data_fft, (length, length))
    t_e = 25
    f_e = 1 / t_e
    freq = np.linspace(-f_e / 2, f_e / 2, length)
    return freq, np.fft.fftshift(np.fft.fftshift(fft_2d))
