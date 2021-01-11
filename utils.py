# -*- coding: utf-8 -*-(
import numpy as np

# import gdspy
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


# def create_GDS(data, L, ax, ay, pitch):
#     ptS = np.ceil(data * L / pitch).astype('int')
#     ax = np.ceil(ax / pitch).astype('int')
#     ay = np.ceil(ay / pitch).astype('int')
#
#     myGds = gdspy.GdsLibrary(unit=1e-06)
#     cell = myGds.new_cell('NP_layer')
#     print(ptS.shape)
#     for pt in ptS:
#         print(pt[0], ', ', pt[1])
#         NP = gdspy.Round((pt[0], pt[1]), ax, tolerance=0.1)
#         cell.add(NP)
#     filename = os.path.splitext(datafile)[0]
#     fullname = "%s_pitch_%d.gds" % (filename, pitch)
#     myGds.write_gds(fullname)
#     return fullname
#
#
def fft2d(data):
    h, w = data.shape
    L = 2 * 2 ** (np.ceil(np.log2(h)).astype('int'))
    C = 2 * 2 ** (np.ceil(np.log2(w)).astype('int'))
    my_arr = data[0:h // 2, 0:w // 2]

    fft2d = np.fft.fft2(my_arr, (L, C))
    Te = 25
    Fe = 1 / Te
    freq = np.linspace(-Fe / 2, Fe / 2, L)
    return freq, np.fft.fftshift(np.abs(fft2d))
