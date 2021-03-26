#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

MIN_FLOAT = np.finfo('float').eps

pil_img = np.array(Image.open('img_256.png'))
mp_img = mpimg.imread('img_256.png')

IJ_FFT = np.genfromtxt('FFT_of_img_256.dat', delimiter='\t')
IJ_FHT = np.genfromtxt('FHT_of_img_256.dat', delimiter='\t')

N = 256
P = 256

# Normalisation grayscale de ImageJ
# IJ_FHT[IJ_FHT == 0] = MIN_FLOAT

IJ_FHT_NORM = np.log(np.abs(IJ_FHT))

MIN_FHT_IJ = np.min(np.min(IJ_FHT))
MAX_FHT_IJ = np.max(np.max(IJ_FHT))

if MIN_FHT_IJ < 1.0:
    MIN_FHT_IJ = 0.0
else:
    MIN_FHT_IJ = np.log(MIN_FHT_IJ)
MAX_FHT_IJ = np.log(MAX_FHT_IJ)
Scale = 253.0 / (MAX_FHT_IJ - MIN_FHT_IJ)

#IJ_FHT[IJ_FHT < 1.0] = 0.0
#IJ_FHT[IJ_FHT >= 1.0] = 0.0
IJ_FHT_NORM = np.log(IJ_FHT)
IJ_FHT_NORM = (IJ_FHT_NORM - MIN_FHT_IJ) * Scale + 0.5 + 1

a = 1.0 / (MAX_FHT_IJ - MIN_FHT_IJ)
b = -MIN_FHT_IJ / (MAX_FHT_IJ - MIN_FHT_IJ)

IJ_FHT_NORM = (a * IJ_FHT_NORM + b) * 253

fig_IJ = plt.figure()
axe_IJ = fig_IJ.subplots(2, 2).ravel()

map_IJ_FFT = axe_IJ[0].pcolor(IJ_FFT, cmap='gray')
axe_IJ[0].set_title('FFT ImageJ')
fig_IJ.colorbar(map_IJ_FFT, ax=axe_IJ[0])

map_IJ_FHT = axe_IJ[1].pcolor(IJ_FHT_NORM, cmap='gray')
axe_IJ[0].set_title('FHT ImageJ')
fig_IJ.colorbar(map_IJ_FHT, ax=axe_IJ[1])

# axe_IJ[2].hist(IJ_FFT.ravel(), bins=255)
# axe_IJ[2].set_title('FFT ImageJ')
# axe_IJ[3].hist(IJ_FHT_NORM.ravel(), bins=255)
# axe_IJ[3].set_title('FHT NORM ImageJ')
# fig_FHT_HISTO = plt.figure()
# axe_FHT_HISTO = fig_FHT_HISTO.subplots(2, 2).ravel()
# h_FHT_PY = axe_FHT_HISTO[0].hist( FHT_PY.ravel(), bins=255 )
# h_FHT_IJ = axe_FHT_HISTO[1].hist( FHT_IJ.ravel(), bins=255 )
# h_PS_IJ  = axe_FHT_HISTO[2].hist(NORM_FHT_IJ.ravel(), bins=255)
# legend('FHT Matlab', 'FHT ImageJ', 'PS ImageJ');

plt.show()
