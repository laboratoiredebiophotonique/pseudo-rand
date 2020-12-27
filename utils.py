# -*- coding: utf-8 -*-(
import os
import numpy as np 

import gdspy
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt


def create_PNG(datafile, L, ax, ay, pitch):

    data  = np.genfromtxt(datafile, delimiter=',')
    ptS   = np.ceil(data*L/pitch).astype('int') 
    W = H = np.ceil(L/pitch).astype('int')
    ax    = np.ceil(ax/pitch).astype('int')
    ay    = np.ceil(ay/pitch).astype('int')

    imgPng  = Image.new("L", (W, H), "white")
    pngdraw = ImageDraw.Draw(imgPng)

    for pt in ptS:
        x_UpperLeft   = pt[0] - ax
        y_UpperLeft   = pt[1] - ay
        x_BottomRight = pt[0] + ax
        y_BottomRight = pt[1] + ay

        pngdraw.ellipse( [x_UpperLeft, y_UpperLeft, x_BottomRight, y_BottomRight], fill="black")

    filename = os.path.splitext(datafile)[0]
    fullname = "%s_pitch_%d.png" % (filename, pitch)
    imgPng.save(fullname, "PNG")
    return fullname
    
def create_GDS(datafile, L, ax, ay, pitch):    
    data  = np.genfromtxt(datafile, delimiter=',')
    ptS   = np.ceil(data*L/pitch).astype('int') 
    ax    = np.ceil(ax/pitch).astype('int')
    ay    = np.ceil(ay/pitch).astype('int')

    myGds = gdspy.GdsLibrary(unit=1e-06)
    cell  = myGds.new_cell('NP_layer')
    print(ptS.shape)
    for pt in ptS:
        print(pt[0],', ', pt[1])
        NP = gdspy.Round((pt[0], pt[1]), ax, tolerance=0.1)
        cell.add(NP)
    filename = os.path.splitext(datafile)[0]
    fullname = "%s_pitch_%d.gds" % (filename, pitch)
    myGds.write_gds(fullname)
    return fullname


def fft2_PG(imgfile):

    imgPng = Image.open(imgfile)

    DataImg = np.array(imgPng)
    l, c = DataImg.shape
    myArr = DataImg[0:l//2, 0:c//2]

    L = 2*2**(np.ceil(np.log2(l)).astype('int'))
    C = 2*2**(np.ceil(np.log2(c)).astype('int'))

    IMGG=np.fft.fft2(myArr, (L,C))
    Te = 25#pitch
    Fe=1/Te;
    freq=np.linspace(-Fe/2, Fe/2, L)

#    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,3))
#    plt.plot( freq, np.fft.fftshift( np.abs(IMGG[0:L,0]) ) )
    plt.pcolor( np.fft.fftshift(np.abs(IMGG)) )
#    plt.xlim( (-6e-3, 6e-3) )
    plt.show()
    
    
