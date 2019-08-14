import ctypes
import numpy as np
import os
from ptu_wrapper import *
import matplotlib.pyplot as plt
from PIL import Image


def genImABfromptu(fname, uselines = np.ones(1, dtype = np.int), xbinning = 1, ybinning= 1, gate = 3228):
    NumRecords = ptuHeader_wrap (fname)
    eventN, tac, t, can = ptu_wrap(fname, NumRecords)
    root, file = os.path.split(fname)
    name, _ = os.path.splitext(file)
    header_name = os.path.join(root, b"header", name + b".txt")
    print('number of records is ' + str(NumRecords))

    dimX, dimY, dwelltime, counttime = read_header(header_name)
    dimX = int(dimX / xbinning)
    dimY = int(dimY / ybinning)
    dwelltime *= xbinning
    imA, imB = SplitOnTacs_wrap(eventN, tac, t, can, dimX, dimY, dwelltime, counttime, 
                                NumRecords, gate = gate, uselines = uselines)
    print("total image intensity is "  + str(np.sum(imA+imB)))

    try:
        os.mkdir(os.path.join(root, file[:-4]))
        print('storing result in new folder')
    except:
        print('overwriting result in existing folder')

    im = Image.fromarray(imA)
    outname = str(os.path.join(root, file[:-4] + b'\\imA.tiff'), encoding='utf-8')
    im.save(outname)
    im = Image.fromarray(imB)
    outname = str(os.path.join(root, file[:-4] + b'\\imB.tiff'), encoding='utf-8')
    im.save(outname)
    return imA, imB
    
def initLifetime(fname, uselines = np.array([1,0]), Gchan = np.array([0,2]), \
    Rchan = np.array([1,3]), Ychan = np.array([1,3]), ntacs = 256):
    """initialisation routine for lifetime image manipulation class.
    Loads the .ptu file located at fname.
    Uselines codes for the used line steps. 
        0 denotes that the line is ignored
        1 denotes that the line is FRET sensitized
        2 denotes that it is PIE sensitized
    Uselines should not be used for binning. A separate functions exists for 
        that.
    Gchan, Rchan and Ychan indicates the detection channels for the 
        corresponding colors. No distinction between P and S 
        polarisation is made. Generally R and Y have the same channels.
    ntacs is the amount of bins used for the tac histogram. Decrease 
        to safe memory usage. Computational efficiency is minimal
        effected."""
    uselines = uselines.astype(np.ubyte)
    Gchan = Gchan.astype(np.ubyte)
    Rchan = Rchan.astype(np.ubyte)
    Ychan = Ychan.astype(np.ubyte)
    root, file = os.path.split(fname)
    NumRecords = ptuHeader_wrap (fname)
    eventN, tac, t, can = ptu_wrap(fname, NumRecords)
    root, file = os.path.split(fname)
    name, _ = os.path.splitext(file)
    header_name = os.path.join(root, b"header", name + b".txt")
    print('number of records is ' + str(NumRecords))

    dimX, dimY, dwelltime, counttime = read_header(header_name)
    imG, imR, imY = genGRYLifetimeWrap(eventN, tac, t, can, dimX, dimY, ntacs, \
        dwelltime, counttime, NumRecords, uselines, Gchan, Rchan, Ychan)
    print('total intensity of imG is %i' % imG.sum())
    return imG, imR, imY



#wdir = b'N:\\Singlem\\singlem19-2\\April\\18_Calibrations_NV\\'
#file = b'10micros\\10micros_1.ptu'
#fname = os.path.join(wdir,file)
#genImABfromptu(fname, uselines = np.ones(2, dtype = np.byte), xbinning = 2, ybinning = 2)

#file = b'20micros\\20micros_1.ptu'
#fname = os.path.join(wdir,file)
#genImABfromptu(fname, uselines = np.ones(1, dtype = np.int))
wdir = b'N:\\Singlem\\singlem19-3\\August\\05_Origami_NV_AB_JHB\\completely labelled'
file = b'Origami_TLMR_6_ROXS.ptu'
fname = os.path.join(wdir, file)
initLifetime(fname, uselines = np.array([1,1,1,2]))