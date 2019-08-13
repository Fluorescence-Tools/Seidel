
import ctypes
import numpy as np
import os
from ptu_wrapper import *
import matplotlib.pyplot as plt
from PIL import Image


def genImABfromptu(fname, uselines = np.ones(1), gate = 3228):
    NumRecords = ptuHeader_wrap (fname)
    eventN, tac, t, can = ptu_wrap(fname, NumRecords)
    root, file = os.path.split(fname)
    name, _ = os.path.splitext(file)
    header_name = os.path.join(root, b"header", name + b".txt")
    print('number of records is ' + str(NumRecords))

    dimX, dimY, dwelltime, counttime = read_header(header_name)
    #uselines = np.ones(1).astype(np.ubyte)
    imA, imB = SplitOnTacs_wrap(eventN, tac, t, can, dimX, dimY, dwelltime, counttime, 
                                NumRecords, gate = gate, uselines = uselines)
    print("total image intensity is "  + str(np.sum(imA+imB)))

    im = Image.fromarray(imA)
    outname = str(os.path.join(wdir, file[:-4] + b'_imA.tiff'), encoding='utf-8')
    im.save(outname)
    im = Image.fromarray(imB)
    outname = str(os.path.join(wdir, file[:-4] + b'_imB.tiff'), encoding='utf-8')
    im.save(outname)
    

wdir = b'N:\\Singlem\\singlem19-1\\March\\29_20nmCrimson_NV'
file = b'2019032920nmCrimsonExc640_4-3STED775_100_1-25W.ptu'
fname = os.path.join(wdir,file)
genImABfromptu(fname, uselines = np.ones(1, dtype = np.int))
#wdir = b"N:\\Singlem\\singlem19-1\\March\\06_FRC_20nmCrimson_NV\\STEDPowerSeries\\"
#file = b"Crimson20nm_Exc_640_1perc_STED010perc_0016AU.ptu"
#fname = os.path.join(wdir,file)
#genImABfromptu(fname, uselines = np.ones(10, dtype = np.int))