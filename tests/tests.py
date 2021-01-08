# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 12:48:04 2021

@author: voort
"""
#%% test read ptu_wrap
import cpp_wrappers
import numpy as np
import os
testfiledir =  r'C:\Users\voort\Documents\201111_lightdataset_dev\freedye'

def test_imread():
    fname = os.path.join(testfiledir, 'LSM_Rh110_test_pulsed3.ptu').encode()
    uselines = np.array([1]).astype(np.ubyte)
    Gchan = np.array([1,1]).astype(np.ubyte)
    Rchan = np.array([0,0]).astype(np.ubyte)
    Ychan = np.array([0]).astype(np.ubyte)
    ntacs = 1024
    TAC_range = 32768
    dwelltime = 10e-6
    framestop = -1
    root, file = os.path.split(fname)
    NumRecords = cpp_wrappers.ptuHeader_wrap (fname)
    eventN, tac, t, can = cpp_wrappers.ptu_wrap(fname, NumRecords)
    root, file = os.path.split(fname)
    name, _ = os.path.splitext(file)
    header_name = os.path.join(root, b"header", name + b".txt")
    print('number of records is ' + str(NumRecords))

    dimX, dimY, _dwelltime, counttime = cpp_wrappers.read_header(header_name)
    if _dwelltime == -1:
        print('dwelltime not found in header, using user-set value')
    else:
        dwelltime = _dwelltime
    imG, imR, imY = cpp_wrappers.genGRYLifetimeWrap(
        eventN, tac, t, can, dimX, dimY, ntacs, TAC_range, dwelltime, 
        counttime, NumRecords, uselines, Gchan, Rchan, Ychan, framestop
        )
    assert np.sum(imG) != 0
    assert np.sum(imR) != 0
    return imG, imR, imY
imG, imR, imY = test_imread()