# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 12:48:04 2021

@author: voort
"""
#%% test read ptu_wrap
import cpp_wrappers
import numpy as np
import os
try:
    dir_path = os.path.dirname(os.path.realpath(__file__))
    testdatadir = os.path.join(dir_path, 'data') #breaks if executing from cell
except NameError: 
    testdatadir = r'K:\vanderVoortN\FRC\tests\data'

#%%
def test_imread():
    fname = os.path.join(testdatadir, 'freedye','LSM_Rh110_test_pulsed3.ptu').encode()
    uselines = np.array([1]).astype(np.ubyte)
    Gchan = np.array([1,1]).astype(np.ubyte)
    Rchan = np.array([0,0]).astype(np.ubyte)
    Ychan = np.array([0,0]).astype(np.ubyte)
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
    dwelltime = 10e-6 #not in header for LSM data
    imG, imR, imY = cpp_wrappers.genGRYLifetimeWrap(
        eventN, tac, t, can, dimX, dimY, ntacs, TAC_range, dwelltime, 
        counttime, NumRecords, uselines, Gchan, Rchan, Ychan, framestop
        )
    assert np.sum(imG) != 0
    assert np.sum(imR) != 0
    print('imread test successful')
    return imG, imR, imY
imG, imR, imY = test_imread()

#%%
import fitDA
def test_fitD0():
    fname = os.path.join(testdatadir, 'D0TAC', 'Donly_#1516_1_G_PS.dat')
    D0dat = np.genfromtxt(fname)[30:300]
    *_, chi2red = fitDA.fitDonly(D0dat, dtime = 0.064)
    assert chi2red < 1.5
    #write intelligent tests
    print('D0fit test successful')
    return chi2red
test_fitD0()

#%%
import fitDA
def test_fitDA1lt():
    datarange = [30, 300]
    D0name = os.path.join(testdatadir, 'D0TAC', 'Donly_#1516_10_G_PS.dat')
    DAname = os.path.join(testdatadir, 'DATAC', 'D1A10_#1516#1605_1_G_PS.dat')
    D0dat = np.genfromtxt(D0name)[datarange[0]:datarange[1]]
    DAdat = np.genfromtxt(DAname)[datarange[0]:datarange[1]]
    popt,*_, chi2red = fitDA.fitDA1lt(DAdat, D0dat, dtime = 0.064)
    #write intelligent tests
    assert chi2red < 1.5
    return chi2red
test_fitDA1lt()

#%%
import fitDA
def test_fitDA2lt():
    datarange = [30, 450]
    D0name = os.path.join(testdatadir, 'D0TAC', 'Donly_#1516_10_G_PS.dat')
    file = 'D1A10_#1516#1605_1_G_PS.dat'
    DAfolder = os.path.join(testdatadir, 'DATAC')
    DAname = os.path.join(DAfolder, 'D1A10_#1516#1605_1_G_PS.dat')
    D0dat = np.genfromtxt(D0name)[datarange[0]:datarange[1]]
    DAdat = np.genfromtxt(DAname)[datarange[0]:datarange[1]]
    popt, _, DAmodel, chi2red = fitDA.fitDA2lt(DAdat, D0dat, dtime = 0.064)
    print(popt)
    _, _, _, D0model, chi2redD0 = fitDA.fitDonly(D0dat)
    #write intelligent tests
    if chi2red > 1.5:
        fitDA.pltDA_eps(DAdat, D0dat, DAmodel, D0model, file, popt, chi2red,
                        chi2redD0, DAfolder)
    return chi2red
test_fitDA2lt()
