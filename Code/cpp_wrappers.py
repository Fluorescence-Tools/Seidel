# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 11:53:04 2019

@author: voort
"""

import numpy as np
import ctypes
import os

readPTU = ctypes.WinDLL (r"S:\64bit dll 's\PQ_PTU_sf\release\PQ_PTU.dll")
#_SplitOnTacs = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\ProcessPhotonStream.dll').SplitOnTacs
debug = False
if debug:
    _SplitOnTacs = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\x64\Debug\ProcessPhotonStream.dll').SplitOnTacs
else:
    _SplitOnTacs = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\x64\Release\ProcessPhotonStream.dll').SplitOnTacs
    _genGRYlifetime = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\x64\Release\ProcessPhotonStream.dll').genGRYlifetime


def ptuHeader_wrap (fname):
    """wrapper for PQ_ptuHeader_sf function. 
    Returns number of entries and writes  the header in a subdirectory named header."""
    fpin = ctypes.create_string_buffer(fname)
    #strip file name for header location
    root, file = os.path.split(fname)
    name, _ = os.path.splitext(file)
    
    try:
        os.mkdir(os.path.join(root, b"header"))
        print(b'creating new directory:' + os.path.join(root, b"header"))
    except:
        print("header dir already exists")
    outfile =  os.path.join(root, b"header", name + b".txt")
    fpout = ctypes.create_string_buffer(outfile)
    print(outfile)
    return readPTU.PQ_ptuHeader_sf(fpin,fpout)

def ptu_wrap(fname, NumRecords):
    #initialize variables in memory for the c routine to write in
    c_longlong_p = ctypes.POINTER(ctypes.c_longlong) #init class for long long pointer
    c_ubyte_p = ctypes.POINTER(ctypes.c_ubyte) #init class for unsigned char pointer
    c_int_p = ctypes.POINTER(ctypes.c_int) #init class for int pointer



    length = ctypes.c_longlong(NumRecords)
    fpin = ctypes.create_string_buffer(fname)

    eventN = np.zeros(NumRecords).astype(np.int64)
    eventN_p = eventN.ctypes.data_as(c_longlong_p)

    tac = np.zeros(NumRecords).astype(np.int32)
    tac_p = tac.ctypes.data_as(c_int_p)

    t = np.zeros(NumRecords).astype(np.int64)
    t_p = t.ctypes.data_as(c_longlong_p)

    can = np.zeros(NumRecords).astype(np.uint8)
    can_p = can.ctypes.data_as(c_ubyte_p)

    j = ctypes.c_longlong()
    ov_in = ctypes.c_longlong()
    stage = ctypes.c_int()

    readPTU.PQ_ptu_sf(fpin, length, eventN_p, tac_p, t_p, can_p, ctypes.byref(j), ctypes.byref(ov_in), ctypes.byref(stage))
    
    return eventN, tac, t, can

import logging
logger = logging.getLogger('readptu')
def read_header(header_name):
    header = np.genfromtxt(header_name.decode(), delimiter = '\n', dtype = str)
    for el in header:
        if "ImgHdr_PixX" in el:
            dimX = int(el[40:])
        if "ImgHdr_PixY" in el:
            dimY = int(el[40:])
        if "ImgHdr_TimePerPixel" in el:
            dwelltime = float(el[40:]) * 1e-3
        if "MeasDesc_GlobalResolution" in el:
            counttime = float(el[40:])
    try:
        return dimX, dimY, dwelltime, counttime
    except NameError:
        logger.error("not all needed variables were found")
        raise

def SplitOnTacs_wrap(eventN, tac, t, can, dimX, dimY, dwelltime, counttime, NumRecords, uselines, gate = 3328):
    c_longlong_p = ctypes.POINTER(ctypes.c_longlong) #init class for long long pointer
    c_ubyte_p = ctypes.POINTER(ctypes.c_ubyte) #init class for unsigned char pointer
    c_int_p = ctypes.POINTER(ctypes.c_int) #init class for int pointer
    
    eventN_p = eventN.ctypes.data_as(c_longlong_p)
    tac_p = tac.ctypes.data_as(c_int_p)
    t_p = t.ctypes.data_as(c_longlong_p)
    can_p = can.ctypes.data_as(c_ubyte_p)
    C_dimX = ctypes.c_int(dimX)
    C_dimY = ctypes.c_int(dimY)
    C_dwelltime = ctypes.c_float(dwelltime)
    C_counttime = ctypes.c_float(counttime)
    C_NumRecords = ctypes.c_int(NumRecords)
    uselines_p = uselines.ctypes.data_as(c_ubyte_p)
    C_gate = ctypes.c_int(gate)
    C_nlines = ctypes.c_int(uselines.shape[0])
    imA = np.zeros(dimX * dimY).astype(np.int)
    imB = np.zeros(dimX * dimY).astype(np.int)
    imA_p = imA.ctypes.data_as(c_int_p)
    imB_p = imB.ctypes.data_as(c_int_p)
    #initialize empty arrays and handles for imA, imB
    _SplitOnTacs(eventN_p, tac_p, t_p, can_p, C_dimX, C_dimY, C_dwelltime, C_counttime, C_NumRecords, C_gate, 
                 C_nlines, uselines_p, imA_p, imB_p)
    imA = imA.reshape((dimX, dimY))
    imB = imB.reshape((dimX, dimY))
    return imA, imB

def fit2DGaussian_wrap(params0, a, im, debug = False):
    """wrapper to fit2DGaussian function from cpp used in Ani Fitting routine
    params0: initial guess for parameters. They are ordered according to Mortensen et al.
    im: 2D Gauss image to be fitted
    a: pixel size"""
    
    #ISSUE: class declaration occurs inside of function, this disables use outside of current
    #function. It cannot be changed easily as imsize must be known for class initialisation.
    
    if debug:
        fit2DGaussian = ctypes.WinDLL(r"K:\vanderVoortN\FRC\dev\Fit2DGaussian\x64\Debug\Fit2DGaussian.dll").fit2DGaussian
    else:
        fit2DGaussian = ctypes.WinDLL(r"K:\vanderVoortN\FRC\Code\Fit2DGaussian.dll").fit2DGaussian

    c_double_p = ctypes.POINTER(ctypes.c_double)
    im = im.flatten()
    imsize = im.shape[0]
    DOUBLEARRAY = ctypes.c_double * imsize

    class LVDoubleArray(ctypes.Structure):
        _fields_ = [
            ('length', ctypes.c_int),#length of array i.e. imsize
            ('data', DOUBLEARRAY) #image array
        ]

    class MGPARAM(ctypes.Structure):
        _fields_ = [
            ('subimage', ctypes.POINTER(ctypes.POINTER(LVDoubleArray))),
            ('osize', ctypes.c_int), #number or rows in squared image
            ('M', ctypes.POINTER(ctypes.POINTER(LVDoubleArray)))
        ]
    #give paramaters in correct format
    variables = np.zeros(17, dtype = np.double).ctypes.data_as(c_double_p)
    variables[0] = params0[0] / a#ux
    variables[1] = params0[1] / a#uy
    variables[3] = params0[2] / a#sx = s
    variables[2] = params0[4] /(variables[3]**2 * 2*np.pi)#N
    variables[4] = params0[5] #ellipticity, for circular 1
    variables[5] = params0[3]**2#b
    variables[6] = params0[6]
    variables[7] = params0[7]
    variables[8] = params0[8]
    variables[9] = params0[9]
    variables[10] = params0[10]
    variables[11] = params0[11]
    variables[12] = 0#info from optimization algorithm
    variables[13] = 1#bool fr weight or no weight: Ask Suren for correct setting: debree, to be deleted
    variables[14] = params[14]#bool if background is fitted: 0 -> background is fitted
    variables[15] = params[15]#bool for ellipticity: 1 means circular, 0 means elliptical
    variables[16] = params[16]# var to select model. 0: 1x 2DGauss, 1: 2x 2DGauss, 2: 3x 2DGauss
    
#    for i in range(9):
#        print('variables element %i has value %f'%(i, variables[i]))

    c_double_p = ctypes.POINTER(ctypes.c_double) #looks like a double declaration. Remove?
    fit2DGaussian.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(MGPARAM)]
    fit2DGaussian.restypes = ctypes.c_double

    c_im = DOUBLEARRAY()
    c_M = DOUBLEARRAY()

    for i, el in enumerate(im):
        c_im[i] = el
        c_M[i] = 1

    c_imsize = ctypes.c_int(imsize)
    c_osize = ctypes.c_int(int(np.sqrt(imsize)))


    subimage = ctypes.pointer(LVDoubleArray(c_imsize, c_im))
    subM = ctypes.pointer(LVDoubleArray(c_imsize, c_M))
    mgparam = MGPARAM(ctypes.pointer(subimage), c_osize, ctypes.pointer(subM))
    fit2DGaussian(variables, mgparam)
    
    params = np.zeros(12)
    params[0] = (variables[0]+0.5) * a #x0 in nm
    params[1] = (variables[1]+0.5) *a #x1 in nm
    params[2] = variables[3] * a #sigma in nm
    params[3] = np.sqrt(variables[5]) #sqrt of bg
    #variables[2] reports fitted amplitude, not total intensity
    #change by multiplying with 2*pi*sigma**2 ?
    params[4] = variables[2] #A0
    params[5] = variables[4] #ellipticity
    params[6] = (variables[6] + 0.5) * a #x1 in nm
    params[7] = (variables[7] + 0.5) * a #y1 in nm
    params[8] = variables[8] #A1
    params[9] = (variables[9] + 0.5) * a # x2 in nm
    params[10] = (variables[10] + 0.5) * a #y2 in nm
    params[11] = variables[11] #A2
    params[12] = variables[12] #info from optimisation algorithm
    
    return params

def genGRYLifetimeWrap(eventN, tac, t, can, dimX, dimY, ntacs, dwelltime, counttime, NumRecords, 
                       uselines, Gchan, Rchan, Ychan):
    """c code wrapper to create tac histogram image. To be used in conjunction
        with PQ_ptuHeader, Read_header, PQ_ptu_sf_wrapper."""
    c_longlong_p = ctypes.POINTER(ctypes.c_longlong) #init class for long long pointer
    c_ubyte_p = ctypes.POINTER(ctypes.c_ubyte) #init class for unsigned char pointer
    c_ushort_p = ctypes.POINTER(ctypes.c_ushort)
    c_int_p = ctypes.POINTER(ctypes.c_int) #init class for int pointer
    
    eventN_p = eventN.ctypes.data_as(c_longlong_p)
    tac_p = tac.ctypes.data_as(c_int_p)
    t_p = t.ctypes.data_as(c_longlong_p)
    can_p = can.ctypes.data_as(c_ubyte_p)
    C_dimX = ctypes.c_int(dimX)
    C_dimY = ctypes.c_int(dimY)
    C_ntacs = ctypes.c_int(ntacs)
    C_dwelltime = ctypes.c_float(dwelltime)
    C_counttime = ctypes.c_float(counttime)
    C_NumRecords = ctypes.c_int(NumRecords)
    uselines_p = uselines.ctypes.data_as(c_ubyte_p)
    Gchan_p = Gchan.ctypes.data_as(c_ushort_p)
    Rchan_p = Rchan.ctypes.data_as(c_ushort_p)
    Ychan_p = Ychan.ctypes.data_as(c_ushort_p)
    C_nlines = ctypes.c_int(uselines.shape[0])
    imG = np.zeros(dimX * dimY * ntacs).astype(np.int)
    imR = np.zeros(dimX * dimY * ntacs).astype(np.int)
    imY = np.zeros(dimX * dimY * ntacs).astype(np.int)
    imG_p = imG.ctypes.data_as(c_int_p)
    imR_p = imR.ctypes.data_as(c_int_p)
    imY_p = imY.ctypes.data_as(c_int_p)
    _genGRYlifetime(eventN_p, tac_p, t_p, can_p, C_dimX, C_dimY, C_ntacs, C_dwelltime, C_counttime, C_NumRecords,
                   C_nlines, uselines_p, Gchan_p, Rchan_p, Ychan_p, imG_p, imR_p, imY_p)
    imG = imG.reshape((dimX, dimY, ntacs))
    imR = imR.reshape((dimX, dimY, ntacs))
    imY = imY.reshape((dimX, dimY, ntacs))
    return imG, imR, imY