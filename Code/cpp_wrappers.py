# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 11:53:04 2019

@author: voort
"""

import numpy as np
import ctypes
import os
wdir = os.path.dirname(__file__)
readPTU = ctypes.WinDLL (os.path.join(wdir, "wrapped", "PQ_PTU.dll"))

#_SplitOnTacs = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\ProcessPhotonStream.dll').SplitOnTacs
debug = False
if debug:
    _SplitOnTacs = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\x64\Debug\ProcessPhotonStream.dll').SplitOnTacs
else:
    wdir = os.path.dirname(__file__)
    _SplitOnTacs = ctypes.CDLL(os.path.join(wdir, "wrapped", 'ProcessPhotonStream.dll')).SplitOnTacs
    _genGRYlifetime = ctypes.CDLL(os.path.join(wdir, "wrapped", 'ProcessPhotonStream.dll')).genGRYlifetime


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
        #print("header dir already exists")
        pass
    outfile =  os.path.join(root, b"header", name + b".txt")
    fpout = ctypes.create_string_buffer(outfile)
    #print(outfile)
    return readPTU.PQ_ptuHeader_sf(fpin,fpout)

def ptu_wrap(fname, NumRecords):
    #bug: when NumRecords is not close to the actual number of photons in the file
        #an error is raised. However it is desirable that a partial dataset is 
        #returned instead.
	#issue: for new imspy reading lists are more convenient and this
	#is more close to ctypes anyway
	#need to write a function that returns list types
    #initialize variables in memory for the c routine to write in
    print('loading file %s' % os.path.split(fname.decode())[1])
    c_longlong_p = ctypes.POINTER(ctypes.c_longlong) #init class for long long pointer
    c_ubyte_p = ctypes.POINTER(ctypes.c_ubyte) #init class for unsigned char pointer
    c_int_p = ctypes.POINTER(ctypes.c_int) #init class for int pointer



    length = ctypes.c_longlong(NumRecords)
    fpin = ctypes.create_string_buffer(fname)

    eventN = np.zeros(NumRecords).astype(np.int64)
    eventN_p = eventN.ctypes.data_as(c_longlong_p)

        #tac can be caught in 16 bit short, but my PQ_ptu_sf function expects an int.
		#potentially this already occurs in Surens code, but maybe also my wrapper is at fault.
		#I will except this inefficiency for now.
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
        dwelltime
    except UnboundLocalError:
        dwelltime = -1
    try:
        return dimX, dimY, dwelltime, counttime
    except NameError:
        logger.error("not all needed variables were found")
        raise

def read_header_v2(header_name):
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
        if  "TTResultFormat_DTimeMask" in el:
            TAC_range = int(el[40:])
    try:
        dwelltime
    except UnboundLocalError:
        dwelltime = -1
    try:
        return dimX, dimY, dwelltime, counttime, TAC_range
    except NameError:
        logger.error("not all needed variables were found")
        raise
###this function was superceded by wrapping using pywrap functionality,
#Producing a .pyd file.
#Left here for now as an example how to deal with LVArray structs.
# def fit2DGaussian_wrap(params0, im, debug = False):
    # """wrapper to fit2DGaussian function from cpp used in Ani Fitting routine
    # params0: initial guess for parameters. They are ordered according to Mortensen et al.
    # im: 2D Gauss image to be fitted
    # a: pixel size"""
    
    # #ISSUE: class declaration occurs inside of function, this disables use outside of current
    # #function. It cannot be changed easily as imsize must be known for class initialisation.
    
    # if debug:
        # dllfile = r"K:\vanderVoortN\FRC\dev\Fit2DGaussian_fresh\build\Debug\Fit2DGaussian.dll"
    # else:
        # dllfile = r"K:\vanderVoortN\FRC\Code\wrapped\Fit2DGaussian.dll"
    # fit2DGaussian = ctypes.WinDLL(dllfile).fit2DGaussian

    # c_double_p = ctypes.POINTER(ctypes.c_double)
    # xshape, yshape = im.shape
    # im = im.flatten()
    # imsize = im.shape[0]
    # DOUBLEARRAY = ctypes.c_double * imsize

    # class LVDoubleArray(ctypes.Structure):
        # _fields_ = [
            # ('length', ctypes.c_int),#length of array i.e. imsize
            # ('data', DOUBLEARRAY) #image array
        # ]

    # class MGPARAM(ctypes.Structure):
        # _fields_ = [
            # ('subimage', ctypes.POINTER(ctypes.POINTER(LVDoubleArray))),
            # ('osize', ctypes.c_int), #number or rows in squared image
            # ('M', ctypes.POINTER(ctypes.POINTER(LVDoubleArray)))
        # ]
    # #give paramaters in correct format
    # params_c = params0.astype(np.double).ctypes.data_as(c_double_p)
    # fit2DGaussian.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(MGPARAM)]
    # fit2DGaussian.restype = ctypes.c_double

    # c_im = DOUBLEARRAY()
    # c_M = DOUBLEARRAY()

    # for i, el in enumerate(im):
        # c_im[i] = el
        # c_M[i] = 1

    # c_imsize = ctypes.c_int(imsize)
    # c_osize = ctypes.c_int(int(np.sqrt(imsize)))


    # subimage = ctypes.pointer(LVDoubleArray(c_imsize, c_im))
    # subM = ctypes.pointer(LVDoubleArray(c_imsize, c_M))
    # mgparam = MGPARAM(ctypes.pointer(subimage), c_osize, ctypes.pointer(subM))
    # Istar = fit2DGaussian(params_c, mgparam)

    # #conversion from ctypes.c_double_p to numpy arrav
    # #there is some bug concerning names that are pointing to the same object.
    # #uncared for TODO
    # paramlength = params0.shape[0]
    # params = []
    # for i in range(paramlength):
        # params.append(params_c[i])
    # params = np.array(params)
    
    # model = []
    # for i in range(im.size):
        # model.append(mgparam.M.data[i])
    # model = np.array(model).reshape(xshape, yshape)

    # del (c_im, c_M, c_imsize, c_osize, subimage, subM, mgparam, params_c)

    # return params, Istar, model

def genGRYLifetimeWrap(eventN, tac, t, can, dimX, dimY, ntacs, TAC_range, 
                        dwelltime, counttime, NumRecords, uselines, 
                        Gchan, Rchan, Ychan, framestop):
    """c code wrapper to create tac histogram image. To be used in conjunction
        with PQ_ptuHeader, Read_header, PQ_ptu_sf_wrapper."""
        #I wish this function would return the number of frames. 
        #this would negate the necessity for the user to set it manually.
    assert len(Gchan) == 2 and len(Rchan) == 2 and len(Ychan) == 2,\
        'in current implementations each Chan must be two ints long,' +\
        'use the same channel number twice to access single channel readout'
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
    C_TAC_range = ctypes.c_int(TAC_range)
    C_dwelltime = ctypes.c_float(dwelltime)
    C_counttime = ctypes.c_float(counttime)
    C_NumRecords = ctypes.c_longlong(NumRecords)
    uselines_p = uselines.ctypes.data_as(c_ubyte_p)
    Gchan_p = Gchan.ctypes.data_as(c_ushort_p)
    Rchan_p = Rchan.ctypes.data_as(c_ushort_p)
    Ychan_p = Ychan.ctypes.data_as(c_ushort_p)
    C_nlines = ctypes.c_int(uselines.shape[0])
    #if there are N frames, there are N+1 frame markers. 
    #the +1 is a hotfix. Should be moved to c code.
    C_framestop = ctypes.c_int(framestop + 1)
    imG = np.zeros(dimX * dimY * ntacs).astype(np.int32)
    imR = np.zeros(dimX * dimY * ntacs).astype(np.int32)
    imY = np.zeros(dimX * dimY * ntacs).astype(np.int32)
    imG_p = imG.ctypes.data_as(c_int_p)
    imR_p = imR.ctypes.data_as(c_int_p)
    imY_p = imY.ctypes.data_as(c_int_p)
    _genGRYlifetime(eventN_p, tac_p, t_p, can_p, C_dimX, C_dimY, C_ntacs, 
                    C_TAC_range, C_dwelltime, C_counttime, C_NumRecords,
                    C_nlines, C_framestop, uselines_p, 
                    Gchan_p, Rchan_p, Ychan_p, imG_p, imR_p, imY_p)
    imG = imG.reshape((dimX, dimY, ntacs))
    imR = imR.reshape((dimX, dimY, ntacs))
    imY = imY.reshape((dimX, dimY, ntacs))
    return imG, imR, imY
#######################Stanislav's Fit Functions################################
# fitdll = ctypes.WinDLL(r"K:\vanderVoortN\FRC\dev\wrapFits\fits2x.dll")

# #DOUBLEARRAY = ctypes.c_double * length
# #INTARRAY = ctypes.c_int * length
# c_short_p = ctypes.POINTER(ctypes.c_short)
# c_int_p = ctypes.POINTER(ctypes.c_int)
# c_double_p = ctypes.POINTER(ctypes.c_double)
# p = ctypes.pointer #shorthand p is a function that returns a c pointer of object
# class LVDoubleArray(ctypes.Structure):
    # _fields_ = [
        # ('length', ctypes.c_int),#length of array i.e. imsize
        # ('data', c_double_p) #image array
    # ]
# class LVI32Array(ctypes.Structure):
    # _fields_ = [
        # ('length', ctypes.c_int),#length of array i.e. imsize
        # ('data', c_int_p) #image array
    # ]
# class MPARAM(ctypes.Structure):
    # _fields_ = [
        # ('expdata', ctypes.POINTER(ctypes.POINTER(LVI32Array))),
        # ('irf', ctypes.POINTER(ctypes.POINTER(LVDoubleArray))),
        # ('bg', ctypes.POINTER(ctypes.POINTER(LVDoubleArray))),
        # ('dt', ctypes.c_double),
        # ('corrections', ctypes.POINTER(ctypes.POINTER(LVDoubleArray))),
        # ('M', ctypes.POINTER(ctypes.POINTER(LVDoubleArray)))
        # ]
# MPARAM_p = ctypes.POINTER(MPARAM)


# def fit23_wrap(x, fixed, expdata, irf, bg, dt):
    # fit23 = fitdll.fit23
    # #maybe we need to do something similar to this.
    # fit23.argtypes = [c_double_p, c_short_p, MPARAM_p]
    # fit23.restype = ctypes.c_double
    # length = expdata.size
    # corrections = np.zeros(expdata.size)
    # M = np.zeros(expdata.size)
    # length_c = ctypes.c_int(length)
    # x_c = x.astype(np.double).ctypes.data_as(c_double_p)
    # fixed_c = fixed.astype(np.short).ctypes.data_as(c_short_p)
    # expdata_c = expdata.astype(np.intc).ctypes.data_as(c_int_p)
    # irf_c = irf.astype(np.double).ctypes.data_as(c_double_p)
    # bg_c = bg.astype(np.double).ctypes.data_as(c_double_p)
    # dt_c = ctypes.c_double(dt)
    # corrections_c = corrections.astype(np.double).ctypes.data_as(c_double_p)
    # M_c = M.astype(np.double).ctypes.data_as(c_double_p)
    # mparam = MPARAM(
        # p(p(LVI32Array(length_c, expdata_c))),
        # p(p(LVDoubleArray(length_c, irf_c))),
        # p(p(LVDoubleArray(length_c, bg_c))),
        # dt_c,
        # p(p(LVDoubleArray(length_c, corrections_c))),
        # p(p(LVDoubleArray(length_c, M_c)))
        # )
    # tIstar = fit23(x_c, fixed_c, p(mparam))
    # return tIstar
    
############################functions rarely used###############################
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
    C_NumRecords = ctypes.c_longlong(NumRecords)
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