import numpy as np
import ctypes
import os


readPTU = ctypes.WinDLL (r"S:\64bit dll 's\PQ_PTU_sf\release\PQ_PTU.dll")
#_SplitOnTacs = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\ProcessPhotonStream.dll').SplitOnTacs
debug = True
if debug:
    _SplitOnTacs = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\x64\Debug\ProcessPhotonStream.dll').SplitOnTacs
    _genGRYlifetime = ctypes.CDLL(r'K:\vanderVoortN\FRC\dev\readPTU\x64\Debug\ProcessPhotonStream.dll').genGRYlifetime
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
#    header = np.genfromtxt(r"header\\Crimson20nm_Exc_640_1perc_STED100perc_0016AU.txt", delimiter = '\n', dtype = str)
    print(header_name.decode())
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
    _SplitOnTacs(eventN_p, tac_p, t_p, can_p, C_dimX, C_dimY, C_dwelltime, C_counttime, C_NumRecords, C_gate, \
                 C_nlines, uselines_p, imA_p, imB_p)
    imA = imA.reshape((dimX, dimY))
    imB = imB.reshape((dimX, dimY))
    return imA, imB

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