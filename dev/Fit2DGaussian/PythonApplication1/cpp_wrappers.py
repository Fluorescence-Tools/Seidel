import numpy as np
import ctypes
import os

def fit2DGaussian_wrap(params0, a, im):
    """wrapper to fit2DGaussian function from cpp used in Ani Fitting routine
    params0: initial guess for parameters. They are ordered according to Mortensen et al.
    im: 2D Gauss image to be fitted
    a: pixel size"""
    
    #ISSUE: class declaration occurs inside of function, this disables use outside of current
    #function. It cannot be changed easily as imsize must be known for class initialisation.
    
    fit2DGaussian = ctypes.WinDLL(r"K:\vanderVoortN\FRC\dev\Fit2DGaussian\x64\Debug\Fit2DGaussian.dll").fit2DGaussian

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
    variables = np.zeros(9, dtype = np.double).ctypes.data_as(c_double_p)
    variables[0] = params0[0] / a#ux
    variables[1] = params0[1] / a#uy
    variables[2] = params0[4]#N
    variables[3] = params0[2] / a#sx = s
    variables[4] = params0[2] / a#sy = s
    variables[5] = params0[3]**2#b
    variables[6] = 0#info from optimization algorithm
    variables[7] = 0#bool fr weight or no weight: Ask Suren for correct setting
    variables[8] = 0#bool if background is fitted: 0 -> background is fitted

    c_double_p = ctypes.POINTER(ctypes.c_double)
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
    
    params = np.zeros(5)
    params[0] = (variables[0]) * a
    params[1] = (variables[1]) *a
    params[2] = np.sqrt(variables[3]**2+variables[4]**2) * a
    params[3] = variables[5]
    params[4] = variables[2]
    return params
