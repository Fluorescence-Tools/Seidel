import ctypes
from ctypes import *
import numpy as np

def twoDGauss(x, y, ux, uy, s):
    """2D circular Gauss function with integral normalised to one. """
    p = 1 / (2 * np.pi * s**2) * np.exp(- ( (x-ux)**2 + (y-uy)**2 )/ (2*s**2))
    return p
    
def Expected(x, y, params, a):
    """calculate the expectation value in point x,y given the gaussian parameters params and pixel size a"""
    #an alteration to the original MatLab code has been made to account for b being defined as
    #background density / nm. This was missing in the original code.
    E = params[4] * a**2 * twoDGauss(x, y, params[0], params[1], params[2]) + params[3]**2
    return E
    
def createGaussImg(dim, params, a):
    """Relies on the function Expected to create an expectation value image. Then takes Poisson statistics."""
    x = (np.arange(dim)+0.5) * a
    y = x
    x, y = np.meshgrid(x,y)
    return np.random.poisson(Expected(x, y, params, a))

def findVar(params, a):
    """find variance according to eq.54 suppl. Mortensen et. al."""
    b2 = params[3]**2
    N = params[4]
    sa = params[2]
    F = lambda t: np.log(t) / (1 + (N * a**2 * t / (2 * np.pi * sa**2 * b2)))
    integral, _ = quad(F, 0, 1)
    print('integral in variance has value %f' % integral)
    return sa**2 / N  /(1+integral)
    
    
def logLikelihood(params0, a, im):
    """the - log of the Fisher Likelihood function for a Gaussian with Poisson statistics.
    see eq. 49 of suppl Mortensen."""
    x, y = im.shape
    x = (np.arange(x)+0.5)*a
    y = (np.arange(y)+0.5)*a
    x, y = np.meshgrid(x, y)
    datafun = lambda params: Expected(x, y, params, a).sum(axis = (0,1)) - (
        im * np.log(Expected(x, y, params, a))).sum(axis = (0,1))
    optimizedResult = minimize(datafun, params0, method = 'Nelder-Mead', 
                           options = {'maxiter':10000,'maxfev':10000, 'fatol':1e-5})
    return optimizedResult.x



params0 = [500,500,100,0,1000]
a = 50
im = createGaussImg(20, params0, a).flatten().astype(dtype = np.double)
#im = np.arange(400)
#_____above is input of function


fit2DGaussian = WinDLL(r"K:\vanderVoortN\FRC\dev\Fit2DGaussian\x64\Release\Fit2DGaussian.dll").fit2DGaussian
#fit2DGaussian = WinDLL(r"K:\vanderVoortN\FRC\dev\sandbox\Fit2DGaussian.dll")._Z13fit2DGaussianPdP7MGParam


c_double_p = POINTER(c_double)
imsize = im.shape[0]
DOUBLEARRAY = c_double * imsize

class LVDoubleArray(ctypes.Structure):
    _fields_ = [
        ('length', c_int),#number of channels in image, i.e. 1
        ('data', DOUBLEARRAY) #image array
    ]

class MGPARAM(Structure):
    _fields_ = [
        ('subimage', POINTER(POINTER(LVDoubleArray))),
        ('osize', c_int), #number or rows in squared image
        ('M', POINTER(POINTER(LVDoubleArray)))
    ]
    
variables = np.array([10,10,1000,2,2,0,0,0,0], dtype = np.double).ctypes.data_as(c_double_p)
    
c_double_p = POINTER(c_double)
fit2DGaussian.argtypes = [POINTER(c_double), POINTER(MGPARAM)]
fit2DGaussian.restypes = c_double

M = np.zeros(imsize)
c_im = DOUBLEARRAY()
c_M = DOUBLEARRAY()

for i, el in enumerate(im):
    c_im[i] = el
    c_M[i] = 1

c_imsize = c_int(imsize)
c_osize = c_int(int(np.sqrt(imsize)))


subimage = pointer(LVDoubleArray(c_imsize, c_im))
subM = pointer(LVDoubleArray(c_imsize, c_M))
mgparam = MGPARAM(pointer(subimage), c_osize, pointer(subM))

fit2DGaussian(variables, mgparam)

for i in range(400):
    print("Python reads %f at pos %i" %(mgparam.M[0][0].data[i], i))
for i in range(8):
    print('value of var %i is %f' %(i, variables[i]))