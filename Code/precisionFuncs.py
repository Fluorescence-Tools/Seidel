# -*- coding: utf-8 -*-
"""

Code adapted from Mortensen et. al.
Nat Methods. 2010 May ; 7(5): 377–381. doi:10.1038/nmeth.1447.

% Maximum likelihood fit to a 2D Gaussian with a constant background
%
%   N* 1/(2*pi*s^2) * exp (-( (x-ux).^2+(y-uy).^2 ) / (2*s^2)) + b^2
%
% Input Parameters:
% Required:
% data -- the image of the isolated probe.
% params0 -- The user's initial guess of the parameters to be fit:
%           [ux, uy, s, b, N]
%           ux, uy and s should be specified in nanometers.
%           note: as the fit minimizes only localy it is important for this
%           inital guess to be fairly close to the true value.
            note 2: this fit optimizes the adjusted s as per eq. 6 suppl. To obtain the FWHM of the PSF,
            it must be corrected for the pixelsize.
            ux, uy: center position in nm
            s: cilindrical sigma in nm
            b: b^2 background density (cnts / px)
            N: Total number of signal photons in all pixels
% a -- pixel size in nanometers
%
% Optional:
% plot_on -- a binary parameter that determines whether the outcome is
% plotted.
% find_variance -- a binary parameter that determines whether the variance
% should be calculated. If find_variance is 0, then -1 is returned as the
% variance.
% noise_floor -- a constant offset added by the camera before the pixel is
% read out. This is necessary if an accurate spot amplitude is desired.
%
% Output parameter:
% paramsF -- the result of the fit.
%
% Optional:
% variance -- the variance of the localization based on eq.(5) of Mortensen
% et al.
Added feature: ability to make a single bead image.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad

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

def findVar(params, a, verbose = False):
    """find variance according to eq.54 suppl. Mortensen et. al."""
    b2 = params[3]**2
    N = params[4]
    sa = np.sqrt(params[2]**2 + a**2 / 12)
    F = lambda t: np.log(t) / (1 + (N * a**2 * t / (2 * np.pi * sa**2 * b2)))
    #quad function does not take array, makes array calculations slow
    integral, _ = quad(F, 0, 1)
    if verbose:
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
