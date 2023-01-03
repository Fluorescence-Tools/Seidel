#library to do simple Donor only calibrated FRET sensitized lifetime fits
#this library should be superceded by automated scripting of the Chisurf programm
#Chisurf has much more extensive and accurate fitting tools.
#author: Nicolaas van der Voort
#AG Seidel
#April 14, 2020

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import lmfit
import warnings

#issue dtime value is shared over all functions. Better to make it a class 
#variable

def CheckFitData(func):
    """decorator to catch common error when bad data causes a fit to fail"""
    def inner(data, *args, **kwargs):
        assert type(data) == np.ndarray
        #I suspect this condition does not work properly. Why?
        if (data == 0).any(): 
            data[data==0] = 1
            warnings.warn('data contains zero values, ' + \
                          'gaussian error estimate cannot handle this, '+\
                          'setting 0 values to 1')
        try: return func(data, *args, **kwargs)
        except RuntimeError:
            plt.plot(data)
            plt.title('data used for fitting')
            plt.yscale('log')
            plt.show()
            raise RuntimeError \
                ('A runtimeError Occurred, did you attempt fit with bad data?')
    return inner

def Donly(t, A, B, tau1, tau2, bg):
    return A * np.exp(-t/tau1) + B * np.exp(-t/tau2) + bg
    
#def Donly_p(t, p):
#    return(Donly(t, p['A'], p['B'], p['tau1'], p['tau2'], p['bg']))

def DA1lt(t, A, x0, kfret, bg, D0):
    return A * D0 * ( (1-x0) * np.exp(- t * kfret) + x0 ) + bg
def DA2lt(t, A, x1, x2, kfret1, kfret2, bg, D0):
    """x1 and tau_fr1 denote the fraction and rate for FRET fraction 1
    x2 and tau_fr2 denote the fraction and rate for FRET fraction 2
    The no FRET fraction is calculated from (1 - x1 - x2)"""
    return A * D0 * ( (1-x1 -x2) \
                   + x1 * np.exp(- t * kfret1)  \
                   + x2 * np.exp(- t * kfret2)) \
                   + bg
DA = DA1lt #legacy name
    
@CheckFitData
def fitDonly(D0dat, dtime = 0.064):
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    p0 = [np.max(D0dat)/2, np.max(D0dat)/2, 1, 3, 100]
    popt, pcov = curve_fit(Donly, fittime, D0dat, p0 = p0, sigma = np.sqrt(D0dat))
    Donlymodel = Donly(fittime, popt[0], popt[1], popt[2], popt[3], popt[4])
    Donly_base = Donly(fittime, popt[0] / ( popt[0] + popt[1] ), \
                            popt[1]/ ( popt[0] + popt[1]), popt[2], popt[3], 0)
    #print(popt)
    chi2red = np.sum( (D0dat-Donlymodel)**2 / Donlymodel) / (Npoints - 5)
    #print('chi2 reduced is %.2f' % chi2red)
    return popt, pcov, Donly_base, Donlymodel, chi2red

@CheckFitData
def fitDA1lt (DAdat, D0dat, dtime = 0.064):
    _, _, Donly_base, _, _ = fitDonly(D0dat, dtime = dtime)
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    p0 = [np.max(DAdat), 0.8, 10, 1000]
    #zero values in data crashes the fitting routine
    for dat in DAdat, D0dat:
        dat[dat==0]=1
    popt, pcov = curve_fit( lambda t, A, x0, kf, bg: DA1lt(t, A, x0, kf, bg, Donly_base), \
                           fittime, DAdat, p0 = p0, sigma = np.sqrt(DAdat),
                           bounds = ([0, 0, 0, 0], [np.inf, 1, 1e4, np.inf]))
    DAmodel = DA1lt(fittime, *popt, Donly_base)
    chi2red = np.sum( (DAdat-DAmodel)**2 / DAmodel) / (Npoints - 3)
    print('chi2 reduced is %.2f' % chi2red)
    return popt, pcov, DAmodel, chi2red

@CheckFitData
def fitDA2lt (DAdat, D0dat, dtime = 0.064):
    _, _, Donly_base, _, _ = fitDonly(D0dat, dtime = dtime)
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    p0 = [np.max(DAdat), 0.2, 0.2, 0.1, 0.2, 0]
    popt, pcov = curve_fit( lambda t, A, x0, x1, kf1, kf2, bg: \
                                DA2lt(t, A, x0, x1, kf1, kf2, bg, Donly_base), \
                           fittime, DAdat, p0 = p0, sigma = np.sqrt(DAdat),
                           bounds = ([0, 0, 0, 0, 0, 0], 
                                     [np.inf, 1, 1, 1e4, 1e4, np.inf]),
                                     method = 'dogbox')
    DAmodel = DA2lt(fittime, *popt, Donly_base)
    chi2red = np.sum( (DAdat-DAmodel)**2 / DAmodel) / (Npoints - len(p0))
    print(chi2red)
    return popt, pcov, DAmodel, chi2red
fitDA = fitDA1lt #legacy name
    

#def get_chi2red(params, func, xdata, ydata, sign = 1, 
#    modelargs = (), modelkwargs = {}):
#    """return chi2red for 1D binned data functions.
#    minimize this function to obtain most likely fit"""
#    #want to modify function such that it just takes modelkwargs
#    #and not xdata, ydata, then accordingly fitfunctions and surface functions
#    #can be generalised
#    model = func(xdata, params, *modelargs, **modelkwargs)
#    nparams = len(xdata) - len(params)
#    chi2red = np.sum((model - ydata)**2 / ydata) / nparams
#    print(chi2red)
#    return chi2red * sign


def plteps(ax, DAdat, D0dat, DAmodel, Donlymodel, 
            makeplot = True, dtime = 0.064, **kwargs):
    #define time axis
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    D0bg = np.mean(Donlymodel[-5:])
    DAbg = np.mean(DAmodel[-5:])
    #calc epsilon
    norm = max(Donlymodel - D0bg) / max(DAmodel - DAbg)
    epsdat = (DAdat - DAbg) / (Donlymodel - D0bg) * norm
    epsmod = (DAmodel - DAbg) / (Donlymodel - D0bg) * norm
   
    #plot
    if makeplot:
        ax.plot(fittime, epsdat, label = '\u03B5(t)')
        ax.plot(fittime, epsmod, label = '\u03B5(t) fit')
        ax.grid()
        ax.set_xlim(0,20)
        ax.set_ylim(0.1,1.1)
        ax.set_xlabel('time(ns)')
        ax.set_ylabel('\u03B5D (t)')
        ax.legend()
        
def pltD0(D0dat, Donlymodel, name, outdir, dtime = 0.064):
    #define time pltis
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    plt.plot(fittime, D0dat, label = 'D(0)')
    plt.plot(fittime, Donlymodel, 'c--', label = 'D(0) fit')
    plt.yscale('log')
    plt.legend()
    plt.xlabel('time (ns)')
    plt.ylabel('cnts')
    plt.xlim(0, 20)
    plt.savefig(os.path.join(outdir,name[:-4]+'.png'), dpi = 300, bbox_inches = 'tight')
    plt.show()
        
def pltDA(ax, DAdat, D0dat, DAmodel, Donlymodel, popt, chi2red, chi2red_D0, dtime = 0.064):
    #define time axis
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    DAmax = max(DAmodel)
    D0max = max(Donlymodel)
    ax.plot(fittime, DAdat / DAmax, label = 'D(A)')
    ax.plot(fittime, D0dat / D0max, label = 'D(0)')
    ax.plot(fittime, DAmodel / DAmax, 'r--', label = 'D(A) fit')
    ax.plot(fittime, Donlymodel / D0max, 'c--', label = 'D(0) fit')
    ax.set_yscale('log')
    ax.tick_params(direction='in', top=True, right=True)
    ax.tick_params(direction='in', labelbottom=False)

    ax.legend()
    #ax.set_xlabel('time (ns)')
    ax.set_ylabel('cnts')
    ax.set_xlim(0, 20)

    ax.text(0.5,0.01, 'x0: %.2f\nk_fret: %.2f (1/ns) \n\u03C72 D(A): %.2f\n\u03C72 D(0): %.2f'\
             % (popt[1], popt[2], chi2red, chi2red_D0), fontsize = 11)
             
def pltDA_eps(DAdat, D0dat, DAmodel, Donlymodel, name, popt, chi2red, chi2red_D0, outdir, **kwargs):
    #consider making some args kwargs.
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, figsize=(7, 6))
    fig.subplots_adjust(hspace=0)
    ax1 = plt.subplot(2,1,1)
    plt.title('FRET induced donor decay for %s' % name)
    pltDA(ax1, DAdat, D0dat, DAmodel, Donlymodel, popt, chi2red, chi2red_D0)
    plteps(ax2, DAdat, D0dat, DAmodel, Donlymodel, **kwargs)
    plt.savefig(os.path.join(outdir,name+'.png'), dpi = 300, bbox_inches = 'tight')
    plt.show()
    return 