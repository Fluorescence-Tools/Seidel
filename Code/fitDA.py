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

#issue dtime value is shared over all functions. Better to make it a class 
#variable

def Donly(t, A, B, tau1, tau2, bg):
    return A * np.exp(-t/tau1) + B * np.exp(-t/tau2) + bg

def DA1lt(t, A, x0, kfret, bg, D0):
    return A * D0 * ( (1-x0) * np.exp(- t * kfret) + x0 ) + bg
def DA2lt(t, A, x1, x2, kfret1, kfret2, bg, D0):
    """x1 and tau_fr1 denote the fraction and rate for FRET fraction 1
    x2 and tau_fr2 denote the fraction and rate for FRET fraction 2
    The no FRET fraction is calculated from (1 - x1 - x2)"""
    return A * D0 * ( (1-x1 -x2) \
                   + x1 * np.exp(- t * kfret1)  \
                   + x2 * np.exp(- t * kfret1)) \
                   + bg
DA = DA1lt #legacy name
def eps(t, x0, tau_fret):
    return (1-x0) * np.exp(-t / tau_fret) + x0
    
def fitDonly(D0dat, dtime = 0.064):
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    p0 = [np.max(D0dat)/2, np.max(D0dat)/2, 1, 3, 100]
    popt, pcov = curve_fit(Donly, fittime, D0dat, p0 = p0, sigma = np.sqrt(D0dat))
    Donlymodel = Donly(fittime, popt[0], popt[1], popt[2], popt[3], popt[4])
    Donly_base = Donly(fittime, popt[0] / ( popt[0] + popt[1] ), \
                            popt[1]/ ( popt[0] + popt[1]), popt[2], popt[3], 0)
    chi2red = np.sum( (D0dat-Donlymodel)**2 / D0dat) / (Npoints - 5)
    #print('chi2 reduced is %.2f' % chi2red)
    return popt, pcov, Donly_base, Donlymodel, chi2red

def fitDA1lt (DAdat, D0dat, dtime = 0.064):
    _, _, Donly_base, _, _ = fitDonly(D0dat, dtime = dtime)
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    p0 = [np.max(DAdat), 0.8, 10, 50]
    popt, pcov = curve_fit( lambda t, A, x0, kf, bg: DA1lt(t, A, x0, kf, bg, Donly_base), \
                           fittime, DAdat, p0 = p0, sigma = np.sqrt(DAdat),
                           bounds = ([0, 0, 0, 0], [np.inf, 1,1e4, np.inf]))
    DAmodel = DA1lt(fittime, *popt, Donly_base)
    chi2red = np.sum( (DAdat-DAmodel)**2 / DAdat) / (Npoints - 3)
    print('chi2 reduced is %.2f' % chi2red)
    return popt, pcov, DAmodel, chi2red
def fitDA2lt (DAdat, D0dat, dtime = 0.064):
    _, _, Donly_base, _, _ = fitDonly(D0dat, dtime = dtime)
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    p0 = [np.max(DAdat), 0.2, 0.2, 1, 10, 50]
    popt, pcov = curve_fit( lambda t, A, x0, x1, kf1, kf2, bg: \
                                DA2lt(t, A, x0, x1, kf1, kf2, bg, Donly_base), \
                           fittime, DAdat, p0 = p0, sigma = np.sqrt(DAdat),
                           bounds = ([0, 0, 0, 0, 0, 0], 
                                     [np.inf, 1, 1, 1e4, 1e4, np.inf]))
    DAmodel = DA2lt(fittime, *popt, Donly_base)
    chi2red = np.sum( (DAdat-DAmodel)**2 / DAdat) / (Npoints - 3)
    print('chi2 reduced is %.2f' % chi2red)
    return popt, pcov, DAmodel, chi2red
fitDA = fitDA1lt #legacy name
    
def plteps(ax, DAdat, D0dat, x0, tau_fret, bgrange = [320,420], makeplot = True, dtime = 0.064):
    #calc backgrounds
    bgest_DA = np.mean(DAdat[bgrange[0]:bgrange[1]])
    bgest_D0 = np.mean(D0dat[bgrange[0]:bgrange[1]])
    #define time axis
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    #calc epsilon
    epsdat = (DAdat - bgest_DA) / (D0dat- bgest_D0) * \
            max(D0dat-bgest_D0) / max(DAdat - bgest_DA)
    epsmod = eps(fittime, x0, tau_fret)
   
    #plot
    if makeplot:
        ax.plot(fittime, epsdat, label = '\u03B5(t)')
        ax.plot(fittime, epsmod, label = '\u03B5(t) fit')
        ax.set_xlim(0,20)
        ax.set_ylim(0.1,1.1)
        ax.set_xlabel('time(ns)')
        ax.set_ylabel('\u03B5D (t)')
        ax.legend()
        
def pltDA(ax, DAdat, D0dat, DAmodel, Donlymodel, file, popt, chi2red, chi2red_D0, dtime = 0.064):
    #define time axis
    Npoints = D0dat.shape[0]
    fittime = np.arange(Npoints) * dtime
    ax.plot(fittime, DAdat, label = 'D(A)')
    ax.plot(fittime, D0dat, label = 'D(0)')
    ax.plot(fittime, DAmodel, 'r--', label = 'D(A) fit')
    ax.plot(fittime, Donlymodel, 'c--', label = 'D(0) fit')
    ax.set_yscale('log')
    ax.tick_params(direction='in', top=True, right=True)
    ax.tick_params(direction='in', labelbottom=False)

    ax.legend()
    #ax.set_xlabel('time (ns)')
    ax.set_ylabel('cnts')
    ax.set_xlim(0, 20)

    ax.text(0.5,100, 'x0: %.2f\n1/k_fret: %.2f ns \n\u03C72 D(A): %.2f\n\u03C72 D(0): %.2f'\
             % (popt[1], popt[2], chi2red, chi2red_D0), fontsize = 11)
             
def pltDA_eps(DAdat, D0dat, DAmodel, Donlymodel, file, popt, chi2red, chi2red_D0, outdir):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, figsize=(7, 6))
    fig.subplots_adjust(hspace=0)
    ax1 = plt.subplot(2,1,1)
    plt.title('FRET induced donor decay for %s' % file)
    pltDA(ax1, DAdat, D0dat, DAmodel, Donlymodel, file, popt, chi2red, chi2red_D0)
    plteps(ax2, DAdat, D0dat, popt[1], popt[2])
    plt.savefig(os.path.join(outdir,file[:-4]+'.png'), dpi = 300, bbox_inches = 'tight')
    plt.show()
    return 