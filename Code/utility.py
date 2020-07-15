#A collection of utility functions
#functions have some common use, but are not intended to be maintained
#author: Nicolaas van der Voort
#AG Seidel, HHU
#14 April 2020

import ImageManipulation as IM
import numpy as np
import os
import fitDA
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import lmfit
import developmental_functions as df

def exportLSMTAC(fname, outdir, dwelltime, pulsetime, uselines = np.array([1]), 
                 Gchan = np.array([0,1]), Rchan = np.array([4,5]), Ychan = np.array([4,5]),
                ntacs = 1024, TAC_range = 4096, PIE_gate = 440):
    """utility function to export GRY (P+S) TAC histogram for LSM microscope"""
    #issue: when Rchan and Ychan are identical (i.e. true PIE), then Y channel
    #remains empty. Workaround implemented. ugly
    _, file = os.path.split(fname)
    data = IM.processLifetimeImage(fname.encode(), uselines = uselines, Gchan = Gchan, 
                                   Rchan = Rchan, 
                                   Ychan = Ychan, ntacs = ntacs, TAC_range = TAC_range, \
                                   pulsetime = pulsetime, dwelltime = dwelltime)
    data.loadLifetime()
    if (Rchan == Ychan).all(): # workaround in case of PIE
       data.workLifetime.Y = copy.deepcopy(data.workLifetime.R)
    data.gate(0,PIE_gate, channel = 'R')
    data.gate(PIE_gate, -1, channel = 'Y')
    GTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'G'))).transpose()
    RTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'R'))).transpose()
    YTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'Y'))).transpose()
    np.savetxt(os.path.join(outdir, file[:-4] + '_G.dat'), GTACS, fmt = '%i')
    np.savetxt(os.path.join(outdir, file[:-4] + '_R.dat'), RTACS, fmt = '%i')
    np.savetxt(os.path.join(outdir, file[:-4] + '_Y.dat'), YTACS, fmt = '%i')
    
def genGYfnames(wdir):
    """utility function to get all matching G and Y TAC decays in a list
    as exported from exportLSMTAC"""
    Gfiles = []
    Yfiles = []
    files = os.listdir(wdir)
    for file in files:
        if file[-6:] == '_G.dat' and file[:5] != 'Donly':
            Gfiles.append(os.path.join(wdir, file))
            #check for Y files that have the same identifier
            Yfile = file[:-6]+'_Y.dat'
            if Yfile in files:
                Yfiles.append(os.path.join(wdir, Yfile))
            else:
                raise FileNotFoundError('not all matching yellow files found')
    return Gfiles, Yfiles
def batchAnalyseDA(Gfiles, Yfiles, D0_fname, plotoutDir, paramout, dtime = 0.064, times = []):
    """utility function to batch-fit a set of Donor only calibrated
    Donor Acceptor fits. 
    function assumes files were exported using exportLSMTAC naming convention"""
    #issue: passing wdir and files is not really clean
    #better would be to pass list of G files with matching list of Yellow files
    #have left it like to so for the time being.
    #settings
    mpl.rcParams.update({'font.size': 16})

    #fit Donly
    D0dat = np.genfromtxt(D0_fname)[30:600,1]
    popt_D0,_,Donly_base, Donlymodel, chi2red_D0 = fitDA.fitDonly(D0dat)

    names = []
    S_y = []
    x0 = []
    tau_fret = []
    chi2red_lst = []
    if times == []:
        times = [0] * len(Gfiles)

    for Gfile, Yfile in zip(Gfiles, Yfiles):
        names.append(Gfile[:-6])
        #get brightness of Yellow signal
        S_y.append(np.sum(np.genfromtxt(Yfile)[30:600,1]))
        #fit D(A)
        print('analyzing file %s' % Gfile)
        DAdat = np.genfromtxt(Gfile)[30:600,1]
        popt, pcov, DAmodel, chi2red = fitDA.fitDA(DAdat, D0dat, dtime)
        x0.append(popt[1])
        tau_fret.append(popt[2])
        chi2red_lst.append(chi2red)
        #plot and store
        print('Donly fraction is %.2f and FRET lifetime is %.2f ns' % (popt[1], popt[2]))
        fitDA.pltDA_eps(DAdat, D0dat, DAmodel, Donlymodel, Gfile, popt, chi2red, chi2red_D0, plotoutDir)
    df = pd.DataFrame ({
        'name' : names,
        'time' : times,
        'x_fret' : 1 - np.array(x0),
        'tau_fret (=1/k_fret)' : tau_fret,
        'k_fret' : 1 / np.array(tau_fret),
        'chi2_red' :  chi2red_lst,
        'total number of yellow photons' : S_y
    })
    df.to_csv(paramout)
    return names, S_y, x0, tau_fret, chi2red_lst
    

def analyzeRuler1D1A(locLst, resdir, identifier, ntacs = 256):
    """set of analysis parameters that should be applied to all Rulers"""
    outname = os.path.join(resdir, identifier + '_an.spots')
    TACout = os.path.join(resdir, identifier + '_1D1A_PS.dat')
    statsout = os.path.join(resdir, identifier + '.pg4')

    locLst_s = df.selectSpotOccurence(locLst, [1], [1])
    locLst_s_an = df.analyseLocLst(locLst_s, Igate = [0,0,0], ltgate = [25, 25, 25],
                                    verbose = False, framestop = 20, outname = outname,
                                    ntacs = ntacs, bgphotons = [1.1, 5.0, 5.6])
    #an data contains always the full TCSPC decay
    df.subensembleTAC(locLst_s_an, ntacs = ntacs, outfile = TACout)
    stats = df.genStats(locLst_s_an, outfile = statsout)
    return locLst_s_an, stats