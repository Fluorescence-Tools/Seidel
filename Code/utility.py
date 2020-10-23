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
import GaussAnalysisPipeline as GAP
import aid_functions as aid


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
    

def analyzeRuler1D1A(locLst, resdir, identifier, ntacs = 256, 
    select11stoic = True):
    """set of analysis parameters that should be applied to all Rulers"""
    outname = os.path.join(resdir, identifier + '_an.spots')
    TACout = os.path.join(resdir, identifier + '_1D1A_PS.dat')
    statsout = os.path.join(resdir, identifier + '.pg4')

    if select11stoic:
        locLst = df.selectSpotOccurence(locLst, [1], [1])
    locLst_an = df.analyseLocLst(locLst, Igate = [0,0,0], ltgate = [25, 25, 25],
                                    verbose = False, framestop = 20, outname = outname,
                                    ntacs = ntacs, bgphotons = [1.1, 5.0, 5.6])
    posdir = os.path.join(resdir, 'positions')
    df.export_position(locLst, posdir)
    #an data contains always the full TCSPC decay
    df.subensembleTAC(locLst_an, ntacs = ntacs, outfile = TACout)
    stats = df.genStats(locLst_an, outfile = statsout, isforMargarita = True)
    return locLst_an, stats
    
def analyzeOrigami(locLst, resdir, identifier, ntacs = 256):
    """set of analysis parameters that should be applied to all Origamis"""
    outname = os.path.join(resdir, identifier + '_an.spots')
    TACout = os.path.join(resdir, identifier + '_all_PS.dat')
    statsout = os.path.join(resdir, identifier + '.pg4')

    locLst_an = df.analyseLocLst(locLst, Igate = [0,0,0], ltgate = [29, 29, 29],
                                    verbose = False, framestop = 20, outname = outname,
                                    ntacs = ntacs, bgphotons = [1.1, 5.0, 5.6])
    posdir = os.path.join(resdir, 'positions')
    df.export_position(locLst, posdir)
    #an data contains always the full TCSPC decay
    df.subensembleTAC(locLst_an, ntacs = ntacs, outfile = TACout)
    stats = df.genStats(locLst_an, outfile = statsout, isforMargarita = True)
    return locLst_an, stats
    
def exportImageAndFit(loc, xstart, ystart, size = 30, Nspots = 2, outdir = None, verbose = False):
    """utility function to export attractive images for display in paper.
    Loc is a data struct to store localisation and FRET information
    The df.analyzeloclst routine already gates the data. Therefore the user should know what gate was used
        to generate these images.
    xstart, ystart and size are used to define the snapshot.
    outdir defines an out directory. the filename is automatically inferred from from the fpath parameter
    in the loc struct. 
    verbose dumps what could also be saved to disc."""
    for C in ['G', 'Y']:

        #get image
        snip = loc[C].bitmap[xstart: xstart + size, ystart: ystart + size]

        #fit image
        FitOpts = GAP.optionsCluster(fitbg = 0, setmodel = Nspots -1)
        param_est = GAP.genParamEstimate(snip)
        FitOpts.transferOptions(param_est)
        param_opt = GAP.GaussFits.Fit2DGauss(param_est, snip)
        model = np.zeros(snip.shape)
        if Nspots ==1: 
            model = np.array(GAP.GaussFits.model2DGaussian(param_opt, model))
        if Nspots == 2:
            model = np.array(GAP.GaussFits.modelTwo2DGaussian(param_opt, model))
        if Nspots == 3:
            model = np.array(GAP.GaussFits.modelThree2DGaussian(param_opt, model))
        if verbose:
            plt.imshow(loc[C].bitmap, cmap = 'hot')
            plt.show()
            plt.imshow(snip, cmap = 'hot')
            plt.show()
            plt.imshow(model, cmap = 'hot')
            plt.colorbar()
            plt.show()
            print(param_opt)

        #save fit and image
        if outdir is not None:
            df.aid.trymkdir(outdir)
            fid = loc['filepath'][-20:]
            outpath = os.path.join(outdir, fid + '_data' + C + '.csv')
            np.savetxt(outpath, snip, delimiter = ',')
            outpath = os.path.join(outdir, fid + '_model' + C + '.csv')
            np.savetxt(outpath, model, fmt = '%.5e', delimiter = ',')
            outpath = os.path.join(outdir, fid + '_model_perc' + C + '.csv')
            np.savetxt(outpath, model / np.max(model)*100, fmt = '%.5e', delimiter = ',')
            outtext = os.path.join(outdir, 'README.txt')
            with open(outtext, 'w') as f:
                f.write('original file path of dataset is:\n')
                f.write(loc['filepath']+ '\n')
                f.write('xstart: %i\n' % xstart)
                f.write('ystart: %i\n' % ystart)
            pickleout = os.path.join(outdir, 'loc' +fid+ '.spots')
            with open(pickleout, 'wb') as output:
                df.pickle.dump(loc, output, 1)
                
def analyzeLSMCells(ffiles, TACoutdir, statsOut, ntacs, pulsetime, 
    dwelltime, Nframes, threshold = 50, TAC_range = 4096):
    """
    This script is intended to automate image analysis for cellular data to avoid 
    tim-consuming manual work in AnI.
    To work, this script neads a functioncal copy of Seidel in the pythonpath
    The processLifetimeImage is not build for Anisotropy, but it can if one 
        mis-uses the channels. I.e. processlifetimeImage takes up to 3 channels 
        labelled Green, Red, Yellow. Now we will abuse by doing:
            Green = parallel
            Red = perpendicular
            Yellow = unused
            Then repeating for both channels
    This workaround should hold for the forseeable future, but should ultimately be 
        replaced.
    output:
        TACPS channel for Green (could add Yellow also)
        intensity, surface, brightness, Countrate for Green, Red+Yellow
    input:
        ffiles: full file path of all cells to be analysed
        TACoutdir: location where TAC PS decays are saved
        ntacs: number of TAC channels
        pulsetime: inverse of laser repetition rate, in ns
        dwelltime: pixel dwell time in seconds
        Nframes: number taken in imreading and for calculating total 
            illumination time
        threshold: all pixels below this threshold are set to zero
        TAC_range: set in hydraharp
        
    """
    #init
    df = pd.DataFrame(columns = ['intensity_G',
                                 'surface_G', 
                                 'brightness_G',
                                 'countrate_G',
                                 'intensity_Y',
                                 'surface_Y',
                                 'brightness_Y',
                                 'countrate_Y'])
    aid.trymkdir(TACoutdir)
    uselines = np.array([1]) 
    #loop over each cell
    for index, ffile in enumerate(ffiles):
        ffile = ffiles[index]
        #load GRY object
        imD = IM.processLifetimeImage(
                ffile, 
                uselines = uselines, 
                Gchan = np.array([1]),
                Rchan = np.array([0]),
                ntacs = ntacs,
                TAC_range = TAC_range,
                pulsetime = pulsetime,
                dwelltime = dwelltime,
                framestop = int(Nframes)) #some bug changed type of Nframes to float
        imA = IM.processLifetimeImage(
                ffile, 
                uselines = uselines, 
                Gchan = np.array([5]),
                Rchan = np.array([4]),
                ntacs = ntacs,
                TAC_range = TAC_range,
                pulsetime = pulsetime,
                dwelltime = dwelltime,
                framestop = int(Nframes))
        
        #need to build a dict row for adding to DataFrame
        pdrow = {}
        TACout = os.path.join(TACoutdir, 'cell%i_G_PS.dat' % index)
        #load each channel, calculate parameters and store in dict
        for image, label in zip([imD, imA], ['_G', '_Y']):
            image.loadLifetime()
            image.loadIntensity()
            mask = image.buildMaskFromIntensityThreshold(threshold = 50, sumchannels = ['G', 'R'])
            image.mask(mask)
            image.loadIntensity()
            #extract lifetime decays
            TACP, TACS, _ = image.getTACS()
            TACPS = np.zeros(ntacs*2)
            TACPS[:ntacs] = TACP
            TACPS[ntacs:] = TACS
            #extract intensity
            Nphot = np.sum(TACP +TACS)
            #extract surface
            surface = np.sum(mask)
            #calculated derived variables intensity and countrate
            brightness = Nphot / surface
            countrate = brightness / (dwelltime * Nframes)
            pdrow.update({'intensity'+label: Nphot,
                       'surface'+label: surface,
                       'brightness'+label: brightness,
                       'countrate'+label: countrate}
                )
        #transfer dict to dataframe
        df = df.append(pdrow, ignore_index = True)
        np.savetxt(TACout, TACPS, fmt = '%i')
        print('finished with file %s\n' % ffile[-20:])
    df.to_csv(statsOut)