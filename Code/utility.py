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
import batchplot



    
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
                
class sampleSet():
    """collection of attributes specific to either Donor only, or acceptor only"""
    def __init__(self, wdir):
        self.wdir = wdir
        self.TACdir = os.path.join(wdir, 'TAC')
        aid.trymkdir(self.TACdir)
        self.ptufiles = bp.appendOnPattern(wdir, 'ptu')
        self.getTACfiles()
        self.TACs = []
    def getTACfiles(self):
        self.TACfiles = bp.appendOnPattern(self.TACdir, '_G_PS.dat')
    def appendTAC(self, TAC):
        assert type(TAC) == np.ndarray
        self.TACs.append(TAC)

class lsm_analysis():
    """contains some shared stuff for LSM analysis functions"""
    def __init__(self, 
        wdir, 
        ntacs = 1024, 
        pulsetime = 50, 
        dwelltime = 20e-6,
        TACrange = 4096):
        """        
        input:
            ntacs: number of TAC channels
            pulsetime: inverse of laser repetition rate, in ns
            dwelltime: pixel dwell time in seconds
            Nframes: number taken in imreading and for calculating total 
                illumination time
            threshold: all pixels below this threshold are set to zero
            TAC_range: set in hydraharp
        """
        self.imStatsHeader = ['name', 'I_G', 'surface_G', 'Br_G', 'rate_G', 'I_Y', 'surface_Y', 'Br_Y', 'rate_Y']
        self.wdir = wdir
        D0dir = os.path.join(wdir, 'D0')
        DAdir = os.path.join(wdir, 'DA')
        self.D0 = sampleSet(D0dir)
        self.DA = sampleSet(DAdir)
        self.ntacs = ntacs
        self.pulsetime = pulsetime
        self.dwelltime = dwelltime
        self.TACrange = TACrange
        self.resdir = os.path.join(wdir, 'results')
        aid.trymkdir(self.resdir)
    def genOldStyleHeader(self):
        self.imStatsHeader =  ['name', 'surface', 'I_G', 'Br_G', 'rate_G', 'I_Y', 'Br_Y', 'rate_Y']
        
    def analyzeGYLSMCells(self,
        sampleSetID,
        imStatsOut = None,
        threshold = 50,
        Nframes = -1):
        """
        This script is intended to automate image analysis for cellular data to avoid 
        time-consuming manual work in AnI.
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
        input: (partly given in sampleSet object)
            ffiles: full file path of all cells to be analysed
            ntacs: number of TAC channels
            pulsetime: inverse of laser repetition rate, in ns
            dwelltime: pixel dwell time in seconds
            Nframes: number taken in imreading and for calculating total 
                illumination time
            threshold: all pixels below this threshold are set to zero
            TAC_range: set in hydraharp
            sampleSetID is e.g. 'D0' or 'DA' for self.D0 and self.DA respectively
        """
        
        #init
        setID = sampleSetID # shorthand
        files = getattr(self, setID).ptufiles
        ffu
        df = pd.DataFrame(columns = self.imStatsHeader)
        uselines = np.array([1]) 
        #loop over each cell
        for index, file in enumerate(files):
            ffile = os.path.join(getattr(self, setID).wdir, file)
            #load GRY object
            try: ffile.encode()
            except: pass
            imD = IM.processLifetimeImage(
                    ffile, 
                    uselines = uselines, 
                    Gchan = np.array([1]),
                    Rchan = np.array([0]),
                    ntacs = self.ntacs,
                    TAC_range = self.TAC_range,
                    pulsetime = self.pulsetime,
                    dwelltime = self.dwelltime,
                    framestop = int(Nframes)) #some bug changed type of Nframes to float
            imA = IM.processLifetimeImage(
                    ffile, 
                    uselines = uselines, 
                    Gchan = np.array([5]),
                    Rchan = np.array([4]),
                    ntacs = self.ntacs,
                    TAC_range = self.TAC_range,
                    pulsetime = self.pulsetime,
                    dwelltime = self.dwelltime,
                    framestop = int(Nframes))
            
            #need to build a dict row for adding to DataFrame
            dfrow = pd.Series()
            TACout = os.path.join(getattr(self, setID).TACdir, file[:-4] + '_G_PS.dat')
            #load each channel, calculate parameters and store in dict
            for image, label in zip([imD, imA], ['_G', '_Y']):
                image.loadLifetime()
                image.loadIntensity()
                mask = image.buildMaskFromIntensityThreshold(threshold = threshold, sumchannels = ['G', 'R'])
                image.mask(mask)
                image.loadIntensity()
                #extract lifetime decays
                TACP, TACS, _ = image.getTACS()
                TACPS = np.zeros(self.ntacs*2)
                TACPS[:self.ntacs] = TACP
                TACPS[self.ntacs:] = TACS
                #extract intensity
                Nphot = np.sum(TACP +TACS)
                #extract surface
                surface = np.sum(mask)
                #calculated derived variables intensity and countrate
                brightness = Nphot / surface
                countrate = brightness / (dwelltime * Nframes)
                dfrow = dfrow.append(pd.Series({'I'+label: Nphot,
                           'surface'+label: surface,
                           'Br'+label: brightness,
                           'rate'+label: countrate})
                    )
            #transfer series to dataframe
            #I don't like the pandas ApI, it loses the Series name when another is
            #added
            dfrow.name = file[:-4]
            df = df.append(dfrow)
            np.savetxt(TACout, TACPS, fmt = '%i')
            print('finished with file %s\n' % file[-20:])
        try:
            df['Sg/Sy'] = df['Br_G'] / df['Br_Y']
        except KeyError:
            print('keywords not found in header, skipping calculated variables')
        if imStatsOut: df.to_csv(statsOut)
        return df
    
    def batchFitD0DA(self, 
                      D0dat,
                      identifier,
                      DAID = 'DA',
                      fitrange = (0, -1),
                      dataselect = (0, -1)):
        """makes simple Donor Only calibrated Donor Acceptor fits
        Reads existing csv data and returns the result.
        overrideHeader should be true to calculate derived variables"""
        #read all DA decays
        DATACs = getattr(self, DAID).TACs[dataselect[0]: dataselect[1]]
        #fit and plot DA
        xFRET = []
        kFRET = []
        chi2redLst = []
        names = []
        
        plotout = os.path.join(self.resdir, identifier + 'D0DAplots')
        aid.trymkdir(plotout)
        for i, DATAC in enumerate(DATACs):
            D0snip = D0dat[fitrange[0]:fitrange[1]]
            DAsnip = DATAC[fitrange[0]:fitrange[1]]
            _, _, _, Donlymodel, chi2red_D0 = fitDA.fitDonly(D0snip)
            popt, pcov, DAmodel, chi2red = fitDA.fitDA (DAsnip, D0snip)
            fitDA.pltDA_eps(DAsnip, D0snip, DAmodel, Donlymodel, DAfiles[i], popt, 
                            chi2red, chi2red_D0, plotout)
            xFRET.append(1-popt[1])
            kFRET.append(1/popt[2])
            chi2redLst.append(chi2red)
            names.append(DAfiles[i])
            
        df = pd.DataFrame(index = names)
        df['xFRET'] = xFRET
        df['kFRET'] = kFRET
        df['chi2red'] = chi2redLst
        outname = os.path.join(self.resdir, identifier + 'D0DAFitData.csv')
        df.to_csv(outname)
        return df
        
    def batchFitD0(self,
                    identifier,
                    D0ID = 'D0'
                    fitrange = (0, -1)):
        """batch fit D0 data assuming twop lifetimes"""

        D0TACs = getattr(self, DOID).TACs
        
        #prep empty arrays
        x0Lst = []
        x1Lst = []
        tau0Lst = []
        tau1Lst = []
        chi2redLst = []
        names = []
        bgLst = []
        
        #fit all data
        for i, D0TAC in enumerate(D0TACs):
            D0snip = D0TAC[fitrange[0]:fitrange[1]]
            popt, _, _, Donlymodel, chi2red = fitDA.fitDonly(D0snip)
            #append to Lists
            for var, val in zip([x0Lst, x1Lst, tau0Lst, tau1Lst, bgLst], popt):
                var.append(val)
            chi2redLst.append(chi2red)
            names.append(D0files[i][:-len('_G_PS.dat')])
            
        #calc derived vars
        taufLst = [(x0 * tau0**2 + x1 * tau1**2) / (x0 * tau0 + x1 * tau1) 
                    for x0, x1, tau0, tau1 in zip(x0Lst, x1Lst, tau0Lst, tau1Lst)]
        tauxLst = [(x0 * tau0 + x1 * tau1) / (x0 + x1)
                    for x0, x1, tau0, tau1 in zip(x0Lst, x1Lst, tau0Lst, tau1Lst)]
        #add all lists to DataFrame
        df = pd.DataFrame({'tauf' : taufLst,
                            'taux' : tauxLst, 
                            'x0' : x0Lst, 
                            'x1' : x1Lst, 
                            'tau0' : tau0Lst, 
                            'tau1' : tau1Lst, 
                            'chi2red' : chi2redLst}, index = names)
        outname = os.path.join(self.resdir, identifier + 'D0FitData.csv')
        df.to_csv(outname)
        return df



#def exportLSMTAC(fname, outdir, dwelltime, pulsetime, uselines = np.array([1]), 
#                 Gchan = np.array([0,1]), Rchan = np.array([4,5]), Ychan = np.array([4,5]),
#                ntacs = 1024, TAC_range = 4096, PIE_gate = 440):
#    """utility function to export GRY (P+S) TAC histogram for LSM microscope,
#    
#    deprecated. Function was written some time in the past and then forgotten about.
#    superceded by analyzeLSMCells
#    """
    ## issue: when Rchan and Ychan are identical (i.e. true PIE), then Y channel
    ## remains empty. Workaround implemented. ugly
    # _, file = os.path.split(fname)
    # data = IM.processLifetimeImage(fname.encode(), uselines = uselines, Gchan = Gchan, 
                                   # Rchan = Rchan, 
                                   # Ychan = Ychan, ntacs = ntacs, TAC_range = TAC_range, \
                                   # pulsetime = pulsetime, dwelltime = dwelltime)
    # data.loadLifetime()
    # if (Rchan == Ychan).all(): # workaround in case of PIE
       # data.workLifetime.Y = copy.deepcopy(data.workLifetime.R)
    # data.gate(0,PIE_gate, channel = 'R')
    # data.gate(PIE_gate, -1, channel = 'Y')
    # GTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'G'))).transpose()
    # RTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'R'))).transpose()
    # YTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'Y'))).transpose()
    # np.savetxt(os.path.join(outdir, file[:-4] + '_G.dat'), GTACS, fmt = '%i')
    # np.savetxt(os.path.join(outdir, file[:-4] + '_R.dat'), RTACS, fmt = '%i')
    # np.savetxt(os.path.join(outdir, file[:-4] + '_Y.dat'), YTACS, fmt = '%i')