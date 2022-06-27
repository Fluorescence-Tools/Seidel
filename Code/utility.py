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
import precisionFuncs as pF
import LSManalysis as LSMan
import batchplot as bp




    
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
                
                
def reportLocStats(locLst, outname = None):
    """prints some statistics of a locLST dataset"""
    s = ''
    stats = df.genStats(locLst)
    # #below section prints yields for cut, but the cut is hardcoded
    # s += 'filter yields:\n'
    # stats = df.genStats(locLst)
    # s += 'total images taken: %i\n' % len(locLst)
    # scut = df.filterStats(stats, 'stoichiometry', 0.2, 0.8)
    # s += 'total FRET pairs: %i\n' % len(scut['posxG'])
    # NF = df.filterStats(scut, 'proxRatio', 0, 0.6)
    # s += 'NF pairs: %i\n' % len(NF['posxG'])
    # stats = df.genStats(locLst)
    # scut = df.filterStats(stats, 'stoichiometry', 0.2, 0.8)
    # HF = df.filterStats(scut, 'proxRatio', 0.75, 1)
    # s += 'HF pairs: %i\n' % len(HF['posxG'])
    # s += '\n' + '%'*80 + '\n\n'
    
    s += 'loc Accury indicators:\n'
    locNG = np.array(stats['AG']) * np.array(stats['sigmaG'])**2 *2 * np.pi
    locNY = np.array(stats['AY']) * np.array(stats['sigmaY'])**2 *2 * np.pi
    locNGmean = np.mean(locNG); locNGstd =np.std(locNG)
    s += 'Green photons is %.0f +- %.0f\n' % (locNGmean, locNGstd)
    sigmaG = np.array(stats['sigmaG'])*10; sigmaY = np.array(stats['sigmaY'])*10
    sigmaGmean = np.mean(sigmaG); sigmaGstd = np.std(sigmaG)
    s += 'green spot sigma is %.2f +- %.2f\n' % (sigmaGmean, sigmaGstd)
    s+= 'green spot FWHM is %.2f +- %.2f\n' % (sigmaGmean * 2.355, sigmaGstd * 2.355)
    bgG = stats['bgG']; bgY = stats['bgY']
    bgGmean = np.mean(bgG); bgGstd = np.std(bgG)
    s += 'Green background is %.2f +- %.2f\n' % (bgGmean, bgGstd)
    locNYmean = np.mean(locNY); locNYstd = np.std(locNY)
    s += 'Yellow photons is %.0f +- %.0f\n' % (locNYmean, locNYstd)
    sigmaYmean = np.mean(sigmaY); sigmaYstd = np.std(sigmaY)
    s += 'Yellow spot sigma is %.2f +- %.2f\n' % (sigmaYmean, sigmaYstd)
    s += 'Yellow spot FWHM is %.2f +- %.2f\n' % (sigmaYmean * 2.355, sigmaYstd * 2.355)
    bgYmean = np.mean(bgY); bgYstd = np.std(bgY)
    s += 'Yellow background is %.2f +- %.2f\n' % (bgYmean, bgYstd)
    #slow, arrays are much faster, quad function in findVar does not take arr
    Gprecision = [np.sqrt(pF.findVar([0,0,sigmaG, bgG, locNG], 10)) \
        for sigmaG, bgG, locNG in zip(sigmaG, bgG, locNG)]
    GprecisionMean = np.sqrt(pF.findVar([0,0,sigmaGmean, bgGmean, locNGmean], 10))
    s += 'standard deviation of G Channel is %.2f nm\n' % GprecisionMean
    Yprecision = [np.sqrt(pF.findVar([0,0,sigmaY, bgY, locNY], 10)) \
        for sigmaY, bgY, locNY in zip(sigmaY, bgY, locNY)]
    YprecisionMean = np.sqrt(pF.findVar([0,0,sigmaYmean, bgYmean, locNYmean], 10))
    s += 'standard deviation of Y Channel is %.2f nm\n' % YprecisionMean
    posprecision = 0
    s += 'uncertainty in dye position is assumed to be %.2f nm \n' % posprecision
    totprecision = [np.sqrt(Gprecision**2 + Yprecision**2) \
        for Gprecision, Yprecision in zip(Gprecision, Yprecision)]
    chiSigma = np.sqrt (GprecisionMean**2 + YprecisionMean**2 + posprecision**2)
    s += 'expected sigma of chi distribution is %.2f nm\n' % chiSigma
    
    def bindata(data, bin_edges, handle):
        binned_data, *_ = plt.hist(data, bins = bin_edges)
        bincenters = (bin_edges[1:]+bin_edges[:-1])/2
        dfrm[handle+'bins'] = bincenters
        dfrm[handle] = binned_data
        plt.clf()
        return binned_data
    dfrm = pd.DataFrame()
    precisionbins = np.linspace(0,20,51)
    bindata(totprecision, precisionbins, 'presision')
    locRateGbins = np.linspace(0, 150, 51)
    locRateYbins = locRateGbins
    integrationTime = 60*7*7*5e-6 #60 frames
    bindata(locNG / integrationTime / 1000, locRateGbins, 'locRateG') #kHz
    bindata(locNY / integrationTime / 1000, locRateYbins, 'locRateY') #kHz
    bgGbins = np.linspace(0,2,51)
    bgYbins = bgGbins 
    bindata(bgG, bgGbins, 'bgG') 
    bindata(bgY, bgYbins, 'bgY')
    sigmaGbins = np.linspace(0,100, 51)
    sigmaYbins = sigmaGbins
    bindata(sigmaG, sigmaGbins, 'sigmaG')
    bindata(sigmaY, sigmaYbins, 'sigmaY')
    dfrmout = outname[:-4]+'_data.csv'
    dfrm.to_csv(dfrmout, index = False)
    
    if outname:
        f = open(outname, 'w')
        f.write(s)
        f.close()
    print (s)
    
    
#####################LSM analysis utility#######################################
def loadLSMUtility(wdir, identifier, load = False, **kwargs):
    """utility function for collecting a bunch of functions often used together"""
    picklepath = os.path.join(wdir, 'results', identifier + '.pickle')
    if load:
        SampleSet = aid.loadpickle(picklepath)
    else:
        SampleSet = LSMan.sampleSet(wdir,
                                    **kwargs)
        SampleSet.analyzeDir(identifier, **kwargs)
        aid.savepickle(SampleSet, picklepath)
    return SampleSet
    
def loadLSMUtilitywMasks(wdir, identifier, load = False, **kwargs):
    """utility function for collecting a bunch of functions often used together"""
    picklepath = os.path.join(wdir, 'results', identifier + '.pickle')
    if load:
        SampleSet = aid.loadpickle(picklepath)
    else:
        SampleSet = LSMan.sampleSet(wdir,
                                    **kwargs)
        SampleSet.analyzeDirwMasks(identifier, **kwargs)
        aid.savepickle(SampleSet, picklepath)
    return SampleSet
    
def FitPlotLSMUtility(SampleSet, normImageG, normImageY, identifier,
                      fitfunc = 'batchFit2lt',
                      fitkwargs = LSMan.genDefaultFitKwargs(),
                      colorby = 'rateY_LPC',
                      **kwargs
                      ):
    SampleSet.batchgenNormDecay(normImageG, normImageY, **kwargs)
    if fitfunc is not None:
        getattr(SampleSet, fitfunc)(identifier, **fitkwargs)
    colorcoding = SampleSet.imstats[colorby]
    bp.pltRelativeDecays(SampleSet, identifier, decaytype = fitkwargs['decaytype'],
                         colorcoding = colorcoding,
                         resdir = SampleSet.resdir, **kwargs)
                         
def prepareForCreatingManualMask(wdir):
    #list all ptu files in folder
    imagefiles = [file for file in os.listdir(wdir)\
                  if file.endswith('.ptu')]
    #make a target directory for saving images for masking
    imagesoutdir = os.path.join(wdir, 'imagesForMasking')
    aid.trymkdir(imagesoutdir)
    for file in imagefiles:
        #generate a directory for manually placing masks
        maskdir = file[:-4] + '_Masks'
        aid.trymkdir(os.path.join(wdir, maskdir))
        #export an image for creating masks
        ffile = os.path.join(wdir, file)
        ltimg = IM.processLifetimeImage(ffile.encode())
        ltimg.loadLifetime()
        ltimg.loadIntensity()
        preposition = file[:-4]
        ltimg.saveWorkIntensityToTiff(imagesoutdir, preposition = preposition)