import shutil
import pandas as pd
import copy
import numpy as np
import matplotlib.pyplot as plt
import ImageManipulation as IM
import os
import findPeaksLib
import GaussFits
import pickle
from GaussAnalysisPipeline import *
import precisionFuncs as pF
from scipy.optimize import minimize
from scipy.special import factorial
import copy
import lmfit

class FRETind:
    def __init__(self):
        self.no = None
        self.GTAC = None
        self.RTAC = None
        self.YTAC = None
        self.tauG = None
        self.tauR = None
        self.tauY = None
        self.NG = None
        self.NR = None
        self.NY = None
        self.Gbg = None
#        self.Rbg = None
        self.Ybg = None
        self.proxRatio = None
        self.stoichiometry = None
        self.distx = None
        self.disty = None
        self.dist = None

def analyseDir(
    options, wdir, files, Ggate = 0, Rgate = 0, Ygate = 0, 
    DTwoIstar = 0.03, garbageBrightness = 20, junkIstar = 0.30,
    outname = '', framestop = -1, ntacs = 256, ROIsize = 20,
    rebin = None, verbose = False, saveplot = False):
    """
    Analyse all ptu files listed in files. A single ROI is found in each file.
    In this ROI, either one, two or Three Gaussians are fitted.
    Each successfull fit is added as an entrie in locLst.
    input:
    locLst      target list to save all fits
    options     class object that contains fitting options
    wdir        filepath to working directory
    files       list of files in wdir
    Ggate       tac channel from which gating startswith for G channel
    Rgate       for Red channel
    YGate       for Yellow channel
    DtwoIstar   minimal 2I* improvement for a more complicated model to be
                considered better
    garbageBrightness
                minimal brightness for a peak to be accepted
    junkIstar   minimal absolute 2I* value for fit to be accepted
                noise images generate typically very low 2I* values
    outname     absolute path. if given, locLst is saved in outname
    framestop   how many frames of ptu file to collected
                -1 means collect all frames
    ntacs       number of channels for tac decay
    ROIsize     size in pixels of ROI side. ROI is square.
    rebin       integer that specifies rebinning factor in x and y_coord
    verbose     if True, plotfits are printed in the terminal
    saveplot    if True, plotfits are saved to disc
    
    """
    #make sure all variables are passed correctly
    locLst = []
    for i, file in enumerate(files):
        if file[-4:] != '.ptu':
            print('not a .ptu file, skipping')
            continue
        print('analysing image no. %i' %i)
        ffile = os.path.join(wdir, file)

        CLR = loadGRYimage(ffile, Ggate = Ggate, Rgate = Rgate,
            Ygate = Ygate, ntacs = ntacs, framestop = framestop,
            rebin = rebin)

        #get ROI based on Y signal
        ROI = getROI(CLR.workIntensity.Y, ROIsize)    
        
        #check that ROI is not touching borders:
        try:
            crop(CLR.workIntensity.G, ROI)
        except IndexError:
            print('ROI touches image borders, skipping')
            continue
            
        locLst.append({})
        locLst[-1]['filepath'] = ffile
        locLst[-1]['ROI'] = ROI
        #loop over G and Y channels
        for color, bitmap in CLR.workIntensity.__dict__.items():
            if color in ['G', 'R', 'Y']:
                if saveplot:
                    outdir = ffile[:-4] + '_' + color
                else:
                    outdir = ''

                snip = crop(bitmap, ROI)
                #fits 1, 2 or 3 gauss spots and determines which one is best
                #returns optimised parameter array
                bestfit, twoIstar, _ = fitNGauss (
                    snip, options, 
                    DTwoIstar = DTwoIstar, garbageBrightness = garbageBrightness,
                    junkIstar = junkIstar, verbose = verbose, outdir = outdir
                    )

                #build array containing all the fit results
                locLst[-1][color] = Channel(snip)
                locLst[-1][color].fillSpotLst(bestfit)
    
    if outname:
        with open(outname, 'wb') as output:
            pickle.dump(locLst, output, 1)
        
    return locLst


def analyseLoc(
    loc, cntr, Ggate = 0, Rgate = 0, Ygate = 0, 
    winSigma = 3, framestop = 20, rebin = None,
    TACCal = 0.0128, ntacs = 256,
    verbose = False, pxSize = 10):
    """finds all closest spot pairs
    for each pair the following are calculate:
        pair distance dist
        proximity ratio E 
        tau_f in Green channel
        tau_f in Red Channel
        tau_f in Yellow channel
    A dictionary entry 'FRETind' is added to input
    The entry points to an object containing the 
    FRET indicators
    input: 
        loc dictionary with 'G' and 'Y' entries
        Ggate: start of Ggate
        Ygate: start of Ygate
        ROIsize: size of ROI, should be phased out
        winSigma: half-side of area around spot center
            around which photons are collected
            e.g. winSigma = 3: 7x7 image"""
    #write all relevant parameters in function call
    #put calculation of FRET indicators in a separate function
    outloc = copy.deepcopy(loc)
    ROI = loc['ROI'] #shorthand
    
    #reload ptu file
    CLR = loadGRYimage(outloc['filepath'], Ggate = Ggate, Rgate = Rgate,
        Ygate = Ygate, ntacs = ntacs, framestop = framestop,
        rebin = rebin)
    #override existing snips with new snips
    for color, bitmap in CLR.workIntensity.__dict__.items():
        assert color in ['G', 'R', 'Y']
        snip = crop(bitmap, ROI)
        try:
            outloc[color].snip = snip
        except KeyError:
            outloc[color] = Channel(snip)
        if verbose:
            plt.imshow(snip)
            plt.show()
    #add pairwise localisation of loc to outloc
    sortSpots(outloc)    

    cntr = calcFRETind(CLR, outloc, winSigma, cntr, verbose, 
        Ggate, Rgate, Ygate, pxSize)
    return outloc, cntr
    
def analyseLocLst(locLst, Ggate = 0, Rgate = 0, Ygate = 0, 
    winSigma = 3, framestop = 20, rebin = None,
    TACCal = 0.0128, ntacs = 256,
    verbose = False, pxSize = 10, outname = ''):
    #add saving functionality
    outLst = []
    cntr = 0
    for loc in locLst:
        outloc, cntr = analyseLoc(loc, cntr, Ggate = Ggate, Rgate = Rgate, 
            Ygate = Ygate, winSigma = winSigma, framestop = framestop,
            rebin = rebin, verbose = verbose)
        outLst.append(outloc)
    if outname:
        with open(outname, 'wb') as output:
            pickle.dump(outLst, output, 1)
    return outLst
    
    
def loadGRYimage(
    ffile, Ggate = 0, Rgate = 0, Ygate = 0, 
    uselines = np.array([1,2]), ntacs = 256, framestop = -1, rebin = None):
    #load image, gate
    CLR = IM.processLifetimeImage(
        ffile.encode(), uselines = uselines, ntacs = ntacs,
        framestop = framestop)
    CLR.loadLifetime()
    if rebin:
        CLR.rebin(rebin, rebin)
    CLR.gate(Ggate, 150, channel = 'G')
    CLR.gate(Rgate, 150, channel = 'R')
    CLR.gate(Ygate, 150, channel = 'Y')
    CLR.loadIntensity()
    return CLR
    

def kickvector(dat, maxval):
    """
    filter all entries in dat whose norm is larger than maxval
    dat must be 2D array
    """
    return dat[np.linalg.norm(dat, axis = 1) < maxval]
    
def filterCoords(coords, maxval = 50):
    """kick outlier and center coordinates of LocLst"""
    coords = kickvector(coords, maxval)
    return coords - np.mean(coords, axis = 0)
    
def plotSinglePair(locLst, pxSize = 10):
    #split into fitting function and filtering function
    #filtering function goes into different function
    coords = np.zeros([len(locLst),2])
    for i, loc in enumerate(locLst):
        coords[i] = loc['G'].spotLst[0].coord - loc['Y'].spotLst[0].coord
    coords *= pxSize
    filterCoords(coords, 50)
    plt.plot(coords[:,0], coords[:,1], '.')
    plt.plot([0,0], [-100,+100], 'r--')
    plt.plot([-100,+100], [0,0], 'r--')
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlim([-50, 50])
    plt.ylim([-50, 50])
    stdx, stdy = np.std(coords, axis = 0)
    print('standard deviation of spots is %.2f and %.2f' % (stdx, stdy))
    plt.show()
    return coords

def plotOccurence(locLst):
    """plots a 2D histogram of how many spots have been fitted
    in the Green and red channel.
    Takes as argument a spotLst Lst"""
    occurence = np.zeros([4,4], np.int)
    for loc in locLst:
        Gspots = len(loc['G'].spotLst)
        Yspots = len(loc['Y'].spotLst)
        occurence[Gspots, Yspots] += 1
    plt.figure(figsize = [4,4])
    x, y = np.meshgrid(np.arange(4),np.arange(4))
    plt.scatter(x,y, s = occurence)
    ax = plt.gca()
    ax.set_xticks([0,1,2,3])
    ax.set_yticks([0,1,2,3])
    plt.ylabel ('# Green spots identified')
    plt.xlabel ('# red spots identified')
    plt.show
    return occurence
    
def selectSpotOccurence(locLst, Gspots, Yspots):
    """selects from locLst those localisations that have
    Gspots and Yspots. Gspots and Spots must be lists.
    returns: cleaned locLst"""
    i=0
    locLst_copy = locLst.copy()
    while i < len(locLst_copy):
        if (len(locLst_copy[i]['G'].spotLst) in Gspots and
                len(locLst_copy[i]['Y'].spotLst) in Yspots):
            i += 1
        else:
            locLst_copy.pop(i)
    return locLst_copy


def get_logLikelihood1DPoisson(params, func, xdata, ydata, sign = 1):
    """return log-likelihood. for 1D binned data function
    maximize this function to obtain most likely fit"""
    model = func(xdata, params)
    l = np.sum(-model + ydata * np.log(model) - np.log(factorial(ydata)))
    return l * sign
def get_logLikelihood(params, func, observations, sign = 1):
    """
    returns the log-likelihood that a series of observations are described
    by a model func with parameters params.
    observations are e.g. a unbinned arrival times of photons in a burst
    """
    loglikelihood = np.log(np.prod(func(observations, params)))
    return loglikelihood * sign
def ncChidistr(r, mu, sig, A, offset):
    """calculates non-centered Chi Distribution"""
    return A * r / sig**2 * np.exp(- (mu**2 + r**2) / (2 * sig**2)) \
        * np.i0( r * mu / (sig**2)) + offset
def NncChidistr(x, p):
    """takes lmfit Parameter object and generates as many ncChiDistr as there 
    are mu%i"""
    model = np.zeros(x.shape)
    v = p.valuesdict()
    N = 0
    while "mu%i" % N in v:
        model += ncChidistr(x, v['mu%i' % N], v['sig%i' % N], v['A%i' % N], 0)
        N += 1
    model += v['bg']
    return model
    
def genPeakEst(Npeaks, x, y):
    """
    returns lmfit parameter object for Npeaks ncChidistr
    """
    p = lmfit.Parameters()
    bg = np.mean(y) / 10
    mu = np.sum(x*y) / np.sum(y) #center of mass
    sig = np.sqrt(np.sum(y * (x - mu) ** 2) / np.sum(y)) # sigma
    A_max = np.sum(y)*10
    A = np.sum(y)
    p.add('bg', bg, True, 0, A_max)
    for i in range(Npeaks):
        p.add_many(('mu%i' % i, mu, True, 0, np.max(x)),
                   ('sig%i' % i, sig, True, 0, np.inf),
                   ('A%i' % i, A, True, 0, A_max))
    return p
    
#def ncChiTwodistr(r, mu1, sig1, A1, mu2, sig2, A2, offset):
#    """calculates to non-centered Chi distributions"""
#    distr1 = ncChidistr(r, mu1, sig1, A1, offset)
#    distr2 = ncChidistr(r, mu2, sig2, A2, 0)
#    return distr1 + distr2

def expDecay(r, tau, A):
    return A * np.exp( - r / tau )
def get_AIC(nFitVars, logLikelihood):
    """calculate Aikaike Information citerion"""
    return 2 * nFitVars - 2 * logLikelihood
def get_AICc(nFitVars, logLikelihood, samplesize):
    """corrected AIC, use for small samplesize"""
    k = nFitVars #shorthand
    n = samplesize #shorthand
    AIC = get_AIC(nFitVars, logLikelihood)
    return AIC + (2 * k ** 2 + 2 * k) / (n - k - 1)
def get_BIC(nFitVars, logLikelihood, samplesize):
    """calculate Bayesian information Criterion
    use for large sample sizes"""
    return np.log(samplesize) * nFitVars - 2 * logLikelihood
    
def estChiSigma(Gsigma, Ysigma, Gphotons, Yphotons, Gbg, Ybg, pxSize, posprecision):
    Gprecision = np.sqrt(pF.findVar([0, 0, Gsigma, Gbg, Gphotons], pxSize))
    Yprecision = np.sqrt(pF.findVar([0,0, Ysigma, Ybg, Yphotons], pxSize))
    chiSigma = np.sqrt(Gprecision**2 + Yprecision**2 + posprecision**2)
    return chiSigma

def fitNChidistr(dist, N = 1, outfile = '', binwidth = 2, 
                maxbin = 60, plotshow = True):
    """fits a chi distribution to a set of distances
    params are ordered: mu, sigma, amplitude, offset"""
    counts, bin_edges, _ = plt.hist(dist, bins = np.arange(0, maxbin, binwidth))
    plt.clf()
    Nbins = bin_edges.shape[0] - 1
    bins = np.zeros(Nbins)
    for i in range(Nbins):
        bins[i] = (bin_edges[i] + bin_edges[i + 1]) / 2
    p = genPeakEst(N, bins, counts)
    fitres = lmfit.minimize(get_logLikelihood1DPoisson, p, method = 'nelder',
        args = (NncChidistr, bins, counts, -1))
    logLikelihood = get_logLikelihood1DPoisson(fitres.params, NncChidistr, 
                                                bins, counts)
    AIC = get_AIC(fitres.nvarys, logLikelihood)
    AICc = get_AICc(fitres.nvarys, logLikelihood, len(dist))
    BIC = get_BIC(fitres.nvarys, logLikelihood, len(dist))
    return fitres, AIC, AICc, BIC
def whichChiIsBest(dist, verbose = False):
    fitresL = []
    AICL = []
    BICL = []
    for N in range(4):
        try:
            fitres, AIC, _, BIC = fitNChidistr(dist, N = N)
        except ValueError:
            print('fit for %i peaks failed' %N)
            break
        fitresL.append(fitres)
        AICL.append(AIC)
        BICL.append(BIC)
    
    bestfit = np.argmin(BICL)
    if verbose:
        for i in range(len(AICL)-1):
            relativeLikelihood = np.exp(AICL[i] - AICL[i+1])
            print('fitting %i peaks is %.2e times more likely than fitting %i peaks according to AIC'
                  % (i+1, relativeLikelihood, i))
        for i in range(len(BICL)-1):
            relativeLikelihood = np.exp(BICL[i] - BICL[i+1])
            print('fitting %i peaks is %.2e times more likely than fitting %i peaks according to BIC'
                  % (i+1, relativeLikelihood, i))
        #plot
        x = np.arange(0,60,.1)
        p = fitresL[bestfit].params
        model = NncChidistr(x, p)
        plt.plot(x, model)
        plt.plot()
        plt.hist(dist, bins = np.arange(0, 60, 2))
        plt.show()
    return fitresL[bestfit]

def fitChiDistr(dist, params0, bounds = None, outfile = '', binwidth = 2, maxbin = 60, plotshow = True):
    """fits a chi distribution to a set of distances
    params are ordered: mu, sigma, amplitude, offset"""
    counts, bin_edges, _ = plt.hist(dist, bins = np.arange(0, maxbin, binwidth))
    Nbins = bin_edges.shape[0] - 1
    bins = np.zeros(Nbins)
    for i in range(Nbins):
        bins[i] = (bin_edges[i] + bin_edges[i + 1]) / 2
    
    fitres = minimize(get_logLikelihood1DPoisson, params0, args = (ncChidistr, bins, counts, -1), 
      method = 'SLSQP', bounds = bounds, 
      options = {'maxiter': 1000, 'ftol': 1e-6})
    
    xgrid = np.arange(0,max(bins),0.1)
    fit = ncChidistr(xgrid, *fitres.x)
    plt.plot(xgrid, fit, label = 'MLE sigma fixed')
    if plotshow:
        plt.show()
        
    if outfile:
        savearray = np.array([xgrid, fit]).transpose()
        header = 'xgrid\t fit'
        np.savetxt(outfile, 
           savearray,
           fmt = '%.4f', 
           header = header, 
           delimiter = '\t')
        print(os.path.splitext(outfile)[0])
        fpath = os.path.splitext(outfile)[0] + '_fit_parameters.txt'
        with open(fpath, 'wt') as f:
            f.write('freefit\n')
            f.write('mu, sigma, A, bg\n')
            f.write(str(fitres.x) + '\n')
            f.write(str(bounds) + '\n')
    return fitres

def sortSpots(loc):
    """
    selects from unorganised Green and Yellow localisations pairs of the closest
    distance.
    First the closest pair is taken and removed from the list. 
    Then, the second closest (if present).
    if the localisations are not equally long, those matching a pair are taken.
    """
    loc_copy = copy.deepcopy(loc)
    NG = len(loc['G'].spotLst)
    NY = len(loc['Y'].spotLst)
    Ndist = min (NG, NY)
    Gcoords = np.zeros([NG, 2])
    Ycoords = np.zeros([NY, 2])
    alldist = np.zeros([NG, NY])
    for i, el in enumerate(loc_copy['G'].spotLst):
        Gcoords[i] = el.coord
    for i, el in enumerate(loc_copy['Y'].spotLst):
        Ycoords[i] = el.coord
    
    for i in range(NG):
        for j in range(NY):
            alldist[i, j] = np.linalg.norm(Gcoords[i] - Ycoords[j])
    loc['G'].spotLst = []
    loc['Y'].spotLst = []
    for i in range(Ndist):
        NG, NY = alldist.shape
        Gpos, Ypos = [alldist.argmin() // NY , alldist.argmin() % NY]
        #re-build channel objects with paired spots
        loc['G'].spotLst.append(loc_copy['G'].spotLst[Gpos])
        loc['Y'].spotLst.append(loc_copy['Y'].spotLst[Ypos])
        alldist[Gpos] = 1e6
        alldist[:, Ypos] = 1e6


def cropSpot(spotcenter, bitmap, winSigma):
    """crops 2D or 3D data in 2D
    spotcenter defines the coordinate of the spot in pixels in bitmap
    The total width of cropwindow is 2*winSigma + 1
    """
    #convert between float peak to pixel position
    xcenter, ycenter = (np.round(spotcenter)).astype(np.int)
    xstart, ystart, xstop, ystop = [xcenter - winSigma,
                                   ycenter-winSigma,
                                   xcenter + winSigma,
                                   ycenter + winSigma]
    return bitmap[xstart : xstop + 1, ystart : ystop + 1]
    
def fitTau(TAC, TACCal = 0.128, verbose = False, params0 = [1, 2]):
        """simple fitting and plotting function that fits 
        a monoexponential decay to a TAC decay using MLE optimization
        No boundaries are implemented.
        """
        tactimes = np.arange(len(TAC)) * TACCal
        # add TACcal/2 because average arrival time in first bin is TACcal / 2
        meantac = np.mean(tactimes * TAC) + TACCal / 2 
        if verbose:
            plt.plot(tactimes, TAC)
            plt.plot(tactimes, expDecay(tactimes, meantac, np.sum(TAC)))
            plt.show()
        return meantac

def calcFRETind(CLR, loc, winSigma, cntr, verbose, Ggate, Rgate, Ygate, pxSize):
    """calc FRET indicators based on intensity and lifetime information
    Takes a loc dict and edits properties of the 'FRETind' entry"""
    ROI = loc['ROI'] #shorthand
    npairs = len(loc['G'].spotLst)
    loc['FRETind'] = []
    for i in range(npairs):
        loc['FRETind'].append(FRETind())
    for i in range(npairs):
        #get No. Green photons
        loc['FRETind'][i].no = cntr
        Gsnip = cropSpot(loc['G'].spotLst[i].coord, 
                         loc['G'].snip, winSigma)
        Gphotons = np.sum(Gsnip)
        #get No. of Red photons using Y localisation
        Rsnip = cropSpot(loc['G'].spotLst[i].coord, 
                         loc['R'].snip, winSigma)
        Rphotons = np.sum(Rsnip)
        Ysnip = cropSpot(loc['Y'].spotLst[i].coord, 
                         loc['Y'].snip, winSigma)
        Yphotons = np.sum(Ysnip)
        loc['FRETind'][i].NG = Gphotons
        loc['FRETind'][i].NR = Rphotons
        loc['FRETind'][i].NY = Yphotons
        loc['FRETind'][i].proxRatio = Rphotons / (Gphotons + Rphotons)
        loc['FRETind'][i].stoichiometry = \
            (Gphotons + Rphotons) / (Gphotons + Rphotons + Yphotons)
            
        if verbose:
            fig, (ax1, ax2, ax3) = plt.subplots(1,3)
            ax1.imshow(Gsnip)
            ax1.set_title('green channel')
            ax2.imshow(Rsnip)
            ax2.set_title('red channel')
            ax3.imshow(Ysnip)
            ax3.set_title('Yellow channel')
            plt.show()

        loc['FRETind'][i].Gbg = loc['G'].spotLst[i].bg
        loc['FRETind'][i].Ybg = loc['Y'].spotLst[i].bg
        loc['FRETind'][i].GTAC = \
            np.sum(cropSpot(loc['G'].spotLst[i].coord + ROI[[0,1]],
            CLR.workLifetime.G, winSigma), axis = (0,1))[Ggate:150]
        loc['FRETind'][i].tauG = fitTau(loc['FRETind'][i].GTAC, 0.128, verbose)
        loc['FRETind'][i].RTAC = \
            np.sum(cropSpot(loc['Y'].spotLst[i].coord + ROI[[0,1]],
            CLR.workLifetime.R, winSigma), axis = (0,1))[Rgate:150]
        loc['FRETind'][i].tauR = fitTau(loc['FRETind'][i].RTAC, 0.128, verbose)        
        loc['FRETind'][i].YTAC = \
            np.sum(cropSpot(loc['Y'].spotLst[i].coord + ROI[[0,1]],
            CLR.workLifetime.Y, winSigma), axis = (0,1))[Ygate:150]
        loc['FRETind'][i].tauY = fitTau(loc['FRETind'][i].YTAC, 0.128, verbose)
        distvec = loc['G'].spotLst[i].coord - loc['Y'].spotLst[i].coord
        loc['FRETind'][i].distx = distvec[0] * pxSize
        loc['FRETind'][i].disty = distvec[1] * pxSize
        loc['FRETind'][i].dist = np.linalg.norm(distvec) * pxSize

        cntr += 1
    return cntr

def getFRETnames(locLst):
    """helper function that returns the names of FRET indicators"""
    names = []
    i = 0
    while len(locLst[i]['FRETind']) == 0:
        i += 1
        if i == 100000:
            raise RuntimeError
    for name, value in locLst[i]['FRETind'][0].__dict__.items():
        if type(value) != np.ndarray and name[0] != '_':
            names.append(name)
    return names

def genFRETstats(locLst, outfile = ''):
    """
    saves al FRET indicators of locLst in a Margarita-readable file.
    If there are more spots in a single localisation, they are treated as 
    independant.
    All non array value properties that are in the first non-empty entry 
    of locLst are saved
    """
    FRETind = {}
    names = getFRETnames(locLst)
    for name in names:
        FRETind[name] = []
    
    for loc in locLst:
        for spot in loc['FRETind']:
            for name, value in spot.__dict__.items():
                if name in names:
                    FRETind[name].append(value)
                    
    if outfile:
        print('saving FRET indicators to disc for Margarita')
        header = ''
        values = []
        for name, value in FRETind.items():
            header += name + '\t'
            values.append(value)
        np.savetxt(outfile,
                    np.array(values).transpose(),
                    delimiter = '\t',
                    header = header,
                    fmt = '%.3e'
                   )
        
    return FRETind, names
    
def loadpickle(outname):
    with open(outname, 'rb') as f:
        return pickle.load(f)

def simplePTUmerge(infolder, outfolder = ''):
    """copies all ptu files in subfolders of infolder to outfolder"""
    if outfolder == '':
        outfolder = os.path.join(infolder, 'all')
    print(outfolder)
    try:
        os.mkdir(outfolder)
        print('creating folder')
    except FileExistsError:
        pass
    folders = [f for f in os.listdir(infolder) if \
        os.path.isdir(os.path.join(infolder, f))]
    for folder in folders:
        subfolder = os.path.join(infolder, folder)
        files = os.listdir(subfolder)
        for file in files:
            if file.endswith('.ptu'):
                source = os.path.join(subfolder, file)
                destination = os.path.join(outfolder, file)
                shutil.copyfile(source, destination)

def filterFRETind(locLst, indicator, vmin, vmax):
    locLst_copy = copy.deepcopy(locLst)
    for loc in locLst_copy:
        for j in range(len(loc['FRETind']) -1, -1, -1): #loop last to first to avoid shifts
            val = getattr(loc['FRETind'][j], indicator)
            if val < vmin or val > vmax:
                loc['FRETind'].pop(j)
                loc['G'].spotLst.pop(j)
                loc['Y'].spotLst.pop(j)
                try:
                    loc['R'].spotLst.pop(j)
                except IndexError:
                    pass
    return locLst_copy

def AnIExport(locLst, outfolder, rebin = 1, sizex = 100):
    """export analysed locLst to AnI compatible format
    To calculate pixelid, rebinning factor and the unbinned xsize is needed
    e.g. a 100x100 10 nm CF image is rebinned to 20x20 pixels: 
        rebin = 5, imsize = 100
    """
    try:
        os.mkdir(outfolder)
    except FileExistsError:
        pass
        
    FRETnames = getFRETnames(locLst)
    headerbase = 'image_ID_cI\tPixel Number_cI'
    
    for loc in locLst:
        locname = os.path.split(loc['filepath'])[1][:-4]
        for color in ['G', 'R', 'Y']:
            if color == 'G':
                prepend = '_Green '
                append = '_cI'
                ext = '.0I4'
            elif color == 'R':
                prepend = '_Red '
                append = '_cI'
                ext = '.II4'
            elif color == 'Y':
                prepend = '_Yellow '
                append = '_cI'
                ext = '.II4'
            
            spotLst = loc[color].spotLst
            #naming convention from Suren
            outfile = os.path.join(outfolder, locname + prepend + \
                    'Photons_Fit results' + ext)
            header = headerbase
            if len(spotLst) != 0:
                for name, value in spotLst[0].__dict__.items():
                    if name == 'coord':
                        header += '\tposx' + append + '\tposy' + append
                    else:
                        header += '\t' + name + append

                if color == 'G':
                    for name in FRETnames:
                        header += '\t' + name
            header += '\n'
            
            with open(outfile, 'wt') as f:
                f.write(header)
                for i, spot in enumerate(spotLst):
                    xpos, ypos = loc['ROI'][:2] * rebin + \
                    np.round(spot.coord * rebin)
                    pixid = ypos * sizex + xpos +1
                    line = '0\t%d' % pixid
                    for name, value in spot.__dict__.items():
                        if name == 'coord':
                            line += '\t%.3f\t%.3f' % \
                                (value[0] * rebin, value[1] * rebin)
                        else:
                            line += "\t%.3f" % value
                    if color == 'G':
                        FRETind = loc['FRETind'][i]
                        for name, value in FRETind.__dict__.items():
                            if name in FRETnames:
                                line += "\t%.3f" % value
                    line += '\n'
                    f.write(line)
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:35:29 2019

@author: buddeja
"""


def mergePTUfiles(path, savepath):
    files = os.listdir(path)
    file_folder = []


    if os.path.exists(savepath) == True:
         print('savepath_exists')
    else:
        os.makedirs(savepath)  

    for i in files:
        if i.startswith('Overview'):
            file_folder.append(i)

    for ii in file_folder:
        print(ii)

        files_ptu = []
        files_xcoord = []
        files_ycoord = []
        x_relative = []
        y_relative = []
        folderpath = os.path.join(path, ii)#'{}{}{}'.format(path, ii, '/')
        files_folderpath = os.listdir(folderpath)
        for i in files_folderpath:
            if i.endswith('.ptu'):
                files_ptu.append(i)
            elif i.endswith('x_transfer.dat.dat'):
                x_coord = i
            elif i.endswith('y_transfer.dat.dat'):
                y_coord = i
        if len(files_ptu) >0 :

            filex = pd.read_csv('{}{}'.format(folderpath,x_coord),sep='\t')
            filey = pd.read_csv('{}{}'.format(folderpath,y_coord),sep='\t')
            filex = [i[0] for i in filex.values.tolist()]
            filey = [i[0] for i in filey.values.tolist()]
            x_relative.extend(filex) 
            y_relative.extend(filey) 
            for i in range(len(files_ptu)):
                name_new = '{}{}{}{}{}{}'.format(files_ptu[i][:-4],'_x_',x_relative[i],'_y_',y_relative[i],'.ptu')
                shutil.copy('{}{}'.format(folderpath, files_ptu[i]), os.path.join(savepath, name_new))
