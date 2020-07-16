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
import GaussAnalysisPipeline as GAP
import precisionFuncs as pF
from scipy.optimize import minimize
from scipy.special import factorial
import copy
import lmfit
import aid_functions as aid
from skimage import feature
from itertools import product

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
    rebin = None, verbose = False, saveplot = False,
    min_distance = 15, ROI_threshold_rel = 0.3, ROI_threshold_abs = 1,
    gateStop = 150):
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
    min_distance
                minimal distance of two ROIs
    ROI_threshold
                threshold for peak identification to determine ROIs
    
    """
    if ROIsize > 2 * min_distance:
        print('warning: ROIsize larger than 2 * min_distance, ROIs may overlap')
    locLst = []
    for i, file in enumerate(files):
        if file[-4:] != '.ptu':
            print('not a .ptu file, skipping')
            continue
        print('analysing image no. %i' %i)
        ffile = os.path.join(wdir, file)

        CLR = loadGRYimage(ffile, Ggate = Ggate, Rgate = Rgate,
            Ygate = Ygate, ntacs = ntacs, framestop = framestop,
            rebin = rebin, gateStop = gateStop)
        #make smooth intensity image to select ROIs
        CLR_smooth = copy.deepcopy(CLR)
        CLR_smooth.smoothIntensity(sigma = 2)
        loc = {}
        loc['filepath'] = ffile
        #loop over G, R and Y channels
        for color in ['G', 'R', 'Y']:
            bitmap = getattr(CLR.workIntensity, color)
            loc[color] = GAP.Channel(bitmap)
            if saveplot:
                outdir = ffile[:-4] + '_' + color
            else:
                outdir = ''
            peakimg = getattr(CLR_smooth.workIntensity, color)
            #find ROIs in image
            pos = feature.peak_local_max(peakimg, \
                min_distance = min_distance, \
                threshold_rel = ROI_threshold_rel,
                threshold_abs = ROI_threshold_abs)
            ROIs = aid.pos2ROI(pos[:,0], pos[:,1], ROIsize / 2)
            
            for ROI in ROIs:

                #check that ROI is not touching borders:
                try:
                    aid.crop(bitmap, ROI)
                except IndexError:
                    print('ROI touches image borders, skipping')
                    continue
                snip = aid.crop(bitmap, ROI)
                #fits 1, 2 or 3 gauss spots and determines which one is best
                #returns optimised parameter array
                bestfit, twoIstar, _ = GAP.fitNGauss (
                    snip, options, 
                    DTwoIstar = DTwoIstar, garbageBrightness = garbageBrightness,
                    junkIstar = junkIstar, verbose = verbose, outdir = outdir
                    )
                #this function appends each time it is called
                loc[color].fillSpotLst(bestfit, ROI)
            if verbose: 
                aid.plotBitmapROI(bitmap, loc[color].spotLst, \
                    title = color + ' channel:' + file)
        locLst.append(loc)
    
    if outname:
        with open(outname, 'wb') as output:
            pickle.dump(locLst, output, 1)
        
    return locLst


def analyseLoc(
    loc, cntr, Igate = [0,0,0], ltgate = [0,0,0], 
    winSigma = 3, framestop = 20, rebin = None,
    TACCal = 0.0128, ntacs = 256,
    verbose = False, pxSize = 10, bgphotons = [0,0,0]):
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
    
    #reload ptu file
    CLR = loadGRYimage(outloc['filepath'], Ggate = Igate[0], Rgate = Igate[1],
        Ygate = Igate[2], ntacs = ntacs, framestop = framestop,
        rebin = rebin)
    #override existing snips with new snips
    for color, bitmap in CLR.workIntensity.__dict__.items():
        assert color in ['G', 'R', 'Y']
        #replace if channel already exists, otherwise create
        try:
            outloc[color].bitmap = bitmap
        except KeyError:
            outloc[color] = GAP.Channel(bitmap)
        if verbose:
            plt.imshow(bitmap)
            plt.show()
    #add pairwise localisation of loc to outloc
    sortSpots(outloc)    

    cntr = calcFRETind(CLR, outloc, winSigma, cntr, verbose, 
        Igate, ltgate, pxSize, rebin, bgphotons)
    return outloc, cntr
    
def analyseLocLst(locLst, Igate = [0,0,0], ltgate = [0,0,0],
    winSigma = 3, framestop = 20, rebin = None,
    TACCal = 0.0128, ntacs = 256,
    verbose = False, pxSize = 10, outname = '', bgphotons = [0,0,0]):
    #add saving functionality
    outLst = []
    cntr = 0
    for loc in locLst:
        outloc, cntr = analyseLoc(loc, cntr, Igate = Igate, ltgate = ltgate, 
            winSigma = winSigma, framestop = framestop,
            rebin = rebin, verbose = verbose, bgphotons = bgphotons)
        print('analysing localisation %i' %cntr)
        outLst.append(outloc)
    if outname:
        with open(outname, 'wb') as output:
            pickle.dump(outLst, output, 1)
    return outLst
    
    
def loadGRYimage(
    ffile, Ggate = 0, Rgate = 0, Ygate = 0, 
    uselines = np.array([1,2]), ntacs = 256, framestop = -1, rebin = None,
    gateStop = 150):
    #load image, gate
    CLR = IM.processLifetimeImage(
        ffile.encode(), uselines = uselines, ntacs = ntacs,
        framestop = framestop)
    CLR.loadLifetime()
    if rebin:
        CLR.rebin(rebin, rebin)
    CLR.gate(Ggate, gateStop, channel = 'G')
    CLR.gate(Rgate, gateStop, channel = 'R')
    CLR.gate(Ygate, gateStop, channel = 'Y')
    CLR.loadIntensity()
    return CLR
    

def kickvector(v, maxval):
    """
    filter all entries in v whose norm is larger than maxval
    v must be 2D array
    """
    return v[np.linalg.norm(v, axis = 1) < maxval]
    
def filterVec(v, maxval = 50, verbose = True, center = None):
    """kick outlier and center coordinates of LocLst"""
    #issue: give maxval a more descriptive name. max_dist?
    v = kickvector(v, maxval)
    if not center: center = np.mean(v, axis = 0)
    if verbose:
        print('%.1f in x and %.1f in y subtracted' % (center[0], center [1]))
    return v - center
    
def plotSinglePair(locLst, pxSize = 10):
    #split into fitting function and filtering function
    #filtering function goes into different function
    posx = np.zeros(len(locLst))
    posy = np.zeros(len(locLst))
    for i, loc in enumerate(locLst):
        posx[i] = loc['G'].spotLst[0].posx - loc['Y'].spotLst[0].posx
        posy[i] = loc['G'].spotLst[0].posy - loc['Y'].spotLst[0].posy
    posx *= pxSize
    posy *= pxSize
    filterVec(np.array([posx, posy]).transpose(), 50)
    plt.plot(posx, posy, '.')
    plt.plot([0,0], [-100,+100], 'r--')
    plt.plot([-100,+100], [0,0], 'r--')
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlim([-50, 50])
    plt.ylim([-50, 50])
    stdx = np.std(posx)
    stdy = np.std(posy)
    print('standard deviation of spots is %.2f and %.2f' % (stdx, stdy))
    plt.show()
    return posx, posy

def plotOccurence(locLst):
    """plots a 2D histogram of how many spots have been fitted
    in the Green and red channel.
    Takes as argument a spotLst Lst"""
    #find max occurence
    Gmax = 0
    Ymax = 0
    for loc in locLst:
        Gspots = len(loc['G'].spotLst)
        Yspots = len(loc['Y'].spotLst)
        if Gspots + 1> Gmax: Gmax = Gspots + 1#+1 to store zero entries
        if Yspots + 1 > Ymax: Ymax = Yspots + 1
    occurence = np.zeros([Gmax,Ymax], np.int)
    for loc in locLst:
        Gspots = len(loc['G'].spotLst)
        Yspots = len(loc['Y'].spotLst)
        occurence[Gspots, Yspots] += 1
    plt.figure(figsize = [Ymax,Gmax])
    x, y = np.meshgrid(np.arange(Ymax),np.arange(Gmax))
    plt.scatter(x,y, s = occurence)
    ax = plt.gca()
    ax.set_xticks(np.arange(Ymax))
    ax.set_yticks(np.arange(Gmax))
    plt.ylabel ('# Green spots identified')
    plt.xlabel ('# red spots identified')
    plt.show()
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
        p.add_many(('mu%i' % i, mu, True, 0, np.amax(y)),
                   ('A%i' % i, A, True, 0, A_max),
                   ('sig%i' % i, sig, True, 0, np.inf)
                   )
    return p

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
    Gprecision = np.sqrt(pF.findVar([0, 0, Gsigma, np.sqrt(Gbg), Gphotons], pxSize))
    Yprecision = np.sqrt(pF.findVar([0,0, Ysigma, np.sqrt(Ybg), Yphotons], pxSize))
    chiSigma = np.sqrt(Gprecision**2 + Yprecision**2 + posprecision**2)
    return chiSigma

def fitNChidistr(p, bins, counts):
    """fits a chi distribution to a set of distances
    params are ordered: mu, sigma, amplitude, offset"""
    #possible to make function generic to fit any distribution?
    #even better: create object that contains the data, binwidth, maxbins etc.
    fitres = lmfit.minimize(get_logLikelihood1DPoisson, p, method = 'nelder',
        args = (NncChidistr, bins, counts, -1))
    logLikelihood = get_logLikelihood1DPoisson(fitres.params, NncChidistr, 
                                                bins, counts)
    samplesize = np.sum(counts)
    AIC = get_AIC(fitres.nvarys, logLikelihood)
    AICc = get_AICc(fitres.nvarys, logLikelihood, samplesize)
    BIC = get_BIC(fitres.nvarys, logLikelihood, samplesize)
    return fitres, AIC, AICc, BIC, logLikelihood
    
def scanLikelihoodSurface(param_ranges, p, bins, counts, verbose = False):
    param_ticks=[]
    param_names=param_ranges.keys()
    for param_name in param_names:
        param_ticks.append(param_ranges[param_name])
        
    nDshape = [len(sublist) for sublist in param_ticks]
    length = np.product(nDshape)
    logLikelihoodSurface = np.zeros(length)

    for i, point in enumerate(product(*param_ticks)):
        for j, name in enumerate(param_names):
            p[name].set(value = point[j], vary = False)
        try:
            *_, logLikelihood = fitNChidistr(p, bins, counts)
        except ValueError:
            if verbose: print ('fit failed for ' + str(point))
            logLikelihood = np.nan
        logLikelihoodSurface[i] = logLikelihood
        if(i % 10 == 0):
            if verbose: print ('calculating point ' + str(point))
    return logLikelihoodSurface.reshape(nDshape)
    
def plotLikelihoodSurface(surface, param_ranges, skip = 2, outname = '', 
    title = ''):
    """utility function, works only in 2D"""
    keyy, keyx = param_ranges.keys()
    bestfit = np.nanmax(surface.flatten())
    contourlevels = [bestfit -2, bestfit - 1, bestfit-0.5, bestfit]
    fig, ax = plt.subplots(1,1)
    
    #img = ax.imshow(surface, cmap = 'hot')
    vmin = bestfit-5
    plt.imshow(surface, cmap = 'hot', vmin = vmin, vmax = bestfit)
    #plt.imshow(surface, cmap = 'hot')
    CS = ax.contour(surface, levels = contourlevels)
    labels = ['13.5% as likely','36% as likely', '60% as likely', 'most likely']
    fmt = {}
    for l, s in zip(CS.levels, labels):
        fmt[l] = s
    ax.clabel(CS, fmt = fmt, fontsize = 6)
    x_label_list = param_ranges[keyx][::skip]
    ax.set_xticks(np.arange(len(param_ranges[keyx]))[::skip])
    ax.set_xticklabels(['%.1f' % x for x in x_label_list])
    ax.set_xlabel(keyx)
    y_label_list = param_ranges[keyy][::skip]
    ax.set_yticks(np.arange(len(param_ranges[keyy]))[::skip])
    ax.set_yticklabels(['%.1f' % x for x in y_label_list])
    ax.set_ylabel(keyy)
    plt.title(title)
    plt.text
    plt.colorbar()
    if outname:
        plt.savefig(outname, dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    
def whichChiIsBest(bins, counts, verbose = False):
    """utility function"""
    binwidth = bins[1]-bins[0]
    fitresL = []
    AICL = []
    BICL = []
    for N in range(4):
        try:
            p = genPeakEst(N, bins, counts)
            fitres, AIC, _, BIC, _ = fitNChidistr(p, bins, counts)
        except ValueError:
            print('fit for %i peaks failed' %N)
            break
        fitresL.append(fitres)
        AICL.append(AIC)
        BICL.append(BIC)
    
    bestfit = np.argmin(AICL)
    if verbose:
        for i, AIC in enumerate(AICL):
            print('AIC is %.1f for %i peaks' % (AIC, i))
        for i, BIC in enumerate(BICL):
            print('BIC is %.1f for %i peaks' % (BIC, i))
    return fitresL[bestfit]
    
def plotdistr(dist, bins, fit = None, title = '', modelout = '', 
    plotout = ''):
    """plots a histogrammed disttribution with bin position bins and 
    counts in each bin. fit is an lmfit.MinimizerResult object, if it is given,
    a model is plotted too.
    plot and data export modalities are integrated."""
    binwidth = bins[1] - bins[0]
    if fit:
        x = np.arange(0,max(bins),.1)
        p = fit.params
        model = NncChidistr(x, p)
        plt.plot(x, model, 'k')
        i = 0
        while True:
            try:
                model = ncChidistr(x, p['mu%i' % i], p['sig%i' % i], p['A%i' % i], 0)
                plt.plot(x,model, '--')
                i+=1
            except KeyError:
                bg = np.ones(len(x))*p['bg']
                plt.plot(x, bg, '--')
                break

        
    plt.hist(dist, bins = bins, color = 'c')
    plt.xlabel('distance (nm)')
    plt.ylabel('localisation events / %.0f nm' % binwidth)
    plt.title(title)
    if plotout:
        plt.savefig(plotout, dpi = 300, bbox_inches = 'tight')
    if modelout:
        dictout = {}
        dictout['grid'] = x
        dictout['model'] = model
        saveDict(dictout, modelout)
        fpath = os.path.splitext(modelout)[0] + '_fit_parameters.txt'
        with open(fpath, 'wt') as f:
            f.write(lmfit.fit_report(fit))
    plt.show()

def exportChiComponents(outname, bins, p_all):
    """utility function"""
    dataFrame = pd.DataFrame({'xgrid':bins})
    for i in range(10):#no more than ten populations are exported
        try:
            p = lmfit.Parameters()
            p.add('bg', 0, False, 0, 10)
            p['mu0'] = p_all['mu' + str(i)]
            p['sig0'] = p_all['sig' + str(i)]
            p['A0'] = p_all['A' + str(i)]
            model = NncChidistr(bins, p)
            dataFrame['ncChi_pop'+str(i)] = model
            print(p)
            #gen fit and export
        except KeyError:
            break
    dataFrame.to_csv(outname)
    return dataFrame

def fitChiDistr(dist, params0, bounds = None, outfile = '', binwidth = 2, 
                maxbin = 60, plotshow = True):
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
    NG = len(loc['G'].spotLst)
    NY = len(loc['Y'].spotLst)
    unsortedG = [x for x in range(NG)]
    unsortedY = [x for x in range(NY)]
    sortedG = []
    sortedY = []
    Ndist = min (NG, NY)
    Gcoords = np.zeros([NG, 2])
    Ycoords = np.zeros([NY, 2])
    alldist = np.zeros([NG, NY])

    for i, el in enumerate(loc['G'].spotLst):
        Gcoords[i] = np.array([el.posx, el.posy])
    for i, el in enumerate(loc['Y'].spotLst):
        Ycoords[i] = np.array([el.posx, el.posy])
    
    for i in range(NG):
        for j in range(NY):
            alldist[i, j] = np.linalg.norm(Gcoords[i] - Ycoords[j])

    for i in range(Ndist):
        Gpos, Ypos = [alldist.argmin() // NY , alldist.argmin() % NY]
        #do bookkeeping
        sortedG.append(Gpos)
        unsortedG.remove(Gpos)
        sortedY.append(Ypos)
        unsortedY.remove(Ypos)
        alldist[Gpos] = 1e6
        alldist[:, Ypos] = 1e6
    #unpaired spots are appended
    sortedG += unsortedG
    sortedY += unsortedY
    #re-order
    loc['G'].spotLst[:] = [loc['G'].spotLst[i] for i in sortedG]
    loc['Y'].spotLst[:] = [loc['Y'].spotLst[i] for i in sortedY]


#def cropSpot(xcenter, ycenter, bitmap, winSigma):
#    """crops 2D or 3D data in 2D
#    spotcenter defines the coordinate of the spot in pixels in bitmap
#    The total width of cropwindow is 2*winSigma + 1
#    """
#    #convert between float peak to pixel position
#    xcenter = np.round(xcenter).astype(np.int)
#    ycenter = np.round(ycenter).astype(np.int)
#    xstart, ystart, xstop, ystop = [xcenter - winSigma,
#                                  ycenter-winSigma,
#                                   xcenter + winSigma,
#                                   ycenter + winSigma]
#    return bitmap[xstart : xstop + 1, ystart : ystop + 1]
    
def fitTau(TAC, TACCal = 0.128, verbose = False, params0 = [1, 2], bgphotons = 0):
        """simple fitting and plotting function that fits 
        a monoexponential decay to a TAC decay using MLE optimization
        No boundaries are implemented.
        """
        tactimes = np.arange(len(TAC)) * TACCal
        bgArrivalTime = len(TAC) * TACCal / 2
        # add TACcal/2 because average arrival time in first bin is TACcal / 2
        #meantac = np.sum(tactimes * TAC) / np.sum(TAC) + TACCal / 2 
        bgcorTAC = (np.sum(tactimes * TAC) - bgphotons * bgArrivalTime) / \
                    (np.sum(TAC) - bgphotons)
        if verbose:
            plt.plot(tactimes, TAC)
            plt.plot(tactimes, expDecay(tactimes, bgcorTAC, \
                    np.sum(TAC) / bgcorTAC * TACCal))
            plt.plot(tactimes, np.ones(len(TAC)) * bgphotons / len(TAC))
            plt.show()
        return bgcorTAC

def calcFRETind(CLR, loc, winSigma, cntr, verbose, Igate, ltgate, pxSize,
                rebin, bgphotons = [0,0,0]):
    """calc FRET indicators based on intensity and lifetime information
    Takes a loc dict and edits properties of the 'FRETind' entry
    2 winSigma + 1 is the siye of the integration area"""
    
    npairs = min(len(loc['G'].spotLst), len(loc['Y'].spotLst))
    loc['FRETind'] = []
    #reload ungated data for lifetime information
    CLR.loadLifetime()
    if rebin:
        CLR.rebin(rebin, rebin)
    for i in range(npairs):
        loc['FRETind'].append(FRETind())
    for i in range(npairs):
        loc['FRETind'][i].no = cntr
        GROI = aid.pos2ROI(loc['G'].spotLst[i].posx,
                         loc['G'].spotLst[i].posy,
                         winSigma)
        YROI = aid.pos2ROI(loc['Y'].spotLst[i].posx,
                         loc['Y'].spotLst[i].posy,
                         winSigma)
        RROI = YROI         #get No. of Red photons using Y localisation
        try:
            Gsnip = aid.crop(loc['G'].bitmap, GROI)
            Rsnip = aid.crop(loc['R'].bitmap, RROI)
            Ysnip = aid.crop(loc['Y'].bitmap, YROI)
        except IndexError:
            print('ROI touches image borders, skipping')
            continue
        GlifetimeSnip = aid.crop(CLR.workLifetime.G, GROI)
        RlifetimeSnip = aid.crop(CLR.workLifetime.R, RROI)
        YlifetimeSnip = aid.crop(CLR.workLifetime.Y, YROI)
        Gsnip = np.sum(GlifetimeSnip[:,:,Igate[0]:150], axis = 2)
        Rsnip = np.sum(RlifetimeSnip[:,:,Igate[1]:150], axis = 2)
        Ysnip = np.sum(YlifetimeSnip[:,:,Igate[2]:150], axis = 2)
        #intensity based indicators
        loc['FRETind'][i].Gbg = loc['G'].spotLst[i].bg
        loc['FRETind'][i].Ybg = loc['Y'].spotLst[i].bg
        Gphotons = np.sum(GlifetimeSnip[:,:,Igate[0]:150])
        Rphotons = np.sum(RlifetimeSnip[:,:,Igate[1]:150])
        Yphotons = np.sum(YlifetimeSnip[:,:,Igate[2]:150])
        loc['FRETind'][i].NG = Gphotons
        loc['FRETind'][i].NR = Rphotons
        loc['FRETind'][i].NY = Yphotons
        loc['FRETind'][i].proxRatio = Rphotons / (Gphotons + Rphotons)
        loc['FRETind'][i].stoichiometry = \
            (Gphotons + Rphotons) / (Gphotons + Rphotons + Yphotons)
        loc['FRETind'][i].GTAC = np.sum(GlifetimeSnip, axis = (0,1))
        loc['FRETind'][i].RTAC = np.sum(RlifetimeSnip, axis = (0,1))
        loc['FRETind'][i].YTAC = np.sum(YlifetimeSnip, axis = (0,1))
        loc['FRETind'][i].tauG = fitTau(loc['FRETind'][i].GTAC[ltgate[0]:150], \
                                        0.128, verbose, bgphotons = bgphotons[0])
        loc['FRETind'][i].tauR = fitTau(loc['FRETind'][i].RTAC[ltgate[1]:150], \
                                        0.128, verbose, bgphotons = bgphotons[1])
        loc['FRETind'][i].tauY = fitTau(loc['FRETind'][i].YTAC[ltgate[2]:150], \
                                        0.128, verbose, bgphotons = bgphotons[2])
        distx = (loc['G'].spotLst[i].posx - loc['Y'].spotLst[i].posx )* pxSize
        disty = (loc['G'].spotLst[i].posy - loc['Y'].spotLst[i].posy )* pxSize
        loc['FRETind'][i].distx = distx
        loc['FRETind'][i].disty = disty
        loc['FRETind'][i].dist = np.linalg.norm([distx, disty])
        if verbose: aid.plotBitmapROI( loc['G'].bitmap, loc['G'].spotLst)
        if verbose:
            fig, (ax1, ax2, ax3) = plt.subplots(1,3)
            ax1.imshow(Gsnip)
            ax1.set_title('green channel')
            ax2.imshow(Rsnip)
            ax2.set_title('red channel')
            ax3.imshow(Ysnip)
            ax3.set_title('Yellow channel')
            plt.show()
            

        cntr += 1
    return cntr

def getFRETnames(locLst):
    """aider function that returns the names of FRET indicators"""
    names = []
    i = 0
    #return empty list if FRETind key does not exist
    try:
        locLst[i]['FRETind']
    except KeyError:
        return []
        
    while len(locLst[i]['FRETind']) == 0:
        i += 1
        if i == 10000:
            raise RuntimeError
    
    for name, value in locLst[i]['FRETind'][0].__dict__.items():
        if type(value) != np.ndarray and name[0] != '_':
            names.append(name)
    return names
def genSpotNames(locLst):
    names = []
    i = 0
    while len(locLst[i]['G'].spotLst) == 0:
        i += 1
        if i == 100000:
            raise RuntimeError
    for name, value in locLst[i]['G'].spotLst[0].__dict__.items():
        #for color in ['G', 'R', 'Y']:
        for color in ['G', 'Y']:
            names.append(name + color)
    return names
    
def genStats(locLst, outfile = ''):
    """
    saves al FRET indicators of locLst in a Margarita-readable file.
    If there are more spots in a single localisation, they are treated as 
    independant.
    If unpaired spots exist in a localisation, all unreachable parameters are
    set to 0.
    If the locLst has paired Green and Yellow entries, each line represents a pair
    If it is unordered, columns are also unordered.
    """
    #issue:fnames should also be included
    statsdict = {}
    spotNames = genSpotNames(locLst)
    #FRETnames contains only single variable FRET indicators
    #empty when locLst['FRETind'] does not exist
    FRETnames = getFRETnames(locLst)
    #ROInames = ['ROI_xstart', 'ROI_ystart', 'ROI_xstop', 'ROI_ystop']
    names = spotNames + FRETnames #+ ROInames# + ['filepath']
    for name in names:
        statsdict[name] = []
    
    for loc in locLst:
        maxdyes = max(len(loc['G'].spotLst),len(loc['Y'].spotLst))
        for i in range(maxdyes):
            for spotName in spotNames:
                attr = spotName[:-1]
                color = spotName[-1]
                #need to take care of all the possible errors when object does not exist
                try: attr = getattr(loc[color].spotLst[i], attr)
                except: attr = 0
                statsdict[spotName].append(attr)
            for FRETname in FRETnames:
                try: attr = getattr(loc['FRETind'][i], FRETname)
                except: attr = 0
                statsdict[FRETname].append( attr )
                    
    if outfile:
        print('saving spectroscopic parameters to disc for Margarita')
        statsDataFrame = pd.DataFrame(statsdict)
        statsDataFrame.to_csv(outfile, sep = '\t', float_format = '%.3f')
        
    return statsdict
    
def tryGetattr(obj, attrName):
    try: 
        attr = getattr(obj, attrName)
    except AttributeError:
        attr = 0
    return attr
    
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
def filterStats(stats, indicator, vmin, vmax):
    """
    stats must be dict with equal length entries
    select rows of stats[indicator] that are between vmin and vmax
    """
    #ISSUE: pop utility deletes original list. Unexpected behaviour.
    #loop from last to first
    for i in range(len(stats[indicator]) -1, -1, -1):
        if stats[indicator][i] < vmin or stats[indicator][i] > vmax:
            for name, value in stats.items():
                stats[name].pop(i)
    return stats
def saveDict(dictionary, outfile):
    """save dict column-wise to file."""
    keys = dictionary.keys()
    header = ''
    for key in keys:
        header += key + '\t'
    header += '\n'
    length = len(dictionary[key])
    with open(outfile, 'wt') as f:
        f.write(header)
        for i in range(length):
            line = ''
            for key in keys:
                line += str(dictionary[key][i]) + '\t'
            line += '\n'
            f.write(line)
            
def subensembleTAC(locLst, ntacs = None, outfile = None):
    """sum all TAC decays to do subensemblefitting"""
    #get length of tac array
    if not ntacs:
        for loc in locLst:
            try:
                ntacs = len(loc['FRETind'][0].GTAC)
                break
            except:
                continue
    eGTAC = np.zeros(ntacs)
    eRTAC = np.zeros(ntacs)
    eYTAC = np.zeros(ntacs)
    dummy_IRF = np.zeros(2 * ntacs)
    
    for loc in locLst:
        for spot in loc['FRETind']:
            eGTAC += spot.GTAC
            eRTAC += spot.RTAC
            eYTAC += spot.YTAC
   
    eGTAC = np.concatenate((eGTAC, np.zeros(ntacs)))
    eRTAC = np.concatenate((eRTAC, np.zeros(ntacs)))
    eYTAC = np.concatenate((eYTAC, np.zeros(ntacs)))
    dummy_IRF[0] = 1
    if outfile:
        outG = outfile[:-4] + '_eGTAC.txt'
        outR = outfile[:-4] + '_eRTAC.txt'
        outY = outfile[:-4] + '_eYTAC.txt'
        outIRF = outfile[:-4] + '_dummy_IRF.txt'
        np.savetxt(outG, eGTAC, fmt = '%i', delimiter = '\t')
        np.savetxt(outR, eRTAC, fmt = '%i', delimiter = '\t')
        np.savetxt(outY, eYTAC, fmt = '%i', delimiter = '\t')
        np.savetxt(outIRF, dummy_IRF, fmt = '%i', delimiter = '\t')
    return eGTAC, eRTAC, eYTAC
            
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
        
    
    headerbaseG = 'image_ID_cI\tPixel Number_cI\t'
    headerbaseRY = 'image_ID_cII\tPixel Number_cII\t'

    i = 0
    while len(locLst[i]['G'].spotLst) == 0:
        i+=1
        if i> 100000:
            break
    for key in locLst[i]['G'].spotLst[0].__dict__.keys():
        headerbaseG += key + '_cI\t'
    for key in locLst[i]['G'].spotLst[0].__dict__.keys():
        headerbaseRY += key + '_cII\t'
    FRETnames = getFRETnames(locLst)
    headerFRET = ''
    for name in FRETnames:
        headerFRET += name + '\t'
    
    for loc in locLst:
        locname = os.path.split(loc['filepath'])[1][:-4]
        for color in ['G', 'R', 'Y']:
            if color == 'G':
                prepend = '_Green '
                append = '_cI'
                ext = '.0I4'
            elif color == 'R':
                prepend = '_Red '
                append = '_cII'
                ext = '.II4'
            elif color == 'Y':
                prepend = '_Yellow '
                append = '_cII'
                ext = '.II4'
            
            spotLst = loc[color].spotLst
            #naming convention from Suren
            outfile = os.path.join(outfolder, locname + prepend + \
                    'Photons_Fit results' + ext)
            if color == 'G':
                header = headerbaseG + headerFRET
            elif color in ['R', 'Y']:
                header = headerbaseRY
            header += '\n'
            
            with open(outfile, 'wt') as f:
                f.write(header)
                for i, spot in enumerate(spotLst):
                    posx = np.round(spot.posx * rebin)
                    posy = np.round(spot.posy * rebin)
                    pixid = posy * sizex + posx +1
                    line = '0\t%d' % pixid
                    for name, value in spot.__dict__.items():
                        line += "\t%.3f" % value
                    #add FRET observables to G file
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
