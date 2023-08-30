import aid_functions as aid
import batchplot as bp
import os
import pandas as pd
import ImageManipulation as IM
import numpy as np
import fitDA
import gc
import tiffile #pip install tiffile if missing
import copy
try:
    import vvvh_support
except ModuleNotFoundError:
    import sys
    sys.path.append(r'K:\vanderVoortN\redRobin\src')
    import vvvh_support

import warnings
from scipy.ndimage import gaussian_filter
#note name df is blocked for dataframe
debug = False
if debug:
    import matplotlib.pyplot as plt

warnings.simplefilter("default")



def saveNpAsImage(array, outname):
    image = Image.fromarray(array)
    image.save(outname)

def createSeriesHiLoMasks(seriesdir):
    """assumes an existing structure of masks and cell images.
    Finds each pair of mask and cellimg from pre-existing file structure
    calls function that works to create additional masks
    Saves these additional masks to disc as .tiff"""
    #get the names of the tiff 
    cellimgdir = 'imagesForMasking'
    cellimgfilesp = os.listdir(os.path.join(seriesdir, cellimgdir))
    cellimgfiles = [file for file in cellimgfilesp if file.endswith('imG.tif')]
    
    #find the Masks folders
    for entry in os.listdir(seriesdir):
        if entry.endswith('_Masks'):
            #get the tiff file that works with this folder
            basename = entry[:-6]
            for cellimgfile in cellimgfiles:
                if (cellimgfile[:-7] == basename):
                    fcellname = os.path.join(seriesdir, cellimgdir, cellimgfile)
                    cellarr = np.array(Image.open(fcellname))
                    #there is a legacy bug that where some values are 
                    # truncated in 8 bit tiffs: check for clipping of data
                    maxval = max(cellarr.flatten())
                    if maxval == 255 or maxval == 2**16-1:
                        with warnings.catch_warnings():
                            warnings.simplefilter("always")
                            warnings.warn('saturation detected for %s' % cellimgfile, )
                    break
            #get the masks 
            maskfiles = os.listdir(os.path.join(seriesdir, entry))
            maskfiles = [file for file in maskfiles if file.endswith('.tif')]
            maskfiles = [file for file in maskfiles if 
                         not (file.endswith('_lo.tif') or file.endswith('_hi.tif'))]
            maskffiles = [os.path.join(seriesdir, entry, maskfile) 
                          for maskfile in maskfiles]
            for maskffile in maskffiles:
                maskarr = np.array(Image.open(maskffile))
                #make sure the cell area is 1, other values are 0
                maskarr = maskarr == 0
                #function that does the calculation
                himask, lomask = createHiLoMasks(cellarr, maskarr)
                hiOutname = maskffile[:-4] + '_hi.tif'
                loOutname = maskffile[:-4] + '_lo.tif'
                saveNpAsImage(himask, hiOutname)
                saveNpAsImage(lomask, loOutname)
                
def deleteHiLoMasks(seriesdir):
    """in '_Masks' subfolders delete all .tiff files ending on '_hi.tif' or 
    '_lo.tif'.
    """
    #find the Masks folders
    for entry in os.listdir(seriesdir):
        fsubdir = os.path.join(seriesdir, entry)
        if entry.endswith('_Masks'):
            for file in os.listdir(fsubdir):
                if file.endswith('_hi.tif') or file.endswith('_lo.tif'):
                    print('deleting mask %s' % file)
                    ffile = os.path.join(fsubdir, file)
                    os.remove(ffile)


def createHiLoMasks(cellarr, maskarr, verbose = False):
    #smooth against shot noise
    cellarr = gaussian_filter(cellarr, 1, output = float)
    #apply mask
    maskedcell = cellarr * maskarr
    total = np.sum(maskedcell)
    #make an ordering of all smoothed values
    sortedcell = np.sort(maskedcell.flatten())
    cumsum = np.cumsum(sortedcell)
    spliti = np.argmax(cumsum > total / 2)
    #get the threshold value
    splitval = sortedcell[spliti]
    lomask = np.logical_and(maskedcell > 0, maskedcell < splitval)
    himask = maskedcell >= splitval
    if verbose:
        #checks
        lowersum = np.sum(maskedcell * lomask)
        uppersum = np.sum(maskedcell * himask)
        print('lower is up to %.2f\nlowersum is %i\nuppersum is %i' %
              (splitval, lowersum, uppersum))
        print('total sum check: %r' % (np.isclose(lowersum + uppersum, total)))
    return himask, lomask

################################################################################
#This class structure is disadvantageous because everytime I test a bugfix,
#I need to reanalyze all data to check it, considering the files are very large
#this takes a long time, solution: make functional

def PS2PandS(TACPS):
    assert len(TACPS) %2 == 0 and len(TACPS.shape) == 1 ,\
        "array not divisible by two or not 1D."
    ntacs = len(TACPS) / 2
    TACP = TACPS[:ntacs]
    TACS = TACPS[ntacs:]
    return TACP, TACS

def PandS2PS(TACP, TACS):
    assert TACP.shape == TACS.shape and len(TACS.shape) ==1, \
        "arrays dimensions mismatch or are not 1D."
    ntacs = len(TACP)
    TACPS = np.zeros(ntacs*2)
    TACPS[:ntacs] = TACP
    TACPS[ntacs:] = TACS
    return TACPS
def calculateDerivedVariables(df, 
                              integrationtime = 1, 
                              Gpower = 1, 
                              Ypower = 1):
    """integration time is dwelltime * Nframes"""
    print(df)
    df['surfaceMax'] = df[['surfaceG', 'surfaceY']].max(axis = 1)
    for label, power in zip(['G', 'R', 'Y'], [Gpower, Ypower, Ypower]):
        df['Br' + label] = df['N' + label+'-tot'] / df['surfaceMax']
        df['rate'+label] = df['Br' + label] / integrationtime
        #LPC stands for LaserPowerCorrected
        df['rate' + label + '_LPC'] = df['rate'+label] / power
    df['rateGrateY'] = df['rateG'] / df['rateY']
    df['rateGrateR'] = df['rateG'] / df['rateR']
    df['rateRrateY'] = df['rateR'] / df['rateY']
    df['rateGrateY_LPC'] = \
        df['rateG_LPC'] / df['rateY_LPC']
    df['rateGrateR_LPC'] = \
        df['rateG_LPC'] / df['rateR_LPC']
    df['rateRrateY_LPC'] = \
        df['rateR_LPC'] / df['rateY_LPC']
    return df

def cleanImage(image):
    del(image.baseLifetime)
    del(image.workLifetime)
    return 1

def saveTACs(image, TACdir, label):
    TACPS = PandS2PS(image.P, image.S)
    TACout = os.path.join(TACdir, image.name + label +'_PS.dat')
    VMout = os.path.join(TACdir, image.name  + label + '_VM.dat')
    rout = os.path.join(TACdir, image.name  + label +'_r.dat')
    np.savetxt(TACout, TACPS, fmt = '%i')
    np.savetxt(VMout, image.VM, fmt = '%i')
    np.savetxt(rout, image.r, fmt = '%.5e')
    return 1

def genDefaultFitKwargs():
    return { #some dummy mono exponential
            'D0dat' : np.exp(-np.arange(0,25, 0.064) / 2.5),
            'decaytype' : 'VM',
            'fitrange' : (30, 380)}
            
def getMask(fname):
    mask = tiffile.imread(fname)
    #we need masks to be binary
    mask [ mask != 0 ] = 1
    return mask #mask is np array
    
def trySaveToCSV(dfrm, outname):
    try:
        dfrm.to_csv(outname)
    except PermissionError:
        print('could not save stats, is the file open in Excel?')

## having the below two functions functional, such that I don't need to re-
# analyze all data
def wrap_OO_fitDO(SampleSet, identifier, data_irf, data_af, af_G_norm, tauxD0,
                    bsel = None):
    #get vvvh format decays
    decays_vv = np.array(SampleSet.getDecay(decaytype = 'P'))
    decays_vh = np.array(SampleSet.getDecay(decaytype = 'S'))
    Ndecays = decays_vv.shape[0]
    length = decays_vv.shape[1]
    decays_vvvh = np.zeros([Ndecays, 2, length])
    decays_vvvh[:,0, :] = decays_vv
    decays_vvvh[:,1, :] = decays_vh
    #path for saving result
    csvout = os.path.join(SampleSet.resdir, identifier + '_2lt_OOFitData.csv')
    dfrm = vvvh_support.fitDO(decays_vvvh, data_irf, data_af, 
                        SampleSet.imstats, af_G_norm, tauxD0, bsel = None,
                        csvout = csvout)
    #add as class variable
    SampleSet.D0FitOODfrm = dfrm
    
def wrap_OO_fitDA(SampleSet, identifier, data_irf, data_af, af_G_norm, tauxD0,
                    bsel = None):
    #get vvvh format decays
    decays_vv = np.array(SampleSet.getDecay(decaytype = 'P'))
    decays_vh = np.array(SampleSet.getDecay(decaytype = 'S'))
    Ndecays = decays_vv.shape[0]
    length = decays_vv.shape[1]
    decays_vvvh = np.zeros([Ndecays, 2, length])
    decays_vvvh[:, 0, :] = decays_vv
    decays_vvvh[:, 1, :] = decays_vh
    #path for saving result
    csvout = os.path.join(SampleSet.resdir, identifier + '_D0DA2lt_OOFitData.csv')
    dfrm = vvvh_support.fitDA(decays_vvvh, data_irf, data_af, 
                        SampleSet.imstats, af_G_norm, tauxD0, bsel = None,
                        csvout = csvout)
    #add as class variable
    SampleSet.DAFitOODfrm = dfrm
        
class sampleSet():
    """This class is intended to automate image analysis for cellular data to
        avoid time-consuming manual work in AnI.
        To work, this script neads a functioncal copy of Seidel in the pythonpath
        """
    #set to False after default settings have been applied once
    isApplyDefaultSettings = True 
    def __init__(self,
                 wdir,
                 **settings
                 ):
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
        #issue: if I want to chance one imreadkwarg, have to specify all of them
        #potential workaround is to init a dict-like class, but this is cumbersome
        #put all this in a separate function, so that it can be updated separately.
        self.setDefaultSettings(wdir)
        self.setUserSettings(**settings)
        return
        
    def setDefaultSettings(self, wdir):
        self.wdir = wdir
        self.imreadkwargs =  {'ntacs' : 1024,
                    'pulsetime' : 50,
                    'dwelltime': 20e-6,
                    'TAC_range': 4096}
        self.TACdir = os.path.join(wdir, 'TAC')
        self.resdir = os.path.join(wdir, 'results')
        self.imdir = os.path.join(wdir, 'images')
        self.images = {'G': [], 'Y': []}#dict entries are channels,
        self.g_factor = 1
        self.dt_glob = 0.064
        self.Gpower = 1
        self.Ypower = 1
        self.FRETPIETACranges = [[0,380], [0, 380], [380,800]]
        self.Nframes = -1
        self.dataselect = (0, None)
        self.PSchannels = [[1,0], [5,4]]
        self.uselinesLst = [np.array([1]), np.array([1])]
        self.PSshift = 0
    
    def setUserSettings(self, **settings):
        for setting, settingvalue in zip(settings, settings.values()):
            setattr(self, setting, settingvalue)
        self.completeSetting()
    def completeSetting(self):
        aid.trymkdir(self.TACdir)
        aid.trymkdir(self.resdir)
        aid.trymkdir(self.imdir)
        self.ptufiles = bp.appendOnPattern(self.wdir, 'ptu')\
            [self.dataselect[0]: self.dataselect[1]]

    def normTACs(self, normdecay,
                 decaymode = 'VM',
                 channel = 'G',
                 bgrange = [0, 5],
                 normshift = 0):
        #consider making normbydecay arguments kwargs
        assert decaymode in ['VM', 'P', 'S'], "decaymode must be VM, P or S"
        assert channel in ['G', 'Y'], "channel must be G or Y"
        for GRYimage in self.images[channel]:
            #normbydecay takes list as first argument, inconvenient
            GRYimage.VMnorm = bp.normbydecay(\
                                              [getattr(GRYimage, decaymode)],
                                              normdecay,
                                              normshift,
                                              bgrange = bgrange)[0]

    def analyzeDir(self, identifier, **kwargs):
        """analyzes all ptu files in a directory into GRY image objects

        kwargs:
        isSave:         if True, TAC decays and tiff images are stored to
                        self.TACdir and self.imdir
        isCleanImage:   if True, 3D lifetime arrays are discarded, freeing
                        up space.
        threshold:      all pixels with an intensity value lower than
                        threshold are set to 0.
        """
        for ptufile in self.ptufiles:
            self.analyzeFile(ptufile, **kwargs)
        self.genImstatsdf(identifier, **kwargs)
        return 1

    def analyzeDirwMasks(self, identifier, maskdirs, timeList, **kwargs):
        """analyzes a set of N files each with M masks, totalling
        NxM image objects
        
        maskdirs is the directory containing the masks. It can have values:
            'automatic':    Each ptufile gets an own set of masks
                            the names are automatically inferred from the names
                            of the ptufiles. E.g. 
                            cell1Masks, cell2Masks 
                            for cell1.ptu, cell2.ptu
            other strings:  One set of masks is applied to all files. 
                            The string is the directory containing the masks.
                            
        time in seconds
        
        See analyzeDir for list of kwargs"""
        #passing identifier each time is pretty cumbersome, 
        #make it a class property?
        if maskdirs == 'automatic':
            maskdirs = [os.path.join( self.wdir, ptufile[:-4] + '_Masks') 
                for ptufile in self.ptufiles]
        else:
            maskdirs = [maskdirs] * len(self.ptufiles)
        if timeList == 'automatic':
            timeList = np.zeros(len(self.ptufiles))
        for ptufile, maskdir, time in \
                zip(self.ptufiles, maskdirs, timeList):
            maskfiles = [os.path.join(maskdir, file)\
                for file in os.listdir(maskdir) if file.endswith('tif')]
            self.analyzeFilewMasks(ptufile, maskfiles, time, **kwargs)
        self.genImstatsdf(identifier, **kwargs)
        
    def analyzeFilewMasks(self, ptufile, maskfiles, time, **kwargs):
        ffile = os.path.join(self.wdir, ptufile).encode()
        for PSchannel, uselines, label in zip(self.PSchannels, \
                self.uselinesLst, ['_G', '_Y']):
            PChan, SChan = PSchannel
            image = self.loadPSImage(ffile, PChan, SChan, uselines, **kwargs)
            for maskfile in maskfiles:
                mask = getMask(maskfile)
                #there is some overhead here, but we optimize development time
                imageCopy = copy.deepcopy(image)
                #give unique name to each image, ugly
                maskid = os.path.split(os.path.splitext(maskfile)[0])[1]
                imageCopy.name += maskid
                imageCopy.maskid = maskid
                imageCopy.time = time
                self.procesPSimage(imageCopy, label, usermask = mask,
                    **kwargs)
                self.images[label[1]].append(imageCopy)
                print('finished applying mask %s' %maskid)
        
    def analyzeFile(self, ptufile, **kwargs):
        #load Donor and acceptor channels
        #P is now called G, S is called R (nuisance)
        ffile = os.path.join(self.wdir, ptufile).encode()
        for PSchannel, uselines, label in zip(self.PSchannels, \
                self.uselinesLst, ['_G', '_Y']):
            PChan, SChan = PSchannel
            image = self.loadPSImage(ffile, PChan, SChan, uselines, **kwargs)
            self.procesPSimage(image, label, **kwargs)
            self.images[label[1]].append(image)
            

    def genImstatsdf(self, identifier,
                   channels = ['G', 'Y', 'Y'],
                   #idea: consider using ['Donor', 'FRET', 'PIE'] nomenclature instead
                   labels = ['G', 'R', 'Y'], **kwargs):
        df = pd.DataFrame()
        for channel, TACrange, label in \
                zip(channels, self.FRETPIETACranges, labels):
            images = self.images[channel]
            for image in images:
                Np = np.sum(image.P[TACrange[0]:TACrange[1]])
                Ns = np.sum(image.S[TACrange[0]:TACrange[1]])
                Ntot = Np + Ns
                #need to change this into Ny-tot, Ny-p and Ny-s
                df.at[image.name, 'N'+label+'-p'] = Np
                df.at[image.name, 'N'+label+'-s'] = Ns
                df.at[image.name, 'N'+label+'-tot'] = Ntot
                #data was masked previously
                surface = np.sum(image.workIntensity.G > 0)
                df.at[image.name, 'surface'+label] = surface
                #these parameters are not dependant on the channel and thus
                #are overwritten each time. It is a bit stupid, but also 
                # a cheap mistake.
                for attr in ['time', 'maskid']:
                    try:
                        df.at[image.name, attr] = getattr(image, attr)
                    except AttributeError:
                        pass
        if self.Nframes == -1:
            print('number of frames not given,' \
                  + 'cannot calculate integration time and derived variables')
        else:
            integrationtime = self.imreadkwargs['dwelltime'] * self.Nframes
            df = calculateDerivedVariables(df, integrationtime,
                                           self.Gpower, self.Ypower)
        #save
        outname = os.path.join(self.resdir, identifier + 'imstats.csv')
        trySaveToCSV(df, outname)
        self.imstats = df
        return 0


    def procesPSimage(self, image, label, intensityThreshold = 0, isSave = True,
        isCleanImage = True, usermask = None, **kwargs):
        """only does work on image, has a lot of dependencies, unwanted
        The processLifetimeImage is not build for Anisotropy, but it can if one
            mis-uses the channels. I.e. processlifetimeImage takes up to 3
            channels labelled Green, Red, Yellow. Now we will abuse by doing:
                Green = parallel
                Red = perpendicular
                Yellow = unused
                Then repeating for both channels
        This workaround should hold for the forseeable future, but should
            ultimately be replaced."""
            
        image.loadLifetime()
        image.loadIntensity()
        #the intensity mask is determined based on workIntensity
        intMask = image.buildMaskFromIntensityThreshold(
                threshold = intensityThreshold, sumchannels = ['G', 'R'])
        #the masking happens on the 3D lifetime arrays
        image.mask(intMask, mode = 'lifetime')
        if usermask is not None:
            image.mask(usermask, mode = 'lifetime')
        #the TAC decays are determined from the worklifetime image
        self.genPSfromGRYImage(image)
        self.genDerivedFromPSDecays(image)
        #reload the intensity, such that it matches the lifetime image
        image.loadIntensity()
        if isSave:
            saveTACs(image, self.TACdir, label)
            #underlying routine is limited to 8 bit Tiff, problem
            #potential solution is to use the tiffile lib, need to test
            image.saveWorkIntensityToTiff(self.imdir, image.name + label)
        if isCleanImage: #free memory intensive 3D array
            cleanImage(image)
        gc.collect()
        return 1

    def loadPSImage(self, ffile, PChan, SChan, uselines, **kwargs):
        GRYim = IM.processLifetimeImage(
                    ffile,
                    uselines = uselines,
                    Gchan = np.array([PChan, PChan]), # duplicity works around bug
                    Rchan = np.array([SChan, SChan]),
                    #Gchan = np.array([2, 0]),
                    #Rchan = np.array([3, 1]),
                    **self.imreadkwargs,
                    framestop = int(self.Nframes))
        return GRYim

    def genPSfromGRYImage (self, image):
        assert hasattr(image, 'workLifetime'), \
            'image object must have workLifetime property'
        #extract lifetime decays
        image.sumLifetime() # generate decays
        image.P, image.S, _ = image.getTACS()
        del(image.decay) #redundant property

    def genDerivedFromPSDecays(self, image):
        TACPS = PandS2PS(image.P, image.S)
        #add VM and r variables
        VM, r = bp.genFr(TACPS, self.g_factor, shift = self.PSshift)
        image.VM = VM
        image.r = r
        return 1

    def genNormDecay(self, image, normimage,
                     decaytypes = ['VM', 'P', 'S'],
                     shift = 0,
                     bgrange = None):
        """implemented on image object level, because it is the least hassle
        has a drawbrack of being inflexible though"""
        from scipy.ndimage import gaussian_filter
        #subtract background if specified
        for decaytype in decaytypes:
            decay = getattr(image, decaytype)
            rawnormdecay = getattr(normimage, decaytype)
            # smooth normdecay to reduce shot noise.
            normdecay = gaussian_filter(rawnormdecay, 2)
            if bgrange is not None:
                decay = decay - np.mean(decay[bgrange[0]:bgrange[1]])
                normdecay = normdecay - np.mean(normdecay[bgrange[0]:bgrange[1]])
            decay, normdecay = bp.intshift(shift, decay, normdecay)
            decay = decay / normdecay
            setattr(image, decaytype + 'norm', decay)
        return 0

    def batchgenNormDecay(self, normimageG, normimageY,
                     decaytypes = ['VM', 'P', 'S'],
                     shift = 0,
                     bgrange = None,
                     **kwargs):
        """applies genNormDecay to each image in sampleSet"""
        for color in ['G', 'Y']:
            if color == 'G':
                normimage = normimageG
            elif color == 'Y':
                normimage = normimageY
            for image in self.images[color]:
                self.genNormDecay(image, normimage,
                                  decaytypes,
                                  shift,
                                  bgrange)
        return 0

    def batchFit1ltD0DA(self,
                      identifier,
                      D0dat = None,
                      fitrange = (25, 380),
                      decaytype = 'VM',
                      **kwargs):
        """makes simple Donor Only calibrated Donor Acceptor fits
        """
         #ugly workaround
        assert D0dat is not None, 'must give a Donor only decay'
        #read all DA decays
        DATACs = self.getDecay(decaytype)
        names = self.getPropertyList('name')
        dfrm = pd.DataFrame()
        #fit and plot DA
        plotout = os.path.join(self.resdir, identifier + 'D0DA1ltplots')
        aid.trymkdir(plotout)
        
        D0snip = D0dat[fitrange[0]:fitrange[1]]
        _, _, _, Donlymodel, chi2red_D0 = fitDA.fitDonly(D0snip, self.dt_glob)
        for name, DATAC in zip(names, DATACs):
            DAsnip = DATAC[fitrange[0]:fitrange[1]]
            popt, pcov, DAmodel, chi2red = \
                fitDA.fitDA1lt (DAsnip, D0snip, self.dt_glob)
            fitDA.pltDA_eps(DAsnip, D0snip, DAmodel, Donlymodel, name, popt,
                            chi2red, chi2red_D0, plotout, **kwargs)
            dfrm.at[name, 'xFRET'] = 1-popt[1]
            dfrm.at[name, 'kFRET'] = popt[2]
            dfrm.at[name, 'chi2red'] = chi2red
        outname = os.path.join(self.resdir, identifier + 'D0DAFitData.csv')
        dfrm.to_csv(outname)
        self.D0DA1ltdfrm = dfrm
        return dfrm

    def batchFit2ltD0DA(self,
                      identifier,
                      D0dat = None,
                      fitrange = (25, 380),
                      decaytype = 'VM',
                      **kwargs):
        """makes Donor Only calibrated (2lt) Donor Acceptor (2lt) fits
        """
        #ugly workaround
        assert D0dat is not None, 'must give a Donor only decay'
        #TODO split tauf, taux, E calculation in generic function
        #init and get data from object
        DATACs = self.getDecay(decaytype)
        names = self.getPropertyList('name')
        pnames = ['A_DA', 'xFRET1', 'xFRET2', 'kFRET1', 'kFRET2', 'bg']
        dfrm = pd.DataFrame()
        #fit and plot DA

        plotout = os.path.join(self.resdir, identifier + 'D0DA2ltplots')
        aid.trymkdir(plotout)
        D0snip = D0dat[fitrange[0]:fitrange[1]]
        poptD0, _, _, Donlymodel, chi2red_D0 = fitDA.fitDonly(D0snip, \
            self.dt_glob)
        x1, x2, tau1, tau2, _ = poptD0
        x1, x2 = [x1 / (x1 + x2), x2 / (x1 + x2)]
        k1, k2 = [1/tau1, 1 / tau2]
        tauxD0 = x1 * tau1 + x2 * tau2
        for name, DATAC in zip(names, DATACs):
            DAsnip = DATAC[fitrange[0]:fitrange[1]]
            popt, pcov, DAmodel, chi2red = \
                fitDA.fitDA2lt (DAsnip, D0snip, self.dt_glob)
            fitDA.pltDA_eps(DAsnip, D0snip, DAmodel, Donlymodel, name, popt,
                            chi2red, chi2red_D0, plotout, **kwargs)
            for p, pname in zip (popt, pnames):
                dfrm.at[name, pname] = p
            dfrm.at[name, 'chi2red'] = chi2red
            #     | kDA1  kDA2
            #___________________
            #kDO1 | x11   x12
            #kDO2 | x21   x22
            #tau_ij = (kD0i + kDAj)^-1
            #sum over all species species / fluorescence weighted
            #tau_x = SUMIJ xij * tau_ij
            taux = 0
            tauf = 0
            for xDA, kDA in zip(popt[[1,2]], popt[[3,4]]):
                for xD0, kD0 in zip ([x1, x2], [k1, k2]):
                    taux += xDA * xD0 * (1 / (kDA + kD0))
                    tauf += xDA * xD0 * (1 / (kDA + kD0))**2
            tauf = tauf / taux
            dfrm.at[name, 'taux'] = taux
            dfrm.at[name, 'tauf'] = tauf
            dfrm.at[name, 'E'] = 1-taux / tauxD0

        outname = os.path.join(self.resdir, identifier + 'D0DAFitData.csv')
        dfrm.to_csv(outname)
        self.D0DA2ltdfrm = dfrm
        return dfrm

    def batchFit2lt(self,
                    identifier,
                    fitrange = (20, 380),
                    decaytype = 'VM',
                    **kwargs):
        """batch fit D0 data assuming two lifetimes
        commonly for D0"""
        #**kwargs can take arguments that are not used
        #TODO split tauf, taux, E calculation in generic function
        pnames = ['x0', 'x1', 'tau0', 'tau1', 'bg']
        names = self.getPropertyList('name')
        dfrm = pd.DataFrame()
        plotout = os.path.join(self.resdir, identifier + 'D02ltplots')
        aid.trymkdir(plotout)
        #fit all data
        TACs = self.getDecay(decaytype)
        print(fitrange)
        assert len(TACs) != 0, 'TACs is empty'
        for TAC, name in zip(TACs, names):
            D0snip = TAC[fitrange[0]:fitrange[1]]
            popt, _, _, Donlymodel, chi2red = fitDA.fitDonly(D0snip, \
                self.dt_glob)
            fitDA.pltD0(D0snip, Donlymodel, name, plotout, dtime = self.dt_glob)
            #fill dataframe row
            for pname, p in zip(pnames, popt):
                dfrm.at[name, pname] = p
            x0, x1, tau0, tau1, bg = popt
            #calc derived vars
            tauf = (x0 * tau0**2 + x1 * tau1**2) / (x0 * tau0 + x1 * tau1)
            taux = (x0 * tau0 + x1 * tau1) / (x0 + x1)
            dfrm.at[name, 'tauf'] = tauf
            dfrm.at[name, 'taux'] = taux
            dfrm.at[name, 'chi2red'] = chi2red
            print('finished fitting with 2lt set %s' % name)
        outname = os.path.join(self.resdir, identifier + '2ltFitData.csv')
        dfrm.to_csv(outname)
        self.fit2ltdfrmrm = dfrm
        return dfrm

    def getDecay(self, decaytype = 'VM'):
        assert decaytype in ['VM', 'P', 'S'], \
            '%s is not a valid decay type' % decaytype
        return self.getPropertyList(decaytype)

    def getPropertyList(self, propertyName, channel = 'G'):
        """function scans all images in channel for PropertyName and
        returns a list of properties"""
        List = [getattr(GRYimage, propertyName) for GRYimage in \
                self.images[channel]]
        return List

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