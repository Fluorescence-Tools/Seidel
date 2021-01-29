import aid_functions as aid
import batchplot as bp
import os
import pandas as pd
import ImageManipulation as IM
import numpy as np
import fitDA
import gc

debug = False
if debug:
    import matplotlib.pyplot as plt

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
def calculateDerivedVariables(df, integrationtime = 1):
    """integration time is dwelltime * Nframes"""
    print(df)
    df['surfaceMax'] = df[['surfaceG', 'surfaceY']].max(axis = 1)
    for label in ['G', 'R', 'Y']:
        df['Br' + label] = df['N' + label+'-tot'] / df['surfaceMax']
        df['rate'+label] = df['Br' + label] / integrationtime
    df['BrG/BrY'] = df['BrG'] / df['BrY']
    df['BrG/BrR'] = df['BrG'] / df['BrR']
    df['BrR/BrY'] = df['BrR'] / df['BrY']
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


    
class sampleSet():
    """collection of attributes specific to either Donor only, or acceptor only"""
    def __init__(self, 
                 wdir,
                 relTACdir = 'TAC',
                 relresdir = 'results',
                 relimdir = 'images',
                 imreadkwargs = {'ntacs' : 1024,
                        'pulsetime' : 50,
                        'dwelltime': 20e-6,
                        'TAC_range': 4096},
                  g_factor = 1,
                  dataselect = (0, None),
                  dt_glob = 0.064,
                  FRETPIETACranges = [[0,380], [0, 380], [380,800]]
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
        self.wdir = wdir
        self.imreadkwargs = imreadkwargs
        self.TACdir = os.path.join(wdir, relTACdir)
        aid.trymkdir(self.TACdir)
        self.resdir = os.path.join(wdir, relresdir)
        aid.trymkdir(self.resdir)
        self.imdir = os.path.join(wdir, relimdir)
        aid.trymkdir(self.imdir)
        self.ptufiles = bp.appendOnPattern(wdir, 'ptu')[dataselect[0]: dataselect[1]]
        #self.names = [os.path.splitext(os.path.split(ptufile)[1])[0] \
        #              for ptufile in self.ptufiles]
        self.images = {'G': [], 'Y': []}#dict entries are channels, 
        self.g_factor = g_factor
        self.dt_glob = dt_glob
        self.FRETPIETACranges = FRETPIETACranges
        
    def getTACfileNames(self, ext = '_G_PS.dat'):
        self.TACfiles = bp.appendOnPattern(self.TACdir, ext)
        
#    def appendTAC(self, TAC):
#        assert type(TAC) == np.ndarray
#        self.TACs.append(TAC)
        
#    def loadTACfiles(self, append = False):
#        if not append: self.TACs = []
#        for TACfile in self.TACfiles:
#            fTACfile = os.path.join(self.TACdir, TACfile)
#            self.TACs.append(np.genfromtxt(fTACfile))
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


                 
    def analyzeLSMCells(self,
        identifier,
        threshold = 0,
        Nframes = -1,
        isSave = True,
        isCleanImage = True):
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
        """
        #init
        df = pd.DataFrame()
        #loop over each cell
        for index, file in enumerate(self.ptufiles):
            ffile = os.path.join(self.wdir, file).encode()
            #load Donor and acceptor channels
            #P is now called G, S is called R (nuisance)
            PSchannels = [[1,0], [5,4]]
            for PSchannel, label in zip(PSchannels, ['_G', '_Y']):
                PChan, SChan = PSchannel
                image = self.loadPSImage(ffile, PChan, SChan, Nframes)
                self.procesPSimage(image, threshold, isSave, 
                                   isCleanImage, label)
                self.images[label[1]].append(image)
        #assert self.names == [image.name for image in self.images['G']],\
        #    "two sets of naming variables should be identical"
        self.genImstatsdf(identifier, Nframes)
        return df
    
    def genImstatsdf(self, identifier, Nframes, 
                   channels = ['G', 'Y', 'Y'], 
                   #idea: consider using ['Donor', 'FRET', 'PIE'] nomenclature instead
                   labels = ['G', 'R', 'Y']):
        df = pd.DataFrame()
        for channel, TACrange, label in \
                zip(channels, self.FRETPIETACranges, labels):
            images = self.images[channel]
            for image in images:
                Np = np.sum(image.P[TACrange[0]:TACrange[1]])
                Ns = np.sum(image.S[TACrange[0]:TACrange[1]])
                Ntot = Np + self.g_factor * 2 * Ns
                #need to change this into Ny-tot, Ny-p and Ny-s
                df.at[image.name, 'N'+label+'-p'] = Np
                df.at[image.name, 'N'+label+'-s'] = Ns
                df.at[image.name, 'N'+label+'-tot'] = Ntot
                #data was masked previously
                surface = np.sum(image.workIntensity.G > 0) 
                df.at[image.name, 'surface'+label] = surface
        if Nframes == -1:
            print('number of frames not given,' \
                  + 'cannot calculate integration time and derived variables')
        else:
            integrationtime = self.imreadkwargs['dwelltime'] * Nframes
            df = calculateDerivedVariables(df, integrationtime)
        
        #save
        outdir = os.path.join(self.resdir, identifier + 'imstats.csv')
        df.to_csv(outdir)
        self.imstats = df
        return 0
            
    
    def procesPSimage(self, image, threshold, isSave, isCleanImage, label):
        """only does work on image, has a lot of dependencies, unwanted"""
        #both image object and sampleSet object keep names.
        #delete the first to avoid trouble
        #del image.name
        image.loadLifetime()
        image.loadIntensity()
        mask = image.buildMaskFromIntensityThreshold(
                threshold = threshold, sumchannels = ['G', 'R'])
        image.mask(mask, mode = 'intensity')
        self.genPSfromGRYImage(image)
        self.genDerivedFromPSDecays(image)
        if isSave:
            saveTACs(image, self.TACdir, label)
            image.saveWorkIntensityToTiff(self.imdir, image.name + 
                                      '_wmask' +label)
        if isCleanImage: #free memory intensive 3D array
            cleanImage(image)
        gc.collect() 
        return 1
    
    def loadPSImage(self, ffile, PChan, SChan, Nframes):
        GRYim = IM.processLifetimeImage(
                    ffile, 
                    uselines = np.array([1]), 
                    Gchan = np.array([PChan, PChan]), # duplicity works around bug
                    Rchan = np.array([SChan, SChan]), 
                    **self.imreadkwargs,
                    framestop = int(Nframes))
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
        VM, r = bp.genFr(TACPS, self.g_factor, shift = 0)
        image.VM = VM
        image.r = r
        return 1
    
    def genNormDecay(self, image, normimage, 
                     decaytypes = ['VM', 'P', 'S'],
                     shift = 0,
                     bgrange = None):
        """implemented on image object level, because it is the least hassle
        has a drawbrack of being inflexible though"""
        #subtract background if specified
        for decaytype in decaytypes:
            if bgrange is not None:
                decay = getattr(image, decaytype)
                normdecay = getattr(normimage, decaytype)
                decay = decay - np.mean(decay[bgrange[0]:bgrange[1]])
                normdecay = normdecay - np.mean(normdecay[bgrange[0]:bgrange[1]])
            shiftdecay, shiftnorm = bp.intshift(shift, decay, normdecay)
            setattr(image, decaytype + 'norm', shiftdecay/shiftnorm)
        return 0
    
    def batchgenNormDecay(self, normimageG, normimageY, 
                     decaytypes = ['VM', 'P', 'S'],
                     shift = 0,
                     bgrange = None):
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
                      decaytype = 'VM'):
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
        _, _, _, Donlymodel, chi2red_D0 = fitDA.fitDonly(D0snip)
        for name, DATAC in zip(names, DATACs):
            DAsnip = DATAC[fitrange[0]:fitrange[1]]
            popt, pcov, DAmodel, chi2red = \
                fitDA.fitDA1lt (DAsnip, D0snip, self.dt_glob)
            fitDA.pltDA_eps(DAsnip, D0snip, DAmodel, Donlymodel, name, popt, 
                            chi2red, chi2red_D0, plotout)
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
                      decaytype = 'VM'):
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
        poptD0, _, _, Donlymodel, chi2red_D0 = fitDA.fitDonly(D0snip)
        x1, x2, tau1, tau2, _ = poptD0
        x1, x2 = [x1 / (x1 + x2), x2 / (x1 + x2)]
        k1, k2 = [1/tau1, 1 / tau2]
        tauxD0 = x1 * tau1 + x2 * tau2
        for name, DATAC in zip(names, DATACs):            
            DAsnip = DATAC[fitrange[0]:fitrange[1]]
            popt, pcov, DAmodel, chi2red = \
                fitDA.fitDA2lt (DAsnip, D0snip, self.dt_glob)
            fitDA.pltDA_eps(DAsnip, D0snip, DAmodel, Donlymodel, name, popt, 
                            chi2red, chi2red_D0, plotout)
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
        self.D0DA2ltdfrmrm = dfrm
        return dfrm
    
    def batchFit2lt(self,
                    identifier,
                    D0dat = None,
                    fitrange = (20, 380),
                    decaytype = 'VM'):
        """batch fit D0 data assuming two lifetimes
        commonly for D0"""
        #D0dat is given as a keyword argument to give all fit functions
        #the same kwargs
        #TODO split tauf, taux, E calculation in generic function
        pnames = ['x0', 'x1', 'tau0', 'tau1', 'bg']
        names = self.getPropertyList('name')
        dfrm = pd.DataFrame()
        plotout = os.path.join(self.resdir, identifier + 'D02ltplots')
        aid.trymkdir(plotout)
        #fit all data
        TACs = self.getDecay(decaytype)
        assert len(TACs) != 0, 'TACs is empty'
        for TAC, name in zip(TACs, names):
            D0snip = TAC[fitrange[0]:fitrange[1]]
            popt, _, _, Donlymodel, chi2red = fitDA.fitDonly(D0snip)
            fitDA.pltD0(D0snip, Donlymodel, name, plotout)
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