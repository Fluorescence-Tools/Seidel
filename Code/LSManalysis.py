import aid_functions as aid
import batchplot as bp
import os
import pandas as pd
import ImageManipulation as IM
import numpy as np
import fitDA
import gc
import pickle 

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
    for label in ['_G', '_Y']:
        df['Br' + label] = df['I' + label] / df['surface' + label]
        df['rate'+label] = df['Br' + label] / integrationtime
    df['Sg/Sy'] = df['Br_G'] / df['Br_Y']
    return df

def cleanImage(image):
    del(image.baseLifetime)
    del(image.workLifetime)
    return 1

def saveTACs(image, TACdir):
    TACPS = PandS2PS(image.P, image.S)
    TACout = os.path.join(TACdir, image.name +'_PS.dat')
    VMout = os.path.join(TACdir, image.name  +'_VM.dat')
    rout = os.path.join(TACdir, image.name  +'_r.dat')
    np.savetxt(TACout, TACPS, fmt = '%i')
    np.savetxt(VMout, image.VM, fmt = '%i')
    np.savetxt(rout, image.r, fmt = '%.5e')
    return 1

def loadpickle(outname):
    with open(outname, 'rb') as f:
        return pickle.load(f)
    
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
                 imStatsHeader = ['I_G', 
                                  'surface_G', 
                                  'Br_G', 
                                  'rate_G', 
                                  'I_Y', 
                                  'surface_Y', 
                                  'Br_Y', 
                                  'rate_Y'],
                  g_factor = 1,
                  dataselect = (0, None),
                  dt_glob = 0.064
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
        self.imStatsHeader = imStatsHeader
        self.TACdir = os.path.join(wdir, relTACdir)
        aid.trymkdir(self.TACdir)
        self.resdir = os.path.join(wdir, relresdir)
        aid.trymkdir(self.resdir)
        self.imdir = os.path.join(wdir, relimdir)
        aid.trymkdir(self.imdir)
        self.ptufiles = bp.appendOnPattern(wdir, 'ptu')[dataselect[0]: dataselect[1]]
        self.names = [os.path.splitext(os.path.split(ptufile)[1])[0] \
                      for ptufile in self.ptufiles]
        self.images = {'G': [], 'Y': []}#dict entries are channels, 
        self.g_factor = g_factor
        self.dt_glob = dt_glob
        
    def getTACfileNames(self, ext = '_G_PS.dat'):
        self.TACfiles = bp.appendOnPattern(self.TACdir, ext)
        
    def genOldStyleHeader(self):
        self.imStatsHeader =  ['surface', 
                               'I_G', 
                               'Br_G', 
                               'rate_G', 
                               'I_Y', 
                               'Br_Y', 
                               'rate_Y']
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

    def getPropertyList(self, propertyName, channel = 'G'):
        List = [getattr(GRYimage, propertyName) for GRYimage in \
                self.images[channel]]
        return List

    def savepickle(self, outname):
        with open(outname, 'wb') as output:
            pickle.dump(self, output, 1)
                 
    def analyzeLSMCells(self,
        identifier,
        threshold = 50,
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
        df = pd.DataFrame(columns = self.imStatsHeader, index = self.names)
        #loop over each cell
        for index, file in enumerate(self.ptufiles):
            ffile = os.path.join(self.wdir, file).encode()
            #load Donor and acceptor channels
            #P is now called G, S is called R (nuisance)
            channels = [[1,0], [5,4]]
            for channel, label in zip(channels, ['_G', '_Y']):
                PChan, SChan = channel
                image = self.loadPSImage(ffile, PChan, SChan, Nframes)
                self.procesPSimage(image, threshold, isSave, 
                                   isCleanImage, label)
                #extract intensity image variables
                df['I'+label][image.name] = np.sum(image.P + image.S)
                surface = np.sum(image.P[image.P > 0])
                df['surface'+label][image.name] = surface
                #append lt image to struct
                self.images[label[1]].append(image)
            print('finished with file %s\n' % file[-20:])
        assert self.names == [image.name for image in self.images['G']],\
            "two sets of naming variables should be identical"
        df = calculateDerivedVariables(df)
        outdir = os.path.join(self.resdir, identifier + 'imstats.csv')
        df.to_csv(outdir)
        return df
    
    def procesPSimage(self, image, threshold, isSave, isCleanImage, label):
        """only does work on image, has a lot of dependencies, unwanted"""
        image.loadLifetime()
        image.loadIntensity()
        mask = image.buildMaskFromIntensityThreshold(
                threshold = threshold, sumchannels = ['G', 'R'])
        image.mask(mask, mode = 'intensity')
        self.genPSfromGRYImage(image)
        self.genDerivedFromPSDecays(image)
        if isSave:
            saveTACs(image, self.TACdir)
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
                    Gchan = np.array([PChan]),
                    Rchan = np.array([SChan]),
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
                      D0dat,
                      identifier,
                      fitrange = (25, 380),
                      decaytype = 'VM'):
        """makes simple Donor Only calibrated Donor Acceptor fits
        """
        #read all DA decays
        DATACs = self.getDecay(decaytype)
        #fit and plot DA
        xFRET = []
        kFRET = []
        chi2redLst = []
        
        plotout = os.path.join(self.resdir, identifier + 'D0DAplots')
        aid.trymkdir(plotout)
        for i, DATAC in enumerate(DATACs):
            D0snip = D0dat[fitrange[0]:fitrange[1]]
            DAsnip = DATAC[fitrange[0]:fitrange[1]]
            _, _, _, Donlymodel, chi2red_D0 = fitDA.fitDonly(D0snip)
            popt, pcov, DAmodel, chi2red = \
                fitDA.fitDA1lt (DAsnip, D0snip, self.dt_glob)
            fitDA.pltDA_eps(DAsnip, D0snip, DAmodel, Donlymodel, self.names[i], popt, 
                            chi2red, chi2red_D0, plotout)
            xFRET.append(1-popt[1])
            kFRET.append(popt[2])
            chi2redLst.append(chi2red)
            
        dfrm = pd.DataFrame(index = self.names)
        dfrm['xFRET'] = xFRET
        dfrm['kFRET'] = kFRET
        dfrm['chi2red'] = chi2redLst
        outname = os.path.join(self.resdir, identifier + 'D0DAFitData.csv')
        dfrm.to_csv(outname)
        self.D0DA1ltdfrm = dfrm
        return dfrm

    def batchFit2ltD0DA(self, 
                      D0dat,
                      identifier,
                      fitrange = (25, 380),
                      decaytype = 'VM'):
        """makes Donor Only calibrated (2lt) Donor Acceptor (2lt) fits
        """
        #read all DA decays
        DATACs = self.getDecay(decaytype)
        #fit and plot DA
        xFRET1 = []
        xFRET2 = []
        kFRET1 = []
        kFRET2 = []
        taufLst = []
        tauxLst = []
        ELst = []
        chi2redLst = []
        
        plotout = os.path.join(self.resdir, identifier + 'D0DAplots')
        aid.trymkdir(plotout)
        D0snip = D0dat[fitrange[0]:fitrange[1]]
        poptD0, _, _, Donlymodel, chi2red_D0 = fitDA.fitDonly(D0snip)
        x1, x2, tau1, tau2, _ = poptD0
        x1, x2 = [x1 / (x1 + x2), x2 / (x1 + x2)]
        k1, k2 = [1/tau1, 1 / tau2]
        tauxD0 = x1 * tau1 + x2 * tau2
        for i, DATAC in enumerate(DATACs):            
            DAsnip = DATAC[fitrange[0]:fitrange[1]]
            popt, pcov, DAmodel, chi2red = \
                fitDA.fitDA2lt (DAsnip, D0snip, self.dt_glob)
            print(popt)
            fitDA.pltDA_eps(DAsnip, D0snip, DAmodel, Donlymodel, self.names[i], popt, 
                            chi2red, chi2red_D0, plotout)
            xFRET1.append(popt[1])
            xFRET2.append(popt[2])
            kFRET1.append(popt[3])
            kFRET2.append(popt[4])
            chi2redLst.append(chi2red)
            

        #     | kDA1  kDA2
        #___________________
        #kDO1 | x11   x12
        #kDO2 | x21   x22
        #tau_ij = (kD0i + kDAj)^-1
        #sum over all species species / fluorescence weighted
        #tau_x = SUMIJ xij * tau_ij
        for n in range(len(xFRET1)):
            taux = 0
            tauf = 0
            for xDA, kDA in zip([xFRET1[n], xFRET2[n]], [kFRET1[n], kFRET2[n]]):
                for xD0, kD0 in zip ([x1, x2], [k1, k2]):
                    taux += xDA * xD0 * (1 / (kDA + kD0))
                    tauf += xDA * xD0 * (1 / (kDA + kD0))**2
            tauf = tauf / taux
            taufLst.append(tauf)
            tauxLst.append(taux)
            ELst.append(1-taux / tauxD0)
        dfrm = pd.DataFrame(index = self.names)
        dfrm['xFRET1'] = xFRET1
        dfrm['xFRET2'] = xFRET2
        dfrm['kFRET1'] = kFRET1
        dfrm['kFRET2'] = kFRET2
        dfrm['tau_x'] = tauxLst
        dfrm['tau_f'] = taufLst
        dfrm['E'] = ELst
        dfrm['chi2red'] = chi2redLst
        outname = os.path.join(self.resdir, identifier + 'D0DAFitData.csv')
        dfrm.to_csv(outname)
        self.D0DA2ltdfrm = dfrm
        return dfrm
    
    def batchFit2lt(self,
                    identifier,
                    fitrange = (20, 380),
                    decaytype = 'VM'):
        """batch fit D0 data assuming two lifetimes
        commonly for D0"""
         #prep empty arrays
        x0Lst = []
        x1Lst = []
        tau0Lst = []
        tau1Lst = []
        chi2redLst = []
        bgLst = []
        
        #fit all data
        TACs = self.getDecay(decaytype)
        assert len(TACs) != 0, 'TACs is empty'
        for i, TAC in enumerate(TACs):
            D0snip = TAC[fitrange[0]:fitrange[1]]
            popt, _, _, Donlymodel, chi2red = fitDA.fitDonly(D0snip)
            #append to Lists
            for var, val in zip([x0Lst, x1Lst, tau0Lst, tau1Lst, bgLst], popt):
                var.append(val)
            chi2redLst.append(chi2red)
            print('finished fitting with 2lt set %i' % i)
            
        #calc derived vars
        taufLst = [(x0 * tau0**2 + x1 * tau1**2) / (x0 * tau0 + x1 * tau1) 
                    for x0, x1, tau0, tau1 in zip(x0Lst, x1Lst, tau0Lst, tau1Lst)]
        tauxLst = [(x0 * tau0 + x1 * tau1) / (x0 + x1)
                    for x0, x1, tau0, tau1 in zip(x0Lst, x1Lst, tau0Lst, tau1Lst)]
        #add all lists to DataFrame
        dfrm = pd.DataFrame({'tauf' : taufLst,
                            'taux' : tauxLst, 
                            'x0' : x0Lst, 
                            'x1' : x1Lst, 
                            'tau0' : tau0Lst, 
                            'tau1' : tau1Lst, 
                            'chi2red' : chi2redLst}, 
                            index = self.names)
        outname = os.path.join(self.resdir, identifier + '2ltFitData.csv')
        dfrm.to_csv(outname)
        self.fit2ltdfrm = dfrm
        return dfrm
    
    def getDecay(self, decaytype = 'VM'):
        assert decaytype in ['VM', 'P', 'S'], \
            '%s is not a valid decay type' % decaytype
        return self.getPropertyList(decaytype)

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