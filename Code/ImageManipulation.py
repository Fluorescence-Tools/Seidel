import cpp_wrappers
import FRCfuncs
import numpy as np
import os
import copy

class GRYLifetime():
    def __init__(self, imG, imR, imY):
        assert (len(imG.shape) == 3 and len(imR.shape) == 3 and len(imY.shape) == 3), \
            "arrays must be 3D"
        self.G = imG
        self.R = imR
        self.Y = imY
    
class GRYIntensity():
    def __init__(self, imG, imR, imY):
        assert (len(imG.shape) == 2 and len(imR.shape) == 2 and len(imY.shape) == 2), \
            "arrays must be 2D"
        self.G = imG
        self.R = imR
        self.Y = imY
    
class processLifetimeImage:
    def __init__(self, fname, uselines = np.array([1,0]), Gchan = np.array([0,2]), \
        Rchan = np.array([1,3]), Ychan = np.array([1,3]), ntacs = 256, pulsetime = 25):
        self.baseLifetime = None #immutable after initialization
        self.baseIntensity = None #immutable after initialization
        self.workLifetime = None #mutable
        self.workIntensity = None #mutable
        self.ntacs = ntacs
        self.tac2time = pulsetime / ntacs #pulsetime in ns
        assert( type(uselines) == np.ndarray and \
               type(Gchan) == np.ndarray and\
               type(Rchan) == np.ndarray and \
               type(Ychan) == np.ndarray), \
            "uselines, Gchan, Ychan and Rchan must be numpy array type!"
        #initialize baseLifetimeObject
        self._makeLifetime(fname, uselines, Gchan, Rchan, Ychan, ntacs)
        #initialize baseIntensityObject
        self._makeIntensity()
        
    def gate(self, mingate, maxgate):
        self.workLifetime.G[:,:,:mingate] = 0
        self.workLifetime.G[:,:,maxgate:] = 0
        self.workLifetime.R[:,:,:mingate] = 0
        self.workLifetime.R[:,:,maxgate:] = 0
        self.workLifetime.Y[:,:,:mingate] = 0
        self.workLifetime.Y[:,:,maxgate:] = 0
        return 0
    
    def rebin(self, xfactor, yfactor):
        """reshape by xfactor and yfactor"""
        shape = self.workLifetime.G.shape
        xshape = shape[0]
        yshape = shape[1]
        nbins = shape[2]
        self.workLifetime.G = self.workLifetime.G.reshape(
            xshape // xfactor, xfactor, yshape // yfactor, yfactor, \
            nbins).sum(axis = 3).sum(axis = 1)
        self.workLifetime.R = self.workLifetime.R.reshape(
            xshape // xfactor, xfactor, yshape // yfactor, yfactor, \
            nbins).sum(axis = 3).sum(axis = 1)
        self.workLifetime.Y = self.workLifetime.Y.reshape(
            xshape // xfactor, xfactor, yshape // yfactor, yfactor, \
            nbins).sum(axis = 3).sum(axis = 1)
        return 0
    
    def filterLifetime(self, mode = 'xyz', window = 3):
        if mode == 'smear_lifetime':
            # even sized window, weighted toward + direction
            if window % 2 == 0: 
                lwindow = window // 2 - 1
                rwindow = window // 2 + 1
            #odd sized window, distributed evenly left and right
            elif window % 2 == 1: 
                lwindow = window // 2 - 1
                rwindow = window // 2 + 2
            
            ltG = self.workLifetime.G
            ltR = self.workLifetime.R
            ltY = self.workLifetime.Y
            resG = np.zeros(ltG.shape)
            resR = np.zeros(ltR.shape)
            resY = np.zeros(ltY.shape)
            xshape, yshape, _ = ltG.shape
            for i in range(xshape):
                #check boundaries
                while i - lwindow < 0:
                    lwindow -= 1
                while i + lwindow > xshape:
                    lwindow -= 1
                for j in range(yshape):
                    while j - lwindow < 0:
                        lwindow -= 1
                    while j + lwindow > xshape:
                        lwindow -= 1
                    resG[i,j] = ltG[i - lwindow : i + rwindow, 
                               j - lwindow: j + rwindow,:].sum(axis = (0,1))
                    resR[i,j] = ltR[i - lwindow : i + rwindow, 
                               j - lwindow: j + rwindow,:].sum(axis = (0,1))
                    resY[i,j] = ltY[i - lwindow : i + rwindow, 
                               j - lwindow: j + rwindow,:].sum(axis = (0,1))
            self.workLifetime.G = resG
            self.workLifetime.R = resR
            self.workLifetime.Y = resY
        return 0
            
    def meanArrivalTime(self, tacstart, threshold = 0):
        """calculates mean arrival time of workLifetime
        ignores all photons arriving bebore tacstart.
        Sets lifetime to 0 for all pixels with less than
        threshold photons."""
        weights = (np.arange(self.ntacs) - tacstart)*self.tac2time
        #ensure no negative weights. Also ensure that photons before
        #tacstart are ignored.
        for i in range(self.ntacs):
            if weights[i] < 0 :
                weights[i] = 0
        imG = np.zeros(self.workLifetime.G.shape[:2])
        imR = np.zeros(self.workLifetime.R.shape[:2])
        imY = np.zeros(self.workLifetime.Y.shape[:2])
        for i in range(imG.shape[0]):
            for j in range(imG.shape[1]):
                pixelsum = np.sum(self.workLifetime.G[i,j,:])
                if pixelsum < threshold:
                    imG[i,j] = 0
                else:
                    imG[i,j] = self.workLifetime.G[i,j,:].dot(weights) / pixelsum
                pixelsum = np.sum(self.workLifetime.R[i,j,:])
                if pixelsum < threshold:
                    imR[i,j] = 0
                else:
                    imR[i,j] = self.workLifetime.R[i,j,:].dot(weights) / pixelsum
                pixelsum = np.sum(self.workLifetime.Y[i,j,:])
                if pixelsum < threshold:
                    imY[i,j] = 0
                else:
                    imY[i,j] = self.workLifetime.Y[i,j,:].dot(weights) / pixelsum
        self.workIntensity = GRYIntensity(imG, imR, imY)
        return 0
        
    
    def filterIntensity(self, mode = 'xyz'):
        pass
    
    def lifetimeFit(self):
        pass
    
    def loadLifetime(self):
        self.workLifetime = copy.deepcopy(self.baseLifetime)
        return 0
        
    def loadIntensity(self):
        imG = self.workLifetime.G.sum(axis=2)
        imR = self.workLifetime.R.sum(axis=2)
        imY = self.workLifetime.Y.sum(axis=2)
        self.workIntensity = GRYIntensity(imG, imR, imY)
        return 0
    
    def getLifetime(self):
        return self.workLifetime
    
    def getIntensity(self):
        return self.workIntensity
    
    def setLifetime(self, LifetimeObject):
        assert(type(LifetimeObject = GRYLifetime)),'not a GRY lifetime object'
        self.workLifetime = LifetimeObject
        return 0
        
    def setIntensity(self, IntensityObject):
        assert(type(IntensityObject = GRYIntensity)),'not a GRY intensity object'
        self.workIntensity = IntensityObject
        return 0
        
    def _makeIntensity(self):
        """sum tac histograms to obtain intensity image"""
        imG = self.baseLifetime.G.sum(axis=2)
        imR = self.baseLifetime.R.sum(axis=2)
        imY = self.baseLifetime.Y.sum(axis=2)
        self.baseIntensity = GRYIntensity(imG, imR, imY)
        return 0
    
    
    
    def _makeLifetime(self, fname, uselines, Gchan, Rchan, Ychan, ntacs):
        """initialisation routine for lifetime image manipulation class.
        Loads the .ptu file located at fname.
        Uselines codes for the used line steps. 
            0 denotes that the line is ignored
            1 denotes that the line is FRET sensitized
            2 denotes that it is PIE sensitized
        Uselines should not be used for binning. A separate functions exists for 
            that.
        Gchan, Rchan and Ychan indicates the detection channels for the 
            corresponding colors. No distinction between P and S 
            polarisation is made. Generally R and Y have the same channels.
        ntacs is the amount of bins used for the tac histogram. Decrease 
            to safe memory usage. Computational efficiency is minimal
            effected.
        returns imG, imR, imY"""
        uselines = uselines.astype(np.ubyte)
        Gchan = Gchan.astype(np.ubyte)
        Rchan = Rchan.astype(np.ubyte)
        Ychan = Ychan.astype(np.ubyte)
        root, file = os.path.split(fname)
        NumRecords = cpp_wrappers.ptuHeader_wrap (fname)
        eventN, tac, t, can = cpp_wrappers.ptu_wrap(fname, NumRecords)
        root, file = os.path.split(fname)
        name, _ = os.path.splitext(file)
        header_name = os.path.join(root, b"header", name + b".txt")
        print('number of records is ' + str(NumRecords))

        dimX, dimY, dwelltime, counttime = cpp_wrappers.read_header(header_name)
        imG, imR, imY = cpp_wrappers.genGRYLifetimeWrap(eventN, tac, t, can, dimX, dimY, ntacs, \
            dwelltime, counttime, NumRecords, uselines, Gchan, Rchan, Ychan)
        self.baseLifetime = GRYLifetime(imG, imR, imY)
        return 0