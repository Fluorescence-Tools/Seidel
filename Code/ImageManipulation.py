import cpp_wrappers
import FRCfuncs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from scipy.optimize import curve_fit
import os
import copy
from PIL import Image

debug = False

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
        
class GRYDecay():
    def __init__(self, decG, decR, decY):
        assert (len(decG.shape) == 1 and len(decR.shape) == 1 and len(decY.shape) == 1), \
            "arrays must be 1D"
        self.G = decG
        self.R = decR
        self.Y = decY
    
class processLifetimeImage:
    def __init__(self, fname, uselines = np.array([1,0]), Gchan = np.array([0,2]), \
        Rchan = np.array([1,3]), Ychan = np.array([1,3]), ntacs = 256, pulsetime = 25):
        """Lifetime Image class
        fname must be a .ptu file path denoted in bytes.
                Uselines codes for the used line steps. 
            0 denotes that the line is ignored
            1 denotes that the line is FRET sensitized
            2 denotes that it is PIE sensitized
        Gchan, Rchan and Ychan code the channel numbers
        ntacs denote the number of bins in the tac decay histogram
            should be a multiple of 2, max tacs 2**15 = 32768
        pulsetime in ns."""
        self.baseLifetime = None #immutable after initialization
        self.baseIntensity = None #immutable after initialization
        self.workLifetime = None #mutable
        self.workIntensity = None #mutable
        self.decay = None #mutable
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
        
    def mask(self, mask, Channel = 'all', mode = 'lifetime'):
    #consider moving this operation to c, as it is very time-consuming
        """mask worklifetime of workIntensity image.
        mask should be interger 0 or 1.
        Channel should be 'all', 'G', 'R' or 'Y'
        mode should be 'lifetime' or intensity'"""
        assert (np.unique(mask) == np.array([0,1])).all(), 'mask has values other than 0 and 1'
        if mode == 'lifetime':
            #loop over all tacs. This avoids saving 3D mask to disk
            for tac in range(self.ntacs):
                if Channel == 'all':
                    self.workLifetime.G[:,:,tac] = self.workLifetime.G[:,:,tac] * mask
                    self.workLifetime.R[:,:,tac] = self.workLifetime.R[:,:,tac] * mask
                    self.workLifetime.Y[:,:,tac] = self.workLifetime.Y[:,:,tac] * mask
                elif Channel == 'G':
                    self.workLifetime.G[:,:,tac] = self.workLifetime.G[:,:,tac] * mask
                elif Channel == 'R':
                    self.workLifetime.R[:,:,tac] = self.workLifetime.R[:,:,tac] * mask
                elif Channel == 'Y':
                    self.workLifetime.Y[:,:,tac] = self.workLifetime.Y[:,:,tac] * mask
        elif mode == 'intensity':
            if Channel == 'all':
                self.workIntensity.G = self.workIntensity.G * mask
                self.workIntensity.R = self.workIntensity.R * mask
                self.workIntensity.Y = self.workIntensity.Y * mask
            elif Channel == 'G':
                self.workIntensity.G = self.workIntensity.G * mask
            elif Channel == 'R':
                self.workIntensity.R = self.workIntensity.R * mask
            elif Channel == 'Y':
                self.workIntensity.Y = self.workIntensity.Y * mask
        return 0
    
    def sumLifetime(self):
        """Integrate over the pixels in worklifetime to generate
            a TAC decay. Operation runs over all channels and is
            stored in self.decay."""
        decG = np.sum(self.workLifetime.G, axis = (0,1) )
        decR = np.sum(self.workLifetime.R, axis = (0,1) )
        decY = np.sum(self.workLifetime.Y, axis = (0,1) )
        self.decay = GRYDecay(decG, decR, decY)
        return 0
        
    def buildMaskFromROI(self, ROIS):
        """the format of the ROIS is inherited from AnI.
            Ani has the reversed axis convention from python.
            i.e. x in Ani corresponds to y in python.
        The logical format for python is therefore:
            ycorner, yrange, xcorner, xrange, ysubstitute, x substitute. With:
            ycorner = ycorner_fromAni + x_offset_fromPython."""
        assert np.issubdtype(ROIS.dtype, np.integer), 'ROIS is not a integer array'
        xshape, yshape = self.workIntensity.G.shape
        mask = np.zeros((xshape, yshape), dtype = np.int)
        for ROI in ROIS:
            if ROI[2] < 0 or ROI[0] < 0 or ROI[2] + ROI[3] >= xshape or \
                    ROI[0] + ROI[1] >= yshape:
                print('ROI is touching image borders, skipping')
                continue
            mask[ROI[2] : ROI[2] + ROI[3], ROI[0] : ROI[0] + ROI[1] ] = 1
        return mask
        
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
        
    def saveGRYTif(GRYobj, outfolder, preposition = '', xmin = 0, xmax = -1, ymin = 0, ymax = -1): 
        """function takes GRY object and saves to outfolder in tif format in 32bit float .tif format
            preposition allows adding identifier to default filenames.
            xmin. xmax, ymin, ymax allow saving snip of np array"""
        #convert all strings to bytes
        try: outfolder = outfolder.encode()
        except: pass
        try: preposition = preposition.encode()
        except: pass
        try:
            os.mkdir(outfolder)
            print('creating new folder %s\n' %os.path.join(outfolder))
        except: 'folder already exists'
        im = Image.fromarray(GRYobj.G[xmin:xmax, ymin:ymax])
        outname = os.path.join(outfolder, preposition + b'imG.tif').decode()
        im.save(outname)
        im = Image.fromarray(GRYobj.R[xmin:xmax, ymin:ymax])
        outname = os.path.join(outfolder, preposition + b'imR.tif').decode()
        im.save(outname)
        im = Image.fromarray(GRYobj.Y[xmin:xmax, ymin:ymax])
        outname = os.path.join(outfolder, preposition + b'imY.tif').decode()
        im.save(outname)
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

        assert type(fname) == bytes, 'Error, fname is not bytes type'
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
        
class fitImage:
    """Class containing fitNGauss and supporting functions.
    Class interfaces with cppwrappers and c.
    refer to fitNGauss documentation."""
    
    @classmethod
    def plotNGauss (cls, image, params, identifier, savefig, showfig, model, outdir):
        """routine for plotting two Gaussian centers on an image
        sides of the box is three sigma.
        params should be ordered according to Gaussian c code convention.
        savefig and showfig are bools to trigger saving to outdir and display, respectively.
        model determines whether one, two or three gaussians are fitted
        identifier is a string used for graph title and filename.
        Currently only 2 Gauss implemented"""
        #draw image
        plt.imshow(image)
        plt.colorbar()
        ax = plt.gca()
        plt.title(identifier, fontsize = 20)
        
        #Draw rectangles to indicate ROIS
        cls._drawrect(params[0], params[1], params[3], ax)
        if params[16] >= 1 :
            cls._drawrect(params[6], params[7], params[3], ax)
        if params[16] >= 2 :
            cls._drawrect(params[9], params[10], params[3], ax)

        if savefig:
            try:
                os.mkdir (outdir)
            except FileExistsError:
                pass
            if outdir == '':
                raise ValueError ('No output dir to save figure specified, aborting.')
            outname = os.path.join( outdir, identifier + r'.png')
            plt.savefig(outname, dpi = 300, bbox_inches = 'tight')
        if showfig:
            plt.show()
        else:
            plt.clf()
        return 0

    def modelTwo2DGaussians(coord, x0, y0, A0, sigma, ellipticity, bg, x1, y1, A1):
        """models two 2D Gaussian functions. Coordinate order is taken from c code.
        ellipticity remains unused."""
        x = coord[0]
        y = coord[1]  
        g = bg + A0 * np.exp( - ((x-x0)**2 + (y-y0)**2) / (2 * sigma**2)) + \
            A1 * np.exp( - ((x-x1)**2 + (y-y1)**2) / (2 * sigma**2))
        return g.ravel()
        
    @classmethod
    def fitTwo2DGaussian(cls, image, params0):
        """fit image with two 2D Gaussians using python code and optimalisation"""
        # Create x and y indices
        x = np.arange(image.shape[0])
        y = np.arange(image.shape[1])
        coord = np.meshgrid(x, y)
        
        popt, pcov = curve_fit(cls.modelTwo2DGaussians, coord, image.ravel(), p0 = params0)
        
        return popt, pcov
        
    def IstarThres(image):
        """calculate an accepted value of Istar assuming that the model
        describes the object very well
        Addition factor 0.5 is chosen arbitrarily"""
        
        thres = 0
        for el in image.ravel():
            if el == 0:
                pass
            else:
                thres += el - el * np.log(el)
        return thres / image.ravel().shape[0] + 0.8
    
    @classmethod
    def tryfit(cls, image, model, identifier):
        Ntries = 50
        Istar = None
        for counter in range(Ntries):
            #generate params0, add offset to sigma to avoid local minima of small rois
            xshape, yshape = image.shape
            Nmax = image.max()
            params0 = np.random.random(17) * \
                        np.array([xshape, yshape, Nmax, 5, 1, 2, xshape, yshape, Nmax, xshape, yshape, Nmax, \
                                 0, 0, 0, 0, 0]) + \
                        np.array([0, 0, 0, 2, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0])
            params = params0

            #decide which fit to use
            if model == 'one2DGaussian_py':
                raise NotImplemented
            elif model == 'two2DGaussian_py':
                try:
                    params[:9] , _ = cls.fitTwo2DGaussian(image, params0[:9])
                except RuntimeError:
                    print ('%s attempt %i: python fit did not converge, trying again' % (identifier, counter))
                    continue
            elif model == 'three2DGaussian_py':
                raise NotImplemented
            elif model == 'one2DGaussian_c':
                params0[16] = 0
                params, Istar = cpp_wrappers.fit2DGaussian_wrap(params0, image, debug = debug)
            elif model == 'two2DGaussian_c':
                params0[16] = 1
                params, Istar = cpp_wrappers.fit2DGaussian_wrap(params0, image, debug = debug)
            elif model == 'three2DGaussian_c':
                params0[16] = 2.01
                params, Istar = cpp_wrappers.fit2DGaussian_wrap(params0, image, debug = debug)
            else:
                raise ValueError('invalid model')
                
            #check if c routine returned convergence error or bad Istar value
            thresh = cls.IstarThres(image)
            if params[12] == -1 or params[12] == 1 or params[12] == 2:
                print('%s attempt %i: c routine returned info %.0f, trying again' % 
                      (identifier, counter, params[12]))
                continue
            elif Istar != None and Istar > thresh:
                print('%s attempt %i: value of Istar is %.2f and Istarthreshold is %.2f, trying again' %
                              (identifier, counter, Istar, thresh) )
                continue
            print('%s:value of Istar is %.2f and Istarthreshold is %.2f, success!' %
                              (identifier, Istar, thresh) )
            break
            
        #signal a failed fit
        if counter == Ntries - 1:
            params[12] = -1

        return params
       
    @classmethod
    def fitNGauss(cls, image, model, identifier, savefig, showfig, outdir = ''):
        """fits N Gaussian model to 2D image. 
        internally params is used among the python and c subroutines.
        params = [0: x0, 1: y0, 2: A0, 3: sigma, 4: ellipticity, 5: bg, 6: x1, 7: y1, 8: A1, 9: x2, 10: y2, 11: A2, \
                12: info, 13: wi_nowi, 14: fit_bg, 15: ellipt_circ, 16: model]
        info contains information from the fitting algorithm
        wi_nowi contains weights or no_wights, outdated debree
        fit_bg asks if background is fitted. 0 -> bg is fitted
        ellipt_circ  determines if elliptical fits are allowed. 1 means circular, 0 elliptical
            Only relevant for model2DGaussian
        model determines the model to be used:
            one2DGaussian_py: not implemented
            two2DGaussian_py
            three2DGaussian_py: not implemented
            one2DGaussian_c
            two2DGaussian_c
            three2DGaussian_c
        savefig and showfig are bools to determine saving and showing of fits
        
        returns: params"""
        #run tryfit
        params = cls.tryfit (image, model, identifier)
        
        #run plotNGauss
        cls.plotNGauss (image, params, identifier, savefig, showfig, model, outdir)
        
        return params

    def _drawrect(centerx, centery, sigma, ax):
        """subroutine to draw rectangle on ax object
        sides of object are 3 sigma"""
        x0, y0, w, h = [centerx - 1.5 * sigma, centery - 1.5 * sigma, 3 * sigma, 3 * sigma]
        rect0 = patches.Rectangle((x0, y0),w, h, linewidth=1,edgecolor='r',facecolor='none')
        ax.add_patch(rect0)
        ax.scatter(int(centerx + 0.5), int(centery + 0.5), s = 3, c = 'r')
        return 0