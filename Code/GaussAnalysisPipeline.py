import numpy as np
import matplotlib.pyplot as plt
import os
import findPeaksLib
import GaussFits
import aid_functions as aid
#############################axes declarations#######################
#the x naming refers to the first dimension of the image
#the y naming refers to the second dimension of the image
#the plt.imshow has the following convention:
#the first dimension selects the rows, corresponding to y
#the second dimension selects columns, corresponding to x
######################class declarations#######################

class optionsCluster:
    def __init__(self, fitbg = 1, setmodel = 0, ellipt_circ = 1, setbg = 0.2):
        self.fitbg = fitbg #1 to fix, 0 to fit
        self.setmodel = setmodel #0: 1 Gauss, 1: 2 Gauss, 2: 3 Gauss
        self.ellipt_circ = ellipt_circ # 1: circular, 0 elliptical
        self.setbg = setbg #guess for background
        
    def transferOptions(self, params):
        assert type(params) == np.ndarray and params.shape == (18, )
        params[5] = self.setbg
        params[14] = self.fitbg
        params[15] = self.ellipt_circ
        params[16] = self.setmodel

class GaussSpot:
    def __init__(self, posx, posy, Amplitude, sigma, epsilon, background):
        self.posx = posx
        self.posy = posy
        self.A = Amplitude
        self.sigma = sigma
        self.eps = epsilon
        self.bg = background
        
    def setROI(self, xstart, ystart, xstop, ystop):
        self.xstart = xstart
        self.ystart = ystart
        self.xstop = xstop
        self.ystop = ystop
    def getROI(self):
        return np.array([self.xstart, self.ystart, self.xstop, self.ystop])
        
class Channel:
    def __init__(self, bitmap):
        self.bitmap = bitmap
        self.spotLst = []
        
    def fillSpotLst(self, params, ROI):
        sigma = params[3]
        eps = params[4]
        bg = params[5]
        for i in range(params[16].astype(np.int) + 1):
            if i == 0:
                posx = params[0] + ROI[0]
                posy = params[1] + ROI[1]
                A = params[2]
                self.spotLst.append(GaussSpot(posx, posy, A, sigma, eps, bg))
                self.spotLst[-1].setROI(*ROI)
            if i == 1:
                posx = params[6] + ROI[0]
                posy = params[7] + ROI[1]
                A = params[8]
                self.spotLst.append(GaussSpot(posx, posy, A, sigma, eps, bg))
                self.spotLst[-1].setROI(*ROI)
            if i == 2:
                posx = params[9] + ROI[0]
                posy = params[10] + ROI[1]
                A = params[11]
                self.spotLst.append(GaussSpot(posx, posy, A, sigma, eps, bg))
                self.spotLst[-1].setROI(*ROI)
                
                
################## Gauss fitting and judging functions ############

def fitNGauss(image, OptionsCluster,
                  DTwoIstar = 0.03, garbageBrightness = 50, junkIstar = 0.4,
                  verbose = False, outdir = None):
    """Analyse an image for up to three spots.
    First a parameter estimate is generated in genParamEstimate.
    Second, one, two and three Gaussians are fitted.
    Third, the best fit is chosen based in chooseBestfit.
    input: 
    image: a square 2D numpy array
    OptionsCluster: class object that carries fit options
    DTwoIstar: minimal difference above which two variables are considered 
        the same
    garbageBrightness: minimal integrated photons in spot for it not to be 
        considered garbage
    junkIstar: minimal absolute Istar value for photons not be junk
    verbose: if True, print spot fits in terminal
    outdir: if not empty, save spot fits to disc. Saving to disc can slow
        fitting significantly."""
    param_est = None
    param1Gauss = None
    param2Gauss = None
    param3Gauss = None
    param_best = None
    
    brightness = np.zeros([3,3])
    
    #obtain paramater estimates
    param_est = genParamEstimate(image)

    #set fitting options
    OptionsCluster.transferOptions(param_est)
    
    #fit 1 Gauss
    param1Gauss = param_est.copy()
    OptionsCluster.setmodel = 0
    OptionsCluster.transferOptions(param1Gauss)
    param1Gauss = GaussFits.Fit2DGauss(param1Gauss, image)
    
    #fit 2 Gauss
    param2Gauss = param_est.copy()
    OptionsCluster.setmodel = 1
    OptionsCluster.transferOptions(param2Gauss)
    param2Gauss = GaussFits.Fit2DGauss(param2Gauss, image)
    
    #fit 3 Gauss
    param3Gauss = param_est.copy()
    OptionsCluster.setmodel = 2
    OptionsCluster.transferOptions(param3Gauss)
    param3Gauss = GaussFits.Fit2DGauss(param3Gauss, image)
    
    if verbose or outdir:
        outfile = os.path.join(outdir, 'contourfits.png')
        fig, axs = plt.subplots(1,3, figsize = (12,3))
        drawSpotFit(axs[0], image, param1Gauss, '1 Gauss Fit')
       # pltFitResiduals(image, param1Gauss, '1 Gauss Fit', 
       #     verbose = verbose, outdir = outdir)
        drawSpotFit(axs[1], image, param2Gauss, '2 Gauss Fit')
       # pltFitResiduals(image, param2Gauss, '2 Gauss Fit', 
       #     verbose = verbose, outdir = outdir)
        drawSpotFit(axs[2], image, param3Gauss, '3 Gauss Fit')
       # pltFitResiduals(image, param3Gauss, '3 Gauss Fit', 
       #     verbose = verbose, outdir = outdir)
        plt.show()
        if outdir:
           aid.createPath(outfile)
           plt.savefig(outfile, dpi = 300, bbox_inches = 'tight')
    
    #choose best fit
    bestfit, twoIstar, brightness = chooseBestfit(param1Gauss, param2Gauss, param3Gauss, 
                                                 DTwoIstar = DTwoIstar, 
                                                 garbageBrightness = garbageBrightness, 
                                                 junkIstar = junkIstar,
                                                 verbose = verbose,
                                                 outdir = outdir)
    #return best parameters
    return bestfit, twoIstar, brightness 
    
def chooseBestfit(param1Gauss, param2Gauss, param3Gauss,
                  DTwoIstar = 0.03, garbageBrightness = 50, junkIstar = 0.4,
                 verbose = False, outdir = None):
    """
    input: optimesed parameters for 1, 2 and 3 Gauss fits
    
    isSignificantlyLower:
    check that the 2I* value is at least DTwoIstar lower than all simpler models
    
    isNoGarbagePeaks:
    check that no peaks have brightness less than garbageBrightness
    these peaks are considered to be noisepeaks
    
    isNoJunkIstar:
    check that Istar values are above junkIstar
    Istar value of noise is often much lower than regular data
    
    returns: most complicated fit that fulfills all conditions."""
    
    brightness = np.zeros([3,3])

    isSignificantlyLower = None 
    isNoJunkIstar = None
    isNoGarbagePeaks = None
    fullfillsAll = None
    
    brightness[0] = getSpotBrightness(param1Gauss)
    brightness[1] = getSpotBrightness(param2Gauss)
    brightness[2] = getSpotBrightness(param3Gauss)
    twoIstar = np.array([param1Gauss[17], param2Gauss[17], param3Gauss[17]])
    
    isNoJunkIstar = twoIstar > junkIstar
    
    #1 Gauss is simplest model
    isSignificantlyLower = np.array([True, False, False])
    for i in [1, 2]:
        isSignificantlyLower[i] = (twoIstar[i] + DTwoIstar < twoIstar[: i]).all()
    
    isNoGarbagePeaks = np.array([False, False, False])
    for i in range(3):
        isNoGarbagePeaks[i] = (brightness[i,:i + 1] > garbageBrightness).all()
         
    #combine conditons
    fullfillsAll = np.logical_and(isSignificantlyLower, isNoJunkIstar)
    fullfillsAll = np.logical_and(fullfillsAll, isNoGarbagePeaks)
    
    np.set_printoptions(precision=3)
    msgLst = []
    msgLst.append('optimised 1 Gauss fit is: \n' + str(param1Gauss) + '\n')
    msgLst.append('optimised 2 Gauss fit is: \n' + str(param2Gauss) + '\n')
    msgLst.append('optimised 3 Gauss fit is: \n' + str(param3Gauss) + '\n')
    msgLst.append('isSignificantlyLower : ' + str(isSignificantlyLower))
    msgLst.append('isNoJunkIstar : ' + str(isNoJunkIstar))
    msgLst.append('isNoGarbagePeaks : ' + str(isNoGarbagePeaks))
    msgLst.append('2I* : ' + str(twoIstar))
    msgLst.append('fullfills all conditions' + str(fullfillsAll))
    
    if verbose:
        for msg in msgLst:
            print(msg)
    if outdir:
        try:
            os.mkdir(outdir)
        except:
            pass
        f = open(os.path.join(outdir, 'fitinfo.txt'), 'w')
        for msg in msgLst:
            f.write(msg)
        f.close()

    
    #choose most complex model that fullfills all conditions
    for i in [2, 1, 0]:
        if fullfillsAll[i]:
            if i == 2:
                bestfit = param3Gauss
                break
            elif i == 1:
                bestfit = param2Gauss
                break
            elif i == 0:
                bestfit = param1Gauss
                break
        #if no suitable candidate is found, model param is set -1
        bestfit = np.zeros(18)
        bestfit[16] = -1
    
    return bestfit, twoIstar, brightness
    
def getSpotBrightness(params):
    """get the total number of photons contained in one spot
    parameter model is used to determine how many spots are calculated
    uncalculated spots are filled with zero.
    returns length 3 array."""
    brightness = np.zeros(3)
    factor = params[3]**2 * 2 * np.pi
    for i in range(params[16].astype(np.int)+1):
        if i == 0:
            brightness[i] = params[2] * factor
        if i == 1:
            brightness[i] = params[8] * factor
        if i == 2:
            brightness[i] = params[11] * factor
    return brightness
    
def genParamEstimate(image):
#write docstring why this format is chosen
    peaks = findPeaksLib.findPeaks(image)
    params_est = np.array([peaks[0,0], #x0
                        peaks[0,1],#y0
                        image[peaks[0,0], peaks[0,1]], #A0
                        1, #sigma
                        1,#eps
                        0, #bg, set using OptionsCluster
                        peaks[1,0], # x1
                        peaks[1,1], #y1
                        image[peaks[1,0], peaks[1,1]], #A1
                        peaks[2,0], #x2
                        peaks[2,1], #y2
                        image[peaks[2,0], peaks[2,1]], #A2
                        0, #info
                        0, #wi_nowi
                        0, #fitbg
                        0, #ellipt_circ
                        0, #model
                        0]) # two Istar
    return params_est
    
    


###################plotting functions #########################

def drawSpotFit(ax, image, fit, title):
#add docstring
    model= np.zeros(image.shape)
    if fit[16] == 0:
        model = GaussFits.model2DGaussian(fit, model)
    if fit[16] == 1:
        model = GaussFits.modelTwo2DGaussian(fit, model)
    if fit[16] == 2:
        model = GaussFits.modelTwo2DGaussian(fit, model)
    ax.contour(model, levels = np.array([0.01, 0.5, 1, 3, 5, 10, 20, 30, 40], \
        dtype = np.double))
    im = ax.imshow(image, cmap = 'hot')
    ax.set_title(title)
    plt.colorbar(im, ax = ax)
    return True

def pltFitResiduals(image, fit, title, verbose = True, outdir = None):
#add docstring
    model= np.zeros(image.shape)
    if fit[16] == 0:
        model = GaussFits.model2DGaussian(fit, model)
    if fit[16] == 1:
        model = GaussFits.modelTwo2DGaussian(fit, model)
    if fit[16] == 2:
        model = GaussFits.modelTwo2DGaussian(fit, model)
    #plt.imshow(image - model, cmap = 'hot')
    plt.imshow(model - image *( 1 + np.log(model / image)), cmap = 'hot')
    plt.colorbar()

    if outdir:
        try:
            os.mkdir(outdir)
        except:
            pass
        plt.savefig(os.path.join(outdir, 'FitResiduals_%iSpots.png' % fit[16]), dpi = 300, bbox_inches = 'tight')
    if verbose:
        plt.show()
    plt.clf()