import numpy as np
import ImageManipulation as IM
import os
import pickle
from GaussAnalysisPipeline import *

def analyseCFdir(locLst, options, wdir, 
                DTwoIstar = 0.03, garbageBrightness = 20, junkIstar = 0.30, outname = ''):
    ntacs = 256
    for i, file in enumerate(files[:nfiles]):
        if file[-4:] != '.ptu':
            print('not a .ptu file, skipping')
            continue
        print('analysing image no. %i' %i)
        ffile = os.path.join(wdir, file)
        #froifile = os.path.join(roidir, roifile)

        #load image, gate
        CLR = IM.processLifetimeImage(ffile.encode(), uselines = np.array([1,2]), ntacs = ntacs)
        CLR.loadLifetime()
        CLR.rebin(5, 5)# rebin 100x100px image into 20x20px, 50nm /px
        CLR.loadIntensity()
        locLst.append({})
        #loop over G and Y channels
        for color, bitmap in CLR.workIntensity.__dict__.items():
            if color in ['G', 'Y']:
                outdir = ffile[:-4] + '_' + color
                #fits 1, 2 or 3 gauss spots and determines which one is best
                #returns optimised parameter array
                bestfit, twoIstar, _ = fitNGauss (bitmap, options, 
                                                  DTwoIstar = DTwoIstar, garbageBrightness = garbageBrightness,
                                                  junkIstar = junkIstar,
                                                  verbose = False, outdir = outdir)
                #print(bestfit)
                #build array containing all the fit results
                locLst[-1][color] = Channel(bitmap)
                locLst[-1][color].fillSpotLst(bestfit)
    fpath = os.path.join(wdir, outname)
    with open(fpath, 'wb') as output:
        pickle.dump(locLst, output, 1)
    return locLst
    
def analyseSTEDdir(locLst, options, wdir, files, Ggate = 0, Ygate = 0, 
                   DTwoIstar = 0.03, garbageBrightness = 20, junkIstar = 0.30,
                   outname = ''):
    ntacs = 256
    ROIsize = 20

    for i, file in enumerate(files):
        if file[-4:] != '.ptu':
            print('not a .ptu file, skipping')
            continue
        print('analysing image no. %i' %i)
        ffile = os.path.join(wdir, file)

        #load image, gate
        CLR = IM.processLifetimeImage(ffile.encode(), uselines = np.array([1,2]), ntacs = ntacs)
        CLR.loadLifetime()
        CLR.gate(Ggate, 150, channel = 'G')
        CLR.gate(Ygate, 150, channel = 'Y')
        CLR.loadIntensity()

        #get ROI based on Y signal
        ROI = getROI(CLR.workIntensity.Y, ROIsize)    
        
        #check that ROI is not touching borders:
        try:
            crop(CLR.workIntensity.G, ROI)
        except IndexError:
            print('ROI touches image borders, skipping')
            continue
            
        locLst.append({})
        #loop over G and Y channels
        for color, bitmap in CLR.workIntensity.__dict__.items():
            if color in ['G', 'Y']:
                outdir = ffile[:-4] + '_' + color

                ROIsnip = crop(bitmap, ROI)
                #fits 1, 2 or 3 gauss spots and determines which one is best
                #returns optimised parameter array
                bestfit, twoIstar, _ = fitNGauss (ROIsnip, options, 
                                                  DTwoIstar = DTwoIstar, garbageBrightness = garbageBrightness,
                                                  junkIstar = junkIstar,
                                                  verbose = False, outdir = outdir)

                #build array containing all the fit results
                locLst[-1][color] = Channel(bitmap)
                locLst[-1][color].fillSpotLst(bestfit)
    
    fpath = os.path.join(wdir, outname)
    
    with open(fpath, 'wb') as output:
        pickle.dump(locLst, output, 1)
        
    return locLst
    
    
wdir = r'N:\Singlem\singlem20-1\January\path\to\your\CF\dir'
files = os.listdir(wdir)[:5] #for testing purposes we take only 5
options = optionsCluster(fitbg = 0, setbg = 0.2, ellipt_circ = 0)
RulerLst = []
analyseCFdir(RulerLst, options, DTwoIstar = 0.03, garbageBrightness = 20, junkIstar = 0.30, outname = 'somename.spots')

#to load data from disk
wdir = r'N:\Singlem\singlem20-1\January\path\to\your\CF\dir'
fpath = os.path.join(wdir, 'somename.spots')
with open(fpath, 'rb') as spotsfile:
    RulerLst_loaded = pickle.load(spotsfile)