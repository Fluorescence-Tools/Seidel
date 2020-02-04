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

class FRETind:
    def __init__(self):
        self.pseudoE = None
        self.tauG = None
        self.tauR = None
        self.tauY = None
        self.dist = None

def analyseCFdir(
    locLst, options, wdir, nfiles = -1, outname = 'SpotsObject.spots',
    framestop = -1, verbose = False, saveplot = False
    ):
    ntacs = 256
    for i, file in enumerate(files[:nfiles]):
        if file[-4:] != '.ptu':
            print('not a .ptu file, skipping')
            continue
        print('analysing image no. %i' %i)
        ffile = os.path.join(wdir, file)
        #froifile = os.path.join(roidir, roifile)

        #load image, gate
        CLR = IM.processLifetimeImage(
            ffile.encode(), uselines = np.array([1,2]), ntacs = ntacs,
            framestop = framestop)
        CLR.loadLifetime()
        CLR.rebin(5, 5)# rebin 100x100px image into 20x20px, 50nm /px
        CLR.loadIntensity()
        locLst.append({})
        locLst[-1]['filepath'] = ffile
        locLst[-1]['ROI'] = ROI
        #loop over G and Y channels
        for color, bitmap in CLR.workIntensity.__dict__.items():
            if color in ['G', 'Y']:
                if saveplot:
                    outdir = ffile[:-4] + '_' + color
                else:
                    outdir = ''
                #fits 1, 2 or 3 gauss spots and determines which one is best
                #returns optimised parameter array
                bestfit, twoIstar, _ = fitNGauss (
                    bitmap, options, verbose = False, outdir = outdir
                    )
                #print(bestfit)
                #build array containing all the fit results
                locLst[-1][color] = Channel(bitmap)
                locLst[-1][color].fillSpotLst(bestfit)
    fpath = os.path.join(wdir, outname)
    with open(fpath, 'wb') as output:
        pickle.dump(locLst, output, 1)
    return locLst

def analyseSTEDdir(
    locLst, options, wdir, files, Ggate = 0, Ygate = 0, 
    DTwoIstar = 0.03, garbageBrightness = 20, junkIstar = 0.30,
    outname = '', framestop = -1, verbose = False, saveplot = False):
    ntacs = 256
    ROIsize = 30

    for i, file in enumerate(files):
        if file[-4:] != '.ptu':
            print('not a .ptu file, skipping')
            continue
        print('analysing image no. %i' %i)
        ffile = os.path.join(wdir, file)

        #load image, gate
        CLR = IM.processLifetimeImage(
            ffile.encode(), uselines = np.array([1,2]), ntacs = ntacs,
            framestop = framestop)
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
        locLst[-1]['filepath'] = ffile
        locLst[-1]['ROI'] = ROI
        #loop over G and Y channels
        for color, bitmap in CLR.workIntensity.__dict__.items():
            if color in ['G', 'Y']:
                if saveplot:
                    outdir = ffile[:-4] + '_' + color
                else:
                    outdir = ''

                ROIsnip = crop(bitmap, ROI)
                #fits 1, 2 or 3 gauss spots and determines which one is best
                #returns optimised parameter array
                bestfit, twoIstar, _ = fitNGauss (
                    ROIsnip, options, 
                    DTwoIstar = DTwoIstar, garbageBrightness = garbageBrightness,
                    junkIstar = junkIstar, verbose = False, outdir = outdir
                    )

                #build array containing all the fit results
                locLst[-1][color] = Channel(bitmap)
                locLst[-1][color].fillSpotLst(bestfit)
    
    fpath = os.path.join(wdir, outname)
    
    with open(fpath, 'wb') as output:
        pickle.dump(locLst, output, 1)
        
    return locLst

def plotSinglePair(locLst, pxSize = 10):
    coords = np.zeros([len(locLst),2])
    for i, loc in enumerate(locLst):
        coords[i] = loc['G'].spotLst[0].coord - loc['Y'].spotLst[0].coord
    coords *= pxSize
    coords = kickvector(coords, 50)
    coords = coords - np.mean(coords, axis = 0)
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

def kickvector(dat, maxval):
    return dat[np.linalg.norm(dat, axis = 1) < maxval]

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
    
def plotSinglePair(locLst, pxSize = 10):
    coords = np.zeros([len(locLst),2])
    for i, loc in enumerate(locLst):
        coords[i] = loc['G'].spotLst[0].coord - loc['Y'].spotLst[0].coord
    coords *= pxSize
    coords = coords - np.mean(coords, axis = 0)
    coords = kickvector(coords, 50)
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

def kickvector(dat, maxval):
    return dat[np.linalg.norm(dat, axis = 1) < maxval]

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


def logLikelihood1D(params, func, xdata, ydata):
    return np.sum(func(xdata, *params) - ydata * np.log(func(xdata, *params)))
def ncChidistr(r, mu, sig, A, offset):
    return A * r / sig**2 * np.exp(- (mu**2 + r**2) / (2 * sig**2)) * np.i0( r * mu / (sig**2)) + offset
def expDecay(r, tau, A):
    return A * np.exp( - r / tau )

def sortSpots(loc_in, loc_out):
    
    NG = len(loc_in['G'].spotLst)
    NY = len(loc_in['Y'].spotLst)
    Ndist = min (NG, NY)
    Gcoords = np.zeros([NG, 2])
    Ycoords = np.zeros([NY, 2])
    alldist = np.zeros([NG, NY])
    for i, el in enumerate(loc_in['G'].spotLst):
        Gcoords[i] = el.coord
    for i, el in enumerate(loc_in['Y'].spotLst):
        Ycoords[i] = el.coord
    
    for i in range(NG):
        for j in range(NY):
            alldist[i, j] = np.linalg.norm(Gcoords[i] - Ycoords[j])
    for i in range(Ndist):
        NG, NY = alldist.shape
        Gpos, Ypos = [alldist.argmin() // NY , alldist.argmin() % NY]
        #re-build channel objects with paired spots
        loc_out['G'].spotLst.append(loc_in['G'].spotLst[Gpos])
        loc_out['Y'].spotLst.append(loc_in['Y'].spotLst[Ypos])
        alldist[Gpos] = 1e6
        alldist[:, Ypos] = 1e6
        #alldist = np.delete(alldist, Gpos, 0)
        #alldist = np.delete(alldist, Ypos, 1)


def cropSpot(spotcenter, data, winSigma, ROI):
    """takes ROI at position spotNumber from ROIchannel
    a crop is the taken from datachannel
    ROIchannel and datachannel can be different if e.g. 'Y' ROI
    is used to crop 'R' photons
    crops 2D or 3D data in 2D"""
    #convert between float peak to pixel position
    #check that convertion is correct
    xcenter, ycenter = (np.round(spotcenter) + \
                        ROI[[0,1]]).astype(np.int)
    xstart, ystart, xstop, ystop = [xcenter - winSigma,
                                   ycenter-winSigma,
                                   xcenter + winSigma,
                                   ycenter + winSigma]
    return data[xstart : xstop + 1, ystart : ystop + 1]
    
def fitTau(TAC, TACCal, verbose):
        params0 = [1, 2]
        tactimes = np.arange(len(TAC)) * TACCal
        fitres = minimize(logLikelihood1D, params0, args = (expDecay, tactimes, TAC), 
                          method = 'SLSQP')
        if verbose:
            plt.plot(tactimes, TAC)
            plt.plot(tactimes, expDecay(tactimes, *fitres.x))
            plt.show()
        return fitres.x[0]

def analyseLoc(
    loc, Ggate = 0, Rgate = 0, Ygate = 0, 
    winSigma = 3, framestop = 20, 
    TACCal = 0.0128, verbose = False):
    """finds all closest spot pairs
    for each pair the following are calculate:
        pair distance dist
        pseudoEfficiency E 
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
    outloc = {}
    
    ntacs = 256
    outloc['filepath'] = loc['filepath']
    
    #load image, gate
    CLR = IM.processLifetimeImage(
        outloc['filepath'].encode(), uselines = np.array([1,2]), ntacs = ntacs,
        framestop = framestop
    )
    CLR.loadLifetime()
    CLR.gate(Ggate, 150, channel = 'G')
    CLR.gate(Rgate, 150, channel = 'R')
    CLR.gate(Ygate, 150, channel = 'Y')
    CLR.loadIntensity()


    #get ROI based on Y signal
    #this part should be saved in object at first analysis
    #outloc['ROI'] = getROI(CLR.workIntensity.Y, ROIsize)
    ROI = loc['ROI'] #shorthand
    outloc['ROI'] = loc['ROI']
    
    for color, bitmap in CLR.workIntensity.__dict__.items():
            assert color in ['G', 'R', 'Y']
            ROIsnip = crop(bitmap, outloc['ROI'])
            outloc[color] = Channel(bitmap)
            if verbose:
                plt.imshow(ROIsnip)
                plt.show()
    #add pairwise localisation of loc to outloc
    sortSpots(loc, outloc)
    npairs = len(outloc['G'].spotLst)
    
    outloc['FRETind'] = []
    for i in range(npairs):
        outloc['FRETind'].append(FRETind())
    for i in range(npairs):
        #get No. Green photons
        Gsnip = cropSpot(outloc['G'].spotLst[i].coord, 
                         outloc['G'].bitmap, winSigma, ROI)
        if verbose:
            plt.imshow(Gsnip)
            plt.show()
        Gphotons = np.sum(Gsnip)
        #get No. of Red photons using Y localisation
        Rsnip = cropSpot(outloc['G'].spotLst[i].coord, 
                         outloc['R'].bitmap, winSigma, ROI)
        if verbose:
            plt.imshow(Rsnip)
            plt.show()
        Rphotons = np.sum(Rsnip)
        outloc['FRETind'][i].pseudoE = Rphotons / (Gphotons + Rphotons)
        GTAC = np.sum(cropSpot(outloc['G'].spotLst[i].coord,
                       CLR.workLifetime.G, winSigma, ROI), axis = (0,1))[Ggate:]
        outloc['FRETind'][i].tauG = fitTau(GTAC, 0.128, verbose)
        RTAC = np.sum(cropSpot(outloc['Y'].spotLst[i].coord,
                       CLR.workLifetime.R, winSigma, ROI), axis = (0,1))[Rgate:]
        outloc['FRETind'][i].tauR = fitTau(RTAC, 0.128, verbose)        
        YTAC = np.sum(cropSpot(outloc['Y'].spotLst[i].coord,
                       CLR.workLifetime.Y, winSigma, ROI), axis = (0,1))[Ygate:]
        outloc['FRETind'][i].tauY = fitTau(YTAC, 0.128, verbose)
        Gcoord = outloc['G'].spotLst[i].coord
        Ycoord = outloc['Y'].spotLst[i].coord
        outloc['FRETind'][i].dist = Gcoord - Ycoord
    return outloc
    

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
        folderpath = '{}{}{}'.format(path, ii, '/')
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
                shutil.copy('{}{}'.format(folderpath, files_ptu[i]), '{}{}'.format(savepath, name_new))
