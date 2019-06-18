# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 17:31:02 2019

@author: voort
"""
import numpy as np
import numpy.fft as ft
import matplotlib.pyplot as plt
from sys import float_info
import os
from scipy.optimize import curve_fit
from cpp_wrappers import *
import ctypes
from PIL import Image

def readtxt(fileA, fileB):
    imA = np.genfromtxt(fileA)
    imB = np.genfromtxt(fileB)
    return imA, imB
	
def readfolder(filedir):
    fnames = os.listdir (filedir)
    frame = np.genfromtxt(os.path.join(filedir, fnames[0]))
    #merge two tuples to generate the desired img shape
    imgshape = ( (len(fnames),) + frame.shape)
    img = np.zeros(imgshape)
    for i, name in enumerate(fnames):
        frame = np.genfromtxt(os.path.join(filedir, name))
        img[i] = frame
    return img
	
def readraw(fileA, fileB, w, h):
    f = open(fileA, "rb")
    imA = np.fromfile(f, dtype = np.uint16).reshape(400,400)
    f.close()
    f = open(fileB, "rb")
    imB = np.fromfile(f, dtype = np.uint16).reshape(400,400)
    f.close()
    return imA, imB

def Hamming(w,h):
    alpha = 0.54
    beta = 1-alpha
    xv = alpha - beta * np.cos(2*np.pi / (w-1) * np.arange(w))
    yv = alpha - beta * np.cos(2*np.pi / (h-1) * np.arange(h))
    hamming = np.zeros([w,h])
    for i in range(h):
        hamming[i] = xv * yv[i]
    return hamming

def perform_FRC(fftA, fftB, h, w, theta):
    #int cast floors values
    xc = int(w / 2)
    yc = int(h / 2)
    rMax = min(w - xc, h - yc)
    #init return values
    smallAngles = np.zeros(rMax)
    largeAngles = np.zeros(rMax)
    threeSigma = np.zeros(rMax)
    fiveSigma = np.zeros (rMax)
    
    #init local vars
    corr_largeAngles = np.zeros(rMax, dtype = 'complex128')
    absA_largeAngles = np.zeros(rMax)
    absB_largeAngles = np.zeros(rMax)
    corr_smallAngles = np.zeros(rMax, dtype = 'complex128')
    absA_smallAngles = np.zeros(rMax)
    absB_smallAngles = np.zeros(rMax)
    
    #loop over all pixels and update local vars
    for x in range(w):
        for y in range(h):
            #r is the ring number that the pixel is in
            r =  int(round(np.sqrt((x-xc)**2+(y-yc)**2)))
            #we include only full rings in the image, no truncated rings
            if (r < rMax):
                threeSigma[r] += 1
                fiveSigma[r] += 1
               # print("value of theta is %f" %theta)
             #   print("angle of data point equals: %f" % np.arctan(abs(y-yc)/ (abs(x - xc) + float_info.epsilon)))
                if (theta == 0 or np.arctan(abs(y-yc)/ (abs(x - xc) + float_info.epsilon)) > theta ):
                    corr_largeAngles[r] += fftA[y, x] * np.conj(fftB[y, x])
                    absA_largeAngles[r] += abs(fftA[y,x]**2)
                    absB_largeAngles[r] += abs(fftB[y,x]**2)
                else:
                    corr_smallAngles[r] += fftA[y, x] * np.conj(fftB[y, x])
                    absA_smallAngles[r] += abs(fftA[y,x]**2)
                    absB_smallAngles[r] += abs(fftB[y,x]**2)
    
    #compute class vars
    largeAngles = abs(corr_largeAngles) / np.sqrt(absA_largeAngles*absB_largeAngles + float_info.epsilon)
    smallAngles = abs(corr_smallAngles) / np.sqrt(absA_smallAngles*absB_smallAngles + float_info.epsilon)
    threeSigma = 3 / np.sqrt(threeSigma / 2)
    fiveSigma = 5 / np.sqrt(fiveSigma / 2)
    
    #set values higher than 1 to 1.
    for i in range(rMax):
        if (threeSigma[i] > 1):
            threeSigma[i] = 1
        if (fiveSigma[i] > 1):
            fiveSigma[i] = 1

    return smallAngles, largeAngles, threeSigma, fiveSigma

def meanFilter(data, width):
    d = data.shape[0]
    smoothData = np.zeros(d)
    for i in range(d):
        if ( i <= width ):
            smoothData[i] = np.mean(data[:i + width + 1])
        elif (i > width and i < d - width):
            smoothData[i] = np.mean(data[i - width : i + width + 1])
        elif ( i >= d - width):
            smoothData[i] = np.mean(data[i - width:])
        else:
            print("error found in meanFilter program")
    return smoothData




def pltFRC(pixelSize, title, threeSigma, fiveSigma, largeAngles, smallAngles, saveas = None):
    d = threeSigma.shape[0]
    #plot smallAngles if array was filled
    if (smallAngles[0] is not None):
        plt.plot(smallAngles, '.-', label = 'small angle')
    plt.plot(largeAngles, '.-', label = 'large Angle')
    plt.plot(threeSigma, label = '3 sigma')
    plt.plot(fiveSigma, label = '5 sigma')
    plt.plot(np.ones(d) / 7, label = '1/7')
    plt.legend()
    plt.title('FRC for ' + title)
    plt.ylabel('correlation')
    plt.xlabel('resolution (nm)')
    plt.gca().set_xlim(left = 0)
    #set custom x ticks
    locs, labels = plt.xticks()
    
    for i in range(locs.shape[0]):
        if locs[i] == 0:
            labels[i] = 'inf'
        else:
            labels[i] = (2 * d * pixelSize / locs[i]).astype(np.int)
    plt.xticks(locs, labels)
    if (saveas is not None):
        plt.savefig(saveas, dpi = 300)
    plt.show()
    
def checkArguments(imA, imB, theta):
    assert imA.shape == imB.shape, "images of unequal size"
    assert theta <= np.pi / 2 and theta >= 0,"invalid theta: theta < 0 or theta > pi / 2"
    
def findIntercept(data, threshold, pixelSize):
    assert( data.shape == threshold.shape), 'data and threshold of unequal length'
    d = data.shape[0]
    intercept = -1
    resolution = -1
    for i in range(3,d-1):
        if (data[i] < threshold[i] and max(data[i-1], data [i-2], data[i-3]) > threshold[i]):
            resolution = 2 * d * pixelSize / (i)
            intercept = i + 1
            break
    return resolution, intercept

def FRCAnalysis(imA, imB, pixelSize, theta = 0, meanFilterwidth = 3, title = 'someTitle', correctDrift = True):
    checkArguments(imA, imB, theta)
    (w,h) = imA.shape
    hamming = Hamming(w,h)
    fftA = ft.fftshift(ft.fft2(imA*hamming))
    fftB = ft.fftshift(ft.fft2(imB*hamming))
    
    if correctDrift:
        #detect drift of imB
        drift, _ = DriftDetect(imA, imB, pixelSize)
        #correct drift
        fftB = fftShift(fftB, drift)
        print("applied drift correction of %f nm in x and and %f nm in y.\n" 
              % (drift[0] * pixelSize, drift[1] * pixelSize))
        
    smallAngles, largeAngles, threeSigma, fiveSigma = perform_FRC(fftA, fftB, h, w, theta)
    sSmallAngles = meanFilter(smallAngles, meanFilterwidth)
    sLargeAngles = meanFilter(largeAngles, meanFilterwidth)
    pltFRC(pixelSize, title, threeSigma, fiveSigma, sLargeAngles, sSmallAngles)
    
    LargeAnglesResolution = np.ones([3,2]) * -1
    smallAnglesResolution = np.ones([3,2]) * -1
    fixedthreshold = np.ones(smallAngles.shape[0])/7.
    
    LargeAnglesResolution[0] = findIntercept(sLargeAngles, fixedthreshold, pixelSize)
    LargeAnglesResolution[1] = findIntercept(sLargeAngles, threeSigma, pixelSize)
    LargeAnglesResolution[2] = findIntercept(sLargeAngles, fiveSigma, pixelSize)
    
    smallAnglesResolution[0] = findIntercept(sSmallAngles, fixedthreshold, pixelSize)
    smallAnglesResolution[1] = findIntercept(sSmallAngles, threeSigma, pixelSize)
    smallAnglesResolution[2] = findIntercept(sSmallAngles, fiveSigma, pixelSize)
    
    print("retrieved resolution for " + title + " is: " + f"{LargeAnglesResolution[0,0]:.0f}nm\n" )
    
    return sSmallAngles, sLargeAngles, threeSigma, fiveSigma, smallAnglesResolution, LargeAnglesResolution

def fftShift(fftImage, shift):
    """multiply fft Image with complex exponential such that image is shifted by shift.
    
    Usage Warning: do not apply to raw ffts. Either apply Hamming Window or pad with zeros."""
    #this function needs some work to enable shifting a stacks of images.
    (w,h) = fftImage.shape
    x = np.arange(- np.floor( w/2), - np.floor(w/2) + w)
    y = np.arange(- np.floor( h/2), - np.floor(h/2) + h)
    [xF, yF] = np.meshgrid(x, y)
    
    return fftImage * np.exp(1j * 2 * np.pi * (shift[0] * xF / w + shift[1] * yF / h))

	
def expDecay(x, a, b, c):
    return a * np.exp(-x / b) + c

def twoD_Gaussian(coord, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x = coord[0]
    y = coord[1]
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def plotBleaching(frameint, popt, pcov, frameN):
    plt.plot(frameint, label = 'data')
    fit = np.zeros(len(frameint))
    for i in frameN:
        fit[i] = expDecay(i, popt[0], popt[1], popt[2] )
    plt.plot(fit, label = 'fit')
    plt.xlabel('frame number')
    plt.ylabel('frame intensity')
    plt.title('Photobleaching at 100% STED power')
    plt.legend()
    plt.show()
    print (u'decay half-life is {0:.2f} \u00B1 {1:.2f} frames'.format(np.log(2) * popt[1],np.sqrt(pcov[1,1]) * np.log(2)))

def fitBleaching(frameint, p0 = [10000, 1, 30000], showFit = False ):
    frameN = np.arange(len(frameint))
    popt, pcov = curve_fit(expDecay, frameN, frameint, p0 = p0)
    if showFit: plotBleaching(frameint, popt, pcov, frameN)
    return popt, pcov
	
def DriftDetect(imA, imB, pixelSize, p0 = (-1, -1, -1, -1, -1, -1, -1), verbose = False, fitWindow =10):
    """Detect drift between imA and imB. 
     Programm uses autocorrelation and least squares fitting from scipy's curve_fit func.
     p0 is used as initial guess. autocorrelation landscape shows local minima when beads are used.
     Therefore, the initial guess must be accurate, especially the center position guess.
     If verbose is True, the autocorrelation, fitted gaussian center
     and zero shift center are plotted. Also the shift will be printed in nm.
     Returns: center, centersigma in number of pixels."""
    assert (imA.shape ==imB.shape and imA.shape[0]==imA.shape[1]), "images non-square or non-equal sized"
    imcenter = int(imA.shape[0]/2)

    
    #perform autocorrelation
    corr = np.fft.fftshift(np.fft.ifft2(np.fft.fft2(imA) * np.conj(np.fft.fft2(imB))))
    #take center region for Gaussian fitting
    corr = corr[imcenter - fitWindow : imcenter + fitWindow, imcenter - fitWindow : imcenter + fitWindow]
    
    # Create x and y indices
    x = np.arange(corr.shape[0])
    y = np.arange(corr.shape[0])
    coord = np.meshgrid(x, y)
    
    #if parameter guess is unset, make intelligent estimate.
    if p0[0] == -1:
        p0 = (np.max(abs(corr)), fitWindow, fitWindow, 1, 1, 0, 0)
    
    
    #perform 2D gauss fit
    popt, pcov = curve_fit(twoD_Gaussian, coord, abs(corr).ravel(), p0=p0)
    
    #take relevant parameters
    drift = popt[1:3] - [fitWindow, fitWindow]
    driftSigma = np.sqrt([pcov[1,1], pcov[2,2]])
    
    #plt and print
    if verbose:
        plt.imshow(abs(corr))
        plt.plot([0,fitWindow*2],[fitWindow,fitWindow], 'r', linewidth = 0.5)
        plt.plot([fitWindow, fitWindow], [0, fitWindow * 2], 'r', linewidth = 0.5)
        plt.colorbar()
        plt.scatter(popt[1], popt[2], s=10)
        plt.show()
        
        print ("using initial parameters (Amplitude, xpos, ypos, xsigma, ysigma, theta, background) equal to " + str(p0))

        print (u"shift in x is  {:.3f} \u00B1 {:.3f}  pixels\n".format(drift[0], driftSigma[0])
               + u"shift in y is  {:.3f} \u00B1 {:.3f} pixels\n".format(drift[1], driftSigma[1]) 
               + u"total shift is {:.3f} \u00B1 {:.3f} nm for pixel size of {:d} nm \n".format(
                   np.linalg.norm(drift-[fitWindow, fitWindow])*pixelSize,
                   np.linalg.norm(driftSigma*pixelSize),
                   pixelSize
               ) 
              )
    return drift, driftSigma

def stackDrift(img, pixelSize, dwelltime, verbose = False):
    """calculates cumulative drift for frames in stack. 
    Drift is calculated wrt the preceding frame. Then cumulative drift is taken to represent total drift.
    pixelSize in nm
    dwelltime in microseconds
    if verbose is True, a plot is generated."""
    frames = img.shape[0]
    #first element of driftArr is 0
    drift = np.zeros((frames,2))
    for i in range(1, frames):
        drift[i], _ = DriftDetect(img[i-1], img[i], pixelSize)
    cumDrift = np.cumsum(drift, axis = 0)*pixelSize
    if verbose:
        xdim = img.shape[1]
        ydim = img.shape[2]
        #the factor 2.5 accounts coarsely for microscope dead time. 
        colors = xdim* ydim* dwelltime* 2.5* np.arange(0,frames,1) / 1e6
        plt.figure(figsize=(6,4.5))
        plt.scatter(cumDrift[:,0], cumDrift[:,1], marker = '.', c = colors, cmap = 'hot')
        plt.ylabel("drift in y (nm)")
        plt.xlabel("drift in x (nm)")
        plt.title("stack drift opver")
        cbar = plt.colorbar()
        cbar.set_label("time in seconds")
        plt.show()
    return cumDrift

def createBeadImg(size, spotLst, spotSigma, spotAmplitude, background, UsePoisson = True):
    """build an artificial bead image. 
    The PSF is modelled with a symmetric 2D gaussian with uniform sigma and amplitude.
    Spots will located at spotLst. Then background and poisson noise is added.
    Parameters:
    size : number of pixels in x and y. non-square images not seported. int
    spotLst : list of spot center. Can be generated using createSpotLst. Nx2 np array.
    Nspot : number of spots, int
    spotSigma : uniform standard deviation of gaussian spots, float
    spotAmplitude: uniform amplitude of gaussian spot. float
    background : constant background added to image. float"""

    #create image canvas
    im = np.zeros((size, size))

    #add spots
    x = np.arange(size)
    y = np.arange(size)
    coord = np.meshgrid(x, y)
    theta = 0 #in degrees
    offset = 0
    for spot in spotLst:
        spot = twoD_Gaussian(coord, spotAmplitude, spot[0], spot[1], spotSigma, spotSigma, theta, offset)
        im += spot.reshape(size,size)
    #add backgound
    im += background
    #add Poisson noise
    if UsePoisson:
        im = np.random.poisson (im)

    return im.astype(np.int)

def createSpotLst(size, Nspots, spotSigma, avoidEdges = False, avoidNeighbour = False):
    """create gaussian spot center list. 
    avoidEdges : if True, all spotcenters are at least 4 * spotSigma from image edge. bool
    avoidNeighbour : if True, all spotcenters are at least 8 * spotSigma apart. bool"""
    #create spot center list
    spotLst = np.zeros((Nspots,2))
    i=0
    counter = 0
    while i < Nspots:
        spotLst[i,0] = np.random.random()*size
        spotLst[i,1] = np.random.random()*size

        #check if bead is close to edge
        if avoidEdges:
            legalEdge = (spotLst[i,0] > 4 * spotSigma and spotLst[i,0] < size - 4 * spotSigma and
                spotLst[i,1] > 4 * spotSigma and spotLst[i,1] < size - 4 * spotSigma)
        else: legalEdge = True

        #check if bead is too close to other bead
        if avoidNeighbour:
            legalNeighbour = True
            for j in range(i):
                dsq = (spotLst[i,0] - spotLst[j,0])**2 + (spotLst[i,1] - spotLst[j,1])**2
                if dsq < (8 * spotSigma)**2:
                    legalNeighbour = False
                    break
        else: legalNeighbour = True

        #Repeat untill both conditions are met. Raise error after 100 false placements

        if (legalEdge and legalNeighbour): 
            i+=1
            counter = 0
        if counter >100: raise Exception('image too crowded: no space for new spots')
        else: counter +=1
    return spotLst

def genImABfromptu(fname, uselines = np.ones(1, dtype = np.int), xbinning = 1, ybinning= 1, gate = 3228):
    NumRecords = ptuHeader_wrap (fname)
    eventN, tac, t, can = ptu_wrap(fname, NumRecords)
    root, file = os.path.split(fname)
    name, _ = os.path.splitext(file)
    header_name = os.path.join(root, b"header", name + b".txt")
    print('number of records is ' + str(NumRecords))

    dimX, dimY, dwelltime, counttime = read_header(header_name)
    dimX = int(dimX / xbinning)
    dimY = int(dimY / ybinning)
    dwelltime *= xbinning
    imA, imB = SplitOnTacs_wrap(eventN, tac, t, can, dimX, dimY, dwelltime, counttime, 
                                NumRecords, gate = gate, uselines = uselines)
    print("total image intensity is "  + str(np.sum(imA+imB)))

    try:
        os.mkdir(os.path.join(root, file[:-4]))
        print('storing result in new folder')
    except:
        print('overwriting result in existing folder')

    im = Image.fromarray(imA)
    outname = str(os.path.join(root, file[:-4] + b'\\imA.tiff'), encoding='utf-8')
    im.save(outname)
    im = Image.fromarray(imB)
    outname = str(os.path.join(root, file[:-4] + b'\\imB.tiff'), encoding='utf-8')
    im.save(outname)
    return imA, imB

def genArrfromTif(fname):
    return np.array(Image.open(fname))