#helper functions for developmental functions and Gauss fitting
#Nicolaas van der Voort
#AG Seidel, 4 March, 2020

import findPeaksLib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

def pos2ROI(xpos, ypos, winSigma):
    """convert center position and width to ROI"""
    xstart = xpos -winSigma
    xstop = xpos + winSigma +1
    ystart = ypos - winSigma
    ystop = ypos + winSigma + 1
    #ROI is column major
    ROI = np.array([ystart, xstart, ystop, xstop], dtype = np.int).transpose()
    return ROI
    
##########helper functions from Gauss Analysis Pipeline#########################
def cropAnI(image, ROI, ROISize, ROIpad = 0):
    """crop ROI from image
    ROI: roi parameters as taken from AnI, corner positions are taken
    ROISize: length of square ROI from Ani
    ROIpad (optional): enlarge ROI from AnI on all sides
    side of ROI has length: ROISize + 2 * ROIpad"""
    xshape, yshape = image.shape
    cornery, cornerx = ROI[[0, 2]].astype(np.int) - np.array([ROIpad, ROIpad])
    ROISize = ROISize + 2 * ROIpad
    if cornerx < 0 or cornery < 0 or cornerx + ROISize > xshape or cornery + ROISize > yshape:
        raise IndexError
    ROIsnip = image[cornerx: cornerx + ROISize, cornery: cornery + ROISize]
    return ROIsnip
    
def matchfiles(files, roifiles, ext = '_Red Photons.roi'):
    """matches a set of roifiles to their original ptu files based on names.
    Returns a list of filename pairs"""
    filepairs = []
    #copy files to avoid popping original list
    filescopy = files[:]
    for roifile in roifiles:
        for i in range(len(filescopy)):
            if filescopy[i][:-4] == roifile[:-len(ext)]:
                filepairs.append([filescopy.pop(i),roifile])
                break
    return filepairs
    

def genROI(im, ROIsize):
    """get single ROI around intensity max of image"""
    assert ROIsize % 2 == 0
    ROIside = ROIsize / 2
    smooth_im = findPeaksLib.smooth_image(im)
    xlen, ylen = im.shape
    maxpos = smooth_im.argmax()
    xpos, ypos = [maxpos // xlen, maxpos % xlen]
    ROI = np.array([xpos - ROIside, ypos - ROIside, xpos + ROIside, ypos + ROIside]).astype(np.int)
    return ROI

def crop(image, ROI):
    """crop ROI from image
    ROI: roi parameters as taken from getROI
    returns cropped image"""
    xshape, yshape, *_ = image.shape #*= is wildcard in case of lifetime image
    if ROI[0] < 0 or ROI[1] < 0 or ROI[2] >= xshape or ROI[3] >= yshape:
        raise IndexError
    return image[ROI[0]: ROI[2], ROI[1]: ROI[3]]
    
    
########################plotting functions#####################################
def plotBitmapROI(bitmap, spotLst):
    fig,ax = plt.subplots(1)
    ax.imshow(bitmap)
    for spot in spotLst:
        rect = patches.Rectangle((spot.ystart, spot.xstart), 
                                spot.ystop - spot.ystart,
                                spot.xstop - spot.xstart,
                                linewidth = 1, edgecolor = 'r', facecolor='none')
        ax.add_patch(rect)
        ax.plot(spot.posx, spot.posy, 'r.')
    plt.show()