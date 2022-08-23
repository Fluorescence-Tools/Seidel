#Function set to find likely positions of peaks in noisy data
#peak estimate serves as initial parameter for fitting with mul-
#tiple 2D Gaussian.
#Code is single-purpose.
#Tests indicate good performance under noisy conditions
#Under high photon count conditions, the peak finding fails
#if the peaks are not reparated: The raindrop approach collects all
#raindrops in the global maxima.
#use findPeaks function as a single-use function to get peaks from image
#author: Nicolaas van der Voort
#institute: Molecular Physical Chemistry, AG Seidel
#date: January 10, 2020

from scipy.ndimage import gaussian_filter

import numpy as np

def findMaxima(data):
    """finds all local peaks according to 'raindrop' model.
    All pixel get a raindrop. The rain flows uphill in neighbouring pixels
    Priority is given to x-direction
    input: A double gaussian-smoothed input image
    returns: a list of where the raindrops have converged"""
    
    xlen = data.shape[0]
    ylen = data.shape[1]
    points = np.zeros([xlen*ylen, 2], dtype = int)
    eps_step = 1
    for x in range(xlen):
        for y in range(ylen):
            points[x*ylen + y, 0] = x
            points[x*ylen + y, 1] = y
    points_copy = points.copy()
    padded_data = np.pad(data, 1, mode = 'constant')

    while(eps_step != 0):
        for index, point in enumerate(points):
            x,y = point
            xset = padded_data[x: x + 3, y + 1]
            if xset[0] > xset[1] and xset[0] > xset[2]: #move point left
                points_copy[index, 0] -= 1
                continue
            elif xset[2] > xset[1] and xset[2] > xset[0]: # move point right
                points_copy[index, 0] += 1
                continue
            yset = padded_data[ x + 1, y: y + 3]
            if yset[0] > yset[1] and yset[0] > yset [2]: #move up
                points_copy[index, 1] -= 1
            if yset[2] > yset[1] and yset[2] > yset [1]: #move down
                points_copy[index, 1] += 1
        eps_step = np.linalg.norm(points-points_copy)
        points = points_copy.copy()
        
    #peaks is array, each row contains xcoord, ycoord and intensity
    #from finding algorithm, peaks contains many duplicates.
    peak_intensities = np.zeros(xlen * ylen)
    for i in range(xlen * ylen):
        peak_intensities[i] = data[points[i,0], points[i,1]]
    peaks = np.array([points[:, 0], points[:, 1]]).transpose()
    #sort according to intensity
    peaks = peaks[peak_intensities.argsort(-1)][-1::-1]
    return peaks
    
#def findUniquePoints(image, points):
#    """
#    selects unique entries in numpy array points.
#    returns an array of unique points and their corresponding smoothed image
#    intensity. The array is sorted with highed intensity on top"""
#    point_list = list(points)
#    peaks = []
#    while (len(point_list) > 0):
#        xpoint, ypoint = point_list.pop(0) #pop new unique value
#        #add unique value + image intensity to list
#        peaks.append(np.array([xpoint, ypoint, image[xpoint, ypoint]]))
#        #check remaining pointlist for duplicates
#        i = 0
#        while i < len(point_list):
#            #if duplicate: pop
#            if (point_list[i][0] == xpoint and point_list[i][1] == ypoint):
#                point_list.pop(i)
#            #else inspect next element
#            else: i+=1
#        #print(len(point_list))
#    peak_arr = np.array(peaks)
#    peak_arr = peak_arr[peak_arr[:,2].argsort(-1)][-1::-1]
#    return peak_arr
    
def sortPeaks(peaks, xlen, ylen, mindiff):
    """finds the most likely real starting locations of peaks
    First the brightest is taken. Then the second brightest is taken
    that is at least mindiff euclidian norm from the first point. etc.
    If there are many identical peaks, duolicates are filtered out via the 
        mindiff criterium.
    peaks: sorted array of peaks according to intensity
    xlen: xshape of image data
    ylen: yshape of image data
    mindiff: minimal euclidian norm from between peaks.
    returns: 3 peaks"""
    if mindiff == 0:
        print('must set mindiff >0, setting mindigg = 0.1')
        mindiff = 0.1
    i = 0
    goodpeak_cntr = 0
    goodpeaks = np.ones([3,2]) * -10
    #first peak is always the highest
    goodpeaks[0] = peaks[0]
    goodpeak_cntr += 1
    i += 1
    while i < peaks.shape[0]:
        #check if next peak is closer than any known peak
        if (np.linalg.norm(goodpeaks - peaks[i,0:2], axis = 1) < mindiff).any():
            i+=1
        #else another goodpeak is found
        else:
            goodpeaks[goodpeak_cntr] = peaks[i,0:2]
            i+= 1
            goodpeak_cntr += 1
        #stop after 3 peaks
        if goodpeak_cntr >= 3:
            break
    #if for some reason not enough peaks are found, supplement
    while goodpeak_cntr < 3:
        print ('not enough peaks found, try decreasing mindiff \n setting remaining peaks at random')
        goodpeaks[goodpeak_cntr] = np.random.random(2) * [xlen, ylen]
        goodpeak_cntr += 1
    return goodpeaks.astype(int)
    
def findPeaks(data, mindiff = 2, smooth_sigma = 1):
    """ finds three initial estimates for spot locations.
    First the image is smoothed.
    Second all local maxima are found
    Third the highest 3 maxima seperated by mindiff are selected
    returns: 3 x,y entries containing the peak coordinated
    """
    xlen, ylen = data.shape
    data = smooth_image(data, sigma = smooth_sigma)
    peaks = findMaxima(data)
    goodpeaks = sortPeaks(peaks, xlen, ylen, mindiff)
    return goodpeaks
    

def smooth_image(im, sigma = 1):
    return gaussian_filter(im.astype(np.double), sigma)