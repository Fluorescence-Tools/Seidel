import numpy as np
import matplotlib.pyplot as plt
import os
import copy

#%%
def appendOnPattern( wdir, pattern):
    fileList = []
    for file in os.listdir(wdir):
        if file.endswith(pattern):
            fileList.append(file)
    return fileList
    #fnames = [x for x in os.listdir(measdir) if x.endswith('.ptu')]
    
def convertToFullPath(wdir, filelist):
    ffileList = []
    for file in filelist:
        ffileList.append(os.path.join(wdir, file))
    return ffileList
    #ffiles = [os.path.join(measdir, fname).encode() for fname in fnames]

def loadList(ffiles):
    outlist = []
    for ffile in ffiles:
        outlist.append(np.genfromtxt(ffile))
    return outlist

def plotdecay(decay,
              dt, 
              plotkwargs, 
              isnorm = False, 
              normrange = (0, None), 
              isbgsubtract = False,
              bgrange = [0,10],
              ax = plt.gca()):
    """pass a **decaykwargs dict to control background correction and 
    normalisation e.g.:
        decaykwargs = {'isnorm' : True, 'normrange' : (20, 25)}"""
    decay = copy.deepcopy(decay)
    if isbgsubtract:
        avgbg = np.mean(decay[bgrange[0]:bgrange[1]])
        decay -= avgbg
    if isnorm:
        decay = decay / max(decay[normrange[0]:normrange[1]])
    xdat = np.arange(len(decay)) * dt
    ax.plot(xdat, decay, **plotkwargs)
            
#deprecated
def plotList(data_list, dt, plotrange = (0, -1), norm = False, labels = None,
             pattern = '-', normrange = (0, -1), ax = None, alpha = 1,
             clist = [0]):
    #ugly work around, code misplaced, should pass only plotkwargs in this 
    #function, or have separate set of kwargs for bg subtraction etc.
    print('warning! this function is deprecated, and superceded by plotdecay')
    if len(clist) == 1: clist = ['#1f77b4']*len(data_list)
    for i, data in enumerate(data_list):
        datasnip = data[plotrange[0]:plotrange[1]]
        if norm:
            datasnip = datasnip / max(data[normrange[0]:normrange[1]])
        if labels is not None:
            label = labels[i]
        else: label = None
        xdat = np.arange(len(datasnip)) * dt
        if ax:
            ax.plot(xdat, datasnip, pattern, label = label, alpha = alpha,
                    c = clist[i])
        else:
            plt.plot(xdat, datasnip, pattern, label = label, alpha = alpha,
                    c = clist[i])

def normbydecay(data_list, norm_decay, shift = 0, bgrange = None):
    norm_list = []
    #subtract background if specified
    if bgrange is not None:
        norm_decay = norm_decay - np.mean(norm_decay[bgrange[0]:bgrange[1]])
        data_list = [data - np.mean(data[bgrange[0]:bgrange[1]])
            for data in data_list]
    for data in data_list:
        shiftdata, shiftnorm = intshift(shift, data, norm_decay)
        norm_list.append(shiftdata/shiftnorm)
#        if shift <= 0:
#            norm_list.append(data[0:shift-1]/norm_decay[-shift:-1])
#        elif shift > 0:
#            norm_list.append(data[shift:-1]/norm_decay[0:-shift-1])
    return norm_list

def genFr(PS, g_factor, shift = 0, bgrange = None):
    """generate magic angle Fluorescence decay F
    generate anisotropy decay r
    g-factor definition according Lakowicz (Peulen)
    """
    assert len(PS) % 2 == 0
    ntacs = int(len(PS) / 2)
    P = PS[: ntacs]
    S = PS[ntacs:]
    if bgrange is not None:
        [P, S] = [Ch - np.mean(Ch[bgrange[0]: bgrange[1]]) for Ch in [P, S]]
    P, S = intshift(shift, P, S)
    F = P + g_factor * 2 * S
    Fr = P - g_factor * S
    r = Fr / F
    return F, r

def genFrList(data_list, shift = 0, g_factor = 1, bgrange = None):
    rout = []
    Fout = []
    for data in data_list:
        F, r = genFr(data, 1, shift, bgrange = bgrange)
        rout.append(r)
        Fout.append(F)
    return Fout, rout

#####################generic functions########################################

def intshift(shift, data):
    if shift < 0:
        return data[:shift]
    elif shift >= 0:
        return data[shift:]
def intshift(shift, data1, data2):
    if shift < 0:
        return data1[:shift], data2[-shift:]
    elif shift > 0:
        return data1[shift:], data2[:-shift]
    elif shift == 0:
        return data1, data2