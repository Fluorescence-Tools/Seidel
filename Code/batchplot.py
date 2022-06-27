import numpy as np
import matplotlib.pyplot as plt
import os
import copy
import matplotlib as mpl
from scipy.ndimage import gaussian_filter

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
              dtime = 0.064, 
              plotkwargs = {}, 
              normrange = None, 
              bgrange = None,
              ax = None,
              xlim = None):
    """pass a **decaykwargs dict to control background correction and 
    normalisation e.g.:
        decaykwargs = {'isnorm' : True, 'normrange' : (20, 25)}"""
    if ax == None:
        ax = plt.gca()
    decay = copy.deepcopy(decay)
    if bgrange is not None:
        avgbg = np.mean(decay[bgrange[0]:bgrange[1]])
        decay -= avgbg
    if normrange is not None:
        decay = decay / max(decay[normrange[0]:normrange[1]])
    xdat = np.arange(len(decay)) * dtime
    ax.plot(xdat, decay, **plotkwargs)
    if xlim:
        ax.set_xlim(xlim)

def plotdecayList(decaylist, plotdecaykwargs = {}, plotkwargs = {}, ax = None,
            plotout = None):
    """plotdecaykwargs and plotkwargs can be dict of list of dict.
    If dict, all curves are plotted with the same settings
    If list of dict, list must have equal lengths as set1 and each line has a
    set of kwargs"""
    assert type(plotdecaykwargs) == type(plotkwargs) and \
            type(plotkwargs) in [dict, list], \
            "both kwarg arguments must either list or dict"
    if ax == None: 
        fig, ax = plt.subplots(figsize = (11, 8))
    if type(plotkwargs) == dict:
        for decay in decaylist:
            plotdecay(decay, **plotdecaykwargs, plotkwargs = plotkwargs)
    if type(plotkwargs) == list:
        for i in range(len(decaylist)):
            plotdecay(decaylist[i], **plotdecaykwargs[i], 
                      plotkwargs = plotkwargs[i])
    return ax            
#deprecated
#def plotList(data_list, dtime, plotrange = (0, -1), norm = False, labels = None,
#             pattern = '-', normrange = (0, -1), ax = None, alpha = 1,
#             clist = [0]):
#    #ugly work around, code misplaced, should pass only plotkwargs in this 
#    #function, or have separate set of kwargs for bg subtraction etc.
#    print('warning! this function is deprecated, and superceded by plotdecay')
#    if len(clist) == 1: clist = ['#1f77b4']*len(data_list)
#    for i, data in enumerate(data_list):
#        datasnip = data[plotrange[0]:plotrange[1]]
#        if norm:
#            datasnip = datasnip / max(data[normrange[0]:normrange[1]])
#        if labels is not None:
#            label = labels[i]
#        else: label = None
#        xdat = np.arange(len(datasnip)) * dtime
#        if ax:
#            ax.plot(xdat, datasnip, pattern, label = label, alpha = alpha,
#                    c = clist[i])
#        else:
#            plt.plot(xdat, datasnip, pattern, label = label, alpha = alpha,
#                    c = clist[i])
#
#################functionality moved to sampleSet object######################
#def normbydecay(data_list, norm_decay, shift = 0, bgrange = None):
#    norm_list = []
#    #subtract background if specified
#    if bgrange is not None:
#        norm_decay = norm_decay - np.mean(norm_decay[bgrange[0]:bgrange[1]])
#        data_list = [data - np.mean(data[bgrange[0]:bgrange[1]])
#            for data in data_list]
#    for data in data_list:
#        shiftdata, shiftnorm = intshift(shift, data, norm_decay)
#        norm_list.append(shiftdata/shiftnorm)
#    return norm_list

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
#function overloading does not work in python
#def intshift(shift, data):
#    if shift < 0:
#        return data[:shift]
#    elif shift >= 0:
#        return data[shift:]
def intshift(shift, data1, data2):
    """shifts data1 by integer amount w.r.t data2"""
    if shift < 0:
        return data1[:shift], data2[-shift:]
    elif shift > 0:
        return data1[shift:], data2[:-shift]
    elif shift == 0:
        return data1, data2
    
################specific functions taking sampleSet objects ##################
#transform below four functions to take sampleSET and to be part of batchplot

def pltD0DArisetherms(D0Set, DASet, identifier, resdir = None,
                      decaytype = 'VM',
                      plotkwargs = {}, 
                      plotdecaykwargs = {'xlim' : (0, 2),
                                         'normrange' : np.array([0, 30])}):
    fig, ax = plt.subplots(figsize=(11,8))
    D0TACs = D0Set.getDecay(decaytype)
    DATACs = DASet.getDecay(decaytype)
    
    plotdecaykwargs['dtime'] = D0Set.dt_glob
    plotkwargs['c'] = 'g'
    plotdecayList(D0TACs, plotdecaykwargs, plotkwargs, ax = ax)
    plotkwargs['c'] = 'r'
    plotdecayList(DATACs, plotdecaykwargs, plotkwargs, ax = ax)
    
    plt.yscale("log")
    plt.grid()
    plt.xlabel('time (ns)')
    plt.ylabel('normalised intensity')
    plt.title(identifier + ' risetherm')
    plt.plot(0,0, 'g', label = 'Donly (%i cells)' % len(D0TACs))
    plt.plot(0,0, 'r', label = 'D(A) (%i cells' % len(DATACs))
    plt.ylim(1e-3, 1)
    plt.legend()
    if resdir:
        plt.savefig(os.path.join(resdir, identifier + '_risetherms.png'),
                dpi = 300, bbox_inches = 'tight')
    return ax
    
def pltPSrisetherms(sampleSet, identifier, resdir = None,
                    plotkwargs = {}, 
                    plotdecaykwargs = {'xlim' : (0, 2),
                                       'normrange' : np.array([0, 30])}):
    fig, ax = plt.subplots(figsize=(11,8))
    PTACs = sampleSet.getPropertyList('P')
    STACs = sampleSet.getPropertyList('S')
    plt.figure(figsize=(11,8))
    plotdecaykwargs['dtime'] = sampleSet.dt_glob
    plotkwargs['c'] = 'm'
    plotdecayList(PTACs, plotdecaykwargs, plotkwargs, ax = ax)
    plotkwargs['c'] = 'c'
    plotdecayList(STACs, plotdecaykwargs, plotkwargs, ax = ax)
    
    
    plt.yscale("log")
    plt.grid()
    plt.xlabel('time (ns)')
    plt.ylabel('normalised intensity')
    plt.title(identifier + ' PS risetherm')
    plt.plot(0,0, 'm', label = 'Donly P (%i cells)' % len(PTACs))
    plt.plot(0,0, 'c', label = 'Donly S (%i cells)' % len(STACs))
    plt.legend()
    if resdir:
        plt.savefig(os.path.join(resdir, identifier +'_PS_risetherms.png'),
                dpi = 300, bbox_inches = 'tight')

def pltRelativeDecays(sampleSet, identifier, 
                    resdir = None,
                    plotkwargs_base = {'linewidth' : 0.3}, 
                    plotdecaykwargs = {'xlim' : (0, 20),
                                       'normrange' : (23, 30)},
                    decaytype = 'VM',
                    colorcoding = None,
                    vmin = 1,
                    vmax = None, 
                    figsize = (10,7),
                    **kwargs):
    """plots all normal and normalised decays in a sampleset    """
    fig, ax1 = plt.subplots(figsize = figsize)
    TACnorms = sampleSet.getPropertyList(decaytype + 'norm')
    TACnorms = [gaussian_filter(TACnorm, 1) for TACnorm in TACnorms]
    TACs = sampleSet.getPropertyList(decaytype)
    if vmax == None: vmax = max(colorcoding)
    #make a list of colorcoding and set brightness accordingly
    cmap = mpl.cm.get_cmap('plasma')
    cmap.set_over(color = 'green')
    cmap.set_under(color = 'black')
    plotkwargLst = []
    for i in range(len(TACs)):
        #different objects
        plotkwargs = copy.deepcopy(plotkwargs_base)
        if colorcoding is not None:
            c = cmap((np.log(colorcoding[i]) - np.log(vmin)) \
                     / (np.log(vmax)- np.log(vmin)))
            plotkwargs['c'] = c
        plotkwargLst.append(plotkwargs)
    #all pointers to same object
    plotdecaykwargLst = [plotdecaykwargs for i in range(len(TACs))]
    
    for el in plotkwargLst:
        el['ls'] = '--'
        
    plotdecayList(TACnorms, plotdecaykwargLst, plotkwargLst, ax = ax1)
    ax1.plot(0,0, **plotkwargLst[0], label = identifier + ' normalised')
    ax1.set_ylabel('\u03B5 (\u03C4) \'......\'')
    ax1.set_xlabel('time (ns)')
    ax1.grid()
    
    cax = fig.add_axes([0.27, 0.2, 0.5, 0.05])
    norm = mpl.colors.LogNorm( vmin = vmin, vmax = vmax)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')

    #twinx needs to come after plotting to ax1, otherwise it is also 
    #set to logarithmic
    ax2 = ax1.twinx() 
    
    for el in plotkwargLst:
        el['ls'] = '-'
    plotdecayList(TACs, plotdecaykwargLst, plotkwargLst, ax = ax2)
    ax2.plot(0,0, **plotkwargLst[0], 
             label = identifier + ' \u03B5 (\u03C4) - %i cells' % 
             len(TACs))
    ax2.set_ylabel('normalised intensity (a.u.)')
    ax2.set_yscale('log')
    ax2.set_ylim([1e-3,3.8])
    ax1.set_ylim([0,1.2])
    plt.grid()
    plt.legend()
    plt.title('%s %s CD95' % (decaytype, identifier))
    if resdir:
        outname = os.path.join(resdir, identifier + '_epstau+normalised.png')
        plt.savefig(outname, dpi = 300, bbox_inches = 'tight')
    return 1


def pltAnisotropies(sampleSet, identifier, resdir = None,
                    plotkwargs = {}, 
                    plotdecaykwargs = {'xlim' : np.array([0, 25])}
                    ):
    fig, ax = plt.figure(figsize=(11,8))
    rs = sampleSet.getPropertyList('r')
    plt.figure(figsize=(11,8))
    plotdecaykwargs['dtime'] = sampleSet.imreadkwargs.dt_glob
    plotdecayList(rs, plotdecaykwargs, plotkwargs, ax = ax)
    plt.plot([1.6,1.6], [0,0.45], 'k--')
    plt.ylim(0,1)
    plt.ylabel('anisotropy')
    plt.xlabel('time (ns)')
    plt.title(identifier + 'anisotropy')
    plt.grid()
    plt.plot(0,0, **plotkwargs, label = 'anisotropy (%i cells)' % len(rs))
    plt.legend()
    if resdir:
        plt.savefig(os.path.join(resdir, identifier +'_anisotropy.png'),
                dpi = 300, bbox_inches = 'tight')
    return fig, ax