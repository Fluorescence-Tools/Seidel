# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 17:24:43 2021

@author: voort
"""

#%% test analyzing a ptufile with a set of masks
import LSManalysis as LSMan
import os
import utility
from importlib import reload
reload(LSMan)
#%%
ntacs = 1024
TAC_range = 4096
pulsetime = 50
dwelltime = 20e-6
Nframes = 10 #checked. NOte: pixel size smaller!
g_factor = 1.17 #checked
shift = 0 #to check
decaytype = 'VM'
rootdir = r'Z:\CD95\data\20-4\201112_D1A1_D1A5'
Gpower = 0.85
Ypower = 1.27
intensityThreshold = 20
imreadkwargs = {'ntacs' : 1024,
                        'pulsetime' : 50,
                        'dwelltime': 20e-6,#checked
                        'TAC_range': 4096}

D0select = 1
dataselect = None
#%%
identifier = 'cellClump'
wdir = os.path.join(rootdir, identifier)
cellClump = LSMan.sampleSet(wdir,
                                    g_factor = g_factor,
                                    #dataselect = dataselect,
                                    imreadkwargs = imreadkwargs,
                                    Gpower = Gpower,
                                    Ypower = Ypower,
                                    Nframes = 160)
#%%
maskdirs = os.path.join(wdir, 'masks')
cellClump.analyzeDirwMasks(identifier, maskdirs, [0])
#%%
identifier = 'cellClump'
normimageG = cellClump.images['G'][D0select]
normimageY = cellClump.images['Y'][D0select]
utility.FitPlotLSMUtility(cellClump, normimageG, normimageY, identifier,
                  fitfunc = 'batchFit2lt')
