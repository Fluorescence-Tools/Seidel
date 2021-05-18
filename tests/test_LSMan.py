# -*- coding: utf-8 -*-
"""


Created on Fri Oct 23 10:00:46 2020

@author: voort
"""

#%%imports
import numpy as np
import matplotlib.pyplot as plt
import os
import LSManalysis as LSMan
import utility

#this test is not good in it's current state, too much unnneeded loading.
#%% measurement settings
ntacs = 1024
TAC_range = 4096
pulsetime = 50
dwelltime = 20e-6
Nframes = 150
g_factor = 1.18 #checked
shift = 0 #to check
decaytype = 'VM'
rootdir = r'K:\vanderVoortN\FRC\tests\data'    
D0select = 1
intensity_threshold = 0
#%% proces images
identifier = 'D0'
D0dir = os.path.join(rootdir, identifier)
load = False
D0 = utility.loadLSMUtility(D0dir, identifier, g_factor = g_factor,
                    Nframes = Nframes, load = load,
                    intensity_threshold = 20)
#%% fit data

#this norm image business has to go, it is not needed in principle
identifier = 'D0'
normimageG = D0.images['G'][D0select]
normimageY = D0.images['Y'][D0select]
utility.FitPlotLSMUtility(D0, normimageG, normimageY, identifier,
                  fitfunc = 'batchFit2lt')
#%% proces images
identifier = 'DA'
load = False
DAdir = os.path.join(rootdir, identifier)
D1A1 = utility.loadLSMUtility(DAdir, identifier, g_factor = g_factor,
                    Nframes = Nframes, load = load,
                    intensity_threshold = 20)
#%% fit data
identifier = 'DA'
normimageG = D0.images['G'][D0select]
normimageY = D0.images['Y'][D0select]
fitkwargs = LSMan.genDefaultFitKwargs()
fitkwargs['D0dat'] = getattr(D0.images['G'][D0select], decaytype)
fitkwargs['fitrange'] = (30, 380)
fitkwargs['decaytype'] = decaytype
utility.FitPlotLSMUtility(D1A1, normimageG, normimageY, identifier,
                  fitfunc = 'batchFit1ltD0DA', fitkwargs = fitkwargs)

