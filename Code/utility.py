#A collection of utility functions
#functions have some common use, but are not intended to be maintained
#author: Nicolaas van der Voort
#AG Seidel, HHU
#14 April 2020

import ImageManipulation as IM
import numpy as np
import os


def exportLSMTAC(fname, outdir, dwelltime, pulsetime, uselines = np.array([1]), 
                 Gchan = np.array([0,1]), Rchan = np.array([4,5]), Ychan = np.array([4,5]),
                ntacs = 1024, TAC_range = 4096, PIE_gate = 440):
    """utility function to export GRY (P+S) TAC histogram for LSM microscope"""
    _, file = os.path.split(fname)
    data = IM.processLifetimeImage(fname.encode(), uselines = uselines, Gchan = Gchan, 
                                   Rchan = Rchan, 
                                   Ychan = Ychan, ntacs = ntacs, TAC_range = TAC_range, \
                                   pulsetime = pulsetime, dwelltime = dwelltime)
    data.loadLifetime()
    GTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'G'))).transpose()
    data.gate(0,PIE_gate, channel = 'R')
    data.gate(PIE_gate, -1, channel = 'Y')
    RTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'R'))).transpose()
    YTACS = np.stack((np.arange(ntacs), data.getTACS(mode = 'Y'))).transpose()
    np.savetxt(os.path.join(outdir, file[:-4] + '_G.dat'), GTACS, fmt = '%i')
    np.savetxt(os.path.join(outdir, file[:-4] + '_R.dat'), RTACS, fmt = '%i')
    np.savetxt(os.path.join(outdir, file[:-4] + '_Y.dat'), YTACS, fmt = '%i')