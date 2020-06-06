# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 15:42:39 2020

@author: voort
"""

#%% imports

import cpp_wrappers
import numpy as np
import os
from importlib import reload

import sys.path
sys.path.append("K:\vanderVoortN\FRC\dev\imspy\build\Debug\")
import imspy as spy
#%%
reload(spy)

#%%build channels
fname = r'K:\vanderVoortN\FRC\dev\imspy\test\images\PQSpcm_2019-02-01_15-59-23.ptu'.encode()
#%% gen header and build channels
assert type(fname) == bytes, 'Error, fname is not bytes type'

NumRecords = cpp_wrappers.ptuHeader_wrap (fname)
eventN, tac, t, can = cpp_wrappers.ptu_wrap(fname, NumRecords)
root, file = os.path.split(fname)
name, _ = os.path.splitext(file)
header_name = os.path.join(root, b"header", name + b".txt")
print('number of records is ' + str(NumRecords))

dimX, dimY, dwelltime, counttime = cpp_wrappers.read_header(header_name)
    
ImOpts = spy.imOpts()
ImOpts.line_ids = [1, 2]
ImOpts.dwelltime = dwelltime
ImOpts.counttime = counttime
ImOpts.NumRecords = 100
ImOpts.linestep = 25e-9
ImOpts.pxsize = 50

G = spy.imChannel()
G.line_id = 1
G.can = [0, 3]

#%% execute func
spy.ProcessPhotonStream(tac, t, can, ImOpts, [G])
print(G.phstream)

#%%
a = np.ones(2)
b = np.ones(2)
c = np.ones(2)
Ph1 = spy.Eigen_array(a, b, c)
print(Ph1.tac)