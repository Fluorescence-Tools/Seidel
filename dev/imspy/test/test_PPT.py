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
import matplotlib.pyplot as plt

import sys
root = r"C:\Users\Niek\work\FRC\dev\imspy"
sys.path.append(os.path.join(root, r"build\Debug"))
import imspy as spy
#%% data file
fname = os.path.join(root, r'test\images\PQSpcm_2019-02-01_15-59-23.ptu').encode()
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
#ImOpts.NumRecords = NumRecords
ImOpts.NumRecords = 1000
ImOpts.linestep = 25e-9
ImOpts.pxsize = 50



GP = spy.imChannel()
GP.line_id = 1
GP.can = [0]
GS = spy.imChannel()
GS.line_id = 1
GS.can = [2]
G = spy.imChannel()
G.line_id = 1
G.can = [0,2]
Channels = [GP, GS, G]
#pass vars as list


#%%
test_image = spy.imspy(list(tac), list(t), list(can), ImOpts)
#%% execute func
#Channels = spy.ProcessPhotonStream(tac, t, can, ImOpts, Channels)
test_image.ProcessPhotonStream()


#%%

test_image.Channels[0].gentacdecay(257)
print(test_image.Channels[0].tacdecay)
