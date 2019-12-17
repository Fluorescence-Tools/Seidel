import ctypes
from ctypes import *
import numpy as np
from precisionFuncs import *
from cpp_wrappers import fit2DGaussian_wrap
params0 = [525,525,100,1e-6,100000]
a = 50
imsize = 1050

im = createGaussImg(imsize / a, params0, a)
#im = np.arange(441)
randparams = []
#for i in range(10000000):
#    im = createGaussImg(imsize / a, params0, a)

for el in params0:
    randparams.append( el*(0.98 +np.random.rand()*0.04))
print(randparams)
fitparams, tIstar = fit2DGaussian_wrap(randparams, a, im)
print(fitparams)
print ('tIstar equals %f'%tIstar)
plt.imshow(im)
plt.colorbar
plt.show()