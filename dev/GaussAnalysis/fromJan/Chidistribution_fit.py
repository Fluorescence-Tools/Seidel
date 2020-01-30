# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 13:53:55 2019

@author: buddeja
"""

import os
import shutil
import matplotlib.pyplot as pp
from PIL import Image
import numpy as np
from tkinter import*
from PIL import Image,ImageTk
from tkinter import ttk
from natsort import natsorted


import numpy as np
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as pp
import math
import os
from scipy.stats import stats
from skimage import img_as_float
from PIL import Image
from scipy import ndimage
from skimage import feature
import numpy.ma as ma
from skimage.feature import peak_local_max
import shutil
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg#, NavigationToolbar2TkAgg
from skimage import color
import matplotlib.cm as cm
from sklearn.neighbors import NearestNeighbors
from astropy.stats import RipleysKEstimator


import numpy as np
import csv
import pandas as pd
import io
from io import StringIO 
import scipy.stats
from scipy.spatial import distance
from astropy.stats import RipleysKEstimator
from sklearn.neighbors import NearestNeighbors
from functools import partial
import time
import seaborn as sns
import math 

import scipy
import scipy.stats
from scipy.optimize import curve_fit
#from scipy.stats import chi2
from scipy.stats import chisquare
from lmfit import Model
from lmfit import conf_interval2d
from mpl_toolkits import mplot3d
from matplotlib import ticker, cm
import matplotlib.colors as colors
from numpy import inf
import lmfit
import emcee
import corner




filepath = 'N:/Singlem/singlem20-1/January/28_dBlowFRET_JHB/3rd_run/Analysis_JHB/export_JHB/Imagewise Colocalization Fit results/'


ending_01 = '.0I4'
ending_02 = '.II4'

#Fit paramter:
Fit_only =                 False
Calc_uncertainty =         False
Plot_init =                False
Parameter_scan =           False
Probability_scan =         False
Simuluation =              False
Gaussian_indiv_plot =      True
Plot_Overview_statistics = False
Add_SgSrSyStoichiometry =  False

# Parameterscan :

Parameter_scan1       = 'R1' # type False for not parameter scan else type i.e. 'mu1'
scan_value1           = np.linspace(0.01,20,20)
Parameter_scan2       = 'sig1' # type False for not parameter scan else type i.e. 'mu1'
scan_value2           = np.linspace(0.01,10,20)
Value_displayed       = 'redchi'   # chi² = chisqr, chi²_red = redchi
Value_displayed_max   = 10
Fitfunction_used      = 'n.c.gaussian2_2D'

Binning = 31
pixelsize = 10 # in nm

shift_x = 0.7441059381107424 # pixelcoorection in x
shift_y = 0.397 # pixelcoorection in y

Fitrange =[0,50]
weighting = 1#  'Poisson'
Method = 'leastsq'#'leastsq'

Init_value = [    10,         4.66,           100,        0.00001,               4,              0,           0.0001]
Bound_min  = [     0,         0.5,              0,              0,             0.5,              0,                0]
Bound_max  = [    20,          20,           1000,              5,              10,           1000,                1]
Value_fit  = [  True,        False,          True,          False,            True,          False,            False]
Value_exp =  [  None,        None,           None,           None,            None,           None,             None]




filter_0I4 = ["data_green['ellipticity'] > 0",
              "data_green['ellipticity'] < 2",
              "data_green['lm_message_cI'] > 0",
              "data_green['lm_message_cI'] <= 4"]


filter_II4 =["data_yellow['ellipticity'] > 0.5",
             "data_yellow['ellipticity'] < 1.5",
             "data_yellow['lm_message_cII'] > 0",
             "data_yellow['lm_message_cII'] <= 4"]
 



###############################################################         





files = os.listdir(filepath)
files = natsorted(files)
file_folder_0I4 = []
file_folder_II4 = []
file_merge = []

for i in files:
    if i.endswith(ending_01):
        file_folder_0I4.append(i)
    elif i.endswith(ending_02):
        file_folder_II4.append(i)
file_folder_0I4 = natsorted(file_folder_0I4)
file_folder_II4 = natsorted(file_folder_II4)
##############################################################################


if Fitfunction_used == 'n.c.gaussian2_2D':
    def Fitfunction(x, R1, sig1, A1, R2, sig2, A2, offset):
        return  A1*x/(sig1*sig1)*  np.exp(-(x*x + R1*R1)/(2 * sig1 * sig1))*np.i0(x*R1/(sig1*sig1)) + A2*x/(sig2*sig1)*  np.exp(-(x*x + R2*R1)/(2 * sig2 * sig2))*np.i0(x*R2/(sig2*sig2)) + offset

elif Fitfunction_used == 'n.c.gaussian3_2D':
    def Fitfunction(x, R1, sig1, A1, R2, sig2, A2, R3, sig3, A3, offset):
        return  A1*x/(sig1*sig1)*  np.exp(-(x*x + R1*R1)/(2 * sig1 * sig1))*np.i0(x*R1/(sig1*sig1)) + A2*x/(sig2*sig1)*  np.exp(-(x*x + R2*R1)/(2 * sig2 * sig2))*np.i0(x*R2/(sig2*sig2)) + offset  

elif Fitfunction_used == 'n.c.gaussian2_3D':
    def Fitfunction(x, R1, sig1, A1, R2, sig2, A2, offset):
        r'f(x) = $\sum_{i=1}^2\frac{A_ix}{2.5066*R_i \sigma_i}\{e^{-(\frac{x-R_i}{2\sigma_i})^{2}} + e^{-(\frac{x+R_i}{2\sigma_i})^{2}}\} + offset$'
        return ((A1*x)/(np.sqrt(2*math.pi)*R1*sig1))*  (np.exp(-np.power(x - R1, 2.) / (2 * np.power(sig1, 2.)))+ np.exp(-np.power(x + R1, 2.) / (2 * np.power(sig1, 2.)))) + ((A2*x)/(np.sqrt(2*math.pi)*R2*sig2))*  (np.exp(-np.power(x - R2, 2.) / (2 * np.power(sig2, 2.)))+ np.exp(-np.power(x + R2, 2.) / (2 * np.power(sig2, 2.)))) + offset
    

def getdata(file_folder): 
    
    file_merge = []
    file_name = []
    x_coord = []
    y_coord = []
    y_offset = []
    
    for i in range(len(file_folder)):
        file = pd.read_csv('{}{}{}'.format(filepath, '/', file_folder[i]),sep='\t')
        
        if len(file)== 0 or len(file)> 1:
            file_values = np.empty(len(list(file)))
            file_values[:] = np.nan
            file_merge.append(file_values)
            file_name.append(file_folder[i])
        elif len(file)== 1:
            for ii in range(len(file)):
                file_values = file.values[ii]
                file_merge.append(file_values)
                file_name.append(file_folder[i])
        s = file_folder[i]
        x_coord.append(s[s.find('_x_')+len('_x_'):s.rfind('_y_')])
        y_coord.append(s[s.find('_y_')+len('_y_'):s.rfind('_G')])
        y_offset.append(s[s.find('s_y')+len('s_f'):s.rfind('_s')])
        
                
                
    data = pd.DataFrame(np.asarray(file_merge), columns = list(file))

    data['data_name'] =  file_name
    data['x_corr']=  x_coord
    data['y_corr']=  y_coord
    data['y_offset']=  y_offset
    
    data= data.dropna()
    header = list(file)
    header.append('data_name')
    header.append('x_corr')
    header.append('y_corr')
    header.append('y_offset')

    return data, header, file_name

def calculateDistance(x1, x2, y1, y2):
    dist_min = []
    d_test = []
    x_min1 = []
    y_min1 = []
    x_min2 = []
    y_min2 = []
    
    for i in range(len(x1)):
        d = []
        for ii in range(len(x2)):
            #dist = np.sqrt((x2[ii] - x1[i])**2 + (y2[ii] - y1[i])**2)
            d.append(np.sqrt(((x2[ii]-x1[i])**2)+((y2[ii]-y1[i])**2))) 
            
        dist_min.append(np.amin(d)*pixelsize)

       # d_test.append(np.amin(d)*pixelsize)
        d_test.append(d[np.where(d == np.amin(d))[0][0]])
        x_min1.append(x1[i])
        y_min1.append(y1[i])
        x_min2.append(x2[np.where(d == np.amin(d))[0][0]])
        y_min2.append(y2[np.where(d == np.amin(d))[0][0]])
        coord_merge =np.vstack((x_min1,y_min1, x_min2, y_min2, dist_min))
        
    return dist_min, coord_merge


###############################################################################  
# Imreading data from filepath:
if Fit_only == False:
    data_0I4 =   getdata(file_folder_0I4)[0]
    header_0I4_old = list(data_0I4)
    data_II4 =   getdata(file_folder_II4)[0]
    header_II4_old = list(data_II4)
    
    
    header_0I4 =[]
    header_II4 =[]
    for i in range(len(header_0I4_old)):
        n = header_0I4_old[i].replace(' ','')
        n = n.replace('(','_')
        n = n.replace(')','_')
        n = n.replace('<','')
        n = n.replace('>','')
        n = n.replace(',','_')
        header_0I4.append(n) 
        
    for i in range(len(header_II4_old)):
        n = header_II4_old[i].replace(' ','')
        n = n.replace('(','_')
        n = n.replace(')','_')
        n = n.replace('<','')
        n = n.replace('>','')
        n = n.replace(',','_')
        header_II4.append(n) 
    
    data_0I4 = pd.DataFrame((data_0I4).values, columns = header_0I4)
    data_II4 = pd.DataFrame((data_II4).values, columns = header_II4)
    
    x1 = np.asarray(data_0I4.peak_x_cI)
    y1 = np.asarray(data_0I4.peak_y_cI)
    x2 = np.asarray(data_II4.peak_x_cII)
    y2 = np.asarray(data_II4.peak_y_cII)

    d_x =[]
    d_y =[]
    for i in range(len(x1)):
        d = []
        for ii in range(len(x2)):
            d.append(np.sqrt(((x2[ii]-x1[i])**2)+((y2[ii]-y1[i])**2)))
        d_x.append(x1[i]-  (x2[np.where(d == np.amin(d))[0][0]]))
        d_y.append(y1[i] - (y2[np.where(d == np.amin(d))[0][0]]))
    fig1 = pp.figure()
    pp.scatter(d_x,d_y)
    d_x_new = list(np.array(d_x)[np.where(np.sqrt(np.array(d_x)**2+ np.array(d_y)**2) < 10 )])
    d_y_new = list(np.array(d_y)[np.where(np.sqrt(np.array(d_x)**2+ np.array(d_y)**2) < 10 )])
    
    fig2 = pp.figure()
    pp.scatter(d_x_new,d_y_new)
    print('Mean_x|Median_x :', np.nanmean(d_x_new), np.median(d_x_new) , 'Mean_y|Median_y :', np.nanmean(d_y_new), np.median(d_y_new))
    
##############################################################################
# Distribution Fitting #
else:    
    data_green = data_0I4
    data_green['ellipticity'] = data_green['sigmax_cI']/data_green['sigmay_cI']
    for i in filter_0I4:
        data_green = (data_green.where(eval(i))).dropna()
        
    data_yellow = data_II4 
    data_yellow['ellipticity'] = data_yellow['sigmax_cII']/data_yellow['sigmay_cII']
    for i in filter_II4:
        data_yellow= (data_yellow.where(eval(i))).dropna()
        
    x_green  =  (np.asarray(data_green.peak_x_cI))
    y_green  =  (np.asarray(data_green.peak_y_cI))
    x_yellow  = (np.asarray(data_yellow.peak_x_cII)+ shift_x)
    y_yellow  = (np.asarray(data_yellow.peak_y_cII)+ shift_y)
    
    c_orig = calculateDistance(x_green, x_yellow, y_green, y_yellow)[1][4]
    c = np.sort(c_orig)
    c_merge = calculateDistance(x_green, x_yellow, y_green, y_yellow)[1]
    d =calculateDistance((np.asarray(data_0I4.peak_x_cI)), (np.asarray(data_II4.peak_x_cII)+ shift_x), (np.asarray(data_0I4.peak_y_cI)), (np.asarray(data_II4.peak_y_cII)+ shift_y))[1][4]
    
    dist = calculateDistance(x_green, x_yellow, y_green, y_yellow)[1][4]
    r =[]
    for i in range(len(c_merge[1])):
        for ii in data_green.peak_y_cI:
            if c_merge[1][i] == ii:
                r.append(c_orig[i])
    data_green['r'] = r
    
    sigx_g = list((data_green.sigmax_cI).values)/np.sqrt(list((data_green.Ncounts_cI).values))
    sigy_g = list((data_green.sigmay_cI).values)/np.sqrt(list((data_green.Ncounts_cI).values))
    sigx_y = list((data_yellow.sigmax_cII).values)/np.sqrt(list((data_yellow.Ncounts_cII).values))
    sigy_y = list((data_yellow.sigmax_cII).values)/np.sqrt(list((data_yellow.Ncounts_cII).values))
      

    
    
    ###########
    
    

    a = np.histogram(c, bins = Binning, range=Fitrange) 
    A = []
    y = []
    
    A.append(0)
    y.append(0)
    for i in a[1]:
        corr = (a[1][1]-a[1][0])
        A.append(i+(corr/2))
        
    for i in a[0]:
        y.append(i)
    x = A[:-1]
    #y = a[0]

    if weighting == 'Poisson':
        W_orig = 1/np.array(y)
        W = 1/np.array(y)
        W[W == inf] = 1
        W[W == -inf] = 1
        W[W == np.nan] = 1
             
    else:
        W =1
        
    model = Model(Fitfunction)
    params = model.make_params()
    P = model.param_names
    for i in range(len(P)):
        params[P[i]].value= Init_value[i]
        params[P[i]].min= Bound_min[i]
        params[P[i]].max= Bound_max[i]
        params[P[i]].vary= Value_fit[i]
        params[P[i]].expr= Value_exp[i]
        
    fig1 = pp.figure(figsize = (7,7))
    pp.subplot(221)
    pp.hist2d(data_green.sigmax_cI, data_green.sigmay_cI, bins=100, cmap='hot')
    pp.xlabel('sigmax_cI')
    pp.ylabel('sigmay_cI')
    
    pp.subplot(222)
    pp.hist2d(data_yellow.sigmax_cII, data_yellow.sigmay_cII, bins=100, cmap='hot')
    pp.xlabel('sigmax_cII')
    pp.ylabel('sigmaY_cII')
    
#    pp.subplot(223)
#    pp.hist2d(data_green.chi_2_cI, data_green.intensity_cI, bins=100, cmap='hot')
#    pp.xlabel('chi_2_cI')
#    pp.ylabel('intensity_cI')
#    
#    pp.subplot(224)
#    pp.hist2d(data_yellow.chi_2_cII, data_yellow.intensity_cII, bins=100, cmap='hot')
#    pp.xlabel('chi_2_cII')
#    pp.ylabel('intensity_cII')  
    
    pp.subplot(223)
    pp.hist2d(data_green.sigmax_cI, data_green.ellipticity, range = [[min(data_green.sigmax_cI), max(data_green.sigmax_cI)], [min(data_green.ellipticity), 1.4]],bins=(100,100), cmap='hot')
    pp.xlabel('sigmax_cI')
    pp.ylabel('Ellipticity_green') 
    
    pp.subplot(224)
    pp.hist2d(data_yellow.sigmax_cII, data_yellow.ellipticity, range = [[min(data_yellow.sigmax_cII), max(data_yellow.sigmax_cII)], [min(data_yellow.ellipticity), 1.4]], bins=(100,100), cmap='hot')
    pp.xlabel('sigmax_cII')
    pp.ylabel('Ellipticity_yellow') 
    
    
    
    pp.savefig('{}{}'.format(filepath,'2D_distribution_chi2.png'))

    if Parameter_scan == False:
        fig4 = pp.figure()
        
        if Simuluation == True:
            x = np.linspace(Fitrange[0],Fitrange[1],10* Fitrange[1])
            pmr = np.array(params)
            y = Fitfunction(x, *pmr)
            y = y + np.random.randn(x.size)
            result = model.fit(y, params, x=x, weights = W, method = Method )
            pp.plot(x, y, 'blue', label='Simulation')
            pp.plot(x, result.best_fit, 'r-', label='best fit')
        
        else: 
            result = model.fit(y, params, x=x, weights = W, method = Method)

            pp.hist(c, bins=Binning, color='blue',alpha= 0.5, label='Cleaned_data', range=Fitrange)
            pp.hist(d, bins=Binning, color='blue',alpha= 0.1, label='Original_data', range=Fitrange)       
            pp.plot(x, result.best_fit, 'r-', label='best fit')
            pp.grid(b = True, linestyle=':')
            pp.xticks(np.arange(Fitrange[0], Fitrange[1], Fitrange[1]/10))
            
            if Gaussian_indiv_plot == True:
                pmr_new = []
                for i in result.best_values.values():
                    pmr_new.append(i)
                pmr_new = np.array(pmr_new)
                null = 0.0000
                y1 = Fitfunction(np.array(x), *(pmr_new*[1,1,1,1,1,null,1]))
                y2 = Fitfunction(np.array(x), *(pmr_new*[1,1,null,1,1,1,1]))
                pp.plot(x, y1, 'gold',       linestyle = '--', label='Gaussian_01')
                pp.plot(x, y2, 'darkorange', linestyle = '--', label='Gaussian_02')

            
        
        if Plot_init == True:
            pp.plot(x, result.init_fit, 'k--', label='initial fit')
        
        
        if Calc_uncertainty == True:
            dely = result.eval_uncertainty(sigma=2)  # 
            pp.fill_between(x, result.best_fit-dely, result.best_fit+dely, color="#ABABAB",label='2-$\sigma$ uncertainty band')
        
        pp.xlim(Fitrange)
        pp.xlabel('Distance RDA [nm]')
        pp.ylabel('counts')
        legend = pp.legend(loc='best', title= '{}{}'.format('Chi²_red: ', np.around(result.redchi,2) ))
        legend.get_title().set_fontsize(13)
        legend._legend_box.align = "left"
        pp.title(Fitfunction.__doc__, loc= 'left')
        pp.savefig('{}{}'.format(filepath,'1D_distribution.png')) 
        pp.show()
        
#        fig7 = pp.figure()
#        result.plot()

        #print(result.fit_report())        
        print('Chi²_red:', np.around( result.redchi, decimals = 2))
        print('#molecules_green:',x_green.size)
        print('#molecules_yellow:',x_yellow.size)
        print('Fit result:',(result.best_values))

        sigma_exp_g = np.sqrt(np.nanmean(sigx_g)**2+np.nanmean(sigy_g)**2)*pixelsize
        sigma_exp_y = np.sqrt(np.nanmean(sigx_y)**2+np.nanmean(sigy_y)**2)*pixelsize
        sigma_all_x = np.sqrt(np.nanmean(sigx_g)**2+np.nanmean(sigy_g)**2+np.std(sigx_g)**2+np.std(sigy_g)**2)*pixelsize
        sigma_all_y = np.sqrt(np.nanmean(sigx_y)**2+np.nanmean(sigy_y)**2+np.std(sigx_y)**2+np.std(sigy_y)**2)*pixelsize
        sigma_all   = np.sqrt(((np.nanmean(sigx_g)**2+np.nanmean(sigy_g)**2)/2 + (np.nanmean(sigx_y)**2+np.nanmean(sigy_y)**2)/2) +np.std(sigx_g)**2+np.std(sigy_g)**2 +np.std(sigx_y)**2+np.std(sigy_y)**2)*pixelsize
        
        #sigma_test= np.mean((np.sqrt((np.sqrt(np.sqrt(sigx_g**2 + sigy_g**2) * np.sqrt(sigx_y**2 + sigy_y**2)))**2))) + np.sqrt((np.sqrt(np.std(sigx_g)**2+ np.std(sigy_g)**2))**2 +(np.sqrt(np.std(sigx_y)**2+ np.std(sigy_y)**2))**2)
        print('sigma_g_expected:', np.around(sigma_exp_g,decimals = 2), 
              'sigma_y_expected:', np.around(sigma_exp_y,decimals = 2),
              'sigma_all_x:', np.around(sigma_all_x,decimals = 2),
              'sigma_all_y:', np.around(sigma_all_y,decimals = 2),
              'sigma_all:', np.around(sigma_all,decimals = 2))
        
        if Gaussian_indiv_plot == True:
            print('frac1: ', np.around(sum(y1)/(sum(y1)+sum(y2)),decimals = 2) )
            print('frac2: ', np.around(sum(y2)/(sum(y1)+sum(y2)),decimals = 2))
        
        if Probability_scan == True:
            emcee_kws = dict(steps=1000, burn=300, thin=20, is_weighted=False,progress=False)
            emcee_params = result.params.copy()
            #emcee_params.add('__lnsigma', value=np.log(0.1), min=np.log(0.001), max=np.log(2.0))
            result_emcee = model.fit(data=y, x=x, params=emcee_params, method='emcee',nan_policy='omit')#, fit_kws=emcee_kws)
            #emcee_corner = corner.corner(result_emcee.flatchain, labels= result_emcee.var_names, truths=list(result_emcee.params.valuesdict().values()))
            #lmfit.report_fit(result_emcee)
            emcee_plot = corner.corner(result_emcee.flatchain, labels= result_emcee.var_names, truths=list(result.params.valuesdict().values()))
            #pp.savefig('{}{}'.format(filepath,'2D_parameter_scan.png')) 
        
    else:
        Chi_red= []
        value_s = []
        Z = np.zeros((len(scan_value1), (len(scan_value2))))
        X = scan_value1
        Y = scan_value2
        for ii in range(len(scan_value1)):
            params[Parameter_scan1].value = scan_value1[ii]
            params[Parameter_scan1].min = -np.inf
            params[Parameter_scan1].max = np.inf
            params[Parameter_scan1].vary= False
            for i in range(len(scan_value2)):
                params[Parameter_scan2].value = scan_value2[i]
                params[Parameter_scan2].min = -np.inf
                params[Parameter_scan2].max = np.inf
                params[Parameter_scan2].vary= False
                result = model.fit(y, params, x=x, weights = W, method = Method)
                Chi_red.append(result.redchi)
                Z[i,ii]= eval('{}{}'.format('result.',Value_displayed))#result.redchi
                value_s.append(i)
                print('Chi²_red:',result.redchi)
        Z_corr = np.where(Z < Value_displayed_max, Z, Value_displayed_max)
        
        fig4 = pp.figure(figsize = (6,5))
        bounds = np.linspace(Z_corr.min(), Value_displayed_max, 1000)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        pp.contourf(X, Y, Z_corr, 500, norm= norm , cmap= 'hot_r')
        pp.colorbar()
        pp.xlabel(Parameter_scan1)
        pp.ylabel(Parameter_scan2)
        pp.savefig('{}{}'.format(filepath,'2D_parameter_scan.png')) 
        

if Plot_Overview_statistics == True:
    filter = ["data_green['r'] > 0","data_green['r'] < 500"]
    for i in filter:
        data_green = (data_green.where(eval(i))).dropna()
    
    xdata = (np.asarray(data_green.x_corr))
    xdata =  np.asarray([float(i) for i in xdata])
    ydata = (np.asarray(data_green.y_corr))
    ydata =  np.asarray([float(i) for i in ydata])
    zdata = (np.asarray(data_green.r))
    zdata =  np.asarray([float(i) for i in zdata])
    
    
    
    
    coord = np.asarray(data_green.y_offset)
    coord_list = list(set(coord))
    coord_list.sort()
    co = coord_list
    co =  np.around([float(i) for i in coord_list], decimals= 6)
    co2 =  np.around([float(i) for i in coord], decimals= 6)
    
    RDA= []
    Npoints= []
    RDA_mean = []
    for ii in coord_list:
        y = data_green.loc[data_green['y_offset'] == ii]
        RDA.append(y.r)
        RDA_mean.append(np.mean(y.r))
        Npoints.append(len(y.r))
    dt = pd.DataFrame(RDA)
    df = pd.DataFrame(np.transpose(dt.values), columns = co)
    
    fig1 = pp.figure(figsize = (8,6.5))
    pp.scatter(xdata, ydata, c=zdata, cmap='hot_r')
    cbar = pp.colorbar()
    pp.clim(0, 25)
    pp.xlabel('x-coord.')
    pp.ylabel('y-coord.')
    pp.title('RDA-distance map [all Overview] ', fontsize=16)
    cbar.set_label('RDA-distance [nm]', rotation=90)
    pp.savefig('{}{}'.format(filepath,'RDA-distance map [all Overview].png'))
    
    fig2, axs = pp.subplots(3, sharex=True, sharey=False, gridspec_kw={'hspace': 0}, figsize = (20,8))
    fig2.suptitle('Mean RDA-distance vs. Overview y_offset', fontsize=16)
    axs[0].plot(coord_list, Npoints, 'bo')
    axs[0].set_ylabel('# spots')
    axs[1].plot(coord_list, RDA_mean, 'ro')
    axs[1].set_ylabel('RDA_mean [nm]')
    sns.violinplot(data=df)
    pp.xticks(rotation='vertical')
    pp.xlabel('Overview y_offset [m]')
    pp.ylabel('RDA_distance [nm]') 
    pp.savefig('{}{}'.format(filepath,'Mean RDA-distance vs. Overview y_offset.png'), bbox_inches='tight') 
    
    fig3 = pp.figure(figsize = (8,6))
    ax = pp.axes(projection='3d')
    im = ax.scatter3D(xdata, ydata, co2, c=zdata, cmap='hot_r')
    #pp.xlabel('x-coord.')
    ax.set_zlabel('test')
    pp.ylabel('y-coord.')
    ax.view_init(45, 45)
    fig3.colorbar(im)
    pp.savefig('{}{}'.format(filepath,'3D-scatter vs. Overview y_offset.png'), bbox_inches='tight')
    

if Add_SgSrSyStoichiometry == True and Fit_only== False :
    location =os.path.dirname( os.path.dirname(os.path.dirname(filepath)))
    #path = '{}{}'.format(location, '/')
    path = 'N:/Singlem/singlem19-4/December/20_dB_noFRET_sb_hFRET_JHB\merged/export_JHB/'
    
    starting1 = 'Overview'   
    files = os.listdir(path)
     
    file_folder = []
    for i in files:
        if i.startswith(starting1):
            file_folder.append(i)
    
    COLOUR = ['Green', 'Red', 'Yellow']
    dic={'max_photons_Green':[],'max_photons_Red':[],'max_photons_Yellow':[]}
    name = []
    
    def getintensity(colours):
        global c
        
        c = []
        sum = 0       
        for ii in file_folder:
            sum += (len(ii)/len(ii))
            print((np.around(sum/len(file_folder)*100,decimals=1), '%'))
            folderpath = '{}{}{}'.format(path, ii, '/')
            Path_data = '{}{}{}{}'.format(folderpath, colours,' Photons', '/')
            data = '{}{}'.format(Path_data,os.listdir(Path_data)[0])
            name.append(ii)
            
            if os.stat(data).st_size == 0 :            
                c.append(0)
            else:
                c.append(np.sum(np.loadtxt(data, delimiter='\t')))

     
        
    
    for colour in COLOUR:
        getintensity(colour)
        dic['{}{}'.format('max_photons_', colour)] = c
        max_photons_green = dic['max_photons_Green']
        max_photons_yellow = dic['max_photons_Yellow']
        max_photons_red = dic['max_photons_Red']
    
    name = name[:int(len(name)/len(COLOUR))]
        
    pp.figure(figsize=(5,5))
    pp.subplot(331)
    pp.plot(max_photons_green)
    pp.plot(max_photons_yellow)
    pp.plot(max_photons_red)
    pp.subplot(332)   
    n =  np.true_divide(max_photons_green, max_photons_red)
    pp.plot(n)
    
    MEAN = np.mean(max_photons_green)
    STD = np.std(max_photons_green)
    Stoich = list(np.array(max_photons_green) + np.array(max_photons_red))/(np.array(max_photons_green) + np.array(max_photons_red) + np.array(max_photons_yellow))
    pp.subplot(333)
    boxdata = [max_photons_green,max_photons_yellow]
    #pp.boxplot(boxdata)
    coord = []
    for i in name:
        coord.append(i[14:24])
    
    coord_list = list(set(coord))
    coord_list.sort()
    
    ###
    
    STORE = pd.DataFrame(name, columns= ['filename'])
    STORE['Sg'] =  (max_photons_green)
    STORE['Sr'] =  (max_photons_red)
    STORE['Sy'] =  (max_photons_yellow)
    STORE['Stoichiometry'] =  Stoich
    STORE['Fr'] =  STORE['Sr']-0.5*STORE['Sg']
        
    box_add = []
    for ii in coord_list:
        box = []
        indices = [i for i, x in enumerate(coord) if x == ii]
        for iii in indices:
            box.append(max_photons_yellow[iii])
        box_add.append(box)
        
    pp.figure(figsize=(5,5))
    #pp.boxplot(box_add)
    
    dt = pd.DataFrame(box_add)
    df = pd.DataFrame(np.transpose(dt.values), columns = coord_list)
    #sns.boxplot(data=df)
    
    pp.figure(figsize=(5,5))
    sns.violinplot(data=df)
    pp.xticks(rotation='vertical')
    pp.xlabel('y-position')
    pp.ylabel('#photons-green') 
    
    
    fig =pp.figure(figsize= (5,5))
    x_ = STORE['Sg']
    y_ = STORE['Fr']
    pp.hist2d(x_, y_, bins=500, cmap=pp.cm.hot)
    pp.xlim(min(x_), 500)
    pp.ylim(min(y_), 500)
    pp.xlabel(x_.name)
    pp.ylabel(y_.name)
    
    
    Sg=[]    
    Sr=[] 
    Sy=[] 
    Stoichiometry=[] 
    for ii in range(len(data_0I4['data_name'])):
        for i in range(len(STORE['filename'])):
            if STORE['filename'][i] in data_0I4['data_name'][ii]:  
                Sg.append(STORE['Sg'][i])
                Sr.append(STORE['Sr'][i])
                Sy.append(STORE['Sy'][i])
                Stoichiometry.append(STORE['Stoichiometry'][i])
    data_0I4['Sg'] = Sg
    data_0I4['Sr'] = Sr
    data_0I4['Sy'] = Sy
    data_0I4['Stoich'] = Stoichiometry

            
            
    

        














