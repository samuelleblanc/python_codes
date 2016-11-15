
# coding: utf-8

# # Info
# Name:  
# 
#     SEAC4RS_cld_compare_v5_small
# 
# Purpose:  
# 
#     Python script to simplify the SEAC4RS_compare_v5 script. 
#     Used for making the figures in the paper:
#         Comparing Cloud properties and radiative effect estimated from airborne measurements of transmitted and reflected light
#         LeBlanc et al., JGR 
# 
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures and save files...
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - run_kisq_retrieval.py : for the retrieval functions
#     - load_utils.py : for loading modis files
#     - matplotlib
#     - mpltools
#     - numpy
#     - plotly : optional
#     - scipy : for saving and reading
#     - math
#     - os
#     - gc
#     - pdb
#     - datetime
#     - pyhdf
#     - mpl_toolkits
#     - gdal (from osgeo)
#     - plotting_utils (user defined plotting routines)
#     - map_utils, dependent on geopy
#     - hdf5storage
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - sp_v3_20130913_4STAR.out: modeled spectra output for SEAC4RS at sza 17.9, idl save file
#   - %%20130219starzen_rad.mat : special zenith radiance 4star matlab file 
#   - ict files from 20130913
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2016-11-11
#              ported from SEAC4RS_compare_v5

# # Load the required python modules

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color

import numpy as np
import scipy.io as sio
import math
import os
import Sp_parameters as Sp

import hdf5storage as hs
from load_utils import mat2py_time, toutc, load_ict

from Sp_parameters import smooth
import cPickle as pickle


# In[2]:

get_ipython().magic(u'matplotlib notebook')


# In[3]:

# set the basic directory path
fp='C:/Users/sleblan2/Research/SEAC4RS/'


# # Load some files for easier processing

# In[4]:

m_dict = hs.loadmat(fp+'20130913_retrieval_output.mat')


# In[5]:

m_dict.keys()


# In[6]:

if not 'emas_tau_full' in vars():
    print 'not defined, loading from file'
    emas_tau_full = m_dict['emas'][1]; emas_ref_full = m_dict['emas'][3]; emas_utc_full = m_dict['emas'][5]
    modis_tau = m_dict['modis'][1]; modis_ref = m_dict['modis'][3]
    ssfr_tau = m_dict['ssfr'][1]; ssfr_ref = m_dict['ssfr'][3]; ssfr_utc = m_dict['ssfr'][5]
    rsp_tau = m_dict['rsp'][1]; rsp_ref = m_dict['rsp'][3]; rsp_utc = m_dict['rsp'][5]
    star_tau = m_dict['star'][1]; star_ref = m_dict['star'][3]
    goes_tau = m_dict['goes'][1]; goes_ref = m_dict['goes'][3]
    goes_utc = m_dict['goes'][5]; star_utc = m_dict['star'][5]


# # Now make a plot of the time series for easier referencing

# In[7]:

plt.figure(figsize=(9,7))
ax1 = plt.subplot(211)
ax1.plot(star_utc,smooth(modis_tau,6),'m->',label='MODIS',markeredgecolor='none')
ax1.plot(emas_utc_full,smooth(emas_tau_full,60),'ko',label='eMAS',markeredgecolor='none')
ax1.plot(ssfr_utc,smooth(ssfr_tau,2),'g-x',label='SSFR')
ax1.plot(rsp_utc,smooth(rsp_tau,70),'c-+',label='RSP')
ax1.plot(goes_utc,smooth(goes_tau,2),'b-*',label='GOES',markeredgecolor='none')
ax1.plot(star_utc,smooth(star_tau,40),'r-s',label='4STAR',markeredgecolor='none')

ax1.errorbar(star_utc,smooth(modis_tau,6),yerr=modis_tau_std*2.0,color='m')
ax1.errorbar(ssfr_utc,smooth(ssfr_tau,2),yerr=ssfr_tau_std*2.0,color='g')
ax1.errorbar(rsp_utc,smooth(rsp_tau,70),yerr=rsp_tau_std*2.0,color='c')
ax1.errorbar(goes_utc,smooth(goes_tau,2),yerr=goes_tau_std*2.0,color='b')
ax1.errorbar(star_utc,smooth(star_tau,40),yerr=star_tau_std*2.0,color='r')

ax1.legend(frameon=False,numpoints=1)
ax1.grid()
#ax1.set_xlabel('UTC [H]')
ax1.set_ylabel('$\\tau$')
ax1.set_ylim([0,100])

ax2 = plt.subplot(212,sharex=ax1)
ax2.plot(star_utc,smooth(modis_ref,6),'m->',label='MODIS',markeredgecolor='none')
ax2.plot(emas_utc_full,smooth(emas_ref_full,60),'k-o',label='eMAS',markeredgecolor='none')
ax2.plot(ssfr_utc,smooth(ssfr_ref,2),'g-x',label='SSFR')
ax2.plot(rsp_utc,smooth(rsp_ref,70),'c-+',label='RSP')
ax2.plot(goes_utc,smooth(goes_ref,2),'b-*',label='GOES',markeredgecolor='none')
ax2.plot(star_utc,smooth(star_ref,40),'r-s',label='4STAR',markeredgecolor='none')

ax2.errorbar(star_utc,smooth(modis_ref,6),yerr=modis_ref_std*2.0,color='m')
ax2.errorbar(ssfr_utc,smooth(ssfr_ref,2),yerr=ssfr_ref_std*2.0,color='g')
ax2.errorbar(rsp_utc,smooth(rsp_ref,70),yerr=rsp_ref_std*2.0,color='c')
ax2.errorbar(goes_utc,smooth(goes_ref,2),yerr=goes_ref_std*2.0,color='b')
ax2.errorbar(star_utc,smooth(star_ref,40),yerr=star_ref_std*2.0,color='r')

#ax2.legend(frameon=False,numpoints=1)
ax2.grid()
ax2.set_ylim([0,60])
ax2.set_xlabel('UTC [H]')
ax2.set_ylabel('r$_{eff}$ [$\\mu$m]')

plt.savefig(fp+'plots/20130911_retrieved_horz_var.png',dpi=600,transparent=True)


# In[ ]:



