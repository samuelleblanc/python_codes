
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


# In[ ]:



