
# coding: utf-8

# # Intro
# Name:  
# 
#     MAP_proposal_profiles
# 
# Purpose:  
# 
#     Python script for plotting AOD profile and extinctions for various places
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
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - os
#     - datetime
#     - mpl_toolkits
#     - plotting_utils (user defined plotting routines)
#     - map_utils, dependent on geopy
#     - Basemap
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, NASA Ames, 2016-06-08

# # Import the required modules and do the setup

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
get_ipython().magic(u'matplotlib notebook')
from mpltools import color
import numpy as np
import scipy.io as sio
import hdf5storage as hs
import math
import os
import Sp_parameters as Sp


# In[2]:

from load_utils import mat2py_time, toutc


# In[3]:

# set the basic directory path
fp='C:\\Users\\sleblan2\\Research\\SEAC4RS\\'


# # Load some required files

# In[4]:

star = sio.loadmat(fp+'dc8\\20130816\\20130816starsun_R2.mat',variable_names=('w','tau_aero','t','Alt','Lat','Lon'))


# In[5]:

star['tt'] = mat2py_time(star['t'])
star['utc'] = toutc(star['tt'])


# In[6]:

star_cl = sio.loadmat(fp+'starsun\\20130913starsun_R2_tauaero.mat',variable_names=('w','tau_aero','t','Alt','Lat','Lon'))


# In[ ]:



