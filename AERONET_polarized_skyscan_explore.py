#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     To go through and process and explore the AERONET polarized skyscans for the Ames' aeronet
# 
# Input:
# 
#     None
# 
# Output:
# 
#     Figure and save files
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - load_utils.py
#     - matplotlib
#     - numpy
#     - Sp_parameters
#     - write_utils
#     - path_utils
#     - hdf5storage
#     - scipy
#     - pandas
# 
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - output from the command: wget --no-check-certificate -q -O NASA_Ames_PPP_ALP_HYP_alltime.dat "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_raw_sky_v3?site=NASA_Ames&year=2013&month=6&day=1&year2=2021&month2=1&day2=14&HYP00=1&PPP00=1&ALP00=1&AVG=10"
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-02-18
#     Modified:
# 

# # Prepare python environment

# In[3]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import os
import pandas as pd


# In[8]:


name = 'rooftop'
vv = 'v1'
fp_dat = getpath(name)
fp = getpath('aeronet')
fp_pl = fp+'plot/'


# # Load files

# In[16]:


d = pd.read_csv(fp_dat+'NASA_Ames_PPP_ALP_HYP_alltime.dat',header=7)


# In[17]:


d.keys()


# # Plot out data

# In[ ]:




