#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     COSR_AOD_extinction
# 
# Purpose:  
# 
#     Analysis of COSR data. Focus on case studies, and define the changes in AOD. Obtain the extinction coefficient from case     studies
#   
# Input:
# 
#     None
#   
# Output:
# 
#     Figures
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - load_utils.py : for loading OMI HDF5 files
#     - matplotlib
#     - numpy
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - COSR starsun.mat, and flag files
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-09-25
#     Modified: 

# # Prepare python environment

# In[1]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
import os
matplotlib.rc_file(os.path.join(os.getcwd(),'file.rc'))
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
import load_utils as lu
import plotting_utils as pu
from path_utils import getpath
import hdf5storage as hs
from datetime import datetime
from scipy.interpolate import UnivariateSpline
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
import scipy.stats as st


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[6]:


fp =getpath('COSR')


# # Load files

# ## Load the starsun

# In[ ]:





# ## Load the flag files

# In[ ]:





# # Run analysis and prepare variables
# Do some of the calculations to the data here

# # Plotting
# Present some fo the early plots here

# In[ ]:




