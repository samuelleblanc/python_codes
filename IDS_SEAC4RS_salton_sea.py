
# coding: utf-8

# # Info
# Name:  
# 
#     ## insert name here
# 
# Purpose:  
# 
#     ## add description
#   
# Input:
# 
#     ## inputs
#   
# Output:
# 
#     ##variables, figures and save files...
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
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-05-18
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


# In[4]:


fp =getpath('SEAC4RS')


# # Load files

# In[5]:


s = lu.load_ict(fp+'data/SEAC4RS-4STAR-AOD-CWV_DC8_20130923_R2.ict')


# In[6]:


plt.figure()
plt.plot(s['Start_UTC'],s['GPS_alt'],'.')


# In[7]:


plt.figure()
plt.plot(s['Start_UTC'],s['AOD0501'],'.')


# In[9]:


plt.figure()
plt.plot(s['Longitude'],s['Latitude'],'.')


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# # Plotting
# Present some fo the early plots here
