
# coding: utf-8

# # Info
# Name:  
# 
#     ATDM Proposal dust heating rate
# 
# Purpose:  
# 
#     Make some figures for the ATMD proposal
#     Focus on the heating rate profiles obtained during ORACLES 2018 transit back, near Cabo Verde, with SSFR
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
#     - 
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - for_Sam_20181025.out
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-08-23
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
import scipy.io as sio


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


fp =getpath('ORACLES')


# # Load files

# ## Load the SSFR files

# In[10]:


ssfr = sio.idl.readsav(fp+'data_other_2018/SSFR/for_Sam_20181025_SSFR.out')


# In[11]:


ssfr.keys()


# ## Load the 4STAR file

# In[14]:


s = lu.load_ict(fp+'aod_ict_2018/4STAR-AOD_P3_20181025_R1.ict')


# In[17]:


sp = sio.loadmat(fp+'data_2018/4STAR_20181025starsun.mat')


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# Get the location of the dust profile

# In[21]:


plt.figure()
plt.plot(s['Start_UTC'],s['GPS_Alt'],'.')


# # Plotting
# Present some fo the early plots here
