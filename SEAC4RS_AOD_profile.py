
# coding: utf-8

# Name:  
# 
#     SEAC4RS_AOD_profile
# 
# Purpose:  
# 
#     Python script for plotting boundary layer AOD profile
# 
# Calling Sequence:
# 
#     python SEAC4RS_AOD_profile
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
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, NASA Ames, 2015-10-09

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
import math
import os
import Sp_parameters as Sp


# In[6]:

from load_modis import mat2py_time, toutc


# In[2]:

# set the basic directory path
fp='C:\\Users\\sleblan2\\Research\\SEAC4RS\\'


# # Load the 4STAR starsun file

# In[3]:

star = sio.loadmat(fp+'dc8\\20130816\\20130816starsun_R2.mat',variable_names=('w','tau_aero','t','Alt','Lat','Lon'))


# In[4]:

star


# In[5]:

star.keys()


# In[7]:

star['tt'] = mat2py_time(star['t'])
star['utc'] = toutc(star['tt'])


# In[8]:

plt.plot(star['utc'],star['Alt'])


# In[11]:

plt.plot(star['t'],star['Alt'])


# In[13]:

plt.figure()
plt.plot(star['t'],star['Alt'])


# In[ ]:



