#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Describe the details ...
# 
# Input:
# 
#     arguments
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
# 
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-08-27
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


# In[5]:


name = 'ARCTAS'
vv = 'v1'
fp = getpath(name)


# # Load files

# In[11]:


dat,dat_d = lu.load_netcdf('/data/sunsat/ARCTAS_2008/data_archival/ARCTAS_AATS14_2008_R2.nc',everything=True)


# In[13]:


dat.keys()


# In[14]:


dat[b'wavelength']


# In[15]:


dat[b'AOD'].shape


# In[17]:


dat[b'AOD'][dat[b'AOD']<=-10] = np.nan


# In[26]:


dat_d[b'time']


# In[24]:


dat[b'time']


# # Plot out data

# In[32]:


plt.figure(figsize=(6,2))

plt.hist(dat[b'AOD'][:,3],bins=50,range=[0,1])
plt.xlabel('AOD$_{{500}}$')
plt.title('Measured AOD distribution from ARCTAS')
plt.tight_layout()


# In[ ]:




