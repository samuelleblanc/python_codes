#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Work with Dan Lubin et al. for cloud remote sensing using AWARE/Antarctic and other cloud retrievals.
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-07-30
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


# In[11]:


import pandas as pd


# In[9]:


name = 'AWARE'
vv = 'v1'
fp = getpath(name)


# # Load files

# In[12]:


fl = os.listdir(fp+'data/')


# In[41]:


fl


# ## Run through and read the values

# In[70]:


d = []
for f in fl:
    d_tmp = pd.read_csv(fp+'data/'+f,delimiter='\t',header=25)
    v_tmp = d_tmp.keys()
    par_tmp = Sp.param(d_tmp[v_tmp[1]].to_numpy(), d_tmp['Wavelength'].to_numpy())
    print(v_tmp[1])
    for i,p in enumerate(par_tmp): print('   eta_{} = {:1.4f}'.format(i+1,p))
    d.append({'wvl':d_tmp['Wavelength'].to_numpy(),'rad':d_tmp[v_tmp[1]].to_numpy(),'name':v_tmp[1],'par':par_tmp})


# # Plot out data

# In[ ]:




