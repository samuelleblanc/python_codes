#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     To look into the autocorrelation distances for ORACLES SSFR cloud properties retrievals
#     - For support on the MODIS/TASNPP 
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-02-01
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


name = 'ORACLES'
vv = 'v4'
fp = getpath(name)


# # Load files

# In[6]:


s = hs.loadmat(fp+'ORACLES_DARE_{}.mat'.format(vv))


# In[7]:


s.keys()


# # Special functions

# In[ ]:





# # Plot out data

# In[ ]:




