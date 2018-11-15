
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_aerosol_above_cloud_from_insitu
# 
# Purpose:  
# 
#     Get an average vertical profile of aerosol and cloud in situ properties from ORACLES 2016
#   
# Input:
# 
#     none
# 
# Output:
#    
#     plots
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - matplotlib
#     - scipy
# 
#   
# Needed Files:
# 
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc,Santa Cruz, CA, 2018-08-17
#     

# # Load python environment modules, and paths

# In[1]:


import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


import load_utils as lu
from Sp_parameters import smooth

import hdf5storage as hs
from mpltools import color
from path_utils import getpath
from write_utils import nearest_neighbor
from tqdm import tqdm_notebook as tqdm


# In[4]:


fp = getpath('ORACLES')
fp


# # Load the desired files
