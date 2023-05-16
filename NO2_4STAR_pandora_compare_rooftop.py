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
#   
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-08-26
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


# In[10]:


name = 'rooftop'
vv = 'v1'
fp = getpath(name)


# # Load files

# ## Load pandora

# In[20]:


i = 0
for a in [1,2,3,4,5,6,7,8,9]:
    if i>2: break
    if a>2: i = i+1
    print(a)


# In[22]:


l = 'Data file version: rout0p1-7'


# In[24]:


linesplit = l.strip().split(':')


# In[26]:


l0 = linesplit[0]


# In[27]:


l0.strip().replace(' ','_')


# In[37]:


filepath = fp+'Pandora34s1_MountainViewCA_L2Tot_rout0p1-7.txt'

header = {}
i = 0
n = 0 
f = open(filepath,'r')
for l in f:
    if i>1: break
    n = n+1
    if l.strip().startswith('--'): 
        i=i+1
        continue
    linesplit = l.strip().split(':')
    header[linesplit[0].strip().replace(' ','_')] = linesplit[1].strip()
f.close()


# In[39]:


n


# In[38]:


header


# In[16]:


p = pd.read_csv(filepath,delimiter=' ',header=n)


# In[17]:


p


# In[ ]:


g =


# # Plot out data

# In[ ]:




