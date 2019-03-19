
# coding: utf-8

# # Info
# Name:  
# 
#     ARM_SaSza_spetra_test
# 
# Purpose:  
# 
#     First go through of SaS-ze spectral data for testing application to cloud properties
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
#     - load_utils.py : for loading netcdf files
#     - matplotlib
#     - numpy
#     - scipy : for saving and reading
#     - os
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, NOAA-Boulder, 2019-02-26
#     Modified: 

# # Load the python modules and environment

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import load_utils as lu
from path_utils import getpath


# In[31]:


import os
import Sp_parameters as Sp


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[7]:


fp = getpath('LASIC',make_path=True)
fp


# In[8]:


daystr = '20160815'


# # Load the files

# In[10]:


fl = os.listdir(fp+'data/')


# In[22]:


for g in fl: 
    if daystr in g:
        if 'vis' in g:
            glv = g
        if 'nir' in g:
            gln = g


# In[23]:


glv,gln


# In[24]:


vis,vish = lu.load_netcdf(fp+'data/'+glv,everything=True)


# In[25]:


nir,nirh = lu.load_netcdf(fp+'data/'+gln,everything=True)


# In[29]:


vish['zenith_radiance'], vish['wavelength']


# In[30]:


nirh['zenith_radiance'], nirh['wavelength']


# # Plot out some data

# In[32]:


help(Sp.Sp)

