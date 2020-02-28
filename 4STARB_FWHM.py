#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     4STARB_FWHM
# 
# Purpose:  
# 
#     Plot out the Full Width Half Max slit function variability of 4STARB's spectrometers
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
#   - load_utils.py
#   - matplotlib
#   - numpy
#   - path_utils
#   - write_utils
#   
# Needed Files:
#   - visFWHM.txt
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, NASA Ames Research Center., CA, 2019-10-04
#     Modified: 
# 

# # Prepare python environment

# In[7]:


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


# In[8]:


name = '4STARB'
vv = 'v1'
fp = getpath('4STAR_data')


# In[ ]:





# # Load files

# In[9]:


vis = np.genfromtxt(fp+'visFWHM.txt')


# In[10]:


vis


# # Plot out data

# In[17]:


plt.figure()
plt.plot(vis[:,0],vis[:,1],'.')
plt.xlabel('Wavelength [nm]')
plt.ylabel('FWHM [nm]')
plt.grid()
plt.title('4STAR FWHM from GSFC 2013-06-05')
plt.savefig(fp+'../../4STAR/4STAR_visFWHM_20130605.png',dpi=600,transparent=True)


# In[ ]:




