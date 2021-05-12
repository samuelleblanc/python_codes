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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-04-07
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


# In[8]:


name = 'cloud_retrieval'
vv = 'v1'
fp = getpath(name)


# # Load Model files for SEAC4RS

# In[9]:


lut = sio.idl.readsav(getpath('SEAC4RS')+'model/sp_v5_20130913_raw.out')


# In[12]:


lut.keys()


# In[19]:


sea_lut = Sp.Sp(lut,irrad=True)


# In[24]:


sea_lut['sp'].shape, sea_lut['sp_irrdn'].shape


# In[26]:


# switch out the radiance for the downwelling irradiance
sea_lut.sp = sea_lut['sp_irrdn']


# ## Calculate the params 

# In[27]:


sea_lut.params()


# In[28]:


sea_lut.param_hires()


# ## Plot out the lut

# In[29]:


lut['sza']


# In[76]:


sea_lut['par'][0,10,20,:],sea_lut['par'][1,10,20,:]


# In[30]:


figl,axl = Sp.plot_lut_vs_tau(sea_lut)


# # Load the NAAMES models

# In[41]:


lutb = hs.loadmat('/data/sam/NAAMES/lut/v2_NAAMES_lut.mat')


# In[42]:


lutb.keys()


# In[51]:


lutb['rad'].shape


# In[77]:


lutn = []
for s in xrange(len(lutb['sza'])):
    sptemp = {}
    sptemp['tau'] = lutb['tau']
    sptemp['ref'] = lutb['ref']
    sptemp['zout'] = lutb['zout']
    sptemp['sza'] = lutb['sza']
    sptemp['phase'] = lutb['phase']
    sptemp['irr_dn_diff'] = lutb['irr_dn_diff'][:,:,:,:,:,s]
    sptemp['irr_dn'] = lutb['irr_dn'][:,:,:,:,:,s]
    sptemp['irr_up'] = lutb['irr_up'][:,:,:,:,:,s]
    sptemp['wvl'] = lutb['wvl']
    sptemp['rad'] = lutb['irr_dn'][:,:,:,:,:,s] #switch out for irradiance
    ltemp = Sp.Sp(sptemp,verbose=False)

    ltemp.params(ice_only=False)
    ltemp.param_hires()
    lutn.append(ltemp)


# In[66]:


lutb['sza']


# ## plot the lut

# In[88]:


lutb['rad'][0,500,0,10,10,5],lutb['rad'][1,500,0,10,10,5]


# In[71]:


lutn[5].par[0,10,20,:]


# In[72]:


lutn[5].par[1,10,20,:]


# In[94]:


# for sza at 61deg
figl,axl = Sp.plot_lut_vs_tau(lutn[5],forceliq=False,forceice=False)
plt.suptitle('Irradiance derived modeled parameters for North Atlantic sza: {}'.format(lutb['sza'][5]))
plt.savefig(fp+'model/NAAMES_irr_lut_sza61.png',dpi=600,transparent=True)


# In[93]:


figl,axl = Sp.plot_lut_vs_tau(lutn[10],forceliq=False,forceice=False)
plt.suptitle('Irradiance derived modeled parameters for North Atlantic sza: {}'.format(lutb['sza'][10]))
plt.savefig(fp+'model/NAAMES_irr_lut_sza71.png',dpi=600,transparent=True)


# In[ ]:




