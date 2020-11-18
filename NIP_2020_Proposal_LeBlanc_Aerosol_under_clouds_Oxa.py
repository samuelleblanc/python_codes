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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2020-09-17
#     Modified:
# 

# # Prepare python environment

# In[1]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.colors as colors
get_ipython().magic(u'matplotlib notebook')
import os


# In[2]:


name = 'KORUS'
vv = 'v1'
fp = getpath(name)


# # Load files

# ## Load DIAL

# In[3]:


da, da_dict = lu.load_hdf(fp+'data_other/korusaq-DIAL_DC8_20160519_R1.h5',all_values=True)


# In[19]:


np.nanmax(da['//Data_Products/ext_532nm_prfl'])


# # Plot out data

# In[4]:


plt.figure()
plt.pcolor(da['//Nav_Data/Midtime'],da['//Nav_Data/Altitudes'],da['//Data_Products/ext_532nm_prfl'],
           vmin=0,vmax=0.1)
plt.colorbar(extend='both',label='Aerosol Extinction [Mm$^{-1}$]')


# In[5]:


np.nanmin(da['//Data_Products/bsc_532nm_prfl']),np.nanmean(da['//Data_Products/bsc_532nm_prfl'])


# In[6]:


norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max())


# In[10]:


plt.figure()
plt.pcolor(da['//Nav_Data/Midtime']/60.0/60.0,da['//Nav_Data/Altitudes'],da['//Data_Products/bsc_532nm_prfl'],
            norm=colors.LogNorm(vmin=0.00001, vmax=0.01))
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')
plt.colorbar(extend='both',label='Aerosol Backscatter 532 nm [Mm$^{-1}$sr$^{-1}$]')
plt.scatter(da['//Nav_Data/Midtime'][cld_mask==1]/60.0/60.0,da['//Nav_Data/Altitudes'][cld_mask==1],'.r')
#plt.pcolor(da['//Nav_Data/Midtime']/60.0/60.0,da['//Nav_Data/Altitudes'],cld_mask,cmap=plt.cm.autumn)


# In[11]:


cld_mask.shape


# In[12]:


da['//Nav_Data/Midtime'].shape


# In[13]:


da['//Nav_Data/Altitudes'].shape


# In[18]:


nx,ny = cld_mask.shape
alts, time = [],[]
for i in xrange(nx):
    for j in xrange (ny):
        if cld_mask[i,j] ==1.0:
            alts.append(da['//Nav_Data/Altitudes'][i])
            time.append(da['//Nav_Data/Midtime'][0,j])
time = np.array(time)
alts = np.array(alts)


# In[45]:


nx,ny = cld_mask.shape
alts_depol, time_depol = [],[]
for i in xrange(nx):
    for j in xrange (ny):
        if (da['//Data_Products/depol_532nm_prfl'][i,j] > 0.1) &           (da['//Data_Products/bsc_532nm_prfl'][i,j]>0.0001) &           (da['//Nav_Data/Altitudes'][i] > 7500.0):
            alts_depol.append(da['//Nav_Data/Altitudes'][i])
            time_depol.append(da['//Nav_Data/Midtime'][0,j])
time_depol = np.array(time_depol)
alts_depol = np.array(alts_depol)


# In[9]:


cld_mask = da['//Data_Products/Cloud_Mask_prfl']
cld_mask[cld_mask==0] = np.nan


# In[25]:


plt.figure()
plt.pcolor(da['//Nav_Data/Midtime']/60.0/60.0,da['//Nav_Data/Altitudes'],da['//Data_Products/bsc_532nm_prfl'],
            norm=colors.LogNorm(vmin=0.00001, vmax=0.01))
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')
plt.colorbar(extend='both',label='Aerosol Backscatter 532 nm [Mm$^{-1}$sr$^{-1}$]')
plt.plot(time/3600.0,alts,'.r')
plt.xlim(23.1,24.7)


# In[68]:


plt.figure(figsize=(6,1.8))
plt.pcolor(da['//Nav_Data/Midtime']/60.0/60.0,da['//Nav_Data/Altitudes']/1000.0,da['//Data_Products/bsc_532nm_prfl'],
            norm=colors.LogNorm(vmin=0.00001, vmax=0.01))
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [km]')
plt.colorbar(extend='both',label='Aerosol Backscatter\n 532 nm [Mm$^{-1}$sr$^{-1}$]')
plt.plot(time_depol/3600.0,alts_depol/1000.0,'.r',label='cloud')
plt.legend(frameon=True)
plt.xlim(23.3,24.6)
plt.ylim(0,17.500)
plt.subplots_adjust(bottom=0.25)
plt.savefig(fp+'plot/KORUS_DIAL_20180519_aer_bcsc_cloud.png',dpi=600,transparent=True)


# In[57]:


plt.figure()
plt.contourf(da['//Nav_Data/Midtime'].flatten()/60.0/60.0,da['//Nav_Data/Altitudes'].flatten(),da['//Data_Products/bsc_532nm_prfl'],
            norm=colors.LogNorm(vmin=0.00001, vmax=0.01),levels=15)
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')
plt.colorbar(extend='both',label='Aerosol Backscatter 532 nm [Mm$^{-1}$sr$^{-1}$]')
plt.plot(time_depol/3600.0,alts_depol,'.r')
plt.xlim(23.3,25.3)


# In[50]:


plt.figure()
plt.pcolor(da['//Nav_Data/Midtime'].flatten()/60.0/60.0,da['//Nav_Data/Altitudes'].flatten(),da['//Data_Products/depol_532nm_prfl'],
          vmin=0.00001, vmax=0.5)
plt.colorbar()
plt.xlim(23.3,25.3)


# In[58]:


plt.figure()
plt.tricontour(da['//Nav_Data/Midtime'].flatten()/60.0/60.0,da['//Nav_Data/Altitudes'].flatten(),da['//Data_Products/depol_532nm_prfl'],
          vmin=0.00001, vmax=0.5,levels=4, linewidths=0.5, colors='k')
#plt.colorbar()
plt.xlim(23.3,25.3)


# In[33]:


plt.figure()
plt.pcolor(da['//Nav_Data/Midtime']/60.0/60.0,da['//Nav_Data/Altitudes'],da['//Data_Products/bsc_532nm_prfl'],
            norm=colors.LogNorm(vmin=0.00001, vmax=0.01))
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')
plt.colorbar(extend='both',label='Aerosol Backscatter 532 nm [Mm$^{-1}$sr$^{-1}$]')
plt.contour(da['//Nav_Data/Midtime'].flatten()/60.0/60.0,da['//Nav_Data/Altitudes'].flatten(),da['//Data_Products/depol_532nm_prfl'])
plt.xlim(23.1,24.7)


# In[34]:


plt.figure()
plt.contour(da['//Nav_Data/Midtime'].flatten()/60.0/60.0,da['//Nav_Data/Altitudes'].flatten(),da['//Data_Products/depol_532nm_prfl'])


# In[38]:


plt.figure()
plt.pcolor(da['//Nav_Data/Midtime'].flatten()/60.0/60.0,da['//Nav_Data/Altitudes'].flatten(),da['//Data_Products/depol_532nm_prfl'],
          vmin=0.00001, vmax=0.5)
plt.colorbar()


# In[7]:


plt.figure()
plt.pcolor(da['//Nav_Data/Midtime'],da['//Nav_Data/Altitudes'],da['//Data_Products/Cloud_Mask_prfl'])
plt.colorbar(extend='both')


# In[ ]:




