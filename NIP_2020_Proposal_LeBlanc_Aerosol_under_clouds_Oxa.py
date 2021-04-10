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


# # Load rtm model data

# In[3]:


import load_utils as lu


# In[5]:


rtm,rtmd = lu.load_netcdf('/home/sam/unl-vrtm/unl-vrtm-2.1/run/O3band.unlvrtm.nc',everything=True)


# In[62]:


rtmd['LinPar']


# In[59]:


rtmd['tauAER_WFS']


# In[65]:


rtm['LinPar']


# In[30]:


rtm['Flux_Direct']


# In[55]:


plt.figure()
plt.plot(rtm['Lamdas'],rtm['omegaAER'].T[:,0,0],label='Aerosol type 1')
plt.ylabel('SSA of the aerosol')
plt.xlabel('Wavelength [nm]')
plt.title('SSA for aerosol')


plt.plot(rtm['Lamdas'],rtm['omegaAER'].T[:,0,1],label='Aerosol type 2')
plt.legend()


# In[43]:


plt.figure()
plt.plot(rtm['Lamdas'],rtm['omegaIOP'].T)
plt.ylabel('SSA of the atmosphere')
plt.xlabel('Wavelength [nm]')
plt.title('SSA at various atmospheric layers')


# In[37]:


plt.figure()
plt.plot(rtm['Lamdas'],rtm['Flux_Direct'][0,0,0,:])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Direct Irradiance [W/m$^2$/nm]')
plt.title('Ox-A modeled direct beam irradiance UNL-VRTM')


# ## For 2 different ssa

# In[70]:


rtm85_2,rtmd85_2 = lu.load_netcdf('/home/sam/unl-vrtm/unl-vrtm-2.1/run/OxA_ssa085.unlvrtm.nc',everything=True)


# In[67]:


rtm97,rtmd97 = lu.load_netcdf('/home/sam/unl-vrtm/unl-vrtm-2.1/run/OxA_ssa097.unlvrtm.nc',everything=True)


# In[68]:


plt.figure()
plt.plot(rtm97['Lamdas'],rtm97['Flux_Direct'][0,0,0,:],label='SSA=0.97')
plt.plot(rtm85['Lamdas'],rtm85['Flux_Direct'][0,0,0,:],label='SSA=0.85')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Direct Irradiance [W/m$^2$/nm]')
plt.title('Ox-A modeled direct beam irradiance UNL-VRTM (AOD=0.4)')
plt.legend()


# In[71]:


plt.figure()
plt.plot(rtm97['Lamdas'],rtm97['Flux_Direct'][0,0,0,:]/rtm97['Flux_Direct'][0,0,0,0],label='SSA=0.97')
plt.plot(rtm85['Lamdas'],rtm85['Flux_Direct'][0,0,0,:]/rtm85['Flux_Direct'][0,0,0,0],label='SSA=0.85')
plt.plot(rtm85['Lamdas'],rtm85_2['Flux_Direct'][0,0,0,:]/rtm85_2['Flux_Direct'][0,0,0,0],label='SSA=0.85')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalized Direct Irradiance')
plt.title('Ox-A modeled direct beam irradiance UNL-VRTM (AOD=0.4)')
plt.legend()


# In[72]:


rtm85_2['Flux_Direct'].shape


# In[73]:


rtm85['Flux_Direct'].shape


# In[77]:


rtm85_d,rtmd85_d = lu.load_netcdf('/home/sam/unl-vrtm/unl-vrtm-2.1/run/OxA_ssa085_deg30.unlvrtm.nc',everything=True)
rtm97_d,rtmd97_d = lu.load_netcdf('/home/sam/unl-vrtm/unl-vrtm-2.1/run/OxA_ssa098_deg30.unlvrtm.nc',everything=True)


# In[78]:


plt.figure()
plt.plot(rtm97_d['Lamdas'],rtm97_d['Flux_Direct'][0,0,0,:]/rtm97_d['Flux_Direct'][0,0,0,0],label='SSA=0.97')
plt.plot(rtm85_d['Lamdas'],rtm85_d['Flux_Direct'][0,0,0,:]/rtm85_d['Flux_Direct'][0,0,0,0],label='SSA=0.85')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalized Direct Irradiance')
plt.title('Ox-A modeled direct beam irradiance UNL-VRTM (AOD=0.4)')
plt.legend()


# In[79]:


plt.figure()
plt.plot(rtm97_d['Lamdas'],rtm97_d['Flux_Direct'][0,0,0,:],label='SSA=0.97')
plt.plot(rtm85_d['Lamdas'],rtm85_d['Flux_Direct'][0,0,0,:],label='SSA=0.85')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalized Direct Irradiance')
plt.title('Ox-A modeled direct beam irradiance UNL-VRTM (AOD=0.4)')
plt.legend()


# In[80]:


import pyunlvrtm as pum


# In[102]:


a = pum.read_unlvrtm('/home/sam/unl-vrtm/unl-vrtm-2.1/run/OxA_ssa085_deg30.unlvrtm.nc')
b = pum.read_unlvrtm('/home/sam/unl-vrtm/unl-vrtm-2.1/run/OxA_ssa098_deg30.unlvrtm.nc')


# In[103]:


a['Stokes']


# In[105]:


plt.figure()
plt.plot(a['Lamdas'],a['Stokes'][0,:])
plt.plot(b['Lamdas'],b['Stokes'][0,:])


# In[104]:


plt.figure()
plt.plot(a['Lamdas'],a['Stokes'][0,:]/a['Stokes'][0,0])
plt.plot(b['Lamdas'],b['Stokes'][0,:]/b['Stokes'][0,0])


# In[ ]:




