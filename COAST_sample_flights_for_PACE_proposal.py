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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2023-04-27
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
get_ipython().run_line_magic('matplotlib', 'notebook')
import os


# In[8]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature


# In[5]:


name = 'sunsat'
vv = 'v1'
fp = getpath(name)


# In[6]:


fp = fp+'COAST_2011/'


# In[16]:


fpp = fp+'plots/'


# # Load files

# In[7]:


a,adict = lu.load_netcdf(fp+'data_archival/COAST_AATS14_2011_R0.nc',everything=True)


# In[21]:


adict[b'qual_flag']


# # Plot out data

# In[11]:


a[b'Latitude'].shape


# In[13]:


a[b'AOD'].shape


# In[14]:


a[b'wavelength']


# In[15]:


i500 = 3


# In[19]:


plt.figure()
plt.plot(a[b'time'],a[b'GPS_Alt'])


# In[30]:


plt.figure()
plt.plot(a[b'time'],a[b'qual_flag'])


# In[31]:


fl = (a[b'GPS_Alt']<1000.0)# & (a[b'qual_flag']<1)


# In[28]:


a[b'qual_flag']


# In[36]:


adict[b'base_time']


# In[70]:


fig = plt.figure(figsize=(4,3))
proj = ccrs.PlateCarree()#(-120, 45)#(central_latitude=np.nanmean(nav['FMS_LAT']), central_longitude=np.nanmean(nav['FMS_LON']))

ax = fig.add_subplot(111,projection=proj)
ax.plot(a[b'Longitude'],a[b'Latitude'],'k-',lw=0.2)
axs = ax.scatter(a[b'Longitude'][fl],a[b'Latitude'][fl],c=a[b'AOD'][fl,i500],marker='o',s=6,zorder=10)
#ax.scatter(nav['FMS_LON'][ifl],nav['FMS_LAT'][ifl],c='r'',marker='x',s=12)
ax.plot(-122.01721311779664,36.95717375523506,'^',color='tab:orange')
ax.text(-122.01721311779664-0.05,36.95717375523506+0.01,'Santa Cruz Wharf')

ax.coastlines(resolution='10m')
gl = ax.gridlines(draw_labels=True,auto_update=True)
gl.xlabels_top = False
gl.ylabels_right = False


cbar_ax = fig.add_axes([0.1, 0.1, 0.1, 0.1])
fig.subplots_adjust(hspace=0.0, wspace=0, top=1, left=0.09)
posn = ax.get_position()
ax.set_position([posn.x0+0.1, posn.y0+0.15, posn.width-0.1, posn.height-0.3])
posn = ax.get_position()
cbar_ax.set_position([posn.x0 , posn.y0-0.09,
                          posn.width, 0.02])
cbar = plt.colorbar(axs,cax=cbar_ax,label='AOD [500 nm]',orientation='horizontal')

ax.set_title('COAST - AATS-14 on Twin Otter')
plt.savefig(fpp+'CAOST_TO_AOD_map.png',dpi=400,transparent=True)


# In[71]:


np.nanmean(a[b'GPS_Alt'][fl])


# In[72]:


np.nanmedian(a[b'GPS_Alt'][fl])


# In[ ]:


ax.gridlines()

