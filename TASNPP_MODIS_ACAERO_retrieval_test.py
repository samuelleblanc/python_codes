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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-12-05
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
get_ipython().run_line_magic('matplotlib', 'notebook')
import os
from pathlib import Path


# In[9]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature


# In[2]:


name = 'TASNPP'
vv = 'v1'
fp = getpath(name)


# # Load files for ORACLES

# ## Load the MODIS ACAERO

# In[3]:


flist_modis = [str(p) for p in Path(fp+'2016/').rglob('*.nc')]
flist_modis.sort()
flist_modis


# In[4]:


f


# In[53]:


import importlib
importlib.reload(load_utils)


# In[5]:


mod,modd = [],[]
for f in flist_modis:
    print('loading file: '+f)
    mtmp,mtmp_dict = lu.load_hdf(f,all_values=True,verbose=False)
    mod.append(mtmp)
    modd.append(mtmp_dict)


# In[14]:


f


# In[13]:


mod[0]


# In[ ]:





# ## Load the 4STAR for 2016

# ## Load the HSRL for 2016

# # Plot out data

# ## Make maps of MODIS ACAERO

# In[15]:


modd[0]['/geophysical_data/Above_Cloud_AOD']


# In[36]:


def make_acaod(m,md):
    acaod = np.ma.masked_array(m['/geophysical_data/Above_Cloud_AOD'],m['/geophysical_data/Above_Cloud_AOD']==float(md['/geophysical_data/Above_Cloud_AOD']['/geophysical_data/Above_Cloud_AOD#_FillValue']))
    acaod = acaod*float(md['/geophysical_data/Above_Cloud_AOD']['/geophysical_data/Above_Cloud_AOD#scale_factor'])
    return acaod


# In[44]:


fig = plt.figure(figsize=(9,7))
proj = ccrs.PlateCarree()#(-120, 45)#(central_latitude=np.nanmean(nav['FMS_LAT']), central_longitude=np.nanmean(nav['FMS_LON']))
ax = fig.add_subplot(111,projection=proj)
for i in range(16):
    axs = ax.scatter(mod[i]['/geolocation_data/longitude'],mod[i]['/geolocation_data/latitude'],c=make_acaod(mod[i],modd[i]),marker='o',s=6,vmin=0.0,vmax=1.0)
    #ax.scatter(nav['FMS_LON'][ifl],nav['FMS_LAT'][ifl],c='r',marker='x',s=12)

ax.coastlines(resolution='50m')
ax.gridlines(draw_labels=True,auto_update=True)

cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
fig.subplots_adjust(hspace=0.0, wspace=0, top=0.925, left=0.06)
posn = ax.get_position()
ax.set_position([posn.x0, posn.y0, posn.width-0.1, posn.height])
posn = ax.get_position()
cbar_ax.set_position([posn.x0 + posn.width + 0.07, posn.y0,
                          0.02, posn.height])

cbar = plt.colorbar(axs,cax=cbar_ax,label='ACAERO AOD')

   # ax.set_title('ORACLES - NASA P-3 - {}'.format(flist_modis[i]))
   # if i>16: break
plt.savefig(fp+'ORACLES_MODIS_ACAERO_AOD1_map_doy240.png',dpi=400,transparent=True)


# In[ ]:




