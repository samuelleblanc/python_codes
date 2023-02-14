#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     ACCDAM project to look into the trends for GOME, MOPITT, OMI, and TOMS (O3, NO2, CO)
# 
# Input:
# 
#     None
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
#     - write_utils
#     - path_utils
#     - hdf5storage
#     - scipy
# 
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - L3 
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2023-02-08
#     Modified:
# 

# # Prepare python environment

# In[30]:


import numpy as np
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'notebook')
import os
import pandas as pd
from datetime import datetime 


# In[50]:


from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf


# In[6]:


name = 'ACCDAM'
vv = 'v1'
fp = getpath(name)


# Notes on the seasonal decomposition  
# 
# 
# from https://www.statsmodels.org/dev/generated/statsmodels.tsa.seasonal.seasonal_decompose.html  
# 
# Notes  
# 
# This is a naive decomposition. More sophisticated methods should be preferred.  
# 
# The additive model is Y[t] = T[t] + S[t] + e[t]  
# 
# The multiplicative model is Y[t] = T[t] * S[t] * e[t]  
# 
# The results are obtained by first estimating the trend by applying a convolution filter to the data. The trend is then removed from the series and the average of this de-trended series for each period is the returned seasonal component.  

# # Load files

# Load the different L3 files  
# Send from Matt Johnson on Jan 25th, 2023

# ## Load GOME Tropos NO2

# In[11]:


gome,gome_dict = lu.load_netcdf(fp+'GOME_SCIAMACHY_GOME2_NO2_L3/GOME_SCIAMACHY_GOME2ab_TroposNO2_v2.3_041996-092017_temis.nc',everything=True)


# In[28]:


for k in gome: 
    print(k,gome[k].shape)


# In[14]:


gome[b'lon']


# In[15]:


gome[b'time']


# In[17]:


gome[b'TroposNO2'].mean()


# In[44]:


gometime = [datetime((t/100.0).astype(int),(((t/100.0)%1)*100.0).astype(int),15) for t in gome[b'time']]


# ### Convert GOME time series to pandas time for one point to start

# Beijing: 39.9042° N, 116.4074° E  
# San Francisco: 37.7749° N, 122.4194° W

# In[21]:


points = {'Beijing':{'lon':116.40,'lat':39.9},
          'San Francisco': {'lon':-122.42,'lat':37.77}}


# In[26]:


ilon = np.argmin(np.abs(gome[b'lon']-points['Beijing']['lon']))
ilat = np.argmin(np.abs(gome[b'lat']-points['Beijing']['lat']))
ilon,ilat


# In[48]:


gome_NO2_pd = pd.DataFrame(data=gome[b'TroposNO2'][:,ilat,ilon])
gome_NO2_pd['time'] = pd.to_datetime(gometime)
gome_NO2_pd['trop_NO2'] = gome[b'TroposNO2'][:,ilat,ilon]
gome_NO2_pd.set_index('time', inplace = True)
gome_NO2_pd.dropna(inplace=True)


# In[52]:


gno2 = gome_NO2_pd.resample('M').mean()
gno2.interpolate(inplace=True)


# ## Load MOPITT
# Should have 235 files

# In[62]:


fp


# In[74]:


fp_mol = os.listdir(fp+'MOPITT/')
fp_mo = [f for f in fp_mol if '.he5' in f]
fp_mo.sort()


# In[75]:


fp_mo


# In[133]:


mops,mop_dicts = [],[]
for f in fp_mo:
    print('Opening file: ', f)
    mop,mop_dict = lu.load_hdf(fp+'MOPITT/'+f,values=(('COday',20),('COnight',23)),verbose=False)
    mops.append(mop)
    mop_dicts.append(mop_dict)


# In[130]:


import h5py
f5 = h5py.File(fp+'MOPITT/'+f)
mop_lat = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Latitude'])
mop_lon = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Longitude'])
mop_pre = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Pressure'])


# In[80]:


mop_dict


# # Set the different regions

# From Matt Johnson on 8 février 2023 11:51  
# 
# Hi Sam,
# 
#  
# 
# Here are the lat/lon bounds we discussed on our meeting today.  I’ve cc’d others from the project which will be interested in this.  I have taken a first pass on defining lat/lon bounds for the regions of: Southeast Asia, China, South Asia, Siberia, and trop Pacific Ocean.  The first three I used HTAP countries to define these bounds.  This is preliminary but at least we can start looking at O3, NO2, and CO trends from the satellite data.  We can easily adjust the lat/lon bounds based on future discussions.
# 
#  
# 
# Regions based on HTAP2 countries:  
# Southeast Asia: 12°S–18°N, 95–140°E  
# China: 21–46°N, 75–127°E  
# South Asia: 6°N–31°N, 67–91°E  
# Regions not based off HTAP2:  
# Siberia: 50–75°N, 70–180°E  
# Tropical Pacific Ocean: 5–35°N, 180–130°W  
# 
# US regions from proposal:  
#  
# 
# |     | Color in Fig. 2 | Lat (°N) | Lon (°W) |  
# | --- | --------------- | -------- | -------- |  
# |Entire Domain | Base Map | 25 – 55 | 130 – 90 |  
# |Coastal | Teal | 32 – 55 | 130 – 113 |  
# | Northwest NA| Purple | 48 – 54 | 130 - 115 |
# | Northern California | Orange | 35 – 41 | 126 – 120 | 
# | Southern California | Olive Green | 35 – 32 | 120 – 113 | 
# | Great Basin | Pink | 36 – 41 | 120 – 115 | 
# 
#    
#    
#    ![HTAP_regions_v01.jpg](attachment:HTAP_regions_v01.jpg)

# ## Definition of the regions

# In[83]:


rgs = {'Southeast Asia':[[-12,95],[18,140]],
       'China':[[21,75],[46,127]],
       'South Asia':[[6,67],[31,91]],
       'Siberia':[[50,70],[75,180]],
       'Tropical Pacific Ocean':[[5,-180],[35,-130]]     
      } #lower left [lat lon], upper right [lat lon]


# In[87]:


rgs['China'][0]


# In[144]:


def multi_stats_pd(data,time,name='trop_NO2',axis=1):
    'to get a dict for the different dataframes (mean, median, min, max)'
    pds = {}
    pds['mean'] = build_pd(np.nanmean(np.nanmean(data,axis=axis),axis=axis),time,name='mean_'+name)
    pds['median'] = build_pd(np.nanmedian(np.nanmedian(data,axis=axis),axis=axis),time,name='median_'+name)
    pds['min'] = build_pd(np.nanmin(np.nanmin(data,axis=axis),axis=axis),time,name='min_'+name)
    pds['max'] = build_pd(np.nanmax(np.nanmax(data,axis=axis),axis=axis),time,name='max_'+name)
    pds['std'] = build_pd(np.nanstd(np.nanstd(data,axis=axis),axis=axis),time,name='std_'+name)
    
    return pds
    


# In[88]:


def build_pd(data,time,name='mean_trop_NO2'):
    'To prepare the dataframe with time settings and regular interpolation, ready for seasonal decomposition'
    data_pd = pd.DataFrame(data=data)
    data_pd['time'] = pd.to_datetime(time)
    data_pd[name] = data
    data_pd.set_index('time', inplace = True)
    data_pd.dropna(inplace=True)
    dat = data_pd.resample('M').mean()
    dat.interpolate(inplace=True)
    return dat


# ## Subset for GOME NO2

# In[145]:


gome_rg = {}
for rg in rgs: 
    print(rg)
    ill_lon = np.argmin(np.abs(gome[b'lon']-rgs[rg][0][1]))
    ill_lat = np.argmin(np.abs(gome[b'lat']-rgs[rg][0][0]))
    iur_lon = np.argmin(np.abs(gome[b'lon']-rgs[rg][1][1]))
    iur_lat = np.argmin(np.abs(gome[b'lat']-rgs[rg][1][0]))
    
    gome_rg[rg] = multi_stats_pd(gome[b'TroposNO2'][:,ill_lat:iur_lat,ill_lon:iur_lon],gometime,name='GOME_tropNO2')
    


# In[148]:


gome_rg['China'].keys()


# In[151]:


gome_rg['China']['mean']['mean_GOME_tropNO2']


# # Plot out data

# ## GOME Tropospheric NO2 one location

# In[49]:


gome_NO2_pd.plot()


# Run the stats models

# In[58]:


result_gome = seasonal_decompose(gno2['trop_NO2'])
p = result_gome.plot()


# In[57]:


fig, ax = plt.subplots(figsize = (4,2))
plot_acf(gno2['trop_NO2'], ax = ax)
plt.show()


# ## GOME NO2 regions

# In[161]:


for rg in gome_rg:
    for typ in ['mean', 'median', 'min', 'max', 'std']:
        result_gome = seasonal_decompose(gome_rg[rg][typ][typ+'_GOME_tropNO2'])
        p = result_gome.plot()
        p.get_axes()[0].set_title(rg+ ': '+typ+'_GOME_tropNO2')
        p.get_axes()[0].set_ylabel('All')


# In[ ]:




