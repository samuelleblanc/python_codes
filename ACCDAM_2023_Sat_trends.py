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

# In[1]:


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
from datetime import datetime, timedelta


# In[2]:


from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf


# In[3]:


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


# In[166]:


len(mops)


# In[179]:


mop_COday = np.array([m['COday']  for m in mops])
mop_COnight = np.array([m['COnight']  for m in mops])
mop_COday.shape, mop_COnight.shape


# In[199]:


mop_COday[mop_COday == -9999] = np.nan
mop_COnight[mop_COnight == -9999] = np.nan


# In[176]:


mop_time = [datetime(int(f.split('-')[1][0:4]),int(f.split('-')[1][4:6]),15)  for f in fp_mo]


# In[130]:


import h5py
f5 = h5py.File(fp+'MOPITT/'+f)
mop_lat = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Latitude'])
mop_lon = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Longitude'])
mop_pre = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Pressure'])


# In[193]:


mop_dict


# ## Load TOMS O3

# In[5]:


fp


# In[6]:


toms_pd = pd.read_excel(fp+'TOMS_O3_L3/TOMS_monthly_trop_O3_L3_v01.xlsx')


# In[19]:


toms_ar = np.array(toms_pd)


# In[25]:


toms_ar[0:73,0]


# In[56]:


start_time = datetime(1978,12,31)
end_time = datetime(2005,12,15)
toms_time = pd.date_range(start=start_time,end=end_time,freq='MS') #jan1979_to_dec2005


# In[78]:


ntime = 27*12
nlon = 72
nlat = 12
toms_O3 = np.zeros((nlat,nlon,ntime)) #lat, lon, time
for n in range(ntime):
    toms_O3[:,:,n] = toms_ar[(n*(nlon+2)):(n*(nlon+2))+nlon,1:].T


# In[82]:


toms_O3[toms_O3>900] = np.nan


# In[85]:


toms_lat = np.arange(-27.5,27.6,5)
toms_lon = np.arange(-177.5,177.6,5)


# ## Load OMI

# ### Load MLS O3

# In[12]:


fp
fpp_mlso3 = fp+'OMI_MLS_O3_L3/OMI_MLS_O3_L3_ncfiles/'


# In[13]:


fp_mls_o3 = os.listdir(fpp_mlso3)
fp_mls_o3 = [f for f in fp_mls_o3 if '.nc' in f]
fp_mls_o3.sort()


# In[10]:


fp_mls_o3


# In[14]:


mlso3_tmp,mlso3_tmp_dict = lu.load_netcdf(fpp_mlso3+fp_mls_o3[0],everything=True)


# In[20]:


mlso3_tmp[b'tropo3'].shape, mlso3_tmp[b'lon'].shape[0], mlso3_tmp[b'lat'].shape[0]


# In[16]:


ntime = len(fp_mls_o3)


# In[24]:


mls_o3 = {'lat':mlso3_tmp[b'lat'],'lon':mlso3_tmp[b'lon'],'time':[],'tropo3':np.zeros((ntime,mlso3_tmp[b'lat'].shape[0],mlso3_tmp[b'lon'].shape[0]))}


# In[43]:


mls_o3['time'] = []
for i,f in list(enumerate(fp_mls_o3)):
    mlso3_tmp,mlso3_tmp_dict = lu.load_netcdf(fpp_mlso3+f,everything=True)
    monthstr = f.split('_')[2]
    mls_o3['time'].append(datetime(int(monthstr[0:4]),int(monthstr[4:6]),15))
    mls_o3['tropo3'][i,:,:] = mlso3_tmp[b'tropo3']


# In[44]:


mls_o3_time = np.array(mls_o3['time'])


# In[51]:


mlso3_tmp_dict[b'tropo3']


# In[54]:


mls_o3['tropo3'] = mls_o3['tropo3']/10.0


# ### ** Load OMI NO2 L2

# In[57]:


fp
fpp_mlsno2 = fp+'OMI_NO2/'


# In[58]:


fp_mls_no2 = os.listdir(fpp_mlsno2)
fp_mls_no2 = [f for f in fp_mls_no2 if '.mat' in f]
fp_mls_no2.sort()


# In[59]:


fp_mls_no2


# In[61]:


mlsno2_tmp = sio.loadmat(fpp_mlsno2+fp_mls_no2[0])


# In[63]:


mlsno2_tmp.keys()


# In[65]:


mlsno2_tmp['NO2_monthly_avg'].shape, mlsno2_tmp['LAT'].shape, mlsno2_tmp['LON'].shape


# In[66]:


ntime = len(fp_mls_no2)*12


# In[72]:


mls_no2 = {'lat':mlsno2_tmp['LAT'],'lon':mlsno2_tmp['LON'],'time':[],'no2':np.zeros((ntime,mlsno2_tmp['LAT'].shape[0],mlsno2_tmp['LON'].shape[1]))}


# In[74]:


mls_no2['time'] = []
for i,f in list(enumerate(fp_mls_no2)):
    mlsno2_tmp = sio.loadmat(fpp_mlsno2+f)
    yearstr = f.split('_')[4]
    [mls_no2['time'].append(datetime(int(yearstr[0:4]),m,15)) for m in range(1,13)]
    mls_no2['no2'][i*12:(i+1)*12,:,:] = mlsno2_tmp['NO2_monthly_avg']


# In[44]:


mls_no2_time = np.array(mls_no2['time'])


# In[51]:


mlsno2_tmp_dict[b'tropno2']


# In[54]:


mls_no2['tropno2'] = mls_no2['tropno2']/10.0


# ## Load CEDS NOx

# ### Aircraft

# In[112]:


ceds_1980,ceds_1980_dict = lu.load_netcdf(fp+'CEDS/CEDS_NOx_aircraft_monthly_1980_1999.nc',everything=True)


# In[113]:


ceds_2000,ceds_2000_dict = lu.load_netcdf(fp+'CEDS/CEDS_NOx_aircraft_monthly_2000_2019.nc',everything=True)


# In[126]:


ceds_2000_dict[b'Time']


# In[114]:


ceds_2000[b'CEDS_NOx_aircraft_emission'].shape


# In[115]:


ceds_2000[b'Lat'].shape, ceds_2000[b'Lon'].shape, ceds_2000[b'Time'].shape


# In[116]:


ceds_1980[b'Lat'].shape, ceds_1980[b'Lon'].shape, ceds_1980[b'Time'].shape


# In[136]:


ceds_AC_NOx = np.vstack((ceds_1980[b'CEDS_NOx_aircraft_emission'],ceds_2000[b'CEDS_NOx_aircraft_emission']))
ceds_time_days = np.hstack((ceds_1980[b'Time'],ceds_2000[b'Time']))


# In[123]:


ceds_AC_NOx.shape


# In[137]:


ceds_time = np.array([datetime(1750,1,1)+timedelta(days=int(d)) for d in ceds_time_days])


# In[140]:


ceds_lon = ceds_2000[b'Lon']
ceds_lat = ceds_2000[b'Lat']


# ### anthropogenic

# In[138]:


ceds_1980a,ceds_1980a_dict = lu.load_netcdf(fp+'CEDS/CEDS_NOx_anthro_monthly_1980_1999.nc',everything=True)
ceds_2000a,ceds_2000a_dict = lu.load_netcdf(fp+'CEDS/CEDS_NOx_anthro_monthly_2000_2019.nc',everything=True)


# In[139]:


ceds_Ant_NOx = np.vstack((ceds_1980a[b'CEDS_NOx_anthro_emission'],ceds_2000a[b'CEDS_NOx_anthro_emission']))


# ## Load TCR

# In[4]:


fp


# In[5]:


fp_tcrl = os.listdir(fp+'TCR/TCR_2_data/')
fp_tcr = [f for f in fp_tcrl if '.nc' in f]
fp_tcr.sort()


# In[6]:


fp_tcrl


# In[17]:


tcr_vals = ['_'.join([f.split('_')[1],f.split('_')[3]]) for f in fp_tcrl]


# In[15]:


tcr_tmp,tcr_tmp_dict = lu.load_netcdf(fp+'TCR/TCR_2_data/'+fp_tcrl[0],everything=True)


# In[16]:


tcr = {'lat':tcr_tmp[b'lat'],'lon':tcr_tmp[b'lon'],'time':tcr_tmp[b'time']}


# In[18]:


tcr_vals


# In[21]:


for i,f in list(enumerate(fp_tcrl)):
    tcr_tmp,tcr_tmp_dict = lu.load_netcdf(fp+'TCR/TCR_2_data/'+f,everything=True)
    tcr[tcr_vals[i]] = tcr_tmp[tcr_vals[i].split('_')[0].encode('utf-8')]


# In[24]:


tcr_tmp_dict[b'time']


# In[38]:


tcr_time = np.array([datetime(2005+int(m//12),1+int(m%12),15) for m in tcr['time']])


# In[63]:


# convert the 0 to 360 into a -180 to 180 longitude
tcr['lon'][tcr['lon']>=180.0] = tcr['lon'][tcr['lon']>=180.0] - 360.0


# In[81]:


tcr_tmp_dict[b'co']


# In[ ]:





# In[84]:


plt.figure()
plt.hist(tcr['co_bio'].flatten(),bins=30)
plt.yscale('log')


# In[85]:


plt.figure()
plt.hist(tcr['co_anth'].flatten(),bins=30)
plt.yscale('log')


# ## ***Load GFED CO

# In[ ]:





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

# In[37]:


rgs = {'Southeast Asia':[[-12,95],[18,140]],
       'China':[[21,75],[46,127]],
       'South Asia':[[6,67],[31,91]],
       'Siberia':[[50,70],[75,180]],
       'Tropical Pacific Ocean':[[5,-180],[35,-130]]     
      } #lower left [lat lon], upper right [lat lon]


# In[38]:


rgs['China'][0]


# In[39]:


def multi_stats_pd(data,time,name='trop_NO2',axis=1):
    'to get a dict for the different dataframes (mean, median, min, max)'
    pds = {}
    pds['mean'] = build_pd(np.nanmean(np.nanmean(data,axis=axis),axis=axis),time,name='mean_'+name)
    pds['median'] = build_pd(np.nanmedian(np.nanmedian(data,axis=axis),axis=axis),time,name='median_'+name)
    pds['min'] = build_pd(np.nanmin(np.nanmin(data,axis=axis),axis=axis),time,name='min_'+name)
    pds['max'] = build_pd(np.nanmax(np.nanmax(data,axis=axis),axis=axis),time,name='max_'+name)
    pds['std'] = build_pd(np.nanstd(np.nanstd(data,axis=axis),axis=axis),time,name='std_'+name)
    
    return pds
    


# In[40]:


def build_pd(data,time,name='mean_trop_NO2'):
    'To prepare the dataframe with time settings and regular interpolation, ready for seasonal decomposition'
    data_pd = pd.DataFrame(data=data)
    data_pd['time'] = pd.to_datetime(time)
    data_pd[name] = data
    data_pd.set_index('time', inplace = True)
    data_pd.dropna(inplace=True)
    # return how many bad data in here
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


# ## Subset for MOPITT (pressure level and regions)

# Should change timeline to start at 2002:  
# 
# Rebecca bucholtz, et al., 2021, https://doi.org/10.1016/j.rse.2020.112275  
# 

# In[200]:


mop_pre


# In[201]:


mop_lon = np.array(mop_lon)
mop_lat = np.array(mop_lat)
mop_lon.shape, mop_lat.shape, mop_COday.shape


# In[208]:


mopitt_day_rg = {}
mopitt_night_rg = {}
for rg in rgs: 
    print(rg)
    ill_lon = np.argmin(np.abs(mop_lon-rgs[rg][0][1]))
    ill_lat = np.argmin(np.abs(mop_lat-rgs[rg][0][0]))
    iur_lon = np.argmin(np.abs(mop_lon-rgs[rg][1][1]))
    iur_lat = np.argmin(np.abs(mop_lat-rgs[rg][1][0]))
    mopitt_day_rg[rg] = {}
    mopitt_night_rg[rg] = {}
    for ip,pr in list(enumerate(mop_pre)):
        print(pr)
        mopitt_day_rg[rg]['{:3.0f}'.format(pr)] = multi_stats_pd(mop_COday[:,ill_lon:iur_lon,ill_lat:iur_lat,ip],mop_time,name='MOPITT_dayCO_{:3.0f}'.format(pr))
        mopitt_night_rg[rg]['{:3.0f}'.format(pr)] = multi_stats_pd(mop_COnight[:,ill_lon:iur_lon,ill_lat:iur_lat,ip],mop_time,name='MOPITT_nightCO_{:3.0f}'.format(pr))


# In[203]:


mopitt_day_rg['China']


# ## Subset for TOMS O3

# In[107]:


toms_rg = {}
for rg in rgs: 
    if rg == 'Siberia': continue # no TOMS data for Siberia
    print(rg)
    
    ill_lon = np.argmin(np.abs(toms_lon-rgs[rg][0][1]))
    ill_lat = np.argmin(np.abs(toms_lat-rgs[rg][0][0]))
    iur_lon = np.argmin(np.abs(toms_lon-rgs[rg][1][1]))
    iur_lat = np.argmin(np.abs(toms_lat-rgs[rg][1][0]))
    
    toms_rg[rg] = multi_stats_pd(toms_O3[ill_lat:iur_lat,ill_lon:iur_lon,:],toms_time,axis=0,name='TOMS_O3')
    


# ## Subset for OMI MLS

# ### Ozone

# In[55]:


mlso3_rg = {}
for rg in rgs: 
    #if rg == 'Siberia': continue # no TOMS data for Siberia
    print(rg)
    
    ill_lon = np.argmin(np.abs(mls_o3['lon']-rgs[rg][0][1]))
    ill_lat = np.argmin(np.abs(mls_o3['lat']-rgs[rg][0][0]))
    iur_lon = np.argmin(np.abs(mls_o3['lon']-rgs[rg][1][1]))
    iur_lat = np.argmin(np.abs(mls_o3['lat']-rgs[rg][1][0]))
    
    mlso3_rg[rg] = multi_stats_pd(mls_o3['tropo3'][:,ill_lat:iur_lat,ill_lon:iur_lon],mls_o3_time,axis=1,name='OMI_MLS_O3')
    


# ## Subset for CEDS NOx

# In[141]:


ceds_ac_rg = {}
ceds_ant_rg = {}
for rg in rgs: 
    print(rg)
    
    ill_lon = np.argmin(np.abs(ceds_lon-rgs[rg][0][1]))
    ill_lat = np.argmin(np.abs(ceds_lat-rgs[rg][0][0]))
    iur_lon = np.argmin(np.abs(ceds_lon-rgs[rg][1][1]))
    iur_lat = np.argmin(np.abs(ceds_lat-rgs[rg][1][0]))
    
    ceds_ac_rg[rg] = multi_stats_pd(ceds_AC_NOx[:,ill_lat:iur_lat,ill_lon:iur_lon],ceds_time,axis=1,name='CEDS_AC_NOx')
    ceds_ant_rg[rg] = multi_stats_pd(ceds_Ant_NOx[:,ill_lat:iur_lat,ill_lon:iur_lon],ceds_time,axis=1,name='CEDS_Ant_NOx')
    


# ## Subset for TCR NOx/CO

# In[43]:


tcr_vals


# In[46]:


tcr['nox_anth'].shape, tcr['lat'].shape, tcr['lon'].shape, tcr['time'].shape


# In[50]:


ill_lon, iur_lon, ill_lat, iur_lat


# In[51]:


tcr[val][:,ill_lat:iur_lat,ill_lon:iur_lon]


# In[52]:


val


# In[65]:


tcr_rg = {}
first = True
for rg in rgs: 
    print(rg)
    
    ill_lon = np.argmin(np.abs(tcr['lon']-rgs[rg][0][1]))
    ill_lat = np.argmin(np.abs(tcr['lat']-rgs[rg][0][0]))
    iur_lon = np.argmin(np.abs(tcr['lon']-rgs[rg][1][1]))
    iur_lat = np.argmin(np.abs(tcr['lat']-rgs[rg][1][0]))
    
    for val in tcr_vals:
        if first: tcr_rg[val] = {}
        tcr_rg[val][rg] = multi_stats_pd(tcr[val][:,ill_lat:iur_lat,ill_lon:iur_lon],tcr_time,axis=1,name='TCR'+val)
    first = False
    


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


# ## MOPITT CO regions

# In[213]:


mop_pre_lbl = ['{:3.0f}'.format(p) for p in mop_pre ]


# In[214]:


for pr in mop_pre_lbl:
    for rg in mopitt_day_rg:
        for typ in ['median']:
            result_mop = seasonal_decompose(mopitt_day_rg[rg][pr][typ][typ+'_MOPITT_dayCO_'+pr])
            p = result_mop.plot()
            p.get_axes()[0].set_title(rg+ ': '+typ+'_MOPITT_dayCO_'+pr)
            p.get_axes()[0].set_ylabel('All')


# In[209]:


for rg in mopitt_day_rg:
    for typ in ['mean', 'median', 'min', 'max', 'std']:
        result_mop = seasonal_decompose(mopitt_day_rg[rg]['800'][typ][typ+'_MOPITT_dayCO_800'])
        p = result_mop.plot()
        p.get_axes()[0].set_title(rg+ ': '+typ+'_MOPITT_dayCO_800')
        p.get_axes()[0].set_ylabel('All')


# In[207]:


mopitt_day_rg['China']


# In[206]:


for rg in mopitt_day_rg:
    for typ in ['mean', 'median', 'min', 'max', 'std']:
        result_mop = seasonal_decompose(mopitt_day_rg[rg][typ][typ+'_MOPITT_dayCO_800'])
        p = result_mop.plot()
        p.get_axes()[0].set_title(rg+ ': '+typ+'_MOPITT_dayCO_800')
        p.get_axes()[0].set_ylabel('All')


# ## TOMS O3 

# In[109]:


for rg in toms_rg:
    for typ in ['mean', 'median', 'std']:
        result_toms = seasonal_decompose(toms_rg[rg][typ][typ+'_TOMS_O3'])
        p = result_toms.plot()
        p.get_axes()[0].set_title(rg+ ': '+typ+'_TOMS_O3')
        p.get_axes()[0].set_ylabel('All')


# ## OMI MLS

# ### O3

# In[79]:


type(result_mlso3.trend)


# In[83]:


mlso3_rg['Southeast Asia']['mean']


# In[56]:


for rg in mlso3_rg:
    for typ in ['mean', 'median', 'std']:
        result_mlso3 = seasonal_decompose(mlso3_rg[rg][typ][typ+'_OMI_MLS_O3'])
        p = result_mlso3.plot()
        p.get_axes()[0].set_title(rg+ ': '+typ+'_OMI_MLS_O3')
        p.get_axes()[0].set_ylabel('All')
        mlso3_rg[rg][typ]['trend'] = p.trend


# ## CEDS NOx

# ### Aircraft

# In[142]:


for rg in ceds_ac_rg:
    for typ in ['mean', 'median', 'std']:
        result_ceds_ac = seasonal_decompose(ceds_ac_rg[rg][typ][typ+'_CEDS_AC_NOx'])
        p = result_ceds_ac.plot()
        p.get_axes()[0].set_title(rg+ ': '+typ+'_CEDS_AC_NOx')
        p.get_axes()[0].set_ylabel('All')


# ### Anthropogenic

# In[143]:


for rg in ceds_ant_rg:
    for typ in ['mean', 'median', 'std']:
        result_ceds_ant = seasonal_decompose(ceds_ant_rg[rg][typ][typ+'_CEDS_Ant_NOx'])
        p = result_ceds_ant.plot()
        p.get_axes()[0].set_title(rg+ ': '+typ+'_CEDS_Ant_NOx')
        p.get_axes()[0].set_ylabel('All')


# ## TCR

# In[66]:


tcr_vals


# In[71]:


tcr_rg['nox_anth']['Southeast Asia']['mean']


# In[ ]:


for val in tcr_vals:
    for rg in tcr_rg[val]:
        for typ in ['mean', 'median', 'std']:
            result_tcr = seasonal_decompose(tcr_rg[val][rg][typ][typ+'_TCR'+val])
            p = result_tcr.plot()
            p.get_axes()[0].set_title(rg+ ': '+typ+'_TCR'+val)
            p.get_axes()[0].set_ylabel('All')


# In[ ]:





# # Get the time series for the different percentile ranges

# In[86]:


pcts = [5.0,33.0,50.0,66.0,95.0]


# In[ ]:


datp = np.percentile(dat,pcts)


# In[ ]:





# # Make figures of the trends for each points

# ## Load the functions

# In[136]:


import cartopy.crs as ccrs
import statsmodels.api as sm
import cartopy


# In[221]:


def trends_and_pval_per_lat_lon(data,lat,lon,time,name='mean_trop_NO2'):
    "Find the linear trend and its pval for each lat/lon point"
    nlat = len(lat)
    nlon = len(lon)
    
    # Create an empty array to hold the trend and p-value for each lat and lon point
    trend = np.empty((nlat, nlon))
    trend_pval = np.empty((nlat, nlon))
    for i in range(nlat):
        for j in range(nlon):
            # Extract the time series for this lat and lon point
            ts = build_pd(data[:, i, j],time,name=name)
            ts.dropna()

            # Perform seasonal decomposition on the time series
            result = seasonal_decompose(ts[name]) #, model='additive', period=12)

            # Extract the trend component from the decomposition
            trend_comp = result.trend
            trend_comp.dropna()

            # Calculate the straight trend using the OLS method
            fitl = sm.OLS(trend_comp, sm.add_constant(range(len(trend_comp))),missing='drop').fit()
            trend[i, j] = fitl.params[1]

            # Calculate the p-value for the trend using a t-test
            trend_pval[i, j] = fitl.pvalues[1]

    return trend,trend_pval


# In[233]:


def plot_trend_and_pval(trend,trend_pval,lon,lat,name='',cax_name='Trend',figsize=(10,4)):
    
    pr = ccrs.PlateCarree()
    # Set up the plot
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=pr)

    # Define the colormap
    cmap = plt.cm.get_cmap('RdBu_r', 12)

    # Create a filled contour plot of the trend
    contour_levels = np.arange(-0.025, 0.04, 0.005)
    contour = ax.contourf(lon, lat, trend, contour_levels,
                          cmap=cmap, transform=pr, extend='both')

    # Add a colorbar
    cbar = plt.colorbar(contour, shrink=0.6, pad=0.02)
    cbar.ax.set_ylabel(cax_name, fontsize=16)

    # Add a land mask
    ax.add_feature(cartopy.feature.LAND, edgecolor='black', facecolor='None', alpha=0.3,zorder=-100)

    # Add a grid
    gl = ax.gridlines(color='gray', linestyle='--',draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False

    # Add a '+' symbol for statistically significant trends
    significant = np.where(trend_pval < 0.05)
    ax.scatter(lon[significant[1]][0::5], lat[significant[0]][0::5], marker='+', color='k', transform=pr,s=0.5)

    # Add a title
    plt.title('Trend for '+name, fontsize=18)
    ax.add_feature(cartopy.feature.COASTLINE)
    return fig


# In[223]:


def span_pd_lat_lon(data,time,lat,lon,name='OMI_MLS_O3',axis=1):
    "build the pandas time series for each point in lat lon"
    nlat = len(lat)
    nlon = len(lon)
    ntime  = len(time)
    dat = np.empty((ntime,nlat,nlon))
    for i in range(nlat):
        for j in range(nlon):
            dat[:,i,j] = build_pd(data[:,i,j],time,name=name)
    return dat


# ## Plot the GOME NO2

# In[ ]:


gome[b'lon']
gome[b'lat']

gome_rg[rg] = multi_stats_pd(gome[b'TroposNO2'][:,ilat,ilon],gometime,name='GOME_tropNO2')


# ##  Plot the OMI MLS

# ### O3

# In[224]:


mls_o3_trend,mls_o3_trend_pval = trends_and_pval_per_lat_lon(mls_o3['tropo3'],mls_o3['lat'],mls_o3['lon'],mls_o3_time,name='OMI_MLS_O3')


# In[228]:


mls_o3_time[0],mls_o3_time[-1]


# In[237]:


fig = plot_trend_and_pval(mls_o3_trend,mls_o3_trend_pval,mls_o3['lon'],mls_o3['lat'],
                          name=name+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(mls_o3_time[0],mls_o3_time[-1]),
                          cax_name='O3 trend [DU/month]',figsize=(10,4))
fig.tight_layout()


# ### NO2

# In[ ]:




