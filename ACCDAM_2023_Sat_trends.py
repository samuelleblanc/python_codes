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
import pickle


# In[2]:


from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf


# In[3]:


import cartopy.crs as ccrs
import statsmodels.api as sm
import cartopy


# In[4]:


import traceback
from tqdm.notebook import tqdm 
from functools import partial
from IPython.core.debugger import set_trace


# In[5]:


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

# In[6]:


gome,gome_dict = lu.load_netcdf(fp+'GOME_SCIAMACHY_GOME2_NO2_L3/GOME_SCIAMACHY_GOME2ab_TroposNO2_v2.3_041996-092017_temis.nc',everything=True)


# In[7]:


for k in gome: 
    print(k,gome[k].shape)


# In[8]:


gome[b'lon']


# In[9]:


gome[b'time']


# In[10]:


gome[b'TroposNO2'].mean()


# In[11]:


gometime = [datetime((t/100.0).astype(int),(((t/100.0)%1)*100.0).astype(int),15) for t in gome[b'time']]


# ### Convert GOME time series to pandas time for one point to start

# Beijing: 39.9042° N, 116.4074° E  
# San Francisco: 37.7749° N, 122.4194° W

# In[12]:


points = {'Beijing':{'lon':116.40,'lat':39.9},
          'San Francisco': {'lon':-122.42,'lat':37.77}}


# In[13]:


ilon = np.argmin(np.abs(gome[b'lon']-points['Beijing']['lon']))
ilat = np.argmin(np.abs(gome[b'lat']-points['Beijing']['lat']))
ilon,ilat


# In[14]:


gome_NO2_pd = pd.DataFrame(data=gome[b'TroposNO2'][:,ilat,ilon])
gome_NO2_pd['time'] = pd.to_datetime(gometime)
gome_NO2_pd['trop_NO2'] = gome[b'TroposNO2'][:,ilat,ilon]
gome_NO2_pd.set_index('time', inplace = True)
gome_NO2_pd.dropna(inplace=True)


# In[15]:


gno2 = gome_NO2_pd.resample('M').mean()
gno2.interpolate(inplace=True)


# ## Load MOPITT
# Should have 235 files

# In[107]:


fp


# In[108]:


fp_mol = os.listdir(fp+'MOPITT/')
fp_mo = [f for f in fp_mol if '.he5' in f]
fp_mo.sort()


# In[109]:


fp_mo


# In[110]:


mop,mop_dict = lu.load_hdf(fp+'MOPITT/'+fp_mo[0],verbose=True,values=(('COday',20),('COnight',23)))


# In[111]:


f5 = h5py.File(fp+'MOPITT/'+fp_mo[0])
mop_COday_tmp = np.array(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['RetrievedCOMixingRatioProfileDay'])
mop_COnight_tmp = np.array(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['RetrievedCOMixingRatioProfileNight'])
f5.close()


# In[112]:


mop_COday_tmp.shape


# In[113]:


nlon_mop,nlat_mop,npres_mop = mop_COday_tmp.shape


# In[114]:


ntime_moppit = len(fp_mo)


# In[115]:


mop_COday = np.zeros((ntime_moppit,nlon_mop,nlat_mop,npres_mop))+np.nan
mop_COnight = np.zeros((ntime_moppit,nlon_mop,nlat_mop,npres_mop))+np.nan
for i,f in list(enumerate(fp_mo)):
    print('Opening file: ', f)
    f5 = h5py.File(fp+'MOPITT/'+f)
    mop_COday_tmp = np.array(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['RetrievedCOMixingRatioProfileDay'])
    mop_COnight_tmp = np.array(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['RetrievedCOMixingRatioProfileNight'])
    f5.close()
    mop_COday[i,:,:,:] = mop_COday_tmp
    mop_COnight[i,:,:,:] = mop_COnight_tmp


# In[116]:


mop_COday[mop_COday == -9999] = np.nan
mop_COnight[mop_COnight == -9999] = np.nan


# In[117]:


mop_time = [datetime(int(f.split('-')[1][0:4]),int(f.split('-')[1][4:6]),15)  for f in fp_mo]


# In[118]:


import h5py
f5 = h5py.File(fp+'MOPITT/'+f)
mop_lat = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Latitude'])
mop_lon = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Longitude'])
mop_pre = list(f5['HDFEOS']['GRIDS']['MOP03']['Data Fields']['Pressure'])
f5.close()


# ## Load TOMS O3

# In[6]:


fp


# In[9]:


toms_pd = pd.read_excel(fp+'TOMS_O3_L3/TOMS_monthly_trop_O3_L3_v01.xlsx')


# In[10]:


toms_ar = np.array(toms_pd)


# In[11]:


toms_ar[0:73,0]


# In[12]:


start_time = datetime(1978,12,31)
end_time = datetime(2005,12,15)
toms_time = pd.date_range(start=start_time,end=end_time,freq='MS') #jan1979_to_dec2005


# In[13]:


ntime = 27*12
nlon = 72
nlat = 12
toms_O3 = np.zeros((nlat,nlon,ntime)) #lat, lon, time
for n in range(ntime):
    toms_O3[:,:,n] = toms_ar[(n*(nlon+2)):(n*(nlon+2))+nlon,1:].T


# In[14]:


toms_O3[toms_O3>900] = np.nan


# In[15]:


toms_lat = np.arange(-27.5,27.6,5)
toms_lon = np.arange(-177.5,177.6,5)


# ## Load OMI

# ### Load MLS O3

# In[6]:


fp
fpp_mlso3 = fp+'OMI_MLS_O3_L3/OMI_MLS_O3_L3_ncfiles/'


# In[7]:


fp_mls_o3 = os.listdir(fpp_mlso3)
fp_mls_o3 = [f for f in fp_mls_o3 if '.nc' in f]
fp_mls_o3.sort()


# In[8]:


fp_mls_o3


# In[9]:


mlso3_tmp,mlso3_tmp_dict = lu.load_netcdf(fpp_mlso3+fp_mls_o3[0],everything=True)


# In[10]:


mlso3_tmp[b'tropo3'].shape, mlso3_tmp[b'lon'].shape[0], mlso3_tmp[b'lat'].shape[0]


# In[11]:


ntime = len(fp_mls_o3)


# In[12]:


mls_o3 = {'lat':mlso3_tmp[b'lat'],'lon':mlso3_tmp[b'lon'],'time':[],'tropo3':np.zeros((ntime,mlso3_tmp[b'lat'].shape[0],mlso3_tmp[b'lon'].shape[0]))}


# In[13]:


mls_o3['time'] = []
for i,f in list(enumerate(fp_mls_o3)):
    mlso3_tmp,mlso3_tmp_dict = lu.load_netcdf(fpp_mlso3+f,everything=True)
    monthstr = f.split('_')[2]
    mls_o3['time'].append(datetime(int(monthstr[0:4]),int(monthstr[4:6]),15))
    mls_o3['tropo3'][i,:,:] = mlso3_tmp[b'tropo3']


# In[14]:


mls_o3_time = np.array(mls_o3['time'])


# In[15]:


mlso3_tmp_dict[b'tropo3']


# In[16]:


mls_o3['tropo3'] = mls_o3['tropo3']/10.0


# ### Load OMI NO2 L2

# In[62]:


fp
fpp_mlsno2 = fp+'OMI_NO2/'


# In[63]:


fp_mls_no2 = os.listdir(fpp_mlsno2)
fp_mls_no2 = [f for f in fp_mls_no2 if '.mat' in f]
fp_mls_no2.sort()


# In[64]:


fp_mls_no2


# In[65]:


mlsno2_tmp = sio.loadmat(fpp_mlsno2+fp_mls_no2[0])


# In[66]:


mlsno2_tmp.keys()


# In[67]:


mlsno2_tmp['NO2_monthly_avg'].shape, mlsno2_tmp['LAT'].shape, mlsno2_tmp['LON'].shape


# In[68]:


ntime = len(fp_mls_no2)*12


# In[69]:


mls_no2 = {'lat':mlsno2_tmp['LAT'],'lon':mlsno2_tmp['LON'],'time':[],'no2':np.zeros((ntime,mlsno2_tmp['LAT'].shape[0],mlsno2_tmp['LON'].shape[1]))}


# In[70]:


mls_no2['time'] = []
for i,f in list(enumerate(fp_mls_no2)):
    mlsno2_tmp = sio.loadmat(fpp_mlsno2+f)
    yearstr = f.split('_')[4]
    [mls_no2['time'].append(datetime(int(yearstr[0:4]),m,15)) for m in range(1,13)]
    mls_no2['no2'][i*12:(i+1)*12,:,:] = mlsno2_tmp['NO2_monthly_avg']


# In[71]:


mls_no2_time = np.array(mls_no2['time'])


# In[72]:


mls_no2.keys()


# In[64]:


## Need to reform the mls_no2 lat and lon into 1d arrays


# In[73]:


def regrid(xin,yin,zin,xskip=1,yskip=1):
    from scipy.interpolate import griddata
    #xskip, yskip = 25,25
    xin = xin[::yskip,::xskip]
    yin = yin[::yskip,::xskip]
    zin = zin[::yskip,::xskip]
    # --------------------------
    # let us take some info from original coordinates:
    x0,x1,dx = np.min(xin),np.max(xin),np.abs(np.mean(np.diff(xin)))
    y0,y1,dy = np.min(yin),np.max(yin),np.abs(np.mean(np.diff(yin.T)))
    # --------------------------
    # let us make new (regular) coordinates:
    xout = np.arange(x0,x1+dx,dx)
    yout = np.arange(y0,y1+dx,dy)
    # --------------------------
    xm,ym = np.meshgrid(xout,yout)
    zo = np.griddata((xin.flatten(),yin.flatten()),zin.flatten(),(xm,ym),'nearest')
    return zo,xout,yout


# In[74]:


mls_no2['no2'].shape, mls_no2['lat'].shape, mls_no2['lon'].shape


# In[81]:


mls_no2_lat = mls_no2['lat'][:,0]


# In[82]:


mls_no2_lon = mls_no2['lon'][0,:]


# In[84]:


mls_no2_lon.shape


# ## Load CEDS NOx

# ### Aircraft

# In[9]:


ceds_1980,ceds_1980_dict = lu.load_netcdf(fp+'CEDS/CEDS_NOx_aircraft_monthly_1980_1999.nc',everything=True)


# In[10]:


ceds_2000,ceds_2000_dict = lu.load_netcdf(fp+'CEDS/CEDS_NOx_aircraft_monthly_2000_2019.nc',everything=True)


# In[11]:


ceds_2000_dict[b'Time']


# In[12]:


ceds_2000[b'CEDS_NOx_aircraft_emission'].shape


# In[13]:


ceds_2000[b'Lat'].shape, ceds_2000[b'Lon'].shape, ceds_2000[b'Time'].shape


# In[14]:


ceds_1980[b'Lat'].shape, ceds_1980[b'Lon'].shape, ceds_1980[b'Time'].shape


# In[15]:


ceds_AC_NOx = np.vstack((ceds_1980[b'CEDS_NOx_aircraft_emission'],ceds_2000[b'CEDS_NOx_aircraft_emission']))
ceds_time_days = np.hstack((ceds_1980[b'Time'],ceds_2000[b'Time']))


# In[16]:


ceds_AC_NOx.shape


# In[17]:


ceds_time = np.array([datetime(1750,1,1)+timedelta(days=int(d)) for d in ceds_time_days])


# In[18]:


ceds_lon = ceds_2000[b'Lon']
ceds_lat = ceds_2000[b'Lat']


# ### anthropogenic

# In[19]:


ceds_1980a,ceds_1980a_dict = lu.load_netcdf(fp+'CEDS/CEDS_NOx_anthro_monthly_1980_1999.nc',everything=True)
ceds_2000a,ceds_2000a_dict = lu.load_netcdf(fp+'CEDS/CEDS_NOx_anthro_monthly_2000_2019.nc',everything=True)


# In[20]:


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


# In[7]:


tcr_vals = ['_'.join([f.split('_')[1],f.split('_')[3]]) for f in fp_tcrl]


# In[8]:


tcr_tmp,tcr_tmp_dict = lu.load_netcdf(fp+'TCR/TCR_2_data/'+fp_tcrl[0],everything=True)


# In[9]:


tcr = {'lat':tcr_tmp[b'lat'],'lon':tcr_tmp[b'lon'],'time':tcr_tmp[b'time']}


# In[10]:


tcr_vals


# In[11]:


for i,f in list(enumerate(fp_tcrl)):
    tcr_tmp,tcr_tmp_dict = lu.load_netcdf(fp+'TCR/TCR_2_data/'+f,everything=True)
    tcr[tcr_vals[i]] = tcr_tmp[tcr_vals[i].split('_')[0].encode('utf-8')]


# In[12]:


tcr_tmp_dict[b'time']


# In[13]:


tcr_time = np.array([datetime(2005+int(m//12),1+int(m%12),15) for m in tcr['time']])


# In[14]:


# convert the 0 to 360 into a -180 to 180 longitude
tcr['lon'][tcr['lon']>=180.0] = tcr['lon'][tcr['lon']>=180.0] - 360.0


# In[15]:


plt.figure()
plt.hist(tcr['co_bio'].flatten(),bins=30)
plt.yscale('log')


# In[16]:


plt.figure()
plt.hist(tcr['co_anth'].flatten(),bins=30)
plt.yscale('log')


# ## Load GFED CO

# In[94]:


fp


# In[96]:


fp_gfedl = os.listdir(fp+'GFED/GFED_4.1/')
fp_gfeds = [f for f in fp_gfedl if '.hdf5' in f]
fp_gfeds.sort()


# In[97]:


fp_gfeds


# In[133]:


gfed_tmp,gfed_tmp_dict = lu.load_hdf(fp+'GFED/GFED_4.1/'+fp_gfeds[-1],all_values=True)


# In[118]:


nlon = len(gfed_tmp['//lon'][0,:])
nlat = len(gfed_tmp['//lat'][:,0])
nyears_gfed = len(fp_gfeds)


# In[139]:


gfed = {'bio_BB':np.zeros((nyears_gfed*12,nlat,nlon))+np.nan,
        'bio_NPP':np.zeros((nyears_gfed*12,nlat,nlon))+np.nan,
        'bio_Rh':np.zeros((nyears_gfed*12,nlat,nlon))+np.nan,
        'burn_fraction':np.zeros((nyears_gfed*12,nlat,nlon))+np.nan,
        'burn_source':np.zeros((nyears_gfed*12,nlat,nlon))+np.nan,
        'emi_C':np.zeros((nyears_gfed*12,nlat,nlon))+np.nan,
        'emi_DM':np.zeros((nyears_gfed*12,nlat,nlon))+np.nan,
        'lat':gfed_tmp['//lat'][:,0],
        'lon':gfed_tmp['//lon'][0,:]}


# In[140]:


gfed_tmp_vals = gfed_tmp.keys()


# In[141]:


gfed_vals = {'bio_BB':'//biosphere/{:02.0f}/BB',
             'bio_NPP':'//biosphere/{:02.0f}/NPP',
             'bio_Rh':'//biosphere/{:02.0f}/Rh',
             'burn_fraction':'//burned_area/{:02.0f}/burned_fraction',
             'burn_source':'//burned_area/{:02.0f}/source',
             'emi_C':'//emissions/{:02.0f}/C',
             'emi_DM':'//emissions/{:02.0f}/DM'}


# In[142]:


gfed['time'] = []


# In[143]:


for i,f in list(enumerate(fp_gfeds)):
    gfed_tmp,gfed_tmp_dict = lu.load_hdf(fp+'GFED/GFED_4.1/'+f,all_values=True,verbose=False)
    print('loaded file: {}'.format(f))
    for j in range(1,13):
        gfed['time'].append(datetime(int(f.split('.')[1].split('_')[1]),j,15))
        for k in gfed_vals:
            gfed[k][j-1+i*12,:,:] = gfed_tmp.get(gfed_vals[k].format(j),np.zeros((nlat,nlon)))


# In[144]:


gfed.keys()


# In[146]:


gfed['emi_C'].mean()


# ## Load MERRA2 - GMI

# ### Lightning NO2

# In[6]:


fp
fpp_gmi_lno = fp+'MERRA2/GMI_LNO_processed/'


# In[8]:


fp_gmi_lno = os.listdir(fpp_gmi_lno)
fp_gmi_lno = [f for f in fp_gmi_lno if '.nc' in f]
fp_gmi_lno.sort()


# In[9]:


fp_gmi_lno


# In[11]:


gmi_lno_tmp,gmi_lno_dict_tmp = lu.load_netcdf(fpp_gmi_lno+fp_gmi_lno[0],everything=True)


# In[43]:


gmi_lno_dict_tmp


# In[12]:


gmi_lno_tmp.keys()


# In[14]:


gmi_lno_tmp[b'M2GMI_NO_lightning_emission'].shape, gmi_lno_tmp[b'Lat'].shape, gmi_lno_tmp[b'Lon'].shape


# In[15]:


ntime = len(fp_gmi_lno)*12


# In[17]:


gmi_lno = {'lat':gmi_lno_tmp[b'Lat'],
           'lon':gmi_lno_tmp[b'Lon'],
           'time':[],
           'no_lightning_emission':np.zeros((ntime,gmi_lno_tmp[b'Lat'].shape[0],gmi_lno_tmp[b'Lon'].shape[0]))}


# In[18]:


gmi_lno['time'] = []
for i,f in list(enumerate(fp_gmi_lno)):
    gmi_lno_tmp,gmi_lno_dict_tmp = lu.load_netcdf(fpp_gmi_lno+f,everything=True)
    yearstr = f.split('_')[-1]
    [gmi_lno['time'].append(datetime(int(yearstr[0:4]),m,15)) for m in range(1,13)]
    gmi_lno['no_lightning_emission'][i*12:(i+1)*12,:,:] = gmi_lno_tmp[b'M2GMI_NO_lightning_emission']


# In[19]:


gmi_lno_time = np.array(gmi_lno['time'])


# In[20]:


gmi_lno.keys()


# ### GMI O3 production

# In[27]:


fpp_gmi_o3prod = fp+'MERRA2/GMI_O3_prod_processed/'


# In[29]:


fp_gmi_o3prod = os.listdir(fpp_gmi_o3prod)
fp_gmi_o3prod = [f for f in fp_gmi_o3prod if '.nc' in f]
fp_gmi_o3prod.sort()


# In[30]:


fp_gmi_o3prod


# In[32]:


gmi_o3prod_tmp,gmi_o3prod_dict_tmp = lu.load_netcdf(fpp_gmi_o3prod+fp_gmi_o3prod[0],everything=True)


# In[42]:


gmi_o3prod_dict_tmp


# In[33]:


gmi_o3prod_tmp[b'Time'].shape,gmi_o3prod_tmp[b'M2GMI_O3_net_prod_PBL'].shape


# In[34]:


gmi = gmi_lno


# In[35]:


gmi.keys()


# In[36]:


gmi['O3prod_PBL'] = np.zeros((ntime,gmi_o3prod_tmp[b'Lat'].shape[0],gmi_o3prod_tmp[b'Lon'].shape[0]))
gmi['O3prod_TROP'] = np.zeros((ntime,gmi_o3prod_tmp[b'Lat'].shape[0],gmi_o3prod_tmp[b'Lon'].shape[0]))


# In[37]:


for i,f in list(enumerate(fp_gmi_o3prod)):
    gmi_tmp,gmi_dict_tmp = lu.load_netcdf(fpp_gmi_o3prod+f,everything=True)
    gmi['O3prod_PBL'][i*12:(i+1)*12,:,:] = gmi_tmp[b'M2GMI_O3_net_prod_PBL']
    gmi['O3prod_TROP'][i*12:(i+1)*12,:,:] = gmi_tmp[b'M2GMI_O3_net_prod_TROP']


# ### GMI Stratospheric O3

# In[38]:


fpp_gmi_o3strat = fp+'MERRA2/GMI_StratO3_processed/'


# In[39]:


fp_gmi_o3strat = os.listdir(fpp_gmi_o3strat)
fp_gmi_o3strat = [f for f in fp_gmi_o3strat if '.nc' in f]
fp_gmi_o3strat.sort()
fp_gmi_o3strat


# In[40]:


gmi_o3strat_tmp,gmi_o3strat_dict_tmp = lu.load_netcdf(fpp_gmi_o3strat+fp_gmi_o3strat[0],everything=True)


# In[41]:


gmi_o3strat_dict_tmp


# In[44]:


gmi_o3strat_tmp[b'Time'].shape,gmi_o3strat_tmp[b'M2GMI_Strat_O3_conc_PBL'].shape


# In[45]:


gmi['O3strat_PBL'] = np.zeros((ntime,gmi_o3strat_tmp[b'Lat'].shape[0],gmi_o3strat_tmp[b'Lon'].shape[0]))
gmi['O3strat_TROP'] = np.zeros((ntime,gmi_o3strat_tmp[b'Lat'].shape[0],gmi_o3strat_tmp[b'Lon'].shape[0]))


# In[46]:


for i,f in list(enumerate(fp_gmi_o3strat)):
    gmi_tmp,gmi_dict_tmp = lu.load_netcdf(fpp_gmi_o3strat+f,everything=True)
    gmi['O3strat_PBL'][i*12:(i+1)*12,:,:] = gmi_tmp[b'M2GMI_Strat_O3_conc_PBL']
    gmi['O3strat_TROP'][i*12:(i+1)*12,:,:] = gmi_tmp[b'M2GMI_Strat_O3_conc_TROP']


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

# ## Definition of the regions and functions

# In[22]:


rgs = {'Southeast Asia':[[-12,95],[18,140]],
       'China':[[21,75],[46,127]],
       'South Asia':[[6,67],[31,91]],
       'Siberia':[[50,70],[75,180]],
       'Tropical Pacific Ocean':[[5,-180],[35,-130]]     
      } #lower left [lat lon], upper right [lat lon]


# In[23]:


rgs['China'][0]


# In[24]:


def multi_stats_pd(data,time,name='trop_NO2',axis=1):
    'to get a dict for the different dataframes (mean, median, min, max)'
    pds = {}
    pds['mean'] = build_pd(np.nanmean(np.nanmean(data,axis=axis),axis=axis),time,name='mean_'+name)
    pds['median'] = build_pd(np.nanmedian(np.nanmedian(data,axis=axis),axis=axis),time,name='median_'+name)
    pds['min'] = build_pd(np.nanmin(np.nanmin(data,axis=axis),axis=axis),time,name='min_'+name)
    pds['max'] = build_pd(np.nanmax(np.nanmax(data,axis=axis),axis=axis),time,name='max_'+name)
    pds['std'] = build_pd(np.nanstd(np.nanstd(data,axis=axis),axis=axis),time,name='std_'+name)
    
    return pds
    


# In[25]:


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
    


# ## Subset for GMI

# In[49]:


gmi.keys()


# In[57]:


gmi['no_lightning_emission'].shape


# In[66]:


gmi_rg = {}
first = True
for rg in rgs: 
    print(rg)
    
    ill_lon = np.argmin(np.abs(gmi['lon']-rgs[rg][0][1]))
    ill_lat = np.argmin(np.abs(gmi['lat']-rgs[rg][0][0]))
    iur_lon = np.argmin(np.abs(gmi['lon']-rgs[rg][1][1]))
    iur_lat = np.argmin(np.abs(gmi['lat']-rgs[rg][1][0]))
    
    if iur_lat<ill_lat: 
        iur_lat,ill_lat = ill_lat,iur_lat
    
    for val in ['no_lightning_emission', 'O3prod_PBL', 'O3prod_TROP', 'O3strat_PBL', 'O3strat_TROP']:
        if first: gmi_rg[val] = {}
        gmi_rg[val][rg] = multi_stats_pd(gmi[val][:,ill_lat:iur_lat,ill_lon:iur_lon],gmi_lno_time,axis=1,name='GMI'+val)
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


# ## GMI

# In[77]:


gmi_rg[val][rg][typ].keys()


# In[ ]:





# In[80]:


for val in ['no_lightning_emission', 'O3prod_PBL', 'O3prod_TROP', 'O3strat_PBL', 'O3strat_TROP']:
    for rg in gmi_rg[val]:
        for typ in ['mean', 'median', 'std']:
            result_gmi_nol = seasonal_decompose(gmi_rg[val][rg][typ][typ+'_GMI'+val])
            p = result_gmi_nol.plot()
            p.get_axes()[0].set_title(rg+ ': '+typ+'_GMI_'+val)
            p.get_axes()[0].set_ylabel('All')
            
            p.savefig(os.path.join(fp,'MERRA2','GMI_Rgs_deseason_'+val+'_'+rg+'_'+typ+'.png'),transparent=True,dpi=600)


# In[72]:


gmi_rg.keys()


# # Make figures of the trends for each points

# ## Load the functions

# In[83]:


def run_the_seasonal_decomp(data,time,name,nlon,nlat,i):
    trend_tmp = np.empty((nlon))
    trend_pval_tmp = np.empty((nlon))
    for j in range(nlon):
        # Extract the time series for this lat and lon point
        ts = build_pd(data[:, i, j],time,name=name)
        ts.dropna()

        # Perform seasonal decomposition on the time series
        result = seasonal_decompose(ts[name],extrapolate_trend='freq') #, model='additive', period=12)

        # Extract the trend component from the decomposition
        trend_comp = result.trend
        trend_comp.dropna()

        # Calculate the straight trend using the OLS method
        fitl = sm.OLS(trend_comp, sm.add_constant(range(len(trend_comp))),missing='drop').fit()
        trend_tmp[j] = fitl.params[1]

        # Calculate the p-value for the trend using a t-test
        trend_pval_tmp[j] = fitl.pvalues[1]
    return trend_pval_tmp,trend_tmp


# In[84]:


def trends_and_pval_per_lat_lon(data,lat,lon,time,name='mean_trop_NO2',parallel=False):
    "Find the linear trend and its pval for each lat/lon point"
    nlat = len(lat)
    nlon = len(lon)
    
    # Create an empty array to hold the trend and p-value for each lat and lon point
    trend = np.empty((nlat, nlon))
    trend_pval = np.empty((nlat, nlon))
    
    

     
    if not parallel:
        rtsd = partial(run_the_seasonal_decomp,data,time,name,nlon,nlat)
        for i in tqdm(range(nlat)):
            trend[i,:],trend_pval[i,:] = rtsd(i)
    else:
        from pqdm.processes import pqdm
        rtsd = partial(run_the_seasonal_decomp,data,time,name,nlon,nlat)
        trend = pqdm(range(nlat), rtsd, n_jobs=16)


    return trend,trend_pval



# In[85]:


def trends_and_pval_per_lat_lon_single(data,lat,lon,time,name='mean_trop_NO2'):
    "Find the linear trend and its pval for each lat/lon point"
    nlat = len(lat)
    nlon = len(lon)
    
    if np.ma.isMaskedArray(data): 
        ma = True
    else:
        ma = False
    
    # Create an empty array to hold the trend and p-value for each lat and lon point
    trend = np.empty((nlat, nlon))
    trend_pval = np.empty((nlat, nlon))
    trend_rmse = np.empty((nlat, nlon))
    
    for i in tqdm(range(nlat)):
        for j in range(nlon):
            # Extract the time series for this lat and lon point
            if ma:
                if data[:,i,j].mask.all():
                    #print('...All NaN for i={}, j={}'.format(i,j))
                    trend[i, j] = np.nan
                    trend_pval[i, j] = np.nan
                    continue
            else:
                if sum(np.isfinite(data[:,i,j]))<2:
                    trend[i, j] = np.nan
                    trend_pval[i, j] = np.nan
                    continue
            ts = build_pd(data[:, i, j],time,name=name)
            ts.dropna(inplace=True)
            
            if len(ts[name])<24:
                #print('...Not enough points for i={}, j={}'.format(i,j))
                trend[i, j] = np.nan
                trend_pval[i, j] = np.nan
                continue

            # Perform seasonal decomposition on the time series
            try:
                result = seasonal_decompose(ts[name],extrapolate_trend='freq') #, model='additive', period=12)
            except Exception:
                print(traceback.format_exc())
                set_trace()

            # Extract the trend component from the decomposition
            trend_comp = result.trend
            trend_comp.dropna(inplace=True)

            # Calculate the straight trend using the OLS method
            fitl = sm.OLS(trend_comp, sm.add_constant(range(len(trend_comp))),missing='drop').fit()
            trend[i, j] = fitl.params[1]
            trend_rmse[i,j] = np.sqrt(fitl.mse_total)

            # Calculate the p-value for the trend using a t-test
            trend_pval[i, j] = fitl.pvalues[1]

    return trend,trend_pval,trend_rmse


# In[86]:


def convert_monthly_to_decadal_trend(trend,time):
    'convert the trends in per month to per decade'
    def diff_month(d1, d2):
        return (d1.year - d2.year) * 12 + d1.month - d2.month
    
    dt_months = diff_month(time[0],time[-1])
    dt_decades = dt_months/12.0/10.0
    
    return trend*dt_months/dt_decades
    


# In[87]:


def plot_trend_and_pval(trend,trend_pval,lon,lat,name='',cax_name='Trend',figsize=(10,4),
                        clevels=np.arange(-0.025, 0.04, 0.005),vmin=-0.025,vmax=0.025,ncolors=12,npval_skip=5,msize=0.5):
    
    pr = ccrs.PlateCarree()
    # Set up the plot
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=pr)

    # Define the colormap
    cmap = plt.cm.get_cmap('RdBu_r', ncolors)

    # Create a filled contour plot of the trend
    contour_levels = clevels
    
    #contour = ax.contourf(lon, lat, trend, contour_levels,
    #                      cmap=cmap, transform=pr, extend='both')
    contour = ax.pcolormesh(lon, lat, trend, transform=pr,cmap=cmap,shading='auto',vmin=vmin,vmax=vmax)

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
    try:
        ax.scatter(lon[significant[1]][0::npval_skip], lat[significant[0]][0::npval_skip], marker='+', color='k',
               transform=pr,s=msize)
    except:
        print('** Error plotting the pvalues')

    # Add a title
    plt.title('Trend for '+name, fontsize=18)
    ax.add_feature(cartopy.feature.COASTLINE)
    return fig


# In[88]:


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


# In[89]:


def trend_and_plot(data,lat,lon,time,name,fp=fp,clevels=np.arange(-0.01,0.01,0.0025),vv=vv,units='DU',npvals=10):
    'combine the trends and ploting'
    trend,trend_pval,trend_rmse = trends_and_pval_per_lat_lon_single(data,lat,lon,time,name=name)
    #trend,trend_pval,trend_rmse = trend_for_multi_lon(data,lat,lon,time,name=name)
    print('saving to the pickle: '+fp+name+'_trend_output.{}.p'.format(vv))
    pickle.dump({name+'_trend':trend,name+'_trend_pval':trend_pval,name+'_trend_rmse':trend_rmse,
             name+'_time':time,name:data},open(fp+name+'_trend_output.{}.p'.format(vv),"wb"))
    
    vmin = np.percentile(convert_monthly_to_decadal_trend(trend,time),0.2)
    
    
    fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(trend,time),trend_pval,lon,lat,
                          name=name+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(time[0],time[-1]),
                          cax_name=name+'\n trend ['+units+'/decade]',figsize=(10,4),
                          clevels=clevels,vmin=vmin,vmax=np.abs(vmin),npval_skip=npvals)
    fig.tight_layout()
    fig.savefig(fp+name+'_trend_output.{}.png'.format(vv),dpi=600,transparent=True)
    print('figure saved to :'+fp+name+'_trend_output.{}.png'.format(vv))
    


# In[90]:


def plot_avg_and_rmse(avg,rmse,trend,lon,lat,name='',cax_name='Average',figsize=(10,4),
                        clevels=np.arange(-0.025, 0.04, 0.005),vmin=None,vmax=None,ncolors=12,nlevels=12):
    
    pr = ccrs.PlateCarree()
    # Set up the plot
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=pr)

    # Define the colormap
    cmap = plt.cm.get_cmap('plasma', ncolors)
    
    if np.ma.is_masked(avg):
        if not vmin: vmin = np.percentile(avg.compressed(),0.5)
        if not vmax: vmax = np.percentile(avg.compressed(),99.5)
    else:
        if not vmin: vmin = np.percentile(avg,0.5)
        if not vmax: vmax = np.percentile(avg,99.5)
    

    # Create a filled contour plot of the trend
    contour_levels = clevels
    
    #contour = ax.contourf(lon, lat, trend, contour_levels,
    #                      cmap=cmap, transform=pr, extend='both')
    contour = ax.pcolormesh(lon, lat, avg, transform=pr,cmap=cmap,shading='auto',vmin=vmin,vmax=vmax)

    # Add a colorbar
    cbar = plt.colorbar(contour, shrink=0.6, pad=0.02)
    cbar.ax.set_ylabel(cax_name, fontsize=16)

    # Add a land mask
    ax.add_feature(cartopy.feature.LAND, edgecolor='black', facecolor='None', alpha=0.3,zorder=-100)

    # Add a grid
    gl = ax.gridlines(color='gray', linestyle='--',draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False

    # Add contour lines where the rmse is greater than the trend
    isub = trend>rmse
    rms = np.ma.masked_array(rmse,mask=isub)
    cl = ax.contourf(lon, lat, rms, 4, hatches=['//','.','..','oo'], colors='none')
    #cl = ax.contour(lon[isub], lat[isub], rmse[isub], nlevels, color='lightgrey',linewidths=0.5, transform=pr)
    #ax.clabel(cl, inline=True, fontsize=10)
    

    # Add a title
    plt.title('Averages for '+name, fontsize=18)
    ax.add_feature(cartopy.feature.COASTLINE)
    return fig


# ### Make multiprocessing

# In[91]:


from multiprocessing import Pool, cpu_count
import signal


# In[92]:


class KeyboardInterruptError(Exception): pass


# In[93]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print('signal: %s' % signal_num)
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[94]:


def single_trend(data,time,name,ma):

    if ma:
        if data.mask.all():
            return np.nan,np.nan,np.nan
    else:
        if sum(np.isfinite(data))<2:
            return np.nan,np.nan,np.nan

    ts = build_pd(data,time,name=name)
    ts.dropna(inplace=True)

    if len(ts[name])<24:
        return np.nan,np.nan,np.nan

    # Perform seasonal decomposition on the time series
    result = seasonal_decompose(ts[name],extrapolate_trend='freq') #, model='additive', period=12)

    # Extract the trend component from the decomposition
    trend_comp = result.trend
    trend_comp.dropna(inplace=True)

    # Calculate the straight trend using the OLS method
    fitl = sm.OLS(trend_comp, sm.add_constant(range(len(trend_comp))),missing='drop').fit()
    trend_tmp = fitl.params[1]
    trend_tmp_rmse = np.sqrt(fitl.mse_total)

    # Calculate the p-value for the trend using a t-test
    trend_tmp_pval = fitl.pvalues[1]

    return trend_tmp,trend_tmp_pval,trend_tmp_rmse


# In[95]:


def trend_for_multi_lon(data,lat,lon,time,name='mean_trop_NO2',subn=None):
    "Find the linear trend and its pval for each lat/lon point"
    nlat = len(lat)
    nlon = len(lon)
    
    if np.ma.isMaskedArray(data): 
        ma = True
    else:
        ma = False
    
    # Create an empty array to hold the trend and p-value for each lat and lon point
    trend = np.empty((nlat, nlon))
    trend_pval = np.empty((nlat, nlon))
    trend_rmse = np.empty((nlat, nlon))
    
    ntotal = nlat*nlon
    if not subn:
        subn=ntotal
    with Pool(cpu_count()-2) as p:#,worker_init)
        time_t = list(map(list, zip(*[time for n in range(ntotal)])))
        inputs = zip(data.reshape(len(time),-1),time_t,[name]*ntotal,[ma]*subn)
        results = p.starmap(single_trend, tqdm(inputs, total=ntotal))
    for k,outs in enumerate(results):
        i,j = np.unravel_index(k,data[0,:,:].shape)
        trend[i,j] = outs[0]
        trend_pval[i,j] = outs[1]
        trend_rmse[i,j] = outs[2]
            
    return trend,trend_pval,trend_rmse


# ## Plot the GOME NO2

# In[116]:


data.reshape(len(time),-1).shape


# In[140]:


trend_and_plot(gome[b'TroposNO2'],gome[b'lat'],gome[b'lon'],gometime,'GOME_NO2_tropos',
               fp=fp+'GOME_SCIAMACHY_GOME2_NO2_L3/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# In[16]:


gome_dict = pd.read_pickle('/data2/ACCDAM_low_ozone/GOME_SCIAMACHY_GOME2_NO2_L3/GOME_NO2_tropos_trend_output.v1.p')


# In[17]:


gome_dict.keys()


# In[36]:


plt.rcParams['text.usetex'] = False


# In[37]:


fig = plot_trend_and_pval(gome_dict['GOME_NO2_tropos_trend'],gome_dict['GOME_NO2_tropos_trend_pval'],
                          gome[b'lon'],gome[b'lat'],
                          name='GOME_NO2_tropos'+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(gome_dict['GOME_NO2_tropos_time'][0],gome_dict['GOME_NO2_tropos_time'][-1]),
                          cax_name='GOME_NO2_tropos\n'+'trend [#/cm^2 /decade]',figsize=(10,4),
                          vmin=-0.08,vmax=0.08,ncolors=25,npval_skip=200,msize=0.3)
fig.tight_layout()
fig.savefig(fp+'GOME_NO2_tropos'+'_trend_output.{}.png'.format(vv),dpi=600,transparent=True)


# In[31]:


fig.show()


# In[46]:


isub = gome_dict['GOME_NO2_tropos_trend_rmse']>gome_dict['GOME_NO2_tropos_trend']


# In[80]:


figa = plot_avg_and_rmse(np.nanmean(gome_dict['GOME_NO2_tropos'],axis=0),gome_dict['GOME_NO2_tropos_trend_rmse'],
                         gome_dict['GOME_NO2_tropos_trend'],gome[b'lon'],gome[b'lat'],
                         name='GOME_NO2_tropos'+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(gome_dict['GOME_NO2_tropos_time'][0],gome_dict['GOME_NO2_tropos_time'][-1]),
                         cax_name='GOME_NO2_tropos\n'+'averages [#/cm^2]',figsize=(10,4),ncolors=25)
figa.tight_layout()
figa.savefig(fp+'GOME_NO2_tropos'+'_averages_output.{}.png'.format(vv),dpi=600,transparent=True)


# In[75]:


gome[b'lon'].shape,gome[b'lat'].shape


# In[62]:


isub[:,0].shape


# ## Plot the MOPITT

# ### MOPITT CO day

# Units: ppbv, name: Retrieved CO Mixing Ratio Profile Day

# In[119]:


mop_pre


# In[120]:


mop_COday.shape, len(mop_lat),len(mop_lon),len(mop_time)


# In[121]:


ip = 0


# In[122]:


mop_COday[:,ip,:,:].shape


# In[123]:


np.nanmean(mop_COday[:,ip,:,:],axis=0).shape


# In[124]:


for ip in range(9):
    plt.figure()
    ctf = plt.contourf(mop_lon,mop_lat,np.nanmean(mop_COday[:,:,:,ip],axis=0).T,cmap=plt.cm.get_cmap('plasma', 25),shading='nearest')

    plt.title('MOPITT CO day Averaged {} mb'.format(mop_pre[ip]))
    plt.ylabel('Latitude [deg]')
    plt.xlabel('Longitude [deg]')
    plt.colorbar(ctf,label='Retrieved CO Mixing Ratio Profile Day [ppbv]')

    plt.savefig(fp+'MOPITT/'+'MOPITT_CO_day_Averaged_{}_v2.png'.format(mop_pre[ip]),dpi=600,transparent=True)


# In[125]:


np.swapaxes(mop_COday[:,:,:,ip],2,1).shape


# In[126]:


len(mop_lat)


# In[127]:


for ip,pre in enumerate(mop_pre):
    trend_and_plot(np.swapaxes(mop_COday[:,:,:,ip],2,1),mop_lat,mop_lon,mop_time,'MOPITT_CO_day_{:3.0f}mb'.format(pre),
                   fp=fp+'MOPITT/',vv=vv,npvals=5,units='ppbv')


# ### Redo figures

# In[103]:


mop_day_900 = pd.read_pickle(fp+'MOPITT/MOPITT_CO_day_900mb_trend_output.v1.p')
mop_day_900.keys()



# In[104]:


mop_lat,mop_lon


# ### MOPITT CO night

# In[250]:


for ip in range(9):
    plt.figure()
    ctf = plt.contourf(mop_lon,mop_lat,np.nanmean(mop_COnight[:,:,:,ip],axis=0).T,cmap=plt.cm.get_cmap('plasma', 25),shading='nearest')

    plt.title('MOPITT CO night Averaged {} mb'.format(mop_pre[ip]))
    plt.ylabel('Latitude [deg]')
    plt.xlabel('Longitude [deg]')
    plt.colorbar(ctf,label='Retrieved CO Mixing Ratio Profile Night [ppbv]')

    plt.savefig(fp+'MOPITT/'+'MOPITT_CO_night_Averaged_{}_v2.png'.format(mop_pre[ip]),dpi=600,transparent=True)


# In[249]:


for ip,pre in enumerate(mop_pre):
    trend_and_plot(np.swapaxes(mop_COnight[:,:,:,ip],2,1),mop_lat,mop_lon,mop_time,'MOPITT_CO_night_{:3.0f}mb'.format(pre),
                   fp=fp+'MOPITT/',vv=vv,npvals=5,units='ppbv')


# In[24]:


mop_night_900 = pd.read_pickle('/data2/ACCDAM_low_ozone/MOPITT/MOPITT_CO_night_900mb_trend_output.v1.p')
mop_night_900.keys()


# In[ ]:


mop_night_900['MOPITT_CO_night_900mb_trend']


# ## Plot the TOMS O3

# In[43]:


np.moveaxis(toms_O3,2,0).shape


# In[35]:


toms_time.shape,toms_lat.shape,toms_lon.shape


# In[44]:


trend_and_plot(np.moveaxis(toms_O3,2,0),toms_lat,toms_lon,toms_time,'TOMS_O3',
               fp=fp+'TOMS_O3_L3/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# ##  Plot the OMI MLS

# ### O3

# In[35]:


mls_o3_trend,mls_o3_trend_pval,mls_o3_trend_rmse = trends_and_pval_per_lat_lon_single(mls_o3['tropo3'],
                                    mls_o3['lat'],mls_o3['lon'],mls_o3_time,name='OMI_MLS_O3')


# In[36]:


mls_o3_time[0],mls_o3_time[-1]


# In[39]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(mls_o3_trend,mls_o3_time).flatten(),bins=50)


# In[40]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(mls_o3_trend,mls_o3_time),
                          mls_o3_trend_pval,mls_o3['lon'],mls_o3['lat'],
                          name='OMI_MLS_O3'+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(mls_o3_time[0],mls_o3_time[-1]),
                          cax_name='O3 trend [DU/decade]',figsize=(10,4),vmin=-4,vmax=4)
fig.tight_layout()
fig.savefig(fp+'OMI_MLS_O3_L3/OMI_MLS_O3_trend_v2.png')


# ### NO2

# In[60]:


mls_no2.keys()


# In[61]:


mls_no2['no2'].shape


# In[63]:


len(mls_no2['time']),mls_no2['lat'].shape,mls_no2['lon'].shape


# In[85]:


trend_and_plot(mls_no2['no2'],mls_no2_lat,mls_no2_lon,mls_no2_time,'MLS_NO2',
               fp=fp+'OMI_NO2/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# In[86]:


mls_no2_dict = pd.read_pickle(fp+'OMI_NO2/MLS_NO2_trend_output.v1.p')


# In[88]:


mls_no2_dict.keys()


# In[89]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(mls_no2_dict['MLS_NO2_trend'],mls_no2_dict['MLS_NO2_time']).flatten(),
         bins=50)


# In[ ]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(mls_no2_dict['MLS_NO2_trend'],mls_no2_dict['MLS_NO2_time']),
                          mls_no2_dict['MLS_NO2_trend_pval'],mls_no2_lon,mls_no2_lat,
                          name='MLS_NO2'+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(mls_no2_dict['MLS_NO2_time'][0],mls_no2_dict['MLS_NO2_time'][-1]),
                          cax_name='MLS_NO2 trend [DU/decade]',figsize=(10,4),
                          vmin=-1e15,vmax=1e15,ncolors=50,npval_skip=200,msize=0.3)
fig.tight_layout()
fig.savefig(fp+'OMI_NO2/MLS_NO2_trend_output.v1.png',dpi=300,transparent=True)


# ## Plot CEDS NOx

# ### NOx Aircraft emissions

# In[39]:


ceds_AC_NOx.shape, ceds_time.shape


# In[49]:


ceds_lat.shape, ceds_lon.shape


# In[40]:


trend_and_plot(ceds_AC_NOx,ceds_lat,ceds_lon,ceds_time,'CED_NOx_AC',
               fp=fp+'CEDS/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# In[46]:


ceds_ac_dict = pd.read_pickle(fp+'CEDS/CED_NOx_AC_trend_output.v1.p')


# In[47]:


ceds_ac_dict.keys()


# In[48]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(ceds_ac_dict['CED_NOx_AC_trend'],ceds_ac_dict['CED_NOx_AC_time']).flatten(),
         bins=50)


# In[51]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(ceds_ac_dict['CED_NOx_AC_trend'],ceds_ac_dict['CED_NOx_AC_time']),
                          ceds_ac_dict['CED_NOx_AC_trend_pval'],ceds_lon,ceds_lat,
                          name='CED_NOx_AC'+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(ceds_ac_dict['CED_NOx_AC_time'][0],ceds_ac_dict['CED_NOx_AC_time'][-1]),
                          cax_name='CED_NOx_AC trend [DU/decade]',figsize=(10,4),
                          vmin=-0.2e-11,vmax=0.2e-11,ncolors=50,npval_skip=50,msize=0.3)
fig.tight_layout()
fig.savefig(fp+'CEDS/CED_NOx_AC_trend_output.v1.png',dpi=300,transparent=True)


# ### NOx Anthropogenic

# In[41]:


ceds_Ant_NOx.shape


# In[42]:


ceds_time.shape


# In[44]:


trend_and_plot(ceds_Ant_NOx,ceds_lat,ceds_lon,ceds_time,'CED_NOx_Ant',
               fp=fp+'CEDS/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# In[52]:


ceds_ant_dict = pd.read_pickle(fp+'CEDS/CED_NOx_Ant_trend_output.v1.p')


# In[53]:


ceds_ant_dict.keys()


# In[55]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(ceds_ant_dict['CED_NOx_Ant_trend'],ceds_ant_dict['CED_NOx_Ant_time']).flatten(),
         bins=50)


# In[58]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(ceds_ant_dict['CED_NOx_Ant_trend'],ceds_ant_dict['CED_NOx_Ant_time']),
                          ceds_ant_dict['CED_NOx_Ant_trend_pval'],ceds_lon,ceds_lat,
                          name='CED_NOx_Ant'+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(ceds_ant_dict['CED_NOx_Ant_time'][0],ceds_ant_dict['CED_NOx_Ant_time'][-1]),
                          cax_name='CED_NOx_Ant trend [DU/decade]',figsize=(10,4),
                          vmin=-0.15e-9,vmax=0.15e-9,ncolors=50,npval_skip=50,msize=0.3)
fig.tight_layout()
fig.savefig(fp+'CEDS/CED_NOx_Ant_trend_output.v1.png',dpi=300,transparent=True)


# ## Plot TCR

# ### TCR Nox Anth

# In[45]:


tcr.keys()


# In[26]:


tcr_noxanth_trend,tcr_noxanth_trend_pval = trends_and_pval_per_lat_lon(tcr['nox_anth'],tcr['lat'],tcr['lon'],tcr_time,name='TCR_NOx_Anth',parallel=False)


# In[35]:


tcr_noxanth_trend,tcr_noxanth_trend_pval = trends_and_pval_per_lat_lon_single(tcr['nox_anth'],tcr['lat'],tcr['lon'],tcr_time,name='TCR_NOx_Anth')


# #### Save the TCR values before plot

# In[42]:


vv_tcr = 'v2'


# In[37]:


pickle.dump({'tcr_noxanth_trend':tcr_noxanth_trend,'tcr_noxanth_trend_pval':tcr_noxanth_trend_pval,
             'tcr_time':tcr_time,'tcr':tcr},open(fp+'TCR/'+'TCR_EMI_NOx_anth_trend_output.{}.p'.format(vv_tcr),"wb"))


# In[38]:


tcr_time[0],tcr_time[-1]


# #### Load from the pickle

# In[43]:


tcr_dict = pd.read_pickle(fp+'TCR/'+'TCR_EMI_NOx_anth_trend_output.{}.p'.format(vv_tcr))


# In[44]:


for k in tcr_dict:
    locals()[k] = tcr_dict[k]


# #### Make plot

# In[47]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(tcr_noxanth_trend,tcr_time).flatten(),bins=50)


# In[50]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(tcr_noxanth_trend,tcr_time),
                          tcr_noxanth_trend_pval,tcr['lon'],tcr['lat'],
                          name='TCR_NOx_Anth'+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(tcr_time[0],tcr_time[-1]),
                          cax_name='TCR_NOx_Anth trend [DU/decade]',figsize=(10,4),
                          clevels=np.arange(-0.000001, 0.01, 0.0000001),vmin=-0.5e-11,vmax=0.5e-11,ncolors=50)
fig.tight_layout()
fig.savefig(fp+'TCR/TCR_EMI_NOx_anth_trend.png',dpi=300,transparent=True)


# ### TCR CO Surface

# In[34]:


tcr.keys()


# In[35]:


tcr_co_sfc_trend,tcr_co_sfc_trend_pval = trends_and_pval_per_lat_lon_single(tcr['co_sfc'],
                                                                              tcr['lat'],tcr['lon'],
                                                                              tcr_time,name='TCR_CO_SFC')


# In[36]:


pickle.dump({'tcr_co_sfc_trend':tcr_co_sfc_trend,'tcr_co_sfc_trend_pval':tcr_co_sfc_trend_pval,
             'tcr_time':tcr_time,'tcr':tcr},open(fp+'TCR/'+'TCR_EMI_CO_SFC_trend_output.{}.p'.format(vv_tcr),"wb"))


# In[51]:


tcr_dict_cosfc = pd.read_pickle(fp+'TCR/'+'TCR_EMI_CO_SFC_trend_output.{}.p'.format(vv_tcr))
for k in tcr_dict_cosfc:
    locals()[k] = tcr_dict_cosfc[k]


# In[53]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(tcr_co_sfc_trend,tcr_time).flatten(),bins=50)


# In[55]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(tcr_co_sfc_trend,tcr_time),
                          tcr_co_sfc_trend_pval,tcr['lon'],tcr['lat'],
                          name='TCR_CO_SFC'+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(tcr_time[0],tcr_time[-1]),
                          cax_name='TCR_CO_SFC trend [DU/decade]',figsize=(10,4),
                          clevels=np.arange(-0.01, 0.01, 0.0025),vmin=-0.5e-9,vmax=0.5e-9)
fig.tight_layout()
fig.savefig(fp+'TCR/'+'TCR_EMI_CO_SFC_trend_output.{}.png'.format(vv_tcr),dpi=600,transparent=True)


# ### TCR co_bio

# In[44]:


nm = 'TCR_CO_BIO'
tcr_co_bio_trend,tcr_co_bio_trend_pval = trends_and_pval_per_lat_lon_single(tcr['co_bio'],
                                                                              tcr['lat'],tcr['lon'],
                                                                              tcr_time,name=nm)


# In[45]:


pickle.dump({'tcr_co_bio_trend':tcr_co_bio_trend,'tcr_co_bio_trend_pval':tcr_co_bio_trend_pval,
             'tcr_time':tcr_time,'tcr':tcr},open(fp+'TCR/'+nm+'_trend_output.{}.p'.format(vv_tcr),"wb"))


# In[56]:


tcr_dict_cosfc = pd.read_pickle(fp+'TCR/'+'TCR_CO_BIO_trend_output.{}.p'.format(vv_tcr))
for k in tcr_dict_cosfc:
    locals()[k] = tcr_dict_cosfc[k]


# In[57]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(tcr_co_bio_trend,tcr_time).flatten(),bins=50)


# In[59]:


nm = 'TCR_CO_BIO'
fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(tcr_co_bio_trend,tcr_time),
                          tcr_co_bio_trend_pval,tcr['lon'],tcr['lat'],
                          name=nm+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(tcr_time[0],tcr_time[-1]),
                          cax_name=nm+' trend [DU/decade]',figsize=(10,4),
                          clevels=np.arange(-0.01, 0.01, 0.0025),vmin=-0.2e-9,vmax=0.2e-9)
fig.tight_layout()
fig.savefig(fp+'TCR/'+nm+'_trend_output.{}.png'.format(vv_tcr),dpi=600,transparent=True)


# ### TCR co anth

# In[60]:


nm = 'TCR_CO_ANTH'


# In[48]:


tcr_co_anth_trend,tcr_co_anth_trend_pval = trends_and_pval_per_lat_lon_single(tcr['co_anth'],
                                                                              tcr['lat'],tcr['lon'],
                                                                              tcr_time,name=nm)


# In[49]:


pickle.dump({'tcr_co_anth_trend':tcr_co_anth_trend,'tcr_co_anth_trend_pval':tcr_co_anth_trend_pval,
             'tcr_time':tcr_time,'tcr':tcr},open(fp+'TCR/'+nm+'_trend_output.{}.p'.format(vv_tcr),"wb"))


# In[62]:


tcr_dict_coanth = pd.read_pickle(fp+'TCR/'+'{}_trend_output.{}.p'.format(nm,vv_tcr))
for k in tcr_dict_coanth:
    locals()[k] = tcr_dict_coanth[k]


# In[63]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(tcr_co_anth_trend,tcr_time).flatten(),bins=50)


# In[64]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(tcr_co_anth_trend,tcr_time),
                          tcr_co_anth_trend_pval,tcr['lon'],tcr['lat'],
                          name=nm+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(tcr_time[0],tcr_time[-1]),
                          cax_name=nm+' trend [DU/decade]',figsize=(10,4),
                          clevels=np.arange(-0.01, 0.01, 0.0025),vmin=-0.2e-9,vmax=0.2e-9)
fig.tight_layout()
fig.savefig(fp+'TCR/'+nm+'_trend_output.{}.png'.format(vv_tcr),dpi=600,transparent=True)


# ### TCR NOX SFC

# In[65]:


nm = 'NOX_SFC'


# In[51]:


tcr_nox_sfc_trend,tcr_nox_sfc_trend_pval = trends_and_pval_per_lat_lon_single(tcr['nox_sfc'],
                                                                              tcr['lat'],tcr['lon'],
                                                                              tcr_time,name=nm)


# In[52]:


pickle.dump({'tcr_nox_sfc_trend':tcr_nox_sfc_trend,'tcr_nox_sfc_trend_pval':tcr_nox_sfc_trend_pval,
             'tcr_time':tcr_time,'tcr':tcr},open(fp+'TCR/'+nm+'_trend_output.{}.p'.format(vv_tcr),"wb"))


# In[66]:


tcr_dict_tmp = pd.read_pickle(fp+'TCR/'+'{}_trend_output.{}.p'.format(nm,vv_tcr))
for k in tcr_dict_tmp:
    locals()[k] = tcr_dict_tmp[k]


# In[67]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(tcr_nox_sfc_trend,tcr_time).flatten(),bins=50)


# In[137]:


np.percentile(convert_monthly_to_decadal_trend(tcr_nox_sfc_trend,tcr_time),0.2)


# In[69]:


fig = plot_trend_and_pval(tcr_nox_sfc_trend,tcr_nox_sfc_trend_pval,tcr['lon'],tcr['lat'],
                          name=nm+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(tcr_time[0],tcr_time[-1]),
                          cax_name=nm+' trend [DU/month]',figsize=(10,4),
                          clevels=np.arange(-0.01, 0.01, 0.0025),vmin=-0.15e-10,vmax=0.15e-10)
fig.tight_layout()
fig.savefig(fp+'TCR/'+nm+'_trend_output.{}.png'.format(vv_tcr),dpi=600,transparent=True)


# ### TCR NOX LIGHT

# In[70]:


nm = 'TCR_NOX_LIGHT'


# In[55]:


tcr_nox_light_trend,tcr_nox_light_trend_pval = trends_and_pval_per_lat_lon_single(tcr['nox_light'],
                                                                              tcr['lat'],tcr['lon'],
                                                                              tcr_time,name=nm)


# In[56]:


pickle.dump({'tcr_nox_light_trend':tcr_nox_light_trend,'tcr_nox_light_trend_pval':tcr_nox_light_trend_pval,
             'tcr_time':tcr_time,'tcr':tcr},open(fp+'TCR/'+nm+'_trend_output.{}.p'.format(vv_tcr),"wb"))


# In[71]:


tcr_dict_tmp = pd.read_pickle(fp+'TCR/'+'{}_trend_output.{}.p'.format(nm,vv_tcr))
for k in tcr_dict_tmp:
    locals()[k] = tcr_dict_tmp[k]


# In[72]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(tcr_nox_light_trend,tcr_time).flatten(),bins=50)


# In[73]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(tcr_nox_light_trend,tcr_time),
                          tcr_nox_light_trend_pval,tcr['lon'],tcr['lat'],
                          name=nm+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(tcr_time[0],tcr_time[-1]),
                          cax_name=nm+' trend [DU/decade]',figsize=(10,4),
                          clevels=np.arange(-0.01, 0.01, 0.0025),vmin=-0.6e-12,vmax=0.5e-12)
fig.tight_layout()
fig.savefig(fp+'TCR/'+nm+'_trend_output.{}.png'.format(vv_tcr),dpi=600,transparent=True)


# ### TCR NOX SOIL

# In[74]:


nm = 'TCR_NOX_SOIL'


# In[58]:


tcr_nox_soil_trend,tcr_nox_soil_trend_pval = trends_and_pval_per_lat_lon_single(tcr['nox_soil'],
                                                                              tcr['lat'],tcr['lon'],
                                                                              tcr_time,name=nm)


# In[59]:


pickle.dump({'tcr_nox_soil_trend':tcr_nox_soil_trend,'tcr_nox_soil_trend_pval':tcr_nox_soil_trend_pval,
             'tcr_time':tcr_time,'tcr':tcr},open(fp+'TCR/'+nm+'_trend_output.{}.p'.format(vv_tcr),"wb"))


# In[75]:


tcr_dict_tmp = pd.read_pickle(fp+'TCR/'+'{}_trend_output.{}.p'.format(nm,vv_tcr))
for k in tcr_dict_tmp:
    locals()[k] = tcr_dict_tmp[k]


# In[76]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(tcr_nox_soil_trend,tcr_time).flatten(),bins=50)


# In[77]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(tcr_nox_soil_trend,tcr_time),
                          tcr_nox_soil_trend_pval,tcr['lon'],tcr['lat'],
                          name=nm+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(tcr_time[0],tcr_time[-1]),
                          cax_name=nm+' trend [DU/decade]',figsize=(10,4),
                          clevels=np.arange(-0.01, 0.01, 0.0025),vmin=-0.6e-12,vmax=0.6e-12)
fig.tight_layout()
fig.savefig(fp+'TCR/'+nm+'_trend_output.{}.png'.format(vv_tcr),dpi=600,transparent=True)


# ## Plot GFED

# In[161]:


gfed.keys()


# ### GFED C emissions

# In[175]:


nm = 'GFED_EMI_C'


# In[170]:


trend_and_plot(gfed['emi_C'],gfed['lat'],gfed['lon'],gfed['time'],'GFED_EMI_C',
               fp=fp+'GFED/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# In[171]:


gfed_dict = pd.read_pickle('/data2/ACCDAM_low_ozone/GFED/GFED_EMI_C_trend_output.v1.p')


# In[172]:


gfed_dict.keys()


# In[174]:


plt.figure()
plt.hist(convert_monthly_to_decadal_trend(gfed_dict['GFED_EMI_C_trend'],gfed_dict['GFED_EMI_C_time']).flatten(),bins=100)


# In[181]:


fig = plot_trend_and_pval(convert_monthly_to_decadal_trend(gfed_dict['GFED_EMI_C_trend'],gfed_dict['GFED_EMI_C_time']),
                          gfed_dict['GFED_EMI_C_trend_pval'],gfed['lon'],gfed['lat'],
                          name=nm+' De-seasonalized [{:%Y/%m}-{:%Y/%m}]'.format(gfed['time'][0],gfed['time'][-1]),
                          cax_name=nm+' trend [g/kg /decade]',figsize=(10,4),npval_skip=80,
                          clevels=np.arange(-0.01, 0.01, 0.0025),vmin=-5,vmax=5)
fig.tight_layout()
fig.savefig(fp+'GFED/'+nm+'_trend_output.{}.png'.format(vv),dpi=600,transparent=True)


# In[182]:


np.nanmean(gfed['emi_C'],axis=0).shape


# In[184]:


fig = plot_avg_and_rmse(np.nanmean(gfed['emi_C'],axis=0),gfed_dict['GFED_EMI_C_trend_rmse'],gfed_dict['GFED_EMI_C_trend'],
                        gfed['lon'],gfed['lat'],name='GFED_EMI_C',cax_name='Average '+'GFED_EMI_C [g/kg]',figsize=(10,4),
                        clevels=np.arange(-0.025, 0.04, 0.005),vmin=None,vmax=None,ncolors=12,nlevels=12)

fig.savefig(fp+'GFED/'+nm+'_averages.{}.png'.format(vv),dpi=600,transparent=True)


# ## Plot GMI / MERRA2

# In[81]:


gmi.keys()


# In[97]:


nm = 'no_lightning_emission'


# In[99]:


trend_and_plot(gmi[nm],gmi['lat'],gmi['lon'],gmi['time'],nm,
               fp=fp+'MERRA2/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# In[100]:


nm = 'O3prod_PBL'
trend_and_plot(gmi[nm],gmi['lat'],gmi['lon'],gmi['time'],nm,
               fp=fp+'MERRA2/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# In[101]:


for nm in ['O3prod_TROP', 'O3strat_PBL', 'O3strat_TROP']:
    trend_and_plot(gmi[nm],gmi['lat'],gmi['lon'],gmi['time'],nm,
               fp=fp+'MERRA2/',clevels=np.arange(-0.01,0.01,0.0025),vv=vv)


# In[ ]:




