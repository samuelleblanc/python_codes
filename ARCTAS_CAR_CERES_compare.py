#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Check out the BRDF comparison of CAR to CERES from ARCTAS
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-01-24
#     Modified:
# 

# # Prepare python environment

# In[3]:


import numpy as np
import load_utils as lu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import os


# In[26]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature


# In[121]:


from datetime import datetime, timedelta


# In[4]:


name = 'ARCTAS'
vv = 'v1'
fp = getpath(name)


# # Load files

# In[111]:


os.listdir(fp+'nasa/20080409/')


# In[8]:


fpp = fp+'nasa/20080409/'


# ## Load P3 flight path

# In[23]:


daystr = '20080409'


# In[12]:


nav,nav_dict = lu.load_ict(fpp+'pds_p3b_20080409_r2.ict',return_header=True)


# In[164]:


lat_rg = [79.7,80.1]
lon_rg = [-101.1,-100.2]


# In[80]:


ifl = np.where((nav['FMS_LON']<-100.2)&(nav['FMS_LON']>-101.1)&(nav['FMS_LAT']<80.1)&(nav['FMS_LAT']>79.7))[0]


# In[81]:


ifl


# In[199]:


fig = plt.figure(figsize=(9,2.5))
proj = ccrs.PlateCarree()#(-120, 45)#(central_latitude=np.nanmean(nav['FMS_LAT']), central_longitude=np.nanmean(nav['FMS_LON']))
ax = fig.add_subplot(111,projection=proj)
axs = ax.scatter(nav['FMS_LON'],nav['FMS_LAT'],c=nav['GPS_ALT'],marker='o',s=6)
ax.scatter(nav['FMS_LON'][ifl],nav['FMS_LAT'][ifl],c='r',marker='x',s=12)

ax.coastlines(resolution='50m')
ax.gridlines(draw_labels=True,auto_update=True)

cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])
fig.subplots_adjust(hspace=0.0, wspace=0, top=0.925, left=0.06)
posn = ax.get_position()
ax.set_position([posn.x0, posn.y0, posn.width-0.1, posn.height])
posn = ax.get_position()
cbar_ax.set_position([posn.x0 + posn.width + 0.07, posn.y0,
                          0.02, posn.height])

cbar = plt.colorbar(axs,cax=cbar_ax,label='GPS Alt [feet]')

ax.set_title('ARCTAS - NASA P-3 - {}'.format(daystr))
plt.savefig(fpp+'ARCTAS_P3_map_{}.png'.format(daystr),dpi=400,transparent=True)


# In[90]:


nav['Time'][ifl]/3600.0


# In[91]:


time = [nav['Time'][ifl[0]]/3600.0,nav['Time'][ifl[-1]]/3600.0]


# ## Plot out the altitude

# In[208]:


plt.figure()
plt.plot(nav['Time']/3600.0,nav['GPS_ALT'],'.')
plt.axvspan(time[0],time[1],alpha=0.4,color='grey',label='Spiral')
plt.xlabel('Time [h]')
plt.ylabel('Altitude [feet]')


# In[214]:


plt.figure()
plt.plot(nav['Time']/3600.0,nav['FMS_LAT'],'.')
plt.axvspan(time[0],time[1],alpha=0.4,color='grey',label='Spiral')
plt.xlabel('time [h]')
plt.ylabel('Latitude [$^\circ$]')
plt.xlim(18.3,19.1)
plt.ylim(79.0,80.25)


# ## Load CAR

# In[ ]:





# ## Load CERES for Arctic April 2008

# In[114]:


ceres,ceres_dict = lu.load_netcdf(fpp+'CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20080401-20080531.nc',everything=True)


# In[122]:


ceres_dict[b'aux_snow_1h']


# In[117]:


ceres[b'adj_atmos_sw_up_all_surface_1h'].shape


# In[197]:


ceres_dict[b'ini_albedo_1h']


# ### Convert to useable time units

# In[120]:


ceres_dict[b'time']


# In[118]:


ceres[b'time']


# In[158]:


time2date = lambda x : datetime(2000,3,1)+timedelta(days=1)*x


# In[159]:


ceres_date = time2date(ceres[b'time'])


# In[160]:


ceres_date[0:3]


# ### find region subset

# In[162]:


ceres[b'lon']


# In[163]:


ceres[b'lat']


# In[168]:


lat_rg,lon_rg[0]+360.0,lon_rg[1]+360.0


# In[171]:


ceres[b'lat'][9],ceres[b'lon'][39]


# In[172]:


ilat = 9
ilon = 39


# ### plot out time series

# In[192]:


time_dt = [datetime(2008,4,9)+timedelta(hours=time[0]),datetime(2008,4,9)+timedelta(hours=time[1])]


# In[193]:


fig = plt.figure()
plt.plot(ceres_date,ceres[b'adj_atmos_sw_up_all_surface_1h'][:,ilat,ilon],'.',label='surf up')
plt.plot(ceres_date,ceres[b'adj_atmos_sw_down_all_surface_1h'][:,ilat,ilon],'.',label='surf down')
plt.ylabel('Shortwave Irradiance [W/m$^2$]')
plt.xlim([datetime(2008,4,7),datetime(2008,4,12)])
fig.autofmt_xdate()
plt.legend(frameon=False)
plt.title('CERES adjusted irradiance over ARCTAS P3 - {}'.format(daystr))
plt.axvspan(time_dt[0],time_dt[1],alpha=0.4,color='grey',label='Spiral')


# In[194]:


plt.figure()
plt.plot(ceres_date,ceres[b'ini_albedo_1h'][:,ilat,ilon],'.')
plt.ylabel('Albedo init')
plt.xlim([datetime(2008,4,7),datetime(2008,4,12)])
fig.autofmt_xdate()
plt.title('CERES alebdo init over ARCTAS P3 - {}'.format(daystr))
plt.axvspan(time_dt[0],time_dt[1],alpha=0.4,color='grey',label='Spiral')


# In[196]:


plt.figure()
plt.plot(ceres_date,ceres[b'ini_aod55_1h'][:,ilat,ilon],'.')
plt.ylabel('AOD init at 550 nm')
plt.xlim([datetime(2008,4,7),datetime(2008,4,12)])
fig.autofmt_xdate()
plt.title('CERES AOD init over ARCTAS P3 - {}'.format(daystr))
plt.axvspan(time_dt[0],time_dt[1],alpha=0.4,color='grey',label='Spiral')


# In[198]:


plt.figure()
plt.plot(ceres_date,ceres[b'aux_ocean_1h'][:,ilat,ilon],'.',label='ocean')
plt.plot(ceres_date,ceres[b'aux_snow_1h'][:,ilat,ilon],'.',label='snow')
plt.ylabel('Fractional surface condition')
plt.xlim([datetime(2008,4,7),datetime(2008,4,12)])
fig.autofmt_xdate()
plt.title('CERES surface condition over ARCTAS P3 - {}'.format(daystr))
plt.axvspan(time_dt[0],time_dt[1],alpha=0.4,color='grey',label='Spiral')
plt.legend()


# ## Load SSFR surface albedo ict file

# In[85]:


ssfr, ssfr_dict = lu.load_ict(fpp+'SSFR_P3B_20080409_R0.ict',return_header=True)


# ### plot out time series

# In[96]:


plt.figure()
plt.plot(ssfr['UTC'],ssfr['UP_BB350_700']/ssfr['DN_BB350_700'],'.',label='350-700 nm')
plt.plot(ssfr['UTC'],ssfr['UP_BB350_2150']/ssfr['DN_BB350_2150'],'.',label='350-2150 nm')
plt.plot(ssfr['UTC'],ssfr['UP380']/ssfr['DN380'],'.',label='380 nm')
plt.plot(ssfr['UTC'],ssfr['UP500']/ssfr['DN500'],'.',label='500 nm')
plt.plot(ssfr['UTC'],ssfr['UP1200']/ssfr['DN1200'],'.',label='1200 nm')

plt.axvspan(time[0],time[1],alpha=0.4,color='grey',label='Spiral')

plt.xlabel('UTC [h]')
plt.ylabel('Albedo')
plt.ylim(0,1)
plt.legend(frameon=True)
plt.title('ARCTAS - SSFR albedo - {}'.format(daystr))


# ## Load BBR

# In[103]:


bbr, bbr_dict = lu.load_ict(fpp+'BBR_P3B_20080409_R0.ict',return_header=True)


# In[105]:


bbr_dict


# ### plot time series

# In[108]:


plt.figure()
plt.plot(bbr['Time_UTC_secs']/3600.0,bbr['Upwelling_SOLAR_Irradiance_Wm2']/bbr['Downwelling_SOLAR_Irradiance_Wm2'],'.',label='solar')
plt.plot(bbr['Time_UTC_secs']/3600.0,bbr['Upwelling_IR_Irradiance_Wm2']/bbr['Downwelling_IR_Irradiance_Wm2'],'.',label='IR')

plt.axvspan(time[0],time[1],alpha=0.4,color='grey',label='Spiral')

plt.xlabel('UTC [h]')
plt.ylabel('Broad Band Albedo (up/dn)')
plt.ylim(0,4.0)
plt.grid()
plt.legend(frameon=True)
plt.title('ARCTAS - BBR - {}'.format(daystr))


# In[110]:


plt.figure()
plt.plot(bbr['Time_UTC_secs']/3600.0,bbr['Upwelling_SOLAR_Irradiance_Wm2'],'.',label='solar up')
plt.plot(bbr['Time_UTC_secs']/3600.0,bbr['Downwelling_SOLAR_Irradiance_Wm2'],'.',label='solar down')
plt.plot(bbr['Time_UTC_secs']/3600.0,bbr['Upwelling_IR_Irradiance_Wm2'],'.',label='IR up')
plt.plot(bbr['Time_UTC_secs']/3600.0,bbr['Downwelling_IR_Irradiance_Wm2'],'.',label='IR down')

plt.axvspan(time[0],time[1],alpha=0.4,color='grey',label='Spiral')

plt.xlabel('UTC [h]')
plt.ylabel('Broad Band Irradiance [W/m$^2$]')
#plt.ylim(0,4.0)
plt.grid()
plt.legend(frameon=True)
plt.title('ARCTAS - BBR - {}'.format(daystr))


# ## Load AATS-14 AOD

# In[97]:


aats,aats_dict = lu.load_ict(fpp+'AATS14_P3B_20080409_R2.ict',return_header=True)


# ### plot time series

# In[100]:


plt.figure()
fqa = aats['cld_scr']>0
plt.plot(aats['Start_UTC'][fqa],aats['AOD354'][fqa],'.',label='354 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD380'][fqa],'.',label='380 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD453'][fqa],'.',label='453 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD499'][fqa],'.',label='499 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD519'][fqa],'.',label='519 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD606'][fqa],'.',label='606 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD675'][fqa],'.',label='675 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD779'][fqa],'.',label='779 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD865'][fqa],'.',label='865 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD1019'][fqa],'.',label='1019 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD1241'][fqa],'.',label='1241 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD1559'][fqa],'.',label='1559 nm')
plt.plot(aats['Start_UTC'][fqa],aats['AOD2139'][fqa],'.',label='2139 nm')

plt.axvspan(time[0],time[1],alpha=0.4,color='grey',label='Spiral')

plt.xlabel('UTC [h]')
plt.ylabel('AOD')
plt.ylim(0,0.256)
plt.grid()
plt.legend(frameon=True)
plt.title('ARCTAS - AATS-14 AOD - {}'.format(daystr))


# # Add notes from the flight

# ![image.png](attachment:image.png)

# ![image.png](attachment:image.png)

# The Sea Ice Concentration (L4, MUR25) layer is a Group for High Resolution Sea Surface Temperature (GHRSST) Level 4 sea surface temperature analysis. It is produced as a retrospective dataset (four day latency) and near-real-time dataset (one day latency) at the JPL Physical Oceanography DAAC using wavelets as basis functions in an optimal interpolation approach on a global 0.25 degree grid. The version 4 Multiscale Ultrahigh Resolution (MUR) L4 analysis is based upon nighttime GHRSST L2P skin and subskin SST observations from several instruments including the NASA Advanced Microwave Scanning Radiometer-EOS (AMSR-E), the JAXA Advanced Microwave Scanning Radiometer 2 on GCOM-W1, the Moderate Resolution Imaging Spectroradiometers (MODIS) on the NASA Aqua and Terra platforms, the US Navy microwave WindSat radiometer, the Advanced Very High Resolution Radiometer (AVHRR) on several NOAA satellites, and in situ SST observations from the NOAA iQuam project. The ice concentration data are from the archives at the EUMETSAT Ocean and Sea Ice Satellite Application Facility (OSI SAF) High Latitude Processing Center and are also used for an improved SST parameterization for the high-latitudes. The dataset also contains an additional SST anomaly variable derived from a MUR climatology.
# 
# This dataset is funded by the NASA MEaSUREs program, and created by a team led by Dr. Toshio M. Chin from JPL. It adheres to the GHRSST Data Processing Specification (GDS) version 2 format specifications.
# 
# References: MUR25-JPL-L4-GLOB-v04.2 doi:10.5067/GHM25-4FJ42

# CAR circle image
# ![ARCTAS2008-car_p3b_20080409184734_680.png](attachment:ARCTAS2008-car_p3b_20080409184734_680.png)

# # Subset data for the same time periods

# In[ ]:





# # Plot out data

# In[ ]:




