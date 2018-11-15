
# coding: utf-8

# # Info
# Name:  
# 
#     KORUS_cloud_aerosol_profile
# 
# Purpose:  
# 
#     Comparison of aerosol number and cloud profiles
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures and save files...
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - load_utils.py : for loading ict and hdf5 files
#     - matplotlib
#     - numpy
#     - scipy : for saving and reading
#     - os
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, NASA Ames Research Center, 2018-08-16
#     Modified: 

# # Load the python modules and setup the environment

# In[1]:


import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.io as sio
import Sp_parameters as Sp
import tables
import load_utils as lm
from mpl_toolkits.basemap import Basemap,cm
from Sp_parameters import deriv, smooth
from path_utils import getpath
import hdf5storage as hs


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


fp = getpath('KORUS')
fp


# # Load the files for plotting

# ## Load the merge files

# In[4]:


mrg,mrg_info = lm.load_ict(fp+'data_other/korusaq-mrg01-dc8_merge_20160512_R4.ict',return_header=True)


# In[5]:


mrg_info


# ## Load DIAL h5 data

# In[6]:


val = (('AOT',0),('alt_prof',1),('cld_b',3),('cld_t',5),('cld_m',4),('ext532',19),('Alt',21),
       ('GPS_alt',24),('GPS_lat',25),('GPS_lon',26),('time',23))
dial,dial_info = lm.load_hdf(fp+'data_other/korusaq-DIAL_DC8_20160512_R1.h5',values=val)


# In[7]:


dial['time'][0,280]/3600.0


# In[8]:


dial_info['Alt']


# In[9]:


dial['cld_t'][:,280],dial['cld_b'][:,280]


# In[10]:


plt.figure()
plt.plot(dial['cld_m'][:,280],dial['Alt'],'.')


# In[89]:


plt.figure()
plt.plot(dial['ext532'][:,273:276]+0.02,dial['Alt'],'.')


# In[59]:


ic = 293


# In[61]:


dial['cld_t'][:,280:300]


# In[54]:


dial['cld_t'][:,ic],dial['cld_b'][:,ic]


# In[92]:


plt.figure()
plt.plot(dial['ext532'][:,273]+0.02,dial['Alt'],'.',label='DIAL aerosol extinction')
plt.axhspan(0.0,dial['cld_t'][0,293],color='grey',alpha=0.1,label='DIAL Cloud')

for i in xrange(1,100+1):
    top = dial['cld_t'][0,293]
    bot = dial['cld_t'][0,293]-dial['cld_t'][0,293]/100.0*(i-1)
    plt.axhspan(bot,top,color='grey',alpha=0.005-0.005/100.0*i)

#plt.axhspan(dial['cld_t'][0,293]*0.6,dial['cld_t'][0,293]*0.79,color='grey',alpha=0.15)
#plt.axhspan(dial['cld_t'][0,293]*0.4,dial['cld_t'][0,293]*0.59,color='grey',alpha=0.1)
#plt.axhspan(dial['cld_t'][0,293]*0.2,dial['cld_t'][0,293]*0.19,color='grey',alpha=0.05)
#plt.axhspan(0.0,dial['cld_t'][0,293]*0.19,color='grey',alpha=0.0)
plt.legend(frameon=False)
plt.xlim(0,0.1)
plt.xlabel('Aerosol extinction [km$^{{-1}}$]')
plt.ylabel('Altitude [m]')


# ## Load the MERRA-2 files

# In[12]:


fpma = fp+'data_other/MERRA2/'


# In[13]:


fpma


# In[14]:


day = '20160512'


# In[15]:


mcld,mcld_info = lm.load_netcdf(fpma+'MERRA2_400.tavg3_3d_cld_Nv.{}.SUB.nc'.format(day),everything=True)


# In[16]:


maero,maero_info = lm.load_netcdf(fpma+'MERRA2_400.inst3_3d_aer_Nv.{}.SUB.nc'.format(day),everything=True)


# In[17]:


ll=[]
for k in maero_info.keys():
    print k,maero_info[k].long_name, maero_info[k].units
    if maero_info[k].units=='kg kg-1': ll.append(k)


# In[18]:


ll


# In[19]:


maero['tot_aero'] = 0
for l in ll:
    maero['tot_aero'] = maero['tot_aero']+maero[l]


# In[20]:


maero['scat'] = maero['tot_aero']*maero['AIRDENS']*1000.0*2.18*1000.0*1000.0


# In[21]:


maero['scat'].shape


# In[22]:


matm,matm_info = lm.load_netcdf(fpma+'MERRA2_400.tavg3_3d_asm_Nv.{}.SUB.nc'.format(day),everything=True)


# In[23]:


matm['H'].shape


# In[24]:


matm['lat'].shape,matm['lon'].shape


# In[25]:


matm['time'][7]/60.0


# In[19]:


mrg['LATITUDE'][ip][0],mrg['LONGITUDE'][ip][0]


# In[26]:


lat,lon = matm['lat'][7],matm['lon'][2]


# In[27]:


maero['isaero'] = maero['scat']>50.0
mcld['iscld'] = mcld['QL']>0.00001


# In[28]:


dz0 = (matm['H'][7,:,7,2][np.where(mcld['iscld'][7,:,7,2])[0][0]]+matm['H'][7,:,7,2][np.where(mcld['iscld'][7,:,7,2])[0][0]-1])/2.0
dz1 = (matm['H'][7,:,7,2][np.where(mcld['iscld'][7,:,7,2])[0][0]+1]+matm['H'][7,:,7,2][np.where(mcld['iscld'][7,:,7,2])[0][0]])/2.0


# In[29]:


plt.figure()
plt.plot(matm['H'][7,1:,7,2],np.diff(matm['H'][7,:,7,2]),'x-')
plt.plot(matm['H'][7,1:,10,2],np.diff(matm['H'][7,:,10,2]),'x-')
plt.xlabel('Altitude [m]')
plt.ylabel('Altitude bin difference [m]')


# In[32]:


ilow = matm['H'][7,1:,7,2]<7000.0
iloww = matm['H'][7,1:,7,2]<3000.0


# In[33]:


np.nanmean(np.diff(matm['H'][7,:,7,2])[ilow]),np.nanmean(np.diff(matm['H'][7,:,7,2])[iloww])


# In[34]:


np.diff(matm['H'][7,:,7,2])[ilow]


# In[35]:


plt.figure()
plt.plot(maero['scat'][7,:,7,2]/1000.0,matm['H'][7,:,7,2],'.',label='MERRA-2 aerosol scattering')
plt.axhspan(dz0,dz1,color='grey',alpha=0.2,label='MERRA-2 Cloud')
plt.ylim(0,12000)
plt.legend(frameon=False)
plt.xlabel('Scattering coefficient [km$^{{-1}}$]')
plt.ylabel('Altitude [m]')


# # Plot some data

# In[23]:


plt.figure()
plt.plot(mrg['UTC'],mrg['HAE_GPS_ALT'],'.')


# In[93]:


ip = (mrg['UTC']>23.75)&(mrg['UTC']<23.94)


# In[25]:


plt.figure()
plt.plot(mrg['UTC'],mrg['HAE_GPS_ALT'],'.')
plt.plot(mrg['UTC'][ip],mrg['HAE_GPS_ALT'][ip],'r.')


# In[21]:


plt.figure()
plt.plot(mrg['UTC'],mrg['CWC'],'.')


# In[26]:


plt.figure()
plt.plot(mrg['SCAT550nmamb_total_stdPT_LARGE'][ip]/1000.0,mrg['HAE_GPS_ALT'][ip],'.')
plt.plot(mrg['CWC'][ip],mrg['HAE_GPS_ALT'][ip],'r.')


# In[94]:


mrg['HAE_GPS_ALT'][ip][np.where(mrg['CWC'][ip]>0.01)[0][0]]


# In[95]:


plt.figure()
plt.plot(mrg['SCAT550nmamb_total_stdPT_LARGE'][ip]/1000.0,mrg['HAE_GPS_ALT'][ip],'.',label='In situ 550 nm scattering')
plt.axhspan(mrg['HAE_GPS_ALT'][ip][np.where(mrg['CWC'][ip]>0.01)[0][-1]],
             mrg['HAE_GPS_ALT'][ip][np.where(mrg['CWC'][ip]>0.01)[0][0]],color='grey',alpha=0.2,label='In situ Cloud')
plt.legend(frameon=False)
plt.xlabel('Aerosol scattering [km$^{{-1}}$]')
plt.xlim(0,0.05)


# In[58]:


mrg['LATITUDE'][ip][0],mrg['LONGITUDE'][ip][0]


# ## Combine plots

# In[96]:


plt.figure(figsize=(10,7))
ax1 = plt.subplot(1,3,1)

plt.plot(mrg['SCAT550nmamb_total_stdPT_LARGE'][ip]/1000.0,mrg['HAE_GPS_ALT'][ip]*1000.0,'.',label='In situ 550 nm scattering')
plt.axhspan(mrg['HAE_GPS_ALT'][ip][np.where(mrg['CWC'][ip]>0.01)[0][-1]]*1000.0,
             mrg['HAE_GPS_ALT'][ip][np.where(mrg['CWC'][ip]>0.01)[0][0]]*1000.0,color='grey',alpha=0.2,label='In situ Cloud')
#plt.legend(frameon=False)
plt.xlabel('Scattering\ncoefficient [km$^{{-1}}$]')
plt.xlim(0,0.05)
plt.ylabel('Altitude [m]')
plt.title('In Situ')

ax2 = plt.subplot(1,3,2,sharey=ax1)
plt.plot(dial['ext532'][:,280],dial['Alt'],'.',label='Aerosol')
plt.axhspan(0.0,dial['cld_t'][0,280],color='grey',alpha=0.2,label='Cloud')
plt.legend(frameon=False)
plt.xlim(0,0.05)
plt.xlabel('Extinction\ncoefficient [km$^{{-1}}$]')
plt.title('DIAL')
#plt.ylabel('Altitude [m]')
plt.setp(ax2.get_yticklabels(), visible=False)

ax3 = plt.subplot(1,3,3,sharey=ax1)
plt.plot(maero['scat'][7,:,7,2]/1000.0,matm['H'][7,:,7,2],'.',label='MERRA-2 aerosol scattering')
plt.axhspan(dz0,dz1,color='grey',alpha=0.2,label='MERRA-2 Cloud')
plt.ylim(0,12000)
#plt.legend(frameon=False)
plt.xlim(0,0.05)
plt.xlabel('Scattering\ncoefficient [km$^{{-1}}$]')
plt.title('MERRA-2')
plt.setp(ax3.get_yticklabels(), visible=False)
#plt.ylabel('Altitude [m]')
plt.suptitle('KORUS-AQ aerosol cloud distribution on {}, {}$^{{\circ}}$N, {}$^{{\circ}}$E'.format(day,lat,lon))
plt.savefig(fpma+'KORUS_aerosol_cloud_insitu_dial_merra_{}.png'.format(day),transparent=True,di=600)


# In[106]:


plt.figure(figsize=(8,5.8))
ax1 = plt.subplot(1,3,1)

plt.plot(mrg['SCAT550nmamb_total_stdPT_LARGE'][ip]/1000.0,mrg['HAE_GPS_ALT'][ip]*1000.0,'.',label='Aerosol')
plt.axhspan(mrg['HAE_GPS_ALT'][ip][np.where(mrg['CWC'][ip]>0.01)[0][-1]]*1000.0,
             mrg['HAE_GPS_ALT'][ip][np.where(mrg['CWC'][ip]>0.01)[0][0]]*1000.0,color='grey',alpha=0.2,label='Cloud')
plt.legend(frameon=False)
plt.xlabel('Scattering [km$^{{-1}}$]')
plt.xlim(0,0.05)
plt.ylabel('Altitude [m]')
plt.title('In Situ')

ax2 = plt.subplot(1,3,3,sharey=ax1)
plt.plot(dial['ext532'][:,273]+0.02,dial['Alt'],'.',label='Aerosol')
plt.axhspan(0.0,dial['cld_t'][0,293],color='grey',alpha=0.1,label='Cloud')

for i in xrange(1,100+1):
    top = dial['cld_t'][0,293]
    bot = dial['cld_t'][0,293]-dial['cld_t'][0,293]/100.0*(i-1)
    plt.axhspan(bot,top,color='grey',alpha=0.005-0.005/100.0*i)
#plt.legend(frameon=False)
plt.xlim(0,0.1)
plt.xlabel('Extinction [km$^{{-1}}$]')
plt.title('nearest DIAL')
#plt.ylabel('Altitude [m]')
plt.setp(ax2.get_yticklabels(), visible=False)

ax3 = plt.subplot(1,3,2,sharey=ax1)
plt.plot(maero['scat'][7,:,7,2]/1000.0,matm['H'][7,:,7,2],'.',label='Aerosol')
plt.axhspan(dz0,dz1,color='grey',alpha=0.2,label='Cloud')
plt.ylim(0,12000)
#plt.legend(frameon=False)
plt.xlim(0,0.05)
plt.xlabel('Scattering [km$^{{-1}}$]')
plt.title('nearest MERRA-2')
plt.setp(ax3.get_yticklabels(), visible=False)
#plt.ylabel('Altitude [m]')
plt.suptitle('KORUS-AQ aerosol cloud distribution on {}, {}$^{{\circ}}$N, {}$^{{\circ}}$E'.format(day,lat,lon))
plt.savefig(fpma+'KORUS_aerosol_cloud_insitu_dial_merra_{}_v2.png'.format(day),transparent=True,di=600)

