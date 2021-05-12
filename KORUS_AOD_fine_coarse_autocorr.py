#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     KORUS_AOD_fine_coarse_autocorr
#     
# Purpose:  
# 
#     Analyse some of the AOD values from KORUS AQ
#     Split up between fine mode and coarse mode AOD
#     Subset for level legs only
#         - interpolate within the level legs
#     Run autocorrelation values for the distance/time travelled 
#   
# Input:
# 
#     None at command line
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
#     - hdf5 python loader
#     - 
#     - matplotlib
#     - numpy
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - '/aod_ict/all_aod_KORUS_R2_ict.mat'
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-05-18
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2020-01-14
#                 Added plotting of maps statistics

# # Prepare python environment

# In[1]:


#%config InlineBackend.rc = {}
import matplotlib 
import os
#matplotlib.rc_file(os.path.join(os.getcwd(),'file.rc'))
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
import load_utils as lu
import plotting_utils as pu
from path_utils import getpath
import hdf5storage as hs
from datetime import datetime
from scipy.interpolate import UnivariateSpline
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
import scipy.stats as st
import scipy.io as sio
import pandas as pd


# In[2]:


import map_utils as mu
import write_utils as wu
from scipy import interpolate
import math
import sys


# In[3]:


from linfit import linfit
import Sun_utils as su


# In[4]:


get_ipython().magic(u'matplotlib notebook')


# In[5]:


fp = getpath('KORUS')


# # Load files

# Load the KORUS 4STAR AOD ict files for better handling

# ## Load 4STAR AOD ict files

# In[6]:


ar = hs.loadmat(fp+'/aod_ict/all_aod_KORUS_R2_ict.mat')


# In[7]:


ka = ar.keys()


# In[8]:


ka.sort()


# In[9]:


ka


# In[10]:


ar['qual_flag']


# In[11]:


ar['fl_QA']


# In[12]:


ar['qual_flag'].sum()/float(len(ar['qual_flag']))


# In[13]:


fl_ci = (ar['qual_flag']==1.) & (ar['AOD0501']<1.0)


# In[14]:


fl_ci.sum()


# In[15]:


len(fl_ci)*1.0/60.0/60.0


# In[16]:


fl_ci.sum()/float(len(fl_ci))


# In[17]:


ar['fl_QA'].sum()*1.0/60.0/60.0


# In[18]:


ar['qual_flag'].sum()


# In[19]:


fl_ci


# In[20]:


nwl = ka[0:17]


# In[21]:


nwl


# In[22]:


nm = [380.0,452.0,501.0,520.0,532.0,550.0,606.0,620.0,675.0,781.0,865.0,1020.0,1040.0,1064.0,1236.0,1559.0,1627.0]


# In[23]:


days = ['20160501','20160503','20160504','20160506','20160510','20160511',
        '20160512','20160516','20160517','20160519','20160521','20160524',
        '20160526','20160529','20160601','20160602','20160604','20160608',
        '20160609','20160614']


# In[24]:


doys = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in days]


# In[25]:


doys


# In[26]:


fdoys = np.array(doys)


# In[27]:


ar['doys'] = fdoys[ar['days'].astype(int)]+ar['Start_UTC']/24.0


# In[28]:


ar['doys']


# ## Load GOCI AOD files

# In[29]:


gl = os.listdir(fp+'data_other/GOCI')


# In[30]:


gl.sort()


# In[31]:


goci = []


# In[ ]:


if false:
    for d in days:
        for q in xrange(8):
            print 'loading : GOCI_YAER_V2_AHICLD_AOPs_{d}0{q}.hdf'.format(d=d,q=q)
            g_tmp,g_dict = lu.load_hdf(fp+'data_other/GOCI/GOCI_YAER_V2_AHICLD_AOPs_{d}0{q}.hdf'.format(d=d,q=q),
                                       values=(('lat',1),('lon',0),('aod',2),('fmf',3),('ssa',4),('AE',6),('CF',10),('t',7)),
                                       verbose=False)
            goci.append(g_tmp)  


# In[32]:


for l in gl:
    print 'loading file: {}'.format(l)
    g_tmp,g_dict = lu.load_hdf(fp+'data_other/GOCI/{}'.format(l),
                                   values=(('lat',1),('lon',0),('aod',2),('fmf',3),('ssa',4),('AE',6),('CF',10),('t',7)),
                                   verbose=False)
    goci.append(g_tmp)  


# In[33]:


len(goci)


# In[34]:


np.nanmax(goci[0]['t']),np.nanmean(goci[0]['t']),np.nanmin(goci[0]['t'])


# In[35]:


g_dict


# In[36]:


gl[0][-14:-10],gl[0][-10:-8],gl[0][-8:-6]


# In[37]:


goci_doy = []
for i,l in enumerate(gl):
    goci[i]['doy'] = datetime(int(l[-14:-10]),int(l[-10:-8]),int(l[-8:-6])).timetuple().tm_yday+(int(l[-6:-4])+0.5)/24.0
    #print goci[i]['doy']
    goci_doy.append(goci[i]['doy'])
    goci[i]['doys'] = goci[i]['doy']+(int(l[-6:-4])+goci[i]['t']/60.0)/24.0 
goci_doy = np.array(goci_doy)


# In[38]:


for og in goci:
    og['aod_f'] = og['fmf']*og['aod']
    og['aod_c'] = (1.0-og['fmf'])*og['aod']


# ### Build collcation of goci data to 4STAR

# In[39]:


ig4 = Sp.find_closest(goci_doy,ar['doys'])


# In[40]:


ig4


# In[47]:


plt.figure()
plt.plot(goci_doy[ig4],'.',label='goci')
plt.plot(ar['doys'],'x',label='4star')
plt.legend()


# In[41]:


bad_ig4 = abs(goci_doy[ig4]-ar['doys'])>(1.0/24.0)


# In[42]:


sum(bad_ig4)/float(len(ig4))


# In[43]:


len(np.where(bad_ig4)[0])


# In[44]:


len(ig4)


# In[45]:


ar['doys']%1.0


# In[46]:


goci[0].keys()


# In[47]:


goci2ar = {'doys':[],'lat':[],'lon':[],'aod':[],'AE':[],'fmf':[],'CF':[],'aod_f':[],'aod_c':[]}
for ii in np.unique(ig4):
    print ii
    imeas = np.where(ig4==ii)[0]
    igoci = mu.map_ind(goci[ii]['lon'],goci[ii]['lat'],ar['Longitude'][imeas],ar['Latitude'][imeas])
    if len(igoci)<1:
        goci2ar['lat'] = np.append(goci2ar['lat'],ar['Latitude'][imeas])
        goci2ar['lon'] = np.append(goci2ar['lon'],ar['Longitude'][imeas])
        goci2ar['aod'] = np.append(goci2ar['aod'],ar['Latitude'][imeas]+np.nan)
        goci2ar['aod_c'] = np.append(goci2ar['aod_c'],ar['Latitude'][imeas]+np.nan)
        goci2ar['aod_f'] = np.append(goci2ar['aod_f'],ar['Latitude'][imeas]+np.nan)
        goci2ar['AE'] = np.append(goci2ar['AE'],ar['Latitude'][imeas]+np.nan)
        goci2ar['fmf'] = np.append(goci2ar['fmf'],ar['Latitude'][imeas]+np.nan)
        goci2ar['CF'] = np.append(goci2ar['CF'],ar['Latitude'][imeas]+np.nan)
        goci2ar['doys'] = np.append(goci2ar['doys'],ar['doys'][imeas])
    else:
        goci2ar['lat'] = np.append(goci2ar['lat'],goci[ii]['lat'][igoci[0],igoci[1]])
        goci2ar['lon'] = np.append(goci2ar['lon'],goci[ii]['lon'][igoci[0],igoci[1]])
        goci2ar['aod'] = np.append(goci2ar['aod'],goci[ii]['aod'][igoci[0],igoci[1]])
        goci2ar['aod_c'] = np.append(goci2ar['aod_c'],goci[ii]['aod_c'][igoci[0],igoci[1]])
        goci2ar['aod_f'] = np.append(goci2ar['aod_f'],goci[ii]['aod_f'][igoci[0],igoci[1]])
        goci2ar['AE'] = np.append(goci2ar['AE'],goci[ii]['AE'][igoci[0],igoci[1]])
        goci2ar['fmf'] = np.append(goci2ar['fmf'],goci[ii]['fmf'][igoci[0],igoci[1]])
        goci2ar['CF'] = np.append(goci2ar['CF'],goci[ii]['CF'][igoci[0],igoci[1]])
        goci2ar['doys'] = np.append(goci2ar['doys'],goci[ii]['doys'][igoci[0],igoci[1]])

for k in goci2ar.keys():
    goci2ar[k] = np.array(goci2ar[k])


# In[48]:


igoci


# In[49]:


len(goci2ar['aod'])


# In[50]:


goci2ar['aod'][2000]


# In[51]:


np.unique(ig4)


# In[52]:


imeas = np.where(ig4==8)[0]


# In[53]:


len(imeas)


# In[54]:


np.argmin((goci[8]['lon']-ar['Longitude'][imeas[0]])**2.0+(goci[8]['lat']-ar['Latitude'][imeas[0]])**2.0)


# In[55]:


goci[8]['lon'].flatten()[94639],goci[8]['lat'].flatten()[94639]


# In[56]:


ar['Longitude'][imeas[0]],ar['Latitude'][imeas[0]]


# In[483]:


plt.figure()
plt.plot(goci2ar['aod'])
plt.figure()
plt.plot(goci2ar['AE'])


# In[57]:


goci2ar['aod'][bad_ig4] = np.nan
goci2ar['aod_c'][bad_ig4] = np.nan
goci2ar['aod_f'][bad_ig4] = np.nan
goci2ar['AE'][bad_ig4] = np.nan
goci2ar['fmf'][bad_ig4] = np.nan
goci2ar['CF'][bad_ig4] = np.nan


# In[58]:


len(goci2ar['aod_c'])


# In[59]:


len(goci2ar['aod'])


# ### Make daily regional averages to compare to flight subsets

# In[66]:


goci[0]['lat'].min(), goci[0]['lat'].max()


# In[67]:


goci[0]['lon'].min(), goci[0]['lon'].max()


# In[68]:


np.nanmin(ar['Latitude']),np.nanmax(ar['Latitude'])


# In[69]:


np.nanmin(ar['Longitude']),np.nanmax(ar['Longitude'])


# In[70]:


ar['doys']


# In[71]:


itransit = ar['days']!=19.0


# In[72]:


for da in np.unique(ar['days']):
    itr = ar['days']==da
    print 'days: {}, lat range: {}, {}, lon range: {}, {}'.format(da,
            np.nanmin(ar['Latitude'][itr]),np.nanmax(ar['Latitude'][itr]),
            np.nanmin(ar['Longitude'][itr]),np.nanmax(ar['Longitude'][itr]))


# In[545]:


itr = (ar['days']!=19.0) &  np.isfinite(ar['Latitude']) &  np.isfinite(ar['Longitude'])
print 'days: {}, lat range: {}, {}, lon range: {}, {}'.format(da,
        np.nanmin(ar['Latitude'][itr]),np.nanmax(ar['Latitude'][itr]),
        np.nanmin(ar['Longitude'][itr]),np.nanmax(ar['Longitude'][itr]))
np.percentile(ar['Latitude'][itr],[5,10,25,75,90,95]),np.percentile(ar['Longitude'][itr],[5,10,25,75,90,95])


# In[546]:


j = 104
'GOCI: lat range: {}, {}, lon range: {}, {}'.format(
    goci[j]['lat'].min(), goci[j]['lat'].max(),goci[j]['lon'].min(), goci[j]['lon'].max())


# In[547]:


lat_rg = [33.8,37.6] # from percentiles of 4STAR data
lon_rg = [124.3,129.4]


# In[548]:


goci[0]['lat'].shape


# In[549]:


igoci_rg = (goci[0]['lat']>lat_rg[0]) & (goci[0]['lat']<lat_rg[1]) & (goci[0]['lon']>lon_rg[0]) & (goci[0]['lon']<lon_rg[1])


# In[550]:


igoci_rg.any()


# In[551]:


goci[0]['lat'][igoci_rg].shape


# In[552]:


goci_time = {}
goci_time['aod_mean'] = np.array([np.nanmean(g['aod'][igoci_rg]) for g in goci])
goci_time['aod_median'] = np.array([np.nanmedian(g['aod'][igoci_rg]) for g in goci])
goci_time['aod_std'] = np.array([np.nanstd(g['aod'][igoci_rg]) for g in goci])
goci_time['ae_mean'] = np.array([np.nanmean(g['AE'][igoci_rg]) for g in goci])
goci_time['ae_median'] = np.array([np.nanmedian(g['AE'][igoci_rg]) for g in goci])
goci_time['ae_std'] = np.array([np.nanstd(g['AE'][igoci_rg]) for g in goci])
goci_time['fmf_mean'] = np.array([np.nanmean(g['fmf'][igoci_rg]) for g in goci])
goci_time['fmf_median'] = np.array([np.nanmedian(g['fmf'][igoci_rg]) for g in goci])
goci_time['fmf_std'] = np.array([np.nanstd(g['fmf'][igoci_rg]) for g in goci])
goci_time['doys'] = np.array([np.nanmean(g['doys'][igoci_rg]) for g in goci])


# ### Make sampling analysis and combined for easier plotting

# In[542]:


gociar_time = {'aod_mean':[],'aod_median':[],'aod_std':[],
               'ae_mean':[],'ae_median':[],'ae_std':[],'fmf_mean':[],'fmf_median':[],'fmf_std':[],'doys':[]}
ar_time = {'aod_mean':[],'aod_median':[],'aod_std':[],
               'ae_mean':[],'ae_median':[],'ae_std':[],'fmf_mean':[],'fmf_median':[],'fmf_std':[],'doys':[]}
for da in np.unique(goci2ar['doys'].astype(int)):
    ida = goci2ar['doys'].astype(int)==da
    gociar_time['aod_mean'].append(np.nanmean(goci2ar['aod'][ida]))
    gociar_time['aod_median'].append(np.nanmedian(goci2ar['aod'][ida]))
    gociar_time['aod_std'].append(np.nanstd(goci2ar['aod'][ida]))
    gociar_time['ae_mean'].append(np.nanmean(goci2ar['AE'][ida]))
    gociar_time['ae_median'].append(np.nanmedian(goci2ar['AE'][ida]))
    gociar_time['ae_std'].append(np.nanstd(goci2ar['AE'][ida]))
    gociar_time['fmf_mean'].append(np.nanmean(goci2ar['fmf'][ida]))
    gociar_time['fmf_median'].append(np.nanmedian(goci2ar['fmf'][ida]))
    gociar_time['fmf_std'].append(np.nanstd(goci2ar['fmf'][ida]))
    gociar_time['doys'].append(np.nanmean(goci2ar['doys'][ida]))
    
    iaa = (ar['doys'].astype(int)==da) & np.isfinite(ar['AOD0501']) & ar['fl_QA'] & (ar['GPS_Alt']<1000.0)
    ar_time['aod_mean'].append(np.nanmean(ar['AOD0501'][iaa]))
    ar_time['aod_median'].append(np.nanmedian(ar['AOD0501'][iaa]))
    ar_time['aod_std'].append(np.nanstd(ar['AOD0501'][iaa]))
    ar_time['ae_mean'].append(np.nanmean(angs[iaa]))
    ar_time['ae_median'].append(np.nanmedian(angs[iaa]))
    ar_time['ae_std'].append(np.nanstd(angs[iaa]))
    ar_time['fmf_mean'].append(np.nanmean(fmf['eta'][iaa]))
    ar_time['fmf_median'].append(np.nanmedian(fmf['eta'][iaa]))
    ar_time['fmf_std'].append(np.nanstd(fmf['eta'][iaa]))
    ar_time['doys'].append(np.nanmean(ar['doys'][iaa]))

for k in gociar_time.keys():
    gociar_time[k] = np.array(gociar_time[k])
    ar_time[k] = np.array(ar_time[k])


# In[543]:


goci2ar.keys()


# In[544]:


goci2ar['doys'].astype(int)


# In[117]:


fig, ax = plt.subplots(2,1,sharex=True)

ax[0].plot(goci_time['doys'],goci_time['aod_mean'],'.-',markersize=0.1,label='GOCI regional average')
ax[0].errorbar(goci_time['doys'],goci_time['aod_mean'],yerr=goci_time['aod_std'],markersize=0.1,marker='.',c='tab:blue',alpha=0.2)
ax[0].plot(gociar_time['doys'],gociar_time['aod_mean'],'.-',label='GOCI flight average')
ax[0].errorbar(gociar_time['doys'],gociar_time['aod_mean'],yerr=gociar_time['aod_std'],marker='.',c='tab:orange',alpha=0.2)
ax[0].plot(ar_time['doys'],ar_time['aod_mean'],'.-',label='4STAR average')
ax[0].errorbar(ar_time['doys'],ar_time['aod_mean'],yerr=ar_time['aod_std'],marker='.',c='tab:green',alpha=0.2)
ax[1].set_xlabel('DOY')
ax[0].set_ylabel('AOD$_{{500}}$')
ax[0].set_title('Daily average over Korea')
ax[0].legend()

ax[1].plot(goci_time['doys'],goci_time['fmf_mean'],'.-',label='GOCI regional average')
ax[1].errorbar(goci_time['doys'],goci_time['fmf_mean'],yerr=goci_time['fmf_std'],marker='.',c='tab:blue',alpha=0.2)
ax[1].plot(gociar_time['doys'],gociar_time['fmf_mean'],'.-',label='GOCI flight average')
ax[1].errorbar(gociar_time['doys'],gociar_time['fmf_mean'],yerr=gociar_time['fmf_std'],marker='.',c='tab:orange',alpha=0.2)
ax[1].plot(ar_time['doys'],ar_time['fmf_mean'],'.-',label='4STAR average')
ax[1].errorbar(ar_time['doys'],ar_time['fmf_mean'],yerr=ar_time['fmf_std'],marker='.',c='tab:green',alpha=0.2)
ax[1].set_ylabel('Fine Mode Fraction')


plt.savefig(fp+'plot/KORUS_GOCI_AOD_daily_avgs.png',dpi=600,transparent=True)


# In[142]:


plt.figure()
plt.hist([goci_time['aod_mean'],gociar_time['aod_mean'],ar_time['aod_mean']],label=['GOCI Regional','GOCI flight','4STAR'],
         range=[0,1.5],normed=True)
plt.legend()

plt.axvline(np.nanmean(goci_time['aod_mean']),ls='-',alpha=0.5,c='tab:blue')
plt.axvline(np.nanmedian(goci_time['aod_mean']),ls='--',alpha=0.5,c='tab:blue')
plt.axvline(np.nanmean(gociar_time['aod_mean']),ls='-',alpha=0.5,c='tab:orange')
plt.axvline(np.nanmedian(gociar_time['aod_mean']),ls='--',alpha=0.5,c='tab:orange')
plt.axvline(np.nanmean(ar_time['aod_mean']),ls='-',alpha=0.5,c='tab:green')
plt.axvline(np.nanmedian(ar_time['aod_mean']),ls='--',alpha=0.5,c='tab:green')
plt.xlabel('AOD$_{{500}}$')


# In[555]:


np.nanmean(goci_time['aod_mean']), np.nanmean(gociar_time['aod_mean']), np.nanmean(ar_time['aod_mean'])


# In[143]:


plt.figure()
plt.hist([goci_time['fmf_mean'],gociar_time['fmf_mean'],ar_time['fmf_mean']],label=['GOCI Regional','GOCI flight','4STAR'],
         range=[0,1.0],normed=True)
plt.legend()

plt.axvline(np.nanmean(goci_time['fmf_mean']),ls='-',alpha=0.5,c='tab:blue')
plt.axvline(np.nanmedian(goci_time['fmf_mean']),ls='--',alpha=0.5,c='tab:blue')
plt.axvline(np.nanmean(gociar_time['fmf_mean']),ls='-',alpha=0.5,c='tab:orange')
plt.axvline(np.nanmedian(gociar_time['fmf_mean']),ls='--',alpha=0.5,c='tab:orange')
plt.axvline(np.nanmean(ar_time['fmf_mean']),ls='-',alpha=0.5,c='tab:green')
plt.axvline(np.nanmedian(ar_time['fmf_mean']),ls='--',alpha=0.5,c='tab:green')
plt.xlabel('Fine Mode Fraction')


# In[556]:


np.nanmean(goci_time['fmf_mean']), np.nanmean(gociar_time['fmf_mean']), np.nanmean(ar_time['fmf_mean'])


# In[553]:


fig, ax = plt.subplots(2,1)
ax[0].hist([goci_time['aod_mean'],gociar_time['aod_mean'],ar_time['aod_mean']],
           label=['GOCI Regional','GOCI flight','4STAR'],range=[0,1.5],normed=True)


ax[0].axvline(np.nanmean(goci_time['aod_mean']),ls='-',alpha=0.5,c='tab:blue')
ax[0].axvline(np.nanmedian(goci_time['aod_mean']),ls='--',alpha=0.5,c='tab:blue')
ax[0].axvline(np.nanmean(gociar_time['aod_mean']),ls='-',alpha=0.5,c='tab:orange')
ax[0].axvline(np.nanmedian(gociar_time['aod_mean']),ls='--',alpha=0.5,c='tab:orange')
ax[0].axvline(np.nanmean(ar_time['aod_mean']),ls='-',alpha=0.5,c='tab:green')
ax[0].axvline(np.nanmedian(ar_time['aod_mean']),ls='--',alpha=0.5,c='tab:green')
ax[0].axvline([np.nan],ls='-',alpha=0.5,c='k',label='mean')
ax[0].axvline([np.nan],ls='--',alpha=0.5,c='k',label='median')
ax[0].legend()
ax[0].set_xlabel('AOD$_{{500}}$')
ax[0].set_ylabel('Normalized counts')

ax[1].hist([goci_time['fmf_mean'],gociar_time['fmf_mean'],ar_time['fmf_mean']],
           label=['GOCI Regional','GOCI flight','4STAR'],range=[0,1.0],normed=True)
ax[1].legend()

ax[1].axvline(np.nanmean(goci_time['fmf_mean']),ls='-',alpha=0.5,c='tab:blue')
ax[1].axvline(np.nanmedian(goci_time['fmf_mean']),ls='--',alpha=0.5,c='tab:blue')
ax[1].axvline(np.nanmean(gociar_time['fmf_mean']),ls='-',alpha=0.5,c='tab:orange')
ax[1].axvline(np.nanmedian(gociar_time['fmf_mean']),ls='--',alpha=0.5,c='tab:orange')
ax[1].axvline(np.nanmean(ar_time['fmf_mean']),ls='-',alpha=0.5,c='tab:green')
ax[1].axvline(np.nanmedian(ar_time['fmf_mean']),ls='--',alpha=0.5,c='tab:green')
ax[1].set_xlabel('Fine Mode Fraction')
ax[1].set_ylabel('Normalized counts')
plt.tight_layout()

plt.savefig(fp+'plot/KORUS_GOCI_hist_daily_avgs.png',dpi=600,transparent=True)


# In[84]:


angs.shape, fmf['eta'].shape


# In[85]:


ar_time['delta_aod_gociar'] = np.zeros_like(ar_time['doys'])+np.nan
gociar_time['delta_aod_goci'] = np.zeros_like(ar_time['doys'])+np.nan
gociar_time['std_aod_goci'] = np.zeros_like(ar_time['doys'])+np.nan
ar_time['delta_ae_gociar'] = np.zeros_like(ar_time['doys'])+np.nan
gociar_time['delta_ae_goci'] = np.zeros_like(ar_time['doys'])+np.nan
gociar_time['std_ae_goci'] = np.zeros_like(ar_time['doys'])+np.nan
ar_time['delta_fmf_gociar'] = np.zeros_like(ar_time['doys'])+np.nan
gociar_time['delta_fmf_goci'] = np.zeros_like(ar_time['doys'])+np.nan
gociar_time['std_fmf_goci'] = np.zeros_like(ar_time['doys'])+np.nan
for ii, da in enumerate(ar_time['doys'].astype(int)):
    igat = gociar_time['doys'].astype(int) == da
    igt  = goci_time['doys'].astype(int) == da
    ar_time['delta_aod_gociar'][ii] = ar_time['aod_mean'][ii] - np.nanmean(gociar_time['aod_mean'][igat])
    gociar_time['delta_aod_goci'][ii] = np.nanmean(gociar_time['aod_mean'][igat]) - np.nanmean(goci_time['aod_mean'][igt])
    gociar_time['std_aod_goci'][ii] = np.nanmean(goci_time['aod_std'][igt])
    ar_time['delta_ae_gociar'][ii] = ar_time['ae_mean'][ii] - np.nanmean(gociar_time['ae_mean'][igat])
    gociar_time['delta_ae_goci'][ii] = np.nanmean(gociar_time['ae_mean'][igat]) - np.nanmean(goci_time['ae_mean'][igt])
    gociar_time['std_ae_goci'][ii] = np.nanmean(goci_time['ae_std'][igt])
    ar_time['delta_fmf_gociar'][ii] = ar_time['fmf_mean'][ii] - np.nanmean(gociar_time['fmf_mean'][igat])
    gociar_time['delta_fmf_goci'][ii] = np.nanmean(gociar_time['fmf_mean'][igat]) - np.nanmean(goci_time['fmf_mean'][igt])
    gociar_time['std_fmf_goci'][ii] = np.nanmean(goci_time['fmf_std'][igt])


# In[86]:


gociar_time['aod_mean'][igat]


# In[87]:


gociar_time['aod_mean'][igat] - goci_time['aod_mean'][igt]


# In[88]:


gociar_time['aod_mean'][igat]


# In[89]:


ar_time['delta_aod_gociar']


# In[90]:


gociar_time['delta_aod_goci']


# In[91]:


goci_time['doys'][igt]


# In[92]:


np.isfinite(gociar_time['std_aod_goci'])


# In[122]:


istd = np.isfinite(gociar_time['std_aod_goci'])

plt.figure()
plt.plot(ar_time['doys'],ar_time['delta_aod_gociar'],'x',label='4STAR - GOCI flights')
iout = abs(gociar_time['delta_aod_goci'])>gociar_time['std_aod_goci']
plt.plot(ar_time['doys'][~iout],gociar_time['delta_aod_goci'][~iout],'+',label='GOCI flights - GOCI regional',mew=2)

plt.plot(ar_time['doys'][iout],gociar_time['delta_aod_goci'][iout],'+',c='tab:orange',mew=0.5,label='Sampling not representative')
plt.axhline(0,ls='--',c='grey',alpha=0.2,zorder=-10)
plt.fill_between(ar_time['doys'][istd],0.0-gociar_time['std_aod_goci'][istd],gociar_time['std_aod_goci'][istd],
                 color='grey',alpha=0.2,label='regional variation',lw=0.0)
plt.legend()
plt.xlabel('DOY')
plt.ylabel('Difference in AOD daily averages')
plt.savefig(fp+'plot/KORUS_GOCI_AOD_sampling_representative.png',dpi=600,transparent=True)


# In[131]:


plt.figure()
plt.hist(abs(gociar_time['delta_aod_goci']),range=[0,0.4],bins=15,label='GOCI')
plt.xlabel('Difference in AOD daily averages (regional to flight sampling)')


# ## Load MERRA2

# ### Old version
# 
# 
# For MERRA2 data
# Downloaded from https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NXGAS.5.12.4/doc/MERRA2.README.pdf
# https://disc.gsfc.nasa.gov/
#     
# on 2020-06-26
# 
# // global attributes:
#                 :History = "Original file generated: Thu May 19 23:35:04 2016 GMT" ;
#                 :Comment = "GMAO filename: d5124_m2_jan10.inst3_2d_gas_Nx.20160501.nc4" ;
#                 :Filename = "MERRA2_400.inst3_2d_gas_Nx.20160501.nc4" ;
#                 :Conventions = "CF-1" ;
#                 :Institution = "NASA Global Modeling and Assimilation Office" ;
#                 :References = "http://gmao.gsfc.nasa.gov" ;
#                 :Format = "NetCDF-4/HDF-5" ;
#                 :SpatialCoverage = "global" ;
#                 :VersionID = "5.12.4" ;
#                 :TemporalRange = "1980-01-01 -> 2016-12-31" ;
#                 :identifier_product_doi_authority = "http://dx.doi.org/" ;
#                 :ShortName = "M2I3NXGAS" ;
#                 :GranuleID = "MERRA2_400.inst3_2d_gas_Nx.20160501.nc4" ;
#                 :ProductionDateTime = "Original file generated: Thu May 19 23:35:04 2016 GMT" ;
#                 :LongName = "MERRA2 inst3_2d_gas_Nx: 2d,3-Hourly,Instantaneous,Single-Level,Assimilation,Aerosol Optical Depth Analysis" ;
#                 :Title = "MERRA2 inst3_2d_gas_Nx: 2d,3-Hourly,Instantaneous,Single-Level,Assimilation,Aerosol Optical Depth Analysis" ;
#                 :SouthernmostLatitude = "-90.0" ;
#                 :NorthernmostLatitude = "90.0" ;
#                 :WesternmostLongitude = "-180.0" ;
#                 :EasternmostLongitude = "179.375" ;
#                 :LatitudeResolution = "0.5" ;
#                 :LongitudeResolution = "0.625" ;
#                 :DataResolution = "0.5 x 0.625" ;
#                 :Source = "CVS tag: GEOSadas-5_12_4_p5" ;
#                 :Contact = "http://gmao.gsfc.nasa.gov" ;
#                 :identifier_product_doi = "10.5067/HNGA0EWW0R09" ;
#                 :RangeBeginningDate = "2016-05-01" ;
#                 :RangeBeginningTime = "00:00:00.000000" ;
#                 :RangeEndingDate = "2016-05-01" ;
#                 :RangeEndingTime = "21:00:00.000000" ;
#                 :DODS_EXTRA.Unlimited_Dimension = "time" ;
#                 :history = "2020-06-26 22:21:20 GMT Hyrax-1.15.1 https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NXGAS.5.12.4/2016/05/MERRA2_400.inst3_2d_gas_Nx.20160501.nc4.nc4?AODANA[0:7][243:259][484:499],time,lat[243:259],lon[484:499]" ;

# In[59]:


ml = os.listdir(fp+'data_other/MERRA2')


# In[60]:


ml = [ m for m in ml if m.endswith('nc4')] 
ml.sort()


# In[61]:


merra = []
for m in ml:
    print 'Loading file: {}'.format(m)
    mtmp, mdict = lu.load_hdf_h5py(fp+'data_other/MERRA2/'+m,all_values=True,verbose=False)
    merra.append(mtmp)


# In[62]:


len(merra),len(ml)


# In[201]:


mdict.keys()


# In[63]:


mdict['AODINC']


# In[64]:


merra[0]['time']


# In[65]:


merra[0]['time']/60.0/24.0


# In[62]:


(1-0.875)/2


# #### Collocate to 4STAR

# In[66]:


merra_doy = [datetime(int(l[-16:-12]),int(l[-12:-10]),int(l[-10:-8])).timetuple().tm_yday for l in ml]


# In[67]:


if False:
    for m in merra:
        m['aod'] = m['AODANA']+m['AODINC']
else:
    for m in merra:
        m['aod'] = m['AODANA']


# In[68]:


im4 = Sp.find_closest(np.array(merra_doy),(ar['doys']+0.0625).astype(int))


# In[69]:


np.unique(im4)


# In[70]:


nmerra = len(merra)
nlat = len(merra[0]['lat'])
nlon = len(merra[0]['lon'])
ntime = len(merra[0]['time'])


# In[71]:


merra2ar = {'doys':[],'lat':[],'lon':[],'aod':[],'ind':[]}
for ii in np.unique(im4):
    print ii
    imeas = np.where(im4==ii)[0]
    itime = Sp.find_closest(merra[ii]['time']/60.0/24.0+merra_doy[ii],ar['doys'][imeas])
    ilat = Sp.find_closest(merra[0]['lat'],ar['Latitude'][imeas])
    ilon = Sp.find_closest(merra[0]['lon'],ar['Longitude'][imeas])
    
    if len(imeas)<1:
        merra2ar['lat'] = np.append(merra2ar['lat'],ar['Latitude'][imeas])
        merra2ar['lon'] = np.append(merra2ar['lon'],ar['Longitude'][imeas])
        merra2ar['aod'] = np.append(merra2ar['aod'],ar['Latitude'][imeas]+np.nan)
        merra2ar['doys'] = np.append(merra2ar['doys'],ar['doys'][imeas])
        merra2ar['ind'] = np.append(merra2ar['ind'],ar['doys'][imeas].astype(int)*0) #pixel index
    else:
        merra2ar['lat'] = np.append(merra2ar['lat'],merra[ii]['lat'][ilat])
        merra2ar['lon'] = np.append(merra2ar['lon'],merra[ii]['lon'][ilon])
        merra2ar['aod'] = np.append(merra2ar['aod'],merra[ii]['aod'][itime,ilat,ilon])
        merra2ar['doys'] = np.append(merra2ar['doys'],merra[ii]['time'][itime]/60.0/24.0+merra_doy[ii])
        merra2ar['ind'] = np.append(merra2ar['ind'],
                                    np.ravel_multi_index((imeas*0+ii,itime,ilat,ilon),(nmerra,ntime,nlat,nlon)).astype(int))

for k in merra2ar.keys():
    merra2ar[k] = np.array(merra2ar[k])


# In[72]:


merra2ar['aod'].shape


# In[73]:


merra2ar['ind']


# In[71]:


plt.figure()
plt.plot(merra2ar['doys'],merra2ar['aod'],'.')
plt.plot(ar['doys'],ar['AOD0501'],'.')


# ### New Version with aerosol diagnostic

# From: M2T1NXAER: MERRA-2 tavg1_2d_aer_Nx: 2d,1-Hourly,Time-averaged,Single-Level,Assimilation,Aerosol Diagnostics V5.12.4
# 
# Downloaded at: https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXAER.5.12.4/2016/05/
# 
# DOI:10.5067/KLICLTZ8EM9D
# 
# 
# 
# Global Modeling and Assimilation Office (GMAO) (2015), MERRA-2 tavg1_2d_aer_Nx: 2d,1-Hourly,Time-averaged,Single-Level,Assimilation,Aerosol Diagnostics V5.12.4, Greenbelt, MD, USA, Goddard Earth Sciences Data and Information Services Center (GES DISC), Accessed: 2020-09-17, 10.5067/KLICLTZ8EM9D

# In[60]:


ml = os.listdir(fp+'data_other/MERRA2/aer')


# In[61]:


ml = [ m for m in ml if m.endswith('nc4')] 
ml.sort()


# In[62]:


mtmp, mdict = lu.load_hdf_h5py(fp+'data_other/MERRA2/aer/'+ml[0])


# In[63]:


merra = []
for m in ml:
    print 'Loading file: {}'.format(m)
    mtmp, mdict = lu.load_hdf_h5py(fp+'data_other/MERRA2/aer/'+m,
                                   values=(('lat',50),('lon',51),('time',52),('TOTEXTTAU',48),('TOTSCATAU',49),('TOTANGSTR',47),
                                           ('BCEXTTAU',2),('DUEXTTAU',13),('OCEXTTAU',22),('SSEXTTAU',35),('SUEXTTAU',43)),
                                   verbose=False)
    merra.append(mtmp)


# In[64]:


len(merra),len(ml)


# In[65]:


mdict.keys()


# In[66]:


mdict['TOTANGSTR']


# In[87]:


mdict['TOTEXTTAU']


# In[88]:


mdict['SSEXTTAU'],mdict['BCEXTTAU'],mdict['OCEXTTAU'],mdict['SUEXTTAU'],mdict['DUEXTTAU']


# In[67]:


merra[0]['time']/60.0/24.0


# In[68]:


merra[0]['time']


# In[69]:


delta_time = (1.0-merra[0]['time'][-1]/60.0/24.0)/2.0


# #### Collocate new version to 4STAR

# In[70]:


merra_doy = [datetime(int(l[-12:-8]),int(l[-8:-6]),int(l[-6:-4])).timetuple().tm_yday for l in ml]


# In[71]:


for m in merra:
    m['aod'] = m['TOTEXTTAU']
    m['ae'] = m['TOTANGSTR']
    m['aod_dust'] = m['DUEXTTAU']
    m['aod_sea'] = m['SSEXTTAU']
    m['aod_bc'] = m['BCEXTTAU']
    m['aod_oc'] = m['OCEXTTAU']
    m['aod_sulf'] = m['SUEXTTAU']


# In[72]:


im4 = Sp.find_closest(np.array(merra_doy),(ar['doys']+delta_time).astype(int))


# In[73]:


np.unique(im4)


# In[74]:


nmerra = len(merra)
nlat = len(merra[0]['lat'])
nlon = len(merra[0]['lon'])
ntime = len(merra[0]['time'])


# In[75]:


merra2ar = {'doys':[],'lat':[],'lon':[],'aod':[],'ae':[],'ind':[],
            'aod_sulf':[],'aod_dust':[],'aod_sea':[],'aod_bc':[],'aod_oc':[]}
for ii in np.unique(im4):
    print ii
    imeas = np.where(im4==ii)[0]
    itime = Sp.find_closest(merra[ii]['time']/60.0/24.0+merra_doy[ii],ar['doys'][imeas])
    ilat = Sp.find_closest(merra[0]['lat'],ar['Latitude'][imeas])
    ilon = Sp.find_closest(merra[0]['lon'],ar['Longitude'][imeas])
    
    if len(imeas)<1:
        merra2ar['lat'] = np.append(merra2ar['lat'],ar['Latitude'][imeas])
        merra2ar['lon'] = np.append(merra2ar['lon'],ar['Longitude'][imeas])
        merra2ar['aod'] = np.append(merra2ar['aod'],ar['Latitude'][imeas]+np.nan)
        merra2ar['aod_dust'] = np.append(merra2ar['aod_dust'],ar['Latitude'][imeas]+np.nan)
        merra2ar['aod_sea'] = np.append(merra2ar['aod_sea'],ar['Latitude'][imeas]+np.nan)
        merra2ar['aod_bc'] = np.append(merra2ar['aod_bc'],ar['Latitude'][imeas]+np.nan)
        merra2ar['aod_oc'] = np.append(merra2ar['aod_oc'],ar['Latitude'][imeas]+np.nan)
        merra2ar['aod_sulf'] = np.append(merra2ar['aod_sulf'],ar['Latitude'][imeas]+np.nan)
        merra2ar['ae'] = np.append(merra2ar['ae'],ar['Latitude'][imeas]+np.nan)
        merra2ar['doys'] = np.append(merra2ar['doys'],ar['doys'][imeas])
        merra2ar['ind'] = np.append(merra2ar['ind'],ar['doys'][imeas].astype(int)*0) #pixel index
    else:
        merra2ar['lat'] = np.append(merra2ar['lat'],merra[ii]['lat'][ilat])
        merra2ar['lon'] = np.append(merra2ar['lon'],merra[ii]['lon'][ilon])
        merra2ar['aod'] = np.append(merra2ar['aod'],merra[ii]['aod'][itime,ilat,ilon])
        merra2ar['aod_dust'] = np.append(merra2ar['aod_dust'],merra[ii]['aod_dust'][itime,ilat,ilon])
        merra2ar['aod_sea'] = np.append(merra2ar['aod_sea'],merra[ii]['aod_sea'][itime,ilat,ilon])
        merra2ar['aod_bc'] = np.append(merra2ar['aod_bc'],merra[ii]['aod_bc'][itime,ilat,ilon])
        merra2ar['aod_oc'] = np.append(merra2ar['aod_oc'],merra[ii]['aod_oc'][itime,ilat,ilon])
        merra2ar['aod_sulf'] = np.append(merra2ar['aod_sulf'],merra[ii]['aod_sulf'][itime,ilat,ilon])
        merra2ar['ae'] = np.append(merra2ar['ae'],merra[ii]['ae'][itime,ilat,ilon])
        merra2ar['doys'] = np.append(merra2ar['doys'],merra[ii]['time'][itime]/60.0/24.0+merra_doy[ii])
        merra2ar['ind'] = np.append(merra2ar['ind'],
                                    np.ravel_multi_index((imeas*0+ii,itime,ilat,ilon),(nmerra,ntime,nlat,nlon)).astype(int))

for k in merra2ar.keys():
    merra2ar[k] = np.array(merra2ar[k])


# In[110]:


plt.figure()
plt.plot(merra2ar['doys'],merra2ar['aod'],'.')
plt.plot(ar['doys'],ar['AOD0501'],'.')


# ### Make daily regional averages

# In[528]:


j = 2
'MERRA: lat range: {}, {}, lon range: {}, {}'.format(
    merra[j]['lat'].min(), merra[j]['lat'].max(),merra[j]['lon'].min(), merra[j]['lon'].max())


# In[529]:


lat_rg = [33.8,37.6] # from percentiles of 4STAR data
lon_rg = [124.3,129.4]


# In[530]:


imerra_rg_lat = (merra[0]['lat']>lat_rg[0]) & (merra[0]['lat']<lat_rg[1])
imerra_rg_lon = (merra[0]['lon']>lon_rg[0]) & (merra[0]['lon']<lon_rg[1])


# In[531]:


imerra_rg_lat.shape, imerra_rg_lon.shape


# In[532]:


merra_time = {}
merra_time['aod_mean'] = np.array([np.nanmean(g['aod'][:,imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
merra_time['aod_median'] = np.array([np.nanmedian(g['aod'][:,imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
merra_time['aod_std'] = np.array([np.nanstd(g['aod'][:,imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
merra_time['ae_mean'] = np.array([np.nanmean(g['ae'][:,imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
merra_time['ae_median'] = np.array([np.nanmedian(g['ae'][:,imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
merra_time['ae_std'] = np.array([np.nanstd(g['ae'][:,imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
#merra_time['fmf_mean'] = np.array([np.nanmean(g['fmf'][imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
#merra_time['fmf_median'] = np.array([np.nanmedian(g['fmf'][imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
#merra_time['fmf_std'] = np.array([np.nanstd(g['fmf'][imerra_rg_lat,:][:,:,imerra_rg_lon]) for g in merra])
merra_time['doys'] = np.array(merra_doy)+0.5


# In[533]:


merraar_time = {'aod_mean':[],'aod_median':[],'aod_std':[],
               'ae_mean':[],'ae_median':[],'ae_std':[],'fmf_mean':[],'fmf_median':[],'fmf_std':[],'doys':[]}
for da in np.unique(merra2ar['doys'].astype(int)):
    ida = merra2ar['doys'].astype(int)==da
    merraar_time['aod_mean'].append(np.nanmean(merra2ar['aod'][ida]))
    merraar_time['aod_median'].append(np.nanmedian(merra2ar['aod'][ida]))
    merraar_time['aod_std'].append(np.nanstd(merra2ar['aod'][ida]))
    merraar_time['ae_mean'].append(np.nanmean(merra2ar['ae'][ida]))
    merraar_time['ae_median'].append(np.nanmedian(merra2ar['ae'][ida]))
    merraar_time['ae_std'].append(np.nanstd(merra2ar['ae'][ida]))
    #merraar_time['fmf_mean'].append(np.nanmean(merra2ar['fmf'][ida]))
    #merraar_time['fmf_median'].append(np.nanmedian(merra2ar['fmf'][ida]))
    #merraar_time['fmf_std'].append(np.nanstd(merra2ar['fmf'][ida]))
    merraar_time['doys'].append(np.nanmean(merra2ar['doys'][ida]))

for k in merraar_time.keys():
    merraar_time[k] = np.array(merraar_time[k])


# In[540]:


fig, ax = plt.subplots(2,1,sharex=True)

ax[0].plot(merra_time['doys'],merra_time['aod_mean'],'.-',label='MERRA-2 regional average')
ax[0].errorbar(merra_time['doys'],merra_time['aod_mean'],yerr=merra_time['aod_std'],marker='.',c='tab:blue',alpha=0.2)
ax[0].plot(merraar_time['doys'],merraar_time['aod_mean'],'.-',label='MERRA-2 flight average')
ax[0].errorbar(merraar_time['doys'],merraar_time['aod_mean'],yerr=merraar_time['aod_std'],marker='.',c='tab:orange',alpha=0.2)
ax[0].plot(ar_time['doys'],ar_time['aod_mean'],'.-',label='4STAR average')
ax[0].errorbar(ar_time['doys'],ar_time['aod_mean'],yerr=ar_time['aod_std'],marker='.',c='tab:green',alpha=0.2)
ax[1].set_xlabel('DOY')
ax[0].set_ylabel('AOD$_{{500}}$')
ax[0].set_title('Daily average over Korea')
ax[0].legend()

ax[1].plot(merra_time['doys'],merra_time['ae_mean'],'.-',label='merra regional average')
ax[1].errorbar(merra_time['doys'],merra_time['ae_mean'],yerr=merra_time['ae_std'],marker='.',c='tab:blue',alpha=0.2)
ax[1].plot(merraar_time['doys'],merraar_time['ae_mean'],'.-',label='merra flight average')
ax[1].errorbar(merraar_time['doys'],merraar_time['ae_mean'],yerr=merraar_time['ae_std'],marker='.',c='tab:orange',alpha=0.2)
ax[1].plot(ar_time['doys'],ar_time['ae_mean'],'.-',label='4STAR average')
ax[1].errorbar(ar_time['doys'],ar_time['ae_mean'],yerr=ar_time['ae_std'],marker='.',c='tab:green',alpha=0.2)
ax[1].set_ylabel('Fine Mode Fraction')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_MERRA_AOD_daily_avgs.png',dpi=600,transparent=True)


# In[538]:


fig, ax = plt.subplots(2,1)
ax[0].hist([merra_time['aod_mean'],merraar_time['aod_mean'],ar_time['aod_mean']],
           label=['MERRA Regional','MERRA flight','4STAR'],range=[0,1.5],normed=True)


ax[0].axvline(np.nanmean(merra_time['aod_mean']),ls='-',alpha=0.5,c='tab:blue')
ax[0].axvline(np.nanmedian(merra_time['aod_mean']),ls='--',alpha=0.5,c='tab:blue')
ax[0].axvline(np.nanmean(merraar_time['aod_mean']),ls='-',alpha=0.5,c='tab:orange')
ax[0].axvline(np.nanmedian(merraar_time['aod_mean']),ls='--',alpha=0.5,c='tab:orange')
ax[0].axvline(np.nanmean(ar_time['aod_mean']),ls='-',alpha=0.5,c='tab:green')
ax[0].axvline(np.nanmedian(ar_time['aod_mean']),ls='--',alpha=0.5,c='tab:green')
ax[0].axvline([np.nan],ls='-',alpha=0.5,c='k',label='mean')
ax[0].axvline([np.nan],ls='--',alpha=0.5,c='k',label='median')
ax[0].legend()
ax[0].set_xlabel('AOD$_{{500}}$')
ax[0].set_ylabel('Normalized counts')

ax[1].hist([merra_time['ae_mean'],merraar_time['ae_mean'],ar_time['ae_mean']],
           label=['MERRA Regional','MERRA flight','4STAR'],range=[0,2.2],normed=True)
ax[1].legend()

ax[1].axvline(np.nanmean(merra_time['ae_mean']),ls='-',alpha=0.5,c='tab:blue')
ax[1].axvline(np.nanmedian(merra_time['ae_mean']),ls='--',alpha=0.5,c='tab:blue')
ax[1].axvline(np.nanmean(merraar_time['ae_mean']),ls='-',alpha=0.5,c='tab:orange')
ax[1].axvline(np.nanmedian(merraar_time['ae_mean']),ls='--',alpha=0.5,c='tab:orange')
ax[1].axvline(np.nanmean(ar_time['ae_mean']),ls='-',alpha=0.5,c='tab:green')
ax[1].axvline(np.nanmedian(ar_time['ae_mean']),ls='--',alpha=0.5,c='tab:green')
ax[1].set_xlabel('Angstrom Exponent')
ax[1].set_ylabel('Normalized counts')
plt.tight_layout()

plt.savefig(fp+'plot/KORUS_MERRA_hist_daily_avgs.png',dpi=600,transparent=True)


# ## Load in situ extinction

# In[103]:


lrgl = os.listdir(fp+'data_other/LARGE')


# In[104]:


lrgl.sort()


# In[105]:


lrgl


# In[106]:


large = []
for g in lrgl:
    print 'Loading file: {}'.format(g)
    gtmp, gdict = lu.load_ict(fp+'data_other/LARGE/'+g,return_header=True)
    large.append(gtmp)


# In[491]:


np.diff(gtmp['UTC_mid'])*3600.0


# In[557]:


gdict


# ### Match to 4STAR

# In[107]:


large[2]['UTC_mid'],large[3]['UTC_mid'] 


# In[108]:


ar['doys']


# In[109]:


large_doys = [datetime(int(l[26:30]),int(l[30:32]),int(l[32:34])).timetuple().tm_yday+np.nanmean(large[i]['UTC_mid'])/24.0 for i,l in enumerate(lrgl)]


# In[110]:


large_doys


# In[111]:


il4 = Sp.find_closest(np.array(large_doys),ar['doys'])


# In[112]:


lg2ar = {'ext_532':[],'SSA_550':[],'scat_450':[],'scat_550':[],'scat_700':[],'AE':[]}
for lll in np.unique(il4):
    imeas = np.where(il4==lll)[0]
    if len(imeas) < 1:
        lg2ar['ext_532'] = np.append(lg2ar['ext_532'],ar['aod'][imeas]+np.nan)
        lg2ar['SSA_550'] = np.append(lg2ar['SSA_550'],ar['aod'][imeas]+np.nan)
        lg2ar['scat_450'] = np.append(lg2ar['scat_450'],ar['aod'][imeas]+np.nan)
        lg2ar['scat_550'] = np.append(lg2ar['scat_550'],ar['aod'][imeas]+np.nan)
        lg2ar['scat_700'] = np.append(lg2ar['scat_700'],ar['aod'][imeas]+np.nan)
        lg2ar['AE'] = np.append(lg2ar['AE'],ar['aod'][imeas]+np.nan)
        
    else: 
        lg2ar['ext_532'] = np.append(lg2ar['ext_532'],
                                     wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambEXT_532_stdPT'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['SSA_550'] = np.append(lg2ar['SSA_550'], 
                                     wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambSSA_550'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['scat_450'] = np.append(lg2ar['scat_450'], 
                                      wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambSC450_stdPT'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['scat_550'] = np.append(lg2ar['scat_550'], 
                                      wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambSC550_stdPT'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['scat_700'] = np.append(lg2ar['scat_700'], 
                                      wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambSC700_stdPT'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['AE'] = np.append(lg2ar['AE'], 
                                wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['dryAEscat_550to700'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
  
    


# ## Load the MODIS reflectances

# In[73]:


fp


# In[74]:


fpd = fp+'data_other/MOD021KM.A2016140.0140.061.2017326011557.hdf'


# In[75]:


fpd


# In[9]:


mod,mod_dict = lu.load_hdf(fpd,values=(('refsb',0),('urefsb',1),('emis',2),('uemis',3),('refsb25',4),('refsw',7)))


# In[10]:


mod1,mod1_dict = lu.load_hdf(fp+'data_other/MOD02QKM.A2016140.0140.061.2017326011557.hdf',values=(('lat',2),('lon',3)))


# In[11]:


moda,moda_dict = lu.load_hdf(fp+'data_other/MOD04_L2.A2016140.0140.061.2017326145822.hdf',values=(('aod',135),('lat',71),('lon',70)))


# In[27]:


moda_dict['aod']


# In[12]:


modc,modc_dict = lu.load_hdf(fp+'data_other/MOD06_L2.A2016140.0140.061.2017326153640.hdf',values=(('cod',72),('cir',106)))


# In[11]:


mod1_dict['lat'].keys()


# In[13]:


refl_off = 316.9721985
refl_scal5 = 2.380933438e-06
refl_scal10 = 3.092909765e-06
refl_scal2 = 4.241572969e-05


# In[14]:


ref660 = (mod['refsb'][5,:,:]+refl_off)*refl_scal5
ref860 = (mod['refsb'][10,:,:]+refl_off)*refl_scal10
ref1240 = (mod['refsw'][2,:,:])*refl_scal2


# In[15]:


mod['refsw'].shape


# In[16]:


mod_dict['refsb']['band_names']


# ### Plot some MODIS data

# In[18]:


def make_map1(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=128,lat_0=36.0,
            llcrnrlon=125.0, llcrnrlat=33.0,
            urcrnrlon=130.0, urcrnrlat=38,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(123,133,11),labels=[0,0,0,1])
    m.drawparallels(np.linspace(31,39,17),labels=[1,0,0,0])
    return m


# In[18]:


fig,ax = plt.subplots(1,3)
m1 = make_map1(ax[2])
m1.pcolor(mod1['lon'],mod1['lat'],mod['refsw'][2,:,:],latlon=True)
ax[2].set_title('1240 nm')


m2 = make_map1(ax[0])
m2.pcolor(mod1['lon'],mod1['lat'],mod['refsb'][5,:,:],latlon=True)
ax[0].set_title('660 nm')

m3 = make_map1(ax[1])
m3.pcolor(mod1['lon'],mod1['lat'],mod['refsb'][10,:,:],latlon=True)
ax[1].set_title('860 nm')


# In[33]:


fig,ax = plt.subplots(1,1)
m0 = make_map1(ax)
mta = m0.pcolor(moda['lon'],moda['lat'],moda['aod'],latlon=True,cmap=plt.cm.YlOrRd,vmin=0,vmax=0.8)
mtc = m0.pcolor(mod1['lon'],mod1['lat'],modc['cod'],latlon=True,cmap=plt.cm.gist_earth,vmin=0,vmax=20.0)



cbxr = plt.gcf().add_axes([0.83, 0.2, 0.02, 0.6])
cbxe = plt.gcf().add_axes([0.80, 0.2, 0.02, 0.6])
#cbg.set_ticks([0,6,12,16,18,20])
#cbb.set_ticks([0,6,12,16,18,20]),cbb.set_ticklabels(['','','','',''])
cbxe.xaxis.set_ticks_position('top')#,cbaxesbl.yaxis.set_ticks_position('left')
#cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

bxr = plt.colorbar(mta,label='AOD',cax=cbxr,orientation='vertical',shrink=0.5)
bxe = plt.colorbar(mtc,cax=cbxe,orientation='vertical',shrink=0.5)
cbxe.yaxis.set_ticks_position('left')
#cbxe.set_label(' ')
cbxe.text(-4.0,0.5,'COD',rotation=90,verticalalignment='center')
plt.savefig(fp+'plot/KORUS_MODIS_AOD_COD_20160519.png',dpi=600,transparent=True)


# In[ ]:


fig,ax = plt.subplots(1,4)
m0 = make_map1(ax[0])
mta = m0.pcolor(moda['lon'],moda['lat'],moda['aod'],latlon=True,cmap=plt.cm.plasma)
mtc = m0.pcolor(mod1['lon'],mod1['lat'],modc['cod'],latlon=True,cmap=plt.cm.gist_earth)

m1 = make_map1(ax[3])
mp = m1.pcolor(mod1['lon'],mod1['lat'],ref1240,latlon=True,vmin=0,vmax=0.5)
ax[2].set_title('1240 nm')

m2 = make_map1(ax[1])
m2.pcolor(mod1['lon'],mod1['lat'],ref660,latlon=True,vmin=0,vmax=0.5)
ax[0].set_title('660 nm')

m3 = make_map1(ax[2])
m3.pcolor(mod1['lon'],mod1['lat'],ref860,latlon=True,vmin=0,vmax=0.5)
ax[1].set_title('860 nm')
plt.colorbar(mp,label='Reflectance',ax=ax[1:],shrink=0.6,location='bottom')
#plt.tight_layout()
#plt.savefig(fp+'plot/KORUS_MODIS_660_860_1240_reflect_20160519.png',dpi=600,transparent=True)


# In[38]:


fig,ax = plt.subplots(1,1)
m2 = make_map1(ax)
mp = m2.pcolor(mod1['lon'],mod1['lat'],ref660/ref860,latlon=True,vmin=0,vmax=4)
ax.set_title('660 nm / 860 nm')
plt.colorbar(mp)#label='Radiance ratio 660 nm / 860 nm')


#m3 = make_map1(ax[1])
#m3.pcolor(mod1['lon'],mod1['lat'],mod['refsb'][10,:,:],latlon=True)
#ax[1].set_title('860 nm')


# In[39]:


ref660/ref860


# In[23]:


mod['refsb'][5,:,:]/mod['refsb'][10,:,:]


# ## Load OCO-2 data

# In[59]:


oco, oco_dict = lu.load_hdf(fp+'data_other/acos_L2s_160511_06_B9200_PolB_190817214259.h5')


# In[60]:


lu.load_hdf_h5py(fp+'data_other/acos_L2s_160511_06_B9200_PolB_190817214259.h5')


# In[65]:


import h5py


# In[66]:


foco = h5py.File(fp+'data_other/acos_L2s_160511_06_B9200_PolB_190817214259.h5','r')


# In[81]:


ko = foco.keys()


# In[85]:


for k in ko:
    jk = foco[k].keys()
    print [(k,i) for i in jk if 'lat' in i]
    


# In[68]:


foco['AlbedoResults'].keys()


# In[86]:


foco['Shapes'].keys()


# In[114]:


foco['Shapes']['InputPtr_Array'].attrs.items()


# In[97]:


ra = foco['Shapes']['Retrieval_Array']


# In[102]:


ra.attrs.items()


# In[76]:


foco['AlbedoResults']['albedo_o2_fph'].attrs.items()


# In[120]:


aoo = foco['AlbedoResults']['albedo_o2_fph']


# In[122]:


aoo.value.shape


# In[125]:


foco['SoundingGeometry']['sounding_latitude'].value


# In[132]:


foco['RetrievalResults']['retrieved_o2_column'].value


# In[129]:


plt.figure()
plt.scatter(foco['SoundingGeometry']['sounding_longitude'].value,
          foco['SoundingGeometry']['sounding_latitude'].value,c=foco['AlbedoResults']['albedo_o2_fph'].value)


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# In[76]:


fl1 = ar['days']==ar['days'][0]


# In[77]:


fl1.shape


# In[78]:


fl = (ar['qual_flag']==0) & (np.isfinite(ar['AOD0501'])) 


# In[79]:


fl1.shape


# ## Calculate the Angstrom Exponent

# In[80]:


nwl,nm


# In[81]:


aodrr = np.array([ar[n] for n in nwl])


# In[82]:


aodrr.shape


# In[83]:


angs = su.calc_angs(ar['Start_UTC'],np.array(nm[1:11]),aodrr[1:11,:])


# ## Calculate the fine mode fraction

# In[84]:


fmf = su.sda(aodrr[1:13,:],np.array(nm[1:13])/1000.0)


# In[85]:


fmf.keys()


# In[86]:


fmf['tauc'].shape, ar['GPS_Alt'].shape


# ## Subset the level legs

# In[87]:


def running_std(x,n):
    'Function to do a running standard deviation on array (x) with window size (n)'
    q = x**2
    q = np.convolve(q, np.ones((n, )), mode="same")
    s = np.convolve(x, np.ones((n, )), mode="same")
    o = (q-s**2/n)/float(n-1)
    return o 


# In[88]:


nbox = 20


# In[89]:


std_alt = running_std(ar['GPS_Alt'][fl],nbox)


# In[90]:


std_alt.shape


# In[91]:


ar['GPS_Alt'][fl].shape


# In[92]:


f_level = np.where(std_alt<5.0)[0]


# In[93]:


std_alt1 = running_std(ar['GPS_Alt'][fl1],nbox)


# In[94]:


f_level1 = np.where(std_alt1<5.0)[0]


# In[95]:


ar['Start_UTC'][fl1][f_level1]


# In[149]:


plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(ar['Start_UTC'][fl1],ar['GPS_Alt'][fl1],'.')
plt.plot(ar['Start_UTC'][fl1][f_level1],ar['GPS_Alt'][fl1][f_level1],'r.')


ax2 = plt.subplot(2,1,2,sharex=ax1)
plt.plot(ar['Start_UTC'][fl1],std_alt1,'.')
plt.plot(ar['Start_UTC'][fl1][f_level1],std_alt[f_level1],'r.')
plt.ylim(0,100)


# In[118]:


plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(ar['Start_UTC'][fl],ar['GPS_Alt'][fl],'.')
plt.plot(ar['Start_UTC'][fl][f_level],ar['GPS_Alt'][fl][f_level],'r.')


ax2 = plt.subplot(2,1,2,sharex=ax1)
plt.plot(ar['Start_UTC'][fl],std_alt,'.')
plt.plot(ar['Start_UTC'][fl][f_level],std_alt[[f_level]],'r.')
plt.ylim(0,100)


# ## Seperate each of the level legs into distinct segments

# In[96]:


def get_segments(index,vals_dict,nsep=150,set_nan=True):
    'Function to seperate continuous segments (within a distance of nsep) based on a prior index'
    disc_flacaod_long = np.where(np.diff(index,1)>nsep)[0]
    
    discontinuity_istart_long =  index[np.append(0,disc_flacaod_long[:-1]+1)]
    discontinuity_iend_long =  index[disc_flacaod_long]
    
    kv = vals_dict.keys()
    d = {k:[] for k in kv}
    for i,start in enumerate(discontinuity_istart_long): # loop through discontinuities 
        if discontinuity_iend_long[i]-start < 2: continue
        for k in kv: # loop through keys
            try:
                d[k].append(vals_dict[k][start:discontinuity_iend_long[i]])
            except:
                print start, discontinuity_iend_long[i]
                continue
                #d[k].append([np.nan])
    
    for k in kv:
        d[k] = np.array(d[k])
        
    return d


# In[97]:


np.unique(ar['days'])


# In[99]:


def get_segments_by_time(index,doys,vals_dict,tsep=5.0/24.0/60.0/60.0,set_nan=True):
    'Function to seperate continuous segments (within a distance in doys of tsep) based on a prior index (default for 5 seconds in doy)'
    disc_flacaod_long = np.where(np.diff(doys[index],1)>tsep)[0]
    
    discontinuity_istart_long =  index[np.append(0,disc_flacaod_long[:-1]+1)]
    discontinuity_iend_long =  index[disc_flacaod_long]
    
    kv = vals_dict.keys()
    d = {k:[] for k in kv}
    for i,start in enumerate(discontinuity_istart_long): # loop through discontinuities 
        if discontinuity_iend_long[i]-start < 2: continue
        for k in kv: # loop through keys
            try:
                d[k].append(vals_dict[k][start:discontinuity_iend_long[i]])
            except:
                print start, discontinuity_iend_long[i]
                continue
                #d[k].append([np.nan])
    
    for k in kv:
        d[k] = np.array(d[k])
        
    return d


# In[101]:


angs[angs>5.0] = np.nan


# In[113]:


lg2ar.keys()


# In[114]:


vals = {'utc':ar['Start_UTC'][fl],'alt':ar['GPS_Alt'][fl],'lat':ar['Latitude'][fl],'lon':ar['Longitude'][fl],
        'aod0500':ar['AOD0501'][fl],'aod1040':ar['AOD1040'][fl],'AE':angs[fl],'doys':ar['doys'][fl],
        'aod_fine':fmf['tauf'][fl],'aod_coarse':fmf['tauc'][fl],'fmf':fmf['eta'][fl],
        'GOCI_AOD':goci2ar['aod'][fl],'GOCI_AE':goci2ar['AE'][fl],'GOCI_fmf':goci2ar['fmf'][fl],
        'GOCI_AOD_f':goci2ar['aod_f'][fl],'GOCI_AOD_c':goci2ar['aod_c'][fl],
        'MERRA_AOD':merra2ar['aod'][fl], 'MERRA_AE':merra2ar['ae'][fl],
        'MERRA_AOD_dust':merra2ar['aod_dust'][fl], 'MERRA_AOD_sea':merra2ar['aod_sea'][fl], 'MERRA_AOD_bc':merra2ar['aod_bc'][fl],
        'MERRA_AOD_oc':merra2ar['aod_oc'][fl], 'MERRA_AOD_sulf':merra2ar['aod_sulf'][fl],
        'situ_ext':lg2ar['ext_532'][fl],'situ_ssa':lg2ar['SSA_550'][fl],'situ_ae':lg2ar['AE'][fl]}


# In[115]:


vals.keys()


# In[116]:


dvals = get_segments(f_level,vals,nsep=100)


# In[117]:


dvalst = get_segments_by_time(f_level,vals['doys'],vals,tsep=200.0/24.0/60.0/60.0)


# In[118]:


dvals.keys()


# In[119]:


len(dvals['AE']),len(dvalst['AE'])


# In[120]:


dvals['doys'][0]


# In[121]:


dvalst['doys'][0]


# In[122]:


dvalst['doys'][1]


# In[123]:


dvals['len_minutes'] = []
dvals['len_doys'] = []
for i,n in enumerate(dvals['utc']):
    try:
        len_min = (n[-1]-n[0])*60.0
        if len_min < 0.0: 
            len_min = (n[-1]-(n[0]-24.0))*60.0
        print i,len_min
        dvals['len_minutes'].append(len_min)
        dvals['len_doys'].append(dvals['doys'][i][-1]-dvals['doys'][i][0])
    except:
        print np.nan
        dvals['len_minutes'].append(np.nan)
        dvals['len_doys'].append(np.nan)


# In[124]:


dvalst['len_minutes'] = []
dvalst['len_doys'] = []
for i,n in enumerate(dvalst['utc']):
    try:
        len_min = (n[-1]-n[0])*60.0
        print i,len_min
        dvalst['len_minutes'].append(len_min)
        dvalst['len_doys'].append(dvalst['doys'][i][-1]-dvalst['doys'][i][0])
    except:
        print np.nan
        dvalst['len_minutes'].append(np.nan)
        dvalst['len_doys'].append(np.nan)


# In[209]:


plt.figure()
plt.hist([dvals['len_doys'],dvalst['len_doys']],bins=50,label=['number based','time based'])
plt.legend()
plt.yscale('log')
plt.ylabel('Number of segments')
plt.xlabel('Fractional day length of segment')


# In[125]:


dvals = dvalst # set the default to use the time seperated segments


# ### discrete colorbar

# In[188]:


def discrete_matshow(data,cmapname='RdBu'):
    ' plotting function for a discrete colormap'
    cmap = plt.get_cmap(cmapname, np.nanmax(data)-np.nanmin(data)+1)
    # set limits .5 outside true range
    scalarmap = plt.cm.ScalarMappable(cmap=cmapname)
    scalarmap.set_array(data)
    #mat = plt.matshow(data,cmap=cmap,vmin = np.min(data)-.5, vmax = np.max(data)+.5)
    #tell the colorbar to tick at integers
    cax = plt.colorbar(scalarmap, ticks=np.arange(np.min(data),np.max(data)+1))
    return cax


# In[ ]:


for q in np.unique(ar['days']):
    flq = ar['days'][fl]==q
    flql = ar['days'][fl][f_level]==q
    plt.figure()
    plt.plot(ar['Start_UTC'][fl][flq],ar['GPS_Alt'][fl][flq],'.')
    plt.plot(ar['Start_UTC'][fl][f_level][flql],ar['GPS_Alt'][fl][f_level][flql],'r.')
    ax = plt.gca()

    ax.set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, len(dvals['utc'])+1)])

    for i,n in enumerate(dvals['utc']):
        plt.plot(n,dvals['alt'][i],'s-',markeredgecolor='None')

    plt.xlabel('UTC [h from midnight]')
    plt.ylabel('Altitude [m]')
    plt.title('Days: {}'.format(q))

    #scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('gist_ncar'))
    #scalarmap.set_array(range(len(dvals['utc'])+1))
    #cb = plt.colorbar(scalarmap)
    cb = discrete_matshow(range(len(dvals['utc'])+1),cmapname='gist_ncar')
    cb.set_label('Level leg number')
#plt.plot(ar['Start_UTC'][fl1][f_level],ar['GPS_Alt'][fl1][f_level],'r.')


# ## Now calculate the distances travelled within each segments

# In[126]:


def get_distances(seg_dict):
    'Function that calculates the cumulative distance and instantaneous change between each point for a set of segments'
    seg_dict['dist'],seg_dict['cumdist'] = [],[]
    for i,l in enumerate(seg_dict['lat']):
        try:
            ckm,km = [],[]
            pos1 = [seg_dict['lat'][i][0],seg_dict['lon'][i][0]] 
            for j,a in enumerate(seg_dict['lat'][i]):
                d = mu.spherical_dist(pos1,[seg_dict['lat'][i][j],seg_dict['lon'][i][j]])
                ckm.append(abs(d))
                km.append(abs(d))
        except:
            cckm,dkm = [np.nan],[np.nan]

        iu = np.where(np.isfinite(ckm))[0]
        try:
            fckm = interpolate.interp1d(seg_dict['utc'][i][iu],np.array(ckm)[iu])
            fkm = interpolate.interp1d(seg_dict['utc'][i][iu],np.array(km)[iu])
            cckm = fckm(seg_dict['utc'][i])
            dkm = fkm(seg_dict['utc'][i])
            seg_dict['cumdist'].append(np.array(cckm))
            seg_dict['dist'].append(np.array(np.diff(dkm)))
        except:
            seg_dict['cumdist'].append(np.array(np.nan))
            seg_dict['dist'].append(np.array(np.nan))

    return seg_dict


# In[127]:


def sort_by_cumdist(dd):
    'function to sort all the values in the dvals dict to be strictly increasing by cumulative distance cumdist'
    ke = dd.keys()
    for i,c in enumerate(dd['cumdist']):
        ic = np.argsort(c)
        for k in ke:
            try:
                if len(dd[k][i])==len(ic):
                    dd[k][i] = dd[k][i][ic]
            except TypeError:
                continue
        #recalculate the dist
        dd['dist'][i] = np.diff(c[ic])
    return dd


# In[128]:


ddv = get_distances(dvals)


# In[129]:


sort_by_cumdist(dvals)


# In[130]:


dvals['dist'][0]


# In[131]:


dvals['cumdist'][0:2]


# In[132]:


len(dvals['cumdist'])


# ## Save to file for easier reloading

# In[156]:


fp


# In[157]:


vv = 'v2'


# In[158]:


hs.savemat(fp+'KORUS_autocorr_dvals_{}.mat'.format(vv),dvals)


# ## Load from file

# In[ ]:


dvals = hs.loadmat(fp+'KORUS_autocorr_dvals.mat')


# In[ ]:





# In[575]:


dvals.keys()


# In[9]:


len(dvals['dist'])


# In[10]:


dvals['dist'][0].shape


# # Calculate the autocorrelation of AOD with respect to distance

# **From Shinozuka and Redemann, 2011, Horizontal variability of aerosol optical depth observed during the ARCTAS airborne experiment, ACP**
# 
# Autocorrelation is the correlation coefficient among all
# data pairs xj and xj+k that exist at a separation, or lag, of k. That is,
# 
# ![image.png](attachment:image.png)
# 
# where k indicates the spatial lag (or distance), m+k and std+k denote the mean and standard deviation, respectively, of all data points that are located a distance of +k away from an- other data point, and m−k and std−k are the corresponding quantities for data points located a distance of −k away from another data point (Redemann et al., 2006; Anderson et al., 2003).
# Figure 1c shows pairs of 499nm AOD measured 20km (±0.2 km) away from each other in the Canada phase. The correlation coefficient, r, is 0.37. This is the autocorrelation for 20km.

# Define the different time periods from Met. From: 
#     Peterson, D.A., Hyer, E.J., Han, S.-O., Crawford, J.H., Park, R.J., Holz, R., Kuehn, R.E., Eloranta, E., Knote, C., Jordan, C.E. and Lefer, B.L., 2019. Meteorology influencing springtime air quality, pollution transport, and visibility in Korea. Elem Sci Anth, 7(1), p.57. DOI: http://doi.org/10.1525/elementa.395
# 
#     - Dynamic meteorology and complex aerosol vertical profiles (01–16 May);
#     - Stagnation under a persistent anticyclone (17–22 May);
#     - Dynamic meteorology, low-level transport, and haze development (25–31 May); (Extreme pollution)
#     - Blocking pattern (01–07 June).
# 
# 
# ![image.png](attachment:image.png)

# The distribution and source of pm 2.5 during different met periods
# 
# ![image.png](attachment:image.png)

# ### Build the limits of the autocorrelation

# In[133]:


# for met times
dt1 = ['20160501','20160516']
dt2 = ['20160517','20160522']
dt3 = ['20160523','20160531']
dt4 = ['20160601','20160607']


# In[134]:


t1 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt1]
t2 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt2]
t3 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt3]
t4 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt4]


# In[135]:


# limits of DOY for each of the met times
t1,t2,t3,t4


# In[136]:


#altitude limits in m
z1 = [0.0, 1000.0] 
z2 = [1000.0, 3000.0]
z3 = [3000.0, 15000.0]


# ### Test out Shinozuka & Redemann autocorrelation 

# In[137]:


dvals['cumdist'][2]


# In[138]:


dvals.keys()


# In[178]:


# method for making one autocorrelation per segment. Old not recommended
if False:
    cr = []
    for i,cd in enumerate(dvals['cumdist']):
        #cd = dvals['cumdist'][i]
        corr = {'aod1040':[],'aod0500':[],'AE':[]}
        corr_ks =[0.1,0.25,0.5,0.75,1.0,1.5,2.0,3.0,5.0,7.5,10.0,12.5,15.0,20.0,
                  25.0,30.0,35.0,40.0,50.0,60.0,75.0,100.0,150.0,200.0] 
        for ik, k in enumerate(corr_ks):

        #k = 5.0 # for 5km distance
            if k>np.nanmax(cd):
                [corr[val].append(np.nan) for val in corr.keys()]
                continue
            ipk = np.argmin(abs(cd-k)) #ipk:
            imk = np.argmin(abs(cd-(cd[-1]-k))) #0:imk
            N = len(cd)
            #c = np.sqrt(2.0/(N-1))*math.gamma(N/2.0)/math.gamma((N-1.0)/2.0)

            for val in dvals.keys():
                if val in ['lon','utc','lat','cumdist','cdist_n','dist','alt','autocor','aod1040_r','AE_r','aod_n']: continue
                #print val, len(dvals[val][i])
                mpk = np.nanmean(dvals[val][i][ipk:]) #mean +k
                mmk = np.nanmean(dvals[val][i][0:imk]) #mean -k
                spk = np.nanstd(dvals[val][i][ipk:]) #std +k
                smk = np.nanstd(dvals[val][i][0:imk]) #std -k
                top = [(dvals[val][i][j]-mpk)*(dvals[val][i][j+ipk]-mmk) for j in xrange(N-ipk-1)]
                #dvals[val+'_r'] = []
                corr[val].append(np.sum(top)/((N-1)*spk*smk))
                if (corr[val][-1]>1.0) | (corr[val][-1]<0.0):
                    print '{} has bad corr: {:2.2f} val for key {}: std+k:{:2.2f}, std-k:{:2.2f}, m+k:{:2.2f}, m-k:{:2.2f} '.format(i,
                        corr[val][-1],val,spk,smk,mpk,mmk)

        for val in corr.keys():
            corr[val] = np.array(corr[val])
        cr.append(corr)


# In[139]:


len(cr)


# In[140]:


min_num = 100 # minimum number of points to be considered valid


# In[141]:


types = ['all','t1','t2','t3','t4','z1','z2','z3']


# ### Alternate calculations (lower mem)

# In[142]:


min_num = 100 # minimum number of points to be considered valid
types = ['all','t1','t2','t3','t4','z1','z2','z3']
corr_ks = [0.08,0.12, 0.18,0.27,0.4,0.6,0.9,1.35,2.0,3.0,5.0,7.5, 10.0, 15.0, 25.0, 35.0, 65.0, 100.0, 160.0,250.0, 380.0]

#[0.08,0.12,0.25,0.5,1.0,1.5,3.0,5.0,10.0,15.0,25.0,
#          35.0,65.0,100.0,150.0,200.0,350.0] 


# In[ ]:


if False:
    dv = 0.20
    for i,cd in enumerate(dvals['cumdist']):
        mat_dist = np.array([abs(cd-d) for d in cd])
        sys.stdout.write( '*{}*'.format(i))
        for ik,k in enumerate(corr_ks):
            #sys.stdout.write("{} ".format(ik))
            iimg,iipg = np.where((mat_dist>(k*(1.0-dv)))&(mat_dist<(k*(1.0+dv))))
            if len(iimg)==0:
                continue
            N = len(cd)
            for val in corr_vals:
                #all
                corr_alln[0][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[0][1][ik][val].append(dvals[val][i][iipg])
                #type 1
                if (dvals['doys'][i][0]> t1[0]) & (dvals['doys'][i][0]< t1[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_alln[1][0][ik][val].append(dvals[val][i][iimg])
                    corr_alln[1][1][ik][val].append(dvals[val][i][iipg])
                #type 2
                if (dvals['doys'][i][0]> t2[0]) & (dvals['doys'][i][0]< t2[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_alln[2][0][ik][val].append(dvals[val][i][iimg])
                    corr_alln[2][1][ik][val].append(dvals[val][i][iipg])
                #type 3
                if (dvals['doys'][i][0]> t3[0]) & (dvals['doys'][i][0]< t3[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_alln[3][0][ik][val].append(dvals[val][i][iimg])
                    corr_alln[3][1][ik][val].append(dvals[val][i][iipg])
                #type 4
                if (dvals['doys'][i][0]> t4[0]) & (dvals['doys'][i][0]< t4[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_alln[4][0][ik][val].append(dvals[val][i][iimg])
                    corr_alln[4][1][ik][val].append(dvals[val][i][iipg])
                #type 5
                if (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z1[1]):
                    corr_alln[5][0][ik][val].append(dvals[val][i][iimg])
                    corr_alln[5][1][ik][val].append(dvals[val][i][iipg])
                #type 6
                if (dvals['alt'][i][0]> z2[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_alln[6][0][ik][val].append(dvals[val][i][iimg])
                    corr_alln[6][1][ik][val].append(dvals[val][i][iipg])
                #type 7
                if (dvals['alt'][i][0]> z3[0]) & (dvals['alt'][i][0]< z3[1]):
                    corr_alln[7][0][ik][val].append(dvals[val][i][iimg])
                    corr_alln[7][1][ik][val].append(dvals[val][i][iipg])


# In[143]:


dvals['doys'][0][0], dvals['alt'][0][0]


# In[144]:


dv = 0.20
mat_dist = [None]*len(dvals['cumdist'])
iis = []
for i,cd in enumerate(dvals['cumdist']):
    mat_dist[i] = np.array([cd-d for d in cd])
    iis.append([np.where((mat_dist[i]>(k*(1.0-dv)))&(mat_dist[i]<(k*(1.0+dv)))) for ik,k in enumerate(corr_ks)])


# In[145]:


N = [len(cd) for i,cd in enumerate(dvals['cumdist'])]


# In[146]:


dvals.keys()


# In[147]:


alldist = np.hstack([dist for dist in dvals['dist']])


# In[23]:


plt.figure()
plt.hist(alldist,bins=150,range=[0,20])
plt.gca().set_yscale('log')
plt.xlabel('Distance from previous point [km]')
plt.ylabel('Numbers of point')


# In[148]:


dvals['dist'][0]


# In[149]:


mat_dist[0]


# In[150]:


iis[0]


# In[151]:


dvals['doys_n'] = np.array([dvals['doys'][i][0] for i,cd in enumerate(dvals['cumdist'])])
dvals['alt_n'] = np.array([dvals['alt'][i][0] for i,cd in enumerate(dvals['cumdist'])])


# In[384]:


itypes = [None]*8
#all data under alt <3km
itypes[0], = np.where((dvals['doys_n']>0.0) &
                        (dvals['alt_n'] > z1[0]) & (dvals['alt_n'] < z2[1]))
#time types for alt <3km
itypes[1], = np.where((dvals['doys_n'] > t1[0]) & (dvals['doys_n'] < t1[1]) &
                        (dvals['alt_n'] > z1[0]) & (dvals['alt_n'] < z2[1]))
itypes[2], = np.where((dvals['doys_n'] > t2[0]) & (dvals['doys_n'] < t2[1]) &
                        (dvals['alt_n'] > z1[0]) & (dvals['alt_n'] < z2[1]))
itypes[3], = np.where((dvals['doys_n'] > t3[0]) & (dvals['doys_n'] < t3[1]) &
                        (dvals['alt_n'] > z1[0]) & (dvals['alt_n'] < z2[1]))
itypes[4], = np.where((dvals['doys_n'] > t4[0]) & (dvals['doys_n'] < t4[1]) &
                        (dvals['alt_n'] > z1[0]) & (dvals['alt_n'] < z2[1]))
#alt types
itypes[5], = np.where((dvals['alt_n'] > z1[0]) & (dvals['alt_n'] < z1[1]))
itypes[6], = np.where((dvals['alt_n'] > z2[0]) & (dvals['alt_n'] < z2[1]))
itypes[7], = np.where((dvals['alt_n'] > z3[0]) & (dvals['alt_n'] < z3[1]))


# In[385]:


z1,z2


# In[386]:


np.array(corr_ks)*(1.0-dv),np.array(corr_ks)*(1.0+dv),corr_ks


# In[396]:


# Spot check the distance indices represent the actual distance values in the each segment
j = 0
for ik,k in enumerate(corr_ks): 
    print k, dvals['cumdist'][j][iis[j][ik][0]]-dvals['cumdist'][j][iis[j][ik][1]]


# In[397]:


val = 'AE'
ik = 1
j = 0
corrp = np.hstack([dvals[val][i][iis[i][ik][1]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
corrm = np.hstack([dvals[val][i][iis[i][ik][0]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])


# In[398]:


dvals[val][0][3]


# In[399]:


iis[0][ik][1]


# In[390]:


corrp


# In[391]:


len(corrp),len(corrm)


# In[395]:


corr_vals = ['aod1040','aod0500','AE','aod_fine','aod_coarse','fmf',
             'GOCI_AOD','GOCI_AE','GOCI_fmf',
             'MERRA_AOD','MERRA_AE',
             'situ_ext','situ_ae','situ_ssa']


# In[393]:


range_vals = {'aod1040':[0.0,3.0],'aod0500':[0.0,5.0],'AE':[-1.0,5.0],'aod_fine':[0.0,5.0],'aod_coarse':[0.0,5.0],
              'fmf':[0.0,1.0],'GOCI_AOD':[0.0,5.0],'GOCI_AE':[-1.0,5.0],'GOCI_fmf':[0.0,1.0],
              'MERRA_AOD':[0.0,5.0],'MERRA_AE':[-1.0,5.0],'situ_ext':[0.0,1000.0],'situ_ae':[-1.0,5.0],'situ_ssa':[0.2,1.0]}


# In[394]:


corr_vals = ['AE','aod1040','fmf','aod0500']


# In[163]:


len(dvals['cumdist'])


# In[242]:


autocorr = {}
autocorr_len = {}

for val in corr_vals:
    autocorr[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_len[val] = np.zeros((len(types),len(corr_ks)))
    print val
    for ik,k in enumerate(corr_ks):
        sys.stdout.write("{} ".format(ik))
        for j,jt in enumerate(types):
            corrp = np.hstack([dvals[val][i][iis[i][ik][1]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
            corrm = np.hstack([dvals[val][i][iis[i][ik][0]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
            badc = (corrp>range_vals[val][1]) | (corrp<range_vals[val][0])
            if any(badc): 
                corrp[badc] = np.nan
                corrm[badc] = np.nan
            
            badc = (corrm>range_vals[val][1]) | (corrm<range_vals[val][0])
            if any(badc): 
                corrp[badc] = np.nan
                corrm[badc] = np.nan
            
            mmk = np.nanmean(corrm)
            mpk = np.nanmean(corrp)
            smk = np.nanstd(corrm)
            spk = np.nanstd(corrp)
            top = (corrm-mpk)*(corrp-mmk) #[(v-mpk)*(corrp[iv]-mmk) for iv,v in enumerate(corrm)]
            #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
            autocorr[val][j,ik] = np.nansum(top)/((len(corrm)-1)*spk*smk)
            autocorr_len[val][j,ik] = len(corrm)
            if autocorr_len[val][j,ik]<min_num:
                autocorr[val][j,ik] = np.nan


# In[194]:


ix


# In[196]:


# For calculating monte carlo over the selection of horizontal legs, and not the samples, first testing
autocorr_mc = {}
autocorr_len_mc = {}
subsamp_ratio = 0.30
N_mc = 50


for val in ['AE']:
    autocorr_mc[val] = np.zeros((len(types),len(corr_ks),N_mc))+np.nan
    autocorr_len_mc[val] = np.zeros((len(types),len(corr_ks),N_mc))
    for ix in range(N_mc):
        print ix, val
        irand = np.random.randint(len(dvals['cumdist']),size=int(subsamp_ratio*len(dvals['cumdist'])))
        for ik,k in enumerate(corr_ks):
            sys.stdout.write("{} ".format(ik))
            for j,jt in enumerate(types):        
                corrp = np.hstack([dvals[val][i][iis[i][ik][1]] for i in irand if i in itypes[j]])
                corrm = np.hstack([dvals[val][i][iis[i][ik][0]] for i in irand if i in itypes[j]])
                badc = (corrp>range_vals[val][1]) | (corrp<range_vals[val][0])
                if any(badc): 
                    corrp[badc] = np.nan
                    corrm[badc] = np.nan

                badc = (corrm>range_vals[val][1]) | (corrm<range_vals[val][0])
                if any(badc): 
                    corrp[badc] = np.nan
                    corrm[badc] = np.nan

                mmk = np.nanmean(corrm)
                mpk = np.nanmean(corrp)
                smk = np.nanstd(corrm)
                spk = np.nanstd(corrp)
                top = (corrm-mpk)*(corrp-mmk) #[(v-mpk)*(corrp[iv]-mmk) for iv,v in enumerate(corrm)]
                #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
                autocorr_mc[val][j,ik,ix] = np.nansum(top)/((len(corrm)-1)*spk*smk)
                autocorr_len_mc[val][j,ik,ix] = len(corrm)
                if autocorr_len_mc[val][j,ik,ix]<min_num:
                    autocorr_mc[val][j,ik,ix] = np.nan


# In[197]:


autocorr_mc['AE'].shape


# In[198]:


plt.figure()
plt.errorbar(corr_ks,np.nanmean(autocorr_mc['AE'][0,:,:],axis=1),yerr=np.nanstd(autocorr_mc['AE'][0,:,:],axis=1))


# In[199]:


np.nanmean(autocorr_mc['AE'][0,:,:],axis=1)


# In[202]:


autocorr_len_mc['AE'][0,6,:]


# ### Functionalize for using multiprocessing

# In[400]:


from multiprocessing import Pool, cpu_count
from tqdm.notebook import tqdm 
import signal


# In[401]:


class KeyboardInterruptError(Exception): pass


# In[402]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[403]:


autocorr = {}
autocorr_len = {}


# In[404]:


corr_vals = ['aod1040','aod0500','AE','aod_fine','aod_coarse','fmf',
             'GOCI_AOD','GOCI_AE','GOCI_fmf',
             'MERRA_AOD','MERRA_AE',
             'situ_ext','situ_ae','situ_ssa']


# In[44]:


def calc_autocorr(val,types=types,corr_ks=corr_ks,dvals=dvals,iis=iis,itypes=itypes,range_vals=range_vals):
    dat = {'c':{},'l':{}}
    dat['c'][val] = np.zeros((len(types),len(corr_ks)))+np.nan
    dat['l'][val] = np.zeros((len(types),len(corr_ks)))
    print val
    for ik,k in enumerate(corr_ks):
        sys.stdout.write("{} ".format(ik))
        for j,jt in enumerate(types):
            corrp = np.hstack([dvals[val][i][iis[i][ik][1]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
            corrm = np.hstack([dvals[val][i][iis[i][ik][0]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
            badc = (corrp>range_vals[val][1]) | (corrp<range_vals[val][0])
            if any(badc): 
                corrp[badc] = np.nan
                corrm[badc] = np.nan
            
            badc = (corrm>range_vals[val][1]) | (corrm<range_vals[val][0])
            if any(badc): 
                corrp[badc] = np.nan
                corrm[badc] = np.nan
            
            mmk = np.nanmean(corrm)
            mpk = np.nanmean(corrp)
            smk = np.nanstd(corrm)
            spk = np.nanstd(corrp)
            top = (corrm-mpk)*(corrp-mmk) #[(v-mpk)*(corrp[iv]-mmk) for iv,v in enumerate(corrm)]
            #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
            dat['c'][val][j,ik] = np.nansum(top)/((len(corrm)-1)*spk*smk)
            dat['l'][val][j,ik] = len(corrm)
            if dat['l'][val][j,ik]<min_num:
                dat['c'][val][j,ik] = np.nan
                
    return dat


# In[405]:


subsamp_ratio = 0.30
N_mc = 50


# In[406]:


if (vv == 'v2') | (vv == 'v3'):
    corr_vals = ['aod_fine','aod_coarse','GOCI_AOD_c','GOCI_AOD_f',
                 'MERRA_AOD_dust','MERRA_AOD_sea','MERRA_AOD_bc','MERRA_AOD_oc','MERRA_AOD_sulf',
                 'aod1040','aod0500','AE','fmf',
                 'GOCI_AOD','GOCI_AE','GOCI_fmf',
                 'MERRA_AOD','MERRA_AE',
                 'situ_ext','situ_ae','situ_ssa']
    range_vals = {'aod_fine':[0.0,5.0],'aod_coarse':[0.0,5.0],'GOCI_AOD_c':[0.0,5.0],'GOCI_AOD_f':[0.0,5.0],
                  'MERRA_AOD_dust':[0.0,5.0],'MERRA_AOD_sea':[0.0,5.0],'MERRA_AOD_bc':[0.0,5.0],
                  'MERRA_AOD_oc':[0.0,5.0],'MERRA_AOD_sulf':[0.0,5.0],
                  'aod1040':[0.0,3.0],'aod0500':[0.0,5.0],'AE':[-1.0,5.0],
                  'fmf':[0.0,1.0],'GOCI_AOD':[0.0,5.0],'GOCI_AE':[-1.0,5.0],'GOCI_fmf':[0.0,1.0],
                  'MERRA_AOD':[0.0,5.0],'MERRA_AE':[-1.0,5.0],'situ_ext':[0.0,1000.0],
                  'situ_ae':[-1.0,5.0],'situ_ssa':[0.2,1.0]}
    


# In[407]:


def calc_autocorr_mc(val,types=types,corr_ks=corr_ks,dvals=dvals,iis=iis,itypes=itypes,range_vals=range_vals,
                     subsamp_ratio=subsamp_ratio,N_mc=N_mc):
    dat = {'c':{},'l':{}}
    dat['c'][val] = np.zeros((len(types),len(corr_ks),N_mc))+np.nan
    dat['l'][val] = np.zeros((len(types),len(corr_ks),N_mc))
    for ix in range(N_mc):
        print ix,val
        irand = np.random.randint(len(dvals['cumdist']),size=int(subsamp_ratio*len(dvals['cumdist'])))
        for ik,k in enumerate(corr_ks):
            sys.stdout.write("{} ".format(ik))
            for j,jt in enumerate(types):
                corrp = np.hstack([dvals[val][i][iis[i][ik][1]] for i in irand if i in itypes[j]])
                corrm = np.hstack([dvals[val][i][iis[i][ik][0]] for i in irand if i in itypes[j]])
                badc = (corrp>range_vals[val][1]) | (corrp<range_vals[val][0])
                if any(badc): 
                    corrp[badc] = np.nan
                    corrm[badc] = np.nan

                badc = (corrm>range_vals[val][1]) | (corrm<range_vals[val][0])
                if any(badc): 
                    corrp[badc] = np.nan
                    corrm[badc] = np.nan

                mmk = np.nanmean(corrm)
                mpk = np.nanmean(corrp)
                smk = np.nanstd(corrm)
                spk = np.nanstd(corrp)
                top = (corrm-mpk)*(corrp-mmk) #[(v-mpk)*(corrp[iv]-mmk) for iv,v in enumerate(corrm)]
                #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
                dat['c'][val][j,ik,ix] = np.nansum(top)/((len(corrm)-1)*spk*smk)
                dat['l'][val][j,ik,ix] = len(corrm)
                if dat['l'][val][j,ik,ix]<min_num:
                    dat['c'][val][j,ik,ix] = np.nan
                
    return dat


# In[408]:


p = Pool(7,worker_init)


# In[409]:


autocorr = {}
autocorr_len = {}
with tqdm(total=len(corr_vals)) as pbar:
    for i, outs in tqdm(enumerate(p.imap_unordered(calc_autocorr, corr_vals))):
        pbar.update()
        k = outs['c'].keys()[0]
        autocorr[k] = outs['c'][k]
        autocorr_len[k] = outs['l'][k]


# In[410]:


autocorr_mc = {}
autocorr_len_mc = {}
with tqdm(total=len(corr_vals)) as pbar:
    for i, outs in tqdm(enumerate(p.imap_unordered(calc_autocorr_mc, corr_vals))):
        pbar.update()
        k = outs['c'].keys()[0]
        autocorr_mc[k] = outs['c'][k]
        autocorr_len_mc[k] = outs['l'][k]


# In[411]:


autocorr_mc['GOCI_fmf'].shape


# ### Save

# In[412]:


vv = 'v3'


# In[413]:


dat_c = {u'autocorr_mc':autocorr_mc,'autocorr_len_mc':autocorr_len_mc}


# In[414]:


import write_utils as wu
dat_u = wu.iterate_dict_unicode(dat_c)


# In[415]:


hs.savemat(fp+'KORUS_fine_coarse_autocorr_mc_{}.mat'.format(vv),dat_u)


# ### Load the models

# ## Old calculation methods (no longer used)

# ### Autocorraltion calculation with all vars

# In[67]:


corr_ks =[0.08,0.1,0.25,0.5,0.75,1.0,1.5,2.0,3.0,5.0,7.5,10.0,12.5,15.0,20.0,
          25.0,30.0,35.0,40.0,50.0,60.0,75.0,100.0,150.0,200.0] 
corr_ks =[0.08,0.1,0.25,0.5,1.0,1.5,3.0,5.0,10.0,15.0,
          35.0,65.0,100.0,150.0,200.0] 
corr_all = [[[{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[],'fmf':[],'GOCI_AOD':[],'GOCI_AE':[],'GOCI_fmf':[],'MERRA_AOD':[],'MERRA_AE':[],'situ_ext':[],'situ_ae':[],'situ_ssa':[]} for i,k in enumerate(corr_ks)],
             [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[],'fmf':[],'GOCI_AOD':[],'GOCI_AE':[],'GOCI_fmf':[],'MERRA_AOD':[],'MERRA_AE':[],'situ_ext':[],'situ_ae':[],'situ_ssa':[]} for i,k in enumerate(corr_ks)]] \
            for j in types] 
corr_alln = [[[{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[],'fmf':[],'GOCI_AOD':[],'GOCI_AE':[],'GOCI_fmf':[],'MERRA_AOD':[],'MERRA_AE':[],'situ_ext':[],'situ_ae':[],'situ_ssa':[]} for i,k in enumerate(corr_ks)],
             [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[],'fmf':[],'GOCI_AOD':[],'GOCI_AE':[],'GOCI_fmf':[],'MERRA_AOD':[],'MERRA_AE':[],'situ_ext':[],'situ_ae':[],'situ_ssa':[]} for i,k in enumerate(corr_ks)]] \
            for j in types] 
#corr_all = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_t1 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_t2 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_t3 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_t4 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_z1 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_z2 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_z3 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
corr_vals = corr_all[0][0][0].keys()
corr_vals.remove('k')
corr_valsn = corr_alln[0][0][0].keys()
corr_valsn.remove('k')


# In[68]:


corr_all[0][0][0],corr_all[0][0][0]


# In[69]:


np.array(corr_all).shape #type, [minusk,plusk], distance


# In[70]:


len(dvals['cumdist'])


# In[71]:


dv = 0.05
for i,cd in enumerate(dvals['cumdist']):
    mat_dist = np.array([abs(cd-d) for d in cd])
    sys.stdout.write( '*{}*'.format(i))
    for ik,k in enumerate(corr_ks):
        #sys.stdout.write("{} ".format(ik))
        iimg,iipg = np.where((mat_dist>(k*(1.0-dv)))&(mat_dist<(k*(1.0+dv))))
        if len(iimg)==0:
            continue
        N = len(cd)
        for val in corr_vals:
            #all
            corr_alln[0][0][ik][val].append(dvals[val][i][iimg])
            corr_alln[0][1][ik][val].append(dvals[val][i][iipg])
            #type 1
            if (dvals['doys'][i][0]> t1[0]) & (dvals['doys'][i][0]< t1[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[1][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[1][1][ik][val].append(dvals[val][i][iipg])
            #type 2
            if (dvals['doys'][i][0]> t2[0]) & (dvals['doys'][i][0]< t2[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[2][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[2][1][ik][val].append(dvals[val][i][iipg])
            #type 3
            if (dvals['doys'][i][0]> t3[0]) & (dvals['doys'][i][0]< t3[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[3][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[3][1][ik][val].append(dvals[val][i][iipg])
            #type 4
            if (dvals['doys'][i][0]> t4[0]) & (dvals['doys'][i][0]< t4[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[4][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[4][1][ik][val].append(dvals[val][i][iipg])
            #type 5
            if (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z1[1]):
                corr_alln[5][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[5][1][ik][val].append(dvals[val][i][iipg])
            #type 6
            if (dvals['alt'][i][0]> z2[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[6][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[6][1][ik][val].append(dvals[val][i][iipg])
            #type 7
            if (dvals['alt'][i][0]> z3[0]) & (dvals['alt'][i][0]< z3[1]):
                corr_alln[7][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[7][1][ik][val].append(dvals[val][i][iipg])


# In[316]:


if False:
    dv = 0.15 #30% leeway in values
    for ik, k in enumerate(corr_ks):
        for i,cd in enumerate(dvals['cumdist']):
            if k>np.nanmax(cd):
                continue
            if ik==0:
                ipk = 1
                imk = len(cd)-2
                iip = range(1,len(cd)-1)
                iipg = iip
                iimg = range(0,imk)
            else:
                ipk = np.argmin(abs(cd-k)) #ipk:
                imk = np.argmin(abs(cd-(cd[-1]-k))) #0:imk
                if imk<2:
                    continue
                iip = Sp.find_closest(cd,cd[0:imk]+k)
                iip_good = (abs(cd[iip]-cd[0:imk])<k*(1.0+dv)) & (abs(cd[iip]-cd[0:imk])>k*(1.0-dv))
                if not any(iip_good):
                    continue
                iipg = iip[iip_good]
                iimg = np.arange(0,imk)[iip_good]
            N = len(cd)
            for val in corr_vals:
                #all
                corr_all[0][0][ik][val] = np.append(corr_all[0][0][ik][val],dvals[val][i][iimg])
                corr_all[0][1][ik][val] = np.append(corr_all[0][1][ik][val],dvals[val][i][iipg])
                #type 1
                if (dvals['doys'][i][0]> t1[0]) & (dvals['doys'][i][0]< t1[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[1][0][ik][val] = np.append(corr_all[1][0][ik][val],dvals[val][i][iimg])
                    corr_all[1][1][ik][val] = np.append(corr_all[1][1][ik][val],dvals[val][i][iipg])
                #type 2
                if (dvals['doys'][i][0]> t2[0]) & (dvals['doys'][i][0]< t2[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[2][0][ik][val] = np.append(corr_all[2][0][ik][val],dvals[val][i][iimg])
                    corr_all[2][1][ik][val] = np.append(corr_all[2][1][ik][val],dvals[val][i][iipg])
                #type 3
                if (dvals['doys'][i][0]> t3[0]) & (dvals['doys'][i][0]< t3[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[3][0][ik][val] = np.append(corr_all[3][0][ik][val],dvals[val][i][iimg])
                    corr_all[3][1][ik][val] = np.append(corr_all[3][1][ik][val],dvals[val][i][iipg])
                #type 4
                if (dvals['doys'][i][0]> t4[0]) & (dvals['doys'][i][0]< t4[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[4][0][ik][val] = np.append(corr_all[4][0][ik][val],dvals[val][i][iimg])
                    corr_all[4][1][ik][val] = np.append(corr_all[4][1][ik][val],dvals[val][i][iipg])
                #type 5
                if (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z1[1]):
                    corr_all[5][0][ik][val] = np.append(corr_all[5][0][ik][val],dvals[val][i][iimg])
                    corr_all[5][1][ik][val] = np.append(corr_all[5][1][ik][val],dvals[val][i][iipg])
                #type 6
                if (dvals['alt'][i][0]> z2[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[6][0][ik][val] = np.append(corr_all[6][0][ik][val],dvals[val][i][iimg])
                    corr_all[6][1][ik][val] = np.append(corr_all[6][1][ik][val],dvals[val][i][iipg])
                #type 7
                if (dvals['alt'][i][0]> z3[0]) & (dvals['alt'][i][0]< z3[1]):
                    corr_all[7][0][ik][val] = np.append(corr_all[7][0][ik][val],dvals[val][i][iimg])
                    corr_all[7][1][ik][val] = np.append(corr_all[7][1][ik][val],dvals[val][i][iipg])


# In[72]:


len(corr_alln[0][0][0]['AE'])


# In[73]:


corr_all[0][0][0]['AE'] = np.hstack(corr_alln[0][0][0]['AE'])


# In[74]:


for j,jt in enumerate(types):
    corr_all[j][1][0]['AE'] = np.hstack(corr_alln[j][1][0]['AE'])
    print corr_all[j][1][0]['AE'].shape


# In[75]:


len(corr_all)


# In[76]:


corr_vals


# In[77]:


types


# In[ ]:


if False:
    autocorr = {}
    autocorr_len = {}

    for val in corr_vals:
        autocorr[val] = np.zeros((len(types),len(corr_ks)))+np.nan
        autocorr_len[val] = np.zeros((len(types),len(corr_ks)))
        print val
        for ik,k in enumerate(corr_ks):
            sys.stdout.write("{} ".format(ik))
            for j,jt in enumerate(types):
                corr_all[j][0][ik][val] = np.hstack(corr_alln[j][0][ik][val])
                corr_all[j][1][ik][val] = np.hstack(corr_alln[j][1][ik][val])
                if False: #val is 'AE':
                    igd = np.where(corr_all[2][1][3]['AE']<5.0)[0] #np.isfinite(corr_all[j][0][ik][val])
                    mmk = np.nanmean(corr_all[j][0][ik][val][igd])
                    mpk = np.nanmean(corr_all[j][1][ik][val][igd])
                    smk = np.nanstd(corr_all[j][0][ik][val][igd])
                    spk = np.nanstd(corr_all[j][1][ik][val][igd])
                    top = [(v-mpk)*(corr_all[j][1][ik][val][igd][iv]-mmk) for iv,v in enumerate(corr_all[j][0][ik][val][igd])]
                    #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
                    autocorr[val][j,ik] = np.nansum(top)/((len(corr_all[j][0][ik][val][igd])-1)*spk*smk)
                else:
                    mmk = np.nanmean(corr_all[j][0][ik][val])
                    mpk = np.nanmean(corr_all[j][1][ik][val])
                    smk = np.nanstd(corr_all[j][0][ik][val])
                    spk = np.nanstd(corr_all[j][1][ik][val])
                    top = [(v-mpk)*(corr_all[j][1][ik][val][iv]-mmk) for iv,v in enumerate(corr_all[j][0][ik][val])]
                    #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
                    autocorr[val][j,ik] = np.nansum(top)/((len(corr_all[j][0][ik][val])-1)*spk*smk)
                    autocorr_len[val][j,ik] = len(corr_all[j][0][ik][val])
                    if autocorr_len[val][j,ik]<min_num:
                        autocorr[val][j,ik] = np.nan


# In[ ]:





# ### Save for offline calculations

# In[100]:


autocorr_dat = {u'corr_vals':corr_vals,u'types':types,u'corr_ks':corr_ks,
                u'corr_all':corr_all,u'corr_alln':corr_alln,u'min_num':min_num}


# In[96]:


import write_utils as wu
wu.dict_keys_to_unicode(autocorr_dat)


# In[ ]:


np.save(fp+'autocorr_dat_{}.npy'.foramt(vv),autocorr_dat,allow_pickle=True)


# ### Run offline

# In[106]:


get_ipython().system(u'python run_autocorr.py')


# In[ ]:


np.load


# In[ ]:


autocorr = {}
autocorr_len = {}


# In[79]:


val = 'aod_coarse'


# In[80]:


def run_autocorr(val,types=types,corr_ks=corr_ks,corr_all=corr_all,corr_alln=corr_alln,min_num=min_num):
    'to run the autocorrelate, and to help with garbage collection'
    print val
    autocorr = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_len = np.zeros((len(types),len(corr_ks)))
    for ik,k in enumerate(corr_ks):
        sys.stdout.write("{} ".format(ik))
        for j,jt in enumerate(types):
            corr_all[j][0][ik][val] = np.hstack(corr_alln[j][0][ik][val])
            corr_all[j][1][ik][val] = np.hstack(corr_alln[j][1][ik][val])

            mmk = np.nanmean(corr_all[j][0][ik][val])
            mpk = np.nanmean(corr_all[j][1][ik][val])
            smk = np.nanstd(corr_all[j][0][ik][val])
            spk = np.nanstd(corr_all[j][1][ik][val])
            top = [(v-mpk)*(corr_all[j][1][ik][val][iv]-mmk) for iv,v in enumerate(corr_all[j][0][ik][val])]
            #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
            autocorr[j,ik] = np.nansum(top)/((len(corr_all[j][0][ik][val])-1)*spk*smk)
            autocorr_len[j,ik] = len(corr_all[j][0][ik][val])
            if autocorr_len[j,ik]<min_num:
                autocorr[j,ik] = np.nan
    return autocorr, autocorr_len


# In[82]:


autocorr = {}
autocorr_len = {}
for val in corr_vals:
    autocorr[val],autocorr_len[val] = run_autocorr(val)


# In[2]:


autocorr_len[val].shape


# In[170]:


np.hstack(corr_all[j+1][0][ik][val])


# In[168]:


np.concatenate(corr_all[j+1][0][ik][val])


# In[ ]:





# In[318]:


autocorr['aod0500'].shape


# In[319]:


autocorr_len['aod0500'].shape


# ## Run monte Carlo sub-sampling to calculate variation in autocorrelation

# ### functions for monte carlo

# In[167]:


def fx_autocor(plus_dist,minus_dist,min_num=100):
    'Function to calculate the autocorrelation from 2 populations of the same distribution - plus and minus'
    # minus_dist = corr_all[j][0][ik][val]
    # plus_dist = corr_all[j][1][ik][val]
    mmk = np.nanmean(minus_dist)
    mpk = np.nanmean(plus_dist)
    smk = np.nanstd(minus_dist)
    spk = np.nanstd(plus_dist)
    top = (minus_dist-mpk)*(plus_dist-mmk) #[(v-mpk)*(plus_dist[iv]-mmk) for iv,v in enumerate(minus_dist)]
    #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
    autocorr = np.nansum(top)/((len(minus_dist)-1)*spk*smk)
    autocorr_len = len(minus_dist)
    if autocorr_len<min_num:
        autocorr = np.nan
    return autocorr, autocorr_len


# In[168]:


def rand_autocor(plus_dist,minus_dist,min_num=100,subsamp_ratio=0.5):
    'fx_autocor, but called with random subsampling at a ratio defined by subsamp_ratio [0-1]'
    irand = np.random.randint(len(plus_dist),size=int(subsamp_ratio*len(plus_dist)))
    return fx_autocor(plus_dist[irand],minus_dist[irand],min_num=min_num)


# In[169]:


j,ik,val


# In[55]:


len(corr_all[j][0][ik][val])


# In[56]:


autocorr[val][j,ik],autocorr_len[val][j,ik]


# In[325]:


fx_autocor(corr_all[j][0][ik][val],corr_all[j][1][ik][val])


# In[326]:


rand_autocor(corr_all[j][0][ik][val],corr_all[j][1][ik][val],subsamp_ratio=0.5)


# In[ ]:





# In[327]:


autocorr_std = {}
autocorrs = {}
num_for_std = 50

for val in corr_vals:
    autocorrs[val] = np.zeros((len(types),len(corr_ks),num_for_std))+np.nan
    autocorr_std[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    print val
    
    for ik,k in enumerate(corr_ks):
        for j,jt in enumerate(types):
            if not len(corr_all[j][0][ik][val]):
                continue
            for u in xrange(num_for_std):
                autocorrs[val][j,ik,u],l = rand_autocor(corr_all[j][0][ik][val],corr_all[j][1][ik][val],subsamp_ratio=0.5)
            autocorr_std[val][j,ik] = np.nanstd(autocorrs[val][j,ik,:])


# ### Reality check of matching monte carlo calcs to standards

# In[57]:


corr_ks


# In[58]:


val = 'AE'
ik = 12
k = corr_ks[ik]
j = 0
jt = types[j]
corrp = np.hstack([dvals[val][i][iis[i][ik][1]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
corrm = np.hstack([dvals[val][i][iis[i][ik][0]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
badc = (corrp>range_vals[val][1]) | (corrp<range_vals[val][0])
if any(badc): 
    corrp[badc] = np.nan
    corrm[badc] = np.nan

badc = (corrm>range_vals[val][1]) | (corrm<range_vals[val][0])
if any(badc): 
    corrp[badc] = np.nan
    corrm[badc] = np.nan


# In[59]:


autocorr[val][j,ik]


# In[61]:


rand_autocor(corrm,corrp,subsamp_ratio=0.4)


# ### Apply serial processing

# In[316]:


autocorr_std = {}
autocorrs = {}
num_for_std = 50

for val in corr_vals:
    autocorrs[val] = np.zeros((len(types),len(corr_ks),num_for_std))+np.nan
    autocorr_std[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    print val
    
    for ik,k in enumerate(corr_ks):
        for j,jt in enumerate(types):
            
            corrp = np.hstack([dvals[val][i][iis[i][ik][1]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
            corrm = np.hstack([dvals[val][i][iis[i][ik][0]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
            badc = (corrp>range_vals[val][1]) | (corrp<range_vals[val][0])
            if any(badc): 
                corrp[badc] = np.nan
                corrm[badc] = np.nan
            
            badc = (corrm>range_vals[val][1]) | (corrm<range_vals[val][0])
            if any(badc): 
                corrp[badc] = np.nan
                corrm[badc] = np.nan
            
            if not len(corrp):
                continue
            for u in xrange(num_for_std):
                autocorrs[val][j,ik,u],l = rand_autocor(corrm,corrp,subsamp_ratio=0.3)
            autocorr_std[val][j,ik] = np.nanstd(autocorrs[val][j,ik,:])


# ### Apply multi-processing

# In[269]:


if (vv == 'v2') | (vv == 'v3'):
    corr_vals = ['aod_fine','aod_coarse','GOCI_AOD_c','GOCI_AOD_f',
                 'MERRA_AOD_dust','MERRA_AOD_sea','MERRA_AOD_bc','MERRA_AOD_oc','MERRA_AOD_sulf']
    range_vals = {'aod_fine':[0.0,5.0],'aod_coarse':[0.0,5.0],'GOCI_AOD_c':[0.0,5.0],'GOCI_AOD_f':[0.0,5.0],
                  'MERRA_AOD_dust':[0.0,5.0],'MERRA_AOD_sea':[0.0,5.0],'MERRA_AOD_bc':[0.0,5.0],
                  'MERRA_AOD_oc':[0.0,5.0],'MERRA_AOD_sulf':[0.0,5.0]}


# In[270]:


autocorr_std = {}
autocorrs = {}
num_for_std = 50

#for val in corr_vals:


# In[184]:


def montecarlo_autocorr(val,dvals=dvals,types=types,num_for_std=num_for_std,corr_ks=corr_ks,
                        range_vals=range_vals,iis=iis,itypes=itypes):
    'Function to enable multiprocessing of montecarlo'
    autocorr_std = {}
    autocorrs = {}
    
    autocorrs[val] = np.zeros((len(types),len(corr_ks),num_for_std))+np.nan
    autocorr_std[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    print val
    
    for ik,k in enumerate(corr_ks):
        for j,jt in enumerate(types):
            
            corrp = np.hstack([dvals[val][i][iis[i][ik][1]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
            corrm = np.hstack([dvals[val][i][iis[i][ik][0]] for i,cd in enumerate(dvals['cumdist']) if i in itypes[j]])
            badc = (corrp>range_vals[val][1]) | (corrp<range_vals[val][0])
            if any(badc): 
                corrp[badc] = np.nan
                corrm[badc] = np.nan
            
            badc = (corrm>range_vals[val][1]) | (corrm<range_vals[val][0])
            if any(badc): 
                corrp[badc] = np.nan
                corrm[badc] = np.nan
            
            if not len(corrp):
                continue
            for u in xrange(num_for_std):
                autocorrs[val][j,ik,u],l = rand_autocor(corrm,corrp,subsamp_ratio=0.3)
            autocorr_std[val][j,ik] = np.nanstd(autocorrs[val][j,ik,:])
    
    return val,autocorrs[val],autocorr_std[val]


# In[185]:


from multiprocessing import Pool, cpu_count
from tqdm.notebook import tqdm 
import signal


# In[186]:


class KeyboardInterruptError(Exception): pass


# In[187]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[188]:


p = Pool(5,worker_init)


# In[189]:


with tqdm(total=len(corr_vals)) as pbar:
    for i, outs in tqdm(enumerate(p.imap_unordered(montecarlo_autocorr, corr_vals))):
        pbar.update()
        k = outs[0]
        autocorrs[k] = outs[1]
        autocorr_std[k] = outs[2]


# ### combine results and calc stats

# In[190]:


autocorr_min = {}
autocorr_max = {}
autocorr_d = {}
autocorr_mean = {}

for val in corr_vals:
    autocorr_min[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_max[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_d[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_mean[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    
    for ik,k in enumerate(corr_ks):
        for j,jt in enumerate(types):
            #if not len(corr_all[j][0][ik][val]):
            #    continue
            autocorr_min[val][j,ik] = np.nanmin(autocorrs[val][j,ik,:])
            autocorr_max[val][j,ik] = np.nanmax(autocorrs[val][j,ik,:])
            autocorr_d[val][j,ik] = autocorr_max[val][j,ik] - autocorr_min[val][j,ik]
            autocorr_mean[val][j,ik] = np.nanmean(autocorrs[val][j,ik,:])


# ### Save autocorr and monte carlo calculations to file

# In[191]:


vv = 'v2'


# In[193]:


dat_c = {u'autocorrs':autocorrs, u'autocorr_std':autocorr_std,u'autocorr_min':autocorr_min,u'autocorr_max':autocorr_max,
        u'autocorr_d':autocorr_d,u'autocorr_mean':autocorr_mean,
        u'corr_vals':corr_vals,u'corr_ks':corr_ks,u'types':types,u'num_for_std':num_for_std}


# In[192]:


dat_c = {u'autocorrs':autocorrs, u'autocorr_std':autocorr_std,u'autocorr_min':autocorr_min,u'autocorr_max':autocorr_max,
        u'autocorr_d':autocorr_d,u'autocorr_mean':autocorr_mean,
        u'corr_vals':corr_vals,u'corr_ks':corr_ks,u'types':types,u'num_for_std':num_for_std,
        u'autocorr':autocorr,'autocorr_len':autocorr_len}


# In[194]:


import write_utils as wu
dat_u = wu.iterate_dict_unicode(dat_c)


# In[195]:


hs.savemat(fp+'KORUS_fine_coarse_autocorr_montecarlo_{}.mat'.format(vv),dat_u)


# ### Load autocorr and monte carlo

# In[196]:


import write_utils as wu


# In[204]:


vv = 'v1'
dat_u = hs.loadmat(fp+'KORUS_fine_coarse_autocorr_montecarlo_{}.mat'.format(vv))


# In[207]:


for k in dat_u.keys():
    if (k in locals()) & (isinstance(dat_u[k],dict)):
        locals()[k] = wu.merge_dicts(dat_u[k],locals()[k])
        print 'merging dicts: '+k
    else:
        locals()[k] = dat_u[k]
        print 'saving variable: '+k


# In[209]:


autocorrs.keys()


# ## Internal old testing

# ### Integrated autocorrelation

# In[193]:


def autocorr(x, t=1):
    return np.corrcoef(np.array([x[:-t], x[t:]]))


# In[194]:


def autocorr2(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]/result.max()


# In[195]:


def autocorr5(x):
    '''numpy.correlate, non partial'''
    n = len(x)
    lags = range(n)
    #mean=x.mean()
    var=np.nanvar(x)
    xp=x-np.nanmean(x)
    corr=np.correlate(xp,xp,'full')[n-1:]/var/n

    return corr[:n]


# In[196]:


authcor = autocorr(dvals['aod0500'][1])
authcor2 = autocorr2(dvals['aod0500'][1])
authcor3 = autocorr5(dvals['aod0500'][1])


# In[197]:


len(authcor2)


# In[237]:


len(dvals['aod0500'][1])


# In[238]:


dvals['dist'][1]


# In[ ]:


[(dvals['dist'][i].mean(),np.size(dvals['dist'][i])) for i in xrange(len(dvals['dist']))]


# ### interpolate AODs to a constant distance grid

# In[204]:


def interp_dist(d,dist=0.12,verbose=False):
    'function to insterpolate the AOD from the dict to an even grid spacing accroding to distance (default 0.12 km)'
    d['cdist_n'],d['aod_n'] = [],[]
    for i,cd in enumerate(d['cumdist']):
        if verbose:
            print i, cd.min(),cd.max(), np.nanmin(cd),np.nanmax(cd)
            if not np.isfinite(cd.min()): print cd
        d['cdist_n'].append(np.arange(cd.min(),cd.max(),dist))
        try:
            fcd = interpolate.interp1d(cd,d['aod0500'][i])
            d['aod_n'].append(fcd(d['cdist_n'][i]))
        except TypeError:
            d['aod_n'].append(np.array(np.nan))


# In[ ]:


dvals['aod0500']


# In[205]:


interp_dist(dvals)


# In[206]:


dvals['autocor'] = [] 
for i,a in enumerate(dvals['aod_n']):
    try:
        dvals['autocor'].append(autocorr5(a))
    except:
        dvals['autocor'].append(np.array(np.nan))


# ### Autocorrelation plots (old)

# In[207]:


plt.figure()
plt.plot(dvals['cdist_n'][1],dvals['autocor'][1],'.')
plt.xlabel('Lag Distance [km]')
plt.ylabel('Correlation')


# In[208]:


plt.figure()
for i,j in enumerate(dvals['cdist_n']):
    try:
        plt.plot(j,dvals['aod_n'][i],'.')
    except:
        pass


# In[209]:


plt.figure()
for i,j in enumerate(dvals['cdist_n']):
    try:
        plt.plot(j,dvals['autocor'][i],'.')
    except:
        pass
plt.ylabel('Correlation Coefficient')
plt.xlabel('Distance [km]')
plt.xscale('log')


# ### Autocorrelation plot with Anderson method

# In[373]:


key_list = ['aod0500','aod1040','AE','aod_fine','aod_coarse']
key_list2 = ['aod0500','aod1040','AE','fmf']
legend_list = ['All','Dynamic','Stagnation','Extreme pollution','Blocking','0-1 km','1-3 km','3+ km']
cl_list = ['k','tab:red','tab:blue','tab:orange','tab:green','tab:olive','tab:cyan','tab:purple']
m_list = ['.','o','s','v','^','*','+','x']
tit = ['AOD$_{{500}}$','AOD$_{{1040}}$','AE','AOD$_{{fine}}$','AOD$_{{coarse}}$']
tit2 = ['AOD$_{{500}}$','AOD$_{{1040}}$','AE','fine-mode fraction']


# In[705]:


fig, ax = plt.subplots(5,3,figsize=(12,9))
for i,k in enumerate(key_list):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks,autocorr[k][j,:],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks,autocorr[k][j,:],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
    if i==1:
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.4])
    if i==4:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_all.png',dpi=600,transparent=True)


# In[240]:


key_list


# In[329]:


fig, ax = plt.subplots(5,3,figsize=(12,9))
for i,k in enumerate(key_list):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
    if i==1:
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.4])
    if i==4:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)


# ## Load the Autocorrelations from Shinozuka & Redemann

# In[246]:


SR_corr_ks = [0.45,1.0,3.0,6.0,10.0,20.0,34.2]
SR_aod_corr_long = [0.998,0.997,0.995,0.985,0.981,0.946,0.917]
SR_aod_corr_loc = [0.975,0.941,0.830,0.712,0.584,0.365,0.328]
SR_AE_corr_long = [1.000,0.977,0.975,0.960,0.944,0.913,0.919]
SR_AE_corr_loc = [0.975,0.956,0.919,0.831,0.747,0.519,0.366]


# In[342]:


#AOD local

0,4376084276355722; 0,9757725145572111
0,9955912496957824; 0,9408298503197231
2,993269282069332; 0,8303668774336445
5,997701705897239; 0,711850105383489
9,981245014250291; 0,5844939806380168
20,005524337736603; 0,36462079805665715
35,166804844012965; 0,32825206301575405


# In[343]:


#AE local
0,43836153494904556; 0,9757721573250459
0,9972623352215253; 0,9556260493694856
2,9925078192756667; 0,9191462151252102
5,985354161428386; 0,8309627406851716
9,993780925398312; 0,7465159146929592
19,996703813167056; 0,5185049833887045
35,22351663341235; 0,36598292430250434


# In[344]:


#AOD long
0,43833272405336876; 0,9987068195620337
0,9971397077970326; 0,9985360625870756
2,9918543848699146; 0,9953484799771375
6,0033122296001995; 0,985586039366985
9,969908496905848; 0,9810416889936774
19,972213345966814; 0,9461254599364131
35,10743016861416; 0,9178948308505701


# In[154]:


#AE long
0,43757781072798063; 1,0001868324223917
0,9989150464154728; 0,9778205265602119
2,992025659700471; 0,9753731289965352
5,9831337971357215; 0,9604326081520385
9,970944322675367; 0,9447901261029547
19,97407661698948; 0,9135730361161722
35,10728129935923; 0,9193744864787629


# In[244]:


note = [['a)','b)'],['c)','d)'],['e)','f)'],['g)','h)'],['i)','j)']]


# In[718]:


fig, ax = plt.subplots(5,3,figsize=(12,9))
for i,k in enumerate(key_list):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks,autocorr[k][j,:],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks,autocorr[k][j,:],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.1])
    
    if i==3:
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==4:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_all_with_SR2011.png',dpi=600,transparent=True)


# In[244]:


key_list


# In[245]:


autocorr.keys()


# In[1003]:


fig, ax = plt.subplots(4,3,figsize=(12,9))
for i,k in enumerate(key_list2):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'AE':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'AE':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit2[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.1])
    
    if i==3:
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==3:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_rel_all_with_SR2011_with_GOCI.png',dpi=600,transparent=True)


# ### Add MERRA

# In[1237]:


fig, ax = plt.subplots(4,3,figsize=(12,9))
for i,k in enumerate(key_list2):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
        if k is 'AE':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
        if k is 'AE':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit2[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.1])
    
    if i==3:
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==3:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_rel_all_with_SR2011_GOCI_MERRA.png',dpi=600,transparent=True)


# ### Add In situ

# In[227]:


fig, ax = plt.subplots(4,3,figsize=(12,9))
for i,k in enumerate(key_list2):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
                        color=cl_list[j],ls='-.',lw=1)
        if k is 'AE':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'AE':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit2[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].plot([],[],'-.',c='grey',lw=1,label='In situ Ext.')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.1])
    
    if i==3:
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==3:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_rel_all_with_SR2011_GOCI_MERRA_insitu.png',dpi=600,transparent=True)


# In[228]:


fig, ax = plt.subplots(4,3,figsize=(12,9))
for i,k in enumerate(key_list2):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
                        color=cl_list[j],ls='-.',lw=1)
        if k is 'AE':
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'AE':
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(-1,1)
    ax[i,1].set_ylim(-1,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('diff r({})'.format(tit2[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        auk_norm = wu.nearest_neighbor(np.array(corr_ks[1:]),autocorr[k][0,1:]/autocorr[k][0,1],np.array(SR_corr_ks),dist=30.0)
        
        ax[i,0].plot(SR_corr_ks,auk_norm-SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,auk_norm-SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].plot([],[],'-.',c='grey',lw=1,label='In situ Ext.')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.6,0.1])
    
    if i==3:
        auk_norm = wu.nearest_neighbor(np.array(corr_ks[1:]),autocorr[k][0,1:]/autocorr[k][0,1],np.array(SR_corr_ks),dist=30.0)
        ax[i,0].plot(SR_corr_ks,auk_norm-SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,auk_norm-SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==3:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_diff_rel_all_with_SR2011_GOCI_MERRA_insitu.png',dpi=600,transparent=True)


# ## Only plot subset of autocorrelation extensive vs intensive

# In[15]:


autocorr_len.keys()


# In[418]:


key_list = ['aod0500','aod1040','AE','aod_fine','aod_coarse']
key_list2 = ['aod0500','aod1040','AE','fmf']
legend_list = ['All','Dynamic','Stagnation','Extreme pollution','Blocking','0-1 km','1-3 km','3+ km']
cl_list = ['k','tab:red','tab:blue','tab:orange','tab:green','tab:olive','tab:cyan','tab:purple']
m_list = ['.','o','s','v','^','*','+','x']
tit = ['AOD$_{{500}}$','AOD$_{{1040}}$','AE','AOD$_{{fine}}$','AOD$_{{coarse}}$']
tit2 = ['AOD$_{{500}}$','AOD$_{{1040}}$','AE','fine-mode fraction']


# In[17]:


fig, ax = plt.subplots(2,2,figsize=(6,4))
i = 0
for j in [0,1,2,3,4]:
    ax[i,0].plot(corr_ks[1:],autocorr_len['aod0500'][j,1:],label=legend_list[j],color=cl_list[j],marker=m_list[j])

    ax[i,0].plot(corr_ks[1:],autocorr_len['GOCI_AOD'][j,1:],
                 color=cl_list[j],ls=':',lw=1)
    ax[i,0].plot(corr_ks[1:],autocorr_len['MERRA_AOD'][j,1:],
                color=cl_list[j],ls='--',lw=1)
    ax[i,0].plot(corr_ks[1:],autocorr_len['situ_ext'][j,1:],
                 color=cl_list[j],ls='-.',lw=1)
for j in [0,5,6,7]:    
    ax[i,1].plot(corr_ks[1:],autocorr_len['aod0500'][j,1:],label=legend_list[j],color=cl_list[j],marker=m_list[j])

    ax[i,1].plot(corr_ks[1:],autocorr_len['GOCI_AOD'][j,1:],
                 color=cl_list[j],ls=':',lw=1)
    ax[i,1].plot(corr_ks[1:],autocorr_len['MERRA_AOD'][j,1:],
                 color=cl_list[j],ls='--',lw=1)
    ax[i,1].plot(corr_ks[1:],autocorr_len['situ_ext'][j,1:],
                 color=cl_list[j],ls='-.',lw=1)
#ax[i,0].set_ylim(0,1)
#ax[i,1].set_ylim(0,1)
ax[i,0].set_xscale('log')
ax[i,1].set_xscale('log')
ax[i,0].set_yscale('log')
ax[i,1].set_yscale('log')
ax[i,0].grid()
ax[i,1].grid()

ax[1,0].set_visible(False)
ax[1,1].set_visible(False)


#print 'r({})'.format(k)
ax[i,0].set_ylabel('Number of points'.format(tit2[i]))
pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)

ax[i,0].set_title('Meteorology')
ax[i,1].set_title('Altitude')
ax[i,0].set_xlabel('Distance [km]')
ax[i,1].set_xlabel('Distance [km]')

ax[i,0].legend(frameon=False,bbox_to_anchor=[0.85,-0.25])
ax[i,1].legend(frameon=False,bbox_to_anchor=[0.95,-0.25])
    

plt.savefig(fp+'plot/KORUS_Autocorr_number_{}.png'.format(vv),dpi=600,transparent=True)


# In[18]:


key_list2


# In[19]:


tit2


# In[235]:


key_list3 = ['aod0500','fmf']
tit3 = ['AOD$_{{500}}$','fine-mode fraction']


# In[302]:


fig, ax = plt.subplots(2,3,figsize=(12,3.5))
for i,k in enumerate(key_list3):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
       #     ax[i,0].plot(corr_ks[1:],autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
       #                 color=cl_list[j],ls='-.',lw=1)
       # if k is 'AE':
       #     ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
       #                  color=cl_list[j],ls=':',lw=1)
       #     ax[i,0].plot(corr_ks[1:],autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
       #                  color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
        #    ax[i,1].plot(corr_ks[1:],autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
        #                 color=cl_list[j],ls='-.',lw=1)
        #if k is 'AE':
        #    ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
        #                 color=cl_list[j],ls=':',lw=1)
        #    ax[i,1].plot(corr_ks[1:],autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
        #                 color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,1].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit3[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=False,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=False,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==0:
    #    ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
    #    ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
    #    ax[i,0].plot([],[],'-.',c='grey',lw=1,label='In situ Ext.')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.1])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[2.3,1.1])

    #if i==1:
    #    ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
    #    ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==1:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.subplots_adjust(left=None, bottom=0.15, right=None, top=None,
                wspace=None, hspace=None)
        
plt.savefig(fp+'plot/KORUS_Autocorr_rel_subset_with_SR2011_GOCI_MERRA_insitu.png',dpi=600,transparent=True)


# ### Add the standard deviation subset

# In[282]:


corr_ks = np.array(corr_ks)


# In[354]:


fig, ax = plt.subplots(2,3,figsize=(14,4.5))
for i,k in enumerate(key_list3):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],
                         yerr=autocorr_std[k][j,1:]/autocorr[k][j,1],color=cl_list[j],
                     marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
            ax[i,0].errorbar(corr_ks[1:]*1.02,autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1], 
                             yerr=autocorr_std['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],capsize=0,ls=':',lw=1)
            ax[i,0].errorbar(corr_ks[1:]*0.98,autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1], 
                             yerr=autocorr_std['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],capsize=0,ls='--',lw=1)

        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].errorbar(corr_ks[1:]*1.02,autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                             yerr=autocorr_std['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1,capsize=0)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,1].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],
                         yerr=autocorr_d[k][j,1:]/autocorr[k][j,1],color=cl_list[j],
                     marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
            ax[i,1].errorbar(corr_ks[1:]*1.02,autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1], 
                             yerr=autocorr_d['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],capsize=0,ls=':',lw=1)
            ax[i,1].errorbar(corr_ks[1:]*0.98,autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1], 
                             yerr=autocorr_d['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],capsize=0,ls='--',lw=1)

        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].errorbar(corr_ks[1:]*1.02,autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                             yerr=autocorr_d['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                            color=cl_list[j],ls=':',lw=1,capsize=0)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,1].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit3[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==0:
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.1])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[2.3,1.1])

    if i==1:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.subplots_adjust(left=None, bottom=0.15, right=None, top=None,
                wspace=None, hspace=None)
        
plt.savefig(fp+'plot/KORUS_Autocorr_std_rel_subset_with_SR2011_GOCI_MERRA_{}.png'.format(vv),dpi=600,transparent=True)


# In[350]:


fig, ax = plt.subplots(2,3,figsize=(12,3.5))
for i,k in enumerate(key_list3):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],yerr=autocorr_d[k][j,1:]/autocorr[k][j,1],
                         color=cl_list[j],marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)

        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],yerr=autocorr_d[k][j,1:]/autocorr[k][j,1],
                         color=cl_list[j],marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)

        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,1].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit3[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=False,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=False,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==0:
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.1])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[2.3,1.1])

    if i==1:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.subplots_adjust(left=None, bottom=0.15, right=None, top=None,
                wspace=None, hspace=None)
        
#plt.savefig(fp+'plot/KORUS_Autocorr_dif_rel_subset_with_SR2011_GOCI_MERRA.png',dpi=600,transparent=True)


# In[437]:


fig, ax = plt.subplots(2,3,figsize=(12,3.5))
for i,k in enumerate(key_list3):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],yerr=autocorr_d[k][j,1:]/autocorr[k][j,1],
                         color=cl_list[j],marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)

        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],yerr=autocorr_d[k][j,1:]/autocorr[k][j,1],
                         color=cl_list[j],marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)

        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1.1)
    ax[i,1].set_ylim(0,1.1)
    ax[i,0].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,1].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('Normalized\nr({})'.format(tit3[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==0:
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[2.7,1.25],title='Meteorology')
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.055,-0.55],title='Altitude')

    if i==1:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.subplots_adjust(left=0.07, bottom=0.15, right=1.25, top=0.94,
                wspace=None, hspace=None)
        
plt.savefig(fp+'plot/KORUS_Autocorr_dif_rel_subset_with_SR2011_GOCI_MERRA_{}.png'.format(vv).format(vv),dpi=600,transparent=True)


# ### With the v3 monte carlo based on segments instead of samples

# In[438]:


autocorr_mm = {}
autocorr_dm = {}
for k in autocorr_mc.keys():
    autocorr_mm[k] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_dm[k] = np.zeros((len(types),len(corr_ks)))+np.nan
    for j,jt in enumerate(types):
        autocorr_mm[k][j,:] = np.nanmean(autocorr_mc[k][j,:,:],axis=1)
        autocorr_dm[k][j,:] = np.nanstd(autocorr_mc[k][j,:,:],axis=1)


# In[439]:


cl_list = ['k','tab:red','tab:blue','tab:orange','tab:green','tab:olive','tab:cyan','tab:purple']
key_list3 = ['aod0500','fmf']
tit3 = ['AOD$_{{500}}$','fine-mode fraction']


# In[444]:


fig, ax = plt.subplots(2,3,figsize=(12,3.5))
for i,k in enumerate(key_list3):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr_mm[k][j,1:]/autocorr_mm[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:]*(0.94+j/50.0),autocorr_mm[k][j,1:]/autocorr_mm[k][j,1],yerr=autocorr_dm[k][j,1:]/autocorr_mm[k][j,1],
                         color=cl_list[j],marker=None,ls='none',alpha=0.6,capsize=1.0,elinewidth=0.6)
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr_mm['GOCI_AOD'][j,1:]/autocorr_mm['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr_mm['MERRA_AOD'][j,1:]/autocorr_mm['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)

        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr_mm['GOCI_fmf'][j,1:]/autocorr_mm['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr_mm[k][j,1:]/autocorr_mm[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,1].errorbar(corr_ks[1:]*(0.92+j/50.0),autocorr_mm[k][j,1:]/autocorr_mm[k][j,1],yerr=autocorr_dm[k][j,1:]/autocorr_mm[k][j,1],
                         color=cl_list[j],marker=None,ls='none',alpha=0.6,capsize=1.0,elinewidth=0.6)
        #if k is 'aod0500':
        #    ax[i,1].plot(corr_ks[1:],autocorr_mm['GOCI_AOD'][j,1:]/autocorr_mm['GOCI_AOD'][j,1],
        #                 color=cl_list[j],ls=':',lw=1)
        #    ax[i,1].plot(corr_ks[1:],autocorr_mm['MERRA_AOD'][j,1:]/autocorr_mm['MERRA_AOD'][j,1],
        #                 color=cl_list[j],ls='--',lw=1)

        #if k is 'fmf':
        #    ax[i,1].plot(corr_ks[1:],autocorr_mm['GOCI_fmf'][j,1:]/autocorr_mm['GOCI_fmf'][j,1],
        #                 color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1.1)
    ax[i,1].set_ylim(0,1.1)
    ax[i,0].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,1].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('Normalized\nr({})'.format(tit3[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==0:
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[2.7,1.25],title='Meteorology')
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.055,-0.55],title='Altitude')

    if i==1:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.subplots_adjust(left=0.07, bottom=0.15, right=1.25, top=0.94,
                wspace=None, hspace=None)
        
plt.savefig(fp+'plot/KORUS_Autocorr_dif_rel_subset_with_SR2011_GOCI_MERRA_{}.png'.format(vv).format(vv),dpi=600,transparent=True)


# ## Autocorrelation of extrinsic vs intrisic properties

# In[285]:


autocorr.keys()


# In[286]:


autocorrs['AE'][j,jj:]


# In[27]:


autocorr_std['AE'][j,jj:]


# In[214]:


autocorrs.keys()


# In[218]:


plt.figure()
jj = 0

plt.plot(corr_ks,autocorrs['aod_coarse'][j,jj:,:])
plt.xlabel('Distance [km]')


# ### Looking into speciation

# In[360]:


vv = 'v2'


# In[222]:


fig, ax = plt.subplots(1,1,sharey=True,figsize=(8,4))
j = 0
jj = 0
ax = [ax]
#ax[0].plot(corr_ks[jj:],autocorr_mean['AE'][j,jj:],ls='--',marker='.',label='AE',color='k')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['aod_fine'][j,jj:],yerr=autocorr_d['aod_fine'][j,jj:],
               ls='--',marker='.',label='4STAR fine')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['aod_coarse'][j,jj:],yerr=autocorr_d['aod_coarse'][j,jj:],
               ls='--',marker='.',label='4STAR coarse')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD_sulf'][j,jj:],yerr=autocorr_d['MERRA_AOD_sulf'][j,jj:],
               ls='--',marker='.',label='MERRA_AOD_sulf')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD_oc'][j,jj:],yerr=autocorr_d['MERRA_AOD_oc'][j,jj:],
               ls='--',marker='.',label='MERRA_AOD_oc')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD_bc'][j,jj:],yerr=autocorr_d['MERRA_AOD_bc'][j,jj:],
               ls='--',marker='.',label='MERRA_AOD_bc')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD_dust'][j,jj:],yerr=autocorr_d['MERRA_AOD_dust'][j,jj:],
               ls='--',marker='.',label='MERRA_AOD_dust')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD_sea'][j,jj:],yerr=autocorr_d['MERRA_AOD_sea'][j,jj:],
               ls='--',marker='.',label='MERRA_AOD_sea')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['GOCI_AOD_f'][j,jj:],yerr=autocorr_d['GOCI_AOD_f'][j,jj:],
               ls='--',marker='.',label='GOCI fine')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['GOCI_AOD_c'][j,jj:],yerr=autocorr_d['GOCI_AOD_c'][j,jj:],
               ls='--',marker='.',label='GOCI coarse')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['situ_ae'][j,jj:],yerr=autocorr_d['situ_ae'][j,jj:],
               ls='--',marker='.',label='Point / in situ')
ax[0].legend(frameon=False)
ax[0].set_ylim(0,1)
ax[0].set_yticks([0,0.25,0.5,0.75,1.0])
ax[0].set_xscale('log')


# In[260]:


e_folding_autocorr = lambda x: corr_ks[np.argmin(abs(x[0]*np.exp(-1.0)-x))]
percentile_autocorr = lambda x,p: corr_ks[np.argmin(abs(x[0]*p-x))]


# In[261]:


ka = autocorr_mean.keys()
ka.sort()


# In[262]:


kas = ['GOCI_AOD_c','GOCI_AOD_f','MERRA_AOD_bc','MERRA_AOD_dust','MERRA_AOD_oc','MERRA_AOD_sea','MERRA_AOD_sulf',
       'aod_coarse','aod_fine']


# In[264]:


legend_list = ['All','Dynamic','Stagnation','Extreme\npollution','Blocking','0-1 km','1-3 km','3+ km']


# In[290]:


e_fold = {}
perc = {}
for k in types:
    e_fold[k] = []
    perc[k] = []
jj = 0
for j,t in enumerate(types):
    print 'Type {}: label {}'.format(t,legend_list[j])
    for k in kas:
        print '..{}: efolding={}, 90%={}, 85%={}'.format(k,
        e_folding_autocorr(autocorr_mean[k][j,jj:]),
        percentile_autocorr(autocorr_mean[k][j,jj:],0.9),percentile_autocorr(autocorr_mean[k][j,jj:],0.85))
        e_fold[t].append(e_folding_autocorr(autocorr_mean[k][j,jj:]))
        perc[t].append(percentile_autocorr(autocorr_mean[k][j,jj:],0.85))


# In[291]:


e_foldr = [e_fold[k] for k in types]
percr = [perc[k] for k in types]


# In[292]:


for i,e in enumerate(e_foldr):
    e.insert(0,legend_list[i])
    percr[i].insert(0,legend_list[i])


# In[293]:


pd_fold = pd.DataFrame(percr,columns=['types'] + kas) 


# In[431]:


cl_list = plt.cm.tab20(range(12)*21)


# In[327]:


pd_fold


# In[357]:


fig, ax = plt.subplots(2,1,sharex=True,figsize=(10,3))

pd_fold.plot(x='types',
        kind='bar',
        stacked=False,
        title='Autocorrelation distance at 85th percentile for speciated AOD',
        ax=ax[0],grid=True,logy=True,rot=0,color=cl_list[1::2],
        y=['aod_coarse','GOCI_AOD_c','MERRA_AOD_dust','MERRA_AOD_sea'])
ax[0].legend(['4STAR Coarse','GOCI Coarse','MERRA Dust','MERRA Sea salt'],
             frameon=False,loc='center left', bbox_to_anchor=(1.0, 0.5))
ax[0].set_yticks([0.1,1,10,100])
ax[0].set_ylabel('Distance\n[km]')

pd_fold.plot(x='types',
        kind='bar',
        stacked=False,
        ax=ax[1],grid=True,logy=True,rot=0,color=cl_list[0::2],
        y=['aod_fine','GOCI_AOD_f','MERRA_AOD_bc','MERRA_AOD_oc','MERRA_AOD_sulf'])
ax[1].legend(['4STAR Fine','GOCI Fine','MERRA Black Carbon','MERRA Org. Carbon','MERRA Sulfate'],
             frameon=False,loc='center left', bbox_to_anchor=(1.0, 0.5))
ax[1].set_yticks([0.1,1,10,100])
ax[1].set_ylabel('Distance\n[km]')

#ax[0].legend(['GOCI Coarse','GOCI Fine','MERRA Black Carbon','MERRA Dust',
#           'MERRA Org. Carbon','MERRA Sea salt','MERRA Sulfate','4STAR Coarse','4STAR Fine'],
#          frameon=False,loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.tight_layout(h_pad=0.5,w_pad=-10)
plt.savefig(fp+'plot/KORUS_autocorr_percentile_speciated_AOD_{}.png'.format(vv),dpi=600,transparent=True)


# In[358]:


lbls = ['GOCI Coarse','GOCI Fine','MERRA Black Carbon','MERRA Dust',
           'MERRA Org. Carbon','MERRA Sea salt','MERRA Sulfate','4STAR Coarse','4STAR Fine']


# ### Add the average AODs for each species

# In[420]:


mean_vals = {}
median_vals = {}
std_vals = {}
for kv in vals.keys():
    mean_vals[kv],median_vals[kv],std_vals[kv] = [],[],[]
    for n,tp in enumerate(itypes):
        dv_temp = np.hstack(dvals[kv][tp])
        mean_vals[kv].append(np.nanmean(dv_temp))
        median_vals[kv].append(np.nanmedian(dv_temp))
        std_vals[kv].append(np.nanstd(dv_temp))


# In[421]:


vals.keys()


# In[259]:


def match_ygrid(ax1,ax2,ticks):
    'function to match the grid ticks to a dual y axis plot, matching limits of ax2 such that the ticks line up with ax1 grid'
    y0,y1 = ax1.get_ybound()
    ti = ax1.get_yticks()
    ax2.set_yticks(ticks)
    if ax1.get_yscale() =='log':
        a = (np.log10(ti[1])-np.log10(ti[0]))/(ticks[1]-ticks[0])
        dy = a*ticks[0]-np.log10(ti[0])
        ax2.set_ylim((np.log10(y0)+dy)/a,(np.log10(y1)+dy)/a)
    else:
        a = (ti[1]-ti[0])/(ticks[1]-ticks[0])
        dy = a*ticks[0]-ti[0]
        ax2.set_ylim((y0+dy)/a,(y1+dy)/a)

    


# In[714]:


fig, ax = plt.subplots(2,1,sharex=True,figsize=(10,3))

y1_lbl = ['aod_coarse','GOCI_AOD_c','MERRA_AOD_dust','MERRA_AOD_sea']
b1 = pd_fold.plot(x='types',
        kind='bar',
        stacked=False,
        title='Autocorrelation distance at 85th percentile for speciated AOD',
        ax=ax[0],grid=True,logy=True,rot=0,color=cl_list[1::2],
        y=y1_lbl)

b1x = np.reshape([rect.get_x()+rect.get_width()/2.0 for rect in b1.patches],(len(y1_lbl),len(types)))
b1c = [rect.get_facecolor() for ir,rect in enumerate(b1.patches) if ir%len(types)==0]
b1x2 = b1.twinx()
for iy, lbl in enumerate(y1_lbl):
    #b1x2.plot(b1x[iy,:],mean_vals[lbl],'^',color='k',markersize=8)
    b1x2.errorbar(b1x[iy,:],mean_vals[lbl],yerr=std_vals[lbl],marker='+',color='k',alpha=0.5,ls='None')
    
b1x2.set_ylabel('AOD')
match_ygrid(b1,b1x2,np.linspace(0,0.45,4))

ax[0].legend(['4STAR Coarse','GOCI Coarse','MERRA Dust','MERRA Sea salt'],
             frameon=False,loc='center left', bbox_to_anchor=(1.085, 0.5))
ax[0].set_yticks([0.1,1,10,100])
ax[0].set_ylabel('Distance\n[km]')

y2_lbl = ['aod_fine','GOCI_AOD_f','MERRA_AOD_bc','MERRA_AOD_oc','MERRA_AOD_sulf']
b2 = pd_fold.plot(x='types',
        kind='bar',
        stacked=False,
        ax=ax[1],grid=True,logy=True,rot=0,color=cl_list[0::2],
        y=y2_lbl)


b2x = np.reshape([rect.get_x()+rect.get_width()/2.0 for rect in b2.patches],(len(y2_lbl),len(types)))
b2c = [rect.get_facecolor() for ir,rect in enumerate(b2.patches) if ir%len(types)==0]
b2x2 = b2.twinx()
for iy, lbl in enumerate(y2_lbl):
    b2x2.errorbar(b2x[iy,:],mean_vals[lbl],yerr=std_vals[lbl],marker='+',color='k',alpha=0.8,ls='None')
    
b2x2.set_ylabel('AOD')
match_ygrid(b2,b2x2,np.linspace(0,0.45,4))
ax[1].legend(['4STAR Fine','GOCI Fine','MERRA Black Carbon','MERRA Org. Carbon','MERRA Sulfate'],
             frameon=False,loc='center left', bbox_to_anchor=(1.085, 0.5))
ax[1].set_yticks([0.1,1,10,100])
ax[1].set_ylabel('Distance\n[km]')
ax[1].set_xlabel('')
#ax[0].legend(['GOCI Coarse','GOCI Fine','MERRA Black Carbon','MERRA Dust',
#           'MERRA Org. Carbon','MERRA Sea salt','MERRA Sulfate','4STAR Coarse','4STAR Fine'],
#          frameon=False,loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.tight_layout(h_pad=0.5,w_pad=-12)
match_ygrid(b2,b2x2,np.linspace(0,0.45,4))
match_ygrid(b1,b1x2,np.linspace(0,0.45,4))
plt.savefig(fp+'plot/KORUS_autocorr_percentile_speciated_AOD_withAODavgs_{}.png'.format(vv),dpi=600,transparent=True)


# In[722]:


np.array(mean_vals['MERRA_AOD_bc'])+np.array(mean_vals['MERRA_AOD_oc'])+np.array(mean_vals['MERRA_AOD_sulf']),np.array(mean_vals['MERRA_AOD_bc'])+np.array(mean_vals['MERRA_AOD_oc'])+np.array(mean_vals['MERRA_AOD_sulf'])          -np.array(mean_vals['aod_fine'])


# In[721]:


np.array(mean_vals['MERRA_AOD_dust'])+np.array(mean_vals['MERRA_AOD_sea']),np.array(mean_vals['MERRA_AOD_dust'])+np.array(mean_vals['MERRA_AOD_sea'])-np.array(mean_vals['aod_coarse'])


# In[723]:


pd_fold


# ### Redo 85th percentil plot with segment monte carlo (v3)

# In[422]:


ka = autocorr_mc.keys()
ka.sort()


# In[423]:


mc = True


# In[424]:


e_fold = {}
perc = {}
for k in types:
    e_fold[k] = []
    perc[k] = []
jj = 0
for j,t in enumerate(types):
    print 'Type {}: label {}'.format(t,legend_list[j])
    for k in kas:
        print '..{}: efolding={}, 90%={}, 85%={}'.format(k,
        e_folding_autocorr(np.nanmean(autocorr_mc[k][j,jj:,:],axis=1)),
        percentile_autocorr(np.nanmean(autocorr_mc[k][j,jj:,:],axis=1),0.9),percentile_autocorr(np.nanmean(autocorr_mc[k][j,jj:,:],axis=1),0.85))
        if mc:
            e_fold[t].append([e_folding_autocorr(a) for a in autocorr_mc[k][j,jj:,:].T])
            perc[t].append([percentile_autocorr(a,0.85) for a in autocorr_mc[k][j,jj:,:].T])
        else:
            e_fold[t].append(e_folding_autocorr(np.nanmean(autocorr_mc[k][j,jj:,:],axis=1)))
            perc[t].append(percentile_autocorr(np.nanmean(autocorr_mc[k][j,jj:,:],axis=1),0.85))


# In[425]:


if mc:
    e_foldr = [list(np.nanmedian(e_fold[k],axis=1)) for k in types]
    percr = [list(np.nanmedian(perc[k],axis=1)) for k in types]
else:
    e_foldr = [e_fold[k] for k in types]
    percr = [perc[k] for k in types]


# In[426]:


for i,e in enumerate(e_foldr):
    e.insert(0,legend_list[i])
    percr[i].insert(0,legend_list[i])


# In[427]:


pd_fold = pd.DataFrame(percr,columns=['types'] + kas) 


# In[428]:


autocorr_mc.keys()


# In[429]:


pd_fold['types'][3] = 'Extreme\npollution'


# In[303]:


fig, ax = plt.subplots(2,1,sharex=True,figsize=(10,3))

y1_lbl = ['aod_coarse','GOCI_AOD_c','MERRA_AOD_dust','MERRA_AOD_sea']
b1 = pd_fold.plot(x='types',
        kind='bar',
        stacked=False,
        title='Autocorrelation distance at 85th percentile for speciated AOD',
        ax=ax[0],grid=True,logy=True,rot=0,color=cl_list[1::2],
        y=y1_lbl)

b1x = np.reshape([rect.get_x()+rect.get_width()/2.0 for rect in b1.patches],(len(y1_lbl),len(types)))
b1c = [rect.get_facecolor() for ir,rect in enumerate(b1.patches) if ir%len(types)==0]
b1x2 = b1.twinx()
for iy, lbl in enumerate(y1_lbl):
    #b1x2.plot(b1x[iy,:],mean_vals[lbl],'^',color='k',markersize=8)
    b1x2.errorbar(b1x[iy,:],mean_vals[lbl],yerr=std_vals[lbl],marker='+',color='k',alpha=0.5,ls='None')
    
b1x2.set_ylabel('Coarse AOD')
match_ygrid(b1,b1x2,np.linspace(0,0.45,4))

ax[0].legend(['4STAR Coarse','GOCI Coarse','MERRA Dust','MERRA Sea salt'],
             frameon=False,loc='center left', bbox_to_anchor=(1.085, 0.5))
ax[0].set_yticks([0.1,1,10,100])
ax[0].set_ylabel('Distance\n[km]')

y2_lbl = ['aod_fine','GOCI_AOD_f','MERRA_AOD_bc','MERRA_AOD_oc','MERRA_AOD_sulf']
b2 = pd_fold.plot(x='types',
        kind='bar',
        stacked=False,
        ax=ax[1],grid=True,logy=True,rot=0,color=cl_list[0::2],
        y=y2_lbl)


b2x = np.reshape([rect.get_x()+rect.get_width()/2.0 for rect in b2.patches],(len(y2_lbl),len(types)))
b2c = [rect.get_facecolor() for ir,rect in enumerate(b2.patches) if ir%len(types)==0]
b2x2 = b2.twinx()
for iy, lbl in enumerate(y2_lbl):
    b2x2.errorbar(b2x[iy,:],mean_vals[lbl],yerr=std_vals[lbl],marker='+',color='k',alpha=0.8,ls='None')
    
b2x2.set_ylabel('Fine AOD')
match_ygrid(b2,b2x2,np.linspace(0,0.45,4))
ax[1].legend(['4STAR Fine','GOCI Fine','MERRA Black Carbon','MERRA Org. Carbon','MERRA Sulfate'],
             frameon=False,loc='center left', bbox_to_anchor=(1.085, 0.5))
ax[1].set_yticks([0.1,1,10,100])
ax[1].set_ylabel('Distance\n[km]')
ax[1].set_xlabel('')
#ax[0].legend(['GOCI Coarse','GOCI Fine','MERRA Black Carbon','MERRA Dust',
#           'MERRA Org. Carbon','MERRA Sea salt','MERRA Sulfate','4STAR Coarse','4STAR Fine'],
#          frameon=False,loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.tight_layout(h_pad=0.5,w_pad=-12)
match_ygrid(b2,b2x2,np.linspace(0,0.45,4))
match_ygrid(b1,b1x2,np.linspace(0,0.45,4))
plt.savefig(fp+'plot/KORUS_autocorr_percentile_speciated_AOD_withAODavgs_{}.png'.format(vv),dpi=600,transparent=True)


# In[432]:


fig, ax = plt.subplots(2,1,sharex=True,figsize=(8,3))

y1_lbl = ['aod_coarse','GOCI_AOD_c','MERRA_AOD_dust','MERRA_AOD_sea']
b1 = pd_fold.iloc[0:5].plot(x='types',
        kind='bar',
        stacked=False,
        title='Autocorrelation distance at 85th percentile for speciated AOD',
        ax=ax[0],grid=True,logy=True,rot=0,color=cl_list[1::2],
        y=y1_lbl)

b1x = np.reshape([rect.get_x()+rect.get_width()/2.0 for rect in b1.patches],(len(y1_lbl),len(types[0:5])))
b1c = [rect.get_facecolor() for ir,rect in enumerate(b1.patches) if ir%len(types)==0]
b1x2 = b1.twinx()
for iy, lbl in enumerate(y1_lbl):
    #b1x2.plot(b1x[iy,:],mean_vals[lbl],'^',color='k',markersize=8)
    b1x2.errorbar(b1x[iy,:],mean_vals[lbl][0:5],yerr=std_vals[lbl][0:5],marker='+',color='k',alpha=0.5,ls='None')
    
b1x2.set_ylabel('Coarse AOD')
match_ygrid(b1,b1x2,np.linspace(0,0.45,4))

ax[0].legend(['4STAR Coarse','GOCI Coarse','MERRA Dust','MERRA Sea salt'],
             frameon=False,loc='center left', bbox_to_anchor=(1.12, 0.5))
ax[0].set_yticks([1,10,100])
ax[0].set_ylabel('Distance\n[km]')

y2_lbl = ['aod_fine','GOCI_AOD_f','MERRA_AOD_bc','MERRA_AOD_oc','MERRA_AOD_sulf']
b2 = pd_fold.iloc[0:5].plot(x='types',
        kind='bar',
        stacked=False,
        ax=ax[1],grid=True,logy=True,rot=0,color=cl_list[0::2],
        y=y2_lbl)


b2x = np.reshape([rect.get_x()+rect.get_width()/2.0 for rect in b2.patches],(len(y2_lbl),len(types[0:5])))
b2c = [rect.get_facecolor() for ir,rect in enumerate(b2.patches) if ir%len(types)==0]
b2x2 = b2.twinx()
for iy, lbl in enumerate(y2_lbl):
    b2x2.errorbar(b2x[iy,:],mean_vals[lbl][0:5],yerr=std_vals[lbl][0:5],marker='+',color='k',alpha=0.8,ls='None')
    
b2x2.set_ylabel('Fine AOD')
match_ygrid(b2,b2x2,np.linspace(0,0.45,4))
ax[1].legend(['4STAR Fine','GOCI Fine','MERRA Black Carbon','MERRA Org. Carbon','MERRA Sulfate'],
             frameon=False,loc='center left', bbox_to_anchor=(1.12, 0.5))
ax[1].set_yticks([1,10,100])
ax[1].set_ylabel('Distance\n[km]')
ax[1].set_xlabel('')
#ax[0].legend(['GOCI Coarse','GOCI Fine','MERRA Black Carbon','MERRA Dust',
#           'MERRA Org. Carbon','MERRA Sea salt','MERRA Sulfate','4STAR Coarse','4STAR Fine'],
#          frameon=False,loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.tight_layout(h_pad=0.5,w_pad=-12)
match_ygrid(b2,b2x2,np.linspace(0,0.45,4))
match_ygrid(b1,b1x2,np.linspace(0,0.45,4))
plt.savefig(fp+'plot/KORUS_autocorr_percentile_speciated_AOD_withAODavgs_mc_{}.png'.format(vv),dpi=600,transparent=True)


# In[304]:


vv


# ### Save autocorr and monte carlo calculations to file

# In[191]:


vv = 'v2'


# In[382]:


dat_c = { u'mean_vals':mean_vals,u'std_vals':std_vals,u'dvals':dvals,
        u'vals':vals,u'types':types,
        u'itypes':itypes,u'corr_ks':corr_ks,u'autocorr_mean':autocorr_mean}


# In[737]:


import write_utils as wu
dat_u = wu.iterate_dict_unicode(dat_c)


# In[738]:


hs.savemat(fp+'KORUS_fine_coarse_autocorr_dvals_{}.mat'.format(vv),dat_u)


# ### intrinsic vs. extrinsic subplots

# In[368]:


fig, ax = plt.subplots(1,2,sharey=True,figsize=(8,4))
j = 0
jj = 0
cl = ['tab:blue','tab:orange','tab:green','tab:red']
#ax[0].plot(corr_ks[jj:],autocorr_mean['AE'][j,jj:],ls='--',marker='.',label='AE',color='k')
ax[1].errorbar(corr_ks[jj:],autocorr_mean['AE'][j,jj:],yerr=autocorr_d['AE'][j,jj:],
               ls='--',marker='.',label='4STAR',color=cl[0])
ax[1].axvline(percentile_autocorr(autocorr_mean['AE'][j,jj:],0.85),ls=':',color=cl[0])
ax[1].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AE'][j,jj:],yerr=autocorr_d['MERRA_AE'][j,jj:],
               ls='--',marker='.',label='MERRA',color=cl[1])
ax[1].axvline(percentile_autocorr(autocorr_mean['MERRA_AE'][j,jj:],0.85),ls=':',color=cl[1])
ax[1].errorbar(corr_ks[jj:],autocorr_mean['GOCI_AE'][j,jj:],yerr=autocorr_d['GOCI_AE'][j,jj:],
               ls='--',marker='.',label='GOCI',color=cl[2])
ax[1].axvline(percentile_autocorr(autocorr_mean['GOCI_AE'][j,jj:],0.85),ls=':',color=cl[2])
ax[1].errorbar(corr_ks[jj:],autocorr_mean['situ_ae'][j,jj:],yerr=autocorr_d['situ_ae'][j,jj:],
               ls='--',marker='.',label='Point / in situ',color=cl[3])
ax[1].axvline(percentile_autocorr(autocorr_mean['situ_ae'][j,jj:],0.85),ls=':',color='k',label='85th percentile')
ax[1].axvline(percentile_autocorr(autocorr_mean['situ_ae'][j,jj:],0.85),ls=':',color=cl[3])
ax[1].set_title('AE')

#ax[1].plot(corr_ks[jj:],autocorr_mean['aod0500'][j,jj:],ls='-',marker='.',label='AOD',color='k')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['aod0500'][j,jj:],yerr=autocorr_d['aod0500'][j,jj:],
               ls='-',marker='.',label='4STAR')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD'][j,jj:],yerr=autocorr_d['MERRA_AOD'][j,jj:],
               ls='-',marker='.',label='MERRA')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['GOCI_AOD'][j,jj:],yerr=autocorr_d['GOCI_AOD'][j,jj:],
               ls='-',marker='.',label='GOCI')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['situ_ext'][j,jj:],yerr=autocorr_d['situ_ext'][j,jj:],
               ls='-',marker='.',label='Point / in situ')
ax[0].axvline(percentile_autocorr(autocorr_mean['aod0500'][j,jj:],0.85),ls=':',color=cl[0],lw=1.2)
ax[0].axvline(percentile_autocorr(autocorr_mean['MERRA_AOD'][j,jj:],0.85),ls=':',color=cl[1])
ax[0].axvline(percentile_autocorr(autocorr_mean['GOCI_AOD'][j,jj:],0.85),ls=':',color=cl[2])
ax[0].axvline(percentile_autocorr(autocorr_mean['situ_ext'][j,jj:],0.85),ls=':',color='k',label='85th percentile')
ax[0].axvline(percentile_autocorr(autocorr_mean['situ_ext'][j,jj:],0.85),ls=':',color=cl[3])
ax[0].set_title('AOD$_{{500}}$')

if False:
    #ax[1].plot(corr_ks[jj:],autocorr_mean['aod0500'][j,jj:],ls='-',marker='.',label='AOD',color='k')
    ax[2].errorbar(corr_ks[jj:],autocorr_mean['fmf'][j,jj:],yerr=autocorr_d['fmf'][j,jj:],ls='-',marker='.',label='4STAR',color='tab:blue')
    #ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD'][j,jj:],yerr=autocorr_d['MERRA_AOD'][j,jj:],ls='-',marker='.',label='MERRA')
    ax[2].errorbar(corr_ks[jj:],autocorr_mean['GOCI_fmf'][j,jj:],yerr=autocorr_d['GOCI_fmf'][j,jj:],ls='-',marker='.',label='GOCI',color='tab:green')
    #ax[0].errorbar(corr_ks[jj:],autocorr_mean['situ_ext'][j,jj:],yerr=autocorr_d['situ_ext'][j,jj:],ls='-',marker='.',label='Point / in situ')
    ax[2].set_title('FMF')

ax[0].legend(frameon=False)

ax[0].set_ylim(0,1)
ax[0].set_yticks([0,0.25,0.5,0.75,1.0])
ax[0].set_xscale('log')
ax[1].set_xscale('log')
#ax[2].set_xscale('log')
ax[0].grid()
ax[1].grid()
#ax[2].grid()

plt.savefig(fp+'plot/KORUS_autocorr_AOD_AE_{}.png'.format(vv),dpi=600,transparent=True)


# In[65]:


fig, ax = plt.subplots(1,2,sharey=True,figsize=(8,4))
j = 6
jj = 0
#ax[0].plot(corr_ks[jj:],autocorr_mean['AE'][j,jj:],ls='--',marker='.',label='AE',color='k')
ax[1].errorbar(corr_ks[jj:],autocorr_mean['AE'][j,jj:],yerr=autocorr_d['AE'][j,jj:],ls='--',marker='.',label='4STAR')
ax[1].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AE'][j,jj:],yerr=autocorr_d['MERRA_AE'][j,jj:],ls='--',marker='.',label='MERRA')
ax[1].errorbar(corr_ks[jj:],autocorr_mean['GOCI_AE'][j,jj:],yerr=autocorr_d['GOCI_AE'][j,jj:],ls='--',marker='.',label='GOCI')
ax[1].errorbar(corr_ks[jj:],autocorr_mean['situ_ae'][j,jj:],yerr=autocorr_d['situ_ae'][j,jj:],ls='--',marker='.',label='Point / in situ')
ax[1].set_title('AE')

#ax[1].plot(corr_ks[jj:],autocorr_mean['aod0500'][j,jj:],ls='-',marker='.',label='AOD',color='k')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['aod0500'][j,jj:],yerr=autocorr_d['aod0500'][j,jj:],ls='-',marker='.',label='4STAR')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD'][j,jj:],yerr=autocorr_d['MERRA_AOD'][j,jj:],ls='-',marker='.',label='MERRA')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['GOCI_AOD'][j,jj:],yerr=autocorr_d['GOCI_AOD'][j,jj:],ls='-',marker='.',label='GOCI')
ax[0].errorbar(corr_ks[jj:],autocorr_mean['situ_ext'][j,jj:],yerr=autocorr_d['situ_ext'][j,jj:],ls='-',marker='.',label='Point / in situ')
ax[0].set_title('AOD$_{{500}}$')

if False:
    #ax[1].plot(corr_ks[jj:],autocorr_mean['aod0500'][j,jj:],ls='-',marker='.',label='AOD',color='k')
    ax[2].errorbar(corr_ks[jj:],autocorr_mean['fmf'][j,jj:],yerr=autocorr_d['fmf'][j,jj:],ls='-',marker='.',label='4STAR',color='tab:blue')
    #ax[0].errorbar(corr_ks[jj:],autocorr_mean['MERRA_AOD'][j,jj:],yerr=autocorr_d['MERRA_AOD'][j,jj:],ls='-',marker='.',label='MERRA')
    ax[2].errorbar(corr_ks[jj:],autocorr_mean['GOCI_fmf'][j,jj:],yerr=autocorr_d['GOCI_fmf'][j,jj:],ls='-',marker='.',label='GOCI',color='tab:green')
    #ax[0].errorbar(corr_ks[jj:],autocorr_mean['situ_ext'][j,jj:],yerr=autocorr_d['situ_ext'][j,jj:],ls='-',marker='.',label='Point / in situ')
    ax[2].set_title('FMF')

ax[0].legend(frameon=False)

ax[0].set_ylim(0,1)
ax[0].set_yticks([0,0.25,0.5,0.75,1.0])
ax[0].set_xscale('log')
ax[1].set_xscale('log')
#ax[2].set_xscale('log')
ax[0].grid()
ax[1].grid()
#ax[2].grid()

#plt.savefig(fp+'plot/KORUS_autocorr_AOD_AE_{}.png'.format(vv),dpi=600,transparent=True)


# ### For v3 instrisic vs. extrinsic (monte carlo on segments instead of samples)

# In[225]:


corr_ks = np.array(corr_ks)


# In[467]:


corr_ks[11:]


# In[477]:


fig, ax = plt.subplots(1,2,sharey=True,figsize=(8,4))
j = 0
jj = 0
jjm = 10
cl = ['tab:blue','tab:orange','tab:green','tab:red']
#ax[0].plot(corr_ks[jj:],np.nanmean(autocorr_mc['AE'][j,jj:],ls='--',marker='.',label='AE',color='k')
ax[1].errorbar(corr_ks[jj:],np.nanmean(autocorr_mc['AE'][j,jj:,:],axis=1),yerr=np.nanstd(autocorr_mc['AE'][j,jj:,:],axis=1),
               elinewidth=0.51,ls='--',marker='.',label='4STAR',color=cl[0])
ax[1].axvline(percentile_autocorr(np.nanmean(autocorr_mc['AE'][j,jj:,:],axis=1),0.85)*1.03,ls=':',color=cl[0])
ax[1].plot(corr_ks[jj:]*1.02,np.nanmean(autocorr_mc['MERRA_AE'][j,jj:,:],axis=1),ls='--',marker='.',alpha=0.3,color=cl[1])
ax[1].errorbar(corr_ks[jjm:]*1.02,np.nanmean(autocorr_mc['MERRA_AE'][j,jjm:,:],axis=1),
               yerr=np.nanstd(autocorr_mc['MERRA_AE'][j,jjm:,:],axis=1),elinewidth=0.51,
               ls='--',marker='.',label='MERRA',color=cl[1])
ax[1].axvline(percentile_autocorr(np.nanmean(autocorr_mc['MERRA_AE'][j,jj:,:],axis=1),0.85),ls=':',color=cl[1])
ax[1].plot(corr_ks[jj:],np.nanmean(autocorr_mc['GOCI_AE'][j,jj:,:],axis=1),ls='--',marker='.',alpha=0.3,color=cl[2])
ax[1].errorbar(corr_ks[jjm:],np.nanmean(autocorr_mc['GOCI_AE'][j,jjm:,:],axis=1),
               yerr=np.nanstd(autocorr_mc['GOCI_AE'][j,jjm:,:],axis=1),elinewidth=0.51,
               ls='--',marker='.',label='GOCI',color=cl[2])
ax[1].axvline(percentile_autocorr(np.nanmean(autocorr_mc['GOCI_AE'][j,jj:,:],axis=1),0.85)*0.97,ls=':',color=cl[2])
ax[1].errorbar(corr_ks[jj:]*1.02,np.nanmean(autocorr_mc['situ_ae'][j,jj:,:],axis=1),
               yerr=np.nanstd(autocorr_mc['situ_ae'][j,jj:,:],axis=1),elinewidth=0.51,
               ls='--',marker='.',label='Point / in situ',color=cl[3])
ax[1].axvline(percentile_autocorr(np.nanmean(autocorr_mc['situ_ae'][j,jj:,:],axis=1),0.85),ls=':',color='k',label='85th percentile')
ax[1].axvline(percentile_autocorr(np.nanmean(autocorr_mc['situ_ae'][j,jj:,:],axis=1),0.85),ls=':',color=cl[3])
ax[1].set_title('AE')

#ax[1].plot(corr_ks[jj:],np.nanmean(autocorr_mc['aod0500'][j,jj:],ls='-',marker='.',label='AOD',color='k')
ax[0].errorbar(corr_ks[jj:],np.nanmean(autocorr_mc['aod0500'][j,jj:,:],axis=1),
               yerr=np.nanstd(autocorr_mc['aod0500'][j,jj:,:],axis=1),elinewidth=0.51,
               ls='-',marker='.',label='4STAR',color=cl[0])
ax[0].plot(corr_ks[jj:]*1.02,np.nanmean(autocorr_mc['MERRA_AOD'][j,jj:,:],axis=1),ls='-',marker='.',alpha=0.3,color=cl[1])
ax[0].errorbar(corr_ks[jjm:]*1.02,np.nanmean(autocorr_mc['MERRA_AOD'][j,jjm:,:],axis=1),
               yerr=np.nanstd(autocorr_mc['MERRA_AOD'][j,jjm:,:],axis=1),elinewidth=0.51,
               ls='-',marker='.',label='MERRA',color=cl[1])
ax[0].plot(corr_ks[jj:],np.nanmean(autocorr_mc['GOCI_AOD'][j,jj:,:],axis=1),ls='-',marker='.',alpha=0.3,color=cl[2])
ax[0].errorbar(corr_ks[jjm:],np.nanmean(autocorr_mc['GOCI_AOD'][j,jjm:,:],axis=1),
               yerr=np.nanstd(autocorr_mc['GOCI_AOD'][j,jjm:,:],axis=1),elinewidth=0.51,
               ls='-',marker='.',label='GOCI',color=cl[2])
ax[0].errorbar(corr_ks[jj:]*1.02,np.nanmean(autocorr_mc['situ_ext'][j,jj:,:],axis=1),
               yerr=np.nanstd(autocorr_mc['situ_ext'][j,jj:,:],axis=1),elinewidth=0.51,
               ls='-',marker='.',label='Point / in situ',color=cl[3])
ax[0].axvline(percentile_autocorr(np.nanmean(autocorr_mc['aod0500'][j,jj:,:],axis=1),0.85),ls=':',color=cl[0],lw=1.2)
ax[0].axvline(percentile_autocorr(np.nanmean(autocorr_mc['MERRA_AOD'][j,jj:,:],axis=1),0.85)*1.03,ls=':',color=cl[1])
ax[0].axvline(percentile_autocorr(np.nanmean(autocorr_mc['GOCI_AOD'][j,jj:,:],axis=1),0.85),ls=':',color=cl[2])
ax[0].axvline(percentile_autocorr(np.nanmean(autocorr_mc['situ_ext'][j,jj:,:],axis=1),0.85),ls=':',color='k',label='85th percentile')
ax[0].axvline(percentile_autocorr(np.nanmean(autocorr_mc['situ_ext'][j,jj:,:],axis=1),0.85),ls=':',color=cl[3])
ax[0].set_title('AOD$_{{500}}$')

if False:
    #ax[1].plot(corr_ks[jj:],np.nanmean(autocorr_mc['aod0500'][j,jj:],ls='-',marker='.',label='AOD',color='k')
    ax[2].errorbar(corr_ks[jj:],np.nanmean(autocorr_mc['fmf'][j,jj:,:],axis=1),yerr=np.nanstd(autocorr_mc['fmf'][j,jj:,:],axis=1),ls='-',marker='.',label='4STAR',color='tab:blue')
    #ax[0].errorbar(corr_ks[jj:],np.nanmean(autocorr_mc['MERRA_AOD'][j,jj:],yerr=np.nanstd(autocorr_mc['MERRA_AOD'][j,jj:],ls='-',marker='.',label='MERRA')
    ax[2].errorbar(corr_ks[jj:],np.nanmean(autocorr_mc['GOCI_fmf'][j,jj:,:],axis=1),yerr=np.nanstd(autocorr_mc['GOCI_fmf'][j,jj:,:],axis=1),ls='-',marker='.',label='GOCI',color='tab:green')
    #ax[0].errorbar(corr_ks[jj:],np.nanmean(autocorr_mc['situ_ext'][j,jj:],yerr=np.nanstd(autocorr_mc['situ_ext'][j,jj:],ls='-',marker='.',label='Point / in situ')
    ax[2].set_title('FMF')

ax[0].legend(frameon=False)

ax[0].set_ylim(0,1)
ax[0].set_yticks([0,0.25,0.5,0.75,1.0])
ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[0].set_xlabel('distance [km]')
ax[1].set_xlabel('distance [km]')
ax[0].set_ylabel('Autocorrelation')
ax[0].grid()
ax[1].grid()
#ax[2].grid()

plt.savefig(fp+'plot/KORUS_autocorr_AOD_AE_{}_sub.png'.format(vv),dpi=600,transparent=True)


# In[474]:


percentile_autocorr(np.nanmean(autocorr_mc['GOCI_AOD'][j,jj:,:],axis=1),0.85)


# In[315]:


j=0
jj=0
for k in ['GOCI_AOD','MERRA_AOD','aod0500','situ_ext']:
    print k, np.nanmean(np.nanstd(autocorr_mc[k][j,jj:-1,:],axis=1))
for k in ['GOCI_AE','MERRA_AE','AE','situ_ae']:
    print k, np.nanmean(np.nanstd(autocorr_mc[k][j,jj:-1,:],axis=1))


# ### Larger plots

# In[359]:


plt.figure(figsize=(6,6))
j = 1
jj = 0
plt.plot(corr_ks[jj:],autocorr['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],label='Point/In situ AE',color='tab:blue',ls='--')
plt.plot(corr_ks[jj:],autocorr['AE'][j,jj:]/autocorr['AE'][j,jj],label='Column/4STAR AE',color='k',ls='--')
plt.plot(corr_ks[jj:],autocorr['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],label='Column/GOCI AE',color='tab:red',ls='--')
plt.plot(corr_ks[jj:],autocorr['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],label='Column/MERRA AE',color='tab:purple',ls='--')

plt.plot(corr_ks[jj:],autocorr['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],label='Point/In situ SSA',color='tab:orange',ls=':')
plt.plot(corr_ks[jj:],autocorr['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],label='Point/In situ Ext.',color='tab:blue',ls='-')

plt.plot(corr_ks[jj:],autocorr['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],label='Column/4STAR AOD',color='k')
plt.plot(corr_ks[jj:],autocorr['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],label='Column/GOCI AOD',color='tab:red')
plt.plot(corr_ks[jj:],autocorr['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],label='Column/MERRA AOD',color='tab:purple')


ax = plt.gca()

ax.errorbar(corr_ks[jj:],autocorr['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],
                 yerr=autocorr_d['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],color='tab:orange',ls=':',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],
                 yerr=autocorr_d['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],color='tab:blue',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],
                 yerr=autocorr_d['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],color='tab:blue',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['AE'][j,jj:]/autocorr['AE'][j,jj],
                 yerr=autocorr_d['AE'][j,jj:]/autocorr['AE'][j,jj],color='k',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],
                 yerr=autocorr_d['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],color='tab:red',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],
                 yerr=autocorr_d['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],color='tab:purple',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],
                 yerr=autocorr_d['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],color='k',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],
                 yerr=autocorr_d['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],color='tab:red',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],
                 yerr=autocorr_d['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],color='tab:purple',
                 alpha=0.6,capsize=0)

ax.set_ylim(0,1)
ax.set_yticks([0,0.25,0.5,0.75,1.0])
ax.set_xscale('linear')
ax.grid()
plt.legend(frameon=True,ncol=2,loc=(0,1.04))
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Autocorrelation')
plt.subplots_adjust(top=0.7)
plt.savefig(fp+'plot/KORUS_Autocorr_intrisic_vs_extrinsic_{}.png'.format(vv),dpi=600,transparent=True)


# In[355]:


plt.figure(figsize=(6,6))
j = 0
jj = 1
plt.plot(corr_ks[jj:],autocorr_mean['situ_ae'][j,jj:]/autocorr_mean['situ_ae'][j,jj],label='Point/In situ AE',color='tab:blue',ls='--')
plt.plot(corr_ks[jj:],autocorr_mean['AE'][j,jj:]/autocorr_mean['AE'][j,jj],label='Column/4STAR AE',color='k',ls='--')
plt.plot(corr_ks[jj:],autocorr_mean['GOCI_AE'][j,jj:]/autocorr_mean['GOCI_AE'][j,jj],label='Column/GOCI AE',color='tab:red',ls='--')
plt.plot(corr_ks[jj:],autocorr_mean['MERRA_AE'][j,jj:]/autocorr_mean['MERRA_AE'][j,jj],label='Column/MERRA AE',color='tab:purple',ls='--')

plt.plot(corr_ks[jj:],autocorr_mean['situ_ssa'][j,jj:]/autocorr_mean['situ_ssa'][j,jj],label='Point/In situ SSA',color='tab:orange',ls=':')
plt.plot(corr_ks[jj:],autocorr_mean['situ_ext'][j,jj:]/autocorr_mean['situ_ext'][j,jj],label='Point/In situ Ext.',color='tab:blue',ls='-')

plt.plot(corr_ks[jj:],autocorr_mean['aod0500'][j,jj:]/autocorr_mean['aod0500'][j,jj],label='Column/4STAR AOD',color='k')
plt.plot(corr_ks[jj:],autocorr_mean['GOCI_AOD'][j,jj:]/autocorr_mean['GOCI_AOD'][j,jj],label='Column/GOCI AOD',color='tab:red')
plt.plot(corr_ks[jj:],autocorr_mean['MERRA_AOD'][j,jj:]/autocorr_mean['MERRA_AOD'][j,jj],label='Column/MERRA AOD',color='tab:purple')


ax = plt.gca()
if False:
    ax.errorbar(corr_ks[jj:],autocorr['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],
                     yerr=autocorr_d['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],color='tab:orange',ls=':',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],
                     yerr=autocorr_d['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],color='tab:blue',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],
                     yerr=autocorr_d['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],color='tab:blue',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['AE'][j,jj:]/autocorr['AE'][j,jj],
                     yerr=autocorr_d['AE'][j,jj:]/autocorr['AE'][j,jj],color='k',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],
                     yerr=autocorr_d['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],color='tab:red',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],
                     yerr=autocorr_d['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],color='tab:purple',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],
                     yerr=autocorr_d['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],color='k',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],
                     yerr=autocorr_d['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],color='tab:red',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],
                     yerr=autocorr_d['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],color='tab:purple',
                     alpha=0.6,capsize=0)

ax.set_ylim(0,1)
ax.set_yticks([0,0.25,0.5,0.75,1.0])
ax.set_xscale('linear')
ax.grid()
plt.legend(frameon=True,ncol=2,loc=(0,1.04))
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Autocorrelation')
plt.subplots_adjust(top=0.7)
#plt.savefig(fp+'plot/KORUS_Autocorr_intrisic_vs_extrinsic.png',dpi=600,transparent=True)


# # Make supporting plots of AOD comparisons and spatial maps

# ## Plot the comparison of GOCI to 4STAR AOD

# In[258]:


plt.figure()
plt.hist([ar['GPS_Alt'],ar['GPS_Alt'][fl]],label=['All data','QA data'],bins=50)
plt.legend()
plt.xlabel('Altitude [m]')
plt.ylabel('number of samples')
plt.title('KORUS-AQ 4STAR sampling')


# In[259]:


plt.figure()
fla = ar['fl_QA'] & (ar['GPS_Alt']<500.0)
flan = ar['fl_QA'] & (ar['GPS_Alt']<500.0) & np.isfinite(ar['AOD0501']) & np.isfinite(goci2ar['aod'])
r = np.corrcoef(ar['AOD0501'][flan],goci2ar['aod'][flan])[0,1]**2.0
plt.plot(ar['AOD0501'][fla],goci2ar['aod'][fla],'.',label='R$^2$ = {:1.3f}'.format(r))
plt.xlim(0,1.5)
plt.ylim(0,1.5)
plt.xlabel('4STAR AOD$_{{500}}$')
plt.ylabel('GOCI AOD')
plt.plot([0,1.5],[0,1.5],'--k',label='1:1')
pu.plot_lin(ar['AOD0501'][fla],goci2ar['aod'][fla],x_err=ar['UNCAOD0501'][fla],labels=True,shaded_ci=True,ci=95)

plt.legend()
plt.title('All 4STAR samples below 0.5km with nearby GOCI')


# In[175]:


len(np.unique(goci2ar['lon'][flan]))


# In[260]:


nsub = len(np.unique(goci2ar['lon'][flan]))
aod_star = np.zeros((nsub))
aod_goci = np.zeros((nsub))
ae_star = np.zeros((nsub))
ae_goci = np.zeros((nsub))
fmf_star = np.zeros((nsub))
fmf_goci = np.zeros((nsub))
aod_star_std = np.zeros((nsub))
ae_star_std = np.zeros((nsub))
fmf_star_std = np.zeros((nsub))
for j,la in enumerate(np.unique(goci2ar['lat'][flan])):
    ipixel = np.where(goci2ar['lat'][flan]==la)[0]
    aod_goci[j] = np.mean(goci2ar['aod'][flan][ipixel])
    ae_goci[j] = np.mean(goci2ar['AE'][flan][ipixel])
    fmf_goci[j] = np.mean(goci2ar['fmf'][flan][ipixel])
    aod_star[j] = np.mean(ar['AOD0501'][flan][ipixel])
    ae_star[j] = np.mean(angs[flan][ipixel])
    fmf_star[j] = np.mean(fmf['eta'][flan][ipixel])
    aod_star_std[j] = np.std(ar['AOD0501'][flan][ipixel])
    ae_star_std[j] = np.std(angs[flan][ipixel])
    fmf_star_std[j] = np.std(fmf['eta'][flan][ipixel])
    


# In[261]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
rbin = np.corrcoef(aod_star,aod_goci)[0,1]**2.0
ax[0].plot(aod_star,aod_goci,'.',label='R$^2$ = {:1.3f}'.format(rbin))
ax[0].plot([0,1.5],[0,1.5],'--k',label='1:1')
pu.plot_lin(aod_star,aod_goci,labels=True,shaded_ci=True,ci=95,ax=ax[0])
ax[0].legend()
ax[0].set_xlim(0,1.5)
ax[0].set_ylim(0,1.5)
ax[0].set_xlabel('4STAR AOD averaged within GOCI pixel')
ax[0].set_ylabel('GOCI AOD')
ax[0].set_title('KORUS-AQ Average AOD below 0.5 km')

flae = np.isfinite(ae_star) & np.isfinite(ae_goci)
rbinae = np.corrcoef(ae_star[flae],ae_goci[flae])[0,1]**2.0
ax[1].plot(ae_star,ae_goci,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
ax[1].plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(ae_star,ae_goci,labels=True,shaded_ci=True,ci=95,ax=ax[1])
ax[1].legend()
ax[1].set_ylim(0.2,1.8)
ax[1].set_xlim(0.2,1.8)
ax[1].set_xlabel('4STAR AE averaged within GOCI pixel')
ax[1].set_ylabel('GOCI AE')
ax[1].set_title('KORUS-AQ Average AE below 0.5 km')

plt.savefig(fp+'plot/KORUS_GOCI_vs_4STAR_AOD_AE.png',dpi=600,transparent=True)


# In[262]:


plt.figure()
flae = np.isfinite(ae_star) & np.isfinite(ae_goci)
rbinae = np.corrcoef(ae_star[flae],ae_goci[flae])[0,1]**2.0
plt.plot(ae_star,ae_goci,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
plt.plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(ae_star,ae_goci,labels=True,shaded_ci=True,ci=95)
plt.legend()
plt.ylim(0.2,1.8)
plt.xlim(0.2,1.8)
plt.xlabel('4STAR AE averaged within GOCI pixel')
plt.ylabel('GOCI AE')
plt.title('Average AE below 0.5 km')


# In[288]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax[0].hist(aod_star-aod_goci,bins=50,normed=True)
ax[0].axvline(0,ls='-',color='k',alpha=1,lw=0.5)
ax[0].axvline(np.nanmean(aod_star-aod_goci),ls='-',color='darkblue',alpha=0.6,lw=3,
            label='mean={:2.2f}'.format(np.nanmean(aod_star-aod_goci)))
ax[0].axvline(np.nanmedian(aod_star-aod_goci),ls='--',color='darkblue',alpha=0.6,lw=3,
            label='median={:2.2f}'.format(np.nanmedian(aod_star-aod_goci)))
ax[0].legend()
flae = np.isfinite(ae_star) & np.isfinite(ae_goci)
ax[1].hist(ae_star[flae]-ae_goci[flae],bins=50,normed=True)
ax[1].axvline(0,ls='-',color='k',alpha=1,lw=0.5)
ax[1].axvline(np.nanmean(ae_star[flae]-ae_goci[flae]),ls='-',color='darkblue',alpha=0.6,lw=3,
            label='mean={:2.2f}'.format(np.nanmean(ae_star[flae]-ae_goci[flae])))
ax[1].axvline(np.nanmedian(ae_star[flae]-ae_goci[flae]),ls='--',color='darkblue',alpha=0.6,lw=3,
            label='median={:2.2f}'.format(np.nanmedian(ae_star[flae]-ae_goci[flae])))
ax[1].legend()
ax[0].set_xlabel('AOD difference (4STAR-GOCI)')
ax[0].set_ylabel('Normalized counts')
ax[1].set_xlabel('AE difference (4STAR-GOCI)')
plt.savefig(fp+'plot/KORUS_hist_diff_GOCI_vs_4STAR_AOD_AE.png',dpi=600,transparent=True)


# ## Compare MERRA2 AOD to 4STAR

# In[264]:


plt.figure()
fla = ar['fl_QA'] & (ar['GPS_Alt']<500.0)
flan = ar['fl_QA'] & (ar['GPS_Alt']<500.0) & np.isfinite(ar['AOD0501']) & np.isfinite(merra2ar['aod'])
r = np.corrcoef(ar['AOD0501'][flan],merra2ar['aod'][flan])[0,1]**2.0
plt.plot(ar['AOD0501'][fla],merra2ar['aod'][fla],'.',label='R$^2$ = {:1.3f}'.format(r))
plt.xlim(0,1.5)
plt.ylim(0,1.5)
plt.xlabel('4STAR AOD$_{{500}}$')
plt.ylabel('MERRA2 AOD')
plt.plot([0,1.5],[0,1.5],'--k',label='1:1')
pu.plot_lin(ar['AOD0501'][fla],merra2ar['aod'][fla],x_err=ar['UNCAOD0501'][fla],labels=True,shaded_ci=True,ci=95)

plt.legend()
plt.title('All 4STAR samples below 0.5km with nearby MERRA2 AOD')


# In[273]:


merra2ar.keys()


# In[274]:


nsub = len(np.unique(merra2ar['ind'][flan]))
aod_star_m = np.zeros((nsub))
aod_merra = np.zeros((nsub))
aod_star_m_std = np.zeros((nsub))
ae_star_m = np.zeros((nsub))
ae_merra = np.zeros((nsub))
ae_star_m_std = np.zeros((nsub))
for j,la in enumerate(np.unique(merra2ar['ind'][flan])):
    ipixel = np.where(merra2ar['ind'][flan]==la)[0]
    aod_merra[j] = np.mean(merra2ar['aod'][flan][ipixel])
    aod_star_m[j] = np.mean(ar['AOD0501'][flan][ipixel])
    aod_star_m_std[j] = np.std(ar['AOD0501'][flan][ipixel])
    ae_merra[j] = np.mean(merra2ar['ae'][flan][ipixel])
    ae_star_m[j] = np.mean(angs[flan][ipixel])
    ae_star_m_std[j] = np.std(angs[flan][ipixel])


# In[267]:


plt.figure(figsize=(7,6))
flae = np.isfinite(aod_star_m) & np.isfinite(aod_merra)
rbinae = np.corrcoef(aod_star_m[flae],aod_merra[flae])[0,1]**2.0
plt.plot(aod_star_m,aod_merra,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
plt.errorbar(aod_star_m,aod_merra,xerr=aod_star_m_std,color='b',marker='.',ls='None',elinewidth=0.4,label='std dev')
plt.plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(aod_star_m,aod_merra,labels=True,shaded_ci=True,ci=95)
plt.legend()
plt.ylim(0.0,1.5)
plt.xlim(0.0,1.5)
plt.xlabel('4STAR AOD averaged within MERRA pixel')
plt.ylabel('MERRA2 AOD')
plt.title('Average AOD below 0.5 km during KORUS-AQ (May-June 2016)')
plt.savefig(fp+'plot/KORUS_MERRA2_vs_4STAR_AOD.png',dpi=600,transparent=True)


# In[277]:


plt.figure(figsize=(7,6))
flae = np.isfinite(ae_star_m) & np.isfinite(ae_merra)
rbinae = np.corrcoef(ae_star_m[flae],ae_merra[flae])[0,1]**2.0
plt.plot(ae_star_m,ae_merra,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
plt.errorbar(ae_star_m,ae_merra,xerr=ae_star_m_std,color='b',marker='.',ls='None',elinewidth=0.4,label='std dev')
plt.plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(ae_star_m,ae_merra,labels=True,shaded_ci=True,ci=95)
plt.legend()
plt.ylim(0.2,2.0)
plt.xlim(0.2,2.0)
plt.xlabel('4STAR AE averaged within MERRA pixel')
plt.ylabel('MERRA2 AE')
plt.title('Average AE below 0.5 km during KORUS-AQ (May-June 2016)')
plt.savefig(fp+'plot/KORUS_MERRA2_vs_4STAR_Ae.png',dpi=600,transparent=True)


# In[287]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
flae = np.isfinite(aod_star_m) & np.isfinite(aod_merra)
rbinae = np.corrcoef(aod_star_m[flae],aod_merra[flae])[0,1]**2.0
ax[0].plot(aod_star_m,aod_merra,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
ax[0].errorbar(aod_star_m,aod_merra,xerr=aod_star_m_std,color='b',marker='.',ls='None',elinewidth=0.4,label='std dev')
ax[0].plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(aod_star_m,aod_merra,labels=True,shaded_ci=True,ci=95,ax=ax[0])
ax[0].legend()
ax[0].set_ylim(0.0,1.5)
ax[0].set_xlim(0.0,1.5)
ax[0].set_xlabel('4STAR AOD averaged within MERRA pixel')
ax[0].set_ylabel('MERRA2 AOD')
fig.suptitle('Average aerosol below 0.5 km during KORUS-AQ (May-June 2016)')

flae = np.isfinite(ae_star_m) & np.isfinite(ae_merra)
rbinae = np.corrcoef(ae_star_m[flae],ae_merra[flae])[0,1]**2.0
ax[1].plot(ae_star_m,ae_merra,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
ax[1].errorbar(ae_star_m,ae_merra,xerr=ae_star_m_std,color='b',marker='.',ls='None',elinewidth=0.4,label='std dev')
ax[1].plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(ae_star_m,ae_merra,labels=True,shaded_ci=True,ci=95,ax=ax[1])
ax[1].legend()
ax[1].set_ylim(0.2,2.0)
ax[1].set_xlim(0.2,2.0)
ax[1].set_xlabel('4STAR AE averaged within MERRA pixel')
ax[1].set_ylabel('MERRA2 AE')

plt.savefig(fp+'plot/KORUS_MERRA2_vs_4STAR_AOD_AE_splitplot.png',dpi=600,transparent=True)


# In[281]:


plt.figure()
flae = np.isfinite(aod_star_m) & np.isfinite(aod_merra)
plt.hist(aod_star_m[flae]-aod_merra[flae],bins=50,normed=True)
plt.axvline(0,ls='-',color='k',alpha=1,lw=0.5)
plt.axvline(np.nanmean(aod_star_m[flae]-aod_merra[flae]),ls='-',color='darkblue',alpha=0.6,
            label='mean={:2.2f}'.format(np.nanmean(aod_star_m[flae]-aod_merra[flae])),lw=3)
plt.axvline(np.nanmedian(aod_star_m[flae]-aod_merra[flae]),ls='--',color='darkblue',alpha=0.6,
            label='median={:2.2f}'.format(np.nanmedian(aod_star_m[flae]-aod_merra[flae])),lw=3)
plt.legend()
plt.xlabel('AOD difference (4STAR-MERRA2)')
plt.ylabel('Normalized counts')
plt.savefig(fp+'plot/KORUS_hist_diff_MERRA2_vs_4STAR_AOD.png',dpi=600,transparent=True)


# In[282]:


plt.figure()
flae = np.isfinite(ae_star_m) & np.isfinite(ae_merra)
plt.hist(ae_star_m[flae]-ae_merra[flae],bins=50,normed=True)
plt.axvline(0,ls='-',color='k',alpha=1,lw=0.5)
plt.axvline(np.nanmean(ae_star_m[flae]-ae_merra[flae]),ls='-',color='darkblue',alpha=0.6,lw=3,
            label='mean={:2.2f}'.format(np.nanmean(ae_star_m[flae]-ae_merra[flae])))
plt.axvline(np.nanmedian(ae_star_m[flae]-ae_merra[flae]),ls='--',color='darkblue',alpha=0.6,lw=3,
            label='median={:2.2f}'.format(np.nanmedian(ae_star_m[flae]-ae_merra[flae])))
plt.legend()
plt.xlabel('AE difference (4STAR-MERRA2)')
plt.ylabel('Normalized counts')
plt.savefig(fp+'plot/KORUS_hist_diff_MERRA2_vs_4STAR_AE.png',dpi=600,transparent=True)


# In[289]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
flae = np.isfinite(aod_star_m) & np.isfinite(aod_merra)
ax[0].hist(aod_star_m[flae]-aod_merra[flae],bins=50,normed=True)
ax[0].axvline(0,ls='-',color='k',alpha=1,lw=0.5)
ax[0].axvline(np.nanmean(aod_star_m[flae]-aod_merra[flae]),ls='-',color='darkblue',alpha=0.6,
            label='mean={:2.2f}'.format(np.nanmean(aod_star_m[flae]-aod_merra[flae])),lw=3)
ax[0].axvline(np.nanmedian(aod_star_m[flae]-aod_merra[flae]),ls='--',color='darkblue',alpha=0.6,
            label='median={:2.2f}'.format(np.nanmedian(aod_star_m[flae]-aod_merra[flae])),lw=3)
ax[0].legend()
ax[0].set_xlabel('AOD difference (4STAR-MERRA2)')
ax[0].set_ylabel('Normalized counts')

flae = np.isfinite(ae_star_m) & np.isfinite(ae_merra)
ax[1].hist(ae_star_m[flae]-ae_merra[flae],bins=50,normed=True)
ax[1].axvline(0,ls='-',color='k',alpha=1,lw=0.5)
ax[1].axvline(np.nanmean(ae_star_m[flae]-ae_merra[flae]),ls='-',color='darkblue',alpha=0.6,lw=3,
            label='mean={:2.2f}'.format(np.nanmean(ae_star_m[flae]-ae_merra[flae])))
ax[1].axvline(np.nanmedian(ae_star_m[flae]-ae_merra[flae]),ls='--',color='darkblue',alpha=0.6,lw=3,
            label='median={:2.2f}'.format(np.nanmedian(ae_star_m[flae]-ae_merra[flae])))
ax[1].legend()
ax[1].set_xlabel('AE difference (4STAR-MERRA2)')
ax[1].set_ylabel('Normalized counts')
plt.savefig(fp+'plot/KORUS_hist_diff_MERRA2_vs_4STAR_AOD_AE_2plt.png',dpi=600,transparent=True)


# ## Now get the angstrom exponent and plot it vertically

# In[212]:


nwl,nm


# In[213]:


aodrr = np.array([ar[n] for n in nwl])


# In[214]:


aodrr.shape


# In[215]:


angs = su.calc_angs(ar['Start_UTC'],np.array(nm[1:11]),aodrr[1:11,:])


# In[282]:


def make_bined_alt(x,alt,days,fl,n=70,rg=None):
    'Function to create binned data for a set range, usually for altitude'
    binned_ang,binned_alt,binned_num,binned_ndays = [],[],[],[]
    if rg:
        dz = (rg[1]-rg[0])/n
    else:
        dz = np.nanmax(alt[fl])/n
        rg = [0.0,np.nanmax(alt[fl])]
    print np.nanmax(alt[fl]),dz
    for i in xrange(n):
        flaa = (alt[fl]>=(i*dz)+rg[0]) & (alt[fl]<((i+1.0)*dz)+rg[0])
        binned_ang.append(x[fl][flaa])
        binned_alt.append(np.mean([(i*dz)+rg[0],((i+1.0)*dz)+rg[0]]))
        binned_num.append(len(x[fl][flaa]))
        binned_ndays.append(len(np.unique(days[fl][flaa])))
    return binned_ang,binned_alt,binned_num,binned_ndays


# ### Plotting of the angstrom vertical dependence over Seoul

# In[217]:


ar['fl_QA_angs'] = ar['fl'] & (ar['AOD0501']>0.05) 


# In[218]:


ar['fl_QA_angs_seoul'] = ar['fl'] & (ar['AOD0501']>0.05) & (ar['Latitude']<37.75) &                        (ar['Latitude']>36.9) & (ar['Longitude']<127.30) & (ar['Longitude']>126.60)


# In[219]:


any(ar['fl_QA_angs_seoul'])


# In[220]:


bang,balt,bnum,bndays = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs'],n=90)


# In[221]:


bangs,balts,bnums,bndayss = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_seoul'],n=90)


# In[222]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(bang,positions=np.array(balt)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom from fit between 452 nm and 865 nm')
plt.ylabel('Altitude [m]')
gr = plt.cm.RdPu
bl = plt.cm.Blues
pu.set_box_whisker_color(gr,bp,bndays)
    
bpc =plt.boxplot(bangs,positions=np.array(balts)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(bl,bpc,bndayss)
bpc['boxes'][0].set_color('grey')

ax = plt.gca()
plt.title('KORUS-AQ Angstrom Exponent')
plt.ylim(0,8000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000,8000])
plt.xlim(-0.1,2.3)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['All data','Near Seoul','Mean','Median','25% - 75%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bndays)
cbaxesgr = plt.gcf().add_axes([0.83, 0.35, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.85, 0.35, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_Angstrom_fit_vertical.png',
            transparent=True,dpi=500)


# ### Plot vertical angstrom exponent for different Met. regimes

# In[226]:


ar['doys']


# In[228]:


ar['fl_QA_angs'] = ar['fl'] & (ar['AOD0501']>0.05) 
ar['fl_QA_angs_met1'] = ar['fl'] & (ar['AOD0501']>0.05) & (ar['doys']> t1[0]) & (ar['doys']< t1[1])
ar['fl_QA_angs_met2'] = ar['fl'] & (ar['AOD0501']>0.05)  & (ar['doys']> t2[0]) & (ar['doys']< t2[1])
ar['fl_QA_angs_met3'] = ar['fl'] & (ar['AOD0501']>0.05)  & (ar['doys']> t3[0]) & (ar['doys']< t3[1])
ar['fl_QA_angs_met4'] = ar['fl'] & (ar['AOD0501']>0.05)  & (ar['doys']> t4[0]) & (ar['doys']< t4[1])

#ar['fl_QA_angs_seoul'] = ar['fl'] & (ar['AOD0501']>0.05) & (ar['Latitude']<37.75) &\
#                        (ar['Latitude']>36.9) & (ar['Longitude']<127.30) & (ar['Longitude']>126.60)


# In[287]:


bang,balt,bnum,bndays = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs'],n=30,rg=[0.0,8000.0])
bangm1,baltm1,bnumm1,bndaysm1 = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_met1'],n=30,rg=[0.0,8000.0])
bangm2,baltm2,bnumm2,bndaysm2 = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_met2'],n=30,rg=[0.0,8000.0])
bangm3,baltm3,bnumm3,bndaysm3 = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_met3'],n=30,rg=[0.0,8000.0])
bangm4,baltm4,bnumm4,bndaysm4 = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_met4'],n=30,rg=[0.0,8000.0])


# In[291]:


plt.figure(figsize=(5.5,7.5))
bp =plt.boxplot(bang,positions=np.array(balt)-25.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom from fit between 452 nm and 865 nm')
plt.ylabel('Altitude [m]')

rd = plt.cm.RdPu
bl = plt.cm.Blues
og = plt.cm.YlOrBr
gr = plt.cm.Greens
k = plt.cm.Greys

pu.set_box_whisker_color(k,bp,bndays,mean_color='grey',median_color='grey')
    
bp1 = plt.boxplot(bangm1,positions=np.array(baltm1)+00.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(rd,bp1,bndaysm1,mean_color='tab:red',median_color='r')
bp1['boxes'][0].set_color('grey')
bp2 = plt.boxplot(bangm2,positions=np.array(baltm2)+25.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(bl,bp2,bndaysm2,mean_color='tab:blue',median_color='b')
bp3 = plt.boxplot(bangm3,positions=np.array(baltm3)+50.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(og,bp3,bndaysm3,mean_color='tab:orange',median_color='y')
bp4 = plt.boxplot(bangm4,positions=np.array(baltm4)+75.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(gr,bp4,bndaysm4,mean_color='tab:green',median_color='g')

ax = plt.gca()
plt.title('KORUS-AQ Angstrom Exponent')
plt.ylim(0,8000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000,8000])
plt.xlim(-0.1,2.5)
plt.grid()
plt.legend([bp['boxes'][5],bp1['boxes'][15],bp2['boxes'][16],bp3['boxes'][7],bp4['boxes'][-4],
            bp['means'][0],bp['medians'][0],bp['boxes'][-2],bp['whiskers'][0]],
           ['All data','Dynamic','Stagnation','Extreme\npollution','Blocking','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bndays)
scalarmapk = plt.cm.ScalarMappable(cmap=k)
scalarmapk.set_array(bndays)
scalarmaprd = plt.cm.ScalarMappable(cmap=rd)
scalarmaprd.set_array(bndays)
scalarmapog = plt.cm.ScalarMappable(cmap=og)
scalarmapog.set_array(bndays)
cbaxesgr = plt.gcf().add_axes([0.83, 0.15, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.87, 0.15, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbaxesk = plt.gcf().add_axes([0.91, 0.15, 0.015, 0.3])
cbk = plt.colorbar(scalarmapk,cax=cbaxesk)
cbaxesog = plt.gcf().add_axes([0.85, 0.15, 0.015, 0.3])
cbo = plt.colorbar(scalarmapog,cax=cbaxesog)
cbaxesrd = plt.gcf().add_axes([0.89, 0.15, 0.015, 0.3])
cbr = plt.colorbar(scalarmaprd,cax=cbaxesrd)
cbg.set_ticks([0,3,6,9,12,15,18])
cbb.set_ticks([0,3,6,9,12,15,18]),cbb.set_ticklabels(['','','','','',''])
cbk.set_ticks([0,3,6,9,12,15,18]),cbk.set_ticklabels(['','','','','',''])
cbo.set_ticks([0,3,6,9,12,15,18]),cbo.set_ticklabels(['','','','','',''])
cbr.set_ticks([0,3,6,9,12,15,18]),cbr.set_ticklabels(['','','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_Angstrom_fit_vertical_met.png',
            transparent=True,dpi=500)


# ## Analyse the Fine mode fraction

# In[50]:


ar['fl_QA_low'] = ar['fl_QA'] & (ar['GPS_Alt']<500.0)
ar['fl_QA_mid'] = ar['fl_QA'] & (ar['GPS_Alt']>2000.0) & (ar['GPS_Alt']<5000.0) 


# In[51]:


ar['fl_QA_fmf'] = ar['fl_QA'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc']))


# In[52]:


bfaod,baltf,bnumf,bndaysf = make_bined_alt(fmf['tauf'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)
bcaod,baltc,bnumc,bndaysc = make_bined_alt(fmf['tauc'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)
beta,balte,bnume,bndayse = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)


# In[53]:


blat,baltl,bnuml,bndaysl = make_bined_alt(ar['Latitude'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)
blon,baltlo,bnumlo,bndayslo = make_bined_alt(ar['Longitude'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)


# In[54]:


blats = [np.nanmedian(ll) for ll in blat]
blons = [np.nanmedian(ll) for ll in blon]


# In[90]:


blons


# ### Plot the fine mode fraction distribution

# In[87]:


plt.figure()
plt.plot(aodrr[2,:],label='measurement AOD 500 nm')
plt.plot(fmf['tau'],'.k',label='fit AOD 500 nm')
plt.plot(fmf['tauf'], 'ob',label='Fine mode fraction AOD')
plt.plot(fmf['tauc'],'sr',label='Coarse mode fraction AOD')
plt.legend()
plt.ylim(0,1.5)


# In[96]:


plt.figure()
plt.hist(fmf['tauc'][ar['fl_QA']]+fmf['tauf'][ar['fl_QA']],range=[0,1.5],bins=50,label='total')
plt.hist(fmf['tauc'][ar['fl_QA']],range=[0,1.5],bins=50,label='Coarse mode')
plt.legend(frameon=False)


# In[106]:


any(ar['fl_QA_mid'])


# ### Plot the histogram distribution of the fine mode fraction

# In[61]:


plt.figure(figsize=(6,7))
ax1 = plt.subplot(3,1,1)
plt.hist([fmf['tauc'][ar['fl_QA']],fmf['tauf'][ar['fl_QA']]],color=['r','b'],histtype='bar',
            bins=50,range=[0.0,1.5],label=['Coarse','Fine'],edgecolor='None',alpha=0.75,normed=False,stacked=True)
plt.legend(frameon=True,loc=1)
plt.title('KORUS-AQ AOD fine/coarse mode - all data')
#plt.xlabel('AOD 500 nm')
plt.ylabel('Counts')
plt.yscale('log'),plt.xscale('log')
plt.ylim(5,250000),plt.xlim(0.01,1.5)
plt.xticks([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.5])
ax1.set_xticklabels([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,'',1.0,'',1.5])

ax2 = plt.subplot(3,1,2,sharex=ax1)
plt.hist([fmf['tauc'][ar['fl_QA_mid']],fmf['tauf'][ar['fl_QA_mid']]],color=['r','b'],histtype='bar',
            bins=50,range=[0.0,1.5],label=['Coarse','Fine'],edgecolor='None',alpha=0.75,normed=False,stacked=True)
#plt.legend(frameon=False)
plt.title('Between 2 and 5 km')
plt.ylabel('Counts')
plt.yscale('log'),plt.xscale('log')
plt.ylim(5,250000),plt.xlim(0.01,1.5)
plt.xticks([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.5])
ax2.set_xticklabels([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,'',1.0,'',1.5])

ax3 = plt.subplot(3,1,3,sharex=ax2)
plt.hist([fmf['tauc'][ar['fl_QA_low']],fmf['tauf'][ar['fl_QA_low']]],color=['r','b'],histtype='bar',
            bins=50,range=[0.0,1.5],label=['Coarse','Fine'],edgecolor='None',alpha=0.75,normed=False,stacked=True)
#plt.legend(frameon=False)
plt.title('Below 0.5 km')
#plt.xlabel('AOD 500 nm')
plt.ylabel('Counts')
plt.yscale('log'),plt.xscale('log')
plt.ylim(5,250000),plt.xlim(0.01,1.5)
plt.xlabel('AOD 500 nm')
plt.xticks([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.5])
ax3.set_xticklabels([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,'',1.0,'',1.5])

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_fine_mode_hist.png',
            transparent=True,dpi=500)


# ### Plot the vertical dependence of the fine mode fraction

# In[110]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(bfaod,positions=np.array(baltf)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('AOD 500 nm')
plt.ylabel('Altitude [m]')
bl = plt.cm.YlOrRd
gr = plt.cm.Blues
pu.set_box_whisker_color(gr,bp,bndaysf)
    
bpc =plt.boxplot(bcaod,positions=np.array(baltc)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(bl,bpc,bndaysc)
bpc['boxes'][-1].set_color('grey')

ax = plt.gca()
plt.title('KORUS-AQ Fine/Coarse mode AOD')
plt.ylim(0,8000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000,8000])
plt.xlim(0.0,0.65)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][-1],bpc['whiskers'][0]],
           ['Fine mode','Coarse mode','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bndaysf)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bndaysc)
cbaxesgr = plt.gcf().add_axes([0.83, 0.35, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.85, 0.35, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,6,12,16,18,20])
cbb.set_ticks([0,6,12,16,18,20]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_fine_mode_AOD_vertical.png',
            transparent=True,dpi=500)


# In[56]:


blats[0]=36.2


# In[57]:


bndm = np.nanmax(blats)*1.0
bndm


# In[100]:


cl = gr
for j,q in enumerate(blats):
    print j, q, cl(blats[j]*1.0/bndm)


# In[111]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(beta,positions=np.array(balte),vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('fine mode fraction')
plt.ylabel('Altitude [m]')
bl = plt.cm.YlOrRd
gr = plt.cm.Blues
pu.set_box_whisker_color(gr,bp,blats,color_not_start_at_zero=True)
    
#bpc =plt.boxplot(bcaod,positions=np.array(baltc)+10.0,vert=False,
#                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
#pu.set_box_whisker_color(bl,bpc,bndaysc)
bp['boxes'][-1].set_color('grey')

ax = plt.gca()
plt.title('KORUS-AQ fine mode fraction')
plt.ylim(0,8000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000,8000])
plt.xlim(0.0,1.0)
plt.grid()
plt.legend([bp['boxes'][5],bp['means'][5],bp['medians'][5],bp['boxes'][-1],bp['whiskers'][5]],
           ['AOD fmf','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(blats)
#scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
#scalarmapbl.set_array(bndays)
cbaxesgr = plt.gcf().add_axes([0.88, 0.35, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
#cbaxesbl = plt.gcf().add_axes([0.85, 0.35, 0.015, 0.3])
#cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
#cbg.set_ticks([0,6,12,15,18])
#cbb.set_ticks([0,6,12,15,18]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left')#,cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-9.0,0.5,'Median Latitude',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_fine_mode_AOD_vertical.png',
            transparent=True,dpi=500)


# ### Split fmf for met

# In[293]:


ar['fl_QA_fmf'] = ar['fl_QA'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc']))
ar['fl_QA_fmf_met1'] = ar['fl'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc'])) & (ar['doys']> t1[0]) & (ar['doys']< t1[1])
ar['fl_QA_fmf_met2'] = ar['fl'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc'])) & (ar['doys']> t2[0]) & (ar['doys']< t2[1])
ar['fl_QA_fmf_met3'] = ar['fl'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc'])) & (ar['doys']> t3[0]) & (ar['doys']< t3[1])
ar['fl_QA_fmf_met4'] = ar['fl'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc'])) & (ar['doys']> t4[0]) & (ar['doys']< t4[1])


# In[ ]:


beta,balte,bnume,bndayse = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)
beta1,balte1,bnume1,bndayse1 = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf_met1'],n=90)
beta2,balte2,bnume2,bndayse2 = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf_met2'],n=90)
beta3,balte3,bnume3,bndayse3 = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf_met3'],n=90)
beta4,balte4,bnume4,bndayse4 = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf_met4'],n=90)


# ## Calculate the autocorrelation of the fine and coarse mode AOD

# In[55]:


fvals = {'utc':ar['Start_UTC'][fl],'alt':ar['GPS_Alt'][fl],'lat':ar['Latitude'][fl],'lon':ar['Longitude'][fl],
        'aod0500':ar['AOD0501'][fl],'aod1040':ar['AOD1040'][fl],'aodf':fmf['tauf'][fl],'aodc':fmf['tauc'][fl],'eta':fmf['eta'][fl]}


# In[56]:


dfvals = get_segments(f_level,fvals,nsep=100)


# In[57]:


ddfv = get_distances(dfvals)


# Now the segments are identified and the cumulative distances are quantified, we must interpolate over the segments, to remove any missing data.

# In[58]:


def interp_dist_fmf(d,dist=0.12):
    'function to insterpolate the AOD from the dict to an even grid spacing accroding to distance (default 0.12 km)'
    d['cdist_n'],d['aod_nf'],d['aod_nc'],d['eta_n'] = [],[],[],[]
    for i,cd in enumerate(d['cumdist']):        
        d['cdist_n'].append(np.arange(cd.min(),cd.max(),dist))
        if np.sum(np.isfinite(d['aodf'][i]))/float(len(d['aodf'][i])) < 0.75: # check if at least 75% of the segment is valid
            af = np.array(np.nan)
            ac = np.array(np.nan)
            et = np.array(np.nan)
        else:
            try:
                fcdf = interpolate.interp1d(cd,d['aodf'][i])
                af = fcdf(d['cdist_n'][i])
                fcdc = interpolate.interp1d(cd,d['aodc'][i])
                ac = fcdc(d['cdist_n'][i])
                fcde = interpolate.interp1d(cd,d['eta'][i])
                et = fcde(d['cdist_n'][i])
            except TypeError:
                af = np.array(np.nan)
                ac = np.array(np.nan)
                et = np.array(np.nan)
        d['aod_nf'].append(af)
        d['aod_nc'].append(ac)
        d['eta_n'].append(et)


# In[59]:


interp_dist_fmf(dfvals)


# In[79]:


dfvals['autocor_f'],dfvals['autocor_c'],dfvals['autocor_e'] = [] ,[],[]
for i,a in enumerate(dfvals['aod_nf']):
    #auf,auc,eut = [],[],[]
    try:
        auf = autocorr5(a)
        auc = autocorr5(dfvals['aod_nc'][i])
        eut = autocorr5(dfvals['eta_n'][i])
    except:
        auf = np.array([np.nan])
        auc = np.array([np.nan])
        eut = np.array([np.nan])
    dfvals['autocor_f'].append(auf[:])
    dfvals['autocor_c'].append(auc[:])
    dfvals['autocor_e'].append(eut[:])


# In[61]:


len(dfvals['aod_nf'])


# In[63]:


dfvals['aod_nf'][i]


# In[64]:


dfvals['cdist_n'][i]


# In[80]:


dfvals['autocor_f'][i]


# In[111]:


mc = np.max([len(m) for m in dfvals['autocor_c']])
imc = np.argmax([len(m) for m in dfvals['autocor_c']])


# In[113]:


cdist = dfvals['cdist_n'][imc]


# In[123]:


autocor_c = np.zeros((len(dfvals['autocor_c']),mc))+np.nan
autocor_f = np.zeros((len(dfvals['autocor_f']),mc))+np.nan
autocor_e = np.zeros((len(dfvals['autocor_e']),mc))+np.nan


# In[124]:


for i,c in enumerate(dfvals['autocor_c']): autocor_c[i,:len(c)]=c
for i,c in enumerate(dfvals['autocor_f']): autocor_f[i,:len(c)]=c
for i,c in enumerate(dfvals['autocor_e']): autocor_e[i,:len(c)]=c


# ### Plot out the autocorrelation 

# In[83]:


plt.figure()
for i,j in enumerate(dfvals['cdist_n']):
    try:
        if len(j)<1: continue
        plt.plot(j,dfvals['autocor_f'][i],'.')
    except:
        continue
plt.ylabel('Correlation Coefficient for fine mode')
plt.xlabel('Distance [km]')
plt.xscale('log')
plt.xlim(0.1,500)


# In[84]:


plt.figure()
for i,j in enumerate(dfvals['cdist_n']):
    try:
        plt.plot(j,dfvals['autocor_c'][i],'.')
    except:
        pass
plt.ylabel('Correlation Coefficient for coarse mode')
plt.xlabel('Distance [km]')
plt.xscale('log')


# In[110]:


dfvals['cdist_n'][1:3]


# In[133]:


autocor_c_ma = np.ma.masked_array(autocor_c,mask=np.isnan(autocor_c))


# In[157]:


def make_binned(x,alt,fl,bn,flb):
    'Function to create binned data for a set range, usually for altitude'
    import numpy as np
    binned_ang,binned_alt,binned_num = [],[],[]
    for i,b in enumerate(bn[:-1]):
        flaa = (alt[flb]>=b) & (alt[flb]<bn[i+1])
        binned_ang.append(x[:,flb][flaa])
        binned_alt.append(np.mean([b,bn[i+1]]))
        binned_num.append(len(x[fl][:,flaa]))
    return binned_ang,binned_alt,binned_num,binned_ndays


# In[154]:


bnc = np.logspace(0.1,3.0)


# In[155]:


bnc


# In[163]:


auc = make_binned(autocor_c,cdist,np.isfinite(autocor_c),bnc)


# In[161]:


b = bnc[0]
i = 0
fl = np.isfinite(autocor_c)


# In[162]:


flaa = (cdist[fl]>=b) & (cdist[fl]<bnc[i+1])


# In[156]:


autocor_c.shape


# In[135]:


plt.figure()
bp = plt.boxplot(autocor_c_ma[:,0:20],positions=cdist[0:20],vert=True,
            showfliers=False,widths=90,showmeans=True,patch_artist=True) 


# In[118]:


autocor_c.shape


# In[134]:


autocor_c_ma[:,0:20]


# In[130]:


np.isfinite(autocor_c)


# In[136]:


autocor_c_ma[:,0]


# ## Map out the level legs and numbers

# In[181]:



#set up a easy plotting function
def make_map(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=128,lat_0=36.0,
            llcrnrlon=123.0, llcrnrlat=32.0,
            urcrnrlon=132.0, urcrnrlat=39,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(123,133,11),labels=[0,0,0,1])
    m.drawparallels(np.linspace(31,39,17),labels=[1,0,0,0])
    return m


# In[303]:


fig,ax = plt.subplots(1,1)
m = make_map(ax)
m.plot(ar['Longitude'],ar['Latitude'],'.',markersize=0.2,color='tab:blue',latlon=True,label='All data')
m.plot(ar['Longitude'][fl][f_level],ar['Latitude'][fl][f_level],'.',markersize=0.5,color='tab:red',latlon=True,label='level legs')
#m.plot(s['Lon'][it],s['Lat'][it],'r+',latlon=True)
plt.legend(markerscale=20)
plt.savefig(fp+'plot/KORUS_map_QA.png',dpi=600,transparent=True)


# ## Plot the AOD histogram by met

# In[497]:


ar['fl_QA_met1'] = ar['fl'] & (ar['doys']> t1[0]) & (ar['doys']< t1[1])
ar['fl_QA_met2'] = ar['fl'] & (ar['doys']> t2[0]) & (ar['doys']< t2[1])
ar['fl_QA_met3'] = ar['fl'] & (ar['doys']> t3[0]) & (ar['doys']< t3[1])
ar['fl_QA_met4'] = ar['fl'] & (ar['doys']> t4[0]) & (ar['doys']< t4[1])


# In[337]:


fig = plt.figure()
n=plt.hist([ar['AOD0501'][ar['fl']],
            ar['AOD0501'][ar['fl_QA_met1']],
          ar['AOD0501'][ar['fl_QA_met2']],
          ar['AOD0501'][ar['fl_QA_met3']],ar['AOD0501'][ar['fl_QA_met4']]],
           bins=20,range=(0,1.4),normed=True,edgecolor='None',alpha=0.7,
         label=['All data','Dynamic','Stagnation','Extreme\npollution','Blocking'])
#y = [(nn+n[1][j+1])/2.0 for j,nn in enumerate(n[1][:-1])]
#for i,p in enumerate(n[-1]):
#    plt.plot(y,n[0][i],'-',color=p[0].get_facecolor(),lw=3)
plt.legend(frameon=False)
plt.grid()
plt.xlim(0,0.5)
plt.xlabel('AOD @ 501 nm')
plt.ylabel('Normalized counts')
#plt.title('Histogram of 4STAR AOD from KORUS-AQ subsetted by altitude')

left, bottom, width, height = [0.63, 0.3, 0.35, 0.2]
ax2 = fig.add_axes([left, bottom, width, height])
n = ax2.hist([ar['AOD0501'][ar['fl']],
            ar['AOD0501'][ar['fl_QA_met1']],
          ar['AOD0501'][ar['fl_QA_met2']],
          ar['AOD0501'][ar['fl_QA_met3']],ar['AOD0501'][ar['fl_QA_met4']]],
             bins=20,range=(0,1.4),normed=True,edgecolor='None',alpha=0.7)
ax2.set_xlim(0.4,1.4)
ax2.set_ylim(0,0.6)
ax2.grid()
ax2.set_xlabel('AOD @ 501 nm')

plt.savefig(fp+'plot/AOD_hist_met_KORUS.png',dpi=600,transparent=True)


# ## Plot the AOD spectra by met

# In[493]:


aod_names = sorted([a for a in ar.keys() if ('AOD' in a) and not ('UNC' in a)])


# In[494]:


aod_names


# In[495]:


wvl = np.array([380,452,501,520,532,550,606,620,675,781,865,1020,1040,1064,1236,1559,1627])
wvl_bins = np.append(wvl[0]-10,wvl+10)


# In[498]:


ars = []
wvs = []
fls = {'fl':[],'fl_QA_met1':[],'fl_QA_met2':[],'fl_QA_met3':[],'fl_QA_met4':[]}
for i,a in enumerate(aod_names):
    ars.append(ar[a])
    wvs.append(ar[a]*0.0+wvl[i])
    for ff in fls.keys():
        fls[ff].append(ar[ff])
ars = np.array(ars)
wvs = np.array(wvs)
for ff in fls.keys():
    fls[ff] = np.array(fls[ff])
    fls[ff] = fls[ff].reshape(fls[ff].size)
arsn = ars.reshape(ars.size)
wvsn = wvs.reshape(wvs.size)


# In[508]:


np.nanmean(ars[2,:])


# In[505]:


wvs[2,0]


# In[502]:


plt.figure()
plt.plot(wvsn[fls['fl']],arsn[fls['fl']],'.',alpha=0)
plt.yscale('log')
plt.xscale('log')
plt.xticks([350,400,500,600,700,800,900,1000,1200,1400,1700])
plt.xlim(350,1700)
plt.ylim([0.01,2.0])

pu.make_boxplot(arsn[fls['fl']],wvsn[fls['fl']],wvl_bins,wvl,y=1,alpha=0.5,label='All points',
                fliers_off=True,color='k')
pu.make_boxplot(arsn[fls['fl_QA_met1']],wvsn[fls['fl_QA_met1']],wvl_bins,wvl+2,y=1,alpha=0.5,label='Dynamic',
                fliers_off=True,color='tab:orange')
pu.make_boxplot(arsn[fls['fl_QA_met2']],wvsn[fls['fl_QA_met2']],wvl_bins,wvl+4,y=1,alpha=0.5,label='Stagnation',
                fliers_off=True,color='tab:green')
pu.make_boxplot(arsn[fls['fl_QA_met3']],wvsn[fls['fl_QA_met3']],wvl_bins,wvl+6,y=1,alpha=0.5,label='Extreme\npollution',
                fliers_off=True,color='tab:red')
pu.make_boxplot(arsn[fls['fl_QA_met4']],wvsn[fls['fl_QA_met4']],wvl_bins,wvl+8,y=1,alpha=0.5,label='Blocking',
                fliers_off=True,color='tab:purple')

plt.legend(frameon=False)

plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.grid()
plt.title('Average AOD spectra')

#plt.savefig(fp+'plot/KORUS_AOD_wvl_loglog_bymet.png',dpi=600,transparent=True)


# # Spatial maps and time traces

# In[509]:


def stats_2d(lat,lon,x,fl=[],bins=26,rg=[[-25,-8],[0,16]],days=[],verbose=True):
    "Combined Statistics function to get the mean, median, number, ndays, and std from a 2d dataset"
    import scipy.stats as st
    import numpy as np
    
    stat = {}
    if not len(fl)>0: fl = np.isfinite(x)
        
    stat['mean'],stat['xm'],stat['ym'],stat['bin'] =           st.binned_statistic_2d(lat[fl],lon[fl],x[fl],bins=bins,range=rg,statistic=np.nanmean)
    stat['mean'] = np.ma.masked_array(stat['mean'],np.isnan(stat['mean']))
    
    stat['median'],stat['xe'],stat['ye'],stat['bine'] =           st.binned_statistic_2d(lat[fl],lon[fl],x[fl],bins=bins,range=rg,statistic=np.nanmedian)
    stat['median'] = np.ma.masked_array(stat['median'],np.isnan(stat['median']))

    stat['std'],stat['xs'],stat['ys'],stat['bins'] =           st.binned_statistic_2d(lat[fl],lon[fl],x[fl],bins=bins,range=rg,statistic=np.nanstd)
    stat['std'] = np.ma.masked_array(stat['std'],np.isnan(stat['std']))
    
    stat['cnt'],stat['xn'],stat['yn'],stat['binn'] =           st.binned_statistic_2d(lat[fl],lon[fl],x[fl],bins=bins,range=rg,statistic='count')
    stat['cnt'] = np.ma.masked_array(stat['cnt'],np.isnan(stat['cnt']))

    if len(days)>0:
        uniq_cnt = lambda x: len(np.unique(x))
        stat['dcnt'],stat['xd'],stat['yd'],stat['bind'] =           st.binned_statistic_2d(lat[fl],lon[fl],days[fl],bins=bins,range=rg,statistic=uniq_cnt)
        stat['dcnt'] = np.ma.masked_array(stat['dcnt'],np.isnan(stat['dcnt']))
    else:
        stat['dcnt'] = stat['cnt']*0.0
    
    if verbose:
        print 'Mean values: mean={}, median={}, std={}, num={}, day={}'.format(                    np.nanmean(stat['mean']),np.nanmean(stat['median']),np.nanmean(stat['std']),
                    np.nanmean(stat['cnt']),np.nanmean(stat['dcnt']))
        print 'Median values: mean={}, median={}, std={}, num={}, day={}'.format(                    np.nanmedian(stat['mean']),np.nanmedian(stat['median']),np.nanmedian(stat['std']),
                    np.nanmedian(stat['cnt']),np.nanmedian(stat['dcnt']))
        print 'STD values: mean={}, median={}, std={}, num={}, day={}'.format(                    np.nanstd(stat['mean']),np.nanstd(stat['median']),np.nanstd(stat['std']),
                    np.nanstd(stat['cnt']),np.nanstd(stat['dcnt']))
    return stat


# In[510]:


rg = [[32.5,38.5],[123.5,131.5]]
nbins = 18


# ## Make a map of 4STAR AOD

# In[511]:


flalt = ar['fl'] & (ar['GPS_Alt']<1000.0)
astat_aod = stats_2d(ar['Latitude'],ar['Longitude'],ar['AOD0501'],fl=flalt,days=ar['days'],
                     bins=nbins,rg=rg,verbose=True)


# In[512]:


astat_aod.keys()


# In[513]:


iao = np.where((astat_aod['cnt'].data>0.0) & (astat_aod['std'].data<1.0))


# In[517]:


astat_aod['xm'][15]


# In[520]:


astat_aod['ym'][8]


# In[522]:


astat_aod['mean'][15,8] # for Seoul


# In[524]:


astat_aod['mean'][14,8] # for South of Seoul


# In[526]:


astat_aod['mean'][13,7] # for South of Incheon/Suwon


# In[556]:


#fig = plt.figure(figsize=(11,4.5))
fig,ax = plt.subplots(1,2,figsize=(11,4.5))
ax1,ax2 = ax[0],ax[1]

m = make_map(ax=ax1)
m.shadedrelief(alpha=0.4)

ym,xm = m(astat_aod['ym'],astat_aod['xm'])
p = ax1.pcolor(ym,xm,astat_aod['mean'],vmin=0.0,vmax=0.8,cmap='magma')

ax1.set_title('4STAR - Mean AOD$_{{501}}$')
cb = plt.colorbar(p,extend='both',ax=ax1,pad=0.02)

m2 = make_map(ax=ax2)
m2.shadedrelief(alpha=0.4)

y2,x2 = m2(astat_aod['ys'][iao[1]],astat_aod['xs'][iao[0]])
p2 = ax2.scatter(y2,x2,8.0+(astat_aod['cnt'].data[iao[0],iao[1]].flatten()/25.0),
                c=astat_aod['std'].data[iao[0],iao[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0,vmax=0.35)
ax2.set_title('4STAR - Standard Deviation AOD$_{{501}}$')
cb2 = plt.colorbar(p2,extend='max',ax=ax2,pad=0.02)

sizes = [10,100,500,1500]
labels = ['N={0}'.format(z) for z in sizes]
points = [ax2.scatter([], [], s=8.0+(z/25.0), c='grey',marker='s',edgecolor='None') for z in sizes]
ax2.legend(points, labels, scatterpoints=1,frameon=True,
           framealpha=0.8,handletextpad=0.1,labelspacing=0.1,borderpad=0.1,loc='lower right')
plt.tight_layout(pad=1.12,h_pad=1.8,w_pad=3.0,rect=(0.05,0,1,1))
plt.savefig(fp+'plot/KORUS_4STAR_AOD_map_{}.png'.format(vv),dpi=600,transparent=True)


# ## Expand map to include MERRA-2

# In[559]:


merra2ar.keys()


# In[560]:


mstat_aod = stats_2d(merra2ar['lat'],merra2ar['lon'],merra2ar['aod'],fl=flalt,days=ar['days'],
                     bins=nbins,rg=rg,verbose=True)


# In[561]:


iao = np.where((mstat_aod['cnt'].data>0.0) & (mstat_aod['std'].data<1.0))


# In[562]:


#fig = plt.figure(figsize=(11,4.5))
fig,ax = plt.subplots(1,2,figsize=(11,4.5))
ax1,ax2 = ax[0],ax[1]

m = make_map(ax=ax1)
m.shadedrelief(alpha=0.4)

ym,xm = m(mstat_aod['ym'],mstat_aod['xm'])
p = ax1.pcolor(ym,xm,mstat_aod['mean'],vmin=0.0,vmax=0.8,cmap='magma')

ax1.set_title('MERRA2 - Mean AOD$_{{501}}$')
cb = plt.colorbar(p,extend='both',ax=ax1,pad=0.02)

m2 = make_map(ax=ax2)
m2.shadedrelief(alpha=0.4)

y2,x2 = m2(mstat_aod['ys'][iao[1]],mstat_aod['xs'][iao[0]])
p2 = ax2.scatter(y2,x2,8.0+(mstat_aod['cnt'].data[iao[0],iao[1]].flatten()/25.0),
                c=mstat_aod['std'].data[iao[0],iao[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0,vmax=0.35)
ax2.set_title('MERRA2 - Standard Deviation AOD$_{{501}}$')
cb2 = plt.colorbar(p2,extend='max',ax=ax2,pad=0.02)

sizes = [10,100,500,1500]
labels = ['N={0}'.format(z) for z in sizes]
points = [ax2.scatter([], [], s=8.0+(z/25.0), c='grey',marker='s',edgecolor='None') for z in sizes]
ax2.legend(points, labels, scatterpoints=1,frameon=True,
           framealpha=0.8,handletextpad=0.1,labelspacing=0.1,borderpad=0.1,loc='lower right')
plt.tight_layout(pad=1.12,h_pad=1.8,w_pad=3.0,rect=(0.05,0,1,1))
plt.savefig(fp+'plot/KORUS_MERRA2_AOD_map_{}.png'.format(vv),dpi=600,transparent=True)


# ## Expand map to include GOCI

# In[563]:


goci2ar.keys()


# In[564]:


gstat_aod = stats_2d(goci2ar['lat'],goci2ar['lon'],goci2ar['aod'],fl=flalt,days=ar['days'],
                     bins=nbins,rg=rg,verbose=True)


# In[565]:


iao = np.where((gstat_aod['cnt'].data>0.0) & (gstat_aod['std'].data<1.0))


# In[566]:


#fig = plt.figure(figsize=(11,4.5))
fig,ax = plt.subplots(1,2,figsize=(11,4.5))
ax1,ax2 = ax[0],ax[1]

m = make_map(ax=ax1)
m.shadedrelief(alpha=0.4)

ym,xm = m(gstat_aod['ym'],gstat_aod['xm'])
p = ax1.pcolor(ym,xm,gstat_aod['mean'],vmin=0.0,vmax=0.8,cmap='magma')

ax1.set_title('GOCI YAER - Mean AOD$_{{501}}$')
cb = plt.colorbar(p,extend='both',ax=ax1,pad=0.02)

m2 = make_map(ax=ax2)
m2.shadedrelief(alpha=0.4)

y2,x2 = m2(gstat_aod['ys'][iao[1]],gstat_aod['xs'][iao[0]])
p2 = ax2.scatter(y2,x2,8.0+(gstat_aod['cnt'].data[iao[0],iao[1]].flatten()/25.0),
                c=gstat_aod['std'].data[iao[0],iao[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0,vmax=0.35)
ax2.set_title('GOCI YAER - Standard Deviation AOD$_{{501}}$')
cb2 = plt.colorbar(p2,extend='max',ax=ax2,pad=0.02)

sizes = [10,100,500,1500]
labels = ['N={0}'.format(z) for z in sizes]
points = [ax2.scatter([], [], s=8.0+(z/25.0), c='grey',marker='s',edgecolor='None') for z in sizes]
ax2.legend(points, labels, scatterpoints=1,frameon=True,
           framealpha=0.8,handletextpad=0.1,labelspacing=0.1,borderpad=0.1,loc='lower right')
plt.tight_layout(pad=1.12,h_pad=1.8,w_pad=3.0,rect=(0.05,0,1,1))
plt.savefig(fp+'plot/KORUS_GOCI_AOD_map_{}.png'.format(vv),dpi=600,transparent=True)


# In[568]:


gstat_aod.keys()


# In[573]:


np.diff(gstat_aod['xs']),np.diff(gstat_aod['ys'])


# In[ ]:




