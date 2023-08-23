#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     MODACAERO project. To look at the results from the improved above cloud aod retrievals from MODIS during ORACLES
#     Compares to 4STAR and HSRL over the same region as MODIS
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
#   - MODACAERO updated retrieval files
#   - 4STAR ORACLES netcdf
#   - HSRL AOD retrievals for ORACLES
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2023-08-11
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
from tqdm.notebook import tqdm 
import warnings
import map_utils as mu
import plotting_utils as pu


# In[2]:


name = 'ACAERO'
vv = 'v1'
fp = getpath(name)


# # Load files

# ## Load 4STAR ORACLES

# In[3]:


fps = getpath('sunsat')


# In[4]:


fpo = getpath('ORACLES')


# ### 2016

# In[5]:


ar6 = hs.loadmat(fpo+'/aod_ict/R4/all_aod_ict_R4_2016.mat')


# In[6]:


ar6['flac'] = (ar6['qual_flag']==0)&(ar6['flag_acaod']==1)
ar6['flacr'] = (ar6['qual_flag']==0)&(ar6['flag_acaod']==1)&(ar6['fl_routine'])
ar6['flaco'] = (ar6['qual_flag']==0)&(ar6['flag_acaod']==1)&~(ar6['fl_routine'])
ar6['flr'] = (ar6['qual_flag']==0) & (ar6['fl_routine'])
ar6['flo'] = (ar6['qual_flag']==0) & ~(ar6['fl_routine'])
ar6['fl'] = (ar6['qual_flag']==0)


# ### 2017

# In[7]:


ar7 = hs.loadmat(fpo+'/aod_ict_2017/R1/all_aod_ict_R1_2017.mat')


# In[8]:


ar7['flac'] = (ar7['qual_flag']==0)&(ar7['flag_acaod']==1)
ar7['flacr'] = (ar7['qual_flag']==0)&(ar7['flag_acaod']==1)&(ar7['fl_routine'])
ar7['flaco'] = (ar7['qual_flag']==0)&(ar7['flag_acaod']==1)&~(ar7['fl_routine'])
ar7['flr'] = (ar7['qual_flag']==0) & (ar7['fl_routine'])
ar7['flo'] = (ar7['qual_flag']==0) & ~(ar7['fl_routine'])
ar7['fl'] = (ar7['qual_flag']==0)


# ### 2018

# In[9]:


ar8 = hs.loadmat(fpo+'/aod_ict_2018/{vv}/all_aod_ict_{vv}_2018.mat'.format(vv='R1'))


# In[10]:


ar8['flac'] = (ar8['qual_flag']==0) & (ar8['flag_acaod']==1)  
ar8['flacr'] = (ar8['qual_flag']==0) & (ar8['flag_acaod']==1)&(ar8['fl_routine']) 
ar8['flaco'] = (ar8['qual_flag']==0) & (ar8['flag_acaod']==1)&~(ar8['fl_routine']) 
ar8['flr'] = (ar8['qual_flag']==0) & (ar8['fl_routine'])
ar8['flo'] = (ar8['qual_flag']==0) & (ar8['fl_routine']==False)
ar8['fl'] = (ar8['qual_flag']==0)


# ### Combined and set days

# In[11]:


days6 = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']
days7 = ['20170801','20170802','20170807','20170809', '20170812','20170813','20170815','20170817','20170818','20170819','20170821',
        '20170824','20170826','20170828','20170830','20170831','20170902','20170903','20170904']
days8 = ['20180921','20180922','20180924','20180927','20180930','20181002','20181003','20181005','20181007','20181010','20181012',
        '20181015','20181017','20181019','20181021','20181023','20181025','20181026','20181027']


# In[12]:


ar6['daysd'] = [days6[i] for i in ar6['days'].astype(int)]
ar7['daysd'] = [days7[i] for i in ar7['days'].astype(int)]
ar8['daysd'] = [days8[i] for i in ar8['days'].astype(int)]


# In[13]:


ar6['ndtime2'] = np.array([datetime(2018,int(d[4:6]),int(d[6:8]),int(ar6['Start_UTC'][i]),
                          int((ar6['Start_UTC'][i]-float(int(ar6['Start_UTC'][i])))*60)) for i,d in enumerate(ar6['daysd'])])
ar7['ndtime2'] = np.array([datetime(2018,int(d[4:6]),int(d[6:8]),int(ar7['Start_UTC'][i]),
                          int((ar7['Start_UTC'][i]-float(int(ar7['Start_UTC'][i])))*60)) for i,d in enumerate(ar7['daysd'])])
ar8['ndtime2'] = np.array([datetime(2018,int(d[4:6]),int(d[6:8]),int(ar8['Start_UTC'][i]),
                          int((ar8['Start_UTC'][i]-float(int(ar8['Start_UTC'][i])))*60)) for i,d in enumerate(ar8['daysd'])])


# ## Load HSRL ORACLES

# ### 2016

# In[14]:


vv_2016 = 'R9'


# In[15]:


fp6 = fpo+'data_other/HSRL/{vv}/'.format(vv=vv_2016)
f6_hsrl = os.listdir(fp6)
f6_hsrl = [ff for ff in f6_hsrl if (vv_2016 in ff) & ff.endswith('h5') ]
f6_hsrl.sort()
f6_hsrl


# In[16]:


hsrl_days6 = [f.split('_')[2] for f in f6_hsrl]
hsrl_doys6 = [datetime.strptime(d,'%Y%m%d').timetuple().tm_yday   for d in hsrl_days6]


# In[17]:


s6 = {}
for i,f in enumerate(f6_hsrl):
    print('Reading file: '+fp6+f)
    hf = hs.h5py.File(fp6+f,'r+')
    h = {}
    h[u'lat'] = hf[u'Nav_Data']['Latitude'][:]
    h[u'lon'] = hf[u'Nav_Data']['Longitude'][:]
    h[u'alt'] = hf[u'Nav_Data']['gps_alt'][:]
    h[u'time'] = hf[u'Nav_Data']['gps_time'][:]
    h[u'header'] = hf['000_Readme'][:]
    h[u'acaod_355'] = hf['DataProducts']['355_AOT_above_cloud'][:]
    h[u'acaod_532'] = hf['DataProducts']['532_AOT_above_cloud'][:]
    h[u'acaod_355_unc'] = hf['DataUncertainty']['355_AOT_above_cloud_unc'][:]
    h[u'acaod_532_unc'] = hf['DataUncertainty']['532_AOT_above_cloud_unc'][:]
    h[u'cloud_top_height'] = hf['DataProducts']['cloud_top_height'][:]
    h[u'date'] = hf['header']['date'][:][0][0]
    h[u'filename'] = f
    alt = hf['DataProducts']['Altitude'][:][0,:]
    bsc_cl = hf['DataProducts']['532_bsc_cloud_screened'][:]
    h[u'aero_bot_height'] = np.array([alt[bks>0.00025][0] if (any(bks) & any(bks>0.00025)) else np.nan for bks in bsc_cl])
    h[u'aero_bot_height'][h[u'aero_bot_height']<200.0] = np.nan 
    h[u'aero_bot_height'][h[u'aero_bot_height']>7000.0] = np.nan 
    h[u'gap_dist'] = h[u'aero_bot_height'] - h[u'cloud_top_height'][:,0]
    s6[u's{:08.0f}'.format(h['date'])] = h


# In[18]:


s6.keys()


# In[19]:


# combine into one array
hr6 = {}
k6 = s6.keys()
for n in s6[list(k6)[0]].keys():
    try:
        hr6[n] = np.array([],dtype=s6[list(k6)[0]][n].dtype)
    except:
        hr6[n] = np.array([])


# In[20]:


for i,d in enumerate(s6.keys()):
    hr6['date'] = np.append(hr6['date'],np.zeros_like(s6[d]['time'][:-1])+s6[d]['date'])
    print( len(np.zeros_like(s6[d]['time'])+s6[d]['date']),len(s6[d]['lat']))
    for n in s6[list(k6)[0]].keys():
        hr6[n] = np.append(hr6[n],s6[d][n])


# In[21]:


hr6.keys() 


# ### 2017

# In[62]:


vv_2017 = 'R1'
fp7 = fpo+'data_other_2017/HSRL/'
f7_hsrl = os.listdir(fp7)
f7_hsrl = [ff for ff in f7_hsrl if ff.endswith('h5') ]
f7_hsrl.sort()
f7_hsrl


# In[63]:


s7 = {}
for i,f in enumerate(f7_hsrl):
    print('Reading file: '+fp7+f)
    hf = hs.h5py.File(fp7+f,'r+')
    h = {}
    h[u'header'] = hf['000_Readme'][:]
    h[u'acaod_355'] = hf['DataProducts']['355_AOT_above_cloud'][:]
    h[u'acaod_532'] = hf['DataProducts']['532_AOT_above_cloud'][:]
  #  h[u'acaod_355_unc'] = hf['DataUncertainty']['355_AOT_above_cloud_unc'][:]
  #  h[u'acaod_532_unc'] = hf['DataUncertainty']['532_AOT_above_cloud_unc'][:]
    h[u'cloud_top_height'] = hf['DataProducts']['cloud_top_height'][:]
    h[u'date'] = hf['header']['date'][:][0][0]
    h[u'filename'] = f
    try:
        h[u'lat'] = hf['Nav_Data']['gps_lat'][:]
        h[u'lon'] = hf['Nav_Data']['gps_lon'][:]
        h[u'alt'] = hf['Nav_Data']['gps_alt'][:]
        h[u'time'] = hf['Nav_Data']['gps_time'][:]
    except:  
        h[u'lat'] = hf['ApplanixIMU']['gps_lat'][:]
        h[u'lon'] = hf['ApplanixIMU']['gps_lon'][:]
        h[u'alt'] = hf['ApplanixIMU']['gps_alt'][:]
        h[u'time'] = hf['ApplanixIMU']['gps_time'][:]
    h[u'fl'] = h['alt']>5000.0
    alt = hf['DataProducts']['Altitude'][:][0,:]
    bsc_cl = hf['DataProducts']['532_bsc_cloud_screened'][:]
    h[u'aero_bot_height'] = np.array([alt[bks>0.00025][0] if (any(bks) & any(bks>0.00025)) else np.nan for bks in bsc_cl])
    h[u'aero_bot_height'][h[u'aero_bot_height']<200.0] = np.nan 
    h[u'aero_bot_height'][h[u'aero_bot_height']>7000.0] = np.nan 
    h[u'gap_dist'] = h[u'aero_bot_height'] - h[u'cloud_top_height'][:,0]
    s7[u's{:08.0f}'.format(h['date'])] = h


# In[75]:


hr7 = {}
k7 = s7.keys()
for n in s7[list(k7)[0]].keys():
    try:
        hr7[n] = np.array([],dtype=s7[list(k7)[0]][n].dtype)    
    except:
        hr7[n] = np.array([])


# In[76]:


for i,d in enumerate(s7.keys()):
    hr7['date'] = np.append(hr7['date'],np.zeros_like(s7[d]['time'][:-1])+s7[d]['date'])
    print( len(np.zeros_like(s7[d]['time'])+s7[d]['date']),len(s7[d]['lat']))
    for n in s7[list(k7)[0]].keys():
        hr7[n] = np.append(hr7[n],s7[d][n])


# In[77]:


hr7.keys()


# ### 2018

# In[64]:


vv_2018 = 'R2'
fp8 = fpo+'data_other_2018/HSRL/'
f8_hsrl = []
for f8l in os.listdir(fp8):
    if f8l.endswith('.h5') and vv_2018 in f8l:
        f8_hsrl.append(f8l)
f8_hsrl.sort()
f8_hsrl


# In[65]:


s8 = {}
for i,f in enumerate(f8_hsrl):
    print( 'Reading file: '+fp8+f)
    hf = hs.h5py.File(fp8+f,'r+') 
    h = {}
    h[u'header'] = hf['000_Readme'][:]
    h[u'acaod_355'] = hf['DataProducts']['355_AOT_above_cloud'][:]
    h[u'acaod_532'] = hf['DataProducts']['532_AOT_above_cloud'][:]
  #  h[u'acaod_355_unc'] = hf['DataUncertainty']['355_AOT_above_cloud_unc'][:]
  #  h[u'acaod_532_unc'] = hf['DataUncertainty']['532_AOT_above_cloud_unc'][:]
    h[u'cloud_top_height'] = hf['DataProducts']['cloud_top_height'][:]
    h[u'date'] = hf['header']['date'][:][0][0]
    h[u'filename'] = f
    try:
        h[u'lat'] = hf['Nav_Data']['gps_lat'][:]
        h[u'lon'] = hf['Nav_Data']['gps_lon'][:]
        h[u'alt'] = hf['Nav_Data']['gps_alt'][:]
        h[u'time'] = hf['Nav_Data']['gps_time'][:]
    except:
        h[u'lat'] = hf['ApplanixIMU']['gps_lat'][:]
        h[u'lon'] = hf['ApplanixIMU']['gps_lon'][:]
        h[u'alt'] = hf['ApplanixIMU']['gps_alt'][:]
        h[u'time'] = hf['ApplanixIMU']['gps_time'][:]
    h[u'fl'] = h['alt'] >5000.0
    s8[u's{:08.0f}'.format(h['date'])] = h


# In[78]:


hr8 = {}
k8 = list(s8.keys())
for n in s8[k8[0]].keys():
    try:
        hr8[n] = np.array([],dtype=s8[k8[0]][n].dtype)
    except:
        hr8[n] = np.array([])


# In[79]:


for i,d in enumerate(s8.keys()):
    hr8['date'] = np.append(hr8['date'],np.zeros_like(s8[d]['time'][:-1])+s8[d]['date'])
    for n in s8[k8[0]].keys():
        hr8[n] = np.append(hr8[n],s8[d][n])


# In[80]:


hr8.keys()


# In[ ]:





# ## Load the AQUA

# In[22]:


fpa = fp+'AQUA'


# In[23]:


files_a = os.listdir(fpa)
files_a.sort()


# In[24]:


len(files_a)


# In[25]:


doys6 = [datetime.strptime(d,'%Y%m%d').timetuple().tm_yday   for d in days6]
doys7 = [datetime.strptime(d,'%Y%m%d').timetuple().tm_yday   for d in days7]
doys8 = [datetime.strptime(d,'%Y%m%d').timetuple().tm_yday   for d in days8]


# In[26]:


files_a6_p = [val for val in files_a if any(str(x) in val.split('.')[1] for x in doys6) ]
files_a6_e = [val for val in files_a if any(str(x) in val.split('.')[1] for x in hsrl_doys6) ]


# In[27]:


files_a6_p.sort()
files_a6_e.sort()


# In[28]:


a1, a1_dict = lu.load_hdf(fpa+'/'+files_a6_p[0],values=(('lat',0),('lon',1),('AAOD',10),('AAOD_UNC',11),
                                                        ('AAOD_2',16),('AAOD2_UNC',17)))


# In[29]:


a1.keys()


# In[30]:


a1['AAOD']


# In[31]:


nx,ny = a1['AAOD'].shape
nt = len(files_a6_p)


# In[32]:


f[:37]


# In[33]:


datetime.strptime(f[:37],'CLDACAERO_L2_MODIS_Aqua.A%Y%j.%H%M')


# In[34]:


rad = mu.spherical_dist([a1['lat'][0,0],a1['lon'][0,0]],[a1['lat'][1,0],a1['lon'][1,0]])*1000.0


# In[35]:


a1['lat'].shape


# In[36]:


a1_dict['AAOD_2']


# In[37]:


a1_dict['lat']


# In[38]:


warnings.filterwarnings('ignore')


# In[104]:


aaod = []
aaod_unc = []
aaod2 = []
aaod2_unc = []
lat = []
lon = []
a_time = []
a_daystr = []
for f in tqdm(files_a6_p):
    #print('loading file: {}'.format(f))
    a1, a1_dict = lu.load_hdf(fpa+'/'+f,values=(('lat',0),('lon',1),('AAOD',10),('AAOD_UNC',11),
                                                        ('AAOD_2',16),('AAOD2_UNC',17)),verbose=False)
    aaod.append(a1['AAOD'].flatten()*float(a1_dict['AAOD']['/geophysical_data/Above_Cloud_AOD#scale_factor']))
    aaod2.append( a1['AAOD_2'].flatten()*float(a1_dict['AAOD_2']['/geophysical_data/Above_Cloud_AOD_Secondary#scale_factor']))
    aaod_unc.append(a1['AAOD_UNC'].flatten()*float(a1_dict['AAOD_UNC']['/geophysical_data/Above_Cloud_AOD_Uncertainty#scale_factor']))
    aaod2_unc.append( a1['AAOD2_UNC'].flatten()*float(a1_dict['AAOD2_UNC']['/geophysical_data/Above_Cloud_AOD_Secondary_Uncertainty#scale_factor']))
    lat.append(a1['lat'].flatten())
    lon.append(a1['lon'].flatten())
    a_time.append(datetime.strptime(f[:37],'CLDACAERO_L2_MODIS_Aqua.A%Y%j.%H%M'))
    a_daystr.append(a_time[-1].strftime('%Y%m%d'))


# In[105]:


aaod_e = []
aaod_e_unc = []
aaod2_e = []
aaod2_e_unc = []
lat_e = []
lon_e = []
a_time_e = []
for f in tqdm(files_a6_e):
    #print('loading file: {}'.format(f))
    a1, a1_dict = lu.load_hdf(fpa+'/'+f,values=(('lat',0),('lon',1),('AAOD',10),('AAOD_UNC',11),
                                                        ('AAOD_2',16),('AAOD2_UNC',17)),verbose=False)
    aaod_e.append(a1['AAOD']*float(a1_dict['AAOD']['/geophysical_data/Above_Cloud_AOD#scale_factor']))
    aaod2_e.append( a1['AAOD_2']*float(a1_dict['AAOD_2']['/geophysical_data/Above_Cloud_AOD_Secondary#scale_factor']))
    aaod_e_unc.append(a1['AAOD_UNC']*float(a1_dict['AAOD_UNC']['/geophysical_data/Above_Cloud_AOD_Uncertainty#scale_factor']))
    aaod2_e_unc.append( a1['AAOD2_UNC']*float(a1_dict['AAOD2_UNC']['/geophysical_data/Above_Cloud_AOD_Secondary_Uncertainty#scale_factor']))
    lat_e.append(a1['lat'])
    lon_e.append(a1['lon'])
    a_time_e.append(datetime.strptime(f[:37],'CLDACAERO_L2_MODIS_Aqua.A%Y%j.%H%M'))


# In[106]:


a_daystr = [[a.strftime('%Y%m%d')]*len(lat[i]) for i,a in enumerate(a_time)]
a_daystr_e = [[a.strftime('%Y%m%d')]*len(lat_e[i]) for i,a in enumerate(a_time_e)]


# In[107]:


a_daystr = np.concatenate(a_daystr)
a_daystr_e = np.concatenate(a_daystr_e)


# In[108]:


aaod = np.concatenate(aaod)
aaod2 = np.concatenate(aaod2)
aaod_unc = np.concatenate(aaod_unc)
aaod2_unc = np.concatenate(aaod2_unc)
lat = np.concatenate(lat)
lon = np.concatenate(lon) 


# In[109]:


aaod[aaod<-0.1] = np.nan
aaod2[aaod2<-0.1] = np.nan
aaod_unc[aaod_unc<-0.1] = np.nan
aaod2_unc[aaod_unc<-0.1] = np.nan


# In[110]:


aaod_e = np.concatenate(aaod_e)
aaod2_e = np.concatenate(aaod2_e)
aaod_e_unc = np.concatenate(aaod_e_unc)
aaod2_e_unc = np.concatenate(aaod2_e_unc)
lat_e = np.concatenate(lat_e)
lon_e = np.concatenate(lon_e) 
aaod_e[aaod_e<-0.1] = np.nan
aaod2_e[aaod2_e<-0.1] = np.nan
aaod_e_unc[aaod_e_unc<-0.1] = np.nan
aaod2_e_unc[aaod2_e_unc<-0.1] = np.nan


# In[111]:


len(aaod_e)


# In[112]:


len(a_daystr_e)


# ## Load the VIIRS

# In[ ]:





# In[ ]:





# In[ ]:





# # Plot out the comparisons

# ## Get the overlaps and matches

# ### For  4STAR

# In[71]:


ar6.keys()


# In[72]:


ar6['mod_AAOD']=np.zeros_like(ar6['AOD0501'])
ar6['mod_AAOD2']=np.zeros_like(ar6['AOD0501'])
ar6['mod_AAOD_UNC']=np.zeros_like(ar6['AOD0501'])
ar6['mod_AAOD2_UNC']=np.zeros_like(ar6['AOD0501'])


# In[ ]:


for d in tqdm(days6,desc='days'):
    i_daymod = np.where(np.array(a_daystr)==d)[0]
    i_dayar = np.where((np.array(ar6['daysd'])==d) & np.isfinite(ar6['Latitude']))[0]
    print(d,len(i_daymod))
    for i in tqdm(i_daymod,desc='is'):
        if (lat[i].max()>np.nanmin(ar6['Latitude'][i_dayar])) & (lat[i].min()<np.nanmax(ar6['Latitude'][i_dayar])) &\
         (lon[i].max()>np.nanmin(ar6['Longitude'][i_dayar])) & (lon[i].min()<np.nanmax(ar6['Longitude'][i_dayar])):
            
            out = mu.stats_within_radius(ar6['Latitude'][i_dayar],ar6['Longitude'][i_dayar],
                                      lat[i],lon[i],aaod[i],rad,subset=False)
            ar6['mod_AAOD'][i_dayar] = out['mean']
            ar6['mod_AAOD2'][i_dayar] = np.array([np.nanmean(aaod2[i].flatten()[io]) for io in out['index']])
            ar6['mod_AAOD_UNC'][i_dayar] = np.array([np.nanmean(aaod_unc[i].flatten()[io]) for io in out['index']])
            ar6['mod_AAOD2_UNC'][i_dayar] = np.array([np.nanmean(aaod2_unc[i].flatten()[io]) for io in out['index']])


# In[82]:


## Save output of colocation


# In[84]:


dat = {'mod_AAOD':ar6['mod_AAOD'],
       'mod_AAOD2':ar6['mod_AAOD2'],
       'mod_AAOD_UNC':ar6['mod_AAOD_UNC'],
       'mod_AAOD2_UNC':ar6['mod_AAOD2_UNC'],
       'lat':ar6['Latitude'],
       'lon':ar6['Longitude'],
       'daystr':ar6['daysd']}
np.save(fpo+'mod_ACAERO_match_4star.npy',dat,allow_pickle=True)


# In[89]:


hs.savemat(fpo+'mod_ACAERO_match_4star.mat',dat)


# In[45]:


## Load output of colocation


# In[40]:


dat = np.load(fpo+'mod_ACAERO_match_4star.npy',allow_pickle=True)


# In[42]:


da = hs.loadmat(fpo+'mod_ACAERO_match_4star.mat')


# In[43]:


da


# ### For HSRL

# In[61]:


for sk in s6.keys():
    plt.figure()
    plt.plot(s6[sk]['time'],s6[sk]['acaod_355'],'b.',label='355 nm')
    plt.plot(s6[sk]['time'],s6[sk]['acaod_532'],'g.',label='532 nm')
    plt.title(sk)
    plt.legend(frameon=False,numpoints=1)
    plt.ylabel('AOD')
    plt.xlabel('Time [UTC]')
    ax = plt.gca()
    axy = ax.twinx()
    axy.plot(s6[sk]['time'],s6[sk]['alt'],'.',color='lightgrey')
    axy.set_ylabel('Altitude [m]')


# In[81]:


for sk in s7.keys():
    plt.figure()
    plt.plot(s7[sk]['time'],s7[sk]['acaod_355'],'b.',label='355 nm')
    plt.plot(s7[sk]['time'],s7[sk]['acaod_532'],'g.',label='532 nm')
    plt.title(sk)
    plt.legend(frameon=False,numpoints=1)
    plt.ylabel('AOD')
    plt.xlabel('Time [UTC]')
    ax = plt.gca()
    axy = ax.twinx()
    axy.plot(s7[sk]['time'],s7[sk]['alt'],'.',color='lightgrey')
    axy.set_ylabel('Altitude [m]')


# In[82]:


for sk in s8.keys():
    plt.figure()
    plt.plot(s8[sk]['time'],s8[sk]['acaod_355'],'b.',label='355 nm')
    plt.plot(s8[sk]['time'],s8[sk]['acaod_532'],'g.',label='532 nm')
    plt.title(sk)
    plt.legend(frameon=False,numpoints=1)
    plt.ylabel('AOD')
    plt.xlabel('Time [UTC]')
    ax = plt.gca()
    axy = ax.twinx()
    axy.plot(s8[sk]['time'],s8[sk]['alt'],'.',color='lightgrey')
    axy.plot(s8[sk]['time'][s8[sk]['fl']],s8[sk]['alt'][s8[sk]['fl']],'.',color='grey')
    axy.set_ylabel('Altitude [m]')


# In[86]:


hr6.keys()


# In[89]:


len(hr6['date'])


# In[94]:


hr6['daysd'] = ['{:8.0f}'.format(d) for d in hr6['date']]


# In[97]:


len(a_daystr_e)


# In[95]:


for d in tqdm(a_daystr_e,desc='days'):
    i_daymod = np.where(np.array(a_daystr_e)==d)[0]
    i_dayar = np.where((np.array(hr6['daysd'])==d) & np.isfinite(hr6['lat']))[0]
    print(d,len(i_daymod))
    for i in tqdm(i_daymod,desc='is'):
        if (lat_e[i].max()>np.nanmin(hr6['lat'][i_dayar])) & (lat_e[i].min()<np.nanmax(hr6['lat'][i_dayar])) &\
         (lon_e[i].max()>np.nanmin(hr6['lon'][i_dayar])) & (lon_e[i].min()<np.nanmax(hr6['lon'][i_dayar])):
            
            out = mu.stats_within_radius(hr6['lat'][i_dayar],hr6['lon'][i_dayar],
                                      lat_e[i],lon_e[i],aaod_e[i],rad,subset=False)
            hr6['mod_AAOD'][i_dayar] = out['mean']
            hr6['mod_AAOD2'][i_dayar] = np.array([np.nanmean(aaod2_e[i].flatten()[io]) for io in out['index']])
            hr6['mod_AAOD_UNC'][i_dayar] = np.array([np.nanmean(aaod_e_unc[i].flatten()[io]) for io in out['index']])
            hr6['mod_AAOD2_UNC'][i_dayar] = np.array([np.nanmean(aaod2_e_unc[i].flatten()[io]) for io in out['index']])


# ## Plot the scatter plots

# ### 4STAR

# In[53]:


ar6.keys()


# In[44]:


fl = (ar6['flag_acaod']==1) & ar6['fl']
flo = (ar6['flag_acaod']==1) & ar6['fl'] & (da['mod_AAOD']>-0.2) & np.isfinite(da['mod_AAOD']) & np.isfinite(ar6['AOD0550'])


# In[58]:


plt.figure()
plt.hist([da['mod_AAOD'][flo],ar6['AOD0550'][flo]],bins=50,label=['MOD ACAERO','4STAR AAOD'])
plt.legend()
plt.xlabel('AOD [550 nm]')
plt.savefig(fpo+'MODACAERO_vs_4STAR_hist.png',dpi=600,transparent=True)


# In[80]:


ar6.keys()


# In[81]:


da.keys()


# In[47]:


fig,ax = plt.subplots(1,1)
plt.plot(ar6['AOD0550'][flo],da['mod_AAOD'][flo],'.')
plt.errorbar(ar6['AOD0550'][flo],da['mod_AAOD'][flo],xerr=ar6['UNCAOD0550'][flo],yerr=da['mod_AAOD_UNC'][flo],
             marker='.',color='tab:blue',linestyle='None',alpha=0.2)
plt.xlabel('4STAR ACAOD [550 nm]')
plt.ylabel('MODIS ACAERO (v1) [550m]')
pu.plot_lin(ar6['AOD0550'][flo],da['mod_AAOD'][flo],
            x_err=ar6['UNCAOD0550'][flo],y_err=da['mod_AAOD_UNC'][flo],
            labels=True,shaded_ci=True,ci=95,ax=ax,
            use_method='york',lblfmt='2.3f',label_prefix='Bivariate fit: [York et al., 2004]\n')
plt.title('ORACLES 2016 Above Cloud AOD (MODIS with 4STAR abs. vs. 4STAR)')
plt.legend()
plt.savefig(fpo+'MODACAERO_vs_4STAR_ORACLES2016.png',dpi=600,transparent=True)


# In[61]:


fpo


# In[ ]:




