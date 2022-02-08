#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Explore some MODIS cloud fraction and subsequent fire counts
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
#     - Sp_parameters
#     - write_utils
#     - path_utils
#     - hdf5storage
#     - scipy
# 
# Needed Files:
#   - MYD06 and MYD14 hdf files
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-05-11
#     Modified:
# 

# # Prepare python environment

# In[37]:


nomap = True
import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt

from datetime import datetime, timedelta

import os
if not nomap:
    from mpl_toolkits.basemap import Basemap
    get_ipython().magic(u'matplotlib notebook')


# In[38]:


name = 'FOG2FIRE'
vv = 'v3'
fp = getpath(name)
yy = '2020'


# In[116]:


import argparse
long_description = """    Quantify the cloud fraction and or the number of fire counts in the files MYD06 and MYD14. 
     Subsets and saves to a npy file for clouds or fires"""
parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-c','--cloud',help='if set, calc the cloud',action='store_true')
parser.add_argument('-f','--fires',help='if set, run the fires calcs',action='store_true')
parser.add_argument('-s','--smap',help='if set, run the smap soil moisture calcs',action='store_true')
parser.add_argument('-y','--year',nargs='?',help='year',default='2020')
parser.add_argument('-m','--modis',nargs='?',help='MYD or MOD for Terra or Aqua',default='MOD')


# In[123]:


in_ = vars(parser.parse_known_args()[0])
run_cloud = in_.get('cloud',False)
run_fire = in_.get('fires',False)
run_smap = in_.get('smap',False)
yy = in_.get('year').strip()
modis = in_.get('modis').strip()


# In[124]:


print(in_)


# Web link to download the MODIS cloud files  
# https://ladsweb.modaps.eosdis.nasa.gov/search/order/3/MOD06_L2--61/2004-01-01..2004-12-31/D/-125,47,-115,32.5

# In[119]:


from osgeo import gdal
if gdal.__version__ > '3.2':
    load_special = True
else:
    load_special = False


# # Plan out the regions

# In[7]:


if not nomap:
    def make_map(ax=plt.gca()):
        m = Basemap(projection='stere',lon_0=-122.0,lat_0=38.0,
                llcrnrlon=-131.0, llcrnrlat=32.0,
                urcrnrlon=-108.0, urcrnrlat=48,resolution='i',ax=ax)
        m.drawcoastlines()
        #m.fillcontinents(color='#AAAAAA')
        m.drawstates()
        m.drawcountries()
        m.drawmeridians(np.linspace(-131,-108,6),labels=[0,0,0,1])
        m.drawparallels(np.linspace(31,48,9),labels=[1,0,0,0])
        def format_coord(x, y):
            return 'x=%.4f, y=%.4f'%(m(x, y, inverse = True))
        ax.format_coord = format_coord
        return m
    make_map()


# In[120]:


rgs = [[[32.5,-121.5],[35.5,-117.0]],
       [[32.5,-117.0],[35.5,-114.0]],
       [[35.5,-123.5],[38.5,-120.8]],
       [[35.5,-120.8],[38.5,-115.0]],
       [[38.5,-125.0],[42.0,-122.0]],
       [[38.5,-122.0],[42.0,-118.0]],
       [[42.0,-125.0],[47.0,-122.0]],
       [[42.0,-122.0],[47.0,-115.0]]]
       #lower left [lat lon], upper right [lat lon]
lbls = ['Socal Coast','Socal land','Central coast','Central Sierras',
        'Norcal coast','Northern Sierras','Oregon Coast','Oregon mountains']


# In[121]:


regions = {'ocean':[[[32.5,-131],[35.5,-121.5]],[[35.5,-131.0],[38.5,-123.5]], [[38.5,-131.0],[42.0,-125.0]],[[42.0,-131.0],[47.0,-125.0]]],
           'coast':[[[32.5,-121.5],[35.5,-117.0]],[[35.5,-123.5],[38.5,-120.8]], [[38.5,-125.0],[42.0,-122.0]],[[42.0,-125.0],[47.0,-122.0]]],
           'land':[[[32.5,-117.0],[35.5,-114.0]],[[35.5,-120.8],[38.5,-115.0]],[[38.5,-122.0],[42.0,-118.0]],[[42.0,-122.0],[47.0,-115.0]]],
           'points':[[[32.5,-117.22],[33.6,-116.7]],[[34.02,-118.18],[34.5,-116.63]],[[34.35,-120.7],[35.08,-118.7]],
                     [[35.37,-121.97],[36.6,-120.69]],[[35.6,-120.01],[38.43,-118.23]],[[36.9,-122.56],[37.68,-121.56]],
                     [[38.24,-121.25],[40.09,-119.92]],[[38.45,-123.98],[40.39,-122.23]],[[40.0,-122.57],[41.73,-120.37]],
                     [[40.77,-124.0],[42.0,-122.54]],[[42.05,-124.51],[43.68,-123.10]],[[43.95,-122.92],[42.82,-120.82]],
                     [[44.04,-124.15],[46.10,-122.92]],[[46.97,-124.87],[48.40,-122.58]],[[46.52,-122.69],[48.59,-120.65]]]
          }

lbls_rg = {'ocean':['SoCal','Central','NorCal','Oregon'],
           'coast':['SoCal','Central','NorCal','Oregon'],
           'land':['SoCal','Central Sierras','Northern Sierras','Oregon mountains'],
           'points':['San Diego','San Bernardino','Los Padres',
                     'Big Sur','Sierra','Santa Cruz',
                     'Eldorado','Mendocino','Lassen-Shasta',
                     'Klamath','Rogue','Mt Hood',
                     'Portland','Olympic','Mt. Baker-Rainier']
          }
ls = {'ocean':':','coast':'-','land':'--','points':'-.'}


# In[32]:


if not nomap:
    fig, ax = plt.subplots(1,2,figsize=(9,4))
    m = make_map(ax[0])
    ax[0].set_prop_cycle(color=[plt.cm.gist_ncar(k) for k in np.linspace(0, 1, len(i_rg[0]))])
    for re in regions:
        for i,r in enumerate(regions[re]):
            m.plot([r[0][1],r[1][1],r[1][1],r[0][1],r[0][1]],[r[1][0],r[1][0],r[0][0],r[0][0],r[1][0]],
                   latlon=True,label=re+'-'+lbls_rg[re][i],lw=2,ls=ls[re])
    ax[0].legend(bbox_to_anchor=[1.0,0.95],ncol=2,loc=2)
    #plt.tight_layout(rect=[0.1,-0.4,0.95,1.5])
    ax[1].set_visible(False)

    plt.savefig(fp+'plots/Map_regions_{}.png'.format(vv),dpi=400,transparent=True)


# # Load files

# ## Load the Cloud files

# In[107]:


if run_cloud:
    print('listdir')
    lc_all = os.listdir(fp+'{}06/{}/'.format(modis,yy))
    lc = [o for o in lc_all if 'A'+yy in o]
    lc.sort()
    nfiles = len(lc)
    ifile = 0
    
    if load_special:
        vals = (('CF','HDF4_SDS:UNKNOWN:"{}":36'),
             ('CF_night','HDF4_SDS:UNKNOWN:"{}":38'),
             ('CF_day','HDF4_SDS:UNKNOWN:"{}":40'),
             ('cld_top','HDF4_SDS:UNKNOWN:"{}":19'),
             ('scan_time','HDF4_SDS:UNKNOWN:"{}":2'),
             ('sza','HDF4_SDS:UNKNOWN:"{}":3'),
             ('surf_temp','HDF4_SDS:UNKNOWN:"{}":16'),
             ('lat','HDF4_SDS:UNKNOWN:"{}":0'),
             ('lon','HDF4_SDS:UNKNOWN:"{}":1'),
             ('QA','HDF4_SDS:UNKNOWN:"{}":111'),
             ('cld_mask','HDF4_SDS:UNKNOWN:"{}":110'))
    else:
        vals = (('CF',160),('CF_night',162),('CF_day',164),('cld_top',143),('scan_time',126),
                ('sza',127),('surf_temp',140),('lat',124),('lon',125),('QA',235),('cld_mask',234))


# In[96]:


if run_cloud:
    print('loading of cloud files : {}'.format(nfiles))
    clds = []
    for i,l in enumerate(lc):
        if not l.endswith('hdf'): continue
        print('loading file: {}/{}, {}'.format(i,nfiles,l))
        try:
            if load_special:
                cld,cld_dict = lu.load_hdf_withsub(fp+'{}06/{}/'.format(modis,yy)+l,values=vals,verbose=False)
            else:
                cld,cld_dict = lu.load_hdf(fp+'{}06/{}/'.format(modis,yy)+l,values=vals,verbose=False)
            cld['surf_temp'] = (cld['surf_temp']+15000)*0.00999999977648258
            clds.append(cld)
        except:
            pass


# In[80]:


if not nomap:
    fig, ax = plt.subplots(1,1)
    m = make_map(ax)
    mcf = m.pcolor(cld['lon'],cld['lat'],cld['CF'],latlon=True)
    plt.colorbar(mcf,label='CF')


# ## Load the fire files

# In[346]:


if run_fire:
    print('listdir')
    lc_all = os.listdir(fp+'MYD14/')
    lc = [o for o in lc_all if 'A'+yy in o]
    lc.sort()
    nfiles = len(lc)
    ifile = 0


# In[347]:


if run_fire:
    fires,time_fires = [],[]
    print('loading of fire files : {}'.format(nfiles))
    for i,l in enumerate(lc):
        if not l.endswith('hdf'): continue
        print('loading file: {}/{}, {}'.format(i,nfiles,l))
        try: 
            fir,fir_dict = lu.load_hdf_sd(fp+'MYD14/'+l,vals=['FP_longitude','FP_latitude','FP_power'],verbose=False)
        except ValueError: 
            fir = {'FP_latitude':[],'FP_longitude':[],'FP_power':[]}
        tm = datetime.strptime(l[7:19],'%Y%j.%H%M')
        fires.append(fir)
        time_fires.append(tm)


# In[348]:


if run_fire:
    fires


# ## Load the SMAP soil moisture files

# Loading files from SMAP:  
# https://nsidc.org/data/SPL2SMP/versions/8
# 
# Chan, S., R. Bindlish, P. E. O'Neill, E. G. Njoku, T. Jackson, A. Colliander, F. Chen, M. Burgin, S. Dunbar, J. R. Piepmeier, S. Yueh, D. Entekhabi, M. Cosh, T. Caldwell, J. Walker, A. Berg, T. Rowlandson, A. Pacheco, H. McNairn, M. Thibeault, J. Martinez-Fernandez, A. González-Zamora, D. Bosch, P. Starks, D. Goodrich, J. Prueger, M. Palecki, E. E. Small, M. Zreda, J. Calvet, W. T. Crow, and Y. Kerr. 2016. Assessment of the SMAP passive soil moisture product, IEEE Transactions on Geoscience and Remote Sensing. 54. 4994–5007. https://doi.org/10.1109/TGRS.2016.2561938
#   
# and   
#   
# O'Neill, P. E., S. Chan, E. G. Njoku, T. Jackson, R. Bindlish, and J. Chaubell. 2021. SMAP L2 Radiometer Half-Orbit 36 km EASE-Grid Soil Moisture, Version 8. [Indicate subset used]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. doi: https://doi.org/10.5067/LPJ8F0TAK6E0. [Accessed 2022-02-08]. (SMAP_L2_SM_P_36065_D)

# In[161]:


print('listdir')
lc_all = os.listdir(fp+'SMAP/{}/'.format(yy))
lc = [o for o in lc_all if (yy in o) and (o.endswith('h5'))]
lc.sort()
nfiles = len(lc)
ifile = 0


# In[162]:


if run_smap:
    vals = (('lat',20),('lon',22),('QA',24),('soil_moist',32),('veg_wat',48),('time',36))
    smaps = []
    smaps_time = []
    print('loading of smap files : {}'.format(nfiles))
    for i,l in enumerate(lc):
        print('loading file: {}/{}, {}'.format(i,nfiles,l))
        sma,sma_dict = lu.load_hdf(fp+'SMAP/{}/'.format(yy)+l,verbose=False,values=vals,i_subdata=0)
        sma['soil_moist'][sma['soil_moist']<-9990.0] = np.nan
        smaps_time.append(datetime(2000,1,1)+timedelta(seconds=sma['time'][0,0]))
        smaps.append(sma)
    smaps_time = np.array(smaps_time)


# ## Define stats and get from regions

# In[154]:


if run_cloud:
    time = []
    for i,c in enumerate(clds):
        time.append(datetime(1993,1,1,0,0,0)+timedelta(seconds=c['scan_time'][0,0]))
    time = np.array(time)


# In[163]:


def stats(lon,lat,data,rg):
    'returns array of mean, median, std, and num for define area'
    if len(lon)>0:
        i = (lon >= rg[0][1]) & (lon <= rg[1][1]) & (lat >= rg[0][0]) & (lon <= rg[1][0])
        me = np.nanmean(data[i])
        md = np.nanmedian(data[i])
        st = np.nanstd(data[i])
        nu = len(np.isfinite(data[i]))
    else:
        me,md,st,nu = 0.0,0.0,0.0,0
    return [me,md,st,nu]


# In[196]:


if run_cloud:
    print('doing CF')
    cf = {}
    for re in regions:
        nre = len(regions[re])
        cf[re] = {u'mean':np.zeros((nfiles,nre))+np.nan,
                  u'median':np.zeros((nfiles,nre))+np.nan,
                  u'std':np.zeros((nfiles,nre))+np.nan,
                  u'num':np.zeros((nfiles,nre))+np.nan}
        for i,r in enumerate(regions[re]):
            for ifile in range(nfiles):
                cf[re]['mean'][ifile,i],cf[re]['median'][ifile,i],cf[re]['std'][ifile,i],cf[re]['num'][ifile,i] =                   stats(clds[ifile]['lon'],clds[ifile]['lat'],clds[ifile]['CF'],r)


# In[197]:


if run_cloud:
    print('doing cloud top')
    cld_top = {}
    for re in regions:
        nre = len(regions[re])
        cld_top[re] = {u'mean':np.zeros((nfiles,nre))+np.nan,
                  u'median':np.zeros((nfiles,nre))+np.nan,
                  u'std':np.zeros((nfiles,nre))+np.nan,
                  u'num':np.zeros((nfiles,nre))+np.nan}
        for i,r in enumerate(regions[re]):
            for ifile in range(nfiles):
                igood = (clds[ifile]['cld_top']>0.0)
                cld_top[re]['mean'][ifile,i],cld_top[re]['median'][ifile,i],                cld_top[re]['std'][ifile,i],cld_top[re]['num'][ifile,i] =                   stats(clds[ifile]['lon'][igood],clds[ifile]['lat'][igood],clds[ifile]['cld_top'][igood],r)


# In[198]:


if run_cloud:
    print('doing surface temp')
    surf_temp = {}
    for re in regions:
        nre = len(regions[re])
        surf_temp[re] = {u'mean':np.zeros((nfiles,nre))+np.nan,
                  u'median':np.zeros((nfiles,nre))+np.nan,
                  u'std':np.zeros((nfiles,nre))+np.nan,
                  u'num':np.zeros((nfiles,nre))+np.nan}
        for i,r in enumerate(regions[re]):
            for ifile in range(nfiles):
                surf_temp[re]['mean'][ifile,i],surf_temp[re]['median'][ifile,i],                surf_temp[re]['std'][ifile,i],surf_temp[re]['num'][ifile,i] =                   stats(clds[ifile]['lon'],clds[ifile]['lat'],clds[ifile]['surf_temp'],r)


# In[161]:


if run_cloud:
    print('doing sza')
    sza = {}
    for re in regions:
        nre = len(regions[re])
        sza[re] = np.zeros((nfiles,nre))+np.nan
        for i,r in enumerate(regions[re]):
            for ifile in range(nfiles):
                sza[re][ifile,i],md,st,nu =                   stats(clds[ifile]['lon'],clds[ifile]['lat'],clds[ifile]['sza'],r)


# In[183]:


if run_cloud:
    surf_temp


# ## Define stats for fire

# In[350]:


if run_fire:
    print('doing fire counts')
    fire_counts = {}
    for re in regions:
        nre = len(regions[re])
        fire_counts[re] = {u'mean':np.zeros((nfiles,nre))+np.nan,
                  u'median':np.zeros((nfiles,nre))+np.nan,
                  u'std':np.zeros((nfiles,nre))+np.nan,
                  u'num':np.zeros((nfiles,nre))+np.nan}
        for i,r in enumerate(regions[re]):
            for ifile in range(nfiles):
                fire_counts[re]['mean'][ifile,i],fire_counts[re]['median'][ifile,i],                fire_counts[re]['std'][ifile,i],fire_counts[re]['num'][ifile,i] =                   stats(fires[ifile]['FP_longitude'],fires[ifile]['FP_latitude'],fires[ifile]['FP_power'],r)


# ## Define stats for SMAP

# In[165]:


if run_smap:
    print('doing SMAP soil moisture')
    soil_moist = {}
    for re in regions:
        nre = len(regions[re])
        soil_moist[re] = {u'mean':np.zeros((nfiles,nre))+np.nan,
                  u'median':np.zeros((nfiles,nre))+np.nan,
                  u'std':np.zeros((nfiles,nre))+np.nan,
                  u'num':np.zeros((nfiles,nre))+np.nan}
        for i,r in enumerate(regions[re]):
            for ifile in range(nfiles):
                soil_moist[re]['mean'][ifile,i],soil_moist[re]['median'][ifile,i],                soil_moist[re]['std'][ifile,i],soil_moist[re]['num'][ifile,i] =                   stats(smaps[ifile]['lon'],smaps[ifile]['lat'],smaps[ifile]['soil_moist'],r)


# # Prep for saving

# In[166]:


if run_cloud:
    data = {'CF':cf,'surf_temp':surf_temp,'sza':sza,'cld_top':cld_top,'time':time,'regions':regions,'lbls_rg':lbls_rg}
if run_fire:
    data = {'FP':fire_counts,'time':time_fires,'regions':regions,'lbls_rg':lbls_rg}
if run_smap:
    data = {'SM':soil_moist,'time':smaps_time,'regions':regions,'lbls_rg':lbls_rg}


# In[167]:


import write_utils as wu
data = wu.iterate_dict_unicode(data)


# In[354]:


if run_cloud: fsuff = '{}06'.format(modis)
if run_fire: fsuff = 'MYD14' 
if run_smap: fsuff = 'SMAP'
print('Saving file to '+fp+'{}_{}_{}.npy'.format(fsuff,yy,vv))
np.save(fp+'{}_{}_{}.npy'.format(fsuff,yy,vv),data,allow_pickle=True)


# # For determining arbitrary if point is in arbitrary polygon

# In[ ]:


if False:
    import fiona
    from shapely.geometry import MultiPoint, Point, Polygon,shape
    from shapely.geometry.polygon import Polygon

    multipol = fiona.open(r"C:\Users\Jordi\Downloads\ESP_adm_shp\ESP_adm0.shp")
    multi = next(iter(multipol))

    point = Point(0,42)
    point.within(shape(multi['geometry']))

