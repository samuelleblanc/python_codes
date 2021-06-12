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

# In[230]:


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


# In[232]:


name = 'FOG2FIRE'
vv = 'v2'
fp = getpath(name)
yy = '2020'


# In[314]:


import argparse
long_description = """    Quantify the cloud fraction and or the number of fire counts in the files MYD06 and MYD14. 
     Subsets and saves to a npy file for clouds or fires"""
parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-c','--cloud',help='if set, calc the cloud',action='store_true')
parser.add_argument('-f','--fires',help='if set, run the fires calcs',action='store_true')
parser.add_argument('-y','--year',nargs='?',help='year',default='2020')


# In[344]:


in_ = vars(parser.parse_known_args(['-y','2020','-f'])[0])
run_cloud = in_.get('cloud',False)
run_fire = in_.get('fires',False)
yy = in_.get('year').strip()


# In[345]:


print(in_)


# # Plan out the regions

# In[41]:


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


# In[43]:


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


# In[49]:


regions = {'ocean':[[[32.5,-131],[35.5,-121.5]],[[35.5,-131.0],[38.5,-123.5]], [[38.5,-131.0],[42.0,-125.0]],[[42.0,-131.0],[47.0,-125.0]]],
           'coast':[[[32.5,-121.5],[35.5,-117.0]],[[35.5,-123.5],[38.5,-120.8]], [[38.5,-125.0],[42.0,-122.0]],[[42.0,-125.0],[47.0,-122.0]]],
           'land':[[[32.5,-117.0],[35.5,-114.0]],[[35.5,-120.8],[38.5,-115.0]],[[38.5,-122.0],[42.0,-118.0]],[[42.0,-122.0],[47.0,-115.0]]]
          }

lbls_rg = {'ocean':['SoCal','Central','NorCal','Oregon'],
           'coast':['SoCal','Central','NorCal','Oregon'],
           'land':['SoCal','Central Sierras','Northern Sierras','Oregon mountains']
          }
ls = {'ocean':':','coast':'-','land':'--'}


# In[56]:


if not nomap:
    fig, ax = plt.subplots(1,1)
    m = make_map(ax)
    for re in regions:
        for i,r in enumerate(regions[re]):
            m.plot([r[0][1],r[1][1],r[1][1],r[0][1],r[0][1]],[r[1][0],r[1][0],r[0][0],r[0][0],r[1][0]],
                   latlon=True,label=re+'-'+lbls_rg[re][i],lw=4,ls=ls[re])
    plt.legend(bbox_to_anchor=[1.0,0.7])
    plt.tight_layout(rect=[0.1,-0.4,0.95,1.5])

    plt.savefig(fp+'plots/Map_regions.png',dpi=400,transparent=True)


# # Load files

# ## Load the Cloud files

# In[62]:


if run_cloud:
    print('listdir')
    lc_all = os.listdir(fp+'MYD06/')
    lc = [o for o in lc_all if 'A'+yy in o]
    lc.sort()
    nfiles = len(lc)
    ifile = 0
    vals = (('CF',160),('CF_night',162),('CF_day',164),('cld_top',143),('scan_time',126),
            ('sza',127),('surf_temp',140),('lat',124),('lon',125),('QA',235),('cld_mask',234))


# In[136]:


if run_cloud:
    print('loading of cloud files : {}'.format(nfiles))
    clds = []
    for i,l in enumerate(lc):
        if not l.endswith('hdf'): continue
        print('loading file: {}/{}, {}'.format(i,nfiles,l))
        try:
            cld,cld_dict = lu.load_hdf(fp+'MYD06/'+l,values=vals,verbose=False)
            cld['surf_temp'] = (cld['surf_temp']+15000)*0.00999999977648258
            clds.append(cld)
        except:
            pass


# In[73]:


if run_cloud:
    cld_dict['CF']


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


# ## Define stats and get from regions

# In[154]:


if run_cloud:
    time = []
    for i,c in enumerate(clds):
        time.append(datetime(1993,1,1,0,0,0)+timedelta(seconds=c['scan_time'][0,0]))
    time = np.array(time)


# In[349]:


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


# # Prep for saving

# In[351]:


if run_cloud:
    data = {'CF':cf,'surf_temp':surf_temp,'sza':sza,'cld_top':cld_top,'time':time,'regions':regions,'lbls_rg':lbls_rg}
if run_fire:
    data = {'FP':fire_counts,'time':time_fires,'regions':regions,'lbls_rg':lbls_rg}


# In[352]:


import write_utils as wu
data = wu.iterate_dict_unicode(data)


# In[354]:


if run_cloud: fsuff = 'MYD06'
if run_fire: fsuff = 'MYD14' 
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

