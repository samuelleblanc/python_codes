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

# In[1]:


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
from mpl_toolkits.basemap import Basemap
get_ipython().magic(u'matplotlib notebook')


# In[11]:


get_ipython().magic(u'matplotlib notebook')


# In[2]:


name = 'FOG2FIRE'
vv = 'v1'
fp = getpath(name)


# In[227]:


nomap = True


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


lc = os.listdir(fp+'MYD06/')


# In[228]:


lc.sort()


# In[125]:


nfiles = len(lc)
ifile = 0


# In[134]:


vals = (('CF',160),('CF_night',162),('CF_day',164),('cld_top',143),('scan_time',126),
('sza',127),('surf_temp',140),('lat',124),('lon',125),('QA',235),('cld_mask',234))


# In[136]:


clds = []
for i,l in enumerate(lc):
    cld,cld_dict = lu.load_hdf(fp+'MYD06/'+l,values=vals,verbose=False)
    clds.append(cld)


# In[181]:


for c in clds:
    c['surf_temp'] = (c['surf_temp']+15000)*0.00999999977648258


# In[73]:


cld_dict['CF']


# In[74]:


cld['CF']


# In[80]:


if not nomap:
    fig, ax = plt.subplots(1,1)
    m = make_map(ax)
    mcf = m.pcolor(cld['lon'],cld['lat'],cld['CF'],latlon=True)
    plt.colorbar(mcf,label='CF')


# ## Define stats and get from regions

# In[154]:


time = []
for i,c in enumerate(clds):
    time.append(datetime(1993,1,1,0,0,0)+timedelta(seconds=c['scan_time'][0,0]))
time = np.array(time)


# In[127]:


def stats(lon,lat,data,rg):
    'returns array of mean, median, std, and num for define area'
    i = (lon >= rg[0][1]) & (lon <= rg[1][1]) & (lat >= rg[0][0]) & (lon <= rg[1][0])
    me = np.nanmean(data[i])
    md = np.nanmedian(data[i])
    st = np.nanstd(data[i])
    nu = len(np.isfinite(data[i]))
    return [me,md,st,nu]


# In[196]:


cf = {}
for re in regions:
    nre = len(regions[re])
    cf[re] = {u'mean':np.zeros((nfiles,nre))+np.nan,
              u'median':np.zeros((nfiles,nre))+np.nan,
              u'std':np.zeros((nfiles,nre))+np.nan,
              u'num':np.zeros((nfiles,nre))+np.nan}
    for i,r in enumerate(regions[re]):
        for ifile in range(nfiles):
            cf[re]['mean'][ifile,i],cf[re]['median'][ifile,i],cf[re]['std'][ifile,i],cf[re]['num'][ifile,i] =               stats(clds[ifile]['lon'],clds[ifile]['lat'],clds[ifile]['CF'],r)


# In[138]:


cf


# In[197]:


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
            cld_top[re]['mean'][ifile,i],cld_top[re]['median'][ifile,i],            cld_top[re]['std'][ifile,i],cld_top[re]['num'][ifile,i] =               stats(clds[ifile]['lon'][igood],clds[ifile]['lat'][igood],clds[ifile]['cld_top'][igood],r)


# In[198]:


surf_temp = {}
for re in regions:
    nre = len(regions[re])
    surf_temp[re] = {u'mean':np.zeros((nfiles,nre))+np.nan,
              u'median':np.zeros((nfiles,nre))+np.nan,
              u'std':np.zeros((nfiles,nre))+np.nan,
              u'num':np.zeros((nfiles,nre))+np.nan}
    for i,r in enumerate(regions[re]):
        for ifile in range(nfiles):
            surf_temp[re]['mean'][ifile,i],surf_temp[re]['median'][ifile,i],            surf_temp[re]['std'][ifile,i],surf_temp[re]['num'][ifile,i] =               stats(clds[ifile]['lon'],clds[ifile]['lat'],clds[ifile]['surf_temp'],r)


# In[161]:


sza = {}
for re in regions:
    nre = len(regions[re])
    sza[re] = np.zeros((nfiles,nre))+np.nan
    for i,r in enumerate(regions[re]):
        for ifile in range(nfiles):
            sza[re][ifile,i],md,st,nu =               stats(clds[ifile]['lon'],clds[ifile]['lat'],clds[ifile]['sza'],r)


# In[183]:


surf_temp


# # Prep for saving

# In[199]:


data = {'CF':cf,'surf_temp':surf_temp,'sza':sza,'cld_top':cld_top,'time':time,'regions':regions,'lbls_rg':lbls_rg}


# In[200]:


import write_utils as wu
data = wu.iterate_dict_unicode(data)


# In[207]:


np.save(fp+'MYD06_{}.npy'.format(vv),data,allow_pickle=True)

