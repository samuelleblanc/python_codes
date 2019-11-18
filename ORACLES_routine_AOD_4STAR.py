#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_routine_4STAR
# 
# Purpose:  
# 
#     Compare the routine flights AOD from 4STAR as compared to climatology
#   
# Input:
# 
#     none
# 
# Output:
#    
#     plots
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - matplotlib
#     - scipy
# 
#   
# Needed Files:
# 
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc,Santa Cruz, CA, 2016-11-08
#     

# # Prepare the python environment
# 

# In[1]:


import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth

import hdf5storage as hs
from mpltools import color
from path_utils import getpath
from write_utils import nearest_neighbor
from tqdm import tqdm_notebook as tqdm


# In[4]:


import plotting_utils as pu


# In[5]:


import scipy.stats as st
import Sun_utils as su


# In[6]:


from mpl_toolkits.basemap import Basemap


# In[7]:


fp = getpath('ORACLES')
fp


# In[8]:


vv = 'R3'


# # Load files

# ## 4STAR ict

# In[9]:


s = hs.loadmat(fp+'/aod_ict/v8/{vv}/all_aod_ict_{vv}_2016.mat'.format(vv=vv))


# In[10]:


s.keys()


# In[10]:


days = np.unique(s['days'])


# In[11]:


len(s['fl_QA'])


# In[12]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[13]:


days = [float(d) for d in days]


# ### Load the acaod flags

# In[14]:


fdat = getpath('4STAR_data')
fdat


# In[15]:


flag_acaod = []
flag_acaod_t = []

for j in days:
    ff = ''
    with open (fdat+'starinfo_{:8.0f}.m'.format(j), 'rt') as in_file:
        for line in in_file:
            if 'flagacaod' in line:
                ff = line.split("'")[1]
    if ff:
        flaod = sio.loadmat(fdat+ff)
        flag_acaod.append(flaod['flag_acaod'])
        flag_acaod_t.append(toutc(mat2py_time(flaod['t'])))
    else:
        flag_acaod.append([])
        flag_acaod_t.append([])


# In[17]:


s['fl_acaod'] = np.zeros_like(s['fl_QA'])
s['fl_acaod_noQA'] = np.zeros_like(s['fl_QA'])

for i,j in enumerate(np.unique(s['days'])):
    if len(flag_acaod[i]):
        iu = s['days']==j
        ii = np.where((flag_acaod_t[i]>=s['Start_UTC'][iu][0]) & (flag_acaod_t[i]<=s['Start_UTC'][iu][-1]))[0]
        fa = flag_acaod[i][ii]
        try:
            s['fl_acaod'][iu] = (fa[:,0]) & (s['fl_QA'][iu])
            s['fl_acaod_noQA'][iu] = (fa[:,0])
        except ValueError:
            ii = np.append(ii[0]-1,ii)
            fa = flag_acaod[i][ii]
            try:
                s['fl_acaod'][iu] = (fa[:,0]) & (s['fl_QA'][iu])
                s['fl_acaod_noQA'][iu] = (fa[:,0])
            except ValueError:
                fa = nearest_neighbor(flag_acaod_t[i],flag_acaod[i],s['Start_UTC'][iu],dist=1.0/3600.0)
                s['fl_acaod'][iu] = (fa[:,0]==1) & (s['fl_QA'][iu])
                s['fl_acaod_noQA'][iu] = (fa[:,0]==1)


# In[16]:


len(s['flag_acaod'])


# In[18]:


s['fl_acaod'] = s['flag_acaod']


# In[19]:


s['fl_acaod'] = (s['fl_acaod']==1) & (s['AOD0501']<4.0)


# In[20]:


len(s['fl_QA'])


# In[21]:


s['fl_below5'] = (s['fl_QA']) & (s['GPS_Alt']<5000.0)


# ## MODIS climatology

# In[20]:


fp


# From Rob Wood, email on 2016-11-08, 15:35, 
# These are just the all-years means (the black lines), for the month of September (DOY 244 to 273)

# In[21]:


m = sio.netcdf_file(fp+'data_other/climatology/mombl_oracles_routine_flight_NW-SE.nc')


# In[22]:


m.variables


# In[23]:


m2 = sio.netcdf_file(fp+'data_other/climatology/mombl_oracles_routine_flight_NW-SE_all.nc')


# In[24]:


m2.variables


# In[25]:


m2.variables['AOD_YRMEAN'].data.shape


# In[26]:


m2.variables['AOD_YRMEAN'].data


# ## MODIS AAC from Meyer for each day, subset by Meloe

# Collocation method:
# 1.     routine flight track is made of 100 points  (i.e. straight line from 23ºS, 13ºE to 10ºS, 0ºE)
# 2.     I select the closest MODIS AAC data point (within 15km) to each of the N=100 pt on the routine flight track.
# 3.     MODIS AAC AODs are then averaged in each of Sam's longitude bins
# 
# 
# Here’s a short description of the dataset (Kerry, please correct me if I’m wrong):
# Combination of MODIS and TERRA AAC AOD, high-resolution 0.1 x 0.1º. Suggested quality filtering is applied. Parameter is "Above_Cloud_AOT_MOD04Abs" i.e retrievals assuming the MOD04 absorbing aerosol model rather than the Haywood/SAFARI-2000 model, which is what is used to create the near-real time retrieval imagery in the field for ORACLES.

# In[27]:


import datetime
import load_utils as lu


# In[28]:


#aac = sio.loadmat(fp+'data_other/MODIS_AAC/MODIS_AAC_per_bin_all_flights_20160801_20161031.mat')
aac = sio.loadmat(fp+'data_other/MODIS_AAC/MODIS_AAC_per_bin_20180705.mat')


# In[29]:


aac['data_aac'][0,0].dtype.names


# In[30]:


daac = []


# In[31]:


for i in xrange(len(aac['data_aac'][0,:])):
    print aac['data_aac'][0,i]['FlightDay'][0][0][0]
    daac.append(lu.recarray_to_dict(aac['data_aac'][0,i]))


# In[32]:


# simplify the dict
for i in xrange(len(daac)):
    for k in daac[i].keys():
        daac[i][k] = daac[i][k][0,0]


# In[33]:


i_fltt,i_augl,i_sepl = [],[],[]
for i in xrange(len(daac)):
    daac[i]['datetime'] = datetime.datetime.strptime(daac[i]['FlightDay'][0], 'A%Y%j')
    daac[i]['date'] = daac[i]['datetime'].strftime('%Y%m%d')
    print i,daac[i]['date'], aac['data_aac'][0,i]['FlightDay'][0][0][0]
    if daac[i]['date'] in ['20160831','20160904','20160908','20160910','20160912','20160925']:
        i_fltt.append(i)
    if daac[i]['date'][0:6] in ['201608']:
        i_augl.append(i)
    if daac[i]['date'][0:6] in ['201609']:
        i_sepl.append(i)


# In[34]:


i_fltt, i_augl,i_sepl


# In[35]:


len(daac)


# In[36]:


i_aug = i_augl# range(0,31)
i_sep = i_sepl #range(31,61)
i_oct = range(61,92)
i_flt = i_fltt#[30,34,38,40,42,55]


# days in the binned for flights only
# 
#  - 0 20160831 A2016244
#  - 1 20160904 A2016248
#  - 2 20160908 A2016252
#  - 3 20160910 A2016254
#  - 4 20160912 A2016256
#  - 5 20160925 A2016269

# # Subset for routine flight

# In[22]:


d_rtn = ['20160831','20160904','20160908','20160910','20160912','20160925']


# In[23]:


d_rtnf = [20160831.0,20160904.0,20160908.0,20160910.0,20160912.0,20160925.0]


# In[24]:


d_irtn = [2.0,4.0,6.0,7.0,8.0,13.0]


# In[25]:


s['days'][0]


# In[26]:


ff


# In[27]:


ff = []
for d in d_irtn:
    ff.append(s['days']==d)


# In[28]:


fl = ((s['days']==d_rtnf[0]) | (s['days']==d_rtnf[1]) | (s['days']==d_rtnf[2]) | 
      (s['days']==d_rtnf[3]) | (s['days']==d_rtnf[4]) | (s['days']==d_rtnf[5]))


# In[29]:


s['fl_rtn'] = fl & s['fl']
s['fl_rtn']


# In[30]:


s['fl_rtna'] = fl & s['fl_acaod']
s['fl_rtna']


# In[31]:


np.nanmean(s['UNCAOD0501'][s['fl_acaod']]), np.nanmean(s['UNCAOD0501']), np.nanmean(s['UNCAOD0501'])/np.nanmean(s['AOD0501']) 


# In[143]:


np.nanmedian(s['UNCAOD0501'][s['fl_acaod']]), np.nanmedian(s['UNCAOD0501']), np.nanmedian(s['UNCAOD0501'])/np.nanmedian(s['AOD0501'])


# In[106]:


np.nanmean(s['UNCAOD0501'][s['fl_acaod']])/np.nanmean(s['AOD0501'][s['fl_acaod']])


# In[146]:


np.nanstd(s['UNCAOD0501'][s['fl_acaod']])


# In[107]:


np.nanmean(s['UNCAOD1020'][s['fl_acaod']]), np.nanmean(s['UNCAOD1020']), np.nanmean(s['UNCAOD1020'])/np.nanmean(s['AOD1020']) 


# In[145]:


np.nanmedian(s['UNCAOD1020'][s['fl_acaod']]),np.nanstd(s['UNCAOD1020'][s['fl_acaod']])


# # Now plot the routine flight aods

# ## all points

# In[168]:


plt.figure()
plt.plot(s['Longitude'][s['fl_rtn']],s['AOD0501'][s['fl_rtn']],'.',label='4STAR [600m - 1800m]',zorder=100)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.5,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.5,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlim(-2,16)
plt.xlabel('Longitude [$^\\circ$]')
plt.legend(numpoints=1,frameon=False)
plt.title('AOD above Clouds along routine flight path')
plt.savefig(fp+'plot/Climat_AAC_4STAR_all_MODIS.png',transparent=True,dpi=600)


# In[169]:


plt.figure()
plt.plot(s['Longitude'][s['fl_rtna']],s['AOD0501'][s['fl_rtna']],'.',label='4STAR [AAC]',zorder=100)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.5,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.5,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlim(-2,16)
plt.xlabel('Longitude [$^\\circ$]')
plt.legend(numpoints=1,frameon=False)
plt.title('AOD above Clouds along routine flight path')
plt.savefig(fp+'plot/Climat_AAC_flaged_4STAR_all_MODIS.png',transparent=True,dpi=600)


# ## heat map of points

# In[65]:


plt.figure()
cc = plt.hist2d(s['Longitude'][s['fl_rtn']],s['AOD0501'][s['fl_rtn']],bins=25,cmap=plt.cm.hot_r)
cb = plt.colorbar()
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.5,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.5,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.legend(numpoints=1,frameon=False)
plt.title('AOD above Clouds along routine flight path')
cb.set_label('4STAR AOD number between 600-1800m')
plt.savefig(fp+'plot\\Climat_AAC_4STAR_hist_MODIS.png',transparent=True,dpi=600)


# ## Use bins from heat map to make box whisker plot

# In[131]:


bins = []
pos = []
for i,c in enumerate(cc[1][0:-2]):
    lon_fl = (s['Longitude'][s['fl_rtn']]>=c)&(s['Longitude'][s['fl_rtn']]<cc[1][i+1])
    bins.append(s['AOD0501'][s['fl_rtn']][lon_fl])
    pos.append((c+cc[1][i+1])/2.0)


# In[107]:


def color_box(bp, color):

    # Define the elements to color. You can also add medians, fliers and means
    elements = ['boxes','caps','whiskers','medians','means','fliers']

    # Iterate over each of the elements changing the color
    for elem in elements:
        [plt.setp(bp[elem][idx], color=color) for idx in xrange(len(bp[elem]))]
    return


# In[151]:


plt.figure()
#plt.plot(s['Longitude'][s['fl_rtn']],s['AOD0501'][s['fl_rtn']],'.',label='4STAR [600m - 1800m]',
#         zorder=60,alpha=0.02,color='green')
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.5,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.3,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
bo = plt.boxplot(bins,0,'g.',showmeans=True,positions=pos)
color_box(bo,'green')
#[plt.setp(bo['fliers'][idx],alpha=0.5,)for idx in xrange(len(bo['fliers']))]
plt.plot(pos,[a.get_ydata()[0] for a in bo['means']],'s-',zorder=100,color='green',label='Mean binned 4STAR [600m-1800m]',lw=2)
plt.legend(numpoints=1,frameon=False)
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
plt.savefig(fp+'plot\\Climat_AAC_4STAR_box_MODIS.png',transparent=True,dpi=600)


# In[156]:


plt.figure()
#plt.plot(s['Longitude'][s['fl_rtn']],s['AOD0501'][s['fl_rtn']],'.',label='4STAR [600m - 1800m]',
#         zorder=60,alpha=0.02,color='green')
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
bo = plt.boxplot(bins,0,'g.',showmeans=True,positions=pos)
color_box(bo,'green')
[plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
plt.plot(pos,[a.get_ydata()[0] for a in bo['means']],'s-',zorder=100,color='green',label='Mean binned 4STAR [600m-1800m]',lw=2.5)
plt.legend(numpoints=1,frameon=False)
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
plt.savefig(fp+'plot\\Climat_AAC_4STAR_box_simpler_MODIS.png',transparent=True,dpi=600)


# ## Make another binning product, with less bins

# In[67]:


plt.figure()
cc2 = plt.hist2d(s['Longitude'][s['fl_rtn']],s['AOD0501'][s['fl_rtn']],bins=15,cmap=plt.cm.hot_r)
cb = plt.colorbar()


# In[38]:


bins2 = []
pos2 = []
for i,c in enumerate(cc2[1][0:-2]):
    lon_fl = (s['Longitude'][s['fl_rtn']]>=c)&(s['Longitude'][s['fl_rtn']]<cc2[1][i+1])
    bins2.append(s['AOD0501'][s['fl_rtn']][lon_fl])
    pos2.append((c+cc2[1][i+1])/2.0)


# In[162]:


plt.figure()
#plt.plot(s['Longitude'][s['fl_rtn']],s['AOD0501'][s['fl_rtn']],'.',label='4STAR [600m - 1800m]',
#         zorder=60,alpha=0.02,color='green')
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
bo = plt.boxplot(bins2,0,'g.',showmeans=True,positions=pos2)
color_box(bo,'green')
[plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
plt.plot(pos2,[a.get_ydata()[0] for a in bo['means']],'s-',zorder=100,color='green',label='Mean binned 4STAR [600m-1800m]',lw=2.5)
plt.legend(numpoints=1,frameon=False)
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
plt.savefig(fp+'plot\\Climat_AAC_4STAR_box2_simpler_MODIS.png',transparent=True,dpi=600)


# ## Now try with same bins as MODIS

# In[81]:


pos3 = m.variables['LONGITUDE'].data[0,:]


# In[82]:


pos3


# In[83]:


lims3 = pos3-0.5
lims3= np.append(lims3,pos3[-1]+0.5)


# In[84]:


lims3


# In[85]:


bins3 = []
for i,c in enumerate(lims3[0:-1]):
    lon_fl = (s['Longitude'][s['fl_rtn']]>=c)&(s['Longitude'][s['fl_rtn']]<lims3[i+1])
    bins3.append(s['AOD0501'][s['fl_rtn']][lon_fl])


# In[86]:


bins3a = []
for i,c in enumerate(lims3[0:-1]):
    lon_fl = (s['Longitude'][s['fl_rtna']]>=c)&(s['Longitude'][s['fl_rtna']]<lims3[i+1])
    bins3a.append(s['AOD0501'][s['fl_rtna']][lon_fl])


# In[87]:


len(lims3)


# In[88]:


len(bins3)


# In[89]:


len(pos3)


# In[54]:


plt.figure()
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
bo = plt.boxplot(bins3,0,'g.',showmeans=True,positions=pos3)
color_box(bo,'green')
[plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],'s-',zorder=100,color='green',label='Mean binned 4STAR [600m-1800m]',lw=2.5)
plt.legend(numpoints=1,frameon=False)
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
plt.savefig(fp+'plot/Climat_AAC_4STAR_box3_simpler_MODIS.png',transparent=True,dpi=600)


# In[180]:


plt.figure()
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
bo = plt.boxplot(bins3a,0,'g.',showmeans=True,positions=pos3)
color_box(bo,'green')
[plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],'s-',zorder=100,color='green',label='Mean binned 4STAR [AAC]',lw=2.5)
plt.legend(numpoints=1,frameon=False)
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
plt.savefig(fp+'plot/Climat_AAC_flagged_4STAR_box3_simpler_MODIS.png',transparent=True,dpi=600)


# In[187]:


plt.figure()
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AODFM_YRMEAN'].data[0,:,:],axis=1),
         '*-',color='r',zorder=50,lw=2.5,label='MODIS Fine mode median')
#plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
#         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)

plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
bo = plt.boxplot(bins3a,0,'g.',showmeans=True,positions=pos3)
color_box(bo,'green')
[plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
plt.plot(pos3,[a.get_ydata()[0] for a in bo['medians']],'x-',zorder=100,color='green',label='Median binned 4STAR [AAC]',lw=2.5)
plt.legend(numpoints=1,frameon=False)
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
plt.savefig(fp+'plot/Climat_AAC_flagged_4STAR_box3_simpler_MODIS_median.png',transparent=True,dpi=600)


# In[46]:


s['fl_alt12'] = (s['GPS_Alt']>=600)&(s['GPS_Alt']<1200)
s['fl_alt16'] = (s['GPS_Alt']>=500)&(s['GPS_Alt']<1600)


# In[47]:


s['fl_alt12']


# In[48]:


s['fl_rtn_12'] = s['fl_alt12'] & s['fl_QA'] & s['fl_rtn']
s['fl_rtn_16'] = s['fl_alt16'] & s['fl_QA'] & s['fl_rtn']


# In[49]:


s['fl_rtn_12']


# In[94]:


bins3_12 = []
for i,c in enumerate(lims3[0:-1]):
    lon_fl = (s['Longitude'][s['fl_rtn_12']]>=c)&(s['Longitude'][s['fl_rtn_12']]<lims3[i+1])
    bins3_12.append(s['AOD0501'][s['fl_rtn_12']][lon_fl])


# In[62]:


plt.figure()
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')

bo = plt.boxplot(bins3,0,'g.',showmeans=True,positions=pos3)
color_box(bo,'green')
[plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],'s-',zorder=100,color='green',label='Mean binned 4STAR [600m-1800m]',lw=2.5)

bo2 = plt.boxplot(bins3_12,0,'b.',showmeans=True,positions=pos3)
color_box(bo2,'blue')
[plt.setp(bo2['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
plt.plot(pos3,[a.get_ydata()[0] for a in bo2['means']],'s-',zorder=100,color='blue',label='Mean binned 4STAR [600m-1200m]',lw=2.5)

plt.legend(numpoints=1,frameon=False)
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
plt.savefig(fp+'plot\\Climat_AAC_4STAR_box3_simpler_2alt_MODIS.png',transparent=True,dpi=600)


# ## make a plot day by day

# In[32]:


d_irtn = [2.0,4.0,6.0,7.0,8.0,13.0]


# In[33]:


s['days']==2.0


# In[34]:


flr2 = s['fl_rtn_16']&(s['days']==d_rtnf[0])
flr4 = s['fl_rtn_16']&(s['days']==d_rtnf[1])
flr6 = s['fl_rtn_16']&(s['days']==d_rtnf[2])
flr7 = s['fl_rtn_16']&(s['days']==d_rtnf[3])
flr8 = s['fl_rtn_16']&(s['days']==d_rtnf[4])
flr13 = s['fl_rtn_16']&(s['days']==d_rtnf[5])


# In[35]:


fr2 = s['fl_rtn']&(s['days']==d_rtnf[0])
fr4 = s['fl_rtn']&(s['days']==d_rtnf[1])
fr6 = s['fl_rtn']&(s['days']==d_rtnf[2])
fr7 = s['fl_rtn']&(s['days']==d_rtnf[3])
fr8 = s['fl_rtn']&(s['days']==d_rtnf[4])
fr13 = s['fl_rtn']&(s['days']==d_rtnf[5])


# In[38]:


fr = [fr2,fr4,fr6,fr7,fr8,fr13]
flr = [flr2,flr4,flr6,flr7,flr8,flr13]


# In[37]:


flr2


# In[118]:


flr2a = s['fl_rtna']&(s['days']==d_rtnf[0])
flr4a = s['fl_rtna']&(s['days']==d_rtnf[1])
flr6a = s['fl_rtna']&(s['days']==d_rtnf[2])
flr7a = s['fl_rtna']&(s['days']==d_rtnf[3])
flr8a = s['fl_rtna']&(s['days']==d_rtnf[4])
flr13a = s['fl_rtna']&(s['days']==d_rtnf[5])


# In[119]:


fr2a = s['fl_rtna']&(s['days']==d_rtnf[0])
fr4a = s['fl_rtna']&(s['days']==d_rtnf[1])
fr6a = s['fl_rtna']&(s['days']==d_rtnf[2])
fr7a = s['fl_rtna']&(s['days']==d_rtnf[3])
fr8a = s['fl_rtna']&(s['days']==d_rtnf[4])
fr13a = s['fl_rtna']&(s['days']==d_rtnf[5])


# In[120]:


fra = [fr2a,fr4a,fr6a,fr7a,fr8a,fr13a]
flra = [flr2a,flr4a,flr6a,flr7a,flr8a,flr13a]


# In[39]:


cls = ['green','blue','yellow','cyan','magenta','orange']


# In[40]:


d_rtn


# In[51]:


plt.figure()
plt.plot(s['Longitude'][f],s['AOD0501'][f],'.')


# In[206]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')

for j,f in enumerate(fr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [0.6-1.8km]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),lw=2.5)

plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.4,1.0))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.savefig(fp+'plot/Climat_AAC_4STAR_box3_days_MODIS.png',transparent=True,dpi=600)


# In[209]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')

for j,f in enumerate(fra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [AAC]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),lw=2.5)

plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.4,1.0))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.savefig(fp+'plot/Climat_AAC_flagged_4STAR_box3_days_MODIS.png',transparent=True,dpi=600)


# In[211]:


plt.figure(figsize=(11,6))
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AODFM_YRMEAN'].data[0,:,:],axis=1),
         '*-',color='r',zorder=50,lw=2.5,label='MODIS Fine mode median')
#plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
#         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')

for j,f in enumerate(fra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['medians']],
             'x-',zorder=100,color=cls[j],label='{}/{} 4STAR [AAC]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),lw=2.5)

plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.4,1.0))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.savefig(fp+'plot/Climat_AAC_flagged_4STAR_box3_days_MODIS_median.png',transparent=True,dpi=600)


# ## Make an average from daily averages

# In[213]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
means = []

for j,f in enumerate(fr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [0.6-1.8km]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
plt.text(0.5, 0.5, 'Preliminary',
        verticalalignment='bottom', horizontalalignment='center',
        transform=plt.gca().transAxes,
        color='k', fontsize=18,zorder=1,alpha=0.3)
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',label='4STAR mean [0.6-1.8km]')
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.4,1.0))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.savefig(fp+'plot/Climat_AAC_4STAR_box3_days_avg_MODIS.png',transparent=True,dpi=600)


# In[215]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
means = []

for j,f in enumerate(fra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [AAC]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
#plt.text(0.5, 0.5, 'Preliminary',
#        verticalalignment='bottom', horizontalalignment='center',
#        transform=plt.gca().transAxes,
#        color='k', fontsize=18,zorder=1,alpha=0.3)
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',label='4STAR mean [AAC]')
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.4,1.0))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.savefig(fp+'plot_v2/Climat_AAC_flagged_4STAR_box3_days_avg_MODIS.png',transparent=True,dpi=600)


# In[217]:


plt.figure(figsize=(11,6))
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AODFM_YRMEAN'].data[0,:,:],axis=1),
         '*-',color='r',zorder=50,lw=2.5,label='MODIS Fine mode median')
#plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
#         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
means = []

for j,f in enumerate(fra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['medians'][idx],alpha=0.05)for idx in xrange(len(bo['medians']))]
    means.append([a.get_ydata()[0] for a in bo['medians']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['medians']],
             'x-',zorder=100,color=cls[j],label='{}/{} 4STAR [AAC]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
#plt.text(0.5, 0.5, 'Preliminary',
#        verticalalignment='bottom', horizontalalignment='center',
#        transform=plt.gca().transAxes,
#        color='k', fontsize=18,zorder=1,alpha=0.3)
plt.plot(pos3,np.nanmedian(np.array(means),axis=0),'s-k',label='4STAR median [AAC]')
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.4,1.0))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.savefig(fp+'plot_v2/Climat_AAC_flagged_4STAR_box3_days_avg_MODIS_median.png',transparent=True,dpi=600)


# In[218]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
means = []

for j,f in enumerate(flr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [0.6-1.2km]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',label='4STAR mean [0.6-1.2km]')
plt.text(0.5, 0.5, 'Preliminary',
        verticalalignment='bottom', horizontalalignment='center',
        transform=plt.gca().transAxes,
        color='k', fontsize=18,zorder=1,alpha=0.3)
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.4,1.0))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.savefig(fp+'plot_v2/Climat_AAC_4STAR_box3_days_avg_12_MODIS.png',transparent=True,dpi=600)


# ## Redo plots for entire column 

# In[50]:


s['fl_alt_6'] = (s['GPS_Alt']<=600)
s['fl_rtn_6'] = s['fl_alt_6'] & s['fl_QA'] & fl
flrr2 = s['fl_rtn_6']&(s['days']==2.0)
flrr4 = s['fl_rtn_6']&(s['days']==4.0)
flrr6 = s['fl_rtn_6']&(s['days']==6.0)
flrr7 = s['fl_rtn_6']&(s['days']==7.0)
flrr8 = s['fl_rtn_6']&(s['days']==8.0)
flrr13 = s['fl_rtn_6']&(s['days']==13.0)
flrr = [flrr2,flrr4,flrr6,flrr7,flrr8,flrr13]


# In[116]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above Clouds along routine flight path')
means = []

for j,f in enumerate(flrr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [0-0.6 km]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',label='4STAR mean [0-0.6 km]')
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.4,1.0))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.savefig(fp+'plot_v2/Climat_AAC_4STAR_box3_days_avg_surf_MODIS.png',transparent=True,dpi=600)


# ## Redo plots but with the MODIS AAC

# In[117]:


daac[0].keys()


# In[161]:


a['meanAODperbin']


# In[120]:


plt.figure()

for i,a in enumerate(daac):
   # print a['BinCenter'].shape,a['meanAODperbin'][0:-1].shape
   # print a['meanAODperbin']
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][0:-1],'x-',lw=1,color='orange',alpha=0.2)


# In[118]:


np.nanmean(ac,axis=0).shape


# In[101]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         '*-',color='b',label='MODIS Total AOD climatology',zorder=51,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above clouds along routine flight path')
means = []

for j,f in enumerate(flr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [0.6-1.2km]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color=cls[i],alpha=0.4)
    ac.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='darkgreen',label='MODIS AAC [Meyer] mean')
    
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,label='4STAR mean [0.6-1.2km]')
plt.text(0.5, 0.5, 'Preliminary',
        verticalalignment='bottom', horizontalalignment='center',
        transform=plt.gca().transAxes,
        color='k', fontsize=18,zorder=1,alpha=0.3)
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.45,1.05))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.78, box.height])

plt.savefig(fp+'plot_v2/MODIS_Climatology_vs_AAC_4STAR_and_Meyer_AAC_v2.png',transparent=True,dpi=600)


# In[119]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         '*-',color='b',label='MODIS Total AOD climatology',zorder=51,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above clouds along routine flight path')
means = []

for j,f in enumerate(flra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [AAC]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color=cls[i],alpha=0.4)
    ac.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='darkgreen',label='MODIS AAC [Meyer] mean')
    
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,label='4STAR mean [AAC]')
#plt.text(0.5, 0.5, 'Preliminary',
#        verticalalignment='bottom', horizontalalignment='center',
#        transform=plt.gca().transAxes,
#        color='k', fontsize=18,zorder=1,alpha=0.3)
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.45,1.05))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.78, box.height])

plt.savefig(fp+'plot_v3/MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC.png',transparent=True,dpi=600)


# In[71]:


plt.figure(figsize=(11,6))
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AODFM_YRMEAN'].data[0,:,:],axis=1),
         '*-',color='r',zorder=50,lw=2.5,label='MODIS Fine mode median')
#plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
#         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         '*-',color='b',label='MODIS Total AOD climatology',zorder=51,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above clouds along routine flight path')
means = []

for j,f in enumerate(flra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['medians'][idx],alpha=0.05)for idx in xrange(len(bo['medians']))]
    means.append([a.get_ydata()[0] for a in bo['medians']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['medians']],
             'x-',zorder=100,color=cls[j],label='{}/{} 4STAR [AAC]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color=cls[i],alpha=0.4)
    ac.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
plt.plot(a['BinCenter'][0,:],np.nanmedian(ac,axis=0),'^-',lw=3,color='darkgreen',label='MODIS AAC [Meyer] median')
    
plt.plot(pos3,np.nanmedian(np.array(means),axis=0),'s-k',lw=3.5,label='4STAR median [AAC]')
#plt.text(0.5, 0.5, 'Preliminary',
#        verticalalignment='bottom', horizontalalignment='center',
#        transform=plt.gca().transAxes,
#        color='k', fontsize=18,zorder=1,alpha=0.3)
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.45,1.05))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.78, box.height])

plt.savefig(fp+'plot_v2/MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_median.png',transparent=True,dpi=600)


# In[52]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         '*-',color='b',label='MODIS Total AOD climatology',zorder=51,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above clouds along routine flight path')
means = []

for j,f in enumerate(flr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [0.6-1.2km]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
ac = []
for i,a in enumerate(daac):
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color=cls[i],alpha=0.4)
    ac.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='darkgreen',label='MODIS AAC [Meyer] mean')
    
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,label='4STAR mean [0.6-1.2km]')
plt.text(0.5, 0.5, 'Preliminary',
        verticalalignment='bottom', horizontalalignment='center',
        transform=plt.gca().transAxes,
        color='k', fontsize=18,zorder=1,alpha=0.3)
#plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.45,1.05))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.78, box.height])

#plt.savefig(fp+'plot\\MODIS_Climatology_vs_AAC_4STAR_and_Meyer_AAC_noleg.png',transparent=True,dpi=600)


# In[72]:


plt.figure(figsize=(11,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine mode climatology',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         '*-',color='b',label='MODIS Total AOD climatology',zorder=51,lw=2.5)
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above clouds along routine flight path')
means = []

for j,f in enumerate(flra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR [AAC]'.format(d_rtn[j][4:6],d_rtn[j][6:8]),
             lw=2.5,alpha=0.2)    
    
ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color=cls[i],alpha=0.4)
    ac.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='darkgreen',label='MODIS AAC [Meyer] mean')
    
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,label='4STAR mean [AAC]')
#plt.text(0.5, 0.5, 'Preliminary',
#        verticalalignment='bottom', horizontalalignment='center',
#        transform=plt.gca().transAxes,
#        color='k', fontsize=18,zorder=1,alpha=0.3)
#plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.45,1.05))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.78, box.height])

plt.savefig(fp+'plot_v2/MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_nolegend.png',transparent=True,dpi=600)


# ## Plot out the previous plot comparison, but with MODIS AAC per monthly

# In[81]:


plt.figure(figsize=(12,6))
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine AOD (Sept. 2001-2013)',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         'v-',color='b',label='MODIS Total AOD (Sept. 2001-2013)',zorder=51,lw=0.5)
#plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,0],
#         'x-',color='grey',alpha=0.1,zorder=10,label='MODIS Yearly averages')
#plt.plot(m2.variables['LONGITUDE'].data[0,:],m2.variables['AODFM_YRMEAN'].data[0,:,:],'x-',color='grey',alpha=0.1,zorder=10)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)
plt.xlabel('Longitude [$^\\circ$]')
plt.title('AOD above clouds along routine flight path')
means = []

for j,f in enumerate(flr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's-',zorder=100,color=cls[j],label='{}/{} 4STAR {vv} [0.5-1.6km]'.format(d_rtn[j][4:6],d_rtn[j][6:8],vv=vv),
             lw=2.5,alpha=0.2)    
    
ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
#    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color=cls[i],alpha=0.4)
    ac.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',label='MODIS AAC for 2016 routine flights')

ac_aug = []
for i,a in enumerate([daac[j] for j in i_aug]):
    ac_aug.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MODIS AAC for Aug. 2016')

ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MODIS AAC for Sept. 2016')

ac_oct = []
for i,a in enumerate([daac[j] for j in i_oct]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'o-',lw=0.5,color='lightblue',label='MODIS AAC for Oct. 2016')


plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR {vv} mean [0.5-1.6km]'.format(vv=vv))
#plt.text(0.5, 0.6, 'Preliminary',
#        verticalalignment='bottom', horizontalalignment='center',
#        transform=plt.gca().transAxes,
#        color='k', fontsize=18,zorder=1,alpha=0.4)
plt.legend(numpoints=1,frameon=True,bbox_to_anchor=(1.68,1.05))
ti = plt.gca().set_xticks([0,2,4,6,8,10,12,14])
tl = plt.gca().set_xticklabels([0,2,4,6,8,10,12,14])
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.67, box.height])

plt.savefig(fp+'plot/MODIS_Climatology_vs_AAC_4STAR_{vv}_and_Meyer_AAC_and_monthlyavg.png'.format(vv=vv),transparent=True,dpi=600)


# In[163]:


fp


# # Prepare the split plot for publication

# In[108]:


plt.figure()
for j,f in enumerate(flr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    plt.plot(s['Longitude'][f],s['AOD0501'][f],'o',color=cls[j])
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],color=cls[j])for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],marker='.')for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's',zorder=100,color=cls[j],label='{}/{} 4STAR'.format(d_rtn[j][4:6],d_rtn[j][6:8],vv=vv),
             lw=2.5,alpha=0.2)    
    break


# In[232]:


plt.figure(figsize=(12,6))

ax = plt.subplot(2,1,1)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)

plt.title('AOD above clouds along 2016 routine flight paths')
means = []

for j,f in enumerate(flr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],color=cls[j])for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],marker='.')for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's',zorder=100,color=cls[j],label='{}/{} 4STAR'.format(d_rtn[j][4:6],d_rtn[j][6:8],vv=vv),
             lw=2.5,alpha=0.2)    
    
#plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR mean')

ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x-',lw=2,color=cls[i],alpha=0.4,
             label='{}/{} MODIS AAC'.format(d_rtn[i][4:6],d_rtn[i][6:8]))
    ac.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',label='MODIS AAC mean')

plt.legend(numpoints=1,frameon=True, ncol=2, bbox_to_anchor=(1.005,0.98),loc=2)



ax2 = plt.subplot(2,1,2,sharex=ax)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine AOD (Sept. 2001-2013)',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         'v-',color='b',label='MODIS Total AOD (Sept. 2001-2013)',zorder=51,lw=0.5)

ac_aug = []
for i,a in enumerate([daac[j] for j in i_aug]):
    ac_aug.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MODIS AAC mean for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MODIS AAC mean for all Sept. 2016')
#ac_oct = []
#for i,a in enumerate([daac[j] for j in i_oct]):
#    ac_sep.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'o-',lw=0.5,color='lightblue',label='MODIS AAC for Oct. 2016')

plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MODIS AAC mean for routine flights')

plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR mean [0.5-1.6km]')
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('AOD 500 nm')

plt.legend(numpoints=1,frameon=True, bbox_to_anchor=(1.005,0.98),loc=2)
ti = ax2.set_xticks([0,2,4,6,8,10,12,14])
tl = ax2.set_xticklabels([0,2,4,6,8,10,12,14])

ti = ax.set_xticks([0,2,4,6,8,10,12,14])
tl = ax.set_xticklabels([0,2,4,6,8,10,12,14])

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
box2 = ax2.get_position()
ax2.set_position([box2.x0, box2.y0, box2.width * 0.6, box2.height])



plt.savefig(fp+'plot_v2/ORACLES2016_MODIS_Climatology_vs_AAC_4STAR_and_Meyer_AAC_split.png',
            transparent=True,dpi=600)


# In[73]:


plt.figure(figsize=(12,6))

ax = plt.subplot(2,1,1)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)

plt.title('AOD above clouds along 2016 routine flight paths')
means = []

for j,f in enumerate(flra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],color=cls[j])for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],marker='.')for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's',zorder=100,color=cls[j],label='{}/{} 4STAR'.format(d_rtn[j][4:6],d_rtn[j][6:8],vv=vv),
             lw=2.5,alpha=0.2)    
    
#plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR mean')

ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x-',lw=2,color=cls[i],alpha=0.4,
             label='{}/{} MODIS AAC'.format(d_rtn[i][4:6],d_rtn[i][6:8]))
    ac.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',label='MODIS AAC mean')

plt.legend(numpoints=1,frameon=True, ncol=2, bbox_to_anchor=(1.005,0.98),loc=2)



ax2 = plt.subplot(2,1,2,sharex=ax)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine AOD (Sept. 2001-2013)',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         'v-',color='b',label='MODIS Total AOD (Sept. 2001-2013)',zorder=51,lw=0.5)

ac_aug = []
for i,a in enumerate([daac[j] for j in i_aug]):
    ac_aug.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MODIS AAC mean for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MODIS AAC mean for all Sept. 2016')
#ac_oct = []
#for i,a in enumerate([daac[j] for j in i_oct]):
#    ac_sep.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'o-',lw=0.5,color='lightblue',label='MODIS AAC for Oct. 2016')

plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MODIS AAC mean for routine flights')

plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR mean [AAC]')
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('AOD 500 nm')

plt.legend(numpoints=1,frameon=True, bbox_to_anchor=(1.005,0.98),loc=2)
ti = ax2.set_xticks([0,2,4,6,8,10,12,14])
tl = ax2.set_xticklabels([0,2,4,6,8,10,12,14])

ti = ax.set_xticks([0,2,4,6,8,10,12,14])
tl = ax.set_xticklabels([0,2,4,6,8,10,12,14])

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
box2 = ax2.get_position()
ax2.set_position([box2.x0, box2.y0, box2.width * 0.6, box2.height])



plt.savefig(fp+'plot_v2/ORACLES2016_MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_split.png',
            transparent=True,dpi=600)


# In[74]:


plt.figure(figsize=(12,6))

ax = plt.subplot(2,1,1)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)

plt.title('AOD above clouds along 2016 routine flight paths')
means = []

for j,f in enumerate(flra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],color=cls[j])for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],marker='.')for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['medians']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's',zorder=100,color=cls[j],label='{}/{} 4STAR'.format(d_rtn[j][4:6],d_rtn[j][6:8],vv=vv),
             lw=2.5,alpha=0.2)    
    
#plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR mean')

ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x-',lw=2,color=cls[i],alpha=0.4,
             label='{}/{} MODIS AAC'.format(d_rtn[i][4:6],d_rtn[i][6:8]))
    ac.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',label='MODIS AAC mean')

plt.legend(numpoints=1,frameon=True, ncol=2, bbox_to_anchor=(1.005,0.98),loc=2)



ax2 = plt.subplot(2,1,2,sharex=ax)
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AODFM_YRMEAN'].data[0,:,:],axis=1),
         '*-',color='r',zorder=50,lw=2.5,label='MODIS Fine AOD (median Sept. 2001-2013)')
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AOD_YRMEAN'].data[0,:,:],axis=1),
         'v-',color='b',zorder=51,lw=0.5,label='MODIS Total AOD (median Sept. 2001-2013)')

#plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
#         '*-',color='r',label='MODIS Fine AOD (Sept. 2001-2013)',zorder=50,lw=2.5)
#plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
#         'v-',color='b',label='MODIS Total AOD (Sept. 2001-2013)',zorder=51,lw=0.5)

ac_aug = []
for i,a in enumerate([daac[j] for j in i_aug]):
    ac_aug.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmedian(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MODIS AAC median for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmedian(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MODIS AAC median for all Sept. 2016')
#ac_oct = []
#for i,a in enumerate([daac[j] for j in i_oct]):
#    ac_sep.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'o-',lw=0.5,color='lightblue',label='MODIS AAC for Oct. 2016')

plt.plot(a['BinCenter'][0,:],np.nanmedian(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MODIS AAC median for routine flights')

plt.plot(pos3,np.nanmedian(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR median [AAC]')
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('AOD 500 nm')

plt.legend(numpoints=1,frameon=True, bbox_to_anchor=(1.005,0.98),loc=2)
ti = ax2.set_xticks([0,2,4,6,8,10,12,14])
tl = ax2.set_xticklabels([0,2,4,6,8,10,12,14])

ti = ax.set_xticks([0,2,4,6,8,10,12,14])
tl = ax.set_xticklabels([0,2,4,6,8,10,12,14])

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.55, box.height])
box2 = ax2.get_position()
ax2.set_position([box2.x0, box2.y0, box2.width * 0.55, box2.height])



plt.savefig(fp+'plot_v2/ORACLES2016_MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_split_median.png',
            transparent=True,dpi=600)


# In[121]:


plt.figure(figsize=(12,8))

ax = plt.subplot(3,1,1)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)

plt.title('AOD above clouds along 2016 routine flight paths')

means = []
medians = []
for j,f in enumerate(flra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],color=cls[j])for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],marker='.')for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    medians.append([a.get_ydata()[0] for a in bo['medians']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's',zorder=100,color=cls[j],label='{}/{} 4STAR'.format(d_rtn[j][4:6],d_rtn[j][6:8],vv=vv),
             lw=2.5,alpha=0.2)   
    
meansr = []
mediansr = []
for j,f in enumerate(flr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    #plt.plot(s['Longitude'][f],s['AOD0501'][f],'.',color=cls[j],alpha=0.02)
    bn = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)

    meansr.append([a.get_ydata()[0] for a in bn['means']])
    mediansr.append([a.get_ydata()[0] for a in bn['medians']])
    [plt.setp(bn['fliers'][idx],alpha=0.0)for idx in xrange(len(bn['fliers']))]
    [plt.setp(bn['means'][idx],alpha=0.0)for idx in xrange(len(bn['means']))]
    [plt.setp(bn['boxes'][idx],alpha=0.0)for idx in xrange(len(bn['boxes']))]
    [plt.setp(bn['medians'][idx],alpha=0.0)for idx in xrange(len(bn['medians']))]
    [plt.setp(bn['whiskers'][idx],alpha=0.0)for idx in xrange(len(bn['whiskers']))]
    [plt.setp(bn['caps'][idx],alpha=0.0)for idx in xrange(len(bn['caps']))]
    
#plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR mean')

ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    
    plt.plot(a['BinCenter_longitude'][0,:],a['meanAODperbin'][1:,0],'x-',lw=2,color=cls[i],alpha=0.4,
             label='{}/{} MODIS AAC'.format(d_rtn[i][4:6],d_rtn[i][6:8]))
    ac.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x--',lw=1,color='k',alpha=0.4,label='MODIS AAC [Meyer] daily avg.')
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',label='MODIS AAC mean')

plt.legend(numpoints=1,frameon=True, ncol=2, bbox_to_anchor=(1.005,1.0),loc=2,handletextpad=0.2,columnspacing=0.7)



ax2 = plt.subplot(3,1,2,sharex=ax)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine AOD (Sept. 2001-2013)',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         'v-',color='b',label='MODIS Total AOD (Sept. 2001-2013)',zorder=51,lw=0.5)

ac_aug = []
for i,a in enumerate([daac[j] for j in i_aug]):
    ac_aug.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter_longitude'][0,:],np.nanmean(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MODIS AAC for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter_longitude'][0,:],np.nanmean(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MODIS AAC for all Sept. 2016')
#ac_oct = []
#for i,a in enumerate([daac[j] for j in i_oct]):
#    ac_sep.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'o-',lw=0.5,color='lightblue',label='MODIS AAC for Oct. 2016')

plt.plot(a['BinCenter_longitude'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MODIS AAC for routine flights')

plt.plot(pos3,np.nanmean(np.array(meansr),axis=0),'x-',color='grey',lw=0.5,zorder=180,label='4STAR [0.5-1.6 km]')
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR AAC')
#plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('Mean AOD 500 nm')
#plt.title('Mean')

#plt.legend(numpoints=1,frameon=True, bbox_to_anchor=(1.005,0.98),loc=2)


ax3 = plt.subplot(3,1,3,sharex=ax)
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AODFM_YRMEAN'].data[0,:,:],axis=1),
         '*-',color='r',zorder=50,lw=2.5,label='MODIS Fine AOD (Sept. 2001-2013)')
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AOD_YRMEAN'].data[0,:,:],axis=1),
         'v-',color='b',zorder=51,lw=0.5,label='MODIS Total AOD (Sept. 2001-2013)')

plt.ylim(0,0.8)
#plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
#         '*-',color='r',label='MODIS Fine AOD (Sept. 2001-2013)',zorder=50,lw=2.5)
#plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
#         'v-',color='b',label='MODIS Total AOD (Sept. 2001-2013)',zorder=51,lw=0.5)

ac_aug = []
for i,a in enumerate([daac[j] for j in i_aug]):
    ac_aug.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter_longitude'][0,:],np.nanmedian(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MODIS AAC for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter_longitude'][0,:],np.nanmedian(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MODIS AAC for all Sept. 2016')
#ac_oct = []
#for i,a in enumerate([daac[j] for j in i_oct]):
#    ac_sep.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'o-',lw=0.5,color='lightblue',label='MODIS AAC for Oct. 2016')

plt.plot(a['BinCenter_longitude'][0,:],np.nanmedian(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MODIS AAC for routine flights')

plt.plot(pos3,np.nanmedian(np.array(mediansr),axis=0),'x-',color='grey',lw=0.5,zorder=180,label='4STAR [0.5-1.6 km]')
plt.plot(pos3,np.nanmedian(np.array(medians),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR AAC')
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('Median AOD 500 nm')
#plt.title('Median')

plt.legend(numpoints=1,frameon=True, bbox_to_anchor=(1.005,1.75),loc=2)
ti = ax3.set_xticks([0,2,4,6,8,10,12,14])
tl = ax3.set_xticklabels([0,2,4,6,8,10,12,14])

ti = ax2.set_xticks([0,2,4,6,8,10,12,14])
tl = ax2.set_xticklabels([0,2,4,6,8,10,12,14])

ti = ax.set_xticks([0,2,4,6,8,10,12,14])
tl = ax.set_xticklabels([0,2,4,6,8,10,12,14])

box = ax.get_position()
ax.set_position([box.x0-0.06, box.y0, box.width * 0.7, box.height])
box2 = ax2.get_position()
ax2.set_position([box2.x0-0.06, box2.y0, box2.width * 0.7, box2.height])
box3 = ax3.get_position()
ax3.set_position([box3.x0-0.06, box3.y0, box3.width * 0.7, box3.height])


plt.savefig(fp+'plot_v3/ORACLES2016_MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_3_split_v2.png',
            transparent=True,dpi=500)


# In[122]:


maa,mad = np.nanmean(np.nanmean(np.array(means),axis=0)), np.nanmedian(np.nanmedian(np.array(medians),axis=0))


# In[123]:


np.nanmean(m.variables['AODFM_CLIMOMEAN'].data[0,:]), np.nanmean(np.nanmean(ac_aug,axis=0)),np.nanmean(np.nanmean(ac_sep,axis=0)),np.nanmean(np.nanmean(np.array(means),axis=0))


# In[124]:


np.median(np.median(m2.variables['AODFM_YRMEAN'].data[0,:,:],axis=1)), np.nanmedian(np.nanmedian(ac_aug,axis=0)),np.nanmedian(np.nanmedian(ac_sep,axis=0)), np.nanmedian(np.nanmedian(np.array(medians),axis=0))


# In[125]:


np.nanmean(m.variables['AODFM_CLIMOMEAN'].data[0,:])/maa*100.0, np.nanmean(np.nanmean(ac_aug,axis=0))/maa*100.0,np.nanmean(np.nanmean(ac_sep,axis=0))/maa*100.0


# In[10]:


1.0-0.30186/0.35921


# In[11]:


0.35921-0.30186


# # Plot Histograms of measured AOD and angstrom

# ## Plot histograms of ACAOD sampled at various altitudes

# In[32]:


plt.figure()
plt.hist(s['AOD0501'][s['fl_QA']],bins=30,range=(0,1.2),edgecolor='None',label='All AOD',alpha=0.3,normed=True)
plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.xlabel('AOD 501 nm')
plt.ylabel('Normalized frequency counts')
plt.title('ORACLES 2016 4STAR AOD')
plt.legend(frameon=False)


# In[22]:


s['fl6'] = s['fl_alt_6'] & s['fl_QA'] 
s['fl6b']  = s['fl_alt_6'] & s['fl_QA'] & (s['AOD0501']<4.0)


# In[23]:


sum(s['fl6']),sum(s['fl6b'])


# In[26]:


sum(s['fl_acaod'])


# In[59]:


np.nanmin(s['AOD0501'][s['fl6b']]),np.nanmax(s['AOD0501'][s['fl6b']])


# In[84]:


plt.figure(figsize=(7,3))
plt.hist([s['AOD0501'][s['fl_QA']],s['AOD0501'][s['fl_acaod']],s['AOD0501'][s['fl6']]], 
         bins=15,range=(0,1.0),edgecolor='None',label=['All AOD','AOD above clouds','AOD below 0.6 km'],alpha=0.6,normed=True)
#plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
#plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_QA']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_QA']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')


plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_QA']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_QA']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2,ls='--')

plt.xlabel('AOD 501 nm')
plt.ylabel('Normalized frequency counts')
plt.title('ORACLES 2016 4STAR AOD')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,4,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False)
#plt.legend(frameon=False)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_AOD501_histogram.png',
            transparent=True,dpi=500)


# ### Plot histogram per day

# In[43]:


np.unique(s['days'])


# In[44]:


iid = s['days']==s['days'][0]


# In[45]:


s['AOD0501'][iid][s['fl6'][iid]]


# In[28]:


for d in np.unique(s['days']):
    
    iid = s['days']==d
    fig = plt.figure(figsize=(7,6))
    plt.subplot(2,1,1)
    try:
        plt.hist([s['AOD0501'][iid][s['fl6'][iid]],s['AOD0501'][iid][s['fl_acaod'][iid]]], bins=15,range=(0,1.0),color=['lightcoral','b'],histtype='bar',
             edgecolor='None',label=['Full Column','Only above clouds'],alpha=0.75,normed=False,stacked=True)
    except:
        continue

    plt.axvline(x=np.nanmean(s['AOD0501'][iid][s['fl_both'][iid]]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
    plt.axvline(x=np.nanmedian(s['AOD0501'][iid][s['fl_both'][iid]]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')


    plt.axvline(x=np.nanmean(s['AOD0501'][iid][s['fl_both'][iid]]),color='k',ymin=0, ymax=10,lw=2)
    plt.axvline(x=np.nanmean(s['AOD0501'][iid][s['fl_acaod'][iid]]),color='b',ymin=0, ymax=10,lw=2)
    plt.axvline(x=np.nanmean(s['AOD0501'][iid][s['fl6'][iid]]),color='lightcoral',ymin=0, ymax=10,lw=2)

    plt.axvline(x=np.nanmedian(s['AOD0501'][iid][s['fl_both'][iid]]),color='k',ymin=0, ymax=10,lw=2,ls='--')
    plt.axvline(x=np.nanmedian(s['AOD0501'][iid][s['fl_acaod'][iid]]),color='b',ymin=0, ymax=10,lw=2,ls='--')
    plt.axvline(x=np.nanmedian(s['AOD0501'][iid][s['fl6'][iid]]),color='lightcoral',ymin=0, ymax=10,lw=2,ls='--')

    plt.xlabel('AOD 501 nm')
    plt.ylabel('Samples')
    plt.title('ORACLES 4STAR AOD for {:8.0f}'.format(d))
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [2,3,0,1]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False)
    #plt.legend(frameon=False)

    fig.subplots_adjust(hspace=.3)
    plt.subplot(2,1,2)
    plt.hist([s['AOD1020'][iid][s['fl6'][iid]],s['AOD1020'][iid][s['fl_acaod'][iid]]],color=['lightcoral','b'],histtype='bar',
            bins=15,range=[0.0,0.3],label=['Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=False,stacked=True)
    plt.xlim(0.0,0.3)
    plt.ylabel('Samples')
    plt.xlabel('AOD 1020 nm')

    plt.axvline(x=np.nanmean(s['AOD1020'][iid][s['fl_both'][iid]]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
    plt.axvline(x=np.nanmedian(s['AOD1020'][iid][s['fl_both'][iid]]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

    plt.axvline(x=np.nanmean(s['AOD1020'][iid][s['fl_both'][iid]]),color='k',ymin=0, ymax=10,lw=2)
    plt.axvline(x=np.nanmean(s['AOD1020'][iid][s['fl_acaod'][iid]]),color='b',ymin=0, ymax=10,lw=2)
    plt.axvline(x=np.nanmean(s['AOD1020'][iid][s['fl6'][iid]]),color='lightcoral',ymin=0, ymax=10,lw=2)

    plt.axvline(x=np.nanmedian(s['AOD1020'][iid][s['fl_both'][iid]]),color='k',ymin=0, ymax=10,lw=2,ls='--')
    plt.axvline(x=np.nanmedian(s['AOD1020'][iid][s['fl_acaod'][iid]]),color='b',ymin=0, ymax=10,lw=2,ls='--')
    plt.axvline(x=np.nanmedian(s['AOD1020'][iid][s['fl6'][iid]]),color='lightcoral',ymin=0, ymax=10,lw=2,ls='--')

    plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_AOD_2wvl_histogram_stack_{:8.0f}.png'.format(d),
                transparent=True,dpi=500)


# ## Now get the angstrom

# ### Using the polynomial fit method

# In[62]:


nn = s.keys()
nn.sort()


# In[63]:


wvl = np.array([float(i[-4:]) for i in nn[0:24]])


# In[64]:


nn[0:24]


# In[65]:


aods = np.array([s[i] for i in nn[0:24]]).T


# In[66]:


aods.shape


# In[67]:


uncaods = np.array([s[i] for i in nn[28:28+24]]).T


# In[113]:


s['polyaod'] = []
s['polylogaod'] = []


# In[114]:


for i,u in enumerate(s['Start_UTC']):
    pp = su.aod_polyfit(wvl,aods[i,:],polynum=3)
    s['polyaod'].append(pp)
    pl = su.logaod_polyfit(wvl,aods[i,:],polynum=3)
    s['polylogaod'].append(pl)
s['polyaod'] = np.array(s['polyaod'])
s['polylogaod'] = np.array(s['polylogaod'])


# In[115]:


s['angs'] = su.angstrom_from_logpoly(s['polylogaod'],[380.0,470.0,500.0,530.0,660.0,865.0,1250.0],polynum=3)


# In[116]:


s['angs'].shape


# In[117]:


awvl = [380.0,470.0,500.0,530.0,660.0,865.0,1250.0]


# In[62]:


plt.figure()
plt.hist([s['angs'][s['fl_acaod'],0],s['angs'][s['fl_acaod'],1],s['angs'][s['fl_acaod'],2]],
        bins=15,range=[-0.5,3.0],label=['380','470','500'],edgecolor='None',alpha=0.75,normed=True)
plt.legend()


# ### Using the typical 2 wavelength method

# In[68]:


wvl


# In[69]:


ja,je = 3,15
wvl[ja],wvl[je]


# In[70]:


s['angs_470_865'] = []
for i,u in enumerate(s['Start_UTC']):
    to = np.log(aods[i,je])-np.log(aods[i,ja])
    bo = np.log(wvl[je])-np.log(wvl[ja])
    s['angs_470_865'].append(to/bo*-1.0)
s['angs_470_865'] = np.array(s['angs_470_865'])


# ## Plot the histogram of angstrom exponent

# In[376]:


plt.figure(figsize=(7,3))
plt.hist([s['angs'][s['fl_QA'],2],s['angs'][s['fl_acaod'],2],s['angs'][s['fl6'],2]],
        bins=15,range=[-0.5,3.0],label=['All data','Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=True)
plt.ylabel('Normalized frequency counts')
plt.xlabel('Angstrom at 500 nm')
plt.title('ORACLES 4STAR column Angstrom exponent')

plt.axvline(x=np.nanmean(s['angs'][s['fl_QA'],2]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_QA'],2]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['angs'][s['fl_QA'],2]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl_acaod'],2]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl6'],2]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['angs'][s['fl_QA'],2]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_acaod'],2]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl6'],2]),color='r',ymin=0, ymax=10,lw=2,ls='--')

handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,4,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False,loc=0)
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_angstrom_histogram.png',
            transparent=True,dpi=500)


# ## Combine both histogram

# In[66]:


ia = 4


# In[104]:


fig = plt.figure(figsize=(7,6))
plt.subplot(2,1,1)
plt.hist([s['AOD0501'][s['fl_QA']],s['AOD0501'][s['fl_acaod']],s['AOD0501'][s['fl6']]], 
         bins=15,range=(0,1.0),edgecolor='None',label=['All data','Only above clouds','Below 0.6 km'],alpha=0.75,normed=True)
#plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
#plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_QA']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_QA']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')


plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_QA']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_QA']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2,ls='--')

plt.xlabel('AOD 501 nm')
plt.ylabel('Normalized frequency counts')
#plt.title('ORACLES 2016 4STAR AOD')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,4,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False)
#plt.legend(frameon=False)

fig.subplots_adjust(hspace=.3)
plt.subplot(2,1,2)
plt.hist([s['angs'][s['fl_QA'],ia],s['angs'][s['fl_acaod'],ia],s['angs'][s['fl6'],ia]],
        bins=15,range=[-0.2,2.5],label=['All data','Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=True)
plt.xlim(-0.2,2.5)
plt.ylabel('Normalized frequency counts')
plt.xlabel('Angstrom at {:4.0f} nm'.format(awvl[ia]))
#plt.title('ORACLES 4STAR column Angstrom exponent')

plt.axvline(x=np.nanmean(s['angs'][s['fl_QA'],ia]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_QA'],ia]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['angs'][s['fl_QA'],ia]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl_acaod'],ia]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl6'],ia]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['angs'][s['fl_QA'],ia]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_acaod'],ia]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl6'],ia]),color='r',ymin=0, ymax=10,lw=2,ls='--')

#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,4,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False,loc=0)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_AOD_angstrom_histogram.png',
            transparent=True,dpi=500)


# In[100]:


fig = plt.figure(figsize=(7,6))
plt.subplot(2,1,1)
plt.hist([s['AOD0501'][s['fl_below5']],s['AOD0501'][s['fl_acaod']],s['AOD0501'][s['fl6']]], 
         bins=15,range=(0,1.0),edgecolor='None',label=['All data below 5km','Only above clouds','Below 0.6 km'],alpha=0.75,normed=True)
#plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
#plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_below5']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_below5']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')


plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_below5']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_below5']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2,ls='--')

plt.xlabel('AOD 501 nm')
plt.ylabel('Normalized frequency counts')
#plt.title('ORACLES 2016 4STAR AOD')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,4,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False)
#plt.legend(frameon=False)

fig.subplots_adjust(hspace=.3)
plt.subplot(2,1,2)
plt.hist([s['angs'][s['fl_below5'],ia],s['angs'][s['fl_acaod'],ia],s['angs'][s['fl6'],ia]],
        bins=15,range=[-0.2,2.5],label=['All data below 5km','Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=True)
plt.xlim(-0.2,2.5)
plt.ylabel('Normalized frequency counts')
plt.xlabel('Angstrom at {:4.0f} nm'.format(awvl[ia]))
#plt.title('ORACLES 4STAR column Angstrom exponent')

plt.axvline(x=np.nanmean(s['angs'][s['fl_below5'],ia]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_below5'],ia]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['angs'][s['fl_below5'],ia]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl_acaod'],ia]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl6'],ia]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['angs'][s['fl_below5'],ia]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_acaod'],ia]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl6'],ia]),color='r',ymin=0, ymax=10,lw=2,ls='--')

#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,4,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False,loc=0)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_AOD_angstrom_histogram_sub.png',
            transparent=True,dpi=500)


# In[139]:


s['fl_both'] = s['fl_acaod'] | s['fl6']
s['fl_bothb'] = s['fl_acaod'] | s['fl6b']


# In[137]:


fig = plt.figure(figsize=(7,6))
plt.subplot(2,1,1)
plt.hist([s['AOD0501'][s['fl_both']],s['AOD0501'][s['fl_acaod']],s['AOD0501'][s['fl6']]], 
         bins=15,range=(0,1.0),edgecolor='None',label=['Combined','Only above clouds','Below 0.6 km'],alpha=0.75,normed=True)
#plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
#plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')


plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_both']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_both']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2,ls='--')

plt.xlabel('AOD 501 nm')
plt.ylabel('Normalized frequency counts')
#plt.title('ORACLES 2016 4STAR AOD')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,4,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False)
#plt.legend(frameon=False)

fig.subplots_adjust(hspace=.3)
plt.subplot(2,1,2)
plt.hist([s['angs'][s['fl_both'],ia],s['angs'][s['fl_acaod'],ia],s['angs'][s['fl6'],ia]],
        bins=15,range=[-0.2,2.5],label=['Combined','Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=True)
plt.xlim(-0.2,2.5)
plt.ylabel('Normalized frequency counts')
plt.xlabel('Angstrom at {:4.0f} nm'.format(awvl[ia]))
#plt.title('ORACLES 4STAR column Angstrom exponent')

plt.axvline(x=np.nanmean(s['angs'][s['fl_both'],ia]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_both'],ia]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['angs'][s['fl_both'],ia]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl_acaod'],ia]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl6'],ia]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['angs'][s['fl_both'],ia]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_acaod'],ia]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl6'],ia]),color='r',ymin=0, ymax=10,lw=2,ls='--')

#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,4,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False,loc=0)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_AOD_angstrom_histogram_comb.png',
            transparent=True,dpi=500)


# Print out the mean, median, std, and number of samples for 501 nm the - acaod, below 0.6 km, and combine below cloud and acaod

# In[125]:


np.nanmean(s['AOD0501'][s['fl_acaod']]), np.nanmedian(s['AOD0501'][s['fl_acaod']]), np.nanstd(s['AOD0501'][s['fl_acaod']]),len(s['AOD0501'][s['fl_acaod']])


# In[71]:


s['fl6b'] = s['fl6'] & (s['AOD1020']<1.5)


# In[138]:


np.nanmean(s['AOD0501'][s['fl6']]), np.nanmedian(s['AOD0501'][s['fl6']]), np.nanstd(s['AOD0501'][s['fl6b']]),len(s['AOD0501'][s['fl6b']])


# In[141]:


np.nanmean(s['AOD0501'][s['fl_both']]), np.nanmedian(s['AOD0501'][s['fl_both']]), np.nanstd(s['AOD0501'][s['fl_bothb']]),len(s['AOD0501'][s['fl_both']])


# In[96]:


plt.figure()
plt.plot(s['AOD0501'][s['fl_acaod']])
plt.plot(s['AOD0501'][s['fl_both']])
plt.plot(s['AOD0501'][s['fl6']])



# In[80]:


plt.figure()
plt.plot(np.where(s['fl_both'])[0])
plt.plot(np.where(s['fl_acaod'])[0])
plt.plot(np.where(s['fl6'])[0])


# In[71]:


sum(s['fl_QA'])


# In[65]:


(sum(s['fl_acaod'])+sum(s['fl6']))/3600.0


# In[66]:


sum(s['fl_acaod'])/3600.0


# ### Combine histogram of AOD at 2 wavelengths

# In[174]:


fig = plt.figure(figsize=(7,6))
plt.subplot(2,1,1)
plt.hist([s['AOD0501'][s['fl_both']],s['AOD0501'][s['fl_acaod']],s['AOD0501'][s['fl6']]], 
         bins=15,range=(0,1.0),edgecolor='None',label=['Combined','Only above clouds','Below 0.6 km'],alpha=0.75,normed=True)
#plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
#plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')


plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_both']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_both']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2,ls='--')

plt.xlabel('AOD 501 nm')
plt.ylabel('Normalized frequency counts')
#plt.title('ORACLES 2016 4STAR AOD')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,4,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False)
#plt.legend(frameon=False)

fig.subplots_adjust(hspace=.3)
plt.subplot(2,1,2)
plt.hist([s['AOD1020'][s['fl_both']],s['AOD1020'][s['fl_acaod']],s['AOD1020'][s['fl6']]],
        bins=15,range=[0.0,0.3],label=['Combined','Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=True)
plt.xlim(0.0,0.3)
plt.ylabel('Normalized frequency counts')
plt.xlabel('AOD 1020 nm')
#plt.title('ORACLES 4STAR column Angstrom exponent')

plt.axvline(x=np.nanmean(s['AOD1020'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['AOD1020'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['AOD1020'][s['fl_both']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD1020'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD1020'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['AOD1020'][s['fl_both']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD1020'][s['fl_acaod']]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD1020'][s['fl6']]),color='r',ymin=0, ymax=10,lw=2,ls='--')

#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,4,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False,loc=0)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_AOD_2wvl_histogram_comb.png',
            transparent=True,dpi=500)


# Print out the mean, median, std, and number of samples for 1020 nm the - acaod, below 0.6 km, and combine below cloud and acaod

# In[126]:


np.nanmean(s['AOD1020'][s['fl_acaod']]), np.nanmedian(s['AOD1020'][s['fl_acaod']]), np.nanstd(s['AOD1020'][s['fl_acaod']]),len(s['AOD1020'][s['fl_acaod']])


# In[132]:


s['fl6b'] = s['fl6'] & (s['AOD1020']<1.5)


# In[135]:


np.nanmean(s['AOD1020'][s['fl6b']]), np.nanmedian(s['AOD1020'][s['fl6b']]), np.nanstd(s['AOD1020'][s['fl6b']]),len(s['AOD1020'][s['fl6b']])


# In[128]:


np.nanmean(s['AOD1020'][s['fl_both']]), np.nanmedian(s['AOD1020'][s['fl_both']]), np.nanstd(s['AOD1020'][s['fl_both']]),len(s['AOD1020'][s['fl_both']])


# ### Combined histogram at 2 wvl, in stacked bars

# In[193]:


fig = plt.figure(figsize=(7,6))
plt.subplot(2,1,1)
#plt.hist([s['AOD0501'][s['fl_both']],s['AOD0501'][s['fl_acaod']],s['AOD0501'][s['fl6']]], 
#         bins=15,range=(0,1.0),edgecolor='None',label=['Combined','Only above clouds','Below 0.6 km'],alpha=0.75,normed=True)
#plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
#plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)
plt.hist([s['AOD0501'][s['fl6']],s['AOD0501'][s['fl_acaod']]], bins=15,range=(0,1.0),color=['lightcoral','b'],histtype='bar',
         edgecolor='None',label=['Full Column','Only above clouds'],alpha=0.75,normed=False,stacked=True)


plt.axvline(x=np.nanmean(s['AOD0501'][s['fl6']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl6']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')


#plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl_acaod']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD0501'][s['fl6']]),color='lightcoral',ymin=0, ymax=10,lw=2)

#plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl_acaod']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD0501'][s['fl6']]),color='lightcoral',ymin=0, ymax=10,lw=2,ls='--')

plt.xlabel('AOD 501 nm')
plt.ylabel('Samples')
pu.sub_note('a)')
#plt.title('ORACLES 2016 4STAR AOD')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False)
#plt.legend(frameon=False)

fig.subplots_adjust(hspace=.3)
plt.subplot(2,1,2)
#plt.hist([s['AOD1020'][s['fl_both']],s['AOD1020'][s['fl_acaod']],s['AOD1020'][s['fl6']]],
#        bins=15,range=[0.0,0.3],label=['Combined','Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=True)
plt.hist([s['AOD1020'][s['fl6']],s['AOD1020'][s['fl_acaod']]],color=['lightcoral','b'],histtype='bar',
        bins=15,range=[0.0,0.3],label=['Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=False,stacked=True)
plt.xlim(0.0,0.3)
plt.ylabel('Samples')
plt.xlabel('AOD 1020 nm')
#plt.title('ORACLES 4STAR column Angstrom exponent')

plt.axvline(x=np.nanmean(s['AOD1020'][s['fl6']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['AOD1020'][s['fl6']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

#plt.axvline(x=np.nanmean(s['AOD1020'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD1020'][s['fl_acaod']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['AOD1020'][s['fl6']]),color='lightcoral',ymin=0, ymax=10,lw=2)

#plt.axvline(x=np.nanmedian(s['AOD1020'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD1020'][s['fl_acaod']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['AOD1020'][s['fl6']]),color='lightcoral',ymin=0, ymax=10,lw=2,ls='--')
pu.sub_note('b)')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,4,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False,loc=0)

plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_AOD_2wvl_histogram_stack.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_AOD_2wvl_histogram_stack.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_AOD_2wvl_histogram_stack.eps')


# In[41]:


np.nanmax(s['AOD0501'][s['fl6b']]),np.nanmax(s['AOD0501'][s['fl_acaod']])


# In[42]:


np.nanmin(s['AOD0501'][s['fl6b']]),np.nanmin(s['AOD0501'][s['fl_acaod']])


# In[45]:


np.nanmin(s['AOD1020'][s['fl_acaodb']]), np.nanmax(s['AOD1020'][s['fl_acaodb']])


# In[44]:


s['fl_acaodb'] = s['fl_acaod'] & (s['AOD1020']>0.0) & (s['AOD1020']< 1.2)


# ## Create a histogram of Angstrom at 2 wavelengths

# In[72]:


ia = 4 #660 nm
ib = 2 #500 nm


# In[73]:


awvl


# In[202]:


fig = plt.figure(figsize=(7,6))
plt.subplot(2,1,1)
plt.hist([s['angs'][s['fl_both'],ib],s['angs'][s['fl_acaod'],ib],s['angs'][s['fl6'],ib]],
        bins=15,range=[-0.2,2.5],label=['Combined','Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=True)
plt.xlim(-0.2,2.5)
#plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
#plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.axvline(x=np.nanmean(s['angs'][s['fl_both'],ib]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_both'],ib]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['angs'][s['fl_both'],ib]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl_acaod'],ib]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl6'],ib]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['angs'][s['fl_both'],ib]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_acaod'],ib]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl6'],ib]),color='r',ymin=0, ymax=10,lw=2,ls='--')


plt.xlabel('Angstrom at {:4.0f} nm'.format(awvl[ib]))
plt.ylabel('Normalized frequency counts')
#plt.title('ORACLES 2016 4STAR AOD')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,4,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False,loc=0)
#plt.legend(frameon=False)

fig.subplots_adjust(hspace=.3)
plt.subplot(2,1,2)
plt.hist([s['angs'][s['fl_both'],ia],s['angs'][s['fl_acaod'],ia],s['angs'][s['fl6'],ia]],
        bins=15,range=[-0.2,2.5],label=['Combined','Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=True)
plt.xlim(-0.2,2.5)
plt.ylabel('Normalized frequency counts')
plt.xlabel('Angstrom at {:4.0f} nm'.format(awvl[ia]))
#plt.title('ORACLES 4STAR column Angstrom exponent')

plt.axvline(x=np.nanmean(s['angs'][s['fl_both'],ia]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_both'],ia]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['angs'][s['fl_both'],ia]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl_acaod'],ia]),color='g',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl6'],ia]),color='r',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['angs'][s['fl_both'],ia]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_acaod'],ia]),color='g',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl6'],ia]),color='r',ymin=0, ymax=10,lw=2,ls='--')


plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_Angstrom_2wvl_histogram_comb.png',
            transparent=True,dpi=500)


# In[ ]:





# ### Remake the angstrom histogram as stacked bars , one with old method, and one with polyfit

# In[191]:


fig = plt.figure(figsize=(7,6))
plt.subplot(2,1,1)
plt.hist([s['angs'][s['fl6'],ib],s['angs'][s['fl_acaod'],ib]],stacked=True,color=['lightcoral','b'],
        bins=20,range=[-0.2,2.5],label=['Full Column','Only above cloud'],edgecolor='None',alpha=0.75,normed=False)
plt.xlim(-0.2,2.5)
#plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
#plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.axvline(x=np.nanmean(s['angs'][s['fl_both'],ib]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_both'],ib]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['angs'][s['fl_both'],ib]),color='k',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl_acaod'],ib]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs'][s['fl6'],ib]),color='lightcoral',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['angs'][s['fl_both'],ib]),color='k',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl_acaod'],ib]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs'][s['fl6'],ib]),color='lightcoral',ymin=0, ymax=10,lw=2,ls='--')


plt.xlabel('Angstrom evaluated at {:4.0f} nm'.format(awvl[ib]))
plt.ylabel('Samples')
pu.sub_note('a)')
#plt.title('ORACLES 2016 4STAR AOD')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False,loc=6)
#plt.legend(frameon=False)

fig.subplots_adjust(hspace=.3)
plt.subplot(2,1,2)
plt.hist([s['angs_470_865'][s['fl6']],s['angs_470_865'][s['fl_acaod']]],color=['lightcoral','b'],stacked=True,
        bins=20,range=[-0.2,2.5],label=['Only above cloud','Below 0.6 km'],edgecolor='None',alpha=0.75,normed=False)
plt.xlim(-0.2,2.5)
plt.ylabel('Samples')
pu.sub_note('b)')
plt.xlabel('Angstrom of 470 nm / 865 nm'.format(awvl[ia]))
#plt.title('ORACLES 4STAR column Angstrom exponent')

plt.axvline(x=np.nanmean(s['angs_470_865'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,label='Mean')
plt.axvline(x=np.nanmedian(s['angs_470_865'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,ls='--',label='Median')

plt.axvline(x=np.nanmean(s['angs_470_865'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs_470_865'][s['fl_acaod']]),color='b',ymin=0, ymax=10,lw=2)
plt.axvline(x=np.nanmean(s['angs_470_865'][s['fl6']]),color='lightcoral',ymin=0, ymax=10,lw=2)

plt.axvline(x=np.nanmedian(s['angs_470_865'][s['fl_both']]),color='k',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs_470_865'][s['fl_acaod']]),color='b',ymin=0, ymax=10,lw=2,ls='--')
plt.axvline(x=np.nanmedian(s['angs_470_865'][s['fl6']]),color='lightcoral',ymin=0, ymax=10,lw=2,ls='--')


plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_Angstrom_2wvl_histogram_stacked.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_Angstrom_2wvl_histogram_stacked.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_Angstrom_2wvl_histogram_stacked.eps')


# In[193]:


np.nanmean(s['angs_470_865'][s['fl_acaod']])


# In[194]:


np.nanmedian(s['angs_470_865'][s['fl_acaod']])


# In[195]:


np.nanstd(s['angs_470_865'][s['fl_acaod']])


# In[196]:


np.nanmean(s['angs_470_865'][s['fl6']]),np.nanmedian(s['angs_470_865'][s['fl6']]),np.nanstd(s['angs_470_865'][s['fl6']])


# In[155]:


np.nanmean(s['angs'][s['fl_acaod'],ib]),np.nanmean(s['angs'][s['fl6'],ib])


# In[197]:


np.nanmedian(s['angs'][s['fl_acaod'],ib]),np.nanmedian(s['angs'][s['fl6'],ib])


# In[198]:


np.nanstd(s['angs'][s['fl_acaod'],ib]),np.nanstd(s['angs'][s['fl6'],ib])


# In[2]:


1.25-1.08


# In[84]:


np.nanmin(s['AOD0501'][s['fl6b']]),np.nanmax(s['AOD0501'][s['fl6b']]) 


# ## Get vertical dependence of Angstrom exponent

# In[74]:


s['fl_QA_angs'] = s['fl_QA'] & (s['AOD0501']>0.1)


# In[75]:


s['fl_QA_angs_aca'] = s['fl_acaod'] & (s['AOD0501']>0.1) & (s['GPS_Alt']>300.0)


# In[81]:


s['fl_QA_angs_nowvb'] = s['fl_QA'] & (s['AOD0501']>0.1) & (s['Longitude']<13.5)


# In[47]:


plt.figure()
plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.')


# In[76]:


binned_ang,binned_alt,binned_num,binned_ndays = [],[],[],[]
for i in xrange(70):
    flaa = (s['GPS_Alt'][s['fl_QA_angs']]>=i*100.0) & (s['GPS_Alt'][s['fl_QA_angs']]<(i+1.0)*100.0)
    binned_ang.append(s['angs_470_865'][s['fl_QA_angs']][flaa])
    binned_alt.append(np.mean([i*100.0,(i+1.0)*100.0]))
    binned_num.append(len(s['angs_470_865'][s['fl_QA_angs']][flaa]))
    binned_ndays.append(len(np.unique(s['days'][s['fl_QA_angs']][flaa])))


# In[77]:


binned_angc,binned_altc,binned_numc,binned_ndaysc = [],[],[],[]
for i in xrange(70):
    flaa = (s['GPS_Alt'][s['fl_QA_angs_aca']]>=i*100.0) & (s['GPS_Alt'][s['fl_QA_angs_aca']]<(i+1.0)*100.0)
    binned_angc.append(s['angs_470_865'][s['fl_QA_angs_aca']][flaa])
    binned_altc.append(np.mean([i*100.0,(i+1.0)*100.0]))
    binned_numc.append(len(s['angs_470_865'][s['fl_QA_angs_aca']][flaa]))
    binned_ndaysc.append(len(np.unique(s['days'][s['fl_QA_angs_aca']][flaa])))


# In[82]:


binned_angnw,binned_altnw,binned_numnw,binned_ndaysnw = [],[],[],[]
for i in xrange(70):
    flaa = (s['GPS_Alt'][s['fl_QA_angs_nowvb']]>=i*100.0) & (s['GPS_Alt'][s['fl_QA_angs_nowvb']]<(i+1.0)*100.0)
    binned_angnw.append(s['angs_470_865'][s['fl_QA_angs_nowvb']][flaa])
    binned_altnw.append(np.mean([i*100.0,(i+1.0)*100.0]))
    binned_numnw.append(len(s['angs_470_865'][s['fl_QA_angs_nowvb']][flaa]))
    binned_ndaysnw.append(len(np.unique(s['days'][s['fl_QA_angs_nowvb']][flaa])))


# In[78]:


plt.figure(figsize=(6,7))
bp =plt.boxplot(binned_ang,positions=binned_alt,vert=False,showfliers=False,widths=100,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom of 470 nm / 865 nm')
plt.ylabel('Altitude [m]')
#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
for b in bp['boxes']:
    b.set_facecolor('green')
    b.set_edgecolor('green')
    b.set_alpha(0.4)
for b in bp['means']:
    b.set_marker('o')
    b.set_color('firebrick')
    b.set_alpha(0.4)
for b in bp['whiskers']:
    b.set_linestyle('-')
    b.set_color('green')
    b.set_alpha(0.4)
for b in bp['caps']:
    b.set_alpha(0.4)
    b.set_color('green')
for b in bp['medians']:
    b.set_linewidth(4)
    b.set_color('gold')
    b.set_alpha(0.4)
    
bpc =plt.boxplot(binned_angc,positions=binned_altc,vert=False,showfliers=False,widths=100,showmeans=True,patch_artist=True)
for b in bpc['boxes']:
    b.set_facecolor('blue')
    b.set_edgecolor('blue')
for b in bpc['means']:
    b.set_marker('o')
    b.set_color('firebrick')
for b in bpc['whiskers']:
    b.set_linestyle('-')
for b in bpc['medians']:
    b.set_linewidth(4)
    b.set_color('gold')
bpc['boxes'][0].set_color('grey')
ax = plt.gca()
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,2.2)
plt.grid()
for j,nn in enumerate(binned_ndays): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])*0.75,binned_alt[j],'{:2.0f}'.format(nn),alpha=0.4,color='g',fontsize=7,verticalalignment='center')
for j,nn in enumerate(binned_ndaysc): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])*0.8,binned_altc[j],'{:2.0f}'.format(nn),color='b',fontsize=7,verticalalignment='center')

plt.legend([bp['boxes'][0],bpc['boxes'][1],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['All data','Above cloud','Mean','Median','25% - 75%','min-max'],
           frameon=False,loc=2,numpoints=1)
plt.tight_layout()

#plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_Angstrom_2wvl_vertical.png',
#            transparent=True,dpi=500)


# In[83]:


plt.figure(figsize=(6,7))
bp =plt.boxplot(binned_ang,positions=binned_alt,vert=False,showfliers=False,widths=100,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom of 470 nm / 865 nm')
plt.ylabel('Altitude [m]')
#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
for b in bp['boxes']:
    b.set_facecolor('green')
    b.set_edgecolor('green')
    b.set_alpha(0.4)
for b in bp['means']:
    b.set_marker('o')
    b.set_color('firebrick')
    b.set_alpha(0.4)
for b in bp['whiskers']:
    b.set_linestyle('-')
    b.set_color('green')
    b.set_alpha(0.4)
for b in bp['caps']:
    b.set_alpha(0.4)
    b.set_color('green')
for b in bp['medians']:
    b.set_linewidth(4)
    b.set_color('gold')
    b.set_alpha(0.4)
    
bpc =plt.boxplot(binned_angnw,positions=binned_altnw,vert=False,showfliers=False,widths=100,showmeans=True,patch_artist=True)
for b in bpc['boxes']:
    b.set_facecolor('blue')
    b.set_edgecolor('blue')
for b in bpc['means']:
    b.set_marker('o')
    b.set_color('firebrick')
for b in bpc['whiskers']:
    b.set_linestyle('-')
for b in bpc['medians']:
    b.set_linewidth(4)
    b.set_color('gold')
bpc['boxes'][0].set_color('grey')
ax = plt.gca()
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,2.2)
plt.grid()
for j,nn in enumerate(binned_ndays): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])*0.75,binned_alt[j],'{:2.0f}'.format(nn),alpha=0.4,color='g',fontsize=7,verticalalignment='center')
for j,nn in enumerate(binned_ndaysnw): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])*0.8,binned_altnw[j],'{:2.0f}'.format(nn),color='b',fontsize=7,verticalalignment='center')

plt.legend([bp['boxes'][0],bpc['boxes'][1],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['All data','No WVB','Mean','Median','25% - 75%','min-max'],
           frameon=False,loc=2,numpoints=1)
plt.tight_layout()

#plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_Angstrom_2wvl_vertical.png',
#            transparent=True,dpi=500)


# In[198]:


plt.figure(figsize=(6,7))
bp =plt.boxplot(binned_ang,positions=np.array(binned_alt)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom of 470 nm / 865 nm')
plt.ylabel('Altitude [m]')
#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
gr = plt.cm.RdPu
bl = plt.cm.Blues
bndm = np.nanmax(binned_ndays)*1.0
for j,b in enumerate(bp['boxes']):
    b.set_facecolor(gr(binned_ndays[j]*1.0/bndm))
    b.set_edgecolor(gr(binned_ndays[j]*1.0/bndm))
    #b.set_alpha(0.4)
for j,b in enumerate(bp['means']):
    b.set_marker('.')
    b.set_color('None')
    b.set_markerfacecolor('darkgreen')
    b.set_markeredgecolor('darkgreen')
    b.set_alpha(0.6)
for j,b in enumerate(bp['whiskers']):
    b.set_linestyle('-')
    b.set_color('pink') #gr(binned_ndays[j]*1.0/bndm))
    b.set_alpha(0.7)
for j,b in enumerate(bp['caps']):
    b.set_alpha(0.7)
    b.set_color('pink')#gr(binned_ndays[j]*1.0/bndm))
for j,b in enumerate( bp['medians']):
    b.set_linewidth(4)
    b.set_color('gold')
    b.set_alpha(0.4)
    
bpc =plt.boxplot(binned_angc,positions=np.array(binned_altc)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
for j,b in enumerate(bpc['boxes']):
    b.set_facecolor(bl(binned_ndaysc[j]*1.0/bndm)) #'blue')
    b.set_edgecolor(bl(binned_ndaysc[j]*1.0/bndm))
for b in bpc['means']:
    b.set_marker('.')
    b.set_color('None')
    b.set_markerfacecolor('darkgreen')
    b.set_markeredgecolor('darkgreen')
    b.set_alpha(0.6)
for b in bpc['whiskers']:
    b.set_linestyle('-')
for b in bpc['medians']:
    b.set_linewidth(4)
    b.set_color('gold')
    b.set_alpha(0.6)
bpc['boxes'][0].set_color('grey')
ax = plt.gca()
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,2.2)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['All data','Above cloud','Mean','Median','25% - 75%','min-max'],
           frameon=False,loc=2,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(binned_ndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(binned_ndays)
cbaxesgr = plt.gcf().add_axes([0.25, 0.2, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.27, 0.2, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
#cbg.ax.yaxis.set_ticks_position('right')
#cbg.set_ticklabels(['0','3','6','9','12','15'],horizontalalignment='left')
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])#,cbg.set_label('Days sampled',verticalalignment='top')
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-4.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_Angstrom_2wvl_vertical_cb.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_Angstrom_2wvl_vertical_cb.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_Angstrom_2wvl_vertical_cb.eps')


# # Make some figures of the map statistic distribution

# ## plot a map of the mean, median, and std of the AOD

# In[46]:


a,xe,ye,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['AOD0501'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]],expand_binnumbers=True)
a = np.ma.masked_array(a,np.isnan(a))


# In[47]:


bn.shape


# In[48]:


np.nanmean(a)


# In[49]:


a.shape


# In[183]:


xe


# In[193]:


ye


# In[210]:


astd[17,5]


# In[50]:


aaa = (s['Latitude'][s['fl_acaod']]<xe[17]) & (s['Latitude'][s['fl_acaod']]>xe[16]) &      (s['Longitude'][s['fl_acaod']]>ye[5]) & (s['Longitude'][s['fl_acaod']]<ye[6])


# In[223]:


xe[17],xe[16]


# In[224]:


ye[5],ye[6]


# In[226]:


any(aaa)


# In[215]:


s['AOD0501'][s['fl_acaod']][aaa]


# In[53]:


xe[1]-xe[0]


# In[54]:


ye[1]-ye[0]


# In[186]:


plt.figure()
p = plt.pcolor(ye,xe,a,vmin=0.0,vmax=1.0,cmap='plasma')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')
plt.title('Mean AOD 501 nm')

cb = plt.colorbar(p,extend='both')
cb.set_label('Mean Above cloud AOD 501 nm')
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_meanAOD_map.png',
            transparent=True,dpi=500)


# In[68]:


astd,xe,ye,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['AOD0501'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]],statistic=np.std)
astd = np.ma.masked_array(astd,np.isnan(astd))
np.nanmean(astd)


# In[185]:


plt.figure()
p = plt.pcolor(ye,xe,astd,vmin=0.0,vmax=0.2,cmap='viridis')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')
plt.title('STD AOD 501 nm')

cb = plt.colorbar(p,extend='both')
cb.set_label('Standard deivation of Above cloud AOD 501 nm')
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_stdAOD_map.png',
            transparent=True,dpi=500)


# In[56]:


med,xe,ye,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['AOD0501'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]],statistic='median')
med = np.ma.masked_array(med,np.isnan(med))
np.nanmedian(med)


# In[184]:


plt.figure()
p = plt.pcolor(ye,xe,med,vmin=0.0,vmax=1.0,cmap='gist_rainbow')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')

plt.title('Median AOD 501 nm')

cb = plt.colorbar(p,extend='both')
cb.set_label('Median Above cloud AOD 501 nm')
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_medianAOD_map.png',
            transparent=True,dpi=500)


# ### Do the uncertainty at 501 and 1020

# In[51]:



stata05uc = {}


# In[52]:


stata05uc['mean'],stata05uc['x'],stata05uc['y'],stata05uc['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['UNCAOD0501'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]])
stata05uc['mean'] = np.ma.masked_array(stata05uc['mean'],np.isnan(stata05uc['mean'])|(stata05tc['mean']>1.5))

stata05uc['median'],stata05uc['xe'],stata05uc['ye'],stata05uc['bine'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['UNCAOD0501'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic='median')
stata05uc['median'] = np.ma.masked_array(stata05uc['median'],np.isnan(stata05uc['median'])|(stata05tc['mean']>1.5))

stata05uc['std'],stata05uc['xs'],stata05uc['ys'],stata05uc['bins'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['UNCAOD0501'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic=np.nanstd)
stata05uc['std'] = np.ma.masked_array(stata05uc['std'],np.isnan(stata05uc['std'])|(stata05tc['mean']>1.5))

np.nanmean(stata05uc['mean']),np.nanmedian(stata05uc['median']), np.nanmean(stata05uc['std'])


# In[53]:


stata10uc = {}


# In[54]:


stata10uc['mean'],stata10uc['x'],stata10uc['y'],stata10uc['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['UNCAOD1020'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]])
stata10uc['mean'] = np.ma.masked_array(stata10uc['mean'],np.isnan(stata10uc['mean'])|(stata05tc['mean']>1.5))

stata10uc['median'],stata10uc['xe'],stata10uc['ye'],stata10uc['bine'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['UNCAOD1020'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic='median')
stata10uc['median'] = np.ma.masked_array(stata05uc['median'],np.isnan(stata05uc['median'])|(stata05tc['mean']>1.5))

stata10uc['std'],stata10uc['xs'],stata10uc['ys'],stata10uc['bins'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['UNCAOD1020'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic=np.nanstd)
stata10uc['std'] = np.ma.masked_array(stata10uc['std'],np.isnan(stata10uc['std'])|(stata05tc['mean']>1.5))

np.nanmean(stata10uc['mean']),np.nanmedian(stata10uc['median']), np.nanmean(stata10uc['std'])


# ### Redo but for ACAOD 1020

# In[55]:


stata10 = {}


# In[152]:


stata10['mean'],stata10['x'],stata10['y'],stata10['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['AOD1020'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]])
stata10['mean'] = np.ma.masked_array(stata10['mean'],np.isnan(stata10['mean']))
np.nanmean(stata10['mean'])


# In[153]:


stata10['median'],stata10['xe'],stata10['ye'],stata10['bine'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['AOD1020'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic='median')
stata10['median'] = np.ma.masked_array(stata10['median'],np.isnan(stata10['median']))
np.nanmedian(stata10['median'])


# In[154]:


stata10['std'],stata10['xs'],stata10['ys'],stata10['bins'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['AOD1020'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic=np.nanstd)
stata10['std'] = np.ma.masked_array(stata10['std'],np.isnan(stata10['std']))
np.nanmean(stata10['std'])


# ### Redo for total column AOD at 501 and 1020

# In[56]:


stata10tc = {}


# In[57]:


stata10tc['mean'],stata10tc['x'],stata10tc['y'],stata10tc['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['AOD1020'][s['fl6']],
                                                                          bins=26,range=[[-25,-8],[0,16]])
stata10tc['mean'] = np.ma.masked_array(stata10tc['mean'],np.isnan(stata10tc['mean'])|(stata10tc['mean']>1.5))

stata10tc['median'],stata10tc['xe'],stata10tc['ye'],stata10tc['bine'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['AOD1020'][s['fl6']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic='median')
stata10tc['median'] = np.ma.masked_array(stata10tc['median'],np.isnan(stata10tc['median'])|(stata10tc['mean']>1.5))

stata10tc['std'],stata10tc['xs'],stata10tc['ys'],stata10tc['bins'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['AOD1020'][s['fl6']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic=np.nanstd)
stata10tc['std'] = np.ma.masked_array(stata10tc['std'],np.isnan(stata10tc['std'])|(stata10tc['mean']>1.5))

np.nanmean(stata10tc['mean']),np.nanmedian(stata10tc['median']), np.nanmean(stata10tc['std'])


# In[58]:


stata05uc = {}


# In[59]:


stata05uc['mean'],stata05uc['x'],stata05uc['y'],stata05uc['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6b']],
                                                                          s['Longitude'][s['fl6b']],
                                                                          s['AOD0501'][s['fl6b']],
                                                                          bins=26,range=[[-25,-8],[0,16]])
stata05uc['mean'] = np.ma.masked_array(stata05uc['mean'],np.isnan(stata05uc['mean'])|(stata05uc['mean']>1.5))

stata05uc['median'],stata05uc['xe'],stata05uc['ye'],stata05uc['bine'] = st.binned_statistic_2d(s['Latitude'][s['fl6b']],
                                                                          s['Longitude'][s['fl6b']],
                                                                          s['AOD0501'][s['fl6b']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic='median')
stata05uc['median'] = np.ma.masked_array(stata05uc['median'],np.isnan(stata05uc['median'])|(stata05uc['mean']>1.5))

stata05uc['std'],stata05uc['xs'],stata05uc['ys'],stata05uc['bins'] = st.binned_statistic_2d(s['Latitude'][s['fl6b']],
                                                                          s['Longitude'][s['fl6b']],
                                                                          s['AOD0501'][s['fl6b']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                          statistic=np.nanstd)
stata05uc['std'] = np.ma.masked_array(stata05uc['std'],np.isnan(stata05uc['std'])|(stata05uc['mean']>1.5))

np.nanmean(stata05uc['mean']),np.nanmedian(stata05uc['median']), np.nanmean(stata05uc['std'])


# In[88]:


stata05uc['dcnt'],stata05uc['dxs'],stata05uc['dys'],stata05uc['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6b']],
                                                                          s['Longitude'][s['fl6b']],
                                                                          s['days'][s['fl6b']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic=uniq_cnt)
stata05uc['dcnt'] = np.ma.masked_array(stata05uc['dcnt'],np.isnan(stata05uc['dcnt']))


# ### combine plots on one figure

# In[222]:


fig = plt.figure(figsize=(15,4.5))
ax1 = plt.subplot(1,3,1)

p = plt.pcolor(ye,xe,a,vmin=0.0,vmax=1.0,cmap='plasma')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')
plt.title('Mean')

cb = plt.colorbar(p,extend='both')
#cb.set_label('Mean Above cloud AOD 501 nm')


ax2 = plt.subplot(1,3,2)
p = plt.pcolor(ye,xe,med,vmin=0.0,vmax=1.0,cmap='plasma')

plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')

plt.title('Median')

cb = plt.colorbar(p,extend='both')
#cb.set_label('Median Above cloud AOD 501 nm')

ax3 = plt.subplot(1,3,3)
p = plt.pcolor(ye,xe,astd,vmin=0.0,vmax=0.2,cmap='viridis')

plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
plt.title('STD')

cb = plt.colorbar(p,extend='both')
#cb.set_label('Standard deivation of Above cloud AOD 501 nm')

#plt.tight_layout()
plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.06)
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAOD_map.png',
            transparent=True,dpi=500)


# In[75]:


fig = plt.figure(figsize=(15,4.5))
ax1 = plt.subplot(1,3,1)

p = plt.pcolor(ye,xe,a,vmin=0.0,vmax=1.0,cmap='plasma')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')
plt.title('Mean ACAOD')

cb = plt.colorbar(p,extend='both')
#cb.set_label('Mean Above cloud AOD 501 nm')


ax2 = plt.subplot(1,3,2)
stata05uc['mean'],stata05uc['x'],stata05uc['y']
p = plt.pcolor(stata05uc['y'],stata05uc['x'],stata05uc['mean'],vmin=0.0,vmax=1.0,cmap='plasma')

plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')

plt.title('Full column')

cb = plt.colorbar(p,extend='both')
#cb.set_label('Median Above cloud AOD 501 nm')

ax3 = plt.subplot(1,3,3)
p = plt.pcolor(ye,xe,astd,vmin=0.0,vmax=0.2,cmap='viridis')

plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
plt.title('STD')

cb = plt.colorbar(p,extend='both')
#cb.set_label('Standard deivation of Above cloud AOD 501 nm')

#plt.tight_layout()
plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.06)
#plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAOD_map.png',
#            transparent=True,dpi=500)


# In[61]:


cnt,xe,ye,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['AOD0501'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]],statistic='count')
cnt = np.ma.masked_array(cnt,np.isnan(cnt))


# In[66]:


dcnt,xed,yed,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['days'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]],statistic=uniq_cnt)
dcnt = np.ma.masked_array(dcnt,np.isnan(dcnt))


# In[69]:


io = np.where((cnt.data>0.0) & (astd.data<1.0))


# In[70]:


cnt.shape


# In[81]:


yy = [(y+ye[i+1])/2.0 for i,y in enumerate(ye[0:-1])]
yys = np.array(yy*26).reshape(26,-1)


# In[82]:


xx = [(x+xe[i+1])/2.0 for i,x in enumerate(xe[0:-1])]
xxs = np.array(xx*26).reshape(26,-1)


# In[83]:


cnt.data[io[0],io[1]].flatten()


# In[372]:


fig = plt.figure(figsize=(11,4.5))
ax1 = plt.subplot(1,2,1)

p = plt.pcolor(ye,xe,a,vmin=0.0,vmax=1.0,cmap='plasma')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')
plt.title('Mean')

cb = plt.colorbar(p,extend='both')


ax2 = plt.subplot(1,2,2)

p = plt.scatter(np.array(yy)[io[1]],np.array(xx)[io[0]],cnt.data[io[0],io[1]].flatten()/10,
                c=astd.data[io[0],io[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis')
plt.xlim(0,16)
plt.ylim(-26,-8)
plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
plt.title('STD')

cb = plt.colorbar(p,extend='max')
#cb.set_label('Standard deivation of Above cloud AOD 501 nm')

sizes = [10,100,500,1500]
labels = ['{0}'.format(z) for z in sizes]
points = [ax.scatter([], [], s=z/10, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')
    

#plt.tight_layout()
plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.06)
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAOD_2panel_map.png',
            transparent=True,dpi=500)


# Get the mean and std from the binned values, to try to remove the sampling issues

# In[112]:


a.shape


# In[113]:


np.nanmean(a),np.nanmedian(a),np.nanmean(astd.data)


# In[193]:


plt.figure()
plt.hist(a.flatten(),25,range=(0,1.2),edgecolor='None',alpha=0.6)


# ## Plot a map of the distribution of the Angstrom exponent

# In[67]:


astat,astat2 = {},{}
castat,castat2 = {},{}


# In[68]:


ia = 2


# In[203]:


astat['mean'],astat['x'],astat['y'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs'][s['fl_acaod'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]])
astat['mean'] = np.ma.masked_array(astat['mean'],np.isnan(astat['mean']))
np.nanmean(astat['mean'])


# In[207]:


castat['mean'],castat['x'],castat['y'],castat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['angs'][s['fl6'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]])
castat['mean'] = np.ma.masked_array(castat['mean'],np.isnan(castat['mean']))
np.nanmean(castat['mean'])


# In[202]:


astat2['mean'],astat2['x'],astat2['y'],astat2['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs_470_865'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]])
astat2['mean'] = np.ma.masked_array(astat2['mean'],np.isnan(astat2['mean']))


# In[208]:


castat2['mean'],castat2['x'],castat2['y'],castat2['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['angs_470_865'][s['fl6']],
                                                                          bins=26,range=[[-25,-8],[0,16]])
castat2['mean'] = np.ma.masked_array(castat2['mean'],np.isnan(castat2['mean']))
np.nanmean(castat2['mean'])


# In[175]:


plt.figure()
p = plt.pcolor(astat['y'],astat['x'],astat['mean'],vmin=0.0,vmax=2.2,cmap='plasma')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')

plt.title('Mean Angstrom')

cb = plt.colorbar(p,extend='both')
cb.set_label('Mean Angstrom exponent \nof Above cloud AOD at 660 nm')
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_meanAngstrom_map.png',
            transparent=True,dpi=500)


# In[204]:


astat['median'],astat['xe'],astat['ye'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs'][s['fl_acaod'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic='median')
astat['median'] = np.ma.masked_array(astat['median'],np.isnan(astat['median']))
np.nanmedian(astat['median'])


# In[209]:


castat['median'],castat['xe'],castat['ye'],castat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['angs'][s['fl6'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic='median')
castat['median'] = np.ma.masked_array(castat['median'],np.isnan(castat['median']))
np.nanmedian(castat['median'])


# In[110]:


astat2['median'],astat2['xe'],astat2['ye'],astat2['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs_470_865'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic='median')
astat2['median'] = np.ma.masked_array(astat2['median'],np.isnan(astat2['median']))


# In[210]:


castat2['median'],castat2['xe'],castat2['ye'],castat2['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['angs_470_865'][s['fl6']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic='median')
castat2['median'] = np.ma.masked_array(castat2['median'],np.isnan(castat2['median']))
np.nanmedian(castat2['median'])


# In[211]:


plt.figure()
p = plt.pcolor(astat['ye'],astat['xe'],astat['median'],vmin=0.0,vmax=2.2,cmap='hsv')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')
plt.title('Median Angstrom')
cb = plt.colorbar(p,extend='both')
cb.set_label('Median Angstrom exponent \nof Above cloud AOD at 660 nm')
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_medianAngstrom_map.png',
            transparent=True,dpi=500)


# In[205]:


astat['std'],astat['xs'],astat['ys'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs'][s['fl_acaod'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic=np.std)
astat['std'] = np.ma.masked_array(astat['std'],np.isnan(astat['std']))
np.nanmean(astat['std'])


# In[211]:


castat['std'],castat['xs'],castat['ys'],castat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['angs'][s['fl6'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic=np.std)
castat['std'] = np.ma.masked_array(castat['std'],np.isnan(castat['std']))
np.nanmean(castat['std'])


# In[178]:


astat2['std'],astat2['xs'],astat2['ys'],astat2['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs_470_865'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic=np.std)
astat2['std'] = np.ma.masked_array(astat2['std'],np.isnan(astat2['std']))


# In[212]:


castat2['std'],castat2['xs'],castat2['ys'],castat2['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl6']],
                                                                          s['Longitude'][s['fl6']],
                                                                          s['angs_470_865'][s['fl6']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic=np.std)
castat2['std'] = np.ma.masked_array(castat2['std'],np.isnan(castat2['std']))
np.nanmean(castat2['std'])


# In[557]:


plt.figure()
p = plt.pcolor(astat['ys'],astat['xs'],astat['std'],vmin=0.0,vmax=0.2,cmap='viridis')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')
plt.title('std Angstrom')

cb = plt.colorbar(p,extend='both')
cb.set_label('Standard deviation \nAngstrom exponent of Above cloud AOD at 660 nm')
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_stdAngstrom_map.png',
            transparent=True,dpi=500)


# ### Combine stats of Angstrom into single figure

# In[223]:


fig = plt.figure(figsize=(15,4.5))
ax1 = plt.subplot(1,3,1)

p = plt.pcolor(astat['y'],astat['x'],astat['mean'],vmin=0.0,vmax=2.2,cmap='gist_rainbow')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')

plt.title('Mean')

cb = plt.colorbar(p,extend='both')


ax2 = plt.subplot(1,3,2)
p = plt.pcolor(astat['ye'],astat['xe'],astat['median'],vmin=0.0,vmax=2.2,cmap='gist_rainbow')

plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
plt.title('Median')
cb = plt.colorbar(p,extend='both')

ax3 = plt.subplot(1,3,3)
p = plt.pcolor(astat['ys'],astat['xs'],astat['std'],vmin=0.0,vmax=0.2,cmap='viridis')

plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
plt.title('STD')

cb = plt.colorbar(p,extend='both')

plt.suptitle('Above cloud Angstrom exponent at 660 nm\n')
plt.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAngstrom_map.png',
            transparent=True,dpi=500)


# In[72]:


astat['cnt'],astat['xs'],astat['ys'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs'][s['fl_acaod'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic='count')
astat['cnt'] = np.ma.masked_array(astat['cnt'],np.isnan(astat['cnt']))


# In[73]:


astat2['cnt'],astat2['xs'],astat2['ys'],astat2['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs_470_865'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic='count')
astat2['cnt'] = np.ma.masked_array(astat2['cnt'],np.isnan(astat2['cnt']))


# In[115]:


s['days']


# In[65]:


uniq_cnt = lambda x: len(np.unique(x))


# In[89]:


astat['dcnt'],astat['dxs'],astat['dys'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['days'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic=uniq_cnt)
astat['dcnt'] = np.ma.masked_array(astat['dcnt'],np.isnan(astat['dcnt']))


# In[70]:


astat2['dcnt'],astat2['xs'],astat2['ys'],astat2['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['days'][s['fl_acaod']],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic=uniq_cnt)
astat2['dcnt'] = np.ma.masked_array(astat2['dcnt'],np.isnan(astat2['dcnt']))


# In[ ]:





# In[71]:


iao = np.where((astat['cnt'].data>0.0) & (astat['std'].data<1.0))
iao2 = np.where((astat2['cnt'].data>0.0) & (astat2['std'].data<1.0))


# In[409]:


fig = plt.figure(figsize=(11,4.5))
ax1 = plt.subplot(1,2,1)

p = plt.pcolor(astat['y'],astat['x'],astat['mean'],vmin=0.0,vmax=2.2,cmap='gist_rainbow')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')

plt.title('Mean')

cb = plt.colorbar(p,extend='both')


ax2 = plt.subplot(1,2,2)

p = plt.scatter(np.array(yy)[iao[1]],np.array(xx)[iao[0]],astat['cnt'].data[iao[0],iao[1]].flatten()/10,
                c=astat['std'].data[iao[0],iao[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis')
plt.xlim(0,16)
plt.ylim(-26,-8)
plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
plt.title('STD')

cb = plt.colorbar(p,extend='max')
#cb.set_label('Standard deivation of Above cloud AOD 501 nm')

sizes = [10,100,500,1500]
labels = ['{0}'.format(z) for z in sizes]
points = [ax.scatter([], [], s=z/10, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')
    

plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.06)


plt.suptitle('Above cloud Angstrom exponent at 660 nm\n')
plt.subplots_adjust(left=0.06, right=0.98, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAngstrom_2panel_map.png',
            transparent=True,dpi=500)


# ### Add the coastline to the figures

# In[79]:


def mapfig(ax=plt.gca(),no_ylabel=False):
    m = Basemap(projection='merc',llcrnrlat=-26,urcrnrlat=-8,llcrnrlon=0,urcrnrlon=16,resolution='l',ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.linspace(0,16,9),labels=[0,0,0,1],linewidth=0.1)
    if no_ylabel:
        m.drawparallels(np.linspace(-26,-8,10),labels=[0,0,0,0],linewidth=0.1)
    else:
        m.drawparallels(np.linspace(-26,-8,10),labels=[1,0,0,0],linewidth=0.1)
    m.shadedrelief(alpha=0.4)
    return m


# In[142]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax = ax.flatten()
ax1 = ax[0]
#ax1 = plt.subplot(1,2,1)
m = mapfig(ax=ax1)

mx,my = m(astat['y'],astat['x'])
p = ax1.pcolor(mx,my,astat['mean'],vmin=0.8,vmax=1.6,cmap='magma')
ax1.set_title('Mean')
cb = m.colorbar(p,extend='both')

ax2 = ax[1]
#ax2 = plt.subplot(1,2,2)
m2 = mapfig(ax=ax2)
mxx,myy = m2(np.array(yy),np.array(xx))
p2 = ax2.scatter(mxx[iao[1]],myy[iao[0]],astat['cnt'].data[iao[0],iao[1]].flatten()/5,
                c=astat['std'].data[iao[0],iao[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis')
#plt.xlim(0,16)
#plt.ylim(-26,-8)
#plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
ax2.set_title('STD')
cb = m2.colorbar(p2,extend='max')
#cb.set_label('Standard deivation of Above cloud AOD 501 nm')

sizes = [10,100,500,1500]
labels = ['{0}'.format(z) for z in sizes]
points = [ax2.scatter([], [], s=z/5, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')
    
plt.suptitle('Above cloud Angstrom exponent at 500 nm\n')
plt.subplots_adjust(left=0.03, right=0.97, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAngstrom500nm_2panel_actualmap_focuscmap.png',
            transparent=True,dpi=500)


# In[185]:


astat2['yy'] = np.array([(y+astat2['ys'][i+1])/2.0 for i,y in enumerate(astat2['ys'][0:-1])])
astat2['xx'] = np.array([(x+astat2['xs'][i+1])/2.0 for i,x in enumerate(astat2['xs'][0:-1])])


# In[334]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax = ax.flatten()
ax1 = ax[0]
#ax1 = plt.subplot(1,2,1)
m = mapfig(ax=ax1)

mx,my = m(astat2['y'],astat2['x'])
p = ax1.pcolor(mx,my,astat2['mean'],vmin=0.2,vmax=2.0,cmap='gist_ncar')
ax1.set_title('Mean')
cb = m.colorbar(p,extend='both')

ax2 = ax[1]
#ax2 = plt.subplot(1,2,2)
m2 = mapfig(ax=ax2)
mxx,myy = m2(astat2['yy'],astat2['xx'])
p2 = ax2.scatter(mxx[iao2[1]],myy[iao2[0]],astat2['cnt'].data[iao2[0],iao2[1]].flatten()/5,
                c=astat2['std'].data[iao2[0],iao2[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0.0,vmax=0.4)
#plt.xlim(0,16)
#plt.ylim(-26,-8)
#plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
ax2.set_title('STD')
cb = m2.colorbar(p2,extend='max')
#cb.set_label('Standard deivation of Above cloud AOD 501 nm')

sizes = [10,100,500,1500]
labels = ['{0}'.format(z) for z in sizes]
points = [ax2.scatter([], [], s=z/5, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')
    
plt.suptitle('Above cloud Angstrom of 470 - 865 nm\n')
plt.subplots_adjust(left=0.03, right=0.97, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAngstrom_470_865_2panel_actualmap.png',
            transparent=True,dpi=500)


# In[190]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax = ax.flatten()
ax1 = ax[0]
#ax1 = plt.subplot(1,2,1)
m = mapfig(ax=ax1)

mx,my = m(astat2['y'],astat2['x'])
p = ax1.pcolor(mx,my,astat2['mean'],vmin=1.2,vmax=2.0,cmap='magma')
ax1.set_title('Mean')
pu.sub_note('b)')
cb = m.colorbar(p,extend='both')
cb.set_label('Mean Angstrom exponent')

ax2 = ax[1]
#ax2 = plt.subplot(1,2,2)
m2 = mapfig(ax=ax2)
mxx,myy = m2(astat2['yy'],astat2['xx'])
p2 = ax2.scatter(mxx[iao2[1]],myy[iao2[0]],astat2['dcnt'].data[iao2[0],iao2[1]].flatten()*30,
                c=astat2['std'].data[iao2[0],iao2[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0.0,vmax=0.4)
#plt.xlim(0,16)
#plt.ylim(-26,-8)
#plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
ax2.set_title('STD')
pu.sub_note('a)')
cb = m2.colorbar(p2,extend='max')
cb.set_label('Standard deviation of Angstrom exponent')
#cb.set_label('Standard deivation of Above cloud AOD 501 nm')

sizes = [1,2,3,5] #[10,100,500,1500]
labels = ['{0} day{1}'.format(z,'s' if z>1 else '') for z in sizes]
points = [ax2.scatter([], [], s=z*30, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left',handletextpad=0.1)
    
#plt.suptitle('Above cloud Angstrom of 470 nm / 865 nm\n')
plt.subplots_adjust(left=0.03, right=0.94, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAngstrom_470_865_2panel_actualmap_days_focus_cmap.png',
            transparent=True,dpi=500)

plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAngstrom_470_865_2panel_actualmap_days_focus_cmap.eps')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAngstrom_470_865_2panel_actualmap_days_focus_cmap.pdf')


# In[135]:


astat2['mean'].mean(),astat2['median'].mean(),astat2['std'].mean()


# In[165]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax = ax.flatten()
ax1 = ax[0]
m = mapfig(ax=ax1)

mxe,mye = m(ye,xe)
p = ax1.pcolor(mxe,mye,a,vmin=0.0,vmax=1.0,cmap='plasma')
ax1.set_title('Mean')
cb = m.colorbar(p,extend='max')

ax2 = ax[1]
m2 = mapfig(ax=ax2)
mxxe,myye = m2(np.array(yy)[io[1]],np.array(xx)[io[0]])
p2 = ax2.scatter(mxxe,myye,cnt.data[io[0],io[1]].flatten()/5.0,
                c=astd.data[io[0],io[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis')
ax2.set_title('STD')
cb = m2.colorbar(p2,extend='max')

sizes = [10,100,500,1500]
labels = ['{0}'.format(z) for z in sizes]
points = [ax2.scatter([], [], s=z/5.0, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')
    
plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.03, right=0.97, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAOD_2panel_actualmap.png',
            transparent=True,dpi=500)


# Now calculate the averages and std, from the binned statistics

# Calculate the difference in Angstrom exponent methods

# In[150]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax = ax.flatten()
ax1 = ax[0]
#ax1 = plt.subplot(1,2,1)
m = mapfig(ax=ax1)

mx,my = m(astat['y'],astat['x'])
p = ax1.pcolor(mx,my,astat['mean']-astat2['mean'],vmin=-0.4,vmax=-0.1,cmap='magma')
ax1.set_title('Mean')
cb = m.colorbar(p,extend='both')
cb.set_label('difference AE$_{{500}}$-AE$_{{470/865}}$')

ax2 = ax[1]
#ax2 = plt.subplot(1,2,2)
m2 = mapfig(ax=ax2)
mxx,myy = m2(np.array(yy),np.array(xx))
p2 = ax2.scatter(mxx[iao[1]],myy[iao[0]],astat['cnt'].data[iao[0],iao[1]].flatten()/5,
                c=astat['std'].data[iao[0],iao[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis')
#plt.xlim(0,16)
#plt.ylim(-26,-8)
#plt.xlabel(u'Longitude [°]')
#plt.ylabel(u'Latitude [°]')
ax2.set_title('STD')
cb = m2.colorbar(p2,extend='max')
#cb.set_label('Standard deivation of Above cloud AOD 501 nm')

sizes = [10,100,500,1500]
labels = ['{0}'.format(z) for z in sizes]
points = [ax2.scatter([], [], s=z/5, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')
    
plt.suptitle('Above cloud Angstrom exponent at 500 nm\n')
plt.subplots_adjust(left=0.03, right=0.97, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_stats_diffAngstrom_2panel_actualmap.png',
            transparent=True,dpi=500)


# ## Check the variability at 8-10 E and 18-22 S

# In[103]:


flsub = (s['Latitude'][s['fl_acaod']]>-22.0) & (s['Latitude'][s['fl_acaod']]<-18.0) & (s['Longitude'][s['fl_acaod']]>8.0) & (s['Longitude'][s['fl_acaod']]<10.0)


# In[105]:


any(flsub)


# In[106]:


plt.figure()
plt.hist(s['AOD0501'][s['fl_acaod']][flsub])


# In[120]:


plt.figure()
plt.hist(s['angs_470_865'][s['fl_acaod']][flsub])


# In[122]:


plt.figure()
plt.plot(s['AOD0501'][s['fl_acaod']][flsub],s['angs_470_865'][s['fl_acaod']][flsub],'.')
plt.ylabel('AE')

plt.xlabel('AOD 500nm')
plt.title('betwen 8-10E, 18-22S')


# # Load the AERONET sites, and put them on the map

# ## get online AERONET data and save

# In[27]:


import aeronet
reload(aeronet)


# In[68]:


anet = aeronet.get_aeronet(daystr='2016-08-15',daystr2='2016-09-30',lat_range=[-8.0,-26.0],lon_range=[0,16.0],
                           avg=False,lev='SDA20',version='3')


# In[69]:


anet.keys()



# In[40]:


sio.savemat(fp+'/data/aeronet_20160801.mat',anet)


# ## Alternatively load the AERONET from file

# In[28]:


anet = sio.loadmat(fp+'/data/aeronet_20160801.mat')


# In[29]:


kanet = anet.keys()
kanet.sort()
kanet


# In[30]:


anet_bad = np.where(anet == -999.0)


# In[51]:


aod_aeronet = {}
anet_lat,anet_lon,anet_mean_fineaod,anet_mean_coarseaod,anet_mean_totaod,anet_name = [],[],[],[],[],[]
anet_std_fineaod,anet_cnt_fineaod = [],[]
for aa in np.unique(anet['AERONET_Site_Name']):
    if type(aa) is np.ndarray: aa = aa[0]
    f_anet = anet['AERONET_Site_Name']==aa
    aod_aeronet[aa] = {'aod_fine':anet['Fine_Mode_AOD_500nm[tau_f]'][f_anet],
                   'aod_coarse':anet['Coarse_Mode_AOD_500nm[tau_c]'][f_anet],
                   'aod_tot':anet['Total_AOD_500nm[tau_a]'][f_anet],
                   'doy':anet[' Day_of_Year(fraction)'][f_anet],
                   'lat':anet['Site_Latitude(Degrees)'][f_anet][0],
                   'lon':anet['Site_Longitude(Degrees)'][f_anet][0],
                   'fmf':anet['FineModeFraction_500nm[eta]'][f_anet]}
    anet_name.append(aa)
    anet_lon.append(anet['Site_Longitude(Degrees)'][f_anet][0])
    anet_lat.append(anet['Site_Latitude(Degrees)'][f_anet][0])

for aa in anet_name:
    aod_aeronet[aa]['mean_aod_fine'] = np.nanmean(aod_aeronet[aa]['aod_fine'][aod_aeronet[aa]['aod_fine']!=-999.0])
    aod_aeronet[aa]['median_aod_fine'] = np.nanmedian(aod_aeronet[aa]['aod_fine'][aod_aeronet[aa]['aod_fine']!=-999.0])
    aod_aeronet[aa]['std_aod_fine'] = np.nanstd(aod_aeronet[aa]['aod_fine'][aod_aeronet[aa]['aod_fine']!=-999.0])
    aod_aeronet[aa]['mean_aod_coarse'] = np.nanmean(aod_aeronet[aa]['aod_coarse'][aod_aeronet[aa]['aod_coarse']!=-999.0])
    aod_aeronet[aa]['median_aod_coarse'] = np.nanmedian(aod_aeronet[aa]['aod_coarse'][aod_aeronet[aa]['aod_coarse']!=-999.0])
    aod_aeronet[aa]['std_aod_coarse'] = np.nanstd(aod_aeronet[aa]['aod_coarse'][aod_aeronet[aa]['aod_coarse']!=-999.0])
    aod_aeronet[aa]['mean_aod_tot'] = np.nanmean(aod_aeronet[aa]['aod_tot'][aod_aeronet[aa]['aod_tot']!=-999.0])
    aod_aeronet[aa]['median_aod_tot'] = np.nanmedian(aod_aeronet[aa]['aod_tot'][aod_aeronet[aa]['aod_tot']!=-999.0])
    aod_aeronet[aa]['std_aod_tot'] = np.nanstd(aod_aeronet[aa]['aod_tot'][aod_aeronet[aa]['aod_tot']!=-999.0])

    anet_mean_fineaod.append(aod_aeronet[aa]['mean_aod_fine'])
    anet_mean_coarseaod.append(aod_aeronet[aa]['mean_aod_coarse'])
    anet_mean_totaod.append(aod_aeronet[aa]['mean_aod_tot'])
    
    anet_std_fineaod.append(aod_aeronet[aa]['std_aod_fine'])
    anet_cnt_fineaod.append(len(aod_aeronet[aa]['aod_fine'][aod_aeronet[aa]['aod_fine']!=-999.0]))
    
    


# In[32]:


anet_name


# In[33]:


anet_lon


# In[34]:


anet_lat


# In[35]:


aod_aeronet['Walvis_Bay_airport']


# In[55]:


plt.figure()
plt.plot(1-aod_aeronet['Walvis_Bay_airport']['fmf'][aod_aeronet['Walvis_Bay_airport']['fmf']!=-999.0],'.')


# In[56]:


np.nanmean(1-aod_aeronet['Walvis_Bay_airport']['fmf'][aod_aeronet['Walvis_Bay_airport']['fmf']!=-999.0]), np.nanstd(1-aod_aeronet['Walvis_Bay_airport']['fmf'][aod_aeronet['Walvis_Bay_airport']['fmf']!=-999.0])


# In[41]:


plt.figure()
plt.plot(aod_aeronet['Walvis_Bay_airport']['aod_coarse'],'.',label='coarse')
plt.plot(aod_aeronet['Walvis_Bay_airport']['aod_fine'],'.',label='fine')
plt.plot(aod_aeronet['Walvis_Bay_airport']['aod_tot'],'.',label='total')
plt.legend()
plt.ylim(0,0.4)


# In[43]:


aod_aeronet['Walvis_Bay_airport']['aod_coarse'][aod_aeronet['Walvis_Bay_airport']['aod_coarse']==-999.0] = np.nan


# In[44]:


aod_aeronet['Walvis_Bay_airport']['aod_fine'][aod_aeronet['Walvis_Bay_airport']['aod_fine']==-999.0] = np.nan


# In[46]:


aod_aeronet['Walvis_Bay_airport']['aod_tot'][aod_aeronet['Walvis_Bay_airport']['aod_tot']==-999.0] = np.nan


# In[47]:


np.nanstd(aod_aeronet['Walvis_Bay_airport']['aod_coarse']), np.nanstd(aod_aeronet['Walvis_Bay_airport']['aod_fine']),np.nanstd(aod_aeronet['Walvis_Bay_airport']['aod_tot'])


# In[48]:


np.nanmean(aod_aeronet['Walvis_Bay_airport']['aod_coarse']), np.nanmean(aod_aeronet['Walvis_Bay_airport']['aod_fine']),np.nanmean(aod_aeronet['Walvis_Bay_airport']['aod_tot'])


# ## Now plot AERONET values on the map of ACAOD 501

# In[84]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax = ax.flatten()
ax1 = ax[0]
m = mapfig(ax=ax1)
mxe,mye = m(ye,xe)
p = ax1.pcolor(mxe,mye,a,vmin=0.0,vmax=1.0,cmap='plasma')

mxxa,myya = m(np.array(anet_lon),np.array(anet_lat))
pa = ax1.scatter(mxxa,myya,150,c=np.array(anet_mean_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='plasma',vmin=0.0,vmax=1.0)

ax1.set_title('Mean')
pu.sub_note('b)')
cb = m.colorbar(p,extend='max')


ax2 = ax[1]
m2 = mapfig(ax=ax2)
mxxe,myye = m2(np.array(yy)[io[1]],np.array(xx)[io[0]])
p2 = ax2.scatter(mxxe,myye,cnt.data[io[0],io[1]].flatten()/5.0,
                c=astd.data[io[0],io[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0.0)
m2xxa,m2yya = m2(np.array(anet_lon),np.array(anet_lat))
pa2 = ax2.scatter(m2xxa,m2yya,np.array(anet_cnt_fineaod)/10.0,c=np.array(anet_std_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='viridis',vmin=0.0,vmax=p2.get_clim()[1])

ax2.set_title('STD')
pu.sub_note('a)')
cb = m2.colorbar(p2,extend='max')

sizes = [10,100,500,1500]
labels = ['{0}'.format(z) for z in sizes]
points = [ax2.scatter([], [], s=z/5.0, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')
    
plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.03, right=0.97, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_statsAOD_2panel_actualmap_withAERONET.png',
            transparent=True,dpi=500)


# In[494]:


p2.get_clim()


# ## number of days map plot with AERONET

# In[205]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax = ax.flatten()
ax1 = ax[0]
m = mapfig(ax=ax1)
mxe,mye = m(ye,xe)
p = ax1.pcolor(mxe,mye,a,vmin=0.0,vmax=1.0,cmap='plasma')

mxxa,myya = m(np.array(anet_lon),np.array(anet_lat))
pa = ax1.scatter(mxxa,myya,150,c=np.array(anet_mean_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='plasma',vmin=0.0,vmax=1.0)

ax1.set_title('Mean')
pu.sub_note('b)')
cb = m.colorbar(p,extend='max')
cb.set_label('Mean of AOD 501 nm')


ax2 = ax[1]
m2 = mapfig(ax=ax2)
mxxe,myye = m2(np.array(yy)[io[1]],np.array(xx)[io[0]])
p2 = ax2.scatter(mxxe,myye,dcnt.data[io[0],io[1]].flatten()*30.0,
                c=astd.data[io[0],io[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0.0)
m2xxa,m2yya = m2(np.array(anet_lon),np.array(anet_lat))
pa2 = ax2.scatter(m2xxa,m2yya,150,c=np.array(anet_std_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='viridis',vmin=0.0,vmax=p2.get_clim()[1])

ax2.set_title('STD')
pu.sub_note('a)')
cb = m2.colorbar(p2,extend='max')
cb.set_label('Standard deviation of AOD 501 nm')

sizes = [1,2,3,5]#[10,100,500,1500]
labels = ['{0} day{1}'.format(z,'s' if z>1 else '') for z in sizes]
points = [ax2.scatter([], [], s=z*30.0, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left',handletextpad=0.1)
    
#plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.03, right=0.94, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_2panel_actualmap_withAERONET_days.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_2panel_actualmap_withAERONET_days.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_2panel_actualmap_withAERONET_days.eps')


# In[91]:


iof = np.where((stata05uc['dcnt'].data>0.0) & (stata05uc['std'].data<1.0))


# In[134]:


fig,ax = plt.subplots(1,3,figsize=(15,5))
ax = ax.flatten()
ax1 = ax[0]
m = mapfig(ax=ax1)
mxe,mye = m(ye,xe)
p = ax1.pcolor(mxe,mye,a,vmin=0.0,vmax=1.0,cmap='plasma')

mxxa,myya = m(np.array(anet_lon),np.array(anet_lat))
pa = ax1.scatter(mxxa,myya,150,c=np.array(anet_mean_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='plasma',vmin=0.0,vmax=1.0)

ax1.set_title('Mean')
pu.sub_note('c)')
cb = m.colorbar(p,extend='max')
cb.set_label('Mean of ACAOD 501 nm')


ax2 = ax[1]
m2 = mapfig(ax=ax2,no_ylabel=True)
mxxe,myye = m2(np.array(yy)[io[1]],np.array(xx)[io[0]])
p2 = ax2.scatter(mxxe,myye,dcnt.data[io[0],io[1]].flatten()*30.0,
                c=astd.data[io[0],io[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0.0)
m2xxa,m2yya = m2(np.array(anet_lon),np.array(anet_lat))
pa2 = ax2.scatter(m2xxa,m2yya,150,c=np.array(anet_std_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='viridis',vmin=0.0,vmax=p2.get_clim()[1])

ax2.set_title('STD')
ax2.get_yaxis().set_visible(False)
pu.sub_note('a)')
cb = m2.colorbar(p2,extend='max')
cb.set_label('Standard deviation of AOD 501 nm')

sizes = [1,2,3,5]#[10,100,500,1500]
labels = ['{0} day{1}'.format(z,'s' if z>1 else '') for z in sizes]
points = [ax2.scatter([], [], s=z*30.0, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left',handletextpad=0.1)


ax3 = ax[2]
m3 = mapfig(ax=ax3,no_ylabel=True)
mxxe3,myye3 = m3(np.array(stata05uc['y'])[iof[1]],np.array(stata05uc['x'])[iof[0]])
p3 = ax3.scatter(mxxe3,myye3,stata05uc['dcnt'].data[iof[0],iof[1]].flatten()*30.0,
                c=stata05uc['mean'].data[iof[0],iof[1]].flatten(),
               marker='s',edgecolor='None',cmap='plasma',vmin=0.0,vmax=1.0)
m3xxa,m3yya = m3(np.array(anet_lon),np.array(anet_lat))
pa3 = ax3.scatter(m3xxa,m3yya,150,c=np.array(anet_mean_totaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='plasma',vmin=0.0,vmax=p3.get_clim()[1])
ax3.set_title('Full column')
ax3.get_yaxis().set_visible(False)
ax3.axes.get_yaxis().set_ticklabels([])
pu.sub_note('b)')
cb = m3.colorbar(p3,extend='max')
cb.set_label('Full Column AOD 501 nm')

sizes = [1,2,3,5]#[10,100,500,1500]
labels = ['{0} day{1}'.format(z,'s' if z>1 else '') for z in sizes]
points = [ax3.scatter([], [], s=z*30.0, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left',handletextpad=0.1)
    
#plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.03, right=0.94, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_3panel_actualmap_withAERONET_days_full.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_3panel_actualmap_withAERONET_days_full.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_3panel_actualmap_withAERONET_days_full.eps')


# In[ ]:


p2 = ax2.scatter(mxxe,myye,cnt.data[io[0],io[1]].flatten()/5.0,
                c=astd.data[io[0],io[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0.0)
m2xxa,m2yya = m2(np.array(anet_lon),np.array(anet_lat))
pa2 = ax2.scatter(m2xxa,m2yya,np.array(anet_cnt_fineaod)/10.0,c=np.array(anet_std_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='viridis',vmin=0.0,vmax=p2.get_clim()[1])


# In[102]:


fig,ax = plt.subplots(1,3,figsize=(15,5))
ax = ax.flatten()
ax1 = ax[0]
m = mapfig(ax=ax1)
mxe,mye = m(ye,xe)
p = ax1.pcolor(mxe,mye,a,vmin=0.0,vmax=1.0,cmap='plasma')

mxxa,myya = m(np.array(anet_lon),np.array(anet_lat))
pa = ax1.scatter(mxxa,myya,150,c=np.array(anet_mean_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='plasma',vmin=0.0,vmax=1.0,zorder=10)


ax1.scatter(mxxe,myye,cnt.data[io[0],io[1]].flatten()/5.0,
                c=a.data[io[0],io[1]].flatten(),alpha=0.5,
               marker='o',edgecolor='k',linewidth=0.5,cmap='plasma',vmin=0.0,vmax=p.get_clim()[1])
m2gxxa,m2gyya = m(np.array(anet_lon),np.array(anet_lat))
ax1.scatter(m2gxxa,m2gyya,np.array(anet_cnt_fineaod)/10.0,c=np.array(anet_mean_fineaod),alpha=0.5,
               marker='o',edgecolor='k',linewidth=0.5,cmap='plasma',vmin=0.0,vmax=pa.get_clim()[1],zorder=-5)

sizes = [10,100,500,1500,3000]
labels = ['N={0}'.format(z) for z in sizes]
points = [ax1.scatter([], [], s=z/5.0, c='white',marker='o',edgecolor='k',linewidth=0.5,alpha=0.75) for z in sizes]
ax1.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')

ax1.set_title('Mean')
pu.sub_note('c)')
cb = m.colorbar(p,extend='max')
cb.set_label('Mean of ACAOD 501 nm')


ax2 = ax[1]
m2 = mapfig(ax=ax2,no_ylabel=True)
mxxe,myye = m2(np.array(yy)[io[1]],np.array(xx)[io[0]])
p2 = ax2.scatter(mxxe,myye,dcnt.data[io[0],io[1]].flatten()*30.0,
                c=astd.data[io[0],io[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis',vmin=0.0)
m2xxa,m2yya = m2(np.array(anet_lon),np.array(anet_lat))
pa2 = ax2.scatter(m2xxa,m2yya,150,c=np.array(anet_std_fineaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='viridis',vmin=0.0,vmax=p2.get_clim()[1])

ax2.set_title('STD')
ax2.get_yaxis().set_visible(False)
pu.sub_note('a)')
cb = m2.colorbar(p2,extend='max')
cb.set_label('Standard deviation of AOD 501 nm')

sizes = [1,2,3,5]#[10,100,500,1500]
labels = ['{0} day{1}'.format(z,'s' if z>1 else '') for z in sizes]
points = [ax2.scatter([], [], s=z*30.0, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left',handletextpad=0.1)


ax3 = ax[2]
m3 = mapfig(ax=ax3,no_ylabel=True)
mxxe3,myye3 = m3(np.array(stata05uc['y'])[iof[1]],np.array(stata05uc['x'])[iof[0]])
p3 = ax3.scatter(mxxe3,myye3,stata05uc['dcnt'].data[iof[0],iof[1]].flatten()*30.0,
                c=stata05uc['mean'].data[iof[0],iof[1]].flatten(),
               marker='s',edgecolor='None',cmap='plasma',vmin=0.0,vmax=1.0)
m3xxa,m3yya = m3(np.array(anet_lon),np.array(anet_lat))
pa3 = ax3.scatter(m3xxa,m3yya,150,c=np.array(anet_mean_totaod),
               marker='^',edgecolor='k',linewidth=0.5,cmap='plasma',vmin=0.0,vmax=p3.get_clim()[1])
ax3.set_title('Full column')
ax3.get_yaxis().set_visible(False)
ax3.axes.get_yaxis().set_ticklabels([])
pu.sub_note('b)')
cb = m3.colorbar(p3,extend='max')
cb.set_label('Full Column AOD 501 nm')

sizes = [1,2,3,5]#[10,100,500,1500]
labels = ['{0} day{1}'.format(z,'s' if z>1 else '') for z in sizes]
points = [ax3.scatter([], [], s=z*30.0, c='grey',marker='s',edgecolor='None') for z in sizes]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left',handletextpad=0.1)
    
#plt.suptitle('Above cloud AOD 501 nm\n')
plt.subplots_adjust(left=0.03, right=0.94, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_3panel_actualmap_withAERONET_days_full_num.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_3panel_actualmap_withAERONET_days_full_num.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_statsAOD_3panel_actualmap_withAERONET_days_full_num.eps')


# In[135]:


np.nanmean(a.data[iof[0],iof[1]]),np.nanmean(stata05uc['mean'].data[iof[0],iof[1]])


# In[137]:


np.nanmean(anet_mean_fineaod), np.nanmean(anet_mean_totaod)


# In[170]:


a[-8,:]


# ## Get the AERONET AOD on the routine flight comparison plot

# In[111]:


from datetime import datetime


# In[207]:


aod_aeronet['DRAGON_Henties_6']


# In[112]:


doy_rtn = [datetime.strptime(dnt,'%Y%m%d').timetuple().tm_yday for dnt in d_rtn]
doy_rtn
iaero_rtn = [(aod_aeronet['DRAGON_Henties_6']['doy']>=doyr*1.0-0.33) & (aod_aeronet['DRAGON_Henties_6']['doy']<=doyr*1.0+1.0) for doyr in doy_rtn]


# In[113]:


iadr = iaero_rtn[0]
aero_mean = []
aero_median = []
aod_aeronet['DRAGON_Henties_6']['aod_fine'][aod_aeronet['DRAGON_Henties_6']['aod_fine']<-998.0] = np.nan
for iai,iaer in enumerate(iaero_rtn):
    iadr  = iadr | iaer
    aero_mean.append(np.nanmean(aod_aeronet['DRAGON_Henties_6']['aod_fine'][iaer]))
    aero_median.append(np.nanmedian(aod_aeronet['DRAGON_Henties_6']['aod_fine'][iaer]))


# In[114]:


datetime(2016,9,1).timetuple().tm_yday


# In[115]:


iaero_aug = aod_aeronet['DRAGON_Henties_6']['doy']<245.0
iaero_sep = aod_aeronet['DRAGON_Henties_6']['doy']>=245.0


# In[116]:


aero_mean


# In[118]:


plt.figure(figsize=(12,8))

ax = plt.subplot(3,1,1)
plt.ylabel('AOD 500 nm')
plt.ylim(0,0.8)

plt.title('AOD above clouds along 2016 routine flight paths')
pu.sub_note('a)')

means = []
medians = []
for j,f in enumerate(flra):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    bo = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)
    color_box(bo,cls[j])
    [plt.setp(bo['fliers'][idx],alpha=0.05)for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],color=cls[j])for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['fliers'][idx],marker='.')for idx in xrange(len(bo['fliers']))]
    [plt.setp(bo['means'][idx],alpha=0.05)for idx in xrange(len(bo['means']))]
    means.append([a.get_ydata()[0] for a in bo['means']])
    medians.append([a.get_ydata()[0] for a in bo['medians']])
    plt.plot(pos3,[a.get_ydata()[0] for a in bo['means']],
             's',zorder=100,color=cls[j],label='{}/{}'.format(d_rtn[j][4:6],d_rtn[j][6:8],vv=vv),
             lw=2.5,alpha=0.4)   
    
meansr = []
mediansr = []
for j,f in enumerate(flr):
    binsf = []
    for i,c in enumerate(lims3[0:-1]):
        lon_fl = (s['Longitude'][f]>=c)&(s['Longitude'][f]<lims3[i+1])
        binsf.append(s['AOD0501'][f][lon_fl])
    bn = plt.boxplot(binsf,0,'.',showmeans=True,positions=pos3)

    meansr.append([a.get_ydata()[0] for a in bn['means']])
    mediansr.append([a.get_ydata()[0] for a in bn['medians']])
    [plt.setp(bn['fliers'][idx],alpha=0.0)for idx in xrange(len(bn['fliers']))]
    [plt.setp(bn['means'][idx],alpha=0.0)for idx in xrange(len(bn['means']))]
    [plt.setp(bn['boxes'][idx],alpha=0.0)for idx in xrange(len(bn['boxes']))]
    [plt.setp(bn['medians'][idx],alpha=0.0)for idx in xrange(len(bn['medians']))]
    [plt.setp(bn['whiskers'][idx],alpha=0.0)for idx in xrange(len(bn['whiskers']))]
    [plt.setp(bn['caps'][idx],alpha=0.0)for idx in xrange(len(bn['caps']))]

ac = []
for i,a in enumerate([daac[j] for j in i_flt]):
    plt.plot(a['BinCenter_longitude'][0,:],a['meanAODperbin'][1:,0],'x-',lw=2,color=cls[i],alpha=0.4)
             #label='{}/{} MODIS AAC'.format(d_rtn[i][4:6],d_rtn[i][6:8]))
    ac.append(a['meanAODperbin'][1:,0])
    
for i,aerm in enumerate(aero_mean):
    plt.plot(14.0,aerm,'o',color=cls[i],alpha=0.4,lw=2,markersize=15)

plt.plot([],[],'ko',alpha=0,label=' ')
plt.plot([],[],'ks',alpha=0.4,label='4STAR')
plt.plot([],[],'kx-',alpha=0.4,label='MOD06ACAERO')
plt.plot([],[],'ko',alpha=0.4,label='AERONET',markersize=15)
plt.plot([],[],'ko',alpha=0,label=' ')
plt.legend(numpoints=1,frameon=True, ncol=2, bbox_to_anchor=(1.005,1.0),loc=2,handletextpad=0.2,columnspacing=0.7)


ax2 = plt.subplot(3,1,2,sharex=ax)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AODFM_CLIMOMEAN'].data[0,:],
         '*-',color='r',label='MODIS Fine AOD (Sept. 2001-2013)',zorder=50,lw=2.5)
plt.plot(m.variables['LONGITUDE'].data[0,:],m.variables['AOD_CLIMOMEAN'].data[0,:],
         'v-',color='b',label='MODIS Total AOD (Sept. 2001-2013)',zorder=51,lw=0.5)
pu.sub_note('b)')
ac_aug = []
for i,a in enumerate([daac[j] for j in i_aug]):
    ac_aug.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter_longitude'][0,:],np.nanmean(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MOD06ACAERO for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter_longitude'][0,:],np.nanmean(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MOD06ACAERO for all Sept. 2016')


plt.plot(a['BinCenter_longitude'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MOD06ACAERO for routine flights')

plt.plot(14.0,np.nanmean(aod_aeronet['DRAGON_Henties_6']['aod_fine'][iaero_aug]),'o',color='fuchsia',markersize=15,alpha=0.4,
         label='AERONET for Aug. 2016')
plt.plot(14.0,np.nanmean(aod_aeronet['DRAGON_Henties_6']['aod_fine'][iaero_sep]),'o',color='aqua',markersize=15,alpha=0.4,
         label='AERONET for Sept. 2016')
plt.plot(14.0,np.nanmean(aero_mean),'o',markersize=15,color='maroon',label='AERONET for routine flights',alpha=0.4)

plt.plot(pos3,np.nanmean(np.array(meansr),axis=0),'x-',color='grey',lw=0.5,zorder=180,label='4STAR [0.5-1.6 km]')
plt.plot(pos3,np.nanmean(np.array(means),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR AAC')
plt.ylabel('Mean AOD 500 nm')


ax3 = plt.subplot(3,1,3,sharex=ax)
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AODFM_YRMEAN'].data[0,:,:],axis=1),
         '*-',color='r',zorder=50,lw=2.5,label='MODIS Fine AOD (Sept. 2001-2013)')
plt.plot(m2.variables['LONGITUDE'].data[0,:],np.median(m2.variables['AOD_YRMEAN'].data[0,:,:],axis=1),
         'v-',color='b',zorder=51,lw=0.5,label='MODIS Total AOD (Sept. 2001-2013)')

plt.ylim(0,0.8)

ac_aug = []
for i,a in enumerate([daac[j] for j in i_aug]):
    ac_aug.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter_longitude'][0,:],np.nanmedian(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MOD06ACAERO for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter_longitude'][0,:],np.nanmedian(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MOD06ACAERO for all Sept. 2016')
plt.plot(a['BinCenter_longitude'][0,:],np.nanmedian(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MOD06ACAERO for routine flights')

plt.plot(14.0,np.nanmedian(aod_aeronet['DRAGON_Henties_6']['aod_fine'][iaero_aug]),'o',color='fuchsia',markersize=15,alpha=0.4,
         label='AERONET for Aug. 2016')
plt.plot(14.0,np.nanmedian(aod_aeronet['DRAGON_Henties_6']['aod_fine'][iaero_sep]),'o',color='aqua',markersize=15,alpha=0.4,
         label='AERONET for Sept. 2016')
plt.plot(14.0,np.nanmedian(aero_median),'o',markersize=15,color='maroon',label='AERONET for routine flights',alpha=0.4)

plt.plot(pos3,np.nanmedian(np.array(mediansr),axis=0),'x-',color='grey',lw=0.5,zorder=180,label='4STAR [0.5-1.6 km]')
plt.plot(pos3,np.nanmedian(np.array(medians),axis=0),'s-k',lw=3.5,zorder=200,label='4STAR ACAOD')
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('Median AOD 500 nm')
#plt.title('Median')
pu.sub_note('c)')

plt.legend(numpoints=1,frameon=True, bbox_to_anchor=(1.005,1.95),loc=2)
ax3.set_xlim(0,14.4)
ti = ax3.set_xticks([0,2,4,6,8,10,12,14])
tl = ax3.set_xticklabels([0,2,4,6,8,10,12,14])

ti = ax2.set_xticks([0,2,4,6,8,10,12,14])
tl = ax2.set_xticklabels([0,2,4,6,8,10,12,14])

ti = ax.set_xticks([0,2,4,6,8,10,12,14])
tl = ax.set_xticklabels([0,2,4,6,8,10,12,14])

box = ax.get_position()
ax.set_position([box.x0-0.06, box.y0, box.width * 0.7, box.height])
box2 = ax2.get_position()
ax2.set_position([box2.x0-0.06, box2.y0, box2.width * 0.7, box2.height])
box3 = ax3.get_position()
ax3.set_position([box3.x0-0.06, box3.y0, box3.width * 0.7, box3.height])


plt.savefig(fp+'plot_v3/ORACLES2016_MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_3_split_AERONET.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_3_split_AERONET.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_3_split_AERONET.eps')


# # Spectral AOD mean, median, pca, std...

# In[29]:


aods.shape


# In[225]:


s['fl_acaod']


# In[226]:


ns = len(s['Start_UTC'][s['fl_acaod']])


# In[227]:


pbar = tqdm(total=ns)


# In[152]:


plt.figure()
color.cycle_cmap(ns,cmap=plt.cm.gist_ncar,ax=plt.gca())

for i,u in enumerate(s['Start_UTC']):
    if s['fl_acaod'][i]:
        plt.plot(wvl,aods[i,:],'-x',lw=0.2)
        pbar.update(1)

plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.ylim([0,1.2])
plt.title('ACAOD')
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar)
scalarmap.set_array(s['Start_UTC'][s['fl_acaod']])
cba = plt.colorbar(scalarmap)
cba.set_label('UTC [h]')
#plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_all_ACAOD_spectra.png',
#            transparent=True,dpi=500)


# In[153]:


plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_all_ACAOD_spectra.png',
            transparent=True,dpi=500)


# ## ACAOD spectra, mean, median, std

# In[44]:


iw = wvl!=700.0


# In[45]:


meanaod = np.nanmean(aods[s['fl_acaod'],:],axis=0)
medianaod = np.nanmedian(aods[s['fl_acaod'],:],axis=0)
stdaod = np.nanstd(aods[s['fl_acaod'],:],axis=0)


# In[46]:


meanuncaod = np.nanmean(uncaods[s['fl_acaod'],:],axis=0)


# In[250]:


plt.figure()
plt.plot(wvl[iw],meanaod[iw],'k-+',label='Mean')
plt.errorbar(wvl[iw],meanaod[iw],yerr=meanuncaod[iw],color='k')
plt.plot(wvl[iw],meanaod[iw]+stdaod[iw],'--x',color='grey',lw=0.4,label='Mean +/- 1 Standard Deviation')
plt.plot(wvl[iw],stdaod[iw],'-^',lw=0.4,color='lightcoral',label='Standard Deviation',zorder=-2,alpha=0.5)
plt.plot(wvl[iw],meanaod[iw]-stdaod[iw],'--x',color='grey',lw=0.4)
plt.plot(wvl[iw],medianaod[iw],'-o',lw=2.0,color='lightblue',label='Median',zorder=-1)
plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.title('ACAOD average spectra from 4STAR at selected wavelengths')
plt.grid()

plt.legend(frameon=False)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_mean_ACAOD_spectra_less700.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_mean_ACAOD_spectra_less700.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_mean_ACAOD_spectra_less700.eps')


# In[47]:


mean_angs = -1.0*(np.log(meanaod[je])-np.log(meanaod[ja]))/(np.log(wvl[je])-np.log(wvl[ja]))
mean_angs


# In[48]:


wvl[[2,15]]


# In[49]:


meanaod[[2,15]], medianaod[[2,15]],stdaod[[2,15]],meanuncaod[[2,15]]


# ## Log of ACAOD spectra

# In[64]:


fig = plt.figure()
plt.plot(wvl[iw],meanaod[iw],'k-+',label='Mean')
plt.errorbar(wvl[iw],meanaod[iw],yerr=meanuncaod[iw],color='k')
plt.plot(wvl[iw],meanaod[iw]+stdaod[iw],'--x',color='grey',lw=0.4,label='Mean +/- 1 Standard Deviation')
plt.plot(wvl[iw],stdaod[iw],'-^',lw=0.4,color='lightcoral',label='Standard Deviation',zorder=-2,alpha=0.5)
plt.plot(wvl[iw],meanaod[iw]-stdaod[iw],'--x',color='grey',lw=0.4)
plt.plot(wvl[iw],medianaod[iw],'-o',lw=2.0,color='lightblue',label='Median',zorder=-1)
plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.title('ACAOD average spectra from 4STAR at selected wavelengths')
plt.grid()

ax = plt.gca()

ax.set_yscale('log')
ax.set_xscale('log')

plt.xlim(340,1700)
plt.ylim(0.007,1)
plt.xticks([350,400,500,600,700,800,900,1000,1200,1400,1600])
ax.set_xticklabels([350,400,500,600,700,800,900,1000,1200,1400,1600])
plt.yticks([0.01,0.02,0.05,0.1,0.15,0.2,0.3,0.5,0.8,1])
ax.set_yticklabels([0.01,0.02,0.05,0.1,0.15,0.2,0.3,0.5,0.8,1])

plt.legend(frameon=False,loc=3)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_mean_ACAOD_spectra_less700_loglog.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_mean_ACAOD_spectra_less700_loglog.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_mean_ACAOD_spectra_less700_loglog.eps')


# ## Get the Principal components from PCA

# In[238]:


ivalid = np.all(np.isfinite(aods[s['fl_acaod'],:]),axis=1)
aods_valid = aods[s['fl_acaod'],:][ivalid,:]


# In[239]:


from sklearn.decomposition import PCA
pca = PCA(n_components=20)
pca.fit(aods_valid)


# In[240]:


pca.explained_variance_ratio_[0]


# In[241]:


pca.explained_variance_ratio_[:]


# In[242]:


pca.components_.shape


# In[243]:


plt.figure()
plt.plot(wvl,meanaod,'k-+',label='Mean')
plt.errorbar(wvl,meanaod,yerr=meanuncaod,color='k')
plt.plot(wvl,meanaod+stdaod,'--x',color='grey',lw=0.4,label='Mean +/- 1 Standard Deviation')
plt.plot(wvl,meanaod-stdaod,'--x',color='grey',lw=0.4)
plt.plot(wvl,stdaod,'-s',color='orange',lw=0.4,label='Standard Deviation')
plt.plot(wvl,medianaod,'-o',lw=2.0,color='lightblue',label='Median',zorder=-1)
plt.plot(wvl,pca.components_[0,:],'-*',lw=0.8,color='coral',zorder=-2,
         label='Principal Component\nrepresenting {:2.1f}% of variance'.format(pca.explained_variance_ratio_[0]*100.0))
plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.title('ACAOD representative spectra')

plt.legend(frameon=False)
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_mean_pca_ACAOD_spectra.png',
            transparent=True,dpi=500)


# In[244]:


plt.figure()
plt.plot(wvl[iw],pca.components_[0,iw])
plt.plot(wvl[iw],pca.components_[1,iw])
plt.plot(wvl[iw],pca.components_[2,iw])
plt.plot(wvl[iw],pca.components_[3,iw])
plt.plot(wvl[iw],pca.components_[4,iw])


# In[245]:


plt.figure()
plt.plot(wvl[iw],meanaod[iw],'k-+',label='Mean')
plt.errorbar(wvl[iw],meanaod[iw],yerr=meanuncaod[iw],color='k')
plt.plot(wvl[iw],meanaod[iw]+stdaod[iw],'--x',color='grey',lw=0.4,label='Mean +/- 1 Standard Deviation')
plt.plot(wvl[iw],meanaod[iw]-stdaod[iw],'--x',color='grey',lw=0.4)
plt.plot(wvl[iw],medianaod[iw],'-o',lw=2.0,color='lightblue',label='Median',zorder=-1)
plt.plot(wvl[iw],pca.components_[0,iw],'-*',lw=0.8,color='coral',zorder=-2,
         label='Principal Component\nrepresenting {:2.1f}% of variance'.format(pca.explained_variance_ratio_[0]*100.0))
plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.title('ACAOD representative spectra')
plt.grid()


plt.legend(frameon=False)
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_mean_pca_ACAOD_spectra_no700.png',
            transparent=True,dpi=500)


# In[574]:


plt.figure()

plt.plot(wvl[iw],meanaod[iw],'k-+',label='Mean')
plt.errorbar(wvl[iw],meanaod[iw],yerr=meanuncaod[iw],color='k')
plt.plot(wvl[iw],meanaod[iw]+stdaod[iw],'--x',color='grey',lw=0.4,label='Mean +/- 1 Standard Deviation')
plt.plot(wvl[iw],meanaod[iw]-stdaod[iw],'--x',color='grey',lw=0.4)
plt.plot(wvl[iw],medianaod[iw],'-o',lw=2.0,color='lightblue',label='Median',zorder=-1)


plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.title('ACAOD representative spectra')
plt.grid()
plt.plot([1000,1000],[0,0],'-*',lw=0.8,color='coral',
        label='Principal Component\nrepresenting {:2.1f}% of variance'.format(pca.explained_variance_ratio_[0]*100.0))
plt.legend(frameon=False)


ax2 = plt.twinx(plt.gca())
ax2.plot(wvl[iw],pca.components_[0,iw],'-*',lw=0.8,color='coral',zorder=-2,
         label='Principal Component\nrepresenting {:2.1f}% of variance'.format(pca.explained_variance_ratio_[0]*100.0))

ax2.set_ylabel('Principal Component Loading',color='coral')


plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_mean_pca_ACAOD_spectra_no700_twinax.png',
            transparent=True,dpi=500)


# In[254]:


plt.figure()
plt.plot(wvl,meanaod,'k-+',label='Mean')
plt.errorbar(wvl,meanaod,yerr=meanuncaod,color='k')
plt.plot(wvl,meanaod+stdaod,'--x',color='grey',lw=0.4,label='Mean +/- 1 Standard Deviation')
plt.plot(wvl,meanaod-stdaod,'--x',color='grey',lw=0.4)
plt.plot(wvl,medianaod,'-o',lw=2.0,color='lightblue',label='Median',zorder=-1)
plt.plot(wvl,pca.components_[0,:],'-*',lw=0.8,color='coral',zorder=-2,
         label='Principal Component\nrepresenting {:2.1f}% of variance'.format(pca.explained_variance_ratio_[0]*100.0))
plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.title('ACAOD representative spectra')
plt.yscale('log')
plt.xscale('log')
plt.xlim(300,1750)
plt.xticks([350,400,500,660,800,900,1000,1200,1400,1650])
#plt.xticklabel([350,400,500,660,800,900,1000,1200,1400,1650])
plt.yticks([0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
#plt.yticklabel([0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
plt.grid()

plt.legend(frameon=False)
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_mean_pca_ACAOD_spectra_loglog.png',
            transparent=True,dpi=500)


# In[551]:


pca_angs = -1.0*(np.log(pca.components_[0,je])-np.log(pca.components_[0,ja]))/(np.log(wvl[je])-np.log(wvl[ja]))
pca_angs


# # Get the spacing and altitudes of gap from ACAOD altitudes and distances

# In[29]:


ii_flacaod = np.where(s['fl_acaod_noQA'])[0] 
ii_flacaod[0]


# In[30]:


ii_flacaod_0_5 = np.where(s['fl_acaod_noQA'] & (s['Longitude']>=0) & (s['Longitude']<5.0))[0] 
ii_flacaod_0_5[0]


# In[31]:


ii_flacaod_5_10 = np.where(s['fl_acaod_noQA'] & (s['Longitude']>=5) & (s['Longitude']<10.0))[0] 
ii_flacaod_5_10[0]


# In[32]:


ii_flacaod_10_15 = np.where(s['fl_acaod_noQA'] & (s['Longitude']>=10) & (s['Longitude']<15.0))[0] 
ii_flacaod_10_15[0]


# In[33]:


ii_flacaod_10_14 = np.where(s['fl_acaod_noQA'] & (s['Longitude']>=10) & (s['Longitude']<14.0))[0] 
ii_flacaod_10_14[0]


# In[34]:


disc_flacaod = np.where(np.diff(ii_flacaod,1)>1)[0]
disc_flacaod_long = np.where(np.diff(ii_flacaod,1)>150)[0]
disc_flacaod_0_5 = np.where(np.diff(ii_flacaod_0_5,1)>150)[0]
disc_flacaod_5_10 = np.where(np.diff(ii_flacaod_5_10,1)>150)[0]
disc_flacaod_10_15 = np.where(np.diff(ii_flacaod_10_15,1)>150)[0]
disc_flacaod_10_14 = np.where(np.diff(ii_flacaod_10_14,1)>150)[0]


# In[25]:


plt.figure()
plt.plot(ii_flacaod[0:100],'.')
plt.plot(np.arange(100)[disc_flacaod[0:6]],ii_flacaod[0:100][disc_flacaod[0:6]],'rx')
plt.plot(np.arange(100)[disc_flacaod[0:6]+1],ii_flacaod[0:100][disc_flacaod[0:6]+1],'go')


# In[35]:


discontinuity_istart =  ii_flacaod[np.append(0,disc_flacaod[:-1]+1)]
discontinuity_iend =  ii_flacaod[disc_flacaod]


# In[36]:


discontinuity_istart_long =  ii_flacaod[np.append(0,disc_flacaod_long[:-1]+1)]
discontinuity_iend_long =  ii_flacaod[disc_flacaod_long]


# In[37]:


discontinuity_istart_0_5 =  ii_flacaod_0_5[np.append(0,disc_flacaod_0_5[:-1]+1)]
discontinuity_iend_0_5 =  ii_flacaod_0_5[disc_flacaod_0_5]
discontinuity_istart_5_10 =  ii_flacaod_5_10[np.append(0,disc_flacaod_5_10[:-1]+1)]
discontinuity_iend_5_10 =  ii_flacaod_5_10[disc_flacaod_5_10]
discontinuity_istart_10_15 =  ii_flacaod_10_15[np.append(0,disc_flacaod_10_15[:-1]+1)]
discontinuity_iend_10_15 =  ii_flacaod_10_15[disc_flacaod_10_15]
discontinuity_istart_10_14 =  ii_flacaod_10_14[np.append(0,disc_flacaod_10_14[:-1]+1)]
discontinuity_iend_10_14 =  ii_flacaod_10_14[disc_flacaod_10_14]


# In[38]:


delta_alt,delta_lon,delta_lat = [],[],[]
for i,start in enumerate(discontinuity_istart):
    try:
        ma = np.nanmax(s['GPS_Alt'][start:discontinuity_iend[i]])
        mi = np.nanmin(s['GPS_Alt'][start:discontinuity_iend[i]])
        delta_alt.append(ma-mi)
        delta_lon.append(np.nanmean(s['Longitude'][start:discontinuity_iend[i]]))
        delta_lat.append(np.nanmean(s['Latitude'][start:discontinuity_iend[i]]))
    except:
        pass
delta_alt = np.array(delta_alt)
delta_lon = np.array(delta_lon)
delta_lat = np.array(delta_lat)


# In[39]:


ldelta_alt,ldelta_lon,ldelta_lat,ldelta_lon_days,ldelta_lat_days = [],[],[],[],[]
for i,start in enumerate(discontinuity_istart_long):
    try:
        ma = np.nanmax(s['GPS_Alt'][start:discontinuity_iend_long[i]])
        mi = np.nanmin(s['GPS_Alt'][start:discontinuity_iend_long[i]])
        ldelta_alt.append(ma-mi)
        ldelta_lon.append(np.nanmean(s['Longitude'][start:discontinuity_iend_long[i]]))
        ldelta_lat.append(np.nanmean(s['Latitude'][start:discontinuity_iend_long[i]]))
        ldelta_lat_days.append(np.unique(s['days'][start:discontinuity_iend_long[i]]))
        ldelta_lon_days.append(np.unique(s['days'][start:discontinuity_iend_long[i]]))
    except:
        pass
ldelta_alt = np.array(ldelta_alt)
ldelta_lon = np.array(ldelta_lon)
ldelta_lat = np.array(ldelta_lat)
ldelta_lat_days = np.array(ldelta_lat_days)
ldelta_lon_days = np.array(ldelta_lon_days)


# In[40]:


delta_alt1,delta_lon1,delta_lat1 = [],[],[]
for i,start in enumerate(discontinuity_istart_0_5):
    try:
        ma = np.nanmax(s['GPS_Alt'][start:discontinuity_iend_0_5[i]])
        mi = np.nanmin(s['GPS_Alt'][start:discontinuity_iend_0_5[i]])
        delta_alt1.append(ma-mi)
        delta_lon1.append(np.nanmean(s['Longitude'][start:discontinuity_iend_0_5[i]]))
        delta_lat1.append(np.nanmean(s['Latitude'][start:discontinuity_iend_0_5[i]]))
    except:
        pass
delta_alt1 = np.array(delta_alt1)
delta_lon1 = np.array(delta_lon1)
delta_lat1 = np.array(delta_lat1)

delta_alt2,delta_lon2,delta_lat2 = [],[],[]
for i,start in enumerate(discontinuity_istart_5_10):
    try:
        ma = np.nanmax(s['GPS_Alt'][start:discontinuity_iend_5_10[i]])
        mi = np.nanmin(s['GPS_Alt'][start:discontinuity_iend_5_10[i]])
        delta_alt2.append(ma-mi)
        delta_lon2.append(np.nanmean(s['Longitude'][start:discontinuity_iend_5_10[i]]))
        delta_lat2.append(np.nanmean(s['Latitude'][start:discontinuity_iend_5_10[i]]))
    except:
        pass
delta_alt2 = np.array(delta_alt2)
delta_lon2 = np.array(delta_lon2)
delta_lat2 = np.array(delta_lat2)

delta_alt3,delta_lon3,delta_lat3 = [],[],[]
for i,start in enumerate(discontinuity_istart_10_15):
    try:
        ma = np.nanmax(s['GPS_Alt'][start:discontinuity_iend_10_15[i]])
        mi = np.nanmin(s['GPS_Alt'][start:discontinuity_iend_10_15[i]])
        delta_alt3.append(ma-mi)
        delta_lon3.append(np.nanmean(s['Longitude'][start:discontinuity_iend_10_15[i]]))
        delta_lat3.append(np.nanmean(s['Latitude'][start:discontinuity_iend_10_15[i]]))
    except:
        pass
delta_alt3 = np.array(delta_alt3)
delta_lon3 = np.array(delta_lon3)
delta_lat3 = np.array(delta_lat3)

delta_alt3b,delta_lon3b,delta_lat3b = [],[],[]
for i,start in enumerate(discontinuity_istart_10_14):
    try:
        ma = np.nanmax(s['GPS_Alt'][start:discontinuity_iend_10_14[i]])
        mi = np.nanmin(s['GPS_Alt'][start:discontinuity_iend_10_14[i]])
        delta_alt3b.append(ma-mi)
        delta_lon3b.append(np.nanmean(s['Longitude'][start:discontinuity_iend_10_14[i]]))
        delta_lat3b.append(np.nanmean(s['Latitude'][start:discontinuity_iend_10_14[i]]))
    except:
        pass
delta_alt3b = np.array(delta_alt3b)
delta_lon3b = np.array(delta_lon3b)
delta_lat3b = np.array(delta_lat3b)


# In[41]:


np.nanmax(delta_alt),np.nanmax(ldelta_alt),np.nanmax(delta_alt1),np.nanmax(delta_alt2),np.nanmax(delta_alt3)


# In[59]:


float(len(ldelta_alt[ldelta_alt<360.0]))*100.0/len(ldelta_alt)


# In[84]:


delta_alt


# In[264]:


plt.figure()
plt.hist(delta_alt,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
#plt.hist(delta_alt_all,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
plt.yscale('log')
plt.xscale('log')
plt.axvline(np.nanmean(delta_alt_all),label='mean')
plt.axvline(np.nanmedian(delta_alt_all),label='median',ls='--')
plt.legend(frameon=False)


# In[265]:


plt.figure()
plt.hist(ldelta_alt,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
#plt.hist(delta_alt_all,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
#plt.yscale('log')
#plt.xscale('log')
plt.axvline(np.nanmean(ldelta_alt),label='mean')
plt.axvline(np.nanmedian(ldelta_alt),label='median',ls='--')
plt.legend(frameon=False)


# In[519]:


plt.figure()
#plt.hist(ldelta_alt,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3,label='All')
#plt.hist(delta_alt1,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3,label=u'0° - 5°E')
#plt.hist(delta_alt2,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3,label=u'5°E - 10°E')
#plt.hist(delta_alt3,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3,label=u'10°E - 15°E')
plt.hist([delta_alt1,delta_alt2,delta_alt3],bins=30,range=[0,2500],stacked=False,normed=False,edgecolor='None',alpha=0.3,
         label=[u'0° - 5°E',u'5°E - 10°E',u'10°E - 15°E'])

#plt.hist(delta_alt_all,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
#plt.yscale('log')
#plt.xscale('log')
plt.axvline(np.nanmean(ldelta_alt),label='mean')
plt.axvline(np.nanmedian(ldelta_alt),label='median',ls='--')
plt.legend(frameon=False)


# In[34]:


s['days']


# In[35]:


len(ldelta_alt)


# ## Plot the ACAOD altitudes

# In[36]:


np.min(s['GPS_Alt'][s['fl_acaod']][s['Longitude'][s['fl_acaod']]>13.0]), np.max(s['GPS_Alt'][s['fl_acaod']][s['Longitude'][s['fl_acaod']]>13.0])


# In[267]:


plt.figure()
plt.plot(s['Longitude'][s['fl_acaod']],s['GPS_Alt'][s['fl_acaod']],'.')
plt.xlim(13,16)


# ## make a box plot of the altitude differences in gap vs. longitude

# In[79]:


delta_lon


# In[37]:


numlon,binlon = np.histogram(delta_lon,range=(0,16))


# In[38]:


binlon


# In[39]:


bin_alt,bin_lon = [],[]
for i in xrange(16):
    fllon = (delta_lon>=i) & (delta_lon<(i+1.0))
    bin_alt.append(delta_alt[fllon])
    bin_lon.append(np.mean([i,i+1.0]))


# In[16]:


lbin_alt,lbin_lon,lbin_ndays,lbin_daysname,lbin_days,lbin_lat = [],[],[],[],[],[]
for i in xrange(16):
    fllon = (ldelta_lon>=i) & (ldelta_lon<(i+1.0))
    lbin_alt.append(ldelta_alt[fllon])
    lbin_lon.append(np.mean([i,(i+1.0)]))
    try:
        lbin_ndays.append(len(np.unique(np.hstack(ldelta_lon_days[fllon]))))
        lbin_daysname.append(np.unique(np.hstack(ldelta_lon_days[fllon])))
        lbin_days.append(np.hstack(ldelta_lon_days[fllon]))
        lbin_lat.append(np.hstack(ldelta_lat[fllon]))
    except:
        lbin_ndays.append(0)
        lbin_daysname.append('no')
        lbin_days.append('no')
        lbin_lat.append(-99)


# In[17]:


lbin_ndays


# In[20]:


for i in xrange(len(lbin_alt)):
    print lbin_lon[i], len(lbin_alt[i]),lbin_ndays[i]


# In[166]:


[(d,lbin_alt[7][i],lbin_lat[7][i]) for i,d in enumerate(lbin_days[7])]


# In[148]:


lbin_daysname


# In[84]:


bin_alt


# In[85]:


plt.figure()
plt.boxplot(bin_alt,positions=bin_lon)
plt.yscale('log')
plt.ylim(0,3700)
plt.ylabel('Gap distance [m]')
plt.xlabel(u'Longitude [°]')


# In[87]:


plt.figure()
plt.boxplot(lbin_alt,positions=lbin_lon)
plt.plot(ldelta_lon,ldelta_alt,'.',alpha=0.5)
#plt.yscale('log')
plt.ylim(-50,3700)
plt.ylabel('Gap distance [m]')
plt.xlabel(u'Longitude [°]')
plt.title('Distribution of gap altitude per longitude')
#plt.savefig()


# In[420]:


plt.figure(figsize=(12,5))
ax = plt.subplot2grid((1,4),(0,0),colspan=3)

ax.boxplot(lbin_alt,positions=lbin_lon,showfliers=False)
ax.plot(ldelta_lon,ldelta_alt,'k.',alpha=0.2,zorder=-1)
#plt.yscale('log')
ax.set_ylim(-50,3700)
ax.set_ylabel('Gap distance [m]')
ax.set_xlabel(u'Longitude [°]')
#plt.xticks(rotation=45,ha="right") 
ax.set_title('Distribution of gap altitude per longitude')

ax2 = plt.subplot2grid((1,4),(0,3),sharey=ax)

ax2.hist(ldelta_alt,bins=50,range=[0,4000],normed=False,edgecolor='None',alpha=0.3,orientation='horizontal')
ax2.axhline(np.nanmean(ldelta_alt),label='Mean')
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.axhline(np.nanmedian(ldelta_alt),label='Median',ls='--')
ax2.annotate('{:4.0f} m'.format(np.nanmean(ldelta_alt)),(15,np.nanmean(ldelta_alt)*1.05),color='b',alpha=0.3)
ax2.annotate('{:4.0f} m'.format(np.nanmedian(ldelta_alt)),(15,np.nanmedian(ldelta_alt)*1.05),color='b',alpha=0.3)

ax2.set_title('Gap distribution')
#plt.xticks(rotation=45,ha="right")
ax2.set_xlabel('Counts')
plt.legend(frameon=False,handletextpad=0.05)

#plt.tight_layout()
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_gap_distribution.png',
            transparent=True,dpi=500)


# In[299]:


plt.figure(figsize=(12,5))
ax = plt.subplot2grid((1,4),(0,0),colspan=3)

bow = ax.boxplot(lbin_alt,positions=lbin_lon,showfliers=True,sym='.')
ax.plot(ldelta_lon,ldelta_alt,'k.',alpha=0.01,zorder=-1)
for j,nn in enumerate(lbin_ndays):
    if nn>0:
        ax.text(lbin_lon[j],np.nanmax(bow['medians'][j].get_data()[1])+30,'{:2.0f}'.format(nn),color='blue',horizontalalignment='center')
    else:
        ax.text(lbin_lon[j],50,'{:2.0f}'.format(nn),color='blue',horizontalalignment='center')
#plt.yscale('log')
ax.set_ylim(-50,3700)
ax.set_ylabel('Vertical extent [m]')
ax.set_xlabel(u'Longitude [°]')
pu.sub_note('a)')
#plt.xticks(rotation=45,ha="right") 
ax.set_title('Distribution of gap vertical extent per longitude')

ax2 = plt.subplot2grid((1,12),(0,9),sharey=ax)
ax2.hist(delta_alt1,bins=50,range=[0,4000],normed=False,edgecolor='None',alpha=0.3,orientation='horizontal')
ax2.axhline(np.nanmean(delta_alt1),label='Mean')
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.axhline(np.nanmedian(delta_alt1),label='Median',ls='--')
ax2.annotate('{:4.0f} m'.format(np.nanmean(delta_alt1)),(1,np.nanmean(delta_alt1)*1.05),color='b',alpha=0.3)
ax2.annotate('{:4.0f} m'.format(np.nanmedian(delta_alt1)),(1,np.nanmedian(delta_alt1)*1.05),color='b',alpha=0.3)
ax2.annotate(' {:2.0f}%'.format(float(len(delta_alt1[delta_alt1<60.0]))*100.0/len(delta_alt1)),(1.4,2450.0),color='b')
ax2.annotate('({:2.0f}%)'.format(float(len(delta_alt1[delta_alt1<360.0]))*100.0/len(delta_alt1)),(1.4,2300.0),color='b')
ax2.set_title(u'0°-5°')
pu.sub_note('b)')
ax2.set_xticks([0,3,6])
#plt.xticks(rotation=45,ha="right")
#ax2.set_xlabel('Counts')
#plt.legend(frameon=False,handletextpad=0.05)


ax3 = plt.subplot2grid((1,12),(0,10),sharey=ax)
ax3.hist(delta_alt2,bins=50,range=[0,4000],normed=False,edgecolor='None',alpha=0.3,orientation='horizontal')
ax3.axhline(np.nanmean(delta_alt2),label='Mean')
plt.setp(ax3.get_yticklabels(), visible=False)
ax3.axhline(np.nanmedian(delta_alt2),label='Median',ls='--')
ax3.annotate('{:4.0f} m'.format(np.nanmean(delta_alt2)),(1.5,np.nanmean(delta_alt2)*1.05),color='b',alpha=0.3)
ax3.annotate('{:4.0f} m'.format(np.nanmedian(delta_alt2)),(1.5,np.nanmedian(delta_alt2)*1.05),color='b',alpha=0.3)
ax3.annotate(' {:2.0f}%'.format(float(len(delta_alt2[delta_alt2<60.0]))*100.0/len(delta_alt2)),(3.1,2450.0),color='b')
ax3.annotate('({:2.0f}%)'.format(float(len(delta_alt2[delta_alt2<360.0]))*100.0/len(delta_alt2)),(3.1,2300.0),color='b')
ax3.set_title('Gap extent distribution\n'+u'5°-10°')
#plt.xticks(rotation=45,ha="right")
ax3.set_xlabel('Counts')
pu.sub_note('c)')
ax3.set_xticks([0,6,12])
plt.legend(frameon=True,handletextpad=0.02,bbox_to_anchor=(1.01, 0.93))#,ncol=2,columnspacing=0.2,loc=1)

ax4 = plt.subplot2grid((1,12),(0,11),sharey=ax)
ax4.hist(delta_alt3,bins=50,range=[0,4000],normed=False,edgecolor='None',alpha=0.3,orientation='horizontal')
ax4.axhline(np.nanmean(delta_alt3),label='Mean')
plt.setp(ax4.get_yticklabels(), visible=False)
ax4.axhline(np.nanmedian(delta_alt3),label='Median',ls='--')
ax4.annotate('{:4.0f} m'.format(np.nanmean(delta_alt3)),(0.5,np.nanmean(delta_alt3)*1.05),color='b',alpha=0.3)
ax4.annotate('{:4.0f} m'.format(np.nanmedian(delta_alt3)),(0.5,np.nanmedian(delta_alt3)*1.05),color='b',alpha=0.3)
ax4.annotate(' {:2.0f}%'.format(float(len(delta_alt3[delta_alt3<60.0]))*100.0/len(delta_alt3)),(1.2,2450.0),color='b')
ax4.annotate('({:2.0f}%)'.format(float(len(delta_alt3[delta_alt3<360.0]))*100.0/len(delta_alt3)),(1.2,2300.0),color='b')
ax4.set_title(u'10°-15°')
pu.sub_note('d)')
ax4.set_xticks([0,3,6])
#ax4.legend(frameon=True,handletextpad=0.02,loc=1)
#plt.xticks(rotation=45,ha="right")
#ax4.set_xlabel('Counts')
#plt.legend(frameon=False,handletextpad=0.05)

#plt.tight_layout()
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_gap_distribution_3sub.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_gap_distribution_3sub.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_gap_distribution_3sub.eps')


# In[45]:


len(ldelta_alt)


# In[15]:


ceal_bin,short_bin = [],[]
for j in lbin_alt:
    try:
        ceal_bin.append(float(len(j[j<360.0]))*100.0/len(j))
        short_bin.append(float(len(j[j<60.0]))*100.0/len(j))
    except ZeroDivisionError:
        ceal_bin.append(np.nan)
        short_bin.append(np.nan)
ceal_bin = np.array(ceal_bin)
short_bin = np.array(short_bin)


# In[16]:


ceal_bin


# In[17]:


lbin_alt


# In[58]:


lbin_lon


# In[140]:


lbin_ndays


# In[57]:


plt.figure(figsize=(12,5))
axx = plt.subplot2grid((1,4),(0,0),colspan=3)
axx.plot(lbin_lon,ceal_bin,'x-',color='darkgoldenrod',label='extent <60 m')
axx.plot(lbin_lon,short_bin,'+--',color='sienna', label='extent <360 m')
axx.set_ylabel('Cases with small extent [$\%$]')
axx.set_ylim([0,100])
axx.spines["left"].set_position(("axes", -0.12))
axx.spines["left"].set_visible(True)
axx.yaxis.set_label_position('left')
axx.yaxis.set_ticks_position('left')
axx.set_xlim([0,16])
axx.yaxis.label.set_color('goldenrod')
axx.tick_params(axis='y', colors=('goldenrod'))
axx.spines["left"].set_edgecolor('goldenrod')
axx.set_xlabel(u'Longitude [°]')
lgx = axx.legend(frameon=False)
lgx.get_texts()[0].set_color('darkgoldenrod')
lgx.get_texts()[1].set_color('sienna')

ax = axx.twinx()
bow = ax.boxplot(lbin_alt,positions=lbin_lon,showfliers=True,sym='.')
ax.plot(ldelta_lon,ldelta_alt,'k.',alpha=0.01,zorder=-1)
for j,nn in enumerate(lbin_ndays):
    if nn>0:
        ax.text(lbin_lon[j],np.nanmax(bow['medians'][j].get_data()[1])+30,'{:2.0f}'.format(nn),color='blue',horizontalalignment='center')
    else:
        ax.text(lbin_lon[j],50,'{:2.0f}'.format(nn),color='blue',horizontalalignment='center')
ax.set_ylim(-50,3700)
ax.set_ylabel('Vertical extent [m]')
ax.set_xlabel(u'Longitude [°]')
pu.sub_note('a)')
ax.set_title('Distribution of gap vertical extent per longitude')
ax.spines["left"].set_visible(True)
ax.yaxis.set_label_position('left')
ax.yaxis.set_ticks_position('left')

ax2 = plt.subplot2grid((1,12),(0,9),sharey=ax)
ax2.hist(delta_alt1,bins=50,range=[0,4000],normed=False,edgecolor='None',alpha=0.3,orientation='horizontal')
ax2.axhline(np.nanmean(delta_alt1),label='Mean')
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.axhline(np.nanmedian(delta_alt1),label='Median',ls='--')
ax2.annotate('{:4.0f} m'.format(np.nanmean(delta_alt1)),(1,np.nanmean(delta_alt1)*1.05),color='b',alpha=0.3)
ax2.annotate('{:4.0f} m'.format(np.nanmedian(delta_alt1)),(1,np.nanmedian(delta_alt1)*1.05),color='b',alpha=0.3)
ax2.annotate(' {:2.0f}%'.format(float(len(delta_alt1[delta_alt1<60.0]))*100.0/len(delta_alt1)),(1.4,2450.0),color='darkgoldenrod')
ax2.annotate('({:2.0f}%)'.format(float(len(delta_alt1[delta_alt1<360.0]))*100.0/len(delta_alt1)),(1.4,2300.0),color='sienna')
ax2.set_title(u'0°-5°')
pu.sub_note('b)')
ax2.set_xticks([0,3,6])

ax3 = plt.subplot2grid((1,12),(0,10),sharey=ax)
ax3.hist(delta_alt2,bins=50,range=[0,4000],normed=False,edgecolor='None',alpha=0.3,orientation='horizontal')
ax3.axhline(np.nanmean(delta_alt2),label='Mean')
plt.setp(ax3.get_yticklabels(), visible=False)
ax3.axhline(np.nanmedian(delta_alt2),label='Median',ls='--')
ax3.annotate('{:4.0f} m'.format(np.nanmean(delta_alt2)),(1.5,np.nanmean(delta_alt2)*1.05),color='b',alpha=0.3)
ax3.annotate('{:4.0f} m'.format(np.nanmedian(delta_alt2)),(1.5,np.nanmedian(delta_alt2)*1.05),color='b',alpha=0.3)
ax3.annotate(' {:2.0f}%'.format(float(len(delta_alt2[delta_alt2<60.0]))*100.0/len(delta_alt2)),(3.1,2450.0),color='darkgoldenrod')
ax3.annotate('({:2.0f}%)'.format(float(len(delta_alt2[delta_alt2<360.0]))*100.0/len(delta_alt2)),(3.1,2300.0),color='sienna')
ax3.set_title('Gap extent distribution\n'+u'5°-10°')
#plt.xticks(rotation=45,ha="right")
ax3.set_xlabel('Counts')
pu.sub_note('c)')
ax3.set_xticks([0,6,12])
#plt.legend(frameon=True,handletextpad=0.02,bbox_to_anchor=(1.01, 0.93))#,ncol=2,columnspacing=0.2,loc=1)

ax4 = plt.subplot2grid((1,12),(0,11),sharey=ax)
ax4.hist(delta_alt3b,bins=50,range=[0,4000],normed=False,edgecolor='None',alpha=0.3,orientation='horizontal')
ax4.axhline(np.nanmean(delta_alt3b),label='Mean')
plt.setp(ax4.get_yticklabels(), visible=False)
ax4.axhline(np.nanmedian(delta_alt3b),label='Median',ls='--')
ax4.annotate('{:4.0f} m'.format(np.nanmean(delta_alt3b)),(0.5,np.nanmean(delta_alt3b)*1.05),color='b',alpha=0.3)
ax4.annotate('{:4.0f} m'.format(np.nanmedian(delta_alt3b)),(0.5,np.nanmedian(delta_alt3b)*1.05),color='b',alpha=0.3)
ax4.annotate(' {:2.0f}%'.format(float(len(delta_alt3b[delta_alt3b<60.0]))*100.0/len(delta_alt3b)),(1.2,2450.0),color='darkgoldenrod')
ax4.annotate('({:2.0f}%)'.format(float(len(delta_alt3b[delta_alt3b<360.0]))*100.0/len(delta_alt3b)),(1.2,2300.0),color='sienna')
ax4.set_title(u'10°-15°')
pu.sub_note('d)')
ax4.set_xticks([0,3,6])
plt.legend(frameon=True,handletextpad=0.02,bbox_to_anchor=(0.5, 0.93),loc=2)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_gap_distribution_3sub_smallextent.png',
            transparent=True,dpi=500)
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_gap_distribution_3sub_smallextent.pdf')
plt.savefig(fp+'plot_v3/ORACLES2016_4STAR_gap_distribution_3sub_smallextent.eps')


# In[292]:


plt.figure()
plt.hist(delta_alt3,bins=50)
plt.xlim(0,1000)
plt.figure()
plt.hist(delta_alt2,bins=50)
plt.xlim(0,1000)
plt.figure()
plt.hist(delta_alt1,bins=50)
plt.xlim(0,1000)


# ## Make a stats test to check if it can be at 60%

# In[60]:


import scipy.stats as stat


# In[ ]:


float(len(ldelta_alt[ldelta_alt<360.0]))*100.0/len(ldelta_alt)


# In[61]:


ldelta_alt[ldelta_alt<360.0]


# In[81]:


ceal_bin_ma = ceal_bin[np.isfinite(ceal_bin)]


# In[82]:


ceal_bin_ma


# In[83]:


t_val,p_val = stat.ttest_1samp(ceal_bin_ma,60.0)


# In[84]:


t_val,p_val


# Extracted data from Rajapshke et al. 2017 fig. 5c (CATS nightime 1064nm aerosol bottom and cloud top)

# In[64]:


aerbot = [2.0157932886957894,
2.04297524618854,
2.0651235078493,
2.083244812844466,
2.104386335338833,
2.1215009011676003,
2.145662641161156,
2.1668041636555166,
2.188952425316277,
2.2090872086442417,
2.2211680786410213,
2.217141121975427,
2.205060251978651,
2.195999599481066,
2.184925468650686,
2.176871555319501,
2.1949928603146667,
2.210093947810641,
2.228215252805807,
2.2433163403017815,
2.2513702536329667,
2.2443230794681774,
2.2302287311386024,
2.2181478611418264,
2.203046773645852,
2.1788850336522962,
2.162777206989926,
2.1376087278299707,
2.111433509503616,
2.083244812844466,
2.061096551183706,
2.072170682014086,
2.0872717695100604,
2.103379596172431,
2.117473944502006,
2.128548075332386,
2.11042677033722,
2.087271769510060,
2.0641167686829007,
2.0389482895229456,
2.0157932886957894,
2.022840462860575,
2.0308943761917604,
2.040961767855741,
2.0510291595197216,
2.0560628553517155,
2.058076333684511,
2.058076333684511,
2.0500224203533257,
2.047002202854131,
2.04297524618854,
2.006732636198201,
1.9604226345438853,
1.9151193720559654,
1.869816109568049,
1.8194791512481387,
1.8184724120817393,
1.8285398037457234,
1.8376004562433046,
1.8476678479072888,
1.8496813262400806,
1.8647824137360551,
1.8768632837328347,
1.8889441537296108,
1.8980048062271955,
1.9050519803919848,
1.9151193720559654,
1.9251867637199496,
1.9322339378847353,
1.9362608945503297,
1.9392811120495246,
1.9261935028863455,
1.8949845887280041,
1.872584642275644]


# In[65]:


cldtop = [1.8315600212449148,
1.81142523791695,
1.793303932921784,
1.7741758887602188,
1.75504784459865,
1.7359198004370882,
1.7480006704338642,
1.7480006704338642,
1.7480006704338642,
1.7510208879330627,
1.7469939312674683,
1.7480006704338642,
1.7480006704338642,
1.7480006704338642,
1.7439737137682698,
1.7318928437714938,
1.724845669606708,
1.7107513212771295,
1.7016906687795483,
1.685582842117178,
1.665448058789213,
1.6594076237908233,
1.6594076237908233,
1.6594076237908233,
1.6594076237908233,
1.654373927958833,
1.6594076237908233,
1.6483334929604432,
1.6483334929604432,
1.6483334929604432,
1.6392728404628585,
1.6362526229636671,
1.6161178396356988,
1.5990032738069324,
1.5879291429765523,
1.5738347946469773,
1.5547067504854084,
1.5345719671574436,
1.5144371838294823,
1.4953091396679135,
1.4711473996743578,
1.4590665296775818,
1.4429587030152113,
1.4308778330184317,
1.4137632671896654,
1.3976554405272914,
1.393628483861697,
1.3956419621944924,
1.4047026146920771,
1.4047026146920771,
1.3986621796936909,
1.3906082663625021,
1.3704734830345409,
1.352352178039375,
1.3332241338778061,
1.3130893505498378,
1.3010084805530617,
1.2828871755578959,
1.265772609729126,
1.2506715222331515,
1.235570434737177,
1.21946260807481,
1.20838847724443,
1.192280650582059,
1.1791930414188805,
1.1640919539229095,
1.1711391280876953,
1.1872469547500692,
1.2013413030796407,
1.2174491297420111,
1.2295299997387907,
1.234563695570781,
1.2275165214059953,
1.2275165214059953]


# In[70]:


raj_gap = -np.array(cldtop)+np.array(aerbot)


# In[73]:


raj_gap.mean()


# ## Save the delta gap to a file

# In[66]:


gap = {'delta_alt':delta_alt,'delta_lat':delta_lat,'delta_lon':delta_lon,
      'ldelta_alt':ldelta_alt,'ldelta_lat':ldelta_lat,'ldelta_lat_days':ldelta_lat_days,'ldelta_lon':ldelta_lon,
       'ldelta_lon_days':ldelta_lon_days,
       'delta_alt1':delta_alt1,'delta_lat1':delta_lat1,'delta_lon1':delta_lon1,
       'delta_alt2':delta_alt2,'delta_lat2':delta_lat2,'delta_lon2':delta_lon2,
       'delta_alt3':delta_alt3,'delta_lat3':delta_lat3,'delta_lon3':delta_lon3}


# In[39]:


hs.savemat(fp+'ORACLES_2016_gap_distance.mat',gap)


# ### Load the delta gap data

# In[9]:


gap = hs.loadmat(fp+'ORACLES_2016_gap_distance.mat')


# In[10]:


gap.keys()


# In[11]:


locals().update(gap)


# In[15]:


lbin_alt


# # Calculate the fine mode and coarse mode aod

# In[555]:


dd = su.sda(aods,wvl)


# In[556]:


dd.keys()


# In[563]:


plt.figure()
plt.hist(dd['tauc'][s['fl_acaod']],range=[0,1.5],bins=35)


# In[560]:


plt.figure()
plt.hist([dd['tauc'][s['fl_acaod']],dd['tauf'][s['fl_acaod']]],bins=25,range=(0,1.5),stacked=True,label=['Coarse','Fine'])
plt.legend(frameon=False)


# In[558]:


plt.figure()
#plt.hist(dd['tauc'][s['fl_acaod']],bins=25,range=(0,5.0),normed=True)
#plt.hist(dd['tauf'][s['fl_acaod']],bins=25,range=(0,5.0),normed=True)
plt.hist(dd['tau'][s['fl_acaod']],bins=25,range=(0,1.0),normed=True)
plt.hist(aods[s['fl_acaod'],4],bins=25,range=(0,1.0),normed=True)


# In[ ]:




