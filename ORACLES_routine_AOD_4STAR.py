
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

# In[11]:


import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[12]:


get_ipython().magic(u'matplotlib notebook')


# In[13]:


from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth

import hdf5storage as hs
from mpltools import color
from path_utils import getpath
from write_utils import nearest_neighbor
from tqdm import tqdm_notebook as tqdm


# In[93]:


import scipy.stats as st
import Sun_utils as su


# In[ ]:


from mpl_toolkits.basemap import Basemap


# In[15]:


fp = getpath('ORACLES')
fp


# In[16]:


vv = 'R2'


# # Load files

# ## 4STAR ict

# In[17]:


s = hs.loadmat(fp+'/aod_ict/{vv}/all_aod_ict_{vv}.mat'.format(vv=vv))


# In[18]:


s.keys()


# In[19]:


days = np.unique(s['days'])


# In[20]:


len(s['fl_QA'])


# ### Load the acaod flags

# In[21]:


fdat = getpath('4STAR_data',make_path=True,path='/mnt/c/Users/sleblanc/Research/4STAR_codes/data_folder')
fdat


# In[22]:


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


# In[300]:


s['fl_acaod'] = np.zeros_like(s['fl_QA'])
s['fl_acaod_noQA'] = np.zeros_like(s['fl_QA'])

for i,j in enumerate(days):
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


# In[24]:


len(s['fl_acaod'])


# In[25]:


len(s['fl_QA'])


# In[26]:


s['fl_below5'] = (s['fl_QA']) & (s['GPS_Alt']<5000.0)


# ## MODIS climatology

# In[27]:


fp


# From Rob Wood, email on 2016-11-08, 15:35, 
# These are just the all-years means (the black lines), for the month of September (DOY 244 to 273)

# In[28]:


m = sio.netcdf_file(fp+'data_other/climatology/mombl_oracles_routine_flight_NW-SE.nc')


# In[29]:


m.variables


# In[30]:


m2 = sio.netcdf_file(fp+'data_other/climatology/mombl_oracles_routine_flight_NW-SE_all.nc')


# In[31]:


m2.variables


# In[32]:


m2.variables['AOD_YRMEAN'].data.shape


# In[33]:


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

# In[34]:


import datetime
import load_utils as lu


# In[35]:


aac = sio.loadmat(fp+'data_other/MODIS_AAC/MODIS_AAC_per_bin_all_flights_20160801_20161031.mat')


# In[36]:


aac['data_aac'][0,0].dtype.names


# In[37]:


daac = []


# In[38]:


for i in xrange(len(aac['data_aac'][0,:])):
    print aac['data_aac'][0,i]['FlightDay'][0][0][0]
    daac.append(lu.recarray_to_dict(aac['data_aac'][0,i]))


# In[39]:


# simplify the dict
for i in xrange(len(daac)):
    for k in daac[i].keys():
        daac[i][k] = daac[i][k][0,0]


# In[40]:


for i in xrange(len(daac)):
    daac[i]['datetime'] = datetime.datetime.strptime(daac[i]['FlightDay'][0], 'A%Y%j')
    daac[i]['date'] = daac[i]['datetime'].strftime('%Y%m%d')
    print i,daac[i]['date'], aac['data_aac'][0,i]['FlightDay'][0][0][0]


# In[30]:


len(daac)


# In[41]:


i_aug = range(0,31)
i_sep = range(31,61)
i_oct = range(61,92)
i_flt = [30,34,38,40,42,55]


# days in the binned for flights only
# 
#  - 0 20160831 A2016244
#  - 1 20160904 A2016248
#  - 2 20160908 A2016252
#  - 3 20160910 A2016254
#  - 4 20160912 A2016256
#  - 5 20160925 A2016269

# # Subset for routine flight

# In[42]:


d_rtn = ['20160831','20160904','20160908','20160910','20160912','20160925']


# In[43]:


d_rtnf = [20160831.0,20160904.0,20160908.0,20160910.0,20160912.0,20160925.0]


# In[44]:


d_irtn = [2.0,4.0,6.0,7.0,8.0,13.0]


# In[45]:


s['days'][0]


# In[36]:


ff


# In[46]:


ff = []
for d in d_rtnf:
    ff.append(s['days']==d)


# In[47]:


fl = ((s['days']==d_rtnf[0]) | (s['days']==d_rtnf[1]) | (s['days']==d_rtnf[2]) | 
      (s['days']==d_rtnf[3]) | (s['days']==d_rtnf[4]) | (s['days']==d_rtnf[5]))


# In[48]:


s['fl_rtn'] = fl & s['fl']
s['fl_rtn']


# In[49]:


s['fl_rtna'] = fl & s['fl_acaod']
s['fl_rtna']


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

# In[50]:


bins = []
pos = []
for i,c in enumerate(cc[1][0:-2]):
    lon_fl = (s['Longitude'][s['fl_rtn']]>=c)&(s['Longitude'][s['fl_rtn']]<cc[1][i+1])
    bins.append(s['AOD0501'][s['fl_rtn']][lon_fl])
    pos.append((c+cc[1][i+1])/2.0)


# In[51]:


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

# In[52]:


pos3 = m.variables['LONGITUDE'].data[0,:]


# In[53]:


pos3


# In[54]:


lims3 = pos3-0.5
lims3= np.append(lims3,pos3[-1]+0.5)


# In[55]:


lims3


# In[56]:


bins3 = []
for i,c in enumerate(lims3[0:-1]):
    lon_fl = (s['Longitude'][s['fl_rtn']]>=c)&(s['Longitude'][s['fl_rtn']]<lims3[i+1])
    bins3.append(s['AOD0501'][s['fl_rtn']][lon_fl])


# In[57]:


bins3a = []
for i,c in enumerate(lims3[0:-1]):
    lon_fl = (s['Longitude'][s['fl_rtna']]>=c)&(s['Longitude'][s['fl_rtna']]<lims3[i+1])
    bins3a.append(s['AOD0501'][s['fl_rtna']][lon_fl])


# In[58]:


len(lims3)


# In[59]:


len(bins3)


# In[60]:


len(pos3)


# In[179]:


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


# In[61]:


s['fl_alt12'] = (s['GPS_Alt']>=600)&(s['GPS_Alt']<1200)
s['fl_alt16'] = (s['GPS_Alt']>=500)&(s['GPS_Alt']<1600)


# In[62]:


s['fl_alt12']


# In[63]:


s['fl_rtn_12'] = s['fl_alt12'] & s['fl_QA'] & s['fl_rtn']
s['fl_rtn_16'] = s['fl_alt16'] & s['fl_QA'] & s['fl_rtn']


# In[64]:


s['fl_rtn_12']


# In[65]:


bins3_12 = []
for i,c in enumerate(lims3[0:-1]):
    lon_fl = (s['Longitude'][s['fl_rtn_12']]>=c)&(s['Longitude'][s['fl_rtn_12']]<lims3[i+1])
    bins3_12.append(s['AOD0501'][s['fl_rtn_12']][lon_fl])


# In[209]:


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

# In[66]:


d_irtn = [2.0,4.0,6.0,7.0,8.0,13.0]


# In[67]:


s['days']==2.0


# In[68]:


flr2 = s['fl_rtn_16']&(s['days']==d_rtnf[0])
flr4 = s['fl_rtn_16']&(s['days']==d_rtnf[1])
flr6 = s['fl_rtn_16']&(s['days']==d_rtnf[2])
flr7 = s['fl_rtn_16']&(s['days']==d_rtnf[3])
flr8 = s['fl_rtn_16']&(s['days']==d_rtnf[4])
flr13 = s['fl_rtn_16']&(s['days']==d_rtnf[5])


# In[69]:


fr2 = s['fl_rtn']&(s['days']==d_rtnf[0])
fr4 = s['fl_rtn']&(s['days']==d_rtnf[1])
fr6 = s['fl_rtn']&(s['days']==d_rtnf[2])
fr7 = s['fl_rtn']&(s['days']==d_rtnf[3])
fr8 = s['fl_rtn']&(s['days']==d_rtnf[4])
fr13 = s['fl_rtn']&(s['days']==d_rtnf[5])


# In[70]:


fr = [fr2,fr4,fr6,fr7,fr8,fr13]
flr = [flr2,flr4,flr6,flr7,flr8,flr13]


# In[71]:


flr2


# In[72]:


flr2a = s['fl_rtna']&(s['days']==d_rtnf[0])
flr4a = s['fl_rtna']&(s['days']==d_rtnf[1])
flr6a = s['fl_rtna']&(s['days']==d_rtnf[2])
flr7a = s['fl_rtna']&(s['days']==d_rtnf[3])
flr8a = s['fl_rtna']&(s['days']==d_rtnf[4])
flr13a = s['fl_rtna']&(s['days']==d_rtnf[5])


# In[73]:


fr2a = s['fl_rtna']&(s['days']==d_rtnf[0])
fr4a = s['fl_rtna']&(s['days']==d_rtnf[1])
fr6a = s['fl_rtna']&(s['days']==d_rtnf[2])
fr7a = s['fl_rtna']&(s['days']==d_rtnf[3])
fr8a = s['fl_rtna']&(s['days']==d_rtnf[4])
fr13a = s['fl_rtna']&(s['days']==d_rtnf[5])


# In[74]:


fra = [fr2a,fr4a,fr6a,fr7a,fr8a,fr13a]
flra = [flr2a,flr4a,flr6a,flr7a,flr8a,flr13a]


# In[75]:


cls = ['green','blue','yellow','cyan','magenta','orange']


# In[76]:


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

# In[77]:


s['fl_alt_6'] = (s['GPS_Alt']<=600)
s['fl_rtn_6'] = s['fl_alt_6'] & s['fl_QA'] & fl
flrr2 = s['fl_rtn_6']&(s['days']==2.0)
flrr4 = s['fl_rtn_6']&(s['days']==4.0)
flrr6 = s['fl_rtn_6']&(s['days']==6.0)
flrr7 = s['fl_rtn_6']&(s['days']==7.0)
flrr8 = s['fl_rtn_6']&(s['days']==8.0)
flrr13 = s['fl_rtn_6']&(s['days']==13.0)
flrr = [flrr2,flrr4,flrr6,flrr7,flrr8,flrr13]


# In[78]:


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

# In[79]:


daac[0].keys()


# In[80]:


a['meanAODperbin']


# In[222]:


plt.figure()

for i,a in enumerate(daac):
   # print a['BinCenter'].shape,a['meanAODperbin'][0:-1].shape
   # print a['meanAODperbin']
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][0:-1],'x-',lw=1,color='orange',alpha=0.2)


# In[62]:


np.nanmean(ac,axis=0).shape


# In[227]:


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

plt.savefig(fp+'plot_v2/MODIS_Climatology_vs_AAC_4STAR_and_Meyer_AAC.png',transparent=True,dpi=600)


# In[70]:


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

plt.savefig(fp+'plot_v2/MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC.png',transparent=True,dpi=600)


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


# In[82]:


fp


# # Prepare the split plot for publication

# In[ ]:


plt.figure()
plt.plot()


# In[118]:





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


# In[75]:


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
    plt.plot(a['BinCenter'][0,:],a['meanAODperbin'][1:,0],'x-',lw=2,color=cls[i],alpha=0.4,
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
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MODIS AAC for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MODIS AAC for all Sept. 2016')
#ac_oct = []
#for i,a in enumerate([daac[j] for j in i_oct]):
#    ac_sep.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'o-',lw=0.5,color='lightblue',label='MODIS AAC for Oct. 2016')

plt.plot(a['BinCenter'][0,:],np.nanmean(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MODIS AAC for routine flights')

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
plt.plot(a['BinCenter'][0,:],np.nanmedian(ac_aug,axis=0),'<-',lw=0.5,color='purple',label='MODIS AAC for all Aug. 2016')
ac_sep = []
for i,a in enumerate([daac[j] for j in i_sep]):
    ac_sep.append(a['meanAODperbin'][1:,0])
plt.plot(a['BinCenter'][0,:],np.nanmedian(ac_sep,axis=0),'d-',lw=0.5,color='orange',label='MODIS AAC for all Sept. 2016')
#ac_oct = []
#for i,a in enumerate([daac[j] for j in i_oct]):
#    ac_sep.append(a['meanAODperbin'][1:,0])
#plt.plot(a['BinCenter'][0,:],np.nanmean(ac_sep,axis=0),'o-',lw=0.5,color='lightblue',label='MODIS AAC for Oct. 2016')

plt.plot(a['BinCenter'][0,:],np.nanmedian(ac,axis=0),'^-',lw=3,color='green',zorder=60,label='MODIS AAC for routine flights')

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


plt.savefig(fp+'plot_v2/ORACLES2016_MODIS_Climatology_vs_AAC_flagged_4STAR_and_Meyer_AAC_3_split.png',
            transparent=True,dpi=500)


# # Plot the histogram of acaod

# In[299]:


plt.figure()
plt.hist(s['AOD0501'][s['fl_QA']],bins=30,range=(0,1.2),edgecolor='None',label='All AOD',alpha=0.3,normed=True)
plt.hist(s['AOD0501'][s['fl_acaod']],bins=30,range=(0,1.2),edgecolor='None',label='AOD above clouds',alpha=0.3,normed=True)
plt.hist(s['AOD0501'][s['fl_alt_6']],bins=30,range=(0,1.2),edgecolor='None',label='AOD below 0.6 km',alpha=0.3,normed=True)

plt.xlabel('AOD 501 nm')
plt.ylabel('Normalized frequency counts')
plt.title('ORACLES 2016 4STAR AOD')
plt.legend(frameon=False)


# In[83]:


s['fl6'] = s['fl_alt_6'] & s['fl_QA']


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


# # Now get the angstrom

# In[85]:


nn = s.keys()
nn.sort()


# In[86]:


wvl = np.array([float(i[-4:]) for i in nn[0:24]])


# In[87]:


aods = np.array([s[i] for i in nn[0:24]]).T


# In[88]:


aods.shape


# In[165]:


uncaods = np.array([s[i] for i in nn[28:28+24]]).T


# In[89]:


s['polyaod'] = []
s['polylogaod'] = []


# In[94]:


for i,u in enumerate(s['Start_UTC']):
    pp = su.aod_polyfit(wvl,aods[i,:],polynum=3)
    s['polyaod'].append(pp)
    pl = su.logaod_polyfit(wvl,aods[i,:],polynum=3)
    s['polylogaod'].append(pl)
s['polyaod'] = np.array(s['polyaod'])
s['polylogaod'] = np.array(s['polylogaod'])


# In[95]:


s['angs'] = su.angstrom_from_logpoly(s['polylogaod'],[380.0,470.0,500.0,530.0,660.0,865.0,1250.0],polynum=3)


# In[96]:


s['angs'].shape


# In[97]:


awvl = [380.0,470.0,500.0,530.0,660.0,865.0,1250.0]


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


# ### Combine both histogram

# In[99]:


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


# In[132]:


s['fl_both'] = s['fl_acaod'] | s['fl6']


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


# # Make some figures of the map statistic distribution

# ## plot a map of the mean, median, and std of the AOD

# In[429]:


a,xe,ye,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['AOD0501'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]])
a = np.ma.masked_array(a,np.isnan(a))


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


# In[430]:


astd,xe,ye,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['AOD0501'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]],statistic=np.std)
astd = np.ma.masked_array(astd,np.isnan(astd))


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


# In[431]:


med,xe,ye,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['AOD0501'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]],statistic='median')
med = np.ma.masked_array(med,np.isnan(med))


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


# In[432]:


cnt,xe,ye,bn = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],s['Longitude'][s['fl_acaod']],s['AOD0501'][s['fl_acaod']],
                           bins=26,range=[[-25,-8],[0,16]],statistic='count')
cnt = np.ma.masked_array(cnt,np.isnan(cnt))


# In[433]:


io = np.where((cnt.data>0.0) & (astd.data<1.0))


# In[434]:


cnt.shape


# In[435]:


yy = [(y+ye[i+1])/2.0 for i,y in enumerate(ye[0:-1])]
yys = np.array(yy*26).reshape(26,-1)


# In[436]:


xx = [(x+xe[i+1])/2.0 for i,x in enumerate(xe[0:-1])]
xxs = np.array(xx*26).reshape(26,-1)


# In[437]:


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


# ## Plot a map of the distribution of the Angstrom exponent

# In[107]:


astat = {}


# In[108]:


astat['mean'],astat['x'],astat['y'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs'][s['fl_acaod'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]])
astat['mean'] = np.ma.masked_array(astat['mean'],np.isnan(astat['mean']))


# In[182]:


plt.figure()
p = plt.pcolor(astat['y'],astat['x'],astat['mean'],vmin=0.0,vmax=2.2,cmap='plasma')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')

plt.title('Mean Angstrom')

cb = plt.colorbar(p,extend='both')
cb.set_label('Mean Angstrom exponent \nof Above cloud AOD at 660 nm')
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_meanAngstrom_map.png',
            transparent=True,dpi=500)


# In[109]:


astat['median'],astat['xe'],astat['ye'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs'][s['fl_acaod'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic='median')
astat['median'] = np.ma.masked_array(astat['median'],np.isnan(astat['median']))


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


# In[110]:


astat['std'],astat['xs'],astat['ys'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs'][s['fl_acaod'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic=np.std)
astat['std'] = np.ma.masked_array(astat['std'],np.isnan(astat['std']))


# In[213]:


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


# In[111]:


astat['cnt'],astat['xs'],astat['ys'],astat['bin'] = st.binned_statistic_2d(s['Latitude'][s['fl_acaod']],
                                                                          s['Longitude'][s['fl_acaod']],
                                                                          s['angs'][s['fl_acaod'],ia],
                                                                          bins=26,range=[[-25,-8],[0,16]],
                                                                           statistic='count')
astat['cnt'] = np.ma.masked_array(astat['cnt'],np.isnan(astat['cnt']))


# In[112]:


iao = np.where((astat['cnt'].data>0.0) & (astat['std'].data<1.0))


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

# In[442]:


def mapfig(ax=plt.gca()):
    m = Basemap(projection='merc',llcrnrlat=-26,urcrnrlat=-8,llcrnrlon=0,urcrnrlon=16,resolution='l')
    m.drawcoastlines()
    m.drawmeridians(np.linspace(0,16,9),labels=[0,0,0,1])
    m.drawparallels(np.linspace(-26,-8,10),labels=[1,0,0,0])
    return m


# In[443]:


fig = plt.figure(figsize=(11,4.5))
ax1 = plt.subplot(1,2,1)
m = mapfig()

mx,my = m(astat['y'],astat['x'])
p = plt.pcolor(mx,my,astat['mean'],vmin=0.0,vmax=2.2,cmap='gist_rainbow')
plt.title('Mean')
cb = plt.colorbar(p,extend='both')

ax2 = plt.subplot(1,2,2)
m = mapfig(ax=ax2)
mxx,myy = m(np.array(yy),np.array(xx))

p = plt.scatter(mxx[iao[1]],myy[iao[0]],astat['cnt'].data[iao[0],iao[1]].flatten()/10,
                c=astat['std'].data[iao[0],iao[1]].flatten(),
               marker='s',edgecolor='None',cmap='viridis')
plt.xlim(0,16)
plt.ylim(-26,-8)
#plt.xlabel(u'Longitude [°]')
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


# In[ ]:


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


# # Prepare the spectral AOD mean, median, pca, std...

# In[144]:


aods.shape


# In[145]:


s['fl_acaod']


# In[150]:


ns = len(s['Start_UTC'][s['fl_acaod']])


# In[151]:


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


# In[159]:


meanaod = np.nanmean(aods[s['fl_acaod'],:],axis=0)
medianaod = np.nanmedian(aods[s['fl_acaod'],:],axis=0)
stdaod = np.nanstd(aods[s['fl_acaod'],:],axis=0)


# In[166]:


meanuncaod = np.nanmean(uncaods[s['fl_acaod'],:],axis=0)


# In[174]:


plt.figure()
plt.plot(wvl,meanaod,'k-+',label='Mean')
plt.errorbar(wvl,meanaod,yerr=meanuncaod,color='k')
plt.plot(wvl,meanaod+stdaod,'--x',color='grey',lw=0.4,label='Standard Deviation')
plt.plot(wvl,meanaod-stdaod,'--x',color='grey',lw=0.4)
plt.plot(wvl,medianaod,'-o',lw=2.0,color='lightblue',label='Median',zorder=-1)
plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.title('ACAOD average spectra from 4STAR at selected wavelengths')

plt.legend(frameon=False)
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_mean_ACAOD_spectra.png',
            transparent=True,dpi=500)


# ## Get the Principal components from PCA

# In[183]:


ivalid = np.all(np.isfinite(aods[s['fl_acaod'],:]),axis=1)
aods_valid = aods[s['fl_acaod'],:][ivalid,:]


# In[201]:


from sklearn.decomposition import PCA
pca = PCA(n_components=20)
pca.fit(aods_valid)


# In[217]:


pca.explained_variance_ratio_[0]


# In[203]:


pca.components_.shape


# In[221]:


plt.figure()
plt.plot(wvl,meanaod,'k-+',label='Mean')
plt.errorbar(wvl,meanaod,yerr=meanuncaod,color='k')
plt.plot(wvl,meanaod+stdaod,'--x',color='grey',lw=0.4,label='Standard Deviation')
plt.plot(wvl,meanaod-stdaod,'--x',color='grey',lw=0.4)
plt.plot(wvl,medianaod,'-o',lw=2.0,color='lightblue',label='Median',zorder=-1)
plt.plot(wvl,pca.components_[0,:],'-*',lw=0.8,color='coral',zorder=-2,
         label='Principal Component\nrepresenting {:2.1f}% of variance'.format(pca.explained_variance_ratio_[0]*100.0))
plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.title('ACAOD representative spectra')

plt.legend(frameon=False)
plt.savefig(fp+'plot_v2/ORACLES2016_4STAR_mean_pca_ACAOD_spectra.png',
            transparent=True,dpi=500)


# # Get the spacing and altitudes of gap from ACAOD altitudes and distances

# In[301]:


ii_flacaod = np.where(s['fl_acaod_noQA'])[0] 


# In[302]:


ii_flacaod[0]


# In[414]:


disc_flacaod = np.where(np.diff(ii_flacaod,1)>1)[0]
disc_flacaod_long = np.where(np.diff(ii_flacaod,1)>150)[0]


# In[304]:


plt.figure()
plt.plot(ii_flacaod[0:100],'.')
plt.plot(np.arange(100)[disc_flacaod[0:6]],ii_flacaod[0:100][disc_flacaod[0:6]],'rx')
plt.plot(np.arange(100)[disc_flacaod[0:6]+1],ii_flacaod[0:100][disc_flacaod[0:6]+1],'go')


# In[281]:


discontinuity_istart =  ii_flacaod[np.append(0,disc_flacaod[:-1]+1)]
discontinuity_iend =  ii_flacaod[disc_flacaod]


# In[415]:


discontinuity_istart_long =  ii_flacaod[np.append(0,disc_flacaod_long[:-1]+1)]
discontinuity_iend_long =  ii_flacaod[disc_flacaod_long]


# In[317]:


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


# In[416]:


ldelta_alt,ldelta_lon,ldelta_lat = [],[],[]
for i,start in enumerate(discontinuity_istart_long):
    try:
        ma = np.nanmax(s['GPS_Alt'][start:discontinuity_iend_long[i]])
        mi = np.nanmin(s['GPS_Alt'][start:discontinuity_iend_long[i]])
        ldelta_alt.append(ma-mi)
        ldelta_lon.append(np.nanmean(s['Longitude'][start:discontinuity_iend_long[i]]))
        ldelta_lat.append(np.nanmean(s['Latitude'][start:discontinuity_iend_long[i]]))
    except:
        pass
ldelta_alt = np.array(ldelta_alt)
ldelta_lon = np.array(ldelta_lon)
ldelta_lat = np.array(ldelta_lat)


# In[306]:


np.nanmax(delta_alt)


# In[307]:


delta_alt


# In[318]:


plt.figure()
plt.hist(delta_alt,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
#plt.hist(delta_alt_all,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
plt.yscale('log')
plt.xscale('log')
plt.axvline(np.nanmean(delta_alt_all),label='mean')
plt.axvline(np.nanmedian(delta_alt_all),label='median',ls='--')
plt.legend(frameon=False)


# In[417]:


plt.figure()
plt.hist(ldelta_alt,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
#plt.hist(delta_alt_all,bins=50,range=[0,2500],normed=True,edgecolor='None',alpha=0.3)
#plt.yscale('log')
#plt.xscale('log')
plt.axvline(np.nanmean(ldelta_alt),label='mean')
plt.axvline(np.nanmedian(ldelta_alt),label='median',ls='--')
plt.legend(frameon=False)


# ## make a box plot of the altitude differences in gap vs. longitude

# In[322]:


delta_lon


# In[324]:


numlon,binlon = np.histogram(delta_lon,range=(0,16))


# In[325]:


binlon


# In[333]:


bin_alt,bin_lon = [],[]
for i in xrange(16):
    fllon = (delta_lon>=i) & (delta_lon<(i+1.0))
    bin_alt.append(delta_alt[fllon])
    bin_lon.append(np.mean([i,i+1.0]))


# In[418]:


lbin_alt,lbin_lon = [],[]
for i in xrange(16):
    fllon = (ldelta_lon>=i) & (ldelta_lon<(i+1.0))
    lbin_alt.append(ldelta_alt[fllon])
    lbin_lon.append(np.mean([i,(i+1.0)]))


# In[330]:


bin_alt


# In[340]:


plt.figure()
plt.boxplot(bin_alt,positions=bin_lon)
plt.yscale('log')
plt.ylim(0,3700)
plt.ylabel('Gap distance [m]')
plt.xlabel(u'Longitude [°]')


# In[361]:


plt.figure()
plt.boxplot(lbin_alt,positions=lbin_lon)
plt.plot(ldelta_lon,ldelta_alt,'.',alpha=0.5)
#plt.yscale('log')
plt.ylim(-50,3700)
plt.ylabel('Gap distance [m]')
plt.xlabel(u'Longitude [°]')
plt.title('Distribution of gap altitude per longitude')
plt.savefig()


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

