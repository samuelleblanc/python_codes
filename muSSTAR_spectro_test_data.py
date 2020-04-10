#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     muSSTAR_spectro_test_data
# 
# Purpose:  
# 
#     Test out the data from different spectrometers to be used with musstar
#   
# Input:
# 
#     none
#   
# Output:
# 
#     figures
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - load_utils.py : for loading OMI HDF5 files
#     - matplotlib
#     - numpy
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-05-18
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2020-04-01

# # Prepare python environment

# In[1]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
import os
matplotlib.rc_file(os.path.join(os.getcwd(),'file.rc'))
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


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[5]:


fp = getpath('musstar')


# In[6]:


import pandas as pd


# # Load files for the Darks analysis of C11482

# In[7]:


f = fp+'spectro_testdata/C11482GA_Data/'


# In[33]:


fs = os.listdir(f)


# In[49]:


d = []
s = []
for ff in fs:
    print ff
    dat = pd.read_csv(f+ff,header=21)
    d.append(dat)
    try:
        sp = np.array([dat['{}'.format(i+1)][1:] for i in xrange(len(dat.keys())-4)])
    except:
        sp = np.array([np.nan, np.nan])
    
    s.append(sp)
    


# In[51]:


s


# In[38]:


dat.keys()


# In[18]:


dat = pd.read_csv(f+'File_HiGain_50msec.csv',header=21)


# In[19]:


dat


# ## Plot out the darks

# In[24]:


plt.figure()
plt.plot(dat['WaveLength'][1:],dat['1'][1:])


# In[26]:


sp = np.array([dat['{}'.format(i+1)][1:] for i in xrange(5000)])


# In[29]:


dat['WaveLength'][200]


# In[58]:


ff.replace('_','\_')


# In[73]:


w = 200
plt.figure(figsize=(12,5))
for i,ss in enumerate(s):
    try:
        plt.plot(ss[:,w],'.',label=fs[i].replace('_','\_'))
    except:
        pass
plt.xlabel('time [s]')
plt.ylabel('Dark at {} nm'.format(dat['WaveLength'][w]))
plt.title('Spectrometer C11482GA dark test')
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
plt.savefig(f+'Darks_1240nm.png',dpi=600,transparent=True)


# In[77]:


w = 75
plt.figure(figsize=(12,5))
for i,ss in enumerate(s):
    try:
        plt.plot(ss[:,w],'.',label=fs[i].replace('_','\_'))
    except:
        pass
plt.xlabel('time [s]')
plt.ylabel('Dark at {} nm'.format(dat['WaveLength'][w]))
plt.title('Spectrometer C11482GA dark test')
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
plt.savefig(f+'Darks_1020nm.png',dpi=600,transparent=True)


# In[80]:


w = 470
plt.figure(figsize=(12,5))
for i,ss in enumerate(s):
    try:
        plt.plot(ss[:,w],'.',label=fs[i].replace('_','\_'))
    except:
        pass
plt.xlabel('time [s]')
plt.ylabel('Dark at {} nm'.format(dat['WaveLength'][w]))
plt.title('Spectrometer C11482GA dark test')
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
plt.savefig(f+'Darks_1650nm.png',dpi=600,transparent=True)


# In[107]:


fs


# # Load files for the direct beam Hammamatsu data

# In[81]:


fd = fp+'spectro_testdata/Hamamatsu_Data/'


# In[82]:


ffd = os.listdir(fd)


# In[83]:


ffd


# In[84]:


dd = []
sd = []
for ff in ffd:
    print ff
    dat = pd.read_csv(fd+ff,header=21)
    dd.append(dat)
    try:
        sp = np.array([dat['{}'.format(i+1)][1:] for i in xrange(len(dat.keys())-4)])
    except:
        sp = np.array([np.nan, np.nan])
    
    sd.append(sp)
    


# ## Plot the spectra

# In[166]:


w = 75
plt.figure(figsize=(12,5))
mk = '.xsdo^+'
for i,ss in enumerate(sd):
    try:
        plt.plot(ss[:,w],mk[i],label=ffd[i].replace('_','\_'),markeredgecolor='None')
    except:
        pass
plt.xlabel('time [s]')
plt.ylabel('Raw direct beam at {} nm'.format(dat['WaveLength'][w]))
plt.title('Spectrometer C11482GA direct beam')
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.legend(loc='center left',bbox_to_anchor=(1, 0.5))
plt.savefig(fd+'Direct_beam_1020nm.png',dpi=600,transparent=True)


# In[108]:


s[7]


# In[109]:


np.nanstd(s[7],axis=0).shape


# In[115]:


plt.figure()
plt.plot(dat['WaveLength'][1:],np.nanmean(sd[4][6:,:],axis=0),'-',label='direct beam - mean')
plt.plot(dat['WaveLength'][1:],np.nanmean(sd[4][6:,:],axis=0)+np.nanstd(sd[4][6:,:],axis=0),'.k',label='direct beam - std')
plt.plot(dat['WaveLength'][1:],np.nanmean(sd[4][6:,:],axis=0)-np.nanstd(sd[4][6:,:],axis=0),'.k')
plt.plot(dat['WaveLength'][1:],np.nanmean(s[7],axis=0),'-',label='dark - mean')
plt.plot(dat['WaveLength'][1:],np.nanmean(s[7],axis=0)+np.nanstd(s[7],axis=0),'.',color='lightgrey',label='dark - std')
plt.plot(dat['WaveLength'][1:],np.nanmean(s[7],axis=0)-np.nanstd(s[7],axis=0),'.',color='lightgrey')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Raw counts')
plt.title('Hammamatsu InGaAs C11482GA')
plt.savefig(fd+'raw_spectra_dir_darks.png',dpi=600,transparent=True)


# In[158]:


plt.figure()
plt.plot(sd[2][7:,75])


# In[164]:


plt.figure()
nm = np.nanmean(sd[2][7:,:],axis=0)
plt.plot(dat['WaveLength'][1:],np.nanmean(sd[2][7:,:],axis=0)/nm*100.0,'-',label='direct beam (3 ms) - mean')
plt.plot(dat['WaveLength'][1:],(np.nanmean(sd[2][7:,:],axis=0)+np.nanstd(sd[2][7:,:],axis=0))/nm*100.0,'.k',label='direct beam - std')
plt.plot(dat['WaveLength'][1:],(np.nanmean(sd[2][7:,:],axis=0)-np.nanstd(sd[2][7:,:],axis=0))/nm*100.0,'.k')
#plt.plot(dat['WaveLength'][1:],np.nanmean(s[7],axis=0),'-',label='dark - mean')
plt.plot(dat['WaveLength'][1:],(nm+np.nanstd(s[7],axis=0))/nm*100.0,'.',color='lightgrey',label='dark (50 ms) - std')
plt.plot(dat['WaveLength'][1:],(nm-np.nanstd(s[7],axis=0))/nm*100.0,'.',color='lightgrey')
plt.plot(dat['WaveLength'][1:],nm*101.0/nm,'-',color='red',label='1$\%$ variation (0.01 AOD)')
plt.plot(dat['WaveLength'][1:],nm*99.0/nm,'-',color='red')
plt.ylim(99.99,101.1)
plt.grid()
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalzed Raw counts [$\\%$]')
plt.title('Hammamatsu InGaAs C11482GA HiGain')
plt.savefig(fd+'raw_spectra_dir_darks_nomralized_hi50.png',dpi=600,transparent=True)


# In[165]:


plt.figure()
nm = np.nanmean(sd[4][6:,:],axis=0)
plt.plot(dat['WaveLength'][1:],np.nanmean(sd[4][6:,:],axis=0)/nm*100.0,'-',label='direct beam - mean')
plt.plot(dat['WaveLength'][1:],(np.nanmean(sd[4][6:,:],axis=0)+np.nanstd(sd[4][6:,:],axis=0))/nm*100.0,'.k',label='direct beam - std')
plt.plot(dat['WaveLength'][1:],(np.nanmean(sd[4][6:,:],axis=0)-np.nanstd(sd[4][6:,:],axis=0))/nm*100.0,'.k')
#plt.plot(dat['WaveLength'][1:],np.nanmean(s[7],axis=0),'-',label='dark - mean')
plt.plot(dat['WaveLength'][1:],(nm++np.nanstd(s[8],axis=0))/nm*100.0,'.',color='lightgrey',label='dark - std')
plt.plot(dat['WaveLength'][1:],(nm-np.nanstd(s[8],axis=0))/nm*100.0,'.',color='lightgrey')
plt.plot(dat['WaveLength'][1:],nm*101.0/nm,'-',color='red',label='1$\%$ variation (0.01 AOD)')
plt.plot(dat['WaveLength'][1:],nm*99.0/nm,'-',color='red')
plt.ylim(99.99,101.1)
plt.grid()
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalzed Raw counts [$\\%$]')
plt.title('Hammamatsu InGaAs C11482GA LowGain 50msec')
plt.savefig(fd+'raw_spectra_dir_darks_nomralized_lo50.png',dpi=600,transparent=True)


# In[153]:


fs[7]


# In[126]:


ffd[4]


# In[127]:


ffd[2]


# In[147]:


fs[8]


# # Load the files for analysis of C11486

# In[10]:


fu = fp+'spectro_testdata/C11486GA/'


# In[11]:


fsu = os.listdir(fu)


# In[12]:


fsu


# In[ ]:


pd.read_excel()


# In[18]:


du = []
for u in fsu:
    ddu = pd.read_excel(fu+u,header=22)
    du.append(ddu)


# In[23]:


du[2].drop(0)


# In[24]:


du[1]['WaveLength']


# In[35]:


len(du[1].keys())


# In[40]:


sp = np.stack([du[1][i].drop(0) for i in range(1,len(du[1].keys())-4)])


# In[41]:


sp.shape


# In[43]:


wvl = du[1]['WaveLength'].drop(0).to_numpy()


# In[55]:


sp1 = np.stack([du[2][i].drop(0) for i in range(1,len(du[2].keys())-4)])


# In[56]:


wvl1 = du[2]['WaveLength'].drop(0).to_numpy()


# ## Plot out some spectra

# In[48]:


from mpltools import color


# In[75]:


plt.figure()
cmap = 'plasma'
color.cycle_cmap(length=len(sp[35:,0])+1,cmap=cmap,ax=plt.gca())
plt.plot(wvl,sp.T[:,35:])
plt.ylabel('Digital Counts')
plt.xlabel('Wavelength [nm]')
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(cmap))
scalarmap.set_array(range(len(sp[35:,0])))
plt.colorbar(scalarmap)
plt.plot(wvl1,np.nanmean(sp.T[:,35:],axis=1),'-k',label='mean')
plt.plot(wvl1,np.nanmean(sp.T[:,35:],axis=1)-np.nanstd(sp.T[:,35:],axis=1),'--k',
         label='std {:2.2f}\%'.format(np.nanstd(sp.T[:,35:],axis=1)[i1020]/np.nanmean(sp.T[:,35:],axis=1)[i1020]*100.0))
plt.plot(wvl1,np.nanmean(sp.T[:,35:],axis=1)+np.nanstd(sp.T[:,35:],axis=1),'--k')
plt.legend(frameon=False)
plt.title('C11486GA spectrometer Low time integration')


# In[52]:


i1020 = np.argmin(abs(wvl-1020.0))
i1240 = np.argmin(abs(wvl-1240.0))
i1630 = np.argmin(abs(wvl-1630.0))


# In[60]:


plt.figure()
plt.plot(sp[:,i1020],label='1020 nm')
plt.plot(sp[:,i1240],label='1240 nm')
plt.plot(sp[:,i1630],label='1630 nm')
plt.legend()
plt.ylabel('Digital Counts')
plt.xlabel('Measurement number')


# In[74]:


plt.figure()
cmap = 'plasma'
color.cycle_cmap(length=len(sp1[25:,0])+1,cmap=cmap,ax=plt.gca())
plt.plot(wvl1,sp1.T[:,25:])
plt.ylabel('Digital Counts')
plt.xlabel('Wavelength [nm]')
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap(cmap))
scalarmap.set_array(range(len(sp1[25:,0])))
plt.colorbar(scalarmap)
plt.plot(wvl1,np.nanmean(sp1.T[:,25:],axis=1),'-k',label='mean')
plt.plot(wvl1,np.nanmean(sp1.T[:,25:],axis=1)-np.nanstd(sp1.T[:,25:],axis=1),'--k',
         label='std {:2.2f}\%'.format(np.nanstd(sp1.T[:,25:],axis=1)[i1020]/np.nanmean(sp1.T[:,25:],axis=1)[i1020]*100.0))
plt.plot(wvl1,np.nanmean(sp1.T[:,25:],axis=1)+np.nanstd(sp1.T[:,25:],axis=1),'--k')
plt.legend(frameon=False)
plt.title('C11486GA spectrometer High time integration')


# In[58]:


plt.figure()
plt.plot(sp1[:,i1020],label='1020 nm')
plt.plot(sp1[:,i1240],label='1240 nm')
plt.plot(sp1[:,i1630],label='1630 nm')
plt.legend()
plt.ylabel('Digital Counts')
plt.xlabel('Measurement number')


# In[ ]:




