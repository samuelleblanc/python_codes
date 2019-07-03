
# coding: utf-8

# # Info
# Name:  
# 
#     ## insert name here
# 
# Purpose:  
# 
#     ## add description
#   
# Input:
# 
#     ## inputs
#   
# Output:
# 
#     ##variables, figures and save files...
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
#     Modified: 

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


# In[6]:


fp = getpath('musstar')


# In[4]:


import pandas as pd


# # Load files for the Darks

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


# # Load files for the direct beam

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


# # Plotting
# Present some fo the early plots here
