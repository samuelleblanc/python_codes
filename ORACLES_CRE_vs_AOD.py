#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_CRE_vs_AOD
# 
# Purpose:  
# 
#     Plot out results of the CRE vs. AOD for ORACLES
#     Quantify the impact of overlying aerosol on cloud radiative effect
#   
# Input:
# 
#     None
#   
# Output:
# 
#     Figures
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-06-24
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


# In[3]:


fp =getpath('ORACLES')


# # Load files

# ## Load the theory results

# In[4]:


mat = hs.loadmat(fp+'rtm/ORACLES_theory_CRE_v1.mat')


# In[5]:


mat.keys()


# In[16]:


mat['star_aero_CRE']['dn'].shape


# In[17]:


CRE_aero = mat['star_aero_CRE']['up'][:,:,:,:,:,2] -mat['star_aero_CRE_clear']['up'][:,:,:,:,:,2] 
CRE_noaero = mat['star_noaero_CRE']['up'][:,:,:,:,:,2] -mat['star_noaero_CRE_clear']['up'][:,:,:,:,:,2] 


# In[18]:


rCRE_sur_aero = mat['star_aero_C'][:,:,:,:,:,0]/mat['star_aero_CRE_clear']['dn'][:,:,:,:,:,2]*100.0 
rCRE_sur_noaero = mat['star_noaero_C'][:,:,:,:,:,0]/mat['star_aero_CRE_clear']['dn'][:,:,:,:,:,2]*100.0 
rCRE_toa_aero = CRE_aero/mat['star_aero_CRE_clear']['dn'][:,:,:,:,:,2]*100.0 
rCRE_toa_noaero = CRE_noaero/mat['star_aero_CRE_clear']['dn'][:,:,:,:,:,2]*100.0


# In[11]:


cod_arr = [1.0,2.5,5.0,7.5,10.0,12.5,15.0,20.0]
ref_arr = [2.0,5.0,7.5,10.0,12.5,15.0]


# In[12]:


# set the range of ext, ssa, and asy to model at 500 nm
ext_arr = [0.05,0.1,0.15,0.2,0.3]
ssa_arr = [0.75,0.8,0.85,0.875,0.9]
asy_arr = [0.6,0.65,0.7,0.75]


# In[14]:


len(cod_arr), len(ref_arr),len(ext_arr),len(ssa_arr),len(asy_arr)


# In[15]:


mat['star_noaero_CRE_clear']['dn'].shape


# ### Plot out the CRE dependences
# Do some of the calculations to the data here

# In[33]:


plt.figure()
plt.plot(cod_arr,rCRE_toa_aero[:,0,0,0,0],':',label='ext:{}'.format(ext_arr[0]))
plt.plot(cod_arr,rCRE_toa_aero[:,0,0,0,1],'-',label='ext:{}'.format(ext_arr[1]))
plt.plot(cod_arr,rCRE_toa_aero[:,0,0,0,2],'+-',label='ext:{}'.format(ext_arr[2]))
plt.plot(cod_arr,rCRE_toa_aero[:,0,0,0,3],'x-',label='ext:{}'.format(ext_arr[3]))
plt.plot(cod_arr,rCRE_toa_aero[:,0,0,0,0],'d-',label='ext:{}'.format(ext_arr[4]))
plt.xlabel('COD')
plt.ylabel('rCRE_aero')
plt.legend(loc=4)


# ## Load the lut with and without overlying aerosol

# In[35]:


from load_utils import load_from_json


# In[36]:


oa = hs.loadmat(fp+'model/v1_ORACLES_lut.mat')
wa = hs.loadmat(fp+'model/v2_ORACLES_lut.mat')


# In[37]:


aero_oa = load_from_json(fp+'model/aero_save.txt')
aero_wa = load_from_json(fp+'model/aero_save_v2.txt')


# In[38]:


oa.keys(), wa.keys()


# In[39]:


aero_oa


# In[40]:


aero_wa


# In[41]:


oa['sza']


# In[42]:


oa['rad'].shape


# In[48]:


oa['sza'].shape, oa['tau'].shape, oa['zout'], oa['phase'], oa['ref'].shape, oa['wvl'].shape


# In[49]:


oa['ref']


# In[50]:


oa['tau']


# In[69]:


def plot_greys(fig=None,ax=None):
    " Plotting of grey regions that indicates the different wavelenght regions where the parameters are defined. "
    cl = '#DDDDDD'
    plt.axvspan(1000,1077,color=cl) #eta1
    plt.axvspan(1192,1194,color=cl) #eta2
    plt.axvspan(1492,1494,color=cl) #eta3
    plt.axvspan(1197,1199,color=cl); plt.axvspan(1235,1237,color=cl);  #eta4
    plt.axvspan(1248,1270,color=cl) #eta5
    plt.axvspan(1565,1644,color=cl) #eta6
    plt.axvspan(1000,1050,color=cl) #eta7
    plt.axvspan(1493,1600,color=cl) #eta8
    plt.axvspan(1000,1077,color=cl) #eta9
    plt.axvspan(1200,1300,color=cl) #eta10
    plt.axvspan(530 ,610 ,color=cl) #eta11
    plt.axvspan(1039,1041,color=cl) #eta12
    plt.axvspan(999 ,1001,color=cl); plt.axvspan(1064,1066,color=cl);  #eta13
    plt.axvspan(599 ,601 ,color=cl); plt.axvspan(869 ,871 ,color=cl);  #eta14
    plt.axvspan(1565,1634,color=cl); #eta15
    


# In[74]:


plt.figure()
plt.plot(oa['wvl'],oa['rad'][0,:,0,6,17,1],label='$\\tau$={t}, r$_{{eff}}$={r} $\\mu$m, aod$_{{500nm}}$={a:1.1f}'.format(t=oa['tau'][17],r=oa['ref'][6],a=0.6))
plt.plot(wa['wvl'],wa['rad'][0,:,0,6,17,1],label='$\\tau$={t}, r$_{{eff}}$={r} $\\mu$m, aod$_{{500nm}}$={a:1.1f}'.format(t=wa['tau'][17],r=wa['ref'][6],a=1.2))
plt.plot(oa['wvl'],oa['rad'][0,:,0,6,23,1],label='$\\tau$={t}, r$_{{eff}}$={r} $\\mu$m, aod$_{{500nm}}$={a:1.1f}'.format(t=oa['tau'][23],r=oa['ref'][6],a=0.6))
plt.plot(wa['wvl'],wa['rad'][0,:,0,6,23,1],label='$\\tau$={t}, r$_{{eff}}$={r} $\\mu$m, aod$_{{500nm}}$={a:1.1f}'.format(t=wa['tau'][23],r=wa['ref'][6],a=1.2))

plt.legend(loc=1,frameon=False)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Zenith Radiance [W/m$^2$ nm sr]')
plt.title('Changes in cloud transmitted light with overlying aerosol')
plt.axvspan(400,1800,color='white')
plot_greys()
plt.savefig(fp+'plot/Aerosol_above_cloud_model_spectra.png',dpi=600,transparent=True)


# In[75]:


plt.figure()
plt.plot(oa['wvl'],oa['rad'][0,:,0,6,17,1]/max(oa['rad'][0,:,0,6,17,1]),label='$\\tau$={t}, r$_{{eff}}$={r} $\\mu$m, aod$_{{500nm}}$={a:1.1f}'.format(t=oa['tau'][17],r=oa['ref'][6],a=0.6))
plt.plot(wa['wvl'],wa['rad'][0,:,0,6,17,1]/max(wa['rad'][0,:,0,6,17,1]),label='$\\tau$={t}, r$_{{eff}}$={r} $\\mu$m, aod$_{{500nm}}$={a:1.1f}'.format(t=wa['tau'][17],r=wa['ref'][6],a=1.2))
plt.plot(oa['wvl'],oa['rad'][0,:,0,6,23,1]/max(oa['rad'][0,:,0,6,23,1]),label='$\\tau$={t}, r$_{{eff}}$={r} $\\mu$m, aod$_{{500nm}}$={a:1.1f}'.format(t=oa['tau'][23],r=oa['ref'][6],a=0.6))
plt.plot(wa['wvl'],wa['rad'][0,:,0,6,23,1]/max(wa['rad'][0,:,0,6,23,1]),label='$\\tau$={t}, r$_{{eff}}$={r} $\\mu$m, aod$_{{500nm}}$={a:1.1f}'.format(t=wa['tau'][23],r=wa['ref'][6],a=1.2))


plt.legend(loc=1,frameon=False)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalized Zenith Radiance')
plt.title('Changes in cloud transmitted light with overlying aerosol')
plt.axvspan(400,1800,color='white')
plot_greys()
plt.savefig(fp+'plot/Aerosol_above_cloud_model_spectra_norm.png',dpi=600,transparent=True)


# In[ ]:





# # Plotting
# Present some fo the early plots here

# In[ ]:




