
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_SSFR_cloud_retrieval
# 
# Purpose:  
# 
#     Run the cloud retrieval on SSFR data from ORACLES 2016 and 2017, but only on the flagacaod times.
#   
# Input:
# 
#     none
# 
# Output:
#    
#     plots
#     mat save files
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
#     Written: Samuel LeBlanc,Santa Cruz, CA, 2018-06-29

# # Prepare the python environment

# In[1]:


import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


import load_utils as lu
import Sp_parameters as Sp


# In[4]:


import hdf5storage as hs
from path_utils import getpath
from write_utils import nearest_neighbor
from tqdm import tqdm_notebook as tqdm
import math


# In[207]:


from scipy import interpolate


# In[5]:


fp = getpath('ORACLES')
fp


# # Load files

# In[9]:


days = ['20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927']


# ## Load the SSFR ict files for 2016

# In[26]:


ssfr_a, ssfr_ah = [],[]
for d in days:
    sf,sfh = lu.load_ict(fp+'data_other/ssfr/SSFR_P3_{}_R1.ict'.format(d),return_header=True)
    ssfr_a.append(lu.recarray_to_dict(sf))
    ssfr_ah.append(sfh)


# In[12]:


ssfr_ah[0]


# In[13]:


ssfr_a[0]


# ## Load the 4STAR files with flagacaod

# In[25]:


star_a, star_ah = [],[]
for d in days:
    sf,sfh = lu.load_ict(fp+'aod_ict/v8/4STAR-AOD_P3_{}_R3.ict'.format(d),return_header=True)
    star_a.append(lu.recarray_to_dict(sf))
    star_ah.append(sfh)


# In[18]:


ssfr_a[3]['Start_UTC'][100]


# ## Get the flagacaod on the timescale of the ssfr measurements

# In[243]:


for i,d in enumerate(days):
    fa = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['flag_acaod'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['flagacaod'] = fa
    am = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['amass_aer'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['airmass'] = am
    ssfr_a[i]['sza'] = np.arccos(1.0/am)*180.0/np.pi
    aod = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['AOD0501'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['AOD_500'] = aod
    a2 = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['AOD_polycoef_a2'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['a2'] = a2
    a1 = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['AOD_polycoef_a1'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['a1'] = a1
    a0 = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['AOD_polycoef_a0'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['a0'] = a0


# In[220]:


ssfr_a[0]['flagacaod'].shape,ssfr_a[0]['Start_UTC'].shape


# ## Load the LUT for 2wvl reflectance retrieval

# In[30]:


lut = hs.loadmat(fp+'rtm/v5_irr_ORACLES_lut.mat')


# In[31]:


lut.keys()


# ## Combine into one array

# In[249]:


nm = ssfr_a[1].keys()


# In[250]:


ar = {}
for n in ssfr_a[1].keys():
    ar[n] = np.array([])


# In[251]:


ar['days'] = np.array([])


# In[252]:


for i,d in enumerate(days):
    ar['days'] = np.append(ar['days'],np.zeros_like(ssfr_a[i]['Start_UTC'])+i)
    for n in nm:
        try:
            ar[n] = np.append(ar[n],ssfr_a[i][n])
        except:
            print 'problem with :'+n
            ar[n] = np.append(ar[n],ssfr_a[i]['Start_UTC']*0)


# # Format the LUT and data for retrievals

# In[253]:


class so:
    pass


# ## Set up the data

# In[254]:


ar['meas'] = so
ar['meas'].sza = ar['sza']
ar['meas'].Rvis = ar['UP500']/ar['DN500']
ar['meas'].Rnir = ar['UP1650']/ar['DN1650']
ar['meas'].utc = ar['Start_UTC']


# In[255]:


# filter out the bad data. 
bad = (ar['meas'].Rvis > 1.0) & (ar['flagacaod']==0) & (ar['meas'].Rnir > 1.0)
ar['meas'].Rvis[bad] = np.nan
ar['meas'].Rvis[bad] = np.nan


# In[256]:


igood = np.where((np.isfinite(ar['meas'].Rvis)) & (ar['meas'].Rvis > 0.0) & (np.isfinite(ar['meas'].Rnir)) & (ar['meas'].Rnir > 0.0) & (ar['flagacaod']==1))[0]


# ## set up the LUT

# In[257]:


lut.keys()


# In[124]:


lut['tau'].shape, lut['ref'].shape, lut['sza'].shape, lut['irr_dn'].shape, lut['wvl'].shape, lut['zout'], lut['phase']


# In[121]:


lut['Rvis'] = np.zeros([23,27,48])
lut['Rnir'] = np.zeros([23,27,48])


# In[122]:


for ir,r in enumerate(lut['ref']):
    for it,t in enumerate(lut['tau']):
        for iz,s in enumerate(lut['sza']):
            lut['Rvis'][ir,it,iz] = lut['irr_up'][0,0,1,ir,it,iz]/lut['irr_dn'][0,0,1,ir,it,iz]
            lut['Rnir'][ir,it,iz] = lut['irr_up'][0,1,1,ir,it,iz]/lut['irr_dn'][0,1,1,ir,it,iz]


# In[126]:


lut['sza']


# ### Make a hires version of the LUT

# In[206]:


lut['tau_hi'] = np.hstack([np.arange(0,25,0.5),np.arange(25,50,1),np.arange(50,102.5,2.5)])
lut['ref_hi'] = np.hstack([np.arange(0,15,0.25),np.arange(15,30.5,0.5)])


# In[208]:


len(lut['tau_hi']), len(lut['ref_hi'])


# In[209]:


lut['Rvis_hi'] = np.zeros([91,96,48])
lut['Rnir_hi'] = np.zeros([91,96,48])


# In[210]:


for i,z in enumerate(lut['sza']):
    fv = interpolate.RectBivariateSpline(lut['ref'],lut['tau'],lut['Rvis'][:,:,i],kx=1,ky=1)
    lut['Rvis_hi'][:,:,i] = fv(lut['ref_hi'],lut['tau_hi'])
    fn = interpolate.RectBivariateSpline(lut['ref'],lut['tau'],lut['Rnir'][:,:,i],kx=1,ky=1)
    lut['Rnir_hi'][:,:,i] = fn(lut['ref_hi'],lut['tau_hi'])


# # Run the retrieval

# In[258]:


ar['tau'], ar['ref'] = np.zeros_like(ar['sza'])*np.nan,np.zeros_like(ar['sza'])*np.nan


# In[259]:


ar['ki'] = np.zeros_like(ar['sza'])


# In[260]:


pbar = tqdm(total=len(ar['sza']))
for i,s in enumerate(ar['sza']):
    pbar.update()
    if (s>73.0) | (np.isnan(s)):
        continue
    if not i in igood:
        continue
    isza = np.argmin(np.abs(lut['sza']-s))
    ki = (ar['meas'].Rvis[i]-lut['Rvis_hi'][:,:,isza]/ar['meas'].Rvis[i])**2+(ar['meas'].Rnir[i]-lut['Rnir_hi'][:,:,isza]/ar['meas'].Rnir[i])**2
    kimin = np.unravel_index(np.nanargmin(ki),ki.shape)
    ar['ki'][i] = np.nanmin(ki)
    ar['tau'][i],ar['ref'][i] = lut['tau_hi'][kimin[1]],lut['ref_hi'][kimin[0]]


# # Plot the retrieval results

# In[215]:


plt.figure()
plt.plot(ar['tau'],'.')


# In[216]:


plt.figure()
plt.plot(ar['ref'],'.')


# In[261]:


len(np.where(np.isfinite(ar['ref']))[0])


# In[262]:


len(ar['ref'])


# # Save the retrieved output

# In[263]:


out = {}


# In[264]:


ar.keys()


# In[265]:


out['tau'] = ar['tau']
out['ref'] = ar['ref']
out['sza'] = ar['sza']
out['aod'] = ar['AOD_500']
out['days'] = ar['days']
out['utc'] = ar['Start_UTC']
out['lat'] = ar['LAT']
out['lon'] = ar['LON']
out['a0'],out['a1'],out['a2'] = ar['a0'],ar['a1'],ar['a2']


# In[266]:


fp


# In[267]:


hs.savemat(fp+'data_other/ssfr_2016_retrieved_COD.mat',out)


# In[268]:


out['days']

