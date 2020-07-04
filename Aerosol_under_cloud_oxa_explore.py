#!/usr/bin/env python
# coding: utf-8

# # Intro
# Name:  
# 
#     Aerosol_under_cloud_oxa_explore
# 
# Purpose:  
# 
#     Explore KORUS-AQ data to see if there are signals of aerosol to be exploited within the cloudy/cirrusy scenes. Relate to Oxygen-A band
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures and save files...
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - plotting_utils (user defined plotting routines)
#     - hdf5storage
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - starsun for KORUS-AQ
#   
#  Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, Humble Sea, 2019-03-15

# # Load the python modules and setup paths
# 

# In[92]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
import scipy.io as sio
import matplotlib.pyplot as plt
from mpltools import color


# In[93]:


import load_utils as lu
from path_utils import getpath


# In[94]:


get_ipython().magic(u'matplotlib notebook')


# In[95]:


import scipy.stats as st
import Sun_utils as su
import plotting_utils as pu


# In[96]:


import Sp_parameters as sp
from scipy import interpolate


# In[97]:


from linfit import linfit


# In[98]:


fp = getpath('KORUS-AQ')


# In[99]:


daystr = '20160519'


# # Load the starsun files

# ## Load the starsun file

# In[12]:


s = sio.loadmat(fp+'data/{}starsun.mat'.format(daystr))


# In[14]:


kk = s.keys()
kk.sort()
kk


# In[18]:


s['utc'] = lu.toutc(lu.mat2py_time(s['t']))


# In[21]:


plt.plot(s['utc'],s['Alt'],'.')


# In[22]:


tprofile = [23.5,23.8]


# In[43]:


ip = (s['utc']>tprofile[0]) & (s['utc']<tprofile[1])


# # Calculate the angstrom exponent

# In[36]:


s['w'].shape


# In[37]:


ja = np.argmin(abs(s['w']-0.470))


# In[39]:


je = np.argmin(abs(s['w']-0.875))


# In[41]:


s['w'][0,ja],s['w'][0,je]


# In[42]:


s['angs_470_865'] = []
for i,u in enumerate(s['utc']):
    to = np.log(s['tau_aero'][i,je])-np.log(s['tau_aero'][i,ja])
    bo = np.log(s['w'][0,je]*1000.0)-np.log(s['w'][0,ja]*1000.0)
    s['angs_470_865'].append(to/bo*-1.0)
s['angs_470_865'] = np.array(s['angs_470_865'])


# In[44]:


plt.figure()
plt.plot(s['Alt'][ip],s['angs_470_865'][ip])


# # Calculate the oxygen-a band depth

# In[52]:


plt.figure()
color.cycle_cmap(len(s['Alt'][ip]),cmap=plt.cm.viridis,ax=plt.gca())
plt.plot(s['w'].T,s['tau_aero'][ip,:].T)
plt.ylim(0,0.6)
plt.xlim(0.35,0.8)
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
scalarmap.set_array(s['Alt'][ip])
cba = plt.colorbar(scalarmap)
cba.set_label('Altitude [m]')


# In[61]:


plt.figure()
plt.plot(s['w'][0,:],s['tau_aero'][ip,:][0,:],'.-')
plt.ylim(0,0.6)
plt.xlim(0.75,0.775)


# In[62]:


y0 = np.argmin(abs(s['w']-0.756))
y1 = np.argmin(abs(s['w']-0.772))


# In[63]:



y0,y1


# In[72]:


s['oxa_flat'] = []
s['oxa_delta'] = []
for i,u in enumerate(s['utc']):
    fx = interpolate.interp1d(s['w'][0,[y0,y1]],s['tau_aero'][i,[y0,y1]])
    s['oxa_flat'].append(fx(s['w'][0,y0:y1]))
    s['oxa_delta'].append(np.nansum(s['tau_aero'][i,y0:y1]-s['oxa_flat']))


# In[71]:


s['w'][0,y0:y1]


# In[87]:


plt.figure()
plt.scatter(s['oxa_delta'],s['angs_470_865'],c=s['Alt'])
plt.xlim(-50000,200000)
plt.colorbar()


# In[90]:


plt.figure()
plt.plot(s['oxa_delta']/s['Alt'],s['angs_470_865'],'.')
#plt.xlim(-50000,200000)
plt.colorbar()


# In[ ]:




