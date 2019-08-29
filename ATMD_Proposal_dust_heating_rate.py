
# coding: utf-8

# # Info
# Name:  
# 
#     ATDM Proposal dust heating rate
# 
# Purpose:  
# 
#     Make some figures for the ATMD proposal
#     Focus on the heating rate profiles obtained during ORACLES 2018 transit back, near Cabo Verde, with SSFR
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
#     - 
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - for_Sam_20181025.out
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-08-23
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
import scipy.io as sio


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


fp =getpath('ORACLES')


# # Load files

# ## Load the SSFR files

# In[4]:


ssfr = sio.idl.readsav(fp+'data_other_2018/SSFR/for_Sam_20181025_SSFR.out')


# In[5]:


ssfr.keys()


# Now interpolate the nadir spectra to the zenith wavelengths

# In[19]:


def interp_spline(xold, yold, ynew):
    uv = UnivariateSpline(xold,yold,ext=0,k=1)
    return uv(ynew)


# In[21]:


ssfr['nadspectra'].shape


# In[25]:


ssfr['nadlambda'][0:10]


# In[26]:


ssfr['nspectra'] = np.array([interp_spline(ssfr['nadlambda'],ssfr['nadspectra'][i,:],ssfr['zenlambda']) for i in xrange(len(ssfr['utc']))])


# In[27]:


ssfr['nspectra'].shape


# In[28]:


ssfr['f_net'] = ssfr['zenspectra'] - ssfr['nadspectra']


# ## Load the 4STAR file

# In[6]:


s = lu.load_ict(fp+'aod_ict_2018/4STAR-AOD_P3_20181025_R1.ict')


# In[7]:


sp = sio.loadmat(fp+'data_2018/4STAR_20181025starsun.mat')


# ## Load the merge file

# In[13]:


mrg,mrg_head = lu.load_netcdf(fp+'data_other_2018/mrg1_P3_20181025_R13.nc',everything=True)


# # Build up the plots to select the region
# Do some of the calculations to the data here

# Get the location of the dust profile

# In[14]:


plt.figure()
plt.plot(s['Start_UTC'],s['GPS_Alt'],'.')


# In[10]:


plt.figure()
plt.plot(s['Latitude'],s['GPS_Alt'],'.')


# In[15]:


pfl = [13.565, 13.95]


# In[17]:


i_ssfr = (ssfr['utc']>pfl[0]) & (ssfr['utc']<pfl[1])


# In[30]:


ssfr['alt'] = interp_spline(s['Start_UTC'],s['GPS_Alt'],ssfr['utc'])


# In[33]:


plt.figure()
plt.pcolor(ssfr['zenlambda'],ssfr['alt'][i_ssfr],ssfr['f_net'][i_ssfr,:])


# In[36]:


ssfr['alt'][i_ssfr]


# # Plotting
# Present some fo the early plots here
