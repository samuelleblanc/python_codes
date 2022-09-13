#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     To Calculate aerosol optical properties for aerosol above cloud reitrewvals using MODIS and VIIRS
#     Using the wavelengths: 
#     0.47, 0.55, 0.67, 0.86, 1.24, 2.1µm
#     
#     - Using the retrieved size distributions
#     - retrieval results of refractive index (imaginary and real) at wavelengths: 400, 500, 675, 870, 995 nm
# 
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
#     - libradtran
# 
# Needed Files:
#   - netcdf of aeroinv from ORACLES
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-09-12
#     Modified:
# 

# # Prepare python environment

# In[3]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'notebook')
import os


# In[7]:


name = 'sunsat_ORACLES2016'
vv = 'v1'
fp = getpath(name)


# In[17]:


fp_lib = getpath('uvspec_bin')


# In[27]:


fp_rtm0 = getpath('rtm')
fp_rtm = fp_rtm0 +'TASNPP_mie/'


# In[28]:


if not os.path.exists(fp_rtm): 
    os.mkdir(fp_rtm)


# # Load files

# In[8]:


f = fp + 'data_archival/AEROINV_nc_R0/NC_FORMAT_NETCDF4_CLASSIC/4STAR-aeroinv_P3_2016_R0.nc'


# In[9]:


ae,ae_dict = lu.load_netcdf(f,everything=True)


# In[14]:


len(ae[b'time'])


# # Run through each retrieval and make input files

# ## Prep functions for printing size distribution

# From the libradtran documentation for mie size distribution file:  
# >Specify a two column file, r [micron], dn(r)/dr, which describes a size distribution

# In[34]:


def print_size_dist_file(fname,r,dnr):
    with open(fname,'w') as f:
        for ir,rr in list(enumerate(r)):
            f.write('{:3.10f} {:3.10f}\n'.format(rr,dnr[ir]))


# In[20]:


ae_dict[b'radius']


# In[21]:


ae_dict[b'psd']


# In[25]:


def convert_dvlnr_to_dndr(psd,r):
     # All about that conversion from the volume size distribution of dV(r)/dln(r) to number size distribution dN(r)/dr
    Nr = psd/(4.0*np.pi/3.0)
    for i,rr in list(enumerate(r)):
        Nr[i] = Nr[i]/rr**4.0
    return Nr


# In[37]:


print_size_dist_file(fp_rtm+'mie_tester.psd',ae[b'radius'],convert_dvlnr_to_dndr(ae[b'psd'][0,:],ae[b'radius']))


# ## Prep function for printing index of refraction

# In[ ]:




