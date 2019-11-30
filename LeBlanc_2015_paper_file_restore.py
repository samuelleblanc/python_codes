#!/usr/bin/env python
# coding: utf-8

# # Intro
# Name:  
# 
#     LeBlanc_2015_paper_file_restore
# 
# Purpose:  
# 
#     Laod the cloud property retrieval data from the paper, for saving into a shareable format with Kokhanovsky
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     ict file save
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
#   - retrieved_kisq_20120525_v9.out  retrieved_kisq_20120806_v9.out  retrieved_kisq_20130110_v9.out
#   
#  Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, 2019-11-29

# # Import modules

# In[1]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
import scipy.io as sio
from path_utils import getpath
import matplotlib.pyplot as plt


# In[2]:


from load_utils import load_from_json, mat2py_time,toutc


# In[3]:


name = 'SSFR3'


# In[4]:


vv = 'v9'


# In[6]:


fp_dat = getpath('SSFR3')+'data/'


# # Load the files

# In[7]:


days = ['20120525','20120806','20130110']
labels = ['Liquid','Mixed-phase','Ice']


# In[25]:


times = [[15.0,16.0],[22.0,23.0],[17.5,19.5]]


# ## Load the 15 parameters retrievals

# In[8]:


f_liq = sio.idl.readsav(fp_dat+'retrieved_kisq_{d}_{v}.out'.format(d=days[0],v=vv))
f_mix = sio.idl.readsav(fp_dat+'retrieved_kisq_{d}_{v}.out'.format(d=days[1],v=vv))
f_ice = sio.idl.readsav(fp_dat+'retrieved_kisq_{d}_{v}.out'.format(d=days[2],v=vv))


# In[9]:


f_liq.keys()


# ## Load the slope retrievals

# In[10]:


sl_liq = sio.idl.readsav(fp_dat+'{d}_cld_parms3_1600nm.out'.format(d=days[0]))
sl_mix = sio.idl.readsav(fp_dat+'{d}_cld_parms3_1600nm.out'.format(d=days[1]))
sl_ice = sio.idl.readsav(fp_dat+'{d}_cld_parms3_ic_1600nm.out'.format(d=days[2]))


# In[11]:


sl_liq.keys()


# ## Load the 2wvl retrievals

# In[13]:


twv_liq = sio.idl.readsav(fp_dat+'{d}_cld_parms3_2wvl_1600nm.out'.format(d=days[0]))
twv_mix = sio.idl.readsav(fp_dat+'{d}_cld_parms3_2wvl_1600nm.out'.format(d=days[1]))
twv_ice = sio.idl.readsav(fp_dat+'{d}_cld_parms3_ic_2wvl_1600nm.out'.format(d=days[2]))


# In[14]:


twv_liq.keys()


# ## Load the spectra

# In[23]:


sp_liq = sio.idl.readsav(fp_dat+'{d}_calibspcs.out'.format(d=days[0],v=vv))
sp_mix = sio.idl.readsav(fp_dat+'{d}_calibspcs.out'.format(d=days[1],v=vv))
sp_ice = sio.idl.readsav(fp_dat+'{d}_calibspcs.out'.format(d=days[2],v=vv))


# In[24]:


sp_liq.keys()


# In[26]:


sp_liq['zspectra'].shape


# In[28]:


sp_liq['tmhrs'].shape


# # Verify the results and combine

# In[ ]:




