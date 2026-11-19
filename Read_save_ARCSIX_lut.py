#!/usr/bin/env python
# coding: utf-8

# # Intro
# Name:  
# 
#     Read_save_ARCSIX_lut
# 
# Purpose:  
# 
#     Read the output files from the ARCSIX lut after libradran run. 
#     Based on inputs created by Prepare_ARCSIX_lut
# 
# Input:
# 
#     none
# 
# Output:
#    
#     Save matlab file of lut files
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - hdf5storage : for saving and reading
#     - Run_Libradtran
#     - os
#   
# Needed Files:
# 
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, 2026-11-16
#              Ported from Read_save_NAAMES_lut

# # Prepare the python environment

# In[1]:


import numpy as np
import Run_libradtran as RL
import hdf5storage as hs
import os


# In[2]:


fp = '/data2/ARCSIX/data/'
fp_atm = '/home/sam/libradtran/libRadtran-2.0.6/data/atmmod/'
fp_out = '/scratch/rtm/input/ARCSIX_v1/'
fp_rtm = '/scratch/rtm/' #input/ARCSIX_v1/'
fp_uvspec = '/home/sam/libradtran/libRadtran-2.0.6/bin/uvspec'


# # Setup the variables used for lut

# In[3]:


vv = 'v1'
sza = np.arange(60.0,63.0,0.5)
tau = np.array([0.1,0.2,0.3,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.3,2.6,3.0,3.5,4.0,4.5,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,22.0,25.0,27.0,30.0,35.0,40.0,50.0,60.0])
ref = np.append(np.append(np.arange(2,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,5.0)))
zout = [0.1,2.0,100.0]


# In[6]:


fp_out = os.path.join(fp_rtm,'output','ARCSIX_%s'%vv)


# In[7]:


dat = RL.read_lut(fp_out,zout=zout,tau=tau,ref=ref,sza=sza,
                  phase=['wc','ic'],
                  fmt='lut_sza{sza:04.1f}_tau{tau:06.2f}_ref{ref:04.1f}_{phase}_w{iwvl:1d}.dat',
                  split_wvl=True)


# In[9]:


print('Saving matlab file')


# In[12]:


import scipy.io as sio


# In[13]:


sio.savemat(fp+'{}_ARCSIX_lut.mat'.format(vv),dat)


# In[14]:


fp


# In[15]:


fp+'{}_ARCSIX_lut.mat'.format(vv)


# In[ ]:




