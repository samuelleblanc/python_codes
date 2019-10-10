
# coding: utf-8

# # Intro
# Name:  
# 
#     Read_save_ARISE_lut
# 
# Purpose:  
# 
#     Read the output files from the NAAMES lut after libradran run. 
#     Based on inputs created by Prepare_NAAMES_lut
# 
# Calling Sequence:
# 
#     python Read_save_ARISE_lut
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
#     Written: Samuel LeBlanc, NASA Ames, 2016-10-19
#              Ported from Read_save_NAAMES_lut

# # Prepare the python environment

# In[9]:


import numpy as np
import Run_libradtran as RL
import hdf5storage as hs
import os


# In[5]:


if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\ARISE\\'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\ARISE\\rtm\\'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/ARISE/'
    fp_rtm = '/nobackup/sleblan2/rtm/'
else:
    raise Exception


# # Setup the variables used for lut

# In[6]:


vv = 'v3'
mu = np.arange(2.7,4.0,0.1)
sza = np.round(np.arccos(1.0/mu)*180.0/np.pi)
tau = np.array([0.1,0.2,0.3,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.3,2.6,3.0,3.5,4.0,4.5,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,22.0,25.0,27.0,30.0,35.0,40.0,50.0])
ref = np.append(np.append(np.arange(2,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,3.0)))
zout = [0.2,2.0,100.0]


# In[7]:


fp_out = os.path.join(fp_rtm,'output','%s_ARISE'%vv)


# In[ ]:


dat = RL.read_lut(fp_out,zout=zout,tau=tau,ref=ref,sza=sza,
                  phase=['wc','ic'],
                  fmt='lut_sza{sza:02.0f}_tau{tau:06.2f}_ref{ref:04.1f}_{phase}_w{iwvl:1d}.dat',
                  split_wvl=True)


# In[ ]:


print 'Saving matlab file'


# In[10]:


hs.savemat(fp+'{}_ARISE_lut.mat'.format(vv),dat)

