
# coding: utf-8

# Name:  
# 
#     Read_save_ORACLES_lut_reflectance
# 
# Purpose:  
# 
#     Read the output files from the ORACLES lut after libradran run for creating irradiance reflectances. 
#     Based on inputs created by Prepare_ORACLES_lut_reflectance
# 
# Calling Sequence:
# 
#     python Read_save_ORACLES_lut_reflectance
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
#   - .out uvspec files
#     
# History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2018-06-27

# # Prepare the python environment

# In[1]:


import numpy as np
import Run_libradtran as RL
import hdf5storage as hs
import os


# In[2]:


if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\ORACLES\\'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\ORACLES\\rtm\\'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/ORACLES/'
    fp_rtm = '/nobackup/sleblan2/rtm/'
else:
    raise Exception


# # Setup the variables used for lut

# In[3]:


vv = 'v5_irr'


# In[ ]:


# try to read from the saved version file
from load_utils import load_from_json
try:
    d = load_from_json(fp+'ORACLES_lut_{}.txt'.format(vv))
    sza = d['lut']['sza']
    tau = d['lut']['tau']
    ref = d['lut']['ref']
    zout = d['geo']['zout']
    fmt = d['lut']['format']
    mu = np.round(1.0/np.cos(sza*np.pi/180.0))
    use_json = True
except ValueError: # not a json file try old way
    use_json = False
    fmt='lut_irr_sza{sza:02.0f}_tau{tau:06.2f}_ref{ref:04.1f}_{phase}_w{iwvl:1d}.dat'
    if vv=='v1':
        mu = np.arange(1.05,4.0,0.2)
        sza = np.round(np.arccos(1.0/mu)*180.0/np.pi)
        tau = np.array([0.1,0.2,0.3,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0,
               6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0,50.0,
               60.0,80.0,100.0])
        ref = np.append(np.append(np.arange(2,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,2.5)))
        zout = [0.2,1.5,100.0]
    elif vv=='v3':
        mu = np.arange(1.05,4.0,0.2)
        sza = np.round(np.arccos(1.0/mu)*180.0/np.pi)
        tau = np.array([0.1,0.2,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,
               6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0,50.0,
               60.0,80.0,100.0])
        ref = np.append(np.append(np.arange(1,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,2.5)))
        zout = [0.2,1.5,100.0]    


# In[7]:


fp_out = os.path.join(fp_rtm,'output','%s_ORACLES'%vv)


# In[ ]:


dat = RL.read_lut(fp_out,zout=zout,tau=tau,ref=ref,sza=sza,
                  phase=['wc'],
                  fmt=fmt,
                  split_wvl=True)


# In[ ]:


if use_json:
    dat['lut_details'] = d


# In[ ]:


print 'Saving matlab file'


# In[10]:


hs.savemat(fp+'{}_ORACLES_lut.mat'.format(vv),dat)

