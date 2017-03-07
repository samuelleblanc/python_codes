
# coding: utf-8

# Name:  
# 
#     Read_save_KORUS_lut
# 
# Purpose:  
# 
#     Read the output files from the KORUS lut after libradran run. 
#     Based on inputs created by Prepare_KORUS_lut
# 
# Calling Sequence:
# 
#     python Read_save_KORUS_lut
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
#     Written: Samuel LeBlanc, Santa Cruz, 2017-02-03

# # Prepare the python environment

# In[9]:

import numpy as np
import Run_libradtran as RL
import hdf5storage as hs
import os
import scipy.io as sio


# In[ ]:

name = 'KORUS'


# In[5]:

if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\{}\\'.format(name)
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\{}\\rtm\\'.format(name)
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/{}/'.format(name)
    fp_rtm = '/nobackup/sleblan2/rtm/'
else:
    raise Exception


# # Setup the variables used for lut

# In[6]:

vv = 'v3'


# In[ ]:

# try to read from the saved version file
from load_utils import load_from_json, deep_convert_dict
try:
    d = load_from_json(fp+'{name}_lut_{vv}.txt'.format(vv=vv,name=name))
    if 'lut' in d.keys():
        sza = d['lut']['sza']
        tau = d['lut']['tau']
        ref = d['lut']['ref']
        fmt = d['lut']['format']
    elif 'lut_details' in d.keys():
        sza = d['lut_details']['sza']
        tau = d['lut_details']['tau']
        ref = d['lut_details']['ref']
        fmt = d['lut_details']['format']
    else:
        raise ValueError
    zout = d['geo']['zout']
    mu = np.round(1.0/np.cos(sza*np.pi/180.0))
    use_json = True
except ValueError: # not a json file try old way
    print '*** LUT definition file problem! Using predefined lut description'
    use_json = False
    fmt='lut_sza{sza:02.0f}_tau{tau:06.2f}_ref{ref:04.1f}_{phase}_w{iwvl:1d}.dat'
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

fp_out = os.path.join(fp_rtm,'output','{vv}_{name}'.format(vv=vv,name=name))


# In[ ]:

dat = RL.read_lut(fp_out,zout=zout,tau=tau,ref=ref,sza=sza,
                  phase=['wc','ic'],
                  fmt=fmt,
                  split_wvl=True)


# In[ ]:

dat = deep_convert_dict(dat)
if use_json:
    dat['lut_details'] = deep_convert_dict(d)


# In[ ]:

print 'Saving matlab file'


# In[10]:

try:
    try:
        hs.savemat(fp+'{vv}_{name}_lut.mat'.format(vv=vv,name=name),dat)
    except:
        sio.savemat(fp+'{vv}_{name}_lut.mat'.format(vv=vv,name=name),dat)
except:
    import pdb
    pdb.set_trace()
    sio.savemat(fp+'{vv}_{name}_lut.mat'.format(vv=vv,name=name),dat)

