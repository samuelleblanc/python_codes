
# coding: utf-8

# Name:  
# 
#     Read_save_NAAMES_lut
# 
# Purpose:  
# 
#     Read the output files from the NAAMES lut after libradran run. 
#     Based on inputs created by Prepare_NAAMES_lut
# 
# Calling Sequence:
# 
#     python Read_save_NAAMES_lut
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
#     - scipy : for saving and reading
#     - Run_Libradtran
# 
#   
# Needed Files:
# 
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc, NASA Ames, 2015-10-14

# # Prepare the python environment

# In[2]:

import numpy as np
import Run_libradtran as RL
import scipy.io as sio
import os


# In[5]:

if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\NAAMES\\'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\NAAMES\\rtm\\'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/NAAMES/'
    fp_rtm = '/nobackup/sleblan2/rtm/'
else:
    raise Exception


# # Setup the variables used for lut

# In[6]:

vv = 'v2'
mu = np.arange(1.05,4.0,0.2)
sza = np.round(np.arccos(1.0/mu)*180.0/np.pi)
tau = np.array([0.1,0.2,0.3,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0,40.0,50.0,
       60.0,70.0,80.0,100.0])
ref = np.append(np.append(np.arange(2,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,2.5)))
zout = [0.2,3.0,100.0]


# In[7]:

fp_out = os.path.join(fp_rtm,'output','%s_NAAMES'%vv)


# In[ ]:

dat = RL.read_lut(fp_out,zout=zout,tau=tau,ref=ref,sza=sza,
                  phase=['wc','ic'],
                  fmt='lut_sza{sza:02d}_ref{tau:02.1f}_tau{ref:03.1f}_{phase}_w{iwvl:1d}.dat',
                  split_wvl=True)


# In[ ]:

print 'Saving matlab file'


# In[ ]:

sio.savemat(fp+'{}_NAAMES_lut.mat'.format(vv),dat)

