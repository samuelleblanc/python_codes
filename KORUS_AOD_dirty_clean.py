
# coding: utf-8

# # Info
# Name:  
# 
#     KORUS_AOD_dirty_clean
# 
# Purpose:  
# 
#     Try to get to the bottom of comparison in window dirtying events. 
#   
# Input:
# 
#     none
# 
# Output:
#    
#     plots
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
#     Written: Samuel LeBlanc,Santa Cruz, CA, 2016-11-08
#     

# # Prepare the python environment

# In[1]:

import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[7]:

from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth


# In[2]:

get_ipython().magic(u'matplotlib notebook')


# In[3]:

fp = 'C:\\Users\\sleblan2\\Research\\KORUS-AQ\\'


# # Load the ict files

# ## Load the 4STAR files

# In[5]:

days = ['20160501','20160503','20160504','20160506',
        '20160510','20160511','20160512','20160513',
        '20160514','20160515','20160516','20160517',
        '20160519','20160521','20160524','20160526',
        '20160529','20160530','20160601','20160602',
        '20160604','20160608','20160609','20160614']


# In[8]:

outaod_RA = []
outaod_head_RA = []
for d in days:
    fname_aod = fp+'aod_ict/korusaq-4STAR-AOD_DC8_{}_RA.ict'.format(d)
    tt,th = load_ict(fname_aod,return_header=True)
    outaod_RA.append(tt)
    outaod_head_RA.append(th)


# In[10]:

dumps(np.array([0.8,0.9,0.95,0.99]))


# In[ ]:



