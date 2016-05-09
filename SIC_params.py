
# coding: utf-8

# Name:  
# 
#     SIC_params
# 
# Purpose:  
# 
#     Python script to looking at Shannon Information Content of the zenith cloud retrievals parameters
# 
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
#     - matplotlib
#     - scipy
#     - load_utils
#     - math
#     - Sp_parameters
#   
# Needed Files:
# 
#   - retr_day_pss_v4.out
#   - retr_day_sub_v4.out
#     
#  Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, CA, 2015-08-12
#             

# In[2]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
get_ipython().magic(u'matplotlib inline')
import numpy as np, h5py
import scipy.io as sio
import scipy
import math, os, IPython
import Sp_parameters as Sp
import load_utils as lm
IPython.InteractiveShell.cache_size = 0
# set the basic directory path
fp='C:/Users/sleblan2/Research/SSFR3/'


# In[3]:

fppss = fp+'retrieved/cst/retr_day_pss_v4.out'
fpsub = fp+'retrieved/cst/retr_day_sub_v4.out'


# In[4]:

pss = sio.idl.readsav(fppss)
sub = sio.idl.readsav(fpsub)


# In[6]:

pss.keys()


# In[7]:

sub.keys()


# In[9]:

pss['hrt'].shape


# In[ ]:



