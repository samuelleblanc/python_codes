
# coding: utf-8

# # Introduction
# 
# Name:  
# 
#     NAAMES_2016_AATS_check
# 
# Purpose:  
# 
#     Testing out the archived files from AATS
#   
# Input:
# 
#     none at command line
# 
# Output:
#    
#     figures
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - math
#     - pdb
#     - datetime
#     - load_utils
#     - Run_fuliou
#   
# Needed Files:
# 
#   - AATS NAAMES AOD ict files
#   
#   
# Modification History:
# 
#     Wrtten: Samuel LeBlanc, NASA Ames, 2018-05-30
#     Modified:

# # Load the defaults modules and setup

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
from Sp_parameters import smooth
from linfit import linfit
from path_utils import getpath
from plotting_utils import make_boxplot


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[7]:


fp = getpath('NAAMES',path='/mnt/c/Users/sleblanc/Research/NAAMES',make_path=True)
fp


# # Load the AOD ict files

# In[9]:


vv = 'R0'
d = '20160527'


# In[11]:


fname_aod = fp+'aod_ict/NAAMES_2/NAAMES-AATS-AOD_C130_{}_{vv}.ict'.format(d,vv=vv)
tt,th = load_ict(fname_aod,return_header=True)


# # Plot out the data

# In[12]:


plt.figure()
plt.plot(tt['Start_UTC'],tt['AOD0499'])


# In[14]:


prof = np.array([46296,48896])/3600.0
prof

