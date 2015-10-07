
# coding: utf-8

# Name:  
# 
#     Prepare_NAAMES_lut
# 
# Purpose:  
# 
#     Create the input libradtran files for creating a lut of low clouds to be used in NAAMES operational cloud retrievals
# 
# Calling Sequence:
# 
#     python Prepare_NAAMES_lut
#   
# Input:
# 
#     none
# 
# Output:
#    
#     input files for libradtran 2.0 (uvspec) 
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - mplt_toolkits for basemap, map plotting
#     - pdb
#     - datetime
# 
#   
# Needed Files:
# 
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc, NASA Ames, 2015-10-06
#     

# # Prepare the python environment

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
import numpy as np
import scipy.io as sio
import os
from mpl_toolkits.basemap import Basemap


# In[45]:

import Run_libradtran as RL
reload(RL)


# In[15]:

if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\NAAMES\\'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/NAAMES/'
else:
    raise Exception


# # Setup the variables used to create the lut

# In[24]:

vv = 'v1'


# In[16]:

sza = np.arange(40,90,2)


# In[19]:

sza


# In[42]:

tau = [1,2]
ref = [2,3]


# In[43]:

geo = {'lat':47.6212167,
       'lon':52.74245,
       'doy':322,
       'zout':[0.2,3.0,100.0]}
aero = {} # none
cloud = {'ztop':1.5,
         'zbot':0.5,
         'write_moments_file':False}
source = {'wvl_range':[350,1750],
          'source':'solar',
          'integrate_values':False,
          'run_fuliou':False,
          'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/'}
albedo = {'create_albedo_file':False,
          'sea_surface_albedo':True,
          'wind_speed':14.0}


# In[46]:

RL.print_version_details(fp+'NAAMES_lut_%s.txt'%vv,vv,geo=geo,
                         aero=aero,cloud=cloud,source=source,albedo=albedo,tau=tau,ref=ref,sza=sza)


# In[ ]:



