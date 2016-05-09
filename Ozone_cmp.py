
# coding: utf-8

# # Info
# Name:  
# 
#     Ozone_cmp
# 
# Purpose:  
# 
#     To compare the ozone from the transit flight of KORUS AQ to OMI ozone amounts.
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
#     - load_utils.py : for loading OMI HDF5 files
#     - matplotlib
#     - numpy
#     - scipy : for saving and reading
#     - math
#     - os
#     - gc
#     - pdb
#     - datetime
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, OSAN AFB, Korea, 2016-05-03
#     Modified: 

# # Import initial modules and default paths

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.io as sio
import Sp_parameters as Sp
import hdf5storage as hs
import load_utils as lm


# In[16]:

from mpl_toolkits.basemap import Basemap,cm


# In[17]:

get_ipython().magic(u'matplotlib notebook')


# In[18]:

fp = 'C:/Users/sleblan2/Research/KORUS-AQ/'


# # Load the various data

# ## Load the OMI ozone data for the transit flight

# In[19]:

oz = hs.File(fp+'ozone/','OMI-Aura_L3-OMDOAO3e_2016m0426_v003-2016m0428t023351.he5','r')


# In[6]:

import h5py


# In[20]:

import tables


# In[21]:

oz = tables.openFile(fp+'ozone/'+'OMI-Aura_L3-OMDOAO3e_2016m0426_v003-2016m0428t023351.he5',mode='r',root_uep='/HDFEOS/GRIDS')


# In[22]:

oz.get_node


# In[25]:

ooo = oz.getNode('/ColumnAmountO3/Data Fields/ColumnAmountO3')


# In[30]:

plt.figure()
plt.contourf(ooo.read(),clevels=[200,225,250,275,300,325,350,375,400,425,450])


# In[24]:

oz.ColumnAmountO3


# In[11]:

oz.keys()


# In[10]:

oz = h5py.File(fp+'ozone/'+'OMI-Aura_L3-OMDOAO3e_2016m0426_v003-2016m0428t023351.he5','r')


# In[ ]:

oz['HDFEOS']


# In[ ]:



