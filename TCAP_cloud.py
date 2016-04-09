
# coding: utf-8

# # Info
# Name:  
# 
#     TCAP_cloud
# 
# Purpose:  
# 
#     To compare the various cloud properties retrieved via different methods from TCAP.
#     Looking at the Feb, 19th, 2013 case
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
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - load_modis.py : for loading modis files
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
#   - 20151117_zen_cld_retrieved.mat: cloud retrieval file
#   - MYD06_L2.A2015321.1540.006.2015322185040.hdf: MODIS file
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, NASA Ames, Santa Cruz, CA, 2016-04-06
#     Modified: 

# # Import initial modules and default paths
# 

# In[2]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.io as sio
import Sp_parameters as Sp
import hdf5storage as hs
import load_modis as lm
import os


# In[3]:

from mpl_toolkits.basemap import Basemap,cm


# In[4]:

import write_utils as wu
import plotting_utils as pu
import map_utils as mu
import Run_libradtran as RL


# In[5]:

get_ipython().magic(u'matplotlib notebook')


# In[6]:

# set the basic directory path
fp='C://Users//sleblan2//Research//TCAP//'


# # Load the various data

# ## Verify the retrievals from 4STAR

# ### Build the proper lut from the idl out file

# In[7]:

vv = 'v1'


# In[8]:

# load the idl save file containing the modeled radiances
s = sio.idl.readsav(fp+'model/sp_v1_20130219_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp,   wvl,  z,  re,  ta)'


# ### Compare to lut from NAAMES

# In[21]:

fp_lut_mat = 'C:\\Users\\sleblan2\\Research\\NAAMES\\lut\\v2_NAAMES_lut.mat'
print('Loading the lut file:{}'.format(fp_lut_mat))
if not os.path.isfile(fp_lut_mat):
    print('File {} does not exist'.format(fp_lut_mat))
    raise IOError('LUT File not found: {}'.format(fp_lut_mat))
luts = hs.loadmat(fp_lut_mat)


# In[22]:

luts.keys()


# In[23]:

luts['irr_up'].shape


# In[27]:

luts['rad'].shape, s.sp.shape


# In[28]:

s.sp[:,:,:,:,:,np.newaxis].shape


# In[38]:

s.ref


# In[39]:

luts['ref']


# ### reshape TCAP LUT to fit NAAMES format and save new file

# In[54]:

tcap_lut = {u'tau':s.tau,
            u'rad':s.sp[:,:,:,:,:,np.newaxis],
            u'sza':np.array(s.sza)[np.newaxis],
            u'irr_dn_diff':s.sp_irrdn[:,:,:,:,:,np.newaxis]*0.0,
            u'irr_dn':s.sp_irrdn[:,:,:,:,:,np.newaxis],
            u'irr_up':s.sp_irrup[:,:,:,:,:,np.newaxis],
            u'zout':s.zout,
            u'wvl':s.zenlambda,
            u'phase':['wc','ic'],
            u'ref':s.ref}


# In[55]:

hs.savemat(fp+'model\\{}_TCAP_lut.mat'.format(vv),tcap_lut)


# ### Run the retrieval via the run_zen_cld_retrieval command line
# in the python_codes directory:
# 
# > python Run_zen_cld_retrieval.py -lut C:\Users\sleblan2\Research\TCAP\model\v1_TCAP_lut.mat C:\Users\sleblan2\Research\TCAP\20130219starzen.mat -o C:\Users\sleblan2\Research\TCAP\plots\ 

# ## Check the 4STAR data

# In[12]:

ss = sio.loadmat(fp+'20130219starzen.mat')


# In[13]:

ss.keys()


# In[14]:

ss['t'].shape


# In[15]:

ss['rad'].shape


# In[16]:

plt.figure()
plt.plot(ss['t'],ss['rad'][:,400])


# ## Load the retrieved cloud properties from 4STAR

# In[27]:

star = hs.loadmat(fp+'4STAR//20130219_zen_cld_retrieved.mat')


# In[28]:

star.keys()


# In[29]:

plt.figure()
plt.plot(star['tau'])


# In[ ]:



