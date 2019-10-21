
# coding: utf-8

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
import load_modis as lm
import os


# In[2]:


hs.__version__


# In[3]:


from mpl_toolkits.basemap import Basemap,cm


# In[4]:


get_ipython().magic(u'matplotlib notebook')


# In[6]:


fp_lut_mat = 'C:\\Users\\sleblan2\\Research\\NAAMES\\lut\\v2_NAAMES_lut.mat'
print('Loading the lut file:{}'.format(fp_lut_mat))
if not os.path.isfile(fp_lut_mat):
    print('File {} does not exist'.format(fp_lut_mat))
    raise IOError('LUT File not found: {}'.format(fp_lut_mat))
luts = hs.loadmat(fp_lut_mat)

