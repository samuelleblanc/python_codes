#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np, h5py
import plotly.plotly as py
import scipy.io as sio
import math
import os
import warnings
warnings.simplefilter('ignore', np.RankWarning)
import Sp_parameters as Sp
py.sign_in("samuelleblanc", "4y3khh7ld4")
print 'C:\\Users\\sleblan2\\Research\\python_codes\\file.rc'
get_ipython().magic(u'matplotlib inline')
from matplotlib import rc_file
rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
# set the basic directory path
fp='C:\\Users\\sleblan2\\Research\\TCAP\\'
if __name__ == "__main__":
    print('yes')


# In[3]:


# load the idl save file containing the modeled radiances
s=sio.idl.readsav(fp+'model/sp_v1_20130219_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp, wvl, z, re, ta)'


# In[4]:


# create custom key for sorting via wavelength
iwvls = np.argsort(s.zenlambda)
s.zenlambda.sort()


# In[5]:


# load the matlab file containing the measured TCAP radiances
m = sio.loadmat(fp+'4STAR/20130219starzen_rad.mat')
sm = sio.idl.AttrDict(m)
print sm.keys()
print 'Measured radiance Shape: ', sm.rad.shape

print np.nanmax(sm.rad[sm.good[100],:])
sm.good[100]


# In[49]:


import Sp_parameters as Sp
reload(Sp)


# In[50]:


lut = Sp.Sp(s)
warnings.simplefilter('ignore')
lut.sp_hires()


# In[51]:


lut.params()


# In[52]:


pcoef = lut.norm_par()


# In[55]:


print pcoef['coef'].shape


# In[ ]:




