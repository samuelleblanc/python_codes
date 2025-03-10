#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     To test out the isofit v3 Hyperion reprocessing
# 
# 
# Output:
# 
#     Figure
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - hyperion package built by Yohei SHinozuka
#     - hyperion L1R files
# 
# Needed Files:
#   - hyperion L1R hdf files
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2025-02-12
#             - based on test.ipynb file from Yohei Shinozuka
#     Modified:
# 

# # Prepare python environment

# In[1]:


import sys
#sys.path = sys.path[1:]
import os
from glob import glob
from importlib import reload
import platform
from datetime import datetime
import filecmp
from IPython.display import HTML, display
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.ticker as mticker
plt.figure() # Run this before importing isofit.utils (via h2i)
from netCDF4 import Dataset
from spectral.io import envi
from scipy.io import loadmat
from hyperion import test, testmemory, hyperion, hyperion2isofit


# ## set the paths

# In[2]:


from path_utils import getpath


# In[4]:


figdir = getpath('hyperion_fig')
figprefix = os.path.join(figdir, '{}')
figexts = ('png', 'eps', )


# In[5]:


from isofit.utils import apply_oe


# In[6]:


apply_oe


# In[63]:


import importlib
importlib.reload(hyperion2isofit)


# In[ ]:


scene = r'EO1H0440342016184110KF'#'EO1H0100612008087110KF'#'EO1H0090582008189110P0' # 'EO1H0440342016184110KF' # 'EO1H0010742008014110K0' # 'EO1H0440342016184110KF'
paths = hyperion2isofit.readHyperionL2(scene)


# # Test out the make hyp L2

# In[135]:


import importlib; importlib.reload(hyperion2isofit)


# In[ ]:


exit_code = hyperion2isofit.MakeHypL2(scene)


# In[ ]:


exit_code


# In[122]:


paths = {'L1R':'/data2/SBG/Hyperion_data/EO1H0440342016184110KF.L1R','radiance':'/data2/SBG/Hyperion_data/EO1H0440342016184110KFradiance.hdr'}


# In[125]:


from spectral.io import envi
from isofit.core.common import envi_header


# In[130]:


radiance_dataset = envi.open(envi_header(paths['radiance'][:-4]))


# In[131]:


radiance_dataset.metadata.keys()


# In[70]:


paths['L1R'].replace('L1R','AUX')


# In[73]:


file = paths['L1R']
dataset = netCDF4.Dataset(file)


# In[80]:


if 'EO1H0440342016184110KF.L1R' in dataset.variables: print('a')


# In[83]:


os.path.split(paths['L1R'])


# In[78]:


dataset['EO1H0440342016184110KF.L1R']


# # Old comparison 

# In[ ]:


if key=='rfl':
    ok = 90 # band
    sa = np.flip((p-m)[..., ok])
    x = np.arange(sa.shape[1])
    xscale = 'linear'
    xlabel = 'Sample #'
elif key in ('subs_rfl', 'subs_h2o'):
    x = paths['wp'][idx[-1][:lwp]]
    xscale = 'log'
    xlabel = 'Wavelength (um)'
    sa = np.flip((p-m).squeeze()[:, idxf[-1]])
clim = np.percentile(sa, [32.1/2,100-32.1/2])
plt.figure()
plt.pcolormesh(x, np.flip(np.arange(sa.shape[0])), sa, vmin=clim[0], vmax=clim[-1])
plt.xscale(xscale)
plt.xlabel(xlabel)
plt.ylabel('Line #')
plt.colorbar()
plt.show()


# In[66]:


# a file in my directory

key = 'subs_rfl'
mine = paths[key] 
# the counterpart in PT's directory
if scene=='EO1H0010742008014110K0':
    PT = os.path.join('/nobackupp10/hyperion/pipeline/results/task-data/27-45-level1/st-50', mine[mine.find(scene)+23:])
# elif scene=='EO1H0090582008189110P0':  
#     PT = os.path.join('/nobackupp10/hyperion/pipeline/results/task-data/29-47-level1/st-168', mine[mine.find(scene)+23:])
# /nobackupp10/hyperion/pipeline/results/task-data/29-47-level1/st-199/EO1H0100612008087110KF.L2.tgz
    replace = ('', '')
elif scene=='EO1H0440342016184110KF':
    replace = ('isofit/{}'.format(scene), 'isofit/{}.L2_20211130'.format(scene))
else:
    replace = ('isofit/{}'.format(scene), 'isofit/{}.L2'.format(scene))
PT = mine.replace(replace[0],replace[1]) 
# mine = mine.replace('isofit/{}'.format(scene),'isofit/{}.L2_20211116'.format(scene)) # !!!!!!!!
print(mine+'\n'+PT)
if not os.path.exists(mine):
    raise ValueError(mine)
if not os.path.exists(PT):
    raise ValueError(PT)
    
# open the files
p = envi.open(PT)[:,:,:]
m = envi.open(mine)[:,:,:]
lwp = len(paths['wp'])
idx = np.append(np.isfinite(paths['wp']),[True]*(p.shape[-1]-lwp)).astype(bool)
idxf = np.append(np.isfinite(paths['wp']),[False]*(p.shape[-1]-lwp)).astype(bool)
idx = (slice(None), slice(None), idx)
idxf = (slice(None), slice(None), idxf)
# are they identical?
if np.all(p==m):
    print('Identical.')
elif np.all(p[idx]==m[idx]):
    print('Identical outside absorption bands.')
else:
    d = np.abs(p-m)
    imax = np.nanargmax(d)
    print('{:0.2f} max diff, p={:0.2g}, m={:0.2g}'.format(d.flatten()[imax], p.flatten()[imax], m.flatten()[imax]))
    r = d/m
    imax = np.nanargmax(r)
    print('{:0.2g}% max diff, p={:0.2g}, m={:0.2g}'.format(np.nanmax(d), r.flatten()[imax]*100, p.flatten()[imax], m.flatten()[imax]))
    d = np.abs(p[idx]-m[idx])
    imax = np.nanargmax(d)
    print('{:0.2f} max diff, p={:0.2g}, m={:0.2g}'.format(d.flatten()[imax], p[idx].flatten()[imax], m[idx].flatten()[imax]))
    r = d/m[idx]
    imax = np.nanargmax(r)
    print('{:0.2g}% max diff, p={:0.2g}, m={:0.2g}'.format(np.nanmax(d), r.flatten()[imax]*100, p[idx].flatten()[imax], m[idx].flatten()[imax]))


# In[ ]:




