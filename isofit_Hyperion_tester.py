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

# In[7]:


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


# In[8]:


from spectral.io import envi
from isofit.core.common import envi_header


# ## set the paths

# In[9]:


from path_utils import getpath


# In[10]:


figdir = getpath('hyperion_fig')
figprefix = os.path.join(figdir, '{}')
figexts = ('png', 'eps', )


# In[11]:


from isofit.utils import apply_oe


# In[12]:


apply_oe


# In[13]:


import importlib
importlib.reload(hyperion2isofit)


# In[14]:


scene = r'EO1H0440342016184110KF'#'EO1H0100612008087110KF'#'EO1H0090582008189110P0' # 'EO1H0440342016184110KF' # 'EO1H0010742008014110K0' # 'EO1H0440342016184110KF'
scene = r'EO1H0400332013093110K7'
paths = hyperion2isofit.readHyperionL2(scene)


# # Test out the make hyp L2

# In[32]:


import importlib; importlib.reload(hyperion2isofit)


# In[41]:


import importlib; importlib.reload(apply_oe)


# In[42]:


help(apply_oe)


# In[45]:


apply_oe.ray.shutdown()


# In[ ]:


exit_code = hyperion2isofit.MakeHypL2(scene)


# In[47]:


exit_code


# In[ ]:


paths = {'L1R':'/data2/SBG/Hyperion_data/EO1H0440342016184110KF.L1R','radiance':'/data2/SBG/Hyperion_data/EO1H0440342016184110KFradiance.hdr'}


# In[130]:


radiance_dataset = envi.open(envi_header(paths['radiance'][:-4]))


# In[131]:


radiance_dataset.metadata.keys()


# # Plot out the output data from this run

# In[48]:


paths


# ## Use plotting from test modules

# In[49]:


scenes = [scene]


# In[50]:


titles = [datetime.strftime(datetime.strptime(scene[10:17], '%Y%j'), '%B %d, %Y')+' '+scene[:22] for scene in scenes]


# In[51]:


key, ylabel, figsuffix = ('rfl', 'Surface Reflectance', 'Reflectance')
sws=(('VNIR', (27, 20, 9)), ('SWIR', (90, 130, 170)))


# In[52]:


arrss = []; ilocss = []; titles = []
for scene in scenes:
    paths = hyperion2isofit.readHyperionL2(scene)
    ilocss.append(test.indicesofinterest(scene))
    arrss.append([envi.open(paths[key])[:,:,:],])
    td = paths['attrs']['Target_Description'] if 'Target_Description' in paths['attrs'] else '' 
    titles.append(td[:td.rfind('[')-1])
    #if key=='rfl_10nm_cover_class':
    #    classnames = ['algae','coral','mud/sand','seagrass']
    #    legendstrs = ['{}, n = {}'.format(classnames[u], np.sum(arrss[0][0]==u)) for u in range(4)]


# In[53]:


if key in ('rfl', 'radiance'):
    with test.makecollage(arrss, paths['wp'], test.sws, ilocss, figprefix.format(scenes[0]), figsuffix, ylabel, ylims=[], labels_maps=titles, labels_spectra=[], autofit=True) as mc:
        # for ax in mc.axss:
        #     ax.set_ylim(0,1.2)
        plt.gcf().patch.set_facecolor('white')


# ## Load some images 

# In[54]:


paths['rfl']


# In[55]:


p = envi.open(paths['rfl'])


# In[56]:


lwp = len(paths['wp'])
idx = np.append(np.isfinite(paths['wp']),[True]*(p.shape[-1]-lwp)).astype(bool)
idxf = np.append(np.isfinite(paths['wp']),[False]*(p.shape[-1]-lwp)).astype(bool)
idx = (slice(None), slice(None), idx)
idxf = (slice(None), slice(None), idxf)


# In[57]:


p.shape


# In[58]:


test.sws


# In[59]:


test.sws[0][1]


# In[60]:


get_ipython().run_line_magic('matplotlib', 'widget')


# In[61]:


np.nanmin(p[:,:,test.sws[0][1]]),np.nanmean(p[:,:,test.sws[0][1]]),np.nanmedian(p[:,:,test.sws[0][1]]),np.nanmax(p[:,:,test.sws[0][1]]) 


# In[62]:


import matplotlib


# In[63]:


matplotlib.scale.get_scale_names()


# In[64]:


fig,ax = plt.subplots(1)
ax.imshow(p[:,:,test.sws[0][1][0]],norm='log') #,vmin=0.000000000000000000001,vmax=0.0001,norm='log')
#plt.colorbar()


# In[65]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature


# In[66]:


test.sws[0][1]


# In[67]:


p_na = np.array(p[:,:,test.sws[0][1]])
p_na[p_na<0] = np.nan


# In[68]:


fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.pcolormesh(paths['lonlat'][0,:,:],paths['lonlat'][1,:,:],p_na,vmin=0,vmax=0.000000000001)
ax.add_feature(cfeature.COASTLINE, edgecolor='grey')
gl = ax.gridlines(draw_labels=True)


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




