#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Describe the details ...
# 
# Input:
# 
#     arguments
# 
# Output:
# 
#     Figure and save files
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - load_utils.py
#     - matplotlib
#     - numpy
#     - Sp_parameters
#     - write_utils
#     - path_utils
#     - hdf5storage
#     - scipy
# 
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-10-15
#     Modified:
# 

# # Prepare python environment

# In[1]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import os


# In[2]:


name = 'ORACLES'
vv = 'v1'
fp = getpath(name)


# In[3]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature


# In[4]:


import cartopy


# In[72]:


import importlib
importlib.reload(ccrs)


# # Test out the cartopy plotting

# In[5]:


mrc = ccrs.Mercator()


# In[6]:


lower_left = [-92,30]
upper_right = [-65,50]


# In[25]:


proj = ccrs.PlateCarree()
ax = plt.axes(projection=proj)
ax.coastlines()

land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',edgecolor='k',
                                        facecolor=cfeature.COLORS['land']+0.0625)
provinces_50m = cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces_lines','50m',facecolor='none')
gl = ax.gridlines(draw_labels=True,auto_update=True)
#mland = ax.add_feature(land_50m,zorder=-100)
ax.add_feature(cfeature.LAKES, facecolor=[0.69375   , 0.81484375, 0.9828125 ],alpha=0.3)
ax.add_feature(cfeature.RIVERS, facecolor=[0.69375   , 0.81484375, 0.9828125 ],alpha=0.2)
ax.add_feature(provinces_50m)
ax.add_feature(cfeature.BORDERS)

ax.set_extent([lower_left[0], upper_right[0],lower_left[1],upper_right[1]])


# In[89]:


ax.figure.canvas.copy_from_bbox(ax.axes.bbox)


# In[50]:


lo = ax.plot([-85.0,-80,-75],[36,39,42],label='ta')


# In[51]:


leg = ax.legend()


# In[52]:


li = leg.get_lines()


# In[53]:


li[0]


# In[55]:


lo[0].set_visible(False)
ax.get_lines()


# In[56]:


lu = ax.get_lines()


# In[70]:


isVisible


# In[69]:


if li[0].get_alpha() is not None and li[0].get_alpha()>0.3: 
    isVisible=True
else: isVisible=False


# In[86]:


True if lo[0].get_alpha() is None else (True if lo[0].get_alpha() else False)


# In[82]:


True if li[0].get_alpha() else False


# In[85]:


lo[0].set_alpha(0.0)


# # Plot out data

# In[15]:


ax.axes.lines


# In[ ]:





# # Test out url reading

# In[34]:


from urllib.request import urlopen
from bs4 import BeautifulSoup


# In[27]:


url = 'http://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v2_globe?year=2021&month=10&day=20&year2=2021&month2=10&day2=20&LEV10=1&AVG=20&lat1=30.000000&lat2=50.519141&lon1=-92.000000&lon2=-65.000000'


# In[28]:


htm = urlopen(url)


# In[29]:


htm


# In[30]:


html = htm.read()


# In[35]:


soup = BeautifulSoup(html)


# In[36]:


soup


# In[ ]:




