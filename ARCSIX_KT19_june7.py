#!/usr/bin/env python
# coding: utf-8

# # Intro
# 
# Name:  
# 
#     ARCSIX_KT19_june7
# 
# Purpose:  
# 
#     Look at the KT19 data in the ARCSIX, in the june 7 case where there was ice in cloud 
# 
# Input:
# 
#     none at command line
# 
# Output:
# 
#     plots
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - numpy
#     - pandas
# 
# Needed Files:
# 
#   - kt19 ict files
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, 2026-04-28
#     Modified: 

# # Load the required modules

# In[1]:


get_ipython().run_line_magic('matplotlib', 'widget')


# In[3]:


import matplotlib 
import matplotlib.pyplot as plt
import numpy as np
import os
import netCDF4 as nc
from datetime import datetime
import matplotlib.patches as patches
import matplotlib.patheffects as pe
import json
import pandas as pd
import cartopy.crs as ccrs
from path_utils import getpath


# In[4]:


import path_utils as pu
import load_utils as lu


# In[5]:


import plotly.express as px
import plotly.io as pio
pio.renderers.default = "notebook"


# In[6]:


fp = getpath('ARCSIX')#,path=r"/data2/ARCSIX/data",make_path=True)
fp


# # Load the files

# In[11]:


fo = fp+'ARCSIX-MetNav-KT19-10Hz_P3B_20240607_R0.ict'
kt19 = lu.load_ict(fo)


# # Plot some of the data

# In[12]:


fig = plt.figure()
plt.plot(kt19['Time_Start']/3600.0,kt19['IR_Surf_Temp'],'.')
plt.xlim(15.34,15.75)
plt.ylim(-4,0)
plt.show()


# In[13]:


import time


# In[14]:


plt.close('all')
fig,ax = plt.subplots()
plt.plot(kt19['Time_Start']/3600.0,kt19['IR_Surf_Temp'],'.')
plt.xlim(15.34,15.75)
plt.ylim(-4,0)
fig.canvas.draw()
time.sleep(0.5)
display(fig.canvas)


# In[10]:


fig = px.line(x=kt19['Time_Start']/3600.0,y=kt19['IR_Surf_Temp'])
fig.show()


# # Load from earthdata

# In[1]:


import earthaccess


# In[19]:


earthaccess.login()


# In[4]:


results = earthaccess.search_datasets(
    keyword="MERRA-2"
    )


# In[5]:


len(results)


# In[19]:


results[0]['umm']['CollectionCitations'][0]['Title']


# In[25]:


results[0]['meta']['native-id']


# In[27]:


results[0]['umm'].keys()


# In[37]:


results[4]['umm']['ShortName']


# In[38]:


results[4]['meta']['concept-id']


# In[28]:


results[0]['meta'].keys()


# In[58]:


[(r['meta']['concept-id'],r['umm']['CollectionCitations'][0].get('Title','NA')) for r in results if 'CollectionCitations' in r['umm']]


# In[41]:


res = earthaccess.search_data(concept_id='C1276812879-GES_DISC',temporal=("2024-05-17","2024-06-17"))


# In[46]:


res = earthaccess.search_data(concept_id='C1276812879-GES_DISC',temporal=("2024-06-07","2024-06-08")) # 3d


# In[47]:


res


# In[52]:


fileobjects = earthaccess.download([res[0]],local_path='/data2/ARCSIX/data/')


# In[53]:


fileobjects


# In[59]:


res2d = earthaccess.search_data(concept_id='C1276812863-GES_DISC',temporal=("2024-06-07","2024-06-08"))


# In[60]:


res2d


# In[61]:


fileobjects2d = earthaccess.download([res2d[0]],local_path='/data2/ARCSIX/data/')


# ## load the MERRA-2 file

# In[55]:


import xarray


# In[56]:


ds = xarray.open_dataset(fileobjects[0])


# In[62]:


ds2d = xarray.open_dataset(fileobjects2d[0])


# In[65]:


ds2d['TROPQ']

