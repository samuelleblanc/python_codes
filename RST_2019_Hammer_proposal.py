
# coding: utf-8

# # Info
# Name:  
# 
#     RST_2019_Hammer_proposal
# 
# Purpose:  
# 
#     Make a map plot of the location of a skyscan by 4STAR.
#   
# Input:
# 
#     none
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
#     - matplotlib
#     - scipy
# 
#   
# Needed Files:
# 
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc, Over North Dakota on the way to Paris, 2019-03-30
#     

# # Load python modules

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
import hdf5storage as hs
from plotting_utils import prelim
from datetime import datetime
from scipy.interpolate import UnivariateSpline
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


fp =getpath('ORACLES')#'C:/Userds/sleblan2/Research/ORACLES/'
fp


# # Make the map

# In[31]:


np.linspace(-15,-14,6)


# In[43]:


def mapfig(ax=plt.gca()):
    m = Basemap(projection='merc',llcrnrlat=-8.2,urcrnrlat=-7.7,llcrnrlon=-14.6,urcrnrlon=-14,resolution='h',ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.linspace(-15,-14,6),labels=[0,0,0,1],linewidth=0.1)
    m.drawparallels(np.linspace(-8.2,-7.5,8),labels=[1,0,0,0],linewidth=0.1)
    m.shadedrelief(alpha=0.7)
    return m


# In[48]:


fp


# In[49]:


fig,ax = plt.subplots(1,1,figsize=(5,5))
m = mapfig(ax=ax)
m.fillcontinents(color='g')
mx,my = m(-7.967,-14.350)
p = m.plot(-14.350,-7.997,'*y',latlon=True,markersize=24)
mx,my = m(-14.380,-7.957)
ax.text(mx,my,'Ascension Island')
#p = ax1.pcolor(mx,my,ar6['bin_mean'],vmin=0.0,vmax=0.8,cmap='plasma')
plt.savefig(fp+'Ascension_Island_location.png',dpi=600,transparent=True)


# In[40]:





# In[19]:


help(m.plot)


# In[35]:


help(m.fillcontinents)

