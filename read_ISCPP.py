
# coding: utf-8

# Name:  
# 
#     read_ISCPP
# 
# Purpose:  
# 
#     Script to read the iscpp cloud cover data in IEEE floating point 
#     number format
#     taken from the website: http://isccp.giss.nasa.gov/products/browsed2.html
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures...
#   
# Dependencies:
# 
#     - struct
#     - matplotlib
#     - basemap
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - iscpp SQ downloaded file: renamed to ISCPP_D2_CLDCOVER_sq_20151207.ieee
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, NASA Ames, Santa Cruz, CA, 2015-12-07

# In[ ]:

print 'yes'


# In[ ]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt


# In[ ]:


import struct
from mpl_toolkits.basemap import Basemap
get_ipython().magic(u'matplotlib notebook')
import numpy as np


# In[ ]:

fp='C:/Users/sleblan2/Research/sat/iscpp/'


# ## Load the ieee file

# In[ ]:

ff = fp+'ISCPP_D2_CLDCOVER_sq_20151207.ieee'
f = open(ff,'rb')


# In[ ]:

n = 10368
num = struct.unpack('{}f'.format(n),f.read(4*n))


# In[ ]:

f.close()


# In[ ]:

cld = np.array(num)


# In[ ]:

cld.shape


# In[ ]:

lat = np.arange(72)*2.5-90.0+1.25
lon = np.arange(144)*2.5-180.0+1.25


# In[ ]:

lat


# In[ ]:

cv = cld.reshape(72,144)


# In[ ]:

fig,ax = plt.subplots(1)
ax.pcolorfast(lon,lat,cv,cmap=plt.cm.rainbow_r)


# In[ ]:

print 'Average cloud cover at :{}'.format(cv.mean())


# In[ ]:

m = Basemap(projection='moll',lon_0=0,resolution='c')
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,20.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,420.,60.))
mlon, mlat = np.meshgrid(lon,lat,sparse=False)
x,y=m(mlon, mlat) 
cs = m.contourf(x,y,cv,np.linspace(0,100,21))
cbar = plt.colorbar(cs)
cbar.set_label('Cloud Cover [\%]')
plt.savefig(fp+'ISCPP_CLDCOVER_mean_annual.png',dpi=600,transparent=True)


# In[ ]:

import cmaps
cmaps.cmaps()


# In[ ]:



