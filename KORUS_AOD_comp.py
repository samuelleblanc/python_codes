
# coding: utf-8

# # Info
# Name:  
# 
#     KORUS_AOD_comp
# 
# Purpose:  
# 
#     Comparison of AOD from 4STAR along flight track and GOCI aerosol
#     Additional calculations of the aerosol extinction profile
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
#     - pytables
#     - os
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, OSAN AFB, Korea, 2016-05-06
#     Modified: 

# # Import the required modules and set up base

# In[164]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.io as sio
import Sp_parameters as Sp
import tables
import load_utils as lm
import hdf5storage as hs
import os


# In[165]:

from mpl_toolkits.basemap import Basemap,cm
get_ipython().magic(u'matplotlib notebook')
fp = 'C:/Users/sleblan2/Research/KORUS-AQ/'


# # Load the various data

# ## Load the 4STAR starsun

# In[33]:

f_star = fp+'data\\20160504starsun.mat'


# In[34]:

s = sio.loadmat(f_star)


# In[35]:

s.keys()


# In[45]:

s['utc'] = lm.toutc(lm.mat2py_time(s['t']))


# In[31]:

s['tau_aero'].shape


# ## Load the GOCI aerosol products

# In[3]:

f_goci = fp+'sat/GOCI/GOCI_YAER_AOP_20160505041644.hdf'


# In[15]:

gg,gg_head = lm.load_hdf(f_goci,values=(('lon',0),('lat',1),('aod550',2),('fmf550',3),('ssa440',4),('type',5),('ang',6),('QA',7),
                                ('obs_time',8),('cf',9),('turbidI',10),('Land_sea_mask',11)))


# In[57]:

gg_head


# ## Get the AERONET data to overlay on plot

# In[163]:

fp


# In[191]:

reload(lm)


# In[183]:

fa = fp+'aeronet/AOT/LEV10/ALL_POINTS/'
fa_l = os.listdir(fa)


# In[192]:

aero = []
for f in fa_l:
    aero.append(lm.load_aeronet(fa+f))


# In[194]:

aero[0].keys()


# # Start making different plots/maps

# In[17]:

#set up a easy plotting function
def make_map(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=128,lat_0=36.0,
            llcrnrlon=123.0, llcrnrlat=33.5,
            urcrnrlon=132.0, urcrnrlat=39,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(123,133,11),labels=[0,0,0,1])
    m.drawparallels(np.linspace(33,39,13),labels=[1,0,0,0])
    return m


# ## Start with simple map plot of GOCI

# In[41]:

fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg['lon'],gg['lat'])
clevels = np.linspace(0,4,41)

plt.title('GOCI AOD 2016-05-05 04:16:44')
cs1 = m.contourf(x,y,gg['aod550'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('AOD 550 nm')

#xx,yy = m(star['lon'],star['lat'])
#m.scatter(xx,yy,c=star['tau'],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
#          alpha=0.5,edgecolors='k',linewidth=0.65)
plt.savefig(fp+'plot/20160505_GOCI_map_AOD.png',dpi=600,transparent=True)


# In[56]:

fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg['lon'],gg['lat'])
clevels = np.linspace(0,40,41)

plt.title('GOCI AOD 2016-05-05 04:16:44')
cs1 = m.contourf(x,y,gg['obs_time'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('Observation Time')


# In[59]:

fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg['lon'],gg['lat'])
clevels = np.linspace(0,3,4)

plt.title('GOCI AOD 2016-05-05 04:16:44')
cs1 = m.contourf(x,y,gg['QA'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('Quality flag')


# ## Overlay 4STAR values

# In[55]:

plt.figure()
plt.plot(s['utc'],s['tau_aero'][:,450])
plt.ylabel('4STAR AOD at {} nm'.format(s['w'][0][469]*1000.0))
plt.xlabel('UTC [H]')


# In[61]:

ig = gg['QA']==3


# In[62]:

ig.shape


# In[63]:

gg['aod550'].shape


# In[64]:

gg['aod550'][ig]=np.nan


# In[67]:

fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg['lon'],gg['lat'])
clevels = np.linspace(0,4,41)

plt.title('GOCI AOD 2016-05-05 04:16:44')
cs1 = m.contourf(x,y,gg['aod550'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('AOD 550 nm')
m.scatter(x,y,c=gg['aod550'],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],edgecolors='None')


xx,yy = m(s['Lon'],s['Lat'])
m.scatter(xx,yy,c=s['tau_aero'][:,469],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
          alpha=0.5,edgecolors='None')
plt.savefig(fp+'plot/20160505_GOCI_4STAR_map_AOD.png',dpi=600,transparent=True)


# #

# In[69]:

fp


# In[ ]:




# In[70]:

fa = fp+'aeronet/AOT/LEV10/ALL_POINTS/160401_160731_SONET_Shanghai.lev10'


# In[ ]:



