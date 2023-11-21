#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     To load a few different maps of the models used for ARCSIX MSS
#     Ensure that pre and post regridding show the same data
# 
# Input:
# 
#     arguments
# 
# Output:
# 
#     Figures
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
# 
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
#     
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2023-11-17
#     Modified:
# 

# # Prepare python environment

# In[1]:


import numpy as np
import netCDF4 as nc
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.colors as colors
get_ipython().run_line_magic('matplotlib', 'notebook')
import os


# In[2]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature


# In[3]:


vv = 'v1'
fp = '/data2/mss/rasmdata/'
name = 'ARCSIX_MSS'


# # Load files

# ## Load OG ARCSIX - RASM-ESRL

# In[4]:


os.listdir(fp+'old/')


# In[5]:


f_og = fp+'old/ARCSIX_2023-11-15.nc'


# In[6]:


og,og_dict = lu.load_netcdf(f_og,everything=True)


# ## Load regridded ARCSIX  - RASM-ESRL

# In[7]:


os.listdir(fp)


# In[8]:


f_re = fp+'RASM_ESRL_ml.fromARCSIX_2023-11-15.nc'


# In[9]:


re,re_dict = lu.load_netcdf(f_re,everything=True)


# # Plot the different maps

# ## Load map defaults

# In[9]:


def make_map(lower_left=[-100,72],upper_right=[20,73]):
    proj = ccrs.NorthPolarStereo()
    m = plt.figure()
    ax = m.add_subplot(111, projection=proj)
    ax.plot([0],[0])
    ax.coastlines()
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',edgecolor='k',
                                            facecolor=cfeature.COLORS['land']+0.0625)
    provinces_50m = cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces_lines','50m',facecolor='none')
    ax.gridlines(draw_labels=True,auto_update=True)
    ax.add_feature(land_50m,alpha=0.1)
    ax.add_feature(cfeature.LAKES, facecolor=[0.69375   , 0.81484375, 0.9828125 ],alpha=0.3)
    ax.add_feature(cfeature.RIVERS,alpha=0.2)
    ax.add_feature(provinces_50m)
    ax.add_feature(cfeature.BORDERS)

    ax.set_extent([lower_left[0], upper_right[0],lower_left[1],upper_right[1]])
    return m,ax


# In[10]:


fig,ax = make_map()
fig.show()


# ## plots maps

# In[11]:


og.keys()


# In[12]:


tr = ccrs.PlateCarree()


# In[37]:


fig,ax = make_map()
cf = ax.pcolormesh(og[b'alon'],og[b'alat'],og[b'lwp'][0,:,:],cmap=plt.cm.summer_r,transform=tr,vmin=0,vmax=500)
plt.title('RASM_ESRL - Cloud liquid water path - 2023-11-15 - original')
plt.colorbar(cf)


# In[14]:


ilat1, ilon1 = 78.0,13.0
ilat2, ilon2 = 80.0,30.0


# In[15]:


iarg1 = np.argmin((og[b'alon']-ilon1)**2+(og[b'alat']-ilat1)**2)
iarg2 = np.argmin((og[b'alon']-ilon2)**2+(og[b'alat']-ilat2)**2)
og[b'lwp'][0,:,:].flatten()[iarg1],og[b'lwp'][0,:,:].flatten()[iarg2]


# In[16]:


iarg1,iarg2


# In[17]:


og[b'lwp'][0,:,:].flatten()[iarg1],og[b'lwp'][0,:,:].flatten()[iarg2]


# In[20]:


plt.figure()
plt.hist(og[b'lwp'][0,:,:].flatten())


# In[21]:


fig,ax = make_map()
cf = ax.pcolormesh(re[b'lon_grid'],re[b'lat_grid'],re[b'lwp'][0,:,:],cmap=plt.cm.summer_r,transform=tr,vmin=0,vmax=500)
plt.colorbar(cf)


# In[24]:


xx,yy = np.meshgrid(re[b'lon_grid'],re[b'lat_grid'])


# In[25]:


iarg1 = np.argmin((xx-ilon1)**2+(yy-ilat1)**2)
iarg2 = np.argmin((xx-ilon2)**2+(yy-ilat2)**2)
re[b'lwp'][0,:,:].flatten()[iarg1],re[b'lwp'][0,:,:].flatten()[iarg2]


# # Make the lat lon regrid by hand to check status

# In[18]:


import scipy.interpolate as interp


# In[19]:


def regrid_latlon(data_file,lon_name='lon',lat_name='lat',named_variables=['sic_merged'],debug=False):
    'Regrid the latitude and longitude 2d matrices into a 1d lat/lon in regular grids'
    import netCDF4 as nc
    import numpy as np
    import scipy.interpolate as interp
    import datetime
    
    vars = {}
    vars_dict = {}
    
    with nc.Dataset(data_file,'a',clobber=False) as file:
        lon_var = file[lon_name][:]
        lat_var = file[lat_name][:]
        time_var = file['time'][:]
        for n in named_variables:
            vars[n] = file[n][:]
            vars_dict[n] = file.variables[n]
        
        new_lon = np.linspace(start=np.nanmin(lon_var),stop=np.nanmax(lon_var),num=lon_var.shape[1])
        new_lat = np.linspace(start=np.nanmin(lat_var),stop=np.nanmax(lat_var),num=lat_var.shape[0])
        xm,ym = np.meshgrid(new_lon,new_lat)
        
        #build the new vars space
        new_vars = {}
        for n in named_variables:
            if len(vars[n].shape) > 2:
                print(f"...regridding file '{data_file}' variable: '{n}'")
                new_vars[n] = np.zeros_like(vars[n])
                vars[n].mask = ~np.isfinite(vars[n])
                for it,t in enumerate(time_var):
                    xmask,ymask = np.where(~vars[n][it,:,:].mask)
                    try:
                        new_vars[n][it,:,:] = interp.griddata((lat_var[xmask,ymask].flatten(),lon_var[xmask,ymask].flatten()),vars[n][it,xmask,ymask].flatten(),(xm,ym),'nearest')
                    except Exception as e:
                        if debug: import ipdb; ipdb.set_trace()
                        print(f"Error gridding variable '{n}': {str(e)}")


# In[20]:


lon_var = og[b'alon'][:]
lat_var = og[b'alat'][:]


# In[21]:


new_lon = np.linspace(start=np.nanmin(lon_var),stop=np.nanmax(lon_var),num=lon_var.shape[1])
new_lat = np.linspace(start=np.nanmin(lat_var),stop=np.nanmax(lat_var),num=lat_var.shape[0])
xm,ym = np.meshgrid(new_lon,new_lat)


# In[22]:


vars = {}
new_vars = {}


# In[23]:


vars['lwp'] = og[b'lwp'][:]


# In[32]:


plt.figure()
plt.hist(vars[n][0,:,:].flatten())


# In[24]:


n = 'lwp'


# In[25]:


new_vars[n] = np.zeros_like(vars[n])


# In[26]:


vars[n].mask = ~np.isfinite(vars[n])


# In[27]:


it = 0


# In[28]:


xmask,ymask = np.where(~vars[n][it,:,:].mask)


# In[83]:


new_vars[n][it,:,:] = interp.griddata((lat_var[xmask,ymask].flatten(),\
                                       lon_var[xmask,ymask].flatten()),\
                                      vars[n][it,xmask,ymask].flatten(),(xm,ym), method='cubic')


# In[84]:


plt.figure()
plt.hist(new_vars[n][0,:,:].flatten())


# In[29]:


new_lon.shape,new_lat.shape


# In[35]:


iarg1 = np.argmin((xm-ilon1)**2+(ym-ilat1)**2)
iarg2 = np.argmin((xm-ilon2)**2+(ym-ilat2)**2)
new_vars[n][0,:,:].flatten()[iarg1],new_vars[n][0,:,:].flatten()[iarg2]


# In[36]:


new_vars[n][0,:,:]


# ## Run new regridding

# In[30]:


from scipy.interpolate import RectSphereBivariateSpline


# In[31]:


def upscale_field(lons, lats, field, x_scale=2, y_scale=2, is_degrees=True):
    '''
    Takes a field defined on a sphere using lons/lats and returns an upscaled
    version, using cubic spline interpolation.
    '''
    if is_degrees:
        lons = lons * np.pi / 180.
        lats = (90.0-lats) * np.pi / 180.

    d_lon = lons[1] - lons[0]
    d_lat = lats[1] - lats[0]

    new_lon = np.linspace(lons[0], lons[-1], len(lons) * x_scale)
    new_lat = np.linspace(lats[0], lats[-1], len(lats) * x_scale)

    mesh_new_lat, mesh_new_lon = np.meshgrid(new_lat, new_lon)

    if True:
        lut = RectSphereBivariateSpline(lats[1:-1], lons[1:-1], field[1:-1,
                                                                      1:-1])

        interp_field = lut.ev(mesh_new_lat[1:-1, 1:-1].ravel(),
                              mesh_new_lon[1:-1, 1:-1].ravel()).reshape(
                                  mesh_new_lon.shape[0] - 2,
                                  mesh_new_lon.shape[1] - 2).T
    else:
        pass
    
    if is_degrees:
        new_lat = 90-(new_lat*180/np.pi)
        new_lon = new_lon*180/np.pi
    
    return new_lon[1:-1], new_lat[1:-1], interp_field


# In[35]:


(90-new_lat)/180*np.pi


# In[40]:


new_lons2,new_lats2,lwp2 = upscale_field(new_lon,new_lat[::-1],vars[n][it,:,::-1],x_scale=1,y_scale=1)


# In[41]:


xm,ym = np.meshgrid(new_lons2,new_lats2)
iarg1 = np.argmin((xm-ilon1)**2+(ym-ilat1)**2)
iarg2 = np.argmin((xm-ilon2)**2+(ym-ilat2)**2)
lwp2.flatten()[iarg1],lwp2.flatten()[iarg2]


# In[42]:


new_lon.shape,lwp2.shape


# In[44]:


fig,ax = make_map()
cf = ax.pcolormesh(new_lons2,new_lats2,lwp2[::-1,:],cmap=plt.cm.summer_r,transform=tr,vmin=-10,vmax=500)
plt.colorbar(cf)


# In[ ]:




