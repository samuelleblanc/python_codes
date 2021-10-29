#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Read in L1 hyperion files. Geolocate with .csv file, output new file
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-08-03
#     Modified:
# 

# # Parse command line

# In[1]:


import argparse


# In[2]:


long_description = """    Pull in Hyperion L1R file and hyperion_metadata.csv to calculate the geolocation
    The output is a new netcdf file.
        if selected creates new file with just the geolocation data, if not saves again the radiance data"""


# In[8]:


parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-f','--file_name',nargs='?',
                    help='Input filename',
                    default='EO1H0080122015195110K4.L1R')
parser.add_argument('-r','--root_dir',nargs='?',
                    help='full file path of the root directory to read from',
                    default='/data/sam/SBG/data/')
parser.add_argument('-m','--hyperionmeta_dir',nargs='?',
                    help='full file path of the directory which has the hyperion metadata file',
                    default='/data/sam/SBG/data/')
parser.add_argument('-o','--out_dir',nargs='?',
                    help='full file path of the output directory',
                    default='/data/sam/SBG/data/')
parser.add_argument('-g','--only_geo',help='if set, will only save the geolocation information to new file',
                    action='store_true')
parser.add_argument('-q','--quiet',help='if set, quiet the comments',
                    action='store_true')


# In[9]:


in_ = vars(parser.parse_known_args()[0])


# In[10]:


fp = in_.get('root_dir','/data/sam/SBG/data/')
fph = in_.get('hyperionmeta_dir','/data/sam/SBG/data/')
fp_out = in_.get('out_dir','/data/sam/SBG/data/')
only_geo = in_.get('only_geo',False)
fname = in_.get('file_name','EO1H0080122015195110K4.L1R')
verbose = not in_.get('quiet',False)


# # Prepare python environment

# In[94]:


import numpy as np
import pandas as pd
import xarray as xr
import pyproj as pp
from datetime import datetime


# In[307]:


vv = '1.0'


# # Load files

# In[317]:


da = xr.open_dataset(fp+fname)
if verbose: print('loaded file: '+fp+fname)


# In[123]:


ny,nx = da.dims['Along Track'],da.dims['Cross Track']


# ## Load metadata

# In[ ]:


g = pd.read_csv(fp+'Hyperion_attributes.csv')
if verbose: print('loaded metadata file: '+fph+'Hyperion_attributes.csv')


# In[45]:


i = g[g['Entity_ID'].str.contains(fname.split('.')[0])]


# # Interpolate the corners

# In[115]:


geoid = pp.Geod(ellps="WGS84")
x_trackpoints_top = geoid.npts(i['NW_Corne_3'],i['NW_Corne_2'],i['NE_Corne_3'],i['NE_Corne_2'],nx) #lon0,lat0,lon1,lat1
x_trackpoints_bottom = geoid.npts(i['SW_Corne_3'],i['SW_Corne_2'],i['SE_Corne_3'],i['SE_Corne_2'],nx)


# In[116]:


if verbose: print('.. interpolating corners for lat & lon')
lat = np.zeros((ny,nx))
lon = np.zeros((ny,nx))
for j,xt in enumerate(x_trackpoints_top):
    tmp = np.array(geoid.npts(xt[0],xt[1],x_trackpoints_bottom[j][0],x_trackpoints_bottom[j][1],ny))
    lat[:,j] = tmp[:,1]
    lon[:,j] = tmp[:,0]


# In[139]:


attributes = {'Geolocation_CreatedBy':'Samuel LeBlanc',
              'Geolocation_version':vv,
              'Geolocation_CreationDate':str(datetime.now()),
              'Geolocation_method':'Great circle interplation between corners from file: Hyperion_attributes.csv, cross-track first, then along track, using pyproj, WSG84',
              'FieldInfo':'see: https://lta.cr.usgs.gov/DD/EO1.html'}


# In[140]:


da['Latitude'] = xr.DataArray(lat,dims=['Along Track','Cross Track'],attrs=attributes)
da['Longitude'] = xr.DataArray(lon,dims=['Along Track','Cross Track'],attrs=attributes)


# # Get the time of each along track line

# In[243]:


start = pd.to_datetime(i['Scene_Star'],format='%Y:%j:%H:%M:%S.%f').values[0]
stop = pd.to_datetime(i['Scene_Stop'],format='%Y:%j:%H:%M:%S.%f').values[0]
dt = np.linspace(start.astype(int),stop.astype(int),ny)
time = [datetime.utcfromtimestamp(dti*1e-9) for dti in dt]


# In[270]:


time_attrs = {'Geolocation_CreatedBy':'Samuel LeBlanc',
              'Geolocation_version':vv,
              'Geolocation_CreationDate':str(datetime.now()),
              'time_method':'from Scene_start and Scene_stop from file Hyperion_attributes.csv'}


# In[271]:


da['time'] = xr.DataArray(time,dims=['Along Track'],attrs=time_attrs)


# # Get the view angles

# In[193]:


distance_sat_to_earth = 705000.0 #m average, could be better by using the two line element orbit descriptor
#vza = arctan(tan(look_angle)+dist_from_center/distance_sat_to_earth)


# In[192]:


# get the distance from the center point
if verbose: print('.. calculating view angles')
ix = [int(nx/2)] #find the corsstrack center point
vza = np.zeros((ny,nx)) 
vaa = np.zeros((ny,nx)) #view azimuth angle
for iy in range(ny):
    faa_tmp,baa_tmp,d_tmp = geoid.inv(lon[iy,ix*nx],lat[iy,ix*nx],lon[iy,:],lat[iy,:]) #forward az, back azi, dist in m    
    vza[iy,:] = np.rad2deg(np.arctan(np.tan(np.deg2rad(float(i['Look_Angle'])))+d_tmp/distance_sat_to_earth))
    vaa[iy,:] = 90.0-faa_tmp


# In[194]:


view_angles_attrs = {'Geolocation_CreatedBy':'Samuel LeBlanc',
                     'Geolocation_version':vv,
                     'Geolocation_CreationDate':str(datetime.now()),
                     'ViewAngle_method':'using pyproj to calculate differences in look angle from center of crosstrack to each pixel. View azimuth angle calculated from normal of along track',
                     'Additional_info':'see https://www.usgs.gov/centers/eros/look-angles-and-coverage-area'}


# In[196]:


da['ViewZenithAngle'] = xr.DataArray(vza,dims=['Along Track','Cross Track'],attrs=view_angles_attrs)
da['ViewAzimuthAngle'] = xr.DataArray(vaa,dims=['Along Track','Cross Track'],attrs=view_angles_attrs)


# # Get the sun angles

# In[ ]:


def get_sza_azi(lat,lon,datetimet,alt=None,return_sunearthfactor=False,return_sunf_and_dec=False):
    """
    Program wrapper for pyephem.Sun to get the solar zenith angle and the solar azimuth angle
    can use inputs of list or numpy arrays
    require input of lat,lon,datetimet
    optional input of altitutde (in meters)
    optional output of sun earth distance factor if return_sunearthfactor is set to True
    optional output of sun earth distance factor and sun declination if return_sunf_and_dec is set to True
    """
    import ephem
    from numpy import pi,isscalar
    sun = ephem.Sun()
    obs = ephem.Observer()
    if isscalar(lat):
        if isscalar(datetimet):
            lat = [lat]
            lon = [lon]
            datetime = [datetimet]
        else:
            lati = [lat for i in range(len(datetimet))]
            loni = [lon for i in range(len(datetimet))]
            lat,lon = lati,loni
    n = len(lat)
    sza = []
    azi = []
    sunf = []
    dec = []
    for i in range(n):
        obs.lat,obs.lon,obs.date = lat[i]/180.0*pi,lon[i]/180.0*pi,datetimet[i]
        if alt:
            obs.elevation = alt
        sun.compute(obs)
        sza.append(90.0-sun.alt*180/pi)
        azi.append(sun.az*180/pi)
        sunf.append(1.0/(sun.earth_distance**2))
        dec.append(sun.dec*180.0/pi)
    if return_sunf_and_dec:
        return sza,azi,sunf,dec
    elif return_sunearthfactor:
        return sza,azi,sunf
    else:
        return sza,azi


# In[286]:


if verbose: print('.. calculating sun angles')
sza = np.zeros((ny,nx))
azi = np.zeros((ny,nx))
for j in range(nx):
    sza_tmp, azi_tmp = get_sza_azi(lat[:,j],lon[:,j],time)
    sza[:,j] = np.array(sza_tmp)
    azi[:,j] = np.array(azi_tmp)


# In[289]:


sun_attrs = {'Geolocation_CreatedBy':'Samuel LeBlanc',
             'Geolocation_version':vv,
             'Geolocation_CreationDate':str(datetime.now()),
             'SolarAngle_method':'using pyephem to calculate solar zenith and azimuth angle from interpolated lat-lon and time positions'}


# In[290]:


da['SolarZenithAngle'] = xr.DataArray(sza,dims=['Along Track','Cross Track'],attrs=sun_attrs)
da['SolarAzimuthAngle'] = xr.DataArray(azi,dims=['Along Track','Cross Track'],attrs=sun_attrs)


# # Save to file

# In[296]:


if only_geo:
    du = da.drop(['Image','Spectral Center Wavelengths','Spectral Bandwidths','Gain Coefficients','Flag Mask'])
    new_path = fp_out+fname.split('.')[0]+'_'+fname.split('.')[1]+'_geo_only.nc'
    print('Saving to : '+new_path)
    du.to_netcdf(new_path)
else:
    new_path = fp_out+fname.split('.')[0]+'_'+fname.split('.')[1]+'_geo.nc'
    print('Saving to : '+new_path)
    da.to_netcdf(new_path)

