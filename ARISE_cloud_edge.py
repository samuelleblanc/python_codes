# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     ARISE_cloud_edge
# 
# Purpose:  
# 
#     Python scripts that builds the analysis for the ARISE data of the Septemebr 19th, 2014, Flgiht # 13
#     Sets up maps and background information, then proceeds to produce cloud radiative properties from 4STAR transmitted data
# 
# Calling Sequence:
# 
#     python ARISE_cloud_edge.py
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
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - run_kisq_retrieval.py : for the retrieval functions
#     - load_modis.py : for loading modis and other files
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - math
#     - os
#     - gc
#     - pdb : for debugging
#     - datetime
#     - mpl_toolkits
#     - gdal (from osgeo)
#     - plotting_utils (user defined plotting routines)
#     - map_utils, dependent on geopy
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - AMSR-E map files (asi-AMSR2-n6250-20140919-v5.hdf and LongitudeLatitudeGrid-n6250-Arctic.hdf)
#   - MODIS retrievals
#   - CALIPSO overpass file
#   - C130 nav data: arise-C130-Hskping_c130_20140919_RA_Preliminary.ict
#   - Probes data: ARISE-LARGE-PROBES_C130_20140919_R0.ict
#   

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
%matplotlib inline
import numpy as np, h5py
import scipy.io as sio
import math, os, IPython
import Sp_parameters as Sp
IPython.InteractiveShell.cache_size = 0
# set the basic directory path
fp='C:/Users/sleblan2/Research/ARISE/'

# <headingcell level=2>

# Get the AMSR data for 2014-09-19

# <codecell>

import load_modis as lm
famsr = fp+'AMSRE/asi-AMSR2-n6250-20140919-v5.hdf'
fll = fp+'AMSRE/LongitudeLatitudeGrid-n6250-Arctic.hdf'
amsr = lm.load_amsr(famsr,fll)

# <codecell>

def plt_amsr(ax='None'):
    from mpl_toolkits.basemap import Basemap
    if not ax:
        fig = plt.figure()
    m = Basemap(projection='stere',lat_0=57,lon_0=-145,
               llcrnrlon=-170,llcrnrlat=53,
               urcrnrlon=-90,urcrnrlat=80 ,resolution='l')
    m.drawcountries()
    m.fillcontinents(color='grey')
    m.drawmeridians(np.linspace(-90,-200,12),labels=[0,0,0,1])
    m.drawparallels(np.linspace(53,80,10),labels=[1,0,0,0])
    x,y = m(amsr['lon'],amsr['lat'])
    clevels = np.linspace(0,100,21)
    cs = m.contourf(x,y,amsr['ice'],clevels,cmap=plt.cm.gist_earth)
    cbar = m.colorbar(cs)
    cbar.set_label('Ice concentration [\%]')
    return m
m = plt_amsr()

# <codecell>

def plt_amsr_zoom(ax='None',colorbar=True):
    from mpl_toolkits.basemap import Basemap
    if not ax:
        fig = plt.figure()
    m = Basemap(projection='stere',lat_0=72,lon_0=-128,
                llcrnrlon=-136,llcrnrlat=71,
                urcrnrlon=-117.5,urcrnrlat=75,resolution='l')
    m.drawcountries()
    m.fillcontinents(color='grey')
    m.drawmeridians(np.linspace(-115,-145,11),labels=[0,0,0,1])
    m.drawparallels(np.linspace(60,75,16),labels=[1,0,0,0])
    x,y = m(amsr['lon'],amsr['lat'])
    clevels = np.linspace(0,100,21)
    cs = m.contourf(x,y,amsr['ice'],clevels,cmap=plt.cm.gist_earth)
    if colorbar:
        cbar = m.colorbar(cs)
        cbar.set_label('Ice concentration [\%]')
    return m
plt_amsr_zoom()

# <headingcell level=2>

# Load C130 nav data

# <codecell>

fnav = fp+'c130/20140919_ARISE_Flight_13/arise-C130-Hskping_c130_20140919_RA_Preliminary.ict'
nav,nav_header = lm.load_ict(fnav,header=True)

# <codecell>

nav_header

# <codecell>

nav['Longitude'][nav['Longitude']==0.0] = np.NaN
nav['Latitude'][nav['Latitude']<0.0] = np.NaN

# <codecell>

plt.figure()
#plt.plot(nav['Longitude'],nav['Latitude'])
ss = plt.scatter(nav['Longitude'],nav['Latitude'],c=nav['Start_UTC'],edgecolor='None',cmap=cm.gist_ncar)
plt.grid()
cbar = plt.colorbar(ss)
cbar.set_label('UTC [h]')
plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.savefig(fp+'plots/20140919_flightpath.png',dpi=600,transparent=True)

# <codecell>

m = plt_amsr()
m.scatter(nav['Longitude'],nav['Latitude'],latlon=True,zorder=10,s=0.5,edgecolor='r')
plt.savefig(fp+'plots/20140919_map_ice_conc.png',dpi=600,transparent=True)

# <codecell>

flt = np.where((nav['Start_UTC']>19.0) & (nav['Start_UTC']<23.0) & (nav['Longitude']<0.0))[0]

# <codecell>

plt.figure()
plt.plot(nav['Longitude'][flt],nav['GPS_Altitude'][flt])
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
plt.ylim([0,7000])
plt.xlim([-138,-128])
ss = plt.scatter(nav['Longitude'][flt],nav['GPS_Altitude'][flt],c=nav['Start_UTC'][flt],edgecolor='None',cmap=cm.gist_ncar)
cbar = plt.colorbar(ss)
cbar.set_label('UTC [h]')
plt.savefig(fp+'plots/20140919_proile_alt.png',dpi=600,transparent=True)

# <codecell>

plt.figure()
plt.plot(nav['Longitude'][flt],nav['GPS_Altitude'][flt])
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
plt.ylim([0,7000])
plt.xlim([-138,-128])
ss = plt.scatter(nav['Longitude'][flt],nav['GPS_Altitude'][flt],c=nav['Static_Air_Temp'][flt],edgecolor='None',cmap=cm.gist_ncar)
cbar = plt.colorbar(ss)
cbar.set_label('Static Air temp [$^{o}$C]')
plt.savefig(fp+'plots/20140919_proile_alt_temp.png',dpi=600,transparent=True)

# <codecell>

plt.figure()
plt.plot(nav['Longitude'][flt],nav['GPS_Altitude'][flt])
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
plt.ylim([0,7000])
plt.xlim([-138,-128])
ss = plt.scatter(nav['Longitude'][flt],nav['GPS_Altitude'][flt],c=nav['Relative_Humidity'][flt],edgecolor='None',cmap=cm.gist_ncar)
cbar = plt.colorbar(ss)
cbar.set_label('Relative humidity [\\%]')
plt.savefig(fp+'plots/20140919_proile_alt_RH.png',dpi=600,transparent=True)

# <headingcell level=2>

# Load Cloud probe data

# <codecell>

fprobe = fp+'c130/20140919_ARISE_Flight_13/ARISE-LARGE-PROBES_C130_20140919_R0.ict'
probe,prb_header = lm.load_ict(fprobe,header=True)

# <codecell>

prb_header

# <codecell>

flt_prb = np.where((probe['UTC_mid']>19.0) & (probe['UTC_mid']<23.0))
probe['TWC_gm3'][probe['TWC_gm3']<0.0] = np.NaN
feet2meter = 0.3048

# <codecell>

plt.figure()
#plt.plot(nav['Longitude'][flt],nav['GPS_Altitude'][flt])
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
ss = plt.scatter(probe['Longitude_deg'][flt_prb],probe['PressAlt_ft'][flt_prb]*feet2meter,
                 c=probe['nCDP_cm3'][flt_prb],edgecolor='None',
                 cmap=cm.gist_ncar_r)
plt.plot(probe['Longitude_deg'][flt_prb],probe['PressAlt_ft'][flt_prb]*feet2meter,c='k',linewidth=0.3)
plt.ylim([0,7000])
plt.xlim([-138,-128])
plt.title('C130 profile on 2014-09-19')
cbar = plt.colorbar(ss)
cbar.set_label('Drop number concentration [cm$^{-3}$]')
plt.savefig(fp+'plots/20140919_proile_alt_ndrop.png',dpi=600,transparent=True)

# <headingcell level=2>

# Load the MODIS cloud properties

# <markdowncell>

# Get the data from MODIS for the proper day of year

# <codecell>

from datetime import datetime
datetime(2014,9,19).timetuple().tm_yday

# <codecell>

fmodis_aqua = fp+'MODIS/MYD06_L2.A2014262.1955.006.2014282222048.hdf'
fmodis_aqua_geo = fp+'MODIS/MYD03.A2014262.1955.006.2014272205731.hdf'
fmodis_terra = fp+'MODIS/MOD06_L2.A2014262.2115.006.2015077072247.hdf'
fmodis_terra_geo = fp+'MODIS/MOD03.A2014262.2115.006.2014272205802.hdf'
fviirs = fp+'' #not yet

# <codecell>

aqua,aqua_dicts = lm.load_modis(fmodis_aqua_geo,fmodis_aqua)

# <codecell>

terra,terra_dicts = lm.load_modis(fmodis_terra_geo,fmodis_terra)

# <codecell>

flt_aqua = np.where(nav['Start_UTC']<19.9)

# <codecell>

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
m = plt_amsr(ax=ax1)
clevels = np.linspace(0,80,41)
cc = m.contourf(aqua['lon'],aqua['lat'],aqua['tau'],clevels,latlon=True,cmap=plt.cm.rainbow,extend='max',zorder=10)
cbar = plt.colorbar(cc,orientation='horizontal')
cbar.set_label('Aqua $\\tau$')
m.scatter(nav['Longitude'][flt_aqua],nav['Latitude'][flt_aqua],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.title('19:55 UTC')

ax2 = fig.add_subplot(1,2,2)
m2 = plt_amsr(ax=ax2)
clevels2 = np.linspace(0,60,31)
cc = m2.contourf(aqua['lon'],aqua['lat'],aqua['ref'],clevels2,latlon=True,cmap=plt.cm.gnuplot2,extend='max',zorder=10)
cbar = plt.colorbar(cc,orientation='horizontal')
cbar.set_label('Aqua r$_{eff}$ [$\\mu$m]')
m2.scatter(nav['Longitude'][flt_aqua],nav['Latitude'][flt_aqua],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.title('19:55 UTC')

plt.savefig(fp+'plots/20140919_map_aqua.png',dpi=600,transparent=True)

# <codecell>

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
m = plt_amsr_zoom(ax=ax1,colorbar=False)
clevels = np.linspace(0,80,41)
cc = m.contourf(aqua['lon'],aqua['lat'],aqua['tau'],clevels,latlon=True,cmap=plt.cm.rainbow,extend='max',zorder=10)
cbar = plt.colorbar(cc,orientation='horizontal')
cbar.set_label('Aqua $\\tau$')
m.scatter(nav['Longitude'][flt_aqua],nav['Latitude'][flt_aqua],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.title('19:55 UTC')

ax2 = fig.add_subplot(1,2,2)
m2 = plt_amsr_zoom(ax=ax2,colorbar=False)
clevels2 = np.linspace(0,60,31)
cc = m2.contourf(aqua['lon'],aqua['lat'],aqua['ref'],clevels2,latlon=True,cmap=plt.cm.gnuplot2,extend='max',zorder=10)
cbar = plt.colorbar(cc,orientation='horizontal')
cbar.set_label('Aqua r$_{eff}$ [$\\mu$m]')
m2.scatter(nav['Longitude'][flt_aqua],nav['Latitude'][flt_aqua],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.title('19:55 UTC')

plt.savefig(fp+'plots/20140919_map_zoom_aqua.png',dpi=600,transparent=True)

# <codecell>

flt_terra = np.where(nav['Start_UTC']<21.25)

# <codecell>

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
m = plt_amsr(ax=ax1)
clevels = np.linspace(0,80,41)
cc = m.contourf(terra['lon'],terra['lat'],terra['tau'],clevels,latlon=True,cmap=plt.cm.rainbow,extend='max',zorder=10)
cbar = plt.colorbar(cc,orientation='horizontal')
cbar.set_label('Terra $\\tau$')
m.scatter(nav['Longitude'][flt_terra],nav['Latitude'][flt_terra],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.title('21:15 UTC')

ax2 = fig.add_subplot(1,2,2)
m2 = plt_amsr(ax=ax2)
clevels2 = np.linspace(0,60,31)
cc = m2.contourf(terra['lon'],terra['lat'],terra['ref'],clevels2,latlon=True,cmap=plt.cm.gnuplot2,extend='max',zorder=10)
cbar = plt.colorbar(cc,orientation='horizontal')
cbar.set_label('Terra r$_{eff}$ [$\\mu$m]')
m2.scatter(nav['Longitude'][flt_terra],nav['Latitude'][flt_terra],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.title('21:15 UTC')

plt.savefig(fp+'plots/20140919_map_terra.png',dpi=600,transparent=True)

# <codecell>

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
m = plt_amsr_zoom(ax=ax1,colorbar=False)
clevels = np.linspace(0,80,41)
cc = m.contourf(terra['lon'],terra['lat'],terra['tau'],clevels,latlon=True,cmap=plt.cm.rainbow,extend='max',zorder=10)
cbar = plt.colorbar(cc,orientation='horizontal')
cbar.set_label('Terra $\\tau$')
m.scatter(nav['Longitude'][flt_terra],nav['Latitude'][flt_terra],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.title('21:15 UTC')

ax2 = fig.add_subplot(1,2,2)
m2 = plt_amsr_zoom(ax=ax2,colorbar=False)
clevels2 = np.linspace(0,60,31)
cc = m2.contourf(terra['lon'],terra['lat'],terra['ref'],clevels2,latlon=True,cmap=plt.cm.gnuplot2,extend='max',zorder=10)
cbar = plt.colorbar(cc,orientation='horizontal')
cbar.set_label('Terra r$_{eff}$ [$\\mu$m]')
m2.scatter(nav['Longitude'][flt_terra],nav['Latitude'][flt_terra],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.title('21:15 UTC')

plt.savefig(fp+'plots/20140919_map_zoom_terra.png',dpi=600,transparent=True)

# <headingcell level=2>

# Load MERRA Reanalysis

# <codecell>

fmerra = fp+'MERRA\MERRA300.prod.assim.inst6_3d_ana_Nv.20140919.SUB.hdf'
merra,merra_dict = lm.load_hdf_sd(fmerra)

# <codecell>

merra_dict['levels']

# <codecell>

def p2alt(p):
    "convert pressure (in hpa) to altitude (in meter)"
    return (1.0-np.power(p/1013.25,1.0/5.25588))/2.25577E-5

# <codecell>

p2alt(merra['levels'][-10])

# <codecell>

fig,ax = plt.subplots(1,1)
plt.plot(merra['levels'],'b+-')
ax2 = ax.twinx()
ax2.plot(p2alt(merra['levels']),'g+-')
ax2.set_ylabel('Altitude [m]',color='g')
ax.set_ylabel('Pressure [hPa]',color='b')

# <codecell>

plt.figure(figsize=(5,10))
m = plt_amsr()
mlon = merra['latitude'][:,np.newaxis]*0.0+merra['longitude'][np.newaxis,:]
mlat = merra['latitude'][:,np.newaxis]+merra['longitude'][np.newaxis,:]*0.0
x,y = m(mlon,mlat)
q = plt.quiver(x[::4,::4],y[::4,::4],merra['u'][3,-10,::4,::4],merra['v'][3,-10,::4,::4],zorder=10,color='r')
x0,y0 = m(-185,79)
plt.quiverkey(q,x0,y0,10,'10 m/s',coordinates='data',color='r')
m.scatter(nav['Longitude'],nav['Latitude'],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.savefig(fp+'plots/20140919_map_wind.png',dpi=600,transparent=True)

# <codecell>

merra['longitude'][50:78]

# <codecell>

print merra['latitude'].shape
print merra['v'].shape

# <codecell>

plt.figure(figsize=(5,10))
m = plt_amsr_zoom()
mlon = merra['latitude'][40:60,np.newaxis]*0.0+merra['longitude'][np.newaxis,50:78]
mlat = merra['latitude'][40:60,np.newaxis]+merra['longitude'][np.newaxis,50:78]*0.0
x,y = m(mlon,mlat)
q = plt.quiver(x,y,merra['u'][3,-10,40:60,50:78],merra['v'][3,-10,40:60,50:78],zorder=10,color='r')
x0,y0 = m(-133.5,74.6)
plt.quiverkey(q,x0,y0,10,'10 m/s',coordinates='data',color='r')
m.scatter(nav['Longitude'],nav['Latitude'],latlon=True,zorder=10,s=0.5,edgecolor='k')
plt.savefig(fp+'plots/20140919_map_zoom_wind.png',dpi=600,transparent=True)

# <codecell>

plt.quiver(mlon,mlat,merra['u'][3,-10,:,:],merra['v'][3,-10,:,:])
plt.ylabel('Latitude')
plt.xlabel('Longitude')

# <headingcell level=2>

# Plotting of profiles

# <markdowncell>

# First get the indices in ice concentration linked to the flight path

# <codecell>

import map_utils as mu

# <codecell>

points_lat = [73.2268,71.817716]
points_lon = [-135.3189167,-128.407833]

# <codecell>

from scipy import interpolate
flat = interpolate.interp1d([0,100],points_lat)
flon = interpolate.interp1d([0,100],points_lon)
path_lat = flat(np.arange(100))
path_lon = flon(np.arange(100))

# <codecell>

plt.plot(path_lon,path_lat)

# <codecell>

amsr_ind = mu.map_ind(amsr['lon']-180.0,amsr['lat'],path_lon,path_lat)

# <codecell>

ind = np.zeros((2,len(path_lat)), dtype=numpy.int)
for i,x in enumerate(path_lat):
    y = path_lon[i]
    ind[:,i] = np.unravel_index(np.nanargmin(np.square(amsr['lat']-x)+np.square(amsr['lon']-180.0-y)),amsr['lat'].shape)

# <codecell>

ind[:,45]

# <codecell>

plt.plot(path_lon,amsr['ice'][ind[0,:],ind[1,:]])

# <codecell>

ind = np.unravel_index(np.nanargmin(abs(amsr['lat']-x)+abs(amsr['lon']-180.0-x)),amsr['lat'].shape)

# <codecell>

ind[:,0]

# <codecell>

amsr['ice'][ind[0,:],ind[1,:]]

# <codecell>

amsr['lon'].shape

# <codecell>

print amsr['lon'][854,908]-180
print amsr['lat'][854,908]

# <codecell>

plt.plot(amsr['lon'].reshape(amsr['lon'].size)-180,amsr['lat'].reshape(amsr['lon'].size),'g+-')
plt.xlim([-160.5,-160])
plt.ylim([55,55.5])

# <codecell>

plt.figure()
plt.plot(nav['Longitude'][flt],nav['GPS_Altitude'][flt])
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
plt.ylim([0,7000])
plt.xlim([-138,-128])
ss = plt.scatter(nav['Longitude'][flt],nav['GPS_Altitude'][flt],c=nav['Start_UTC'][flt],edgecolor='None',cmap=cm.gist_ncar)
si = plt.scatter(path_lon,path_lon*0.0,c=amsr['ice'][amsr_ind[0,:],amsr_ind[1,:]],edgecolor='None',cmap=cm.gist_earth)
cbar = plt.colorbar(ss)
cbar.set_label('UTC [h]')
#cbar = plt.colorbar(si,orientation='horizontal')
#cbar.set_label('Ice Concentration [\\%%]')

plt.savefig(fp+'plots/20140919_proile_alt_ice.png',dpi=600,transparent=True)

# <codecell>

plt.plot(path_lon,amsr['ice'][amsr_ind[0,:],amsr_ind[1,:]])

# <headingcell level=2>

# Load the 4STAR data

# <codecell>

fstar = fp+'starzen/20140919starzen.mat'
star = sio.loadmat(fstar)
star.keys()

# <codecell>

star['iset']

# <codecell>

star['utc'] = lm.toutc(lm.mat2py_time(star['t']))
star['utcr'] = lm.toutc(lm.mat2py_time(star['t_rad']))

# <codecell>

flt_star_ice = np.where((star['utcr']>19.0) & (star['utcr']<23.0) & (star['Lon'][star['iset']]<-133) & (star['Alt'][star['iset']]<900.0))[0]
flt_star_wat = np.where((star['utcr']>19.0) & (star['utcr']<23.0) & (star['Lon'][star['iset']]>-133) & (star['Alt'][star['iset']]<900.0))[0]

# <codecell>

star['utc'].shape

# <codecell>

flt_star_ice

# <codecell>

fig = plt.figure()
fig.add_subplot(1,2,1)
plt.plot(star['utcr'][flt_star_ice],star['rads'][flt_star_ice,400],'b+',label='Over ice')
plt.plot(star['utcr'][flt_star_wat],star['rads'][flt_star_wat,400],'g+',label='Over water')
plt.ylabel('Zenith Radiance')
plt.xlabel('UTC [h]')
plt.ylim([0,250])
plt.title('at %.3f $\mu$m' % star['w'][0,400])
plt.legend(frameon=False)
fig.add_subplot(1,2,2)
plt.plot(star['utcr'][flt_star_ice],star['rads'][flt_star_ice,1200],'b+',label='Over ice')
plt.plot(star['utcr'][flt_star_wat],star['rads'][flt_star_wat,1200],'g+',label='Over water')
plt.ylabel('Zenith Radiance')
plt.xlabel('UTC [h]')
plt.ylim([0,60])
plt.title('at %.3f $\mu$m' % star['w'][0,1200])
plt.legend(frameon=False)

# <codecell>


