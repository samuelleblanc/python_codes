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
import scipy
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

# <markdowncell>

# To Note from Lee Thornhill about this version of the archived data, received (2015-05-27 5:03 am PST):
# 
# The CDP size distributions are in dN, not dlogD. I make a mistake on that label. Usually I put size distributions into the archive as dNdlogD, but I missed a step in doing that. I will correct that with the next archive file update. 

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

# <markdowncell>

# plotting of drop size distribution

# <codecell>

bin_diameters = np.array([3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50])

# <codecell>

nd_dist = np.vstack((probe['CDP01_dNdlogD'],probe['CDP02_dNdlogD'],probe['CDP03_dNdlogD'],
                     probe['CDP04_dNdlogD'],probe['CDP05_dNdlogD'],probe['CDP06_dNdlogD'],
                     probe['CDP07_dNdlogD'],probe['CDP08_dNdlogD'],probe['CDP09_dNdlogD'],
                     probe['CDP10_dNdlogD'],probe['CDP11_dNdlogD'],probe['CDP12_dNdlogD'],
                     probe['CDP13_dNdlogD'],probe['CDP14_dNdlogD'],probe['CDP15_dNdlogD'],
                     probe['CDP16_dNdlogD'],probe['CDP17_dNdlogD'],probe['CDP18_dNdlogD'],
                     probe['CDP19_dNdlogD'],probe['CDP20_dNdlogD'],probe['CDP21_dNdlogD'],
                     probe['CDP22_dNdlogD'],probe['CDP23_dNdlogD'],probe['CDP24_dNdlogD'],
                     probe['CDP25_dNdlogD'],probe['CDP26_dNdlogD'],probe['CDP27_dNdlogD'],
                     probe['CDP28_dNdlogD'],probe['CDP29_dNdlogD'],probe['CDP30_dNdlogD']))

# <codecell>

nd_dist[nd_dist<0.5] = np.nan

# <codecell>

clevels = np.arange(0,51)
vv = plt.contourf(probe['UTC_mid'],bin_diameters,nd_dist,clevels,cmap=plt.cm.gist_ncar_r,extend='both')
cbar = plt.colorbar(vv)
cbar.set_label('dN/dlogD [cm$^{-3}$]')
plt.xlabel('UTC [H]')
plt.ylabel('$D_{eff}$ [$\\mu$m]')
plt.title('CDP cloud drop size density')
plt.savefig(fp+'plots/20140919_utc_dndlogd.png',dpi=600,transparent=True)

# <codecell>

clevels = np.arange(0,51)
vv = plt.contourf(probe['UTC_mid'],bin_diameters,nd_dist,clevels,cmap=plt.cm.gist_ncar_r,extend='both')
cbar = plt.colorbar(vv)
cbar.set_label('dN/dlogD [cm$^{-3}$]')
plt.xlabel('UTC [H]')
plt.ylabel('$D_{eff}$ [$\\mu$m]')
plt.title('CDP cloud drop size density')
plt.xlim([19.5,23])
plt.savefig(fp+'plots/20140919_utc_dndlogd_zoom.png',dpi=600,transparent=True)

# <headingcell level=2>

# Convert the drop size distribution to $r_{eff}$

# <markdowncell>

# Trying to replicate the $r_{eff}$ equation:
# 
# $$r_{eff}=\frac{\int^{\infty}_{0}Q_{ext}(r)r^{3}n(r)dr}{\int^{\infty}_{0}Q_{ext}(r)r^{2}n(r)dr}$$
# 
# 
# With $Q_{ext}(r)$ representing the extinction efficiency, $r$ the radius of the particle, $n(r)$ the drop size distribution

# <markdowncell>

# Load $Q_{ext}(r)$ for mie

# <markdowncell>

# For now use Q_ext = parameterized values with exponenetial decrease determined from Mie_Calc ipython notebook for lambda at 1.70 $\mu$m 

# <codecell>

nd_dist.shape

# <codecell>

nd_dist[:,1000]

# <codecell>

from Sp_parameters import nanmasked
def calc_ref(nd,diameter):
    " calculate the effective radius using trapezoid integration rule "
    nda,ndi = nanmasked(nd)
    popt = [ 0.87495122 ,0.22795525 ,2.08105225] # from Mie_Calc.py
    def qext(x, a, c, d):
        return a*np.exp(-c*x)+d
    if not any(ndi):
        re = np.NaN
    else:
        re = np.trapz(qext(diameter[ndi]/2.0,*popt)*(diameter[ndi]/2.0)**3.0*nda)/ \
            np.trapz(qext(diameter[ndi]/2.0,*popt)*(diameter[ndi]/2.0)**2.0*nda)
    return re

# <codecell>

print calc_ref(nd_dist[:,14000],bin_diameters)

# <codecell>

ref = np.zeros(len(nd_dist[0,:]))
for i in xrange(len(nd_dist[0,:])):
    ref[i] = calc_ref(nd_dist[:,i],bin_diameters)

# <codecell>

plt.plot(probe['UTC_mid'],ref)
plt.title('Calculated effective radius from cloud probes')
plt.xlabel('UTC [H]')
plt.ylabel('$r_{eff}$ [$\\mu$m]')
plt.savefig(fp+'plots/20140919_utc_probes_ref.png',dpi=600,transparent=True)

# <codecell>

probe_ref = Sp.smooth(ref,10,nan=False)

# <codecell>

plt.figure()
#plt.plot(nav['Longitude'][flt],nav['GPS_Altitude'][flt])
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
ss = plt.scatter(probe['Longitude_deg'][flt_prb],probe['PressAlt_ft'][flt_prb]*feet2meter,
                 c=probe_ref[flt_prb],edgecolor='None',
                 cmap=cm.gist_ncar_r)
plt.plot(probe['Longitude_deg'][flt_prb],probe['PressAlt_ft'][flt_prb]*feet2meter,c='k',linewidth=0.3)
plt.ylim([0,7000])
plt.xlim([-138,-128])
plt.title('Profile of calculated $r_{eff}$ from cloud probes')
cbar = plt.colorbar(ss)
cbar.set_label('$r_{eff}$ [$\\mu$m]')
plt.savefig(fp+'plots/20140919_proile_alt_ref_calc.png',dpi=600,transparent=True)

# <codecell>

plt.figure()
plt.plot(probe['Longitude_deg'][flt_prb],probe_ref[flt_prb],label='Raw')
plt.plot(probe['Longitude_deg'][flt_prb],Sp.smooth(probe_ref[flt_prb],10,nan=False),label='Smoothed 10x',linewidth=2)
plt.plot(probe['Longitude_deg'][flt_prb],Sp.smooth(probe_ref[flt_prb],30,nan=False),label='Smoothed 30x',linewidth=2)
plt.title('Cloud probes calculated $r_{eff}$')
plt.xlabel('Longitude')
plt.ylabel('$r_{eff}$ [$\\mu$m]')
plt.legend(frameon=False)
plt.savefig(fp+'plots/20140919_lon_ref_calc.png',dpi=600,transparent=True)

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

terra_dicts['tau']

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

merra_dict['delp']

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
#m.contour(mlon,mlat,merra['delp'][3,-10,:,:],latlon=True,zorder=10)
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

# <markdowncell>

# Create a lat lon slice to use for correlating the sea ice concentration

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

ind = np.zeros((2,len(path_lat)), dtype=numpy.int)
for i,x in enumerate(path_lat):
    y = path_lon[i]
    ind[:,i] = np.unravel_index(np.nanargmin(np.square(amsr['lat']-x)+np.square(amsr['lon']-360.0-y)),amsr['lat'].shape)

# <codecell>

plt.plot(path_lon,amsr['ice'][ind[0,:],ind[1,:]])
plt.ylabel('Ice Concentration [\%]')
plt.xlabel('Longitude')
plt.title('Below flight path ice concentration')

# <codecell>

plt.plot(amsr['lon'].reshape(amsr['lon'].size)-360,amsr['lat'].reshape(amsr['lon'].size),'g+-')
plt.plot(amsr['lon'][ind[0,:],ind[1,:]]-360.0,amsr['lat'][ind[0,:],ind[1,:]],'r+')
plt.plot(nav['Longitude'],nav['Latitude'],'k')
plt.xlim([-138,-128])
plt.ylim([71,74])
clevels = np.linspace(0,100,21)
plt.contourf(amsr['lon']-360.0,amsr['lat'],amsr['ice'],clevels,cmap=plt.cm.gist_earth)
plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.title('AMSR ice concentration for along C130 flight path')

# <codecell>

plt.figure()
plt.plot(nav['Longitude'][flt],nav['GPS_Altitude'][flt])
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
plt.ylim([-200,7000])
plt.xlim([-138,-128])
ss = plt.scatter(nav['Longitude'][flt],nav['GPS_Altitude'][flt],c=nav['Start_UTC'][flt],edgecolor='None',cmap=cm.gist_ncar)
si = plt.scatter(path_lon,path_lon*0.0-100,c=amsr['ice'][ind[0,:],ind[1,:]],
                 marker='s',edgecolor='None',cmap=cm.gist_earth,s=30,clip_on=False)
cbar = plt.colorbar(ss)
cbar.set_label('UTC [h]')
cbar = plt.colorbar(si,orientation='horizontal')
cbar.set_label('Ice Concentration [\%]')

plt.savefig(fp+'plots/20140919_proile_alt_ice.png',dpi=600,transparent=True)

# <codecell>

plt.figure()
plt.plot(probe['Longitude_deg'][flt_prb],probe['PressAlt_ft'][flt_prb]*feet2meter,'k',linewidth=0.6)
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
plt.ylim([-200,7000])
plt.xlim([-138,-128])
#ss = plt.scatter(nav['Longitude'][flt],nav['GPS_Altitude'][flt],c=nav['Start_UTC'][flt],edgecolor='None',cmap=cm.gist_ncar)
si = plt.scatter(path_lon,path_lon*0.0-100,c=amsr['ice'][ind[0,:],ind[1,:]],
                 marker='s',edgecolor='None',cmap=cm.gist_earth,s=30,clip_on=False)
#cbar = plt.colorbar(ss)
#cbar.set_label('UTC [h]')
ss = plt.scatter(probe['Longitude_deg'][flt_prb],probe['PressAlt_ft'][flt_prb]*feet2meter,
                 c=probe['nCDP_cm3'][flt_prb],edgecolor='None',
                 cmap=cm.gist_ncar_r)
cbar = plt.colorbar(ss)
cbar.set_label('Drop number concentration [cm$^{-3}$]')

cbar = plt.colorbar(si,orientation='horizontal')
cbar.set_label('Ice Concentration [\%]')

plt.savefig(fp+'plots/20140919_profile_alt_cloud_ice.png',dpi=600,transparent=True)

# <codecell>

plt.figure()
plt.plot(probe['Longitude_deg'][flt_prb],probe['PressAlt_ft'][flt_prb]*feet2meter,'k',linewidth=0.6)
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
plt.ylim([-200,2000])
plt.xlim([-136,-128])
#ss = plt.scatter(nav['Longitude'][flt],nav['GPS_Altitude'][flt],c=nav['Start_UTC'][flt],edgecolor='None',cmap=cm.gist_ncar)
si = plt.scatter(path_lon,path_lon*0.0-100,c=amsr['ice'][ind[0,:],ind[1,:]],
                 marker='s',edgecolor='None',cmap=cm.gist_earth,s=30,clip_on=False)
#cbar = plt.colorbar(ss)
#cbar.set_label('UTC [h]')
ss = plt.scatter(probe['Longitude_deg'][flt_prb],probe['PressAlt_ft'][flt_prb]*feet2meter,
                 c=probe['nCDP_cm3'][flt_prb],edgecolor='None',
                 cmap=cm.gist_ncar_r)
cbar = plt.colorbar(ss)
cbar.set_label('Drop number concentration [cm$^{-3}$]')

cbar = plt.colorbar(si,orientation='horizontal')
cbar.set_label('Ice Concentration [\%]')

plt.savefig(fp+'plots/20140919_proile_zoom_alt_cloud_ice.png',dpi=600,transparent=True)

# <markdowncell>

# Now get a MERRA profile of winds on top of this.

# <codecell>


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

# <headingcell level=3>

# Define 4STAR zenith data as an SP object and prepare for calculations

# <codecell>

import Sp_parameters as Sp
reload(Sp)

# <codecell>

star['good'] = np.where((star['utcr']>19.0) & (star['utcr']<23.0) & (star['Alt'][star['iset'][:,0],0]<900.0))[0]
len(star['good'])

# <codecell>

stars = Sp.Sp(star)
stars.params()

# <headingcell level=3>

# Now load the different luts

# <codecell>

sice = sio.idl.readsav(fp+'model/sp_v2_20140919_ice_4STAR.out')
swat = sio.idl.readsav(fp+'model/sp_v2_20140919_wat_4STAR.out')
print sice.keys()
print swat.keys()

# <codecell>

lutice = Sp.Sp(sice)

# <codecell>

lutwat = Sp.Sp(swat)

# <codecell>

lutice.params()
lutice.param_hires()
lutice.sp_hires()

# <codecell>

lutwat.params()
lutwat.param_hires()
lutwat.sp_hires()

# <headingcell level=3>

# Now load and run the retrieval on 4STAR data subset

# <codecell>

import run_kisq_retrieval as rk
reload(rk)

# <codecell>

(stars.icetau,stars.iceref,stars.icephase,stars.iceki) = rk.run_retrieval(stars,lutice,force_liq=True)

# <codecell>

(stars.wattau,stars.watref,stars.watphase,stars.watki) = rk.run_retrieval(stars,lutwat,force_liq=True)

# <headingcell level=3>

# Plot the resulting retrieved properties

# <markdowncell>

# Assure that the retrieved properties make sense first, remove bad data

# <codecell>

stars.icetau[stars.icetau>=80] = np.nan
stars.wattau[stars.wattau>=80] = np.nan

# <codecell>

stars.iceref[stars.iceref>=60] = np.nan
stars.watref[stars.watref>=60] = np.nan
stars.iceref[stars.iceref==17] = np.nan
stars.watref[stars.watref==17] = np.nan
#stars.iceref[stars.iceref==30] = np.nan
#stars.watref[stars.watref==30] = np.nan
#stars.iceref[stars.iceref==40] = np.nan
#stars.watref[stars.watref==40] = np.nan

# <markdowncell>

# Now smooth over data to reduce variability

# <codecell>

stars.icetau = Sp.smooth(stars.icetau[:],5,nan=False)
stars.wattau = Sp.smooth(stars.wattau[:],5,nan=False)

# <codecell>

stars.iceref = Sp.smooth(stars.iceref[:],5,nan=False)
stars.watref = Sp.smooth(stars.watref[:],5,nan=False)

# <markdowncell>

# Plot the resulting figures

# <codecell>

plt.plot(stars.lon,stars.icetau,label='All clouds over ice')
plt.plot(stars.lon,stars.wattau,label='All clouds over water')
plt.plot(stars.lon[stars.icephase==0],stars.icetau[stars.icephase==0],'c+',label='Liq clouds over ice')
plt.plot(stars.lon[stars.watphase==0],stars.wattau[stars.watphase==0],'y+',label='Liq clouds over water')
plt.plot(stars.lon[stars.icephase==1],stars.icetau[stars.icephase==1],'cx',label='Ice clouds over ice')
plt.plot(stars.lon[stars.watphase==1],stars.wattau[stars.watphase==1],'yx',label='Ice clouds over water')
plt.ylabel('$\\tau$')
plt.xlabel('Longitude')
plt.legend(frameon=False,loc=2)
plt.savefig(fp+'plots/20140919_tau_retrieved.png',dpi=600,transparent=True)
#savemetapng(fp+'plots/20140919_tau_retrieved.png',theNotebookName)

# <codecell>

from PIL import Image
im2 = Image.open(fp+'plots/20140919_tau_retrieved.png')
print im2.info

# <codecell>

plt.plot(stars.lon,stars.iceref,label='All clouds over ice')
plt.plot(stars.lon,stars.watref,label='All clouds over water')
plt.plot(stars.lon[stars.icephase==0],stars.iceref[stars.icephase==0],'c+',label='Liq clouds over ice')
plt.plot(stars.lon[stars.watphase==0],stars.watref[stars.watphase==0],'y+',label='Liq clouds over water')
plt.plot(stars.lon[stars.icephase==1],stars.iceref[stars.icephase==1],'cx',label='Ice clouds over ice')
plt.plot(stars.lon[stars.watphase==1],stars.watref[stars.watphase==1],'yx',label='Ice clouds over water')
plt.ylabel('r$_{eff}$ [$\\mu$m]')
plt.xlabel('Longitude')
plt.legend(frameon=False,loc=2)

# <headingcell level=3>

# Compile into statistics over water and over ice

# <codecell>

fltice = np.where(stars.lon<-132.0)[0]
fltwat = np.where(stars.lon>-131.0)[0]

# <codecell>

(stars.wattaumask,stars.iwattau) = Sp.nanmasked(stars.wattau[fltwat])
(stars.icetaumask,stars.iicetau) = Sp.nanmasked(stars.icetau[fltice])

# <codecell>

def data2figpoints(x,dx):
    "function to tranform data locations to relative figure coordinates (in fractions of total figure"
    flen = fig.transFigure.transform([1,1])
    bot = ax1.transAxes.transform([0,0])/flen
    top = ax1.transAxes.transform([1,1])/flen
    
    start = ax1.transData.transform([x,0])/flen
    end = ax1.transData.transform([x+dx,0])/flen
    left = start[0]
    bottom = bot[1]
    width = end[0]-start[0]
    height = top[1]-bot[1] 
    return left,bottom,width,height

# <codecell>

def plot_vert_hist(fig,y,pos,ylim,color='grey',label=None,legend=False,onlyhist=True,loc=2):
    "function to plot a 'bean' like vertical histogram"
    import Sp_parameters as Sp
    (ymask,iy) = Sp.nanmasked(y)
    ax = fig.add_axes(data2figpoints(pos,0.4),frameon=False,ylim=ylim)
    ax.tick_params(axis='both', which='both', labelleft='off', labelright='off',bottom='off',top='off',
               labelbottom='off',labeltop='off',right='off',left='off')
    ax.hist(ymask,orientation='horizontal',normed=True,color=color,edgecolor='None',bins=30,alpha=0.5,label=label)
    if onlyhist:
        label_mean = None
        label_median = None
    else:
        label_mean = 'Mean'
        label_median = 'Median'
    ax.axhline(mean(ymask),color='red',linewidth=2,label=label_mean)
    ax.axhline(median(ymask),color='k',linewidth=2,linestyle='--',label=label_median)
    if legend:
        ax.legend(frameon=False,loc=loc)
    ax = fig.add_axes(data2figpoints(pos+0.01,-0.4),frameon=False,ylim=ylim)
    ax.tick_params(axis='both', which='both', labelleft='off', labelright='off',bottom='off',top='off',
                   labelbottom='off',labeltop='off',right='off',left='off')
    ax.hist(ymask,orientation='horizontal',normed=True,color=color,edgecolor='None',bins=30,alpha=0.5)
    ax.axhline(mean(ymask),color='red',linewidth=2)
    ax.axhline(median(ymask),color='k',linewidth=2,linestyle='--')

# <codecell>

fig = plt.figure()
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,40],xlim=[-1,7])
ax1.set_ylabel('$\\tau$')
ax1.set_xticks([1,5])
ax1.set_xticklabels(['Over ice','Over water'])
plot_vert_hist(fig,stars.icetau[fltice],0,[0,40],color='grey',label='Total')
plot_vert_hist(fig,stars.icetau[fltice[stars.icephase[fltice]==1]],1,[0,40],color='red',label='Ice')
plot_vert_hist(fig,stars.icetau[fltice[stars.icephase[fltice]==0]],2,[0,40],color='blue',label='Liquid')
ax.legend()

plot_vert_hist(fig,stars.wattau[fltwat],4,[0,40],color='grey')
plot_vert_hist(fig,stars.wattau[fltwat[stars.watphase[fltwat]==1]],5,[0,40],color='red')
plot_vert_hist(fig,stars.wattau[fltwat[stars.watphase[fltwat]==0]],6,[0,40],color='blue')
ax1.set_title('Distribution of $\\tau$ over different surfaces')
plt.savefig(fp+'plots/20140919_pdf_surf_tau.png',dpi=600,transparent=True)

# <codecell>

fig = plt.figure()
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,60],xlim=[-1,7])
ax1.set_ylabel('r$_{eff}$ [$\\mu$m]')
ax1.set_xticks([1,5])
ax1.set_xticklabels(['Over ice','Over water'])
plot_vert_hist(fig,stars.iceref[fltice],0,[0,60],color='grey',label='Total')
plot_vert_hist(fig,stars.iceref[fltice[stars.icephase[fltice]==1]],1,[0,60],color='red',label='Ice')
plot_vert_hist(fig,stars.iceref[fltice[stars.icephase[fltice]==0]],2,[0,60],color='blue',label='Liquid')

plot_vert_hist(fig,stars.watref[fltwat],4,[0,60],color='grey')
plot_vert_hist(fig,stars.watref[fltwat[stars.watphase[fltwat]==1]],5,[0,60],color='red')
plot_vert_hist(fig,stars.watref[fltwat[stars.watphase[fltwat]==0]],6,[0,60],color='blue')
ax1.set_title('Distribution of Effective Radius over different surfaces')
plt.savefig(fp+'plots/20140919_pdf_surf_ref.png',dpi=600,transparent=True)

# <codecell>

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,30],xlim=[0,2])
ax1.set_ylabel('r$_{eff}$ [$\\mu$m]')
ax1.set_xticks([0.5,1.5])
ax1.set_xticklabels(['Over ice','Over open ocean'])
plot_vert_hist(fig,stars.iceref[fltice],0.5,[0,30],color='grey',legend=True,onlyhist=False)
plot_vert_hist(fig,stars.watref[fltwat],1.5,[0,30],color='grey')

ax1.set_title('Distribution of r$_{eff}$ over different surfaces')
plt.savefig(fp+'plots/20140919_pdf_surf_tot_ref_forceliq.png',dpi=600,transparent=True)

# <codecell>

print np.nanmean(stars.iceref[fltice])
print np.nanmean(stars.watref[fltwat])
print np.nanmean(stars.watref[fltwat])-np.nanmean(stars.iceref[fltice])
print np.nanmean(stars.watref[fltwat])-np.nanmean(stars.iceref[fltice])/np.nanmean(stars.watref[fltwat])*100.0

# <codecell>

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,60],xlim=[0,2])
ax1.set_ylabel('$\\tau$')
ax1.set_xticks([0.5,1.5])
ax1.set_xticklabels(['Over ice','Over open ocean'])
plot_vert_hist(fig,stars.icetau[fltice],0.5,[0,60],color='grey',legend=True,onlyhist=False)
plot_vert_hist(fig,stars.wattau[fltwat],1.5,[0,60],color='grey')
ax1.set_title('Distribution of $\\tau$ over different surfaces')
plt.savefig(fp+'plots/20140919_pdf_surf_tot_tau_forceliq.png',dpi=600,transparent=True)

# <codecell>

ntau,intau = Sp.nanmasked(stars.icetau[fltice])
plt.figure()
plt.hist(ntau,bins=30,color='grey',edgecolor='None',alpha=0.5)
plt.xlim([0,60])
plt.xlabel('$\\tau$')
plt.ylabel('Counts')
plt.axvline(np.mean(ntau),color='r',label='Mean',linewidth=2)
plt.axvline(np.median(ntau),color='k',linestyle='--',label='Median',linewidth=2)
plt.legend(frameon=False)
plt.savefig(fp+'plots/20140919_hist_tau_forceliq.png',dpi=600,transparent=True)

# <codecell>

print np.nanmean(stars.icetau[fltice])
print np.nanmean(stars.wattau[fltwat])
print np.nanmean(stars.wattau[fltwat])-np.nanmean(stars.icetau[fltice])
print (np.nanmean(stars.wattau[fltwat])-np.nanmean(stars.icetau[fltice]))/np.nanmean(stars.wattau[fltwat])*100.0

# <codecell>

print np.nanmean(stars.iceref[fltice])
print np.nanmean(stars.watref[fltwat])
print np.nanmean(stars.watref[fltwat])-np.nanmean(stars.iceref[fltice])
print (np.nanmean(stars.watref[fltwat])-np.nanmean(stars.iceref[fltice]))/np.nanmean(stars.watref[fltwat])*100.0

# <markdowncell>

# Plotting the over ice, tau and reff together, and also over ocean tau and ref together

# <codecell>

fig = plt.figure(figsize=(3,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,60],xlim=[0.5,1.5])
ax1.set_ylabel('$\\tau$')
ax1.set_xticks([1])
ax1.set_xticklabels(['Over ice'])
plot_vert_hist(fig,stars.icetau[fltice],1,[0,60],color='grey',legend=True,onlyhist=False,loc=1)
#plot_vert_hist(fig,stars.wattau[fltwat],1.5,[0,40],color='grey')
ax1.set_title('Optical thickness')
plt.savefig(fp+'plots/20140919_hist_tau_over_ice_forceliq.png',dpi=600,transparent=True)

# <codecell>

fig = plt.figure(figsize=(3,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,60],xlim=[0.5,1.5])
ax1.set_ylabel('$\\tau$')
ax1.set_xticks([1])
ax1.set_xticklabels(['Over Open ocean'])
#plot_vert_hist(fig,stars.icetau[fltice],1,[0,40],color='grey',legend=True,onlyhist=False,loc=1)
plot_vert_hist(fig,stars.wattau[fltwat],1,[0,60],color='grey',legend=True,onlyhist=False,loc=1)
ax1.set_title('Optical thickness')
plt.savefig(fp+'plots/20140919_hist_tau_over_ocean_forceliq.png',dpi=600,transparent=True)

# <codecell>

fig = plt.figure(figsize=(3,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,30],xlim=[0.5,1.5])
ax2 = ax1.twinx()
ax2.set_ylabel('r$_{eff}$ [$\\mu$m]')
ax2.set_ylim(ax1.get_ylim())
ax1.set_xticks([1])
ax1.set_xticklabels(['Over ice'])
ax1.tick_params(axis='y',which='both',labelleft='off', labelright='off')
ax2.tick_params(axis='both', which='both', labelleft='off', labelright='on',bottom='on',top='off',
               labelbottom='on',labeltop='off',right='on',left='off')
plot_vert_hist(fig,stars.iceref[fltice],1,[0,30],color='green',legend=True,onlyhist=False,loc=1)
#plot_vert_hist(fig,stars.wattau[fltwat],1.5,[0,40],color='grey')
ax1.set_title('Effective radius')
plt.savefig(fp+'plots/20140919_hist_ref_over_ice_forceliq.png',dpi=600,transparent=True)

# <codecell>

fig = plt.figure(figsize=(3,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,30],xlim=[0.5,1.5])
ax2 = ax1.twinx()
ax2.set_ylabel('r$_{eff}$ [$\\mu$m]')
ax2.set_ylim(ax1.get_ylim())
ax1.set_xticks([1])
ax1.set_xticklabels(['Over Open ocean'])
ax1.tick_params(axis='y',which='both',labelleft='off', labelright='off')
ax2.tick_params(axis='both', which='both', labelleft='off', labelright='on',bottom='on',top='off',
               labelbottom='on',labeltop='off',right='on',left='off')
#plot_vert_hist(fig,stars.iceref[fltice],1,[0,60],color='green',legend=True,onlyhist=False,loc=1)
plot_vert_hist(fig,stars.watref[fltwat],1,[0,30],color='green',legend=True,onlyhist=False,loc=1)
ax1.set_title('Effective radius')
plt.savefig(fp+'plots/20140919_hist_ref_over_ocean_forceliq.png',dpi=600,transparent=True)

# <markdowncell>

# Testing other types of plots

# <codecell>

import statsmodels.api as sm

# <codecell>

fig = sm.graphics.beanplot([stars.wattaumask,stars.icetaumask],
                           labels=['Over water','Over ice'],
                           jitter=True,
                           plot_opts={'jitter_marker':'+','jitter_marker_size':2})

# <headingcell level=2>

# Now Match MODIS points to flight path

# <codecell>

aqua.keys()

# <codecell>

aqua['lon'].shape

# <codecell>

ind_aqua = np.zeros((2,len(path_lat)), dtype=numpy.int)
for i,x in enumerate(path_lat):
    y = path_lon[i]
    ind_aqua[:,i] = np.unravel_index(np.nanargmin(np.square(aqua['lat']-x)+np.square(aqua['lon']-y)),aqua['lat'].shape)

# <codecell>

ind_terra = np.zeros((2,len(path_lat)), dtype=numpy.int)
for i,x in enumerate(path_lat):
    y = path_lon[i]
    ind_terra[:,i] = np.unravel_index(np.nanargmin(np.square(terra['lat']-x)+np.square(terra['lon']-y)),terra['lat'].shape)

# <headingcell level=2>

# Plotting resulting MODIS cloud properties

# <codecell>

plt.plot(path_lon,aqua['tau'][ind_aqua[0,:],ind_aqua[1,:]],'b+-',label='Aqua')
plt.plot(path_lon,terra['tau'][ind_terra[0,:],ind_terra[1,:]],'g+-',label='Terra')
plt.legend(frameon=False)
plt.ylabel('$\\tau$')
plt.xlabel('Longitude')
plt.title('MODIS Cloud optical depth along flight path')

# <codecell>

plt.plot(path_lon,aqua['ref'][ind_aqua[0,:],ind_aqua[1,:]],'b+-',label='Aqua')
plt.plot(path_lon,terra['ref'][ind_terra[0,:],ind_terra[1,:]],'g+-',label='Terra')
plt.legend(frameon=False)
plt.ylabel('r$_{eff}$ [$\\mu$m]')
plt.xlabel('Longitude')
plt.title('MODIS Cloud effective radius along flight path')

# <codecell>

pathice = np.where(path_lon<-132.0)[0]
pathwat = np.where(path_lon>-131.0)[0]

# <headingcell level=2>

# Build MODIS and 4STAR comparison plot

# <codecell>

fig = plt.figure(figsize=(7,3))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,60],xlim=[-1,7])
ax1.set_ylabel('$\\tau$')
ax1.set_xticks([1,5])
ax1.set_xticklabels(['Over ice','Over water'])
ax1.tick_params(axis='x', which='both', top='off')
plot_vert_hist(fig,stars.icetau[fltice],0,[0,60],color='grey',label='4STAR',legend=True,onlyhist=False,loc=2)
plot_vert_hist(fig,aqua['tau'][ind_aqua[0,pathice],ind_aqua[1,pathice]],1,[0,60],color='blue')
plot_vert_hist(fig,terra['tau'][ind_terra[0,pathice],ind_terra[1,pathice]],2,[0,60],color='green')

plot_vert_hist(fig,stars.wattau[fltwat],4,[0,60],color='grey')
plot_vert_hist(fig,aqua['tau'][ind_aqua[0,pathwat],ind_aqua[1,pathwat]],5,[0,60],color='blue',label='Aqua',legend=True,loc=1)
plot_vert_hist(fig,terra['tau'][ind_terra[0,pathwat],ind_terra[1,pathwat]],6,[0,60],color='green',label='Terra',legend=True,loc=9)
ax1.set_title('Distribution of $\\tau$ over different surfaces')
plt.savefig(fp+'plots/20140919_comp_modis_tau_forceliq.png',dpi=600,transparent=True)

# <codecell>

fig = plt.figure(figsize=(7,3))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,30],xlim=[-1,7])
ax1.set_ylabel('r$_{eff}$ [$\\mu$m]')
ax1.set_xticks([1,5])
ax1.set_xticklabels(['Over ice','Over water'])
ax1.tick_params(axis='x', which='both', top='off')
plot_vert_hist(fig,stars.iceref[fltice],0,[0,30],color='grey',label='4STAR',legend=True,onlyhist=False,loc=2)
plot_vert_hist(fig,aqua['ref'][ind_aqua[0,pathice],ind_aqua[1,pathice]],1,[0,30],color='blue')
plot_vert_hist(fig,terra['ref'][ind_terra[0,pathice],ind_terra[1,pathice]],2,[0,30],color='green')

plot_vert_hist(fig,stars.watref[fltwat],4,[0,30],color='grey')
plot_vert_hist(fig,aqua['ref'][ind_aqua[0,pathwat],ind_aqua[1,pathwat]],5,[0,30],color='blue',label='Aqua',legend=True,loc=1)
plot_vert_hist(fig,terra['ref'][ind_terra[0,pathwat],ind_terra[1,pathwat]],6,[0,30],color='green',label='Terra',legend=True,loc=9)
ax1.set_title('Distribution of Effective Radius over different surfaces')
plt.savefig(fp+'plots/20140919_comp_modis_ref_forceliq.png',dpi=600,transparent=True)

# <headingcell level=2>

# Compare to $r_{eff}$ derived from Cloud probes

# <codecell>

probe.dtype.names

# <codecell>

flt_probe_ice = np.where((probe['UTC_mid']>19.0) & (probe['UTC_mid']<23.0) & (probe['Longitude_deg']<-133.5))[0]
flt_probe_wat = np.where((probe['UTC_mid']>19.0) & (probe['UTC_mid']<23.0) & (probe['Longitude_deg']>-130.0))[0] 

# <codecell>

fig = plt.figure(figsize=(7,3))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,30],xlim=[-1,5])
ax1.set_ylabel('r$_{eff}$ [$\\mu$m]')
ax1.set_xticks([0.5,3.5])
ax1.set_xticklabels(['Over ice','Over water'])
ax1.tick_params(axis='x', which='both', top='off')
plot_vert_hist(fig,stars.iceref[fltice],0,[0,30],color='grey',label='4STAR',legend=True,onlyhist=False,loc=2)
plot_vert_hist(fig,probe_ref[flt_probe_ice],1,[0,30],color='blue')
#plot_vert_hist(fig,terra['ref'][ind_terra[0,pathice],ind_terra[1,pathice]],2,[0,30],color='green')

plot_vert_hist(fig,stars.watref[fltwat],3,[0,30],color='grey')
plot_vert_hist(fig,probe_ref[flt_probe_wat],4,[0,30],color='blue',label='Probes',legend=True,loc=1)
#plot_vert_hist(fig,aqua['ref'][ind_aqua[0,pathwat],ind_aqua[1,pathwat]],5,[0,30],color='blue',label='Aqua',legend=True,loc=1)
#plot_vert_hist(fig,terra['ref'][ind_terra[0,pathwat],ind_terra[1,pathwat]],6,[0,30],color='green',label='Terra',legend=True,loc=9)
ax1.set_title('Distribution of Effective Radius over different surfaces')
plt.savefig(fp+'plots/20140919_comp_probes_ref_forceliq.png',dpi=600,transparent=True)

# <headingcell level=2>

# Testing different programs for metadata saving for pngs and python files

# <codecell>

def getname(): 
    import IPython
    IPython.display.display(IPython.display.Javascript('IPython.notebook.kernel.execute("theNotebookName = " + "\'"+IPython.notebook.notebook_name+"\'");'))
    IPython.display.display(IPython.display.Javascript('IPython.notebook.kernel.execute("theNotebookPath = " + "\'"+IPython.notebook.notebook_path+"\'");'))

# <codecell>

getname()
thisfilepath = os.getcwd()+os.path.sep+theNotebookPath+theNotebookName
print thisfilepath
fstat = os.stat(thisfilepath)
import datetime
time_of_mod = datetime.datetime.fromtimestamp(fstat.st_mtime).strftime('%Y-%m-%d %H:%M:%S.%f')
import getpass
user = getpass.getuser()

# <codecell>

def savemetapng(filein,Notebookname,Notes=None):
    from PIL import Image
    from PIL import PngImagePlugin
    import datetime
    import getpass
    import os
    
    metadata = dict()
    metadata['Created with'] = 'IPython'+os.sys.version
    metadata['Script Name'] = Notebookname
    metadata['Time created'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
    metadata['Author'] = getpass.getuser()
    if Notes:
        metadata['Notes'] = Notes
    
    im = Image.open(filein)
    meta = PngImagePlugin.PngInfo()

    for x in metadata:
        meta.add_text(x, metadata[x])
    im.save(filein, "png", pnginfo=meta)

# <codecell>


