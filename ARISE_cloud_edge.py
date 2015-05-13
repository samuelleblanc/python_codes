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

def plt_amsr():
    from mpl_toolkits.basemap import Basemap
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

# <headingcell level=2>

# Load C130 nav data

# <codecell>

fnav = fp+'c130/20140919_ARISE_Flight_13/arise-C130-Hskping_c130_20140919_RA_Preliminary.ict'
nav = lm.load_ict(fnav)

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

flt = np.where((nav['Start_UTC']>19.0) & (nav['Start_UTC']<23.0) & (nav['Longitude']<0.0))

# <codecell>

plt.figure()
plt.plot(nav['Longitude'][flt],nav['GPS_Altitude'][flt])
plt.xlabel('Longitude')
plt.ylabel('Altitude [m]')
plt.savefig(fp+'plots/20140919_proile_alt.png',dpi=600,transparent=True)

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

fmodis_aqua = fp+'MODIS/'

