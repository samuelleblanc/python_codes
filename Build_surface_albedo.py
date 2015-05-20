# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     Build_surface_albedo
# 
# Purpose:  
# 
#     Python script that builds surface albedo from SSFR measurements
#     saves in a format that is useful for libradtran inputs
# 
# Calling Sequence:
# 
#     python Build_surface_albedo.py
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures and text files
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - os
#     - pdb : for debugging
#     - datetime
#     - mpl_toolkits
#     - plotting_utils (user defined plotting routines)
#   
# Needed Files:
# 
#   - SSFR calibspcs .out file

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

# Load SSFR file

# <codecell>

daystr = '20140919'

# <codecell>

ssfr = sio.idl.readsav(fp+'c130/'+daystr+'_ARISE_Flight_13/'+daystr+'_calibspcs_attcorr.out')

# <codecell>

plt.plot(ssfr.tmhrs,ssfr.zspectra[:,50],label='Downwelling')
plt.plot(ssfr.tmhrs,ssfr.nspectra[:,50],label='Upwelling')
plt.ylabel('Irradiance [W/m$^{2}$]')
plt.xlabel('UTC [h]')
plt.title('Irradiance on '+daystr+' at %4i nm' % ssfr.zenlambda[50])
plt.legend(frameon=False,loc=2)

# <markdowncell>

# Loop and interpolate nadir spectra to zenith wavelengths

# <codecell>

from scipy.interpolate import interp1d
fn = interp1d(ssfr.nadlambda,ssfr.nspectra,axis=1,bounds_error=False)
ssfr.nnspectra = fn(ssfr.zenlambda)

# <codecell>

wvl = np.arange(350,1701)

# <codecell>

ssfr.nspectra1 = fn(wvl)

# <codecell>

fz = interp1d(ssfr.zenlambda,ssfr.zspectra,axis=1,bounds_error=False)
ssfr.zspectra1 = fz(wvl)

# <headingcell level=2>

# Find the correct times for plotting and outputting

# <codecell>

import Sp_parameters as Sp

# <codecell>

iice = np.argmin(abs(ssfr.tmhrs-21.5828))
iwat = np.argmin(abs(ssfr.tmhrs-20.9261))

# <codecell>

alb = ssfr.nspectra1/ssfr.zspectra1
alb[alb<=0.0] = 0.0
alb[alb>=1.0] = 1.0
alb[np.isnan(alb)] = 0.0

# <codecell>

plt.plot(wvl,alb[iwat,:],'b.')
plt.plot(wvl,Sp.smooth(alb[iwat,:],5),'r')

plt.xlim([350,1700])
plt.ylim([0,1])
plt.ylabel('Albedo')
plt.xlabel('Wavelength [nm]')
plt.title('Surface albedo, above water UTC: %.4f' % ssfr.tmhrs[iwat])
plt.savefig(fp+'c130/'+daystr+'_ARISE_Flight_13/'+daystr+'_surface_albedo_water.png',dpi=600,transparent=True)

# <codecell>

plt.plot(wvl,alb[iice,:],'b.')
plt.plot(wvl,Sp.smooth(alb[iice,:],5),'r')
plt.xlim([350,1700])
plt.ylim([0,1])
plt.ylabel('Albedo')
plt.xlabel('Wavelength [nm]')
plt.title('Surface albedo, above ice UTC: %.4f' % ssfr.tmhrs[iice])
plt.savefig(fp+'c130/'+daystr+'_ARISE_Flight_13/'+daystr+'_surface_albedo_ice.png',dpi=600,transparent=True)

# <codecell>

alb_wat = Sp.smooth(alb[iwat,:],5)
alb_ice = Sp.smooth(alb[iice,:],5)

# <headingcell level=2>

# Save to file

# <codecell>

f = open(fp+daystr+'_surf_alb_water.dat','w')
for i in xrange(len(wvl)):
    f.write('%i %.4f \n' % (wvl[i],alb_wat[i]))
f.close()

# <codecell>

f = open(fp+daystr+'_surf_alb_ice.dat','w')
for i in xrange(len(wvl)):
    f.write('%i %.4f \n' % (wvl[i],alb_ice[i]))
f.close()

