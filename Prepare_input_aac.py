# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     Prepare_input_aac
# 
# Purpose:  
# 
#     Created input files for libradtran based on the calipso, global Aerosol over Clouds
#     For Meloë's study
# 
# Calling Sequence:
# 
#     python Prepare_input_aac
#   
# Input:
# 
#     none
# 
# Output:
#    
#     input files for libradtran 2.0 (uvspec) 
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - mplt_toolkits for basemap, map plotting
#     - pdb
#     - datetime
# 
#   
# Needed Files:
# 
#   - matlab input files: Input_to_DARF_mmm.mat

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap

# <codecell>

import Run_libradtran as RL

# <codecell>

fp = 'C:\Users\sleblan2\Research\Calipso\meloe/'

# <headingcell level=2>

# Read Matlab input files

# <codecell>

input_DJF = sio.loadmat(fp+'Input_to_DARF_DJF.mat',mat_dtype=True)['data_input_darf']

# <codecell>

input_DJF.dtype.names

# <codecell>

input_DJF['MODIS_lat'][0,0].shape

# <codecell>

input_DJF['MODIS_COD_mean'][0,0][input_DJF['MODIS_COD_mean'][0,0]==-9999] = np.nan

# <codecell>

ctr = plt.contourf(input_DJF['MODIS_lon'][0,0][:,0],input_DJF['MODIS_lat'][0,0][:,0],input_DJF['MODIS_COD_mean'][0,0])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(ctr)
cbar.set_label('MODIS MEAN COD')
plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=11,c='k')
plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=10,c=input_DJF['MODIS_COD_mean'][0,0])
plt.xlim([-180,180])
plt.ylim([-90,90])
plt.savefig(fp+'plots/MODIS_mean_COD_flat.png',dpi=600,transparent=True)

# <codecell>

MODIS_lon[20,72]

# <codecell>

input_DJF['MODIS_COD_mean'][0,0][:,72]

# <codecell>

m = Basemap(projection='moll',lon_0=0,resolution='c')
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,20.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,420.,60.))
MODIS_lon, MODIS_lat = np.meshgrid(input_DJF['MODIS_lon'][0,0][:,0],input_DJF['MODIS_lat'][0,0][:,0],sparse=False)
x,y=m(MODIS_lon, MODIS_lat) 
cs = m.contourf(x,y,input_DJF['MODIS_COD_mean'][0,0],np.linspace(0,60,21),extend='max')
cbar = plt.colorbar(cs)
cbar.set_label('MODIS mean COD')
plt.savefig(fp+'plots/MODIS_mean_COD.png',dpi=600,transparent=True)

# <codecell>

input_DJF['MODIS_COD_mean'][0,0][20,20]

# <codecell>

input_DJF['MODIS_effrad_mean'][0,0][20,20]

# <codecell>

input_DJF['MOC_wavelengths'][0,0]

# <codecell>

print MODIS_lon[20,20],MODIS_lat[20,20]

# <codecell>

input_DJF['MOC_asym_mean'][0,0][20,20]

# <codecell>

tt = np.array([1,2,3,4,5.0])

# <codecell>

type(tt)

# <codecell>

to = [1,2,3,4,5.0]

# <codecell>

type(to)

# <codecell>

len(tt)

# <codecell>

len(to)

# <codecell>

ext = [[1,2,3],[4.0,5.0,6.0]]

# <codecell>

len(ext)

# <codecell>

ee = np.array(ext)

# <codecell>

len(ee)

# <codecell>

size(ext)

# <codecell>

size(ee)

# <codecell>

shape(ee)

# <codecell>

shape(ext)

# <codecell>

tr = [  4.85364736e-02,   4.66139195e-02,   4.47312609e-02,
         4.30849589e-02,   4.19923201e-02,   4.03355801e-02,
         3.74764159e-02,   3.45595009e-02,   3.19684762e-02,
         2.94772306e-02,   2.74202103e-02,   2.57334360e-02,
         2.39507641e-02,   2.08731768e-02,   1.67933569e-02,
         1.29016393e-02,   9.04034361e-03,   6.65431703e-03,
         4.35758656e-03,   3.47793084e-03,   2.59552084e-03,
         1.92045503e-03,   1.48977972e-03,   1.14460091e-03,
         7.92407241e-04,   5.10274383e-04,   3.17954425e-04,
         1.69683997e-04,   8.10304392e-05,   3.35441191e-05]

# <codecell>

shape(tr)

# <codecell>


