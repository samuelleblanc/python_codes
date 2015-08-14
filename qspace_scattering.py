# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     qspace_scattering
# 
# Purpose:  
# 
#     Create a scattering phase function but from qspace
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
#     - 
#   
# Needed Files:
# 
#   - baum_2011_ice.out
#     
#  Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, CA, 2015-08-13
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
import load_modis as lm
IPython.InteractiveShell.cache_size = 0
# set the basic directory path
fp='C:/Users/sleblan2/Research/libradtran/ice/'

# <codecell>

baum = sio.idl.readsav(fp+'baum_2011_ice.out')

# <codecell>

baum.keys()

# <codecell>

print baum.phase.shape
print baum.ref.shape
print baum.wvl.shape
print baum.theta.shape

# <codecell>

baum.ref

# <codecell>

baum.wvl[10]

# <codecell>

baum.k = 2.0*np.pi/baum.wvl

# <codecell>

baum.qspace = 2.0*np.sin(baum.theta*np.pi/180.0/2.0)

# <codecell>

baum.qspace.shape

# <codecell>

baum.k.shape

# <codecell>

for i,k in enumerate(baum.k):
    baum.qspace[i,:,:,:] = baum.qspace[i,:,:,:]*baum.k[i]

# <codecell>

import cmaps
cmaps.cmaps()

# <codecell>

for i,r in enumerate(baum.ref):
    if (i%2 == 0):
        plt.plot(baum.qspace[10,i,0,:],baum.phase[10,i,0,:],label='r$_{e}$=%2i$\\mu$m'%r,color=plt.cm.gist_ncar(float(i)/len(baum.ref)))
plt.yscale('log')
plt.xscale('log')
plt.legend(frameon=False)
plt.title('Ice crystal General Habit mixture - 500 nm')
plt.xlabel('qspace')
plt.ylabel('P$_{11}$')
plt.savefig(fp+'baum_2011_qspace_500nm.png',dpi=600,transparent=True)

# <codecell>


