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

import Sp_parameters as sp

# <codecell>

baum.dq = baum.qspace
for i,r in enumerate(baum.ref):
    for j,w in enumerate(baum.wvl):
        baum.dq[j,i,0,:] = sp.deriv(np.log(baum.phase[j,i,0,:]),np.log(baum.qspace[j,i,0,:]))

# <codecell>

for i,r in enumerate(baum.ref):
    if (i%2 == 0):
        plt.plot(baum.dq[10,i,0,:],baum.phase[10,i,0,:],label='r$_{e}$=%2i$\\mu$m'%r,color=plt.cm.gist_ncar(float(i)/len(baum.ref)))
#plt.yscale('log')
#plt.xscale('log')
plt.legend(frameon=False)
plt.title('Ice crystal GHM derivative - 500 nm')
plt.xlabel('qspace')
plt.ylabel('dP$_{11}$')
plt.savefig(fp+'baum_2011_deriv_qspace_500nm.png',dpi=600,transparent=True)

# <codecell>

blogq = np.log(baum.qspace[10,10,0,:])
blogs = np.log(baum.phase[10,10,0,:])

# <codecell>

plt.plot(baum.qspace[10,10,0,:],baum.phase[10,10,0,:])

# <codecell>

plt.plot(blogq,blogs)

# <codecell>

towwrite = np.array(([blogq],[blogs]))[:,0,:]

# <codecell>

np.savetxt(fp+'baum_2011_qspace.dat',np.flipud(towwrite.T))

# <codecell>

np.flipud(towwrite.T)

# <codecell>

np.flipud(towwrite).T.shape

# <headingcell level=2>

# Do legendre expansion of the qspace log log

# <codecell>

qleg = np.polynomial.legendre.legfit(baum.qspace[10,10,0,:],baum.phase[10,10,0,:],50)

# <codecell>

qleg.shape

# <codecell>

qleg

# <codecell>

ll = np.polynomial.legendre.Legendre(qleg)

# <codecell>

x,y = ll.linspace()

# <codecell>

plt.plot(x,y)
#plt.yscale('log')
#plt.xscale('log')

# <codecell>


