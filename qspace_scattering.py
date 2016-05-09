
# coding: utf-8

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

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
get_ipython().magic(u'matplotlib inline')
import numpy as np, h5py
import scipy.io as sio
import scipy
import math, os, IPython
import Sp_parameters as Sp
import load_utils as lm
IPython.InteractiveShell.cache_size = 0
# set the basic directory path
fp='C:/Users/sleblan2/Research/libradtran/ice/'


# In[57]:

baum = sio.idl.readsav(fp+'baum_2011_ice.out')


# In[58]:

baum.keys()


# In[59]:

print baum.phase.shape
print baum.ref.shape
print baum.wvl.shape
print baum.theta.shape


# In[60]:

baum.ref


# In[61]:

baum.wvl[10]


# In[62]:

baum.k = 2.0*np.pi/baum.wvl


# In[63]:

baum.qspace = 2.0*np.sin(baum.theta*np.pi/180.0/2.0)


# In[64]:

baum.qspace.shape


# In[65]:

baum.k.shape


# In[66]:

for i,k in enumerate(baum.k):
    baum.qspace[i,:,:,:] = baum.qspace[i,:,:,:]*baum.k[i]


# In[36]:

import cmaps
cmaps.cmaps()


# In[67]:

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


# In[47]:

import Sp_parameters as sp


# In[52]:

baum.dq = baum.qspace
for i,r in enumerate(baum.ref):
    for j,w in enumerate(baum.wvl):
        baum.dq[j,i,0,:] = sp.deriv(np.log(baum.phase[j,i,0,:]),np.log(baum.qspace[j,i,0,:]))


# In[53]:

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


# In[68]:

blogq = np.log(baum.qspace[10,10,0,:])
blogs = np.log(baum.phase[10,10,0,:])


# In[69]:

plt.plot(baum.qspace[10,10,0,:],baum.phase[10,10,0,:])


# In[70]:

plt.plot(blogq,blogs)


# In[85]:

towwrite = np.array(([blogq],[blogs]))[:,0,:]


# In[97]:

np.savetxt(fp+'baum_2011_qspace.dat',np.flipud(towwrite.T))


# In[96]:

np.flipud(towwrite.T)


# In[92]:

np.flipud(towwrite).T.shape


# ## Do legendre expansion of the qspace log log

# In[102]:

qleg = np.polynomial.legendre.legfit(baum.qspace[10,10,0,:],baum.phase[10,10,0,:],50)


# In[103]:

qleg.shape


# In[104]:

qleg


# In[108]:

ll = np.polynomial.legendre.Legendre(qleg)


# In[111]:

x,y = ll.linspace()


# In[113]:

plt.plot(x,y)
#plt.yscale('log')
#plt.xscale('log')


# In[ ]:



