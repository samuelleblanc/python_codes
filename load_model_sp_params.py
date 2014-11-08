
# coding: utf-8

# In[28]:

import matplotlib.pyplot as plt
import numpy as np, h5py
import plotly.plotly as py
import scipy.io as sio
get_ipython().magic(u'matplotlib inline')
py.sign_in("samuelleblanc", "4y3khh7ld4")


# In[29]:

# set the basic directory path
fp='C:/Users/sleblan2/Research/TCAP/'


# In[51]:

# load the idl save file containing the modeled radiances
s=sio.idl.readsav(fp+'model/sp_v1_20130219_4STAR.out')
print s.keys()
print 'sp', s.sp.shape


# In[52]:

# create custom key for sorting via wavelength
iwvls = np.argsort(s.zenlambda)
s.zenlambda.sort()


# In[107]:

# load the matlab file containing the measured TCAP radiances
m = sio.loadmat(fp+'4STAR/20130219starzen_rad.mat')
sm = sio.idl.AttrDict(m)
print sm.keys()
print 'Measured radiance Shape: ', sm.rad.shape

print np.nanmax(sm.rad[sm.good[100],:])
sm.good[100]


# In[113]:

def nanmasked(x):
    " Build an array with nans masked out and the mask output"
    mask = ~np.isnan(x)
    maskA = x[mask]
    return (maskA,mask)

def closestindex(a,x):
    " Get the index from a of the closest value from x "
    return min(range(len(a)), key=lambda i: abs(a[i]-x))

def norm2max(x):
    " Returns a spectrum, x, that is normalized by its maximum value "
    return x/np.nanmax(x)
    
time_ref=17.22
ii = closestindex(sm.utc,time_ref)
rad,mask = nanmasked(sm.rad[sm.good[ii],:])


# In[121]:

# set up plotting of a few of the zenith radiance spectra
fig = plt.figure(1)
plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
#plt.plot(s.zenlambda,s.sp[0,iwvls,0,10,10])
plt.plot(sm.nm[mask],rad,lw=2, c='k', label='4STAR measured at: '+str(time_ref))
plt.title('Zenith spectra')
plt.ylabel('Radiance [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
plt.xlabel('Wavelength [nm]')
plt.xlim([350,1700])
plt.ylim([0,0.2])
plt.legend(frameon=False)
plt.show()


# In[102]:

#plot_url = py.plot_mpl(fig)


# In[126]:

# now go through and add the different modeled spectra
plt.figure(1)
plt.plot(s.zenlambda,s.sp[0,iwvls,0,10,10])


# In[ ]:



