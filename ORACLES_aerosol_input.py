
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_aerosol_input
# 
# Purpose:  
# 
#     Compare the input of aerosol properties for cloud retrievals. For perparing the aerosol file in Preperation_ORACLES
#   
# Input:
# 
#     none
# 
# Output:
#    
#     plots
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - matplotlib
#     - scipy
# 
#   
# Needed Files:
# 
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc,Santa Cruz, CA, 2016-11-04
#     

# # Prepare the python environment

# In[1]:

import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[2]:

get_ipython().magic(u'matplotlib notebook')


# In[3]:

fp = 'C:\\Users\\sleblan2\\Research\\ORACLES\\'


# # Get the possible aerosol values

# In[21]:

aero = []


# In[22]:

aero.append({'z_arr':[2.0,5.0],
           'ext':np.array([0.6,0.4,0.10,0.04]),
           'ssa':np.array([0.8,0.85,0.9,0.95]),
           'asy':np.array([0.8,0.8,0.8,0.8]),
           'wvl_arr':[400.0,500.0,650.0,940.0],
           'disort_phase':False,
           'expand_hg':True,
           'label':'v1'})


# In[23]:

aero.append({'ssa':np.array([0.867100,0.844300,0.834100,0.832300]),
                  'asy':np.array([0.655889,0.621976,0.621509,0.627778]),
                  'wvl_arr':[439.0,675.0,870.0,1018.0],
             'label':'St_helena'})


# In[24]:

aero.append({'ssa':np.array([0.878400,0.859500,0.846800,0.838000]),
                  'asy':np.array([0.695644,0.654167,0.657314,0.669468]),
                  'wvl_arr':[441.0,675.0,871.0,1022.0],
            'label':'Ascension'})


# In[25]:

aero.append({'ssa':np.array([0.878400,0.859500,0.846800,0.838000]),
                  'asy':np.array([0.699605,0.624855,0.611445,0.631931]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
            'label':'Namibe'})


# In[26]:

aero.append({'ssa':np.array([0.867600,0.832200,0.795800,0.777500]),
                  'asy':np.array([0.672050,0.598835,0.562660,0.556078]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
               'label':'Lubango'})


# In[28]:

aero.append({'ssa':np.array([0.867100,0.844300,0.834100,0.832300]),
                  'asy':np.array([0.699605,0.624855,0.611445,0.631931]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
               'label':'v3'})


# ## Plot the aerosol properties

# In[29]:

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212,sharex=ax1)
for d in aero:
    ax1.plot(d['wvl_arr'],d['ssa'],'x-',label=d['label'])
    ax2.plot(d['wvl_arr'],d['asy'],'o-',label=d['label'])
ax2.set_xlabel('Wavelength [nm]')
ax1.set_ylabel('SSA')
ax2.set_ylabel('ASY')
ax1.legend(frameon=False,numpoints=1)


# ## Make a spline interpolation to extend the properties to longer and horter wavelengths

# In[30]:

from scipy import interpolate


# In[37]:

wvl_new = [350.0,400.0,500.0,650.0,875.0,980.0,1020.0,1240.0,1710.0]


# In[39]:

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212,sharex=ax1)
for d in aero:
    l = ax1.plot(d['wvl_arr'],d['ssa'],'x',label=d['label'])
    f_new = interpolate.InterpolatedUnivariateSpline(d['wvl_arr'],d['ssa'],k=3,ext=0)
    ax1.plot(wvl_new,f_new(wvl_new),'-',color=l[0].get_color())
    d['wvl_new'] = wvl_new
    d['ssa_new'] = f_new(wvl_new)
    l2 = ax2.plot(d['wvl_arr'],d['asy'],'o',label=d['label'])
    f2_new = interpolate.InterpolatedUnivariateSpline(d['wvl_arr'],d['asy'],k=3,ext=0)
    ax2.plot(wvl_new,f2_new(wvl_new),'-',color=l2[0].get_color())
    d['asy_new'] = f2_new(wvl_new)
ax2.set_xlabel('Wavelength [nm]')
ax1.set_ylabel('SSA')
ax2.set_ylabel('ASY')
ax1.legend(frameon=False,numpoints=1)


# ## set the new values

# In[45]:

aero[-1]['ssa_new'] = aero[3]['ssa_new']
aero[-1]['asy_new'] = aero[1]['asy_new']


# # Get the extinction of the layer in question

# In[40]:

import scipy.io as sio


# In[42]:

s = sio.loadmat(fp+'data\\4STAR_20160914starsun.mat')


# In[43]:

s.keys()


# In[ ]:



