
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

# from original first try at aerosol properties
aero.append({'z_arr':[2.0,5.0],
           'ext':np.array([0.6,0.4,0.10,0.04]),
           'ssa':np.array([0.8,0.85,0.9,0.95]),
           'asy':np.array([0.8,0.8,0.8,0.8]),
           'wvl_arr':[400.0,500.0,650.0,940.0],
           'disort_phase':False,
           'expand_hg':True,
           'label':'v1'})


# In[23]:

# from aeeronet at St helena on 2016-09-03
aero.append({'ssa':np.array([0.867100,0.844300,0.834100,0.832300]),
                  'asy':np.array([0.655889,0.621976,0.621509,0.627778]),
                  'wvl_arr':[439.0,675.0,870.0,1018.0],
             'label':'St_helena'})


# In[24]:

# from aeronet at Ascension on 2016-09-08
aero.append({'ssa':np.array([0.878400,0.859500,0.846800,0.838000]),
                  'asy':np.array([0.695644,0.654167,0.657314,0.669468]),
                  'wvl_arr':[441.0,675.0,871.0,1022.0],
            'label':'Ascension'})


# In[25]:

# from aeronet at Namibe on 2016-09-13
aero.append({'ssa':np.array([0.878400,0.859500,0.846800,0.838000]),
                  'asy':np.array([0.699605,0.624855,0.611445,0.631931]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
            'label':'Namibe'})


# In[26]:

# from aeronet at Lubango on 2016-09-14
aero.append({'ssa':np.array([0.867600,0.832200,0.795800,0.777500]),
                  'asy':np.array([0.672050,0.598835,0.562660,0.556078]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
               'label':'Lubango'})


# In[28]:

# saved merge properties from the above
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

# In[46]:

from load_utils import mat2py_time, toutc


# In[40]:

import scipy.io as sio


# In[42]:

s = sio.loadmat(fp+'data\\4STAR_20160914starsun.mat')


# In[43]:

s.keys()


# In[47]:

s['pyt'] = mat2py_time(s['t'])
s['utc'] = toutc(s['pyt'])


# In[48]:

s['tau_aero'].shape


# In[50]:

plt.figure()
plt.plot(s['utc'],s['tau_aero'][:,400])


# In[51]:

i = np.argmin(abs(s['utc']-14.1))


# In[52]:

i


# In[53]:

from Sp_parameters import smooth


# In[65]:

ss = np.nanmean(s['tau_aero'][i-6:i+6,:],axis=0)


# In[69]:

s['w'].shape


# In[80]:

plt.figure()
plt.plot(s['w'][0,:],s['tau_aero'][i,:]/3.0,'.')
plt.plot(s['w'][0,:],ss/3.0)
plt.ylim(0,0.2)
plt.xlim(0.35,1.7)
plt.grid()


# In[77]:

ext_new = []
for w in wvl_new:
    iw = np.argmin(abs(s['w']*1000.0-w))
    ext_new.append(ss[iw]/3.0)


# In[78]:

plt.figure()
plt.plot(wvl_new,ext_new,'x-')


# In[81]:

ext_new[-1] = 0.015
ext_new[5] = 0.035


# In[82]:

ext_new


# In[88]:

plt.figure()
plt.plot(wvl_new,ext_new,'x-')
plt.ylabel('Extinction coefficient [km$^{{-1}}$]')
plt.xlabel('Wavelenght [nm]')


# # Make the dict to be saved and used in next calculations of cloud properties

# In[83]:

aero_save = {'z_arr':[2.0,5.0],
        'ext':np.array([ext_new,ext_new]),
        'ssa':np.array([aero[-1]['ssa_new'],aero[-1]['ssa_new']]),
        'asy':np.array([aero[-1]['asy_new'],aero[-1]['asy_new']]),
        'wvl_arr':wvl_new,
        'disort_phase':False,
        'expand_hg':True}


# In[85]:

aero_save['ext'][1,:] = aero_save['ext'][1,:]*0.0


# In[86]:

aero_save


# ## Save the dict

# In[92]:

import json


# In[89]:

import pickle


# In[91]:

with open(fp+'aero_file.txt', 'wb') as handle:
  pickle.dump(aero_save, handle)


# In[95]:

json.dump(aero_save, open(fp+"aero_save.txt",'w'))


# In[ ]:



