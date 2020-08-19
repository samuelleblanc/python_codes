#!/usr/bin/env python
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


# In[19]:


import plotting_utils as pu
import hdf5storage as hs
import os
import scipy.io as sio
from path_utils import getpath


# In[3]:


fp = 'C:\\Users\\sleblan2\\Research\\ORACLES\\'


# In[4]:


name = 'ORACLES'
vv = 'v3'
vr = 'R3'


# In[5]:


fp = getpath(name)
fp_rtm = getpath('rtm')


# # Get the possible aerosol values

# In[4]:


aero = []


# In[5]:


# from original first try at aerosol properties
aero.append({'z_arr':[2.0,5.0],
           'ext':np.array([0.6,0.4,0.10,0.04]),
           'ssa':np.array([0.8,0.85,0.9,0.95]),
           'asy':np.array([0.8,0.8,0.8,0.8]),
           'wvl_arr':[400.0,500.0,650.0,940.0],
           'disort_phase':False,
           'expand_hg':True,
           'label':'v1'})


# In[6]:


# from aeeronet at St helena on 2016-09-03
aero.append({'ssa':np.array([0.867100,0.844300,0.834100,0.832300]),
                  'asy':np.array([0.655889,0.621976,0.621509,0.627778]),
                  'wvl_arr':[439.0,675.0,870.0,1018.0],
             'label':'St_helena'})


# In[7]:


# from aeronet at Ascension on 2016-09-08
aero.append({'ssa':np.array([0.878400,0.859500,0.846800,0.838000]),
                  'asy':np.array([0.695644,0.654167,0.657314,0.669468]),
                  'wvl_arr':[441.0,675.0,871.0,1022.0],
            'label':'Ascension'})


# In[8]:


# from aeronet at Namibe on 2016-09-13
aero.append({'ssa':np.array([0.878400,0.859500,0.846800,0.838000]),
                  'asy':np.array([0.699605,0.624855,0.611445,0.631931]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
            'label':'Namibe'})


# In[9]:


# from aeronet at Lubango on 2016-09-14
aero.append({'ssa':np.array([0.867600,0.832200,0.795800,0.777500]),
                  'asy':np.array([0.672050,0.598835,0.562660,0.556078]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
               'label':'Lubango'})


# In[10]:


# saved merge properties from the above
aero.append({'ssa':np.array([0.867100,0.844300,0.834100,0.832300]),
                  'asy':np.array([0.699605,0.624855,0.611445,0.631931]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
               'label':'v3'})


# ## Plot the aerosol properties

# In[11]:


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


# In[12]:


fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212,sharex=ax1)
for d in aero:
    if d['label'] in ['v1','v3']:
        break
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

# In[126]:


from load_utils import load_from_json, save_to_json


# In[122]:


save_to_json(fp+'aero_save.txt',aero_save)


# In[123]:


a = ()


# In[127]:


a = load_from_json(fp+'aero_save.txt')


# # Load all the skyscan retirevals for 2016

# In[6]:


try:
    ae,ae_dict = lu.load_netcdf(fp+'aeroinv_2016/netcdf4/4STAR-aeroinv_P3_2016_R0.nc',everything=True)
except:
    import h5py as h5
    f5 = h5.File(fp+'aeroinv_2016/netcdf4/4STAR-aeroinv_P3_2016_R0.nc','r')
    ae5 = {}
    ae5_dict = {}
    for ka,kd in f5.iteritems():
        ae5[ka] = kd.value
        ae5_dict[ka] = {}
        for kdict in kd.attrs.iteritems():
            if type(kdict[1])!=type(np.array([])):
                ae5_dict[ka][kdict[0]] = kdict[1]
    ae = ae5
    ae_dict = ae5_dict


# In[15]:


ka = ae_dict.keys()
ka.sort()
ka


# ## Plot out the g

# In[106]:


plt.figure()
plt.plot(ae['wavelength'],ae['g_fine'].T,'o--r')
plt.plot(ae['wavelength'],ae['g_coarse'].T,'d:b')
plt.plot(ae['wavelength'],ae['g_total'].T,'+-k')

plt.plot(ae['wavelength'],ae['g_fine'].T[:,0],'o--r',label='fine')
plt.plot(ae['wavelength'],ae['g_coarse'].T[:,0],'d:b',label='coarse')
plt.plot(ae['wavelength'],ae['g_total'].T[:,0],'+-k',label='total')



plt.xlabel('Wavelength [nm]')
plt.ylabel('Asymmetry Parameter')
plt.title('ORACLES 2016 4STAR skyscan retrievals')
plt.legend(frameon=False)


plt.savefig(fp+'plot/ORACLES2016_asym.png',dpi=600,transparent=True)


# ## Plot out the refractive index

# In[17]:


fig,ax = plt.subplots(2,1)

ax[0].plot(ae['wavelength'],ae['n_real'].T,'.--')
ax[1].plot(ae['wavelength'],ae['n_imag'].T,'.--')

ax[1].set_xlabel('Wavelength [nm]')
ax[0].set_ylabel('Real Index of ref.')
ax[1].set_ylabel('Imag. Index of ref.')
ax[0].set_title('ORACLES 2016 skyscan retrievals')


# In[105]:


fig,ax = plt.subplots(2,1)

ax[0].boxplot(ae['n_real'],positions=ae['wavelength'],widths=20.0,showmeans=True)
ax[1].boxplot(ae['n_imag'],positions=ae['wavelength'],widths=20.0,showmeans=True)

ax[0].set_xlim(380,1020)
ax[1].set_xlim(380,1020)

ax[1].set_xlabel('Wavelength [nm]')
ax[0].set_ylabel('Real Index of ref.')
ax[1].set_ylabel('Imag. Index of ref.')
ax[0].set_title('ORACLES 2016 skyscan retrievals')

plt.savefig(fp+'plot/ORACLES2016_refractive_index.png',dpi=600,transparent=True)


# ## Plot out the size distribution

# In[40]:


plt.figure()
plt.plot(ae['radius'],ae['psd'].T)
plt.xscale('log')

plt.xlabel('Radius [$\mu$m]')
plt.ylabel('Volume density [dV/dlnR]')


# In[41]:


ae['psd'].shape


# In[43]:


ae['radius'].shape


# In[44]:


ae['radius']


# In[49]:


ae_fine_mode = [ae['radius'][np.argmax(p)] for p in ae['psd'][:,0:9]]
ae_coarse_mode = [ae['radius'][9:][np.argmax(p)] for p in ae['psd'][:,9:]]


# In[53]:


plt.figure()
plt.boxplot([ae_fine_mode,ae_coarse_mode],positions=[0,0.1],vert=False)


# In[107]:


fig = plt.figure()
#
bp = plt.boxplot([ae_fine_mode,ae_coarse_mode],positions=[0.15,0.15],vert=False,showmeans=True,widths=0.02)
plt.plot(ae['radius'],ae['psd'].T)
plt.xscale('log')

[plt.axvline(m.get_data()[0],color='g',ls=':') for m in bp['means']]
[plt.text(m.get_data()[0]+0.01,0.165,'{:1.3f}'.format(m.get_data()[0][0])) for m in bp['means']]

plt.xlabel('Radius [$\mu$m]')
plt.ylabel('Volume density [dV/dlnR]')
plt.ylim(0,0.18)
plt.yticks([0.0,0.05,0.1,0.15])
plt.gca().set_yticklabels([0.0,0.05,0.1,0.15])

plt.title('ORACLES 2016 4STAR retrievals')

plt.savefig(fp+'plot/ORACLES2016_size_dist.png',dpi=600,transparent=True)

