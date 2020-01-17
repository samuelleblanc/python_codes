#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_skyscan_aerosol_prop
# 
# Purpose:  
# 
#     ORACLES skyscan results and compared to local aeronet values
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
#     Written: Samuel LeBlanc,Santa Cruz, CA, 2016-12-08
#     

# # Prepare the python environment

# In[1]:


import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[2]:


import hdf5storage as hs


# In[15]:


import scipy.io as sio


# In[44]:


import plotting_utils as pu


# In[45]:


from map_interactive import build_basemap


# In[3]:


get_ipython().magic(u'matplotlib notebook')


# In[4]:


fp = 'C:\\Users\\sleblan2\\Research\\ORACLES\\'


# # Get the possible aerosol values

# In[117]:


aero = []


# In[118]:


# from aeeronet at St helena on 2016-09-03
aero.append({'ssa':np.array([0.867100,0.844300,0.834100,0.832300]),
                  'asy':np.array([0.655889,0.621976,0.621509,0.627778]),
                  'wvl_arr':[439.0,675.0,870.0,1018.0],
             'label':'St-Helena\n2016-09-03',
             'lat':-15.94194444444444,
             'lon':-5.66694444444444})


# In[119]:


# from aeronet at Ascension on 2016-09-08
aero.append({'ssa':np.array([0.878400,0.859500,0.846800,0.838000]),
             'asy':np.array([0.695644,0.654167,0.657314,0.669468]),
             'wvl_arr':[441.0,675.0,871.0,1022.0],
            'label':'Ascension\n2016-09-08',
             'lat':-7.96694444444444,
             'lon':-14.35})


# In[120]:


# from aeronet at Namibe on 2016-09-13
aero.append({'ssa':np.array([0.853800,0.849600,0.844100,0.840000]),
                  'asy':np.array([0.699605,0.624855,0.611445,0.631931]),
                  'wvl_arr':[440.0,675.0,869.0,1018.0],
            'label':'Namibe\n2016-09-13',
             'lat':-15.158889,
             'lon':12.177778})


# In[121]:


# from aeronet at Lubango on 2016-09-14
aero.append({'ssa':np.array([0.867600,0.832200,0.795800,0.777500]),
             'asy':np.array([0.672050,0.598835,0.562660,0.556078]),
             'wvl_arr':[440.0,675.0,869.0,1018.0],
             'label':'Lubango\n2016-09-14',
             'lat':-14.9577778,
             'lon':13.445})


# ## Read in the matlab retrieved file

# In[12]:


fp


# In[14]:


fp+'skyscan\\4STAR_.4STAR_20160920_033_STARSKYP.created_20161125_005553.mat'


# In[16]:


m = sio.loadmat(fp+'skyscan\\4STAR_.4STAR_20160920_033_STARSKYP.created_20161125_005553.mat')


# In[82]:


m2 = sio.loadmat(fp+'skyscan\\4STAR_.4STAR_20160920_030_STARSKYA.created_20161125_004250.mat')


# In[17]:


m.keys()


# In[19]:


m['Wavelength']


# In[20]:


m['lat'] = -16.5866
m['lon'] = 10.1636


# In[83]:


m2['lat'] = -16.4704
m2['lon'] = 10.5000


# ## Plot the aerosol properties

# In[24]:


m['ssa_total'].shape


# In[25]:


m['Wavelength'].shape


# In[86]:


import matplotlib
matplotlib.rcParams.update({'font.size': 18})


# In[128]:


fig = plt.figure(figsize=(9,6))
ax1 = fig.add_subplot(111)
#ax2 = fig.add_subplot(122,sharex=ax1)
for d in aero:
    ax1.plot(d['wvl_arr'],d['ssa'],'+-',label=d['label'])
    #ax2.plot(d['wvl_arr'],d['asy'],'o-',label=d['label'])
ax1.plot(m['Wavelength'][:,0]*1000.0,m['ssa_total'][0,:],'s-k',zorder=100,label='4STAR\n2016-09-20')
#ax1.plot(m2['Wavelength'][:,0]*1000.0,m2['ssa_total'][0,:],'d-',color='grey',zorder=100,label='4STAR #2')
ax1.set_xlabel('Wavelength [nm]')
ax1.set_ylabel('SSA')
#ax2.set_ylabel('ASY')
pu.prelim()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.7, box.height])

ax1.legend(numpoints=1,bbox_to_anchor=[1.55,1.0])    
plt.savefig(fp+'plot//SSA_4STAR_retrieval.png',transprent=True,dpi=600)


# In[129]:


fig = plt.figure()
ax = fig.add_subplot(111)
mb = build_basemap(lower_left=[-15,-25],upper_right=[18,10],ax=ax,larger=False)

cl = ['blue','green','red','cyan']

for i,d in enumerate(aero):
    xlo,yla = mb(d['lon'],d['lat'])
    ax.text(xlo,yla,'.')
    ax.text(xlo,yla,d['label'],horizontalalignment='left',verticalalignment='bottom',color=cl[i])
    #print d['label'], d['lon'],d['lat']

xlo,yla = mb(m['lon'],m['lat'])
ax.text(xlo,yla,'.')
ax.text(xlo,yla,'4STAR',horizontalalignment='left',verticalalignment='bottom')

xlo2,yla2 = mb(m2['lon'],m2['lat'])
#ax.text(xlo2,yla2,'*')
#ax.text(xlo2,yla2,'4STAR #2',color='grey',horizontalalignment='left',verticalalignment='bottom')

plt.savefig(fp+'plot//SSA_map.png',transparent=True,dpi=600)


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

