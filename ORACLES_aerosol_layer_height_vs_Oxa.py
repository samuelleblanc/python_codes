#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Investigate the Oxygen-A changes when the Aerosol layer varies
#     Based on gap analysis from ORACLES
# 
# Input:
# 
#     None
# 
# Output:
# 
#     Figure and save files
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - load_utils.py
#     - matplotlib
#     - numpy
#     - Sp_parameters
#     - write_utils
#     - path_utils
#     - hdf5storage
#     - scipy
# 
# Needed Files:
#   - starsun.mat files from ORACLES 2016
#   - ORACLES2016_gap_v1.npy:  gap analysis from ORACLES 2016, from ORACLES_AOD_summary_allyears.ipynb
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-02-11
#     Modified:
# 

# # Prepare python environment

# In[2]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import os
from scipy import interpolate


# In[3]:


name = 'ORACLES'
vv = 'v1'
fp = getpath(name)


# In[4]:


fpo = '/data/sunsat/ORACLES_2016/data_processed/starsuns/R3/'


# # Load files

# ## Load gap files

# In[5]:


gap6 = np.load(fp+'ORACLES2016_gap_{}.npy'.format(vv),allow_pickle=True,fix_imports=True,encoding='latin1')
gap6 = gap6.item()


# In[13]:


gap7 = np.load(fp+'ORACLES2017_gap_{}.npy'.format(vv),allow_pickle=True,fix_imports=True,encoding='latin1')
gap7 = gap7.item()


# In[14]:


gap8 = np.load(fp+'ORACLES2018_gap_{}.npy'.format(vv),allow_pickle=True,fix_imports=True,encoding='latin1')
gap8 = gap8.item()


# In[15]:


gap6.keys()


# In[9]:


days6 = [d.strftime('%Y%m%d') for d in gap6['aero_base_day']]
day6 = np.unique(days6)


# In[10]:


day6


# In[16]:


days7 = [d.strftime('%Y%m%d') for d in gap7['aero_base_day']]
day7 = np.unique(days7)


# In[19]:


days8 = [d.strftime('%Y%m%d') for d in gap8['aero_base_day']]
day8 = np.unique(days8)


# ## Load starsun files

# ### Load 2016

# In[27]:


s6 = []
for d in day6:
    ss = hs.loadmat(fpo+'4STAR_{}starsun.mat'.format(d),variable_names=['tau_aero','w','t','Pst'])
    ss['filename'] = fpo+'4STAR_{}starsun.mat'.format(d)
    s6.append(ss)


# In[28]:


len(s6)


# In[29]:


s6[0].keys()


# In[30]:


for ss in s6:
    ss['utc'] = lu.mat2py_time(ss['t'])


# ### Load 2017

# In[20]:


fpo7 = '/data/sunsat/ORACLES_2017/data_processed/starsuns/R0/'


# In[22]:


s7 = []
for d in day7:
    ss = hs.loadmat(fpo7+'4STAR_{}starsun.mat'.format(d),variable_names=['tau_aero','w','t','Pst'])
    ss['filename'] = fpo7+'4STAR_{}starsun.mat'.format(d)
    print(ss['filename'],'/',len(day7))
    s7.append(ss)


# In[31]:


for ss in s7:
    ss['utc'] = lu.mat2py_time(ss['t'])


# ### Load 2018

# In[23]:


fpo8 = '/data/sunsat/ORACLES_2018/data_processed/starsuns/R0/'


# In[24]:


s8 = []
for d in day8:
    ss = hs.loadmat(fpo8+'4STAR_{}starsun.mat'.format(d),variable_names=['tau_aero','w','t','Pst'])
    ss['filename'] = fpo8+'4STAR_{}starsun.mat'.format(d)
    print(ss['filename'])
    s8.append(ss)


# In[32]:


for ss in s8:
    ss['utc'] = lu.mat2py_time(ss['t'])


# ## Calculate the Ox-a depth

# In[33]:


def oxa_depth(w,spec):
    'calculate the oxygen-a depth. w: wavelenght in microns, spec: spectra(tau_aero)'
    y0 = np.argmin(abs(w-0.756))
    y1 = np.argmin(abs(w-0.772))
    oxa_flat,oxa_delta = [],[]
    for i in range(len(spec)):
        fx = interpolate.interp1d(w[0,[y0,y1]],spec[i,[y0,y1]])
        oxa_flat.append(fx(w[0,y0:y1]))
        oxa_delta.append(np.nansum(spec[i,y0:y1]-oxa_flat))
    return np.array(oxa_delta)


# ### 2016 Ox-A

# In[34]:


for j,ss in enumerate(s6):
    print('Doing file {}, which is {}/{}'.format(day6[j],j,len(s6)))
    if 'oxa_delta' in ss.keys(): 
        print('Yas queen already in there')
        continue
    try:
        ss['oxa_delta'] = oxa_depth(ss['w'],ss['tau_aero'])
    except:
        pass


# ### 2017 Ox-A

# In[35]:


for j,ss in enumerate(s7):
    print('Doing file {}, which is {}/{}'.format(day7[j],j,len(s7)))
    if 'oxa_delta' in ss.keys(): 
        print('Yas queen already in there')
        continue
    try:
        ss['oxa_delta'] = oxa_depth(ss['w'],ss['tau_aero'])
    except:
        pass


# ### 2018 Ox-A

# In[36]:


for j,ss in enumerate(s8):
    print('Doing file {}, which is {}/{}'.format(day8[j],j,len(s8)))
    if 'oxa_delta' in ss.keys(): 
        print('Yas queen already in there')
        continue
    try:
        ss['oxa_delta'] = oxa_depth(ss['w'],ss['tau_aero'])
    except:
        pass


# ## Get the oxa depth for the gap times

# In[37]:


dd = gap6['aero_base_day'][1]


# In[38]:


utc_fx = lambda x: x.hour+(x.minute+(x.second/60.0))/60.0


# In[39]:


s6[iu].keys()


# In[40]:


len(gap6['aero_base_alt'])


# ### 2016

# In[45]:


n = len(gap6['aero_base_day'])
oxa_aero_base,oxa_meas,press_meas,aod = np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan 
for i,dd in enumerate(gap6['aero_base_day']):
    iu = np.where(day6==days6[i])[0][0]
    it = np.argmin(abs(dd-s6[iu]['utc']))
    try:
        oxa_aero_base[i] = s6[iu]['oxa_delta'][it]
    except KeyError:
        continue
    it2 = np.argmin(abs(gap6['meas_low_day'][i]-s6[iu]['utc']))
    oxa_meas[i] = s6[iu]['oxa_delta'][it2]
    press_meas[i] = s6[iu]['Pst'][it2]
    aod[i] = s6[iu]['tau_aero'][it2,407]


# In[42]:


s6[iu]['w'][0,407]


# ### 2017

# In[44]:


n = len(gap7['aero_base_day'])
oxa_aero_base7,oxa_meas7,press_meas7,aod7 = np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan 
for i,dd in enumerate(gap7['aero_base_day']):
    iu = np.where(day7==days7[i])[0][0]
    it = np.argmin(abs(dd-s7[iu]['utc']))
    try:
        oxa_aero_base7[i] = s7[iu]['oxa_delta'][it]
    except KeyError:
        continue
    it2 = np.argmin(abs(gap7['meas_low_day'][i]-s7[iu]['utc']))
    oxa_meas7[i] = s7[iu]['oxa_delta'][it2]
    press_meas7[i] = s7[iu]['Pst'][it2]
    aod7[i] = s7[iu]['tau_aero'][it2,407]


# ### 2018

# In[46]:


n = len(gap8['aero_base_day'])
oxa_aero_base8,oxa_meas8,press_meas8,aod8 = np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan 
for i,dd in enumerate(gap8['aero_base_day']):
    iu = np.where(day8==days8[i])[0][0]
    it = np.argmin(abs(dd-s8[iu]['utc']))
    try:
        oxa_aero_base8[i] = s8[iu]['oxa_delta'][it]
    except KeyError:
        continue
    it2 = np.argmin(abs(gap8['meas_low_day'][i]-s8[iu]['utc']))
    oxa_meas8[i] = s8[iu]['oxa_delta'][it2]
    press_meas8[i] = s8[iu]['Pst'][it2]
    aod8[i] = s8[iu]['tau_aero'][it2,407]


# # Plot out data

# In[56]:


plt.figure()
plt.plot(s6[0]['tau_aero'][:,500],s6[0]['oxa_delta'],'.')


# In[47]:


oxa_meas[oxa_meas==0.0] = np.nan


# In[48]:


import plotting_utils as pu


# In[90]:


plt.figure()
i_lowa = (press_meas>880.0) & (press_meas<950.0) & (aod>0.0) & (aod<0.2)
plt.scatter(oxa_meas[i_lowa],gap6['aero_base_alt'][i_lowa],
            label='AOD>0.0 to AOD<0.2, R={:1.2f}'.format(np.corrcoef(gap6['aero_base_alt'][i_lowa],oxa_meas[i_lowa])[1,0]))
pu.plot_lin(oxa_meas[i_lowa],gap6['aero_base_alt'][i_lowa],color='tab:blue',lblfmt='2.4f')

i_lowb = (press_meas>880.0) & (press_meas<950.0)& (aod>0.1) & (aod<0.3)
plt.scatter(oxa_meas[i_lowb],gap6['aero_base_alt'][i_lowb],
            label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(gap6['aero_base_alt'][i_lowb],oxa_meas[i_lowb])[1,0]))
pu.plot_lin(oxa_meas[i_lowb],gap6['aero_base_alt'][i_lowb],color='tab:orange',lblfmt='2.4f')

i_lowc = (press_meas>880.0) & (press_meas<950.0) & (aod>0.4)
plt.scatter(oxa_meas[i_lowc],gap6['aero_base_alt'][i_lowc],
            label='AOD>0.4, R={:1.2f}'.format(np.corrcoef(gap6['aero_base_alt'][i_lowc],oxa_meas[i_lowc])[1,0]))
pu.plot_lin(oxa_meas[i_lowc],gap6['aero_base_alt'][i_lowc],color='tab:green',lblfmt='2.4f')

#i_low = (press_meas>880.0) & (press_meas<950.0) & (aod>0.6)
#plt.scatter(oxa_meas[i_low],gap6['aero_base_alt'][i_low],label='AOD>0.6')
#pu.plot_lin(oxa_meas[i_low],gap6['aero_base_alt'][i_low],color='tab:red')
#plt.colorbar(label='Pressure [mb]')
plt.legend()
plt.xlabel('Depth of Oxygen-A [sum of optical depth deviation]')
plt.ylabel('Aerosol layer base height [m]')
plt.title('ORACLES 2016\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between 880 mb to 950 mb]')
plt.savefig(fp+'ORACLES2016_OXa_vs_aero_alt_{}.png'.format(vv),dpi=600,transparent=True)


# In[62]:


press_meas7


# In[89]:


plt.figure()
i_lowa7 = (press_meas7>850.0) & (press_meas7<920.0) & (aod7>0.0) & (aod7<0.25)
plt.scatter(oxa_meas7[i_lowa7],gap7['aero_base_alt'][i_lowa7],
        label='AOD>0.0 to AOD<0.2, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowa7],oxa_meas7[i_lowa7])[1,0]))
pu.plot_lin(oxa_meas7[i_lowa7],gap7['aero_base_alt'][i_lowa7],color='tab:blue',lblfmt='2.4f')

i_lowb7 = (press_meas7>850.0) & (press_meas7<920.0)& (aod7>0.1) & (aod7<0.3)
plt.scatter(oxa_meas7[i_lowb7],gap7['aero_base_alt'][i_lowb7],
        label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowb7],oxa_meas7[i_lowb7])[1,0]))
pu.plot_lin(oxa_meas7[i_lowb7],gap7['aero_base_alt'][i_lowb7],color='tab:orange',lblfmt='2.4f')

i_lowc7 = (press_meas7>850.0) & (press_meas7<920.0) & (aod7>0.3)
plt.scatter(oxa_meas7[i_lowc7],gap7['aero_base_alt'][i_lowc7],
        label='AOD>0.3, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowc7],oxa_meas7[i_lowc7])[1,0]))
pu.plot_lin(oxa_meas7[i_lowc7],gap7['aero_base_alt'][i_lowc7],color='tab:green',lblfmt='2.4f')

#i_low = (press_meas>880.0) & (press_meas<950.0) & (aod>0.6)
#plt.scatter(oxa_meas[i_low],gap6['aero_base_alt'][i_low],label='AOD>0.6')
#pu.plot_lin(oxa_meas[i_low],gap6['aero_base_alt'][i_low],color='tab:red')
#plt.colorbar(label='Pressure [mb]')
plt.legend()
plt.xlabel('Depth of Oxygen-A [sum of optical depth deviation]')
plt.ylabel('Aerosol layer base height [m]')
plt.title('ORACLES 2017\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between 850 mb to 920 mb]')
plt.savefig(fp+'ORACLES2017_OXa_vs_aero_alt_850mb_{}.png'.format(vv),dpi=600,transparent=True)


# In[88]:


plt.figure()
i_lowa7 = (press_meas7>880.0) & (press_meas7<950.0) & (aod7>0.0) & (aod7<0.25)
plt.scatter(oxa_meas7[i_lowa7],gap7['aero_base_alt'][i_lowa7],
        label='AOD>0.0 to AOD<0.2, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowa7],oxa_meas7[i_lowa7])[1,0]))
pu.plot_lin(oxa_meas7[i_lowa7],gap7['aero_base_alt'][i_lowa7],color='tab:blue',lblfmt='2.4f')

i_lowb7 = (press_meas7>880.0) & (press_meas7<950.0)& (aod7>0.1) & (aod7<0.3)
plt.scatter(oxa_meas7[i_lowb7],gap7['aero_base_alt'][i_lowb7],
        label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowb7],oxa_meas7[i_lowb7])[1,0]))
pu.plot_lin(oxa_meas7[i_lowb7],gap7['aero_base_alt'][i_lowb7],color='tab:orange',lblfmt='2.4f')

i_lowc7 = (press_meas7>880.0) & (press_meas7<950.0) & (aod7>0.3)
plt.scatter(oxa_meas7[i_lowc7],gap7['aero_base_alt'][i_lowc7],
        label='AOD>0.3, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowc7],oxa_meas7[i_lowc7])[1,0]))
pu.plot_lin(oxa_meas7[i_lowc7],gap7['aero_base_alt'][i_lowc7],color='tab:green',lblfmt='2.4f')

#i_low = (press_meas>880.0) & (press_meas<950.0) & (aod>0.6)
#plt.scatter(oxa_meas[i_low],gap6['aero_base_alt'][i_low],label='AOD>0.6')
#pu.plot_lin(oxa_meas[i_low],gap6['aero_base_alt'][i_low],color='tab:red')
#plt.colorbar(label='Pressure [mb]')
plt.legend()
plt.xlabel('Depth of Oxygen-A [sum of optical depth deviation]')
plt.ylabel('Aerosol layer base height [m]')
plt.title('ORACLES 2017\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between 880 mb to 950 mb]')
plt.savefig(fp+'ORACLES2017_OXa_vs_aero_alt_{}.png'.format(vv),dpi=600,transparent=True)


# In[87]:


plt.figure()
i_lowa8 = (press_meas8>850.0) & (press_meas8<930.0) & (aod8>0.0) & (aod8<0.25)
plt.scatter(oxa_meas8[i_lowa8],gap8['aero_base_alt'][i_lowa8],
            label='AOD>0.0 to AOD<0.2, R={:1.2f}'.format(np.corrcoef(gap8['aero_base_alt'][i_lowa8],oxa_meas8[i_lowa8])[1,0]))
pu.plot_lin(oxa_meas8[i_lowa8],gap8['aero_base_alt'][i_lowa8],color='tab:blue',lblfmt='2.4f')

i_lowb8 = (press_meas8>850.0) & (press_meas8<930.0)& (aod8>0.1) & (aod8<0.3)
plt.scatter(oxa_meas8[i_lowb8],gap8['aero_base_alt'][i_lowb8],
            label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(gap8['aero_base_alt'][i_lowb8],oxa_meas8[i_lowb8])[1,0]))
pu.plot_lin(oxa_meas8[i_lowb8],gap8['aero_base_alt'][i_lowb8],color='tab:orange',lblfmt='2.4f')

i_lowc8 = (press_meas8>850.0) & (press_meas8<930.0) & (aod8>0.3)
plt.scatter(oxa_meas8[i_lowc8],gap8['aero_base_alt'][i_lowc8],
            label='AOD>0.3, R={:1.2f}'.format(np.corrcoef(gap8['aero_base_alt'][i_lowc8],oxa_meas8[i_lowc8])[1,0]))
pu.plot_lin(oxa_meas8[i_lowc8],gap8['aero_base_alt'][i_lowc8],color='tab:green',lblfmt='2.4f')

#i_low = (press_meas>880.0) & (press_meas<950.0) & (aod>0.6)
#plt.scatter(oxa_meas[i_low],gap6['aero_base_alt'][i_low],label='AOD>0.6')
#pu.plot_lin(oxa_meas[i_low],gap6['aero_base_alt'][i_low],color='tab:red')
#plt.colorbar(label='Pressure [mb]')
plt.legend()
plt.xlabel('Depth of Oxygen-A [sum of optical depth deviation]')
plt.ylabel('Aerosol layer base height [m]')
plt.title('ORACLES 2018\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between 850 mb to 930 mb]')
plt.savefig(fp+'ORACLES2018_OXa_vs_aero_alt_{}.png'.format(vv),dpi=600,transparent=True)


# In[69]:


gap8['aero_base_alt']


# ## Combined plots

# In[100]:


np.append(np.append(oxa_meas[i_lowa],oxa_meas7[i_lowa7]),oxa_meas8[i_lowa8])


# In[116]:


plt.figure()

plt.scatter(oxa_meas[i_lowa],gap6['aero_base_alt'][i_lowa],marker='.',color='tab:blue')
#            label='AOD>0.0 to AOD<0.2, R={:1.2f}'.format(np.corrcoef(gap6['aero_base_alt'][i_lowa],oxa_meas[i_lowa])[1,0]))
plt.scatter(oxa_meas7[i_lowa7],gap7['aero_base_alt'][i_lowa7],marker='+',color='tab:blue')
#        label='AOD>0.0 to AOD<0.2, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowa7],oxa_meas7[i_lowa7])[1,0]))
#pu.plot_lin(oxa_meas7[i_lowa7],gap7['aero_base_alt'][i_lowa7],color='tab:blue',lblfmt='2.4f')
plt.scatter(oxa_meas8[i_lowa8],gap8['aero_base_alt'][i_lowa8],marker='^',color='tab:blue')
#            label='AOD>0.0 to AOD<0.2, R={:1.2f}'.format(np.corrcoef(gap8['aero_base_alt'][i_lowa8],oxa_meas8[i_lowa8])[1,0]))
#pu.plot_lin(oxa_meas8[i_lowa8],gap8['aero_base_alt'][i_lowa8],color='tab:blue',lblfmt='2.4f')

oxa_a = np.append(np.append(oxa_meas[i_lowa],oxa_meas7[i_lowa7]),oxa_meas8[i_lowa8])
bas_a = np.append(np.append(gap6['aero_base_alt'][i_lowa],gap7['aero_base_alt'][i_lowa7]),gap8['aero_base_alt'][i_lowa8])
plt.scatter(oxa_meas[i_lowa],gap6['aero_base_alt'][i_lowa],marker='.',color='tab:blue',
            label='AOD>0.0 to AOD<0.2, R={:1.2f}'.format(np.corrcoef(bas_a,oxa_a)[1,0]))
pu.plot_lin(oxa_a,bas_a,color='tab:blue',lblfmt='2.4f')

plt.scatter(oxa_meas[i_lowb],gap6['aero_base_alt'][i_lowb],marker='.',color='tab:orange')
#            label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(gap6['aero_base_alt'][i_lowb],oxa_meas[i_lowb])[1,0]))
plt.scatter(oxa_meas7[i_lowb7],gap7['aero_base_alt'][i_lowb7],marker='+',color='tab:orange')
#        label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowb7],oxa_meas7[i_lowb7])[1,0]))
#pu.plot_lin(oxa_meas7[i_lowa7],gap7['aero_base_alt'][i_lowa7],color='tab:blue',lblfmt='2.4f')
plt.scatter(oxa_meas8[i_lowb8],gap8['aero_base_alt'][i_lowb8],marker='^',color='tab:orange')
#            label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(gap8['aero_base_alt'][i_lowb8],oxa_meas8[i_lowb8])[1,0]))
#pu.plot_lin(oxa_meas8[i_lowa8],gap8['aero_base_alt'][i_lowa8],color='tab:blue',lblfmt='2.4f')
plt.scatter(oxa_meas[i_lowb],gap6['aero_base_alt'][i_lowb],marker='.',color='tab:orange', 
            label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(bas_b,oxa_b)[1,0]))

oxa_b = np.append(np.append(oxa_meas[i_lowb],oxa_meas7[i_lowb7]),oxa_meas8[i_lowb8])
bas_b = np.append(np.append(gap6['aero_base_alt'][i_lowb],gap7['aero_base_alt'][i_lowb7]),gap8['aero_base_alt'][i_lowb8])

pu.plot_lin(oxa_b,bas_b,color='tab:orange',lblfmt='2.4f')

plt.scatter(oxa_meas[i_lowc],gap6['aero_base_alt'][i_lowc],marker='.',color='tab:green')
#            label='AOD>0.3, R={:1.2f}'.format(np.corrcoef(gap6['aero_base_alt'][i_lowc],oxa_meas[i_lowc])[1,0]))
plt.scatter(oxa_meas7[i_lowc7],gap7['aero_base_alt'][i_lowc7],marker='+',color='tab:green')
#        label='AOD>0.3, R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowc7],oxa_meas7[i_lowc7])[1,0]))
#pu.plot_lin(oxa_meas7[i_lowa7],gap7['aero_base_alt'][i_lowa7],color='tab:blue',lblfmt='2.4f')
plt.scatter(oxa_meas8[i_lowc8],gap8['aero_base_alt'][i_lowc8],marker='^',color='tab:green')
#            label='AOD>0.3, R={:1.2f}'.format(np.corrcoef(gap8['aero_base_alt'][i_lowc8],oxa_meas8[i_lowc8])[1,0]))
#pu.plot_lin(oxa_meas8[i_lowa8],gap8['aero_base_alt'][i_lowa8],color='tab:blue',lblfmt='2.4f')

oxa_c = np.append(np.append(oxa_meas[i_lowc],oxa_meas7[i_lowc7]),oxa_meas8[i_lowc8])
bas_c = np.append(np.append(gap6['aero_base_alt'][i_lowc],gap7['aero_base_alt'][i_lowc7]),gap8['aero_base_alt'][i_lowc8])
plt.scatter(oxa_meas[i_lowc],gap6['aero_base_alt'][i_lowc],marker='.',color='tab:green',
            label='AOD>0.3, R={:1.2f}'.format(np.corrcoef(bas_c,oxa_c)[1,0]))
pu.plot_lin(oxa_c,bas_c,color='tab:green',lblfmt='2.4f')

plt.legend()
plt.xlabel('Depth of Oxygen-A [sum of optical depth deviation]')
plt.ylabel('Aerosol layer base height [m]')
plt.xlim(0,200000)
plt.title('ORACLES 2016[.], 2017[+], and 2018[^]\nAerosol height and Ox-A depth from 4STAR\n[measured between 850 mb to 930 mb]')
plt.savefig(fp+'ORACLES_all_OXa_vs_aero_alt_{}.png'.format(vv),dpi=600,transparent=True)


# In[ ]:




