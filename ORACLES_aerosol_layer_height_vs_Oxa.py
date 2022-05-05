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

# In[3]:


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
import plotting_utils as pu


# In[4]:


name = 'ORACLES'
vv = 'v1'
fp = getpath(name)


# In[5]:


fpo = '/data/sunsat/ORACLES_2016/data_processed/starsuns/R3/'


# # Load files

# ## Load gap files

# In[6]:


gap6 = np.load(fp+'ORACLES2016_gap_{}.npy'.format(vv),allow_pickle=True,fix_imports=True,encoding='latin1')
gap6 = gap6.item()


# In[7]:


gap7 = np.load(fp+'ORACLES2017_gap_{}.npy'.format(vv),allow_pickle=True,fix_imports=True,encoding='latin1')
gap7 = gap7.item()


# In[8]:


gap8 = np.load(fp+'ORACLES2018_gap_{}.npy'.format(vv),allow_pickle=True,fix_imports=True,encoding='latin1')
gap8 = gap8.item()


# In[9]:


gap6.keys()


# In[10]:


days6 = [d.strftime('%Y%m%d') for d in gap6['aero_base_day']]
day6 = np.unique(days6)


# In[11]:


day6


# In[12]:


days7 = [d.strftime('%Y%m%d') for d in gap7['aero_base_day']]
day7 = np.unique(days7)


# In[13]:


days8 = [d.strftime('%Y%m%d') for d in gap8['aero_base_day']]
day8 = np.unique(days8)


# ## Load starsun files

# ### Load 2016

# In[74]:


s6 = []
for d in day6:
    ss = hs.loadmat(fpo+'4STAR_{}starsun.mat'.format(d),variable_names=['tau_aero','w','t','Pst','m_aero'])
    ss['filename'] = fpo+'4STAR_{}starsun.mat'.format(d)
    s6.append(ss)


# In[75]:


len(s6)


# In[76]:


s6[0].keys()


# In[77]:


for ss in s6:
    ss['utc'] = lu.mat2py_time(ss['t'])


# ### Load 2017

# In[18]:


fpo7 = '/data/sunsat/ORACLES_2017/data_processed/starsuns/R0/'


# In[72]:


s7 = []
for d in day7:
    ss = hs.loadmat(fpo7+'4STAR_{}starsun.mat'.format(d),variable_names=['tau_aero','w','t','Pst','m_aero'])
    ss['filename'] = fpo7+'4STAR_{}starsun.mat'.format(d)
    print(ss['filename'],'/',len(day7))
    s7.append(ss)


# In[73]:


for ss in s7:
    ss['utc'] = lu.mat2py_time(ss['t'])


# ### Load 2018

# In[21]:


fpo8 = '/data/sunsat/ORACLES_2018/data_processed/starsuns/R0/'


# In[67]:


s8 = []
for d in day8:
    ss = hs.loadmat(fpo8+'4STAR_{}starsun.mat'.format(d),variable_names=['tau_aero','w','t','Pst','m_aero'])
    ss['filename'] = fpo8+'4STAR_{}starsun.mat'.format(d)
    print(ss['filename'])
    s8.append(ss)


# In[68]:


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


# In[23]:


def oxa_depth_v2(w,spec):
    'calculate the oxygen-a depth. w: wavelenght in microns, spec: spectra(tau_aero)'
    y0 = np.argmin(abs(w-0.756))
    y1 = np.argmin(abs(w-0.772))
    oxa_flat,oxa_ratio,oxa_ratio2,oxa_delta = [],[],[],[]
    for i in range(len(spec)):
        fx = interpolate.interp1d(w[0,[y0,y1]],spec[i,[y0,y1]])
        oxa_flat.append(fx(w[0,y0:y1]))
        oxa_delta.append(np.nansum(spec[i,y0:y1]-oxa_flat))
        oxa_ratio.append(np.nansum(spec[i,y0:y1]/oxa_flat))
        oxa_ratio2.append(np.nanmean(spec[i,y0:y1]/oxa_flat))
    return np.array(oxa_delta),np.array(oxa_ratio),np.array(oxa_ratio2)


# In[72]:


ss['w'].flatten()


# In[73]:


w = ss['w'].flatten()
spec = ss['tau_aero']


# In[89]:


y0 = np.argmin(abs(w-0.756))
y1 = np.argmin(abs(w-0.772))
ym = np.argmin(abs(w-0.761))
oxa_flat = [(s[y1]-s[y0])/(w[y1]-w[y0])*(w[y0:y1]-w[y0])+s[y0] for s in spec]
oxa_delta = [np.nansum(s[y0:y1]-oxa_flat[i]) for i,s in enumerate(spec)]


# In[113]:


ym-y0


# In[87]:


for s in spec[5111:]:
    print(s.shape)
    print((s[y1]-s[y0])/(w[y1]-w[y0])*(w[y0:y1]-w[y0])+s[y0])
    print((s[y1],s[y0]))
    print((w[y1],w[y0]))
    break


# In[79]:


def oxa_depth_v3(w,spec):
    'calculate the oxygen-a depth. w: wavelenght in microns, spec: spectra(tau_aero)'
    w = w.flatten()
    y0 = np.argmin(abs(w-0.756))
    y1 = np.argmin(abs(w-0.772))
    ym = np.argmin(abs(w-0.761))
    oxa_flat = [(s[y1]-s[y0])/(w[y1]-w[y0])*(w[y0:y1]-w[y0])+s[y0] for s in spec]
    oxa_delta = np.array([np.nansum(s[y0:y1]-oxa_flat[i]) for i,s in enumerate(spec)])
    oxa_ratio = np.array([np.nansum(s[y0:y1]/oxa_flat[i]) for i,s in enumerate(spec)])
    oxa_ratio2 = np.array([np.nanmean(s[y0:y1]/oxa_flat[i]) for i,s in enumerate(spec)])
    oxa_depthr = np.array([s[ym]/oxa_flat[i][ym-y0] for i,s in enumerate(spec)])
    oxa_depthm = np.array([s[ym]-oxa_flat[i][ym-y0] for i,s in enumerate(spec)])
    return oxa_delta,oxa_ratio,oxa_ratio2,oxa_depthr,oxa_depthm


# ### 2016 Ox-A

# In[25]:


plt.figure()
plt.plot(ss['w'][0,:],ss['tau_aero'][6111,:],'.-')
plt.xlim(0.750,0.775)
plt.ylim(0,0.15)


# In[80]:


for j,ss in enumerate(s6):
    print('Doing file {}, which is {}/{}'.format(day6[j],j,len(s6)))
    if 'oxa_delta' in ss.keys(): 
        print('Yas queen already in there')
        #continue
    try:
        ss['oxa_delta'],ss['oxa_ratio'],ss['oxa_ratio2'],ss['oxa_depthr'],ss['oxa_depthm'] = oxa_depth_v3(ss['w'],ss['tau_aero'])
    except Exception as e:
        print(e)


# ### 2017 Ox-A

# In[81]:


for j,ss in enumerate(s7):
    print('Doing file {}, which is {}/{}'.format(day7[j],j,len(s7)))
    if 'oxa_delta' in ss.keys(): 
        print('Yas queen already in there')
        #continue
    try:
        ss['oxa_delta'],ss['oxa_ratio'],ss['oxa_ratio2'],ss['oxa_depthr'],ss['oxa_depthm'] = oxa_depth_v3(ss['w'],ss['tau_aero'])
    except:
        pass


# ### 2018 Ox-A

# In[82]:


for j,ss in enumerate(s8):
    print('Doing file {}, which is {}/{}'.format(day8[j],j,len(s8)))
    if 'oxa_delta' in ss.keys(): 
        print('Yas queen already in there')
        #continue
    try:
        ss['oxa_delta'],ss['oxa_ratio'],ss['oxa_ratio2'],ss['oxa_depthr'],ss['oxa_depthm'] = oxa_depth_v3(ss['w'],ss['tau_aero'])
    except:
        pass


# ## Get the oxa depth for the gap times

# In[83]:


dd = gap6['aero_base_day'][1]


# In[84]:


utc_fx = lambda x: x.hour+(x.minute+(x.second/60.0))/60.0


# In[85]:


iu = 2
s6[iu].keys()


# In[86]:


len(gap6['aero_base_alt'])


# ### 2016

# In[87]:


navg = 3


# In[88]:


n = len(gap6['aero_base_day'])
oxa_aero_base,oxa_meas,press_meas,aod = np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan 
oxa_measr,oxa_measr2 = np.zeros((n))+np.nan,np.zeros((n))+np.nan
oxa_measdr,oxa_measdm = np.zeros((n))+np.nan,np.zeros((n))+np.nan
m_aero = np.zeros((n))+np.nan
for i,dd in enumerate(gap6['aero_base_day']):
    iu = np.where(day6==days6[i])[0][0]
    it = np.argmin(abs(dd-s6[iu]['utc']))
    try:
        oxa_aero_base[i] = s6[iu]['oxa_delta'][it]
    except KeyError:
        continue
    it2 = np.argmin(abs(gap6['meas_low_day'][i]-s6[iu]['utc']))
    oxa_meas[i] = np.nanmean(s6[iu]['oxa_delta'][it2-navg:it2+navg])
    oxa_measr[i] = np.nanmean(s6[iu]['oxa_ratio'][it2-navg:it2+navg])
    oxa_measr2[i] = np.nanmean(s6[iu]['oxa_ratio2'][it2-navg:it2+navg])
    oxa_measdr[i] = np.nanmean(s6[iu]['oxa_depthr'][it2-navg:it2+navg])
    oxa_measdm[i] = np.nanmean(s6[iu]['oxa_depthm'][it2-navg:it2+navg])
    press_meas[i] = np.nanmean(s6[iu]['Pst'][it2-navg:it2+navg])
    m_aero[i] = np.nanmean(s6[iu]['m_aero'][it2-navg:it2+navg])
    aod[i] = s6[iu]['tau_aero'][it2,407]


# In[89]:


s6[iu]['w'][0,407]


# ### 2017

# In[94]:


n = len(gap7['aero_base_day'])
oxa_aero_base7,oxa_meas7,press_meas7,aod7 = np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan 
oxa_measr7,oxa_measr27 = np.zeros((n))+np.nan,np.zeros((n))+np.nan
oxa_measdr7,oxa_measdm7 = np.zeros((n))+np.nan,np.zeros((n))+np.nan
m_aero7 = np.zeros((n))+np.nan
for i,dd in enumerate(gap7['aero_base_day']):
    iu = np.where(day7==days7[i])[0][0]
    it = np.argmin(abs(dd-s7[iu]['utc']))
    try:
        oxa_aero_base7[i] = s7[iu]['oxa_delta'][it]
    except KeyError:
        continue
    it2 = np.argmin(abs(gap7['meas_low_day'][i]-s7[iu]['utc']))
    oxa_meas7[i] = np.nanmean(s7[iu]['oxa_delta'][it2-navg:it2+navg])
    oxa_measr7[i] = np.nanmean(s7[iu]['oxa_ratio'][it2-navg:it2+navg])
    oxa_measr27[i] = np.nanmean(s7[iu]['oxa_ratio2'][it2-navg:it2+navg])
    oxa_measdr7[i] = np.nanmean(s7[iu]['oxa_depthr'][it2-navg:it2+navg])
    oxa_measdm7[i] = np.nanmean(s7[iu]['oxa_depthm'][it2-navg:it2+navg])
    press_meas7[i] = s7[iu]['Pst'][it2]
    m_aero7[i] = s7[iu]['m_aero'][it2,0]
    aod7[i] = s7[iu]['tau_aero'][it2,407]


# ### 2018

# In[95]:


n = len(gap8['aero_base_day'])
oxa_aero_base8,oxa_meas8,press_meas8,aod8 = np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan,np.zeros((n))+np.nan 
oxa_measr8,oxa_measr28 = np.zeros((n))+np.nan,np.zeros((n))+np.nan
oxa_measdr8,oxa_measdm8 = np.zeros((n))+np.nan,np.zeros((n))+np.nan
m_aero8 = np.zeros((n))+np.nan
for i,dd in enumerate(gap8['aero_base_day']):
    iu = np.where(day8==days8[i])[0][0]
    it = np.argmin(abs(dd-s8[iu]['utc']))
    try:
        oxa_aero_base8[i] = s8[iu]['oxa_delta'][it]
    except KeyError:
        continue
    it2 = np.argmin(abs(gap8['meas_low_day'][i]-s8[iu]['utc']))
    oxa_meas8[i] = np.nanmean(s8[iu]['oxa_delta'][it2-navg:it2+navg])
    oxa_measr8[i] = np.nanmean(s8[iu]['oxa_ratio'][it2-navg:it2+navg])
    oxa_measr28[i] = np.nanmean(s8[iu]['oxa_ratio2'][it2-navg:it2+navg])
    oxa_measdr8[i] = np.nanmean(s8[iu]['oxa_depthr'][it2-navg:it2+navg])
    oxa_measdm8[i] = np.nanmean(s8[iu]['oxa_depthm'][it2-navg:it2+navg])
    press_meas8[i] = s8[iu]['Pst'][it2]
    m_aero8[i] = s8[iu]['m_aero'][it2,0]
    aod8[i] = s8[iu]['tau_aero'][it2,407]


# # Plot out data

# In[97]:


plt.figure()
plt.plot(s6[0]['tau_aero'][:,500],s6[0]['oxa_delta'],'.')


# In[96]:


oxa_meas[oxa_meas==0.0] = np.nan


# In[99]:


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


# In[43]:


press_meas7


# In[44]:


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


# In[45]:


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


# In[46]:


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


# In[47]:


gap8['aero_base_alt']


# ## Combined plots

# In[48]:


np.append(np.append(oxa_meas[i_lowa],oxa_meas7[i_lowa7]),oxa_meas8[i_lowa8])


# In[105]:


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


oxa_b = np.append(np.append(oxa_meas[i_lowb],oxa_meas7[i_lowb7]),oxa_meas8[i_lowb8])
bas_b = np.append(np.append(gap6['aero_base_alt'][i_lowb],gap7['aero_base_alt'][i_lowb7]),gap8['aero_base_alt'][i_lowb8])
plt.scatter(oxa_meas[i_lowb],gap6['aero_base_alt'][i_lowb],marker='.',color='tab:orange', 
            label='AOD>0.1 to AOD<0.3, R={:1.2f}'.format(np.corrcoef(bas_b,oxa_b)[1,0]))

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
plt.xlim(0,3.5)
plt.title('ORACLES 2016[.], 2017[+], and 2018[^]\nAerosol height and Ox-A depth from 4STAR\n[measured between 850 mb to 930 mb]')
plt.savefig(fp+'ORACLES_all_OXa_vs_aero_alt_{}.png'.format(vv),dpi=600,transparent=True)


# ## Plot out the different oxa measurements

# In[106]:


oxa_dict6 = [{'o':oxa_meas,'label':'Depth of Oxygen-A [sum of optical depth deviation]'},
             {'o':oxa_measr,'label':'Depth of Oxygen-A [sum of optical depth ratio]'},
             {'o':oxa_measr2,'label':'Depth of Oxygen-A [mean of optical depth ratio]'},
             {'o':oxa_measdr,'label':'Depth of Oxygen-A [optical depth ratio at 761 nm]'},
             {'o':oxa_measdm,'label':'Depth of Oxygen-A [optical depth difference at 761 nm]'}]


# In[107]:


oxa_dict7 = [{'o':oxa_meas7,'label':'Depth of Oxygen-A [sum of optical depth deviation]'},
             {'o':oxa_measr7,'label':'Depth of Oxygen-A [sum of optical depth ratio]'},
             {'o':oxa_measr27,'label':'Depth of Oxygen-A [mean of optical depth ratio]'},
             {'o':oxa_measdr7,'label':'Depth of Oxygen-A [optical depth ratio at 761 nm]'},
             {'o':oxa_measdm7,'label':'Depth of Oxygen-A [optical depth difference at 761 nm]'}]


# In[108]:


oxa_dict8 = [{'o':oxa_meas8,'label':'Depth of Oxygen-A [sum of optical depth deviation]'},
             {'o':oxa_measr8,'label':'Depth of Oxygen-A [sum of optical depth ratio]'},
             {'o':oxa_measr28,'label':'Depth of Oxygen-A [mean of optical depth ratio]'},
             {'o':oxa_measdr8,'label':'Depth of Oxygen-A [optical depth ratio at 761 nm]'},
             {'o':oxa_measdm8,'label':'Depth of Oxygen-A [optical depth difference at 761 nm]'}]


# In[54]:


for ox in oxa_dict6:
    plt.figure()
    i_lowa = (press_meas>820.0) & (press_meas<980.0) & (aod>0.0) & np.isfinite(ox['o']) & (ox['o']!=0.0) # & (aod<0.2)
    plt.scatter(ox['o'][i_lowa],gap6['aero_base_alt'][i_lowa],
                label='R={:1.2f}'.format(np.corrcoef(gap6['aero_base_alt'][i_lowa],ox['o'][i_lowa])[1,0]),c=aod[i_lowa])
    pu.plot_lin(ox['o'][i_lowa],gap6['aero_base_alt'][i_lowa],color='tab:blue',lblfmt='2.4f')
    plt.colorbar(label='AOD')
    plt.legend()
    plt.xlabel(ox['label'])
    plt.ylabel('Aerosol layer base height [m]')
    plt.title('ORACLES 2016\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between 720 mb to 950 mb]')


# In[55]:


for ox in oxa_dict7:
    plt.figure()
    i_lowa = (press_meas7>820.0) & (press_meas7<980.0) & (aod7>0.0) & np.isfinite(ox['o']) & (ox['o']!=0.0)#  & (aod<0.2)
    plt.scatter(ox['o'][i_lowa],gap7['aero_base_alt'][i_lowa],
                label='R={:1.2f}'.format(np.corrcoef(gap7['aero_base_alt'][i_lowa],ox['o'][i_lowa])[1,0]),c=press_meas7[i_lowa])
    pu.plot_lin(ox['o'][i_lowa],gap7['aero_base_alt'][i_lowa],color='tab:blue',lblfmt='2.4f')
    plt.colorbar(label='Press')
    plt.legend()
    plt.xlabel(ox['label'])
    plt.ylabel('Aerosol layer base height [m]')
    plt.title('ORACLES 2017\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between 720 mb to 950 mb]')


# In[56]:


for ox in oxa_dict8:
    plt.figure()
    i_lowa = (press_meas8>820.0) & (press_meas8<980.0) & (aod8>0.0) & np.isfinite(ox['o']) & (ox['o']!=0.0) # & (aod<0.2)
    plt.scatter(ox['o'][i_lowa],gap8['aero_base_alt'][i_lowa],
                label='R={:1.2f}'.format(np.corrcoef(gap8['aero_base_alt'][i_lowa],ox['o'][i_lowa])[1,0]),c=aod8[i_lowa])
    pu.plot_lin(ox['o'][i_lowa],gap8['aero_base_alt'][i_lowa],color='tab:blue',lblfmt='2.4f')
    plt.colorbar(label='AOD')
    plt.legend()
    plt.xlabel(ox['label'])
    plt.ylabel('Aerosol layer base height [m]')
    plt.title('ORACLES 2018\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between 720 mb to 950 mb]')


# In[57]:


for ox in oxa_dict8:
    plt.figure()
    i_lowa = (press_meas8>720.0) & (press_meas8<950.0) & (aod8>0.0) & np.isfinite(ox['o']) & (ox['o']!=0.0) # & (aod<0.2)
    plt.scatter(ox['o'][i_lowa],gap8['aero_base_alt'][i_lowa],
                label='R={:1.2f}'.format(np.corrcoef(gap8['aero_base_alt'][i_lowa],ox['o'][i_lowa])[1,0]),c=press_meas8[i_lowa])
    pu.plot_lin(ox['o'][i_lowa],gap8['aero_base_alt'][i_lowa],color='tab:blue',lblfmt='2.4f')
    plt.colorbar(label='Pressure [hPa]')
    plt.legend()
    plt.xlabel(ox['label'])
    plt.ylabel('Aerosol layer base height [m]')
    plt.title('ORACLES 2018\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between 720 mb to 950 mb]')


# In[58]:


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


# ## Multi year combination plots

# In[109]:


press = np.append(np.append(press_meas,press_meas7),press_meas8)
m = np.append(np.append(m_aero,m_aero7),m_aero8)
aodd = np.append(np.append(aod,aod7),aod8)
aero_base_alt = np.append(np.append(gap6['aero_base_alt'],gap7['aero_base_alt']),gap8['aero_base_alt'])
aero_base_UTC = np.append(np.append(gap6['aero_base_UTC'],gap7['aero_base_UTC']),gap8['aero_base_UTC'])
oxa = np.append(np.append(oxa_measr2,oxa_measr27),oxa_measr28)


# In[115]:


plt.figure()
plt.hist(m,bins=50)


# In[120]:


plt.figure()
p2,p1 = 820.0,980.0
i_lowa = (press>p2) & (press<p1) & (aodd>0.0) & np.isfinite(oxa) & (oxa!=0.0) & (m<1.2) & (oxa<3.0) # & (aod<0.2)
plt.scatter(oxa[i_lowa],aero_base_alt[i_lowa],
            label='R={:1.2f}'.format(np.corrcoef(aero_base_alt[i_lowa],oxa[i_lowa])[1,0]),c=press[i_lowa])
pu.plot_lin(oxa[i_lowa],aero_base_alt[i_lowa],x_err=oxa[i_lowa]*0.03,y_err=aero_base_alt[i_lowa]*0.05,color='tab:blue',lblfmt='2.4f',use_method='york')
plt.colorbar(label='Press [mb]')
plt.legend()
plt.xlabel('OxA depth - mean OD ratio')
plt.ylabel('Aerosol layer base height [m]')
plt.title('ORACLES\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between {:3.0f} mb to {:3.0f} mb]'.format(p1,p2))


# In[119]:


plt.figure()
p2,p1 = 820.0,980.0
i_lowa = (press>p2) & (press<p1) & (aodd>0.0) & np.isfinite(oxa) & (oxa!=0.0) & (m<1.8) & (oxa<3.0) # & (aod<0.2)
plt.scatter(oxa[i_lowa],aero_base_alt[i_lowa],
            label='R={:1.2f}'.format(np.corrcoef(aero_base_alt[i_lowa],oxa[i_lowa])[1,0]),c=m[i_lowa])
pu.plot_lin(oxa[i_lowa],aero_base_alt[i_lowa],x_err=oxa[i_lowa]*0.03,y_err=aero_base_alt[i_lowa]*0.05,color='tab:blue',lblfmt='2.4f',use_method='york')
plt.colorbar(label='Air Mass factor')
plt.legend()
plt.xlabel('OxA depth - mean OD ratio')
plt.ylabel('Aerosol layer base height [m]')
plt.title('ORACLES\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between {:3.0f} mb to {:3.0f} mb]'.format(p1,p2))


# In[129]:


plt.figure()
p2,p1 = 820.0,980.0
i_lowa = (press>p2) & (press<p1) & (aodd>0.0) & np.isfinite(oxa) & (oxa!=0.0) & (m<1.2) & (oxa<3.0) # & (aod<0.2)
plt.scatter(aero_base_alt[i_lowa]/1000.0,oxa[i_lowa],
            label='{}'.format(pu.stats_label(aero_base_alt[i_lowa]/1000.0,oxa[i_lowa])),c=m[i_lowa])
pu.plot_lin(aero_base_alt[i_lowa]/1000.0,oxa[i_lowa],y_err=oxa[i_lowa]*0.03,x_err=aero_base_alt[i_lowa]/1000.0*0.05,color='tab:blue',lblfmt='2.4f')
plt.colorbar(label='Air Mass factor')
plt.legend()
plt.ylabel('OxA depth - mean OD ratio')
plt.xlabel('Aerosol layer base height [km]')
plt.title('ORACLES\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between {:3.0f} mb to {:3.0f} mb]'.format(p1,p2))


# In[121]:


pu.stats_label(oxa[i_lowa],aero_base_alt[i_lowa])


# In[65]:


plt.figure()
p2,p1 = 820.0,980.0
i_lowa = (press>p2) & (press<p1) & (aodd>0.0) & np.isfinite(oxa) & (oxa!=0.0) # & (aod<0.2)
plt.scatter(oxa[i_lowa],aero_base_alt[i_lowa],
            label='R={:1.2f}'.format(np.corrcoef(aero_base_alt[i_lowa],oxa[i_lowa])[1,0]),c=aodd[i_lowa])
pu.plot_lin(oxa[i_lowa],aero_base_alt[i_lowa],color='tab:blue',lblfmt='2.4f')
plt.colorbar(label='AOD')
plt.legend()
plt.xlabel('OxA depth - mean OD ratio')
plt.ylabel('Aerosol layer base height [m]')
plt.title('ORACLES\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between {:3.0f} mb to {:3.0f} mb]'.format(p1,p2))


# In[64]:


plt.figure()
p2,p1 = 880.0,950.0
i_lowa = (press>p2) & (press<p1) & (aodd>0.0) & np.isfinite(oxa) & (oxa!=0.0) & (aero_base_UTC<14.0) #& (aodd<0.6) & (oxa>-0.7)
plt.scatter(oxa[i_lowa],aero_base_alt[i_lowa],
            label='R={:1.2f}'.format(np.corrcoef(aero_base_alt[i_lowa],oxa[i_lowa])[1,0]),c=aero_base_UTC[i_lowa])
pu.plot_lin(oxa[i_lowa],aero_base_alt[i_lowa],color='tab:blue',lblfmt='2.4f')
plt.colorbar(label='UTC')
plt.legend()
plt.xlabel('OxA depth - mean OD ratio')
plt.ylabel('Aerosol layer base height [m]')
plt.title('ORACLES\nAerosol height, Ox-A, and AOD from 4STAR\n[measured between {:3.0f} mb to {:3.0f} mb]'.format(p1,p2))


# In[ ]:




