#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_DARE_explore
# 
# Purpose:  
# 
#     view the calculated DARE values frrom the 4STAR skyscans, aod, and cloud retrieval
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
#   
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-12-03
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2020-04-14, under stay-at-home order for coronavirus
#          - Adding version 3 plotting
# 

# # Prepare python environment

# In[1]:


import numpy as np
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
import os


# In[2]:


import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
import plotting_utils as pu


# In[3]:


get_ipython().magic(u'matplotlib notebook')


# In[4]:


name = 'ORACLES'
vv = 'v4'
fp = getpath(name)


# # Load files

# ## Load DARE calculations

# In[91]:


s = hs.loadmat(fp+'ORACLES_DARE_{}.mat'.format(vv))


# In[92]:


s.keys()


# In[93]:


s['dare'].shape


# In[7]:


ss = hs.loadmat(fp+'ORACLES_DARE_aero_prop_{}.mat'.format(vv))


# In[8]:


s['sza'] = ss['sza']


# ### Load the instant to 24h DARE ratio

# In[9]:


import pandas as pd


# In[10]:


fp


# In[11]:


inst2day = pd.read_csv(fp+'data_other/DARE_instant_over_24h_Redemann_2006.dat')


# In[12]:


i2d = inst2day.to_numpy() #mu, ratio (i/24h)


# In[13]:


i2d


# In[94]:


s['dare'].shape


# In[370]:


s['mu'] = np.cos(np.deg2rad(s['sza']))


# In[371]:


i_fac = [np.argmin(abs(i2d[:,0]-mu)) for mu in s['mu']]


# In[95]:


s['dare_24h'] = s['dare']*2.0+np.nan


# In[96]:


s['dare_24h'][:,0] = s['dare'][:,0]/i2d[i_fac,1]
s['dare_24h'][:,1] = s['dare'][:,1]/i2d[i_fac,1]
s['dare_24h'][:,2] = s['dare'][:,2]/i2d[i_fac,1]


# ## Load the 24h DARE calculations

# In[41]:


sh = hs.loadmat(fp+'ORACLES_DARE_{}_24h.mat'.format(vv))


# In[19]:


sh.keys()


# In[20]:


utcx = np.arange(0,24,0.5)


# In[21]:


ig = np.where(np.isfinite(sh['dare_avg'][:,2]) & (sh['dare_avg'][:,2] != 0.0) & (abs(sh['dare_avg'][:,2]) > 0.25))[0]


# ### Check the daily edges

# In[22]:


from matplotlib import animation, rc


# In[23]:


print sh['dn'][ig[0],:,2] - sh['dn_noa'][ig[0],:,2]


# In[24]:


sh['zout']


# In[25]:


fig = plt.figure()
line = []
line.append(plt.plot(utcx,sh['dare'][ig[0],:,2]*100.0,'.',label='DAREx100')[0])
line.append(plt.plot(utcx,sh['up'][ig[0],:,2],'.',label='Up')[0])
line.append(plt.plot(utcx,sh['dn'][ig[0],:,2],'.',label='Dn')[0])
line.append(plt.plot(utcx,sh['up_noa'][ig[0],:,2],'.',label='Up clear')[0])
line.append(plt.plot(utcx,sh['dn_noa'][ig[0],:,2],'.',label='Dn clear')[0])
plt.xlabel('UTC [h]')
plt.ylabel('Irradiance [W/m$^2$]')
plt.legend()
plt.ylim(-20,1300)


# In[29]:


def init():
    for l in line:
        l.set_data([], [])
    return (line,)


# In[26]:


def animate(i):
    line[0].set_data(utcx, sh['dare'][ig[i],:,2]*100.0)
    line[1].set_data(utcx, sh['up'][ig[i],:,2])
    line[2].set_data(utcx, sh['dn'][ig[i],:,2])
    line[3].set_data(utcx, sh['up_noa'][ig[i],:,2])
    line[4].set_data(utcx, sh['dn_noa'][ig[i],:,2])
    line[0].axes.set_title('{:d}'.format(ig[i]))
    return (line,)


# In[27]:


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(ig), interval=1, 
                               blit=True)


# In[30]:


fig = plt.figure()
line = []
line.append(plt.plot(utcx,sh['dare'][ig[0],:,2]*100.0,'.',label='DAREx100')[0])
line.append(plt.plot(utcx,sh['up'][ig[0],:,2],'.',label='Up')[0])
line.append(plt.plot(utcx,sh['dn'][ig[0],:,2],'.',label='Dn')[0])
line.append(plt.plot(utcx,sh['up_noa'][ig[0],:,2],'.',label='Up clear')[0])
line.append(plt.plot(utcx,sh['dn_noa'][ig[0],:,2],'.',label='Dn clear')[0])
plt.xlabel('UTC [h]')
plt.ylabel('Irradiance [W/m$^2$]')
plt.legend()
plt.ylim(-20,1500)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(ig), interval=1, 
                               blit=True)
anim


# In[31]:


plt.figure()
plt.plot(sh['dare'][ig,10,2],'.',label='edge-1')
plt.plot(sh['dare'][ig,11,2],'.',label='edge')
plt.plot(sh['dare'][ig,12,2],'.',label='next')
plt.legend()


# In[68]:


plt.figure()
plt.plot(sh['dare'][ig,11,2]/sh['dare'][ig,12,2],'.')


# In[32]:


plt.figure()
plt.plot(sh['dare'][ig,33,2],'.',label='edge-1')
plt.plot(sh['dare'][ig,34,2],'.',label='edge')
plt.plot(sh['dare'][ig,35,2],'.',label='next')
plt.legend()


# ### Adjust the calculated DARE for issues with the edges

# In[33]:


bad_edge = abs(sh['dare'][:,11,2])>abs(sh['dare'][:,12,2]) 
bad_edge2 = (abs(sh['dare'][:,34,2])>abs(sh['dare'][:,33,2]))
bad_edge3 = (abs(sh['dare'][:,35,2])>abs(sh['dare'][:,34,2]))


# In[34]:


sh['dare'][bad_edge,11,2] = sh['dare'][bad_edge,12,2]*0.1
sh['dare'][bad_edge2,34,2] = sh['dare'][bad_edge2,33,2]*0.1
sh['dare'][bad_edge3,35,2] = sh['dare'][bad_edge3,34,2]*0.1


# In[35]:


a = np.array([[1,2,3],[2,3,4],[4,5,6],[7,8,9]])


# In[36]:


a.shape


# In[37]:


np.nanmean(a,axis=0)


# In[38]:


sh['dare_avg'][:,2] = np.nanmean(sh['dare'][:,:,2],axis=1)


# ## Load the DARE 24h with diurnal changing clouds

# In[42]:


shd = hs.loadmat(fp+'ORACLES_DARE_{}_24h_diurnal.mat'.format(vv))


# In[43]:


shd.keys()


# ### Check and adjust the daily edges

# In[44]:


plt.figure()
plt.plot(shd['dare'][ig,10,2],'.',label='edge-1')
plt.plot(shd['dare'][ig,11,2],'.',label='edge')
plt.plot(shd['dare'][ig,12,2],'.',label='next')
plt.legend()


# In[45]:


bad_edged = abs(shd['dare'][:,11,2])>abs(shd['dare'][:,12,2]) 
bad_edge2d = (abs(shd['dare'][:,34,2])>abs(shd['dare'][:,33,2]))
bad_edge3d = (abs(shd['dare'][:,35,2])>abs(shd['dare'][:,34,2]))


# In[46]:


shd['dare'][bad_edged,11,2] = shd['dare'][bad_edged,12,2]*0.1
shd['dare'][bad_edge2d,34,2] = shd['dare'][bad_edge2d,33,2]*0.1
shd['dare'][bad_edge3d,35,2] = shd['dare'][bad_edge3d,34,2]*0.1


# In[47]:


shd['dare_avg'][:,2] = np.nanmean(shd['dare'][:,:,2],axis=1)


# ## Load DARE parameters from SARE

# In[130]:


sa = sio.loadmat(fp+'ORACLES_2016_DARE_params_{}_px.mat'.format('v4'))


# In[131]:


sa.keys()


# ## Load DARE from ALADIN

# In[53]:


AL_dare_2016 = np.array([3.68,4.86,4.34,3.79,1.31,0.48,0.26,-0.02])
AL_Lat = np.array([-14.0,-12.0,-10.0,-8.0,-6.0,-4.0,-2.0,0.0])


# # Plot out data

# In[101]:


s['doy']


# In[102]:


s['doys'] = s['doy']+s['utc']/24.0


# In[103]:


ibad = (s['ref']==0.0) | (s['ref']>=25.0) | (s['cod']==0.0)


# In[104]:


s['cod'][ibad] = np.nan
s['ref'][ibad] = np.nan
s['dn'][ibad,:] = np.nan
s['up'][ibad,:] = np.nan
s['dn_noa'][ibad,:] = np.nan
s['up_noa'][ibad,:] = np.nan
s['dare'][ibad,:] = np.nan


# ## Plot the time trace

# In[58]:


plt.figure()
plt.plot(s['doys'],s['cod'])
plt.xlabel('DOY')
plt.ylabel('COD')
plt.grid()
plt.title('ORACLES 2016 SSFR retrieved COD')


# In[59]:


plt.figure()
plt.plot(s['doys'],s['ref'],'.')
plt.xlabel('DOY')
plt.ylabel('REF')
plt.grid()
plt.title('ORACLES 2016 SSFR retrieved REF')


# ## histogram of retrieved properties

# In[60]:


plt.figure()
plt.hist(s['cod'],range=[0,60],bins=30)
plt.xlabel('COD')
plt.title('ORACLES 2016 SSFR COD histogram')


# In[23]:


plt.figure()
plt.hist(s['ref'],range=[0,30],bins=30)
plt.xlabel('REF [${{\\mu}}$m]')
plt.title('ORACLES 2016 SSFR REF histogram')


# In[61]:


plt.figure()
plt.hist(s['ref'],range=[0,60],bins=30,label='R$_{{eff}}$',alpha=0.6,normed=True,color='g')
plt.hist(s['cod'],range=[0,60],bins=30,label='COD',alpha=0.6,normed=True,color='m')

plt.axvline(np.nanmean(s['ref']),ls='-',color='g',label='mean')
plt.axvline(np.nanmedian(s['ref']),ls='--',color='g',label='median')

plt.axvline(np.nanmean(s['cod']),ls='-',color='m')
plt.axvline(np.nanmedian(s['cod']),ls='--',color='m')
plt.legend(frameon=False)
plt.title('ORACLES 2016 SSFR above cloud retrievals')

plt.xlabel('R$_{{eff}}$ [$\mu$m], COD')
plt.savefig(fp+'plot_DARE/ORACLES2016_COD_REF_hist_{}.png'.format(vv),dpi=600,transparent=True)


# In[29]:


icod= s['cod']>0.0
plt.figure(figsize=(4,3))
#plt.hist(s['ref'],range=[0,60],bins=30,label='R$_{{eff}}$',alpha=0.6,normed=True,color='g')
plt.hist(s['cod'][icod],range=[0,60],bins=30,label='Cloud\nOptical Depth',alpha=0.6,normed=True,color='m')

#plt.axvline(np.nanmean(s['ref']),ls='-',color='g',label='mean')
#plt.axvline(np.nanmedian(s['ref']),ls='--',color='g',label='median')

plt.axvline(np.nanmean(s['cod'][icod]),ls='-',color='m',label='mean: {:3.1f}'.format(np.nanmean(s['cod'])))
plt.axvline(np.nanmedian(s['cod'][icod]),ls='--',color='m',label='median: {:3.1f}'.format(np.nanmedian(s['cod'])))
plt.legend(frameon=False)
plt.title('ORACLES 2016 SSFR above cloud retrievals')

plt.xlabel('COD')
plt.tight_layout()

plt.savefig(fp+'plot_DARE/ORACLES2016_COD_SSFR_hist_{}.png'.format(vv),dpi=600,transparent=True)


# In[30]:


iref= s['ref']>0.0
plt.figure(figsize=(4,3))
plt.hist(s['ref'][iref],range=[0,30],bins=30,label='Cloud drop\neffective radius',alpha=0.6,normed=True,color='g')
#plt.hist(s['cod'][icod],range=[0,60],bins=30,label='COD',alpha=0.6,normed=True,color='m')

plt.axvline(np.nanmean(s['ref'][iref]),ls='-',color='g',label='mean: {:3.1f}'.format(np.nanmean(s['ref'])))
plt.axvline(np.nanmedian(s['ref'][iref]),ls='--',color='g',label='median: {:3.1f}'.format(np.nanmedian(s['ref'])))

#plt.axvline(np.nanmean(s['cod'][icod]),ls='-',color='m')
#plt.axvline(np.nanmedian(s['cod'][icod]),ls='--',color='m')
plt.legend(frameon=False)
plt.title('ORACLES 2016 SSFR above cloud retrievals')

plt.xlabel('R$_{{eff}}$ [$\mu$m]')
plt.tight_layout()

plt.savefig(fp+'plot_DARE/ORACLES2016_REF_SSFR_hist_{}.png'.format(vv),dpi=600,transparent=True)


# In[31]:


np.nanmean(s['cod']),np.nanmedian(s['cod'])


# In[32]:


np.nanmean(s['ref']),np.nanmedian(s['ref'])


# In[33]:


s['ssa'].shape


# ## Aerosol products

# In[62]:


plt.figure()
plt.plot(s['ssa'][:,2])


# In[63]:


plt.figure(figsize=(6,3))
plt.hist(s['ssa'][:,2],bins=30,range=[np.nanmin(s['ssa'][:,2]),np.nanmax(s['ssa'][:,2])],normed=True,alpha=0.6)
plt.axvline(np.nanmean(s['ssa'][:,2]),ls='-',color='k',label='mean')
plt.axvline(np.nanmedian(s['ssa'][:,2]),ls='--',color='k',label='median')
plt.legend(frameon=False)
plt.xlabel('SSA at 500 nm')
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_SSA_500_hist_{}.png'.format(vv),dpi=600,transparent=True)


# In[64]:


np.nanmean(s['ssa'][:,2]),np.nanmedian(s['ssa'][:,2])


# In[65]:


plt.figure(figsize=(6,3))
plt.hist(s['asy'][:,2],bins=30,range=[np.nanmin(s['asy'][:,2]),np.nanmax(s['asy'][:,2])],normed=True,alpha=0.6,color='orange')
plt.axvline(np.nanmean(s['asy'][:,2]),ls='-',color='k',label='mean')
plt.axvline(np.nanmedian(s['asy'][:,2]),ls='--',color='k',label='median')
plt.legend(frameon=False)
plt.xlabel('ASY at 500 nm')
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_ASY_500_hist_{}.png'.format(vv),dpi=600,transparent=True)


# In[33]:


s['wvl']


# In[66]:


plt.figure()
plt.hist(s['asy'][:,5],range=[0.2,1],bins=30)


# In[67]:


plt.figure()
plt.hist(s['asy'][:,2],bins=30,range=[np.nanmin(s['asy'][:,2]),np.nanmax(s['asy'][:,2])],normed=True,alpha=0.6,color='orange')
plt.axvline(np.nanmean(s['asy'][:,2]),ls='-',color='k',label='mean')
plt.axvline(np.nanmedian(s['asy'][:,2]),ls='--',color='k',label='median')
plt.legend(frameon=False)
plt.xlabel('ASY at 500 nm')
plt.savefig(fp+'plot_DARE/ORACLES2016_ASY_500_hist_{}.png'.format(vv),dpi=600,transparent=True)


# In[68]:


np.nanmean(s['asy'][:,2]),np.nanmedian(s['asy'][:,2])


# In[69]:


s['ext'].shape


# In[70]:


s['wvl']


# In[71]:


plt.figure()
plt.hist(s['ext'][:,2],range=[np.nanmin(s['ext'][:,2]),np.nanmax(s['ext'][:,2])],bins=30)


# ## DARE

# In[97]:


s.keys()


# In[98]:


plt.figure()
plt.hist2d(s['up'][:,2]/s['dn'][:,2],s['sza'],bins=40,range=[[0,1],[0,90]],cmap=plt.cm.jet)
plt.ylabel('SZA [$^\\circ$]')
plt.xlabel('TOA Albedo (BB)')
plt.title('ORACLES 2016 Albedo from DARE calc {}'.format(vv))

cb = plt.colorbar()
cb.set_label('counts')
pu.prelim()

plt.savefig(fp+'plot/ORACLES_2016_albedo_from_calc_vs_SZA_{}.png'.format(vv),dpi=600,transparent=True)


# In[99]:


plt.figure()
plt.hist2d(s['dare'][:,2],s['sza'],bins=40,range=[[-70,150],[0,90]],cmap=plt.cm.Greens)
plt.ylabel('SZA [$^\\circ$]')
plt.xlabel('DARE [W/m$^2$]')
plt.title('ORACLES 2016 DARE TOA from P3 SSFR and 4STAR - DARE calc {}'.format(vv))

cb = plt.colorbar()
cb.set_label('counts')
pu.prelim()

plt.savefig(fp+'plot/ORACLES_2016_DARE_TOA_from_calc_vs_SZA_{}.png'.format(vv),dpi=600,transparent=True)


# In[70]:


plt.figure()
plt.hist2d(s['dare'][:,0],s['sza'],bins=40,range=[[-170,150],[0,90]],cmap=plt.cm.Greens)
plt.ylabel('SZA [$^\\circ$]')
plt.xlabel('DARE [W/m$^2$]')
plt.title('ORACLES 2016 DARE Surf from P3 SSFR and 4STAR - DARE calc {}'.format(vv))

cb = plt.colorbar()
cb.set_label('counts')
pu.prelim()

plt.savefig(fp+'plot/ORACLES_2016_DARE_surf_from_calc_vs_SZA_{}.png'.format(vv),dpi=600,transparent=True)


# In[105]:


plt.figure()
plt.plot(s['doys'],s['dare'][:,0],'.',label='Surface')
plt.plot(s['doys'],s['dare'][:,1],'.',label='Above cloud')
plt.plot(s['doys'],s['dare'][:,2],'.',label='TOA')
plt.axhline(0,ls=':')
plt.legend(frameon=False)
plt.xlabel('DOY')
plt.ylabel('DARE [W/m$^{{2}}$]')
plt.title('ORACLES 2016 DARE calculations {}'.format(vv))
plt.grid()
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_{}.png'.format(vv),dpi=600,transparent=True)


# In[106]:


plt.figure()
plt.hist(s['dare'][:,2],range=[-30,150],bins=50,label='TOA',alpha=0.6)
plt.hist(s['dare'][:,1],range=[-150,0],bins=50,label='Above cloud',alpha=0.6)
plt.hist(s['dare'][:,0],range=[-150,0],bins=50,label='Surface',alpha=0.6)

plt.axvline(np.nanmean(s['dare'][:,2]),ls='-',color='b',label='Mean')
plt.axvline(np.nanmedian(s['dare'][:,2]),ls='--',color='b',label='Median')
plt.axvline(np.nanmean(s['dare'][:,1]),ls='-',color='orange')
plt.axvline(np.nanmedian(s['dare'][:,1]),ls='--',color='orange')
plt.axvline(np.nanmean(s['dare'][:,0]),ls='-',color='g')
plt.axvline(np.nanmedian(s['dare'][:,0]),ls='--',color='g')

pu.prelim()

plt.xlabel('DARE [W/m$^{{2}}$]')
plt.legend(frameon=False)

plt.title('ORACLES 2016 DARE calculations {}'.format(vv))
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_hist_{}.png'.format(vv),dpi=600,transparent=True)


# In[107]:


plt.figure(figsize=(6,4))

plt.hist(s['dare'][:,0],range=[-160,0],bins=50,label='Surface',alpha=0.6,normed=True)
plt.hist(s['dare'][:,1],range=[-160,0],bins=50,label='Above cloud',alpha=0.6,normed=True)
plt.hist(s['dare'][:,2],range=[-30,150],bins=50,label='Top of Atmosphere',alpha=0.6,normed=True)

plt.axvline(np.nanmean(s['dare'][:,2]),ls='-',color='b',label='Mean')
plt.axvline(np.nanmedian(s['dare'][:,2]),ls='--',color='b',label='Median')
plt.axvline(np.nanmean(s['dare'][:,1]),ls='-',color='orange')
plt.axvline(np.nanmedian(s['dare'][:,1]),ls='--',color='orange')
plt.axvline(np.nanmean(s['dare'][:,0]),ls='-',color='g')
plt.axvline(np.nanmedian(s['dare'][:,0]),ls='--',color='g')

#pu.prelim()

plt.xlabel('Instantaneous DARE [W/m$^{{2}}$]')
plt.legend(frameon=False,loc=2)

plt.title('ORACLES August 2016 DARE calculations')
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_normhist_{}.png'.format(vv),dpi=600,transparent=True)


# In[108]:


sa['dare_px'].shape


# In[109]:


plt.figure()
plt.hist(sa['dare_p'][0,:],range=[-75,150],bins=50,label='Parameterisation, PX(Reflec.,AOD,SSA)',alpha=0.6,normed=True)


# In[110]:


plt.figure(figsize=(5,3))
plt.hist(s['dare'][:,2],range=[-75,200],bins=50,label='Calculations from\ndirect measurements',alpha=0.7,normed=True)
plt.hist(sa['dare_p'][0,:],range=[-75,200],bins=50,label='Parameterisation\nP(Reflec.,AOD)',alpha=0.6,normed=True)
plt.hist(sa['dare_px'][0,:],range=[-75,200],bins=50,label='Parameterisation\nPX(Reflec.,AOD,SSA)',alpha=0.4,normed=True)

plt.plot([0,0],[0,0],ls='-',color='k',alpha=0.6,label='Mean')
plt.plot([0,0],[0,0],ls='--',color='k',alpha=0.6,label='Median')
plt.axvline(np.nanmean(s['dare'][:,2]),ls='-',color='b')
plt.axvline(np.nanmedian(s['dare'][:,2]),ls='--',color='b')
plt.axvline(np.nanmean(sa['dare_p'][0,:]),ls='-',color='orange')
plt.axvline(np.nanmedian(sa['dare_p'][0,:]),ls='--',color='orange')
plt.axvline(np.nanmean(sa['dare_px'][0,:]),ls='-',color='g')
plt.axvline(np.nanmedian(sa['dare_px'][0,:]),ls='--',color='g')

plt.text(np.nanmean(s['dare'][:,2])+3,0.014,'{:3.1f}'.format(np.nanmean(s['dare'][:,2])),color='tab:blue',alpha=0.7)
plt.text(np.nanmean(sa['dare_p'][0,:])-25,0.014,'{:3.1f}'.format(np.nanmean(sa['dare_p'][0,:])),color='tab:orange',alpha=0.6)
plt.text(np.nanmean(sa['dare_px'][0,:])-25,0.013,'{:3.1f}'.format(np.nanmean(sa['dare_px'][0,:])),color='tab:green',alpha=0.4)

plt.xlabel('TOA Instantaneous DARE [W/m$^{{2}}$]')
plt.ylabel('Normalized distribution')
#plt.grid()
plt.legend(frameon=False,loc=1)
plt.tight_layout()
plt.savefig(fp+'plot_DARE/DARE_TOA_hist_with_params_PX_{}.png'.format(vv),transparent=True,dpi=600)


# In[111]:


igood = np.isfinite(s['dare'][:,0])


# In[112]:


plt.figure()
plt.hist(s['dare'][igood,2],bins=30,range=[-30,175])


# In[76]:


np.nanmean(s['dare'][igood,2]),np.nanmedian(s['dare'][igood,2]),np.nanstd(s['dare'][igood,2])


# ## Make the boxes

# The observations and model data are aggregated within horizontal domains of at least 2o by 2o indicated in Fig. 2. One of the three main regions encompasses the routine flight track, with individual grid boxes centered at (14oE, 24oS), (12oE, 22oS), (10oE, 20oS), (8oE, 18oS), (6oE, 16oS), (4oE, 14oS), (2oE, 12oS) and (0oE, 10oS). Another more coastal north-south track has the southernmost grid box centered on 22oS, spanning between 9oE and 11.75oE. Seven grid boxes are located every 2 degrees north of this, with the northernmost grid box centered on 8oS. A third, zonal track covers the larger domain of the ER2 measurements, with individual grid boxes spanning latitudinally between 10oS and 6oS and separated longitudinally at two degree intervals beginning at 3oW to the west and 13oE in the east. The box for St. Helena Island spans between 6.72 oW and 4.72 oW, between 16.933 oS and 14.933 oS.

# In[113]:


boxes_diag = []
boxes_ns = []
boxes_ew = []


# In[114]:


boxes_diag_ct = [[14.0,-24.0], [12.0,-22.0],[10.0,-20.0],[8.0,-18.0],[6.0,-16.0],[4.0,-14.0],[2.0,-12.0],[0.0,-10.0]]
boxes_ns_ct = [[10.5,-22.0],[10.5,-20.0],[10.5,-18.0],[10.5,-16.0],[10.5,-14.0],[10.5,-12.0],[10.5,-10.0],[10.5,-8.0]]
boxes_ew_ct = [[-3.0,-8.0],[-1.0,-8.0],[1.0,-8.0],[3.0,-8.0],[5.0,-8.0],[7.0,-8.0],[9.0,-8.0],[11.0,-8.0],[13.0,-8.0]]


# Corners are [x0,x1,y0,y1]

# In[115]:


boxes_ns = [[9.0,11.75,i[1]-1.0,i[1]+1.0] for i in boxes_ns_ct]


# In[116]:


boxes_ew = [[-10.0,-6.0,i[0]-1.0,i[0]+1.0] for i in boxes_ew_ct]


# In[117]:


boxes_diag = [[i[0]-1.0,i[0]+1,i[1]-1.0,i[1]+1.0] for i in boxes_diag_ct]


# In[118]:


boxes_diag


# In[119]:


boxes_ew


# In[120]:


boxes_ns


# In[121]:


bins_diag = []
for i,b in enumerate(boxes_diag):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) & (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare'][igood,2]))
    bins_diag.append(s['dare'][igood,2][ia])


# In[122]:


bins_diag[0] = bins_diag[1][0:5]


# In[123]:


bins_ns,bins_ns_num = [],[]
for i,b in enumerate(boxes_ns):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) & (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare'][igood,2]))
    bins_ns.append(s['dare'][igood,2][ia])
    bins_ns_num.append(len(s['dare'][igood,2][ia]))


# In[124]:


bins_ns[-1] = bins_ns[-2][0:5]


# In[125]:


bins_ew = []
for i,b in enumerate(boxes_ew):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) & (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare'][igood,2]))
    bins_ew.append(s['dare'][igood,2][ia])


# In[128]:


len(boxes_diag),len(bins_diag)


# In[129]:


[fig,ax] = plt.subplots(1,8,figsize=(13,3))

for i,b in enumerate(boxes_diag_ct):
    ax[i].hist(bins_diag[i],bins=30,edgecolor='None',alpha=0.7,color='g',range=(-30,150),
               zorder=10,normed=True,orientation='horizontal',label='Calc')
    ax[i].hist(sa['bins_diag'][0,i][0],bins=30,edgecolor='None',alpha=0.3,color='m',range=(-30,150),
               zorder=-1,normed=True,orientation='horizontal',label='Param')
    ax[i].axhline(np.nanmean(bins_diag[i]),color='g')#,label='mean')
    ax[i].axhline(np.nanmedian(bins_diag[i]),color='g',linestyle='--')#,label='median')
    xmin, xmax = ax[i].get_xlim()
    ax[i].set_xticks(np.round(np.linspace(xmin, xmax, 3), 2))
    if i>0:
        [ag.set_visible(False) for ag in ax[i].yaxis.get_ticklabels()]
    if i==0:
        ax[i].legend(frameon=False,bbox_to_anchor=[0.08, 0.85])
    ax[i].set_title('{}$^\\circ$ ,{}$^\\circ$'.format(b[0],b[1]))
    ax[i].axhline(0,ls=':',color='k',alpha=0.2)
    #if i%2: pu.prelim(ax[i])
ax[0].set_ylabel('DARE [W/m$^2$]')
fig.suptitle('ORACLES 2016 Routine Diagonal (Lon,Lat) - TOA {}'.format(vv))
fig.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_TOA_calc_diag_boxes_{}.png'.format(vv),dpi=600,transparent=True)


# In[94]:


[fig,ax] = plt.subplots(1,8,figsize=(13,3))

for i,b in enumerate(boxes_ns_ct):
    ax[i].hist(bins_ns[i],bins=30,edgecolor='None',alpha=0.7,color='b',range=(-50,150),zorder=10,normed=True,orientation='horizontal',label='Calc')
    ax[i].hist(sa['bins_ns'][0,i][0],bins=30,edgecolor='None',alpha=0.2,color='orange',range=(-50,150),zorder=-1,normed=True,orientation='horizontal',label='Param')
    ax[i].axhline(np.nanmean(bins_ns[i]),color='b')#,label='mean')
    ax[i].axhline(np.nanmedian(bins_ns[i]),color='b')#,linestyle='--',label='median')
    xmin, xmax = ax[i].get_xlim()
    ax[i].set_xticks(np.round(np.linspace(xmin, xmax, 3), 2))
    if i>0:
        [ag.set_visible(False) for ag in ax[i].yaxis.get_ticklabels()]
    else:
        ax[i].legend(frameon=False,bbox_to_anchor=[0.08,0.85])
    ax[i].set_title('{}$^\\circ$ ,{}$^\\circ$'.format(b[0],b[1]))
    ax[i].axhline(0,ls=':',color='k',alpha=0.2)
    #if i%2: pu.prelim(ax[i])
ax[0].set_ylabel('DARE [W/m$^2$]')
fig.suptitle('ORACLES 2016 North-South (Lon,Lat) - TOA {}'.format(vv))
fig.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_TOA_calc_ns_boxes_{}.png'.format(vv),dpi=600,transparent=True)


# In[93]:


plt.figure()
sca = plt.scatter(s['lon'][igood],s['lat'][igood],c=s['dare'][igood,2],edgecolor='None',s=40,alpha=0.5,cmap=plt.cm.viridis)
plt.grid()
plt.xlim(-1,16)
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('Latitude [$^\\circ$]')
cb = plt.colorbar(sca)
cb.set_label('DARE [W/m$^2$]')
pu.prelim()

for i,b in enumerate(boxes_ns): 
    plt.plot([b[0],b[0],b[1],b[1],b[0]],[b[2],b[3],b[3],b[2],b[2]],'-b')
for i,b in enumerate(boxes_diag): 
    plt.plot([b[0],b[0],b[1],b[1],b[0]],[b[2],b[3],b[3],b[2],b[2]],'-g')

plt.ylim(-25,-7)
plt.title('ORACLES 2016 DARE TOA from calculations of\n4STAR AOD, skyscans and SSFR cloud retrievals')
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_TOA_calc_map_param_{}.png'.format(vv),dpi=600,transparent=True)


# ## Put into box-whisker plots

# In[132]:


ns_ct = np.array(boxes_ns_ct)


# In[133]:


ns_ct[:,1]


# In[134]:


def color_boxes(bp,color):
    for b in bp['boxes']:
        b.set_facecolor(color)
        b.set_edgecolor(color)
        b.set_alpha(0.4)
    for b in bp['means']:
        b.set_marker('o')
        b.set_color('firebrick')
        b.set_markerfacecolor(color)
        b.set_markeredgecolor(color)
        b.set_alpha(0.4)
    for b in bp['whiskers']:
        b.set_linestyle('-')
        b.set_color(color)
        b.set_alpha(0.4)
    for b in bp['caps']:
        b.set_alpha(0.4)
        b.set_color(color)
    for b in bp['medians']:
        b.set_linewidth(4)
        b.set_color('gold')
        b.set_alpha(0.4)
    for b in bp['fliers']:
        b.set_marker('.')
        b.set_markeredgecolor('None')
        b.set_markerfacecolor('lightgrey')
        b.set_alpha(0.5)


# In[136]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('TOA Instantaneous DARE [W/m$^2$]')
plt.ylabel('Latitude [$^{{\\circ}}$S]')

#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
color_boxes(bp,'tab:blue')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(dare_means[:,0],dare_means[:,1],'-',color='tab:blue',alpha=0.6,label='Means')
plt.xlim(-50,200)
#plt.ylim(0,2500)

for j,nn in enumerate(bins_ns_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1],'{:2.0f}'.format(nn),
                 color='tab:blue',fontsize=7,verticalalignment='center',horizontalalignment='left')


plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
           ['Mean','Median','25% - 75%','min-max','outliers'],
           frameon=False,loc=1,numpoints=1)
#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_boxplot_NS_lat_{}.png'.format(vv),dpi=600,transparent=True)


# In[250]:


sa.keys()


# In[137]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns,positions=ns_ct[:,1]-0.5,vert=False,showfliers=True,widths=0.5,showmeans=True,patch_artist=True)

bpp = plt.boxplot(sa['bins_ns_p'][0],positions=ns_ct[:,1]+0.5,vert=False,showfliers=True,widths=0.5,showmeans=True,patch_artist=True)
bpx = plt.boxplot(sa['bins_ns_px'][0],positions=ns_ct[:,1],vert=False,showfliers=True,widths=0.5,showmeans=True,patch_artist=True)


plt.xlabel('TOA Instantaneous DARE [W/m$^2$]')
plt.ylabel('Latitude [$^{{\\circ}}$S]')

color_boxes(bp,'tab:blue')
color_boxes(bpp,'tab:orange')
color_boxes(bpx,'tab:green')

#dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
#plt.plot(dare_means[:,0],dare_means[:,1],'-',color='tab:blue',alpha=0.6,label='Calc.')
#dare_means_p = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bpp['means']])
#plt.plot(dare_means_p[:,0],dare_means_p[:,1],'-',color='tab:orange',alpha=0.6,label='P(R,AOD)')
#dare_means_px = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bpx['means']])
#plt.plot(dare_means_px[:,0],dare_means_px[:,1],'-',color='tab:green',alpha=0.6,label='PX(R,AOD,SSA)')

for j,nn in enumerate(bins_ns_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1]-0.5,'{:2.0f}'.format(nn),
                 color='tab:blue',fontsize=7,verticalalignment='center',horizontalalignment='left')
for j,nn in enumerate(sa['bins_ns_p'][0]):
    if len(nn)>0:
        plt.text(min(bpp['means'][j].get_data()[0])+5,ns_ct[j,1]+0.5,'{:2.0f}'.format(len(nn[0])),
         color='tab:orange',fontsize=7,verticalalignment='center',horizontalalignment='left')
for j,nn in enumerate(sa['bins_ns_px'][0]):
    if len(nn)>0:
        plt.text(min(bpx['means'][j].get_data()[0])+5,ns_ct[j,1],'{:2.0f}'.format(len(nn[0])),
         color='tab:green',fontsize=7,verticalalignment='center',horizontalalignment='left')

plt.xlim(-50,200)
plt.legend([bpp['boxes'][0],bpx['boxes'][0],bp['boxes'][0]],['P(R,AOD)','PX(R,AOD,SSA)','Calc.'],frameon=False,loc=4,numpoints=1)
#plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
#           ['Mean','Median','25% - 75%','min-max','outliers'],
#           frameon=False,loc=1,numpoints=1)
#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_boxplot_NS_lat_withParams_{}.png'.format(vv),dpi=600,transparent=True)


# In[138]:


sa['bins_ns_px'][0]


# ## Make the box-whisker plots for 24h DARE

# In[139]:


np.nanmean(s['dare_24h'][:,2]),np.nanstd(s['dare_24h'][:,2])


# In[140]:


bins_ns_24,bins_ns_24_num = [],[]
for i,b in enumerate(boxes_ns):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) & (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare_24h'][igood,2]))
    bins_ns_24.append(s['dare_24h'][igood,2][ia])
    bins_ns_24_num.append(len(s['dare_24h'][igood,2][ia]))


# In[141]:


bins_ns_24h,bins_ns_24h_num = [],[]
for i,b in enumerate(boxes_ns):
    ia = (sh['lon']>= b[0]) & (sh['lon']<=b[1]) &(sh['lat']>=b[2]) &    (sh['lat']<=b[3]) & (np.isfinite(sh['dare_avg'][:,2]))
    bins_ns_24h.append(sh['dare_avg'][:,2][ia])
    bins_ns_24h_num.append(len(sh['dare_avg'][:,2][ia]))


# In[142]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_24,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('TOA Diurnally averaged DARE [W/m$^2$]')
plt.ylabel('Latitude [$^{{\\circ}}$S]')

#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
color_boxes(bp,'tab:blue')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(dare_means[:,0],dare_means[:,1],'-',color='tab:blue',alpha=0.6,label='Means')
plt.xlim(-25,75)
#plt.ylim(0,2500)

for j,nn in enumerate(bins_ns_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1],'{:2.0f}'.format(nn),
                 color='tab:blue',fontsize=7,verticalalignment='center',horizontalalignment='left')
pu.prelim()

plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
           ['Mean','Median','25% - 75%','min-max','outliers'],
           frameon=False,loc=1,numpoints=1)
#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_24h_boxplot_NS_lat_{}.png'.format(vv),dpi=600,transparent=True)


# In[487]:


np.nanmean(sh['dare_avg'][:,2]),np.nanstd(sh['dare_avg'][:,2])


# In[143]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_24h,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('TOA Diurnally averaged DARE [W/m$^2$]')
plt.ylabel('Latitude [$^{{\\circ}}$S]')

#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
color_boxes(bp,'tab:blue')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(dare_means[:,0],dare_means[:,1],'-',color='tab:blue',alpha=0.6,label='Means')
plt.xlim(-25,60)
#plt.ylim(0,2500)

for j,nn in enumerate(bins_ns_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1],'{:2.0f}'.format(nn),
                 color='tab:blue',fontsize=7,verticalalignment='center',horizontalalignment='left')
#pu.prelim()
plt.axvline(0,ls='--',lw=1,color='k',alpha=0.6)

plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
           ['Mean','Median','25% - 75%','min-max','outliers'],
           frameon=False,loc=1,numpoints=1)
#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_24h_from_calc_boxplot_NS_lat_{}.png'.format(vv),dpi=600,transparent=True)


# In[145]:


dare_means


# In[146]:


np.arange(0,-22,-2)


# In[147]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_24h,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('TOA Diurnally averaged DARE [W/m$^2$]')
plt.ylabel('Latitude [$^{{\\circ}}$S]')

#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
color_boxes(bp,'tab:blue')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(dare_means[:,0],dare_means[:,1],'-',color='tab:blue',alpha=0.6,label='Means')
la = plt.plot(AL_dare_2016, AL_Lat,'s-',color='tab:red',alpha=0.6,label='ALADIN Model')
plt.ylim(-23,1)
plt.yticks(np.arange(-22,2,2))
plt.gca().set_yticklabels(np.arange(-22,2,2))
#plt.yticks([-22.0-20.0,-18.0,-16.0,-14.0,-12.0,-10.0,-8.0,-6.0,-4.0,-2.0])
plt.xlim(-25,60)
#plt.ylim(0,2500)

for j,nn in enumerate(bins_ns_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1],'{:2.0f}'.format(nn),
                 color='tab:blue',fontsize=7,verticalalignment='center',horizontalalignment='left')
#pu.prelim()
plt.axvline(0,ls='--',lw=1,color='k',alpha=0.6)

plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0],la[0]],
           ['Mean','Median','25% - 75%','min-max','outliers','ALADIN\n(August 2016)'],
           frameon=False,loc=1,numpoints=1)
#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_24h_from_calc_with_ALADIN_boxplot_NS_lat_{}.png'.format(vv),dpi=600,transparent=True)


# In[299]:


sunf = lambda x: np.cos(np.deg2rad(x)) if np.cos(np.deg2rad(x))>=0 else 0.0


# In[301]:


u = np.cos(np.radians(np.arange(-180,180)))


# In[302]:


u[u<0] = 0.0


# In[304]:


u/np.mean(u)


# In[306]:


plt.figure()
plt.plot(u,u/np.mean(u))


# In[309]:


ud = np.array([u,u/np.mean(u)])


# In[310]:


ud.shape


# In[313]:


i_fac_s = [np.argmin(abs(ud[0,:]-mu)) for mu in s['mu']]


# In[316]:


s['dare_24hs'] = s['dare']+np.nan
s['dare_24hs'][:,0] = s['dare'][:,0]/ud[1,i_fac_s]
s['dare_24hs'][:,1] = s['dare'][:,1]/ud[1,i_fac_s]
s['dare_24hs'][:,2] = s['dare'][:,2]/ud[1,i_fac_s]


# In[324]:


bins_ns_24,bins_ns_24_num = [],[]
for i,b in enumerate(boxes_ns):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) &           (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare_24hs'][igood,2]))
    bins_ns_24.append(s['dare_24hs'][igood,2][ia])
    bins_ns_24_num.append(len(s['dare_24hs'][igood,2][ia]))


# In[318]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_24,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('TOA Diurnally averaged DARE [W/m$^2$]')
plt.ylabel('Latitude [$^{{\\circ}}$S]')

#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
color_boxes(bp,'tab:blue')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(dare_means[:,0],dare_means[:,1],'-',color='tab:blue',alpha=0.6,label='Means')
plt.xlim(-25,75)
#plt.ylim(0,2500)

for j,nn in enumerate(bins_ns_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1],'{:2.0f}'.format(nn),
                 color='tab:blue',fontsize=7,verticalalignment='center',horizontalalignment='left')


plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
           ['Mean','Median','25% - 75%','min-max','outliers'],
           frameon=False,loc=1,numpoints=1)
#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_24h_boxplot_NS_lat.png',dpi=600,transparent=True)


# In[561]:


plt.figure()
plt.plot(dare_means[:,1],dare_means[:,0],'-',color='tab:blue',alpha=0.6,label='Diurnal Averaged Calc.')
plt.plot(AL_Lat,AL_dare_2016,'-',color='tab:red',alpha=0.6,label='ALADIN Diurnal Averages')


# ## Make box plot for 24h DARE, with clouds changing

# In[ ]:


bins_ns_24,bins_ns_24_num = [],[]
for i,b in enumerate(boxes_ns):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) & (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare_24h'][igood,2]))
    bins_ns_24.append(s['dare_24h'][igood,2][ia])
    bins_ns_24_num.append(len(s['dare_24h'][igood,2][ia]))


# In[ ]:


bins_ns_24h,bins_ns_24h_num = [],[]
for i,b in enumerate(boxes_ns):
    ia = (sh['lon']>= b[0]) & (sh['lon']<=b[1]) &(sh['lat']>=b[2]) &    (sh['lat']<=b[3]) & (np.isfinite(sh['dare_avg'][:,2]))
    bins_ns_24h.append(sh['dare_avg'][:,2][ia])
    bins_ns_24h_num.append(len(sh['dare_avg'][:,2][ia]))


# In[149]:


binsd_ns_24h,binsd_ns_24h_num = [],[]
for i,b in enumerate(boxes_ns):
    ia = (shd['lon']>= b[0]) & (shd['lon']<=b[1]) &(shd['lat']>=b[2]) &    (shd['lat']<=b[3]) & (np.isfinite(shd['dare_avg'][:,2]))
    binsd_ns_24h.append(shd['dare_avg'][:,2][ia])
    binsd_ns_24h_num.append(len(shd['dare_avg'][:,2][ia]))


# In[180]:


plt.figure(figsize=(5,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_24h,positions=ns_ct[:,1]-0.5,vert=False,showfliers=True,widths=0.5,showmeans=True,patch_artist=True)
bpd = plt.boxplot(binsd_ns_24h,positions=ns_ct[:,1],vert=False,showfliers=True,widths=0.5,showmeans=True,patch_artist=True)
plt.xlabel('TOA Diurnally averaged DARE [W/m$^2$]')
plt.ylabel('Latitude [$^{{\\circ}}$S]')

#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
color_boxes(bp,'tab:blue')
color_boxes(bpd,'tab:orange')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(dare_means[:,0],dare_means[:,1],'-',color='tab:blue',alpha=0.6,label='Means')
dare_meansd = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bpd['means']])
plt.plot(dare_meansd[:,0],dare_meansd[:,1],'-',color='tab:orange',alpha=0.6,label='Means')
plt.xlim(-15,70)
#plt.ylim(0,2500)

for j,nn in enumerate(bins_ns_24h_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1]-0.5,'{:2.0f}'.format(nn),
                 color='tab:blue',fontsize=7,verticalalignment='center',horizontalalignment='left')
for j,nn in enumerate(binsd_ns_24h_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1],'{:2.0f}'.format(nn),
                 color='tab:orange',fontsize=7,verticalalignment='center',horizontalalignment='left')
#pu.prelim()
plt.axvline(0,ls='--',lw=1,color='k',alpha=0.6)

plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0],bp['boxes'][0],bpd['boxes'][0]],
           ['Mean','Median','25% - 75%','min-max','outliers','constant','sine'],
           frameon=True,loc=1,numpoints=1,bbox_to_anchor=[1.4,0.7])
#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_24h_from_calc_boxplot_NS_lat_withcl_{}.png'.format(vv),dpi=600,transparent=True)


# In[394]:


plt.figure(figsize=(3.55,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_24h,positions=ns_ct[:,1]-0.5,vert=False,showfliers=True,widths=0.5,showmeans=True,patch_artist=True)
bpd = plt.boxplot(binsd_ns_24h,positions=ns_ct[:,1],vert=False,showfliers=True,widths=0.5,showmeans=True,patch_artist=True)
plt.xlabel('TOA Diurnally averaged DARE [W/m$^2$]')
plt.ylabel('Latitude [$^{{\\circ}}$S]')

#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
color_boxes(bp,'tab:blue')
color_boxes(bpd,'tab:orange')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(dare_means[:,0],dare_means[:,1],'-',color='tab:blue',alpha=0.6,label='Means')
dare_meansd = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bpd['means']])
plt.plot(dare_meansd[:,0],dare_meansd[:,1],'-',color='tab:orange',alpha=0.6,label='Means')
plt.xlim(-15,70)
#plt.ylim(0,2500)

for j,nn in enumerate(bins_ns_24h_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1]-0.5,'{:2.0f}'.format(nn),
                 color='tab:blue',fontsize=7,verticalalignment='center',horizontalalignment='left')
for j,nn in enumerate(binsd_ns_24h_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,ns_ct[j,1],'{:2.0f}'.format(nn),
                 color='tab:orange',fontsize=7,verticalalignment='center',horizontalalignment='left')
#pu.prelim()
plt.axvline(0,ls='--',lw=1,color='k',alpha=0.6)

#plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0],bp['boxes'][0],bpd['boxes'][0]],
#           ['Mean','Median','25% - 75%','min-max','outliers','constant','sine'],
#           frameon=True,loc=1,numpoints=1,bbox_to_anchor=[1.4,0.7])
#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_24h_from_calc_boxplot_NS_lat_withcl_noleg_{}.png'.format(vv),dpi=600,transparent=True)


# In[171]:


sh['dare_avg'][sh['dare_avg'][:,2]==0,2] = np.nan
shd['dare_avg'][shd['dare_avg'][:,2]==0,2] = np.nan


# In[198]:


plt.figure(figsize=(4,3))
plt.hist(sh['dare_avg'][:,2],color='tab:blue',bins=30,range=[-15,75],alpha=0.5,label='Contant')
plt.hist(shd['dare_avg'][:,2],color='tab:orange',bins=30,range=[-15,75],alpha=0.5,label='Sine')

plt.axvline(0,ls=':',c='grey',zorder=-10)
plt.xlabel('TOA Diurnally averaged DARE [W/m$^2$]')


plt.axvline(np.nanmean(sh['dare_avg'][:,2]),ls='-',color='tab:blue',label='Mean')
plt.axvline(np.nanmedian(sh['dare_avg'][:,2]),ls='--',color='tab:blue',label='Median')
plt.text(np.nanmean(sh['dare_avg'][:,2])+3,2000,'{:3.1f}'.format(np.nanmean(sh['dare_avg'][:,2])),
         color='tab:blue',alpha=0.7)

plt.axvline(np.nanmean(shd['dare_avg'][:,2]),ls='-',color='tab:orange')
plt.axvline(np.nanmedian(shd['dare_avg'][:,2]),ls='--',color='tab:orange')
plt.text(np.nanmean(shd['dare_avg'][:,2])-12,1800,'{:3.1f}'.format(np.nanmean(shd['dare_avg'][:,2])),
         color='tab:orange',alpha=0.7)

plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_24h_hist_{}.png'.format(vv),dpi=600,transparent=True)


# In[182]:


np.nanmean(sh['dare_avg'][:,2]),np.nanmean(shd['dare_avg'][:,2]),np.nanstd(sh['dare_avg'][:,2]),np.nanstd(shd['dare_avg'][:,2])


# ## Check for trends in underlying properties

# In[381]:


s.keys()


# In[382]:


s['wvl']


# In[383]:


bins_ns_aod,bins_ns_aodn = [],[]
bins_ns_cod,bins_ns_codn = [],[]
bins_ns_clear,bins_ns_clearn = [],[]
bins_ns_ssa,bins_ns_ssan = [],[]
bins_ns_asy,bins_ns_asyn = [],[]
bins_ns_ref,bins_ns_refn = [],[]
for i,b in enumerate(boxes_ns):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) &           (s['lat'][igood]<=b[3]) & (np.isfinite(s['ext'][igood,2]))
    bins_ns_aod.append(s['ext'][igood,2][ia]*3.0)
    bins_ns_aodn.append(len(s['ext'][igood,2][ia]))
    ib = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) &       (s['lat'][igood]<=b[3]) & (np.isfinite(s['cod'][igood]))
    bins_ns_cod.append(s['cod'][igood][ib])
    bins_ns_codn.append(len(s['cod'][igood][ib]))
    ic = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) &       (s['lat'][igood]<=b[3]) & (np.isfinite(s['ssa'][igood,2]))
    bins_ns_ssa.append(s['ssa'][igood,2][ic])
    bins_ns_ssan.append(len(s['ssa'][igood,2][ic]))
    icc = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) &           (s['lat'][igood]<=b[3]) & (np.isfinite(s['asy'][igood,2]))
    bins_ns_asy.append(s['asy'][igood,2][icc])
    bins_ns_asyn.append(len(s['asy'][igood,2][icc]))
    ie = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) &           (s['lat'][igood]<=b[3]) & (np.isfinite(s['ref'][igood]))
    bins_ns_ref.append(s['ref'][igood][ie])
    bins_ns_refn.append(len(s['ref'][igood][ie]))
    iaf = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) &           (s['lat'][igood]<=b[3]) & (np.isfinite(s['ref'][igood]))
    bins_ns_ref.append(s['ref'][igood][iaf])
    bins_ns_refn.append(len(s['ref'][igood][iaf]))


# In[384]:


#plt.figure(figsize=(4,4))
fig,ax = plt.subplots(1,2)
ax=ax.ravel()
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = ax[0].boxplot(bins_ns_aod,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
ax[0].set_xlabel('AOD')
ax[0].set_ylabel('Latitude [$^{{\\circ}}$S]')

color_boxes(bp,'red')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
ax[0].plot(dare_means[:,0],dare_means[:,1],'-',color='red',alpha=0.6,label='Means')

ax[0].legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
           ['Mean','Median','25% - 75%','min-max','outliers'],
           frameon=False,loc=1,numpoints=1)

bc = ax[1].boxplot(bins_ns_cod,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
ax[1].set_xlabel('COD')
color_boxes(bc,'m')
cod_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bc['means']])
ax[1].plot(cod_means[:,0],cod_means[:,1],'-',color='m',alpha=0.6,label='Means')


#plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
#plt.savefig(fp+'plot_DARE/ORACLES2016_DARE_24h_boxplot_NS_lat.png',dpi=600,transparent=True)


# In[395]:


#plt.figure(figsize=(4,4))
fig,ax = plt.subplots(1,1,figsize=(3.6,4))
#ax=ax.ravel()
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = ax.boxplot(bins_ns_aod,positions=ns_ct[:,1]-0.1,vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
ax.set_xlabel('AOD',color='red')
ax.set_ylabel('Latitude [$^{{\\circ}}$S]')
ax.set_xlim(0,0.8)
color_boxes(bp,'red')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
ax.plot(dare_means[:,0],dare_means[:,1],'-',color='red',alpha=0.6,label='Means')

#ax[0].legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
#           ['Mean','Median','25% - 75%','min-max','outliers'],
#           frameon=False,loc=1,numpoints=1)
ax.tick_params(axis='x', labelcolor='red')

ax1 = ax.twiny()

bc = ax1.boxplot(bins_ns_cod,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
ax1.set_xlabel('COD',color='m')
color_boxes(bc,'m')
ax1.set_xlim(0,50)
cod_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bc['means']])
ax1.plot(cod_means[:,0],cod_means[:,1],'-',color='m',alpha=0.6,label='Means')
ax1.tick_params(axis='x', labelcolor='m')
ax1.legend([bc['boxes'][0],bp['boxes'][0]],['Cloud','Aerosol'],frameon=False,loc=4)

plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES_2016_COD_AOD_double_box_lat_{}.png'.format(vv),dpi=600,transparent=True)


# In[386]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_cod,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('COD')


# In[340]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_ref,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('REF')


# In[338]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_ssa,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('SSA')


# In[339]:


plt.figure(figsize=(4,4))
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = plt.boxplot(bins_ns_asy,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
plt.xlabel('ASY')


# In[396]:


fig,ax = plt.subplots(1,1,figsize=(3.6,4))
#ax=ax.ravel()
#plt.plot(s['dare'][igood,2],s['lat'][igood],'.',alpha=0.05)
bp = ax.boxplot(bins_ns_ssa,positions=ns_ct[:,1]-0.1,vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
ax.set_xlabel('SSA',color='blue')
ax.set_ylabel('Latitude [$^{{\\circ}}$S]')
#ax.set_xlim(0,0.8)
color_boxes(bp,'blue')

dare_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
ax.plot(dare_means[:,0],dare_means[:,1],'-',color='blue',alpha=0.6,label='Means')

#ax[0].legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
#           ['Mean','Median','25% - 75%','min-max','outliers'],
#           frameon=False,loc=1,numpoints=1)
ax.tick_params(axis='x', labelcolor='blue')

ax1 = ax.twiny()

bc = ax1.boxplot(bins_ns_asy,positions=ns_ct[:,1],vert=False,showfliers=True,widths=1,showmeans=True,patch_artist=True)
ax1.set_xlabel('ASY',color='y')
color_boxes(bc,'y')
#ax1.set_xlim(0,50)
cod_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bc['means']])
ax1.plot(cod_means[:,0],cod_means[:,1],'-',color='y',alpha=0.6,label='Means')
ax1.tick_params(axis='x', labelcolor='y')
ax1.legend([bc['boxes'][0],bp['boxes'][0]],['ASY','SSA'],frameon=False,loc=4)

plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES_2016_SSA_ASY_double_box_NS_lat_{}.png'.format(vv),dpi=600,transparent=True)


# ## Cloud optical depth and ref diurnal variations

# In[199]:


shd.keys()


# In[211]:


ishd = np.isfinite(shd['cods'][:,24])


# In[261]:


mod = sio.loadmat(fp+'data_other/MODIS/MODIS_ORACLES2016_match.mat')


# In[259]:


mod['mod']['time'][0,0][0,:].shape


# In[253]:


fig,ax = plt.subplots(3,1,figsize=(6,4),sharex=True)
utcs = np.arange(0,24,0.5)
bpd = ax[0].boxplot(shd['cods'][ishd,:],positions=utcs,showmeans=True,patch_artist=True,showfliers=False,widths=0.4)
color_boxes(bpd,'lightgrey')
bpd_day = ax[0].boxplot(shd['cods'][ishd,14:38],positions=utcs[14:38],showmeans=True,patch_artist=True,showfliers=False,widths=0.4)
color_boxes(bpd_day,'m')
ax[0].set_ylabel('COD')
ax[0].set_xlim([4,22])
ax[0].set_xticks(np.arange(0,24,2))
ax[0].set_xticklabels(np.arange(0,24,2))

bpd = ax[1].boxplot(shd['refs'][ishd,:],positions=utcs,showmeans=True,patch_artist=True,showfliers=False,widths=0.4)
color_boxes(bpd,'lightgrey')
bpd_day = ax[1].boxplot(shd['refs'][ishd,14:38],positions=utcs[14:38],showmeans=True,patch_artist=True,showfliers=False,widths=0.4)
color_boxes(bpd_day,'g')
ax[1].set_ylabel('R$_{{eff}}$ [$\mu$m]')
#ax[1].set_xlabel('UTC [h]')
ax[1].set_xlim([4,22])
ax[1].set_xticks(np.arange(0,24,2))
ax[1].set_xticklabels(np.arange(0,24,2))

ax[2].hist(shd['utc'],bins=48,alpha=0.5,range=[0,24],label='P3',normed=False,color='k')
ax[2].hist(mod['mod']['time'][0,0][0,:],bins=48,alpha=0.5,range=[0,24],label='TERRA',normed=False,color='green')
ax[2].hist(mod['myd']['time'][0,0][0,:],bins=48,alpha=0.5,range=[0,24],label='AQUA',normed=False,color='tab:blue')
ax[2].set_ylabel('# of Samples')

ax[2].legend(frameon=False)
ax[2].set_xlabel('UTC [h]')

ax[2].set_xticks(np.arange(0,24,2))
ax[2].set_xticklabels(np.arange(0,24,2))
ax[2].set_xlim([4,22])

plt.savefig(fp+'plot_DARE/ORACLES_2016_Cloud_diurnal_var.png',dpi=600,transparent=True)


# In[262]:


mod['diurn'].dtype.names


# ## Plot cloud diurnal cycle amplitude

# In[275]:


plt.figure()
plt.scatter(mod['mod']['lon'][0,0][0,mod['diurn']['iaes'][0,0][0,:]],
            mod['mod']['lat'][0,0][0,mod['diurn']['iaes'][0,0][0,:]],
            40.0,
            c=mod['diurn']['A_cod'][0,0][0,:])


# In[276]:


plt.figure()
plt.hist(mod['diurn']['A_cod'][0,0][0,:])


# In[277]:


import scipy.stats as stats


# In[341]:


NX, NY = 50, 30
statistic, xedges, yedges, binnumber = stats.binned_statistic_2d(
    mod['mod']['lon'][0,0][0,mod['diurn']['iaes'][0,0][0,:]],
    mod['mod']['lat'][0,0][0,mod['diurn']['iaes'][0,0][0,:]], 
    values=mod['diurn']['A_cod'][0,0][0,:], statistic='mean', 
    bins=[np.arange(-11.0,15,0.5), np.arange(-23.0,-7.0,0.5)])


# In[319]:


yedges,xedges,len(xedges),len(yedges)


# In[320]:


xmean = [(x+xedges[i+1])/2.0 for i,x in enumerate(xedges[:-1])]
ymean = [(y+yedges[i+1])/2.0 for i,y in enumerate(yedges[:-1])]


# In[363]:


plt.figure(figsize=(4,3.5))
im = plt.imshow(statistic.T,extent=[xmean[0],xmean[-1],ymean[-1],ymean[0]],origin='low',cmap='plasma')
plt.xlabel('Longitude [$^{{\circ}}$]')
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.xlim([-1,15])
plt.colorbar(im,label='COD sine amplitude',fraction=0.046, pad=0.04)
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES_2016_COD_amplitude.png',dpi=600,transparent=True)


# In[353]:


statisticr, xedgesr, yedgesr, binnumber = stats.binned_statistic_2d(
    mod['mod']['lon'][0,0][0,mod['diurn']['iaes'][0,0][0,:]],
    mod['mod']['lat'][0,0][0,mod['diurn']['iaes'][0,0][0,:]], 
    values=mod['diurn']['A_ref'][0,0][0,:], statistic='mean', 
    bins=[np.arange(-11.0,15,0.5), np.arange(-23.0,-7.0,0.5)])


# In[354]:


xmeanr = [(x+xedgesr[i+1])/2.0 for i,x in enumerate(xedgesr[:-1])]
ymeanr = [(y+yedgesr[i+1])/2.0 for i,y in enumerate(yedgesr[:-1])]


# In[366]:


plt.figure(figsize=(4,3.5))
im = plt.imshow(statisticr.T,extent=[xmean[0],xmean[-1],ymean[-1],ymean[0]],origin='low')
plt.xlabel('Longitude [$^{{\circ}}$]')
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.xlim([-1,15])
plt.colorbar(im,label='R$_{{eff}}$ sine amplitude',fraction=0.046, pad=0.04)
plt.tight_layout()
plt.savefig(fp+'plot_DARE/ORACLES_2016_REF_amplitude.png',dpi=600,transparent=True)


# ## Ratio of DARE 24h average to instantaneous

# In[379]:


# need to ratio instantaneous dare from 's' to 24h averages from 'sh'


# In[367]:


utcx = np.arange(0,24,0.5)


# In[368]:


s['sza'].shape


# In[374]:


plt.figure()
plt.plot(s['mu'],s['dare'][:,2]/sh['dare_avg'][:,2],'.')
plt.ylim(-150,150)


# In[557]:


ig = np.where(np.isfinite(sh['dare_avg'][:,2]) & np.isfinite(s['dare'][:,2]) &              (s['dare'][:,2]!=0.0) & (sh['dare_avg'][:,2] != 0.0) & (abs(sh['dare_avg'][:,2]) > 0.25) &               (abs(s['dare'][:,2])>0.25))[0]


# In[375]:


len(ig)


# In[376]:


plt.figure()
plt.plot(s['dare'][ig,2]/sh['dare_avg'][ig,2],'.')


# In[378]:


ii = 14610


# In[379]:


s['dare'][ig[ii],2]


# In[530]:


sh['dare_avg'][ig[ii],2]


# In[ ]:


plt.figure()
plt.plot(s['dare'][ig[ii],2])


# In[534]:


i2d.shape


# In[380]:


plt.figure(figsize=(8,3))
plt.scatter(s['mu'][ig],s['dare'][ig,2]/sh['dare_avg'][ig,2],c=s['cod'][ig],vmin=0,vmax=40,label='ORACLES 2016 calculations')
plt.plot(i2d[:,0],i2d[:,1],label='Redemann et al. 2006',c='tab:orange')
plt.ylim(0,5)

plt.ylabel('DARE$_{{instant}}$ / DARE$_{{diurnal average}}$')
plt.xlabel('$\mu_0$ (cos(SZA))')
plt.legend()
plt.colorbar(label='Cloud Optical Depth',extend='max')
plt.tight_layout()
#plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_ratio_instant_to_24h.png',dpi=600,transparent=True)


# In[ ]:




