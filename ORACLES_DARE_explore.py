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
#     Modified: 
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


# In[30]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


name = 'ORACLES'
vv = 'v1'
fp = getpath(name)


# # Load files

# ## Load DARE calculations

# In[4]:


s = hs.loadmat(fp+'model/ORACLES_DARE_v2.mat')


# In[5]:


s.keys()


# In[8]:


ss = hs.loadmat(fp+'model/ORACLES_DARE_aero_prop_v1.mat')


# In[9]:


s['sza'] = ss['sza']


# ## Load DARE parameters from SARE

# In[11]:


sa = sio.loadmat(fp+'ORACLES_2016_DARE_params_v2.mat')


# In[12]:


sa.keys()


# # Plot out data

# In[13]:


s['doy']


# In[14]:


s['doys'] = s['doy']+s['utc']/24.0


# In[15]:


ibad = (s['ref']==0.0) | (s['ref']>=25.0) | (s['cod']==0.0)


# In[16]:


s['cod'][ibad] = np.nan
s['ref'][ibad] = np.nan
s['dn'][ibad,:] = np.nan
s['up'][ibad,:] = np.nan
s['dn_noa'][ibad,:] = np.nan
s['up_noa'][ibad,:] = np.nan
s['dare'][ibad,:] = np.nan


# ## Plot the time trace

# In[17]:


plt.figure()
plt.plot(s['doys'],s['cod'])
plt.xlabel('DOY')
plt.ylabel('COD')
plt.grid()
plt.title('ORACLES 2016 SSFR retrieved COD')


# In[18]:


plt.figure()
plt.plot(s['doys'],s['ref'],'.')
plt.xlabel('DOY')
plt.ylabel('REF')
plt.grid()
plt.title('ORACLES 2016 SSFR retrieved REF')


# ## histogram of retrieved properties

# In[19]:


plt.figure()
plt.hist(s['cod'],range=[0,60],bins=30)
plt.xlabel('COD')
plt.title('ORACLES 2016 SSFR COD histogram')


# In[20]:


plt.figure()
plt.hist(s['ref'],range=[0,30],bins=30)
plt.xlabel('REF [${{\\mu}}$m]')
plt.title('ORACLES 2016 SSFR REF histogram')


# In[71]:


plt.figure()
plt.hist(s['ref'],range=[0,50],bins=30,label='R$_{{eff}}$',alpha=0.6,normed=True,color='g')
plt.hist(s['cod']*2,range=[0,50],bins=30,label='COD',alpha=0.6,normed=True,color='m')

plt.axvline(np.nanmean(s['ref']),ls='-',color='g',label='mean')
plt.axvline(np.nanmedian(s['ref']),ls='--',color='g',label='median')

plt.axvline(np.nanmean(s['cod']),ls='-',color='m')
plt.axvline(np.nanmedian(s['cod']),ls='--',color='m')
plt.legend(frameon=False)

plt.xlabel('R$_{{eff}}$ [$\mu$m], COD')
plt.savefig(fp+'plot_DARE/ORACLES2016_COD_REF_hist.png',dpi=600,transparent=True)


# In[21]:


np.nanmean(s['cod']*2.0),np.nanmedian(s['cod']*2.0)


# In[22]:


np.nanmean(s['ref']),np.nanmedian(s['ref'])


# In[23]:


s['ssa'].shape


# ## Aerosol products

# In[92]:


plt.figure()
plt.plot(s['ssa'][:,2])


# In[58]:


plt.figure()
plt.hist(s['ssa'][:,2],bins=30,range=[np.nanmin(s['ssa'][:,2]),np.nanmax(s['ssa'][:,2])],normed=True,alpha=0.6)
plt.axvline(np.nanmean(s['ssa'][:,2]),ls='-',color='k',label='mean')
plt.axvline(np.nanmedian(s['ssa'][:,2]),ls='--',color='k',label='median')
plt.legend(frameon=False)
plt.xlabel('SSA at 500 nm')
plt.savefig(fp+'plot_DARE/ORACLES2016_SSA_500_hist.png',dpi=600,transparent=True)


# In[75]:


np.nanmean(s['ssa'][:,2]),np.nanmedian(s['ssa'][:,2])


# In[61]:


plt.figure()
plt.hist(s['asy'][:,2],bins=30,range=[np.nanmin(s['asy'][:,2]),np.nanmax(s['asy'][:,2])],normed=True,alpha=0.6,color='orange')
plt.axvline(np.nanmean(s['asy'][:,2]),ls='-',color='k',label='mean')
plt.axvline(np.nanmedian(s['asy'][:,2]),ls='--',color='k',label='median')
plt.legend(frameon=False)
plt.xlabel('ASY at 500 nm')
plt.savefig(fp+'plot_DARE/ORACLES2016_ASY_500_hist.png',dpi=600,transparent=True)


# In[19]:


s['wvl']


# In[23]:


plt.figure()
plt.hist(s['asy'][:,5],range=[0.2,1],bins=30)


# In[ ]:


plt.figure()
plt.hist(s['asy'][:,2],bins=30,range=[np.nanmin(s['asy'][:,2]),np.nanmax(s['asy'][:,2])],normed=True,alpha=0.6,color='orange')
plt.axvline(np.nanmean(s['asy'][:,2]),ls='-',color='k',label='mean')
plt.axvline(np.nanmedian(s['asy'][:,2]),ls='--',color='k',label='median')
plt.legend(frameon=False)
plt.xlabel('ASY at 500 nm')
plt.savefig(fp+'plot_DARE/ORACLES2016_ASY_500_hist.png',dpi=600,transparent=True)


# In[76]:


np.nanmean(s['asy'][:,2]),np.nanmedian(s['asy'][:,2])


# In[14]:


s['ext'].shape


# In[15]:


s['wvl']


# In[111]:


plt.figure()
plt.hist(s['ext'][:,2],range=[np.nanmin(s['ext'][:,2]),np.nanmax(s['ext'][:,2])],bins=30)


# ## DARE

# In[32]:


s.keys()


# In[33]:


import cmaps


# In[34]:


plt.figure()
cmaps.cmaps()


# In[35]:


plt.figure()
plt.hist2d(s['up'][:,2]/s['dn'][:,2],s['sza'],bins=40,range=[[0,1],[0,90]],cmap=plt.cm.jet)
plt.ylabel('SZA [$^\\circ$]')
plt.xlabel('TOA Albedo (BB)')
plt.title('ORACLES 2016 Albedo from DARE calc')

cb = plt.colorbar()
cb.set_label('counts')
pu.prelim()

plt.savefig(fp+'plot/ORACLES_2016_albedo_from_calc_vs_SZA.png',dpi=600,transparent=True)


# In[31]:


plt.figure()
plt.hist2d(s['dare'][:,2],s['sza'],bins=40,range=[[-70,150],[0,90]],cmap=plt.cm.Greens)
plt.ylabel('SZA [$^\\circ$]')
plt.xlabel('DARE [W/m$^2$]')
plt.title('ORACLES 2016 DARE from P3 SSFR and 4STAR - DARE calc')

cb = plt.colorbar()
cb.set_label('counts')
pu.prelim()

plt.savefig(fp+'plot/ORACLES_2016_DARE_from_calc_vs_SZA.png',dpi=600,transparent=True)


# In[53]:


plt.figure()
plt.plot(s['doys'],s['dare'][:,0],'.',label='Surface')
plt.plot(s['doys'],s['dare'][:,1],'.',label='Above cloud')
plt.plot(s['doys'],s['dare'][:,2],'.',label='TOA')
plt.axhline(0,ls=':')
plt.legend(frameon=False)
plt.xlabel('DOY')
plt.ylabel('DARE [W/m$^{{2}}$]')
plt.title('ORACLES 2016 DARE calculations')
plt.grid()
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE.png',dpi=600,transparent=True)


# In[93]:


plt.figure()
plt.hist(s['dare'][:,2],range=[-30,35],bins=30,label='TOA',alpha=0.6)
plt.hist(s['dare'][:,1],range=[-50,0],bins=30,label='Above cloud',alpha=0.6)
plt.hist(s['dare'][:,0],range=[-50,0],bins=30,label='Surface',alpha=0.6)

plt.axvline(np.nanmean(s['dare'][:,2]),ls='-',color='b',label='Mean')
plt.axvline(np.nanmedian(s['dare'][:,2]),ls='--',color='b',label='Median')
plt.axvline(np.nanmean(s['dare'][:,1]),ls='-',color='orange')
plt.axvline(np.nanmedian(s['dare'][:,1]),ls='--',color='orange')
plt.axvline(np.nanmean(s['dare'][:,0]),ls='-',color='g')
plt.axvline(np.nanmedian(s['dare'][:,0]),ls='--',color='g')

pu.prelim()

plt.xlabel('DARE [W/m$^{{2}}$]')
plt.legend(frameon=False)

plt.title('ORACLES 2016 DARE calculations')
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_hist.png',dpi=600,transparent=True)


# In[34]:


igood = np.isfinite(s['dare'][:,0])


# In[35]:


plt.figure()
plt.hist(s['dare'][igood,2],bins=30,range=[-30,35])


# In[36]:


np.nanmean(s['dare'][igood,2]),np.nanmedian(s['dare'][igood,2]),np.nanstd(s['dare'][igood,2])


# ## Make the boxes

# The observations and model data are aggregated within horizontal domains of at least 2o by 2o indicated in Fig. 2. One of the three main regions encompasses the routine flight track, with individual grid boxes centered at (14oE, 24oS), (12oE, 22oS), (10oE, 20oS), (8oE, 18oS), (6oE, 16oS), (4oE, 14oS), (2oE, 12oS) and (0oE, 10oS). Another more coastal north-south track has the southernmost grid box centered on 22oS, spanning between 9oE and 11.75oE. Seven grid boxes are located every 2 degrees north of this, with the northernmost grid box centered on 8oS. A third, zonal track covers the larger domain of the ER2 measurements, with individual grid boxes spanning latitudinally between 10oS and 6oS and separated longitudinally at two degree intervals beginning at 3oW to the west and 13oE in the east. The box for St. Helena Island spans between 6.72 oW and 4.72 oW, between 16.933 oS and 14.933 oS.

# In[37]:


boxes_diag = []
boxes_ns = []
boxes_ew = []


# In[38]:


boxes_diag_ct = [[14.0,-24.0], [12.0,-22.0],[10.0,-20.0],[8.0,-18.0],[6.0,-16.0],[4.0,-14.0],[2.0,-12.0],[0.0,-10.0]]
boxes_ns_ct = [[10.5,-22.0],[10.5,-20.0],[10.5,-18.0],[10.5,-16.0],[10.5,-14.0],[10.5,-12.0],[10.5,-10.0],[10.5,-8.0]]
boxes_ew_ct = [[-3.0,-8.0],[-1.0,-8.0],[1.0,-8.0],[3.0,-8.0],[5.0,-8.0],[7.0,-8.0],[9.0,-8.0],[11.0,-8.0],[13.0,-8.0]]


# Corners are [x0,x1,y0,y1]

# In[39]:


boxes_ns = [[9.0,11.75,i[1]-1.0,i[1]+1.0] for i in boxes_ns_ct]


# In[40]:


boxes_ew = [[-10.0,-6.0,i[0]-1.0,i[0]+1.0] for i in boxes_ew_ct]


# In[41]:


boxes_diag = [[i[0]-1.0,i[0]+1,i[1]-1.0,i[1]+1.0] for i in boxes_diag_ct]


# In[42]:


boxes_diag


# In[43]:


boxes_ew


# In[44]:


boxes_ns


# In[45]:


bins_diag = []
for i,b in enumerate(boxes_diag):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) & (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare'][igood,2]))
    bins_diag.append(s['dare'][igood,2][ia])


# In[46]:


bins_diag[0] = bins_diag[1][0:5]


# In[47]:


bins_ns = []
for i,b in enumerate(boxes_ns):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) & (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare'][igood,2]))
    bins_ns.append(s['dare'][igood,2][ia])


# In[48]:


bins_ns[-1] = bins_ns[-2][0:5]


# In[49]:


bins_ew = []
for i,b in enumerate(boxes_ew):
    ia = (s['lon'][igood]>= b[0]) & (s['lon'][igood]<=b[1]) &(s['lat'][igood]>=b[2]) & (s['lat'][igood]<=b[3]) & (np.isfinite(s['dare'][igood,2]))
    bins_ew.append(s['dare'][igood,2][ia])


# In[50]:


len(boxes_diag),len(bins_diag)


# In[51]:


[fig,ax] = plt.subplots(1,8,figsize=(13,3))

for i,b in enumerate(boxes_diag_ct):
    ax[i].hist(bins_diag[i],bins=30,edgecolor='None',alpha=0.7,color='g',range=(-20,35),
               zorder=10,normed=True,orientation='horizontal',label='Calculations')
    ax[i].hist(sa['bins_diag'][0,i][0],bins=30,edgecolor='None',alpha=0.3,color='m',range=(-30,100),
               zorder=-1,normed=True,orientation='horizontal',label='Parameterization')
    ax[i].axhline(np.nanmean(bins_diag[i]),color='g',label='mean')
    ax[i].axhline(np.nanmedian(bins_diag[i]),color='g',linestyle='--',label='median')
    xmin, xmax = ax[i].get_xlim()
    ax[i].set_xticks(np.round(np.linspace(xmin, xmax, 3), 2))
    if i>0:
        [ag.set_visible(False) for ag in ax[i].yaxis.get_ticklabels()]
    #else:
        #ax[i].legend(frameon=False)
    ax[i].set_title('{}$^\\circ$ ,{}$^\\circ$'.format(b[0],b[1]))
    ax[i].axhline(0,ls=':',color='k',alpha=0.2)
    if i%2: pu.prelim(ax[i])
ax[0].set_ylabel('DARE [W/m$^2$]')
fig.suptitle('ORACLES 2016 Routine Diagonal (Lon,Lat) - TOA')
fig.tight_layout()
#plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_TOA_calc_diag_boxes.png',dpi=600,transparent=True)


# In[ ]:





# In[53]:


[fig,ax] = plt.subplots(1,8,figsize=(13,3))

for i,b in enumerate(boxes_ns_ct):
    ax[i].hist(bins_ns[i],bins=30,edgecolor='None',alpha=0.7,color='b',range=(-50,100),zorder=10,normed=True,orientation='horizontal')
    ax[i].hist(sa['bins_ns'][0,i][0],bins=30,edgecolor='None',alpha=0.2,color='orange',range=(-50,100),zorder=-1,normed=True,orientation='horizontal')
    ax[i].axhline(np.nanmean(bins_ns[i]),color='b',label='mean')
    ax[i].axhline(np.nanmedian(bins_ns[i]),color='b',linestyle='--',label='median')
    xmin, xmax = ax[i].get_xlim()
    ax[i].set_xticks(np.round(np.linspace(xmin, xmax, 3), 2))
    if i>0:
        [ag.set_visible(False) for ag in ax[i].yaxis.get_ticklabels()]
    ax[i].set_title('{}$^\\circ$ ,{}$^\\circ$'.format(b[0],b[1]))
    ax[i].axhline(0,ls=':',color='k',alpha=0.2)
    if i%2: pu.prelim(ax[i])
ax[0].set_ylabel('DARE [W/m$^2$]')
fig.suptitle('ORACLES 2016 North-South (Lon,Lat) - TOA')
fig.tight_layout()
#plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_TOA_calc_ns_boxes.png',dpi=600,transparent=True)


# In[128]:


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
plt.savefig(fp+'plot_DARE/ORACLES_2016_DARE_TOA_calc_map_param.png',dpi=600,transparent=True)


# In[ ]:




