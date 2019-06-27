
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_SSFR_cloud_retrieval
# 
# Purpose:  
# 
#     Run the cloud retrieval on SSFR data from ORACLES 2016 and 2017, but only on the flagacaod times.
#   
# Input:
# 
#     none
# 
# Output:
#    
#     plots
#     mat save files
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
#     Written: Samuel LeBlanc,Santa Cruz, CA, 2018-06-29

# # Prepare the python environment

# In[1]:


import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


import load_utils as lu
import Sp_parameters as Sp


# In[4]:


import hdf5storage as hs
from path_utils import getpath
from write_utils import nearest_neighbor
from tqdm import tqdm_notebook as tqdm
import math


# In[5]:


from scipy import interpolate


# In[6]:


fp = getpath('ORACLES')
fp


# # Load files

# In[7]:


days = ['20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927']


# ## Load the SSFR ict files for 2016

# In[8]:


ssfr_a, ssfr_ah = [],[]
for d in days:
    sf,sfh = lu.load_ict(fp+'data_other/ssfr/SSFR_P3_{}_R1.ict'.format(d),return_header=True)
    ssfr_a.append(lu.recarray_to_dict(sf))
    ssfr_ah.append(sfh)


# In[10]:


ssfr_ah[0]


# In[11]:


ssfr_a[0]


# ## Load the 4STAR files with flagacaod

# In[12]:


star_a, star_ah = [],[]
for d in days:
    sf,sfh = lu.load_ict(fp+'aod_ict/v8/4STAR-AOD_P3_{}_R3.ict'.format(d),return_header=True)
    star_a.append(lu.recarray_to_dict(sf))
    star_ah.append(sfh)


# In[13]:


ssfr_a[3]['Start_UTC'][100]


# In[14]:


star_ah[4]


# ## Get the flagacaod on the timescale of the ssfr measurements

# In[15]:


for i,d in enumerate(days):
    fa = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['flag_acaod'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['flagacaod'] = fa
    am = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['amass_aer'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['airmass'] = am
    ssfr_a[i]['sza'] = np.arccos(1.0/am)*180.0/np.pi
    aod = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['AOD0501'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['AOD_500'] = aod
    a2 = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['AOD_polycoef_a2'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['a2'] = a2
    a1 = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['AOD_polycoef_a1'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['a1'] = a1
    a0 = nearest_neighbor(star_a[i]['Start_UTC'],star_a[i]['AOD_polycoef_a0'],ssfr_a[i]['Start_UTC'],dist=3.0/3600.0)
    ssfr_a[i]['a0'] = a0


# In[16]:


ssfr_a[0]['flagacaod'].shape,ssfr_a[0]['Start_UTC'].shape


# ## Load the LUT for 2wvl reflectance retrieval

# In[30]:


lut = hs.loadmat(fp+'rtm/v5_irr_ORACLES_lut.mat')


# In[31]:


lut.keys()


# ## Combine into one array

# In[17]:


nm = ssfr_a[1].keys()


# In[18]:


ar = {}
for n in ssfr_a[1].keys():
    ar[n] = np.array([])


# In[19]:


ar['days'] = np.array([])


# In[20]:


for i,d in enumerate(days):
    ar['days'] = np.append(ar['days'],np.zeros_like(ssfr_a[i]['Start_UTC'])+i)
    for n in nm:
        try:
            ar[n] = np.append(ar[n],ssfr_a[i][n])
        except:
            print 'problem with :'+n
            ar[n] = np.append(ar[n],ssfr_a[i]['Start_UTC']*0)


# # Format the LUT and data for retrievals

# In[21]:


class so:
    pass


# ## Set up the data

# In[22]:


ar['meas'] = so
ar['meas'].sza = ar['sza']
ar['meas'].Rvis = ar['UP500']/ar['DN500']
ar['meas'].Rnir = ar['UP1650']/ar['DN1650']
ar['meas'].utc = ar['Start_UTC']


# In[23]:


# filter out the bad data. 
bad = (ar['meas'].Rvis > 1.0) & (ar['flagacaod']==0) & (ar['meas'].Rnir > 1.0)
ar['meas'].Rvis[bad] = np.nan
ar['meas'].Rvis[bad] = np.nan


# In[24]:


igood = np.where((np.isfinite(ar['meas'].Rvis)) & (ar['meas'].Rvis > 0.0) & (np.isfinite(ar['meas'].Rnir)) & (ar['meas'].Rnir > 0.0) & (ar['flagacaod']==1))[0]


# ## Plot the histogram of cloud reflectances

# In[71]:


plt.figure()

plt.hist([ar['meas'].Rvis[igood],ar['meas'].Rnir[igood]],bins=30,edgecolor='None',color=['b','r'],alpha=0.7,normed=True,label=['500 nm','1650 nm'])

plt.ylabel('Normalized counts')
plt.xlabel('Cloud Reflectance')
plt.title('Cloud Albedo under aerosol For ORACLES 2016 from P3')
plt.xlim([0,1])


plt.axvline(np.nanmean(ar['meas'].Rvis[igood]),color='b')
plt.axvline(np.nanmedian(ar['meas'].Rvis[igood]),color='b',linestyle='--')
plt.axvline(np.nanmean(ar['meas'].Rnir[igood]),color='r')
plt.axvline(np.nanmedian(ar['meas'].Rnir[igood]),color='r',linestyle='--')

plt.axvline(-0.1,color='k',alpha=0.7,label='Mean')
plt.axvline(-0.1,color='k',alpha=0.7,linestyle='--',label='Median')

plt.legend(frameon=False)

plt.savefig(fp+'plot/Cloud_reflectance_ORACLES_2016.png',dpi=600,transparent=True)


# In[44]:


plt.figure()
plt.hist2d(ar['meas'].Rvis[igood],ar['meas'].sza[igood],bins=40,range=[[0,1],[0,90]])
plt.ylabel('SZA')
plt.xlabel('Cloud reflectance at 500 nm')
cb = plt.colorbar()
cb.set_label('Counts')
plt.savefig(fp+'plot/ORACLES_2016_2dhist_SZA_vs_cloud_refl500.png',dpi=600,transparent=True)


# In[45]:


plt.figure()
plt.hist2d(ar['meas'].Rnir[igood],ar['meas'].sza[igood],bins=40,range=[[0,1],[0,90]])
plt.ylabel('SZA')
plt.xlabel('Cloud reflectance at 1650 nm')
cb = plt.colorbar()
cb.set_label('Counts')
plt.savefig(fp+'plot/ORACLES_2016_2dhist_SZA_vs_cloud_refl1650.png',dpi=600,transparent=True)


# In[70]:


[fig, ax] = plt.subplots(1,3,figsize=(10,3))
ax[0].hist([ar['meas'].Rvis[igood],ar['meas'].Rnir[igood]],bins=30,edgecolor='None',color=['b','r'],alpha=0.7,label=['500 nm','1650 nm'])
ax[0].set_xlim(0,1)
ax[0].set_xlabel('Cloud Reflectance')
ax[0].set_ylabel('Counts')
ax[0].axvline(np.nanmean(ar['meas'].Rvis[igood]),color='b')
ax[0].axvline(np.nanmedian(ar['meas'].Rvis[igood]),color='b',linestyle='--')
ax[0].axvline(np.nanmean(ar['meas'].Rnir[igood]),color='r')
ax[0].axvline(np.nanmedian(ar['meas'].Rnir[igood]),color='r',linestyle='--')

ax[0].axvline(-0.1,color='k',alpha=0.7,label='Mean')
ax[0].axvline(-0.1,color='k',alpha=0.7,linestyle='--',label='Median')
ax[0].legend(frameon=False)

ax[1].set_title('Cloud Albedo under aerosol For ORACLES 2016 from P3')


# In[ ]:


ax[0].


# In[57]:


plt.figure()
plt.hist2d(ar['meas'].Rvis[igood],ar['AOD_500'][igood],bins=40,range=[[0,1],[0,0.8]],cmap=plt.cm.plasma)
plt.ylabel('AOD$_{{500}}$')
plt.xlabel('Cloud reflectance at 500 nm')
cb = plt.colorbar()
cb.set_label('Counts')


# In[47]:


ar.keys()


# ## set up the LUT

# In[257]:


lut.keys()


# In[124]:


lut['tau'].shape, lut['ref'].shape, lut['sza'].shape, lut['irr_dn'].shape, lut['wvl'].shape, lut['zout'], lut['phase']


# In[121]:


lut['Rvis'] = np.zeros([23,27,48])
lut['Rnir'] = np.zeros([23,27,48])


# In[122]:


for ir,r in enumerate(lut['ref']):
    for it,t in enumerate(lut['tau']):
        for iz,s in enumerate(lut['sza']):
            lut['Rvis'][ir,it,iz] = lut['irr_up'][0,0,1,ir,it,iz]/lut['irr_dn'][0,0,1,ir,it,iz]
            lut['Rnir'][ir,it,iz] = lut['irr_up'][0,1,1,ir,it,iz]/lut['irr_dn'][0,1,1,ir,it,iz]


# In[126]:


lut['sza']


# ### Make a hires version of the LUT

# In[416]:


lut['tau_hi'] = np.hstack([np.arange(1.0,25,0.5),np.arange(25,50,1),np.arange(50,102.5,2.5)])
lut['ref_hi'] = np.hstack([np.arange(0,15,0.25),np.arange(15,30.5,0.5)])


# In[417]:


len(lut['tau_hi']), len(lut['ref_hi'])


# In[418]:


lut['Rvis_hi'] = np.zeros([91,94,48])
lut['Rnir_hi'] = np.zeros([91,94,48])


# In[419]:


for i,z in enumerate(lut['sza']):
    fv = interpolate.RectBivariateSpline(lut['ref'],lut['tau'],lut['Rvis'][:,:,i],kx=1,ky=1)
    lut['Rvis_hi'][:,:,i] = fv(lut['ref_hi'],lut['tau_hi'])
    fn = interpolate.RectBivariateSpline(lut['ref'],lut['tau'],lut['Rnir'][:,:,i],kx=1,ky=1)
    lut['Rnir_hi'][:,:,i] = fn(lut['ref_hi'],lut['tau_hi'])


# In[434]:


plt.figure()
for i,r in enumerate(lut['tau_hi']):
    plt.plot(lut['Rvis_hi'][i,:,0],lut['Rnir_hi'][i,:,0],'x-')


# # Run the retrieval

# In[420]:


ar['tau'], ar['ref'] = np.zeros_like(ar['sza'])*np.nan,np.zeros_like(ar['sza'])*np.nan


# In[421]:


ar['ki'] = np.zeros_like(ar['sza'])


# In[429]:


ar['isza'] = []


# In[430]:


pbar = tqdm(total=len(ar['sza']))
for i,s in enumerate(ar['sza']):
    pbar.update()
    if (s>73.0) | (np.isnan(s)):
        continue
    if not i in igood:
        continue
    isza = np.argmin(np.abs(lut['sza']-s))
    ar['isza'].append(isza)
    ki = (ar['meas'].Rvis[i]-lut['Rvis_hi'][:,:,isza]/ar['meas'].Rvis[i])**2+(ar['meas'].Rnir[i]-lut['Rnir_hi'][:,:,isza]/ar['meas'].Rnir[i])**2
    kimin = np.unravel_index(np.nanargmin(ki),ki.shape)
    ar['ki'][i] = np.nanmin(ki)
    ar['tau'][i],ar['ref'][i] = lut['tau_hi'][kimin[1]],lut['ref_hi'][kimin[0]]


# In[432]:


plt.figure()
plt.hist(ar['isza'],bins=20)


# # Plot the retrieval results

# In[423]:


plt.figure()
plt.plot(ar['tau'],'.')


# In[426]:


np.nanmean(ar['tau'])


# In[427]:


plt.figure()
plt.hist(ar['tau'][np.isfinite(ar['tau'])])


# In[216]:


plt.figure()
plt.plot(ar['ref'],'.')


# In[424]:


len(np.where(np.isfinite(ar['ref']))[0])


# In[425]:


len(ar['ref'])


# In[269]:


plt.figure()
plt.plot(np.where(np.isfinite(ar['ref']))[0])


# # Save the retrieved output

# In[263]:


out = {}


# In[264]:


ar.keys()


# In[265]:


out['tau'] = ar['tau']
out['ref'] = ar['ref']
out['sza'] = ar['sza']
out['aod'] = ar['AOD_500']
out['days'] = ar['days']
out['utc'] = ar['Start_UTC']
out['lat'] = ar['LAT']
out['lon'] = ar['LON']
out['a0'],out['a1'],out['a2'] = ar['a0'],ar['a1'],ar['a2']


# In[266]:


fp


# In[267]:


hs.savemat(fp+'data_other/ssfr_2016_retrieved_COD.mat',out)


# In[268]:


out['days']


# # Load the results of the CRE calculations
# See the ORACLES_cld_CRE_from_SSFR jupyter notebook for details

# In[273]:


c = hs.loadmat(fp+'rtm/ORACLES_CRE_v6_irr.mat')


# In[274]:


c.keys()


# In[276]:


c['ssfr_aero_CRE'].keys()


# In[277]:


CRE_toa = c['ssfr_aero_CRE']['up'][:,2]-c['ssfr_aero_CRE_clear']['up'][:,2]


# In[282]:


flt = np.isfinite(out['aod']) & (CRE_toa>0.0)


# In[294]:


out['days']


# In[295]:


plt.figure()
for i in np.unique(out['days']):
    fi = np.where((flt==1) & (i==out['days']))
    plt.plot(out['aod'][fi],CRE_toa[fi],'.')
plt.xlabel('AOD$_{{500nm}}$')
plt.ylabel('TOA CRE [W/m$^2$]')


# In[284]:


plt.figure()
plt.plot(out['aod'][flt],c['ssfr_aero_C'][flt,0],'.')
plt.xlabel('AOD$_{{500nm}}$')
plt.ylabel('Surface CRE [W/m$^2$]')


# In[293]:


plt.figure()
plt.cm.viridis
plt.hist2d(out['aod'][flt],CRE_toa[flt],bins=30)
plt.colorbar()


# In[300]:


import seaborn as sns
import plotting_utils as pu


# In[299]:


plt.figure()
#plt.hexbin(out['aod'][flt],CRE_toa[flt])
sns.kdeplot(out['aod'][flt],CRE_toa[flt], shade=True)


# In[307]:


plt.figure()
plt.plot(out['aod'][flt],CRE_toa[flt],'.',alpha=0.05)
pu.plot_lin(out['aod'][flt],CRE_toa[flt])
plt.legend(frameon=False)
plt.xlabel('AOD$_{{500nm}}$')
plt.ylabel('TOA CRE [W/m$^2$]')


# In[308]:


plt.figure()
plt.plot(out['aod'][flt],c['ssfr_aero_C'][flt,0],'.',alpha=0.05)
pu.plot_lin(out['aod'][flt],c['ssfr_aero_C'][flt,0])
plt.legend(frameon=False)
plt.xlabel('AOD$_{{500nm}}$')
plt.ylabel('Surface CRE [W/m$^2$]')


# In[324]:


plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(out['aod'][flt],CRE_toa[flt],'.',color='grey',alpha=0.05)
pu.plot_lin(out['aod'][flt],CRE_toa[flt],color='r')
plt.legend(frameon=False,loc=2)
plt.ylabel('TOA CRE [W/m$^2$]')
plt.xlim(0,.9)

ax2 = plt.subplot(2,1,2)
plt.plot(out['aod'][flt],c['ssfr_aero_C'][flt,0],'.',color='blue',alpha=0.05)
pu.plot_lin(out['aod'][flt],c['ssfr_aero_C'][flt,0],color='r')
plt.legend(frameon=False,loc=3)
plt.xlabel('AOD$_{{500nm}}$')
plt.ylabel('Surface CRE [W/m$^2$]')
plt.xlim(0,.9)


# In[325]:


rCRE_toa = CRE_toa/c['ssfr_aero_CRE']['dn'][:,2]*100.0
rCRE_sur = c['ssfr_aero_C'][:,0]/c['ssfr_aero_CRE']['dn'][:,2]*100.0


# In[328]:


plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(out['aod'][flt],rCRE_toa[flt],'.',color='grey',alpha=0.05)
pu.plot_lin(out['aod'][flt],rCRE_toa[flt],color='r',ci=0.99)
plt.legend(frameon=False,loc=1)
plt.ylabel('TOA relative CRE [%]')
plt.xlim(0,.9)

ax2 = plt.subplot(2,1,2)
plt.plot(out['aod'][flt],rCRE_sur[flt],'.',color='blue',alpha=0.05)
pu.plot_lin(out['aod'][flt],rCRE_sur[flt],color='r',ci=0.99)
plt.legend(frameon=False,loc=4)
plt.xlabel('AOD$_{{500nm}}$')
plt.ylabel('Surface relative CRE [%]')
plt.xlim(0,.9)


# In[330]:


help(pu.make_boxplot)


# In[337]:


lims = np.arange(0,1.3,0.1)
pos = np.arange(0.05,1.2,0.1)


# In[339]:


lims,pos,len(lims),len(pos)


# In[362]:


fp


# In[404]:


plt.figure()
ax1 = plt.subplot(2,1,1)
pu.make_boxplot(rCRE_toa[flt],out['aod'][flt],lims,pos,color='k',fliers_off=True,widths=0.09,patch_artist=True,alpha=0.5)
pu.plot_lin(out['aod'][flt],rCRE_toa[flt],color='r',ci=0.99,zorder=200)
plt.xlim(0.0,1.2)
plt.legend(frameon=False,loc=1)
plt.ylabel('TOA relative CRE [%]')

ax2 = plt.subplot(2,1,2)
pu.make_boxplot(rCRE_sur[flt],out['aod'][flt],lims,pos,color='blue',fliers_off=True,widths=0.09,patch_artist=True,alpha=0.6)
pu.plot_lin(out['aod'][flt],rCRE_sur[flt],color='r',ci=0.99,zorder=200)
plt.xlim(0.0,1.2)
plt.legend(frameon=False,loc=4)
plt.ylabel('Surface relative CRE [%]')
plt.xlabel('AOD$_{{500nm}}$')
plt.savefig(fp+'plot/SSFR_CRE_vs_AOD.png',transparent=True,dpi=600)


# In[387]:


flo = (flt==1) & (out['tau']>0)


# In[446]:


plt.figure(figsize=(7,4))
plt.subplot(2,1,1)
plt.hist(out['tau'][flo]*4.0,normed=True,edgecolor='None',color='g',alpha=0.6,bins=100)
plt.xlabel('COD')
plt.xlim(0,80)
plt.ylabel('Normed counts')
plt.grid()

plt.subplot(2,1,2)
plt.hist(out['ref'][flo],normed=True,edgecolor='None',color='grey',alpha=0.6,bins=30)
plt.xlabel('r$_{{eff}}$ [$\mu$m]')
plt.ylabel('Normed counts')
plt.grid()
plt.xlim(0,30)
plt.tight_layout()
plt.savefig(fp+'plot/SSFR_COD_ref_ORACLES2016_flagacaod.png',transparent=True,dpi=600)


# In[445]:


np.nanmean(out['tau'][flo]*4.0)


# In[408]:


import scipy.stats as st


# In[449]:


a,xe,ye,bn = st.binned_statistic_2d(out['lat'][flo],out['lon'][flo],out['tau'][flo]*4.0,
                           bins=26,range=[[-25,-8],[0,16]],statistic='mean')
a = np.ma.masked_array(a,np.isnan(a))


# In[450]:


plt.figure()
p = plt.pcolor(ye,xe,a,vmin=0.0,vmax=1.0,cmap='plasma')

plt.xlabel(u'Longitude [°]')
plt.ylabel(u'Latitude [°]')
plt.title('Mean COD')

cb = plt.colorbar(p,extend='both')
cb.set_label('Mean COD')


# In[405]:


np.nanmean(out['tau'][flo])


# In[406]:


np.nanmean(out['ref'][flo])


# In[394]:


fp

