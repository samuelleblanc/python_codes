
# coding: utf-8

# # Intro
# Name:  
# 
#     ASR_2019_Proposal_MPL_SASZe
# 
# Purpose:  
# 
#     Make some figures to be put into the proposal, showing MPL data of aerosol over clouds, and examples of SASZe data. 
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures and save files...
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - plotting_utils (user defined plotting routines)
#     - hdf5storage
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ARM netcdf files for MPL and SASZe
#   
#  Modification History:
#  
#      Written: by Samuel LeBlanc, Halfway to Tahoe, 2019-03-08

# # Load the python modules and setup paths

# In[1]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
import scipy.io as sio
import matplotlib.pyplot as plt


# In[2]:


import load_utils as lu
from path_utils import getpath


# In[3]:


get_ipython().magic(u'matplotlib notebook')


# In[4]:


import scipy.stats as st
import Sun_utils as su
import plotting_utils as pu


# In[5]:


import Sp_parameters as sp


# In[6]:


fp = getpath('LASIC')
fp


# # Load the files

# ## Load the MPL files

# In[7]:


daystr = '20160814'


# In[9]:


ff = fp+'data/MPL/asi30smplcmask1zwangM1.c1.{}.00000?.cdf'.format(daystr)


# In[8]:


fl = os.listdir(fp+'data/MPL/')


# In[9]:


for i,f in enumerate(fl):
    if daystr in f:
        fo = f


# In[10]:


fo


# In[11]:


mpl,mplh = lu.load_netcdf(fp+'data/MPL/'+fo,everything=True)


# In[14]:


mplh['qc_backscatter']


# In[15]:


mplh['backscatter']


# In[16]:


mplh['background_signal']


# In[20]:


mplh['time']


# In[23]:


mplh['height']


# ## Load the SASZe data

# In[158]:


fs = os.listdir(fp+'data/SASZe/')


# In[159]:


for i,f in enumerate(fs):
    if daystr in f:
        if 'vis' in f:
            fgv = f
        if 'nir' in f:
            fgn = f


# In[160]:


vis,vish = lu.load_netcdf(fp+'data/SASZe/'+fgv,everything=True)


# In[162]:


vish['zenith_radiance']


# In[165]:


vish['wavelength']


# In[163]:


vish['time']


# In[161]:


nir,nirh = lu.load_netcdf(fp+'data/SASZe/'+fgn,everything=True)


# In[166]:


nirh['zenith_radiance']


# ### Build the full spectra from vis and nir

# In[243]:


vis['wavelength'][100:1200]


# In[244]:


nir['wavelength'][15:]


# In[286]:


rad = np.append(vis['zenith_radiance'][:,100:1200],nir['zenith_radiance'][:,15:]*2.0,axis=1)
wvl = np.append(vis['wavelength'][100:1200],nir['wavelength'][15:])


# In[223]:


wvl.shape


# In[224]:


rad.shape


# In[287]:


s = {'zenlambda':wvl,'rad':rad,'utc':vis['time']/3600.0,'good':np.where(vis['solar_zenith']<70.0)[0]}


# In[288]:


sr = sp.Sp(s)


# ## Load the LUT for zenith spectra from ORACLES

# In[255]:


fpo = getpath('ORACLES')


# In[256]:


fpo


# In[263]:


model = hs.loadmat(fpo+'model/v3_ORACLES_lut.mat')


# In[259]:


model.keys()


# In[267]:


1.0/np.cos(model['sza']*np.pi/180.0)


# In[270]:


sptemp = {}
sptemp['tau'] = model['tau']
sptemp['ref'] = model['ref']
sptemp['zout'] = model['zout']
sptemp['sza'] = model['sza']
sptemp['phase'] = model['phase']
sptemp['irr_dn_diff'] = model['irr_dn_diff'][:,:,:,:,:,1]
sptemp['irr_dn'] = model['irr_dn'][:,:,:,:,:,1]
sptemp['irr_up'] = model['irr_up'][:,:,:,:,:,1]
sptemp['wvl'] = [model['wvl']]
sptemp['rad'] = model['rad'][:,:,:,:,:,1]


# In[271]:


lut = sp.Sp(sptemp)


# In[319]:


lut.tau


# In[321]:


lut.ref,lut.ref.shape


# In[323]:


lut.norm[0,:,0,4,16]


# # MPL plotting

# In[35]:


mpl['backscatter'][0,:]


# In[12]:


bsc = mpl['backscatter'].T
bsc[bsc <0.005] = np.nan
bsc[bsc >10.005] = np.nan


# In[13]:


ib = np.where(mpl['cloud_base']>0)[0]
iafter = np.append(ib[np.where(np.diff(ib,1)>1)[0]]+1,ib[-1]+1)
ibefore = np.append(ib[0]-1,ib[np.where(np.diff(ib,1)>1)[0]+1]-1)


# In[14]:


len(iafter),len(ibefore)


# In[140]:


ib


# In[150]:


mpl['time'][iafter]/3600.0


# In[ ]:


bsc


# In[99]:


bsc.shape


# In[105]:


mpl['height'][70]


# In[129]:


x


# In[15]:


#plt.figure()
fig,ax = plt.subplots(1,1,figsize=(12,3))
pc = ax.pcolorfast(mpl['time']/3600.0,mpl['height'],bsc[:-1,:-1],vmax=5.0,vmin=0.0005,cmap='viridis')
ax.plot(mpl['time']/3600.0,mpl['cloud_base'],'.r',label='cloud base')
pre,post = [],[]
for x in np.where(mpl['cloud_base']>0)[0]:
    ax.axvline(mpl['time'][x]/3600.0,color='#DDDDDD',alpha=0.3,linewidth=5.0)
ax.axvline(x,color='#DDDDDD',alpha=0.3,linewidth=5.0,label='cloud event')

for i,y in enumerate(iafter):
    if bsc[70,ibefore[i]]>0.001 and bsc[70,y]>0.001:
        #print i,mpl['time'][ibefore[i]]/3600.0,mpl['time'][y]/3600.0
        ax.axvspan(mpl['time'][ibefore[i]]/3600.0,mpl['time'][y]/3600.0,color='r',alpha=0.2)
        ax.axvline(mpl['time'][ibefore[i]]/3600.0,color='g',alpha=0.5)
        ax.axvline(mpl['time'][y]/3600.0,color='g',alpha=0.5)
ax.axvspan(mpl['time'][ibefore[i]]/3600.0,mpl['time'][y]/3600.0,color='r',alpha=0.2,label='AAC event')
ax.axvline(mpl['time'][y]/3600.0,color='g',alpha=0.5,label='Before/After cloud')
#for x in mpl['time'][bsc[70,:]>0.1]/3600.0:
#     ax.axvline(x,color='r',alpha=0.05,linewidth=5.0)
#ax.axvline(x,color='r',alpha=0.05,linewidth=5.0,label='AAC event')


#ax.plot(mpl['time']/3600.0,mpl['cloud_top'],'.r',label='cloud top')

plt.legend(frameon=True,numpoints=1)
aa = plt.colorbar(pc)
#aa.set_label('Total attenuated backscatter\n[counts/microsecond]')
aa.set_label('Lidar signal')
ax.set_ylabel('Altitude [km]')
ax.set_xlabel('Time UTC [H]')
ax.set_xlim([19,22])
ax.set_ylim([0,4.5])
ax.set_title('MPL - Ascension Island, 2016-08-14')
plt.tight_layout()
plt.savefig(fp+'data/MPL_ASI_AAC_events.png',dpi=600,transparent=True)


# # SASZe plotting

# In[186]:


np.argmin(abs(sr.utc-19.45))


# In[208]:


np.argmin(abs(sr.utc-14.57))


# In[268]:


vis['solar_zenith'][23760]


# In[265]:


1.0/np.cos(vis['solar_zenith'][23760]*np.pi/180.0)


# In[251]:


plt.figure()
plt.plot(sr.utc,sr.sp[:,400])


# In[298]:


import matplotlib.ticker as ticker
from mpltools import color


# In[304]:


def frmt(x,pos):
    return '{:5.2f}'.format(x)


# In[310]:


len(range(23000,25000,50))


# In[356]:


plt.figure(figsize=(10,3))
color.cycle_cmap(41,cmap=plt.cm.plasma,ax=plt.gca())
for j in xrange(23000,25000,50):
    plt.plot(sr.wvl,sr.norm[j,:])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.plasma)
scalarmap.set_array(sr.utc[range(23000,25000,50)])
cba = plt.colorbar(scalarmap,format=ticker.FuncFormatter(frmt))
cba.set_label('UTC [h]')
plt.ylim(0,1)
plt.xlim(340,1750)
plt.title('SASZe normalized spectra, Ascension Island, 2016-08-14')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalized Radiance')
plt.tight_layout()
#plt.plot(lut.wvl,lut.norm[0,:,0,8,20],'k-',lw=3,
#         label='Modeled, $\\tau$={:2.1f}, r$_{{eff}}$={:2.0f}$\\mu$m'.format(lut.tau[20],lut.ref[8]))
plt.plot(lut.wvl,lut.norm[0,:,0,4,16],'-',color='k',lw=3,alpha=0.8,
         label='Modeled, $\\tau$={:2.1f}, r$_{{eff}}$={:2.0f}$\\mu$m'.format(lut.tau[16],lut.ref[4]))
plt.plot(lut.wvl,lut.norm[0,:,0,4,1],'-',color='lightgrey',lw=3,alpha=0.8,
         label='Modeled, $\\tau$={:2.1f}, r$_{{eff}}$={:2.0f}$\\mu$m'.format(lut.tau[1],lut.ref[4]))
plt.legend(handletextpad=0.08,labelspacing=0.05,frameon=True)
plt.savefig(fp+'data/SASZe_sample_spectra_with_model.png',dpi=600,transparent=True)


# In[184]:


fig = sp.plt_zenrad(sr)
plt.ylim(0,0.5)
#plt.figure()


# In[253]:


sp.plt_norm_zenrad(sr)

