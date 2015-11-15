
# coding: utf-8

# Name:  
# 
#     SEAC4RS_compare_v4
# 
# Purpose:  
# 
#     Python script that is used to step through the building and use of the parameters from SEAC4RS
#     Starts with the measured 4STAR zenith radiances, then loads the idl save file containing the modeled lut for that day
#     Regroups all the necessary steps to build all the figures used to analyze the data
#     Runs the ki^2 retrieval with 15 parameters
#     plots the results
#     Compares to eMAS retrievals for that day
#     Also loads the SSFR irradiances from the ER2 and calculates the cloud properties from above
#     Uses the Baum Aggregate Solid Columns (sp_v3_20130913) instead of the General Habit Mixtures 
# 
# Calling Sequence:
# 
#     python SEAC4RS_compare_v3.py
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
#     - run_kisq_retrieval.py : for the retrieval functions
#     - load_modis.py : for loading modis files
#     - matplotlib
#     - mpltools
#     - numpy
#     - plotly : optional
#     - scipy : for saving and reading
#     - math
#     - os
#     - gc
#     - pdb
#     - datetime
#     - pyhdf
#     - mpl_toolkits
#     - gdal (from osgeo)
#     - plotting_utils (user defined plotting routines)
#     - map_utils, dependent on geopy
#     - hdf5storage
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - sp_v3_20130913_4STAR.out: modeled spectra output for SEAC4RS at sza 17.9, idl save file
#   - %%20130219starzen_rad.mat : special zenith radiance 4star matlab file 
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, NASA Ames
#     Modified: Samuel LeBlanc, NASA Ames, Santa Cruz, CA, 2015-10-05
#             - ported from version 3, to 4
#             - new normalization of radiance spectra applied, new SSFR data used
#     Modified: Samuel LeBlanc, NASA Ames, Flying to St-John's, Newfoundland, 2015-11-14
#             - added saving of retrieved properties in matlab format with hdf5storrage

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
#matplotlib.use('nbagg')
import numpy as np, h5py
#import plotly.plotly as py
import scipy.io as sio
import math
import os
import Sp_parameters as Sp
import hdf5storage as hs
#py.sign_in("samuelleblanc", "4y3khh7ld4")
#import mpld3
#mpld3.enable_notbeook()


# In[2]:

get_ipython().magic(u'matplotlib notebook')


# In[3]:

import IPython
IPython.InteractiveShell.cache_size = 0


# In[4]:

# set the basic directory path
fp='C:/Users/sleblan2/Research/SEAC4RS/'


# ## Get the lookup table for the 4STAR data

# In[6]:

# load the idl save file containing the modeled radiances
vv = 'v2'
s=sio.idl.readsav(fp+'model\\sp_'+vv+'_20130913_4STAR.out')#fp+'model/sp_v1.1_20130913_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp, wvl, z, re, ta)'
# create custom key for sorting via wavelength
iwvls = np.argsort(s.zenlambda)
s.wv = np.sort(s.zenlambda)


# In[21]:

if 'Sp' in locals():
    reload(Sp)
if 'lut' in locals():
    del lut
    import gc; gc.collect()


# In[22]:

lut = Sp.Sp(s,irrad=True)
lut.params()
lut.param_hires()


# In[23]:

lut.sp_hires()


# In[24]:

print lut.tau.shape
print lut.ref.shape
print lut.sp.ndim
print lut.par.size
print lut.par.shape
print lut.ref
print lut.par.shape


# In[9]:

fig3,ax3 = plt.subplots(4,4,sharex=True,figsize=(15,8))
ax3 = ax3.ravel()

for i in range(lut.npar-1):
    color.cycle_cmap(len(lut.ref[lut.ref<30]),cmap=plt.cm.RdBu,ax=ax3[i])
    for j in xrange(len(lut.ref)):
        ax3[i].plot(lut.tau,lut.par[0,j,:,i])
    ax3[i].set_title('Parameter '+str(i+1))
    ax3[i].grid()
    ax3[i].set_xlim([0,100])
    ax3[i].set_ylabel('$\\eta_{%i}$' % (i+1))
    if i > 10: 
        ax3[i].set_xlabel('$\\tau$')

for tk in ax3[11].get_xticklabels():
    tk.set_visible(True)
ax3[-1].axis('off')

fig3.tight_layout()
plt.suptitle('Liquid')
plt.subplots_adjust(top=0.93,right=0.93)

cbar_ax = fig3.add_axes([0.95,0.10,0.02,0.8])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.RdBu,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.ref[lut.ref<30])
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
cba.ax.set_ylabel('$r_{eff}$ [$\\mu$m]')
cba.ax.set_yticklabels(np.linspace(lut.ref[0],29,6));

plt.savefig(fp+'plots/20130913_'+vv+'_param_liq_vs_tau.png',dpi=600,transparent=True)
plt.show()


# In[48]:

fig4,ax4 = plt.subplots(4,4,sharex=True,figsize=(15,8))
ax4 = ax4.ravel()

for i in range(lut.npar-1):
    color.cycle_cmap(len(lut.tau),cmap=plt.cm.gist_ncar,ax=ax4[i])
    for j in xrange(len(lut.tau)):
        ax4[i].plot(lut.ref,lut.par[0,:,j,i])
    ax4[i].set_title('Parameter '+str(i+1))
    ax4[i].grid()
    ax4[i].set_xlim([0,30])
    ax4[i].set_ylabel('$\\eta_{%i}$' % (i+1))
    if i > 10: 
        ax4[i].set_xlabel('$r_{eff}$ [$\mu$m]')

for tk in ax4[11].get_xticklabels():
    tk.set_visible(True)
ax4[-1].axis('off')
fig4.tight_layout()
plt.suptitle('Liquid')
plt.subplots_adjust(top=0.93,right=0.93)

cbar_ax = fig4.add_axes([0.95,0.10,0.02,0.8])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.tau)
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
cba.ax.set_ylabel('$\\tau$')
labels = ['%5.1f' %F for F in np.linspace(lut.tau[0],lut.tau[-1],6)]
cba.ax.set_yticklabels(labels);

plt.savefig(fp+'plots/20130913_'+vv+'_param_liq_vs_ref.png',dpi=600,transparent=True)
plt.show()


# In[49]:

fig5,ax5 = plt.subplots(4,4,sharex=True,figsize=(15,8))
ax5 = ax5.ravel()

for i in range(lut.npar-1):
    color.cycle_cmap(len(lut.ref),cmap=plt.cm.RdBu,ax=ax5[i])
    for j in xrange(len(lut.ref)):
        ax5[i].plot(lut.tau,lut.par[1,j,:,i])
    ax5[i].set_title('Parameter '+str(i+1))
    ax5[i].grid()
    ax5[i].set_xlim([0,100])
    ax5[i].set_ylabel('$\\eta_{%i}$' % (i+1))
    if i > 10: 
        ax5[i].set_xlabel('$\\tau$')

for tk in ax5[11].get_xticklabels():
    tk.set_visible(True)
ax5[-1].axis('off')
fig5.tight_layout()
plt.suptitle('Ice')
plt.subplots_adjust(top=0.93,right=0.93)
cbar_ax = fig5.add_axes([0.95,0.10,0.02,0.8])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.RdBu,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.ref)
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
cba.ax.set_ylabel('$r_{eff}$ [$\\mu$m]')
cba.ax.set_yticklabels(np.linspace(lut.ref[0],lut.ref[-1],6));

plt.savefig(fp+'plots/20130913_'+vv+'_param_ice_vs_tau.png',dpi=600,transparent=True)
plt.show()


# In[50]:

fig6,ax6 = plt.subplots(4,4,sharex=True,figsize=(15,8))
ax6 = ax6.ravel()
for i in range(lut.npar-1):
    color.cycle_cmap(len(lut.tau),cmap=plt.cm.gist_ncar,ax=ax6[i])
    for j in xrange(len(lut.tau)):
        ax6[i].plot(lut.ref,lut.par[1,:,j,i])
    ax6[i].set_title('Parameter '+str(i+1))
    ax6[i].grid()
    ax6[i].set_xlim([0,50])
    ax6[i].set_ylabel('$\\eta_{%i}$' % (i+1))
    if i > 10: 
        ax6[i].set_xlabel('$r_{eff}$ [$\mu$m]')

for tk in ax6[11].get_xticklabels():
    tk.set_visible(True)
ax6[-1].axis('off')
fig6.tight_layout()
plt.suptitle('Ice')
plt.subplots_adjust(top=0.93,right=0.93)
cbar_ax = fig6.add_axes([0.95,0.10,0.02,0.8])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.tau)
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
cba.ax.set_ylabel('$\\tau$')
labels = ['%5.1f' %F for F in np.linspace(lut.tau[0],lut.tau[-1],6)]
cba.ax.set_yticklabels(labels);

plt.savefig(fp+'plots/20130913_'+vv+'_param_ice_vs_ref.png',dpi=600,transparent=True)
plt.show()


# ## Get the DC8 nav data

# In[25]:

import load_modis
reload(load_modis)
from load_modis import mat2py_time, toutc, load_ict


# In[26]:

dc8,dc8_header = load_ict(fp+'dc8/20130913/SEAC4RS-MMS-1HZ_DC8_20130913_R0.ict',return_header=True)


# In[8]:

plt.figure()
plt.plot(dc8['TIME_UTC'],dc8['G_ALT'])
plt.xlabel('UTC [h]')
plt.ylabel('GPS Altitude [m]')
plt.title('DC8 Altitude on 2013-09-13')
plt.savefig(fp+'plots/20130913_DC8_alt.png',dpi=600,transparent=True)


# In[27]:

print dc8['TIME_UTC'][12000]/3600.0
print dc8['G_ALT'][12000]


# In[10]:

plt.figure()
plt.plot(dc8['TIME_UTC'],dc8['W'])
plt.xlabel('UTC [h]')
plt.ylabel('Vertical Airspeed [m/s]')
plt.title('DC8 Vertical airspeed on 2013-09-13')
plt.savefig(fp+'plots/20130913_DC8_W.png',dpi=600,transparent=True)


# ## Get the 4STAR data

# In[28]:

# load the matlab file containing the measured TCAP radiances
mea = sio.loadmat(fp+'../4STAR/SEAC4RS/20130913/20130913starzen_3.mat')
mea.keys()


# Go through and get the radiances for good points, and for the time selected

# In[29]:

print mea['t']
tt = mat2py_time(mea['t'])
mea['utc'] = toutc(tt)


# In[30]:

mea['good'] = np.where((mea['utc']>18.5) & (mea['utc']<19.75) & (mea['Str'].flatten()!=0) & (mea['sat_time'].flatten()==0))


# In[31]:

mea['w'][0][1068]


# In[39]:

plt.plot(mea['w'][0])


# In[16]:

plt.figure()
plt.plot(mea['utc'],abs(mea['rad'][:,1024]-mea['rad'][:,1068]))
plt.xlabel('UTC [Hours]')
plt.ylabel('difference in radiance at 1024 nm')


# In[17]:

plt.figure()
plt.plot(mea['utc'],mea['rad'][:,400]/1000.)
plt.xlabel('UTC [hours]')
plt.ylabel('Radiance at 400 nm [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')


# In[32]:

reload(Sp)
if 'meas' in locals():
    del meas
    import gc; gc.collect()
# first convert measurements to Sp class, with inherent parameters defined
meas = Sp.Sp(mea)
meas.params()


# In[33]:

meas.sp.shape


# In[15]:

fig = plt.figure()
color.cycle_cmap(len(meas.utc),cmap=plt.cm.gist_ncar,ax=plt.gca())
for i in range(len(meas.utc)):
    plt.plot(meas.wvl,meas.sp[i,:]/1000.0)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Radiance [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
plt.title('All radiance spectra')
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar)
scalarmap.set_array(meas.utc)
cba = plt.colorbar(scalarmap)
cba.set_label('UTC [h]')


# In[21]:

fig = plt.figure()
n=100
color.cycle_cmap(len(meas.utc)/n,cmap=plt.cm.gist_ncar,ax=plt.gca())
for i in range(len(meas.utc)):
    if not i%n:
        plt.plot(meas.wvl,meas.sp[i,:]/1000.0)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Radiance [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
plt.title('All radiance spectra')
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar)
scalarmap.set_array(meas.utc)
cba = plt.colorbar(scalarmap)
cba.set_label('UTC [h]')
plt.savefig(fp+'plots/20130913_zenrad_spotcheck.png',dpi=600,transparent=True)


# In[34]:

isubwvl = np.where((meas.wvl>315.0)&(meas.wvl<940.0))[0]


# In[35]:

meas.sp.shape


# In[111]:

fig = plt.figure()
color.cycle_cmap(len(meas.utc),cmap=plt.cm.gist_ncar,ax=plt.gca())
for i in range(len(meas.utc)):
    plt.plot(meas.wvl,meas.norm[i,:])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Normalized Radiance')
plt.ylim([0,1])
plt.title('All normalized radiance spectra')


# In[23]:

fig,ax = plt.subplots(1,1,figsize=(8,14))
pco = ax.pcolorfast(meas.wvl,np.where(meas.utc)[0],meas.sp[:-1,:-1]/1000.0,cmap='gist_ncar',vmin=0,vmax=0.8)
cba = plt.colorbar(pco)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Measurement number')
plt.title('All zenith radiance spectra')
cba.set_label('Zenith radiance [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
plt.savefig(fp+'plots/20130913_curtain_zenrad_num.png',dpi=600,transparent=True)


# In[24]:

fig,ax = plt.subplots(1,1,figsize=(8,14))
pco = ax.pcolorfast(meas.wvl,meas.utc,meas.sp[:-1,:-1]/1000.0,cmap='gist_ncar',vmin=0,vmax=0.8)
cba = plt.colorbar(pco)
plt.xlabel('Wavelength [nm]')
plt.ylabel('UTC [h]')
plt.title('All zenith radiance spectra')
cba.set_label('Zenith radiance [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
plt.savefig(fp+'plots/20130913_curtain_zenrad_utc.png',dpi=600,transparent=True)


# In[127]:

fig,ax = plt.subplots(1,1,figsize=(8,14))
pco = ax.pcolorfast(meas.wvl,np.where(meas.utc)[0],meas.norm[:-1,:-1],cmap='gist_ncar',vmin=0,vmax=1)
plt.colorbar(pco)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Measurement number')
plt.title('All normalized radiance spectra')


# ## Run the retrieval on 4STAR data

# In[36]:

from Sp_parameters import smooth


# In[37]:

import run_kisq_retrieval as rk
reload(rk)


# In[38]:

subp = [1,2,4,5,6,10,11,13]
#subp = [2,3,5,6,7,11,12,14]


# In[39]:

print max(meas.good)
print type(mea['good'])
print isinstance(mea['good'],tuple)
print type(mea['good'][0])


# In[40]:

(meas.tau,meas.ref,meas.phase,meas.ki) = rk.run_retrieval(meas,lut)


# In[41]:

print meas.utc.shape
print len(meas.good), max(meas.good)


# In[135]:

fig,ax = plt.subplots(4,sharex=True)
ax[0].set_title('Retrieval results time trace')
ax[0].plot(meas.utc,meas.tau,'rx')
ax[0].plot(meas.utc[meas.good],smooth(meas.tau[meas.good],20),'k')
ax[0].set_ylabel('$\\tau$')
ax[1].plot(meas.utc,meas.ref,'g+')
ax[1].set_ylabel('R$_{ef}$ [$\\mu$m]')
ax[1].plot(meas.utc[meas.good],smooth(meas.ref[meas.good],20),'k')
ax[2].plot(meas.utc,meas.phase,'k.')
ax[2].set_ylabel('Phase')
ax[2].set_ylim([-0.5,1.5])
ax[2].set_yticks([0,1])
ax[2].set_yticklabels(['liq','ice'])
ax[3].plot(meas.utc,meas.ki)
ax[3].set_ylabel('$\\chi^{2}$')
ax[3].set_xlabel('UTC [Hours]')
ax[3].set_xlim([18.5,19.05])
plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results_v4.png',dpi=600)
#plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results_v3.pdf',bbox='tight')


# In[136]:

from Sp_parameters import smooth
fig,ax = plt.subplots(4,sharex=True)
ax[0].set_title('Retrieval results time trace')
ax[0].plot(meas.utc,meas.tau,'rx')
ax[0].plot(meas.utc[meas.good],smooth(meas.tau[meas.good],20),'k')
ax[0].set_ylabel('$\\tau$')
ax[1].plot(meas.utc,meas.ref,'g+')
ax[1].set_ylabel('R$_{ef}$ [$\\mu$m]')
ax[1].plot(meas.utc[meas.good],smooth(meas.ref[meas.good],20),'k')
ax[2].plot(meas.utc,meas.phase,'k.')
ax[2].set_ylabel('Phase')
ax[2].set_ylim([-0.5,1.5])
ax[2].set_yticks([0,1])
ax[2].set_yticklabels(['liq','ice'])
ax[3].plot(meas.utc,meas.ki)
ax[3].set_ylabel('$\\chi^{2}$')
ax[3].set_xlabel('UTC [Hours]')
ax[3].set_xlim([18.5,19.05])
plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results_v4.png',dpi=600)
#plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results_v3.pdf',bbox='tight')


# Smooth out the retrieved 4STAR data

# In[42]:

meas.tau[meas.good] = smooth(meas.tau[meas.good],20)
meas.ref[meas.good] = smooth(meas.ref[meas.good],20)


# ## Get SSFR data from ER2

# In[43]:

import load_modis as lm
if 'lm' in locals():
    reload(lm)
from load_modis import load_ict


# ### First load the ER2 files

# In[44]:

ssfr_er2_file = fp+'er2/20130913/SEAC4RS-SSFR_ER2_20130913_R0.ict'
ssfr_er2 = load_ict(ssfr_er2_file)
print len(ssfr_er2['UTC'])


# In[45]:

nasdat_er2_file = fp+'er2/20130913/seac4rs-nasdat_er2_20130913_r0.ict'
er2 = load_ict(nasdat_er2_file)
print len(er2['Start_UTC'])
print er2['Start_UTC']
print np.min(er2['Solar_Zenith_Angle'])


# In[141]:

plt.figure()
feet2meter=0.3048
plt.plot(er2['Start_UTC'],er2['GPS_Altitude'],label="GPS",color='b')
plt.plot(er2['Start_UTC'],er2['Pressure_Altitude']*feet2meter,label="Pressure",color='g')
plt.grid(True)
plt.xlabel('Time [UTC]')
plt.ylabel('ER-2 Altitude [m]')
plt.legend(frameon=False)
plt.savefig(fp+'plots\\20130913_er2_alt.png',dpi=600,transparent=True)


# In[46]:

print er2['Start_UTC'][12000]
print er2['GPS_Altitude'][12000]
print er2['Pressure_Altitude'][12000]


# In[145]:

plt.figure()
feet2meter=0.3048
plt.plot(er2['Start_UTC'],er2['Roll_Angle'],'b.',label="Roll")
plt.plot(er2['Start_UTC'],er2['Pitch_Angle'],'g.',label="Pitch")
plt.grid(True)
plt.xlabel('Time [UTC]')
plt.ylabel('ER-2 Atitude [$^{\circ}$]')
plt.legend(frameon=False)
plt.savefig(fp+'plots\\20130913_er2_roll_pitch.png',dpi=600,transparent=True)


# ### Now load the SSFR files from the ER2

# In[47]:

ssfr_idl_file = fp+'er2/20130913/20130913_calibspcs.out'
ssfr_idl = sio.idl.readsav(ssfr_idl_file)
print ssfr_idl.keys()
print np.shape(ssfr_idl['zspectra'])
print np.shape(ssfr_idl['sat'])


# In[48]:

class object():
    pass
ssfr = object()
if False:
    ssfr.utc = ssfr_er2['UTC']
    ssfr.Rvis = ssfr_er2['UP500']/ssfr_er2['DN500']
    ssfr.Rnir = ssfr_er2['UP1600']/ssfr_er2['DN1600']
    ssfr.Rvis[ssfr_er2['UP500']==-999] = np.nan
    ssfr.Rnir[ssfr_er2['UP500']==-999] = np.nan
    ssfr.Rnir[ssfr_er2['UP500']<=0] = np.nan
    ssfr.Rvis[ssfr_er2['UP500']<=0] = np.nan
    print np.max(er2['Start_UTC'])
    print np.max(ssfr.utc)
else:
    ssfr.utc = ssfr_idl['tmhrs']
    i500 = np.argmin(abs(ssfr_idl['zenlambda']-500.0))
    i1650 = np.argmin(abs(ssfr_idl['zenlambda']-1700.0))
    ssfr.Rvis = ssfr_idl['zspectra'][:,i500]/ssfr_idl['nspectra'][:,i500]
    ssfr.Rnir = ssfr_idl['zspectra'][:,i1650]/ssfr_idl['nspectra'][:,i1650-1]
    ssfr.Rvis[ssfr.Rvis<0] = np.nan
    ssfr.Rvis[ssfr.Rvis>1] = np.nan
    ssfr.Rnir[ssfr.Rnir<0] = np.nan
    ssfr.Rnir[ssfr.Rnir>1] = np.nan


# In[49]:

print i500, i1650
i1600 = np.argmin(abs(ssfr_idl['zenlambda']-1600.0))


# In[50]:

iu = np.argmin(abs(ssfr_idl['tmhrs']-18.5))
iu2 = np.argmin(abs(ssfr_er2['UTC']-18.5))
print iu,iu2


# In[51]:

print ssfr_er2['UP500'][iu2], ssfr_idl['zspectra'][iu,i500], ssfr_er2['UP500'][iu2]/ssfr_idl['zspectra'][iu,i500]
print ssfr_er2['DN500'][iu2+200], ssfr_idl['nspectra'][iu+200,i500], ssfr_er2['DN500'][iu2+200]/ssfr_idl['nspectra'][iu+200,i500]
print ssfr_er2['UTC'][iu2+200]
print ssfr_idl['tmhrs'][iu+200]
print ssfr_idl['nspectra'][iu+200,i500]
rup500 = ssfr_er2['UP500'][iu2]/ssfr_idl['zspectra'][iu,i500]
rdn500 = ssfr_er2['DN500'][iu2+200]/ssfr_idl['nspectra'][iu+200,i500]


# In[52]:

print ssfr_er2['UP1600'][iu2], ssfr_idl['zspectra'][iu,i1600], ssfr_er2['UP1600'][iu2]/ssfr_idl['zspectra'][iu,i1600]
print ssfr_er2['DN1600'][iu2+200], ssfr_idl['nspectra'][iu+200,i1600], ssfr_er2['DN1600'][iu2+200]/ssfr_idl['nspectra'][iu+200,i1600]
print ssfr_er2['UTC'][iu2+200]
print ssfr_idl['tmhrs'][iu+200]
print ssfr_idl['nspectra'][iu+200,i500]
rup1600 = ssfr_er2['UP1600'][iu2]/ssfr_idl['zspectra'][iu,i1600]
rdn1600 = ssfr_er2['DN1600'][iu2+200]/ssfr_idl['nspectra'][iu+200,i1600]


# In[53]:

print rup500/rdn500


# In[54]:

ssfr.Rvis = ssfr.Rvis*rup500/rdn500
ssfr.Rnir = ssfr.Rnir*rup1600/rdn1600


# In[55]:

print np.nanmin(ssfr.Rvis), np.nanmax(ssfr.Rvis)
print ssfr_idl['nadlambda'][i1650-2]
print i500, i1650
print len(ssfr.utc), len(ssfr.Rvis)


# In[56]:

from scipy import interpolate
sza_fx = interpolate.interp1d(er2['Start_UTC'],er2['Solar_Zenith_Angle'],bounds_error=False)
ssfr.sza = sza_fx(ssfr.utc)
print sza_fx(18.41)
print ssfr.utc[17000], ssfr.sza[17000]


# In[28]:

fig = plt.figure()
fig.add_subplot(1,2,1)
plt.plot(ssfr.utc,ssfr_idl['zspectra'][:,i500],label='Zenith 500 nm')
plt.plot(ssfr.utc,ssfr_idl['nspectra'][:,i500],'r', label='Nadir 500 nm')
plt.title('SSFR Irradiance')
plt.xlabel('UTC [h]')
plt.ylabel('Irradiance [$Wm^{-2}nm^{-1}sr^{-1}$]')
plt.legend(frameon=False)
plt.xlim([17.8,19.2])
fig.add_subplot(1,2,2)
plt.plot(ssfr.utc,ssfr.Rvis,'b',ssfr.utc,ssfr.Rnir,'r')
plt.xlim([17.8,19.2])
#plt.ylim([0.22,0.32])


# ##### Check the lut for reflectance

# In[57]:

lut.sp_hires(doirrad=True)


# In[58]:

lut.reflect.shape


# In[59]:

lut.ref


# In[80]:

figref = plt.figure(15,16)
figref,axref = plt.subplots(10,3,sharex=True,figsize=(15,16))
axref = axref.ravel()
for i in range(lut.ref.size/2-1):
    color.cycle_cmap(len(lut.tau),cmap=plt.cm.gist_ncar,ax=axref[i])
    for j in xrange(len(lut.tau)):
        axref[i].plot(lut.wvl,lut.reflect[1,:,1,2*i,j])
    axref[i].set_title('$r_{eff}$=%2.0f $\\mu$m' % lut.ref[i*2+1])
    axref[i].grid()
    axref[i].set_xlim([350,1700])
    axref[i].set_ylim([0,1])
    if i > 24: 
        axref[i].set_xlabel('Wavelength [nm]')

ididnt = [u for u in range(axref.size) if u>i]
for u in ididnt: 
    axref[u].axis('off')
for u in range(len(ididnt)):
    for tk in axref[i-u-1].get_xticklabels():
        tk.set_visible(True)
    
figref.tight_layout()
plt.suptitle('Ice')
plt.subplots_adjust(top=0.93,right=0.93)
cbar_ax = figref.add_axes([0.95,0.10,0.02,0.8])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.tau)
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
cba.ax.set_ylabel('$\\tau$')
labels = ['%5.1f' %F for F in np.linspace(lut.tau[0],lut.tau[-1],6)]
cba.ax.set_yticklabels(labels);
plt.savefig(fp+'plots/Spectral_reflectance_model_ice.png',dpi=600,transparent=True)
#plt.show()


# In[64]:

figref,axref = plt.subplots(10,3,sharex=True,figsize=(15,16))
axref = axref.ravel()
for i in range(lut.ref.size/2-1):
    color.cycle_cmap(len(lut.tau),cmap=plt.cm.gist_ncar,ax=axref[i])
    for j in xrange(len(lut.tau)):
        axref[i].plot(lut.wvl,lut.reflect[0,:,1,2*i,j])
    axref[i].set_title('$r_{eff}$=%2.0f $\\mu$m' %lut.ref[i*2+1])
    axref[i].grid()
    axref[i].set_xlim([350,1700])
    axref[i].set_ylim([0,1])
    if i > 23: 
        axref[i].set_xlabel('Wavelength [nm]')
        
ididnt = [u for u in range(axref.size) if u>i]
for u in ididnt: 
    axref[u].axis('off')
for u in range(len(ididnt)):
    for tk in axref[i-u-1].get_xticklabels():
        tk.set_visible(True)
        
figref.tight_layout()
plt.suptitle('Liquid')
plt.subplots_adjust(top=0.93,right=0.93)
cbar_ax = figref.add_axes([0.95,0.10,0.02,0.8])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.tau)
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
cba.ax.set_ylabel('$\\tau$')
labels = ['%5.1f' %F for F in np.linspace(lut.tau[0],lut.tau[-1],6)]
cba.ax.set_yticklabels(labels);
plt.savefig(fp+'plots/Spectral_reflectance_model_liq.png',dpi=600,transparent=True)
plt.show()


# ### plot the lut

# In[60]:

w500 = np.argmin(abs(lut.wvl-500.0))
w1700 = np.argmin(abs(lut.wvl-1700.0))


# In[61]:

lut.reflect.shape


# In[79]:

plt.figure()
taus = [1,2,3,4,6,8,10,15,20,30,50,100]
refs = [1,5,10,15,20,25,30,35,40,50,60]
for ir,r in enumerate(lut.ref):
    if r in refs: 
        plt.plot(lut.reflect[1,w500,1,ir,:],lut.reflect[1,w1700,1,ir,:],'k-')
        plt.annotate('%2i $\mu$m'%r,xy=(lut.reflect[1,w500,1,ir,-1]+0.01,lut.reflect[1,w1700,1,ir,-1]-0.01))
for it,t in enumerate(lut.tau):
    if t in taus:
        plt.plot(lut.reflect[1,w500,1,:,it],lut.reflect[1,w1700,1,:,it],'k--')
        if t<6:
            plt.annotate('$\\tau$=%2i'%t,xy=(lut.reflect[1,w500,1,-1,it]-0.02,lut.reflect[1,w1700,1,-1,it]-0.025))
        else:
            plt.annotate('%2i'%t,xy=(lut.reflect[1,w500,1,-1,it]-0.01,lut.reflect[1,w1700,1,-1,it]-0.025))
plt.xlabel('Reflectance at 500 nm') 
plt.ylabel('Reflectance at 1700 nm')
plt.xlim([0.45,1.05])
plt.ylim([0.25,0.8])
plt.savefig(fp+'plots/Reflectance_lut_ice.png',dpi=600,transparent=True)


# ### Now run the retrieval on reflectance from SSFR measurements based on the ER2 platorm

# In[62]:

import run_2wvl_retrieval as rw
if 'rw' in locals():
    reload(rw)


# In[63]:

(ssfr.tau,ssfr.ref,ssfr.ki) = rw.run_2wvl_retrieval(ssfr,lut,wvls=[500.0,1600.0])


# In[64]:

ssfr.tau[ssfr.ref==60] = np.nan
ssfr.ref[ssfr.ref==60] = np.nan


# In[371]:

figsr,axsr = plt.subplots(3,1,sharex=True)
axsr[0].plot(ssfr.utc,ssfr.tau[:,1])
axsr[1].plot(ssfr.utc,ssfr.ref[:,1])
axsr[2].plot(ssfr.utc,ssfr.ki[:,1])
axsr[2].set_xlabel('Wavelength [nm]')
axsr[0].set_ylabel('$\\tau$')
axsr[1].set_ylabel('R$_{e}$ [$\\mu$m]')
axsr[1].set_ylim([0,61])
axsr[2].set_ylabel('$\\chi^{2}$')
axsr[0].set_title('SSFR ER-2 Retrieval results')
axsr[2].set_xlim([18.0,19.2])


# ## Load data from CPL to compare flight profiles of DC8 and ER2 to cloud layers

# In[65]:

cpl_layer_file = fp+'er2\\20130913\\layers_13965_13sep13.txt'
import load_modis as lm
reload(lm)
from load_modis import load_cpl_layers
cpl_layers = load_cpl_layers(cpl_layer_file)


# In[66]:

cpl_layers.dtype


# In[58]:

fig,ax = plt.subplots(1,1)
ax2 = ax.twinx()
ax2.set_yticks([0,1000,2000,3000])
ax2.set_yticklabels(['None','PBL','Aerosol','Cloud'])
ax2.set_ylim([-200,20000])
for i in xrange(10):
    ax.vlines(cpl_layers['utc'],cpl_layers['bot'][:,i],cpl_layers['top'][:,i],lw=0.1,color=plt.cm.gist_rainbow(i/10.0))
    ax2.plot(cpl_layers['utc'],cpl_layers['type'][:,i]*1000,'+',color=plt.cm.gist_rainbow(i/10.0))
plt.title('Cloud layers from CPL')
ax.set_xlabel('UTC [H]')
ax.set_ylabel('Altitude [m]')
#plt.tick_params(axis='y', which='both', labelleft='on', labelright='on')
plt.savefig(fp+'plots/20130913_cpl_layers.png',dpi=600,transparent=True)


# In[69]:

fig,ax = plt.subplots(1,1)
ax2 = ax.twinx()
ax2.set_yticks([0,1000,2000,3000])
ax2.set_yticklabels(['None','PBL','Aerosol','Cloud'])
ax2.set_ylim([-200,20000])
for i in xrange(10):
    ax.vlines(cpl_layers['utc'],cpl_layers['bot'][:,i],cpl_layers['top'][:,i],lw=0.1,color=plt.cm.gist_rainbow(i/10.0))
    ax2.plot(cpl_layers['utc'],cpl_layers['type'][:,i]*1000,'+',color=plt.cm.gist_rainbow(i/10.0))
plt.title('Cloud layers from CPL')
ax.set_xlabel('UTC [H]')
ax.set_ylabel('Altitude [m]')
ax.set_xlim([18.5,19.2])
#plt.tick_params(axis='y', which='both', labelleft='on', labelright='on')
plt.savefig(fp+'plots/20130913_cpl_layers_zoom.png',dpi=600,transparent=True)


# In[48]:

fig,ax = plt.subplots(1,1)
ax2 = ax.twinx()
ax2.set_yticks([0,1000,2000,3000])
ax2.set_yticklabels(['None','PBL','Aerosol','Cloud'])
ax2.set_ylim([-200,20000])
for i in xrange(10):
    ax.vlines(cpl_layers['utc'],cpl_layers['bot'][:,i],cpl_layers['top'][:,i],lw=0.1,color=plt.cm.gist_rainbow(i/10.0))
    ax2.plot(cpl_layers['utc'],cpl_layers['type'][:,i]*1000,'+',color=plt.cm.gist_rainbow(i/10.0))
plt.title('Cloud layers from CPL')
ax.set_xlabel('UTC [H]')
ax.set_ylabel('Altitude [m]')
ax.set_xlim([18.0,19.5])
ax.plot(er2['Start_UTC'],er2['GPS_Altitude'],label="ER-2",color='b')
ax.plot(dc8['TIME_UTC'],dc8['G_ALT'],label="DC8",color='k')
ax.legend(frameon=False)
#plt.tick_params(axis='y', which='both', labelleft='on', labelright='on')
plt.savefig(fp+'plots/20130913_cpl_layers_zoom_flightpath.png',dpi=600,transparent=True)


# In[191]:

fig,ax = plt.subplots(1,1)
#ax2 = ax.twinx()
#ax2.set_yticks([0,1000,2000,3000])
#ax2.set_yticklabels(['None','PBL','Aerosol','Cloud'])
#ax2.set_ylim([-200,20000])
for i in xrange(10):
    ax.vlines(cpl_layers['utc'],cpl_layers['bot'][:,i],cpl_layers['top'][:,i],lw=0.1,color=plt.cm.Greys(200))
    #ax2.plot(cpl_layers['utc'],cpl_layers['type'][:,i]*1000,'+',color=plt.cm.gist_rainbow(i/10.0))
plt.title('Cloud layers from CPL')
ax.set_xlabel('UTC [H]')
ax.set_ylabel('Altitude [m]')
ax.set_xlim([18.0,19.5])
ax.set_ylim([0,25000])
ax.plot(er2['Start_UTC'],er2['GPS_Altitude'],label="ER-2",color='r')
ax.plot(dc8['TIME_UTC'],dc8['G_ALT'],label="DC8",color='b')
ax.legend(frameon=False)
#plt.tick_params(axis='y', which='both', labelleft='on', labelright='on')
plt.savefig(fp+'plots/20130913_cpl_layers_zoom_flightpath_grey.png',dpi=600,transparent=True)


# In[67]:

ia = abs(cpl_layers['utc']-19.03).argmin()
print cpl_layers['top'][ia,0]
print cpl_layers['bot'][ia,:]


# ## Get the data from MODIS to compare

# In[68]:

from mpl_toolkits.basemap import Basemap,cm
myd06_file = fp+'modis\\20130913\\MYD06_L2.A2013256.1910.006.2014267222159.hdf'
myd03_file = fp+'modis\\20130913\\MYD03.A2013256.1910.006.2013257153624.hdf'
print os.path.isfile(myd03_file) #check if it exists
print os.path.isfile(myd06_file)


# In[69]:

import load_modis as lm
reload(lm)
if 'modis' in locals():
    del modis, modis_dicts
    import gc; gc.collect()


# In[70]:

modis,modis_dicts = lm.load_modis(myd03_file,myd06_file)


# Now plot the resulting imagery

# In[75]:

#set up a easy plotting function
def seac_map(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=-89,lat_0=20,
            llcrnrlon=-98, llcrnrlat=18,
            urcrnrlon=-85, urcrnrlat=32,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(-85,-99,8),labels=[0,0,0,1])
    m.drawparallels(np.linspace(18,32,8),labels=[1,0,0,0])
    return m


# In[37]:

figm2,axm2 = plt.subplots(1,2,figsize=(13,10))
m1 = seac_map(axm2[0])
xt,yt = m1(-95.3831,29.7628)
axm2[0].text(xt,yt,'+')
axm2[0].text(xt,yt,'Houston, TX',horizontalalignment='right',verticalalignment='top')
xh,yh = m1(-95.3,19.1)
axm2[0].text(xh,yh,'+')
axm2[0].text(xh,yh,'Tropical Storm Ingrid',horizontalalignment='left',verticalalignment='bottom')
m2 = seac_map(axm2[1])
x,y = m1(modis['lon'],modis['lat'])
clevels = np.linspace(0,80,41)

cs1 = m1.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m1.colorbar(cs1)
cbar.set_label('$\\tau$')
axm2[0].set_title('MODIS - AQUA Cloud optical Thickness')
x1,y1 = m1(meas.lon,meas.lat)
xer2,yer2 = m1(er2['Longitude'],er2['Latitude'])
xdc8,ydc8 = m1(dc8['G_LONG'],dc8['G_LAT'])
#m1.plot(xer2,yer2,'r',lw=2.0)
#m1.plot(xdc8,ydc8,'b',lw=2.0)
xx1,yy1 = m1(-93.8,27.8)
xx2,yy2 = m1(-96.5,28)
#plt.text(xx1,yy1,'ER2',color='r')
#plt.text(xx2,yy2,'DC8',color='b')
m1.scatter(x1,y1,c=meas.tau,cmap=plt.cm.jet,marker='o',vmin=clevels[0],vmax=clevels[-1],alpha=0.5,edgecolors='k',linewidth=0.15)

clevels2 = np.linspace(0,60,31)
cs2 = m2.contourf(x,y,modis['ref'],clevels2,cmap=plt.cm.gist_earth,extend='max')
cbar = m2.colorbar(cs2)
cbar.set_label('R$_{eff}$ [$\\mu$m]')
axm2[1].set_title('MODIS - AQUA Cloud effective radius')
m2.plot(xer2,yer2,'r',lw=2)
m2.plot(xdc8,ydc8,'b',lw=2)
plt.text(xx1,yy1,'ER2',color='r')
plt.text(xx2,yy2,'DC8',color='b')
m2.scatter(x1,y1,c=meas.ref,cmap=plt.cm.gist_earth,marker='o',vmin=clevels2[0],vmax=clevels2[-1],alpha=0.5,edgecolors='k',linewidth=0.15)
figm2.subplots_adjust(wspace=0.3)
plt.savefig(fp+'plots/modis_dc8_tau_ref_comp.png',dpi=600,transparent=True)
#plt.savefig(fp+'plots/modis_dc8_tau_ref_comp.pdf',bbox='tight')
plt.show()


# ## Import eMAS values

# In[71]:

if 'lm' in locals():
    reload(lm)
from load_modis import load_emas, load_hdf


# In[72]:

emas_file = fp+'er2/20130913/EMASL2_13965_13_20130913_1905_1918_V00.hdf'
print os.path.isfile(emas_file)


# In[73]:

emas,emas_dicts = load_hdf(emas_file)


# In[74]:

emas_values = (('lat',0),('lon',1),('tau',15),('ref',23),('phase',58),('layer',59),('qa',68))
emas,emas_dicts = load_hdf(emas_file,values=emas_values)


# In[80]:

plt.figure()
plt.plot(emas['tau'])


# Now Redo the load of emas data, but with the new V01 files, which includes the newest calibration as of 20150122, which is considered final for SEAC4RS. thermal band revisions (3.7 um and higher) may still occur.  
# 
# From Tom Arnold 20150120: While the calibration for the V01 data is considered final, some minor revision may still be made to the themal band (3.7um and higher) radiance data for two or three of the tracks I have given you. Separate from the calibration process, filtering (for a coherent noise problem) has been applied to all the eMAS thermal bands (bands 26-38).  Evaluation of the quality of the filtering has shown for some eMAS tracks, some additional filtering is still required (and will likely affect two or three of the tracks I have given you).  I will make that data available when it is completed for the tracks you need, though for the thick cirrus in the tracks you are interested in,  I don’t expect much impact to the cloud retrievals (most impact will be relatively small adjustment to some of the cloud top property data - such as cloud top temperature or pressure). I expect the re-filtered data to be available in the next couple weeks.

# There is multiple eMAS files, representing each a different time slice. Load all of them.

# In[75]:

emas_file_v1_10 = fp+'emas/20130913/EMASL2_13965_10_20130913_1815_1828_V01.hdf'
emas_file_v1_11 = fp+'emas/20130913/EMASL2_13965_11_20130913_1832_1845_V01.hdf'
emas_file_v1_12 = fp+'emas/20130913/EMASL2_13965_12_20130913_1848_1901_V01.hdf'
emas_file_v1_14 = fp+'emas/20130913/EMASL2_13965_14_20130913_1921_1934_V01.hdf'
emas_file_v1 = fp+'emas/20130913/EMASL2_13965_13_20130913_1905_1918_V01.hdf'
print os.path.isfile(emas_file_v1)


# In[76]:

print fp
print emas_file_v1


# In[77]:

emas_v1,emas_dicts_v1 = load_hdf(emas_file_v1)


# In[78]:

emas_values = (('lat',0),('lon',1),('tau',15),('ref',23),('phase',58),('layer',59),('qa',68))
emas_v1,emas_dicts_v1 = load_hdf(emas_file_v1, values=emas_values)


# In[79]:

emas_v1_10,emas_dicts_v1_10 = load_hdf(emas_file_v1_10, values=emas_values, verbose=False)
emas_v1_11,emas_dicts_v1_11 = load_hdf(emas_file_v1_11, values=emas_values, verbose=False)
emas_v1_12,emas_dicts_v1_12 = load_hdf(emas_file_v1_12, values=emas_values, verbose=False)
emas_v1_14,emas_dicts_v1_14 = load_hdf(emas_file_v1_14, values=emas_values, verbose=False)


# In[80]:

emas_dicts_v1['tau']


# In[81]:

from map_utils import map_ind
dc8_ind = map_ind(emas['lon'],emas['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good'][0])


# Create different good filters for each different time slice

# In[82]:

mea['good_10'] = np.where((mea['utc']>18.15) & (mea['utc']<18.50) & (mea['Str'].flatten()!=0) & (mea['sat_time'].flatten()==0))[0]
mea['good_11'] = np.where((mea['utc']>18.50) & (mea['utc']<18.85) & (mea['Str'].flatten()!=0) & (mea['sat_time'].flatten()==0))[0]
mea['good_12'] = np.where((mea['utc']>18.75) & (mea['utc']<19.10) & (mea['Str'].flatten()!=0) & (mea['sat_time'].flatten()==0))[0]
mea['good_13'] = np.where((mea['utc']>19.00) & (mea['utc']<19.35) & (mea['Str'].flatten()!=0) & (mea['sat_time'].flatten()==0))[0]
mea['good_14'] = np.where((mea['utc']>19.25) & (mea['utc']<19.60) & (mea['Str'].flatten()!=0) & (mea['sat_time'].flatten()==0))[0]


# In[83]:

mea['good_13']


# In[84]:

print mea['Lat'][mea['good_13']].min(), mea['Lat'][mea['good_13']].max()
print mea['Lon'][mea['good_13']].min(), mea['Lon'][mea['good_13']].max()


# In[85]:

print emas_v1['lat'].min(), emas_v1['lat'].max()
print emas_v1['lon'].min(), emas_v1['lon'].max()


# In[86]:

import map_utils as mu
reload(mu)


# In[87]:

dc8_ind_10 = mu.map_ind(emas_v1_10['lon'],emas_v1_10['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good_10'])
dc8_ind_11 = mu.map_ind(emas_v1_11['lon'],emas_v1_11['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good_11'])
dc8_ind_12 = mu.map_ind(emas_v1_12['lon'],emas_v1_12['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good_12'])
dc8_ind_13 = mu.map_ind(emas_v1['lon'],emas_v1['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good_13'])
dc8_ind_14 = mu.map_ind(emas_v1_14['lon'],emas_v1_14['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good_14'])


# In[258]:

print np.shape(dc8_ind)
print dc8_ind_10.shape
print dc8_ind_11.shape
print dc8_ind_12.shape
print dc8_ind_14.shape


# In[88]:

print dc8_ind[0,-1],dc8_ind[1,-1]
print emas['lon'][388,715],emas['lat'][388,715]
print emas['lon'][388,714],emas['lat'][388,714]

sdist= lambda lon1,lat1,lon2,lat2:1000.0 * 3958.75 * np.arccos(np.cos(np.radians(lat1)-np.radians(lat2)) - np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * (1 - np.cos(np.radians(lon1)-np.radians(lon2))))


# In[89]:

print '%20.17f' % sdist(emas['lon'][388,715],emas['lat'][388,715],emas['lon'][385,714],emas['lat'][385,714])


# ### Now combine all eMAS results

# In[90]:

emas_full = dict()
emas_full['lon'] = np.concatenate([emas_v1_10['lon'][dc8_ind_10[0,:],dc8_ind_10[1,:]],
                    emas_v1_11['lon'][dc8_ind_11[0,:],dc8_ind_11[1,:]],
                    emas_v1_12['lon'][dc8_ind_12[0,:],dc8_ind_12[1,:]],
                    emas_v1['lon'][dc8_ind_13[0,:],dc8_ind_13[1,:]],
                    emas_v1_14['lon'][dc8_ind_14[0,:],dc8_ind_14[1,:]]])
emas_full['lat'] = np.concatenate([emas_v1_10['lat'][dc8_ind_10[0,:],dc8_ind_10[1,:]],
                    emas_v1_11['lat'][dc8_ind_11[0,:],dc8_ind_11[1,:]],
                    emas_v1_12['lat'][dc8_ind_12[0,:],dc8_ind_12[1,:]],
                    emas_v1['lat'][dc8_ind_13[0,:],dc8_ind_13[1,:]],
                    emas_v1_14['lat'][dc8_ind_14[0,:],dc8_ind_14[1,:]]])
emas_full['tau'] = np.concatenate([emas_v1_10['tau'][dc8_ind_10[0,:],dc8_ind_10[1,:]],
                    emas_v1_11['tau'][dc8_ind_11[0,:],dc8_ind_11[1,:]],
                    emas_v1_12['tau'][dc8_ind_12[0,:],dc8_ind_12[1,:]],
                    emas_v1['tau'][dc8_ind_13[0,:],dc8_ind_13[1,:]],
                    emas_v1_14['tau'][dc8_ind_14[0,:],dc8_ind_14[1,:]]])
emas_full['ref'] = np.concatenate([emas_v1_10['ref'][dc8_ind_10[0,:],dc8_ind_10[1,:]],
                    emas_v1_11['ref'][dc8_ind_11[0,:],dc8_ind_11[1,:]],
                    emas_v1_12['ref'][dc8_ind_12[0,:],dc8_ind_12[1,:]],
                    emas_v1['ref'][dc8_ind_13[0,:],dc8_ind_13[1,:]],
                    emas_v1_14['ref'][dc8_ind_14[0,:],dc8_ind_14[1,:]]])


# Do the mods calculations

# In[91]:

from map_utils import map_ind
dc8_ind_modis = map_ind(modis['lon'],modis['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good'][0])


# ## Load APR-2 HDF files for DC8 radar images

# In[92]:

fa = fp+'dc8/20130913//SEAC4RS-APR2_DC8_20130913/SEAC4RS-APR2_DC8_20130913'
fe = '_R23.h4'
files = ['180527','181019','182329','184933','190145','192149','194031']
aprfiles = [fa+s+fe for s in files]


# In[93]:

aprfiles


# In[94]:

reload(load_modis)


# In[95]:

from load_modis import load_apr
apr = load_apr(aprfiles)


# In[40]:

levels = np.arange(0,7000,30)
plt.figure()
csz = plt.contourf(apr['lonz'],apr['altflt'],apr['dbz'],levels,cmap=plt.cm.jet)
plt.colorbar(csz)


# In[96]:

apr['utcz'] = apr['lonz']*0.0
for i in xrange(len(apr['utc'])):
    apr['utcz'][:,i] = apr['utc'][i]


# In[97]:

levels = np.arange(0,7000,30)


# In[43]:

plt.figure()
csz = plt.contourf(apr['utcz'],apr['altflt'],apr['dbz'],levels,cmap=plt.cm.Greys)
plt.colorbar(csz)
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')


# In[45]:

plt.figure()
csz = plt.contourf(apr['utcz'],apr['altflt'],apr['dbz'],levels,cmap=plt.cm.jet)
apr_t = [18.1827,18.3147,18.4175,18.4416,18.4749,18.499,18.5524,18.6285,18.6405,18.6633,18.6886,18.7327,
         18.9308,19.0749,19.1495,19.1655,19.1793,19.1976,19.2137,19.271,19.4728,19.4888,19.536,19.5851,19.6085,19.6151,19.7116,19.738,19.7753,19.8368]
plt.colorbar(csz)
for t in apr_t:
    plt.axvline(t)
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')


# In[46]:

plt.figure()
plt.plot(apr['lonz'][11,:],apr['latz'][11,:])
color.cycle_cmap(len(apr_t),cmap=plt.cm.jet,ax=plt.gca())
for t in apr_t:
    ii = np.abs(apr['utc']-t).argmin()
    plt.plot(apr['lonz'][11,ii],apr['latz'][11,ii],'o')
    plt.annotate('%f'%t,(apr['lonz'][11,ii],apr['latz'][11,ii]))


# In[52]:

plt.figure()
plt.plot(dc8['G_LONG'],dc8['G_LAT'],'b-')
plt.plot(apr['lonz'][11,:],apr['latz'][11,:],'k-',linewidth=4)
color.cycle_cmap(len(apr_t),cmap=plt.cm.jet,ax=plt.gca())
for t in apr_t[:-6]:
    ii = np.abs(apr['utc']-t).argmin()
    plt.plot(apr['lonz'][11,ii],apr['latz'][11,ii],'o')
    plt.annotate(t2str(t),(apr['lonz'][11,ii],apr['latz'][11,ii]),fontsize=18,color='r')
plt.xlabel('Longitude [$^{\circ}$]')
plt.ylabel('Latitude [$^{\circ}$]')


# In[98]:

def t2str(t):
    'Simple function to transform decimal time to string of time HH:MM'
    hh = int(t)
    mm = int((t-hh)*60)
    return '%2i:%02i'% (hh,mm)


# In[99]:

apr_t[:-6]


# In[ ]:

for t in apr_t[:-6]:
    print t2str(t)


# In[ ]:

it = np.abs(apr['utc']-19.1).argmin()


# In[67]:

plt.figure()
plt.plot(apr['dbz'][:,it],apr['altflt'][:,it])
plt.xlim([0,8000])
plt.xlabel('dbZ')
plt.ylabel('Altitude [m]')
plt.title('APR2 Zenith radar reflectivity at: %f2.2',apr['utc'][it])


# In[163]:

plt.figure
plt.plot(apr['dbz'][:,0],apr['altflt'][:,0],'+')
plt.ylim([6000,8000])


# In[100]:

inoisezone = apr['dbz'][:,2000].argmax()
inoisezone


# In[101]:

apr['dbz'][:,]


# In[172]:

plt.figure
plt.plot(apr['dbz'][:,1000],apr['altflt'][:,1000],'+')
plt.ylim([6000,10000])


# In[160]:

apr['altflt'][:,0]


# In[213]:

plt.figure
plt.plot(apr['dbz'][11,:,0],apr['altz'][11,:,0]+apr['alt'][11,0])
plt.xlim([0,8000])
plt.ylim([0,20000])


# In[116]:

fig,ax = plt.subplots(1,1)
ax2 = ax.twinx()
ax2.set_yticks([0,1000,2000,3000])
ax2.set_yticklabels(['None','PBL','Aerosol','Cloud'])
ax2.set_ylim([-200,20000])
csz = ax.contourf(apr['utcz'],apr['altflt'],apr['dbz'],levels,cmap=plt.cm.Greys)
for i in xrange(10):
    ax.vlines(cpl_layers['utc'],cpl_layers['bot'][:,i],cpl_layers['top'][:,i],lw=0.1,color=plt.cm.gist_rainbow(i/10.0))
    ax2.plot(cpl_layers['utc'],cpl_layers['type'][:,i]*1000,'+',color=plt.cm.gist_rainbow(i/10.0))
plt.title('Cloud layers from CPL')
ax.set_xlabel('UTC [H]')
ax.set_ylabel('Altitude [m]')
ax.set_xlim([18.0,19.5])
ax.plot(er2['Start_UTC'],er2['GPS_Altitude'],label="ER-2",color='b')
ax.plot(dc8['TIME_UTC'],dc8['G_ALT'],label="DC8",color='k')
ax.legend(frameon=False)


# In[102]:

it = np.abs(apr['utc']-19.19).argmin()
inoisezone = apr['dbz'][:,it].argmax()
i=range(inoisezone-15)+range(inoisezone+15,len(apr['altflt'][:,0]))


# In[103]:

print 0.19*60
print 0.10*60


# In[104]:

apr['utc'][it]


# In[152]:

plt.figure
plt.plot(smooth(apr['dbz'][i,it],10),apr['altflt'][i,it],'k')
plt.xlim([0,8000])
plt.xlabel('dbZ')
plt.ylabel('Altitude [m]')
#plt.title('APR2 Zenith radar reflectivity at: ',apr['utc'][it])


# In[268]:

plt.figure()
for j in xrange(20):
    plt.plot(apr['dbz'][i,it+j],apr['altflt'][i,it+j],label='%f.2'%apr['utc'][it+j])
plt.legend(frameon=False)


# ## Load the cloud probe data

# In[105]:

prb_file = fp+'dc8/20130913/SEAC4RS_20130913_Reff.txt'
probes = np.genfromtxt(prb_file,skip_header=2)


# In[106]:

print probes[:,1]
print probes[:,2]


# In[107]:

print probes.shape
print probes[:,7]
plt.figure()
plt.hist(probes[:,7])


# ## Load 2DS data for effective radius at specific times

# In[108]:

twoDS = load_ict(fp+'dc8/20130913/seac4rs-2DS_DC8_20130913_R0.ict')


# In[61]:

plt.figure
plt.plot(twoDS['Start_UTC'],twoDS['effectiveD']/2.0)
plt.plot(twoDS['Start_UTC'],smooth(twoDS['effectiveD']/2.0,60),'r')
plt.ylim([0,100])
plt.xlabel('UTC [h]')
plt.ylabel('Effective radius [$\\mu$m]')
plt.title('2DS effective radius from 2013-09-13')


# In[62]:

plt.figure
plt.plot(twoDS['effectiveD'][1:]/2.0,dc8['G_ALT'])
plt.xlim([0,100])
plt.xlabel('Effective radius [$\\mu$m]')
plt.ylabel('Altitude [m]')


# In[109]:

# filter the data and only show a part
twoDS['effectiveD'][twoDS['effectiveD']<0] = np.nan
fl = np.where((twoDS['Start_UTC'] > 18.025) & (twoDS['Start_UTC'] < 19.45) & (twoDS['effectiveD'] > 28.0))


# In[100]:

plt.figure
plt.plot(twoDS['effectiveD'][fl]/2.0,dc8['G_ALT'][fl],label='Cloud Probes (2DS)')
plt.xlim([0,100])
plt.xlabel('Effective radius [$\\mu$m]')
plt.ylabel('Altitude [m]')


# In[102]:

plt.figure
plt.plot(smooth(twoDS['effectiveD'][fl]/2.0,10),dc8['G_ALT'][fl],label='Cloud Probes (2DS)')
plt.xlim([0,100])
plt.xlabel('Effective radius [$\\mu$m]')
plt.ylabel('Altitude [m]')


# ## Load CCN results

# In[110]:

ccn,ccnhead = load_ict(fp+'dc8/20130913/SEAC4RS-CCN_DC8_20130913_R0.ict',return_header=True)


# In[673]:

ccnhead


# In[111]:

ccn['Number_Concentration'][ccn['Number_Concentration']==-8888]=np.nan
ccn['Number_Concentration_STP'][ccn['Number_Concentration_STP']==-8888]=np.nan
ccn['Supersaturation'][ccn['Supersaturation']==-999]=np.nan


# In[686]:

fig,ax = plt.subplots(3,1,sharex=True)
ax[0].plot(ccn['UTC_mid'],ccn['Number_Concentration'])
ax[1].plot(ccn['UTC_mid'],ccn['Number_Concentration_STP'])
ax[2].plot(ccn['UTC_mid'],ccn['Supersaturation'])
ax[0].set_title('Number Concentration')
ax[1].set_title('Number Concentration STP')
ax[2].set_title('SuperSaturation')
ax[2].set_xlabel('UTC [h]')
ax[0].set_ylabel('[\#/cm$^{3}$]')
ax[1].set_ylabel('[\#/cm$^{3}$]')
ax[2].set_ylabel('[\%]')


# In[687]:

fig,ax = plt.subplots(3,1,sharex=True)
ax[0].plot(ccn['UTC_mid'],ccn['Number_Concentration'])
ax[1].plot(ccn['UTC_mid'],ccn['Number_Concentration_STP'])
ax[2].plot(ccn['UTC_mid'],ccn['Supersaturation'])
ax[0].set_title('Number Concentration')
ax[1].set_title('Number Concentration STP')
ax[2].set_title('SuperSaturation')
ax[2].set_xlabel('UTC [h]')
ax[2].set_xlim([18.0,19.5])
ax[0].set_ylim([0,200])
ax[1].set_ylim([0,500])
ax[0].set_ylabel('[\#/cm$^{3}$]')
ax[1].set_ylabel('[\#/cm$^{3}$]')
ax[2].set_ylabel('[\%]')


# In[688]:

plt.plot(ccn['UTC_mid'],np.log(ccn['Number_Concentration']))
plt.xlim([18,19.5])


# In[112]:

from Sp_parameters import find_closest
id = find_closest(dc8['TIME_UTC'],ccn['UTC_mid'])


# In[113]:

ccn_good = np.where((ccn['UTC_mid']>18.0)&(ccn['UTC_mid']<19.5)& (np.isfinite(ccn['Number_Concentration'])))[0]


# In[692]:

plt.plot(dc8['G_LAT'][id[ccn_good]],ccn['Number_Concentration'][ccn_good],'b.')
plt.xlabel('Latitude')
plt.ylabel('Number Concentration [\#/cm$^{3}$]')
plt.title('CCN number concentration from DC8')


# ## Load RSP results

# In[114]:

rsp,rsp_header = load_ict(fp+'er2/20130913/SEAC4RS-RSP-ICECLD_ER2_20130913_R2.ict',return_header=True)
print rsp['Phase']
rsp_good = np.where((rsp['UTC']>17.8) & (rsp['UTC']<19.2) & (rsp['Phase']==0) & (rsp['COT']>0) & (rsp['R_eff159']>0))
print np.shape(rsp_good)
print len(rsp_good[0])


# In[302]:

rsp_header


# In[86]:

plt.plot(rsp['UTC'][rsp_good[0]],rsp['R_eff159'][rsp_good[0]],label='1.59 micron')
plt.plot(rsp['UTC'][rsp_good[0]],rsp['R_eff225'][rsp_good[0]],label='2.25 micron')
plt.xlabel('UTC [Hours]')
plt.ylabel('R$_{eff}$ [$\\mu$m]')
plt.title('RSP effective radius retrievals')
plt.legend(frameon=False)
plt.savefig(fp+'plots/20130913_RSP_reff.png',dpi=600,transparent=True)


# In[115]:

print rsp_good[0][-1]


# In[116]:

# get the distance between two points for RSP
from map_utils import spherical_dist
print rsp['Lat'][rsp_good[0][-1]], rsp['Lon'][rsp_good[0][-1]]
print rsp['Lat'][rsp_good[0][-2]], rsp['Lon'][rsp_good[0][-2]]
print spherical_dist(np.array((rsp['Lat'][rsp_good[0][-1]],rsp['Lon'][rsp_good[0][-1]])),
                     np.array((rsp['Lat'][rsp_good[0][-2]],rsp['Lon'][rsp_good[0][-2]])))


# ## Calculate the histogram for each comparisons

# In[117]:

from Sp_parameters import nanmasked
modis_tau,im = nanmasked(modis['tau'][dc8_ind_modis[0,:],dc8_ind_modis[1,:]])
emas_tau,ie = nanmasked(emas['tau'][dc8_ind[0,:],dc8_ind[1,:]])
emas_tau_v1,ie1 = nanmasked(emas_v1['tau'][dc8_ind[0,:],dc8_ind[1,:]])
emas_tau_full,ief = nanmasked(emas_full['tau'])
modis_ref,im = nanmasked(modis['ref'][dc8_ind_modis[0,:],dc8_ind_modis[1,:]])
emas_ref,ie = nanmasked(emas['ref'][dc8_ind[0,:],dc8_ind[1,:]])
emas_ref_v1,ie1 = nanmasked(emas_v1['ref'][dc8_ind[0,:],dc8_ind[1,:]])
emas_ref_full,ief = nanmasked(emas_full['ref'])
star_tau,ist = nanmasked(meas.tau[meas.good])
star_ref,ist = nanmasked(meas.ref[meas.good])
ssfr.good = np.where((ssfr.utc>17.8)&(ssfr.utc<19.2))
ssfr_tau,iss = nanmasked(ssfr.tau[ssfr.good[0],1])
ssfr_ref,iss = nanmasked(ssfr.ref[ssfr.good[0],1])
rsp_tau,irs = nanmasked(rsp['COT'][rsp_good[0]])
rsp_ref,irs = nanmasked(rsp['R_eff159'][rsp_good[0]])


# In[119]:

star_g = np.where((meas.utc>19.0)&(meas.utc<19.2)&np.isfinite(meas.tau))
print '4STAR',meas.tau[star_g[0][0]], meas.ref[star_g[0][0]]
ssfr_g = np.where((ssfr.utc>19.0)&(ssfr.utc<19.2)&np.isfinite(ssfr.tau[:,1]))
print 'ssfr reflect',ssfr.tau[ssfr_g[0][0],1],ssfr.ref[ssfr_g[0][0],1]
rsp_g = np.where((rsp['UTC']>19.0)&(rsp['UTC']<19.2)&np.isfinite(rsp['COT']))
print 'rsp',rsp['COT'][rsp_g[0][0]],rsp['R_eff159'][rsp_g[0][0]]
print 'emas',emas['tau'][dc8_ind[0,0],dc8_ind[1,0]],emas['ref'][dc8_ind[0,0],dc8_ind[1,0]]
print 'emas_v1',emas_v1['tau'][dc8_ind[0,0],dc8_ind[1,0]],emas_v1['ref'][dc8_ind[0,0],dc8_ind[1,0]]
print 'modis',modis['tau'][dc8_ind_modis[0,0],dc8_ind_modis[1,0]],modis['ref'][dc8_ind_modis[0,0],dc8_ind_modis[1,0]] 


# In[91]:

def plot_median_mean(x,lbl=False,color='k'):
    "plot the vertical median and mean over a histogram, if lbl set to true, sets the labels of the lines,default color='k'"
    if lbl:
        llm = 'Mean'
        lld = 'Median'
    else:
        llm = None
        lld = None
    plt.axvline(np.nanmean(x),color='k',label=llm,c=color,lw=2)
    plt.axvline(np.median(x),color='k',linestyle='--',label=lld,c=color,lw=2)


# Now we present the histogram of retrieved cloud properties for each instrument
# The instruments have different cross sectional area at cloud height that they sample
# These are:, with the ratio of smoothing:
#  - MODIS : 500m : 6
#  - emas : 50m : 60
#  - SSFR(reflectance) : 3000m : 1
#  - RSP: 42 m : 70
#  - 4STAR: 70 m : 40

# ## Run through data and convolve to get the same area

# ER-2 Altitude: 19018 m
# DC8 Altitude: 8047 m
# Cloud base height: 13610 m
# Cloud top height: 16489 m
#  - Modis: 500 m
#  - eMAS: 50m or 2.5 mrad FOV ()
#  - RSP: 14 mrad FOV
#  - 4STAR: 2 degree, 0.0349 rads
#  - GOES: 1km
#  - SSFR: 180 degrees, Pi rads

# In[122]:

cld_base = 13610.0
cld_top = 16489.0
er2_hgt = 19018.0
dc8_hgt = 8047.0
r_eMAS = (er2_hgt - cld_top)*np.tan(0.0025)
r_RSP = (er2_hgt - cld_top)*np.tan(0.014)
r_SSFR = (er2_hgt - cld_top)*np.tan(np.pi/3)
r_4STAR = (cld_base - dc8_hgt)*np.tan(0.0349)


# In[123]:

print r_eMAS
print r_RSP
print r_SSFR
print r_4STAR


# ## Save the output of the retrievals and the other retrievals

# In[125]:

import datetime
m_dict = {'time_created':datetime.datetime.today(),
          'emas':{'tau':emas_tau_full,'ref':emas_ref_full},
          'modis':{'tau':modis_tau,'ref':modis_ref},
          'ssfr':{'tau':ssfr_tau,'ref':ssfr_ref},
          'rsp':{'tau':rsp_tau,'ref':rsp_ref},
          'star':{'tau':star_tau,'ref':star_ref,'utc':meas.utc[meas.good][ist]}}


# In[126]:

hs.savemat(fp+'20130913_retrieval_output.mat',m_dict)


# ### Optionally load the retrieval output for easier and faster plotting

# In[ ]:

m_dict = hs.loadmat(fp+'20130913_retrieval_output.mat')


# ## Plot histogram of different tau and ref comparison

# In[275]:

fig,ax = plt.subplots(1,figsize=(9,6))

import matplotlib.transforms as mtransforms
trans1 = mtransforms.blended_transform_factory(ax.transAxes, ax.transAxes)
trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

plt.fill_between([0,1],0,1,transform=trans1,color='#FFFFFF')
plt.fill_between([np.nanmean(smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2))-7.8,np.nanmean(smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2))+7.8],
                 0,1,transform=trans,color='g',edgecolor='g',hatch='x',linewidth=0.2,alpha=0.3)
plt.fill_between([np.nanmean(smooth(modis_tau,6))-2.0,np.nanmean(smooth(modis_tau,6))+2.0],0,1,transform=trans,
                 color='m',edgecolor='m',alpha=0.3,hatch='x',linewidth=0.2)
plt.fill_between([np.nanmean(smooth(rsp_tau,70))-2.0,np.nanmean(smooth(rsp_tau,70))+2.0],0,1,transform=trans,
                 color='c',edgecolor='c',alpha=0.3,hatch='/',linewidth=0.2)
plt.fill_between([np.nanmean(smooth(star_tau,40))-2.0,np.nanmean(smooth(star_tau,40))+2.0],0,1,transform=trans,
                 color='r',edgecolor='r',alpha=0.3,hatch='\\',linewidth=0.2)

plt.hist(smooth(modis_tau,6),bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='Modis (Reflected)',range=(0,40))
#plt.hist(smooth(emas_tau,60),bins=30, histtype='stepfilled', normed=True, color='b',alpha=0.6, label='eMAS (Reflected)',range=(0,40))
#plt.hist(smooth(emas_tau_v1,60),bins=30, histtype='stepfilled', normed=True, color='k',alpha=0.6, label='eMAS (Reflected)',range=(0,40))
plt.hist(smooth(emas_tau_full,60),bins=30, histtype='stepfilled', normed=True, color='k',alpha=0.6, label='eMAS (Reflected)',range=(0,40))
plt.hist(smooth(ssfr_tau,2),bins=20, histtype='stepfilled', normed=True, color='g',alpha=0.6, label='SSFR (Reflected)',range=(5,30))
plt.hist(smooth(rsp_tau,70),bins=30, histtype='stepfilled', normed=True, color='c',alpha=0.6, label='RSP (Reflected)',range=(0,40))
plt.hist(smooth(star_tau,40),bins=30, histtype='stepfilled', normed=True, color='r',alpha=0.6, label='4STAR (Transmitted)',range=(0,40))
plt.title('Optical Thickness histogram')
plt.ylabel('Normed probability')
plt.xlabel('$\\tau$')
plot_median_mean(smooth(modis_tau,6),color='m')
#plot_median_mean(smooth(emas_tau ,60),color='b')
#plot_median_mean(smooth(emas_tau_v1,60),color='k')
plot_median_mean(smooth(emas_tau_full,60),color='k')
plot_median_mean(smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2),color='g')
plot_median_mean(smooth(star_tau,40),color='r',lbl=True)
plot_median_mean(smooth(rsp_tau,70),color='c')

xr = ax.get_xlim()
yr = ax.get_ylim()
ax.add_patch(plt.Rectangle((0,0),0,0,color='none',edgecolor='g',hatch='x',linewidth=0.2,alpha=0.5,label='Horizontal variability'))
ax.set_xlim(xr)
ax.set_ylim(yr)

plt.legend(frameon=False)
plt.xlim([0,60])
plt.savefig(fp+'plots/hist_modis_4star_tau_v4_fill.png',dpi=600,transparent=True)
#plt.savefig(fp+'plots/hist_modis_4star_tau.pdf',bbox='tight')


# In[274]:

fig,ax = plt.subplots(1,figsize=(9,6))
import matplotlib.transforms as mtransforms
trans1 = mtransforms.blended_transform_factory(ax.transAxes, ax.transAxes)
trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

plt.fill_between([0,1],0,1,transform=trans1,color='#FFFFFF')
plt.fill_between([np.nanmean(ssfr_ref)-2.7,np.nanmean(ssfr_ref)+2.7],0,1,transform=trans,
                 color='g',edgecolor='g',hatch='x',linewidth=0.2,alpha=0.3)
plt.fill_between([np.nanmean(modis_ref)-1.7,np.nanmean(modis_ref)+1.7],0,1,transform=trans,
                 color='m',edgecolor='m',alpha=0.3,hatch='x',linewidth=0.2)
plt.fill_between([np.nanmean(rsp_ref)-1.7,np.nanmean(rsp_ref)+1.7],0,1,transform=trans,
                 color='c',edgecolor='c',alpha=0.3,hatch='/',linewidth=0.2)
plt.fill_between([np.nanmean(star_ref)-1.8,np.nanmean(star_ref)+1.8],0,1,transform=trans,
                 color='r',edgecolor='r',alpha=0.3,hatch='\\',linewidth=0.2)

#plt.axvspan(0,80,color='#FFFFFF')
#plt.axvspan(np.nanmean(modis_ref)-1.7,np.nanmean(modis_ref)+1.7,color='m',alpha=0.7)
#plt.axvspan(np.nanmean(ssfr_ref)-2.7,np.nanmean(ssfr_ref)+2.7,color='g',alpha=0.7)
#plt.axvspan(np.nanmean(rsp_ref)-1.7,np.nanmean(rsp_ref)+1.7,color='c',alpha=0.7)
#plt.axvspan(np.nanmean(star_ref)-1.8,np.nanmean(star_ref)+1.8,color='r',alpha=0.7)
plt.hist(smooth(modis_ref,6),bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='Modis (Reflected)',range=(0,59))
#plt.hist(smooth(emas_ref,60),bins=30, histtype='stepfilled', normed=True, color='b',alpha=0.6, label='eMAS (Reflected)',range=(0,59))
#plt.hist(smooth(emas_ref_v1,60),bins=30, histtype='stepfilled', normed=True, color='k',alpha=0.6, label='eMAS (Reflected)',range=(0,59))
plt.hist(smooth(emas_ref_full,60),bins=30, histtype='stepfilled', normed=True, color='k',alpha=0.6, label='eMAS (Reflected)',range=(0,59))
plt.hist(smooth(ssfr_ref,2),bins=30, histtype='stepfilled', normed=True, color='g',alpha=0.6, label='SSFR (Reflected)',range=(0,59))
plt.hist(smooth(rsp_ref,70),bins=30, histtype='stepfilled', normed=True, color='c',alpha=0.6, label='RSP (Reflected)',range=(0,79))
plt.hist(star_ref,bins=30, histtype='stepfilled', normed=True, color='r',alpha=0.6, label='4STAR (Transmitted)',range=(0,59))
plt.hist(probes[:,7],bins=20, histtype='stepfilled', normed=True, color='y',alpha=0.6, label='Cloud probes (In Situ)',range=(0,79))
plt.title('Cloud particle effective radius histogram')
plt.ylabel('Normed probability')
plt.xlabel('R$_{eff}$ [$\\mu$m]')
plot_median_mean(smooth(modis_ref,6),color='m')
#plot_median_mean(smooth(emas_ref,60),color='b')
#plot_median_mean(smooth(emas_ref_v1,60),color='k')
plot_median_mean(smooth(emas_ref_full,60),color='k')
plot_median_mean(smooth(ssfr_ref,2),color='g')
plot_median_mean(smooth(rsp_ref,70),color='c')
plot_median_mean(probes[:,7],color='y')
plot_median_mean(star_ref,lbl=True,color='r')

ax.add_patch(plt.Rectangle((0,0),0,0,color='none',edgecolor='g',hatch='x',linewidth=0.2,alpha=0.5,label='Horizontal variability'))

plt.legend(frameon=False,loc='upper right')
plt.xlim([10,80])
plt.ylim([0,0.3])
plt.savefig(fp+'plots/hist_modis_4star_ref_v4_fill.png',dpi=600,transparent=True)


# In[131]:

plt.figure()
plt.boxplot(smooth(modis_tau,6),vert=False,color='m')


# In[132]:

fig,(ax1,ax2) = plt.subplots(2,figsize=(9,6))

import matplotlib.transforms as mtransforms
trans1 = mtransforms.blended_transform_factory(ax1.transAxes, ax1.transAxes)
trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)

ax1.fill_between([0,1],0,1,transform=trans1,color='#FFFFFF')
ax1.fill_between([np.nanmean(smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2))-7.8,np.nanmean(smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2))+7.8],
                 0,1,transform=trans,color='g',edgecolor='g',hatch='x',linewidth=0.1,alpha=0.3)
ax1.fill_between([np.nanmean(smooth(modis_tau,6))-2.0,np.nanmean(smooth(modis_tau,6))+2.0],0,1,transform=trans,
                 color='m',edgecolor='m',alpha=0.3,hatch='x',linewidth=0.1)
ax1.fill_between([np.nanmean(smooth(rsp_tau,70))-2.0,np.nanmean(smooth(rsp_tau,70))+2.0],0,1,transform=trans,
                 color='c',edgecolor='c',alpha=0.3,hatch='/',linewidth=0.1)
ax1.fill_between([np.nanmean(smooth(star_tau,40))-2.0,np.nanmean(smooth(star_tau,40))+2.0],0,1,transform=trans,
                 color='r',edgecolor='r',alpha=0.3,hatch='\\',linewidth=0.1)

ax1.hist(smooth(modis_tau,6),bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='Modis (Reflected)',range=(0,40))
#plt.hist(smooth(emas_tau,60),bins=30, histtype='stepfilled', normed=True, color='b',alpha=0.6, label='eMAS (Reflected)',range=(0,40))
ax1.hist(smooth(emas_tau_v1,60),bins=30, histtype='stepfilled', normed=True, color='k',alpha=0.6, label='eMAS (Reflected)',range=(0,40))
ax1.hist(smooth(ssfr_tau,2),bins=20, histtype='stepfilled', normed=True, color='g',alpha=0.6, label='SSFR (Reflected)',range=(5,30))
ax1.hist(smooth(rsp_tau,70),bins=30, histtype='stepfilled', normed=True, color='c',alpha=0.6, label='RSP (Reflected)',range=(0,40))
ax1.hist(smooth(star_tau,40),bins=30, histtype='stepfilled', normed=True, color='r',alpha=0.6, label='4STAR (Transmitted)',range=(0,40))
ax1.set_title('Optical Thickness histogram')
ax1.set_ylabel('Normed probability')
ax1.set_xlabel('$\\tau$')
plot_median_mean(smooth(modis_tau,6),color='m')
#plot_median_mean(smooth(emas_tau ,60),color='b')
plot_median_mean(smooth(emas_tau_v1,60),color='k')
plot_median_mean(smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2),color='g')
plot_median_mean(smooth(star_tau,40),color='r',lbl=True)
plot_median_mean(smooth(rsp_tau,70),color='c')

xr = ax1.get_xlim()
yr = ax1.get_ylim()
ax1.add_patch(plt.Rectangle((0,0),0,0,color='none',edgecolor='g',hatch='x',linewidth=0.1,alpha=0.5,label='Horizontal variability'))
ax1.set_xlim(xr)
ax1.set_ylim(yr)

plt.legend(frameon=False)
ax1.set_xlim([0,60])


plt.savefig(fp+'plots/hist_modis_4star_tau_ref.png',dpi=600,transparent=True)


# ## Plot histogram of retrieved properties (tau and ref) in new 'bean' plots

# In[95]:

from plotting_utils import plot_vert_hist


# In[96]:

import plotting_utils
reload(plotting_utils)


# In[917]:

fig = plt.figure(figsize=(7,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,40],xlim=[-0.5,4.5])
ax1.set_ylabel('$\\tau$')
ax1.set_xticks([0,1,2,3,4])
ax1.set_xticklabels(['MODIS\n(reflected)','eMAS\n(reflected)','SSFR\n(reflected)','RSP\n(reflected)','4STAR\n(transmitted)'])
plot_vert_hist(fig,ax1,smooth(modis_tau,6),0,[0,40],legend=True,onlyhist=False,loc=2,color='m')
plot_vert_hist(fig,ax1,smooth(emas_tau_full,60),1,[0,40],legend=True,color='k')
plot_vert_hist(fig,ax1,smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2),2,[0,40],legend=True,color='g')
plot_vert_hist(fig,ax1,smooth(rsp_tau,70),3,[0,40],legend=True,color='c')
plot_vert_hist(fig,ax1,smooth(star_tau,40),4,[0,40],legend=True,color='r')
(_,caps,_) = ax1.errorbar([0,1,2,3,4],
             [np.nanmean(smooth(modis_tau,6)),
              np.nanmean(smooth(emas_tau_full,60)),
              np.nanmean(smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2)),
              np.nanmean(smooth(rsp_tau,70)),
              np.nanmean(smooth(star_tau,40))],yerr=[2,0,7.8,2,2],
             color='r',linestyle='.',linewidth=3,label='Variability within FOV',capsize=7)
for cap in caps:
    cap.set_markeredgewidth(3)
ax1.legend(frameon=False,numpoints=1)
plt.savefig(fp+'plots/vert_hist_tau_v4.png',dpi=600,transparent=True)


# In[922]:

fig = plt.figure(figsize=(8,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,80],xlim=[-0.5,5.5])
ax1.set_ylabel('r$_{eff}$ [$\\mu$m]')
ax1.set_xticks([0,1,2,3,4,5])
ax1.set_xticklabels(['MODIS\n(reflected)','eMAS\n(reflected)','SSFR\n(reflected)','RSP\n(reflected)','4STAR\n(transmitted)','Cloud probes\n(In Situ)'])
plot_vert_hist(fig,ax1,smooth(modis_ref,6),0,[0,80],legend=True,onlyhist=False,loc=2,color='m')
plot_vert_hist(fig,ax1,smooth(emas_ref_full,60),1,[0,80],legend=True,color='k')
plot_vert_hist(fig,ax1,smooth(ssfr_ref,2),2,[0,80],legend=True,color='g')
plot_vert_hist(fig,ax1,smooth(rsp_ref,70),3,[0,80],legend=True,color='c')
plot_vert_hist(fig,ax1,smooth(star_ref,40),4,[0,80],legend=True,color='r')
plot_vert_hist(fig,ax1,probes[:,7],5,[0,80],legend=True,color='y')
(_,caps,_) = ax1.errorbar([0,1,2,3,4],
             [np.nanmean(smooth(modis_ref,6)),
              np.nanmean(smooth(emas_ref_full,60)),
              np.nanmean(smooth(ssfr_ref,2)),
              np.nanmean(smooth(rsp_ref,70)),
              np.nanmean(smooth(star_ref,40))],yerr=[1.7,0,4.7,1.7,1.8],
             color='r',linestyle='.',linewidth=3,label='Variability within FOV',capsize=7)
for cap in caps:
    cap.set_markeredgewidth(3)
ax1.legend(frameon=False,numpoints=1)
plt.savefig(fp+'plots/vert_hist_ref_v4.png',dpi=600,transparent=True)


# ## Scatter plots of the different effective radius

# In[97]:

ids = find_closest(er2['Start_UTC'],ssfr.utc)


# In[280]:

plt.plot(rsp['Lat'][rsp_good],smooth(rsp_ref,70),'c.',label='RSP')
plt.plot(er2['Latitude'][ids[ssfr.good[0][iss]]],smooth(ssfr_ref,2),'g.',label='SSFR')
plt.plot(emas_full['lat'][ief],smooth(emas_ref_full,60),'k.',label='eMAS')
plt.plot(meas.lat[meas.good[ist]],smooth(star_ref,40),'r.',label='4STAR')
plt.plot(modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]],smooth(modis_ref,6),'m.',label='MODIS')
plt.legend(frameon=False)
plt.xlabel('Latitude')
plt.ylabel('r$_{eff}$ [$\\mu$m]')


# In[693]:

plt.plot(dc8['G_LAT'][id[ccn_good]],ccn['Number_Concentration'][ccn_good],'b.')
plt.plot(dc8['G_LAT'][id[ccn_good]],smooth(ccn['Number_Concentration'][ccn_good],5,nan=False),'r.')
plt.xlabel('Latitude')
plt.ylabel('Number Concentration [\#/cm$^{3}$]')


# In[582]:

plt.plot(rsp['Lat'][rsp_good],smooth(rsp_ref,70)*smooth(rsp_tau,70)*2/3/10,'c.',label='RSP')
plt.plot(er2['Latitude'][ids[ssfr.good[0][iss]]],smooth(ssfr_ref,2)*smooth(ssfr_tau,2)*2/3/10,'g.',label='SSFR')
plt.plot(emas_full['lat'][ief],smooth(emas_ref_full,60)*smooth(emas_tau_full,60)*2/3/10,'k.',label='eMAS')
plt.plot(meas.lat[meas.good[ist]],smooth(star_ref,40)*smooth(star_tau,40)*2/3/10,'r.',label='4STAR')
plt.plot(modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]],smooth(modis_ref,6)*smooth(modis_tau,6)*2/3/10,'m.',label='MODIS')
plt.legend(frameon=False)
plt.xlabel('Latitude')
plt.ylabel('IWP [g/m$^{2}$]')


# In[694]:

id_w_up = np.where(dc8['W'][id[ccn_good]]>0.0)[0]


# In[695]:

plt.plot(dc8['G_LAT'][id[ccn_good]],dc8['W'][id[ccn_good]],'b.')
plt.plot(dc8['G_LAT'][id[ccn_good]],smooth(dc8['W'][id[ccn_good]],20),'r')
plt.plot(dc8['G_LAT'][id[ccn_good[id_w_up]]],dc8['W'][id[ccn_good[id_w_up]]],'g.')
plt.axhline(0,color='k')
plt.ylabel('Vertical wind speed [m/s]')
plt.xlabel('Latitude')


# ### Plot CCN vs. ref

# In[696]:

id_ccn_rsp = find_closest(dc8['G_LAT'][id[ccn_good]],rsp['Lat'][rsp_good])
id_ccn_ssfr = find_closest(dc8['G_LAT'][id[ccn_good]],er2['Latitude'][ids[ssfr.good[0][iss]]])
id_ccn_emas = find_closest(dc8['G_LAT'][id[ccn_good]],emas_full['lat'][ief])
id_ccn_star = find_closest(dc8['G_LAT'][id[ccn_good]],meas.lat[meas.good[ist]])
id_ccn_modis = find_closest(dc8['G_LAT'][id[ccn_good]],modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]])


# In[697]:

id_ccn_up_rsp = find_closest(dc8['G_LAT'][id[ccn_good[id_w_up]]],rsp['Lat'][rsp_good])
id_ccn_up_ssfr = find_closest(dc8['G_LAT'][id[ccn_good[id_w_up]]],er2['Latitude'][ids[ssfr.good[0][iss]]])
id_ccn_up_emas = find_closest(dc8['G_LAT'][id[ccn_good[id_w_up]]],emas_full['lat'][ief])
id_ccn_up_star = find_closest(dc8['G_LAT'][id[ccn_good[id_w_up]]],meas.lat[meas.good[ist]])
id_ccn_up_modis = find_closest(dc8['G_LAT'][id[ccn_good[id_w_up]]],modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]])


# In[698]:

iu_rsp = ccn_good[id_w_up[id_ccn_up_rsp]]
iu_ssfr = ccn_good[id_w_up[id_ccn_up_ssfr]]
iu_emas = ccn_good[id_w_up[id_ccn_up_emas]]
iu_star = ccn_good[id_w_up[id_ccn_up_star[:,0]]]
iu_modis = ccn_good[id_w_up[id_ccn_up_modis]]


# In[699]:

import plotting_utils as pu


# In[700]:

plt.plot(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_rsp]),np.log(smooth(rsp_ref,70)),'c.',label='RSP')
pu.plot_lin(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_rsp]),np.log(smooth(rsp_ref,70)),
            x_err=np.log(ccn['Number_Concentration'][ccn_good][id_ccn_rsp]*0.055),y_err=np.log(smooth(rsp_ref,70)*0.15),
            color='c',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_ssfr]),np.log(smooth(ssfr_ref,2)),'g+',label='SSFR')
pu.plot_lin(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_ssfr]),np.log(smooth(ssfr_ref,2)),
            x_err=np.log(ccn['Number_Concentration'][ccn_good][id_ccn_ssfr]*0.055),color='g',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_emas]),np.log(smooth(emas_ref_full,60)),'k*',label='eMAS')
pu.plot_lin(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_emas]),np.log(smooth(emas_ref_full,60)),
            x_err=np.log(ccn['Number_Concentration'][ccn_good][id_ccn_emas]*0.055),color='k',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_star]),np.log(smooth(star_ref,40)),'ro',label='4STAR')
pu.plot_lin(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_star])[:,0],np.log(smooth(star_ref,40)),
            x_err=np.log(ccn['Number_Concentration'][ccn_good][id_ccn_star]*0.055)[:,0],color='r',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_modis]),np.log(smooth(modis_ref,6)),'mx',label='MODIS')
pu.plot_lin(np.log(ccn['Number_Concentration'][ccn_good][id_ccn_modis]),np.log(smooth(modis_ref,6)),
            x_err=np.log(ccn['Number_Concentration'][ccn_good][id_ccn_modis]*0.055),color='m',use_method='statsmodels',ci=75)
plt.xlabel('Log(Number Concentration)')
plt.ylabel('Log(r$_{eff}$)')
plt.xlim([2.2,4.8])
plt.ylim([2.5,4.2])

plt.legend(bbox_to_anchor=[1,1],loc=2)
plt.savefig(fp+'plots/20130913_log_reff_vs_log_ccn.png',dpi=600,transparent=True)


# In[701]:

plt.plot(np.log(ccn['Number_Concentration'][iu_rsp]),np.log(smooth(rsp_ref,70)),'c.',label='RSP')
pu.plot_lin(np.log(ccn['Number_Concentration'][iu_rsp]),np.log(smooth(rsp_ref,70)),
            x_err=np.log(ccn['Number_Concentration'][iu_rsp]*0.055),y_err=np.log(smooth(rsp_ref,70)*0.15),
            color='c',use_method='statsmodels',ci=95)
plt.plot(np.log(ccn['Number_Concentration'][iu_ssfr]),np.log(smooth(ssfr_ref,2)),'g+',label='SSFR')
pu.plot_lin(np.log(ccn['Number_Concentration'][iu_ssfr]),np.log(smooth(ssfr_ref,2)),
            x_err=np.log(ccn['Number_Concentration'][iu_ssfr]*0.055),color='g',use_method='statsmodels',ci=95)
plt.plot(np.log(ccn['Number_Concentration'][iu_emas]),np.log(smooth(emas_ref_full,60)),'k*',label='eMAS')
pu.plot_lin(np.log(ccn['Number_Concentration'][iu_emas]),np.log(smooth(emas_ref_full,60)),
            x_err=np.log(ccn['Number_Concentration'][iu_emas]*0.055),color='k',use_method='statsmodels',ci=95)
plt.plot(np.log(ccn['Number_Concentration'][iu_star]),np.log(smooth(star_ref,40)),'ro',label='4STAR')
pu.plot_lin(np.log(ccn['Number_Concentration'][iu_star]),np.log(smooth(star_ref,40)),
            x_err=np.log(ccn['Number_Concentration'][iu_star]*0.055),color='r',use_method='statsmodels',ci=95)
plt.plot(np.log(ccn['Number_Concentration'][iu_modis]),np.log(smooth(modis_ref,6)),'mx',label='MODIS')
pu.plot_lin(np.log(ccn['Number_Concentration'][iu_modis]),np.log(smooth(modis_ref,6)),
            x_err=np.log(ccn['Number_Concentration'][iu_modis]*0.055),color='m',use_method='statsmodels',ci=95)
plt.xlabel('Log(Number Concentration)')
plt.ylabel('Log(r$_{eff}$)')
plt.title('CCN - r$_{eff}$ relationship - Upwelling only')
plt.xlim([2.2,4.8])
plt.ylim([2.5,4.2])

plt.legend(bbox_to_anchor=[1,1],loc=2)
plt.savefig(fp+'plots/20130913_log_reff_vs_log_ccn_up.png',dpi=600,transparent=True)


# In[291]:

plt.plot(smooth(star_ref,40),smooth(modis_ref,6),'mo')
#plt.plot(smooth(star_ref,40),smooth(emas_ref_v1,60),'ko')
#plt.plot(smooth(star_ref,40),smooth(ssfr_ref,2),'go')
#plt.plot(smooth(star_ref,40),smooth(rsp_ref,70),'co')


# Plotting the ccn vs. ref results in low numbers and repeats of values. This is because we are using a nearest-neighbor approximation. Now we repeat the process but with interpolation instead.

# In[702]:

from scipy import interpolate


# In[753]:

idc8_lat = np.argsort(dc8['G_LAT'][id[ccn_good]])
dc8_sort = np.sort(dc8['G_LAT'][id[ccn_good]])
ccn_sort = ccn['Number_Concentration'][ccn_good][idc8_lat]


# In[763]:

dc8_sort_nonunique,inonunique = np.unique(dc8_sort,return_index=True)
ccn_sort_nonunique = ccn_sort[inonunique]


# In[754]:

plt.plot(dc8_sort,ccn_sort)
plt.plot(dc8['G_LAT'][id[ccn_good]],ccn['Number_Concentration'][ccn_good],'.')


# In[788]:

f_ccn = interpolate.InterpolatedUnivariateSpline(dc8_sort_nonunique,ccn_sort_nonunique,k=1)
ccn_ref_rsp = f_ccn(rsp['Lat'][rsp_good])
ccn_ref_ssfr = f_ccn(er2['Latitude'][ids[ssfr.good[0][iss]]])
ccn_ref_emas = f_ccn(emas_full['lat'][ief])
ccn_ref_star = f_ccn(meas.lat[meas.good[ist],0])
ccn_ref_modis = f_ccn(modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]])


# In[803]:

i_extrap_rsp = np.where((rsp['Lat'][rsp_good] > np.max(dc8_sort_nonunique)) | (rsp['Lat'][rsp_good] < np.min(dc8_sort_nonunique)))
i_extrap_ssfr = np.where((er2['Latitude'][ids[ssfr.good[0][iss]]] > np.max(dc8_sort_nonunique)) | (er2['Latitude'][ids[ssfr.good[0][iss]]] < np.min(dc8_sort_nonunique)))
i_extrap_emas = np.where((emas_full['lat'][ief] > np.max(dc8_sort_nonunique)) | (emas_full['lat'][ief] < np.min(dc8_sort_nonunique)))
i_extrap_star = np.where((meas.lat[meas.good[ist],0] > np.max(dc8_sort_nonunique)) | (meas.lat[meas.good[ist],0] < np.min(dc8_sort_nonunique)))
i_extrap_modis = np.where((modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]] > np.max(dc8_sort_nonunique)) | (modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]] < np.min(dc8_sort_nonunique)))


# In[804]:

ccn_ref_rsp[i_extrap_rsp] = np.nan
ccn_ref_ssfr[i_extrap_ssfr] = np.nan
ccn_ref_emas[i_extrap_emas] = np.nan
ccn_ref_star[i_extrap_star] = np.nan
ccn_ref_modis[i_extrap_modis] = np.nan


# In[805]:

plt.plot(rsp['Lat'][rsp_good],ccn_ref_rsp,'c.',label='CCN interp to RSP')
plt.plot(er2['Latitude'][ids[ssfr.good[0][iss]]],ccn_ref_ssfr,'g+',label='CCN interp to SSFR')
plt.plot(emas_full['lat'][ief],ccn_ref_emas,'k*',label='CCN interp to eMAS')
plt.plot(meas.lat[meas.good[ist],0],ccn_ref_star,'ro',label='CCN interp to 4STAR')
plt.plot(modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]],ccn_ref_modis,'mx',label='CCN interp to MODIS')
plt.plot(rsp['Lat'][rsp_good],smooth(rsp_ref,70),'g.',label='RSP reff')
plt.plot(dc8['G_LAT'][id[ccn_good]],ccn['Number_Concentration'][ccn_good],'r.',label='original CCN')
plt.legend(frameon=False)
plt.xlabel('Latitude')
plt.ylabel('Number concentration [\#/cm$^{3}$] \underline{or} $r_{eff}$ [$\mu$m]')


# In[806]:

plt.plot(np.log(ccn_ref_rsp),np.log(smooth(rsp_ref,70)),'c.',label='RSP')
pu.plot_lin(np.log(ccn_ref_rsp),np.log(smooth(rsp_ref,70)),
            x_err=np.log(ccn_ref_rsp*0.055),y_err=np.log(smooth(rsp_ref,70)*0.15),
            color='c',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn_ref_ssfr),np.log(smooth(ssfr_ref,2)),'g+',label='SSFR')
pu.plot_lin(np.log(ccn_ref_ssfr),np.log(smooth(ssfr_ref,2)),
            x_err=np.log(ccn_ref_ssfr*0.055),color='g',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn_ref_emas),np.log(smooth(emas_ref_full,60)),'k*',label='eMAS')
pu.plot_lin(np.log(ccn_ref_emas),np.log(smooth(emas_ref_full,60)),
            x_err=np.log(ccn_ref_emas*0.055),color='k',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn_ref_star),np.log(smooth(star_ref,40)),'ro',label='4STAR')
pu.plot_lin(np.log(ccn_ref_star),np.log(smooth(star_ref,40)),
            x_err=np.log(ccn_ref_star*0.055),color='r',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn_ref_modis),np.log(smooth(modis_ref,6)),'mx',label='MODIS')
pu.plot_lin(np.log(ccn_ref_modis),np.log(smooth(modis_ref,6)),
            x_err=np.log(ccn_ref_modis*0.055),color='m',use_method='statsmodels',ci=75)
plt.xlabel('Log(Number Concentration)')
plt.ylabel('Log(r$_{eff}$)')
plt.xlim([2.2,4.8])
plt.ylim([2.5,4.2])

plt.legend(bbox_to_anchor=[1,1.05],loc=2)
plt.savefig(fp+'plots/20130913_log_reff_vs_log_ccn_interp.png',dpi=600,transparent=True)


# Now set up for only using the values that have positive vertical velocities as measured by dc8
# 
# First need to interpolate the vertical wind speed, then find the points

# In[833]:

idc8_good = np.where((dc8['TIME_UTC']>18.0)&(dc8['TIME_UTC']<19.5))[0]
idc8_sort_lat = np.argsort(dc8['G_LAT'][idc8_good])
dc8_lat_sort = dc8['G_LAT'][idc8_good[idc8_sort_lat]]
dc8_w_sort = dc8['W'][idc8_good[idc8_sort_lat]]
dc8_lat_nonunique,idc8_lat_nonunique = np.unique(dc8_lat_sort,return_index=True)
dc8_w_sort_nounique = dc8_w_sort[idc8_lat_nonunique]


# In[843]:

f_w_dc8 = interpolate.InterpolatedUnivariateSpline(dc8_lat_nonunique,smooth(dc8_w_sort_nounique,2))
dc8_w_rsp = f_w_dc8(rsp['Lat'][rsp_good])
dc8_w_ssfr = f_w_dc8(er2['Latitude'][ids[ssfr.good[0][iss]]])
dc8_w_emas = f_w_dc8(emas_full['lat'][ief])
dc8_w_star = f_w_dc8(meas.lat[meas.good[ist],0])
dc8_w_modis = f_w_dc8(modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]])


# In[853]:

dc8_w_rsp[dc8_w_rsp>20] = np.nan
dc8_w_ssfr[dc8_w_ssfr>20] = np.nan
dc8_w_emas[dc8_w_emas>20] = np.nan
dc8_w_star[dc8_w_star>20] = np.nan
dc8_w_modis[dc8_w_modis>20] = np.nan


# In[858]:

plt.plot(dc8['G_LAT'][id],dc8['W'][id],'r.')
plt.plot(rsp['Lat'][rsp_good],dc8_w_rsp,'c.')
plt.plot(er2['Latitude'][ids[ssfr.good[0][iss]]],dc8_w_ssfr,'g+')
plt.plot(emas_full['lat'][ief],dc8_w_emas,'k*')
plt.plot(meas.lat[meas.good[ist],0],dc8_w_star,'ro')
plt.plot(modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]],dc8_w_modis,'mx')
plt.plot([20.9,22.3],[0,0],'k')
plt.xlim([20.9,22.3])


# In[863]:

iup_rsp = dc8_w_rsp>0
iup_ssfr = dc8_w_ssfr>0
iup_emas = dc8_w_emas>0
iup_star = dc8_w_star>0
iup_modis = dc8_w_modis>0


# In[865]:

plt.plot(np.log(ccn_ref_rsp[iup_rsp]),np.log(smooth(rsp_ref,70)[iup_rsp]),'c.',label='RSP')
pu.plot_lin(np.log(ccn_ref_rsp[iup_rsp]),np.log(smooth(rsp_ref,70)[iup_rsp]),
            x_err=np.log(ccn_ref_rsp[iup_rsp]*0.055),y_err=np.log(smooth(rsp_ref,70)[iup_rsp]*0.15),
            color='c',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn_ref_ssfr[iup_ssfr]),np.log(smooth(ssfr_ref,2)[iup_ssfr]),'g+',label='SSFR')
pu.plot_lin(np.log(ccn_ref_ssfr[iup_ssfr]),np.log(smooth(ssfr_ref,2)[iup_ssfr]),
            x_err=np.log(ccn_ref_ssfr[iup_ssfr]*0.055),color='g',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn_ref_emas[iup_emas]),np.log(smooth(emas_ref_full,60)[iup_emas]),'k*',label='eMAS')
pu.plot_lin(np.log(ccn_ref_emas[iup_emas]),np.log(smooth(emas_ref_full,60)[iup_emas]),
            x_err=np.log(ccn_ref_emas[iup_emas]*0.055),color='k',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn_ref_star[iup_star]),np.log(smooth(star_ref,40)[iup_star]),'ro',label='4STAR')
pu.plot_lin(np.log(ccn_ref_star[iup_star]),np.log(smooth(star_ref,40)[iup_star]),
            x_err=np.log(ccn_ref_star[iup_star]*0.055),color='r',use_method='statsmodels',ci=75)
plt.plot(np.log(ccn_ref_modis[iup_modis]),np.log(smooth(modis_ref,6)[iup_modis]),'mx',label='MODIS')
pu.plot_lin(np.log(ccn_ref_modis[iup_modis]),np.log(smooth(modis_ref,6)[iup_modis]),
            x_err=np.log(ccn_ref_modis[iup_modis]*0.055),color='m',use_method='statsmodels',ci=75)
plt.xlabel('Log(Number Concentration)')
plt.ylabel('Log(r$_{eff}$)')
plt.xlim([2.2,4.8])
plt.ylim([2.5,4.2])

plt.legend(bbox_to_anchor=[1,1.05],loc=2)
plt.savefig(fp+'plots/20130913_log_reff_vs_log_ccn_interp_up.png',dpi=600,transparent=True)


# Recalculate the reff vs. ccn with derivative

# In[866]:

from Sp_parameters import deriv


# In[882]:

ccn_rsp_dif = smooth(deriv(np.log(smooth(rsp_ref,70)),np.log(ccn_ref_rsp)),5)
ccn_ssfr_dif = smooth(deriv(np.log(smooth(ssfr_ref,2)),np.log(ccn_ref_ssfr)),5)
ccn_emas_dif = smooth(deriv(np.log(smooth(emas_ref_full,60)),np.log(ccn_ref_emas)),5)
ccn_star_dif = smooth(deriv(np.log(smooth(star_ref,40)),np.log(ccn_ref_star)),5)
ccn_modis_dif = smooth(deriv(np.log(smooth(modis_ref,6)),np.log(ccn_ref_modis)),5)


# In[876]:

np.nanmean(ccn_rsp_dif)


# In[887]:

plt.plot(rsp['Lat'][rsp_good],ccn_rsp_dif,'c.')
plt.plot(er2['Latitude'][ids[ssfr.good[0][iss]]],ccn_ssfr_dif,'g+')
plt.plot(emas_full['lat'][ief],ccn_emas_dif,'k*')
plt.plot(meas.lat[meas.good[ist],0],ccn_star_dif,'ro')
plt.plot(modis['lat'][dc8_ind_modis[0,im],dc8_ind_modis[1,im]],ccn_modis_dif,'mx')
plt.plot([20.8,22.4],[0,0],'k')
plt.ylim([-5,5])


# In[895]:

plt.hist(ccn_rsp_dif,bins=50, histtype='stepfilled', normed=True, color='c',alpha=0.6,range=[-20,20])
plt.hist(ccn_ssfr_dif,bins=50, histtype='stepfilled', normed=True, color='g',alpha=0.6,range=(-20,20))
plt.hist(ccn_emas_dif,bins=50, histtype='stepfilled', normed=True, color='k',alpha=0.6,range=(-20,20))
plt.hist(ccn_star_dif,bins=50, histtype='stepfilled', normed=True, color='r',alpha=0.6,range=(-20,20))
plt.hist(ccn_modis_dif,bins=50, histtype='stepfilled', normed=True, color='m',alpha=0.6,range=(-20,20))


# In[898]:

reload(pu)


# In[901]:

fig = plt.figure(figsize=(7,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-2,2],xlim=[-0.5,4.5])
ax1.set_ylabel('ACI$_{r}$')
ax1.set_xticks([0,1,2,3,4])
ax1.set_xticklabels(['MODIS\n(reflected)','eMAS\n(reflected)','SSFR\n(reflected)','RSP\n(reflected)','4STAR\n(transmitted)'])
pu.plot_vert_hist(fig,ax1,ccn_modis_dif,0,[-2,2],legend=True,onlyhist=False,loc=2,color='m')
pu.plot_vert_hist(fig,ax1,ccn_emas_dif,1,[-2,2],legend=True,color='k')
pu.plot_vert_hist(fig,ax1,ccn_ssfr_dif,2,[-2,2],legend=True,color='g')
pu.plot_vert_hist(fig,ax1,ccn_rsp_dif,3,[-2,2],legend=True,color='c')
pu.plot_vert_hist(fig,ax1,ccn_star_dif,4,[-2,2],legend=True,color='r')
plt.savefig(fp+'plots/vert_hist_ACI_r.png',dpi=600,transparent=True)


# ## Combine the vertical information into one figure

# In[133]:

meas.good.shape


# In[134]:

it = np.abs(apr['utc']-19.19).argmin()
inoisezone = apr['dbz'][:,it].argmax()
i=range(inoisezone-15)+range(inoisezone+15,len(apr['altflt'][:,0]))


# In[160]:

plt.figure
plt.plot(smooth(twoDS['effectiveD'][fl]/2.0,10),dc8['G_ALT'][fl],label='Cloud Probes (2DS)')
plt.plot(emas_ref_v1,cpl_layers['top'][dc8_ind[0,ie1],0],'k+',label='eMAS')
plt.plot(rsp_ref,cpl_layers['top'][rsp_good[0][irs],0],'c+',label='RSP')
plt.plot(ssfr_ref,cpl_layers['top'][ssfr.good[0][iss],0],'g+',label='SSFR')
plt.plot(modis_ref,cpl_layers['top'][dc8_ind_modis[0,im],0],'m+',label='MODIS')
plt.plot(star_ref,dc8['G_ALT'][meas.good[ist]],'r+',label='4STAR')
#plt.plot(smooth(apr['dbz'][i,it],10)/50.0,apr['altflt'][i,it],'k',label='APR-2 Reflectivity')
plt.legend(frameon=False)
plt.xlim([0,100])
plt.ylim([4000,18000])
plt.xlabel('Effective radius [$\\mu$m]')
plt.ylabel('Altitude [m]')
plt.title('Vertical profile of Effective radius')
plt.savefig(fp+'plots/ref_profile_seac4rs.png',dpi=600,transparent=True)


# In[159]:

plt.figure
plt.plot(smooth(twoDS['effectiveD'][fl]/2.0,10),dc8['G_ALT'][fl],label='Cloud Probes (2DS)')
plt.plot(emas_ref_v1,cpl_layers['top'][dc8_ind[0,ie1],0],'k+',label='eMAS')
plt.plot(rsp_ref,cpl_layers['top'][rsp_good[0][irs],0],'c+',label='RSP')
plt.plot(ssfr_ref,cpl_layers['top'][ssfr.good[0][iss],0],'g+',label='SSFR')
plt.plot(modis_ref,cpl_layers['top'][dc8_ind_modis[0,im],0],'m+',label='MODIS')
plt.plot(star_ref,dc8['G_ALT'][meas.good[ist]],'r+',label='4STAR')
plt.plot(smooth(apr['dbz'][i,it],10)/50.0,apr['altflt'][i,it],'k',label='APR-2 Reflectivity')
plt.legend(frameon=False)
plt.xlim([0,100])
plt.ylim([4000,18000])
plt.xlabel('Effective radius [$\\mu$m]')
plt.ylabel('Altitude [m]')
plt.title('Vertical profile of Effective radius')
plt.savefig(fp+'plots/ref_profile_seac4rs_radar.png',dpi=600,transparent=True)


# ## eMAS V00 and V01 comparison

# In[73]:

emas_v1.keys()


# In[148]:

plt.figure();
plt.plot(emas_v1['lon'][dc8_ind[0,:],dc8_ind[1,:]],emas_tau_v1,label='V01')
plt.plot(emas['lon'][dc8_ind[0,:],dc8_ind[1,:]], emas_tau,label='V00')



# # 1:1 relationship

# In[75]:

plt.figure()
plt.plot(emas_tau_v1,emas_tau,'+',label=r'eMAS $\tau$')
plt.plot([10,35],[10,35],'k--',label='one-to-one')
from linfit import linfit
emas_fit,z = linfit(emas_tau_v1,emas_tau)
plt.plot(np.linspace(10,39), emas_fit[1]+emas_fit[0]*np.linspace(10,39),'r--',label='Linear fit:\n y='+str('%.3f' % emas_fit[0])+'x+'+str('%.3f' % emas_fit[1]))
plt.title(r'eMAS version comparison along DC8 flight track on 2013-09-13')
plt.xlabel(r'eMAS V01 $\tau$')
plt.ylabel(r'eMAS V00 $\tau$')
plt.legend(frameon=False,loc=4)
plt.savefig(fp+'plots/emas_v00_compare_v01_tau.png',dpi=600,transparent=True)


# In[76]:

print emas_tau_v1, emas_tau
print linfit(emas_tau_v1,emas_tau)


# In[77]:

plt.figure()
plt.hist(emas_tau_v1-emas_tau,bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='eMAS tau difference',range=(-5,5))
plt.xlabel(r'$\tau$ difference')
plt.ylabel('Normed probability')
plt.title(r'eMAS $\tau$ difference (V01 - V00)')
plt.savefig(fp+'plots/emas_diff_histogram_v01_v00_tau.png',dpi=600,transparent=True)
np.max(emas_tau_v1-emas_tau)


# In[78]:

plt.figure()
plt.plot(emas_ref_v1,emas_ref,'+',label='eMAS $r_{ef}$')
plt.plot([22,36],[22,36],'k--', label='one-to-one')
plt.title('eMAS version comparison along DC8 flight track on 2013-09-13')
plt.xlabel('eMAS V01 $r_{ef}$ [$\mu$m]')
plt.ylabel('eMAS V00 $r_{ef}$ [$\mu$m]')
plt.legend(frameon=False,loc=4)
plt.savefig(fp+'plots/emas_v00_compare_v01_ref.png',dpi=600,transparent=True)


# In[79]:

plt.figure()
plt.hist(emas_ref_v1-emas_ref,bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='eMAS Ref difference',range=(-5,5))
plt.xlabel('r$_{ef}$ difference [$\mu$m]')
np.max(emas_ref_v1-emas_ref)


# ## Find the mean tau and ref for each instrument

# In[75]:

print np.nanmean(modis_ref)


# In[76]:

print np.shape(modis_tau)


# In[77]:

means = object()
means.tau = object()
means.ref = object()
means.tau.modis = np.nanmean(modis_tau)
means.ref.modis = np.nanmean(modis_ref)
means.tau.emas = np.nanmean(emas_tau)
means.ref.emas = np.nanmean(emas_ref)
means.tau.ssfr = np.nanmean(ssfr_tau)
means.ref.ssfr = np.nanmean(ssfr_ref)
means.tau.star = np.nanmean(star_tau)
means.ref.star = np.nanmean(star_ref)
means.ref.probes = np.nanmean(probes[:,7])
means.tau.rsp = np.nanmean(rsp_tau)
means.ref.rsp = np.nanmean(rsp_ref)


# In[78]:

sttaumod = np.nanstd(modis_tau)
strefmod = np.nanstd(modis_ref)
sttauemas = np.nanstd(emas_tau)
strefemas = np.nanstd(emas_ref)
sttaussfr = np.nanstd(ssfr_tau)
strefssfr = np.nanstd(ssfr_ref)
sttaursp = np.nanstd(rsp_tau)
strefrsp = np.nanstd(rsp_ref)
sttaustar = np.nanstd(star_tau)
strefstar = np.nanstd(star_ref)


# In[83]:

print lut


# In[84]:

print lut.sp_irrdn.shape
print lut.tausp.shape
print lut.refsp.shape
print lut.wvl.shape


# interpolate the modeled irradiance at z=0 (dc8 altitude) to the retrieved values of tau and ref

# In[192]:

lut.sp_irrdn[1,:,0,4,0]


# In[79]:

from scipy import interpolate
#fx = interpolate.RectBivariateSpline(ref[refranges[ph]],tau,sp[ph,w,z,refranges[ph],:],kx=1,ky=1)
#sp_hires[ph,w,z,:,:] = fx(ref_hires,tau_hires)


# In[81]:

def interpirr(re,ta):
    "interpolate over irradiance to get the spectrum at ta and re"
    sp = np.zeros(lut.wvl.size)
    print re,ta
    for w in xrange(lut.wvl.size):
        irrdnfx = interpolate.RectBivariateSpline(lut.refsp[lut.refsp>5],lut.tausp,lut.sp_irrdn[1,w,0,lut.refsp>5,:],kx=1,ky=1)
        sp[w] = irrdnfx(re,ta)
    return smooth(sp,20)


# In[82]:

def interpirr_up(re,ta):
    "interpolate over irradiance to get the spectrum at ta and re"
    sp = np.zeros(lut.wvl.size)
    print re,ta
    for w in xrange(lut.wvl.size):
        irrupfx = interpolate.RectBivariateSpline(lut.refsp[lut.refsp>5],lut.tausp,lut.sp_irrup[1,w,0,lut.refsp>5,:],kx=1,ky=1)
        sp[w] = irrupfx(re,ta)
    return smooth(sp,20)


# In[83]:

sp_modis = interpirr(means.ref.modis,means.tau.modis)
sp_emas = interpirr(means.ref.emas,means.tau.emas)
sp_ssfr = interpirr(means.ref.ssfr,means.tau.ssfr)
sp_star = interpirr(means.ref.star,means.tau.star)
sp_rsp = interpirr(means.ref.rsp,means.tau.rsp)
sp_modis_no = interpirr(means.ref.modis,means.tau.modis*0.0)
sp_emas_no = interpirr(means.ref.emas,means.tau.emas*0.0)
sp_ssfr_no = interpirr(means.ref.ssfr,means.tau.ssfr*0.0)
sp_star_no = interpirr(means.ref.star,means.tau.star*0.0)
sp_rsp_no = interpirr(means.ref.rsp,means.tau.rsp*0.0)


# In[84]:

sp_modisp = interpirr(means.ref.modis+sttaumod,means.tau.modis+sttaumod)
sp_emasp = interpirr(means.ref.emas+sttauemas,means.tau.emas+sttauemas)
sp_ssfrp = interpirr(means.ref.ssfr+sttaussfr,means.tau.ssfr+sttaussfr)
sp_starp = interpirr(means.ref.star+sttaustar,means.tau.star+sttaustar)
sp_rspp = interpirr(means.ref.rsp+sttaursp,means.tau.rsp+sttaursp)


# In[85]:

sp_modisp_up = interpirr_up(means.ref.modis+sttaumod,means.tau.modis+sttaumod)
sp_emasp_up = interpirr_up(means.ref.emas+sttauemas,means.tau.emas+sttauemas)
sp_ssfrp_up = interpirr_up(means.ref.ssfr+sttaussfr,means.tau.ssfr+sttaussfr)
sp_starp_up = interpirr_up(means.ref.star+sttaustar,means.tau.star+sttaustar)
sp_rspp_up = interpirr_up(means.ref.rsp+sttaursp,means.tau.rsp+sttaursp)


# In[86]:

sp_modis_up = interpirr_up(means.ref.modis,means.tau.modis)
sp_emas_up = interpirr_up(means.ref.emas,means.tau.emas)
sp_ssfr_up = interpirr_up(means.ref.ssfr,means.tau.ssfr)
sp_star_up = interpirr_up(means.ref.star,means.tau.star)
sp_rsp_up = interpirr_up(means.ref.rsp,means.tau.rsp)
sp_modis_no_up = interpirr_up(means.ref.modis,means.tau.modis*0.0)
sp_emas_no_up = interpirr_up(means.ref.emas,means.tau.emas*0.0)
sp_ssfr_no_up = interpirr_up(means.ref.ssfr,means.tau.ssfr*0.0)
sp_star_no_up = interpirr_up(means.ref.star,means.tau.star*0.0)
sp_rsp_no_up = interpirr_up(means.ref.rsp,means.tau.rsp*0.0)


# In[87]:

print list(enumerate(lut.refsp))
print list(enumerate(lut.tausp))
print lut.sp_irrdn.shape


# In[89]:

plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,13,13]*10.0,'m',label='modis')
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,10,12]*10.0,'b',label='emas')
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,29,17]*10.0,'g',label='SSFR')
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,24,11]*10.0,'r',label='RSP')
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,9,13]*10.0,'c',label='4STAR')
#plt.plot(ssfr_dc8[''])
plt.legend(frameon=False)


# In[90]:

print sp_rsp


# Load the SSFR from the DC8

# In[137]:

ssfr_dc8 = sio.idl.readsav(fp+'dc8/20130913/20130913_calibspcs.out')
print ssfr_dc8.keys()
iutc185 = np.nanargmin(abs(ssfr_dc8['tmhrs']-18.5))
iutc189 = np.nanargmin(abs(ssfr_dc8['tmhrs']-18.9))
iutc191 = np.nanargmin(abs(ssfr_dc8['tmhrs']-19.1))
iutc192 = np.nanargmin(abs(ssfr_dc8['tmhrs']-19.2))
print ssfr_dc8['zspectra'].shape
dn = np.nanmean(ssfr_dc8['zspectra'][iutc189:iutc191,:],axis=0)
print dn.shape
up = np.nanmean(ssfr_dc8['nspectra'][iutc189:iutc191,:],axis=0)


# In[138]:

dnstd = np.nanstd(ssfr_dc8['zspectra'][iutc191:iutc192,:],axis=0)


# In[139]:

plt.plot(dnstd)


# In[140]:

dn_no = interpirr(5.0,0)


# In[141]:

plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,4,0])
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,4,1])
plt.plot(lut.wvl,dn_no)
plt.plot(ssfr_dc8['zenlambda'],dn)
plt.plot(ssfr_dc8['zenlambda'],up)
plt.plot((dn-up)-(dn_no-up))


# In[142]:

ssfr_dc8_ict = load_ict(fp+'dc8/20130913/SEAC4RS-SSFR_DC8_20130913_R0.ict')


# In[92]:

iutc185tt = np.nanargmin(abs(ssfr_dc8_ict['UTC']-18.5))
iutc192tt = np.nanargmin(abs(ssfr_dc8_ict['UTC']-19.2))
iutc190tt = np.nanargmin(abs(ssfr_dc8_ict['UTC']-19.0))
ssfr_dc8_ict_good = np.where((ssfr_dc8_ict['UTC']>19.0) & (ssfr_dc8_ict['UTC']<19.2) & (ssfr_dc8_ict['DN500']>0))


# In[93]:

s500 = np.nanmean(ssfr_dc8_ict['DN500'][ssfr_dc8_ict_good[0]])
print s500
plt.plot(ssfr_dc8_ict['UTC'][ssfr_dc8_ict_good[0]],ssfr_dc8_ict['DN500'][ssfr_dc8_ict_good[0]])


# In[143]:

idc500 = np.nanargmin(abs(ssfr_dc8['zenlambda']-500.0))
r500 = s500/dn[idc500]
print r500
fssp = interpolate.interp1d(ssfr_dc8['zenlambda'],dn*r500,bounds_error=False)
ssp = fssp(lut.wvl)


# In[156]:

fssps = interpolate.interp1d(ssfr_dc8['zenlambda'],(dn+dnstd)*r500,bounds_error=False)
ssps = fssps(lut.wvl)


# In[144]:

fspup = interpolate.interp1d(ssfr_dc8['zenlambda'],up*r500,bounds_error=False)
sspup = fspup(lut.wvl)


# In[145]:

plt.plot(lut.wvl,dn_no*sspup/ssp)


# In[157]:

plt.plot(lut.wvl[250:],ssp[250:])
plt.plot(lut.wvl[250:],ssps[250:])


# In[158]:

plt.plot(lut.wvl,sspup,'k')
plt.plot(lut.wvl,ssp,'r')
plt.plot(lut.wvl,dn_no,'b')
plt.plot(lut.wvl,(ssp-sspup)-(dn_no-sspup))
print 'SSFR forcing, measurement based: %f W/m^2' % np.nansum(((ssp[250:]-sspup[250:])-(dn_no[250:]-sspup[250:]))*lut.wvl[250:]/1000.0)
print 'SSFR forcing, std measurement based: %f W/m^2' % np.nansum(((ssps[250:]-sspup[250:])-(dn_no[250:]-sspup[250:]))*lut.wvl[250:]/1000.0)
                                                               


# In[100]:

plt.figure
plt.plot(lut.wvl,(sp_modis-sp_modis_up)-(sp_modis_no-sp_modis_no_up),c='m',label='Modis')
plt.plot(lut.wvl,(sp_emas-sp_emas_up)-(sp_emas_no-sp_emas_no_up),c='b',label='eMAS')
plt.plot(lut.wvl,(sp_ssfr-sp_ssfr_up)-(sp_ssfr_no-sp_ssfr_no_up),c='g',label='SSFR (reflected)')
plt.plot(lut.wvl,(sp_rsp-sp_rsp_up)-(sp_rsp_no-sp_rsp_no_up),c='c',label='RSP')
plt.plot(lut.wvl,(sp_star-sp_star_up)-(sp_star_no-sp_star_no_up),c='r',label='4STAR')
plt.xlim([400,1700])
plt.legend(frameon=False,prop={'size':10},loc=4)
plt.title('Spectral radiative effect')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Change in net irradiance [$Wm^{-2}nm^{-1}$]')
plt.savefig(fp+'plots/sp_rad_effect.png',dpi=600,transparent=True)


# In[109]:

radeff_modis = np.sum(((sp_modis[250:]-sp_modis_up[250:])-(sp_modis_no[250:]-sp_modis_no_up[250:]))*lut.wvl[250:]/1000.0)
radeff_emas = np.sum(((sp_emas[250:]-sp_emas_up[250:])-(sp_emas_no[250:]-sp_emas_no_up[250:]))*lut.wvl[250:]/1000.0)
radeff_ssfr = np.sum(((sp_ssfr[250:]-sp_ssfr_up[250:])-(sp_ssfr_no[250:]-sp_ssfr_no_up[250:]))*lut.wvl[250:]/1000.0)
radeff_rsp = np.sum(((sp_rsp[250:]-sp_rsp_up[250:])-(sp_rsp_no[250:]-sp_rsp_no_up[250:]))*lut.wvl[250:]/1000.0)
radeff_star = np.sum(((sp_star[250:]-sp_star_up[250:])-(sp_star_no[250:]-sp_star_no_up[250:]))*lut.wvl[250:]/1000.0)


# In[110]:

radeff_modisp = np.sum(((sp_modisp-sp_modisp_up)-(sp_modis_no-sp_modis_no_up))*lut.wvl/1000.0)
radeff_emasp = np.sum(((sp_emasp-sp_emasp_up)-(sp_emas_no-sp_emas_no_up))*lut.wvl/1000.0)
radeff_ssfrp = np.sum(((sp_ssfrp-sp_ssfrp_up)-(sp_ssfr_no-sp_ssfr_no_up))*lut.wvl/1000.0)
radeff_rspp = np.sum(((sp_rspp-sp_rspp_up)-(sp_rsp_no-sp_rsp_no_up))*lut.wvl/1000.0)
radeff_starp = np.sum(((sp_starp-sp_starp_up)-(sp_star_no-sp_star_no_up))*lut.wvl/1000.0)


# In[111]:

print 'Modis: ',radeff_modis, ' +/- ', radeff_modisp-radeff_modis
print 'eMAS: ',radeff_emas, ' +/- ', radeff_emasp-radeff_emas
print 'SSFR: ',radeff_ssfr, ' +/- ', radeff_ssfrp-radeff_ssfr
print 'RSP: ',radeff_rsp, ' +/- ', radeff_rspp-radeff_rsp
print '4star: ', radeff_star, ' +/- ', radeff_starp-radeff_star


# In[169]:

plt.figure(figsize=(13,10))
f0,a0 = plt.subplots(2,1,sharex=True)
a0[0].plot(lut.wvl,ssp,label='SSFR Measurement',c='k')
a0[0].axvspan(350,1700,color='#FFFFFF')
a0[1].axvspan(350,1700,color='#FFFFFF')
a0[0].plot(lut.wvl,sp_modis,c='m',label='Modis')
a0[0].plot(lut.wvl,sp_emas,c='b',label='eMAS')
a0[0].plot(lut.wvl,sp_ssfr,c='g',label='SSFR (reflected)')
a0[0].plot(lut.wvl,sp_rsp/10.0,c='c',label='RSP')
a0[0].plot(lut.wvl,sp_star,c='r',label='4STAR')
a0[0].legend(frameon=False,prop={'size':10})
a0[0].set_title('Below cloud downwelling irradiance')
a0[1].set_xlabel('Wavelength [nm]')
a0[0].set_ylabel('Irradiance\n[Wm$^{-2}$nm$^{-1}$]')
a0[0].set_xlim([400,1700])
a0[1].plot(lut.wvl,(ssp-sp_modis)/ssp*100.0,c='m',label='Modis')
a0[1].plot(lut.wvl,(ssp-sp_emas)/ssp*100.0,c='b',label='eMAS')
a0[1].plot(lut.wvl,(ssp-sp_ssfr)/ssp*100.0,c='g',label='SSFR (reflected)')
a0[1].plot(lut.wvl,(ssp-sp_rsp/10.0)/ssp*100.0,c='c',label='RSP')
a0[1].plot(lut.wvl,(ssp-sp_modis)/ssp*100.0,c='r',label='4STAR')
a0[1].axhline(0,linestyle='--',c='k')
a0[1].set_ylabel('Irradiance difference\n[$\%$]')
a0[1].set_ylim([-40,100])
plt.savefig(fp+'plots/sp_compare_mod_meas.png',dpi=600,transparent=True)


# ## Get GOES - Problematic now...

# In[120]:

goes_file = fp+'er2/20130913/G13V04.0.CONUS.2013256.2045.PX.04K.NC'
goes,goes_dicts = load_hdf(goes_file)


# In[133]:

goes_values = (('lat',6),('lon',7),('tau',19),('ref',20))
goes,goes_dicts = load_hdf(goes_file,values=goes_values)


# In[125]:

figmg = plt.figure()
axmg = figm.add_axes([0.1,0.1,0.8,0.8])
mg = Basemap(projection='stere',lon_0=-87,lat_0=20,
        llcrnrlon=-100, llcrnrlat=15,
        urcrnrlon=-85, urcrnrlat=25,resolution='h',ax=axmg)
mg.drawcoastlines()
#m.fillcontinents(color='#AAAAAA')
mg.drawstates()
mg.drawcountries()
figmg.show()
x,y = mg(goes['lon'],goes['lat'])
csg = mg.contourf(x,y,goes['tau'],clevels,cmap=plt.cm.gist_ncar)


# ## Plot out FOV of each instrument on eMAS figure

# In[113]:

from plotting_utils import circles
from map_utils import radius_m2deg, spherical_dist


# In[17]:

plt.figure()
clevels = range(0,50,2)
cf = plt.contourf(emas_v1['lon'],emas_v1['lat'],emas_v1['tau'],clevels,cmap=plt.cm.spectral,extend='max')
ier2 = np.where((er2['Start_UTC']>19.08) & (er2['Start_UTC']<19.3))[0]
idc8 = np.where((dc8['TIME_UTC']>19.08) & (dc8['TIME_UTC']<19.3))[0]
plt.plot(er2['Longitude'][ier2],er2['Latitude'][ier2],'r',label='ER2 flight path')
plt.plot(dc8['G_LONG'][idc8],dc8['G_LAT'][idc8],'b',label='DC8 flight path')
#plt.plot(dc8['G_LONG'][idc8[400]],dc8['G_LAT'][idc8[400]],'k*')

circles(dc8['G_LONG'][idc8[400]],dc8['G_LAT'][idc8[400]],radius_m2deg(dc8['G_LONG'][idc8[400]],dc8['G_LAT'][idc8[400]],r_4STAR),c='r',alpha=0.4)
circles(er2['Longitude'][ier2[400]],er2['Latitude'][ier2[400]],radius_m2deg(er2['Longitude'][ier2[400]],er2['Latitude'][ier2[400]],r_SSFR),c='b',alpha=0.4)
circles(er2['Longitude'][ier2[400]],er2['Latitude'][ier2[400]],radius_m2deg(er2['Longitude'][ier2[400]],er2['Latitude'][ier2[400]],r_eMAS),c='g',alpha=0.4)
circles(er2['Longitude'][ier2[400]],er2['Latitude'][ier2[400]],radius_m2deg(er2['Longitude'][ier2[400]],er2['Latitude'][ier2[400]],r_RSP),c='y',alpha=0.4)

plt.title('FOV comparison')
plt.ylabel('Longitude')
plt.xlabel('Latitude')
plt.legend(loc=2,frameon=False)

cbar = plt.colorbar(cf)
cbar.set_label('$\\tau$')
plt.savefig(fp+'plots/emas_FOV.png',dpi=600,transparent=True)


# In[114]:

ier2 = np.where((er2['Start_UTC']>19.08) & (er2['Start_UTC']<19.3))[0]
idc8 = np.where((dc8['TIME_UTC']>19.08) & (dc8['TIME_UTC']<19.3))[0]
clevels = range(0,50,2)

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(122)
cf = ax.contourf(emas_v1['lon'],emas_v1['lat'],emas_v1['tau'],clevels,cmap=plt.cm.rainbow,extend='max')

ax.plot(er2['Longitude'][ier2],er2['Latitude'][ier2],'r',label='ER2 flight path')
ax.plot(dc8['G_LONG'][idc8],dc8['G_LAT'][idc8],'b--',label='DC8 flight path')

circles(dc8['G_LONG'][idc8[380]],dc8['G_LAT'][idc8[380]],radius_m2deg(dc8['G_LONG'][idc8[380]],dc8['G_LAT'][idc8[380]],r_4STAR),c='r',alpha=0.5,ax=ax,label='4STAR')
circles(er2['Longitude'][ier2[390]],er2['Latitude'][ier2[390]],radius_m2deg(er2['Longitude'][ier2[390]],er2['Latitude'][ier2[390]],r_SSFR),c='g',alpha=0.5,ax=ax,label='SSFR')
circles(er2['Longitude'][ier2[385]],er2['Latitude'][ier2[385]],radius_m2deg(er2['Longitude'][ier2[385]],er2['Latitude'][ier2[385]],r_eMAS),c='k',alpha=0.5,ax=ax,label='eMAS')
circles(er2['Longitude'][ier2[395]],er2['Latitude'][ier2[395]],radius_m2deg(er2['Longitude'][ier2[395]],er2['Latitude'][ier2[395]],r_RSP),c='c',alpha=0.5,ax=ax,label='RSP')
#circles(er2['Longitude'][ier2[375]],er2['Latitude'][ier2[375]],radius_m2deg(er2['Longitude'][ier2[375]],er2['Latitude'][ier2[375]],500.0),c='m',alpha=0.5,ax=ax,label='MODIS')
vo = [[modis['lon'][495,1208],modis['lat'][495,1208]],
      [modis['lon'][495,1209],modis['lat'][495,1209]],
      [modis['lon'][496,1209],modis['lat'][496,1209]],
      [modis['lon'][496,1208],modis['lat'][496,1208]],
      [modis['lon'][495,1208],modis['lat'][495,1208]]]
ax.add_patch(Polygon(vo,closed=True,color='m',alpha=0.5))


ax.text(dc8['G_LONG'][idc8[380]]+0.002,dc8['G_LAT'][idc8[380]],'4STAR',color='r')
ax.text(-93.958,21.53,'SSFR',color='g')
ax.text(er2['Longitude'][ier2[385]]+0.002,er2['Latitude'][ier2[385]],'eMAS',color='k')
ax.text(er2['Longitude'][ier2[395]]+0.002,er2['Latitude'][ier2[395]],'RSP',color='c')
ax.text(er2['Longitude'][ier2[375]]+0.005,er2['Latitude'][ier2[375]],'MODIS',color='m')

#dc8toer2 = spherical_dist(np.array((er2['Latitude'][ier2[393]],er2['Longitude'][ier2[393]])),np.array((dc8['G_LAT'][idc8[380]],dc8['G_LONG'][idc8[380]])))
#ax.annotate('' , (er2['Longitude'][ier2[393]],er2['Latitude'][ier2[393]]),(dc8['G_LONG'][idc8[380]],dc8['G_LAT'][idc8[380]]), arrowprops={'arrowstyle':'<->'})
#ax.text(-93.991,21.515,'%.2f km' % dc8toer2)

ya,yb,xa,xb = 21.46,21.54,-94.03,-93.94
ax.set_ylim([ya,yb])
ax.set_xlim([xa,xb])
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.set_title('Field-Of-View comparison')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

cbar = plt.colorbar(cf,fraction=0.046, pad=0.04)
cbar.set_label('$\\tau$')

ax2 = fig.add_subplot(121)
#ax2.locator_params(nbins=4)
cf = ax2.contourf(emas_v1['lon'],emas_v1['lat'],emas_v1['tau'],clevels,cmap=plt.cm.rainbow,extend='max')
ax2.plot(er2['Longitude'][ier2],er2['Latitude'][ier2],'r',label='ER2 flight path')
ax2.plot(dc8['G_LONG'][idc8],dc8['G_LAT'][idc8],'b--',label='DC8 flight path')
ax2.plot([xa,xa,xb,xb,xa],[ya,yb,yb,ya,ya],'k')
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')
ax2.set_title('eMAS Cloud optical thickness')
ax2.legend(loc=4,frameon=False)
yr = ax2.get_ylim()
xr = ax2.get_xlim()
ax2.plot([xa,-93.02],[yb,22.107],'k',clip_on=False)
ax2.plot([xa,-93.02],[ya,20.912],'k',clip_on=False)
ax2.set_ylim(yr)
ax2.set_xlim(xr)
plt.tight_layout()
plt.savefig(fp+'plots/emas_FOV_comp.png',dpi=600,transparent=True)


# In[63]:

import map_utils
reload(map_utils)
from map_utils import stats_within_radius


# In[64]:

out_ssfrtau = stats_within_radius(er2['Latitude'][ier2[::10]],er2['Longitude'][ier2[::10]],emas_v1['lat'],emas_v1['lon'],emas_v1['tau'],r_SSFR)
out_ssfrref = stats_within_radius(er2['Latitude'][ier2[::10]],er2['Longitude'][ier2[::10]],emas_v1['lat'],emas_v1['lon'],emas_v1['ref'],r_SSFR)


# In[70]:

out_rspref = stats_within_radius(er2['Latitude'][ier2[::10]],er2['Longitude'][ier2[::10]],emas_v1['lat'],emas_v1['lon'],emas_v1['ref'],r_RSP)
out_rsptau = stats_within_radius(er2['Latitude'][ier2[::10]],er2['Longitude'][ier2[::10]],emas_v1['lat'],emas_v1['lon'],emas_v1['tau'],r_RSP)


# In[72]:

out_starref = stats_within_radius(dc8['G_LAT'][idc8],dc8['G_LONG'][idc8],emas_v1['lat'],emas_v1['lon'],emas_v1['ref'],r_4STAR)
out_startau = stats_within_radius(dc8['G_LAT'][idc8],dc8['G_LONG'][idc8],emas_v1['lat'],emas_v1['lon'],emas_v1['tau'],r_4STAR)


# In[74]:

r_MODIS = 500.0
out_modisref = stats_within_radius(dc8['G_LAT'][idc8],dc8['G_LONG'][idc8],emas_v1['lat'],emas_v1['lon'],emas_v1['ref'],r_MODIS)
out_modistau = stats_within_radius(dc8['G_LAT'][idc8],dc8['G_LONG'][idc8],emas_v1['lat'],emas_v1['lon'],emas_v1['tau'],r_MODIS)


# In[78]:

plt.plot(out_ssfrtau['std'],'g',label='SSFR')
plt.plot(out_rsptau['std'],'c',label='RSP')
plt.plot(out_modistau['std'],'m',label='MODIS')
plt.plot(out_startau['std'],'r',label='4STAR')
plt.xlabel('points')
plt.ylabel('$\\tau$')
plt.title('Standard deviation along horizontal area')
plt.legend(frameon=False)


# In[157]:

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(122)
cf = ax.contourf(emas_v1['lon'],emas_v1['lat'],emas_v1['tau'],clevels,cmap=plt.cm.rainbow,extend='max')

ax.plot(er2['Longitude'][ier2],er2['Latitude'][ier2],'r',label='ER2 flight path')
ax.plot(dc8['G_LONG'][idc8],dc8['G_LAT'][idc8],'b--',label='DC8 flight path')

circles(dc8['G_LONG'][idc8[380]],dc8['G_LAT'][idc8[380]],radius_m2deg(dc8['G_LONG'][idc8[380]],dc8['G_LAT'][idc8[380]],r_4STAR),c='r',alpha=0.5,ax=ax,label='4STAR')
circles(er2['Longitude'][ier2[390]],er2['Latitude'][ier2[390]],radius_m2deg(er2['Longitude'][ier2[390]],er2['Latitude'][ier2[390]],r_SSFR),c='g',alpha=0.5,ax=ax,label='SSFR')
circles(er2['Longitude'][ier2[385]],er2['Latitude'][ier2[385]],radius_m2deg(er2['Longitude'][ier2[385]],er2['Latitude'][ier2[385]],r_eMAS),c='k',alpha=0.5,ax=ax,label='eMAS')
circles(er2['Longitude'][ier2[395]],er2['Latitude'][ier2[395]],radius_m2deg(er2['Longitude'][ier2[395]],er2['Latitude'][ier2[395]],r_RSP),c='c',alpha=0.5,ax=ax,label='RSP')
#circles(er2['Longitude'][ier2[375]],er2['Latitude'][ier2[375]],radius_m2deg(er2['Longitude'][ier2[375]],er2['Latitude'][ier2[375]],500.0),c='m',alpha=0.5,ax=ax,label='MODIS')
vo = [[modis['lon'][495,1208],modis['lat'][495,1208]],
      [modis['lon'][495,1209],modis['lat'][495,1209]],
      [modis['lon'][496,1209],modis['lat'][496,1209]],
      [modis['lon'][496,1208],modis['lat'][496,1208]],
      [modis['lon'][495,1208],modis['lat'][495,1208]]]
ax.add_patch(Polygon(vo,closed=True,color='m',alpha=0.5))


#ax.text(dc8['G_LONG'][idc8[380]]+0.002,dc8['G_LAT'][idc8[380]],'4STAR',color='r')
#ax.text(-93.958,21.53,'SSFR',color='g')
#ax.text(er2['Longitude'][ier2[385]]+0.002,er2['Latitude'][ier2[385]],'eMAS',color='k')
#ax.text(er2['Longitude'][ier2[395]]+0.002,er2['Latitude'][ier2[395]],'RSP',color='c')
#ax.text(er2['Longitude'][ier2[375]]+0.005,er2['Latitude'][ier2[375]],'MODIS',color='m')

#dc8toer2 = spherical_dist(np.array((er2['Latitude'][ier2[393]],er2['Longitude'][ier2[393]])),np.array((dc8['G_LAT'][idc8[380]],dc8['G_LONG'][idc8[380]])))
#ax.annotate('' , (er2['Longitude'][ier2[393]],er2['Latitude'][ier2[393]]),(dc8['G_LONG'][idc8[380]],dc8['G_LAT'][idc8[380]]), arrowprops={'arrowstyle':'<->'})
#ax.text(-93.991,21.515,'%.2f km' % dc8toer2)

#ya,yb,xa,xb = 21.46,21.54,-94.03,-93.94
#ax.set_ylim([ya,yb])
#ax.set_xlim([xa,xb])
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.set_title('eMAS $\\tau$')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

cbar = plt.colorbar(cf,fraction=0.046, pad=0.04)
cbar.set_label('$\\tau$')

ax2 = fig.add_subplot(121,sharey=ax)
ax2.plot(out_ssfrtau['std'],er2['Latitude'][ier2[::10]],'g',label='SSFR')
ax2.plot(out_rsptau['std'],er2['Latitude'][ier2[::10]],'c',label='RSP')
ax2.plot(out_modistau['std'],er2['Latitude'][ier2[::10]],'m',label='MODIS')
ax2.plot(out_startau['std'],er2['Latitude'][ier2[::10]],'r',label='4STAR')
ax2.set_ylabel('Lattude')
ax2.set_xlabel('$\\sigma\\tau$')
ax2.set_title('Variability of $\\tau$ in each Field-of-View')
plt.legend(frameon=False)
plt.tight_layout()

plt.savefig(fp+'plots/emas_tau_variations.png',dpi=600,transparent=True)


# In[124]:

print 'SSFR tau: %.2f, %.2f, max, min:%.2f, %.2f ' % (np.nanmean(out_ssfrtau['std']), np.median(out_ssfrtau['std']), np.nanmax(out_ssfrtau['std']), np.nanmin(out_ssfrtau['std']))
print 'RSP tau: %.2f, %.2f, max, min:%.2f, %.2f ' % (np.nanmean(out_rsptau['std']), np.median(out_rsptau['std']),np.nanmax(out_rsptau['std']), np.nanmin(out_rsptau['std']))
print 'MODIS tau: %.2f, %.2f, max, min:%.2f, %.2f ' % (np.nanmean(out_modistau['std']), np.median(out_modistau['std']),np.nanmax(out_modistau['std']), np.nanmin(out_modistau['std']))
print '4STAR tau: %.2f , %.2f, max, min:%.2f, %.2f' % (np.nanmean(out_startau['std']), np.median(out_startau['std']),np.nanmax(out_startau['std']), np.nanmin(out_startau['std']))


# In[125]:

print 'SSFR ref: %.2f, %.2f, max, min: %.2f,%.2f' % (np.nanmean(out_ssfrref['std']),np.median(out_ssfrref['std']),np.nanmax(out_ssfrref['std']),np.nanmin(out_ssfrref['std']))
print 'RSP ref: %.2f, %.2f, max, min: %.2f,%.2f' % (np.nanmean(out_rspref['std']),np.median(out_rspref['std']),np.nanmax(out_rspref['std']),np.nanmin(out_rspref['std']))
print 'MODIS ref: %.2f, %.2f, max, min: %.2f,%.2f' % (np.nanmean(out_modisref['std']),np.median(out_modisref['std']),np.nanmax(out_modisref['std']),np.nanmin(out_modisref['std']))
print '4STAR ref: %.2f, %.2f, max, min: %.2f,%.2f' % (np.nanmean(out_starref['std']),np.median(out_starref['std']),np.nanmax(out_starref['std']),np.nanmin(out_starref['std']))


# In[114]:

fig = plt.figure(figsize=(9,5))
ax = fig.add_subplot(122)

cf = ax.contourf(emas_v1['lon'],emas_v1['lat'],emas_v1['ref'],clevels,cmap=plt.cm.gist_earth,extend='max')

ax.plot(er2['Longitude'][ier2],er2['Latitude'][ier2],'r',label='ER2 flight path')
ax.plot(dc8['G_LONG'][idc8],dc8['G_LAT'][idc8],'b--',label='DC8 flight path')

circles(dc8['G_LONG'][idc8[380]],dc8['G_LAT'][idc8[380]],radius_m2deg(dc8['G_LONG'][idc8[380]],dc8['G_LAT'][idc8[380]],r_4STAR),c='r',alpha=0.5,ax=ax,label='4STAR')
circles(er2['Longitude'][ier2[390]],er2['Latitude'][ier2[390]],radius_m2deg(er2['Longitude'][ier2[390]],er2['Latitude'][ier2[390]],r_SSFR),c='g',alpha=0.5,ax=ax,label='SSFR')
circles(er2['Longitude'][ier2[385]],er2['Latitude'][ier2[385]],radius_m2deg(er2['Longitude'][ier2[385]],er2['Latitude'][ier2[385]],r_eMAS),c='k',alpha=0.5,ax=ax,label='eMAS')
circles(er2['Longitude'][ier2[395]],er2['Latitude'][ier2[395]],radius_m2deg(er2['Longitude'][ier2[395]],er2['Latitude'][ier2[395]],r_RSP),c='c',alpha=0.5,ax=ax,label='RSP')
#circles(er2['Longitude'][ier2[375]],er2['Latitude'][ier2[375]],radius_m2deg(er2['Longitude'][ier2[375]],er2['Latitude'][ier2[375]],500.0),c='m',alpha=0.5,ax=ax,label='MODIS')
vo = [[modis['lon'][495,1208],modis['lat'][495,1208]],
      [modis['lon'][495,1209],modis['lat'][495,1209]],
      [modis['lon'][496,1209],modis['lat'][496,1209]],
      [modis['lon'][496,1208],modis['lat'][496,1208]],
      [modis['lon'][495,1208],modis['lat'][495,1208]]]
ax.add_patch(Polygon(vo,closed=True,color='m',alpha=0.5))

ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.set_title('FOV comparison')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

cbar = plt.colorbar(cf)
cbar.set_label('$R_{eff}$ [$\\mu$m]')

ax2 = fig.add_subplot(121)
ax2.plot(out_ssfrref['std'],range(80),'g',label='SSFR')
ax2.plot(out_rspref['std'],range(80),'c',label='RSP')
ax2.plot(out_modisref['std'],range(80),'m',label='MODIS')
ax2.plot(out_starref['std'],range(80),'r',label='4STAR')
ax2.set_ylabel('points')
ax2.set_xlabel('$\\sigma R_{eff}$')
ax2.set_title('$R_{eff}$ variations in FOV')
plt.legend(frameon=False,loc=4)

plt.savefig(fp+'plots/emas_ref_variations.png',dpi=600,transparent=True)


# In[115]:

plt.hist(out_ssfrtau['std'],bins=30, histtype='stepfilled', normed=True,color='g',label='SSFR',alpha=0.6)
plt.hist(out_rsptau['std'][~isnan(out_rsptau['std'])],bins=30, histtype='stepfilled', normed=True,color='c',label='RSP',alpha=0.6)
plt.hist(out_modistau['std'][~isnan(out_modistau['std'])],bins=30, histtype='stepfilled', normed=True,color='m',label='MODIS',alpha=0.6)
plt.hist(out_startau['std'][~isnan(out_startau['std'])],bins=30, histtype='stepfilled', normed=True,color='r',label='4STAR',alpha=0.6)
plt.xlim([0,3])
plt.legend(frameon=False)

