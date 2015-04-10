# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     SEAC4RS_compare_v3
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
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - sp_v3_20130913_4STAR.out: modeled spectra output for SEAC4RS at sza 17.9, idl save file
#   - %%20130219starzen_rad.mat : special zenith radiance 4star matlab file 

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
%matplotlib inline
import numpy as np, h5py
#import plotly.plotly as py
import scipy.io as sio
import math
import os
import Sp_parameters as Sp
#py.sign_in("samuelleblanc", "4y3khh7ld4")
#import mpld3
#mpld3.enable_notbeook()

# <codecell>

import IPython
IPython.InteractiveShell.cache_size = 0

# <codecell>

# set the basic directory path
fp='C:/Users/sleblan2/Research/SEAC4RS/'

# <headingcell level=2>

# Get the lookup table for the 4STAR data

# <codecell>

# load the idl save file containing the modeled radiances
s=sio.idl.readsav(fp+'model\\sp_v3_20130913_4STAR.out')#fp+'model/sp_v1.1_20130913_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp, wvl, z, re, ta)'
# create custom key for sorting via wavelength
iwvls = np.argsort(s.zenlambda)
s.wv = np.sort(s.zenlambda)

# <codecell>

if 'Sp' in locals():
    reload(Sp)
if 'lut' in locals():
    del lut
    import gc; gc.collect()

# <codecell>

lut = Sp.Sp(s,irrad=True)
lut.params()
lut.param_hires()

# <codecell>

lut.sp_hires()

# <codecell>

print lut.tau.shape
print lut.ref.shape
print lut.sp.ndim
print lut.par.size
print lut.par.shape
print lut.ref
print lut.par.shape

# <codecell>

fig3,ax3 = plt.subplots(5,3,sharex=True,figsize=(15,8))
ax3 = ax3.ravel()

for i in range(lut.npar-1):
    color.cycle_cmap(len(lut.ref[lut.ref<30]),cmap=plt.cm.RdBu,ax=ax3[i])
    for j in xrange(len(lut.ref)):
        ax3[i].plot(lut.tau,lut.par[0,j,:,i])
    ax3[i].set_title('Parameter '+str(i))
    ax3[i].grid()
    ax3[i].set_xlim([0,100])
    if i > 11: 
        ax3[i].set_xlabel('Tau')

fig3.tight_layout()
plt.suptitle('Liquid')
plt.subplots_adjust(top=0.93,right=0.93)

cbar_ax = fig3.add_axes([0.95,0.10,0.02,0.8])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.RdBu,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.ref[lut.ref<30])
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
cba.ax.set_ylabel('R$_{ef}$ [$\\mu$m]')
cba.ax.set_yticklabels(np.linspace(lut.ref[0],29,6));

plt.show()

# <codecell>

fig4,ax4 = plt.subplots(5,3,sharex=True,figsize=(15,8))
ax4 = ax4.ravel()

for i in range(lut.npar-1):
    color.cycle_cmap(len(lut.tau),cmap=plt.cm.gist_ncar,ax=ax4[i])
    for j in xrange(len(lut.tau)):
        ax4[i].plot(lut.ref,lut.par[0,:,j,i])
    ax4[i].set_title('Parameter '+str(i))
    ax4[i].grid()
    ax4[i].set_xlim([0,30])
    if i > 11: 
        ax4[i].set_xlabel('R$_{ef}$ [$\mu$m]')

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

plt.show()

# <codecell>

fig5,ax5 = plt.subplots(5,3,sharex=True,figsize=(15,8))
ax5 = ax5.ravel()

for i in range(lut.npar-1):
    color.cycle_cmap(len(lut.ref),cmap=plt.cm.RdBu,ax=ax5[i])
    for j in xrange(len(lut.ref)):
        ax5[i].plot(lut.tau,lut.par[1,j,:,i])
    ax5[i].set_title('Parameter '+str(i))
    ax5[i].grid()
    ax5[i].set_xlim([0,100])
    if i > 11: 
        ax5[i].set_xlabel('Tau')

fig5.tight_layout()
plt.suptitle('Ice')
plt.subplots_adjust(top=0.93,right=0.93)
cbar_ax = fig5.add_axes([0.95,0.10,0.02,0.8])
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.RdBu,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.ref)
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
cba.ax.set_ylabel('R$_{ef}$ [$\\mu$m]')
cba.ax.set_yticklabels(np.linspace(lut.ref[0],lut.ref[-1],6));
plt.show()

# <codecell>

fig6,ax6 = plt.subplots(5,3,sharex=True,figsize=(15,8))
ax6 = ax6.ravel()
for i in range(lut.npar-1):
    color.cycle_cmap(len(lut.tau),cmap=plt.cm.gist_ncar,ax=ax6[i])
    for j in xrange(len(lut.tau)):
        ax6[i].plot(lut.ref,lut.par[1,:,j,i])
    ax6[i].set_title('Parameter '+str(i))
    ax6[i].grid()
    ax6[i].set_xlim([0,50])
    if i > 11: 
        ax6[i].set_xlabel('R$_{ef}$ [$\mu$m]')

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
plt.show()

# <headingcell level=2>

# Get the 4STAR data

# <codecell>

import load_modis
reload(load_modis)
from load_modis import mat2py_time, toutc, load_ict

# <codecell>

dc8 = load_ict(fp+'dc8/20130913/SEAC4RS-MMS-1HZ_DC8_20130913_RB.ict')

# <codecell>

plt.figure()
plt.plot(dc8['TIME_UTC'],dc8['G_ALT'])
plt.xlabel('UTC [h]')
plt.ylabel('GPS Altitude [m]')
plt.title('DC8 Altitude on 2013-09-13')

# <codecell>

print dc8['TIME_UTC'][12000]/3600.0
print dc8['G_ALT'][12000]

# <codecell>

# load the matlab file containing the measured TCAP radiances
mea = sio.loadmat(fp+'../4STAR/SEAC4RS/20130913/20130913starzen_3.mat')
mea.keys()

# <markdowncell>

# Go through and get the radiances for good points, and for the time selected

# <codecell>

print mea['t']
tt = mat2py_time(mea['t'])
mea['utc'] = toutc(tt)

# <codecell>

mea['good'] = np.where((mea['utc']>18.5) & (mea['utc']<19.75) & (mea['Str'].flatten()!=0) & (mea['sat_time'].flatten()==0))

# <codecell>

plt.figure()
plt.plot(mea['utc'],abs(mea['rad'][:,1024]-mea['rad'][:,1068]))
plt.xlabel('UTC [Hours]')
plt.ylabel('difference in radiance at 1024 nm')

# <codecell>

plt.figure()
plt.plot(mea['utc'],mea['rad'][:,400]/1000.)
plt.xlabel('UTC [hours]')
plt.ylabel('Radiance at 400 nm [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')

# <codecell>

reload(Sp)
if 'meas' in locals():
    del meas
    import gc; gc.collect()
# first convert measurements to Sp class, with inherent parameters defined
meas = Sp.Sp(mea)
meas.params()

# <headingcell level=2>

# Run the retrieval on 4STAR data

# <codecell>

import run_kisq_retrieval as rk
reload(rk)

# <codecell>

subp = [1,2,4,5,6,10,11,13]
subp = [2,3,5,6,7,11,12,14]

# <codecell>

subp = [0,1,3,4,5,6,7,10,11,13,14]

# <codecell>

print max(meas.good)
print type(mea['good'])
print isinstance(mea['good'],tuple)
print type(mea['good'][0])

# <codecell>

(meas.tau,meas.ref,meas.phase,meas.ki) = rk.run_retrieval(meas,lut,subp)

# <codecell>

print meas.utc.shape
print len(meas.good), max(meas.good)

# <codecell>

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
plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results_v3.png',dpi=600)
plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results_v3.pdf',bbox='tight')

# <codecell>

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
plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results_v3.png',dpi=600)
plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results_v3.pdf',bbox='tight')

# <markdowncell>

# Smooth out the retrieved 4STAR data

# <codecell>

meas.tau[meas.good] = smooth(meas.tau[meas.good],20)
meas.ref[meas.good] = smooth(meas.ref[meas.good],20)

# <headingcell level=2>

# Get SSFR data from ER2

# <codecell>

import load_modis as lm
if 'lm' in locals():
    reload(lm)
from load_modis import load_ict
ssfr_er2_file = fp+'er2/20130913/SEAC4RS-SSFR_ER2_20130913_R0.ict'
ssfr_er2 = load_ict(ssfr_er2_file)
print len(ssfr_er2['UTC'])

# <codecell>

nasdat_er2_file = fp+'er2/20130913/seac4rs-nasdat_er2_20130913_r0.ict'
er2 = load_ict(nasdat_er2_file)
print len(er2['Start_UTC'])
print er2['Start_UTC']
print np.min(er2['Solar_Zenith_Angle'])

# <codecell>

plt.figure()
feet2meter=0.3048
plt.plot(er2['Start_UTC'],er2['GPS_Altitude'],label="GPS",color='b')
plt.plot(er2['Start_UTC'],er2['Pressure_Altitude']*feet2meter,label="Pressure",color='g')
plt.grid(True)
plt.xlabel('Time [UTC]')
plt.ylabel('ER-2 Altitude [m]')
plt.legend(frameon=False)

# <codecell>

print er2['Start_UTC'][12000]
print er2['GPS_Altitude'][12000]
print er2['Pressure_Altitude'][12000]

# <codecell>

ssfr_idl_file = fp+'er2/20130913/20130913_calibspcs.out'
ssfr_idl = sio.idl.readsav(ssfr_idl_file)
print ssfr_idl.keys()
print np.shape(ssfr_idl['zspectra'])
print np.shape(ssfr_idl['sat'])

# <codecell>

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

# <codecell>

print i500, i1650
i1600 = np.argmin(abs(ssfr_idl['zenlambda']-1600.0))

# <codecell>

iu = np.argmin(abs(ssfr_idl['tmhrs']-18.5))
iu2 = np.argmin(abs(ssfr_er2['UTC']-18.5))
print iu,iu2

# <codecell>

print ssfr_er2['UP500'][iu2], ssfr_idl['zspectra'][iu,i500], ssfr_er2['UP500'][iu2]/ssfr_idl['zspectra'][iu,i500]
print ssfr_er2['DN500'][iu2+200], ssfr_idl['nspectra'][iu+200,i500], ssfr_er2['DN500'][iu2+200]/ssfr_idl['nspectra'][iu+200,i500]
print ssfr_er2['UTC'][iu2+200]
print ssfr_idl['tmhrs'][iu+200]
print ssfr_idl['nspectra'][iu+200,i500]
rup500 = ssfr_er2['UP500'][iu2]/ssfr_idl['zspectra'][iu,i500]
rdn500 = ssfr_er2['DN500'][iu2+200]/ssfr_idl['nspectra'][iu+200,i500]

# <codecell>

print ssfr_er2['UP1600'][iu2], ssfr_idl['zspectra'][iu,i1600], ssfr_er2['UP1600'][iu2]/ssfr_idl['zspectra'][iu,i1600]
print ssfr_er2['DN1600'][iu2+200], ssfr_idl['nspectra'][iu+200,i1600], ssfr_er2['DN1600'][iu2+200]/ssfr_idl['nspectra'][iu+200,i1600]
print ssfr_er2['UTC'][iu2+200]
print ssfr_idl['tmhrs'][iu+200]
print ssfr_idl['nspectra'][iu+200,i500]
rup1600 = ssfr_er2['UP1600'][iu2]/ssfr_idl['zspectra'][iu,i1600]
rdn1600 = ssfr_er2['DN1600'][iu2+200]/ssfr_idl['nspectra'][iu+200,i1600]

# <codecell>

print rup500/rdn500

# <codecell>

ssfr.Rvis = ssfr.Rvis*rup500/rdn500
ssfr.Rnir = ssfr.Rnir*rup1600/rdn1600

# <codecell>

print np.nanmin(ssfr.Rvis), np.nanmax(ssfr.Rvis)
print ssfr_idl['nadlambda'][i1650-2]
print i500, i1650
print len(ssfr.utc), len(ssfr.Rvis)

# <codecell>

from scipy import interpolate
sza_fx = interpolate.interp1d(er2['Start_UTC'],er2['Solar_Zenith_Angle'],bounds_error=False)
ssfr.sza = sza_fx(ssfr.utc)
print sza_fx(18.41)
print ssfr.utc[17000], ssfr.sza[17000]

# <codecell>

fig = plt.figure()
fig.add_subplot(1,2,1)
plt.plot(ssfr.utc,ssfr_idl['zspectra'][:,i500])
plt.plot(ssfr.utc,ssfr_idl['nspectra'][:,i500],'r')
plt.xlim([17.8,19.2])
fig.add_subplot(1,2,2)
plt.plot(ssfr.utc,ssfr.Rvis,'b',ssfr.utc,ssfr.Rnir,'r')
plt.xlim([17.8,19.2])
#plt.ylim([0.22,0.32])

# <codecell>

rad = np.linspace(0,pi,50)
plt.plot(rad,np.cos(rad))
plt.plot(rad,np.cos(rad).cumsum()/np.cos(rad).cumsum().max())
print rad.cumsum()
print np.sum(np.cos(rad[0:24]))/np.sum(np.cos(rad))

# <headingcell level=5>

# Check the lut for reflectance

# <codecell>

lut.sp_hires(doirrad=True)

# <codecell>

lut.reflect.shape

# <codecell>

figref,axref = plt.subplots(10,3,sharex=True,figsize=(15,16))
axref = axref.ravel()
for i in range(lut.ref.size/2-1):
    color.cycle_cmap(len(lut.tau),cmap=plt.cm.gist_ncar,ax=axref[i])
    for j in xrange(len(lut.tau)):
        axref[i].plot(lut.wvl,lut.reflect[1,:,1,2*i,j])
    axref[i].set_title('Ref i='+str(2*i))
    axref[i].grid()
    axref[i].set_xlim([350,1700])
    axref[i].set_ylim([0,1])
    if i > 26: 
        axref[i].set_xlabel('Wavelength [nm]')
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
plt.show()

# <codecell>

figref,axref = plt.subplots(10,3,sharex=True,figsize=(15,16))
axref = axref.ravel()
for i in range(lut.ref.size/2-1):
    color.cycle_cmap(len(lut.tau),cmap=plt.cm.gist_ncar,ax=axref[i])
    for j in xrange(len(lut.tau)):
        axref[i].plot(lut.wvl,lut.reflect[0,:,1,2*i,j])
    axref[i].set_title('Ref i='+str(2*i))
    axref[i].grid()
    axref[i].set_xlim([350,1700])
    axref[i].set_ylim([0,1])
    if i > 26: 
        axref[i].set_xlabel('Wavelength [nm]')
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
plt.show()

# <headingcell level=3>

# Now run the retrieval on reflectance from SSFR measurements based on the ER2 platorm

# <codecell>

import run_2wvl_retrieval as rw
if 'rw' in locals():
    reload(rw)

# <codecell>

(ssfr.tau,ssfr.ref,ssfr.ki) = rw.run_2wvl_retrieval(ssfr,lut,wvls=[500.0,1600.0])

# <codecell>

ssfr.tau[ssfr.ref==60] = np.nan
ssfr.ref[ssfr.ref==60] = np.nan

# <codecell>

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

# <headingcell level=2>

# Load data from CPL to caompare flight profiles of DC8 and ER2 to cloud layers

# <codecell>

cpl_layer_file = fp+'er2\\20130913\\layers_13965_13sep13.txt'
import load_modis as lm
reload(lm)
from load_modis import load_cpl_layers
cpl_layers = load_cpl_layers(cpl_layer_file)

# <codecell>

cpl_layers.dtype

# <codecell>

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

# <codecell>

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

# <codecell>

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

# <headingcell level=2>

# Get the data from MODIS to compare

# <codecell>

from mpl_toolkits.basemap import Basemap,cm
myd06_file = fp+'modis\\20130913\\MYD06_L2.A2013256.1910.006.2014267222159.hdf'
myd03_file = fp+'modis\\20130913\\MYD03.A2013256.1910.006.2013257153624.hdf'
print os.path.isfile(myd03_file) #check if it exists
print os.path.isfile(myd06_file)

# <codecell>

import load_modis as lm
reload(lm)
if 'modis' in locals():
    del modis, modis_dicts
    import gc; gc.collect()

# <codecell>

modis,modis_dicts = lm.load_modis(myd03_file,myd06_file)

# <markdowncell>

# Now plot the resulting imagery

# <codecell>

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

# <codecell>

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

cs1 = m1.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.gist_ncar,extend='max')
cbar = m1.colorbar(cs1)
cbar.set_label('$\\tau$')
axm2[0].set_title('MODIS - AQUA Cloud optical Thickness')
x1,y1 = m1(meas.lon,meas.lat)
xer2,yer2 = m1(er2['Longitude'],er2['Latitude'])
xdc8,ydc8 = m1(dc8['G_LONG']/100000.0,dc8['G_LAT']/100000.0)
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

# <headingcell level=2>

# Import eMAS values

# <codecell>

if 'lm' in locals():
    reload(lm)
from load_modis import load_emas, load_hdf

# <codecell>

emas_file = fp+'er2/20130913/EMASL2_13965_13_20130913_1905_1918_V00.hdf'
print os.path.isfile(emas_file)

# <codecell>

emas,emas_dicts = load_hdf(emas_file)

# <codecell>

emas_values = (('lat',0),('lon',1),('tau',15),('ref',23),('phase',58),('layer',59),('qa',68))
emas,emas_dicts = load_hdf(emas_file,values=emas_values)

# <codecell>

plt.figure()
plt.plot(emas['tau'])

# <markdowncell>

# Now Redo the load of emas data, but with the new V01 files, which includes the newest calibration as of 20150122, which is considered final for SEAC4RS. thermal band revisions (3.7 um and higher) may still occur.  
# 
# From Tom Arnold 20150120: While the calibration for the V01 data is considered final, some minor revision may still be made to the themal band (3.7um and higher) radiance data for two or three of the tracks I have given you. Separate from the calibration process, filtering (for a coherent noise problem) has been applied to all the eMAS thermal bands (bands 26-38).  Evaluation of the quality of the filtering has shown for some eMAS tracks, some additional filtering is still required (and will likely affect two or three of the tracks I have given you).  I will make that data available when it is completed for the tracks you need, though for the thick cirrus in the tracks you are interested in,  I don’t expect much impact to the cloud retrievals (most impact will be relatively small adjustment to some of the cloud top property data - such as cloud top temperature or pressure). I expect the re-filtered data to be available in the next couple weeks.

# <codecell>

emas_file_v1 = fp+'emas/20130913/EMASL2_13965_13_20130913_1905_1918_V01.hdf'
print os.path.isfile(emas_file_v1)

# <codecell>

print fp
print emas_file_v1

# <codecell>

emas_v1,emas_dicts_v1 = load_hdf(emas_file_v1)

# <codecell>

emas_v1,emas_dicts_v1 = load_hdf(emas_file_v1, values=emas_values)

# <codecell>

emas_dicts_v1['tau']

# <codecell>

reload(lm)
from load_modis import map_ind
dc8_ind = map_ind(emas['lon'],emas['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good'][0])

# <codecell>

print np.shape(dc8_ind)

# <codecell>

print dc8_ind[0,-1],dc8_ind[1,-1]
print emas['lon'][388,715],emas['lat'][388,715]
print emas['lon'][388,714],emas['lat'][388,714]

sdist= lambda lon1,lat1,lon2,lat2:1000.0 * 3958.75 * np.arccos(np.cos(np.radians(lat1)-np.radians(lat2)) - np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * (1 - np.cos(np.radians(lon1)-np.radians(lon2))))

# <codecell>

print '%20.17f' % sdist(emas['lon'][388,715],emas['lat'][388,715],emas['lon'][385,714],emas['lat'][385,714])

# <markdowncell>

# Do the mods calculations

# <codecell>

reload(lm)
from load_modis import map_ind
dc8_ind_modis = map_ind(modis['lon'],modis['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good'][0])

# <headingcell level=2>

# Load the cloud probe data

# <codecell>

prb_file = fp+'dc8/20130913/SEAC4RS_20130913_Reff.txt'
probes = np.genfromtxt(prb_file,skip_header=2)

# <codecell>

print probes[:,1]
print probes[:,2]

# <codecell>

print probes.shape
print probes[:,7]
plt.figure()
plt.hist(probes[:,7])

# <headingcell level=2>

# Load 2DS data for effective radius at specific times

# <codecell>

twoDS = load_ict(fp+'dc8/20130913/seac4rs-2DS_DC8_20130913_R0.ict')

# <codecell>

plt.figure
plt.plot(twoDS['Start_UTC'],twoDS['effectiveD']/2.0)
plt.plot(twoDS['Start_UTC'],smooth(twoDS['effectiveD']/2.0,60),'r')
plt.ylim([0,100])
plt.xlabel('UTC [h]')
plt.ylabel('Effective radius [$\\mu$m]')
plt.title('2DS effective radius from 2013-09-13')

# <codecell>

len(twoDS['Start_UTC'])

# <codecell>

len(dc8['G_ALT'])

# <codecell>

plt.figure
plt.plot(twoDS['effectiveD'][1:]/2.0,dc8['G_ALT'])
plt.xlim([0,100])
plt.xlabel('Effective radius [$\\mu$m]')
plt.ylabel('Altitude [m]')

# <headingcell level=2>

# Load RSP results

# <codecell>

rsp = load_ict(fp+'er2/20130913/SEAC4RS-RSP-ICECLD_ER2_20130913_R2.ict')
print rsp['Phase']
rsp_good = np.where((rsp['UTC']>17.8) & (rsp['UTC']<19.2) & (rsp['Phase']==0) & (rsp['COT']>0) & (rsp['R_eff159']>0))
print np.shape(rsp_good)
print len(rsp_good[0])

# <codecell>

plt.plot(rsp['UTC'][rsp_good[0]],rsp['R_eff159'][rsp_good[0]])

# <codecell>

print rsp_good[0][-1]

# <codecell>

# get the distance between two points for RSP
from load_modis import spherical_dist
print rsp['Lat'][rsp_good[0][-1]], rsp['Lon'][rsp_good[0][-1]]
print rsp['Lat'][rsp_good[0][-2]], rsp['Lon'][rsp_good[0][-2]]
print spherical_dist(np.array(rsp['Lat'][rsp_good[0][-1]],rsp['Lon'][rsp_good[0][-1]]),
                     np.array(rsp['Lat'][rsp_good[0][-2]],rsp['Lon'][rsp_good[0][-2]]))

# <headingcell level=2>

# Calculate the histogram for each comparisons

# <codecell>

from Sp_parameters import nanmasked
modis_tau,im = nanmasked(modis['tau'][dc8_ind_modis[0,:],dc8_ind_modis[1,:]])
emas_tau,ie = nanmasked(emas['tau'][dc8_ind[0,:],dc8_ind[1,:]])
emas_tau_v1,ie1 = nanmasked(emas_v1['tau'][dc8_ind[0,:],dc8_ind[1,:]])
modis_ref,im = nanmasked(modis['ref'][dc8_ind_modis[0,:],dc8_ind_modis[1,:]])
emas_ref,ie = nanmasked(emas['ref'][dc8_ind[0,:],dc8_ind[1,:]])
emas_ref_v1,ie1 = nanmasked(emas_v1['ref'][dc8_ind[0,:],dc8_ind[1,:]])
star_tau,ist = nanmasked(meas.tau[meas.good])
star_ref,ist = nanmasked(meas.ref[meas.good])
ssfr.good = np.where((ssfr.utc>17.8)&(ssfr.utc<19.2))
ssfr_tau,iss = nanmasked(ssfr.tau[ssfr.good[0],1])
ssfr_ref,iss = nanmasked(ssfr.ref[ssfr.good[0],1])
rsp_tau,irs = nanmasked(rsp['COT'][rsp_good[0]])
rsp_ref,irs = nanmasked(rsp['R_eff159'][rsp_good[0]])

# <codecell>

star_g = np.where((meas.utc>19.0)&(meas.utc<19.2)&isfinite(meas.tau))
print '4STAR',meas.tau[star_g[0][0]], meas.ref[star_g[0][0]]
ssfr_g = np.where((ssfr.utc>19.0)&(ssfr.utc<19.2)&isfinite(ssfr.tau[:,1]))
print 'ssfr reflect',ssfr.tau[ssfr_g[0][0],1],ssfr.ref[ssfr_g[0][0],1]
rsp_g = np.where((rsp['UTC']>19.0)&(rsp['UTC']<19.2)&isfinite(rsp['COT']))
print 'rsp',rsp['COT'][rsp_g[0][0]],rsp['R_eff159'][rsp_g[0][0]]
print 'emas',emas['tau'][dc8_ind[0,0],dc8_ind[1,0]],emas['ref'][dc8_ind[0,0],dc8_ind[1,0]]
print 'emas_v1',emas_v1['tau'][dc8_ind[0,0],dc8_ind[1,0]],emas_v1['ref'][dc8_ind[0,0],dc8_ind[1,0]]
print 'modis',modis['tau'][dc8_ind_modis[0,0],dc8_ind_modis[1,0]],modis['ref'][dc8_ind_modis[0,0],dc8_ind_modis[1,0]] 

# <codecell>

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

# <markdowncell>

# Now we present the histogram of retrieved cloud properties for each instrument
# The instruments have different cross sectional area at cloud height that they sample
# These are:, with the ratio of smoothing:
#  - MODIS : 500m : 6
#  - emas : 50m : 60
#  - SSFR(reflectance) : 3000m : 1
#  - RSP: 42 m : 70
#  - 4STAR: 70 m : 40

# <headingcell level=2>

# Run through data and convolve to get the same area

# <markdowncell>

# ER-2 Altitude: 19018 m
# DC8 Altitude: 8047 m
# Cloud base height: 
# Cloud top height: 
#  - Modis: 500 m
#  - eMAS: 50m or 2.5 mrad FOV ()
#  - RSP: 14 mrad FOV
#  - 4STAR: 2 degree, 0.0349 rads
#  - GOES: 1km
#  - SSFR: 180 degrees, Pi rads

# <headingcell level=2>

# Plot histogram of different tau and ref comparison

# <codecell>

plt.figure(figsize=(9,6))
plt.axvspan(0,80,color='#FFFFFF')
plt.hist(smooth(modis_tau,6),bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='Modis (Reflected)',range=(0,40))
#plt.hist(smooth(emas_tau,60),bins=30, histtype='stepfilled', normed=True, color='b',alpha=0.6, label='eMAS (Reflected)',range=(0,40))
plt.hist(smooth(emas_tau_v1,60),bins=30, histtype='stepfilled', normed=True, color='k',alpha=0.6, label='eMAS (Reflected)',range=(0,40))
plt.hist(smooth(ssfr_tau,2),bins=20, histtype='stepfilled', normed=True, color='g',alpha=0.6, label='SSFR (Reflected)',range=(5,30))
plt.hist(smooth(rsp_tau,70),bins=30, histtype='stepfilled', normed=True, color='c',alpha=0.6, label='RSP (Reflected)',range=(0,40))
plt.hist(smooth(star_tau,40),bins=30, histtype='stepfilled', normed=True, color='r',alpha=0.6, label='4STAR (Transmitted)',range=(0,40))
plt.title('Optical Thickness histogram')
plt.ylabel('Normed probability')
plt.xlabel('$\\tau$')
plot_median_mean(smooth(modis_tau,6),color='m')
#plot_median_mean(smooth(emas_tau ,60),color='b')
plot_median_mean(smooth(emas_tau_v1,60),color='k')
plot_median_mean(smooth(ssfr_tau[(ssfr_tau>5)&(ssfr_tau<30)],2),color='g')
plot_median_mean(smooth(star_tau,40),color='r',lbl=True)
plot_median_mean(smooth(rsp_tau,70),color='c')
plt.legend(frameon=False)
plt.xlim([0,60])
plt.savefig(fp+'plots/hist_modis_4star_tau_v3.png',dpi=600,transparent=True)
#plt.savefig(fp+'plots/hist_modis_4star_tau.pdf',bbox='tight')

# <codecell>

plt.figure(figsize=(9,6))
plt.axvspan(0,80,color='#FFFFFF')
plt.hist(smooth(modis_ref,6),bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='Modis (Reflected)',range=(0,59))
#plt.hist(smooth(emas_ref,60),bins=30, histtype='stepfilled', normed=True, color='b',alpha=0.6, label='eMAS (Reflected)',range=(0,59))
plt.hist(smooth(emas_ref_v1,60),bins=30, histtype='stepfilled', normed=True, color='k',alpha=0.6, label='eMAS (Reflected)',range=(0,59))
plt.hist(smooth(ssfr_ref,2),bins=30, histtype='stepfilled', normed=True, color='g',alpha=0.6, label='SSFR (Reflected)',range=(0,59))
plt.hist(smooth(rsp_ref,70),bins=30, histtype='stepfilled', normed=True, color='c',alpha=0.6, label='RSP (Reflected)',range=(0,79))
plt.hist(star_ref,bins=30, histtype='stepfilled', normed=True, color='r',alpha=0.6, label='4STAR (Transmitted)',range=(0,59))
plt.hist(probes[:,7],bins=20, histtype='stepfilled', normed=True, color='y',alpha=0.6, label='Cloud probes (In Situ)',range=(0,79))
plt.title('Cloud particle effective radius histogram')
plt.ylabel('Normed probability')
plt.xlabel('R$_{eff}$ [$\\mu$m]')
plot_median_mean(smooth(modis_ref,6),color='m')
#plot_median_mean(smooth(emas_ref,60),color='b')
plot_median_mean(smooth(emas_ref_v1,60),color='k')
plot_median_mean(smooth(ssfr_ref,2),color='g')
plot_median_mean(smooth(rsp_ref,70),color='c')
plot_median_mean(probes[:,7],color='y')
plot_median_mean(star_ref,lbl=True,color='r')
plt.legend(frameon=False,loc='upper right')
plt.xlim([10,80])
plt.savefig(fp+'plots/hist_modis_4star_ref_v3.png',dpi=600,transparent=True)

# <headingcell level=2>

# eMAS V00 and V01 comparison

# <codecell>

emas_v1.keys()

# <codecell>

plt.figure();
plt.plot(emas_v1['lon'][dc8_ind[0,:],dc8_ind[1,:]],emas_tau_v1,label='V01')
plt.plot(emas['lon'][dc8_ind[0,:],dc8_ind[1,:]], emas_tau,label='V00')


# <headingcell level=1>

# 1:1 relationship

# <codecell>

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

# <codecell>

print emas_tau_v1, emas_tau
print linfit(emas_tau_v1,emas_tau)

# <codecell>

plt.figure()
plt.hist(emas_tau_v1-emas_tau,bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='eMAS tau difference',range=(-5,5))
plt.xlabel(r'$\tau$ difference')
plt.ylabel('Normed probability')
plt.title(r'eMAS $\tau$ difference (V01 - V00)')
plt.savefig(fp+'plots/emas_diff_histogram_v01_v00_tau.png',dpi=600,transparent=True)
np.max(emas_tau_v1-emas_tau)

# <codecell>

plt.figure()
plt.plot(emas_ref_v1,emas_ref,'+',label='eMAS $r_{ef}$')
plt.plot([22,36],[22,36],'k--', label='one-to-one')
plt.title('eMAS version comparison along DC8 flight track on 2013-09-13')
plt.xlabel('eMAS V01 $r_{ef}$ [$\mu$m]')
plt.ylabel('eMAS V00 $r_{ef}$ [$\mu$m]')
plt.legend(frameon=False,loc=4)
plt.savefig(fp+'plots/emas_v00_compare_v01_ref.png',dpi=600,transparent=True)

# <codecell>

plt.figure()
plt.hist(emas_ref_v1-emas_ref,bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.6, label='eMAS Ref difference',range=(-5,5))
plt.xlabel('r$_{ef}$ difference [$\mu$m]')
np.max(emas_ref_v1-emas_ref)

# <headingcell level=2>

# Find the mean tau and ref for each instrument

# <codecell>

print np.nanmean(modis_ref)

# <codecell>

print np.shape(modis_tau)

# <codecell>

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

# <codecell>

print lut

# <codecell>

print lut.sp_irrdn.shape
print lut.tausp.shape
print lut.refsp.shape
print lut.wvl.shape

# <markdowncell>

# interpolate the modeled irradiance at z=0 (dc8 altitude) to the retrieved values of tau and ref

# <codecell>

from scipy import interpolate
#fx = interpolate.RectBivariateSpline(ref[refranges[ph]],tau,sp[ph,w,z,refranges[ph],:],kx=1,ky=1)
#sp_hires[ph,w,z,:,:] = fx(ref_hires,tau_hires)

# <codecell>

def interpirr(re,ta):
    "interpolate over irradiance to get the spectrum at ta and re"
    sp = np.zeros(lut.wvl.size)
    print re,ta
    for w in xrange(lut.wvl.size):
        irrdnfx = interpolate.RectBivariateSpline(lut.refsp[lut.refsp>5],lut.tausp,lut.sp_irrdn[1,w,0,lut.refsp>5,:],kx=1,ky=1)
        sp[w] = irrdnfx(re,ta)
    return smooth(sp,20)

# <codecell>

sp_modis = interpirr(means.ref.modis,means.tau.modis)
sp_emas = interpirr(means.ref.emas,means.tau.emas)
sp_ssfr = interpirr(means.ref.ssfr,means.tau.ssfr)
sp_star = interpirr(means.ref.star,means.tau.star)
sp_rsp = interpirr(means.ref.rsp,means.tau.rsp)

# <codecell>

print list(enumerate(lut.refsp))
print list(enumerate(lut.tausp))
print lut.sp_irrdn.shape

# <codecell>

plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,13,13]*10.0,'m',label='modis')
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,10,12]*10.0,'b',label='emas')
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,29,17]*10.0,'g',label='SSFR')
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,24,11]*10.0,'r',label='RSP')
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,9,13]*10.0,'c',label='4STAR')
#plt.plot(ssfr_dc8[''])
plt.legend(frameon=False)

# <codecell>

print sp_rsp

# <markdowncell>

# Load the SSFR from the DC8

# <codecell>

ssfr_dc8 = sio.idl.readsav(fp+'dc8/20130913/20130913_calibspcs.out')
print ssfr_dc8.keys()
iutc185 = np.nanargmin(abs(ssfr_dc8['tmhrs']-18.5))
iutc192 = np.nanargmin(abs(ssfr_dc8['tmhrs']-19.2))
print ssfr_dc8['zspectra'].shape
dn = np.nanmean(ssfr_dc8['zspectra'][iutc185:iutc192,:],axis=0)
print dn.shape

# <codecell>

ssfr_dc8_ict = load_ict(fp+'dc8/20130913/SEAC4RS-SSFR_DC8_20130913_R0.ict')

# <codecell>

iutc185tt = np.nanargmin(abs(ssfr_dc8_ict['UTC']-18.5))
iutc192tt = np.nanargmin(abs(ssfr_dc8_ict['UTC']-19.2))
iutc190tt = np.nanargmin(abs(ssfr_dc8_ict['UTC']-19.0))
ssfr_dc8_ict_good = np.where((ssfr_dc8_ict['UTC']>19.0) & (ssfr_dc8_ict['UTC']<19.2) & (ssfr_dc8_ict['DN500']>0))

# <codecell>

s500 = np.nanmean(ssfr_dc8_ict['DN500'][ssfr_dc8_ict_good[0]])
print s500
plt.plot(ssfr_dc8_ict['UTC'][ssfr_dc8_ict_good[0]],ssfr_dc8_ict['DN500'][ssfr_dc8_ict_good[0]])

# <codecell>

idc500 = np.nanargmin(abs(ssfr_dc8['zenlambda']-500.0))
r500 = s500/dn[idc500]
print r500
fssp = interpolate.interp1d(ssfr_dc8['zenlambda'],dn*r500,bounds_error=False)
ssp = fssp(lut.wvl)

# <codecell>

plt.figure(figsize=(13,10))
f0,a0 = plt.subplots(2,1,sharex=True)
a0[0].plot(lut.wvl,ssp,label='SSFR Measurement',c='k')
a0[0].axvspan(350,1700,color='#FFFFFF')
a0[1].axvspan(350,1700,color='#FFFFFF')
a0[0].plot(lut.wvl,sp_modis,c='m',label='Modis')
a0[0].plot(lut.wvl,sp_emas,c='b',label='eMAS')
a0[0].plot(lut.wvl,sp_ssfr,c='g',label='SSFR (reflected)')
a0[0].plot(lut.wvl,sp_rsp/10.0,c='c',label='RSP')
a0[0].plot(lut.wvl,sp_modis,c='r',label='4STAR')
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

# <headingcell level=2>

# Get GOES - Problematic now...

# <codecell>

goes_file = fp+'er2/20130913/G13V04.0.CONUS.2013256.2045.PX.04K.NC'
goes,goes_dicts = load_hdf(goes_file)

# <codecell>

goes_values = (('lat',6),('lon',7),('tau',19),('ref',20))
goes,goes_dicts = load_hdf(goes_file,values=goes_values)

# <codecell>

print goes['lat']

# <codecell>

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

