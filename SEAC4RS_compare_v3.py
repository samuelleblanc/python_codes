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
#     - plotting_utils (user defined plotting routines)
#     - map_utils, dependent on geopy
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
s=sio.idl.readsav(fp+'model\\sp_v2_20130913_4STAR.out')#fp+'model/sp_v1.1_20130913_4STAR.out')
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

from Sp_parameters import smooth

# <codecell>

import run_kisq_retrieval as rk
reload(rk)

# <codecell>

subp = [1,2,4,5,6,10,11,13]
#subp = [2,3,5,6,7,11,12,14]

# <codecell>

print max(meas.good)
print type(mea['good'])
print isinstance(mea['good'],tuple)
print type(mea['good'][0])

# <codecell>

(meas.tau,meas.ref,meas.phase,meas.ki) = rk.run_retrieval(meas,lut)

# <codecell>

print meas.utc.shape
print len(meas.good), max(meas.good)

# <codecell>

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

# <codecell>

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

# <codecell>

ia = abs(cpl_layers['utc']-19.03).argmin()
print cpl_layers['top'][ia,0]
print cpl_layers['bot'][ia,:]

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

emas_values = (('lat',0),('lon',1),('tau',15),('ref',23),('phase',58),('layer',59),('qa',68))
emas_v1,emas_dicts_v1 = load_hdf(emas_file_v1, values=emas_values)

# <codecell>

emas_dicts_v1['tau']

# <codecell>

from map_utils import map_ind
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

from map_utils import map_ind
dc8_ind_modis = map_ind(modis['lon'],modis['lat'],mea['Lon'],mea['Lat'],meas_good=mea['good'][0])

# <headingcell level=2>

# Load APR-2 HDF files for DC8 radar images

# <codecell>

fa = fp+'dc8/20130913//SEAC4RS-APR2_DC8_20130913/SEAC4RS-APR2_DC8_20130913'
fe = '_R23.h4'
files = ['180527','181019','182329','184933','190145','192149','194031']
aprfiles = [fa+s+fe for s in files]

# <codecell>

aprfiles

# <codecell>

reload(lm)
from load_modis import load_apr
apr = load_apr(aprfiles)

# <codecell>

levels = np.arange(0,7000,30)
plt.figure
csz = plt.contourf(apr['lonz'],apr['altflt'],apr['dbz'],levels,cmap=plt.cm.jet)
plt.colorbar(csz)

# <codecell>

apr['utcz'] = apr['lonz']*0.0
for i in xrange(len(apr['utc'])):
    apr['utcz'][:,i] = apr['utc'][i]

# <codecell>

levels = np.arange(0,7000,30)

# <codecell>

plt.figure
csz = plt.contourf(apr['utcz'],apr['altflt'],apr['dbz'],levels,cmap=plt.cm.Greys)
plt.colorbar(csz)
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')

# <codecell>

it = np.abs(apr['utc']-19.1).argmin()

# <codecell>

plt.figure
for i in xrange(len()):

# <codecell>

plt.figure
plt.plot(apr['dbz'][:,it],apr['altflt'][:,it])
plt.xlim([0,8000])
plt.xlabel('dbZ')
plt.ylabel('Altitude [m]')
plt.title('APR2 Zenith radar reflectivity at: %f2.2',apr['utc'][it])

# <codecell>

plt.figure
plt.plot(apr['dbz'][:,0],apr['altflt'][:,0],'+')
plt.ylim([6000,8000])

# <codecell>

inoisezone = apr['dbz'][:,2000].argmax()
inoisezone

# <codecell>

apr['dbz'][:,]

# <codecell>

plt.figure
plt.plot(apr['dbz'][:,1000],apr['altflt'][:,1000],'+')
plt.ylim([6000,10000])

# <codecell>

apr['altflt'][:,0]

# <codecell>

plt.figure
plt.plot(apr['dbz'][11,:,0],apr['altz'][11,:,0]+apr['alt'][11,0])
plt.xlim([0,8000])
plt.ylim([0,20000])

# <codecell>

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

# <codecell>

it = np.abs(apr['utc']-19.19).argmin()
inoisezone = apr['dbz'][:,it].argmax()
i=range(inoisezone-15)+range(inoisezone+15,len(apr['altflt'][:,0]))

# <codecell>

print 0.19*60
print 0.10*60

# <codecell>

apr['utc'][it]

# <codecell>

plt.figure
plt.plot(smooth(apr['dbz'][i,it],10),apr['altflt'][i,it],'k')
plt.xlim([0,8000])
plt.xlabel('dbZ')
plt.ylabel('Altitude [m]')
#plt.title('APR2 Zenith radar reflectivity at: ',apr['utc'][it])

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

plt.figure
plt.plot(twoDS['effectiveD'][1:]/2.0,dc8['G_ALT'])
plt.xlim([0,100])
plt.xlabel('Effective radius [$\\mu$m]')
plt.ylabel('Altitude [m]')

# <codecell>

# filter the data and only show a part
twoDS['effectiveD'][twoDS['effectiveD']<0] = np.nan
fl = np.where((twoDS['Start_UTC'] > 18.025) & (twoDS['Start_UTC'] < 19.45) & (twoDS['effectiveD'] > 28.0))

# <codecell>

plt.figure
plt.plot(twoDS['effectiveD'][fl]/2.0,dc8['G_ALT'][fl],label='Cloud Probes (2DS)')
plt.xlim([0,100])
plt.xlabel('Effective radius [$\\mu$m]')
plt.ylabel('Altitude [m]')

# <codecell>

plt.figure
plt.plot(smooth(twoDS['effectiveD'][fl]/2.0,10),dc8['G_ALT'][fl],label='Cloud Probes (2DS)')
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

plt.plot(rsp['UTC'][rsp_good[0]],rsp['R_eff159'][rsp_good[0]],label='1.59 micron')
plt.plot(rsp['UTC'][rsp_good[0]],rsp['R_eff225'][rsp_good[0]],label='2.25 micron')
plt.xlabel('UTC [Hours]')
plt.ylabel('R$_{eff}$ [$\\mu$m]')
plt.title('RSP effective radius retrievals')
plt.legend(frameon=False)

# <codecell>

print rsp_good[0][-1]

# <codecell>

# get the distance between two points for RSP
from map_utils import spherical_dist
print rsp['Lat'][rsp_good[0][-1]], rsp['Lon'][rsp_good[0][-1]]
print rsp['Lat'][rsp_good[0][-2]], rsp['Lon'][rsp_good[0][-2]]
print spherical_dist(np.array((rsp['Lat'][rsp_good[0][-1]],rsp['Lon'][rsp_good[0][-1]])),
                     np.array((rsp['Lat'][rsp_good[0][-2]],rsp['Lon'][rsp_good[0][-2]])))

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
# Cloud base height: 13610 m
# Cloud top height: 16489 m
#  - Modis: 500 m
#  - eMAS: 50m or 2.5 mrad FOV ()
#  - RSP: 14 mrad FOV
#  - 4STAR: 2 degree, 0.0349 rads
#  - GOES: 1km
#  - SSFR: 180 degrees, Pi rads

# <codecell>

cld_base = 13610.0
cld_top = 16489.0
er2_hgt = 19018.0
dc8_hgt = 8047.0
r_eMAS = (er2_hgt - cld_top)*tan(0.0025)
r_RSP = (er2_hgt - cld_top)*tan(0.014)
r_SSFR = (er2_hgt - cld_top)*tan(pi/3)
r_4STAR = (cld_base - dc8_hgt)*tan(0.0349)

# <codecell>

print r_eMAS
print r_RSP
print r_SSFR
print r_4STAR

# <headingcell level=2>

# Plot histogram of different tau and ref comparison

# <codecell>

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

xr = ax.get_xlim()
yr = ax.get_ylim()
ax.add_patch(plt.Rectangle((0,0),0,0,color='none',edgecolor='g',hatch='x',linewidth=0.2,alpha=0.5,label='Horizontal variability'))
ax.set_xlim(xr)
ax.set_ylim(yr)

plt.legend(frameon=False)
plt.xlim([0,60])
plt.savefig(fp+'plots/hist_modis_4star_tau_v3_fill.png',dpi=600,transparent=True)
#plt.savefig(fp+'plots/hist_modis_4star_tau.pdf',bbox='tight')

# <codecell>

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

ax.add_patch(plt.Rectangle((0,0),0,0,color='none',edgecolor='g',hatch='x',linewidth=0.2,alpha=0.5,label='Horizontal variability'))

plt.legend(frameon=False,loc='upper right')
plt.xlim([10,80])
plt.ylim([0,0.3])
plt.savefig(fp+'plots/hist_modis_4star_ref_v3_fill.png',dpi=600,transparent=True)

# <codecell>

plt.figure()
plt.boxplot(smooth(modis_tau,6),vert=False,color='m')

# <codecell>

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

# <headingcell level=2>

# Combine the vertical information into one figure

# <codecell>

meas.good.shape

# <codecell>

it = np.abs(apr['utc']-19.19).argmin()
inoisezone = apr['dbz'][:,it].argmax()
i=range(inoisezone-15)+range(inoisezone+15,len(apr['altflt'][:,0]))

# <codecell>

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

# <codecell>

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

lut.sp_irrdn[1,:,0,4,0]

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

def interpirr_up(re,ta):
    "interpolate over irradiance to get the spectrum at ta and re"
    sp = np.zeros(lut.wvl.size)
    print re,ta
    for w in xrange(lut.wvl.size):
        irrupfx = interpolate.RectBivariateSpline(lut.refsp[lut.refsp>5],lut.tausp,lut.sp_irrup[1,w,0,lut.refsp>5,:],kx=1,ky=1)
        sp[w] = irrupfx(re,ta)
    return smooth(sp,20)

# <codecell>

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

# <codecell>

sp_modisp = interpirr(means.ref.modis+sttaumod,means.tau.modis+sttaumod)
sp_emasp = interpirr(means.ref.emas+sttauemas,means.tau.emas+sttauemas)
sp_ssfrp = interpirr(means.ref.ssfr+sttaussfr,means.tau.ssfr+sttaussfr)
sp_starp = interpirr(means.ref.star+sttaustar,means.tau.star+sttaustar)
sp_rspp = interpirr(means.ref.rsp+sttaursp,means.tau.rsp+sttaursp)

# <codecell>

sp_modisp_up = interpirr_up(means.ref.modis+sttaumod,means.tau.modis+sttaumod)
sp_emasp_up = interpirr_up(means.ref.emas+sttauemas,means.tau.emas+sttauemas)
sp_ssfrp_up = interpirr_up(means.ref.ssfr+sttaussfr,means.tau.ssfr+sttaussfr)
sp_starp_up = interpirr_up(means.ref.star+sttaustar,means.tau.star+sttaustar)
sp_rspp_up = interpirr_up(means.ref.rsp+sttaursp,means.tau.rsp+sttaursp)

# <codecell>

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
iutc189 = np.nanargmin(abs(ssfr_dc8['tmhrs']-18.9))
iutc191 = np.nanargmin(abs(ssfr_dc8['tmhrs']-19.1))
iutc192 = np.nanargmin(abs(ssfr_dc8['tmhrs']-19.2))
print ssfr_dc8['zspectra'].shape
dn = np.nanmean(ssfr_dc8['zspectra'][iutc189:iutc191,:],axis=0)
print dn.shape
up = np.nanmean(ssfr_dc8['nspectra'][iutc189:iutc191,:],axis=0)

# <codecell>

dnstd = np.nanstd(ssfr_dc8['zspectra'][iutc191:iutc192,:],axis=0)

# <codecell>

plt.plot(dnstd)

# <codecell>

dn_no = interpirr(5.0,0)

# <codecell>

plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,4,0])
plt.plot(lut.wvl,lut.sp_irrdn[1,:,0,4,1])
plt.plot(lut.wvl,dn_no)
plt.plot(ssfr_dc8['zenlambda'],dn)
plt.plot(ssfr_dc8['zenlambda'],up)
plt.plot((dn-up)-(dn_no-up))

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

fssps = interpolate.interp1d(ssfr_dc8['zenlambda'],(dn+dnstd)*r500,bounds_error=False)
ssps = fssps(lut.wvl)

# <codecell>

fspup = interpolate.interp1d(ssfr_dc8['zenlambda'],up*r500,bounds_error=False)
sspup = fspup(lut.wvl)

# <codecell>

plt.plot(lut.wvl,dn_no*sspup/ssp)

# <codecell>

plt.plot(lut.wvl[250:],ssp[250:])
plt.plot(lut.wvl[250:],ssps[250:])

# <codecell>

plt.plot(lut.wvl,sspup,'k')
plt.plot(lut.wvl,ssp,'r')
plt.plot(lut.wvl,dn_no,'b')
plt.plot(lut.wvl,(ssp-sspup)-(dn_no-sspup))
print 'SSFR forcing, measurement based: %f W/m^2' % np.nansum(((ssp[250:]-sspup[250:])-(dn_no[250:]-sspup[250:]))*lut.wvl[250:]/1000.0)
print 'SSFR forcing, std measurement based: %f W/m^2' % np.nansum(((ssps[250:]-sspup[250:])-(dn_no[250:]-sspup[250:]))*lut.wvl[250:]/1000.0)
                                                               

# <codecell>

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

# <codecell>

radeff_modis = np.sum(((sp_modis[250:]-sp_modis_up[250:])-(sp_modis_no[250:]-sp_modis_no_up[250:]))*lut.wvl[250:]/1000.0)
radeff_emas = np.sum(((sp_emas[250:]-sp_emas_up[250:])-(sp_emas_no[250:]-sp_emas_no_up[250:]))*lut.wvl[250:]/1000.0)
radeff_ssfr = np.sum(((sp_ssfr[250:]-sp_ssfr_up[250:])-(sp_ssfr_no[250:]-sp_ssfr_no_up[250:]))*lut.wvl[250:]/1000.0)
radeff_rsp = np.sum(((sp_rsp[250:]-sp_rsp_up[250:])-(sp_rsp_no[250:]-sp_rsp_no_up[250:]))*lut.wvl[250:]/1000.0)
radeff_star = np.sum(((sp_star[250:]-sp_star_up[250:])-(sp_star_no[250:]-sp_star_no_up[250:]))*lut.wvl[250:]/1000.0)

# <codecell>

radeff_modisp = np.sum(((sp_modisp-sp_modisp_up)-(sp_modis_no-sp_modis_no_up))*lut.wvl/1000.0)
radeff_emasp = np.sum(((sp_emasp-sp_emasp_up)-(sp_emas_no-sp_emas_no_up))*lut.wvl/1000.0)
radeff_ssfrp = np.sum(((sp_ssfrp-sp_ssfrp_up)-(sp_ssfr_no-sp_ssfr_no_up))*lut.wvl/1000.0)
radeff_rspp = np.sum(((sp_rspp-sp_rspp_up)-(sp_rsp_no-sp_rsp_no_up))*lut.wvl/1000.0)
radeff_starp = np.sum(((sp_starp-sp_starp_up)-(sp_star_no-sp_star_no_up))*lut.wvl/1000.0)

# <codecell>

print 'Modis: ',radeff_modis, ' +/- ', radeff_modisp-radeff_modis
print 'eMAS: ',radeff_emas, ' +/- ', radeff_emasp-radeff_emas
print 'SSFR: ',radeff_ssfr, ' +/- ', radeff_ssfrp-radeff_ssfr
print 'RSP: ',radeff_rsp, ' +/- ', radeff_rspp-radeff_rsp
print '4star: ', radeff_star, ' +/- ', radeff_starp-radeff_star

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

# <headingcell level=2>

# Get GOES - Problematic now...

# <codecell>

goes_file = fp+'er2/20130913/G13V04.0.CONUS.2013256.2045.PX.04K.NC'
goes,goes_dicts = load_hdf(goes_file)

# <codecell>

goes_values = (('lat',6),('lon',7),('tau',19),('ref',20))
goes,goes_dicts = load_hdf(goes_file,values=goes_values)

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

# <headingcell level=2>

# Plot out FOV of each instrument on eMAS figure

# <codecell>

from plotting_utils import circles
from map_utils import radius_m2deg, spherical_dist

# <codecell>

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

# <codecell>

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

# <codecell>

import map_utils
reload(map_utils)
from map_utils import stats_within_radius

# <codecell>

out_ssfrtau = stats_within_radius(er2['Latitude'][ier2[::10]],er2['Longitude'][ier2[::10]],emas_v1['lat'],emas_v1['lon'],emas_v1['tau'],r_SSFR)
out_ssfrref = stats_within_radius(er2['Latitude'][ier2[::10]],er2['Longitude'][ier2[::10]],emas_v1['lat'],emas_v1['lon'],emas_v1['ref'],r_SSFR)

# <codecell>

out_rspref = stats_within_radius(er2['Latitude'][ier2[::10]],er2['Longitude'][ier2[::10]],emas_v1['lat'],emas_v1['lon'],emas_v1['ref'],r_RSP)
out_rsptau = stats_within_radius(er2['Latitude'][ier2[::10]],er2['Longitude'][ier2[::10]],emas_v1['lat'],emas_v1['lon'],emas_v1['tau'],r_RSP)

# <codecell>

out_starref = stats_within_radius(dc8['G_LAT'][idc8],dc8['G_LONG'][idc8],emas_v1['lat'],emas_v1['lon'],emas_v1['ref'],r_4STAR)
out_startau = stats_within_radius(dc8['G_LAT'][idc8],dc8['G_LONG'][idc8],emas_v1['lat'],emas_v1['lon'],emas_v1['tau'],r_4STAR)

# <codecell>

r_MODIS = 500.0
out_modisref = stats_within_radius(dc8['G_LAT'][idc8],dc8['G_LONG'][idc8],emas_v1['lat'],emas_v1['lon'],emas_v1['ref'],r_MODIS)
out_modistau = stats_within_radius(dc8['G_LAT'][idc8],dc8['G_LONG'][idc8],emas_v1['lat'],emas_v1['lon'],emas_v1['tau'],r_MODIS)

# <codecell>

plt.plot(out_ssfrtau['std'],'g',label='SSFR')
plt.plot(out_rsptau['std'],'c',label='RSP')
plt.plot(out_modistau['std'],'m',label='MODIS')
plt.plot(out_startau['std'],'r',label='4STAR')
plt.xlabel('points')
plt.ylabel('$\\tau$')
plt.title('Standard deviation along horizontal area')
plt.legend(frameon=False)

# <codecell>

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

# <codecell>

print 'SSFR tau: %.2f, %.2f, max, min:%.2f, %.2f ' % (np.nanmean(out_ssfrtau['std']), np.median(out_ssfrtau['std']), np.nanmax(out_ssfrtau['std']), np.nanmin(out_ssfrtau['std']))
print 'RSP tau: %.2f, %.2f, max, min:%.2f, %.2f ' % (np.nanmean(out_rsptau['std']), np.median(out_rsptau['std']),np.nanmax(out_rsptau['std']), np.nanmin(out_rsptau['std']))
print 'MODIS tau: %.2f, %.2f, max, min:%.2f, %.2f ' % (np.nanmean(out_modistau['std']), np.median(out_modistau['std']),np.nanmax(out_modistau['std']), np.nanmin(out_modistau['std']))
print '4STAR tau: %.2f , %.2f, max, min:%.2f, %.2f' % (np.nanmean(out_startau['std']), np.median(out_startau['std']),np.nanmax(out_startau['std']), np.nanmin(out_startau['std']))

# <codecell>

print 'SSFR ref: %.2f, %.2f, max, min: %.2f,%.2f' % (np.nanmean(out_ssfrref['std']),np.median(out_ssfrref['std']),np.nanmax(out_ssfrref['std']),np.nanmin(out_ssfrref['std']))
print 'RSP ref: %.2f, %.2f, max, min: %.2f,%.2f' % (np.nanmean(out_rspref['std']),np.median(out_rspref['std']),np.nanmax(out_rspref['std']),np.nanmin(out_rspref['std']))
print 'MODIS ref: %.2f, %.2f, max, min: %.2f,%.2f' % (np.nanmean(out_modisref['std']),np.median(out_modisref['std']),np.nanmax(out_modisref['std']),np.nanmin(out_modisref['std']))
print '4STAR ref: %.2f, %.2f, max, min: %.2f,%.2f' % (np.nanmean(out_starref['std']),np.median(out_starref['std']),np.nanmax(out_starref['std']),np.nanmin(out_starref['std']))

# <codecell>

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

# <codecell>

plt.hist(out_ssfrtau['std'],bins=30, histtype='stepfilled', normed=True,color='g',label='SSFR',alpha=0.6)
plt.hist(out_rsptau['std'][~isnan(out_rsptau['std'])],bins=30, histtype='stepfilled', normed=True,color='c',label='RSP',alpha=0.6)
plt.hist(out_modistau['std'][~isnan(out_modistau['std'])],bins=30, histtype='stepfilled', normed=True,color='m',label='MODIS',alpha=0.6)
plt.hist(out_startau['std'][~isnan(out_startau['std'])],bins=30, histtype='stepfilled', normed=True,color='r',label='4STAR',alpha=0.6)
plt.xlim([0,3])
plt.legend(frameon=False)

# <codecell>


