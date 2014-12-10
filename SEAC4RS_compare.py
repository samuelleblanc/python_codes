# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     SEAC4RS_compare
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
# 
# Calling Sequence:
# 
#     python SEAC4RS_compare.py
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
#   - sp_v2_20130913_4STAR.out: modeled spectra output for SEAC4RS at sza 17.9, idl save file
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
fp='C:\\Users\\sleblan2\\Research\\SEAC4RS\\'

# <headingcell level=2>

# Get the lookup table for the 4STAR data

# <codecell>

# load the idl save file containing the modeled radiances
s=sio.idl.readsav(fp+'model/sp_v2_20130913_4STAR.out')#fp+'model/sp_v1.1_20130913_4STAR.out')
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
from load_modis import mat2py_time, toutc

# <codecell>

# load the matlab file containing the measured TCAP radiances
m = sio.loadmat(fp+'../4STAR/SEAC4RS/20130913/20130913starzen_3.mat')
m.keys()

# <markdowncell>

# Go through and get the radiances for good points, and for the time selected

# <codecell>

print m['t']
tt = mat2py_time(m['t'])
m['utc'] = toutc(tt)

# <codecell>

m['good'] = np.where((m['utc']>18.5) & (m['utc']<19.75) & (m['Str'].flatten()!=0) & (m['sat_time'].flatten()==0))

# <codecell>

plt.figure()
plt.plot(m['utc'],abs(m['rad'][:,1024]-m['rad'][:,1068]))
plt.xlabel('UTC [Hours]')
plt.ylabel('difference in radiance at 1024 nm')

# <codecell>

plt.figure()
plt.plot(m['utc'],m['rad'][:,400]/1000.)
plt.xlabel('UTC [hours]')
plt.ylabel('Radiance at 400 nm [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')

# <codecell>

reload(Sp)
if 'meas' in locals():
    del meas
    import gc; gc.collect()
# first convert measurements to Sp class, with inherent parameters defined
meas = Sp.Sp(m)
meas.params()

# <headingcell level=2>

# Run the retrieval on 4STAR data

# <codecell>

import run_kisq_retrieval as rk
reload(rk)

# <codecell>

print max(meas.good)
print type(m['good'])
print isinstance(m['good'],tuple)
print type(m['good'][0])

# <codecell>

(meas.tau,meas.ref,meas.phase,meas.ki) = rk.run_retrieval(meas,lut)

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
plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results.png',dpi=600)
plt.savefig(fp+'plots\\SEAC4RS_20130913_retri_results.pdf',bbox='tight')

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

ssfr_idl_file = fp+'er2/20130913/20130913_calibspcs.out'
ssfr_idl = sio.idl.readsav(ssfr_idl_file)
print ssfr_idl.keys()
print np.shape(ssfr_idl['zspectra'])
print np.shape(ssfr_idl['sat'])

# <codecell>

class object():
    pass
ssfr = object()
if True:
    ssfr.utc = ssfr_er2['UTC']
    ssfr.Rvis = ssfr_er2['UP500']/ssfr_er2['DN500']
    ssfr.Rnir = ssfr_er2['UP1600']/ssfr_er2['DN1600']
    ssfr.Rvis[ssfr_er2['UP500']==-999] = np.nan
    ssfr.Rnir[ssfr_er2['UP500']==-999] = np.nan
    ssfr.Rnir[ssfr_er2['UP500']<=0] = np.nan
    ssfr.Rvis[ssfr_er2['UP500']<=0] = np.nan
    print np.max(er2['Start_UTC'])
    print np.max(ssfr.utc)

# <codecell>

if False:
    ssfr.utc = ssfr_idl['tmhrs']
    i500 = np.argmin(abs(ssfr_idl['zenlambda']-500.0))
    i1650 = np.argmin(abs(ssfr_idl['zenlambda']-1650.0))
    ssfr.Rvis = ssfr_idl['zspectra'][:,i500]/ssfr_idl['nspectra'][:,i500]
    ssfr.Rnir = ssfr_idl['zspectra'][:,i1650]/ssfr_idl['nspectra'][:,i1650-1]
    ssfr.Rvis[ssfr.Rvis<0] = np.nan
    ssfr.Rvis[ssfr.Rvis>1] = np.nan
    ssfr.Rnir[ssfr.Rnir<0] = np.nan
    ssfr.Rnir[ssfr.Rnir>1] = np.nan

# <codecell>

print np.nanmin(ssfr.Rvis), np.nanmax(ssfr.Rvis)

# <codecell>

print ssfr_idl['nadlambda'][i1650-2]
print i500, i1650

# <codecell>

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
fig.add_subplot(1,2,2)
plt.plot(ssfr.utc,ssfr.Rvis,'b',ssfr.utc,ssfr.Rnir,'r')

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

# Now run the retrieval

# <codecell>

import run_2wvl_retrieval as rw
if 'rw' in locals():
    reload(rw)

# <codecell>

(ssfr.tau,ssfr.ref,ssfr.ki) = rw.run_2wvl_retrieval(ssfr,lut,wvls=[500.0,1650.0])

# <codecell>

plt.figure()
plt.plot(ssfr.tau[:,1])
print ssfr.tau.shape

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
            llcrnrlon=-97, llcrnrlat=18,
            urcrnrlon=-89, urcrnrlat=23,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(-89,-97,9),labels=[0,0,0,1])
    m.drawparallels(np.linspace(18,23,5),labels=[1,0,0,0])
    return m

# <codecell>

figm = plt.figure()
axm = figm.add_axes([0.1,0.1,0.8,0.8])
ma = seac_map(axm)
x,y = ma(modis['lon'],modis['lat'])
clevels = np.linspace(0,80,40)
cs = ma.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.gist_ncar)
cbar = ma.colorbar(cs)
cbar.set_label('$\\tau$')
axm.set_title('MODIS - AQUA Cloud optical Thickness')

# <codecell>

figm2,axm2 = plt.subplots(1,2,figsize=(13,13))
m1 = seac_map(axm2[0])
m2 = seac_map(axm2[1])
x,y = m1(modis['lon'],modis['lat'])

cs1 = m1.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.gist_ncar)
cbar = m1.colorbar(cs1)
cbar.set_label('$\\tau$')
axm2[0].set_title('MODIS - AQUA Cloud optical Thickness')

clevels2 = np.linspace(0,60,30)
cs2 = m2.contourf(x,y,modis['ref'],clevels2,cmap=plt.cm.gist_earth)
cbar = m2.colorbar(cs2)
cbar.set_label('R$_{ef}$ [$\\mu$m]')
axm2[1].set_title('MODIS - AQUA Cloud effective radius')
figm2.subplots_adjust(wspace=0.3)
plt.show()
plt.savefig(fp+'plots/modis_only_tau_ref_comp.png',dpi=600)
plt.savefig(fp+'plots/modis_only_tau_ref_comp.pdf',bbox='tight')

# <codecell>

figm2,axm2 = plt.subplots(1,2,figsize=(13,13))
m1 = seac_map(axm2[0])
m2 = seac_map(axm2[1])
x,y = m1(modis['lon'],modis['lat'])

cs1 = m1.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.gist_ncar)
cbar = m1.colorbar(cs1)
cbar.set_label('$\\tau$')
axm2[0].set_title('MODIS - AQUA Cloud optical Thickness')
x1,y1 = m1(meas.lon,meas.lat)
m1.scatter(x1,y1,c=meas.tau,cmap=plt.cm.gist_ncar,marker='o',vmin=clevels[0],vmax=clevels[-1],alpha=0.5,edgecolors='k',linewidth=0.15)

cs2 = m2.contourf(x,y,modis['ref'],clevels2,cmap=plt.cm.gist_earth)
cbar = m2.colorbar(cs2)
cbar.set_label('R$_{ef}$ [$\\mu$m]')
axm2[1].set_title('MODIS - AQUA Cloud effective radius')
m2.scatter(x1,y1,c=meas.ref,cmap=plt.cm.gist_earth,marker='o',vmin=clevels2[0],vmax=clevels2[-1],alpha=0.5,edgecolors='k',linewidth=0.15)
figm2.subplots_adjust(wspace=0.3)
plt.show()
plt.savefig(fp+'plots/modis_dc8_tau_ref_comp.png',dpi=600)
plt.savefig(fp+'plots/modis_dc8_tau_ref_comp.pdf',bbox='tight')

# <headingcell level=2>

# Import eMAS values

# <codecell>

if 'lm' in locals():
    reload(lm)
from load_modis import load_emas

# <codecell>

emas_file = fp+'er2/20130913/EMASL2_13965_11_20130913_1832_1845_V00.hdf'
print os.path.isfile(emas_file)

# <codecell>

emas,emad_dicts = load_emas(emas_file)

# <codecell>


