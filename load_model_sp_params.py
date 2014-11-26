# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:
#   load_model_sp_params
# 
# Purpose:
#   Python script that is used to step through the building and use of the parameters from TCAP
#   Starts with the measured 4STAR zenith radiances, then loads the idl save file containing the modeled lut for that day
#   Regroups all the necessary steps to build all the figures used
# 
# Calling Sequence:
#   python load_model_sp_params.py
#   
# Input:
#   none at command line
#   
# Output:
#   figures and save files...
#   
# Keywords:
#   none
#   
# Dependencies:
#   see below imports
#   
# Needed Files:
#   - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#   - file.rc : for consistent creation of look of matplotlib figures
#   - sp_v1_20130219_4STAR.out : modeled spectra output for TCAP in idl save file
#   - 20130219starzen_rad.mat : special zenith radiance 4star matlab file 

# <codecell>

%config InlineBackend.rc = {}
import matplotlib
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
%matplotlib inline

# <codecell>

import numpy as np, h5py
import plotly.plotly as py
import scipy.io as sio
import math
import os
import warnings
warnings.simplefilter('ignore', np.RankWarning)
import Sp_parameters as Sp
py.sign_in("samuelleblanc", "4y3khh7ld4")
#import mpld3
#mpld3.enable_notbeook()

# <codecell>

# set the basic directory path
fp='C:\\Users\\sleblan2\\Research\\TCAP\\'

# <codecell>

# load the idl save file containing the modeled radiances
s=sio.idl.readsav(fp+'model/sp_v1_20130219_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp, wvl, z, re, ta)'

# <codecell>

# create custom key for sorting via wavelength
iwvls = np.argsort(s.zenlambda)
s.zenlambda.sort()

# <codecell>

# load the matlab file containing the measured TCAP radiances
m = sio.loadmat(fp+'4STAR/20130219starzen_rad.mat')
sm = sio.idl.AttrDict(m)
print sm.keys()
print 'Measured radiance Shape: ', sm.rad.shape

print np.nanmax(sm.rad[sm.good[100],:])
sm.good[100]

# <headingcell level=4>

# Next section loads a few functions that can be used for typical analysis

# <codecell>

from Sp_parameters import nanmasked, closestindex, norm2max
    
time_ref=17.22
ii = closestindex(sm.utc,time_ref)
rad,mask = nanmasked(sm.rad[sm.good[ii],:])

# <headingcell level=4>

# Plotting functions defined

# <codecell>

# set up plotting of a few of the zenith radiance spectra
def pltzen(fig=None,ax=None, tit='Zenith spectra'):
    "Plotting of zenith measurements in radiance units"
    if ax is None: 
        fig,ax = plt.subplots()
        doaxes = True
    else:
        doaxes = False
    ax.plot(sm.nm[mask],rad,lw=2, c='k', label='4STAR measured at: '+str(time_ref))
    if doaxes:
        plt.title(tit)
        plt.ylabel('Radiance [Wm$^{-2}$nm$^{-1}$sr$^{-1}$]')
        plt.xlabel('Wavelength [nm]')
        plt.xlim([350,1700])
        plt.ylim([0,0.22])
        plt.legend(frameon=False)
    #plot_url = py.plot_mpl(fig)
    return fig,ax

def norm(fig=None,ax=None):
    "Plotting of zenith measurements in normalized radiance"
    if ax is None:
        fig,ax = plt.subplots()
        doaxes = True
    else:
        doaxes = False
    ax.plot(sm.nm[mask],norm2max(rad),lw=2, c='k', label='4STAR measured at: '+str(time_ref))
    if doaxes:
        plt.title('Zenith spectra')
        plt.ylabel('Normalized Radiance')
        plt.xlabel('Wavelength [nm]')
        plt.xlim([350,1700])
        plt.ylim([0,1.0])
        plt.legend(frameon=False)
    #plot_url = py.plot_mpl(fig)
    return fig,ax

def dashlen(dashlength,dashseperation,fig=plt.gcf()):
    """ Build a list of dash length that fits within the current figure or figure denoted by fig, 
        each dash is length dashlength, with its centers at dashseperation """
    totallen = fig.get_figwidth()
    numdash = int(totallen/dashseperation)*2
    f=lambda i: dashlength if i%2==0 else dashseperation-dashlength
    return tuple([f(i) for i in range(numdash)])

def plot_line_gradients(ax,s,names,cmap,iphase,irefs,itau,iwvls,pos,normalize=False):
    """ Make multiple lines on the subplot ax of the spectra s, for the case defined by names with the cmap
      for one particular phase (iphase), range of refs (irefs) and at one itau. Returns the axis handles for the thin and thick ref """
    rf = range(irefs[0],irefs[1])
    colors = plt.cm._generate_cmap(cmap,int(len(rf)*2.25))
    for ir in rf:
        if not(normalize):
            a1 = ax.plot(s.zenlambda,s.sp[iphase,iwvls,0,ir,itau],
                         color=(0.2,0.2,0.2),
                         lw=1.8*ir/irefs[1])
            ax.plot(s.zenlambda,s.sp[iphase,iwvls,0,ir,itau],
                     color=colors(ir),
                     lw=1.7*ir/irefs[1])
            ax.text(pos[0],pos[1],names,color=colors(irefs[1]))
        else:
            a1 = ax.plot(s.zenlambda,norm2max(s.sp[iphase,iwvls,0,ir,itau]),
                         color=(0.2,0.2,0.2),
                         lw=1.8*ir/irefs[1])
            ax.plot(s.zenlambda,norm2max(s.sp[iphase,iwvls,0,ir,itau]),
                     color=colors(ir),
                     lw=1.7*ir/irefs[1])    
            ax.text(pos[0],pos[1]/0.22,names,color=colors(irefs[1]))
        if ir == rf[0]:
            alow = a1
        if ir == rf[-1]:
            ahigh = a1
    return [alow,ahigh]

def plot_greys(fig=None,ax=None):
    " Plotting of grey regions that indicates the different wavelenght regions where the parameters are defined. "
    cl = '#CCCCCC'
    plt.axvspan(1000,1077,color=cl) #eta1
    plt.axvspan(1192,1194,color=cl) #eta2
    plt.axvspan(1492,1494,color=cl) #eta3
    plt.axvspan(1197,1199,color=cl); plt.axvspan(1235,1237,color=cl);  #eta4
    plt.axvspan(1248,1270,color=cl) #eta5
    plt.axvspan(1565,1644,color=cl) #eta6
    plt.axvspan(1000,1050,color=cl) #eta7
    plt.axvspan(1493,1600,color=cl) #eta8
    plt.axvspan(1000,1077,color=cl) #eta9
    plt.axvspan(1200,1300,color=cl) #eta10
    plt.axvspan(530 ,610 ,color=cl) #eta11
    plt.axvspan(1039,1041,color=cl) #eta12
    plt.axvspan(999 ,1001,color=cl); plt.axvspan(1064,1066,color=cl);  #eta13
    plt.axvspan(599 ,601 ,color=cl); plt.axvspan(869 ,871 ,color=cl);  #eta14
    plt.axvspan(1565,1634,color=cl); #eta15
    

# <headingcell level=3>

# Plotting iterations

# <codecell>

fig,ax=pltzen()

# <markdowncell>

# Next figure with modeled spectra

# <codecell>

# now go through and add the different modeled spectra
fig,ax=pltzen()

lines = [('Liquid Cloud Model, COD=0.5','Reds',0,[0,13],1,[420,0.01]),
         ('Ice Cloud Model, COD=0.5','Greens',1,[13,34],1,[380,0.02]),
         ('Liquid Cloud Model, COD=10','Purples',0,[0,13],9,[700,0.16]),
         ('Ice Cloud Model, COD=10','Blues',1,[13,34],9,[750,0.15])]

for names,cmap,iphase,irefs,itau,pos in lines:
    [alow,ahigh] = plot_line_gradients(ax,s,names,cmap,iphase,irefs,itau,iwvls,pos)
    
lbl=["Small R$_{eff}$ (Ice=" + str(s.ref[34]) + " $\mu m$, Liquid=" + str(s.ref[13]) + " $\mu m$)",
     "Large R$_{eff}$ (Ice=" + str(s.ref[13]) + " $\mu m$, Liquid=" + str(s.ref[0]) + " $\mu m$)"]
plt.legend([alow[0],ahigh[0]],
           lbl,
           frameon=False,
           loc=7)
ax.text(600,0.19,'4STAR Measurement')
pltzen(fig,ax)

# <headingcell level=3>

# Next figure with normalized spectra and areas of parameters

# <codecell>

fig,ax=norm()
for names,cmap,iphase,irefs,itau,pos in lines:
    [alow,ahigh] = plot_line_gradients(ax,s,names,cmap,iphase,irefs,itau,iwvls,pos,normalize=True)
plt.legend([alow[0],ahigh[0]],
           lbl,
           frameon=False,
           loc=7)
ax.text(600,0.19/0.22,'4STAR Measurement')
norm(fig,ax)
plot_greys()

# <headingcell level=3>

# Now calculate the parameters for the measured spectra

# <codecell>

# first convert measurements to Sp class, with inherent parameters defined
meas = Sp.Sp(m)
meas.params()

# <markdowncell>

# Plot the parameters for the specified time

# <codecell>

fig2,ax2 = plt.subplots(5,3,sharex=True,figsize=(15,8))
ax2 = ax2.ravel()
for i in range(meas.npar-1):
    ax2[i].plot(meas.utc,Sp.smooth(meas.par[:,i],3))
    ax2[i].set_title('Parameter '+str(i))
    ax2[i].grid()
    ax2[i].set_xlim([17,19])
    if i > 11: 
        ax2[i].set_xlabel('UTC [h]')

fig2.tight_layout()
plt.show()

# <headingcell level=3>

# Prepare the LUT for the modeled spectra

# <codecell>

reload(Sp)

# <codecell>

lut = Sp.Sp(s)
warnings.simplefilter('ignore')
lut.params()

# <codecell>

lut.sp_hires()
lut.params()

# <codecell>

print lut.par.shape

# <markdowncell>

# Now plot the resulting lut of parameters

# <codecell>

fig3,ax3 = plt.subplots(5,3,sharex=True,figsize=(15,8))
ax3 = ax3.ravel()
for i in range(lut.npar-1):
    for j in xrange(len(lut.ref)):
        ax3[i].plot(lut.tau,lut.par[0,j,:,i])
    ax3[i].set_title('Parameter '+str(i))
    ax3[i].grid()
    ax3[i].set_xlim([0,100])
    if i > 11: 
        ax3[i].set_xlabel('Tau')

fig3.tight_layout()
plt.show()

# <codecell>

print lut.sp.shape
print lut.tau[80]

# <codecell>

plt.figure()
for i in xrange(30):
    plt.plot(lut.wvl,lut.sp[0,:,0,i,80])

# <codecell>

all_zeros = not np.any(lut.sp[0,:,0,40,39])
print all_zeros
print np.max(lut.sp[0,:,0,20,80]), np.min(lut.sp[0,:,0,20,80])
print np.any(lut.sp[0,:,0,20,70])
print lut.ref[28]
plt.figure()
for i in xrange(65,75):
    plt.plot(lut.wvl,lut.sp[0,:,0,28,i])

# <codecell>

print not np.any(lut.sp[0,:,0,28,60])
print lut.ref[3]

# <codecell>

import run_kisq_retrieval as rk
reload(rk)
import Sp_parameters as Sp
reload(Sp)
#del lut
#del meas

# <codecell>

 (tau,ref,phase,ki) = rk.run_retrieval(meas,lut)

# <codecell>

print meas.good
print meas.utc.shape

# <codecell>

plt.figure()
plt.plot(meas.utc,tau)

# <codecell>


# <rawcell>

#     num_noise = 200

