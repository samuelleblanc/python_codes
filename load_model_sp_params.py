
# coding: utf-8

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

# In[18]:

import matplotlib
get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np, h5py
import plotly.plotly as py
import scipy.io as sio
import math
import os
import warnings
warnings.simplefilter('ignore', np.RankWarning)
py.sign_in("samuelleblanc", "4y3khh7ld4")
print 'C:\\Users\\sleblan2\\Research\\python_codes\\file.rc'
from matplotlib import rc_file
rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')

#import mpld3
#mpld3.enable_notbeook()


# In[3]:

# set the basic directory path
fp='C:\\Users\\sleblan2\\Research\\TCAP\\'


# In[4]:

# load the idl save file containing the modeled radiances
s=sio.idl.readsav(fp+'model/sp_v1_20130219_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp, wvl, z, re, ta)'


# In[5]:

# create custom key for sorting via wavelength
iwvls = np.argsort(s.zenlambda)
s.zenlambda.sort()


# In[6]:

# load the matlab file containing the measured TCAP radiances
m = sio.loadmat(fp+'4STAR/20130219starzen_rad.mat')
sm = sio.idl.AttrDict(m)
print sm.keys()
print 'Measured radiance Shape: ', sm.rad.shape

print np.nanmax(sm.rad[sm.good[100],:])
sm.good[100]


##### Next section loads a few functions that can be used for typical analysis

# In[7]:

def nanmasked(x):
    " Build an array with nans masked out and the mask output"
    mask = ~np.isnan(x)
    maskA = x[mask]
    return (maskA,mask)

def closestindex(a,x):
    " Get the index from a of the closest value from x "
    return min(range(len(a)), key=lambda i: abs(a[i]-x))

def norm2max(x):
    " Returns a spectrum, x, that is normalized by its maximum value, ignores nans "
    return x/np.nanmax(x)
    
time_ref=17.22
ii = closestindex(sm.utc,time_ref)
rad,mask = nanmasked(sm.rad[sm.good[ii],:])


##### Plotting functions defined

# In[8]:

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
    


#### Plotting iterations

# In[19]:

fig,ax=pltzen()


# Next figure with modeled spectra

# In[20]:

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


#### Next figure with normalized spectra and areas of parameters

# In[21]:

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


#### Now calculate the parameters for the measured spectra

# In[22]:

from Sp_parameters import *
# first convert measurements to Sp class, with inherent parameters defined
meas = Sp(m)
meas.params()


# Plot the parameters for the specified time

# In[23]:

fig2,ax2 = plt.subplots(5,3,sharex=True,figsize=(15,8))
ax2 = ax2.ravel()
for i in range(meas.npar-1):
    ax2[i].plot(meas.utc,smooth(meas.par[:,i],3))
    ax2[i].set_title('Parameter '+str(i))
    ax2[i].grid()
    ax2[i].set_xlim([17,19])
    if i > 11: 
        ax2[i].set_xlabel('UTC [h]')

fig2.tight_layout()
plt.show()


#### Prepare the LUT for the modeled spectra

# In[24]:

lut = Sp(s)
lut.params()


# In[184]:

lut.tau
lut.ref


# In[198]:

r=np.array([[5,10],[8,13],[10,15]])


# In[199]:

print r


# In[230]:

x = np.array([0,3,5])
y = np.array([2,4])
yi = r*0.+y
print yi
xi = np.transpose(r.T*0.+x.T)
print xi


# In[231]:

xx = [0,0.5,0.8,1,2,3,4,5]
yy = [1,2,3,4]


# In[232]:

import scipy.interpolate as interp


# In[234]:

rr = interp.interp2d(xi,yi,r,bounds_error=False,fill_value=None)


# In[272]:

print lut.sp[0,0,0,:,:].shape,lut.ref.shape
ttau = lut.sp[0,0,0,:,:]*0.+lut.tau
rref = np.transpose(lut.sp[0,0,0,:,:].T*0.+lut.ref.T)
print ttau


# In[290]:

print lut.par.shape
po=np.arange(0,16)
print po.shape
dd = np.transpose(lut.par[0,:,0,:,:].T*po.T-po.T)
print dd.shape


# In[289]:

dd[1,20,12]


                
                