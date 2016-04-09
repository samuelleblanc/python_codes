
# coding: utf-8

# Name:  
# 
#     load_model_sp_params
# 
# Purpose:  
# 
#     Python script that is used to step through the building and use of the parameters from TCAP
#     Starts with the measured 4STAR zenith radiances, then loads the idl save file containing the modeled lut for that day
#     Regroups all the necessary steps to build all the figures used to analyze the data
#     Runs the ki^2 retrieval with 15 parameters
#     plots the results
#     Compares to MODIS retrievals for that day
# 
# Calling Sequence:
# 
#     python load_model_sp_params.py
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
#   - sp_v1_20130219_4STAR.out : modeled spectra output for TCAP in idl save file
#   - 20130219starzen_rad.mat : special zenith radiance 4star matlab file 

# In[ ]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
get_ipython().magic(u'matplotlib notebook')
import numpy as np, h5py
#import plotly.plotly as py
import scipy.io as sio
import math
import os
import Sp_parameters as Sp
#py.sign_in("samuelleblanc", "4y3khh7ld4")
#import mpld3
#mpld3.enable_notbeook()


# In[ ]:

import IPython
IPython.InteractiveShell.cache_size = 0


# In[ ]:

# set the basic directory path
fp='C:\\Users\\sleblan2\\Research\\TCAP\\'


# In[4]:

# load the idl save file containing the modeled radiances
s=sio.idl.readsav(fp+'model/sp_v1_20130219_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp, wvl, z, re, ta)'
# create custom key for sorting via wavelength
iwvls = np.argsort(s.zenlambda)
s.wv = np.sort(s.zenlambda)


# In[5]:

# load the matlab file containing the measured TCAP radiances
m = sio.loadmat(fp+'4STAR/20130219starzen_rad.mat')
sm = sio.idl.AttrDict(m)
print sm.keys()
print 'Measured radiance Shape: ', sm.rad.shape

print np.nanmax(sm.rad[sm.good[100],:])
sm.good[100]
print sm.good.shape, sm.good.max()


# #### Next section loads a few functions that can be used for typical analysis

# In[6]:

from Sp_parameters import nanmasked, closestindex, norm2max
    
time_ref=17.22
ii = closestindex(sm.utc,time_ref)
rad,mask = nanmasked(sm.rad[sm.good[ii],:])


# #### Plotting functions defined

# In[9]:

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

def norm(fig=None,ax=None,dolegend=True):
    "Plotting of zenith measurements in normalized radiance"
    if ax is None:
        fig,ax = plt.subplots()
        doaxes = True
    else:
        doaxes = False
    ax.plot(sm.nm[mask],norm2max(rad),lw=2, c='k', label='4STAR measured at: '+str(time_ref))
    if doaxes:
        plt.title('Zenith radiance spectra')
        plt.ylabel('Normalized Radiance')
        plt.xlabel('Wavelength [nm]')
        plt.xlim([350,1700])
        plt.ylim([0,1.0])
        if dolegend:
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
            a1 = ax.plot(s.wv,s.sp[iphase,iwvls,0,ir,itau],
                         color=(0.2,0.2,0.2),
                         lw=1.0+1.4*float(ir)/irefs[1])
            ax.plot(s.wv,s.sp[iphase,iwvls,0,ir,itau],
                     color=colors(ir),
                     lw=1.0+1.3*float(ir)/irefs[1])
            ax.text(pos[0],pos[1],names,color=colors(irefs[1]))
        else:
            a1 = ax.plot(s.wv,norm2max(s.sp[iphase,iwvls,0,ir,itau]),
                         color=(0.2,0.2,0.2),
                         lw=1.0+1.4*float(ir)/irefs[1])
            ax.plot(s.wv,norm2max(s.sp[iphase,iwvls,0,ir,itau]),
                     color=colors(ir),
                     lw=1.0+1.3*float(ir)/irefs[1])    
            ax.text(pos[0],pos[1]/0.22,names,color=colors(irefs[1]))
        if ir == rf[0]:
            alow = a1
        if ir == rf[-1]:
            ahigh = a1
    return [alow,ahigh]

def plot_greys(fig=None,ax=None):
    " Plotting of grey regions that indicates the different wavelenght regions where the parameters are defined. "
    cl = '#DDDDDD'
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
    


# ### Plotting iterations

# In[9]:

fig,ax=pltzen()


# Next figure with modeled spectra

# In[10]:

# now go through and add the different modeled spectra
fig,ax=pltzen()

lines = [('Liquid Cloud Model, COD=0.5','Reds',0,[0,13],1,[420,0.01]),
         ('Ice Cloud Model, COD=0.5','Greens',1,[13,34],1,[380,0.02]),
         ('Liquid Cloud Model, COD=10','RdPu',0,[0,13],9,[700,0.16]),
         ('Ice Cloud Model, COD=10','Blues',1,[13,34],9,[750,0.15])]

for names,cmap,iphase,irefs,itau,pos in lines:
    [alow,ahigh] = plot_line_gradients(ax,s,names,cmap,iphase,irefs,itau,iwvls,pos)
    
lbl=["Small R$_{eff}$ (Ice=" + str(s.ref[13]) + " $\mu m$, Liquid=" + str(s.ref[0]) + " $\mu m$)",
     "Large R$_{eff}$ (Ice=" + str(s.ref[34]) + " $\mu m$, Liquid=" + str(s.ref[13]) + " $\mu m$)"]
plt.legend([alow[0],ahigh[0]],
           lbl,
           frameon=False,loc=7,prop={'size':10})
ax.text(600,0.19,'4STAR Measurement from TCAP')
pltzen(fig,ax)


# ### Next figure with normalized spectra and areas of parameters

# In[13]:

fig,ax = norm(dolegend=False)
plt.axvspan(350,1700,color='#FFFFFF')
ax.text(600,0.19/0.22,'4STAR Measurement')
plt.savefig(fp+'plots/zen_spectra_nomodel.png',dpi=600,transparent=True)


# In[11]:

fig,ax=norm()
lines = [('Liquid Cloud Model, $\\tau$=0.5','Reds',0,[0,13],1,[420,0.01]),
         ('Ice Cloud Model, $\\tau$=0.5','Greens',1,[13,34],1,[380,0.02]),
         ('Liquid Cloud Model, $\\tau$=10','RdPu',0,[0,13],9,[700,0.16]),
         ('Ice Cloud Model, $\\tau$=10','Blues',1,[13,34],9,[750,0.15])]
for names,cmap,iphase,irefs,itau,pos in lines:
    [alow,ahigh] = plot_line_gradients(ax,s,names,cmap,iphase,irefs,itau,iwvls,pos,normalize=True)
plt.legend([alow[0],ahigh[0]],lbl,
           frameon=False,loc=7,prop={'size':10})
ax.text(600,0.19/0.22,'4STAR Measurement')
norm(fig,ax)
plt.axvspan(350,1700,color='#FFFFFF')
plot_greys()
plt.savefig(fp+'plots/zen_spectra_model.png',dpi=600,transparent=True)
#plt.savefig(fp+'plots/zen_spectra_model.eps')


# Plot the same but with squeezed axes

# In[17]:

fp


# In[16]:

fig,ax=norm()
fig.set_figheight(5)
fig.set_figwidth(6)
lines = [('Liquid Cloud Model, $\\tau$=0.5','Reds',0,[0,13],1,[420,0.01]),
         ('Ice Cloud Model, $\\tau$=0.5','Greens',1,[13,34],1,[380,0.02]),
         ('Liquid Cloud Model, $\\tau$=10','RdPu',0,[0,13],9,[700,0.16]),
         ('Ice Cloud Model, $\\tau$=10','Blues',1,[13,34],9,[750,0.15])]
for names,cmap,iphase,irefs,itau,pos in lines:
    [alow,ahigh] = plot_line_gradients(ax,s,names,cmap,iphase,irefs,itau,iwvls,pos,normalize=True)
plt.legend([alow[0],ahigh[0]],lbl,
           frameon=False,loc=7,prop={'size':10})
ax.text(600,0.19/0.22,'4STAR Measurement')
norm(fig,ax)
plt.axvspan(350,1700,color='#FFFFFF')
plot_greys()
plt.savefig(fp+'plots/zen_spectra_model_squeeze.png',dpi=600,transparent=True)


# ### Now calculate the parameters for the measured spectra

# In[12]:

map(lambda x:x*x,[-1,1,24])


# In[13]:

reload(Sp)
if 'meas' in locals():
    del meas
    import gc; gc.collect()


# In[14]:

# first convert measurements to Sp class, with inherent parameters defined
meas = Sp.Sp(m)
meas.params()


# Plot the parameters for the specified time

# In[15]:

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


# ### Prepare the LUT for the modeled spectra

# In[ ]:

reload(Sp)
if 'lut' in locals():
    del lut
    import gc; gc.collect()


# In[16]:

lut = Sp.Sp(s)
lut.params()
lut.param_hires()


# In[17]:

import gc; gc.collect()
import sys
print sys.version


# In[18]:

lut.sp_hires()


# In[19]:

print lut.tau.shape
print lut.ref.shape
print lut.sp.ndim
print lut.par.size
print lut.par.shape


# In[20]:

print lut.ref
print lut.sp[0,400,0,23,10]
print lut.sp[1,400,0,:,10]
print lut.par.shape


# Now plot the resulting lut of parameters

# In[21]:

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


# In[22]:

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


# In[23]:

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


# In[24]:

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


# Now run through a few spectra for double checking the input to make sure everythin matches

# In[25]:

print lut.sp.shape
print lut.tau[80]


# In[26]:

print meas.par.shape


# In[27]:

print meas.par[meas.good[200],13]


# In[28]:

plt.figure()
color.cycle_cmap(31,cmap=plt.cm.gist_ncar)
for i in xrange(30):
    plt.plot(lut.wvl,lut.sp[0,:,0,i,80])
plt.title('Liquid water cloud with varying R$_{ef}$ at $\\tau$ of '+str(lut.tau[80]))
plt.xlabel('Wavelength [nm]')
plt.ylabel('Radiance [Wm$^{-2}$sr$^{-1}$nm$^{-1}$]')
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar
                                  ,norm=plt.Normalize(vmin=0,vmax=1))
scalarmap.set_array(lut.ref[range(30)])
cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6))
cba.ax.set_ylabel('R$_{ef}$ [$\\mu$m]')
cba.ax.set_yticklabels(np.linspace(lut.ref[0],lut.ref[30],6));


# In[29]:

all_zeros = not np.any(lut.sp[0,:,0,40,39])
print all_zeros
print np.max(lut.sp[0,:,0,20,80]), np.min(lut.sp[0,:,0,20,80])
print np.any(lut.sp[0,:,0,20,70])
print lut.ref[28]
plt.figure()
for i in xrange(75,85):
    plt.plot(lut.wvl,lut.sp[0,:,0,28,i], label="tau:"+str(lut.tau[i]))
plt.legend()
plt.title('Liquid water cloud with R$_{ef}$ = '+str(lut.ref[28])+' at varying $\\tau$ between '+str(lut.tau[75])+' to '+str(lut.tau[85]))
plt.xlabel('Wavelength [nm]')
plt.ylabel('Radiance [Wm$^{-2}$sr$^{-1}$nm$^{-1}$]')


# In[30]:

plt.figure()
plt.plot(lut.wvl,lut.sp[0,:,0,28,70])
plt.title('Liquid water cloud with R$_{ef}$ = '+str(lut.ref[28])+' at $\\tau$ '+str(lut.tau[70]))
plt.xlabel('Wavelength [nm]')
plt.ylabel('Radiance [Wm$^{-2}$sr$^{-1}$nm$^{-1}$]')


# In[31]:

print np.any(lut.sp[0,:,0,28,70])
print lut.ref[3]


# In[32]:

ro = (range(1,20),range(4,50))
print ro[1]
print ro[0]


# In[33]:

from scipy import interpolate
print np.shape([lut.tau[69],lut.tau[71]])
print np.shape([lut.sp[0,:,0,28,69],lut.sp[0,:,0,28,69]])
fs = interpolate.interp1d([lut.tau[69],lut.tau[71]],[lut.sp[0,:,0,28,69],lut.sp[0,:,0,28,69]],axis=0)
sss = fs(lut.tau[70])


# In[34]:

print np.shape(sss)
plt.figure()
plt.plot(lut.wvl,sss)
plt.plot(lut.wvl,lut.sp[0,:,0,28,70],'r')
print type(sss)
plt.title('Liquid water cloud with R$_{ef}$ = '+str(lut.ref[28])+' interpolated through $\\tau$='+str(lut.tau[70]))
plt.xlabel('Wavelength [nm]')
plt.ylabel('Radiance [Wm$^{-2}$sr$^{-1}$nm$^{-1}$]')


# ## Now run through the retrieval scheme

# In[35]:

import run_kisq_retrieval as rk
reload(rk)
import Sp_parameters as Sp
reload(Sp)
#del lut
#del meas


# In[36]:

print max(meas.good)


# In[37]:

(tau,ref,phase,ki) = rk.run_retrieval(meas,lut)


# In[38]:

del lut


# In[39]:

print meas.utc.shape
print len(meas.good), max(meas.good)


# In[40]:

reload(Sp)
from Sp_parameters import smooth


# In[41]:

fig,ax = plt.subplots(4,sharex=True)
ax[0].set_title('Retrieval results time trace')
ax[0].plot(meas.utc,tau,'rx')
ax[0].plot(meas.utc[meas.good[:,0]],smooth(tau[meas.good[:,0],0],20),'k')
ax[0].set_ylabel('$\\tau$')
ax[1].plot(meas.utc,ref,'g+')
ax[1].set_ylabel('R$_{ef}$ [$\\mu$m]')
ax[1].plot(meas.utc[meas.good[:,0]],smooth(ref[meas.good[:,0],0],20),'k')
ax[2].plot(meas.utc,phase,'k.')
ax[2].set_ylabel('Phase')
ax[2].set_ylim([-0.5,1.5])
ax[2].set_yticks([0,1])
ax[2].set_yticklabels(['liq','ice'])
ax[3].plot(meas.utc,ki)
ax[3].set_ylabel('$\\chi^{2}$')
ax[3].set_xlabel('UTC [Hours]')
ax[3].set_xlim([17,19.05])
plt.savefig(fp+'plots/TCAP_retri_results.png',dpi=600)
#plt.savefig(fp+'plots/TCAP_retri_results.eps')


# Now save the smoothed values

# In[42]:

print tau.shape


# In[43]:

tau[meas.good[:,0],0] = smooth(tau[meas.good[:,0],0],20)
ref[meas.good[:,0],0] = smooth(ref[meas.good[:,0],0],20)


# ### Now load the results from MODIS to compare

# Check the day of year

# In[44]:

import datetime
doy = datetime.datetime(2013,2,19)
print doy.timetuple().tm_yday
dd = datetime.datetime(2014,9,19)
print dd.timetuple().tm_yday


# Now import the hdf files of MODIS

# In[45]:

from mpl_toolkits.basemap import Basemap,cm


# In[46]:

myd06_file = fp+'MODIS\\MYD06_L2.A2013050.1725.006.2014260074007.hdf'
myd03_file = fp+'MODIS\\MYD03.A2013050.1725.006.2013051163424.hdf'
print os.path.isfile(myd03_file) #check if it exists
print os.path.isfile(myd06_file)


# In[47]:

import load_modis as lm
reload(lm)
if 'modis' in locals():
    del modis, modis_dicts
    import gc; gc.collect()
    


# Load the geolocated data:

# In[48]:

print 'after import'
import gc; gc.collect()
modis,modis_dicts = lm.load_modis(myd03_file,myd06_file)


# In[49]:

print modis_dicts['qa']['description']


# In[50]:

modis['qa'].shape


# In[51]:

bb = modis['qa'][100,100,:]
print bb
bb.dtype
print bb.astype('ubyte')
def bin(x):
    return ''.join(x & (1 << i) and '1' or '0' for i in range(7,-1,-1))
bin(bb[0])


# In[52]:

bin8 = lambda x : ''.join(reversed( [str((x >> i) & 1) for i in range(8)] ) )


# In[53]:

mqa = modis['qa']


# In[54]:

for i in mqa:
    print i


# In[55]:

if bin8(bb[1])[-4]:
    print 'yes'
bin8(bb[1])[-4]


# Now plot the resulting imagery

# In[56]:

#set up a easy plotting function
def tcap_map(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=-70,lat_0=42,
            llcrnrlon=-72, llcrnrlat=40,
            urcrnrlon=-66, urcrnrlat=45,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(-66,-72,6),labels=[0,0,0,1])
    m.drawparallels(np.linspace(40,45,5),labels=[1,0,0,0])
    return m


# In[57]:

figm = plt.figure()
axm = figm.add_axes([0.1,0.1,0.8,0.8])
m = tcap_map(axm)
x,y = m(modis['lon'],modis['lat'])
clevels = np.linspace(0,25,25)
cs = m.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.gist_ncar)
cbar = m.colorbar(cs)
cbar.set_label('$\\tau$')
axm.set_title('MODIS - AQUA Cloud optical Thickness')


# In[58]:

get_ipython().magic(u'pinfo plt.savefig')


# In[59]:

figm2,axm2 = plt.subplots(1,2,figsize=(13,13))
m1 = tcap_map(axm2[0])
m2 = tcap_map(axm2[1])
x,y = m1(modis['lon'],modis['lat'])
clevels = np.linspace(0,25,25)

cs1 = m1.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.gist_ncar)
cbar = m1.colorbar(cs1)
cbar.set_label('$\\tau$')
axm2[0].set_title('MODIS - AQUA Cloud optical Thickness')

clevels2 = np.linspace(0,50,25)
cs2 = m2.contourf(x,y,modis['ref'],clevels2,cmap=plt.cm.gist_earth)
cbar = m2.colorbar(cs2)
cbar.set_label('R$_{ef}$ [$\\mu$m]')
axm2[1].set_title('MODIS - AQUA Cloud effective radius')
figm2.subplots_adjust(wspace=0.3)
plt.show()
plt.savefig(fp+'plots/modis_only_tau_ref_comp.png',dpi=600)
plt.savefig(fp+'plots/modis_only_tau_ref_comp.pdf',bbox='tight')


# Now load the aircraft telemetry onto the plot

# In[60]:

reload(lm)


# In[61]:

# load the ict file and check out the results
iwg = lm.load_ict(fp+'arm-iop/aaf.iwg1001s.g1.TCAP.20130219.145837.a1.dat')
print iwg.dtype.names
print iwg.dtype.names.index('Date_Time')
print iwg['Date_Time'][0]
print type(iwg['Date_Time'][0])
iwg['Lat']


# In[62]:

iwg_utch = np.array([i.hour+i.minute/60.+i.second/3600.+i.microsecond/3600000. for i in iwg['Date_Time']])


# In[86]:

fig = plt.figure()
fig.add_subplot(3,1,1)
ax = plt.plot(iwg_utch,iwg['Lat'])
plt.ylabel('Latitude')
plt.title('flight path')
fig.add_subplot(3,1,2)
plt.plot(iwg_utch,iwg['Lon'])
plt.ylabel('Longitude')
fig.add_subplot(3,1,3)
plt.plot(iwg_utch,iwg['GPS_MSL_Alt'])
plt.xlabel('UTC [hours]')
plt.ylabel('Altitude [m]')


# ## interpolate the lat and lons and alts to retrieved values

# In[64]:

from scipy import interpolate


# In[65]:

flat = interpolate.interp1d(iwg_utch,iwg['Lat'])
meas.lat = flat(meas.utc)
flon = interpolate.interp1d(iwg_utch,iwg['Lon'])
meas.lon = flon(meas.utc)
falt = interpolate.interp1d(iwg_utch,iwg['GPS_MSL_Alt'])
meas.alt = falt(meas.utc)


# Now plot on top of the maps

# In[92]:

figm2,axm2 = plt.subplots(1,2,figsize=(13,13))
m1 = tcap_map(axm2[0])
m2 = tcap_map(axm2[1])
x,y = m1(modis['lon'],modis['lat'])
clevels = np.linspace(0,25,25)

cs1 = m1.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.gist_ncar)
cbar = m1.colorbar(cs1)
cbar.set_label('$\\tau$')
axm2[0].set_title('MODIS - AQUA Cloud optical Thickness')
x1,y1 = m1(meas.lon,meas.lat)
#m1.scatter(x1,y1,c=tau,cmap=plt.cm.gist_ncar,marker='o',vmin=clevels[0],vmax=clevels[-1],alpha=0.5,edgecolors='k',linewidth=0.25)
xt1, yt1 = m1(-70.5,42.2)
m1.plot(x1,y1,'k',lw=2)
axm2[0].text(xt1,yt1,'G-1',color='k')

clevels2 = np.linspace(0,50,25)
cs2 = m2.contourf(x,y,modis['ref'],clevels2,cmap=plt.cm.gist_earth)
cbar = m2.colorbar(cs2)
cbar.set_label('R$_{ef}$ [$\\mu$m]')
axm2[1].set_title('MODIS - AQUA Cloud effective radius')
#m2.scatter(x1,y1,c=ref,cmap=plt.cm.gist_earth,marker='o',vmin=clevels2[0],vmax=clevels2[-1],alpha=0.5,edgecolors='k',linewidth=0.25)
xt2, yt2 = m2(-70.5,42.2)
axm2[1].text(xt2,yt2,'G-1',color='k')
m2.plot(x1,y1,'k',lw=2)
figm2.subplots_adjust(wspace=0.3)
plt.savefig(fp+'plots/modis_g1_tau_ref_path.png',dpi=600,transparent=True)
#plt.savefig(fp+'plots/modis_g1_tau_ref_comp.pdf',bbox='tight')
#plt.show()


# In[94]:

figm2,axm2 = plt.subplots(1,2,figsize=(13,13))
m1 = tcap_map(axm2[0])
m2 = tcap_map(axm2[1])
x,y = m1(modis['lon'],modis['lat'])
clevels = np.linspace(0,25,25)

cs1 = m1.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.gist_ncar)
cbar = m1.colorbar(cs1)
cbar.set_label('$\\tau$')
axm2[0].set_title('MODIS - AQUA Cloud optical Thickness')
x1,y1 = m1(meas.lon,meas.lat)
m1.scatter(x1,y1,c=tau-3.0,cmap=plt.cm.gist_ncar,marker='o',vmin=clevels[0],vmax=clevels[-1],alpha=0.5,edgecolors='k',linewidth=0.25)
xt1, yt1 = m1(-70.5,42.2)
plt.text(xt1,yt1,'G-1',color='k')

clevels2 = np.linspace(0,50,25)
cs2 = m2.contourf(x,y,modis['ref'],clevels2,cmap=plt.cm.gist_earth)
cbar = m2.colorbar(cs2)
cbar.set_label('R$_{ef}$ [$\\mu$m]')
axm2[1].set_title('MODIS - AQUA Cloud effective radius')
m2.scatter(x1,y1,c=ref-3.0,cmap=plt.cm.gist_earth,marker='o',vmin=clevels2[0],vmax=clevels2[-1],alpha=0.5,edgecolors='k',linewidth=0.25)
xt2, yt2 = m2(-70.5,42.2)
plt.text(xt2,yt2,'G-1',color='k')

figm2.subplots_adjust(wspace=0.3)
plt.savefig(fp+'plots/modis_g1_tau_ref_comp.png',dpi=600,transparent=True)
#plt.savefig(fp+'plots/modis_g1_tau_ref_comp.pdf',bbox='tight')
plt.show()


# Now find the modis points along flight track that match the most

# In[67]:

from Sp_parameters import closestindex,startprogress,progress,endprogress
import scipy.spatial


# In[68]:

if max(meas.good) > meas.utc.size:
    meas.good = sm.good
print len(meas.good)
print meas.utc.size
print len(meas.lon)
print meas.sp.size
print meas.good.max()
plt.figure()
plt.plot(meas.good)


# In[69]:

imodis = np.logical_and(np.logical_and(modis['lon']>min(meas.lon[meas.good])-0.02 , modis['lon']<max(meas.lon[meas.good])+0.02),
                        np.logical_and(modis['lat']>min(meas.lat[meas.good])-0.02 , modis['lat']<max(meas.lat[meas.good])+0.02))


# In[70]:

wimodis = np.where(imodis)
print np.shape(wimodis)


# In[71]:

def spherical_dist(pos1, pos2, r=3958.75):
    pos1 = pos1 * np.pi / 180
    pos2 = pos2 * np.pi / 180
    cos_lat1 = np.cos(pos1[..., 0])
    cos_lat2 = np.cos(pos2[..., 0])
    cos_lat_d = np.cos(pos1[..., 0] - pos2[..., 0])
    cos_lon_d = np.cos(pos1[..., 1] - pos2[..., 1])
    return r * np.arccos(cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))


# In[72]:

N1 = modis['lon'][imodis].size
modis_grid = np.hstack([modis['lon'][imodis].reshape((N1,1)),modis['lat'][imodis].reshape((N1,1))])
print N1
#measurement
N2 = len(meas.good)
print N2
meas_grid = np.hstack([np.array(meas.lon[meas.good]).reshape((N2,1)),np.array(meas.lat[meas.good]).reshape((N2,1))])
meas_in = meas_grid.astype(int)
print len(meas_grid[0])


# Test if the spherical dist works for one point along the track

# In[73]:

d = spherical_dist(meas_grid[0],modis_grid)
print d.shape
print np.argmin(d)
print len(wimodis[0])
print len(wimodis[1])
print wimodis[0][np.argmin(d)]
print wimodis[1][np.argmin(d)]


# In[74]:

print meas.lat[0]
print meas.lon[0]
print modis['lon'][292,891]
print modis['lat'][292,891]


# In[75]:

meas.ind = np.array([meas.good.ravel()*0,meas.good.ravel()*0])
print np.shape(meas.ind)
meas.ind[0,0] = 2
meas.ind[1,0] = 3
print meas.ind[:,0]


# The spherical distance works for one point along track, now loop through all values

# In[76]:

startprogress('Running through flight track')
for i in xrange(meas.good.size):
    d = spherical_dist(meas_grid[i],modis_grid)
    meas.ind[0,i] = wimodis[0][np.argmin(d)]
    meas.ind[1,i] = wimodis[1][np.argmin(d)]
    progress(float(i)/len(meas.good)*100)
endprogress()


# In[77]:

print modis['tau'].shape
print np.shape(modis['tau'][meas.ind])
print modis['tau'][meas.ind[0,:],meas.ind[1,:]]
print meas.utc[meas.good].ravel().shape
print meas.good.shape


# ## Now plot MODIS and the retrieval result

# In[78]:

fig = plt.figure()
plt.plot(modis['lat'][meas.ind[0,:],meas.ind[1,:]],modis['lon'][meas.ind[0,:],meas.ind[1,:]])
plt.title('Nearest MODIS points')
plt.ylabel('Latitude')
plt.xlabel('longitude')


# In[79]:

fig = plt.figure()
plt.plot(meas.utc[meas.good].ravel(),modis['tau'][meas.ind[0,:],meas.ind[1,:]])
plt.title('raw tau values')
plt.xlabel('UTC [Hours]')


# In[80]:

modis_dicts['phase']


# In[81]:

fig,ax = plt.subplots(4,sharex=True)
ax[0].set_title('Retrieval results time trace')
ax[0].plot(meas.utc,tau,'rx',label='4STAR')
#ax[0].plot(meas.utc[meas.good[:,0]],smooth(tau[meas.good[:,0],0],20),'k')
ax[0].plot(meas.utc[meas.good].ravel(),modis['tau'][meas.ind[0,:],meas.ind[1,:]],'m+',label='MODIS')
ax[0].set_ylabel('$\\tau$')
ax[0].set_ylim([0,50])
ax[0].legend()
ax[1].plot(meas.utc,ref,'g+',label='4STAR')
ax[1].set_ylabel('R$_{ef}$ [$\\mu$m]')
ax[1].plot(meas.utc[meas.good].ravel(),modis['ref'][meas.ind[0,:],meas.ind[1,:]],'m+',label='MODIS')
ax[1].set_ylim([0,60])
ax[1].legend()
#ax[1].plot(meas.utc[meas.good[:,0]],smooth(ref[meas.good[:,0],0],20),'k')
ax[2].plot(meas.utc,phase,'k.')
ax[2].set_ylabel('Phase')
ax[2].set_ylim([-0.5,1.5])
ax[2].set_yticks([0,1])
ax[2].set_yticklabels(['liq','ice'])
ax[2].plot(meas.utc[meas.good].ravel(),modis['phase'][meas.ind[0,:],meas.ind[1,:]]-1,'m+')
ax[3].plot(meas.utc,ki)
ax[3].set_ylabel('$\\chi^{2}$')
ax[3].set_xlabel('UTC [Hours]')
ax[3].set_xlim([17,19.05])
plt.savefig(fp+'plots/modis_4star_time_comp.png',dpi=600)
plt.savefig(fp+'plots/modis_4star_time_comp.pdf',bbox='tight')


# Redo plot above, but with only tau and ref for comparison

# In[82]:

fig,ax = plt.subplots(2,sharex=True)
ax[0].set_title('Retrieval results time trace')
ax[0].plot(meas.utc,tau,'rx',label='4STAR')
#ax[0].plot(meas.utc[meas.good[:,0]],smooth(tau[meas.good[:,0],0],20),'k')
ax[0].plot(meas.utc[meas.good].ravel(),modis['tau'][meas.ind[0,:],meas.ind[1,:]],'m+',label='MODIS')
ax[0].set_ylabel('$\\tau$')
ax[0].set_ylim([0,15])
ax[0].legend()
ax[1].plot(meas.utc,ref,'g+',label='4STAR')
ax[1].set_ylabel('R$_{ef}$ [$\\mu$m]')
ax[1].plot(meas.utc[meas.good].ravel(),modis['ref'][meas.ind[0,:],meas.ind[1,:]],'m+',label='MODIS')
ax[1].set_ylim([0,60])
ax[1].legend()
#ax[1].plot(meas.utc[meas.good[:,0]],smooth(ref[meas.good[:,0],0],20),'k')
#ax[2].plot(meas.utc,phase,'k.')
#ax[2].set_ylabel('Phase')
#ax[2].set_ylim([-0.5,1.5])
#ax[2].set_yticks([0,1])
#ax[2].set_yticklabels(['liq','ice'])
#ax[2].plot(meas.utc[meas.good].ravel(),modis['phase'][meas.ind[0,:],meas.ind[1,:]]-1,'m+')
#ax[3].plot(meas.utc,ki)
#ax[3].set_ylabel('$\\chi^{2}$')
ax[1].set_xlabel('UTC [Hours]')
ax[1].set_xlim([17.25,19.05])
plt.savefig(fp+'plots/modis_4star_time_comp_tau_ref.png',dpi=600)
#plt.savefig(fp+'plots/modis_4star_time_comp_tau_ref.pdf',bbox='tight')


# plot onto a map the points that were selected

# In[83]:

figs2,axs2 = plt.subplots(1,2,figsize=(13,13))
ma1 = tcap_map(axs2[0])
ma2 = tcap_map(axs2[1])
xr,yr = ma1(modis['lon'][meas.ind[0,:],meas.ind[1,:]],modis['lat'][meas.ind[0,:],meas.ind[1,:]])
xi,yi = ma1(meas.lon[meas.good].ravel(),meas.lat[meas.good].ravel()+0.03)
axs2[0].set_title('Optical thickness')
ma1.scatter(xr,yr,c=modis['tau'][meas.ind[0,:],meas.ind[1,:]],
            cmap=plt.cm.gist_ncar,marker='o',vmin=clevels[0],
            vmax=clevels[-1],alpha=0.5,edgecolors='#999999',linewidth=0.15)
ma1.scatter(xi,yi,c=tau[meas.good],
            cmap=plt.cm.gist_ncar,marker='o',vmin=clevels[0],
            vmax=clevels[-1],alpha=0.5,edgecolors='#FFFFFF',linewidth=0.15)

axs2[1].set_title('Effective Radius')
ma2.scatter(xr,yr,c=modis['ref'][meas.ind[0,:],meas.ind[1,:]],
            cmap=plt.cm.gist_earth,marker='o',vmin=clevels2[0],
            vmax=clevels2[-1],alpha=0.5,edgecolors='#999999',linewidth=0.15)
ma2.scatter(xi,yi,c=ref[meas.good],
            cmap=plt.cm.gist_earth,marker='o',vmin=clevels2[0],
            vmax=clevels2[-1],alpha=0.5,edgecolors='#FFFFFF',linewidth=0.15)


# ## Build the comparison of histograms

# build a filter for only showing 1-1 plots

# In[256]:

utc = meas.utc[meas.good].ravel()
print utc.shape
ut = (utc > 17.25) & (utc < 17.75)


# In[279]:

reload(Sp)
from Sp_parameters import doublenanmask
mtau,meastau = doublenanmask(modis['tau'][meas.ind[0,ut],meas.ind[1,ut]],tau[meas.good[ut]].ravel())


# In[280]:

print len(meastau)
print len(mtau)


# In[291]:

plt.figure()
plt.hist(mtau, bins=10, histtype='stepfilled', normed=True, color='m',alpha=0.7, label='Modis')
plt.hist(meastau, bins=40, histtype='stepfilled', normed=True, color='r',alpha=0.7, label='4STAR')
plt.title('Optical Thickness histogram')
plt.ylabel('Normed probability')
plt.xlabel('$\\tau$')
plt.axvline(np.nanmean(mtau),color='k')
plt.axvline(np.nanmean(meastau),color='k',label='Mean')
plt.axvline(np.median(meastau),color='k',linestyle='--',label='Median')
plt.axvline(np.median(mtau),color='k',linestyle='--')
plt.legend(frameon=False)
plt.xlim([0,15])
plt.savefig(fp+'plots/hist_modis_4star_tau.png',dpi=600)
plt.savefig(fp+'plots/hist_modis_4star_tau.pdf',bbox='tight')


# Now for Ref

# In[282]:

mref,measref = doublenanmask(modis['ref'][meas.ind[0,ut],meas.ind[1,ut]],ref[meas.good[ut]].ravel())


# In[283]:

plt.figure()
plt.hist(mref, bins=30, histtype='stepfilled', normed=True, color='m',alpha=0.7, label='Modis')
plt.hist(measref, bins=30, histtype='stepfilled', normed=True, color='g',alpha=0.7, label='4STAR')
plt.title('Effective Radius histogram')
plt.ylabel('Normed probability')
plt.xlabel('R$_{ef}$ [$\\mu$m]')
plt.axvline(np.nanmean(mref),color='k')
plt.axvline(np.nanmean(measref),color='k',label='Mean')
plt.axvline(np.median(measref),color='k',linestyle='--')
plt.axvline(np.median(mref),color='k',linestyle='--',label='Median')
plt.legend(frameon=False)
plt.savefig(fp+'plots/hist_modis_4star_ref.png',dpi=600)
plt.savefig(fp+'plots/hist_modis_4star_ref.pdf',bbox='tight')

