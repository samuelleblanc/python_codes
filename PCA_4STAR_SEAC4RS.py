
# coding: utf-8

# Name:  
# 
#     PCA_4STAR_SEAC4RS
# 
# Purpose:  
# 
#     Python script that is used to view 4STAR full spectral raw counts
#     Animates the raw over time
# 
# Calling Sequence:
# 
#     python Animate_4STAR_v2.py
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     animated figures 
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - load_utils: for mat2py_time and toutc functions
#     - Sp_parameters: for find_closest function
#     - mpl_toolkits
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - General star.mat: or other 4STAR star file containing raw spectra

# In[ ]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
get_ipython().magic(u'matplotlib inline')
import numpy as np, h5py
#import plotly.plotly as py
import scipy.io as sio
import math
import os
from load_utils import mat2py_time, toutc
from Sp_parameters import find_closest


# In[6]:

# set the basic directory path
fp = 'C:\\Users\\sleblan2\\Research\\4STAR\\SEAC4RS\\'
datestr = '20130816'
filesep = '\\'
print filesep


# In[7]:

file_name = fp+'20130816\\20130816starsun.mat'#'20141002_ARISE_Flight_16\\20141002flt16starsun.mat' #'20141002_ARISE_Flight_16\\20141002flt16star.mat'


# ## Load the 4STAR data and select time points

# In[8]:

star = sio.loadmat(file_name)
star.keys()


# In[9]:

if 't' in star.keys():
    star['utc'] = toutc(mat2py_time(star['t']))
else:
    star['utc'] = toutc(mat2py_time(star['vis_sun']['t'][0][0]))
#star['vis_sun']['t'][0][0]


# In[14]:

figure;
plt.plot(star['utc'],'b.')


# In[15]:

star['good'] = np.where((star['utc']>12.0) & 
                        (star['utc']<23.0) & 
                        (star['Str'].flatten()!=0) & 
                        (star['sat_time'].flatten()==0) & 
                        (np.isfinite(star['tau_aero'][:,319])))[0]
print len(star['good'])
print star['good'][0]


# In[16]:

if 't' in star.keys():
    rate = star['rateaero'][star['good'],:]
    raw = star['rawcorr'][star['good'],:]
    m_aero = star['m_aero'][star['good']]
    print rate.shape
else:
    raw_vis = star['vis_sun']['raw'][0][0][star['good'],:]
    print raw_vis.shape


# ## Plot time traces at select wavelengths

# In[17]:

print star['tau_aero'].shape
print star['w'].shape
iw = find_closest(star['w'][0],np.array([380,400,430,500,515,635,875,1000,1140,1200,1600])/1000.0)
print star['w'][0].shape
print iw


# In[21]:

plt.figure()
color.cycle_cmap(len(iw)+1,cmap=plt.cm.gist_ncar)
for i in iw:
    plt.plot(star['utc'][star['good']],star['tau_aero'][star['good'],i], label=str('%.3f $\mu$m' % star['w'][0,i]))
plt.xlabel('UTC [h]')
plt.ylabel('AOD')
plt.ylim([0,0.35])
plt.legend(bbox_to_anchor=[1,0.1],loc=3)
plt.savefig(fp+datestr+'_AOD_time.png',dpi=600)


# In[23]:

plt.figure()
plt.plot(star['w'].T,star['tau_aero'][star['good'][100],:])
plt.ylim([0,0.5])
plt.xlim([0.3,1.75])
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('AOD')
plt.title('AOD at single time point')
plt.savefig(fp+datestr+'_AOD_spectrum.png',dpi=600)


# In[25]:

fi,ax = plt.subplots(figsize=(8,13))
aodlevels = np.linspace(0,0.3,26)
cs = plt.contourf(star['w'][0,:],star['utc'][star['good']],star['tau_aero'][star['good'],:],aodlevels,cmap=plt.cm.gist_ncar,extend='both')
cbar = plt.colorbar(cs,pad=0.15)
cbar.set_label('AOD')
cbar.ax.set_aspect(30)
ax.set_xlabel('Wavelength [$\mu$m]')
ax.set_ylabel('UTC [h]')
ax.set_title('AOD for SEAC4RS flight ('+datestr+')')
ax2 = ax.twinx()
ax2.set_ylabel('airmass')
ax1Ys = ax.get_yticks()
ax2Ys = []
for Y in ax1Ys:
    ax2Ys.append('%4.1f' % star['m_aero'][(np.abs(star['utc']-Y)).argmin()][0])
ax2.set_yticks(ax1Ys)
ax2.set_ybound(ax.get_ybound())
ax2.set_yticklabels(ax2Ys)
plt.draw()
plt.savefig(fp+datestr+'_AOD_time_spectrum.png',dpi=600)


# In[27]:

plt.figure(figsize=(8,13))
aodlevels = np.linspace(0.005,.035,26)
cs = plt.contourf(star['w'][0,:],star['utc'][star['good']],star['tau_aero'][star['good'],:],aodlevels,cmap=plt.cm.gist_ncar,extend='both')
cbar = plt.colorbar(cs)
cbar.set_label('AOD')
plt.xlabel('Wavelength [$\mu$m]')
plt.xlim([0.4,0.5])
plt.ylabel('UTC [h]')
plt.title('AOD for SEAC4RS flight near 430 nm ('+datestr+')')
plt.savefig(fp+datestr+'_AOD_time_spectrum_zoom.png',dpi=600)


# ## Prepare for testing PCA analysis of raw spectra

# In[28]:

from sklearn.decomposition import PCA


# In[29]:

if 'rate' in locals():
    print rate.shape
    rate_fl = where(np.isfinite(rate[:,0]))
    wvl = star['w'][0]
else:
    print raw_vis.shape
    raw_fl = where(np.isfinite(raw_vis[:,0]))
    wvl = np.genfromtxt('C:\\Users\\sleblan2\\Research\\4STAR\\vis_4STAR_wvl.dat')


# In[30]:

wvl.shape


# Get the log of the rate

# In[31]:

ratel = np.log(rate)
rate_m = ratel/m_aero


# In[32]:

print rate_m.shape
print m_aero.shape
print ratel.shape


# In[33]:

print ratel.shape
print isfinite(ratel[:,400]).all()


# In[34]:

wflrate = [isfinite(ratel).any(0)][0]
print where(~wflrate)
wflrate = range(210,1020)+range(1050,1530)


# In[35]:

flrate = [isfinite(ratel[:,wflrate]).all(1)][0]
print flrate.shape
print where(~flrate)


# In[36]:

ratel1 = ratel[flrate,:]
print ratel1.shape
ratel2 = ratel1[:,wflrate]
print ratel2.shape
if len(wvl) != len(wflrate):
    wvl = star['w'][0][wflrate]


# In[37]:

print isfinite(ratel2).all()


# In[38]:

raw1 = raw[flrate,:]
raw2 = raw1[:,wflrate]
rate_m1 = rate_m[flrate,:]
rate_m2 = rate_m1[:,wflrate]


# In[39]:

if 'rate' in locals():
    ratel1 = ratel[flrate,:]
    pca = PCA(n_components=20)
    pca.fit(ratel2)
    evals = pca.explained_variance_ratio_
    pca_r = PCA(n_components=20)
    pca_r.fit(raw2)
    evals_r = pca_r.explained_variance_ratio_
    pca_m = PCA(n_components=20)
    pca_m.fit(rate_m2)
    evals_m = pca_m.explained_variance_ratio_
else:
    pca = PCA(n_components=20)
    pca.fit(raw_vis)
    evals = pca.explained_variance_ratio_


# In[40]:

plt.plot(evals_r,label='each pca')
plt.plot(evals_r.cumsum(),'r',label='cumulative info')
plt.title('Variability of rate counts explained by number of pca')
plt.legend(frameon=False)


# In[41]:

plt.plot(evals, label='each pca')
plt.plot(evals.cumsum(),'r', label='cumulative info')
plt.title('Vaiability of raw counts in pca')
plt.legend(frameon=False)


# In[42]:

print evals.cumsum()
print evals
coms = pca.components_
print coms.shape
coms_r = pca_r.components_
coms_m = pca_m.components_


# In[43]:

fig2,ax2 = plt.subplots(5,2,figsize=(15,11))
ax2 = ax2.ravel()
for i in range(10):
    if i ==0 :
        ax2[i].plot(wvl, coms_r[i,:]*(-1.0))
    else:
        ax2[i].plot(wvl, coms_r[i,:])
    ax2[i].set_title('pca %d' % i)
    ax2[i].grid()
    #ax3[i].set_xlim([0.3,1.75])
    if i >7:
        ax2[i].set_xlabel('Wavelength [$\mu$m]')
plt.suptitle('Raw counts (rate)')
plt.tight_layout()
plt.savefig(fp+datestr+'_langley_pca_raw.png',dpi=600)


# In[44]:

for i,j in list(enumerate(wvl)): 
    print i,j


# In[45]:

fig3,ax3 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax3 = ax3.ravel()
yli = ([-0.2,0.00001],
       [-0.05,0.02],
       [-0.045,0.015],
       [-0.05,0.03],
       [-0.01,0.01],
       [-0.03,0.04],
       [-0.03,0.02],
       [-0.02,0.02],
       [-0.025,0.009],
       [-0.025,0.04])
for i in range(10):
    ax3[i].plot(wvl, coms[i,:])
    ax3[i].set_title('pca %d' % i)
    ax3[i].grid()
    #ax3[i].set_xlim([0.3,1.75])
    if i < 10: 
        ax3[i].set_ylim(yli[i])
    if i >7:
        ax3[i].set_xlabel('Wavelength [$\mu$m]')
plt.suptitle('Log of raw rate zoomed')
plt.tight_layout()
plt.savefig(fp+datestr+'_langley_pca_log_zoom.png',dpi=600)


# In[46]:

fig4,ax4 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax4 = ax4.ravel()
for i in range(10):
    ax4[i].plot(wvl, coms[i,:])
    ax4[i].set_title('pca %d' % i)
    ax4[i].grid()
    if i >7:
        ax4[i].set_xlabel('Wavelength [$\mu$m]')
plt.suptitle('Log of raw rate')
plt.tight_layout()
plt.savefig(fp+datestr+'_langley_pca_log.png',dpi=600)


# In[47]:

fig5,ax5 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax5 = ax5.ravel()
yli = ([-0.025,0.0],
       [-0.030,-0.020],
       [-0.010,0.01],
       [0.015,0.03],
       [-0.03,0.03],
       [0.00,0.013],
       [-0.01,0.1],
       [-0.015,0.01],
       [-0.005,0.008],
       [-0.005,0.002])
for i in range(10):
    ax5[i].plot(wvl, coms[i,:])
    ax5[i].set_title('pca %d' % i)
    ax5[i].grid()
    ax5[i].set_xlim([0.3,0.65])
    if i < 10: 
        ax5[i].set_ylim(yli[i])
    if i >7:
        ax5[i].set_xlabel('Wavelength [$\mu$m]')
plt.suptitle('Log of raw rate')
plt.tight_layout()
plt.savefig(fp+datestr+'_langley_pca_log_vis_430nm.png',dpi=600)


# In[48]:

print coms[0,:]


# In[49]:

fig6,ax6 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax6 = ax6.ravel()
yli = ([-0.05,0.01],
       [-0.045,0.015],
       [-0.04,0.02],
       [-0.02,0.015],
       [-0.04,0.04],
       [-0.09,0.05],
       [-0.05,0.08],
       [-0.03,0.02],
       [-0.025,0.01],
       [-0.015,0.015])
for i in range(10):
    ax6[i].plot(wvl, coms_m[i,:])
    ax6[i].set_title('pca %d' % i)
    ax6[i].grid()
    #ax3[i].set_xlim([0.3,1.75])
    if i < 10: 
        ax6[i].set_ylim(yli[i])
    if i >7:
        ax6[i].set_xlabel('Wavelength [$\mu$m]')
plt.suptitle('Log of raw rate divided by airmass')
plt.tight_layout()
plt.savefig(fp+datestr+'_langley_pca_log_airmass.png',dpi=600)


# ## Use PCA to decompose and recompose the spectra

# In[50]:

pca_r_low = PCA(n_components=0.999)
pca_r_low.fit(raw2)
evals_r_low = pca_r_low.explained_variance_ratio_
print pca_r_low.n_components


# In[51]:

traw2 = pca_r_low.transform(raw2)
newraw2 = pca_r_low.inverse_transform(traw2)


# In[52]:

print raw2.shape
scores_r = pca_r.transform(raw2)
print scores_r.shape
new_raw2 = pca_r.inverse_transform(scores_r)
print coms_r.shape
print new_raw2.shape


# Run PCA on rate_aero

# In[53]:

rate1 = rate[flrate,:]
rate2 = rate1[:,wflrate[1:]]
print np.isfinite(rate2).shape
print np.isfinite(rate2).all()
print rate2.shape
print np.isfinite(rate2[0,:]).all()


# In[54]:

pp = np.mean(rate2,0)
print pp.shape
print where(~np.isfinite(pp))
print pp


# In[55]:

fll = range(0,1842)+range(1844,2836)+range(2861,5159)+range(5161,6753)+range(6792,len(rate2[:,500]))
r2_vis = rate2[fll,0:809]
r2_nir = rate2[fll,809:]
wvl_vis = wvl[0:809]
wvl_nir = wvl[810:]


# In[56]:

pca_vis = PCA(n_components=0.999)
pca_vis.fit(r2_vis)
tr2_vis = pca_vis.transform(r2_vis)
n2_vis = pca_vis.inverse_transform(tr2_vis)
pca_nir = PCA(n_components=0.995)
pca_nir.fit(r2_nir)
tr2_nir = pca_nir.transform(r2_nir)
n2_nir = pca_nir.inverse_transform(tr2_nir)


# In[57]:

print pca_vis.n_components
print pca_nir.n_components


# In[58]:

pca_rate = PCA(n_components=5)
pca_rate.fit(rate2)
print pca_rate.n_components


# In[59]:

print wvl[800:850]
print wvl[809]


# In[60]:

trate2 = pca_rate.transform(rate2)
newrate2 = pca_rate.inverse_transform(trate2)


# In[61]:

print rate2.shape
print trate2.shape
print newrate2.shape


# In[62]:

plt.plot(rate2[:,500], label='original singal')
plt.plot(newrate2[:,500],'r', label='reconstructed signal from pca')
plt.title('Time trace over time of Rate at %s $\mu$m' % wvl[500])
plt.legend(frameon=False,loc='lower left')


# In[63]:

plt.plot(r2_vis[:,500], label='original signal')
plt.plot(n2_vis[:,500],'r', label='reconstructed signal from pca')
plt.xlabel('Count number')
plt.ylabel('Raw rate')
plt.title('Time trace over time of Rate at %s $\mu$m' % wvl_vis[500])
plt.legend(frameon=False,loc='lower left')


# In[64]:

plt.plot(r2_nir[:,300], label='original signal')
plt.plot(n2_nir[:,300],'r', label='reconstructed signal from pca')
plt.title('Time trace over time of Rate at %s $\mu$m' % wvl_nir[300])
plt.legend(frameon=False,loc='lower left')
plt.xlabel('Count number')
plt.ylabel('Raw rate')


# In[65]:

print wvl_vis.shape
print r2_vis[500,:].shape
print wvl_nir.shape
print r2_nir[500,:].shape


# In[66]:

figw,axw = plt.subplots(1,2,figsize=(15,5))
plt.tight_layout()
axw = axw.ravel()
axw[0].plot(wvl_vis,r2_vis[500,:], linewidth=3,label='Original Signal')
axw[1].plot(wvl_nir,r2_nir[500,:], linewidth=3,label='Original Signal')
axw[0].plot(wvl_vis,n2_vis[500,:],'r', linewidth=2, label='Reconstructed Signal from pca')
axw[1].plot(wvl_nir,n2_nir[500,:],'r', linewidth=2, label='Reconstructed Signal from pca')

for i in range(2):
    axw[i].set_xlabel('Wavelength [$\mu$m]')
    axw[i].set_ylabel('Rate counts')

axw[1].legend(frameon=False)

plt.suptitle('Single spectrum')
plt.tight_layout()
plt.savefig(fp+datestr+'_single_spectrum_reconstructed.png',dpi=600)


# In[67]:

figw,axw = plt.subplots(1,2,figsize=(15,5))
plt.tight_layout()
axw = axw.ravel()
axw[0].plot(wvl_vis,(r2_vis[500,:]-n2_vis[500,:])/r2_vis[500,:]*100.0)
axw[1].plot(wvl_nir,(r2_nir[500,:]-n2_nir[500,:])/r2_nir[500,:]*100.0)

for i in range(2):
    axw[i].set_xlabel('Wavelength [$\mu$m]')
    axw[i].set_ylabel('Rate counts difference [$\%$]')

plt.suptitle('Single spectrum difference between original and reconstructed from PCA')
plt.tight_layout()
plt.savefig(fp+datestr+'_single_spectrum_reconstructed_difference.png',dpi=600)


# In[68]:

ff, axs = plt.subplots(5,figsize=(15,9))
for a in range(5):
    axs[a].plot(coms_r2[a]*trate2[500,a]+pca_rate.mean_)
    axs[a].set_title('component %d' % a)

#plt.plot(coms_r2[0]*trate2[500,0])
print trate2[500,0]


# In[ ]:

print m_aero.shape
print n2_vis[:,300].shape


# In[ ]:

maero = m_aero[flrate]
mm = maero[fll]
print mm.shape


# In[ ]:

figl,axl = plt.subplots(1,2,figsize=(15,5))
plt.tight_layout()
axl = axl.ravel()

axl[0].plot(mm,np.log(r2_vis[:,300]),label='Original signal',color='b',linestyle=':')
axl[0].plot(mm,np.log(n2_vis[:,300]),label='After pca reconstruction with %s componenents' %pca_vis.n_components,color='r',linewidth=0.5)
axl[0].set_title('Langley plot at %s $\mu$m' % wvl_vis[300])
axl[0].set_xlabel('Airmass')
axl[0].set_ylabel('Log (Rate of counts)')
axl[0].legend(frameon=False)

axl[1].plot(mm,np.log(r2_nir[:,300]),label='Original signal',color='b',linestyle=':')
axl[1].plot(mm,np.log(n2_nir[:,300]),label='After pca reconstruction with %s componenents' %pca_nir.n_components,color='r',linewidth=0.5)
axl[1].set_title('Langley plot at %s $\mu$m' % wvl_nir[300])
axl[1].set_xlabel('Airmass')
axl[1].set_ylabel('Log (Rate of counts)')
axl[1].legend(frameon=False)

plt.savefig(fp+datestr+'_langley_smooth.png',dpi=600)


# In[ ]:



