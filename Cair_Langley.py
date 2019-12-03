#!/usr/bin/env python
# coding: utf-8

# # Intro
# 
# Name:  
# 
#     Cair_Lanley
# 
# Purpose:  
# 
#     Create a Langley plot for the Cair tubes attached to 2STAR
#     Assumes that it is tracking at all times, does not filter for bad tracking
# 
# Input:
# 
#     none at command line
# 
# Output:
# 
#     dict and plots 
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - numpy
#     - Pyephem
#     - Sun_utils
#     - Cair_utils
# 
# Needed Files:
# 
#   - Cair raw calibrated comma delimited files (URC)
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, On flight from BWI to SJC, 2016-08-06
#     Modified: 

# # Load the required modules and prepare the paths

# In[1]:


import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import numpy as np


# In[2]:


import map_utils as mu
import Sun_utils as su
import Cair_utils as cu
import plotting_utils as pu


# In[3]:


from mpltools import color


# In[4]:


fp = 'C:\\Users\\sleblan2\\Research\\4STAR\\MLO_2016\\'


# # Read some files

# In[5]:


f = fp+'20160702_MLO5\\CAST_001_160702_090020_URC.csv'


# In[56]:


c = cu.read_Cair(fp+'20160702_MLO5\\CAST_001_160702_090020_URC.csv')


# In[7]:


c.keys()


# In[8]:


c['DateTimeUTC'][0]


# # Generate the bandwith filter functions and a new wavelength range 

# ## Prepare the filter functions

# In[57]:


def Gamma2sigma(Gamma):
    '''Function to convert FWHM (Gamma) to standard deviation (sigma)'''
    import numpy as np
    return Gamma * np.sqrt(2) / ( np.sqrt(2 * np.log(2)) * 2 )


# In[58]:


def gaussian(x_center,fwhm):
    'Function that generates a gaussian distribution and a new x array to fit it onto'
    import numpy as np
    x = np.linspace(x_center-2*fwhm,x_center+2*fwhm,41)
    dg = np.exp( - ((x-x_center)/Gamma2sigma(fwhm))**2 )
    ndg = dg/np.sum(dg)#*np.diff([x,x[-1]+1]))
    return ndg, x


# In[59]:


fwhm = [10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,15,30]


# In[60]:


band_f = []
band_wvl = []
for i,l in enumerate(c['wvl']):
    f,wl = gaussian(l,fwhm[i])
    band_f.append(f)
    band_wvl.append(wl)


# In[61]:


c['band_f'] = band_f
c['band_wvl'] = band_wvl


# ## Plot the filter functions

# In[14]:


fig = plt.figure()
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=ax)
for i,l in enumerate(c['wvl']):
    p = ax.plot(band_wvl[i],band_f[i])
    ax.axvline(l,color=p[0].get_color())
    ax.fill_between(band_wvl[i],band_f[i],0,color=p[0].get_color(),alpha=0.5)
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('Band filter function')
ax.set_title('Synthetic C-air bandwidth')
plt.savefig(fp+'C-air_synthetic_bandwidth.png',dpi=600,transparent=True)


# ## Get the ozone optical depth for each wavelenght band

# ### Get the ozone cross sections

# In[62]:


f_oz = fp+'..\\O3_273K_V3_0.dat'


# In[63]:


oz = np.genfromtxt(f_oz,skip_header=41)


# In[64]:


oz.shape


# In[65]:


gas = {}
gas['o3_coef'] = oz[:,1]
gas['o3_wvl'] = oz[:,0]


# In[66]:


gas['o3_du'] = 261.0
gas['o3_tau'] = gas['o3_coef']*2.54070E19*gas['o3_du']/1000.0 


# In[67]:


plt.figure()
plt.semilogy(gas['o3_wvl'],gas['o3_tau'])
plt.xlabel('Wavelength [nm]')
plt.ylabel('O3 $\\tau$')


# ### calculate the ozone tau for the band wavelengths

# In[68]:


gas['o3_tau_band'] = cu.calc_gas_tau(c['band_wvl'],gas['o3_wvl'],gas['o3_tau'])


# In[69]:


gas['o3_tau_band'].shape


# ### Special case, remove the dependence of Ozone at the lowest wavelength band

# In[70]:


gas['o3_tau_band'][0,:] = 0.0


# # Run analysis and calculate the airmass and rayleigh

# In[71]:


lat, lon, alt = 19.5365,-155.57615,3428.0


# ## Calculate the airmass and sza

# In[72]:


c = su.calc_sza_airmass(c['DateTimeUTC'],lat,lon,alt,c=c)


# In[17]:


c.keys()


# In[18]:


c['Lt'].shape


# ## get the rayleigh tau at center bands

# In[73]:


c['tau_rayleigh'],c['tau_rayleigh_err'] = cu.calc_rayleigh(c,press=680.0)


# In[20]:


c['wvl']


# ### Verify the calculated tau_rayleigh

# In[22]:


import matplotlib.dates as mdates
fmt = mdates.DateFormatter('%H:%M')


# In[22]:


fig = plt.figure()
for i,l in enumerate(c['wvl']):
    plt.plot(c['DateTimeUTC'],c['tau_rayleigh'][i,:],'.',label='{} nm'.format(l))
fig.axes[0].xaxis.set_major_formatter(fmt)
plt.legend(frameon=True)
plt.xlabel('UTC [H]')
plt.ylabel('$\\tau$')


# In[23]:


c['DateTimeUTC'][-1]


# In[24]:


fig = plt.figure()
plt.plot(c['DateTimeUTC'],c['m_ray'],'.')
fig.axes[0].xaxis.set_major_formatter(fmt)
plt.legend(frameon=True)
plt.xlabel('UTC [H]')
plt.ylabel('Rayleigh arimass')


# In[25]:


fig = plt.figure()
plt.plot(c['DateTimeUTC'],c['sunearthf'],'.')
fig.axes[0].xaxis.set_major_formatter(fmt)
plt.legend(frameon=True)
plt.xlabel('UTC [H]')
plt.ylabel('Sun earth distance factor')


# In[26]:


fig = plt.figure()
for i,l in enumerate(c['wvl']):
    plt.plot(c['DateTimeUTC'],np.exp(-np.array(c['m_ray'])*c['tau_rayleigh'])[i,:],'.',label='{} nm'.format(l))
fig.axes[0].xaxis.set_major_formatter(fmt)
plt.legend(frameon=True)
plt.xlabel('UTC [H]')
plt.ylabel('Transmission due to Rayleigh')


# ## Get the Rayleigh tau at wavelengths of filter bands

# In[74]:


c['tau_rayleigh_fl'],c['tau_rayleigh_fl_err'] = cu.calc_rayleigh_filter(c,band_wvl,press=680.0)


# ## Calculate the 'rateaero' or Lt_aero
# Which is the Lt values divided by the impact of rayleigh, and trace gases

# In[75]:


c['Lt_aero'] = c['Lt']/np.array(c['sunearthf'])/np.exp(-np.array(c['m_ray'])*c['tau_rayleigh'])


# In[25]:


c['Lt_aero'].shape


# ### Calculate the filter modified Lt_aero

# In[76]:


c['Lt_aero_fl'] = np.zeros_like(c['Lt_aero'])
c['tr_rayleigh'] = np.zeros_like(c['Lt_aero'])


# In[77]:


for i,l in enumerate(c['wvl']):
    tr = 1.0/np.exp(-np.array(c['m_ray'])[:,np.newaxis]*c['tau_rayleigh_fl'][i,:,:])
    c['tr_rayleigh'][i,:] = np.dot(tr,band_f[i])
    c['Lt_aero_fl'][i,:] = c['Lt'][i,:]/np.array(c['sunearthf'])*c['tr_rayleigh'][i,:]


# ### plot the filter modified Lt_aero

# In[32]:


import matplotlib.dates as mdates
fmt = mdates.DateFormatter('%H:%M')


# In[33]:


fig = plt.figure()
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(c['DateTimeUTC'],1.0/c['tr_rayleigh'][i,:],'.',label='{} nm'.format(l))
fig.axes[0].xaxis.set_major_formatter(fmt)
plt.legend(frameon=True)
plt.xlabel('UTC [H]')
plt.ylabel('Transmission due to filter modified Rayleigh')


# ### Calculate the filter modified and refined wit Ozone Lt_aero
# 

# In[78]:


c['tr_o3'] = np.zeros_like(c['Lt_aero'])


# In[250]:


c['m_o3']


# In[29]:


gas['o3_tau_band'].shape


# In[79]:


for i,l in enumerate(c['wvl']):
    tr = 1.0/np.exp(-np.array(c['m_o3'])[:,np.newaxis]*gas['o3_tau_band'][i,:])
    c['tr_o3'][i,:] = np.dot(tr,c['band_f'][i])
    c['Lt_aero_fl'][i,:] = c['Lt_aero_fl'][i,:]*c['tr_o3'][i,:]


# ### plot the ozone transmission

# In[268]:


fig = plt.figure()
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(c['DateTimeUTC'],1.0/c['tr_o3'][i,:],'.',label='{} nm'.format(l))
fig.axes[0].xaxis.set_major_formatter(fmt)
plt.legend(frameon=True,ncol=2,numpoints=1)
plt.xlabel('UTC [H]')
plt.ylabel('Transmission due to Ozone absorption')


# # Plot the resulting langleys

# ## find the portions of the day for the langley

# In[80]:


from datetime import datetime
import dateutil
from dateutil.tz import tzutc


# In[81]:


fl = np.where((c['DateTimeUTC']>datetime(2016,7,2,16,tzinfo=dateutil.tz.tzutc())) & (c['DateTimeUTC']<datetime(2016,7,2,21,tzinfo=dateutil.tz.tzutc())))[0]


# In[82]:


fl_mu = (np.array(c['m_aero'])[fl]<12.0)&(np.array(c['m_aero'])[fl]>1.5)


# In[141]:


for i,l in enumerate(c['wvl']):
    fig = plt.figure()
    plt.plot(c['DateTimeUTC'],c['Lt'][i,:],'.',color='grey',label='all data')
    plt.plot(c['DateTimeUTC'][fl],c['Lt'][i,fl],'.',color='blue',label='time subset')
    plt.plot(c['DateTimeUTC'][fl[fl_mu]],c['Lt'][i,fl[fl_mu]],'.',color='green',label='airmass subset')
    plt.ylabel('Lt [$\\mu$W/(sr cm$^{2}$ nm)]')
    fig.axes[0].xaxis.set_major_formatter(fmt)
    plt.title('Lt {} nm calibrated values for {:%Y-%m-%d}'.format(c['wvl'][i],c['DateTimeUTC'][0]))
    plt.legend(loc=0)
    plt.savefig(fp+'{:%Y%m%d}_timeseries_Lt_{:04.0f}.png'.format(c['DateTimeUTC'][0],c['wvl'][i]),transparent=True,dpi=600)


# In[195]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(c['DateTimeUTC'],c['Lt'][i,:],'.',color='grey')
    #plt.plot(c['DateTimeUTC'][fl],c['Lt'][i,fl],'.',label='{} nm'.format(l))
    plt.plot(c['DateTimeUTC'][fl[fl_mu]],c['Lt'][i,fl[fl_mu]],'o',label='{} nm'.format(l),markeredgecolor='None')
plt.ylabel('Lt [$\\mu$W/(sr cm$^{2}$ nm)]')
plt.xlabel('UTC [HH:MM]')
fig.axes[0].xaxis.set_major_formatter(fmt)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
plt.title('Lt calibrated values for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.legend(loc='center left',bbox_to_anchor=(1.0,0.5),frameon=True,ncol=2,numpoints=1,columnspacing=0.1)
plt.savefig(fp+'{:%Y%m%d}_timeseries_Lt.png'.format(c['DateTimeUTC'][0]),transparent=True,dpi=600)


# ## Plot the langley for each wavelength

# ### Plot the langleys seperated by channel

# In[158]:


fig,ax = plt.subplots(5,4,sharex=True,figsize=(12,9))
ax = ax.ravel()
for i,l in enumerate(c['wvl']):
    ax[i].plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero'][i,fl]),'.',color='grey',label='bad data')
    ax[i].plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero'][i,fl[fl_mu]]),'.',color='b',label='good data')
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero'][i,fl[fl_mu]]),
                         color='k',ci=0.99,ax=ax[i],labels=False)
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero'][i,fl[fl_mu]]))
    ax[i].text(0.1,0.1,'b={}\nR$^2$={}'.format(p[0],r[0,1]**2),transform=ax[i].transAxes)
    
    ax[i].set_title('{} nm'.format(l))
    ax[i].grid()
    if i > 14: 
        ax[i].set_xlabel('Airmass')
    if i in [0,4,8,12,16]:
        ax[i].set_ylabel('log(Lt)')
    #ax[i].legend(frameon=True)
    ax[i].set_xlim(0,15)
    ax[i].set_ylim(min(np.log(c['Lt_aero'][i,fl[fl_mu]]))*0.98,max(np.log(c['Lt_aero'][i,fl[fl_mu]]))*1.02)
    
for tk in ax[15].get_xticklabels():
    tk.set_visible(True)
ax[i+1].axis('off')
plt.suptitle('C-air Refined Langley (Rayleigh subtracted) for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Langley_refined.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# In[83]:


fig,ax = plt.subplots(5,4,sharex=True,figsize=(12,9))
ax = ax.ravel()
c['Lt_0'] = np.zeros_like(c['wvl'])
c['langley_aod'] = np.zeros_like(c['wvl'])
c['langley_aod_err'] = np.zeros_like(c['wvl'])
c['Lt_0_err'] = np.zeros_like(c['wvl'])
c['Lt_0_r2'] = np.zeros_like(c['wvl'])
for i,l in enumerate(c['wvl']):
    ax[i].plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero_fl'][i,fl]),'.',color='grey',label='bad data')
    ax[i].plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]]),'.',color='b',label='good data')
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]]),
                         color='k',ci=0.99,ax=ax[i],labels=False)
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]]))
    c['Lt_0'][i] = p[0]
    c['Lt_0_err'][i] = perr[0]
    c['langley_aod'][i] = p[1]*-1.0
    c['langley_aod_err'][i] = perr[1]
    c['Lt_0_r2'][i] = r[0,1]**2
    ax[i].text(0.1,0.1,'b={:2.6}\nR$^2$={:1.5}'.format(p[0],r[0,1]**2),transform=ax[i].transAxes)
    
    ax[i].set_title('{} nm'.format(l))
    ax[i].grid()
    if i > 14: 
        ax[i].set_xlabel('Airmass')
    if i in [0,4,8,12,16]:
        ax[i].set_ylabel('log(Lt)')
    #ax[i].legend(frameon=True)
    ax[i].set_xlim(0,15)
    ax[i].set_ylim(min(np.log(c['Lt_aero_fl'][i,fl[fl_mu]]))*0.98,max(np.log(c['Lt_aero_fl'][i,fl[fl_mu]]))*1.02)
    
for tk in ax[15].get_xticklabels():
    tk.set_visible(True)
ax[i+1].axis('off')
plt.suptitle('C-air Refined Langley (Rayleigh subtracted bandwidth adjusted) for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Langley_refined_bandwidth.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# ### Plot the langleys on the same figure

# In[36]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero'][i,fl]),'.',color='grey')
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero'][i,fl[fl_mu]]))
    plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero'][i,fl[fl_mu]]),'.',
             label='{} nm\nr$^2$={:1.3}'.format(l,r[0,1]**2))
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero'][i,fl[fl_mu]]),color='k',ci=0.99,labels=False)
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero'][i,fl[fl_mu]]))
    #ax[i].text(0.1,0.1,'b={}\nR$^2$={}'.format(p[0],r[0,1]**2),transform=ax[i].transAxes)
    #ax[i].set_title('{} nm'.format(l))
plt.grid()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
plt.xlim(0,15)
plt.ylim(9,12)
plt.xlabel('Airmass')
plt.ylabel('log(Lt)')

plt.title('C-air Refined Langley (Rayleigh subtracted) for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Langley_refined_one.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# ### Plot the langleys, but adjusted for bandwidth

# In[39]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero_fl'][i,fl]),'.',color='grey')
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]]))
    plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]]),'.',
             label='{} nm\nr$^2$={:1.3}'.format(l,r[0,1]**2))
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]]),color='k',ci=0.99,labels=False)
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]]))
    #ax[i].text(0.1,0.1,'b={}\nR$^2$={}'.format(p[0],r[0,1]**2),transform=ax[i].transAxes)
    #ax[i].set_title('{} nm'.format(l))
plt.grid()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
plt.xlim(0,15)
plt.ylim(9,12)
plt.xlabel('Airmass')
plt.ylabel('log(Lt)')

plt.title('C-air Refined Langley (Rayleigh subtracted adjusted for bandwidth) for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Langley_refined_bandwitdth_one.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# In[270]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero_fl'][i,fl])/c['Lt_0'][i],'.',color='grey')
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])/c['Lt_0'][i])
    plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])/c['Lt_0'][i],'.',
             label='{} nm\nr$^2$={:1.3}\n$\\tau$={:1.2}'.format(l,r[0,1]**2,c['langley_aod'][i]))
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])/c['Lt_0'][i],
                         color='k',ci=0.99,labels=False)
    print l, p[1]
plt.grid()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
plt.xlim(0,15)
plt.ylim(0.8,1)
plt.xlabel('Airmass')
plt.ylabel('log(Lt)/log(Lt$_0$)')

plt.title('C-air Refined Langley (Rayleigh subtracted adjusted for bandwidth) for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Langley_refined_bandwitdth_norm_one.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# In[271]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero_fl'][i,fl])/c['Lt_0'][i],'.',color='grey')
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])/c['Lt_0'][i])
    plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])/c['Lt_0'][i],'.',
             label='{} nm\nr$^2$={:1.3}\n$\\tau$={:1.2}'.format(l,r[0,1]**2,c['langley_aod'][i]))
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])/c['Lt_0'][i],
                         color='k',ci=0.99,labels=False)
plt.grid()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
plt.xlim(0,15)
plt.ylim(0.95,1)
plt.xlabel('Airmass')
plt.ylabel('log(Lt)/log(Lt$_0$)')

plt.title('C-air Refined Langley (Rayleigh subtracted adjusted for bandwidth) for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Langley_refined_bandwitdth_norm_zoom_one.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# In[272]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero_fl'][i,fl])-c['Lt_0'][i],'.',color='grey')
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i])
    plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],'.',
             label='{} nm\nr$^2$={:1.3}\n$\\tau$={:1.2}'.format(l,r[0,1]**2,c['langley_aod'][i]))
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],
                         color='k',ci=0.99,labels=False)
    print l, p[1]
plt.grid()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
plt.xlim(0,15)
plt.ylim(-1,0)
plt.xlabel('Airmass')
plt.ylabel('log(Lt) - log(Lt$_0$)')

plt.title('C-air Refined Langley (Rayleigh subtracted adjusted for bandwidth) for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Langley_refined_bandwitdth_subtract_zoom_one.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# In[273]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
for i,l in enumerate(c['wvl']):
    plt.plot(np.array(c['m_aero'])[fl],c['Lt_aero_fl'][i,fl]/np.exp(c['Lt_0'][i]),'.',color='grey')
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],c['Lt_aero_fl'][i,fl[fl_mu]]/np.exp(c['Lt_0'][i]))
    plt.plot(np.array(c['m_aero'])[fl][fl_mu],c['Lt_aero_fl'][i,fl[fl_mu]]/np.exp(c['Lt_0'][i]),'.',
             label='{} nm\nr$^2$={:1.3}\n$\\tau$={:1.2}'.format(l,r[0,1]**2,c['langley_aod'][i]))
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],c['Lt_aero_fl'][i,fl[fl_mu]]/np.exp(c['Lt_0'][i]),
                         color='k',ci=0.99,labels=False)
    print l,p[1]
plt.grid()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
plt.xlim(0,15)
#plt.ylim(-1,0)
plt.xlabel('Airmass')
plt.ylabel('Lt/Lt$_0$')

plt.title('C-AIR ratio of raw to Lt$_0$(Rayleigh subtracted adjusted for bandwidth) for {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Lt_Lt0_one.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# ## Plot the AOD spectra

# In[37]:


from matplotlib.ticker import FormatStrFormatter


# In[113]:


plt.figure()
plt.loglog(c['wvl'],c['langley_aod'],'o',label='C-AIR')
plt.errorbar(c['wvl'], c['langley_aod'],marker='o',linestyle='None',yerr=c['langley_aod_err'],color='b')
plt.ylim(10**-4,10**-1)
plt.xlim(300,2200)
plt.ylabel('Aerosol Optical Depth')
plt.xlabel('Wavelength [nm]')
plt.grid(which='both')
plt.gca().xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
plt.gca().xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
plt.title('C-AIR Aerosol optical depth from Langley on {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_Langley_AOD.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# ### Plot the AATS data on top of C-AIR

# In[114]:


xa,ya = pu.plotmatfig(fp+'..\\figs\\AATS20160702amstdev_mult3LangleyAOD.fig')


# In[117]:


aats_wvl,aats_tau = xa[0],ya[0]


# In[125]:


plt.figure()
plt.loglog(c['wvl'],c['langley_aod'],'o',label='C-AIR')
plt.errorbar(c['wvl'], c['langley_aod'],marker='o',linestyle='None',yerr=c['langley_aod_err'],color='b')
plt.loglog(aats_wvl*1000.0,aats_tau,'or',label='AATS')
plt.legend(loc=3,frameon=False,numpoints=1)
plt.ylim(10**-4,10**-1)
plt.xlim(300,2200)
plt.ylabel('Aerosol Optical Depth')
plt.xlabel('Wavelength [nm]')
plt.grid(which='both')
plt.gca().xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
plt.gca().xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
plt.title('Aerosol optical depth from Langley on {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))
plt.savefig(fp+'{:%Y%m%d}_CAIR_AATS_Langley_AOD.png'.format(c['DateTimeUTC'][0]),dpi=600,transparent=True)


# ## Make a plot with AATS data and C-AIR

# In[123]:


pu.plotmatfig(fp+'AATS20160702amstdev_mult3LangleyperV0calc.fig')


# In[39]:


x,y = pu.plotmatfig(fp+'AATS20160702amstdev_mult3LangleyperlogV0calc.fig')


# ### Extract the AATS data from the plot

# In[40]:


len(x)


# In[41]:


aats = {}
aats['wvl'] = [353.3,380.0,451.2,499.4,520.4,605.8,675.1,779.1,864.5,1019.1,1241.3,1558.5,2139.1]


# In[42]:


aats['V'] = np.zeros((len(x)/2,len(y[0])))*np.nan
aats['V_good'] = np.zeros((len(x)/2,len(y[0])))*np.nan
aats['m'] = np.zeros((len(x)/2,len(y[0])))*np.nan
aats['m_good'] = np.zeros((len(x)/2,len(y[0])))*np.nan


# In[43]:


for i in range(0,len(x),2):
    aats['m'][i/2,:] = x[i]
    aats['m_good'][i/2,:len(x[i+1])] = x[i+1]
    aats['V'][i/2,:] = y[i]
    aats['V_good'][i/2,:len(y[i+1])] = y[i+1]


# ### Make the plot of log I - log I0 for cAIR and AATS

# In[287]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(121)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
col = []
for i,l in enumerate(c['wvl']):
    plt.plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero_fl'][i,fl])-c['Lt_0'][i],'.',color='grey')
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i])
    pp = plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],'.',
             label='{:4.0f} nm'.format(l))
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],
                         color='k',ci=0.99,labels=False)
    col.append(pp[0].get_color())
    #print l, p[1]
plt.grid()
#box = plt.gca().get_position()
#plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
#plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
l = plt.legend(frameon=True,loc='bottom left',ncol=4,markerscale=0,handlelength=0,handletextpad=0,columnspacing=0.5)
for i,text in enumerate(l.get_texts()):
    text.set_color(col[i])
plt.xlim(0,15)
plt.ylim(-1,0)
plt.xlabel('Airmass')
plt.ylabel('log(Lt) - log(Lt$_0$)')
ax.set_title('C-AIR')
plt.suptitle('Refined Langley from Mauna Loa on {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))

ax = fig.add_subplot(122)
color.cycle_cmap(len(aats['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
col = []
for i,l in enumerate(aats['wvl']):
    plt.plot(aats['m'][i,:],aats['V'][i,:],'.',color='grey')
    r = np.corrcoef(aats['m_good'][i,:],aats['V_good'][i,:])
    pp = plt.plot(aats['m_good'][i,:],aats['V_good'][i,:],'.',
             label='{:4.0f} nm'.format(l))
    p,perr = pu.plot_lin(aats['m_good'][i,:],aats['V_good'][i,:],
                         color='k',ci=0.99,labels=False)
    col.append(pp[0].get_color())
    #print l, p[1]
plt.grid()
#box = plt.gca().get_position()
#plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
#plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
l = plt.legend(frameon=True,loc='bottom left',ncol=3,markerscale=0,handlelength=0,handletextpad=0,columnspacing=0.5)
for i,text in enumerate(l.get_texts()):
    text.set_color(col[i])
plt.xlim(0,15)
plt.ylim(-1,0)
plt.xlabel('Airmass')
plt.ylabel('log(V) - log(V$_0$)')
ax.set_title('AATS-14')

#plt.title('Refined Langley from Mauna Loa on {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))


plt.savefig(fp+'{:%Y%m%d}_CAIR_AATS_Langley_subtract.png'.format(c['DateTimeUTC'][0]),
            dpi=600,transparent=True)


# In[84]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(121)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
col = []
for i,l in enumerate(c['wvl']):
    plt.plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero_fl'][i,fl])-c['Lt_0'][i],'.',color='grey')
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i])
    pp = plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],'.',
             label='{:4.0f} nm'.format(l))
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],
                         color='k',ci=0.99,labels=False)
    col.append(pp[0].get_color())
    #print l, p[1]
plt.grid()
#box = plt.gca().get_position()
#plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
#plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
l = plt.legend(frameon=True,loc='bottom left',ncol=4,markerscale=0,handlelength=0,handletextpad=0,columnspacing=0.5)
for i,text in enumerate(l.get_texts()):
    text.set_color(col[i])
plt.xlim(0,14)
plt.ylim(-0.3,0)
plt.xlabel('Airmass')
plt.ylabel('log(Lt) - log(Lt$_0$)')
ax.set_title('C-AIR')
plt.suptitle('Refined Langley from Mauna Loa on {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))

ax = fig.add_subplot(122)
color.cycle_cmap(len(aats['wvl'])+1,cmap=plt.cm.gist_ncar,ax=fig.axes[0])
col = []
for i,l in enumerate(aats['wvl']):
    plt.plot(aats['m'][i,:],aats['V'][i,:],'.',color='grey')
    r = np.corrcoef(aats['m_good'][i,:],aats['V_good'][i,:])
    pp = plt.plot(aats['m_good'][i,:],aats['V_good'][i,:],'.',
             label='{:4.0f} nm'.format(l))
    p,perr = pu.plot_lin(aats['m_good'][i,:],aats['V_good'][i,:],
                         color='k',ci=0.99,labels=False)
    col.append(pp[0].get_color())
    #print l, p[1]
plt.grid()
#box = plt.gca().get_position()
#plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
#plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.0,0.5),loc='center left',ncol=2)
l = plt.legend(frameon=True,loc='bottom left',ncol=3,markerscale=0,handlelength=0,handletextpad=0,columnspacing=0.5)
for i,text in enumerate(l.get_texts()):
    text.set_color(col[i])
plt.xlim(0,14)
plt.ylim(-0.15,0)
plt.xlabel('Airmass')
plt.ylabel('log(V) - log(V$_0$)')
ax.set_title('AATS-14')

#plt.title('Refined Langley from Mauna Loa on {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))


plt.savefig(fp+'{:%Y%m%d}_CAIR_AATS_Langley_subtract_zoom_adjusted.png'.format(c['DateTimeUTC'][0]),
            dpi=600,transparent=True)


# In[52]:


i_cair = [2,4,5,6,11,14,15,16,17]
i_aats = [1,2,3,4,6,7,8,9,10]


# In[91]:


c_cair_aats = plt.cm.jet(len(i_cair))


# In[103]:


c_cair_aats = plt.cm.jet(np.arange(len(i_cair))/float(len(i_cair)))


# In[51]:


import cmaps
cmaps.cmaps()


# In[93]:


c_cair_aats


# In[110]:


fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(121)
color.cycle_cmap(len(c['wvl'])+1,cmap=plt.cm.Pastel1,ax=fig.axes[0])
col = []
j = 0
for i,l in enumerate(c['wvl']):
    if i in i_cair:
        co = c_cair_aats[j]
        j+=1
        alpha = 1.0
    else:
        co = None
        alpha = 0.5
    plt.plot(np.array(c['m_aero'])[fl],np.log(c['Lt_aero_fl'][i,fl])-c['Lt_0'][i],'.',color='grey',alpha=alpha)
    r = np.corrcoef(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i])
    pp = plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],'.',
             label='{:4.0f} nm'.format(l),color=co,alpha=alpha)
    p,perr = pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],
                         color='k',ci=0.99,labels=False,alpha=alpha)
    col.append(pp[0].get_color())
j = 0
for i in i_cair:
    plt.plot(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],'.',c=c_cair_aats[j])
    pu.plot_lin(np.array(c['m_aero'])[fl][fl_mu],np.log(c['Lt_aero_fl'][i,fl[fl_mu]])-c['Lt_0'][i],color='k',ci=0.99,labels=False)
    j += 1
plt.grid()
l = plt.legend(frameon=True,loc='bottom left',ncol=4,markerscale=0,handlelength=0,handletextpad=0,columnspacing=0.5)
for i,text in enumerate(l.get_texts()):
    text.set_color(col[i])
plt.xlim(0,14)
plt.ylim(-0.3,0)
plt.xlabel('Airmass')
plt.ylabel('log(Lt) - log(Lt$_0$)')
ax.set_title('C-AIR')
plt.suptitle('Refined Langley from Mauna Loa on {:%Y-%m-%d}'.format(c['DateTimeUTC'][0]))

ax = fig.add_subplot(122)
color.cycle_cmap(len(aats['wvl'])+1,cmap=plt.cm.Pastel1,ax=ax)
col = []
j = 0
for i,l in enumerate(aats['wvl']):
    if i in i_aats:
        co = c_cair_aats[j]
        j += 1
        alpha = 1.0
    else:
        co = None
        alpha = 0.5
    plt.plot(aats['m'][i,:],aats['V'][i,:],'.',color='grey',alpha=alpha)
    r = np.corrcoef(aats['m_good'][i,:],aats['V_good'][i,:])
    pp = plt.plot(aats['m_good'][i,:],aats['V_good'][i,:],'.',
             label='{:4.0f} nm'.format(l),color=co,alpha=alpha)
    p,perr = pu.plot_lin(aats['m_good'][i,:],aats['V_good'][i,:],
                         color='k',ci=0.99,labels=False,alpha=alpha)
    col.append(pp[0].get_color())
j = 0
for i in i_aats:
    plt.plot(aats['m_good'][i,:],aats['V_good'][i,:],'.',c=c_cair_aats[j])
    pu.plot_lin(aats['m_good'][i,:],aats['V_good'][i,:],color='k',ci=0.99,labels=False)
    j += 1
plt.grid()
l = plt.legend(frameon=True,loc='bottom left',ncol=3,markerscale=0,handlelength=0,handletextpad=0,columnspacing=0.5)
for i,text in enumerate(l.get_texts()):
    text.set_color(col[i])
plt.xlim(0,14)
plt.ylim(-0.15,0)
plt.xlabel('Airmass')
plt.ylabel('log(V) - log(V$_0$)')
ax.set_title('AATS-14')
plt.savefig(fp+'{:%Y%m%d}_CAIR_AATS_Langley_highlight.png'.format(c['DateTimeUTC'][0]),
            dpi=600,transparent=True)


# In[ ]:




