#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     COSR_DARE
# 
# Purpose:  
# 
#     To Build the COSR DARE calculations
#   
# Input:
# 
#     arguments
#   
# Output:
# 
#     Figure and save files
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - load_utils.py
#     - matplotlib
#     - numpy
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-10-04
#     Modified: 

# # Prepare python environment

# In[1]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
from datetime import datetime
from scipy.interpolate import UnivariateSpline, interp1d
from scipy import interpolate 
import pandas as pd
import matplotlib.pyplot as plt
import Sun_utils as su


# In[195]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


import Run_libradtran as RL
import os
import write_utils as wu
from tqdm.notebook import tqdm


# In[4]:


name = 'COSR'
vv = 'v3'


# In[5]:


fp =getpath('COSR')


# In[6]:


day = '20180609'
day = '20180618'
#day = '20180624'
#day = '20180625'


# In[7]:


#get the paths
fp_rtm = getpath('rtm')
fpp_in = fp_rtm+'input/{}_{}_CRE_{}/'.format(name,day,vv)
fpp_out = fp_rtm+'output/{}_{}_CRE_{}/'.format(name,day,vv)
fp_uv = getpath('uvspecb')+'uvspec'
fp_librad = getpath('libradtranb')+'data/'


# # Load files

# ## Load the in situ files

# In[8]:


situ = pd.read_csv(fp+'data_other/{}_nephclap.csv'.format(day))


# In[9]:


situ


# In[10]:


insitu = situ.to_dict('list')


# In[11]:


insitu.keys()


# In[12]:


insitu['ssa_500nm'] = np.array(insitu['totScatCalc_500nm'])/np.array(insitu['extCalc500nm'])


# In[13]:


insitu['extCalc500nm'] = np.array(insitu['extCalc500nm'])


# In[14]:


gu = pd.to_datetime(situ['DateTimeUTC']).to_list()
insitu['utc'] = np.array([g.hour+g.minute/60.0+g.second/3600.0 for g in gu])


# In[15]:


plt.figure()
plt.plot(insitu['utc'],insitu['ssa_500nm'],'.',label='calculated from CLAP+neph at 500 nm')
plt.plot(insitu['utc'],Sp.smooth(insitu['ssa_500nm'],60,old=True),'-r',label='smooth over 60s')
plt.legend()
plt.xlabel('UTC [h]')
plt.ylabel('SSA at 500 nm')
plt.title('SSA from in situ calculated on: {}'.format(day))
plt.ylim(0.5,1.0)
plt.grid()
#plt.savefig(fp+'plots/SSA_500_CLAP_neph_smooth_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[634]:


plt.figure(figsize=(7,2.5))
plt.plot(insitu['utc'],insitu['ssa_500nm'],'.',alpha=0.3,label='calculated from CLAP+neph at 500 nm')
plt.plot(insitu['utc'],Sp.smooth(insitu['ssa_500nm'],60,old=True),'-r',label='smooth over 60s')
plt.legend()
plt.xlabel('UTC [h]')
plt.ylabel('SSA at 500 nm')
#plt.title('SSA from in situ calculated on: {}'.format(day))
plt.ylim(0.45,1.0)
plt.xlim(15.4,20.3)
plt.grid()
plt.tight_layout()
plt.savefig(fp+'plots/SSA_500_CLAP_neph_smooth_transp_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[16]:


ssa_insitu = Sp.smooth(insitu['ssa_500nm'],60,old=True)


# In[898]:


len(ssa_insitu)


# ### Apply Altitude change filter

# Need to load the 4STAR spectra, then run the next cell first

# In[42]:


f_alt = interp1d(x=s['utc'],y=s['Alt'][:,0])
insitu['alt'] = f_alt(insitu['utc'])


# In[43]:


plt.figure()
plt.plot(insitu['utc'],Sp.smooth(insitu['ssa_500nm']*2000.0,60,old=True),'-r')
plt.plot(insitu['utc'],insitu['totScat_550nm'],'.g')
plt.plot(insitu['utc'],insitu['alt'])


# In[46]:


def running_std(x,n):
    'Function to do a running standard deviation on array (x) with window size (n)'
    q = x**2
    try: 
        q = np.convolve(q, np.ones((n, )), mode="valid")
    except ValueError:
        x = x.flatten()
        q = x**2
        q = np.convolve(q, np.ones((n, )), mode="valid")        
    s = np.convolve(x, np.ones((n, )), mode="valid")
    o = (q-s**2/n)/float(n-1)
    return o 


# In[49]:


nbox=30


# In[50]:


std_alt = running_std(insitu['alt'],nbox)


# In[51]:


infl = np.where(std_alt<50)[0]
infl_not = np.where(std_alt>50)[0]


# In[52]:


insitu['ssa_500nm'][infl].shape


# In[53]:


insitu['ssa_500nm'].shape


# In[54]:


ssa_in = insitu['ssa_500nm'][:]*2.0
ssa_in[infl_not] = np.nan
ssa_in = ssa_in/2.0


# In[55]:


ssa_in.shape


# In[56]:


#xout = np.convolve(x, np.ones(window)/window, 'same')
window = 15
xmasked, mask = Sp.nanmasked(ssa_in)
fx = interpolate.interp1d(insitu['utc'][mask],xmasked,bounds_error=False)
xinterp = fx(insitu['utc'])
ssaout = np.convolve(xinterp, np.ones(window)/window, 'same')

ssaout[~mask] = np.nan

#ssaout = np.convolve(ssa_in,np.ones(30)/30,'same')


# In[59]:


fig,ax1 = plt.subplots()

ax1.plot(insitu['utc'],insitu['alt'],label='alt all')
ax1.plot(insitu['utc'][infl],insitu['alt'][infl],'ob',label='good alt')

ax2 = ax1.twinx()
ax2.plot(insitu['utc'],insitu['ssa_500nm'],'--k',label='ssa all')
ax2.plot(insitu['utc'][infl],insitu['ssa_500nm'][infl],'xg',label='good ssa')
ax2.plot(insitu['utc'],ssa_in,'.y',label='presmooth')
ax2.plot(insitu['utc'],ssaout,'sc',label='smooth ssa')

plt.legend()


# In[60]:


plt.figure(figsize=(7,2.5))
plt.plot(insitu['utc'],insitu['ssa_500nm'],'.',color='lightgrey',alpha=0.3,label='calculated from CLAP+neph at 500 nm')
plt.plot(insitu['utc'][infl],insitu['ssa_500nm'][infl],'.',alpha=0.3,label='Filtered for altitude changes')
plt.plot(insitu['utc'],ssaout,'-',color='tab:red',label='smooth over 15s')
plt.legend()
plt.xlabel('UTC [h]')
plt.ylabel('SSA at 500 nm')
#plt.title('SSA from in situ calculated on: {}'.format(day))
plt.ylim(0.45,1.0)
plt.xlim(15.0,18.7)
plt.grid()
plt.tight_layout()
plt.savefig(fp+'plots/SSA_500_CLAP_neph_smooth_transp_filtered_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[61]:


ssa_insitu = ssaout


# ## Load the UHSAS files

# In[62]:


uh = sio.loadmat(fp+'data_other/{}_UHSAS_ECCC.mat'.format(day))


# In[63]:


uh.keys()


# In[64]:


uh['utc'] = lu.toutc(lu.mat2py_time(uh['t']))


# In[65]:


uh['utc'].shape


# In[66]:


uh['binDataConc'].sum(axis=1).shape


# In[67]:


uh['nConc'] = uh['binDataConc'].sum(axis=1)


# In[68]:


plt.figure()
plt.plot(uh['utc'],uh['binDataConc'].sum(axis=1),'.')


# ### Build a file for quick mie calculations

# In[462]:


np.argmin(abs(uh['utc']-16.527))


# In[454]:


uh['binDataConc'].shape


# In[463]:


plt.figure()
plt.plot(uh['binFrom'][0,:],uh['binDataConc'][5770,:])
plt.xlabel('Bin size [nm]')
plt.ylabel('#/cc')


# In[466]:


bin_averages = (uh['binFrom'][0,:] +uh['binTo'][0,:])/2.0/1000.0 #in micron


# In[467]:


bin_averages 


# In[468]:


fp


# In[479]:


with open(fp+'data_other/UHSAS_sizedistr_{}.dat'.format(day),'w') as f:
    for j,r in enumerate(bin_averages):
        f.write('{:0.7f} {:3.7f}\n'.format(r,uh['binDataConc'][5770,j]))


# In[480]:


wvl


# ### Read the mie calc

# File input to mie code from libradtran 2.0.2:
#     
#     mie_program MIEV0
#     refrac user 1.60 0.0032
#     size_distribution_file /data/sam/COSR/data_other/UHSAS_sizedistr_20180609.dat
#     wavelength /data/sam/COSR/data_other/mie_cal_wavelengths.dat
#     
# Command:
# 
#     ~/libradtran/libRadtran-2.0.2/bin/mie < mie_aero.inp > mie_aero.out

# In[483]:


mie_out = np.genfromtxt(fp+'data_other/mie_aero.out')


# In[484]:


mie_out


# In[487]:


mout = {}
mout['wvl'] = mie_out[:,0]
mout['refr'] = mie_out[:,1]
mout['refi'] = mie_out[:,2]
mout['qext'] = mie_out[:,3]
mout['ssa'] = mie_out[:,4]
mout['asy'] = mie_out[:,5]


# ## Load the flt table

# In[69]:


fp


# In[70]:


flttable = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format(day))


# In[323]:


flttable


# In[71]:


fromtime = flttable['FromTime'][flttable['FlightType']=='in plume']


# In[72]:


totime = flttable['ToTime'][flttable['FlightType']=='in plume']


# In[73]:


tt = fromtime.to_numpy()


# In[74]:


tt[0].second


# In[75]:


def time_utc(x):
    return np.array([y.hour+y.minute/60.0+y.second/3600.0 for y in x])


# In[76]:


from_utc = time_utc(fromtime.to_numpy())
to_utc = time_utc(totime.to_numpy())


# ## Load the 4STAR AOD

# In[18]:


s = sio.loadmat(fp+'os_data/4STAR_{}starsun.mat'.format(day))


# In[212]:


s.keys()


# In[19]:


s['utc'] = lu.toutc(lu.mat2py_time(s['t']))


# In[423]:


plt.figure()
plt.pcolor(s['w'],s['utc'],s['tau_aero'][:-1],cmap='gist_ncar',vmin=0,vmax=0.8)


# In[ ]:





# ### use the polyfit aod on the wavelength array

# In[20]:


s['tau_aero_polynomial'].shape


# In[21]:


wvl = np.array([0.25,0.35,0.4,0.5,0.675,0.87,0.995,1.2,1.4,1.6,2.1,3.2,4.9])


# In[22]:


s['aod'] = np.zeros((len(s['utc']),len(wvl)))
for i in xrange(len(s['utc'])):
    s['aod'][i,:] = np.exp(np.polyval([s['tau_aero_polynomial'][i,0],s['tau_aero_polynomial'][i,1],
                                       s['tau_aero_polynomial'][i,2]],
                                         np.log(wvl)))


# In[23]:


s['aod'].shape


# In[24]:


s['aod'][10:-1:200,:].shape


# In[25]:


plt.figure()
plt.plot(wvl,s['aod'][10:-1:200,:].T)
plt.ylim(0,1)
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.01,10)
plt.xlim(0.25,5)


# ### Alternative build of tau aero polynomials for AOD

# In[26]:


wvl = np.array([0.25,0.35,0.4,0.5,0.675,0.87,0.995,1.2,1.4,1.6,2.1,3.2,4.9])


# In[27]:


sai = s['aerosolcols'][0,:].astype(int)


# In[28]:


plt.figure()
plt.plot(s['tau_aero'][:,400],'.')


# In[580]:


i=13884


# In[29]:


s['tau_aero'].shape


# In[30]:


s['w'].shape


# In[31]:


plt.figure()
plt.plot(s['w'][0,:],s['tau_aero'][i,:],'-')
plt.plot(s['w'][0,sai],s['tau_aero'][i,sai],'.')
#plt.plot(wvl,s['aod'][i,:],'x-')
plt.ylim(0,0.5)


# In[737]:


np.where(flag)[0]


# In[738]:


5%10


# In[788]:


plt.figure(figsize=(8,3))
jfirst=True
for i in np.where(flag)[0]:
    if i%500 ==0:
        if jfirst:
            plt.plot(s['w'][0,:],s['tau_aero_subtract_all'][i,:],'-',color='grey',alpha=0.4,label='all spectra')
            jfirst=False
        else:
            plt.plot(s['w'][0,:],s['tau_aero_subtract_all'][i,:],'-',color='grey',alpha=0.4)
        plt.plot(s['w'][0,saii],s['tau_aero_subtract_all'][i,saii],'.',label='Good AOD at {:2.2f} UTC'.format(s['utc'][i]))
#plt.plot(wvl,s['aod'][i,:],'x-')
plt.ylim(0,0.5)
plt.xlim(0.3,1.75)
plt.grid()
plt.legend(bbox_to_anchor=(0.98, 1), loc='upper left', ncol=1)

plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('AOD [tau_aero_subtract_all]')
plt.title('AOD spectra on {}'.format(day))
plt.tight_layout()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.96, box.height])
plt.savefig(fp+'plots/Good_AOD_spectra_examples_{}.png'.format(day),dpi=600,transparent=True)


# In[777]:


from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


# In[787]:


plt.figure(figsize=(8,3))
jfirst=True
for i in np.where(flag)[0]:
    if i%500 ==0:
        if jfirst:
            plt.plot(s['w'][0,:]*1000.0,s['tau_aero_subtract_all'][i,:],'-',color='grey',alpha=0.4,label='all spectra')
            jfirst=False
        else:
            plt.plot(s['w'][0,:]*1000.0,s['tau_aero_subtract_all'][i,:],'-',color='grey',alpha=0.4)
        plt.plot(s['w'][0,saii]*1000.0,s['tau_aero_subtract_all'][i,saii],'.',label='Good AOD at {:2.2f} UTC'.format(s['utc'][i]))
#plt.plot(wvl,s['aod'][i,:],'x-')
plt.ylim(0,0.5)
plt.xlim(330,500)
plt.grid(which='major')

plt.gca().xaxis.set_minor_locator(MultipleLocator(5))
plt.gca().tick_params(which='minor', length=4, color='grey')

plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', ncol=1)

plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD [tau_aero_subtract_all]')
plt.title('AOD spectra on {}'.format(day))
plt.tight_layout()
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.96, box.height])
plt.savefig(fp+'plots/Good_AOD_spectra_examples_zoomed_{}.png'.format(day),dpi=600,transparent=True)


# In[752]:


s['w'][0,saii]


# In[583]:


plt.figure()
plt.plot(s['w'][0,:],s['tau_aero'][i,:],'-')
plt.plot(s['w'][0,saii],s['tau_aero'][i,saii],'.')
plt.ylim(0,0.5)


# In[32]:


saii = sai[((s['w'][0,sai]>0.364) & (s['w'][0,sai]<0.370)) | ((s['w'][0,sai]>0.435) & (s['w'][0,sai]<1.566)) | 
            ((s['w'][0,sai]>1.584) & (s['w'][0,sai]<1.597)) ]


# In[33]:


pl = su.logaod_polyfit(np.append(s['w'][0,saii],[2.2,2.7,3.5]),np.append(s['tau_aero'][i,saii],[s['tau_aero'][i,saii][-1]*2.0/3.0,s['tau_aero'][i,saii][-1]/2.0,s['tau_aero'][i,saii][-1]/4.0]),polynum=2)


# In[34]:


plt.figure()
plt.plot(s['w'][0,:],s['tau_aero'][i,:],'-')
plt.plot(s['w'][0,sai],s['tau_aero'][i,sai],'.')
#plt.plot(wvl,s['aod'][i,:],'x-')
plt.plot(wvl,np.exp(np.polyval(pl,np.log(wvl))))

plt.yscale('log')
plt.xscale('log')
plt.ylim(0.01,2)
plt.xlim(0.3,5)


# In[35]:


from scipy import polyfit


# In[36]:


tau_aero_good = np.array([np.append(s['tau_aero_subtract_all'][i,saii],[s['tau_aero_subtract_all'][i,saii][-1]*2.0/3.0,
                                                           s['tau_aero_subtract_all'][i,saii][-1]/2.0,s['tau_aero_subtract_all'][i,saii][-1]/10.0]) \
                          for i in xrange(len(s['t']))])


# In[37]:


poly = np.array([polyfit(np.log(np.append(s['w'][0,saii],[2.2,2.7,3.7])),np.log(aodd),2) for aodd in tau_aero_good])


# In[38]:


poly.shape


# In[39]:


s['paod'] = np.zeros((len(s['utc']),len(wvl)))
for i in xrange(len(s['utc'])):
    s['paod'][i,:] = np.exp(np.polyval(poly[i,:],np.log(wvl)))


# In[40]:


plt.figure()
plt.plot(wvl,s['paod'][10:-1:500,:].T)
plt.plot(s['w'][0,saii],s['tau_aero'][10:-1:500,saii].T,'.')
plt.ylim(0,1)
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.01,10)
plt.xlim(0.25,5)


# In[77]:


s['aod'] = s['paod']


# ### Load the flag files

# In[78]:


fmat = getpath('4STAR_data')


# In[79]:


with open (fmat+'starinfo_{}.m'.format(day), 'rt') as in_file:
    for line in in_file:
        if 'flagfilename ' in line:
            ff = line.split("'")[1]
sf = hs.loadmat(fmat+ff)


# In[80]:


sf.keys()


# In[81]:


flag = sf['manual_flags']['good'][0,:,0]


# In[82]:


flag.shape


# In[83]:


sum(flag)


# ## Load the solar spectrum

# In[615]:


fp_librad


# In[616]:


sol = np.genfromtxt(fp_librad+'solar_flux/kurudz_1.0nm.dat')


# In[617]:


sol.shape


# In[619]:


plt.figure()
plt.plot(sol[:,0],sol[:,1],'.')
plt.ylabel('Irradiance [mW/m$^2$nm]')
plt.title('Kurudz solar spectrum [1nm]')
plt.xlabel('Wavelength [nm]')


# In[620]:


sol[:,1].sum()


# In[623]:


plt.figure()
plt.plot(sol[:,0],sol[:,1].cumsum(),'.')
plt.ylabel('Cumulative Irradiance [mW/m$^2$nm]')
plt.title('Kurudz solar spectrum [1nm]')
plt.xlabel('Wavelength [nm]')


# In[624]:


plt.figure()
plt.plot(sol[:,0],sol[:,1].cumsum()/sol[:,1].sum()*100.0,'.')
plt.ylabel('Cumulative Irradiance [\%]')
plt.title('Kurudz solar spectrum [1nm]')
plt.xlabel('Wavelength [nm]')
plt.grid()


# ## Load the skyscan results

# In[84]:


sk_names = {
    '20180624':'4STAR_20180624_135_SKYP.created_20190329_003621.ppl_lv15.mat',
    '20180625':'4STAR_20180625_026_SKYP.created_20190507_213718.ppl_lv15.mat',
    '20180620':'4STAR_20180620_017_SKYP.created_20190507_225712.ppl_lv15.mat',
    '20180618':'4STAR_20180618_029_SKYA.created_20190507_232752.avg_lv10.mat',
    '20180609':'4STAR_20180609_041_SKYP.created_20190508_003116.ppl_lv15.mat',
    '20180705':'4STAR_20180705_061_SKYP.created_20190508_003930.ppl_lv15.mat'}


# In[85]:


sk_n = {
    '20180624':'135',
    '20180625':'026',
    '20180620':'017',
    '20180618':'029',
    '20180609':'041',
    '20180705':'061'}


# In[86]:


fp_name = sk_names[day]#'4STAR_20180624_135_SKYP.created_20190329_003621.ppl_lv15.mat'


# In[87]:


sky = sio.loadmat(fp+fp_name)
sky.keys()


# In[88]:


sky['refractive_index_real_r']


# In[89]:


sky['refractive_index_imaginary_r']


# # Plot out some data

# ## Plot out the retrieved skyscans

# In[451]:


plt.figure()
plt.plot(sky['radius'],sky['psd'])
plt.title('Size distribution')
plt.xlabel('Radius [micron]')
plt.ylabel('PSD')


# In[440]:


plt.figure()
plt.plot(sky['Wavelength'],sky['ssa_total'][0,:],'x-',label='SSA')
plt.plot(sky['Wavelength'],sky['g_tot'][0,:],'x-',label='ASY')
plt.plot(sky['Wavelength'],sky['sfc_alb'][0,:],'x-',label='Albedo')


plt.legend(frameon=False)
plt.xlabel('Wavelength [micron]')
plt.title('4STAR skyscan results from: \n' + fp_name)
plt.savefig(fp+'plots/4STAR_skyscan_result_{}_{}_{}_SKYP.png'.format(day,sk_n[day],vv),dpi=600,transparent=True)


# In[244]:


sky['g_tot'][-1]


# ## Expand the sky scans results to longer wavelengths

# In[90]:


#wvl = np.array([0.35,0.4,0.5,0.675,0.87,0.995,1.2,1.4,1.6,2.1,4.0])
f_asy = interp1d(np.append(sky['Wavelength'][:,0],[2.0,2.6]),np.append(sky['g_tot'][0,:],[sky['g_tot'][0,-1]+0.008,sky['g_tot'][0,-1]-0.058]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
asy = f_asy(wvl)
f_ssa = interp1d(np.append(sky['Wavelength'][:,0],[1.5,2.0,2.3]),np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]+0.003,sky['ssa_total'][0,-1]+0.0015,sky['ssa_total'][0,-1]-0.002]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
ssa = f_ssa(wvl)
f_alb = interp1d(np.append(np.append([0.25,0.35],sky['Wavelength'][:,0]),[1.65,2.2,3.0,5.0]),
                 np.append(np.append([0.01,0.01],sky['sfc_alb'][0,:]),[sky['sfc_alb'][0,-1]*0.5,sky['sfc_alb'][0,-1]*0.25,0.04,0.02]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
alb = f_alb(wvl)


# In[70]:


# for 20180624 only
if day=='20180624':
    f_asy = interp1d(np.append(sky['Wavelength'][:,0],[1.1,2.4]),np.append(sky['g_tot'][0,:],[sky['g_tot'][0,-1]-0.01,sky['g_tot'][0,-1]-0.13]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    asy = f_asy(wvl)
    f_ssa = interp1d(np.append(sky['Wavelength'][:,0],[1.1,1.2,1.3]),np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.012,sky['ssa_total'][0,-1]-0.022,sky['ssa_total'][0,-1]-0.034]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    ssa = f_ssa(wvl)


# In[162]:


if day=='20180625':
    f_asy = interp1d(np.append(sky['Wavelength'][:,0],[1.1,2.4]),np.append(sky['g_tot'][0,:],[sky['g_tot'][0,-1]-0.005,sky['g_tot'][0,-1]-0.008]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    asy = f_asy(wvl)
    f_ssa = interp1d(np.append(sky['Wavelength'][:,0],[1.1,1.2,1.3]),np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.012,sky['ssa_total'][0,-1]-0.022,sky['ssa_total'][0,-1]-0.034]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    ssa = f_ssa(wvl)


# In[286]:


if day=='20180609':
    f_asy = interp1d(np.append(sky['Wavelength'][:,0],[1.1,2.4]),np.append(sky['g_tot'][0,:],[sky['g_tot'][0,-1]-0.002,sky['g_tot'][0,-1]-0.005]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    asy = f_asy(wvl)
    f_ssa = interp1d(np.append(sky['Wavelength'][:,0],[1.1,1.2,1.3]),np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.008,sky['ssa_total'][0,-1]-0.016,sky['ssa_total'][0,-1]-0.020]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    ssa = f_ssa(wvl)


# In[91]:


np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.008,sky['ssa_total'][0,-1]-0.002])


# In[92]:


sky['Wavelength'][:,0]


# In[93]:


sky['sfc_alb'][0,:]


# In[442]:


plt.figure()
plt.plot(wvl,ssa,'-x',label='SSA')
plt.plot(wvl,asy,'-*',label='ASY')
plt.plot(wvl,alb,'-s',label='Albedo')
plt.xscale('log')
plt.xticks([0.25,0.35,0.5,0.75,1.0,2.0,4.0],[0.25,0.35,0.5,0.75,1.0,2.0,4.0])
plt.grid()
plt.legend()
plt.xlabel('Wavelength [micron]')
plt.title('Extended wavelength properties for DARE calc. {}'.format(day))
plt.savefig(fp+'plots/AERO_prop_for_DARE_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[529]:


plt.figure()

plt.plot(sky['Wavelength'],sky['ssa_total'][0,:],'x-',color='b',label='Retrieved SSA')
plt.plot(sky['Wavelength'],sky['g_tot'][0,:],'*-',color='orange',label='Retrieved ASY')
plt.plot(sky['Wavelength'],sky['sfc_alb'][0,:],'s-',color='g',label='MODIS Albedo')

plt.plot(mout['wvl']/1000.0,mout['asy'],':',color='peru',label='ASY from Mie calculations')
plt.plot(mout['wvl']/1000.0,mout['ssa'],':',color='aqua',label='SSA from Mie calculations')

plt.plot(wvl,ssa,'--x',label='extrapolated SSA',color='b',alpha=0.5)
plt.plot(wvl,asy,'--*',label='extrapolated ASY',color='orange',alpha=0.5)
plt.plot(wvl,alb,'--s',label='extrapolated Albedo',color='g',alpha=0.5)
plt.xscale('log')
plt.xticks([0.25,0.35,0.5,0.75,1.0,1.5,2.0,4.0],[0.25,0.35,0.5,0.75,1.0,1.5,2.0,4.0])
plt.grid()
plt.legend(frameon=False,loc=6)
plt.xlabel('Wavelength [$\\mu$m]')
#plt.title('Extended wavelength properties for DARE calc. {}'.format(day))
plt.savefig(fp+'plots/AERO_prop_extrapolated_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# ### Load the spectral albedo for grass

# In[508]:


albg_wvl = np.genfromtxt(fp+'data_other/splib07a_Wavelengths_BECK_Beckman_0.2-3.0_microns.txt',skip_header=1)
albg_alb = np.genfromtxt(fp+'data_other/splib07a_Lawn_Grass_GDS91_green_BECKa_AREF.txt',skip_header=1)


# In[513]:


albg_alb[albg_alb<0] = 0.0


# In[625]:


plt.figure(figsize=(7,3.5))

plt.plot(sky['Wavelength'],sky['ssa_total'][0,:],'x-',color='b',label='SSA Retrieved')
plt.plot(mout['wvl']/1000.0,mout['ssa'],':',color='aqua',label='SSA from Mie calculations')
plt.plot(wvl,ssa,'--x',label='SSA extrapolated',color='b',alpha=0.4)

plt.plot(sky['Wavelength'],sky['g_tot'][0,:],'*-',color='orange',label='ASY Retrieved')
plt.plot(mout['wvl']/1000.0,mout['asy'],':',color='peru',label='ASY from Mie calculations')
plt.plot(wvl,asy,'--*',label='ASY extrapolated',color='orange',alpha=0.4)

plt.plot(sky['Wavelength'],sky['sfc_alb'][0,:],'s-',color='g',label='Albedo from MODIS') #  [MCD43GF; Schaaf et al., ]
plt.plot(albg_wvl,albg_alb*0.3,':',color='lime',label='Albedo for Grass (scaled)') #[Kokaly et al., 2017]
plt.plot(wvl,alb,'--s',label='Albedo extrapolated',color='g',alpha=0.4)

i = 1201
plt.plot(s['w'][0,saii],s['tau_aero'][flag,:][i,saii],'.',color='red',label='AOD measured')
plt.plot(wvl,s['paod'][flag,:].T[:,i],'--.',color='coral',alpha=0.4,label='AOD polynomial expansion')

plt.xscale('log')
plt.xticks([0.25,0.35,0.5,0.75,1.0,1.5,2.0,4.0],[0.25,0.35,0.5,0.75,1.0,1.5,2.0,4.0])
plt.xlim(0.23,5.5)
plt.grid()

plt.xlabel('Wavelength [$\\mu$m]')
#plt.title('Extended wavelength properties for DARE calc. {}'.format(day))
plt.legend(frameon=True,bbox_to_anchor=(0.93,0.91),loc=2,fontsize='small')

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0+0.1, box.width * 0.76, box.height])

plt.savefig(fp+'plots/AERO_prop_extrapolated_plusalb_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# ## Get the vertical dependence of the extinction

# In[94]:


f_alt = interp1d(x=s['utc'],y=s['Alt'][:,0])
insitu['alt'] = f_alt(insitu['utc'])


# In[253]:


if day=='20180624':
    insitu['extCalc500nm'] = insitu['extCalc500nm']*10.0


# In[258]:


if day=='20180609':
    insitu['extCalc500nm'] = insitu['extCalc500nm']*10.0


# In[95]:


plt.figure()
plt.plot(insitu['extCalc500nm'],insitu['alt'],'.')


# In[96]:


plt.figure()
plt.plot(insitu['utc'],insitu['alt'])


# In[255]:


plt.figure()
plt.scatter(insitu['utc'],insitu['alt'],c=insitu['extCalc500nm'],cmap=plt.cm.magma,s=(insitu['extCalc500nm'])**1.5-20.0)
cb = plt.colorbar()
cb.set_label('Extinction at 500 nm [1/Mm]')
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')
plt.title('In Situ Extinction coefficients on {}'.format(day))
plt.savefig(fp+'plots/Extinction_UTC_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[75]:


np.isfinite(insitu['extCalc500nm'])


# In[97]:


binned_ext,binned_alt,binned_num = [],[],[]
for i in xrange(17):
    flaa = (insitu['alt']>=i*100.0) & (insitu['alt']<(i+1.0)*100.0) & (np.isfinite(insitu['extCalc500nm']))
    if flaa.any():
        binned_ext.append(insitu['extCalc500nm'][flaa])
        binned_alt.append(np.mean([i*100.0,(i+1.0)*100.0]))
        binned_num.append(len(insitu['extCalc500nm'][flaa]))


# In[98]:


plt.figure()
bp =plt.boxplot(binned_ext,positions=binned_alt,vert=False,showfliers=True,widths=100,showmeans=True,patch_artist=True)
plt.xlabel('Extinction Calculated 500 nm [1/Mm]')
plt.ylabel('Altitude [m]')
#plt.plot(s['angs_470_865'][s['fl_QA_angs']],s['GPS_Alt'][s['fl_QA_angs']],'.',alpha=0.005)
for b in bp['boxes']:
    b.set_facecolor('green')
    b.set_edgecolor('green')
    b.set_alpha(0.4)
for b in bp['means']:
    b.set_marker('o')
    b.set_color('firebrick')
    b.set_alpha(0.4)
for b in bp['whiskers']:
    b.set_linestyle('-')
    b.set_color('green')
    b.set_alpha(0.4)
for b in bp['caps']:
    b.set_alpha(0.4)
    b.set_color('green')
for b in bp['medians']:
    b.set_linewidth(4)
    b.set_color('gold')
    b.set_alpha(0.4)
for b in bp['fliers']:
    b.set_marker('.')
    b.set_alpha(0.2)
ext_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(ext_means[:,0],ext_means[:,1],'-k')
plt.xlim(0,100)
plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.savefig(fp+'plots/extinction_vertical_bins_clap_neph_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[99]:


ext_z = ext_means[:,1]/1000.0
ext_ = ext_means[:,0]/1000.0
aod_ = ext_.sum()/10.0


# In[100]:


ext_z = np.append(ext_z,ext_z[-1]+0.1)


# In[101]:


nz = len(ext_z)


# In[102]:


ext_ = np.append(ext_,ext_[-1]*0.0)


# In[103]:


ext_,aod_


# In[104]:


ext_


# # Now build the input files for DARE calculations
# 

# In[105]:


from write_utils import nearest_neighbor


# In[106]:


geo = {'zout':[0,3,100],'year':int(day[0:4]),'month':int(day[4:6]),'day':int(day[6:8]),'hour':12,'minute':0,'second':0}
aero_base = {'z_arr':(ext_z+0.05),'wvl_arr':wvl,'ssa':np.array([ssa,]*nz),'asy':np.array([asy,]*nz)}
source = {'integrate_values':True,'dat_path':fp_librad,
          'run_fuliou':True,'wvl_range':[350,4000]}
albedo = {'create_albedo_file':True,'alb':alb,'alb_wvl':wvl*1000.0}
cloud = {'dummy':None}


# In[107]:


# Old way with no altitude filtering
if False:
    fx_ssa_insitu = interp1d(insitu['utc'],ssa_insitu,bounds_error=False)
    ssa_u = fx_ssa_insitu(s['utc'])
    ssa_u[ssa_u<0.8] = np.nan
else:
    ssa_u = nearest_neighbor(insitu['utc'],ssa_insitu,s['utc'],dist=1/3600.0)


# In[108]:


def expand_ext_vert_and_spect(ext_,ext_z,aod_sp,alt,wvl):
    """
    create a 2d array of extintion (altitude, spectra). 
    Inputs:
        ext_: vertical profile of ext, 
        ext_z: altitudes of profile [km], 
        aod_sp: aod spectra, 
        alt: measured alt [km]
        wvl: wavelength array of the aod_sp [microns]
    """
    iz = np.argmin(abs(ext_z-alt))
    iw = np.argmin(abs(wvl-0.5))
    aod_ = ext_[iz:].sum()/10.0
    
    factor = aod_sp[iw]/aod_
    exts = np.array([aod_sp*e for e in ext_])*factor
    return exts 


# In[109]:


geo.pop('sza')


# ## Write the input files

# In[110]:


file_list = file(fp_rtm+'COSR_DARE_list_file_{d}_{v}.sh'.format(d=day,v=vv),'w')
file_list_clean = file(fp_rtm+'COSR_DARE_list_file_clean_{d}_{v}.sh'.format(d=day,v=vv),'w')
print 'Starting list file'
fpp_out = fp_rtm+'output/COSR_{d}_{v}/'.format(d=day,v=vv)
fpp_in = fp_rtm+'input/COSR_{d}_{v}/'.format(d=day,v=vv)
if not os.path.exists(fpp_out):
    os.mkdir(fpp_out)
if not os.path.exists(fpp_in):
    os.mkdir(fpp_in)

nu = len(s['utc'])
pbar = tqdm(total=nu)
for i,u in enumerate(s['utc']):
    aero = {}
    if flag[i] & np.isfinite(ssa_u[i]):
        aod = s['aod'][i,:]
        ext = expand_ext_vert_and_spect(ext_,ext_z,s['aod'][i,:],s['Alt'][i]/1000.0,wvl)
        iw = np.argmin(abs(aero_base['wvl_arr']-0.5))
        aero['ssa'] = aero_base['ssa']*ssa_u[i]/aero_base['ssa'][0,iw]
        aero['asy'] = aero_base['asy']
        aero['z_arr'] = aero_base['z_arr']
        aero['wvl_arr'] = aero_base['wvl_arr']*1000.0
        
        try: aero['ssa'][aero['ssa']<0.0] = 0.0
        except: pass
        try: aero['ssa'][aero['ssa']>1.0] = 1.0
        except: pass
        try: aero['asy'][aero['asy']<0.0] = 0.0
        except: pass
        try: aero['asy'][aero['asy']>1.0] = 1.0
        except: pass
        try: ext[ext<0.0] = 0.0
        except: pass
        
        #geo['sza'] = s['sza'][i]
        geo['lat'] = s['Lat'][i]
        geo['lon'] = s['Lon'][i]
        
        for ux in np.arange(0,24,0.5):
            geo['utc'] = ux
            geo['hour'] = int(ux)
            geo['minute'] = int((ux-int(ux))*60.0)
            geo['second'] = int((ux-geo['hour']-geo['minute']/60.0)*3600.0)

            aero['ext'] = ext
            fname = 'COSR_DARE_{d}_{v}_{i:06d}_{ux:04.1f}.dat'.format(d=day,v=vv,i=i,ux=ux)
            RL.write_input_aac(fpp_in+fname,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                       verbose=False,make_base=False,set_quiet=True,solver='twostr')

            file_list.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+fname,out=fpp_out+fname))

            fnamec = 'COSR_DARE_{d}_{v}_{i:06d}_{ux:04.1f}_clean.dat'.format(d=day,v=vv,i=i,ux=ux)
            aero.pop('ext',None)
            RL.write_input_aac(fpp_in+fnamec,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                       verbose=False,make_base=False,set_quiet=True,solver='twostr')

            file_list_clean.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+fnamec,out=fpp_out+fnamec))
    pbar.update(1)
    #print '{} / {}'.format(i,nu)
file_list.close()
file_list_clean.close()
print 'done'


# In[111]:


aero_base['z_arr']


# ## Run the files

# In[112]:


fp_file_list = fp_rtm+'COSR_DARE_list_file_{d}_{v}.sh'.format(d=day,v=vv)
fp_file_list_clean = fp_rtm+'COSR_DARE_list_file_clean_{d}_{v}.sh'.format(d=day,v=vv)


# In[113]:


get_ipython().system(u' wc -l $fp_file_list')
get_ipython().system(u' wc -l $fp_file_list_clean')


# In[114]:


fp_file_list_out = fp_file_list+'.out'
fp_file_list_clean_out = fp_file_list_clean+'.out'


# In[115]:


get_ipython().system(u'parallel --jobs=22 --bar --results /scratch/output_dir/out.csv < $fp_file_list  #2> $fp_file_list_out')


# In[175]:


get_ipython().system(u'parallel --jobs=22 --bar --results /scratch/output_dir/out_cl.csv < $fp_file_list_clean ')


# ## Now read the libradtran output files

# In[176]:


out = {'ssa':[],'asy':[],'ext':[],'albedo':alb,'aod':[],'alt':[],
       'sza':[],'utc':[],'lat':[],'lon':[],'wvl':wvl,'z_aero':aero_base['z_arr']}


# In[177]:


nzout = len(geo['zout'])


# In[178]:


fl = flag & (np.isfinite(ssa_u))


# In[179]:


nl = fl.sum()


# In[180]:


nut = len(np.arange(0,24,0.5))


# In[181]:


star_aero = {'dn':np.zeros((nl,nut,nzout))+np.nan,'up':np.zeros((nl,nut,nzout))+np.nan}
star_aero_cl = {'dn':np.zeros((nl,nut,nzout))+np.nan,'up':np.zeros((nl,nut,nzout))+np.nan}
star_aero_C = np.zeros((nl,nut,nzout))+np.nan
star_aero_C_avg = np.zeros((nl,nzout))+np.nan


# In[164]:


import pixiedust


# In[182]:


#%%pixie_debugger
nu = len(s['utc'])
pbar = tqdm(total=nu)
out['utcx'] = np.arange(0,24,0.5)
print 'Reading files'
j = 0
for i,u in enumerate(s['utc']):
    aero = {}
    if flag[i] & np.isfinite(ssa_u[i]):
        aod = s['aod'][i,:]
        aero['ext'] = expand_ext_vert_and_spect(ext_,ext_z,s['aod'][i,:],s['Alt'][i]/1000.0,wvl)
        iw = np.argmin(abs(aero_base['wvl_arr']-0.5))
        aero['ssa'] = aero_base['ssa']*ssa_u[i]/aero_base['ssa'][0,iw]
        aero['asy'] = aero_base['asy']
        aero['z_arr'] = aero_base['z_arr']
        aero['wvl_arr'] = aero_base['wvl_arr']*1000.0
        
        try: aero['ssa'][aero['ssa']<0.0] = 0.0
        except: pass
        try: aero['ssa'][aero['ssa']>1.0] = 1.0
        except: pass
        try: aero['asy'][aero['asy']<0.0] = 0.0
        except: pass
        try: aero['asy'][aero['asy']>1.0] = 1.0
        except: pass
        try: aero['ext'][aero['ext']<0.0] = 0.0
        except: pass
        
        out['utc'].append(u)
        out['sza'].append(s['sza'][i])
        out['lat'].append(s['Lat'][i])
        out['lon'].append(s['Lon'][i])
        out['alt'].append(s['Alt'][i])
        out['ext'].append(aero['ext'])
        out['ssa'].append(aero['ssa'])
        out['asy'].append(aero['asy'])
        out['aod'].append(aod)
        
        
        for ui, ux in enumerate(np.arange(0,24,0.5)):
        
            try:
                fname = 'COSR_DARE_{d}_{v}_{i:06d}_{ux:04.1f}.dat'.format(d=day,v=vv,i=i,ux=ux)
                sa = RL.read_libradtran(fpp_out+fname,zout=geo['zout'])
                star_aero['dn'][j,ui,:] = sa['diffuse_down']+sa['direct_down']
                star_aero['up'][j,ui,:] = sa['diffuse_up']
         
                fnamec = 'COSR_DARE_{d}_{v}_{i:06d}_{ux:04.1f}_clean.dat'.format(d=day,v=vv,i=i,ux=ux)
                sc = RL.read_libradtran(fpp_out+fnamec,zout=geo['zout'])
                star_aero_cl['dn'][j,ui,:] = sc['diffuse_down']+sc['direct_down']
                star_aero_cl['up'][j,ui,:] = sc['diffuse_up']
            except IOError:
                star_aero['dn'][j,ui,:] = np.nan
                star_aero['up'][j,ui,:] = np.nan
                star_aero_cl['dn'][j,ui,:] = np.nan
                star_aero_cl['up'][j,ui,:] = np.nan
                pass
            
            star_aero_C[j,ui,:] = (star_aero['dn'][j,ui,:]-star_aero['up'][j,ui,:]) -                                 (star_aero_cl['dn'][j,ui,:]-star_aero_cl['up'][j,ui,:])
        star_aero_C_avg[j,:] = np.nanmean(star_aero_C[j,:,:],axis=0)
        
        j = j +1
    pbar.update(1)

out['utc'] = np.array(out['utc'])
out['sza'] = np.array(out['sza'])
out['lat'] = np.array(out['lat'])
out['lon'] = np.array(out['lon'])
out['alt'] = np.array(out['alt'])
out['ext'] = np.array(out['ext'])
out['ssa'] = np.array(out['ssa'])
out['asy'] = np.array(out['asy'])
out['aod'] = np.array(out['aod'])
out['dare'] = star_aero_C
out['dare_avg'] = star_aero_C_avg
out['dn_aero'] = star_aero['dn']
out['up_aero'] = star_aero['up']
out['dn_clear'] = star_aero_cl['dn']
out['up_clear'] = star_aero_cl['up']

print 'done'


# In[183]:


out['dare_avg'][out['dare_avg']==0] = np.nan


# In[184]:


#star = wu.iterate_dict_unicode(out)
out['zout'] = geo['zout']
print 'saving file to: '+fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d=day,vv=vv)
hs.savemat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d=day,vv=vv),out)


# # Now plot out the DARE for this flight

# In[186]:


out['dare'].shape


# In[187]:


out['aod'][:-1,3].shape


# In[188]:


out['dare_avg'].shape


# In[189]:


out['sza'][4000]


# In[190]:


out['dare'][4000,:,0]


# In[191]:


out['ext'][4230,:,3]


# In[193]:


plt.figure()
plt.scatter(out['utc'][:],out['dare_avg'][:,0],c=out['aod'][:,3])
plt.xlabel('UTC [h]')
plt.ylabel('Daily average DARE [W/m$^2$]')
ax2 = plt.gca().twinx()
ax2.plot(out['utc'],out['sza'],'.r')
ax2.set_ylabel('SZA [$^{{\\circ}}$]',color='r')
ax2.tick_params(axis='y', labelcolor='r')
cb = plt.colorbar( pad=0.15)
cb.set_label('AOD at 500 nm')

plt.title('Surface DARE from Oil Sands on {}'.format(day))
plt.savefig(fp+'plots/DARE_avg_aod_utc_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[196]:


plt.figure()
plt.scatter(out['utc'][:],out['dare'][:,24,0],c=out['aod'][:,3])
plt.xlabel('UTC [h]')
plt.ylabel('Instanenous DARE [W/m$^2$]')
ax2 = plt.gca().twinx()
ax2.plot(out['utc'],out['sza'],'.r')
ax2.set_ylabel('SZA [$^{{\\circ}}$]',color='r')
ax2.tick_params(axis='y', labelcolor='r')
cb = plt.colorbar( pad=0.15)
cb.set_label('AOD at 500 nm')

plt.title('DARE from Oil Sands on {}'.format(day))
plt.savefig(fp+'plots/DARE_aod_utc_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[197]:


plt.figure()
plt.scatter(out['ssa'][:,0,3],out['dare'][:,0],c=out['utc'])
cb = plt.colorbar()
cb.set_label('UTC')
plt.xlabel('SSA 500nm')
plt.ylabel('Instantaneous DARE [W/m$^2$]')
plt.title('DARE calculations from Oil Sands on {}'.format(day))

plt.savefig(fp+'plots/DARE_vs_SSA_{}.png'.format(day),dpi=600,transparent=True)


# ## Plot compared to UHSAS

# In[198]:


fx_h = interp1d(uh['utc'],uh['nConc'],bounds_error=False)
nConc = fx_h(out['utc'])


# In[199]:


plt.figure()
plt.scatter(nConc,out['dare '][:,0],c=out['ssa'][:,0,3])
cb = plt.colorbar()
cb.set_label('SSA 500 nm')
plt.xlabel('UHSAS number Concentration [#]')
plt.ylabel('Instantaneous DARE [W/m$^2$]')
plt.title('DARE calculations from Oil Sands on {}'.format(day))
plt.savefig(fp+'plots/DARE_vs_nConc_UHSAS_{}.png'.format(day),dpi=600,transparent=True)


# In[550]:


plt.figure()
plt.hist2d(nConc,out['dare'][:,0],range=[[0,7500],[-25,0]],bins=50)
plt.colorbar(label='Number')
plt.xlabel('UHSAS number Concentration [#]')
plt.ylabel('Instantaneous DARE [W/m$^2$]')


# ## Plot filtered for in plume

# In[1048]:


out = out1


# In[1049]:


ipl = []
dare_pl = []
alt_pl = []
dare_out = []
alt_out = []
for ii, fo in enumerate(from_utc):
    pl = (out['utc']>=fo)&(out['utc']<=to_utc[ii])
    if pl.any():
        alt_pl = np.append(alt_pl,out['alt'][pl])
        dare_pl = np.append(dare_pl,out['dare_avg'][pl,0])
        ipl.append(pl)
    dare_out = np.append(dare_out,out['dare_avg'][~pl,0])
    alt_out = np.append(alt_out,out['alt'][~pl,0])


# In[1050]:


plt.figure()
plt.hist([dare_out,dare_pl],color=['r','g'],label=['Background','In plume'],normed=True,bins=30)
plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.title('Diurnally averaged Surface DARE for {}'.format(day))


# In[688]:


ssa_u[1235]


# In[724]:


np.where(flag&np.isfinite(ssa_u))[0][367]


# In[690]:


out['ext'][367,:,3]


# In[692]:


out['z_aero']


# In[693]:


out['ext'][367,:,3].sum()


# In[695]:


out['aod'][367,3]


# In[725]:


i = 4105


# In[726]:


ext = expand_ext_vert_and_spect(ext_,ext_z,s['aod'][i,:],s['Alt'][i]/1000.0,wvl)


# In[727]:


ext.shape


# In[728]:


ext[:,3].sum()


# In[729]:


s['aod'][i,3]


# In[705]:


out1['daystr'] = '20180609'
out2['daystr'] = '20180618'
out3['daystr'] = '20180624'
out4['daystr'] = '20180625'


# In[706]:


for outi in [out1,out2,out3,out4]:
    plt.figure()
    plt.plot(outi['aod'][:,3],outi['ext'][:,:,3].sum(axis=1)/10.0,'.')
    plt.plot([0,1],[0,1],'--k')
    plt.xlabel('AOD')
    plt.ylabel('extinction')
    plt.title('{}'.format(outi['daystr']))


# In[733]:


for outi in [out1,out2,out3,out4]:
    plt.figure()
    plt.plot(outi['utc'],outi['aod'][:,3],'.',label='AOD')
    plt.plot(outi['utc'],outi['ext'][:,:,3].sum(axis=1)/10.0,'.',label='Ext')
    plt.plot(outi['utc'],outi['alt']/2000.0,'.',color='grey')
    #plt.plot([0,1],[0,1],'--k')
    plt.xlabel('UTC')
    plt.ylabel('AOD 500 nm')
    plt.legend()
    plt.title('{}'.format(outi['daystr']))


# ## Plot filtered for level legs

# In[811]:


nbox = 30


# In[812]:


std_alt1 = running_std(out1['alt'],nbox)


# In[813]:


plt.figure()
plt.plot(out1['utc'],out1['alt'],'.',label='alt')
plt.plot(out1['utc'][:1-nbox],std_alt1,'.',label='std box')
plt.legend()


# In[816]:


falt = np.where(std_alt1<100.0)[0]


# In[827]:


fig, ax1 = plt.subplots()

ax1.plot(out1['utc'],out1['alt'],'.',label='all')
ax1.plot(out1['utc'][falt],out1['alt'][falt],'s',label='filetered')
plt.legend()

ax2 = ax1.twinx()
ax2.plot(out1['utc'],out1['ssa'][:,7,3],'ok',label='ssa')
ax2.plot(out1['utc'][falt],out1['ssa'][falt,7,3],'xg',label='good ssa')
plt.legend()


# In[819]:


out1.keys()


# In[820]:


out1['ssa'].shape


# In[823]:


out1['z_aero'].shape


# # Combine the DARE calc for all days

# In[ ]:


day = '20180609'
#day = '20180618'
#day = '20180624'
#day = '20180625'


# In[999]:


out1 = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d='20180609',vv=vv))


# In[204]:


out2 = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d='20180618',vv=vv))


# In[373]:


out3 = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d='20180624',vv=vv))


# In[374]:


out4 = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d='20180625',vv=vv))


# ## Get the time tables 

# ### 20180609

# In[1069]:


flttable1 = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format('20180609'))
fromtime1 = flttable['FromTime'][flttable['FlightType']=='in plume']
totime1 = flttable['ToTime'][flttable['FlightType']=='in plume']
plumeid1 = flttable['PlumeId'][(flttable['PlumeId']=='A') | (flttable['PlumeId']=='B')].to_numpy()


# In[1064]:


def time_utc(x):
    return np.array([y.hour+y.minute/60.0+y.second/3600.0 for y in x])


# In[1065]:


from_utc1 = time_utc(fromtime1.to_numpy())
to_utc1 = time_utc(totime1.to_numpy())


# In[1072]:


ipl1 = []
dare_pl1,dare_pl1a,dare_pl1b = [],[],[]
alt_pl1 = []
dare_out1 = []
alt_out1 = []
dare_pl1_toa,dare_pl1a_toa,dare_pl1b_toa, dare_out1_toa  = [],[],[],[]
for ii, fo in enumerate(from_utc1):
    pl1 = (out1['utc']>=fo)&(out1['utc']<=to_utc1[ii])
    if pl1.any():
        alt_pl1 = np.append(alt_pl1,out1['alt'][pl1])
        dare_pl1 = np.append(dare_pl1,out1['dare_avg'][pl1,0])
        dare_pl1_toa = np.append(dare_pl1_toa,out1['dare_avg'][pl1,2])
        if plumeid1[ii]=='A':
            dare_pl1a = np.append(dare_pl1a,out1['dare_avg'][pl1,0])
            dare_pl1a_toa = np.append(dare_pl1a_toa,out1['dare_avg'][pl1,2])
        else:
            dare_pl1b = np.append(dare_pl1b,out1['dare_avg'][pl1,0])
            dare_pl1b_toa = np.append(dare_pl1b_toa,out1['dare_avg'][pl1,2])
        ipl1.append(pl1)
    dare_out1 = np.append(dare_out1,out1['dare_avg'][~pl1,0])
    dare_out1_toa = np.append(dare_out1_toa,out1['dare_avg'][~pl1,2])
    alt_out1 = np.append(alt_out1,out1['alt'][~pl1,0])


# In[1004]:


out1.keys()


# In[1005]:


ppl1 = np.array(ipl1).flatten()


# In[1010]:


# Plot out the input extinction and aod
plt.figure()
plt.plot(insitu['utc'],insitu['extCalc500nm']/10.0,'.k',label='Ext')
plt.plot(s['utc'],s['tau_aero'][:,400],'.r',label='AOD')
plt.plot(s['utc'],s['Alt']/1000.0,'+',color='grey')
plt.legend()
for ii, fo in enumerate(from_utc1):
    plsitu = (insitu['utc']>=fo)&(insitu['utc']<=to_utc1[ii])
    if plsitu.any():
        plt.plot(insitu['utc'][plsitu],insitu['extCalc500nm'][plsitu]/10.0,'og')
    plaod = (s['utc']>=fo)&(s['utc']<=to_utc1[ii])
    if plaod.any():
        plt.plot(s['utc'][plaod],s['tau_aero'][plaod,400],'ob')
    
    plaodf = (s['utc'][flag]>=fo)&(s['utc'][flag]<=to_utc1[ii])
    if plaod.any():
        plt.plot(s['utc'][flag][plaodf],s['tau_aero'][flag,400][plaodf],'xc')
    
    
        


# In[1009]:


insitu['extCalc500nm'] = np.array(insitu['extCalc500nm'])


# In[1011]:


fig = plt.figure()
ax1 = fig.add_subplot(311)
ax1.plot(out1['utc'],out1['dare_avg'][:,0],'.k',label='all')
for ipp in ipl1:
    ax1.plot(out1['utc'][ipp],out1['dare_avg'][ipp,0],'o')
ax1.set_ylabel('DARE')
    
ax2 = fig.add_subplot(312,sharex = ax1)
ax2.plot(out1['utc'],out1['aod'][:,0],'.k',label='all')
for ipp in ipl1:
    ax2.plot(out1['utc'][ipp],out1['aod'][ipp,0],'o')
ax2.set_ylabel('AOD')
    
ax3 = fig.add_subplot(313,sharex = ax1)
ax3.plot(out1['utc'],out1['ssa'][:,2,3],'.k',label='all')
for ipp in ipl1:
    ax3.plot(out1['utc'][ipp],out1['ssa'][ipp,2,3],'o')
ax3.set_ylabel('SSA')


# In[677]:


out1['ext'].shape


# In[1073]:


plt.figure()
plt.hist([dare_out1,dare_pl1a,dare_pl1b],color=['r','g','b'],label=['Background','In plume A','In plume B'],
         normed=True,bins=30)
plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.title('Diurnally averaged Surface DARE for {}'.format(day))


# In[1074]:


plt.figure()
plt.hist([dare_out1_toa,dare_pl1a_toa,dare_pl1b_toa],color=['r','g','b'],label=['Background','In plume A','In plume B'],
         normed=True,bins=30)
plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.title('Diurnally averaged TOA DARE for {}'.format(day))


# ### 20180618

# In[200]:


flttable2 = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format('20180618'))
fromtime2 = flttable2['FromTime'][flttable2['FlightType']=='in plume']
totime2 = flttable2['ToTime'][flttable2['FlightType']=='in plume']


# In[209]:


plumeid2 = flttable2['PlumeId'][(flttable2['PlumeId']=='A') | (flttable2['PlumeId']=='B')].to_numpy()


# In[201]:


flttable2


# In[202]:


from_utc2 = time_utc(fromtime2.to_numpy())
to_utc2 = time_utc(totime2.to_numpy())


# In[210]:


ipl2 = []
dare_pl2 = []
dare_pl2a,dare_pl2b = [], []
alt_pl2 = []
dare_out2 = []
alt_out2 = []
dare_pl2_toa,dare_pl2a_toa,dare_pl2b_toa, dare_out2_toa  = [],[],[],[]
for ii, fo in enumerate(from_utc2):
    pl2 = (out2['utc']>=fo)&(out2['utc']<=to_utc2[ii])
    if pl2.any():
        alt_pl2 = np.append(alt_pl2,out2['alt'][pl2])
        dare_pl2 = np.append(dare_pl2,out2['dare_avg'][pl2[:],0])
        dare_pl2_toa = np.append(dare_pl2_toa,out2['dare_avg'][pl2[:],2])
        if plumeid2[ii]=='A':
            dare_pl2a = np.append(dare_pl2a,out2['dare_avg'][pl2,0])
            dare_pl2a_toa = np.append(dare_pl2a_toa,out2['dare_avg'][pl2,2])
        else:
            dare_pl2b = np.append(dare_pl2b,out2['dare_avg'][pl2,0])
            dare_pl2b_toa = np.append(dare_pl2b_toa,out2['dare_avg'][pl2,2])
        ipl2.append(pl2)
    dare_out2 = np.append(dare_out2,out2['dare_avg'][~pl2[:],0])
    dare_out2_toa = np.append(dare_out2_toa,out2['dare_avg'][~pl2[:],2])
    alt_out2 = np.append(alt_out2,out2['alt'][~pl2,0])


# In[211]:


plt.figure()
plt.hist([dare_out2,dare_pl2a,dare_pl2b],color=['r','g','b'],label=['Background','In plume A','In plume B'],
         normed=True,bins=30)
plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.title('Diurnally averaged Surface DARE for {}'.format(day))


# ### 20180624

# In[387]:


flttable3 = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format('20180624'))
fromtime3 = flttable3['FromTime'][flttable3['FlightType']=='in plume']
totime3 = flttable3['ToTime'][flttable3['FlightType']=='in plume']


# In[388]:


from_utc3 = time_utc(fromtime3.to_numpy())
to_utc3 = time_utc(totime3.to_numpy())


# In[418]:


ipl3 = []
dare_pl3 = []
alt_pl3 = []
dare_out3 = []
alt_out3 = []
dare_pl3_toa, dare_out3_toa  = [],[]
for ii, fo in enumerate(from_utc3):
    pl3 = (out3['utc']>=fo)&(out3['utc']<=to_utc3[ii])
    if pl3.any():
        alt_pl3 = np.append(alt_pl3,out3['alt'][pl3])
        dare_pl3 = np.append(dare_pl3,out3['dare_avg'][pl3,0])
        dare_pl3_toa = np.append(dare_pl3_toa,out3['dare_avg'][pl3,2])
        ipl3.append(pl3)
    dare_out3 = np.append(dare_out3,out3['dare_avg'][~pl3,0])
    dare_out3_toa = np.append(dare_out3_toa,out3['dare_avg'][~pl3,0])
    alt_out3 = np.append(alt_out3,out3['alt'][~pl3,0])


# ### 20180625

# In[391]:


flttable4 = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format('20180625'))
fromtime4 = flttable4['FromTime'][flttable4['FlightType']=='in plume']
totime4 = flttable4['ToTime'][flttable4['FlightType']=='in plume']


# In[392]:


from_utc4 = time_utc(fromtime4.to_numpy())
to_utc4 = time_utc(totime4.to_numpy())


# In[419]:


ipl4 = []
dare_pl4 = []
alt_pl4 = []
dare_out4 = []
alt_out4 = []
dare_pl4_toa, dare_out4_toa  = [],[]
for ii, fo in enumerate(from_utc4):
    pl4 = (out4['utc']>=fo)&(out4['utc']<=to_utc4[ii])
    if pl4.any():
        alt_pl4 = np.append(alt_pl4,out4['alt'][pl4])
        dare_pl4 = np.append(dare_pl4,out4['dare_avg'][pl4,0])
        dare_pl4_toa = np.append(dare_pl4_toa,out4['dare_avg'][pl4,0])
        ipl4.append(pl4)
    dare_out4 = np.append(dare_out4,out4['dare_avg'][~pl4,0])
    dare_out4_toa = np.append(dare_out4_toa,out4['dare_avg'][~pl4,0])
    alt_out4 = np.append(alt_out4,out4['alt'][~pl4,0])


# ## Plot the DAREs

# In[412]:


fp


# In[414]:


plt.figure()
plt.hist([dare_pl1,dare_pl2,dare_pl3,dare_pl4],color=['b','r','gold','indigo'],
         label=['In Pume 20180609','In Pume 20180618','In Pume 20180624','In Pume 20180625'],normed=True,bins=30,alpha=0.6)
plt.hist([dare_out1,dare_out2,dare_out3,dare_out4],color=['lightblue','coral','yellow','pink'],
         label=['Background 20180609','Background 20180618','Background 20180624','Background 20180625'],
         normed=True,bins=30,alpha=0.6)
plt.xlim(-150,0)
plt.yscale('log')
plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.title('Diurnally averaged Surface DARE')
plt.savefig(fp+'plots/DARE_avg_inplume_out_COSR.png',transparent=True,dpi=600)


# In[421]:


plt.figure()
plt.hist([dare_pl1_toa,dare_pl2_toa,dare_pl3_toa,dare_pl4_toa],color=['b','r','gold','indigo'],
         label=['In Pume 20180609','In Pume 20180618','In Pume 20180624','In Pume 20180625'],normed=True,bins=30,alpha=0.6)
plt.hist([dare_out1_toa,dare_out2_toa,dare_out3_toa,dare_out4_toa],color=['lightblue','coral','yellow','pink'],
         label=['Background 20180609','Background 20180618','Background 20180624','Background 20180625'],
         normed=True,bins=30,alpha=0.6)
plt.xlim(-150,0)
plt.yscale('log')
plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.title('Diurnally averaged TOA DARE')
plt.savefig(fp+'plots/DARE_avg_inplume_out_COSR_TOA.png',transparent=True,dpi=600)


# In[406]:


np.nanmean(dare_pl1),np.nanmean(dare_pl2),np.nanmean(dare_pl3),np.nanmean(dare_pl4)


# In[407]:


np.nanmedian(dare_pl1),np.nanmedian(dare_pl2),np.nanmedian(dare_pl3),np.nanmedian(dare_pl4)


# In[408]:


np.nanstd(dare_pl1),np.nanstd(dare_pl2),np.nanstd(dare_pl3),np.nanstd(dare_pl4)


# In[409]:


np.nanmean(dare_out1),np.nanmean(dare_out2),np.nanmean(dare_out3),np.nanmean(dare_out4)


# In[410]:


np.nanmedian(dare_out1),np.nanmedian(dare_out2),np.nanmedian(dare_out3),np.nanmedian(dare_out4)


# In[411]:


np.nanstd(dare_out1),np.nanstd(dare_out2),np.nanstd(dare_out3),np.nanstd(dare_out4)


# In[ ]:




