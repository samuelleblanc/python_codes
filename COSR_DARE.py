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


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


import Run_libradtran as RL
import os
import write_utils as wu
from tqdm.notebook import tqdm


# In[4]:


name = 'COSR'
vv = 'v3b'


# In[5]:


fp =getpath('COSR')


# In[6]:


day = '20180609'
#day = '20180618'
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

# ## Load the 4STAR AOD

# In[22]:


s = sio.loadmat(fp+'os_data/4STAR_{}starsun.mat'.format(day))


# In[23]:


s.keys()


# In[24]:


s['utc'] = lu.toutc(lu.mat2py_time(s['t']))


# ### Load the flag files

# In[25]:


fmat = getpath('4STAR_data')


# In[26]:


with open (fmat+'starinfo_{}.m'.format(day), 'rt') as in_file:
    for line in in_file:
        if 'flagfilename ' in line:
            ff = line.split("'")[1]
sf = hs.loadmat(fmat+ff)


# In[27]:


sf.keys()


# In[28]:


flag = sf['manual_flags']['good'][0,:,0]


# In[29]:


flag.shape


# In[30]:


sum(flag)


# ### use the polyfit aod on the wavelength array

# In[11]:


s['tau_aero_polynomial'].shape


# In[12]:


wvl = np.array([0.25,0.35,0.4,0.5,0.675,0.87,0.995,1.2,1.4,1.6,2.1,3.2,4.9])


# In[13]:


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

# In[38]:


wvl = np.array([0.25,0.35,0.4,0.5,0.675,0.87,0.995,1.2,1.4,1.6,2.1,3.2,4.9])


# In[37]:


sai = s['aerosolcols'][0,:].astype(int)


# In[298]:


plt.figure()
plt.plot(s['tau_aero'][:,400],'.')


# In[296]:


i=13884


# In[35]:


s['tau_aero'].shape


# In[36]:


s['w'].shape


# In[300]:


plt.figure()
plt.plot(s['w'][0,:],s['tau_aero'][i,:],'-')
plt.plot(s['w'][0,sai],s['tau_aero'][i,sai],'.')
#plt.plot(wvl,s['aod'][i,:],'x-')
plt.ylim(0,0.5)


# In[31]:


np.where(flag)[0]


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


# In[32]:


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


# In[33]:


s['w'][0,saii]


# In[302]:


plt.figure()
plt.plot(s['w'][0,:],s['tau_aero'][i,:],'-')
plt.plot(s['w'][0,saii],s['tau_aero'][i,saii],'.')
plt.ylim(0,0.5)


# In[34]:


saii = sai[((s['w'][0,sai]>0.364) & (s['w'][0,sai]<0.370)) | ((s['w'][0,sai]>0.435) & (s['w'][0,sai]<1.566)) | 
            ((s['w'][0,sai]>1.584) & (s['w'][0,sai]<1.597)) ]


# In[304]:


pl = su.logaod_polyfit(np.append(s['w'][0,saii],[2.2,2.7,3.5]),np.append(s['tau_aero'][i,saii],[s['tau_aero'][i,saii][-1]*2.0/3.0,s['tau_aero'][i,saii][-1]/2.0,s['tau_aero'][i,saii][-1]/4.0]),polynum=2)


# In[305]:


plt.figure()
plt.plot(s['w'][0,:],s['tau_aero'][i,:],'-')
plt.plot(s['w'][0,sai],s['tau_aero'][i,sai],'.')
#plt.plot(wvl,s['aod'][i,:],'x-')
plt.plot(wvl,np.exp(np.polyval(pl,np.log(wvl))))

plt.yscale('log')
plt.xscale('log')
plt.ylim(0.01,2)
plt.xlim(0.3,5)


# In[306]:


from scipy import polyfit


# In[307]:


tau_aero_good = np.array([np.append(s['tau_aero_subtract_all'][i,saii],[s['tau_aero_subtract_all'][i,saii][-1]*2.0/3.0,
                                                           s['tau_aero_subtract_all'][i,saii][-1]/2.0,s['tau_aero_subtract_all'][i,saii][-1]/10.0]) \
                          for i in xrange(len(s['t']))])


# In[308]:


poly = np.array([polyfit(np.log(np.append(s['w'][0,saii],[2.2,2.7,3.7])),np.log(aodd),2) for aodd in tau_aero_good])


# In[309]:


poly.shape


# In[310]:


s['paod'] = np.zeros((len(s['utc']),len(wvl)))
for i in xrange(len(s['utc'])):
    s['paod'][i,:] = np.exp(np.polyval(poly[i,:],np.log(wvl)))


# In[311]:


plt.figure()
plt.plot(wvl,s['paod'][flag,:].T[:,10:-1:500])
plt.plot(s['w'][0,saii],s['tau_aero'][:,saii][flag,:].T[:,10:-1:500],'.')
plt.ylim(0,1)
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.01,10)
plt.xlim(0.25,5)


# In[312]:


s['aod'] = s['paod']


# In[313]:


i=3778


# In[314]:


s['aod'][i,3]


# In[315]:


s['Alt'][i]


# ## Load the in situ files

# In[316]:


situ = pd.read_csv(fp+'data_other/{}_nephclap.csv'.format(day))


# In[ ]:


situ


# In[317]:


insitu = situ.to_dict('list')


# In[318]:


insitu.keys()


# In[319]:


insitu['ssa_500nm'] = np.array(insitu['totScatCalc_500nm'])/np.array(insitu['extCalc500nm'])


# In[320]:


insitu['extCalc500nm'] = np.array(insitu['extCalc500nm'])


# In[321]:


gu = pd.to_datetime(situ['DateTimeUTC']).to_list()
insitu['utc'] = np.array([g.hour+g.minute/60.0+g.second/3600.0 for g in gu])


# In[322]:


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


# In[325]:


plt.figure(figsize=(7,2.5))
plt.plot(insitu['utc'],insitu['ssa_500nm'],'.',alpha=0.3,label='calculated from CLAP+neph at 500 nm')
plt.plot(insitu['utc'],Sp.smooth(insitu['ssa_500nm'],60,old=True),'-r',label='smooth over 60s')
plt.legend()
plt.xlabel('UTC [h]')
plt.ylabel('SSA at 500 nm')
#plt.title('SSA from in situ calculated on: {}'.format(day))
plt.ylim(0.45,1.0)
#plt.xlim(15.4,20.3)
#plt.xlim(14.9,19.3)
plt.xlim(14.7,23.5)
plt.grid()
plt.tight_layout()
plt.savefig(fp+'plots/SSA_500_CLAP_neph_smooth_transp_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[326]:


ssa_insitu = Sp.smooth(insitu['ssa_500nm'],60,old=True)


# In[327]:


len(ssa_insitu)


# ### Apply Altitude change filter

# Need to load the 4STAR spectra, then run the next cell first

# In[328]:


f_alt = interp1d(x=s['utc'],y=s['Alt'][:,0])
insitu['alt'] = f_alt(insitu['utc'])


# In[412]:


plt.figure()
plt.plot(insitu['utc'],Sp.smooth(insitu['ssa_500nm']*2000.0,60,old=True),'-r')
plt.plot(insitu['utc'],insitu['totScatCalc_500nm'],'.g')
plt.plot(insitu['utc'],insitu['totAbsCalc_500nm'],'xy')
plt.plot(insitu['utc'],insitu['alt'])


# In[330]:


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


# In[331]:


nbox=30


# In[332]:


std_alt = running_std(insitu['alt'],nbox)


# In[333]:


infl = np.where(std_alt<50)[0]
infl_not = np.where(std_alt>50)[0]


# In[334]:


insitu['ssa_500nm'][infl].shape


# In[335]:


insitu['ssa_500nm'].shape


# In[336]:


ssa_in = insitu['ssa_500nm'][:]*2.0
ssa_in[infl_not] = np.nan
ssa_in = ssa_in/2.0


# In[337]:


ssa_in.shape


# In[338]:


#xout = np.convolve(x, np.ones(window)/window, 'same')
window = 15
xmasked, mask = Sp.nanmasked(ssa_in)
fx = interpolate.interp1d(insitu['utc'][mask],xmasked,bounds_error=False)
xinterp = fx(insitu['utc'])
ssaout = np.convolve(xinterp, np.ones(window)/window, 'same')

ssaout[~mask] = np.nan

#ssaout = np.convolve(ssa_in,np.ones(30)/30,'same')


# In[339]:


fig,ax1 = plt.subplots()

ax1.plot(insitu['utc'],insitu['alt'],label='alt all')
ax1.plot(insitu['utc'][infl],insitu['alt'][infl],'ob',label='good alt')

ax2 = ax1.twinx()
ax2.plot(insitu['utc'],insitu['ssa_500nm'],'--k',label='ssa all')
ax2.plot(insitu['utc'][infl],insitu['ssa_500nm'][infl],'xg',label='good ssa')
ax2.plot(insitu['utc'],ssa_in,'.y',label='presmooth')
ax2.plot(insitu['utc'],ssaout,'sc',label='smooth ssa')

plt.legend()


# In[340]:


plt.figure(figsize=(7,2.5))
plt.plot(insitu['utc'],insitu['ssa_500nm'],'.',color='lightgrey',alpha=0.3,label='calculated from CLAP+neph at 500 nm')
plt.plot(insitu['utc'][infl],insitu['ssa_500nm'][infl],'.',alpha=0.3,label='Filtered for altitude changes')
plt.plot(insitu['utc'],ssaout,'-',color='tab:red',label='smooth over 15s')
plt.legend()
plt.xlabel('UTC [h]')
plt.ylabel('SSA at 500 nm')
#plt.title('SSA from in situ calculated on: {}'.format(day))
plt.ylim(0.45,1.0)
#plt.xlim(15.0,18.5)
#plt.xlim(14.9,19.3)

plt.xlim(14.7,23.5)
plt.grid()
plt.tight_layout()
plt.savefig(fp+'plots/SSA_500_CLAP_neph_smooth_transp_filtered_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[341]:


ssa_insitu = ssaout


# In[342]:


plt.figure()
plt.plot(ssa_insitu,insitu['alt'],'.')


# ## Load the UHSAS files

# In[71]:


uh = sio.loadmat(fp+'data_other/{}_UHSAS_ECCC.mat'.format(day))


# In[72]:


uh.keys()


# In[73]:


uh['utc'] = lu.toutc(lu.mat2py_time(uh['t']))


# In[74]:


uh['utc'].shape


# In[75]:


uh['binDataConc'].sum(axis=1).shape


# In[76]:


uh['nConc'] = uh['binDataConc'].sum(axis=1)


# In[77]:


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

# In[78]:


fp


# In[79]:


flttable = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format(day))


# In[ ]:


flttable


# In[80]:


fromtime = flttable['FromTime'][flttable['FlightType']=='in plume']


# In[81]:


totime = flttable['ToTime'][flttable['FlightType']=='in plume']


# In[82]:


tt = fromtime.to_numpy()


# In[83]:


tt[0].second


# In[84]:


def time_utc(x):
    return np.array([y.hour+y.minute/60.0+y.second/3600.0 for y in x])


# In[85]:


from_utc = time_utc(fromtime.to_numpy())
to_utc = time_utc(totime.to_numpy())


# ## Load the solar spectrum

# In[443]:


fp_librad


# In[444]:


sol = np.genfromtxt(fp_librad+'solar_flux/kurudz_1.0nm.dat')


# In[445]:


sol.shape


# In[446]:


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


# In[447]:


plt.figure()
plt.plot(sol[:,0],sol[:,1].cumsum()/sol[:,1].sum()*100.0,'.')
plt.ylabel('Cumulative Irradiance [\%]')
plt.title('Kurudz solar spectrum [1nm]')
plt.xlabel('Wavelength [nm]')
plt.grid()


# In[448]:


# relative forcing efficiency to focing efficiency from fresh aerosol milagro LeBlanc 2012 / Schmidt 2010.
wvlm = np.array([400,500,600,680,780,850,1000,2500])
effm = np.array([60.0,40.0,22.0,20.0,16.0,11.0,10.0,1.0])/100.0
feffm = interp1d(wvlm,effm,bounds_error=False,fill_value='extrapolate',kind='linear')
reffm = feffm(sol[:,0])


# In[453]:


plt.figure()
plt.plot(sol[:,0],(sol[:,1]*reffm).cumsum()/sol[:,1].sum()*100.0,'.')
plt.ylabel('Cumulative Irradiance [\%]')
plt.title('Reduction of the Kurudz solar spectrum due to MILAGRO fresh aerosol')
plt.xlabel('Wavelength [nm]')
plt.grid()


# In[450]:


(sol[:,1]*reffm).sum()/sol[:,1].sum()*100.0


# ## Load the skyscan results

# In[343]:


sk_names = {
    '20180624':'4STAR_20180624_135_SKYP.created_20190329_003621.ppl_lv15.mat',
    '20180625':'4STAR_20180625_026_SKYP.created_20190507_213718.ppl_lv15.mat',
    '20180620':'4STAR_20180620_017_SKYP.created_20190507_225712.ppl_lv15.mat',
    '20180618':'4STAR_20180618_029_SKYA.created_20190507_232752.avg_lv10.mat',
    '20180609':'4STAR_20180609_041_SKYP.created_20190508_003116.ppl_lv15.mat',
    '20180705':'4STAR_20180705_061_SKYP.created_20190508_003930.ppl_lv15.mat'}


# In[344]:


sk_n = {
    '20180624':'135',
    '20180625':'026',
    '20180620':'017',
    '20180618':'029',
    '20180609':'041',
    '20180705':'061'}


# In[345]:


fp_name = sk_names[day]#'4STAR_20180624_135_SKYP.created_20190329_003621.ppl_lv15.mat'


# In[346]:


sky = sio.loadmat(fp+fp_name)
sky.keys()


# In[347]:


sky['refractive_index_real_r']


# In[348]:


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


# In[135]:


sky['g_tot'][-1]


# ## Expand the sky scans results to longer wavelengths

# In[349]:


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


# In[213]:


# for 20180624 only
if day=='20180624':
    f_asy = interp1d(np.append(sky['Wavelength'][:,0],[1.1,2.4]),np.append(sky['g_tot'][0,:],[sky['g_tot'][0,-1]-0.01,sky['g_tot'][0,-1]-0.13]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    asy = f_asy(wvl)
    f_ssa = interp1d(np.append(sky['Wavelength'][:,0],[1.1,1.2,1.3]),np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.012,sky['ssa_total'][0,-1]-0.022,sky['ssa_total'][0,-1]-0.034]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    ssa = f_ssa(wvl)


# In[350]:


if day=='20180625':
    f_asy = interp1d(np.append(sky['Wavelength'][:,0],[1.1,2.4]),np.append(sky['g_tot'][0,:],[sky['g_tot'][0,-1]-0.005,sky['g_tot'][0,-1]-0.008]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    asy = f_asy(wvl)
    f_ssa = interp1d(np.append(sky['Wavelength'][:,0],[1.1,1.2,1.3]),np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.012,sky['ssa_total'][0,-1]-0.022,sky['ssa_total'][0,-1]-0.034]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    ssa = f_ssa(wvl)


# In[137]:


if day=='20180609':
    f_asy = interp1d(np.append(sky['Wavelength'][:,0],[1.1,2.4]),np.append(sky['g_tot'][0,:],[sky['g_tot'][0,-1]-0.002,sky['g_tot'][0,-1]-0.005]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    asy = f_asy(wvl)
    f_ssa = interp1d(np.append(sky['Wavelength'][:,0],[1.1,1.2,1.3]),np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.008,sky['ssa_total'][0,-1]-0.016,sky['ssa_total'][0,-1]-0.020]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
    ssa = f_ssa(wvl)


# In[351]:


np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.008,sky['ssa_total'][0,-1]-0.002])


# In[352]:


sky['Wavelength'][:,0]


# In[353]:


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

# In[354]:


f_alt = interp1d(x=s['utc'],y=s['Alt'][:,0])
insitu['alt'] = f_alt(insitu['utc'])


# In[292]:


if day=='20180624':
    insitu['extCalc500nm'] = insitu['extCalc500nm']*10.0


# In[258]:


if day=='20180609':
    insitu['extCalc500nm'] = insitu['extCalc500nm']*10.0


# In[355]:


plt.figure()
plt.plot(insitu['extCalc500nm'][infl],insitu['alt'][infl],'.')


# In[90]:


plt.figure()
plt.plot(insitu['utc'],insitu['alt'])


# In[356]:


plt.figure()
plt.scatter(insitu['utc'],insitu['alt'],c=insitu['extCalc500nm'],cmap=plt.cm.magma,s=(insitu['extCalc500nm'])**1.5-20.0)
cb = plt.colorbar()
cb.set_label('Extinction at 500 nm [1/Mm]')
plt.xlabel('UTC [h]')
plt.ylabel('Altitude [m]')
plt.title('In Situ Extinction coefficients on {}'.format(day))
plt.savefig(fp+'plots/Extinction_UTC_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[357]:


insitu['utc'].shape,insitu['alt'].shape,insitu['extCalc500nm'].shape


# In[358]:


np.isfinite(insitu['extCalc500nm'])


# In[359]:


dz_bin = 200.0 #in meters


# In[360]:


binned_ext,binned_alt,binned_num = [],[],[]
for i in xrange(15):
    flaa = (insitu['alt'][infl]>=i*dz_bin) & (insitu['alt'][infl]<(i+1.0)*dz_bin) & (np.isfinite(insitu['extCalc500nm'][infl]))
    if flaa.any():
        binned_ext.append(insitu['extCalc500nm'][infl][flaa])
        binned_alt.append(np.mean([i*dz_bin,(i+1.0)*dz_bin]))
        binned_num.append(len(insitu['extCalc500nm'][infl][flaa]))


# In[361]:


plt.figure(figsize=(7,3))
bp =plt.boxplot(binned_ext,positions=binned_alt,vert=False,showfliers=True,widths=100,showmeans=True,patch_artist=True)
plt.xlabel('Extinction Calculated at 500 nm [1/Mm]')
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
    b.set_markeredgecolor('None')
    b.set_markerfacecolor('lightgrey')
    b.set_alpha(0.5)
ext_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(ext_means[:,0],ext_means[:,1],'-k',label='Means')
plt.xlim(0,200)
plt.ylim(0,2500)

for j,nn in enumerate(binned_num): 
    if nn>0:
        plt.text(min(bp['means'][j].get_data()[0])+5,binned_alt[j],'{:2.0f}'.format(nn),
                 color='g',fontsize=7,verticalalignment='center',horizontalalignment='left')


plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0],bp['fliers'][0]],
           ['Mean','Median','25% - 75%','min-max','outliers'],
           frameon=False,loc=1,numpoints=1)
plt.title('In situ calculated extinction CLAP+neph: {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plots/extinction_vertical_bins_clap_neph_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[222]:


ext_z = np.arange(binned_alt[0],binned_alt[-1]+2*dz_bin,dz_bin)/1000.0 # in km
fext_ = interp1d(ext_means[:,1]/1000.0,ext_means[:,0]/1000.0,fill_value=0.0,bounds_error=False)
ext_ = fext_(ext_z) # in 1/km


# In[362]:


nz = len(ext_z)


# In[363]:


ext_z


# In[364]:


ext_


# In[365]:


dz_km = ext_z[1]-ext_z[0]


# In[366]:


dz_km = ext_z[1]-ext_z[0]
np.sum(ext_*dz_km)


# ##  Calculate the extinction coefficient from 4STAR
# Abandoned inquiry for now, did not seem feasible - no full column profile within center of plume

# In[228]:


s['w'][0,405]


# In[ ]:


plt.scatter(insitu['utc'],insitu['alt'],c=insitu['extCalc500nm'],cmap=plt.cm.magma,s=(insitu['extCalc500nm'])**1.5-20.0)


# In[341]:


s['Lat'][flag].shape,s['Alt'][flag].shape,s['tau_aero'][flag,405].shape


# In[351]:


plt.figure()
plt.scatter(s['Lat'][fl_ea2][:,0],s['Alt'][fl_ea2][:,0],c=s['tau_aero'][fl_ea2,405][:],cmap=plt.cm.plasma,s=10)


# In[345]:


plt.figure()
plt.plot(s['utc'][flag],s['tau_aero'][flag,405],'.')
plt.plot(s['utc'][flag],s['Alt'][flag]/1000.0,'.')


# In[347]:


plt.figure()
plt.plot(s['utc'],s['Lat'],'.')


# In[325]:


fl_ea = flag & (s['utc']>15.50) & (s['utc']<17.00)


# In[350]:


fl_ea2 = flag & (s['utc']>16.2) & (s['utc']<18.00)


# In[326]:


plt.figure()
plt.plot(s['tau_aero'][fl_ea,405],s['Alt'][fl_ea],'.')


# # Now build the input files for DARE calculations
# 

# In[367]:


from write_utils import nearest_neighbor


# In[368]:


geo = {'zout':[0,3,100],'year':int(day[0:4]),'month':int(day[4:6]),'day':int(day[6:8]),'hour':12,'minute':0,'second':0}
aero_base = {'z_arr':(ext_z+dz_km/2.0),'wvl_arr':wvl,'ssa':np.array([ssa,]*nz),'asy':np.array([asy,]*nz)}
source = {'integrate_values':True,'dat_path':fp_librad,
          'run_fuliou':True,'wvl_range':[350,4000]}
albedo = {'create_albedo_file':True,'alb':alb,'alb_wvl':wvl*1000.0}
cloud = {'dummy':None}


# In[369]:


# Old way with no altitude filtering
if False:
    fx_ssa_insitu = interp1d(insitu['utc'],ssa_insitu,bounds_error=False)
    ssa_u = fx_ssa_insitu(s['utc'])
    ssa_u[ssa_u<0.8] = np.nan
else:
    ssa_u = nearest_neighbor(insitu['utc'],ssa_insitu,s['utc'],dist=1/3600.0)


# In[370]:


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
    dz_km = ext_z[1]-ext_z[0]
    aod_ = np.sum(ext_[iz:]*dz_km)
    #aod_ = ext_[iz:].sum()/10.0
    
    factor = aod_sp[iw]/aod_
    exts = np.array([aod_sp/aod_sp[iw]*factor*e for e in ext_])
    return exts 


# In[371]:


geo.pop('sza')


# ## Write the input files

# In[372]:


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


# ## Run the files

# In[373]:


fp_file_list = fp_rtm+'COSR_DARE_list_file_{d}_{v}.sh'.format(d=day,v=vv)
fp_file_list_clean = fp_rtm+'COSR_DARE_list_file_clean_{d}_{v}.sh'.format(d=day,v=vv)


# In[374]:


get_ipython().system(u' wc -l $fp_file_list')
get_ipython().system(u' wc -l $fp_file_list_clean')


# In[375]:


fp_file_list_out = fp_file_list+'.out'
fp_file_list_clean_out = fp_file_list_clean+'.out'


# In[ ]:


get_ipython().system(u'parallel --jobs=22 --bar --results /scratch/output_dir/out.csv < $fp_file_list  #2> $fp_file_list_out')


# In[ ]:


get_ipython().system(u'parallel --jobs=22 --bar --results /scratch/output_dir/out_cl.csv < $fp_file_list_clean ')


# ## Now read the libradtran output files

# In[378]:


out = {'ssa':[],'asy':[],'ext':[],'albedo':alb,'aod':[],'alt':[],
       'sza':[],'utc':[],'lat':[],'lon':[],'wvl':wvl,'z_aero':aero_base['z_arr']}


# In[379]:


nzout = len(geo['zout'])


# In[380]:


fl = flag & (np.isfinite(ssa_u))


# In[381]:


nl = fl.sum()


# In[382]:


nut = len(np.arange(0,24,0.5))


# In[383]:


star_aero = {'dn':np.zeros((nl,nut,nzout))+np.nan,'up':np.zeros((nl,nut,nzout))+np.nan}
star_aero_cl = {'dn':np.zeros((nl,nut,nzout))+np.nan,'up':np.zeros((nl,nut,nzout))+np.nan}
star_aero_C = np.zeros((nl,nut,nzout))+np.nan
star_aero_C_avg = np.zeros((nl,nzout))+np.nan


# In[ ]:


import pixiedust


# In[255]:


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


# ## Filter for bad data

# ### For too high data

# In[258]:


itoohigh = (out['alt'] >= ext_z[-3]*1000.0)[:,0]


# In[259]:


itoohigh


# In[260]:


out['dare'][itoohigh,:,:] = np.nan
out['dare_avg'][itoohigh,:] = np.nan


# ### For zeros

# In[261]:


out['dare_avg'][out['dare_avg']==0] = np.nan


# ## Save DARE to file

# In[262]:


#star = wu.iterate_dict_unicode(out)
out['zout'] = geo['zout']
print 'saving file to: '+fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d=day,vv=vv)
hs.savemat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d=day,vv=vv),out)


# In[12]:


out = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d=day,vv=vv))


# # Now plot out the DARE for this flight

# In[263]:


out['dare'].shape


# In[264]:


out['aod'][:-1,3].shape


# In[265]:


out['dare_avg'].shape


# ## Plot selection of input/output variables

# In[266]:


plt.figure()
plt.plot(out['utc'],out['dn_aero'][:,24,0],'.',label='dn_aero')
plt.plot(out['utc'],out['up_aero'][:,24,0],'.',label='up_aero')
plt.plot(out['utc'],out['dn_clear'][:,24,0],'.',label='dn_clear')
plt.plot(out['utc'],out['up_clear'][:,24,0],'.',label='up_clear')
plt.plot(out['utc'],out['aod'][:,3]*100.0,'.',label='aod')
plt.plot(out['utc'],out['asy'][:,6,3]*100.0,'x',label='asy')
plt.plot(out['utc'],out['ssa'][:,6,3]*100.0,'o',label='ssa')
plt.plot(out['utc'],out['ext'][:,6,3]*100,'s',label='ext')
plt.plot(out['utc'],out['alt']/10.0,'-',label='Altitude')
leg = plt.legend()
leg.draggable()


# In[304]:


plt.figure()

plt.plot(out['ext'][2730,:,3],ext_z,'.',label='utc: {:2.2f}'.format(out['utc'][2730]))
plt.plot(out['aod'][2730,3],[1.0],'s',label='AOD utc: {:2.2f}'.format(out['utc'][2730]))
plt.plot(out['ext'][2050,:,3],ext_z,'.',label='utc: {:2.2f}'.format(out['utc'][2050]))
plt.plot(out['aod'][2050,3],[1.0],'s',label='AOD utc: {:2.2f}'.format(out['utc'][2050]))
plt.plot(out['ext'][2750,:,3],ext_z,'.',label='utc: {:2.2f}'.format(out['utc'][2750]))
plt.plot(out['aod'][2750,3],[1.0],'s',label='AOD utc: {:2.2f}'.format(out['utc'][2750]))

plt.legend()
plt.xlabel('Extinction [km$^{{-1}}$]')
plt.ylabel('Altitude [km]')


# ## Plots the DARE

# In[267]:


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


# In[131]:


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


# In[268]:


plt.figure()
plt.scatter(out['ssa'][:,0,3],out['dare_avg'][:,0],c=out['aod'][:,3])
cb = plt.colorbar()
cb.set_label('AOD 500 nm')
plt.xlabel('SSA 500nm')
plt.ylabel('Instantaneous DARE [W/m$^2$]')
plt.title('DARE calculations from Oil Sands on {}'.format(day))

plt.savefig(fp+'plots/DARE_24havg_vs_SSA_{}_{}.png'.format(vv,day),dpi=600,transparent=True)


# ## Plot compared to UHSAS

# In[337]:


fx_h = interp1d(uh['utc'],uh['nConc'],bounds_error=False)
nConc = fx_h(out['utc'])


# In[338]:


plt.figure()
plt.scatter(nConc,out['dare_avg'][:,0],c=out['ssa'][:,0,3])
cb = plt.colorbar()
cb.set_label('SSA 500 nm')
plt.xlabel('UHSAS number Concentration [#]')
plt.ylabel('Instantaneous DARE [W/m$^2$]')
plt.title('DARE calculations from Oil Sands on {}'.format(day))
plt.savefig(fp+'plots/DARE_vs_nConc_UHSAS_{}.png'.format(day),dpi=600,transparent=True)


# In[339]:


plt.figure()
plt.hist2d(nConc,out['dare_avg'][:,0],range=[[0,7500],[-100,0]],bins=50)
plt.colorbar(label='Number')
plt.xlabel('UHSAS number Concentration [#]')
plt.ylabel('Instantaneous DARE [W/m$^2$]')


# ## Plot filtered for in plume

# In[341]:


out = out1


# In[342]:


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


# In[343]:


plt.figure()
plt.hist([dare_out,dare_pl],color=['r','g'],label=['Background','In plume'],normed=True,bins=30)
plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.title('Diurnally averaged Surface DARE for {}'.format(day))


# In[344]:


ssa_u[1235]


# In[311]:


np.where(flag&np.isfinite(ssa_u))[0][367]


# In[312]:


out['ext'][367,:,3]


# In[313]:


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

# In[6]:


day = '20180609'
#day = '20180618'
#day = '20180624'
#day = '20180625'


# In[7]:


out1 = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d='20180609',vv=vv))


# In[8]:


out2 = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d='20180618',vv=vv))


# In[9]:


out3 = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d='20180624',vv=vv))


# In[10]:


out4 = hs.loadmat(fp+'{name}_DARE_{d}_{vv}.mat'.format(name=name,d='20180625',vv=vv))


# ## Get the time tables 

# In[13]:


def time_utc(x):
    return np.array([y.hour+y.minute/60.0+y.second/3600.0 for y in x])


# ### 20180609

# In[14]:


flttable1 = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format('20180609'))
fromtime1 = flttable1['FromTime'][flttable1['FlightType']=='in plume']
totime1 = flttable1['ToTime'][flttable1['FlightType']=='in plume']
plumeid1 = flttable1['PlumeId'][(flttable1['PlumeId']=='A') | (flttable1['PlumeId']=='B')].to_numpy()


# In[15]:


from_utc1 = time_utc(fromtime1.to_numpy())
to_utc1 = time_utc(totime1.to_numpy())


# In[23]:


from_utc1


# In[17]:


out1 = out


# In[71]:


# Remove bad data from neph scattering before 16:25
ibad1 = (out1['utc']<16.4) | (out1['utc']>17.6)
out1['dare_avg'][ibad1,:] = np.nan
out1['aod'][ibad1,:] = np.nan
out1['ssa'][ibad1,:] = np.nan


# In[72]:


ipl1 = []
dare_pl1,dare_pl1a,dare_pl1b = [],[],[]
aod_pl1,aod_pl1a,aod_pl1b,aod_out1 = [],[],[],[]
ssa_pl1,ssa_pl1a,ssa_pl1b,ssa_out1 = [],[],[],[]
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
        aod_pl1 = np.append(aod_pl1,out1['aod'][pl1,3])
        ssa_pl1 = np.append(ssa_pl1,out1['ssa'][pl1,2,3])
        if plumeid1[ii]=='A':
            dare_pl1a = np.append(dare_pl1a,out1['dare_avg'][pl1,0])
            dare_pl1a_toa = np.append(dare_pl1a_toa,out1['dare_avg'][pl1,2])
            aod_pl1a = np.append(aod_pl1a,out1['aod'][pl1,3])
            ssa_pl1a = np.append(ssa_pl1a,out1['ssa'][pl1,2,3])
        else:
            dare_pl1b = np.append(dare_pl1b,out1['dare_avg'][pl1,0])
            dare_pl1b_toa = np.append(dare_pl1b_toa,out1['dare_avg'][pl1,2])
            aod_pl1b = np.append(aod_pl1b,out1['aod'][pl1,3])
            ssa_pl1b = np.append(ssa_pl1b,out1['ssa'][pl1,2,3])
        ipl1.append(pl1)
    dare_out1 = np.append(dare_out1,out1['dare_avg'][~pl1,0])
    dare_out1_toa = np.append(dare_out1_toa,out1['dare_avg'][~pl1,2])
    alt_out1 = np.append(alt_out1,out1['alt'][~pl1,0])
    aod_out1 = np.append(aod_out1,out1['aod'][~pl1,3])
    ssa_out1 = np.append(ssa_out1,out1['ssa'][~pl1,2,3])


# In[73]:


out1.keys()


# In[74]:


np.nanmean(aod_pl1a),np.nanmean(aod_pl1b),np.nanmean(aod_out1)


# In[75]:


np.nanmean(ssa_pl1a),np.nanmean(ssa_pl1b),np.nanmean(ssa_out1)


# In[69]:


ppl1 = np.array(ipl1).flatten()


# In[19]:


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


# In[76]:


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


# In[38]:


out1['ext'].shape


# In[77]:


avgs1 = {'bmea':np.nanmean(dare_out1),'bmed':np.nanmedian(dare_out1),'bstd':np.nanstd(dare_out1),
         'pamea':np.nanmean(dare_pl1a),'pamed':np.nanmedian(dare_pl1a),'pastd':np.nanstd(dare_pl1a),
         'pbmea':np.nanmean(dare_pl1b),'pbmed':np.nanmedian(dare_pl1b),'pbstd':np.nanstd(dare_pl1b)}


# In[78]:


np.nanmin(dare_pl1a),np.nanmax(dare_pl1a)


# In[79]:


np.nanmin(dare_pl1b),np.nanmax(dare_pl1b)


# In[80]:


np.nanmin(dare_out1),np.nanmax(dare_out1)


# In[81]:


plt.figure(figsize=(7,3))
plt.hist([dare_pl1a,dare_pl1b,dare_out1],color=['violet','silver','lightgreen'],
         label=['In plume A, {pamea:2.1f} $\\pm$ {pastd:2.1f} W/m$^2$'.format(**avgs1),
                'In plume B, {pbmea:2.1f} $\\pm$ {pbstd:2.1f} W/m$^2$'.format(**avgs1),
                'Out of plume, {bmea:2.1f} $\\pm$ {bstd:2.1f} W/m$^2$'.format(**avgs1)],
         normed=True,bins=30,range=[-90,0])

plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.ylabel('Normalized counts')
plt.title('Diurnally averaged Surface DARE for {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plots/DARE_hist_surface_background_inplume_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[82]:


avgs1_toa = {'bmea':np.nanmean(dare_out1_toa),'bmed':np.nanmedian(dare_out1_toa),'bstd':np.nanstd(dare_out1_toa),
         'pamea':np.nanmean(dare_pl1a_toa),'pamed':np.nanmedian(dare_pl1a_toa),'pastd':np.nanstd(dare_pl1a_toa),
         'pbmea':np.nanmean(dare_pl1b_toa),'pbmed':np.nanmedian(dare_pl1b_toa),'pbstd':np.nanstd(dare_pl1b_toa)}


# In[83]:


avgs1_toa_rg = {'bmin':np.nanmin(dare_out1_toa),'bmax':np.nanmax(dare_out1_toa),
         'pamin':np.nanmin(dare_pl1a_toa),'pamax':np.nanmax(dare_pl1a_toa),
         'pbmin':np.nanmin(dare_pl1b_toa),'pbmax':np.nanmax(dare_pl1b_toa)}
avgs1_toa_rg


# In[84]:


plt.figure(figsize=(7,3))
plt.hist([dare_pl1a_toa,dare_pl1b_toa,dare_out1_toa],color=['violet','silver','lightgreen'],
         label=['In plume A, {pamea:2.1f} $\\pm$ {pastd:2.1f} W/m$^2$'.format(**avgs1_toa),
                'In plume B, {pbmea:2.1f} $\\pm$ {pbstd:2.1f} W/m$^2$'.format(**avgs1_toa),
                'Out of plume, {bmea:2.1f} $\\pm$ {bstd:2.1f} W/m$^2$'.format(**avgs1_toa)],
         normed=True,bins=30,range=[-50,0])

plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.ylabel('Normalized counts')
plt.title('Diurnally averaged TOA DARE for {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plots/DARE_hist_TOA_background_inplume_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[24]:


np.nanmean(dare_pl1a_toa[dare_pl1a_toa>(-15)]), np.nanmean(dare_pl1a_toa[dare_pl1a_toa<(-15)])


# In[85]:


fig = plt.figure()
ax1 = fig.add_subplot(411)
ax1.plot(out1['utc'],out1['dare_avg'][:,2],'.k',label='all')
for ii,ipp in enumerate(ipl1):
    if plumeid1[ii]=='A':
        ax1.plot(out1['utc'][ipp],out1['dare_avg'][ipp,2],'o',color='magenta')
    else:
        ax1.plot(out1['utc'][ipp],out1['dare_avg'][ipp,2],'o',color='grey')
ax1.set_ylabel('DARE')
    
ax2 = fig.add_subplot(412,sharex = ax1)
ax2.plot(out1['utc'],out1['aod'][:,0],'.k',label='all')
for ii,ipp in enumerate(ipl1):
    if plumeid1[ii]=='A':
        ax2.plot(out1['utc'][ipp],out1['aod'][ipp,0],'o',color='magenta')
    else:
        ax2.plot(out1['utc'][ipp],out1['aod'][ipp,0],'o',color='grey')
ax2.set_ylabel('AOD')
    
ax3 = fig.add_subplot(413,sharex = ax1)
ax3.plot(out1['utc'],out1['ssa'][:,2,3],'.k',label='all')
for ii,ipp in enumerate(ipl1):
    if plumeid1[ii]=='A':
        ax3.plot(out1['utc'][ipp],out1['ssa'][ipp,2,3],'o',color='magenta')
    else:
        ax3.plot(out1['utc'][ipp],out1['ssa'][ipp,2,3],'o',color='grey')
ax3.set_ylabel('SSA')

ax4 = fig.add_subplot(414,sharex = ax1)
ax4.plot(out1['utc'],out1['asy'][:,2,3],'.k',label='all')
for ipp in ipl1:
    ax4.plot(out1['utc'][ipp],out1['asy'][ipp,2,3],'o')
ax4.set_ylabel('ASY')


# In[ ]:


fig = plt.figure()
ax1 = fig.add_subplot(411)
ax1.plot(out1['utc'],out1['dare_avg'][:,2],'.k',label='all')
for ipp in ipl1:
    ax1.plot(out1['utc'][ipp],out1['dare_avg'][ipp,2],'o')
ax1.set_ylabel('DARE')
    
ax2 = fig.add_subplot(412,sharex = ax1)
ax2.plot(out1['utc'],out1['aod'][:,0],'.k',label='all')
for ipp in ipl1:
    ax2.plot(out1['utc'][ipp],out1['aod'][ipp,0],'o')
ax2.set_ylabel('AOD')
    
ax3 = fig.add_subplot(413,sharex = ax1)
ax3.plot(out1['utc'],out1['ssa'][:,2,3],'.k',label='all')
for ipp in ipl1:
    ax3.plot(out1['utc'][ipp],out1['ssa'][ipp,2,3],'o')
ax3.set_ylabel('SSA')

ax4 = fig.add_subplot(414,sharex = ax1)
ax4.plot(out1['utc'],out1['asy'][:,2,3],'.k',label='all')
for ipp in ipl1:
    ax4.plot(out1['utc'][ipp],out1['asy'][ipp,2,3],'o')
ax4.set_ylabel('ASY')


# #### Paired T-test to see if differences in plumes vs background

# In[176]:


from scipy import stats


# In[86]:


#[dare_pl1a_toa,dare_pl1b_toa,dare_out1_toa]
tval,pval = stats.ttest_ind(dare_pl1b_toa,dare_out1_toa,nan_policy='omit')
tval,pval


# In[87]:


tval,pval = stats.ttest_ind(dare_pl1a_toa,dare_out1_toa,nan_policy='omit',equal_var=False)
tval,pval


# #### Estimate the uncertainty in DARE

# In[403]:


import plotting_utils as pu


# In[401]:


out1['ssa'].shape


# In[409]:


plt.figure()
plt.plot(out1['ssa'][:,2,3],out1['dare_avg'][:,0],'.',label='All data')
pu.plot_lin(out1['ssa'][:,2,3],out1['dare_avg'][:,0],labels=True)
plt.xlabel('SSA at 500 nm')
plt.ylabel('DARE surface')
plt.legend()


# In[422]:


np.nanmean(out1['ssa'][:,2,3])


# In[425]:


delta_ssa = 0.1*np.nanmean(out1['ssa'][:,2,3]) #assuming 10% uncertainty in SSA
delta_dare_ssa = 40.41*delta_ssa
delta_dare_ssa


# In[410]:


plt.figure()
plt.plot(out1['aod'][:,3],out1['dare_avg'][:,0],'.',label='All data')
pu.plot_lin(out1['aod'][:,3],out1['dare_avg'][:,0],labels=True)
plt.xlabel('AOD at 500 nm')
plt.ylabel('DARE surface')
plt.legend()


# In[442]:


plt.figure()
plt.plot(aod_pl1a,dare_pl1a,'.',color='violet',label='plume A')
pu.plot_lin(aod_pl1a,dare_pl1a,color='violet',labels=True)
plt.plot(aod_pl1b,dare_pl1b,'.',color='silver',label='plume B')
pu.plot_lin(aod_pl1b,dare_pl1b,color='silver',labels=True)
plt.plot(aod_out1,dare_out1,'.',color='lightgreen',label='out of plume')
pu.plot_lin(aod_out1,dare_out1,color='lightgreen',labels=True)
plt.legend()
#color=['violet','silver','lightgreen']


# In[426]:


delta_aod = 0.02 #assuming 0.02 uncertainty in AOD
delta_dare = -88.70*delta_aod
delta_dare


# In[433]:


np.sqrt(delta_dare**2.0+delta_dare_ssa**2.0)


# In[432]:


np.sqrt((delta_dare_ssa*delta_ssa)**2.0+(delta_dare*delta_aod)**2.0)


# In[418]:


ac = 16.0
sc = 500.0
da = 3.2
ds = 100.0
dn = (1.0/(ac+sc)**2)
up = np.sqrt((sc*da)**2+(ac*sc*ds)**2)


# In[419]:


up*dn


# In[454]:


out1['dn_clear'].shape


# In[458]:


out1['dn_clear'][0,:,0]*0.2412


# In[460]:


dare_milagro = np.nanmean(out1['dn_clear'][0,:,0] - out1['dn_clear'][0,:,0]*0.2412)
dare_milagro


# In[456]:


plt.figure()
plt.plot(out1['dn_clear'][0,:,2])


# ### 20180618

# In[135]:


flttable2 = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format('20180618'))
fromtime2 = flttable2['FromTime'][flttable2['FlightType']=='in plume']
totime2 = flttable2['ToTime'][flttable2['FlightType']=='in plume']


# In[136]:


plumeid2 = flttable2['PlumeId'][(flttable2['PlumeId']=='A') | (flttable2['PlumeId']=='B')].to_numpy()


# In[ ]:


flttable2


# In[139]:


from_utc2 = time_utc(fromtime2.to_numpy())
to_utc2 = time_utc(totime2.to_numpy())


# In[140]:


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


# In[144]:


avgs2 = {'bmea':np.nanmean(dare_out2),'bmed':np.nanmedian(dare_out2),'bstd':np.nanstd(dare_out2),
         'pamea':np.nanmean(dare_pl2a),'pamed':np.nanmedian(dare_pl2a),'pastd':np.nanstd(dare_pl2a),
         'pbmea':np.nanmean(dare_pl2b),'pbmed':np.nanmedian(dare_pl2b),'pbstd':np.nanstd(dare_pl2b)}


# In[146]:


plt.figure(figsize=(7,3))
plt.hist([dare_pl2a,dare_pl2b,dare_out2],color=['violet','silver','lightgreen'],
         label=['In plume A, {pamea:2.2f} $\\pm$ {pastd:2.2f} W/m$^2$'.format(**avgs2),
                'In plume B, {pbmea:2.2f} $\\pm$ {pbstd:2.2f} W/m$^2$'.format(**avgs2),
                'Background, {bmea:2.2f} $\\pm$ {bstd:2.2f} W/m$^2$'.format(**avgs2)],
         normed=True,bins=30)

plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.ylabel('Normalized counts')
plt.title('Diurnally averaged Surface DARE for {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plots/DARE_hist_surface_background_inplume_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[147]:


avgs2_toa = {'bmea':np.nanmean(dare_out2_toa),'bmed':np.nanmedian(dare_out2_toa),'bstd':np.nanstd(dare_out2_toa),
         'pamea':np.nanmean(dare_pl2a_toa),'pamed':np.nanmedian(dare_pl2a_toa),'pastd':np.nanstd(dare_pl2a_toa),
         'pbmea':np.nanmean(dare_pl2b_toa),'pbmed':np.nanmedian(dare_pl2b_toa),'pbstd':np.nanstd(dare_pl2b_toa)}


# In[148]:


plt.figure(figsize=(7,3))
plt.hist([dare_pl2a_toa,dare_pl2b_toa,dare_out2_toa],color=['violet','silver','lightgreen'],
         label=['In plume A, {pamea:2.2f} $\\pm$ {pastd:2.2f} W/m$^2$'.format(**avgs2_toa),
                'In plume B, {pbmea:2.2f} $\\pm$ {pbstd:2.2f} W/m$^2$'.format(**avgs2_toa),
                'Background, {bmea:2.2f} $\\pm$ {bstd:2.2f} W/m$^2$'.format(**avgs2_toa)],
         normed=True,bins=30)

plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.ylabel('Normalized counts')
plt.title('Diurnally averaged TOA DARE for {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plots/DARE_hist_TOA_background_inplume_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# ### 20180624

# In[271]:


flttable3 = pd.read_excel(fp+'flt_table/fltable_{}.xlsx'.format('20180624'))
fromtime3 = flttable3['FromTime'][flttable3['FlightType']=='in plume']
totime3 = flttable3['ToTime'][flttable3['FlightType']=='in plume']


# In[272]:


from_utc3 = time_utc(fromtime3.to_numpy())
to_utc3 = time_utc(totime3.to_numpy())


# In[273]:


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


# In[274]:


fig = plt.figure()
ax1 = fig.add_subplot(311)
ax1.plot(out3['utc'],out3['dare_avg'][:,0],'.k',label='all')
for ipp in ipl3:
    ax1.plot(out3['utc'][ipp],out3['dare_avg'][ipp,0],'o')
ax1.set_ylabel('DARE')
    
ax2 = fig.add_subplot(312,sharex = ax1)
ax2.plot(out3['utc'],out3['aod'][:,0],'.k',label='all')
for ipp in ipl3:
    ax2.plot(out3['utc'][ipp],out3['aod'][ipp,0],'o')
ax2.set_ylabel('AOD')
    
ax3 = fig.add_subplot(313,sharex = ax1)
ax3.plot(out3['utc'],out3['ssa'][:,2,3],'.k',label='all')
for ipp in ipl3:
    ax3.plot(out3['utc'][ipp],out3['ssa'][ipp,2,3],'o')
ax3.set_ylabel('SSA')


# In[278]:


avgs3 = {'bmea':np.nanmean(dare_out3),'bmed':np.nanmedian(dare_out3),'bstd':np.nanstd(dare_out3),
         'pamea':np.nanmean(dare_pl3),'pamed':np.nanmedian(dare_pl3),'pastd':np.nanstd(dare_pl3)}


# In[282]:


plt.figure(figsize=(7,3))
plt.hist([dare_pl2,dare_out2],color=['violet','lightgreen'],
         label=['In plume, {pamea:2.2f} $\\pm$ {pastd:2.2f} W/m$^2$'.format(**avgs2),
                'Background, {bmea:2.2f} $\\pm$ {bstd:2.2f} W/m$^2$'.format(**avgs2)],
         normed=True,bins=30)

plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.ylabel('Normalized counts')
plt.title('Diurnally averaged Surface DARE for {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plots/DARE_hist_surface_background_inplume_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


# In[281]:


avgs2_toa = {'bmea':np.nanmean(dare_out2_toa),'bmed':np.nanmedian(dare_out2_toa),'bstd':np.nanstd(dare_out2_toa),
         'pamea':np.nanmean(dare_pl2_toa),'pamed':np.nanmedian(dare_pl2_toa),'pastd':np.nanstd(dare_pl2_toa)}


# In[284]:


plt.figure(figsize=(7,3))
plt.hist([dare_pl2_toa,dare_out2_toa],color=['violet','lightgreen'],
         label=['In plume, {pamea:2.2f} $\\pm$ {pastd:2.2f} W/m$^2$'.format(**avgs2_toa),
                'Background, {bmea:2.2f} $\\pm$ {bstd:2.2f} W/m$^2$'.format(**avgs2_toa)],
         normed=True,bins=30)

plt.legend(frameon=False)
plt.xlabel('DARE [W/m$^2$]')
plt.ylabel('Normalized counts')
plt.title('Diurnally averaged TOA DARE for {}'.format(day))
plt.tight_layout()
plt.savefig(fp+'plots/DARE_hist_TOA_background_inplume_{}_{}.png'.format(day,vv),dpi=600,transparent=True)


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
# Old Versions (Prior to Jan 2020) - Should not use  

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





# # Get the p-value for the delta tau vs. distance of AERONET vs 4STAR

# In[177]:


import plotting_utils as pu
from scipy import stats


# In[9]:


s2aero = pd.read_excel(fp+'starAero.xlsx'.format('20180609'))


# In[10]:


s2aero.keys()


# In[14]:


plt.figure()
plt.plot(s2aero['distKM'],s2aero['deltaaodfine'],'.',label='')
p,perr = pu.plot_lin(s2aero['distKM'],s2aero['deltaaodfine'],y_err=s2aero['deltaaodfine']*0.0+0.03,labels=True,lblfmt='4.5f')
plt.axhline(0,ls='--')
plt.legend()
p,perr
plt.xlabel('Distance from AERONET [km]')
plt.ylabel('$\\tau_{diff}$')
slope, intercept, r_value, p_value, std_err = stats.linregress(s2aero['distKM'],s2aero['deltaaodfine'])
p_value


# ## Subset for different days

# In[168]:


s2aero['day'] = s2aero['Time'].dt.date


# In[169]:


ds = s2aero['day'].unique()


# In[170]:


d = ds[0]


# In[173]:


print '{}'.format(d)


# In[209]:


plt.figure(figsize=(8,4))
for d in s2aero['day'].unique():
    iid = (s2aero['day']==d)
    slope, intercept, r_value, p_value, std_err = stats.linregress(s2aero['distKM'][iid],s2aero['deltaaodfine'][iid])
    pl = plt.plot(s2aero['distKM'][iid],s2aero['deltaaodfine'][iid],'.',label='{},   p={:1.4f}'.format(d,p_value))
    p,perr = pu.plot_lin(s2aero['distKM'][iid],s2aero['deltaaodfine'][iid],y_err=s2aero['deltaaodfine']*0.0+0.03,
                         labels=True,lblfmt='4.5f',color=pl[0].get_c())
    
    
    p,perr
    plt.xlabel('Distance from AERONET [km]')
    plt.ylabel('$\\tau_{diff}$')
    
    print d, p_value
plt.axhline(0,ls='--',color='k')
slope, intercept, r_value, p_value, std_err = stats.linregress(s2aero['distKM'],s2aero['deltaaodfine'])
plt.plot(s2aero['distKM'],s2aero['deltaaodfine'],'.',label='All data,  p={:1.4f}'.format(p_value),zorder=-10,color='grey')
p,perr = pu.plot_lin(s2aero['distKM'],s2aero['deltaaodfine'],y_err=s2aero['deltaaodfine']*0.0+0.03,
                     labels=True,lblfmt='4.3f',color='grey',shaded_ci=False)

plt.legend(bbox_to_anchor=[1.05,0.99])
plt.tight_layout()
plt.savefig(fp+'COSR_distance_to_AERONET_taudiff.png',dpi=600,transparent=True)


# # Add MODIS AOD plots

# In[65]:


from mpl_toolkits.basemap import Basemap
import georaster


# In[17]:


mod,modh = lu.load_hdf(fp+'data_other/MOD04_3K.A2018160.1855.061.2018161080147.hdf',values=(('lat',52),('lon',51),('aod',10)))


# In[112]:


modh['aod']


# In[62]:


def make_map1(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=-111.0,lat_0=57.6,
            llcrnrlon=-113.1, llcrnrlat=56.6,
            urcrnrlon=-109.1, urcrnrlat=58.6,resolution='h',ax=ax)
    #m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(-109,-113.2,7),labels=[0,0,0,1])
    m.drawparallels(np.linspace(56.6,58.6,11),labels=[1,0,0,0])
    return m


# In[54]:


plt.figure()
plt.scatter(mod['lon'],mod['lat'],50,mod['aod'],marker='s')
plt.colorbar()


# In[67]:


import rasterio
from rasterio.plot import show


# In[68]:


src = rasterio.open(fp+'data_other/snapshot-2018-06-09T00_00_00Z.tiff')


# In[136]:


fla = np.where(flag & (s['Alt'][:,0]<1500.0))


# In[143]:


plt.figure()
show(src.read(),transform=src.transform)
plt.xlim(-113.1,-109.8)
plt.ylim(56.3,58.3)
plt.pcolor(mod['lon'],mod['lat'],mod['aod'],vmin=0,vmax=0.5,cmap='plasma',label='MODIS DT TERRA')
plt.colorbar(label='AOD',extend='max',shrink=0.68,pad=0.03)
plt.plot(s['Lon'],s['Lat'],'-',color='k',markersize=0.2,label='flight path')
plt.scatter(s['Lon'][fla,0],s['Lat'][fla,0],50,s['tau_aero'][fla,406],marker='o',
            cmap='plasma',vmin=0,vmax=0.5,zorder=10,label='4STAR',lw=0.1,edgecolor='k',alpha=0.2)
plt.legend()
plt.xlabel('Longitude [$^{{\circ}}$]')
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.savefig(fp+'COSR_20180609_AOD_MODIS_4STAR_Truecolor.png',dpi=600,transparent=True)

