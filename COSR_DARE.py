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

# In[119]:


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


# In[15]:


import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')


# In[3]:


name = 'COSR'
vv = 'v1'


# In[4]:


fp =getpath('COSR')
fp_rtm = '/nobackup/sleblan2/rtm/'
fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
fp_rtmdat = '/nobackup/sleblan2/COSR/rtm/' #'/u/sleblan2/4STAR/rtm_dat/'


# In[5]:


day = '20180624'


# ## Setup command line interface

# In[6]:


import argparse


# In[7]:


long_description = """    Prepare or save the direct Aerosol radiative effect files for calculations. """


# In[32]:


parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-doread','--doread',help='if set, will only read the output, not produce them',
                    action='store_true')
parser.add_argument('-d','--daystr',nargs='?',help='The day string (yyyymmdd) for the desired flight data. Defaults to 20180624')


# In[9]:


in_ = vars(parser.parse_args())
do_read = in_.get('doread',False)
day = in_.get('daystr','20180624')


# # Load files

# ## Load the in situ files

# In[11]:


## neph = sio.loadmat(fp+'20180624_nephclap.csv_20191004_152550.mat.mat')


# In[15]:


## neph.keys()


# In[28]:


## neph['None'][0][4]


# In[33]:


situ = pd.read_csv(fp+'data_other/{}_nephclap.csv'.format(day))


# In[34]:


situ


# In[35]:


insitu = situ.to_dict('list')


# In[36]:


insitu.keys()


# In[37]:


insitu['ssa_500nm'] = np.array(insitu['totScatCalc_500nm'])/np.array(insitu['extCalc500nm'])


# In[38]:


plt.figure()
plt.plot(insitu['ssa_500nm'],'.')
plt.plot(Sp.smooth(insitu['ssa_500nm'],60,old=True),'-r')


# In[314]:


ssa_insitu = Sp.smooth(insitu['ssa_500nm'],60,old=True)


# In[315]:


len(ssa_insitu)


# In[316]:


len(s['utc'])


# ## Load the 4STAR AOD

# In[331]:


import os


# In[332]:


os.path.exists(fp+'os_data/4STAR_{}starsun.mat'.format(day))


# In[334]:


s = sio.loadmat(fp+'os_data/4STAR_{}starsun.mat'.format(day))


# In[335]:


s.keys()


# In[336]:


s['utc'] = lu.toutc(lu.mat2py_time(s['t']))


# ### use the polyfit aod on the wavelength array

# In[337]:


s['tau_aero_polynomial'].shape


# In[338]:


wvl = np.array([0.35,0.4,0.5,0.675,0.87,0.995,1.2,1.4,1.6,2.1,4.0])


# In[356]:


s['aod'] = np.zeros((len(s['utc']),len(wvl)))
for i in xrange(len(s['utc'])):
    s['aod'][i,:] = np.exp(np.polyval([s['tau_aero_polynomial'][i,0],s['tau_aero_polynomial'][i,1],s['tau_aero_polynomial'][i,2]],
                                         np.log(wvl)))


# In[357]:


s['aod'].shape


# In[358]:


s['aod'][10:-1:200,:].shape


# In[360]:


plt.figure()
plt.plot(wvl,s['aod'][10:-1:200,:].T)
plt.ylim(0,1)


# ### Load the flag files

# In[362]:


fmat = getpath('4STAR_data')


# In[365]:


with open (fmat+'starinfo_{}.m'.format(day), 'rt') as in_file:
    for line in in_file:
        if 'flagfilename ' in line:
            ff = line.split("'")[1]
sf = hs.loadmat(fmat+ff)


# In[366]:


sf.keys()


# In[377]:


flag = sf['manual_flags']['good'][0,:,0]


# In[378]:


flag.shape


# In[379]:


sum(flag)


# ## Load the skyscan results

# In[77]:


fp_name = '4STAR_20180624_135_SKYP.created_20190329_003621.ppl_lv15.mat'
sky = sio.loadmat(fp+fp_name)


# In[32]:


sky.keys()


# # Plot out some data

# ## Plot out the retrieved skyscans

# In[204]:


plt.figure()
plt.plot(sky['Wavelength'],sky['ssa_total'][0,:],'x-',label='SSA')
plt.plot(sky['Wavelength'],sky['g_tot'][0,:],'x-',label='ASY')
plt.plot(sky['Wavelength'],sky['sfc_alb'][0,:],'x-',label='Albedo')


plt.legend(frameon=False)
plt.xlabel('Wavelength [micron]')
plt.title('4STAR skyscan results from: \n' + fp_name)
plt.savefig(fp+'plots/4STAR_skyscan_result_20180624_135_SKYP.png',dpi=600,transparent=True)


# In[232]:


sky['g_tot'][-1]


# In[295]:


wvl = np.array([0.35,0.4,0.5,0.675,0.87,0.995,1.2,1.4,1.6,2.1,4.0])
f_asy = interp1d(np.append(sky['Wavelength'][:,0],[1.1,1.2]),np.append(sky['g_tot'][0,:],[sky['g_tot'][0,-1]-0.008,sky['g_tot'][0,-1]-0.011]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
asy = f_asy(wvl)
f_ssa = interp1d(np.append(sky['Wavelength'][:,0],[1.1,1.2,1.3]),np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.008,sky['ssa_total'][0,-1]-0.01,sky['ssa_total'][0,-1]-0.014]),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
ssa = f_ssa(wvl)
f_alb = interp1d(np.append(sky['Wavelength'][:,0],1.1),np.append(sky['sfc_alb'][0,:],sky['sfc_alb'][0,-1]+0.003),
                 bounds_error=False,fill_value='extrapolate',kind='slinear')
alb = f_alb(wvl)


# In[283]:


np.append(sky['ssa_total'][0,:],[sky['ssa_total'][0,-1]-0.008,sky['ssa_total'][0,-1]-0.002])


# In[298]:


plt.figure()
plt.plot(wvl,ssa,'-x',label='SSA')
plt.plot(wvl,asy,'-*',label='ASY')
plt.plot(wvl,alb,'-s',label='Albedo')
plt.legend()
plt.xlabel('Wavelength [micron]')


# ## Get the vertical dependence of the extinction

# In[40]:


gu = pd.to_datetime(situ['DateTimeUTC']).to_list()
insitu['utc'] = np.array([g.hour+g.minute/60.0+g.second/3600.0 for g in gu])


# In[41]:


from scipy.interpolate import interp1d


# In[42]:


f_alt = interp1d(x=s['utc'],y=s['Alt'][:,0])
insitu['alt'] = f_alt(insitu['utc'])


# In[43]:


plt.figure()
plt.plot(insitu['extCalc500nm'],insitu['alt'],'.')


# In[26]:


plt.figure()
plt.plot(insitu['utc'],insitu['alt'])


# In[197]:


plt.figure()
plt.scatter(insitu['utc'],insitu['alt'],c=insitu['extCalc500nm'],cmap=plt.cm.magma,s=(insitu['extCalc500nm'])**1.5-20.0)
plt.colorbar()


# In[67]:


insitu['extCalc500nm'] = np.array(insitu['extCalc500nm'])


# In[71]:


np.isfinite(insitu['extCalc500nm'])


# In[160]:


binned_ext,binned_alt,binned_num = [],[],[]
for i in xrange(14):
    flaa = (insitu['alt']>=i*100.0) & (insitu['alt']<(i+1.0)*100.0) & (np.isfinite(insitu['extCalc500nm']))
    if flaa.any():
        binned_ext.append(insitu['extCalc500nm'][flaa])
        binned_alt.append(np.mean([i*100.0,(i+1.0)*100.0]))
        binned_num.append(len(insitu['extCalc500nm'][flaa]))


# In[203]:


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
plt.title('In situ calculated extinction CLAP+neph: 20180624')
plt.savefig(fp+'plots/extinction_vertical_bins_clap_neph_20180624.png',dpi=600,transparent=True)


# In[166]:


ext_z = ext_means[:,1]/1000.0
ext_ = ext_means[:,0]/1000.0
aod_ = ext_.sum()/10.0


# In[312]:


nz = len(ext_z)


# In[323]:


ext_,aod_


# # Now build the input files for DARE calculations
# 

# In[382]:


import Run_libradtran as RL


# In[467]:


#get the paths
f = open(fp+'rtm/{}_CRE_{}.sh'.format(name,vv),'w')
fpp_in = '/nobackup/sleblan2/rtm/input/{}_CRE_{}/'.format(name,vv)
fpp_out = '/nobackup/sleblan2/rtm/output/{}_CRE_{}/'.format(name,vv)
fp_uv = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
fp_in = fp+'rtm/input/CRE/'


# In[456]:


geo = {'zout':[0,3,100],'year':2018,'month':6,'day':24,'hour':12,'minute':0,'second':0}
aero_base = {'z_arr':(ext_z+50.0)/1000.0,'wvl_arr':wvl,'ssa':np.array([ssa,]*nz),'asy':np.array([asy,]*nz)}
source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/','run_fuliou':True}
albedo = {'create_albedo_file':True,'alb':alb,'alb_wvl':wvl}


# In[457]:


fx_ssa_insitu = interp1d(insitu['utc'],ssa_insitu,bounds_error=False)
ssa_u = fx_ssa_insitu(s['utc'])
ssa_u[ssa_u<0.8] = np.nan


# In[458]:


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


# In[472]:


'{i:06d}'.format(i=5300)


# In[ ]:


try:
    file_list = file(fp_rtm+'COSR_DARE_list_file_{d}_{v}.sh'.format(d=day,v=vv),'w')
except Exception,e:
    print 'Problem with accessing file, return Exception: ',e
    return
print 'Starting list file'
fp_out = fp_rtm+'output/COSR_{d}_{v}/'.format(d=day,v=vv)
fp_in = fp_rtm+'input/COSR_{d}_{v}/'.format(d=day,v=vv)
if not os.path.exists(fp_out):
    os.mkdir(fp_out)
if not os.path.exists(fp_in):
    os.mkdir(fp_in)

fp_base_file = fp_out+'base.inp'
make_base = True
for i,u in enumerate(s['utc']):
    
    if flag[i] & np.isfinite(ssa_u[i]):
        aod = s['aod'][i,:]
        aero['ext'] = expand_ext_vert_and_spect(ext_,ext_z,s['aod'][i,:],s['Alt'][i]/1000.0,wvl)
        iw = np.argmin(abs(aero_base['wvl_arr']-0.5))
        aero['ssa'] = aero_base['ssa']*ssa_u[i]/aero_base['ssa'][0,iw]
        aero['asy'] = aero_base['asy']
        
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
        
        geo['utc'] = u
        geo['sza'] = s['sza'][i]
        geo['lat'] = s['Lat'][i]
        geo['lon'] = s['Lon'][i]
        
        fname = 'COSR_DARE_{d}_{v}_{i:06d}.dat'.format(d=day,v=vv,i=i)
        
        
        if not list_only:
            RL.write_fuliou_input(file_in,geo=geo,aero=aero,albedo=albedo,verbose=False)

        file_list.write(fp_fuliou+' '+file_in+' '+file_out+'\n')
        print mmm,ilat,ilon


# In[ ]:



stdfac_dict = RL.merge_dicts({'ext':0.0,'ssa':0.0,'asym':0.0,
                              'COD':0.0,'ref':0.0},stdfac_dict)
make_base = True

                #sanitize inputs after adding subtracting standard deviations
            try: aero['ssa'][aero['ssa']<0.0] = 0.0
            except: pass
            try: aero['ssa'][aero['ssa']>1.0] = 1.0
            except: pass
            try: aero['asy'][aero['asy']<0.0] = 0.0
            except: pass
            try: aero['asy'][aero['asy']>1.0] = 1.0
            except: pass
            
            aero['ext'] = np.array([aero['ext'],aero['ext']])
            aero['ssa'] = np.array([aero['ssa'],aero['ssa']])
            aero['asy'] = np.array([aero['asy'],aero['asy']])
            
            
            #for HH in xrange(24):
            geo['hour'] = lon/15.0*-1.0+12.0
            form = {'ilat':ilat,'ilon':ilon,'mmm':mmm}
            file_in = fp_out2+'AACFLin_lat{ilat:02.0f}_lon{ilon:02.0f}_{mmm}.datin'.format(**form)
            file_out = fp_output+'AACFLout_lat{ilat:02.0f}_lon{ilon:02.0f}_{mmm}.wrt'.format(**form)
                
            if not list_only:
                RF.write_fuliou_input(file_in,geo=geo,aero=aero,albedo=albedo,verbose=False)
          
            file_list.write(fp_fuliou+' '+file_in+' '+file_out+'\n')
            print mmm,ilat,ilon
    del alb_geo
    del input_mmm
    file_list.close()


# # Plotting
# Present some fo the early plots here

# In[ ]:




