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


# ## Load the 4STAR AOD

# In[28]:


s = hs.loadmat(fp+'data/4STAR_{}starsun_for_starflag.mat'.format(day))


# In[29]:


s.keys()


# In[30]:


s['utc'] = lu.toutc(lu.mat2py_time(s['t']))


# ## Load the skyscan results

# In[77]:


fp_name = '4STAR_20180624_135_SKYP.created_20190329_003621.ppl_lv15.mat'
sky = sio.loadmat(fp+fp_name)


# In[32]:


sky.keys()


# # Plot out some data

# ## Plot out the retrieved skyscans

# In[81]:


plt.figure()
plt.plot(sky['Wavelength'],sky['ssa_total'][0,:],'x-',label='SSA')
plt.plot(sky['Wavelength'],sky['g_tot'][0,:],'x-',label='ASY')

plt.legend(frameon=False)
plt.xlabel('Wavelength [micron]')
plt.title('4STAR skyscan results from: \n' + fp_name)
plt.savefig(fp+'plots/4STAR_skyscan_result_20180624_135_SKYP.png',dpi=600,transparent=True)


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


# In[162]:


plt.figure()
bp =plt.boxplot(binned_ext,positions=binned_alt,vert=False,showfliers=False,widths=100,showmeans=True,patch_artist=True)
plt.xlabel('Extinction Calculated 500 nm')
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
ext_means = np.array([[b.get_data()[0][0],b.get_data()[1][0]] for b in bp['means']])
plt.plot(ext_means[:,0],ext_means[:,1],'-k')
plt.title('In situ calculated extinction CLAP+neph: 20180624')
plt.savefig(fp+'plots/extinction_vertical_bins_clap_neph_20180624.png',dpi=600,transparent=True)


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# # Plotting
# Present some fo the early plots here

# In[ ]:




