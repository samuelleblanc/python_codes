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


# In[14]:


plt.figure()
plt.plot(insitu['ssa_500nm'],'.')
plt.plot(Sp.smooth(insitu['ssa_500nm'],60,old=True),'-r')


# ## Load the 4STAR AOD

# In[16]:


s = hs.loadmat(fp+'data/4STAR_{}starsun_for_starflag.mat'.format(day))


# In[17]:


s.keys()


# In[18]:


s['utc'] = lu.toutc(lu.mat2py_time(s['t']))


# ## Load the skyscan results

# In[19]:


sky = sio.loadmat(fp+'4STAR_20180624_135_SKYP.created_20190329_003621.ppl_lv15.mat')


# In[20]:


sky.keys()


# # Plot out some data

# ## Plot out the retrieved skyscans

# In[21]:


plt.figure()
plt.plot(sky['Wavelength'],sky['ssa_total'][0,:],'x-',label='SSA')
plt.plot(sky['Wavelength'],sky['g_tot'][0,:],'x-',label='ASY')

plt.legend(frameon=False)
plt.xlabel('Wavelength [micron]')


# ## Get the vertical dependence of the extinction

# In[22]:


gu = pd.to_datetime(situ['DateTimeUTC']).to_list()
insitu['utc'] = np.array([g.hour+g.minute/60.0+g.second/3600.0 for g in gu])


# In[23]:


from scipy.interpolate import interp1d


# In[24]:


f_alt = interp1d(x=s['utc'],y=s['Alt'][:,0])
insitu['alt'] = f_alt(insitu['utc'])


# In[25]:


plt.figure()
plt.plot(insitu['extCalc500nm'],insitu['alt'],'.')


# In[26]:


plt.figure()
plt.plot(insitu['utc'],insitu['alt'])


# In[ ]:


binned_ang,binned_alt,binned_num,binned_ndays = [],[],[],[]
for i in xrange(70):
    flaa = (s['GPS_Alt'][s['fl_QA_angs']]>=i*100.0) & (s['GPS_Alt'][s['fl_QA_angs']]<(i+1.0)*100.0)
    binned_ang.append(s['angs_470_865'][s['fl_QA_angs']][flaa])
    binned_alt.append(np.mean([i*100.0,(i+1.0)*100.0]))
    binned_num.append(len(s['angs_470_865'][s['fl_QA_angs']][flaa]))
    binned_ndays.append(len(np.unique(s['days'][s['fl_QA_angs']][flaa])))


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# # Plotting
# Present some fo the early plots here

# In[ ]:




