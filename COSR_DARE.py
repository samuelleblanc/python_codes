
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

# In[13]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
from datetime import datetime
from scipy.interpolate import UnivariateSpline


# In[29]:


import pandas as pd


# In[2]:


import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')


# In[31]:


name = 'COSR'
vv = 'v1'


# In[5]:


fp =getpath('COSR')
fp_rtm = '/nobackup/sleblan2/rtm/'
fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
fp_rtmdat = '/nobackup/sleblan2/COSR/rtm/' #'/u/sleblan2/4STAR/rtm_dat/'


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

# In[14]:


neph = sio.loadmat(fp+'20180624_nephclap.csv_20191004_152550.mat.mat')


# In[15]:


neph.keys()


# In[28]:


neph['None'][0][4]


# In[34]:


situ = pd.read_csv(fp+'{}_neph_clap.csv'.format(day))


# In[36]:


situ


# In[44]:


insitu = situ.to_dict('list')


# In[46]:


insitu.keys()


# In[56]:


[insitu[k]=np.array(insitu[k]) for k in insitu.keys()]


# In[47]:


insitu['ssa_500nm'] = np.array(insitu['totScatCalc_500nm'])/np.array(insitu['extCalc500nm'])


# In[50]:


import Sp_parameters as Sp


# In[55]:


plt.figure()
plt.plot(insitu['ssa_500nm'],'.')
plt.plot(Sp.smooth(insitu['ssa_500nm'],60,old=True),'-r')


# ## Load the skyscan results

# In[17]:


sky = sio.loadmat(fp+'4STAR_20180624_135_SKYP.created_20190329_003621.ppl_lv15.mat')


# In[18]:


sky.keys()


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# # Plotting
# Present some fo the early plots here
