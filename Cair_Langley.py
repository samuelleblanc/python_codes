
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
#     Written: Samuel LeBlanc, On fligth from BWI to SJC, 2016-08-06
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


# In[3]:

fp = 'C:\\Users\\sleblan2\\Research\\4STAR\\MLO_2016\\'


# # Read some files

# In[4]:

f = fp+'20160702_MLO5\\CAST_001_160702_090020_URC.csv'


# In[6]:

reload(cu)


# In[7]:

c = cu.read_Cair(fp+'20160702_MLO5\\CAST_001_160702_090020_URU.csv')


# In[8]:

c.keys()


# # Run analysis and calculate the airmass and rayleigh

# In[10]:

lat, lon, alt = 19.5365,-155.57615,3428.0


# In[26]:

reload(su)


# ## Calculate the airmass and sza

# In[12]:

c = su.calc_sza_airmass(c['DateTimeUTC'],lat,lon,alt,c=c)


# In[13]:

c.keys()


# In[33]:

c['Lt'].shape


# ## get the rayleigh tau

# In[14]:

c['wvl']


# In[34]:

c['tau_rayleigh'],c['tau_rayleigh_err'] = np.zeros_like(c['Lt']),np.zeros_like(c['Lt'])


# In[35]:

for i,d in enumerate(c['DateTimeUTC']):
    c['tau_rayleigh'][:,i],c['tau_rayleigh_err'][:,i] = su.tau_rayleigh(np.array(c['wvl'])/1000.0,680.0,latitude=c['lat'],date=d)


# ## Calculate the 'rateaero' or Lt_aero
# Which is the Lt values divided by the impact of rayleigh, and trace gases

# In[45]:

c['Lt_aero'] = c['Lt']/c['sunearthf']/np.exp(-c['m_aero']*c['tau_rayleigh'])


# In[ ]:



