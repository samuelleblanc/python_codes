#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Redo of muSSTAR radiometer tests, now with python 3.
#     This is the repository of notes and test results
# 
# Input:
# 
#     None
# 
# Output:
# 
#     Figures and analysis
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
#     - write_utils
#     - path_utils
#     - scipy
# 
# Needed Files:
#   - 5STARG_* RADIOMETERS.dat files
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-05-11
#     Modified:
# 

# # Notes

# ## Test from 2021-05-05 with labjack
# Included is the data set from last Wednesday we were discussing. For this test the cold block was set to 20C and the hot block was set to 30C. The TEC NIR drivers were turned off to avoid a thermal overrun. This was with an air cooled system.
# 
# This test features a slightly new data format because we are now using the labjack. All of the normal channel data should be the same but there are now 3 new connectors to plug into. For this test I plugged a temp sensor and other things into J7 which was recorded. J7 has 8 analog inputs. They were used as follows:
# 
#  >J7_1 - 5v Power supply monitoring  
#  J7_2 - Cold block controller temp monitor  
#  J7_3 - Hot block controller temp monitor  
#  J7_4 - Cold block temp sensor  
#  J7_5 - NIR Board temp sensor  
#  J7_6 - Hot block temp sensor  
#  J7_7 - Lid temperature  
#  J7_8 - Ambient Temperature  
# 
# The cold and hot block controller temp monitors (J7_2 and J7_3) are different from our normal temperature sensors. This is the output from the thermistor that the controllers use for feedback control. It runs a constant current through the thermistors making a voltage output. Our current is set at 100 uA. I've included the thermistor data sheet. Unfortunately thermistor response is non linear at the macro level so unless you want to do a lot of ln() calculations there isn't a simple equation to convert from resistance to temperature. I'll list off the relative temp and equivalent voltages.
# 
# > 19C - 1.307V  
#  20C - 1.249V  
#  21C - 1.194V  
#   
# > 29C - 0.840V  
# 30C - 0.805V  
# 31C - 0.772V  
#   
# Probably not a lot of useful data in this test but it should be a good primer for the next data set that I'm working on right now.
# 
# Let me know if you have any questions.
# 
# Thanks,
# 
# Conrad Esch  
#  5STARG_20210505_152031_RADIOMETERS.dat

# ![image.png](attachment:image.png)

# ![image.png](attachment:image.png)

# ![image.png](attachment:image.png)

# ## Test from 2021-05-11 - Long test
# Included is the data set I got yesterday.
# 
# This should be almost the same as the last data set with one minor difference. I accidentally flipped the Lid and Ambient air temperature sensors so now J7_7 is in free air and J7_8 is attached to the lid. Apologies for the inconvenience. The test ran for 7 hours. Below are the labels for J7.
# 
# > J7_1 - 5v Power supply monitoring  
# J7_2 - Cold block controller temp monitor  
# J7_3 - Hot block controller temp monitor  
# J7_4 - Cold block temp sensor  
# J7_5 - NIR Board temp sensor  
# J7_6 - Hot block temp sensor  
# J7_7 - Ambient Temperature  
# J7_8 - Lid Temperature  
# 
#  5STARG_20210511_131357_RADIOMETERS.dat  
# 
# I'm running 2 more tests today. One is a cold start test to see how long we should wait for the system to stabilize (I suspect the temperatures are a bit oscillatory in the beginning.) The other is a second long test to verify the data in this email. Please let me know if you have any questions.
# 
# Thanks,
# 
# Conrad Esch

# # Prepare python environment

# In[13]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import os


# In[17]:


plt.rcParams['figure.dpi'] = 100 # 200 e.g. is really fine, but slower 


# In[16]:


name = 'muSSTAR'
vv = 'v1'
fp = getpath(name)


# In[101]:


import muSSTAR_utils as u
import importlib
importlib.reload(u)


# # Load and plot tests

# ## Test 2021-05-05

# In[58]:


s = u.reader_RADIOMETER(fp+'data/5STARG_20210505_152031_RADIOMETERS.dat',islabjack=True)


# In[59]:


s.keys()


# In[61]:


u.define_long_names(s)


# In[63]:


s.label = 'labjack_first_test'


# ### Plots

# In[71]:


u.clean_up_and_prep_for_plots(s,fp+'plots/')


# In[76]:


u.plt_hskp(s,fp+'plots/')


# ## Test 2021-05-12

# In[102]:


s = u.reader_RADIOMETER(fp+'data/5STARG_20210511_131357_RADIOMETERS.dat',islabjack=True)


# In[103]:


s.keys()


# In[79]:


u.define_long_names(s)


# In[80]:


s.label = 'labjack_long_test'


# ### Plots

# In[81]:


u.clean_up_and_prep_for_plots(s,fp+'plots/')


# In[82]:


u.plt_hskp(s,fp+'plots/')


# In[105]:


plt.figure()
plt.plot(s['UTC'],s['j7_1'])
plt.gca().xaxis.set_major_formatter(myFmt)


# In[100]:


plt.figure()
# 'J7_1' : '5v Power supply monitoring',
#        'J7_2' : 'Cold block controller temp monitor',
#        'J7_3' : 'Hot block controller temp monitor',
#        'J7_4' : 'Cold block temp sensor',
#        'J7_5' : 'NIR Board temp sensor',
#        'J7_6' : 'Hot block temp sensor',
#        'J7_7':'Ambient Temperature',
#        'J7_8':'Lid Temperature'
import matplotlib.dates as mdates
plt.plot(s['UTC'],s['j7_2'],label='cold block ctrl')
plt.plot(s['UTC'],s['j7_3'],label='hot block ctrl')
plt.plot(s['UTC'],s['j7_4'],label='cold block temp')
plt.plot(s['UTC'],s['j7_5'],label='NIR board temp')
plt.plot(s['UTC'],s['j7_6'],label='hot block temp')
plt.plot(s['UTC'],s['j7_7'],label='Ambient temp')
plt.plot(s['UTC'],s['j7_8'],label='Lid temp')
myFmt = mdates.DateFormatter('%H:%M:%S')
plt.gca().xaxis.set_major_formatter(myFmt)
plt.legend()


# In[88]:


u.plt_gains(s,fp+'plots/')


# In[89]:


u.plot_channels(s,fp+'plots/{s.daystr}_{s.label}/'.format(s=s),dpi=200)


# In[96]:


u.plot_all_gainratios(s,fp+'plots/')


# In[94]:


u.plt_corrs(s,fp+'plots/')


# In[ ]:




