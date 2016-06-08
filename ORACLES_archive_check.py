
# coding: utf-8

# # Intro
# Simple Program to load and check the flight planning archive files.
# 
# For RA of ORACLES

# # Load the defaults and imports

# In[2]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth


# In[3]:

get_ipython().magic(u'matplotlib notebook')


# In[15]:

fp ='C:/Users/sleblan2/Research/flight_planning/ORACLES_stm/'


# # Load the files

# In[16]:

fe,feh = load_ict(fp+'ORACLES-Flt-plan_er2_2016-09-01_RA.ict',return_header=True)


# In[17]:

fpp,fph = load_ict(fp+'ORACLES-Flt-plan_p3_2016-09-01_RA.ict',return_header=True)


# In[27]:

plt.figure()
plt.plot(fpp['Start_UTC'],fpp['AZI'])


# In[20]:

plt.figure()
plt.plot(fpp['Start_UTC'],fpp['Latitude'])


# In[8]:

plt.figure()
plt.plot(fp['Start_UTC'],fp['Longitude'])


# In[ ]:



