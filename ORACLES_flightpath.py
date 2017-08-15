
# coding: utf-8

# # Intro 
# 
# Codes to make some simple plots out of the housekeeping ict files

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth


# In[2]:

from linfit import linfit


# In[3]:

get_ipython().magic(u'matplotlib notebook')


# In[4]:

fp ='C:/Users/sleblan2/Research/ORACLES/'


# # Load some ict files

# In[5]:

fname = 'Hskping_P3_20170809_R0.ict'


# In[6]:

f = fp+'data_other_2017/'+fname


# In[11]:

h = load_ict(f)


# In[14]:

h.dtype.names


# # Plot out some figures of altitudes and lat lon

# In[19]:

plt.figure()
plt.plot(h['Start_UTC'],h['MSL_GPS_Altitude'])
plt.xlabel('UTC [h]')
plt.ylabel('GPS Altitude [m]')
plt.title('PRF00Y17')
plt.savefig(fp+'plot/2017/20170809_PRF00Y17_alt.png',dpi=600,transparent=True)


# In[47]:

kmh_2_knts = 1.94384449246 
base_speed = 110.0
speed_per_alt = 0.0075
alt = np.linspace(0,5500,55)


# In[48]:

sp = (base_speed+alt*speed_per_alt)*kmh_2_knts


# In[37]:

plt.figure()
plt.plot(h['True_Air_Speed'],h['MSL_GPS_Altitude'],'.')
plt.xlabel('True Air Speed [knts]')
plt.ylabel('GPS Altitude [m]')
plt.title('PRF00Y17')
plt.plot(sp,alt,'r',lw=7)
#plt.savefig(fp+'plot/2017/20170809_PRF00Y17_alt.png',dpi=600,transparent=True)


# In[43]:

np.diff(h['MSL_GPS_Altitude'])>0


# In[49]:

plt.figure()
plt.plot(h['True_Air_Speed'],h['MSL_GPS_Altitude'],'.')
plt.plot(h['True_Air_Speed'][np.diff(h['MSL_GPS_Altitude'])>0],h['MSL_GPS_Altitude'][np.diff(h['MSL_GPS_Altitude'])>0],'g.')
plt.plot(h['True_Air_Speed'][np.diff(h['MSL_GPS_Altitude'])<0],h['MSL_GPS_Altitude'][np.diff(h['MSL_GPS_Altitude'])<0],'r.')

plt.xlabel('True Air Speed [knts]')
plt.ylabel('GPS Altitude [m]')
plt.title('PRF00Y17')
plt.plot(sp,alt,'k',lw=7)


# In[40]:

plt.figure()
plt.hist(h['True_Air_Speed'],bins=40)


# In[ ]:



