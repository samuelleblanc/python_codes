
# coding: utf-8

# Program to check out previous flights from the P3 to calculate various speeds..

# In[4]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
#matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import numpy as np
#import Pysolar.solar as sol
import datetime
fp='C:/Users/sleblan2/Research/flight_planning/p3_flights/'


# In[7]:


from path_utils import getpath


# In[5]:


import load_utils as lm
reload(lm)


# In[6]:


from Sp_parameters import smooth


# In[9]:


fp = getpath('ORACLES')
fp


# # P3 during ORACLES 2018

# In[12]:


ora,orh = lm.load_ict(fp+'data_other_2018/Hskping_P3_20180927_R0.ict',return_header=True)


# In[29]:


plt.figure()
plt.plot(ora['Ground_Speed'],ora['GPS_Altitude'],'.',label='actual')
alts = np.arange(1000,5400,100)
plt.plot(alts*0.007+110,alts,'g-',lw=6,label='parameterized')
plt.plot(np.arange(5400,7000,100)*0.0+155.0,np.arange(5400,7000,100),'g-',lw=6)
plt.legend(frameon=False,loc=2)
plt.ylabel('Altitude [m]')
plt.xlabel('Ground Speed [m/s]')


# ## P3 during ARCTAS

# In[ ]:


arctas,header = lm.load_ict(fp+'pds_p3b_20080419_r2.ict',return_header=True)
header


# In[41]:


vert_speed_ft = np.diff(arctas['GPS_ALT'])
arctas['GPS_ALT']


# In[45]:


vert_speed = vert_speed_ft*0.3084


# In[55]:


plt.plot(smooth(vert_speed,10),smooth(arctas['GPS_ALT'][1:]*0.3084,10),'b.')
plt.xlabel('Vertical speed [m/s]')
plt.ylabel('Altitude [m]')
plt.title('P-3 vertical speed for ARCTAS')


# ## P3 during DISCOVER-AQ Denver

# In[58]:


discover = lm.load_ict(fp+'discoveraq-pds_p3b_20140807_r1.ict')


# In[59]:


d_vert_speed_ft = np.diff(discover['GPS_ALT'])


# In[60]:


d_vert_speed = d_vert_speed_ft*0.3084


# In[68]:


plt.plot(smooth(d_vert_speed,10),smooth(discover['GPS_ALT'][1:]*0.3084,10),'bo',label='DISCOVER-AQ:\n Boulder')
plt.plot(smooth(vert_speed,10),smooth(arctas['GPS_ALT'][1:]*0.3084,10),'r.',alpha=0.3,label='ARCTAS')
plt.xlabel('Vertical speed [m/s]')
plt.ylabel('Altitude [m]')
plt.title('P-3 vertical speed')
plt.legend(frameon=False,numpoints=1)
plt.savefig(fp+'P3_vert_speed.png',dpi=600,transparent=True)


# In[109]:


vs_up = smooth(vert_speed[(vert_speed>1)&(arctas['GPS_ALT'][1:]*0.3084>6000)],10)


# In[110]:


alt_up = smooth(arctas['GPS_ALT'][1:][(vert_speed>1)&(arctas['GPS_ALT'][1:]*0.3084>6000)],10)


# In[112]:


v,xx = linfit.linfit(alt_up,vs_up)


# In[114]:


p3_slope,p3_interp = v


# In[115]:


p3_slope


# In[118]:


alt1


# In[123]:


alt0 = 5000.0


# In[124]:


speed = 5.0


# In[125]:


climb_time = (alt1-alt0)/speed


# In[126]:


climb_time


# ## For ER2 during SEAC4RS

# In[4]:


er2,header = lm.load_ict(fp+'seac4rs-nasdat_er2_20130821_r0.ict',return_header=True)


# In[5]:


er2_2 = lm.load_ict(fp+'seac4rs-nasdat_er2_20130922_r0.ict',return_header=False)


# In[6]:


header


# In[7]:


er2_vs = np.diff(er2['GPS_Altitude'])


# In[8]:


er22_vs = np.diff(er2_2['GPS_Altitude'])


# In[9]:


import plotting_utils as pu


# In[10]:


plt.plot(smooth(er22_vs,20),smooth(er2_2['GPS_Altitude'][1:],20),'.',color='grey')
plt.plot(smooth(er2_vs,20),smooth(er2['GPS_Altitude'][1:],20),'b.')
plt.plot(smooth(er2_vs[er2_vs<0],20),smooth(er2['GPS_Altitude'][1:][er2_vs<0],20),'r.')
pu.plot_lin(smooth(er2_vs[er2_vs>2],20),smooth(er2['GPS_Altitude'][1:][er2_vs>2],20))
pu.plot_lin(smooth(er2_vs[er2_vs<0],20),smooth(er2['GPS_Altitude'][1:][er2_vs<0],20),color='r')
plt.xlabel('Vertical speed [m/s]')
plt.ylabel('Altitude [m]')
plt.title('ER2 vertical speed from SEAC4RS')
plt.legend(frameon=False)
plt.savefig(fp+'ER2_vert_speed.png',dpi=600,transparent=True)


# Get the inverse relationship for alt to vert speed

# In[11]:


import linfit
v = linfit.linfit(smooth(er2['GPS_Altitude'][1:][er2_vs>2],20),smooth(er2_vs[er2_vs>2],20))


# In[12]:


slope = v[0][0]
intercept = v[0][1]


# In[13]:


slope,intercept


# In[23]:


plt.figure()
plt.plot(smooth(er2['True_Air_Speed'],10),smooth(er2['GPS_Altitude'],10),'.',color='blue')
pu.plot_lin(smooth(er2['True_Air_Speed'],10),smooth(er2['GPS_Altitude'],10),color='blue')
plt.legend(frameon=False,loc=2)
plt.xlabel('True Airspeed [m/s]')
plt.ylabel('Altitude [m]')
plt.title('ER2 Airspeed from SEAC4RS')
plt.savefig(fp+'ER2_airspeed.png',dpi=600,transparent=True)


# In[24]:


vy = linfit.linfit(smooth(er2['GPS_Altitude'],10),smooth(er2['True_Air_Speed'],10))


# In[25]:


slopey,intercepty = vy[0]


# In[26]:


slopey,intercepty


# ## For DC8 during SEAC4RS

# In[127]:


dc8,dc8header = lm.load_ict(fp+'nav_dc8_20080320_r1.ict',return_header=True)


# In[128]:


dc8header


# In[132]:


dc8_vs = np.diff(dc8['GPS_ALT'])


# In[141]:


plt.plot(smooth(dc8_vs,10),smooth(dc8['GPS_ALT'][1:],10),'b.')
plt.plot(15.0-0.001*np.linspace(0,12000),np.linspace(0,12000),label='vs = 15-0.001*alt')
plt.title('DC8 during SEAC4RS')
plt.xlabel('Vertical speed [m/s]')
plt.ylabel('Altitude [m/s]')
plt.legend(frameon=False)
plt.savefig(fp+'DC8_vert_speed.png',dpi=600,transparent=True)


# ## C130 during ARISE

# In[142]:


c130,c130header = lm.load_ict(fp+'arise-C130-Hskping_c130_20140911_RA_Preliminary.ict',return_header=True)


# In[143]:


c130header


# In[147]:


plt.plot(smooth(c130['Vertical_Speed'],10),smooth(c130['GPS_Altitude'],10),'b.')
plt.plot(10-0.001*np.linspace(0,7500),np.linspace(0,7500),label='vs = 10-0.001*alt')
plt.title('C130 vertical speed during ARISE')
plt.xlabel('Vertical Speed [m/s]')
plt.ylabel('Altitude [m]')
plt.legend(frameon=False)
plt.savefig(fp+'C130_vert_speed.png',dpi=600,transparent=True)

