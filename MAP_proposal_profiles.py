
# coding: utf-8

# # Intro
# Name:  
# 
#     MAP_proposal_profiles
# 
# Purpose:  
# 
#     Python script for plotting AOD profile and extinctions for various places
# 
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures and save files...
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - os
#     - datetime
#     - mpl_toolkits
#     - plotting_utils (user defined plotting routines)
#     - map_utils, dependent on geopy
#     - Basemap
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, NASA Ames, 2016-06-08

# # Import the required modules and do the setup

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
get_ipython().magic(u'matplotlib notebook')
from mpltools import color
import numpy as np
import scipy.io as sio
import hdf5storage as hs
import math
import os
import Sp_parameters as Sp


# In[7]:

import load_utils as lu
from load_utils import mat2py_time, toutc


# In[16]:

from mpl_toolkits.basemap import Basemap,cm


# In[3]:

# set the basic directory path
fp='C:\\Users\\sleblan2\\Research\\SEAC4RS\\'


# # Load some required files

# ## Load some matlab files

# ### Load 4STAR data

# In[4]:

star = sio.loadmat(fp+'dc8\\20130816\\20130816starsun_R2.mat',variable_names=('w','tau_aero','t','Alt','Lat','Lon'))


# In[5]:

star['tt'] = mat2py_time(star['t'])
star['utc'] = toutc(star['tt'])


# In[6]:

star_cl = sio.loadmat(fp+'starsun\\20130913starsun_R2_tauaero.mat',variable_names=('w','tau_aero','t','Alt','Lat','Lon'))


# ## Load some ict archived files

# ### 4STAR data

# In[ ]:

istar_cld,istar_cldh = lu.load_ict(fp+'starsun_ict\\SEAC4RS-4STAR-AOD-CWV_DC8_20130913_R2.ict',return_header=True)


# In[35]:

istar_cldh


# In[9]:

istar_bnd,istar_bndh = lu.load_ict(fp+'starsun_ict\\SEAC4RS-4STAR-AOD-CWV_DC8_20130816_R2.ict',return_header=True)


# In[38]:

istar_cld2,istar_cld2h = lu.load_ict(fp+'starsun_ict\\SEAC4RS-4STAR-AOD-CWV_DC8_20130823_R2.ict',return_header=True)


# In[222]:

istar_val,istar_valh = lu.load_ict(fp+'starsun_ict\\SEAC4RS-4STAR-AOD-CWV_DC8_20130806_R2.ict',return_header=True)


# ### insitu data

# In[11]:

iso2_cld,iso2_cldh = lu.load_ict(fp+'dc8\\SEAC4RS-GTCIMS-SO2_DC8_20130913_R1.ict',return_header=True)


# In[12]:

iso2_bnd,iso2_bndh = lu.load_ict(fp+'dc8\\SEAC4RS-GTCIMS-SO2_DC8_20130816_R1.ict',return_header=True)


# In[13]:

iams_bnd,iams_bndh = lu.load_ict(fp+'dc8\\SEAC4RS-AMS_DC8_20130816_R1.ict',return_header=True)


# In[14]:

iams_cld,iams_cldh = lu.load_ict(fp+'dc8\\SEAC4RS-AMS_DC8_20130913_R1.ict',return_header=True)


# In[220]:

iso2_cld,iso2_cldh = lu.load_ict(fp+'dc8\\SEAC4RS-GTCIMS-SO2_DC8_20130823_R1.ict',return_header=True)


# In[221]:

iams_cld,iams_cldh = lu.load_ict(fp+'dc8\\SEAC4RS-AMS_DC8_20130823_R1.ict',return_header=True)


# In[223]:

iso2_val,iso2_valh = lu.load_ict(fp+'dc8\\SEAC4RS-GTCIMS-SO2_DC8_20130806_R1.ict',return_header=True)


# In[ ]:

iams_val,iams_valh = lu.load_ict(fp+'dc8\\SEAC4RS-AMS_DC8_20130806_R1.ict',return_header=True)


# # Now start plotting the files

# In[17]:

import map_interactive as mi


# ## Check out the cld case #1

# In[27]:

fig,ax = plt.subplots(1,1)
m = mi.build_basemap(lower_left=[-100,15],upper_right=[-85,32],ax=ax)
xt,yt = m(-95.3831,29.7628)
ax.text(xt,yt,'+')
ax.text(xt,yt,'Houston, TX',horizontalalignment='right',verticalalignment='top')
cs = m.scatter(istar_cld['Longitude'],istar_cld['Latitude'],marker='o',latlon=True,
               c=istar_cld['Start_UTC'],s=20,cmap=plt.cm.rainbow,edgecolor='none')
cb = plt.colorbar(cs)
cb.set_label('UTC [h]')


# In[30]:

it = (istar_cld['Start_UTC']>19)&(istar_cld['Start_UTC']<22)


# In[31]:

fig,ax = plt.subplots(1,1)
m = mi.build_basemap(lower_left=[-100,15],upper_right=[-85,32],ax=ax)
xt,yt = m(-95.3831,29.7628)
ax.text(xt,yt,'+')
ax.text(xt,yt,'Houston, TX',horizontalalignment='right',verticalalignment='top')
cs = m.scatter(istar_cld['Longitude'],istar_cld['Latitude'],marker='o',latlon=True,
               c=istar_cld['Start_UTC'],s=20,cmap=plt.cm.rainbow,edgecolor='none')
cb = plt.colorbar(cs)
cb.set_label('UTC [h]')
cs = m.scatter(istar_cld['Longitude'][it],istar_cld['Latitude'][it],marker='o',latlon=True,
               c=istar_cld['Start_UTC'][it],s=20,cmap=plt.cm.rainbow)


# In[36]:

plt.figure()
ig = istar_cld['qual_flag'][it]==0
plt.plot(istar_cld['AOD0501'][it],istar_cld['GPS_alt'][it],'.')
plt.plot(istar_cld['AOD0501'][it][ig],istar_cld['GPS_alt'][it][ig],'xr')
plt.xlabel('AOD0501')
plt.ylabel('Alitude [m]')


# In[37]:

plt.figure()
plt.plot(istar_cld['qual_flag'][it])


# ## Check out the cld case #2

# In[44]:

fig,ax = plt.subplots(1,1)
m = mi.build_basemap(lower_left=[-97,28],upper_right=[-90,35],ax=ax)
xt,yt = m(-95.3831,29.7628)
ax.text(xt,yt,'+')
ax.text(xt,yt,'Houston, TX',horizontalalignment='right',verticalalignment='top')
cs = m.scatter(istar_cld2['Longitude'],istar_cld2['Latitude'],marker='o',latlon=True,
               c=istar_cld2['Start_UTC'],s=20,cmap=plt.cm.gist_ncar,edgecolor='none')
cb = plt.colorbar(cs)
cb.set_label('UTC [h]')


# In[48]:

it = (istar_cld2['Start_UTC']>17)&(istar_cld2['Start_UTC']<20.5)


# In[49]:

plt.figure()
ig2 = istar_cld2['qual_flag'][it]==0
plt.plot(istar_cld2['AOD0501'][it],istar_cld2['GPS_alt'][it],'.')
plt.plot(istar_cld2['AOD0501'][it][ig2],istar_cld2['GPS_alt'][it][ig2],'xr')
plt.xlabel('AOD0501')
plt.ylabel('Alitude [m]')


# In[52]:

fig,ax = plt.subplots(1,1)
m = mi.build_basemap(lower_left=[-97,28],upper_right=[-90,35],ax=ax)
xt,yt = m(-95.3831,29.7628)
ax.text(xt,yt,'+')
ax.text(xt,yt,'Houston, TX',horizontalalignment='right',verticalalignment='top')
cs = m.scatter(istar_cld2['Longitude'],istar_cld2['Latitude'],marker='o',latlon=True,
               c=istar_cld2['Start_UTC'],s=20,cmap=plt.cm.gist_ncar,edgecolor='none')
m.scatter(istar_cld2['Longitude'][it][ig2],istar_cld2['Latitude'][it][ig2],marker='x',latlon=True,s=20)
cb = plt.colorbar(cs)
cb.set_label('UTC [h]')


# In[53]:

plt.figure()
ig2 = istar_cld2['qual_flag'][it]==0
plt.plot(istar_cld2['AOD0501'][it][ig2],istar_cld2['GPS_alt'][it][ig2],'xr')
plt.xlabel('AOD0501')
plt.ylabel('Alitude [m]')


# In[54]:

alt = istar_cld2['GPS_alt'][it][ig2]
aod = istar_cld2['AOD0501'][it][ig2]


# In[71]:

plt.figure()
plt.plot(alt,aod,'.')
plt.plot(bins,results,'+r')


# In[213]:

from Sp_parameters import ext_prof,smooth


# In[210]:

ext,bins = ext_prof(istar_cld2['GPS_alt'][it][ig2],istar_cld2['AOD0501'][it][ig2],binsize=50,verbose=True)


# In[219]:

plt.figure()
plt.plot(ext,bins,'.')


# ## Plot figure for the CA Valley case

# In[225]:

fig,ax = plt.subplots(1,1)
m = mi.build_basemap(lower_left=[-127,33],upper_right=[-115,45],ax=ax)
xt,yt = m(-95.3831,29.7628)
ax.text(xt,yt,'+')
ax.text(xt,yt,'Houston, TX',horizontalalignment='right',verticalalignment='top')
cs = m.scatter(istar_val['Longitude'],istar_val['Latitude'],marker='o',latlon=True,
               c=istar_val['Start_UTC'],s=20,cmap=plt.cm.gist_ncar,edgecolor='none')
#m.scatter(istar_cld2['Longitude'][it][ig2],istar_cld2['Latitude'][it][ig2],marker='x',latlon=True,s=20)
cb = plt.colorbar(cs)
cb.set_label('UTC [h]')


# In[ ]:



