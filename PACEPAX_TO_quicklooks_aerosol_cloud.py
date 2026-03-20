#!/usr/bin/env python
# coding: utf-8

# # Intro
# Purpose:
# 
#     Make some plots of the PACE-PAX aerosol and cloud from the Twin Otter (TO)
# 
# Input:
# 
#     None
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
#   - load_utils.py
#   - matplotlib
#   - numpy
#   - write_utils
#   - path_utils
# 
# Needed Files:
#   - NA
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2025-07-29
#     Modified:

# # Prepare python environment

# In[78]:


import numpy as np
#import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import os
import pandas as pd
from datetime import datetime


# In[79]:


name = 'PACEPAX_TO'
vv = 'v0'
fp = getpath(name)


# # Load files

# In[81]:


f = os.listdir(fp)
f.sort()
f


# In[82]:


metnav = {}
for fn in f:
    if (fn.endswith('.ict') & ('MetNav' in fn)): 
        daystr = fn.split('_')[2] 
        if len(fn.split('_'))>4: daystr = daystr +'_'+ fn.split('_')[4].split('.')[0]
        print('doing : '+daystr)

        metnav[daystr] = lu.load_ict(os.path.join(fp,fn))


# In[83]:


aer = {}
for fn in f:
    if (fn.endswith('.ict') & ('MICROPHYSICAL' in fn)): 
        daystr = fn.split('_')[2] 
        if len(fn.split('_'))>4: daystr = daystr +'_'+ fn.split('_')[4].split('.')[0]
        print('doing : '+daystr)

        aer[daystr] = lu.load_ict(os.path.join(fp,fn))


# # Plot out the different values

# ## Interpolate the alt, lat, lon to the aerosol values

# In[84]:


from scipy.interpolate import interp1d
def inter(xold,yold,xnew):
    igood = np.isfinite(xold) & np.isfinite(yold)
    fx = interp1d(xold[igood],yold[igood],fill_value="extrapolate")
    return fx(xnew)


# In[85]:


import numpy.lib.recfunctions as rfn


# In[86]:


aern = {}
for k in list(aer.keys()):
    if not k in metnav: continue
    print(k)
    aern[k] = {}
    for ko in ['Latitude', 'Longitude', 'GPS_Altitude']:
        aern[k][ko] = inter(metnav[k]['Time_Start'],metnav[k][ko],aer[k]['Time_Mid'].data)
        #aer[k] = rfn.append_fields(aer[k],ko,inter(metnav[k]['Time_Start'],metnav[k][ko],aer[k]['Time_Mid'].data),dtypes='<f8')


# ## plot the aerosol number per alt

# In[87]:


aer[k]


# In[88]:


for k in list(aern.keys()):
    plt.figure()
    plt.plot(aer[k]['IntegN_100to1000nm_UHSAS'],aern[k]['GPS_Altitude'],'.',color=aer[k]['Time_Mid'])
    plt.xlabel('UHSAS aerosol #')
    plt.ylabel('GPS Altitude [m]')
    plt.title('TO Aerosol for:'+k)


# ## plot the Cloud LWC per alt against time

# In[89]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[90]:


for k in list(metnav.keys()):
    plt.figure(figsize=(6,2))
    plt.plot(metnav[k]['Time_Start']/3600.0,metnav[k]['GPS_Altitude'],lw=0.4)#'.',color='gray')
    mask = metnav[k]['LWC_Wire'] > 0.05
    p = plt.scatter(metnav[k]['Time_Start'][mask]/3600.0,metnav[k]['GPS_Altitude'][mask],1,marker='o',c=metnav[k]['LWC_Wire'][mask], vmin=0.05, vmax=0.25)
    plt.xlabel('Time UTC [h]')
    plt.ylabel('GPS Altitude [m]')
    plt.title('TO Cloud LWC for:'+k)
    plt.colorbar(label='LWC [g/m^3] from CAPS')
    plt.grid()
    


# In[91]:


for k in list(metnav.keys()):
    plt.figure(figsize=(6,2))
    plt.plot(metnav[k]['Longitude'],metnav[k]['GPS_Altitude'],lw=0.4)#'.',color='gray')
    mask = metnav[k]['LWC_Wire'] > 0.05
    p = plt.scatter(metnav[k]['Longitude'][mask],metnav[k]['GPS_Altitude'][mask],1,marker='o',c=metnav[k]['LWC_Wire'][mask], vmin=0.05, vmax=0.25)
    plt.xlabel('Longitude [°]')
    plt.ylabel('GPS Altitude [m]')
    plt.title('TO Cloud LWC for:'+k)
    plt.colorbar(label='LWC [g/m^3] from CAPS')
    plt.grid()
    


# In[92]:


for k in list(metnav.keys()):
    plt.figure(figsize=(6,2))
    plt.plot(metnav[k]['Latitude'],metnav[k]['GPS_Altitude'],lw=0.4)#'.',color='gray')
    mask = metnav[k]['LWC_Wire'] > 0.05
    p = plt.scatter(metnav[k]['Latitude'][mask],metnav[k]['GPS_Altitude'][mask],1,marker='o',c=metnav[k]['LWC_Wire'][mask], vmin=0.05, vmax=0.25)
    plt.xlabel('Latitude [°]')
    plt.ylabel('GPS Altitude [m]')
    plt.title('TO Cloud LWC for:'+k)
    plt.colorbar(label='LWC [g/m^3] from CAPS')
    plt.grid()


# ## Plot out the cloud and aerosol time trace

# In[95]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[93]:


import matplotlib.colors as mcolors


# In[101]:


for k in list(metnav.keys()):
    plt.figure(figsize=(10,3))
    plt.plot(metnav[k]['Time_Start']/3600.0,metnav[k]['GPS_Altitude'],lw=0.4,zorder=3)#'.',color='gray')
    mask = metnav[k]['LWC_Wire'] > 0.05
    p = plt.scatter(metnav[k]['Time_Start'][mask]/3600.0,metnav[k]['GPS_Altitude'][mask],20,cmap='GnBu',
                    marker='o',c=metnav[k]['LWC_Wire'][mask], vmin=0.05, vmax=0.25,zorder=2)
    plt.colorbar(label='LWC from CAPS [g/m^3]',extend="max")
    plt.xlabel('Time UTC [h]')
    plt.ylabel('GPS Altitude [m]')
    plt.title('TO Cloud LWC and aerosol number for:'+k)
    plt.grid()
    try:
        fl = aer[k]['IntegN_100to1000nm_UHSAS'] > 50
        norm = mcolors.LogNorm(vmin=1, vmax=50000)
        pa = plt.scatter(aer[k]['Time_Mid'][fl]/3600.0,aern[k]['GPS_Altitude'][fl],40,marker='^',
                         c=aer[k]['IntegN_100to1000nm_UHSAS'][fl],cmap='PuRd',zorder=1,vmin=0,vmax=1000)
        plt.colorbar(label='UHSAS aerosol\n100-1000nm [#]',extend="max")
    except KeyError:
        pass
    plt.tight_layout()
    if k=='20240903':
        #plt.xlim(19.5,21)
        break
    
    #plt.savefig(fp+f'PACE-PAX_TO_LWC_aerosolnum_{k}.png',transparent=True,dpi=600)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




