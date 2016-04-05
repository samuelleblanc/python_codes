
# coding: utf-8

# # Info
# Name:  
# 
#     NAAMES_cld_comp
# 
# Purpose:  
# 
#     To compare the various cloud properties retrieved via different methods from NAAMES.
#     Looking at the Nov. 17th, 2015 case
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
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - load_modis.py : for loading modis files
#     - matplotlib
#     - numpy
#     - scipy : for saving and reading
#     - math
#     - os
#     - gc
#     - pdb
#     - datetime
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - 20151117_zen_cld_retrieved.mat: cloud retrieval file
#   - MYD06_L2.A2015321.1540.006.2015322185040.hdf: MODIS file
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, NASA Ames, Santa Cruz, CA, 2016-03-24
#     Modified: 

# # Import initial modules and default paths

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.io as sio
import Sp_parameters as Sp
import hdf5storage as hs
import load_modis as lm


# In[16]:

from mpl_toolkits.basemap import Basemap,cm


# In[64]:

import write_utils as wu


# In[2]:

get_ipython().magic(u'matplotlib notebook')


# In[3]:

# set the basic directory path
fp='C:/Users/sleblan2/Research/NAAMES/'


# # Load the various data

# ## Load the retrieved cloud properties from 4STAR

# In[4]:

star = hs.loadmat(fp+'retrieve/20151117_zen_cld_retrieved.mat')


# In[5]:

star.keys()


# In[131]:

plt.figure()
plt.subplot(3,1,1)
plt.plot(star['utc'],star['tau'],'r.')
plt.plot(star['utc'],Sp.smooth(star['tau'],4,nan=False),'-xr')
plt.ylabel('$\\tau$')
plt.subplot(3,1,2)
plt.plot(star['utc'],star['ref'],'g.')
plt.plot(star['utc'],Sp.smooth(star['ref'],4,nan=False),'-xg')
plt.ylabel('r$_{eff}$ [$\\mu$m]')
plt.subplot(3,1,3)
plt.plot(star['utc'],star['phase'],'x')
plt.ylabel('phase')
plt.xlabel('UTC [h]')


# In[132]:

star['tau'] = Sp.smooth(star['tau'],4,nan=False)
star['ref'] = Sp.smooth(star['ref'],4,nan=False)


# ## Load the MODIS file

# In[13]:

myd3 = fp+'c130/20151117_flt/MYD03.A2015321.1540.006.2015322160400.hdf'
myd6 = fp+'c130/20151117_flt/MYD06_L2.A2015321.1540.006.2015322185040.hdf'


# In[14]:

modis,modis_dicts = lm.load_modis(myd3,myd6)


# ## Load the insitu probe liquid vs ice water content

# In[42]:

help(lm.load_ict)


# In[43]:

wf = fp+'c130/20151117_flt/NAAMES-LARGE-WCM_C130_20151117_RA.ict'
wcm,wcm_head = lm.load_ict(wf,return_header=True)


# In[57]:

plt.figure()
plt.plot(wcm['Date_UTC']/3600,Sp.smooth(wcm['TWC_gm_3'],30),label='TWC')
plt.plot(wcm['Date_UTC']/3600,Sp.smooth(wcm['LWC1_gm_3'],30),label='LWC')
plt.ylabel('Water content [g/m$^3$]')
plt.xlabel('UTC [h]')
plt.xlim([14.7,16.3])
plt.legend(frameon=False,loc=2)


# # Start plotting MODIS and the retrieved values

# In[31]:

#set up a easy plotting function
def naames_map(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=-40,lat_0=50,
            llcrnrlon=-42, llcrnrlat=48,
            urcrnrlon=-37, urcrnrlat=52,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(-35,-45,11),labels=[0,0,0,1])
    m.drawparallels(np.linspace(45,55,11),labels=[1,0,0,0])
    return m


# In[33]:

fig,ax = plt.subplots(1,1,figsize=(11,8))
m = naames_map(ax)
x,y = m(modis['lon'],modis['lat'])
clevels = np.linspace(0,80,41)

cs1 = m.contourf(x,y,modis['tau'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('$\\tau$')

xx,yy = m(star['lon'],star['lat'])
m.scatter(xx,yy,c=star['tau'],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
          alpha=0.5,edgecolors='k',linewidth=0.25)


# In[34]:

fig,ax = plt.subplots(1,1,figsize=(11,8))
m = naames_map(ax)
x,y = m(modis['lon'],modis['lat'])
clevels = np.linspace(0,60,31)

cs1 = m.contourf(x,y,modis['ref'],clevels,cmap=plt.cm.gist_earth,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('r$_{eff}$ [$\\mu$m]')

xx,yy = m(star['lon'],star['lat'])
m.scatter(xx,yy,c=star['ref'],cmap=plt.cm.gist_earth,marker='o',vmin=clevels[0],vmax=clevels[-1],
          alpha=0.5,edgecolors='k',linewidth=0.25)


# In[36]:

modis_dicts['phase']


# In[40]:

fig,ax = plt.subplots(1,1,figsize=(11,8))
m = naames_map(ax)
x,y = m(modis['lon'],modis['lat'])
clevels = np.linspace(0,6,7)

cs1 = m.contourf(x,y,modis['phase'],clevels,cmap=plt.cm.gist_earth,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('phase')

xx,yy = m(star['lon'],star['lat'])
m.scatter(xx,yy,c=star['phase'],cmap=plt.cm.rainbow,marker='o',vmin=0,vmax=1,
          alpha=0.5,edgecolors='k',linewidth=0.25)


# # Now prepare a subsection of retrieved values to be saved in ict

# ## create the dicts for writing to ict

# In[60]:

star.keys()


# In[151]:

d_dict =  {'utc':{'data':star['utc'],'unit':'hours from midnight UTC',
                     'long_description':'Fractional hours starting a midnight, continuous'},
          'COD':{'data':star['tau'],'unit':'None','format':'.1f',
                 'long_description':'Cloud Optical Depth of the cloud above the Aircraft'},
          'REF':{'data':star['ref'],'unit':'microns','format':'.1f',
                 'long_description':'Cloud particle effective radius, pertains to liquid cloud drops and ice crystals'},
          'PHASE':{'data':star['phase'],'unit':'None','format':'.0f',
                   'long_description':'Thermodynamic phase of cloud above, 0: pure liquid cloud, 1: pure ice cloud, mixed phase not retrieved'},
          'LAT':{'data':star['lat'],'unit':'Degrees','format':'.6f',
                 'long_description':'Aircraft position latitude (North positive)'},
          'LON':{'data':star['lon'],'unit':'Degrees','format':'.6f',
                 'long_description':'Aircraft position longitude (East positive)'},
          'ALT':{'data':star['alt'],'unit':'meter','format':'.1f',
                 'long_description':'Aircraft altitude'},
          'SZA':{'data':star['sza'],'unit':'Degrees','format':'.2f',
                 'long_description':'Solar Zenith Angle, angle of the sun between it and zenith'}
          }


# In[156]:

h_dict ={'PI':'Jens Redemann',
         'Institution':'NASA Ames Research Center',
         'Instrument':'4STAR',
         'campaign':'NAAMES #1',
         'special_comments':'Preliminary retrieved cloud properties data',
         'PI_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
         'platform':'C130',
         'location':"based out of St-John's, NL, Canada, actual location of C130 described by lat and lon below",
         'instrument_info':'Retrieved products from the 4STAR zenith radiance measurements',
         'data_info':'For references see LeBlanc et al.(2015) AMT, doi:10.5194/amt-8-1361-2015',
         'time_interval':10.0,
         'uncertainty':'Preliminary 7% in REF and 5% in COD',
         'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
         'project_info':'NAAMES field mission',
         'stipulations':'Prior OK from PI',
         'rev_comments':"""RA: preliminary retrieved values, may be subject to multiple errors
    including due to clouds influencing presumed surface albedo, non-homogeneous clouds, or mixed phase clouds"""
        }


# In[63]:

order=['LAT','LON','ALT','SZA','COD','REF','PHASE']


# ## Verify the input, plot and write the file. Subset only valid time.

# In[152]:

data_dict = wu.prep_data_for_ict(d_dict,Start_UTC=15.04,End_UTC=15.34,time_interval=10.0)


# In[143]:

data_dict.keys()


# In[154]:

wu.make_plots_ict(data_dict,filepath=fp+'plot/',data_id='4STAR-zen-cloud',loc_id='C130',date='20151117',rev='RA',
                  plot_together=['COD','REF','PHASE'],plot_together2=['LAT','LON','SZA'])


# In[157]:

wu.write_ict(h_dict,data_dict,filepath=fp,data_id='4STAR-zen-cloud',loc_id='C130',date='20151117',rev='RA',order=order)

