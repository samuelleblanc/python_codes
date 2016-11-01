
# coding: utf-8

# # Intro
# Name:  
# 
#     Explore_cld_retrieval
# 
# Purpose:  
# 
#     Run throught the retrieved cloud properties and either flag or assure retrieval quality
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
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - plotting_utils (user defined plotting routines)
#     - hdf5storage
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - 4STAR_cloud retrieval .mat files
#   
#  Modification History:
#  
#      Written: by Samuel LeBlanc, NASA Ames, Moffett Field, CA, 2016-10-26

# # Import of modules

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
get_ipython().magic(u'matplotlib notebook')
import numpy as np
import scipy.io as sio
import hdf5storage as hs
import Sp_parameters as Sp


# In[2]:

# set the basic directory path
fp = 'C:/Users/sleblan2/Research/ORACLES/starzen/'
fp_plot = 'C:/Users/sleblan2/Research/ORACLES/plot/'


# # Load the files

# In[20]:

dds = ['20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927']


# In[21]:

rts = []
sps = []


# In[22]:

for daystr in dds:
    print daystr
    rt = hs.loadmat(fp+'{}_zen_cld_retrieved.mat'.format(daystr))
    s = sio.loadmat(fp+'4STAR_{}starzen.mat'.format(daystr))
    sp = Sp.Sp(s)
    rts.append(rt)
    sps.append(sp)


# ## Load the cloud probe incloud flag

# In[105]:

from load_utils import mat2py_time,toutc


# In[23]:

p = sio.netcdf_file(fp+'..//data_other//oracles.cloud.timings.nc','r')


# In[24]:

p.variables


# In[25]:

p.variables['timevec_20160914'].data


# In[17]:

t_0914 = mat2py_time(p.variables['timevec_20160914'].data)


# In[19]:

plt.figure()
plt.plot(t_0914,p.variables['cloud_time_20160914'].data,'x')


# # Start plotting the results

# In[26]:

rt.keys()


# In[27]:

plt.figure()
plt.plot(rt['utc'],rt['tau'])


# In[36]:

rt = rts[9]


# In[33]:

plt.figure()
plt.plot(rts[9]['utc'],rts[9]['tau'],'.')
plt.plot(rts[9]['utc'],rts[9]['utc'],'r+')


# In[35]:

plt.figure()
plt.plot(rts[9]['tau'],rts[9]['ref'],'.')


# In[30]:

igood = rts[9]['tau']>0


# In[38]:

igood[0:10]


# In[39]:

sp = sps[9]


# In[69]:

i=68
i_vis = [1061,1062,1064]
i_nir = [1060,1063]
plt.figure()
plt.plot(sp.wvl,sp.norm[i,:])
#plt.xlim(970,1030)
plt.plot(sp.wvl[i_vis],sp.norm[i,i_vis],'rx')
plt.plot(sp.wvl[i_nir],sp.norm[i,i_nir],'g+')


# In[48]:

np.nanmean(sp.norm[i,iw])


# In[49]:

np.nanmean(sp.norm[i,ii])


# # Now setup filters to weed out bad data

# ## Filter out data points where nir and vis spectrometers don't match

# In[51]:

i_vis = [1061,1062,1064]
i_nir = [1060,1063]


# In[82]:

for i,daystr in enumerate(dds):
    nvis = np.nanmean(sps[i].norm[:,i_vis],axis=1)
    nnir = np.nanmean(sps[i].norm[:,i_nir],axis=1)
    rts[i]['delta'] = abs(nvis-nnir)
    rts[i]['fl_match'] = rts[i]['delta']<0.06
    print daystr,rts[i]['delta'].shape,rts[i]['delta'][rts[i]['fl_match']].shape,        float(rts[i]['delta'][rts[i]['fl_match']].shape[0])/ float(rts[i]['delta'].shape[0])*100.0


# ## Now filter out the times which were at too high altitude

# In[85]:

fl_alt = rt['alt']<1000.0


# In[103]:

for i,daystr in enumerate(dds):
    rts[i]['fl_alt'] = rts[i]['alt'][:,0]<1000.0
    print daystr,rts[i]['utc'].shape,rts[i]['utc'][rts[i]['fl_alt']].shape,        float(rts[i]['utc'][rts[i]['fl_alt']].shape[0])/ float(rts[i]['utc'].shape[0])*100.0


# ## Filter for in cloud

# In[107]:

from write_utils import nearest_neighbor


# In[120]:

for i,daystr in enumerate(dds):
    try:
        p_time = mat2py_time(p.variables['timevec_{}'.format(daystr)].data)
    except KeyError: # no in cloud data, so choose all of them
        rts[i]['fl_incld'] = rts[i]['utc']>0.0
        continue
    putc = toutc(p_time)
    rts[i]['incld'] = nearest_neighbor(putc,p.variables['cloud_time_{}'.format(daystr)].data,rts[i]['utc'],dist=1.0/3600)
    rts[i]['fl_incld'] = rts[i]['incld']==0
    print daystr,rts[i]['utc'].shape,rts[i]['utc'][rts[i]['fl_incld']].shape,        float(rts[i]['utc'][rts[i]['fl_incld']].shape[0])/ float(rts[i]['utc'].shape[0])*100.0


# ## Filter for high ki squared residuas

# In[126]:

for i,daystr in enumerate(dds):
    rts[i]['fl_ki'] = rts[i]['ki']<0.6
    print daystr,rts[i]['utc'].shape,rts[i]['utc'][rts[i]['fl_ki']].shape,        float(rts[i]['utc'][rts[i]['fl_ki']].shape[0])/ float(rts[i]['utc'].shape[0])*100.0


# ## Combine the filters

# In[127]:

for i,daystr in enumerate(dds):
    rts[i]['fl'] = rts[i]['fl_match'] & rts[i]['fl_alt'] & rts[i]['fl_incld'] & rts[i]['fl_ki']
    print daystr,rts[i]['utc'].shape,rts[i]['utc'][rts[i]['fl']].shape,        float(rts[i]['utc'][rts[i]['fl']].shape[0])/ float(rts[i]['utc'].shape[0])*100.0 


# # Now plot each retrieved product, filtered

# In[129]:

from Sp_parameters import smooth


# In[132]:

for i,daystr in enumerate(dds):
    plt.figure()
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212,sharex=ax1)
    ax1.plot(rts[i]['utc'],rts[i]['tau'],'b.')
    ax1.plot(rts[i]['utc'][rts[i]['fl']],rts[i]['tau'][rts[i]['fl']],'g+')
    ax1.plot(rts[i]['utc'][rts[i]['fl']],smooth(rts[i]['tau'][rts[i]['fl']],6),'kx')
    ax1.set_ylabel('tau')
    
    ax2.plot(rts[i]['utc'],rts[i]['ref'],'b.')
    ax2.plot(rts[i]['utc'][rts[i]['fl']],rts[i]['ref'][rts[i]['fl']],'g+')
    ax2.plot(rts[i]['utc'][rts[i]['fl']],smooth(rts[i]['ref'][rts[i]['fl']],6),'kx')
    ax2.set_ylabel('ref')
    ax2.set_xlabel('UTC')
    ax1.set_title(daystr)


# In[ ]:

for i,daystr in enumerate(dds):
    rts[i]['tau_fl'] = smooth(rts[i]['tau'][rts[i]['fl']],6)
    rts[i]['ref_fl'] = smooth(rts[i]['ref'][rts[i]['fl']],6)


# In[134]:

rt.keys()


# # Now write these values to ict file

# In[133]:

import write_utils as wu


# In[ ]:

i=9
d_dict = {'Start_UTC':{'data':rts[i]['uct'][rts[i]['fl']]*3600.0,'unit':'seconds from midnight UTC','long_description':'time keeping'},
      'COD':{'data':rts[i]['tau_fl'],'unit':'None','long_description':'Cloud optical Depth of overlying cloud'},
      'REF':{'data':rts[i]['ref_fl'],'unit':'micrometer','long_description':'Cloud drop effective radius for liquid clouds'},
      'LAT':{'data':rts[i]['lat'][rts[i]['fl']],'unit':'Degrees','long_description':'Latitude of measurement, negative for Southern hemisphere'},
      'LON':{'data':rts[i]['lon'][rts[i]['fl']],'unit':'Degrees','long_description':'Longitude of measurement, East is positive, from -180 to 180'}
      }
print d_dict
hdict = {'PI':'Samuel LeBlanc',
     'Institution':'NASA Ames Research Center',
     'Instrument':'4STAR (Spectrometer for Sky-Scanning, Sun Tracking Atmospheric Research)',
     'campaign':'ORACLES-2016',
     'special_comments':'Retrieved cloud properties',
     'PI_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
     'platform':'NASA P3',
     'location':'based out of Walvis Bay, Namibia, actual location of measurement included in file',
     'instrument_info':'Derived product from 4STAR zenith measurements',
     'data_info':'',
     'uncertainty':'Undefined',
     'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
     'project_info':'ORACLES project deployed during August-September, 2016',
     'stipulations':'None',
     'rev_comments':"""  R0: Preliminary archival of cloud properties retreived from 4STAR sky radiance measurements. Final radiance calibration not yet applied."""
    }
print hdict
order = ['X1','X2','X3']
write_ict(hdict,d_dict,filepath='C:/Users/sleblan2/Research/NAAMES/',
          data_id='4STAR_test',loc_id='C130',date='20160402',rev='RA',order=order)    

