
# coding: utf-8

# # Info
# Name:  
# 
#     Explor_AERONET
# 
# Purpose:  
# 
#     Explore AERONET values for single sites linked to high altitude AOD values
#   
# Input:
# 
#     none at command line
#     see methods of module
# 
# Output:
#    
#     plots
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - math
#     - pdb
#     - datetime
#     - load_utils
#   
# Needed Files:
# 
#     - AEERONET files
#   
#   
# Modification History:
# 
#     Wrtten: Samuel LeBlanc, NASA Ames, from Santa Cruz, 2017-05-19
#     Modified: 

# # Import the required python modules and set paths

# In[32]:

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
get_ipython().magic(u'matplotlib notebook')


# In[33]:

import load_utils as lu


# In[34]:

fp = 'C:\\Users\\sleblan2\\Research\\ORACLES\\data_other\\'


# In[35]:

fp2 = 'C:\\Users\\sleblan2\\Research\\ORACLES\\data_other_2017\\'


# In[36]:

def nat_sort(l):
    n = max([len(j) for j in l])+1
    l.sort(key=lambda x: '{0:0>{n}}'.format(x,n=n).lower())


# # Load the files

# ## Load the Aeronet Mount_Chacaltaya

# In[5]:

a = lu.load_aeronet(fp+'150708_150708_Mount_Chacaltaya.lev20')


# In[6]:

a.keys()


# In[80]:

len(a['Timehhmmss'])


# In[29]:

nm = [float(k[4:]) for k in a.keys() if k.startswith('AOT_')]
nm.sort()


# In[31]:

a['wl'] = nm


# In[38]:

aods = [k for k in a.keys() if k.startswith('AOT_')]


# In[40]:

aods.sort()


# In[43]:

aods.sort(key=lambda x: '{0:0>8}'.format(x).lower())


# In[44]:

aods


# In[86]:

aod = []
for i in xrange(len(a['Timehhmmss'])):
    u = [a[d][i] for d in aods]
    aod.append(u)
a['aod'] = np.array(aod)


# ## Load the MLO AERONET

# In[14]:

reload(lu)


# In[15]:

b = lu.load_aeronet(fp+'20160801_20161031_Mauna_Loa.lev15',version=3)


# In[16]:

b.keys()


# In[35]:

nm = [float(k[4:-2]) for k in b.keys() if (k.startswith('AOD_')&k.endswith('nm'))]
nm.sort()


# In[36]:

nm


# In[37]:

b['wl'] = nm


# In[49]:

baods = [k for k in b.keys() if (k.startswith('AOD_')&k.endswith('nm'))]


# In[64]:

nat_sort(baods)


# In[78]:

baods


# In[88]:

baod = []
for i in xrange(len(b['Timehhmmss'])):
    u = [b[d][i] for d in baods]
    baod.append(u)
b['aod'] = np.array(baod)


# In[91]:

b['aod'].shape


# ## Load the file for Sao Tomé

# In[44]:

c = lu.load_aeronet(fp2+'170801_170930_Sao_Tome.lev10',version=2)


# In[45]:

c.keys()


# In[46]:

from datetime import datetime, timedelta


# In[47]:

list(c['Julian_Day'])


# In[48]:

dt = [datetime(2017,1,1)+timedelta(jd-1) for jd in list(c['Julian_Day'])]


# In[49]:

plt.figure()
plt.plot(dt,c['AOT_500'],'.')
plt.grid()


# ## Load the Ascension island aeronet data (from the airport)

# In[50]:

d = lu.load_aeronet(fp2+'170801_170930_Ascension_Island.lev15',version=2)


# In[51]:

d['dt'] = [datetime(2017,1,1)+timedelta(jd-1) for jd in list(d['Julian_Day'])]


# In[53]:

plt.figure()
plt.plot(d['dt'],d['AOT_500'],'.')
plt.grid()


# In[ ]:




# # Plot out the AOD spectra for each point

# ## Plot out the Mount_Chacaltaya aeronet

# In[169]:

fig,ax = plt.subplots(1)
s_m = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar, norm=plt.normalize(vmin=min(a['Julian_Day']),vmax=max(a['Julian_Day'])))
s_m.set_array([])
for i,j in enumerate(a['Julian_Day']):
    im = np.isfinite(a['aod'][i,:])
    plt.plot(np.array(a['wl'])[im],a['aod'][i,im],'x-',color=s_m.to_rgba(j))
plt.plot(a['wl'],np.nanmean(a['aod'],axis=0),'o-k',zorder=200,lw=4)
plt.plot(a['wl'],np.nanmedian(a['aod'],axis=0),'s-',color='grey',zorder=200,lw=4)
plt.title('Mount_Chacaltaya AERONET')
plt.xlabel('Wavelenght [nm]')
plt.ylabel('AOD')
plt.grid()
cb = plt.colorbar(s_m, format='%3.2f')
cb.set_label('DOY')
plt.savefig(fp+'AERONET_Chacaltaya_high_alt_AOD_log_sept_2016.png',dpi=600,transparent=True)


# In[170]:

fig,ax = plt.subplots(1)
s_m = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar, norm=plt.normalize(vmin=min(a['Julian_Day']),vmax=max(a['Julian_Day'])))
s_m.set_array([])
for i,j in enumerate(a['Julian_Day']):
    im = np.isfinite(a['aod'][i,:])
    plt.plot(np.array(a['wl'])[im],a['aod'][i,im],'x-',color=s_m.to_rgba(j))
plt.plot(a['wl'],np.nanmean(a['aod'],axis=0),'o-k',zorder=200,lw=4)
plt.plot(a['wl'],np.nanmedian(a['aod'],axis=0),'s-',color='grey',zorder=200,lw=4)
plt.title('Mount_Chacaltaya AERONET')
plt.xlabel('Wavelenght [nm]')
plt.ylabel('AOD')
plt.grid()

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.001,0.1)
ax.set_xlim(300.0,1700.0)
plt.xticks([350,400,500,600,800,1000,1200,1400,1650])
plt.yticks([0.001,0.002,0.005,0.01,0.015,0.02,0.03,0.05,0.1])
ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

cb = plt.colorbar(s_m, format='%3.2f')
cb.set_label('DOY')
plt.savefig(fp+'AERONET_Chacaltaya_high_alt_AOD_log_sept_2016.png',dpi=600,transparent=True)


# ## Plot out the MLO AERONET

# In[166]:

fig,ax = plt.subplots(1)
s_m = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar, norm=plt.normalize(vmin=min(b['Day_of_Year']),vmax=max(b['Day_of_Year'])))
s_m.set_array([])
for i,j in enumerate(b['Day_of_Year']):
    b['aod'][i,b['aod'][i,:]==-999] = np.nan
    im = np.isfinite(b['aod'][i,:])
    plt.plot(np.array(b['wl'])[im],b['aod'][i,im],'x-',color=s_m.to_rgba(j))
plt.plot(b['wl'],np.nanmean(b['aod'],axis=0),'o-k',zorder=200,lw=4)
plt.title('MLO AERONET')
plt.xlabel('Wavelenght [nm]')
plt.ylabel('AOD')
cb = plt.colorbar(s_m, format='%3.2f')
cb.set_label('DOY')
plt.savefig(fp+'AERONET_MLO_high_alt_AOD_sept_2016.png',dpi=600,transparent=True)


# In[165]:

fig,ax = plt.subplots(1)
s_m = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar, norm=plt.normalize(vmin=min(b['Day_of_Year']),vmax=max(b['Day_of_Year'])))
s_m.set_array([])
for i,j in enumerate(b['Day_of_Year']):
    b['aod'][i,b['aod'][i,:]==-999] = np.nan
    im = np.isfinite(b['aod'][i,:])
    plt.plot(np.array(b['wl'])[im],b['aod'][i,im],'x-',color=s_m.to_rgba(j))
plt.plot(b['wl'],np.nanmean(b['aod'],axis=0),'o-k',zorder=200,lw=4)
plt.plot(b['wl'],np.nanmedian(b['aod'],axis=0),'s-',color='grey',zorder=200,lw=4)
plt.title('MLO AERONET')
plt.xlabel('Wavelenght [nm]')
plt.ylabel('AOD')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.001,0.1)
ax.set_xlim(300.0,1700.0)
plt.xticks([350,400,500,600,800,1000,1200,1400,1650])
plt.yticks([0.001,0.002,0.005,0.01,0.015,0.02,0.03,0.05,0.1])
ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
plt.grid()

cb = plt.colorbar(s_m, format='%3.2f')
cb.set_label('DOY')
plt.savefig(fp+'AERONET_MLO_high_alt_AOD_log_sept_2016.png',dpi=600,transparent=True)


# In[162]:

np.array(b['wl'])[im],np.nanmean(b['aod'][:,im],axis=0),np.nanmedian(b['aod'][:,im],axis=0)


# In[172]:

plt.figure()
plt.plot(b['Day_of_Year'],b['aod'][:,7],'x')


# In[ ]:



