
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

# In[2]:


import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
get_ipython().magic(u'matplotlib notebook')
#import path


# In[3]:


import load_utils as lu
from path_utils import getpath


# In[4]:


fp = getpath('ORACLES')


# In[3]:


fp = 'C:\\Users\\sleblan2\\Research\\ORACLES\\data_other\\'


# In[4]:


fp2 = 'C:\\Users\\sleblan2\\Research\\ORACLES\\data_other_2017\\'


# In[5]:


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


# ## Load 2018 Misamfu

# In[7]:


fp


# In[8]:


dm = lu.load_aeronet(fp+'data_other_2018/AERONET/20181219_20181219_Misamfu.lev15',version=3)


# In[12]:


km = dm.keys()
km.sort()


# In[13]:


km


# In[94]:


dm['Site_Elevationm'], dm['Site_LatitudeDegrees'], dm['Site_LongitudeDegrees']


# In[15]:


nm = [float(k[4:-2]) for k in km if (k.startswith('AOD_')&k.endswith('nm'))]
nm.sort()


# In[16]:


nm


# In[20]:


maods = [k for k in km if (k.startswith('AOD_')&k.endswith('nm'))]
nat_sort(maods)


# In[22]:


m = {}
m['nm'] = nm
m['keys'] = maods


# In[23]:


maod = []
for i in xrange(len(dm['Timehhmmss'])):
    u = [dm[d][i] for d in maods]
    maod.append(u)
m['aod'] = np.array(maod)


# In[24]:


m['aod'].shape


# In[31]:


m['aod'][m['aod']==-999.0]=np.nan


# In[33]:


plt.figure()
plt.plot(m['nm'],m['aod'].T,'x')


# In[93]:


for i,aa in enumerate(np.nanmean(m['aod'][6:9,:].T,axis=1)):
    if np.isfinite(aa):
        print m['nm'][i],aa


# ## Load the El Farafa 2018

# In[34]:


e = lu.load_aeronet(fp+'data_other_2018/AERONET/20180102_20180102_El_Farafra.lev20',version=3)


# In[63]:


e['Site_Elevationm'],e['Site_LatitudeDegrees'],e['Site_LongitudeDegrees']


# In[62]:


ke


# In[36]:


ke = e.keys()
ke.sort()
nm = [float(k[4:-2]) for k in ke if (k.startswith('AOD_')&k.endswith('nm'))]
nm.sort()


# In[37]:


eaods = [k for k in ke if (k.startswith('AOD_')&k.endswith('nm'))]
nat_sort(eaods)


# In[38]:


em = {}
em['nm'] = nm
em['keys'] = eaods


# In[39]:


eaod = []
for i in xrange(len(e['Timehhmmss'])):
    u = [e[d][i] for d in eaods]
    eaod.append(u)
em['aod'] = np.array(eaod)


# In[40]:


em['aod'][em['aod']==-999.0]=np.nan


# In[43]:


em['aod'].shape


# In[64]:


em['nm']


# In[60]:


em['aod'][-1,:]


# In[61]:


plt.figure()
plt.plot(em['nm'],em['aod'][-5:-1,:].T,'x')


# ## Load from Ragged point 2018

# In[75]:


r = lu.load_aeronet(fp+'data_other_2018/AERONET/20180923_20180923_Ragged_Point.lev15',version=3)


# In[76]:


kr = r.keys()
kr.sort()
nm = [float(k[4:-2]) for k in kr if (k.startswith('AOD_')&k.endswith('nm'))]
nm.sort()


# In[77]:


raods = [k for k in kr if (k.startswith('AOD_')&k.endswith('nm'))]
nat_sort(raods)


# In[78]:


rm = {}
rm['nm'] = nm
rm['keys'] = raods


# In[79]:


raod = []
for i in xrange(len(r['Timehhmmss'])):
    u = [r[d][i] for d in raods]
    raod.append(u)
rm['aod'] = np.array(raod)


# In[80]:


rm['aod'][rm['aod']==-999.0]=np.nan
rm['aod'].shape


# In[83]:


plt.figure()
plt.plot(rm['nm'],rm['aod'][0:5,:].T,'x')


# In[88]:


for i,aa in enumerate(np.nanmean(rm['aod'][0:1,:].T,axis=1)):
    if np.isfinite(aa):
        print rm['nm'][i],aa


# ## Load the La Parguera 2018 (near Barbados)

# In[66]:


p = lu.load_aeronet(fp+'data_other_2018/AERONET/20180926_20180926_La_Parguera.lev20',version=3)


# In[67]:


kp = p.keys()
kp.sort()
nm = [float(k[4:-2]) for k in kp if (k.startswith('AOD_')&k.endswith('nm'))]
nm.sort()


# In[68]:


paods = [k for k in kp if (k.startswith('AOD_')&k.endswith('nm'))]
nat_sort(paods)


# In[69]:


pm = {}
pm['nm'] = nm
pm['keys'] = paods


# In[70]:


paod = []
for i in xrange(len(p['Timehhmmss'])):
    u = [p[d][i] for d in paods]
    paod.append(u)
pm['aod'] = np.array(paod)


# In[71]:


pm['aod'][pm['aod']==-999.0]=np.nan
pm['aod'].shape


# In[74]:


plt.figure()
plt.plot(pm['nm'],pm['aod'][5:10,:].T,'x')


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


# # Test out the aeronet aerosol properties in sao tome for ORACLES 2017

# In[11]:


from json_tricks import dump, dumps, load, loads, strip_comments


# In[12]:


ff = '/mnt/c/Users/sleblanc/Research/ORACLES/aero_file_v2.txt'


# In[13]:


n = load(ff)


# In[14]:


n.keys()


# In[15]:


n['asy'].shape


# In[13]:


plt.figure()
plt.plot(n['wvl_arr'],n['asy'][0,:])


# In[22]:


aero_asy = [0.702247,0.580717,0.520858,0.504244]
aero_wvl = [440.0,674.0,870.0,1020.0]


# In[21]:


new_asy = [ 0.75  ,  0.71,  0.69,  0.645,  0.58,
         0.52,  0.512,  0.51,  0.49 ,  0.45,
         0.427843  ,  0.377843  ]


# In[87]:


plt.figure()
plt.plot(n['wvl_arr'],n['asy'][0,:],label='lut_v3_aero')
plt.plot(aero_wvl,aero_asy,label='Sao Tome Aug 18th, 2017')
plt.plot(n['wvl_arr'],new_asy,label='new lut_v4_aero')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Asymmetry parameter')
plt.savefig('/mnt/c/Users/sleblanc/Research/ORACLES/aero_asy_v4.png',dpi=600,transparent=True)


# In[89]:


plt.figure()
plt.plot(n['wvl_arr'],n['asy'][0,:],'x-',label='2016')
plt.plot(aero_wvl,aero_asy,'x-',label='Sao Tome Aug 18th, 2017')
plt.plot(n['wvl_arr'],new_asy,'x-',label='2017')
plt.gca().set_xscale('log')
plt.legend(frameon=False)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Asymmetry parameter')
plt.grid()
plt.savefig('/mnt/c/Users/sleblanc/Research/ORACLES/aero_asy_v4_yr.png',dpi=600,transparent=True)


# In[20]:


aero_ssa = [0.869100,0.863700,0.849200,0.840400]
aero_ssa2 = [0.899500,0.889300,0.865900,0.843400]


# In[19]:


new_ssa = [ 0.885  ,  0.88,  0.878 ,  0.875,  0.87 ,
        0.851,  0.841,  0.838 ,  0.82,  0.79110621,
        0.761106  ,  0.721106  ]


# In[69]:


plt.figure()
plt.plot(n['wvl_arr'],n['ssa'][0,:],'x-',label='lut_v3_aero')
plt.plot(aero_wvl,aero_ssa,'x-',label='Sao Tome Aug 18th, 2017')
plt.plot(aero_wvl,aero_ssa2,'x-',label='Sao Tome Aug 13th, 2017')
plt.plot(n['wvl_arr'],new_ssa,'x-',label='new lut_v4_aero')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('SSA')
plt.savefig('/mnt/c/Users/sleblanc/Research/ORACLES/aero_ssa_v4.png',dpi=600,transparent=True)


# In[71]:


import numpy as np


# In[17]:


aero_ext = np.array([0.703200,0.357900,0.227200,0.164300])/3.0


# In[18]:


new_ext = n['ext'][0,:]*1.4


# In[76]:


plt.figure()
plt.plot(n['wvl_arr'],n['ext'][0,:],'x-',label='lut_v3_aero')
plt.plot(aero_wvl,aero_ext,'x-',label='Sao Tome Aug 18th, 2017')
#plt.plot(aero_wvl,aero_ssa2,'x-',label='Sao Tome Aug 13th, 2017')
plt.plot(n['wvl_arr'],new_ext,'x-',label='new lut_v4_aero')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Extinction')
plt.savefig('/mnt/c/Users/sleblanc/Research/ORACLES/aero_ext_v4.png',dpi=600,transparent=True)


# In[9]:


ssa_4star = [0.820466265,0.803721687,0.757131325,0.733398795]
ssa_4star_std = [0.031832575,0.032310169,0.040917537,0.046946279]
wvl_4star = [440.0,675.0,870.0,995.0]


# In[37]:


fig,ax = plt.subplots(3,sharex=True,figsize=(7,5))
ax = ax.ravel()


ax[0].plot(n['wvl_arr'],n['asy'][0,:],'x-',label='2016')
ax[0].plot(aero_wvl,aero_asy,'x-',label='Sao Tome Aug 18th, 2017')
ax[0].plot(n['wvl_arr'],new_asy,'x-',label='2017')
ax[0].set_xscale('log')

#ax[0].set_xlabel('Wavelength [nm]')
ax[0].set_ylabel('Asymmetry\nparameter')
ax[0].grid()
ax[0].set_title('Optical properties of aerosol above cloud')

ax[1].plot(n['wvl_arr'],n['ssa'][0,:],'x-')
ax[1].plot(aero_wvl,aero_ssa,'x-')
ax[1].plot(n['wvl_arr'],new_ssa,'x-')
ax[1].errorbar(wvl_4star,np.array(ssa_4star)+0.03,yerr=(np.array(ssa_4star_std)+0.02),label='4STAR SSA from 2016')
#ax[0].set_xlabel('Wavelength [nm]')
ax[1].set_ylabel('Single Scattering\nAlbedo')
ax[1].grid()
ax[1].legend(frameon=False)

ax[2].plot(n['wvl_arr'],n['ext'][0,:],'x-',label='2016')
ax[2].plot(aero_wvl,aero_ext,'x-',label='Sao Tome 2017-08-18')
ax[2].plot(n['wvl_arr'],new_ext,'x-',label='2017')
ax[2].set_ylabel('Extinction\nCoefficient [km$^{{-1}}$]')
ax[2].grid()
ax[2].legend(frameon=False)
ax[2].set_xlabel('Wavelength [nm]')
ax[2].set_xlim(250,5000)
ax[2].set_xticks([320,400,500,600,750,1000,1240,1600,2500,4000])
ax[2].set_xticklabels([320,400,500,600,750,1000,1240,1600,2500,4000])


plt.savefig('/mnt/c/Users/sleblanc/Research/ORACLES/aero_v4_yr.png',dpi=600,transparent=True)


# In[78]:


n['asy'][0,:] = new_asy
n['ssa'][0,:] = new_ssa
n['ext'][0,:] = new_ext


# In[79]:


help(dump)


# In[80]:


dump(n,'/mnt/c/Users/sleblanc/Research/ORACLES/aero_file_v4.txt')

