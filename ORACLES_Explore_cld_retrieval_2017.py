
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


# In[4]:

# set the basic directory path
fp = 'C:/Users/sleblan2/Research/ORACLES/starzen_2017/'
fp_plot = 'C:/Users/sleblan2/Research/ORACLES/starzen_2017/plot/'


# In[5]:

vr = 'R0'


# # Load the files

# In[ ]:

dds = ['20170809','20170812','20170813','20170815','20170817','20170818','20170819','20170821','20170824','20170826']


# In[6]:

dds = ['20170809','20170812','20170813','20170815','20170817','20170818']


# In[7]:

rts = []
sps = []


# In[8]:

for daystr in dds:
    print daystr
    rt = hs.loadmat(fp+'{}_zen_cld_retrieved.mat'.format(daystr))
    s = sio.loadmat(fp+'4STAR_{}starzen.mat'.format(daystr))
    sp = Sp.Sp(s)
    rts.append(rt)
    sps.append(sp)


# ## Load the cloud probe incloud flag

# In[9]:

from load_utils import mat2py_time,toutc


# In[8]:

p = sio.netcdf_file(fp+'..//data_other//oracles.cloud.timings.nc','r')


# In[9]:

p.variables


# In[10]:

p.variables['timevec_20160914'].data


# In[11]:

t_0914 = mat2py_time(p.variables['timevec_20160914'].data)


# In[28]:

plt.figure()
plt.plot(t_0914,p.variables['cloud_time_20160914'].data,'x')


# # Start plotting the results

# In[10]:

rt.keys()


# In[11]:

plt.figure()
plt.plot(rt['utc'],rt['tau'])


# In[12]:

rt = rts[9]


# In[13]:

plt.figure()
plt.plot(rts[9]['utc'],rts[9]['tau'],'.')
plt.plot(rts[9]['utc'],rts[9]['utc'],'r+')


# In[16]:

plt.figure()
plt.plot(rts[9]['tau'],rts[9]['ref'],'.')


# In[12]:

igood = rts[9]['tau']>0


# In[13]:

igood[0:10]


# In[14]:

sp = sps[9]


# In[15]:

i=68
i_vis = [1061,1062,1064]
i_nir = [1060,1063]
plt.figure()
plt.plot(sp.wvl,sp.norm[i,:])
#plt.xlim(970,1030)
plt.plot(sp.wvl[i_vis],sp.norm[i,i_vis],'rx')
plt.plot(sp.wvl[i_nir],sp.norm[i,i_nir],'g+')


# In[16]:

np.nanmean(sp.norm[i,iw])


# In[17]:

np.nanmean(sp.norm[i,ii])


# ## Plot some of the sza for each day to ensure good fitting of lut

# In[14]:

plt.figure()
plt.plot(sps[7].utc,sps[7].sza,'x-')


# # Now setup filters to weed out bad data

# ## Filter out data points where nir and vis spectrometers don't match

# In[15]:

i_vis = [1061,1062,1064]
i_nir = [1060,1063]


# In[16]:

for i,daystr in enumerate(dds):
    nvis = np.nanmean(sps[i].norm[:,i_vis],axis=1)
    nnir = np.nanmean(sps[i].norm[:,i_nir],axis=1)
    rts[i]['delta'] = abs(nvis-nnir)
    rts[i]['fl_match'] = rts[i]['delta']<0.06
    print daystr,rts[i]['delta'].shape,rts[i]['delta'][rts[i]['fl_match']].shape,        float(rts[i]['delta'][rts[i]['fl_match']].shape[0])/ float(rts[i]['delta'].shape[0])*100.0


# ## Now filter out the times which were at too high altitude

# In[17]:

fl_alt = rt['alt']<1000.0


# In[18]:

for i,daystr in enumerate(dds):
    rts[i]['fl_alt'] = rts[i]['alt'][:,0]<1000.0
    print daystr,rts[i]['utc'].shape,rts[i]['utc'][rts[i]['fl_alt']].shape,        float(rts[i]['utc'][rts[i]['fl_alt']].shape[0])/ float(rts[i]['utc'].shape[0])*100.0


# ## Filter for in cloud

# In[19]:

from write_utils import nearest_neighbor


# In[20]:

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

# In[21]:

for i,daystr in enumerate(dds):
    rts[i]['fl_ki'] = rts[i]['ki']<0.6
    print daystr,rts[i]['utc'].shape,rts[i]['utc'][rts[i]['fl_ki']].shape,        float(rts[i]['utc'][rts[i]['fl_ki']].shape[0])/ float(rts[i]['utc'].shape[0])*100.0


# ## Combine the filters

# In[22]:

tot=0
tot_fl=0
for i,daystr in enumerate(dds):
    rts[i]['fl'] = rts[i]['fl_match'] & rts[i]['fl_alt'] & rts[i]['fl_ki'] #& rts[i]['fl_incld']
    print daystr,rts[i]['utc'].shape,rts[i]['utc'][rts[i]['fl']].shape,        float(rts[i]['utc'][rts[i]['fl']].shape[0])/ float(rts[i]['utc'].shape[0])*100.0 
    tot = tot+len(rts[i]['utc'])
    tot_fl = tot_fl+len(rts[i]['utc'][rts[i]['fl']])


# In[23]:

print tot, tot_fl, float(tot_fl)/float(tot)*100.0


# # Now plot each retrieved product, filtered

# In[24]:

from Sp_parameters import smooth


# In[25]:

for i,daystr in enumerate(dds):
    plt.figure()
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212,sharex=ax1)
    ax1.plot(rts[i]['utc'],rts[i]['tau'],'b.')
    ax1.plot(rts[i]['utc'][rts[i]['fl']],rts[i]['tau'][rts[i]['fl']],'g+')
    try:
        ax1.plot(rts[i]['utc'][rts[i]['fl']],smooth(rts[i]['tau'][rts[i]['fl']],6),'kx')
    except:
        pass
    ax1.set_ylabel('tau')
    
    ax2.plot(rts[i]['utc'],rts[i]['ref'],'b.')
    ax2.plot(rts[i]['utc'][rts[i]['fl']],rts[i]['ref'][rts[i]['fl']],'g+')
    try:
        ax2.plot(rts[i]['utc'][rts[i]['fl']],smooth(rts[i]['ref'][rts[i]['fl']],6),'kx')
    except:
        pass
    ax2.set_ylabel('ref')
    ax2.set_xlabel('UTC')
    ax1.set_title(daystr)


# In[26]:

for i,daystr in enumerate(dds):
    try:
        rts[i]['tau_fl'] = smooth(rts[i]['tau'][rts[i]['fl']],6)
        rts[i]['ref_fl'] = smooth(rts[i]['ref'][rts[i]['fl']],6)
    except:
        print 'except',i
        rts[i]['tau_fl'] = rts[i]['tau'][rts[i]['fl']]
        rts[i]['ref_fl'] = rts[i]['ref'][rts[i]['fl']]
    rts[i]['lat_fl'] = rts[i]['lat'][rts[i]['fl']]
    rts[i]['lon_fl'] = rts[i]['lon'][rts[i]['fl']]
    rts[i]['alt_fl'] = rts[i]['alt'][rts[i]['fl']]
    rts[i]['utc_fl'] = rts[i]['utc'][rts[i]['fl']]


# In[27]:

rt.keys()


# # Now write these values to ict file

# In[28]:

import write_utils as wu


# In[29]:

hdict = {'PI':'Jens Redemann',
     'Institution':'NASA Ames Research Center',
     'Instrument':'Spectrometers for Sky-Scanning, Sun-Tracking Atmospheric Research (4STAR)',
     'campaign':'ORACLES 2017',
     'special_comments':'Retrieved cloud properties',
     'PI_contact':'Jens.Redemann-1@nasa.gov',
     'platform':'NASA P3',
     'location':'based out of Sao Tome, actual location of measurement included in file',
     'instrument_info':'Derived product from 4STAR zenith measurements',
     'data_info':'Using the cloud property retrieval method based on spectral transmitted light measurements described by LeBlanc, Pileskie, Schmidt, and Coddington (2015), AMT, modified to include impact of overlying aerosol layer.',
     'uncertainty':'Uncertainty of retrieved properties will be defined in future archival.',
     'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
     'project_info':'ORACLES 2017 deployment; August-September 2017; Sao Tome',
     'stipulations':'Prior OK from PI is requested before using this data',
     'rev_comments':"""R0: Preliminary archival of cloud properties retrieved from 4STAR sky radiance measurements. Final radiance calibration not yet applied. Filtered out bad measurements, and high clouds. 
    Data is subject to uncertainties linked to detector stability, transfer efficiency of light through fiber optic cable, and deposition on the instrument window."""
    }
order = ['LAT','LON','COD','REF']


# In[30]:

for i,daystr in enumerate(dds):
    d_dict = {'Start_UTC':{'data':rts[i]['utc'][rts[i]['fl']]*3600.0,'unit':'seconds from midnight UTC','long_description':'time keeping'},
              'utc':{'data':rts[i]['utc'][rts[i]['fl']],'unit':'seconds from midnight UTC','long_description':'time keeping'},
          'COD':{'data':rts[i]['tau_fl'],'unit':'None','long_description':'Cloud Optical Depth of overlying cloud'},
          'REF':{'data':rts[i]['ref_fl'],'unit':'micrometer','long_description':'Cloud drop effective radius for liquid clouds'},
          'LAT':{'data':rts[i]['lat'][rts[i]['fl']],'unit':'Degrees','long_description':'Latitude of measurement, negative for Southern hemisphere'},
          'LON':{'data':rts[i]['lon'][rts[i]['fl']],'unit':'Degrees','long_description':'Longitude of measurement, East is positive, from -180 to 180'}
          }
    d_dict_out = wu.prep_data_for_ict(d_dict,in_var_name='utc',out_var_name='Start_UTC', in_input=True,time_interval=1.0)
    wu.write_ict(hdict,d_dict_out,filepath=fp+'..//zen_ict_2017/v1/',
              data_id='4STAR_CLD',loc_id='P3',date=daystr,rev='R0',order=order)    


# ## For use of this python, save values to mat files

# In[31]:

rtss = {str(i):rr for i,rr in enumerate(rts)}


# In[32]:

def dict_keys_to_unicode(d):
    out = dict()
    for k, v in d.items():
        out[k.decode()] = v
    return out

for n in rtss.keys():
    if type(rtss[n]) is list:
        print n
        for i,t in enumerate(rtss[n]):
            rtss[n][i] = dict_keys_to_unicode(t)
    else:
        print 'no',n
        rtss[n] = dict_keys_to_unicode(rtss[n])


# In[33]:

hs.savemat(fp+'..//zen_ict_2017/v1/{}_all_retrieved.mat'.format(vr),rtss)


# ## Optionally load the saved mat files

# In[23]:

rtss = hs.loadmat(fp+'..//zen_ict/v3/{}_all_retrieved.mat'.format(vr))


# In[34]:

if not 'rts' in locals():
    rts = []
    for n in sorted([int(u) for u in rtss.keys()]):
        rts.append(rtss[str(n)])
elif not rts:
    for n in sorted([int(u) for u in rtss.keys()]):
        rts.append(rtss[str(n)])


# # Make plots

# ## Read the files as a verification

# In[34]:

vv = 'R0'


# In[35]:

from load_utils import load_ict


# In[36]:

out_RA = []
out_head_RA = []
for d in dds:
    fname_aod = fp+'..//zen_ict_2017/v1/4STAR-CLD_P3_{}_{vr}.ict'.format(d,vr=vr)
    tt,th = load_ict(fname_aod,return_header=True)
    out_RA.append(tt)
    out_head_RA.append(th)


# In[37]:

out_head_RA[0]


# In[38]:

nm = out_RA[0].dtype.names


# In[39]:

nm


# In[40]:

for i,d in enumerate(dds):
    fig,ax = plt.subplots(2,sharex=True,figsize=(9,5))
    ax = ax.ravel()
    ax[0].set_title('Cloud variables {} saved file for flight {}'.format(vv,d),y=1.25)
    #ax[0].set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, len(wl))])
    ax[0].plot(out_RA[i][nm[0]],out_RA[i]['COD'],'.')
    ax[0].set_ylabel('COD')
    ax[0].set_ylim(0,60)
    ax[0].axhline(0,color='k')
    ax[0].grid()
    axy0 = ax[0].twiny()
    axy0.set_xlim(ax[0].get_xlim())
    xt = ax[0].get_xticks()
    xl = []
    for x in xt:
        ii = np.argmin(abs(out_RA[i][nm[0]]-x))
        if np.isfinite(out_RA[i]['LAT'][ii]):
            xl.append('{:2.2f}'.format(out_RA[i]['LAT'][ii]))
        else:
            ia = np.isfinite(out_RA[i]['LAT'][ii-1200:ii+1200])
            if any(ia):
                laa = np.interp([1200],np.arange(2400)[ia],out_RA[i]['LAT'][ii-1200:ii+1200][ia])
                if not np.isfinite(laa[0]):
                    xl.append(' ')
                else:
                    xl.append('{:2.2f}'.format(laa[0]))
            else: xl.append(' ')
    axy0.set_xticks(xt)
    axy0.set_xticklabels(xl)
    axy0.set_xlabel('Latitude [$^\\circ$]')
    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width, box.height*0.88])
    axy0.set_position([box.x0, box.y0, box.width, box.height*0.88])
    
    ax[1].plot(out_RA[i][nm[0]],out_RA[i]['REF'],'g.')
    ax[1].set_ylabel('r$_{{eff}}$ [$\\mu$m]')
    ax[1].set_xlabel('UTC [h]')
    ax[1].grid()
    axy1 = ax[1].twiny()
    axy1.set_xlim(ax[1].get_xlim())
    x1t = ax[1].get_xticks()
    x1l = []
    for x in x1t:
        ii = np.argmin(abs(out_RA[i][nm[0]]-x))
        if np.isfinite(out_RA[i]['LON'][ii]):
            x1l.append('{:2.2f}'.format(out_RA[i]['LON'][ii]))
        else:
            iio = np.isfinite(out_RA[i]['LON'][ii-1200:ii+1200])
            if any(iio):
                loo = np.interp([1200],np.arange(2400)[iio],out_RA[i]['LON'][ii-1200:ii+1200][iio])
                if not np.isfinite(loo[0]):
                    x1l.append(' ')
                else:
                    x1l.append('{:2.2f}'.format(loo[0]))
            else: x1l.append(' ')
    axy1.set_xticks(x1t)
    axy1.set_xticklabels(x1l)
    axy1.set_xlabel('Longitude [$^\\circ$]')
    box = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0, box.width, box.height*0.88])
    axy1.set_position([box.x0, box.y0, box.width, box.height*0.88])
    plt.savefig(fp+'..//zen_ict_2017/v1/{vv}_{}.png'.format(d,vv=vv),dpi=600,transparent=True)


# ## Combine the data into a single array

# In[41]:

ar = {}
for n in rts[0].keys():
    ar[n] = np.array([])


# In[42]:

ar['days'] = np.array([])


# In[43]:

for i,d in enumerate(dds):
    ar['days'] = np.append(ar['days'],np.zeros_like(rts[i]['utc'])+i)
    for n in rts[0].keys():
        ar[n] = np.append(ar[n],rts[i][n])


# ## Save the combined array

# In[44]:

import hdf5storage as hs


# In[45]:

hs.savemat(fp+'..//zen_ict_2017/v1/{}_all_cld_ict.mat'.format(vr),ar)


# ## Optionally load the all ict file

# In[7]:

if not 'ar' in locals():
    ar = hs.loadmat(fp+'..//zen_ict/v3/{}_all_cld_ict.mat'.format(vr))


# ## plot the data on a map

# In[46]:

import plotting_utils as pu


# In[47]:

from map_interactive import build_basemap


# In[48]:

rts[i]['tau_fl']


# In[49]:

for i,daystr in enumerate(dds):
    print rts[i]['lat'][rts[i]['fl']][:,0].shape,rts[i]['lon'][rts[i]['fl']][:,0].shape,rts[i]['tau_fl'].shape


# In[50]:

fig = plt.figure()
ax = plt.subplot(111)
m = build_basemap(lower_left=[-2,-17],upper_right=[15,0],ax=ax,larger=False)
sa = []
for i,daystr in enumerate(dds):
    x,y = m(rts[i]['lon'][rts[i]['fl']][:,0]+i*0.03,rts[i]['lat'][rts[i]['fl']][:,0])
    sca = ax.scatter(x,y,c=rts[i]['tau_fl'],
              s=10,alpha=0.7,vmin=0.0,vmax=60.0,edgecolor='None')
    sa.append(sca)
#pu.prelim()
cb = plt.colorbar(sa[0])
cb.set_label('COD')
plt.savefig(fp+'..//zen_ict_2017/v1/{}_COD_map.png'.format(vr),transparent=True,dpi=600)


# In[51]:

fig = plt.figure()
ax = plt.subplot(111)
m = build_basemap(lower_left=[-2,-17],upper_right=[15,0],ax=ax,larger=False)
sa = []
for i,daystr in enumerate(dds):
    x,y = m(rts[i]['lon'][rts[i]['fl']][:,0]+i*0.03,rts[i]['lat'][rts[i]['fl']][:,0])
    sca = ax.scatter(x,y,c=rts[i]['ref_fl'],
              s=10,alpha=0.7,vmin=0.0,vmax=30.0,edgecolor='None',cmap=plt.cm.gist_earth)
    sa.append(sca)
#pu.prelim()
cb = plt.colorbar(sa[0])
cb.set_label('r$_{{eff}}$ [$\\mu$m]')
plt.savefig(fp+'..//zen_ict_2017/v1/{}_REF_map.png'.format(vr),transparent=True,dpi=600)


# ## Plot out some statistics of all retrievals

# In[52]:

plt.figure()
plt.plot(ar['lat_fl'],ar['tau_fl'],'.',color='grey',alpha=0.1)
plt.hist2d(ar['lat_fl'],ar['tau_fl'],bins=40,normed=True)
plt.xlabel('Latitude [$^\\circ$]')
plt.ylabel('COD')
cb = plt.colorbar()
cb.set_label('Normalized counts')
plt.title('4STAR Cloud optical depth for all ORACLES flights')
plt.savefig(fp+'..//zen_ict_2017/v1/{}_COD_hist_lat.png'.format(vr),transparent=True,dpi=600)


# In[53]:

plt.figure()
plt.plot(ar['lon_fl'],ar['tau_fl'],'.',color='grey',alpha=0.1)
plt.hist2d(ar['lon_fl'],ar['tau_fl'],bins=40,normed=True)
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('COD')
cb = plt.colorbar()
cb.set_label('Normalized counts')
plt.title('4STAR Cloud optical depth for all ORACLES flights')
plt.savefig(fp+'..//zen_ict_2017/v1/COD_hist_lon.png',transparent=True,dpi=600)


# In[54]:

plt.figure()
plt.plot(ar['lon_fl'],ar['ref_fl'],'.',color='grey',alpha=0.1)
plt.hist2d(ar['lon_fl'],ar['ref_fl'],bins=40,normed=True,cmap=plt.cm.gist_earth)
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('r$_{{eff}}$ [$\\mu$m]')
plt.ylim(2,10)
cb = plt.colorbar()
cb.set_label('Normalized counts')
plt.title('4STAR Effective Radius for all ORACLES flights')
plt.savefig(fp+'..//zen_ict_2017/v1/ref_hist_lon.png',transparent=True,dpi=600)


# In[55]:

plt.figure()
plt.plot(ar['lat_fl'],ar['ref_fl'],'.',color='grey',alpha=0.1)
plt.hist2d(ar['lat_fl'],ar['ref_fl'],bins=40,normed=True,cmap=plt.cm.gist_earth)
plt.ylim(2,10)
plt.xlabel('Latitude [$^\\circ$]')
plt.ylabel('r$_{{eff}}$ [$\\mu$m]')
cb = plt.colorbar()
cb.set_label('Normalized counts')
plt.title('4STAR Effective Radius for all ORACLES flights')
plt.savefig(fp+'..//zen_ict_2017/v1/ref_hist_lat.png',transparent=True,dpi=600)


# In[56]:

fig = plt.figure()
plt.hist(ar['tau_fl'],bins=30,edgecolor='None',color='g',alpha=0.7,normed=True,label='filtered')
plt.hist(ar['tau'],bins=30,edgecolor='None',color='b',alpha=0.1,normed=True,range=(0,70),label='All points')
plt.ylabel('Normed counts')
plt.xlabel('COD')
plt.grid()
pu.prelim()
plt.legend(frameon=False)
plt.savefig(fp+'..//zen_ict_2017/v1/cod_hist.png',transparent=True,dpi=600)


# In[57]:

ar.keys()


# In[58]:

aam = ar['utc_fl']<12.0
apm = ar['utc_fl']>12.0


# In[59]:

fig = plt.figure()
plt.hist(ar['tau_fl'],bins=30,edgecolor='None',color='g',alpha=0.3,normed=True,label='filtered')
plt.hist(ar['tau'],bins=30,edgecolor='None',color='b',alpha=0.1,normed=True,range=(0,70),label='All points')
plt.hist(ar['tau_fl'][aam],bins=30,edgecolor='None',color='b',alpha=0.7,normed=True,label='AM')
plt.hist(ar['tau_fl'][apm],bins=30,edgecolor='None',color='r',alpha=0.7,normed=True,label='PM')
plt.ylabel('Normed counts')
plt.xlabel('COD')
plt.grid()
pu.prelim()
plt.legend(frameon=False)
plt.savefig(fp+'..//zen_ict_2017/v1/cod_hist_pm_am.png',transparent=True,dpi=600)


# In[60]:

fig = plt.figure()
plt.hist(ar['tau_fl'],bins=30,edgecolor='None',color='g',alpha=0.7,normed=False,label='filtered')
plt.hist(ar['tau'],bins=30,edgecolor='None',color='b',alpha=0.1,normed=False,range=(0,70),label='All points')
plt.ylabel('Counts')
plt.xlabel('COD')
plt.legend(frameon=False)
plt.savefig(fp+'..//zen_ict_2017/v1/cod_hist_all.png',transparent=True,dpi=600)


# In[61]:

np.nanmean(ar['tau_fl'])


# In[62]:

np.nanmean(ar['ref_fl'])


# In[63]:

fig = plt.figure()
plt.hist(ar['ref_fl'],bins=30,edgecolor='None',color='grey',alpha=0.7,normed=True,label='filtered')
plt.hist(ar['ref'],bins=30,edgecolor='None',color='b',alpha=0.1,normed=True,range=(0,30),label='all points')
plt.ylabel('Normed counts')
plt.xlabel('r$_{{eff}}$ [$\\mu$m]')
plt.grid()
#pu.prelim()
plt.legend(frameon=False)
plt.savefig(fp+'..//zen_ict_2017/v1/{}_ref_hist.png'.format(vr),transparent=True,dpi=600)


# In[65]:

fig = plt.figure()
plt.hist(ar['ref_fl'],bins=30,edgecolor='None',color='grey',alpha=0.7,normed=False,label='filtered')
plt.hist(ar['ref'],bins=30,edgecolor='None',color='b',alpha=0.1,normed=False,range=(0,30),label='all points')
plt.ylabel('Counts')
plt.xlabel('r$_{{eff}}$ [$\\mu$m]')
plt.legend(frameon=False)
plt.savefig(fp+'..//zen_ict_2017/v1/{}_ref_hist_all.png'.format(vr),transparent=True,dpi=600)


# In[66]:

reload(pu)


# In[67]:

fig,ax = plt.subplots(2,1)
ax = ax.ravel()
ax[0].hist(ar['tau_fl'],bins=30,edgecolor='None',color='g',alpha=0.7,normed=True,label='filtered')
ax[0].hist(ar['tau'],bins=30,edgecolor='None',color='b',alpha=0.1,normed=True,range=(0,70),label='all points')
ax[0].set_ylabel('Normed counts')
ax[0].set_xlabel('COD')
ax[0].grid()
pu.prelim(ax=ax[0])
ax[0].legend(frameon=False)

ax[1].hist(ar['ref_fl'],bins=30,edgecolor='None',color='grey',alpha=0.7,normed=True,label='filtered')
ax[1].hist(ar['ref'],bins=30,edgecolor='None',color='b',alpha=0.1,normed=True,range=(0,30),label='all points')
ax[1].set_ylabel('Normed counts')
ax[1].set_xlabel('r$_{{eff}}$ [$\\mu$m]')
plt.grid()
pu.prelim(ax=ax[1])
plt.legend(frameon=False)

plt.tight_layout()

plt.savefig(fp+'..//zen_ict_2017/v1/{}_ref_cod_hist.png'.format(vr),transparent=True,dpi=600)


# # Evaluate the Cloud Radiative Effect (CRE) from calculated retrieved values

# Based on the calculations of CRE found in Link to [ORACLES_cld_CRE](ORACLES_cld_CRE.ipynb)
# 
# After running calculations on Pleaides, results are read in and operated

# ## Load results

# In[37]:

fp


# In[38]:

c = hs.loadmat(fp+'../rtm/ORACLES_CRE_{}.mat'.format('v2'))


# In[39]:

c.keys()


# In[40]:

c['star_aero_C']


# ## Start plotting results of CRE

# In[41]:

import plotting_utils as pu


# In[42]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-1000,0],xlim=[-0.5,1.5])
ax1.set_ylabel('Cloud Radiative Effect [W/m$^2$]')
ax1.set_title('Cloud Radiative Effect from 4STAR retrievals\nSurface')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Aerosols','Without Aerosols'])
pu.plot_vert_hist(fig,ax1,c['star_aero_C'][:,0],0,[-1000,0],legend=True,onlyhist=False,loc=2,color='g',bins=30)
pu.plot_vert_hist(fig,ax1,c['star_noaero_C'][:,0],1,[-1000,0],legend=True,color='r',bins=30)
plt.savefig(fp+'../plot/ORACLES_surface_CRE_4STAR.png',transparent=True,dpi=600)


# In[43]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-1000,0],xlim=[-0.5,1.5])
ax1.set_ylabel('Cloud Radiative Effect [W/m$^2$]')
ax1.set_title('Cloud Radiative Effect from 4STAR retrievals\nTop of Atmosphere')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Aerosols','Without Aerosols'])
pu.plot_vert_hist(fig,ax1,c['star_aero_C'][:,2],0,[-1000,0],legend=True,onlyhist=False,loc=2,color='g',bins=30)
pu.plot_vert_hist(fig,ax1,c['star_noaero_C'][:,2],1,[-1000,0],legend=True,color='r',bins=30)
plt.savefig(fp+'../plot/ORACLES_CRE_toa_4STAR.png',transparent=True,dpi=600)


# In[48]:

print 'Surface CRE'
print 'mean aero: {}, no aero: {}'.format(np.nanmean(c['star_aero_C'][:,0]),np.nanmean(c['star_noaero_C'][:,0]))
print 'median aero: {}, no aero: {}'.format(np.nanmedian(c['star_aero_C'][:,0]),np.nanmedian(c['star_noaero_C'][:,0]))
print 'std aero: {}, no aero: {}'.format(np.nanstd(c['star_aero_C'][:,0]),np.nanstd(c['star_noaero_C'][:,0]))


# In[49]:

print 'TOA CRE'
print 'mean aero: {}, no aero: {}'.format(np.nanmean(c['star_aero_C'][:,2]),np.nanmean(c['star_noaero_C'][:,2]))
print 'median aero: {}, no aero: {}'.format(np.nanmedian(c['star_aero_C'][:,2]),np.nanmedian(c['star_noaero_C'][:,2]))
print 'std aero: {}, no aero: {}'.format(np.nanstd(c['star_aero_C'][:,2]),np.nanstd(c['star_noaero_C'][:,2]))


# In[50]:

plt.figure()
plt.hist(c['star_aero_CRE']['up'][:,2]/c['star_aero_CRE']['dn'][:,2],normed=False,edgecolor='None',color='g',
         alpha=0.6,label='With Aerosol')
plt.hist(c['star_noaero_CRE']['up'][:,2]/c['star_noaero_CRE']['dn'][:,2],normed=False,edgecolor='None',color='r',
         alpha=0.6,label='Without Aerosol')
plt.xlabel('Broadband SW albedo TOA')
plt.legend(frameon=False)
plt.title('TOA albedo')
plt.savefig(fp+'../plot/ORACLES_albedo_toa_4STAR.png',transparent=True,dpi=600)


# ## calculate and plot the relative CRE

# In[62]:

star_aero_rC = np.zeros_like(c['star_aero_C'])
star_noaero_rC = np.zeros_like(c['star_aero_C'])


# In[88]:

star_aero_rC[:,0] = c['star_aero_C'][:,0]/c['star_aero_CRE']['dn'][:,2]*100.0
star_aero_rC[:,1] = c['star_aero_C'][:,1]/c['star_aero_CRE']['dn'][:,2]*100.0
star_aero_rC[:,2] = c['star_aero_C'][:,2]/c['star_aero_CRE']['dn'][:,2]*100.0
star_noaero_rC[:,0] = c['star_noaero_C'][:,0]/c['star_noaero_CRE']['dn'][:,2]*100.0
star_noaero_rC[:,1] = c['star_noaero_C'][:,1]/c['star_noaero_CRE']['dn'][:,2]*100.0
star_noaero_rC[:,2] = c['star_noaero_C'][:,2]/c['star_noaero_CRE']['dn'][:,2]*100.0


# In[91]:

star_aero_rC_abc = np.zeros_like(c['star_aero_C'])
star_noaero_rC_abc = np.zeros_like(c['star_aero_C'])

star_aero_rC_abc[:,0] = c['star_aero_C'][:,0]/c['star_aero_CRE']['dn'][:,1]*100.0
star_noaero_rC_abc[:,0] = c['star_noaero_C'][:,0]/c['star_noaero_CRE']['dn'][:,1]*100.0


# In[89]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-100,0],xlim=[-0.5,1.5])
ax1.set_ylabel('relative Cloud Radiative Effect [\%]')
ax1.set_title('relative Cloud Radiative Effect from 4STAR retrievals\nSurface')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Aerosols','Without Aerosols'])
pu.plot_vert_hist(fig,ax1,star_aero_rC[:,0],0,[-100,0],legend=True,onlyhist=False,loc=4,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star_noaero_rC[:,0],1,[-100,0],legend=True,color='r',bins=50)
plt.savefig(fp+'../plot/ORACLES_rCRE_surface_4STAR.png',transparent=True,dpi=600)


# In[90]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-100,0],xlim=[-0.5,1.5])
ax1.set_ylabel('relative Cloud Radiative Effect [\%]')
ax1.set_title('relative CRE to above cloud from 4STAR retrievals\nSurface')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Aerosols','Without Aerosols'])
pu.plot_vert_hist(fig,ax1,c['star_aero_C'][:,0]/c['star_aero_CRE']['dn'][:,1]*100.0,0,[-100,0],legend=True,onlyhist=False,loc=4,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,c['star_noaero_C'][:,0]/c['star_noaero_CRE']['dn'][:,1]*100.0,1,[-100,0],legend=True,color='r',bins=50)
plt.savefig(fp+'../plot/ORACLES_rCRE_above_cloud_for_surface_4STAR.png',transparent=True,dpi=600)


# In[84]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-100,0],xlim=[-0.5,1.5])
ax1.set_ylabel('relative Cloud Radiative Effect [\%]')
ax1.set_title('relative Cloud Radiative Effect from 4STAR retrievals\nTop of Atmosphere')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Aerosols','Without Aerosols'])
pu.plot_vert_hist(fig,ax1,star_aero_rC[:,2],0,[-100,0],legend=True,onlyhist=False,loc=4,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star_noaero_rC[:,2],1,[-100,0],legend=True,color='r',bins=50)
plt.savefig(fp+'../plot/ORACLES_rCRE_toa_4STAR.png',transparent=True,dpi=600)


# In[85]:

print 'Surface rCRE'
print 'mean aero: {}, no aero: {}'.format(np.nanmean(star_aero_rC[:,0]),np.nanmean(star_noaero_rC[:,0]))
print 'median aero: {}, no aero: {}'.format(np.nanmedian(star_aero_rC[:,0]),np.nanmedian(star_noaero_rC[:,0]))
print 'std aero: {}, no aero: {}'.format(np.nanstd(star_aero_rC[:,0]),np.nanstd(star_noaero_rC[:,0]))


# In[86]:

print 'TOA rCRE'
print 'mean aero: {}, no aero: {}'.format(np.nanmean(star_aero_rC[:,2]),np.nanmean(star_noaero_rC[:,2]))
print 'median aero: {}, no aero: {}'.format(np.nanmedian(star_aero_rC[:,2]),np.nanmedian(star_noaero_rC[:,2]))
print 'std aero: {}, no aero: {}'.format(np.nanstd(star_aero_rC[:,2]),np.nanstd(star_noaero_rC[:,2]))


# ## plot the aerosol forcing

# In[68]:

c.keys()


# In[69]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-100,20],xlim=[-0.5,1.5])
ax1.set_ylabel('Direct Aerosol Radiative Effect [W/m$^2$]')
ax1.set_title('Direct Aerosol Radiative Effect from 4STAR retrievals\nSurface')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Clouds','Without Clouds'])
DAREs = (c['star_aero_CRE']['dn'][:,0]-c['star_aero_CRE']['up'][:,0])-(c['star_noaero_CRE']['dn'][:,0]-c['star_noaero_CRE']['up'][:,0])
DAREs_clear = (c['star_aero_CRE_clear']['dn'][:,0]-c['star_aero_CRE_clear']['up'][:,0])-(c['star_noaero_CRE_clear']['dn'][:,0]-c['star_noaero_CRE_clear']['up'][:,0])
pu.plot_vert_hist(fig,ax1,DAREs,0,[-100,20],legend=True,onlyhist=False,loc=2,color='b',bins=30)
pu.plot_vert_hist(fig,ax1,DAREs_clear,1,[-100,20],legend=True,color='y',bins=30)
plt.savefig(fp+'../plot/ORACLES_DARE_surface_4STAR.png',transparent=True,dpi=600)


# In[70]:

print 'Surface DARE'
print 'mean clouds: {}, no clouds: {}'.format(np.nanmean(DAREs),np.nanmean(DAREs_clear))
print 'median clouds: {}, no clouds: {}'.format(np.nanmedian(DAREs),np.nanmedian(DAREs_clear))
print 'std clouds: {}, no clouds: {}'.format(np.nanstd(DAREs),np.nanstd(DAREs_clear))


# In[71]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-50,100],xlim=[-0.5,1.5])
ax1.set_ylabel('Direct Aerosol Radiative Effect [W/m$^2$]')
ax1.set_title('Direct Aerosol Radiative Effect from 4STAR retrievals\nTop of Atmosphere')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Clouds','Without Clouds'])
DAREt = (c['star_aero_CRE']['dn'][:,2]-c['star_aero_CRE']['up'][:,2])-       (c['star_noaero_CRE']['dn'][:,2]-c['star_noaero_CRE']['up'][:,2])
DAREt_clear = (c['star_aero_CRE_clear']['dn'][:,2]-c['star_aero_CRE_clear']['up'][:,2])-             (c['star_noaero_CRE_clear']['dn'][:,2]-c['star_noaero_CRE_clear']['up'][:,2])
pu.plot_vert_hist(fig,ax1,DAREt,0,[-50,100],legend=True,onlyhist=False,loc=2,color='b',bins=30)
pu.plot_vert_hist(fig,ax1,DAREt_clear,1,[-50,100],legend=True,color='y',bins=30)
plt.savefig(fp+'../plot/ORACLES_DARE_toa_4STAR.png',transparent=True,dpi=600)


# In[72]:

print 'TOA DARE'
print 'mean clouds: {}, no clouds: {}'.format(np.nanmean(DAREt),np.nanmean(DAREt_clear))
print 'median clouds: {}, no clouds: {}'.format(np.nanmedian(DAREt),np.nanmedian(DAREt_clear))
print 'std clouds: {}, no clouds: {}'.format(np.nanstd(DAREt),np.nanstd(DAREt_clear))


# ### Calculate the relative forcing efficiency

# In[73]:

tau_500 = 0.3689


# In[74]:

rfes = DAREs/c['star_aero_CRE']['dn'][:,2]*100.0/tau_500
rfes_clear = DAREs_clear/c['star_aero_CRE']['dn'][:,2]*100.0/tau_500
rfet = DAREt/c['star_aero_CRE']['dn'][:,2]*100.0/tau_500
rfet_clear = DAREt_clear/c['star_aero_CRE']['dn'][:,2]*100.0/tau_500 


# In[75]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-30,5],xlim=[-0.5,1.5])
ax1.set_ylabel('relative DARE efficiency [\%/$\\tau_{{500nm}}$]')
ax1.set_title('relative DARE efficiency from 4STAR retrievals\nSurface')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Clouds','Without Clouds'])
pu.plot_vert_hist(fig,ax1,rfes,0,[-30,5],legend=True,onlyhist=False,loc=4,color='b',bins=50)
pu.plot_vert_hist(fig,ax1,rfes_clear,1,[-30,5],legend=True,color='y',bins=50)
plt.savefig(fp+'../plot/ORACLES_rDAREe_surface_4STAR.png',transparent=True,dpi=600)


# In[76]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-25,25],xlim=[-0.5,1.5])
ax1.set_ylabel('relative DARE efficiency [\%/$\\tau_{{500nm}}$]')
ax1.set_title('relative DARE efficiency from 4STAR retrievals\nTop of Atmosphere')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['With Clouds','Without Clouds'])
pu.plot_vert_hist(fig,ax1,rfet,0,[-25,25],legend=True,onlyhist=False,loc=4,color='b',bins=50)
pu.plot_vert_hist(fig,ax1,rfet_clear,1,[-25,25],legend=True,color='y',bins=50)
plt.savefig(fp+'../plot/ORACLES_rDAREe_TOA_4STAR.png',transparent=True,dpi=600)


# In[77]:

print 'Surface rDAREe'
print 'mean clouds: {}, no clouds: {}'.format(np.nanmean(rfes),np.nanmean(rfes_clear))
print 'median clouds: {}, no clouds: {}'.format(np.nanmedian(rfes),np.nanmedian(rfes_clear))
print 'std clouds: {}, no clouds: {}'.format(np.nanstd(rfes),np.nanstd(rfes_clear))


# In[78]:

print 'TOA rDAREe'
print 'mean clouds: {}, no clouds: {}'.format(np.nanmean(rfet),np.nanmean(rfet_clear))
print 'median clouds: {}, no clouds: {}'.format(np.nanmedian(rfet),np.nanmedian(rfet_clear))
print 'std clouds: {}, no clouds: {}'.format(np.nanstd(rfet),np.nanstd(rfet_clear))


# ## Calculate the impact of aerosol on CRE

# In[138]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-5,15],xlim=[-0.5,1.5])
ax1.set_ylabel('relative Cloud Radiative Effect [\%]')
ax1.set_title('difference in relative CRE due to aerosols')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['Surface','TOA'])
pu.plot_vert_hist(fig,ax1,star_aero_rC[:,0]-star_noaero_rC[:,0],0,[-5,15],legend=True,onlyhist=True,loc=4,color='b',bins=50)
pu.plot_vert_hist(fig,ax1,star_aero_rC[:,2]-star_noaero_rC[:,2],1,[-5,15],legend=True,color='k',bins=50)
ax1.grid()
plt.savefig(fp+'../plot/ORACLES_rCRE_from_aerosol_4STAR.png',transparent=True,dpi=600)


# In[139]:

print 'difference in relative CRE due to aerosol'
print 'mean, surface: {}, toa: {}'.format(np.nanmean(star_aero_rC[:,0]-star_noaero_rC[:,0]),np.nanmean(star_aero_rC[:,2]-star_noaero_rC[:,2]))


# In[127]:

fig = plt.figure(figsize=(7,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-5,15],xlim=[-0.5,2.5])
ax1.set_ylabel('relative Cloud Radiative Effect [\%]')
ax1.set_title('relative Cloud Radiative Effect due to aerosols')
ax1.set_xticks([0,1,2])
ax1.set_xticklabels(['Surface\nrelative to above cloud','Surface','TOA'])
pu.plot_vert_hist(fig,ax1,star_aero_rC_abc[:,0]-star_noaero_rC_abc[:,0],0,[-5,15],legend=True,onlyhist=False,loc=4,color='b',bins=50)
pu.plot_vert_hist(fig,ax1,star_aero_rC[:,0]-star_noaero_rC[:,0],1,[-5,15],legend=True,onlyhist=True,loc=4,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star_aero_rC[:,2]-star_noaero_rC[:,2],2,[-5,15],legend=True,color='r',bins=50)
ax1.grid()
plt.savefig(fp+'../plot/ORACLES_rCRE_from_aerosol_abc_4STAR.png',transparent=True,dpi=600)


# In[98]:

c.keys()


# In[101]:

c['star_aero_C'].shape


# In[107]:

ar.keys()


# In[108]:

ar['tau_fl'].shape


# In[112]:

plt.figure()
plt.hist(np.exp(-1.0*0.36*np.cos(ar['sza'][ar['fl'].astype(bool)]*np.pi/180.0))*100.0,edgecolor='None')
plt.xlabel('Expected reduction due to aerosol optical depth')
#plt.title('TOA albedo')
#plt.savefig(fp+'../plot/ORACLES_albedo_toa_4STAR.png',transparent=True,dpi=600)


# In[115]:

plt.figure()
plt.hist( np.cos(ar['sza'][ar['fl'].astype(bool)]*np.pi/180.0) )


# In[117]:

np.cos(4*np.pi/180)


# In[118]:

1.0/np.cos(45.0*np.pi/180.0)


# Ratio the expected reduction due to aerosol

# In[119]:

aero_tau_ratio = np.exp(-1.0*0.36*np.cos(ar['sza'][ar['fl'].astype(bool)]*np.pi/180.0))


# In[123]:

plt.figure()
plt.hist(aero_tau_ratio)


# In[120]:

print 'Surface CRE'
print 'mean aero: {}, no aero: {}'.format(np.nanmean(c['star_aero_C'][:,0]),np.nanmean(c['star_noaero_C'][:,0]))
print 'median aero: {}, no aero: {}'.format(np.nanmedian(c['star_aero_C'][:,0]),np.nanmedian(c['star_noaero_C'][:,0]))
print 'std aero: {}, no aero: {}'.format(np.nanstd(c['star_aero_C'][:,0]),np.nanstd(c['star_noaero_C'][:,0]))


# In[122]:

print 'Surface CRE corrected for ratio'
print 'mean aero: {}, no aero: {}'.format(np.nanmean(c['star_aero_C'][:,0]),np.nanmean(c['star_noaero_C'][:,0]*aero_tau_ratio))
print 'median aero: {}, no aero: {}'.format(np.nanmedian(c['star_aero_C'][:,0]),np.nanmedian(c['star_noaero_C'][:,0]*aero_tau_ratio))
print 'std aero: {}, no aero: {}'.format(np.nanstd(c['star_aero_C'][:,0]),np.nanstd(c['star_noaero_C'][:,0]*aero_tau_ratio))


# In[125]:

-571.2/-610.3


# In[134]:

np.exp(-0.1)


# In[ ]:



