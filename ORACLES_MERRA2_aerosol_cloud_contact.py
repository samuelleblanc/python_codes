
# coding: utf-8

# # Intro
# Name:  
# 
#     ORACLES_MERRA2_aerosol_cloud_contact
# 
# Purpose:  
# 
#     Prepare some first try analysis of the MERRA2 reanalysis with comparisons to ORACLES campaign for studying when/if aerosol   
#     contact clouds
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
#     - load_utils.py : for loading modis and other files
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - datetime
#     - mpl_toolkits
#     - gdal (from osgeo)
#     - plotting_utils (user defined plotting routines)
#     - map_utils, dependent on geopy
#   
# Needed Files:
# 
#   - MERRA2: AOD files
#   - MERRA2: Cloud files
#   - 4STAR AOD from ORACLES 2016
#   
# Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, CA, 2018-08-08

# # Load the python modules and set paths

# In[1]:


import numpy as np
import scipy.io as sio
import os
import matplotlib.pyplot as plt


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


import load_utils as lu
from Sp_parameters import smooth

import hdf5storage as hs
from mpltools import color
from path_utils import getpath
from write_utils import nearest_neighbor
from tqdm import tqdm_notebook as tqdm


# In[4]:


import scipy.stats as st
import Sun_utils as su
from mpl_toolkits.basemap import Basemap


# In[5]:


fp = getpath('ORACLES')
fp


# # Load the various files

# In[6]:


days = ['']


# In[7]:


day = '20160914'


# ## Load 4STAR gap distance calculations for ORACLES 2016

# In[9]:


gap = hs.loadmat(fp+'ORACLES_2016_gap_distance.mat')


# In[10]:


gap['ldelta_alt'].shape


# In[287]:


gap['days'] = [i[0] for i in gap['ldelta_lat_days']]


# In[303]:


{i:gap['days'].count(i) for i in gap['days']}


# In[304]:


days = ['{:8.0f}'.format(i) for i in np.unique(gap['days'])]
days


# ## Load the MERRA2 2d aerosol files

# In[6]:


fpma = fp+'data_other/MERRA2/'


# In[248]:


merra2_gas2d = lu.load_netcdf(fpma+'MERRA2_400.inst3_2d_gas_Nx.{}.nc4.nc'.format(day),values=(('AODANA',0),('AODINC',1),('lat',2),('lon',3),('time',4)))


# In[8]:


merra2_gas2d


# ## Load the MERRA2 clouds

# In[348]:


mclds = []
for day in days:
    print day
    mcld,mcld_info = lu.load_netcdf(fpma+'clouds_model/MERRA2_400.tavg3_3d_cld_Nv.{}.SUB.nc'.format(day),everything=True)
    mclds.append(mcld)


# In[10]:


mcld_info


# In[11]:


mcld['lev'].shape


# In[217]:


mcld['DELP'].shape


# In[13]:


mcld['DELP'][1,:,50,50]


# ### Load the MERRA2 cloud pressure levels

# In[250]:


mcldp,mcldp_info = lu.load_netcdf(fpma+'clouds_press/MERRA2_400.tavg3_3d_cld_Np.{}.SUB.nc'.format(day),everything=True)


# In[15]:


mcldp['lev'].shape


# ## Load the MERRA2 aerosol

# In[349]:


maeros = []
for day in days:
    print day
    maero,maero_info = lu.load_netcdf(fpma+'aero_model/MERRA2_400.inst3_3d_aer_Nv.{}.nc4'.format(day),everything=True)
    maeros.append(maero)


# In[322]:


maero_info


# In[350]:


ll=[]
for k in maero_info.keys():
    print k,maero_info[k].long_name, maero_info[k].units
    if maero_info[k].units=='kg kg-1': ll.append(k)


# In[351]:


ll


# In[312]:


maero['tot_aero'] = 0
for l in ll:
    maero['tot_aero'] = maero['tot_aero']+maero[l]


# In[352]:


for i,day in enumerate(days):
    maeros[i]['tot_aero'] = 0
    for l in ll:
        maeros[i]['tot_aero'] = maeros[i]['tot_aero']+maeros[i][l]
    maeros[i]['scat'] = maeros[i]['tot_aero']*maeros[i]['AIRDENS']*1000.0*2.18*1000.0*1000.0
    


# In[261]:


maero['tot_aero'].shape


# In[262]:


maero['scat'] = maero['tot_aero']*maero['AIRDENS']*1000.0*2.18*1000.0*1000.0


# In[257]:


maero['AIRDENS'].shape


# ## Load the MERRA2 atmosphere files for altitudes

# In[353]:


matms = []
for day in days:
    print day
    matm,matm_info = lu.load_netcdf(fpma+'atm_model/MERRA2_400.inst3_3d_asm_Nv.{}.SUB.nc'.format(day),everything=True)
    matms.append(matm)


# In[264]:


matm_info['H']


# In[265]:


matm['H'].shape


# In[26]:


matm['H'][5,:,50,50]/1000.0


# # Make a few plots

# In[27]:


maero['lon'].shape


# In[28]:


maero['lat'].shape


# In[29]:


maero['time'].shape


# In[30]:


maero['lev'].shape


# In[31]:


maero['DMS'].shape


# In[32]:


maero['time']


# In[33]:


maero['lev']


# In[227]:


plt.figure()
plt.pcolor(maero['lon'],maero['lat'],maero['tot_aero'][2,10,:,:])


# In[37]:


maero['time']/60.0


# In[43]:


maero['lon'][60],maero['lat'][70]


# In[49]:


mcld_info['QL']


# In[63]:


matm_info['H']


# In[84]:


mcld['cld_n'] = mcld['QL']*maero['AIRDENS']


# In[228]:


mcld['time']/60.0


# In[229]:


maero['time']/60.0


# In[230]:


matm['time']/60.0


# In[231]:


ii,ij = 60,70


# In[232]:


iis = [40,50,60,70]
ijs = [60,70,80,90]


# In[90]:


for ii in iis:
    for ij in ijs:
        plt.figure()
        plt.plot(maero['scat'][1,:,ii,ij],matm['H'][3,:,ii,ij],'k-+',label='Aero')
        plt.legend(frameon=False)

        #plt.xlabel('Total Aerosol Content Mass fraction [kg kg$^{{-1}}$]')
        plt.xlabel('Total Aerosol scattering [1/Mm]')
        ax2 = plt.gca().twiny()
        ax2.plot(mcld['QL'][3,:,ii,ij],matm['H'][3,:,ii,ij],'r-s',label='Cloud')
        ax2.set_ylabel('Altitude [m]')
        ax2.set_xlabel('Mass fraction of Liquid water cloud [kg kg$^{{-1}}$]')
        plt.ylim(0,7000)
        plt.legend(frameon=False,loc=4)
        plt.title('Lat: {}, Lon: {}\n\n'.format(maero['lon'][ii],maero['lat'][ij]))


# # Now combine cloud and aerosol to find vertical distance between them

# ## Test with one case

# In[266]:


maero['isaero'] = maero['scat']>50.0
mcld['iscld'] = mcld['QL']>0.00001


# In[354]:


for i,day in enumerate(days):
    maeros[i]['isaero'] = maeros[i]['scat']>50.0
    mclds[i]['iscld'] = mclds[i]['QL']>0.00001


# In[333]:


day


# In[267]:


maero['isaero'].shape, mcld['iscld'].shape


# In[146]:


plt.figure()
plt.plot(maero['isaero'][1,:,60,72],matm['H'][3,:,60,72],'+')
plt.plot(mcld['iscld'][3,:,60,72],matm['H'][3,:,60,72],'x')
plt.xlim(-0.5,1.5)
plt.ylim(0,7000)


# In[98]:


def start_stop(a, trigger_val, len_thresh=2):
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False,np.equal(a, trigger_val),False]

    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])

    # Get lengths
    lens = idx[1::2] - idx[::2]

    return idx.reshape(-1,2)[lens>len_thresh]-[0,1]


# In[147]:


idx = start_stop(maero['isaero'][1,:,60,72],True,len_thresh=0)
idc = start_stop(mcld['iscld'][3,:,60,72],True,len_thresh=0)
idx,idc


# In[148]:


iwithin = start_stop((maero['isaero'][1,:,60,72] & mcld['iscld'][3,:,60,72]),True,len_thresh=0)
iwithin


# In[149]:


iabove,ibelow = [],[]
for ix in idx:
    for ic in idc:
        if ix[0]>=ic[1]: ibelow.append([ic[1],ix[0]])
        if ix[1]<=ic[0]: iabove.append([ix[1],ic[0]])
iabove,ibelow


# In[153]:


zwithin,zabove,zbelow = [],[],[]
for iw in iwithin: zwithin.append((matm['H'][3,iw[0],60,72]+matm['H'][3,iw[0]-1,60,72])/2.0                                  -(matm['H'][3,iw[1],60,72]+matm['H'][3,iw[1]+1,60,72])/2.0)
if not zwithin: zwithin.append(0.0)
for ia in iabove: zabove.append(matm['H'][3,ia[0],60,72]-matm['H'][3,ia[1],60,72])
if not zabove: zabove.append(0.0)
for ib in ibelow: zbelow.append(matm['H'][3,ib[0],60,72]-matm['H'][3,ib[1],60,72])
if not zbelow: zbelow.append(0.0)
zwithin,zabove,zbelow


# ## Now redo but for entire region

# In[360]:


zwi,zab,zbe = [],[],[]
for i,maeroi in enumerate(maeros):
    matmi = matms[i]
    mcldi = mclds[i]
    print days[i]
    pbar = tqdm(total=len(maeroi['time'])*len(maeroi['lat'])*len(maeroi['lon']))
    zw = np.zeros((len(maeroi['time']),len(maeroi['lat']),len(maeroi['lon']),20))*np.nan
    za =  np.zeros((len(maeroi['time']),len(maeroi['lat']),len(maeroi['lon']),20))*np.nan
    zb =  np.zeros((len(maeroi['time']),len(maeroi['lat']),len(maeroi['lon']),20))*np.nan
    for it in range(len(maeroi['time'])):
        ita = np.argmin(abs(maeroi['time'][it]-matmi['time']))
        itc = np.argmin(abs(maeroi['time'][it]-mcldi['time']))
        for ila, lat in enumerate(maeroi['lat']):
            for ilo, lon in enumerate(maeroi['lon']):
                idx = start_stop(maeroi['isaero'][it,:,ila,ilo],True,len_thresh=0)
                idc = start_stop(mcldi['iscld'][itc,:,ila,ilo],True,len_thresh=0)
                iwithin = start_stop((maeroi['isaero'][it,:,ila,ilo] & mcldi['iscld'][itc,:,ila,ilo]),True,len_thresh=0)
                iabove,ibelow = [],[]
                for ix in idx:
                    for ic in idc:
                        if ix[0]>=ic[1]: ibelow.append([ic[1],ix[0]])
                        if ix[1]<=ic[0]: iabove.append([ix[1],ic[0]])
                zwithin,zabove,zbelow = [],[],[]
                for iw in iwithin: 
                    try:
                        zwithin.append((matmi['H'][ita,iw[0],ila,ilo]+matmi['H'][ita,iw[0]-1,ila,ilo])/2.0-                                       (matmi['H'][ita,iw[1],ila,ilo]+matmi['H'][ita,iw[1]+1,ila,ilo])/2.0)
                    except IndexError:
                        zwithin.append(np.nan)
                if not zwithin: zwithin.append(np.nan)
                for ia in iabove: zabove.append(matmi['H'][ita,ia[0],ila,ilo]-matmi['H'][ita,ia[1],ila,ilo])
                if not zabove: zabove.append(np.nan)
                for ib in ibelow: zbelow.append(matmi['H'][ita,ib[0],ila,ilo]-matmi['H'][ita,ib[1],ila,ilo])
                if not zbelow: zbelow.append(np.nan)
                zw[it,ila,ilo,:len(zwithin)] = zwithin
                za[it,ila,ilo,:len(zabove)] = zabove
                zb[it,ila,ilo,:len(zbelow)] = zbelow
                pbar.update(1)
    zwi.append(np.ma.masked_array(zw,np.isnan(zw)))
    zab.append(np.ma.masked_array(za,np.isnan(za)))
    zbe.append(np.ma.masked_array(zb,np.isnan(zb)))


# In[365]:


size(zab)


# In[361]:


maeroi['time'],mcldi['time'],matm['time']


# In[269]:


zw = np.ma.masked_array(zw,np.isnan(zw))
za = np.ma.masked_array(za,np.isnan(za))
zb = np.ma.masked_array(zb,np.isnan(zb))


# ## plot the resulting MERRA overlap

# In[277]:


plt.figure()
plt.pcolor(maero['lon'],maero['lat'],za[1,:,:,0])
plt.colorbar(label='Vertical distance [m]')
plt.title('MERRA-2 {} Gap between cloud and aerosol layer'.format(day))


# # Combine MERRA and 4STAR ORACLES gap

# In[208]:


gap.keys()


# In[278]:


float(day)


# In[280]:


iday = np.array(gap['days']) == float(day)


# In[286]:


plt.figure()
plt.pcolor(maero['lon'],maero['lat'],za[1,:,:,0])
plt.colorbar(label='Vertical distance [m]')
plt.title('MERRA-2 {} Gap between cloud and aerosol layer'.format(day))

plt.scatter(gap['ldelta_lon'][iday],gap['ldelta_lat'][iday],50,c=gap['ldelta_alt'][iday],
            marker='o')


# # Combine to make an average for all days

# In[363]:


za_avg = np.zeros((len(maeroi['lat']),len(maeroi['lon'])))*np.nan
za_all = np.zeros((len(days),len(maeroi['lat']),len(maeroi['lon'])))*np.nan


# In[367]:


for i,d in enumerate(days):
    za_all[i,:,:] = zab[i][4,:,:,0]


# In[372]:


za_avg = np.nanmean(za_all,axis=0)


# In[374]:


za_avg = np.ma.masked_array(za_avg,np.isnan(za_avg))


# In[373]:


za_avg.shape


# In[417]:


fp


# In[420]:


plt.figure()
plt.contourf(maeroi['lon'],maeroi['lat'],za_avg[:,:],60,vmin=0.0,vmax=3000.0,extend='max')
plt.colorbar(label='Vertical distance [m]')
plt.title('MERRA-2 Gap between cloud and aerosol layer \nORACLES 2016 average')
plt.scatter(gap['ldelta_lon'],gap['ldelta_lat'],100,c=gap['ldelta_alt'],marker='o',edgecolor='grey',vmin=0.0,vmax=3240.0,lw=1)
plt.xlim(-30,20)
plt.ylim(-30,5)
plt.savefig(fp+'MERRA2_ORACLES2016_gap_distance.png',transparent=True,dpi=600)


# ## Add a coastline

# In[400]:


def mapfig(ax=plt.gca()):
    m = Basemap(projection='merc',llcrnrlat=-26,urcrnrlat=-8,llcrnrlon=-18,urcrnrlon=16,resolution='l',ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.linspace(-28,16,9),labels=[0,0,0,1],linewidth=0.1)
    m.drawparallels(np.linspace(-26,-8,10),labels=[1,0,0,0],linewidth=0.1)
    m.shadedrelief(alpha=0.4)
    return m


# In[405]:


lon,lat = np.meshgrid(maeroi['lon'],maeroi['lat'])


# In[416]:


fig,axa = plt.subplots(1,1)
m = mapfig(axa)
lon,lat = np.meshgrid(maeroi['lon'],maeroi['lat'])
mx,my = m(lon,lat)
axa.contourf(mx,my,za_avg[:,:],60,vmin=0.0,vmax=3000.0,extend='max',cmap='jet')
axa.set_title('MERRA-2 Gap between cloud and aerosol layer \nORACLES 2016 average')
plt.colorbar()
mxa,mya = m(gap['ldelta_lon'],gap['ldelta_lat'])
ax1.scatter(mxa,mya,100,c=gap['ldelta_alt'],marker='o',edgecolor='grey',vmin=0.0,vmax=3240.0,lw=1)

