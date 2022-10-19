#!/usr/bin/env python
# coding: utf-8

# Name:  
# 
#     Prepare_input_aac
# 
# Purpose:  
# 
#     Created input files for libradtran based on the calipso, global Aerosol over Clouds
#     For Meloë's study
# 
# Calling Sequence:
# 
#     python Prepare_input_aac
#   
# Input:
# 
#     none
# 
# Output:
#    
#     input files for libradtran 2.0 (uvspec) 
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - mplt_toolkits for basemap, map plotting
#     - pdb
#     - datetime
#     - load_utils
# 
#   
# Needed Files:
# 
#   - matlab input files: Input_to_DARF_mmm.mat

# In[ ]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap


# In[2]:


import Run_libradtran as RL
reload(RL)


# In[3]:


fp = 'C:\Users\sleblan2\Research\Calipso\meloe/'


# ## Read Matlab input files

# In[5]:


input_DJF = sio.loadmat(fp+'Input_to_DARF_DJF.mat',mat_dtype=True)['data_input_darf']


# In[6]:


input_DJF.dtype.names


# In[7]:


enumerate(input_DJF['MODIS_lat'][0,0])


# In[8]:


input_DJF['MODIS_lat'][0,0].shape


# In[9]:


input_DJF['MODIS_COD_mean'][0,0][input_DJF['MODIS_COD_mean'][0,0]==-9999] = np.nan


# In[21]:


ilat,ilon = 5,15


# In[23]:


wvl = np.append(wvl,100000.0)


# In[24]:


wvl


# In[27]:


ext = np.abs(input_DJF['MOC_ext_mean'][0,0][ilat,ilon,:])


# In[22]:


input_DJF['MOC_ext_mean'][0,0][ilat,ilon,:]


# In[11]:


wvl = input_DJF['MOC_wavelengths'][0,0][0,:]*1000.0


# In[91]:


ctr = plt.contourf(input_DJF['MODIS_lon'][0,0][:,0],input_DJF['MODIS_lat'][0,0][:,0],input_DJF['MODIS_COD_mean'][0,0])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(ctr)
cbar.set_label('MODIS MEAN COD')
plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=11,c='k')
plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=10,c=input_DJF['MODIS_COD_mean'][0,0])
plt.xlim([-180,180])
plt.ylim([-90,90])
plt.savefig(fp+'plots/MODIS_mean_COD_flat.png',dpi=600,transparent=True)


# In[105]:


MODIS_lon[20,72]


# In[104]:


input_DJF['MODIS_COD_mean'][0,0][:,72]


# In[93]:


m = Basemap(projection='moll',lon_0=0,resolution='c')
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,20.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,420.,60.))
MODIS_lon, MODIS_lat = np.meshgrid(input_DJF['MODIS_lon'][0,0][:,0],input_DJF['MODIS_lat'][0,0][:,0],sparse=False)
x,y=m(MODIS_lon, MODIS_lat) 
cs = m.contourf(x,y,input_DJF['MODIS_COD_mean'][0,0],np.linspace(0,60,21),extend='max')
cbar = plt.colorbar(cs)
cbar.set_label('MODIS mean COD')
plt.savefig(fp+'plots/MODIS_mean_COD.png',dpi=600,transparent=True)


# In[106]:


input_DJF['MODIS_COD_mean'][0,0][20,20]


# In[107]:


input_DJF['MODIS_effrad_mean'][0,0][20,20]


# In[111]:


input_DJF['MOC_wavelengths'][0,0]


# In[113]:


print MODIS_lon[20,20],MODIS_lat[20,20]


# In[119]:


input_DJF['MOC_asym_mean'][0,0][20,20]


# In[26]:


input_mmm['MOC_wavelengths'][0,0][0,:]*1000.0


# In[30]:


input_mmm['MOC_ext_mean'][0,0].shape


# In[34]:


input_mmm['MODIS_effrad_mean'][0,0].shape


# In[29]:


input_mmm['MODIS_lat'][0,0].shape


# ## Read sample MODIS albedo file

# In[5]:


import load_utils as lm
reload(lm)


# In[6]:


alb_geo,alb_geo_dict = lm.load_hdf_sd('C:\Users\sleblan2\Research\Calipso\meloe\MCD43GF_geo_shortwave_001_2007.hdf')


# In[7]:


alb_geo['MCD43GF_CMG'].shape


# In[8]:


alb_geo['MCD43GF_CMG'][1,0]


# In[76]:


alb_geo_dict['MCD43GF_CMG']['_FillValue']


# In[73]:


input_DJF['MODIS_COD_mean'][0,0].shape


# In[12]:


alb_geo_sub = np.nanmean(np.nanmean(alb_geo['MCD43GF_CMG'].reshape([48,21600/48,75,43200/75]),3),1)


# In[19]:


alb_geo_lat = np.linspace(90,-90,num=48)
alb_geo_lon = np.linspace(-180,180,num=75)


# In[21]:


co = plt.contourf(alb_geo_lon,alb_geo_lat,alb_geo_sub)
cbar = plt.colorbar(co)


# Check netcdf

# In[9]:


fp_rtm='C:/Users/sleblan2/Research/4STAR/rtm_dat/'


# In[10]:


mie = sio.netcdf_file(fp_rtm+'wc_allpmom.sol.mie.cdf','r')


# In[19]:


mie_short = {'wvl':mie.variables['wavelen'].data,
                'ref':mie.variables['reff'].data,
                'ntheta':np.swapaxes(mie.variables['ntheta'].data[:,:,0],0,1),
                'rho':np.swapaxes(mie.variables['rho'].data,0,1),
                'nmom':np.swapaxes(mie.variables['nmom'].data,0,1),
                'ssa':np.swapaxes(mie.variables['ssa'].data,0,1),
                'ext':np.swapaxes(mie.variables['ext'].data,0,1),
                'nim':mie.variables['refim'].data,
                'nre':mie.variables['refre'].data,
                'pmom':np.swapaxes(mie.variables['pmom'].data[:,:,0,:],0,1),
                'phase':np.swapaxes(mie.variables['phase'].data[:,:,0,:],0,1),
                'theta': np.swapaxes(mie.variables['theta'].data[:,:,0,:],0,1)}


# In[57]:


mie.variables['theta'].data.shape


# In[58]:


mie_long.variables['theta'].data.shape


# In[59]:


mie_short['theta'].shape


# In[60]:


pmom = {'wvl':np.append(mie_short['wvl'],mie_long.variables['wavelen'].data[7:]),
        'ref':mie_short['ref'],
        'ntheta':np.concatenate((mie_short['ntheta'],np.swapaxes(mie_long.variables['ntheta'].data[7:,:-5,0],0,1)),axis=1),
        'rho':mie_short['rho'],
        'nmom':np.concatenate((mie_short['nmom'],np.swapaxes(mie_long.variables['nmom'].data[7:,:-5,0],0,1)),axis=1),
        'ssa':np.concatenate((mie_short['ssa'],np.swapaxes(mie_long.variables['ssa'].data[7:,:-5],0,1)),axis=1),
        'ext':np.concatenate((mie_short['ext'],np.swapaxes(mie_long.variables['ext'].data[7:,:-5],0,1)),axis=1),
        'nim':np.append(mie_short['nim'],mie_long.variables['refim'].data[7:]),
        'nre':np.append(mie_short['nre'],mie_long.variables['refre'].data[7:]),
        'pmom':np.concatenate((mie_short['pmom'],np.concatenate((np.swapaxes(mie_long.variables['pmom'].data[7:,:-5,0,:],0,1),
                                                                 np.zeros((25,72,2500))),axis=2)),axis=1),
        'phase':np.concatenate((mie_short['phase'],np.swapaxes(mie_long.variables['phase'].data[7:,:-5,0,:],0,1)),axis=1),
        'theta':np.concatenate((mie_short['theta'],np.swapaxes(mie_long.variables['theta'].data[7:,:-5,0,:],0,1)),axis=1)}


# ## Prepare inputs of aac for libradtran

# Set up the default, not changing properties

# In[36]:


input_mmm = input_DJF
mmm = 'DJF'


# In[40]:


pmom = RL.make_pmom_inputs()


# In[59]:


geo = {'zout':[0,3,100],'year':2007,'month':1,'day':15,'minute':0,'second':0}
aero = {'z_arr':[3.0,4.0]}
cloud = {'ztop':3.0,'zbot':2.0,'phase':'wc','write_moments_file':True,'moms_dict':pmom}
source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/'}
albedo = {'create_albedo_file':False}


# In[2]:


fp = 'C:\Users\sleblan2\Research\Calipso\meloe/'
fp_alb = fp
fp_out = 'C:\Users\sleblan2\Research\Calipso\meloe\input/'


# In[ ]:


RL.build_aac_input(fp,fp_alb,fp_out)


# ## Read the output of DARF

# In[4]:


fp


# In[5]:


mam = sio.loadmat(fp+'DARF/MAM_DARF.mat')


# In[8]:


mam.keys()


# In[9]:


mam['SW_DARF'].shape


# In[10]:


mam['lon'].shape


# In[11]:


ctr = plt.contourf(mam['lon'][:,0],mam['lat'][:,0],mam['SW_DARF'][2,:,:])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(ctr)
cbar.set_label('SW DARF [W/m$^{2}$]')
#plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=11,c='k')
#plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=10,c=input_DJF['MODIS_COD_mean'][0,0])
plt.xlim([-180,180])
plt.ylim([-90,90])
plt.savefig(fp+'plots/MAM_DARF_SW.png',dpi=600,transparent=True)


# In[12]:


ctr = plt.contourf(mam['lon'][:,0],mam['lat'][:,0],mam['LW_DARF'][2,:,:])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(ctr)
cbar.set_label('LW DARF [W/m$^{2}$]')
#plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=11,c='k')
#plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=10,c=input_DJF['MODIS_COD_mean'][0,0])
plt.xlim([-180,180])
plt.ylim([-90,90])
plt.savefig(fp+'plots/MAM_DARF_LW.png',dpi=600,transparent=True)


# In[15]:


ctr = plt.contourf(mam['lon'][:,0],mam['lat'][:,0],mam['LW_DARF'][2,:,:]+mam['SW_DARF'][2,:,:])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
cbar = plt.colorbar(ctr)
cbar.set_label('net DARF [W/m$^{2}$]')
#plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=11,c='k')
#plt.scatter(MODIS_lon, MODIS_lat,marker='x',s=10,c=input_DJF['MODIS_COD_mean'][0,0])
plt.xlim([-180,180])
plt.ylim([-90,90])
plt.savefig(fp+'plots/MAM_DARF_net.png',dpi=600,transparent=True)


# In[54]:


def make_darf_plots(mam,mmm):
    fig,ax = plt.subplots(3,1,figsize=(6,13))
    ax.flatten()
    cl0 = np.linspace(-200,150,31)
    ctr = ax[0].contourf(mam['lon'][:,0],mam['lat'][:,0],mam['SW_DARF'][2,:,:],cl0,extend='both')
    ax[0].set_xlabel('Longitude')
    ax[0].set_ylabel('Latitude')
    cbar = plt.colorbar(ctr,ax=ax[0])
    cbar.set_label('SW DARF [W/m$^{2}$]')
    ax[0].set_xlim([-180,180])
    ax[0].set_ylim([-90,90])

    cl1 = np.linspace(0,150,31)
    ctr1 = ax[1].contourf(mam['lon'][:,0],mam['lat'][:,0],mam['LW_DARF'][2,:,:],cl1,cmap=plt.cm.gist_heat,extend='max')
    ax[1].set_xlabel('Longitude')
    ax[1].set_ylabel('Latitude')
    cbar = plt.colorbar(ctr1,ax=ax[1])
    cbar.set_label('LW DARF [W/m$^{2}$]')
    ax[1].set_xlim([-180,180])
    ax[1].set_ylim([-90,90])

    cl2 = np.linspace(-150,150,41)
    ctr2 = ax[2].contourf(mam['lon'][:,0],mam['lat'][:,0],mam['SW_DARF'][2,:,:]+mam['LW_DARF'][2,:,:],cl2,cmap=plt.cm.RdBu,extend='both')
    ax[2].set_xlabel('Longitude')
    ax[2].set_ylabel('Latitude')
    cbar = plt.colorbar(ctr2,ax=ax[2])
    cbar.set_label('net DARF [W/m$^{2}$]')
    ax[2].set_xlim([-180,180])
    ax[2].set_ylim([-90,90])
    ax[0].set_title(mmm)

    plt.savefig(fp+'plots/'+mmm+'_DARF.png',dpi=600,transparent=True)


# In[59]:


mam = sio.loadmat(fp+'DARF/MAM_DARF.mat')
make_darf_plots(mam,'MAM')


# In[56]:


djf = sio.loadmat(fp+'DARF/DJF_DARF.mat')
make_darf_plots(djf,'DJF')


# In[57]:


jja = sio.loadmat(fp+'DARF/JJA_DARF.mat')
make_darf_plots(jja,'JJA')


# In[58]:


son = sio.loadmat(fp+'DARF/SON_DARF.mat')
make_darf_plots(son,'SON')


# In[62]:


for ilat,lat in enumerate(mam['lon']):
    print ilat,lat


# In[61]:


ctr = plt.contourf(mam['lon'][:,0],mam['lat'][:,0],mam['SW_DARF'][2,:,:])
plt.xlim(25,75)
plt.ylim(-15,-35)


# In[77]:


mam_sw = sio.loadmat(fp+'DARF/AAC_MAM.mat')


# In[78]:


mam_sw.keys()


# In[79]:


ctr = plt.contourf(mam_sw['lon'][:,0],mam_sw['lat'][:,0],mam_sw['SW_irr_dn_avg'][2,:,:],50)
plt.colorbar(ctr)


# In[80]:


ctr = plt.contourf(mam_sw['lon'][:,0],mam_sw['lat'][:,0],mam_sw['SW_irr_up_avg'][2,:,:],50)
plt.colorbar(ctr)


# In[81]:


ctr = plt.contourf(mam_sw['lon'][:,0],mam_sw['lat'][:,0],mam_sw['LW_irr_dn_avg'][2,:,:],50)
plt.colorbar(ctr)


# In[82]:


ctr = plt.contourf(mam_sw['lon'][:,0],mam_sw['lat'][:,0],mam_sw['LW_irr_up_avg'][2,:,:],50)
plt.colorbar(ctr)


# In[72]:


mam_sw_cl= sio.loadmat(fp+'DARF/AAC_MAM_clear.mat')


# In[76]:


ctr = plt.contourf(mam_sw_cl['lon'][:,0],mam_sw_cl['lat'][:,0],mam_sw_cl['LW_irr_dn_avg'][2,:,:],50)
plt.colorbar(ctr)


# ## Check the difference when assuming kato2 vs. fu liou calculations

# In[83]:


djf_fu = sio.loadmat(fp+'DARF/AAC_DJF.mat')
djf_kato = sio.loadmat(fp+'DARF/AAC_DJF_kato.mat')


# In[84]:


djf_fu.keys()


# In[85]:


djf_kato.keys()


# In[89]:


djf_fu['SW_irr_up_avg'].shape


# In[153]:


fig,(ax1,ax2) = plt.subplots(1,2,figsize=(13,5))
ctr = ax1.contourf(djf_fu['lon'][:,0],djf_fu['lat'][:,0],djf_fu['SW_irr_up_avg'][2,:,:])
plt.colorbar(ctr,ax=ax1)
ctr = ax2.contourf(djf_kato['lon'][:,0],djf_kato['lat'][:,0],djf_kato['SW_irr_up_avg'][2,:,:])
plt.colorbar(ctr,ax=ax2)


# In[96]:


djf_dif = djf_fu['SW_irr_up_avg'][2,:,:]-djf_kato['SW_irr_up_avg'][2,:,:]


# In[99]:


ctr = plt.contourf(djf_fu['lon'][:,0],djf_fu['lat'][:,0],djf_dif)
plt.colorbar(ctr)


# In[104]:


dddjf = djf_dif.flatten()


# In[109]:


from Sp_parameters import nanmasked


# In[110]:


dddjf,idjf = nanmasked(djf_dif.flatten())


# In[126]:


dddjf.shape


# In[150]:


plt.plot(djf_fu['SW_irr_up_avg'][2,:,:].flatten())


# In[147]:


plt.hist(dddjf,normed=True,bins=30)
plt.xlabel('SW Irradiance TOA difference [W/m$^{2}$]')
plt.title('Fu liou - Kato')
plt.savefig(fp+'plots/diff_fuliou_kato.png',dpi=600,transparent=True)


# In[151]:


djf_kato_clear = sio.loadmat(fp+'DARF/AAC_DJF_clear_kato.mat')
djf_fu_clear = sio.loadmat(fp+'DARF/AAC_DJF_clear.mat')


# In[154]:


fig,(ax1,ax2) = plt.subplots(1,2,figsize=(13,5))
ctr = ax1.contourf(djf_fu_clear['lon'][:,0],djf_fu_clear['lat'][:,0],djf_fu_clear['SW_irr_up_avg'][2,:,:])
plt.colorbar(ctr,ax=ax1)
ctr = ax2.contourf(djf_kato_clear['lon'][:,0],djf_kato_clear['lat'][:,0],djf_kato_clear['SW_irr_up_avg'][2,:,:])
plt.colorbar(ctr,ax=ax2)


# In[155]:


darf_kato = djf_kato['SW_irr_up_avg'][2,:,:] - djf_kato_clear['SW_irr_up_avg'][2,:,:]
darf_fu = djf_fu['SW_irr_up_avg'][2,:,:] - djf_fu_clear['SW_irr_up_avg'][2,:,:]


# In[172]:


fig,(ax1,ax2) = plt.subplots(1,2,figsize=(13,5))
ctr = ax1.contourf(djf_fu_clear['lon'][:,0],djf_fu_clear['lat'][:,0],darf_fu,40)
cbr = plt.colorbar(ctr,ax=ax1)
cbr.set_label('DARF [W/m$^{2}$]')
ax1.set_title('DJF DARF Fu Liou')
ctr = ax2.contourf(djf_kato_clear['lon'][:,0],djf_kato_clear['lat'][:,0],darf_kato,40)
cbr = plt.colorbar(ctr,ax=ax2)
cbr.set_label('DARF [W/m$^{2}$]')
ax2.set_title('DJF DARF Kato')
plt.savefig(fp+'plots/DARF_fu_vs_kato.png',dpi=600,transparent=True)


# In[157]:


dif_darf = darf_fu-darf_kato


# In[193]:


clevels = np.linspace(-10,10,41)


# In[195]:


clevels[20]


# In[198]:


help(plt.contourf)


# In[201]:



ctr = plt.contourf(djf_fu_clear['lon'][:,0],djf_fu_clear['lat'][:,0],dif_darf,levels=clevels,cmap=plt.cm.bwr)
cbr = plt.colorbar(ctr)
cbr.set_label('DARF difference [W/m$^{2}$]')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('DARF difference Fu Liou - Kato')
plt.savefig(fp+'plots/DARF_diff_lat_lon.png',dpi=600,transparent=True)


# In[183]:


fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
ax1.hist(nanmasked(dif_darf)[0],bins=60,normed=True)
ax1.set_title('DARF difference Fu Liou - Kato')
#ax1.set_xlabel('DARF difference [W/m$^{2}$]')
ax1.set_ylabel('Relative frequency')

ax2.hist(nanmasked(dif_darf)[0],bins=60,normed=True)
#ax2.set_title('DARF difference Fu Liou - Kato')
ax2.set_xlabel('DARF difference [W/m$^{2}$]')
ax2.set_ylabel('Relative frequency')
ax2.set_ylim([0,0.08])
ax2.text(-10,0.05,'Zoomed in')
plt.savefig(fp+'plots/DARF_difference_fu_kato.png',dpi=600,transparent=True)


# In[162]:


plt.hist(nanmasked(dif_darf/darf_kato)[0],bins=60,normed=True)


# In[167]:


rms_darf = np.sqrt(np.nanmean(darf_fu**2-darf_kato**2))


# In[168]:


rms_darf

