# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

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
# 
#   
# Needed Files:
# 
#   - matlab input files: Input_to_DARF_mmm.mat

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap

# <codecell>

import Run_libradtran as RL
reload(RL)

# <codecell>

fp = 'C:\Users\sleblan2\Research\Calipso\meloe/'

# <headingcell level=2>

# Read Matlab input files

# <codecell>

input_DJF = sio.loadmat(fp+'Input_to_DARF_DJF.mat',mat_dtype=True)['data_input_darf']

# <codecell>

input_DJF.dtype.names

# <codecell>

enumerate(input_DJF['MODIS_lat'][0,0])

# <codecell>

input_DJF['MODIS_lat'][0,0].shape

# <codecell>

input_DJF['MODIS_COD_mean'][0,0][input_DJF['MODIS_COD_mean'][0,0]==-9999] = np.nan

# <codecell>

ilat,ilon = 5,15

# <codecell>

wvl = np.append(wvl,100000.0)

# <codecell>

wvl

# <codecell>

ext = np.abs(input_DJF['MOC_ext_mean'][0,0][ilat,ilon,:])

# <codecell>

input_DJF['MOC_ext_mean'][0,0][ilat,ilon,:]

# <codecell>

wvl = input_DJF['MOC_wavelengths'][0,0][0,:]*1000.0

# <codecell>

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

# <codecell>

MODIS_lon[20,72]

# <codecell>

input_DJF['MODIS_COD_mean'][0,0][:,72]

# <codecell>

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

# <codecell>

input_DJF['MODIS_COD_mean'][0,0][20,20]

# <codecell>

input_DJF['MODIS_effrad_mean'][0,0][20,20]

# <codecell>

input_DJF['MOC_wavelengths'][0,0]

# <codecell>

print MODIS_lon[20,20],MODIS_lat[20,20]

# <codecell>

input_DJF['MOC_asym_mean'][0,0][20,20]

# <codecell>

input_mmm['MOC_wavelengths'][0,0][0,:]*1000.0

# <codecell>

input_mmm['MOC_ext_mean'][0,0].shape

# <codecell>

input_mmm['MODIS_effrad_mean'][0,0].shape

# <codecell>

input_mmm['MODIS_lat'][0,0].shape

# <headingcell level=2>

# Read sample MODIS albedo file

# <codecell>

import load_modis as lm
reload(lm)

# <codecell>

alb_geo,alb_geo_dict = lm.load_hdf_sd('C:\Users\sleblan2\Research\Calipso\meloe\MCD43GF_geo_shortwave_001_2007.hdf')

# <codecell>

alb_geo['MCD43GF_CMG'].shape

# <codecell>

alb_geo['MCD43GF_CMG'][1,0]

# <codecell>

alb_geo_dict['MCD43GF_CMG']['_FillValue']

# <codecell>

input_DJF['MODIS_COD_mean'][0,0].shape

# <codecell>

alb_geo_sub = np.nanmean(np.nanmean(alb_geo['MCD43GF_CMG'].reshape([48,21600/48,75,43200/75]),3),1)

# <codecell>

alb_geo_lat = np.linspace(90,-90,num=48)
alb_geo_lon = np.linspace(-180,180,num=75)

# <codecell>

co = plt.contourf(alb_geo_lon,alb_geo_lat,alb_geo_sub)
cbar = plt.colorbar(co)

# <markdowncell>

# Check netcdf

# <codecell>

fp_rtm='C:/Users/sleblan2/Research/4STAR/rtm_dat/'

# <codecell>

mie = sio.netcdf_file(fp_rtm+'wc_allpmom.sol.mie.cdf','r')

# <codecell>

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

# <codecell>

mie.variables['theta'].data.shape

# <codecell>

mie_long.variables['theta'].data.shape

# <codecell>

mie_short['theta'].shape

# <codecell>

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

# <headingcell level=2>

# Prepare inputs of aac for libradtran

# <markdowncell>

# Set up the default, not changing properties

# <codecell>

input_mmm = input_DJF
mmm = 'DJF'

# <codecell>

pmom = RL.make_pmom_inputs()

# <codecell>

geo = {'zout':[0,3,100],'year':2007,'month':1,'day':15,'minute':0,'second':0}
aero = {'z_arr':[3.0,4.0]}
cloud = {'ztop':3.0,'zbot':2.0,'phase':'wc','write_moments_file':True,'moms_dict':pmom}
source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/'}
albedo = {'create_albedo_file':False}

# <codecell>

fp = 'C:\Users\sleblan2\Research\Calipso\meloe/'
fp_alb = fp
fp_out = 'C:\Users\sleblan2\Research\Calipso\meloe\input/'

# <codecell>

RL.build_aac_input(fp,fp_alb,fp_out)

# <headingcell level=2>

# Read the output of DARF

# <codecell>

fp

# <codecell>

mam = sio.loadmat(fp+'DARF/MAM_DARF.mat')

# <codecell>

mam.keys()

# <codecell>

mam['SW_DARF'].shape

# <codecell>

mam['lon'].shape

# <codecell>

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

# <codecell>

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

# <codecell>

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

# <codecell>

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

# <codecell>

mam = sio.loadmat(fp+'DARF/MAM_DARF.mat')
make_darf_plots(mam,'MAM')

# <codecell>

djf = sio.loadmat(fp+'DARF/DJF_DARF.mat')
make_darf_plots(djf,'DJF')

# <codecell>

jja = sio.loadmat(fp+'DARF/JJA_DARF.mat')
make_darf_plots(jja,'JJA')

# <codecell>

son = sio.loadmat(fp+'DARF/SON_DARF.mat')
make_darf_plots(son,'SON')

# <codecell>

for ilat,lat in enumerate(mam['lon']):
    print ilat,lat

# <codecell>

ctr = plt.contourf(mam['lon'][:,0],mam['lat'][:,0],mam['SW_DARF'][2,:,:])
plt.xlim(25,75)
plt.ylim(-15,-35)

# <codecell>

mam_sw = sio.loadmat(fp+'DARF/AAC_MAM.mat')

# <codecell>

mam_sw.keys()

# <codecell>

ctr = plt.contourf(mam_sw['lon'][:,0],mam_sw['lat'][:,0],mam_sw['SW_irr_dn_avg'][2,:,:],50)
plt.colorbar(ctr)

# <codecell>

ctr = plt.contourf(mam_sw['lon'][:,0],mam_sw['lat'][:,0],mam_sw['SW_irr_up_avg'][2,:,:],50)
plt.colorbar(ctr)

# <codecell>

ctr = plt.contourf(mam_sw['lon'][:,0],mam_sw['lat'][:,0],mam_sw['LW_irr_dn_avg'][2,:,:],50)
plt.colorbar(ctr)

# <codecell>

ctr = plt.contourf(mam_sw['lon'][:,0],mam_sw['lat'][:,0],mam_sw['LW_irr_up_avg'][2,:,:],50)
plt.colorbar(ctr)

# <codecell>

mam_sw_cl= sio.loadmat(fp+'DARF/AAC_MAM_clear.mat')

# <codecell>

ctr = plt.contourf(mam_sw_cl['lon'][:,0],mam_sw_cl['lat'][:,0],mam_sw_cl['LW_irr_dn_avg'][2,:,:],50)
plt.colorbar(ctr)

# <headingcell level=2>

# Check the difference when assuming kato2 vs. fu liou calculations

# <codecell>

djf_fu = sio.loadmat(fp+'DARF/AAC_DJF.mat')
djf_kato = sio.loadmat(fp+'DARF/AAC_DJF_kato.mat')

# <codecell>

djf_fu.keys()

# <codecell>

djf_kato.keys()

# <codecell>

djf_fu['SW_irr_up_avg'].shape

# <codecell>

fig,(ax1,ax2) = plt.subplots(1,2)
ctr = ax1.contourf(djf_fu['lon'][:,0],djf_fu['lat'][:,0],djf_fu['SW_irr_up_avg'][2,:,:])
plt.colorbar(ctr,ax=ax1)
ctr = ax2.contourf(djf_kato['lon'][:,0],djf_kato['lat'][:,0],djf_kato['SW_irr_up_avg'][2,:,:])
plt.colorbar(ctr,ax=ax2)

# <codecell>

djf_dif = djf_fu['SW_irr_up_avg'][2,:,:]-djf_kato['SW_irr_up_avg'][2,:,:]

# <codecell>

ctr = plt.contourf(djf_fu['lon'][:,0],djf_fu['lat'][:,0],djf_dif)
plt.colorbar(ctr)

# <codecell>

dddjf = djf_dif.flatten()

# <codecell>

from Sp_parameters import nanmasked

# <codecell>

dddjf,idjf = nanmasked(djf_dif.flatten())

# <codecell>

dddjf.shape

# <codecell>

plt.plot(djf_fu['SW_irr_up_avg'][2,:,:].flatten())

# <codecell>

plt.hist(dddjf,normed=True,bins=30)
plt.xlabel('SW Irradiance TOA difference [W/m$^{2}$]')
plt.title('Fu liou - Kato')
plt.savefig(fp+'plots/diff_fuliou_kato.png',dpi=600,transparent=True)

