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

def build_aac_input(fp,fp_alb,fp_out):
    """
    Program to build the inputs of libradtran for Meloe's AAC study
    """
    import numpy as np
    import scipy.io as sio
    import Run_libradtran as RL
    import load_modis as lm
    pmom = RL.make_pmom_inputs()
    geo = {'zout':[0,3,100],'year':2007,'day':15,'minute':0,'second':0}
    aero = {'z_arr':[3.0,4.0]}
    cloud = {'ztop':3.0,'zbot':2.0,'phase':'wc','write_moments_file':True,'moms_dict':pmom}
    source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/'}
    albedo = {'create_albedo_file':False}
    for mmm in ['DJF','MAM','JJA','SON']:
        fpm = fp+'Input_to_DARF_%s.mat' % mmm
        print 'in %s months, getting mat file: %s' % (mmm,fpm)
        input_mmm = sio.loadmat(fpm,mat_dtype=True)['data_input_darf']
        if mmm=='DJF':
            geo['month'] = 1
            doy = 17
        elif mmm=='MAM':
            geo['month'] = 4
            doy = 137
        elif mmm=='JJA':
            geo['month'] = 7
            doy = 225
        elif mmm=='SON':
            geo['month'] = 10
            doy = 321

        fpa = fp_alb+'MCD43GF_geo_shortwave_%03i_2007.hdf' % doy
        print 'Getting albedo files: '+fpa
        alb_geo,alb_geo_dict = lm.load_hdf_sd(fpa)
        alb_geo_sub = np.nanmean(np.nanmean(alb_geo['MCD43GF_CMG'].reshape([48,21600/48,75,43200/75]),3),1)
        alb_geo_lat = np.linspace(90,-90,num=48)
        alb_geo_lon = np.linspace(-180,180,num=75)

        print 'Running through the files'
        for ilat,lat in enumerate(input_mmm['MODIS_lat'][0,0]):
            for ilon,lon in enumerate(input_mmm['MODIS_lon'][0,0]):
                geo['lat'],geo['lon'] = lat,lon
                # set the aerosol values
                aero['wvl_arr'] = input_mmm['MOC_wavelengths'][0,0][0,:]*1000.0
                aero['ext'] = input_mmm['MOC_ext_mean'][0,0][ilat,ilon,:]
                if np.isnan(aero['ext']).all():
                    print 'skipping lat:%i, lon:%i' % (ilat,ilon)
                    continue
                aero['ssa'] = input_mmm['MOC_ssa_mean'][0,0][ilat,ilon,:]
                aero['asy'] = input_mmm['MOC_asym_mean'][0,0][ilat,ilon,:]
                # set the cloud values
                cloud['tau'] = input_mmm['MODIS_COD_mean'][0,0][ilat,ilon]
                cloud['ref'] = input_mmm['MODIS_effrad_mean'][0,0][ilat,ilon]
                # set the albedo
                alb = alb_geo_sub[np.argmin(abs(alb_geo_lat-lat)),np.argmin(abs(alb_geo_lon-lon))]
                if np.isnan(alb): 
                    albedo['sea_surface_albedo'] = True
                else:
                    albedo['albedo'] = alb

                for HH in xrange(24):
                    geo['hour'] = HH
                    #build the solar input file
                    source['source'] = 'solar'
                    source['wvl_range'] = [250,5600]
                    cloud['write_moments_file'] = True
                    file_out_sol = fp_out+'AAC_input_lat%02i_lon%02i_%s_HH%02i_sol.inp' % (ilat,ilon,mmm,HH)
                    RL.write_input_aac(file_out_sol,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False)
                    #build the thermal input file
                    source['source'] = 'thermal'
                    source['wvl_range'] = [4000,50000]
                    cloud['write_moments_file'] = False
                    file_out_thm = fp_out+'AAC_input_lat%02i_lon%02i_%s_HH%02i_thm.inp' % (ilat,ilon,mmm,HH)
                    RL.write_input_aac(file_out_thm,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False)
                    print ilat,ilon,HH
        del alb_geo
        del input_mmm

# <codecell>

build_aac_input(fp,fp_alb,fp_out)

# <codecell>


