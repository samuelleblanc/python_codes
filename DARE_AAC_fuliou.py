
# coding: utf-8

# Name:  
# 
#     DARE_AAC_fuliou
# 
# Purpose:  
# 
#     Python modules to run the DARE calculations from the AAC paper (Meloë's) using the Fuliou radiative trasnfer algorithm
#     is used to port over the same functions as the Run_libradtran.build_aac_input function, but without the underlying cloud
#     properties.
#     
#     Also houses the reading of the files
#     
# 
# Calling Sequence:
# 
#     import DARE_AAC_fuiou as daf
#   
# Input:
# 
#     none at command line
#     see methods of module
# 
# Output:
#    
#     input files for fuliou
#     read of output of files
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
#     - Run_fuliou
#   
# Needed Files:
# 
#   - ...
#   
#   
# Modification History:
# 
#     Wrtten: Samuel LeBlanc, NASA Ames, 2018-05-29
#     Modified: 

# In[ ]:


if __name__ == '__main__':
    import argparse
    long_description = """    Prepare or read the fuliou DARE AAC calculations"""
    parser = argparse.ArgumentParser(description=long_description)
    parser.add_argument('-doread','--doread',help='if set, will only read the output, not produce them',
                        action='store_true')
    parser.add_argument('-dowrite','--dowrite',help='if set, will write the input and list files for fuliou',
                        action='store_true')
    parser.add_argument('-i','--index',help='Sets the index of the pixel file to use (19374,22135). Default is 0',type=int)
    in_ = vars(parser.parse_args())
    doread = in_.get('doread',False)
    dowrite = in_.get('dowrite',False)
    i = in_.get('index',0)

    if dowrite:
        run_fuliou_linux(i=i)
    elif doread:
        read_fuliou_linux(i=i)
    else:
        run_fuliou_linux(i=i)


# For writing

# In[ ]:


# Program to prepare and call the build_aac_input

import Run_libradtran as RL
import os
import argparse
long_d = """ Program to build the AAC files """
parser = argparse.ArgumentParser(description=long_d)
parser.add_argument('mmm',nargs='?',help='3 letter season designator (DJF, MAM, JJA, SON)')
parser.add_argument('sub',nargs='?',help='string designating the subset (ssam,ssap,...)')
parser.add_argument('-c','--clear',help='if set, will write the input for no aerosol case',
                    action='store_true',default=False)
parser.add_argument('-t','--tmp_folder',nargs='?',help='The path to the temp directory to use')
in_ = vars(parser.parse_args())
doclear = in_.get('clear',False)
if doclear:
    print 'doing the clear no aerosol version'
ms = in_.get('mmm',False)
if not ms:
    ms = ['DJF','MAM','JJA','SON']
else:
    ms = [ms]
print 'doing the season:',ms

if in_.get('tmp_folder'):
    tmp_folder = in_.get('tmp_folder')
    fp_out = tmp_folder+'/AAC_DARF/input/v3b_faac/'
else:
    fp_out = '/nobackup/sleblan2/AAC_DARF/input/v3b_faac/'

if not os.path.exists(fp_out):
    os.makedirs(fp_out)
fp_in = fp_out.replace('input','output')
if not os.path.exists(fp_in):
    os.makedirs(fp_in)

print 'working on folder: '+fp_out

vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'
fp = '/u/sleblan2/meloe_AAC/v3_faac/'
fp_alb = '/nobackup/sleblan2/AAC_DARF/surface_albedo/'
fp_pmom = '/nobackup/sleblan2/AAC_DARF/rtm/'
fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
wvl_thm = '/nobackup/sleblan2/AAC_DARF/rtm/wvl.thm.dat'

stdfac = {'COD':+0.0}
RL.build_aac_input(fp=fp,fp_alb=fp_alb,fp_out=fp_out,fp_pmom=fp_pmom,version=vv,stdfac_dict=stdfac,list_only=False,aero_clear=doclear,mmmlist=ms,deltascale=False)


# For reading

# In[ ]:


import Run_libradtran as RL
import scipy.io as sio

import argparse
long_d = """ Program to read in the AAC files """
parser = argparse.ArgumentParser(description=long_d)
parser.add_argument('mmm',nargs='?',help='3 letter season designator (DJF, MAM, JJA, SON)')
parser.add_argument('sub',nargs='?',help='string designating the subset (ssam,ssap,...)')
parser.add_argument('-t','--tmp_folder',nargs='?',help='The path to the temp directory to use')
in_ = vars(parser.parse_args())
ms = in_.get('mmm',False)
if not ms:
    ms = ['DJF','MAM','JJA','SON']
else:
    ms = [ms]
print 'doing the season:',ms

if in_.get('tmp_folder'):
    tmp_folder = in_.get('tmp_folder')
    fp_out = tmp_folder+'/AAC_DARF/output/v3b_faac/'
else:
    fp_out = '/nobackup/sleblan2/AAC_DARF/output/v3b_faac/'

sub = in_.get('sub',False)
if not sub:
     sub = '' #'_asymp'
     subi = sub
elif sub=='clear':
     sub = '_clear'
     subi = '_clear'
else:
     sub = '_'+sub
     subi = sub
print 'doing the subset: ', sub

import time
time.sleep(0.5)

#vv = 'v3_20160408_CALIOP_AAC_AOD_With_Sa_only'
#vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'
vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'
for mmm in ms:
#     fp_out = '/nobackup/sleblan2/AAC_DARF/output/v3_without/{}/'.format(mmm)
     fp_mat = '/u/sleblan2/meloe_AAC/v3_faac/Input_to_DARF_{m}_{v}.mat'.format(m=mmm,v=vv)
#     aac = RL.read_aac(fp_out,fp_mat,mmm=mmm)
#     sio.savemat('/nobackup/sleblan2/AAC_DARF/mat/v3_without/AAC_{m}_{v}.mat'.format(m=mmm,v=vv),aac)

     fp_outc = fp_out+'{}/'.format(mmm+sub)
     aac_clear = RL.read_aac(fp_outc,fp_mat,mmm=mmm)
     sio.savemat('/nobackup/sleblan2/AAC_DARF/mat/v3c_faac/AAC_{m}_{v}.mat'.format(m=mmm+subi,v=vv),aac_clear)


# In[ ]:


def build_aac_input(fp,fp_alb,fp_out,fp_pmom=None,fp_uvspec='/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec',fp_output=None,
                    wvl_file_sol=None,wvl_file_thm=None,aero_clear=False,version='v1',stdfac_dict={},max_nmom=None,list_only=False,
                    start_lat=None,start_lon=None,mmmlist=['DJF','MAM','JJA','SON']):
    """
    Purpose:
    
        Program to build the inputs of libradtran for Meloë's AAC study
    
    Inputs:
    
        fp: path of directory to matlab input files
        fp_alb: full path to where (without the filename) of the MODIS albedo
        fp_out: full path to where the input files will be saved
        fp_pmom: full path (without the filename) to pmom files
        fp_uvspec: full path to the uvspec program, defaults to : /u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec
        fp_output: path to output of uvspec, if none, the fp_out is used, with the last assumed directory 
                    /input/ to be changed to /output/
        wvl_file_sol: full path to solar wavelength file (wavelengths in nm in second column)
        wvl_file_thm: full path of thermal wavelength file (wavelengths in nm in second column)
        aero_clear: if set to True, then aerosol extinction is set to zero in all cases. (defaults to False) 
        version: (defaults to v1) version number of the files for tracking
        stdfac_dict: Dict that contains the multiplicative factors to the standard deviation stored in the keys:
                     'ext','ssa','asym', 'COD','ref'
        max_nmom: Maximum number of phase function moments to write out. (defaults to none)
        list_only: (default False) When True, only write out the list file (for debugging)
        start_lat: (default None) if set to an integer, uses the index to start the latitutde loops creating the files (for debugging)
        start_lon: (default None) if set to an integer, uses the index to start the longitutde loops creating the files (for debugging)
        mmmlist: (default ['DJF','MAM','JJA','SON']) the list of months to go through
        
    Dependencies:
    
        numpy
        scipy
        load_utils
        Run_libradtran 
        os
        
    Required files:
    
        Input_to_DARF_mmm.mat : input files from Meloë
        surface albedo files from MODIS (MCD43GF_geo_shortwave_doy_2007.hdf)
        pmom files in netcdf (thermal and solar)
        
    Example:
    
        ...
        
    Modification History:
    
        Written: Samuel LeBlanc, 2015-07-01, Happy Canada Day!, Santa Cruz, CA
        Modified: Samuel LeBlanc, 2015-07-07, NASA Ames
                - Modified the calling paths to include fp_pmom
                - Added comments
                - Changed out of Prepare_input_aac to Run_libradtran
        Modified: Samuel LeBlanc, 2015-07-08, Santa Cruz, CA
                - added creator of list file (for running on cluster), creating list file in the path described by fp_out
                - added fp_uvspec 
                - added fp_output for list file creation of uvspec output
                - changed to use the base_file method
                - changed to use wavelength_files
        Modified: Samuel LeBlanc, 2015-07-09, NASA Ames, CA
                - using one set of cloud files per lat lon combinations, not for every hour
                - added aero_clear keyword to define clear value
        Modified: Samuel LeBlanc, 2015-07-10, Santa Cruz, CA
                - changed to have seperate path and list file for each mmm
        MOdified: Samuel LeBlanc, 2017-03-06, Santa Cruz, CA
                - added the version keyword for version tracking
        
    """
    import numpy as np
    import scipy.io as sio
    import Run_libradtran as RL
    import load_utils as lm
    import os
    if fp_pmom:
        pmom_solar = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='solar')
        pmom_thermal = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='thermal')
    else:
        pmom_solar = RL.make_pmom_inputs(source='solar')
        pmom_thermal = RL.make_pmom_inputs(source='thermal')
    if max_nmom:
        pmom_solar['max_nmoms'] = max_nmom
        pmom_thermal['max_nmoms'] = max_nmom
    geo = {'zout':[0,3,100],'year':2007,'day':15,'minute':0,'second':0}
    aero = {'z_arr':[3.0,4.0]}
    cloud = {'ztop':3.0,'zbot':2.0,'phase':'wc','write_moments_file':True}
    source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/','run_fuliou':True}
    albedo = {'create_albedo_file':False}
    if not fp_output:
        change_fp_output = True
    else:
        change_fp_output = False
    
    stdfac_dict = RL.merge_dicts({'ext':0.0,'ssa':0.0,'asym':0.0,
                                  'COD':0.0,'ref':0.0},stdfac_dict)
    std_label = ''
    if aero_clear:
        std_label = '_clear'
    for k in stdfac_dict.keys():
        if stdfac_dict[k] != 0.0:
            if stdfac_dict[k]>0.0: 
                n='p' 
            else: n='m'
            std_label = std_label+'_'+k+nstdfac_dict
    
    for mmm in mmmlist: #['DJF','MAM','JJA','SON']:
        fpm = fp+'Input_to_DARF_{mmm}_{vv}.mat'.format(mmm=mmm,vv=version)
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
            
        try:
            file_list = file(fp_out+'AAC_list_file_{m}_{v}{lbl}.sh'.format(m=mmm,v=version,lbl=std_label),'w')
        except Exception,e:
            print 'Problem with accessing file, return Exception: ',e
            return
        print 'Starting list file'
        fp_out2 = fp_out+mmm+std_label+'/'
        if not os.path.exists(fp_out2):
            os.mkdir(fp_out2)
        if change_fp_output:
            fp_output = fp_out2.replace('input','output')
            if not os.path.exists(fp_output):
                os.mkdir(fp_output)
        fp_base_file = fp_out2+'base.inp'
        make_base = True

        fpa = fp_alb+'MCD43GF_geo_shortwave_%03i_2007.hdf' % doy
        print 'Getting albedo files: '+fpa
        alb_geo,alb_geo_dict = lm.load_hdf_sd(fpa)
        print 'done loading albedo files'
        alb_geo_sub = np.nanmean(np.nanmean(alb_geo['MCD43GF_CMG'].reshape([48,21600/48,75,43200/75]),3),1)
        alb_geo_lat = np.linspace(90,-90,num=48)
        alb_geo_lon = np.linspace(-180,180,num=75)
        
        print 'Running through the files'
        for ilat,lat in enumerate(input_mmm['MODIS_lat'][0,0]):
            if start_lat:
                if not ilat>=start_lat:
                    continue
            for ilon,lon in enumerate(input_mmm['MODIS_lon'][0,0]):
                if start_lon:
                    if not ilon>=start_lon:
                        continue
                geo['lat'],geo['lon'] = lat,lon
                # set the aerosol values
                aero['wvl_arr'] = input_mmm['MOC_wavelengths'][0,0][0,:]*1000.0
                aero['ext'] = np.abs(input_mmm['MOC_ext_mean'][0,0][ilat,ilon,:]) # is in AOD units, not ext, but since it is for 1 km, does not matter
                print 'ext:',stdfac_dict['ext'],aero['ext'][9],stdfac_dict['ext']*np.abs(input_mmm['MOC_ext_std'][0,0][ilat,ilon,9])/1000.0
                aero['ext'] = aero['ext']+stdfac_dict['ext']*np.abs(input_mmm['MOC_ext_std'][0,0][ilat,ilon,:])/1000.0 #convert  per Mm to per km
                aero['ext'][aero['ext']<0.0] = 0.0
                if aero_clear:
                    aero['ext'] = aero['ext']*0.0
                if np.isnan(aero['ext']).all():
                    print 'skipping lat:%i, lon:%i' % (ilat,ilon)
                    continue
                aero['ssa'] = input_mmm['MOC_ssa_mean'][0,0][ilat,ilon,:]
                aero['ssa'] = aero['ssa']+stdfac_dict['ssa']*input_mmm['MOC_ssa_std'][0,0][ilat,ilon,:]
                aero['asy'] = input_mmm['MOC_asym_mean'][0,0][ilat,ilon,:]
                aero['asy'] = aero['asy']+stdfac_dict['asym']*input_mmm['MOC_asym_std'][0,0][ilat,ilon,:]
                
                #sanitize inputs after adding subtracting standard deviations
                try: aero['ssa'][aero['ssa']<0.0] = 0.0
                except: pass
                try: aero['ssa'][aero['ssa']>1.0] = 1.0
                except: pass
                try: aero['asy'][aero['asy']<0.0] = 0.0
                except: pass
                try: aero['asy'][aero['asy']>1.0] = 1.0
                except: pass
                
                if aero['wvl_arr'].max()<100000.0:
                    aero['wvl_arr'] = np.append(aero['wvl_arr'],100000.0)
                    aero['ext'] = np.append(aero['ext'],aero['ext'][-1])
                    aero['ssa'] = np.append(aero['ssa'],aero['ssa'][-1])
                    aero['asy'] = np.append(aero['asy'],aero['asy'][-1])
                # set the cloud values
                cloud['tau'] = input_mmm['MODIS_COD_mean'][0,0][ilat,ilon]
                cloud['tau'] = cloud['tau']+stdfac_dict['COD']*input_mmm['MODIS_COD_std'][0,0][ilat,ilon]
                cloud['ref'] = input_mmm['MODIS_effrad_mean'][0,0][ilat,ilon]
                cloud['ref'] = cloud['ref']+stdfac_dict['ref']*input_mmm['MODIS_Effrad_std'][0,0][ilat,ilon]
                try: cloud['tau'][cloud['tau']<0.0] = 0.0
                except: pass
                try: cloud['ref'][cloud['ref']<2.0] = 2.0
                except: pass
                
                # set the albedo
                alb = alb_geo_sub[np.argmin(abs(alb_geo_lat-lat)),np.argmin(abs(alb_geo_lon-lon))]
                if np.isnan(alb): 
                    albedo['sea_surface_albedo'] = True
                else:
                    albedo['albedo'] = alb
                    
                cloud['link_to_mom_file'] = False
                aero['link_to_mom_file'] = False
                cloud_file_name_sol = fp_out2+'AAC_input_lat%02i_lon%02i_%s_sol.inp_cloud' % (ilat,ilon,mmm)
                cloud_file_name_thm = fp_out2+'AAC_input_lat%02i_lon%02i_%s_thm.inp_cloud' % (ilat,ilon,mmm)
                aero['file_name'] = fp_out2+'AAC_input_lat%02i_lon%02i_%s_thm.inp_aero' % (ilat,ilon,mmm)

                for HH in xrange(24):
                    geo['hour'] = HH
                    #build the solar input file
                    source['source'] = 'solar'
                    if wvl_file_sol:
                        source['wvl_filename'] = wvl_file_sol
                    else:
                        source['wvl_range'] = [250,5600]
                        source['wvl_filename'] = None
                    cloud['moms_dict'] = pmom_solar
                    cloud['file_name'] = cloud_file_name_sol
                    file_out_sol = fp_out2+'AAC_input_lat%02i_lon%02i_%s_HH%02i_sol.inp' % (ilat,ilon,mmm,HH)
                    if not list_only:
                        RL.write_input_aac(file_out_sol,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                                       make_base=make_base,fp_base_file=fp_base_file,set_quiet=True,solver='rodents')
                    if make_base:
                        make_base = False
                    #build the thermal input file
                    source['source'] = 'thermal'
                    if wvl_file_thm:
                        source['wvl_filename'] = wvl_file_thm
                    else:
                        source['wvl_range'] = [4000,50000-1]
                        source['wvl_filename'] = None
                    cloud['moms_dict'] = pmom_thermal
                    cloud['file_name'] = cloud_file_name_thm
                    file_out_thm = fp_out2+'AAC_input_lat%02i_lon%02i_%s_HH%02i_thm.inp' % (ilat,ilon,mmm,HH)
                    
                    if not list_only:
                        RL.write_input_aac(file_out_thm,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                                       make_base=False,fp_base_file=fp_base_file,set_quiet=True,solver='rodents')
                    file_list.write(fp_uvspec+' < '+file_out_sol+' > '+fp_output
                                    +'AAC_input_lat%02i_lon%02i_%s_HH%02i_sol.out\n' % (ilat,ilon,mmm,HH))
                    file_list.write(fp_uvspec+' < '+file_out_thm+' > '+fp_output
                                    +'AAC_input_lat%02i_lon%02i_%s_HH%02i_thm.out\n' % (ilat,ilon,mmm,HH))
                    if not cloud['link_to_mom_file']:
                        cloud['link_to_mom_file'] = True
                    if not aero['link_to_mom_file']:
                        aero['link_to_mom_file'] = True
                    print mmm,ilat,ilon,HH
        del alb_geo
        del input_mmm
        file_list.close()


# In[ ]:


def read_aac(fp_out,fp_mat,mmm=None,read_sol=True,read_thm=True):
    """
    Purpose:
    
        Simple program to read all the output of the libradtran runs for AAC
        Program to read the output of libradtran irradiance files
        Very simple output of 3 zout, one wavelength
    
    Inputs:
    
        fp_out: full path of the directory with the files to read
        fp_mat: full path of mat file with lat and lon to use
        mmm: string with the season defined (can be DJF,MAM,JJA, or SON)
        read_sol: set to True (default) to read the sol files
        read_thm: set to True (default) to read the thm files
        
    Outputs:
    
        out: dictionary with the saved output 
        
    Dependencies:
    
        Run_libradtran (this file)
        os
        numpy
        scipy
        
    Required files:
    
        files of libradtran output
        meloë's .mat files
        
    Example:
    
        >>> import Run_libradtran as RL
        >>> fp_out = '/nobackup/sleblan2/AAC_DARF/output/v2/DJF/'
        >>> fp_mat = '/u/sleblan2/meloe_AAC/Input_to_DARF_DJF.mat'
        >>> mmm = 'DJF'
        >>> DJF = RL.read_aac(fp_out,fp_mat,mmm=mmm)
        
        DJF 0 0
        File not found skip: lat00_lon00_DJF_HH00 
        .
        .
        .
        
        >>> DJF.keys()
        ['UTC', 
         'LW_irr_up_utc', 
         'lon', 
         'SW_irr_dn_avg', 
         'SW_irr_up_utc',
         'LW_irr_dn_avg', 
         'zout', 
         'LW_irr_up_avg', 
         'lat', 
         'SW_irr_up_avg', 
         'SW_irr_dn_utc', 
         'LW_irr_dn_utc']
 
    Modification History:
    
        Written: Samuel LeBlanc, 2015-07-13, Santa Cruz, CA
        Modified: Samuel LeBlanc, 2015-07-15, Santa Cruz, CA
                    - added read_sol and read_thm keywords
                    - added example
        
    """
    import os
    import scipy.io as sio
    import Run_libradtran as RL
    import numpy as np
    
    if not mmm:
        raise NameError('no season string defined')
        return
    
    input_mmm = sio.loadmat(fp_mat,mat_dtype=True)['data_input_darf']
    output = {'lat':input_mmm['MODIS_lat'][0,0],
              'lon':input_mmm['MODIS_lon'][0,0],
              'UTC':range(24),
              'zout':[0,3,100]}
    nlat,nlon,nz = len(output['lat']),len(output['lon']),len(output['zout'])
    output['SW_irr_dn_utc'] = np.zeros((nz,nlat,nlon,24))
    output['SW_irr_up_utc'] = np.zeros((nz,nlat,nlon,24))
    output['LW_irr_dn_utc'] = np.zeros((nz,nlat,nlon,24))
    output['LW_irr_up_utc'] = np.zeros((nz,nlat,nlon,24))
    output['SW_irr_dn_avg'] = np.zeros((nz,nlat,nlon))
    output['SW_irr_up_avg'] = np.zeros((nz,nlat,nlon))
    output['LW_irr_dn_avg'] = np.zeros((nz,nlat,nlon))
    output['LW_irr_up_avg'] = np.zeros((nz,nlat,nlon))    
    for ilat in xrange(nlat):
        for ilon in xrange(nlon):        
            for iutc in output['UTC']:
                file_out_sol = fp_out+'AAC_input_lat%02i_lon%02i_%s_HH%02i_sol.out' % (ilat,ilon,mmm,iutc)
                file_out_thm = fp_out+'AAC_input_lat%02i_lon%02i_%s_HH%02i_thm.out' % (ilat,ilon,mmm,iutc)
                try:
                    if read_sol:
                        sol = RL.read_libradtran(file_out_sol,zout=output['zout'])
                    if read_thm:
                        thm = RL.read_libradtran(file_out_thm,zout=output['zout'])
                except IOError:
                    print 'File not found skip: lat%02i_lon%02i_%s_HH%02i' %(ilat,ilon,mmm,iutc)
                    if iutc==0:
                        print file_out_sol
                    continue
                except ValueError:
                    print 'Problem with file: lat%02i_lon%02i_%s_HH%02i' %(ilat,ilon,mmm,iutc)
                    output['SW_irr_dn_utc'][:,ilat,ilon,iutc] = np.nan
                    output['SW_irr_up_utc'][:,ilat,ilon,iutc] = np.nan
                    output['LW_irr_dn_utc'][:,ilat,ilon,iutc] = np.nan
                    output['LW_irr_up_utc'][:,ilat,ilon,iutc] = np.nan
                    continue
                if read_sol:
                    output['SW_irr_dn_utc'][:,ilat,ilon,iutc] = sol['direct_down']+sol['diffuse_down']
                    output['SW_irr_up_utc'][:,ilat,ilon,iutc] = sol['diffuse_up']
                if read_thm:
                    output['LW_irr_dn_utc'][:,ilat,ilon,iutc] = thm['direct_down']+thm['diffuse_down']
                    output['LW_irr_up_utc'][:,ilat,ilon,iutc] = thm['diffuse_up']
            print mmm,ilat,ilon
            output['SW_irr_dn_avg'][:,ilat,ilon] = np.mean(output['SW_irr_dn_utc'][:,ilat,ilon,:],axis=1)
            output['SW_irr_up_avg'][:,ilat,ilon] = np.mean(output['SW_irr_up_utc'][:,ilat,ilon,:],axis=1)
            output['LW_irr_dn_avg'][:,ilat,ilon] = np.mean(output['LW_irr_dn_utc'][:,ilat,ilon,:],axis=1)
            output['LW_irr_up_avg'][:,ilat,ilon] = np.mean(output['LW_irr_up_utc'][:,ilat,ilon,:],axis=1)
    return output

