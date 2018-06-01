
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

# In[4]:


if __name__ == '__main__':
    import argparse
    import os
    from DARE_AAC_fuliou import write_fuliou_aac, read_fuliou_aac
    long_description = """    Prepare or read the fuliou DARE AAC calculations for Meloe's study"""
    parser = argparse.ArgumentParser(description=long_description)
    parser.add_argument('-doread','--doread',help='if set, will only read the output, not produce them',
                        action='store_true')
    parser.add_argument('-dowrite','--dowrite',help='if set, will write the input and list files for fuliou',
                        action='store_true')
    parser.add_argument('mmm',nargs='?',help='3 letter season designator (DJF, MAM, JJA, SON)')
    parser.add_argument('sub',nargs='?',help='string designating the subset (ssam,ssap,...)')
    parser.add_argument('-c','--clear',help='if set, will write the input for no aerosol case',
                        action='store_true',default=False)
    parser.add_argument('-t','--tmp_folder',nargs='?',help='The path to the temp directory to use')
    parser.add_argument('-v','--vn',nargs='?',help='The version identifier (subfolder names)')
    in_ = vars(parser.parse_args())
    
    doread = in_.get('doread',False)
    dowrite = in_.get('dowrite',False)
    doclear = in_.get('clear',False)
    if doclear:
        print 'doing the clear no aerosol version'
    ms = in_.get('mmm',False)
    if not ms:
        ms = ['DJF','MAM','JJA','SON']
    else:
        ms = [ms]
    print 'doing the season:',ms

    vn = in_.get('vn','v5_fuliou')
    
    if in_.get('tmp_folder'):
        tmp_folder = in_.get('tmp_folder')
        fp_out = tmp_folder+'/AAC_DARF/input/{vn}/'.format(vn=vn)
    else:
        fp_out = '/nobackup/sleblan2/AAC_DARF/input/{vn}/'.format(vn=vn)
    
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
    if not os.path.exists(fp_out):
        os.makedirs(fp_out)
    fp_in = fp_out.replace('input','output')
    if not os.path.exists(fp_in):
        os.makedirs(fp_in)

    if dowrite:
        write_fuliou_aac(fp_out,ms=ms,tmp_folder=tmp_folder,doclear=doclear,vn=vn)
    elif doread:
        read_fuliou_aac(fp_out,ms=ms,tmp_folder=tmp_folder,doclear=doclear,vn=vn)
    else:
        write_fuliou_aac(fp_out,ms=ms,tmp_folder=tmp_folder,doclear=doclear,vn=vn)


# For writing

# In[5]:


# Program to prepare and call the build_aac_input
def write_fuliou_aac(fp_out,ms=['SON'],tmp_folder='',doclear=False,vn='v5_fuliou'):
    from DARE_AAC_fuliou import build_aac_FLinput
    print 'working on folder: '+fp_out

    vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'
    fp = '/u/sleblan2/meloe_AAC/{vn}/'.format(vn=vn)
    fp_alb = '/nobackup/sleblan2/AAC_DARF/surface_albedo/'
    fp_pmom = '/nobackup/sleblan2/AAC_DARF/rtm/'
    fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'

    stdfac = {'COD':+0.0}
    build_aac_FLinput(fp=fp,fp_alb=fp_alb,fp_out=fp_out,fp_pmom=fp_pmom,version=vv,stdfac_dict=stdfac,list_only=False,aero_clear=doclear,mmmlist=ms)


# For reading

# In[6]:


def read_fuliou_aac(fp_out,ms=['SON'],tmp_folder='',doclear=False,vn='v5_fuliou'):
    
    from DARE_AAC_fuliou import read_aac_fl
    import scipy.io as sio

    import time
    time.sleep(0.5)

    #vv = 'v3_20160408_CALIOP_AAC_AOD_With_Sa_only'
    #vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'
    vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'
    for mmm in ms:
    #     fp_out = '/nobackup/sleblan2/AAC_DARF/output/v3_without/{}/'.format(mmm)
         fp_mat = '/u/sleblan2/meloe_AAC/{vn}/Input_to_DARF_{m}_{v}.mat'.format(m=mmm,v=vv,vn=vn)
    #     aac = RL.read_aac(fp_out,fp_mat,mmm=mmm)
    #     sio.savemat('/nobackup/sleblan2/AAC_DARF/mat/v3_without/AAC_{m}_{v}.mat'.format(m=mmm,v=vv),aac)

         fp_outc = fp_out+'{}/'.format(mmm+sub)
         aac_clear = read_aac_fl(fp_outc,fp_mat,mmm=mmm)
         sio.savemat('/nobackup/sleblan2/AAC_DARF/mat/{vn}/AAC_{m}_{v}.mat'.format(m=mmm+subi,v=vv,vn=vn),aac_clear)


# In[7]:


def build_aac_FLinput(fp,fp_alb,fp_out,fp_pmom=None,fp_uvspec='/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec',fp_output=None,
                    wvl_file_sol=None,wvl_file_thm=None,aero_clear=False,version='v1',stdfac_dict={},max_nmom=None,list_only=False,
                    start_lat=None,start_lon=None,mmmlist=['DJF','MAM','JJA','SON']):
    """
    Purpose:
    
        Program to build the inputs of libradtran for Meloë's AAC study 
        Subset for using the Fu Liou codes, without clouds underneath
    
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
    
        Written: Samuel LeBlanc, 2015-05-30, Santa Cruz, CA
                 imported from Run_libradtran and ported to use Run_fuliou
        Modified: 
        
    """
    import numpy as np
    import scipy.io as sio
    import Run_libradtran as RL
    import Run_fuliou as RF
    import load_utils as lm
    import os
    
    geo = {'zout':[0,3,100],'year':2007,'day':15,'minute':0,'second':0}
    aero = {'z_arr':[3.0,4.0]}
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
                
                # set the albedo
                alb = alb_geo_sub[np.argmin(abs(alb_geo_lat-lat)),np.argmin(abs(alb_geo_lon-lon))]
                if np.isnan(alb): 
                    albedo['sea_surface_albedo'] = True
                else:
                    albedo['modis_albedo'] = get_MODIS_surf_albedo(fp_alb,doy,geo['lat'],geo['lon'],year_of_MODIS=2007)
                    albedo['sea_surface_albedo'] = False

                for HH in xrange(24):
                    geo['hour'] = HH

                    file_in = fp_out2+'AAC_input_fuliou_lat%02i_lon%02i_%s_HH%02i.datin' % (ilat,ilon,mmm,HH)
                    file_out = fp_out2+'AAC_input_fuliou_lat%02i_lon%02i_%s_HH%02i.wrt' % (ilat,ilon,mmm,HH)
                    
                    if not list_only:
                        RF.write_fuliou_input(file_in,geo=geo,aero=aero,albedo=albedo,verbose=False)
              
                    file_list.write(fp_uvspec+' '+file_in+' '+fp_output
                                    +file_out+'\n' % (ilat,ilon,mmm,HH))
                    print mmm,ilat,ilon,HH
        del alb_geo
        del input_mmm
        file_list.close()


# In[8]:


def read_aac_fl(fp_out,fp_mat,mmm=None,read_sol=True,read_thm=True):
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

