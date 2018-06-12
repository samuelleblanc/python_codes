
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

# In[9]:


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
        tmp_folder = ''
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
        read_fuliou_aac(fp_out,ms=ms,tmp_folder=tmp_folder,vn=vn)
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
    fp_fuliou = '/home5/sleblan2/fuliou/v20180201/fuliou'

    stdfac = {'COD':+0.0}
    build_aac_FLinput(fp=fp,fp_alb=fp_alb,fp_out=fp_out,fp_fuliou=fp_fuliou,version=vv,stdfac_dict=stdfac,list_only=False,aero_clear=doclear,mmmlist=ms)


# For reading

# In[2]:


def read_fuliou_aac(fp_out,ms=['SON'],tmp_folder='',vn='v5_fuliou'):
    
    from DARE_AAC_fuliou import read_aac_fl
    import scipy.io as sio

    from DARE_AAC_fuliou import build_aac_FLinput
    print 'working on folder: '+fp_out

    vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'
    fp = '/u/sleblan2/meloe_AAC/{vn}/'.format(vn=vn)
    fp_alb = '/nobackup/sleblan2/AAC_DARF/surface_albedo/'
    fp_fuliou = '/home5/sleblan2/fuliou/v20180201/fuliou'
    
    for mmm in ms:
        fp_mat = '/u/sleblan2/meloe_AAC/{vn}/Input_to_DARF_{m}_{v}.mat'.format(m=mmm,v=vv,vn=vn)
            
        fp_outc = fp_out+'{}/'.format(mmm)
        aac = read_aac_fl(fp_outc,fp,mmmlist=mmm)
        sio.savemat('/nobackup/sleblan2/AAC_DARF/mat/{vn}/AAC_DARE_FuLiou_{m}_{v}.mat'.format(m=mmm,v=vv,vn=vn),aac)


# In[7]:


def build_aac_FLinput(fp,fp_alb,fp_out,fp_fuliou='/home5/sleblan2/fuliou/v20180201/fuliou',fp_output=None,
                    aero_clear=False,version='v1',stdfac_dict={},list_only=False,
                    start_lat=None,start_lon=None,mmmlist=['DJF','MAM','JJA','SON']):
    """
    Purpose:
    
        Program to build the inputs of libradtran for Meloë's AAC study 
        Subset for using the Fu Liou codes, without clouds underneath
    
    Inputs:
    
        fp: path of directory to matlab input files
        fp_alb: full path to where (without the filename) of the MODIS albedo
        fp_out: full path to where the input files will be saved
        fp_fuliou: full path to the fuliou program, defaults to : home5/sleblan2/fuliou/v20180201/fuliou
        fp_output: path to output of fuliou, if none, the fp_out is used, with the last assumed directory 
                    /input/ to be changed to /output/
        aero_clear: if set to True, then aerosol extinction is set to zero in all cases. (defaults to False) 
        version: (defaults to v1) version number of the files for tracking
        stdfac_dict: Dict that contains the multiplicative factors to the standard deviation stored in the keys:
                     'ext','ssa','asym', 'COD','ref'
        list_only: (default False) When True, only write out the list file (for debugging)
        start_lat: (default None) if set to an integer, uses the index to start the latitutde loops creating the files (for debugging)
        start_lon: (default None) if set to an integer, uses the index to start the longitutde loops creating the files (for debugging)
        mmmlist: (default ['DJF','MAM','JJA','SON']) the list of months to go through
        
    Dependencies:
    
        numpy
        scipy
        load_utils
        Run_libradtran
        Run_fuliou
        os
        
    Required files:
    
        Input_to_DARF_mmm.mat : input files from Meloë
        surface albedo files from MODIS (MCD43GF_geo_shortwave_doy_2007.hdf)
        
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
                
                aero['ext'] = np.array([aero['ext'],aero['ext']])
                aero['ssa'] = np.array([aero['ssa'],aero['ssa']])
                aero['asy'] = np.array([aero['asy'],aero['asy']])
                
                # set the albedo
                alb = alb_geo_sub[np.argmin(abs(alb_geo_lat-lat)),np.argmin(abs(alb_geo_lon-lon))]
                if np.isnan(alb): 
                    albedo['sea_surface_albedo'] = True
                else:
                    albedo['modis_albedo'] = RF.get_MODIS_surf_albedo(fp_alb,doy,geo['lat'],geo['lon'],year_of_MODIS=2007)
                    albedo['sea_surface_albedo'] = False

                #for HH in xrange(24):
                geo['hour'] = 12.0
                form = {'ilat':ilat,'ilon':ilon,'mmm':mmm}
                file_in = fp_out2+'AACFLin_lat{ilat:02.0f}_lon{ilon:02.0f}_{mmm}.datin'.format(**form)
                file_out = fp_output+'AACFLout_lat{ilat:02.0f}_lon{ilon:02.0f}_{mmm}.wrt'.format(**form)
                    
                if not list_only:
                    RF.write_fuliou_input(file_in,geo=geo,aero=aero,albedo=albedo,verbose=False)
              
                file_list.write(fp_fuliou+' '+file_in+' '+file_out+'\n')
                print mmm,ilat,ilon
        del alb_geo
        del input_mmm
        file_list.close()


# In[1]:


def read_aac_fl(fp_out,fp,mmmlist=['DJF','MAM','JJA','SON'],outstr='AACFLout_lat{ilat:02.0f}_lon{ilon:02.0f}_{mmm}.wrt'):
    """
    Purpose:
    
        Simple program to read all the output of the Fu Liou runs for AAC
        Program to read the output of libradtran irradiance files
        Very simple output of 3 zout, one wavelength
    
    Inputs:
    
        fp_out: full path of the directory with the files to read
        fp: full path of mat file with lat and lon to use and with MOC results
        fp_alb: full path of the  albedo files
        mmmlist: list of strings with the season defined (can be DJF,MAM,JJA, or SON)
        outstr: string format of the file names
        
    Outputs:
    
        out: dictionary with the saved output 
        
    Dependencies:
    
        Run_fuliou
        os
        numpy
        scipy
        load_utils
        
        
    Required files:
    
        files of Fu Liou output
        Input_to_DARF_mmm.mat : input files from Meloë
        albedo files (MODIS gap filled)
        
    Example:
    
 
    Modification History:
    
        Written: Samuel LeBlanc, 2018-06-11, Santa Cruz, CA
                Ported from Run_libradtran to use Fuliou outputs
        Modified: 
        
    """
    import os
    import scipy.io as sio
    import Run_libradtran as RL
    import Run_fuliou as RF
    import numpy as np
    import load_utils as lm
    
    geo = {'zout':[0,3,100],'year':2007,'day':15,'minute':0,'second':0}
    aero = {'z_arr':[3.0,4.0]}
    source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/','run_fuliou':True}
    albedo = {'create_albedo_file':False}
    
    out = {}
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
    
        d = []
        
        print 'Running through the files'
        for ilat,lat in enumerate(input_mmm['MODIS_lat'][0,0]):
            for ilon,lon in enumerate(input_mmm['MODIS_lon'][0,0]):
                # set the aerosol values
                aero['wvl_arr'] = input_mmm['MOC_wavelengths'][0,0][0,:]*1000.0
                aero['ext'] = np.abs(input_mmm['MOC_ext_mean'][0,0][ilat,ilon,:]) # is in AOD units, not ext, but since it is for 1 km, does not matter
                print 'ext:',aero['ext'][9]
                aero['ext'][aero['ext']<0.0] = 0.0
                if np.isnan(aero['ext']).all():
                    print 'skipping lat:%i, lon:%i' % (ilat,ilon)
                    continue
                aero['ssa'] = input_mmm['MOC_ssa_mean'][0,0][ilat,ilon,:]
                aero['asy'] = input_mmm['MOC_asym_mean'][0,0][ilat,ilon,:]
                
                #sanitize inputs after adding subtracting standard deviations
                try: aero['ssa'][aero['ssa']<0.0] = 0.0
                except: pass
                try: aero['ssa'][aero['ssa']>1.0] = 1.0
                except: pass
                try: aero['asy'][aero['asy']<0.0] = 0.0
                except: pass
                try: aero['asy'][aero['asy']>1.0] = 1.0
                except: pass

                #for HH in xrange(24):
                geo['hour'] = 12.0
                form = {'ilat':ilat,'ilon':ilon,'mmm':mmm}
                file_out = fp_out+outstr.format(**form)
    
                d.append(RF.read_fuliou_output(file_out))
                d[i]['ssa'] = np.array(aero['ssa'])
                d[i]['asy'] = np.array(aero['asy'])
                d[i]['ext'] = np.array(aero['ext'])
                d[i]['wvl'] = np.array(aero['wvl_arr'])
                d[i]['lat'],d[i]['lon'] = lat,lon
                RF.analyse_fuliou_output(d[i])
        out[mmm] = d
    return out

