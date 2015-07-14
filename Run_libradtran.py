# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     Run_libradtran
# 
# Purpose:  
# 
#     Python modules that combine the different modules used for writing and reading libradtran files. 
#     Used with libradtran 2.0
#     Same outputs as IDL procedures 'write_inp_mix.pro' but for libradtran 2.0
# 
# Calling Sequence:
# 
#     import Run_libradtran as RL
#   
# Input:
# 
#     none at command line
#     see methods of module
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
#     - math
#     - pdb
#     - datetime
#     - load_modis
#   
# Needed Files:
# 
#   - ...
#   
#   
# Modification History:
# 
#     Wrtten: Samuel LeBlanc, NASA Ames, from Santa Cruz, 2015-06-26

# <codecell>

def write_cloud_file(output_file,tau,ref,zbot,ztop,verbose=False):
    """
    Purpose:

        Program to write out the cloud file profile
        for either ic_file or wc_file in 1D settings
        outputs 3 columns: 
            # z [km]    IWC/LWC [g/m^3]    Reff [um]
            
        ** translates tau and ref to lwc/iwc with the equation LWP = 2/3*tau*ref and LWC=LWP/z_extent*0.001
        ** Writes out file based on layer properties, not level properties. (given value of ref and tau are for a cloud from zbot to ztop)
    
    Input: 
  
        output_file: full path to file to be written
        tau: value of cloud optical thickness
        ref: value of effective radius in microns
        ztop: location in km of cloud top
        zbot: location in km of cloud bottom
    
    Output:

        output_file
                
    
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:

        numpy
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-26, NASA Ames, from Santa Cruz, CA
    """
    import numpy as np
    if (zbot >= ztop):
        raise ValueError('*** Error ztop must be larger than zbot ***')
    if (ztop < 1.0):
        print('ztop is smaller than one, check units, ztop should be in km')
        if verbose:
            print('..file preperations continuing')
    if tau:
        if np.isnan(tau):
            raise ValueError('*** Error tau is NaN ***')
    if ref:
        if np.isnan(ref):
            raise ValueError('*** Error ref is NaN ***')
    lwp = 2.0/3.0*tau*ref
    lwc = lwp/(ztop-zbot)*0.001
    if verbose:
        print('..Cloud water content is: %f' % lwc)
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose:
        print('..printing to file: %s' % output_file)
    output.write('# z [km]    IWC/LWC [g/m^3]    Reff [um] \n')
    output.write('%4.4f\t%4.5f\t%3.2f\n' % (ztop,0,0))
    output.write('%4.4f\t%4.5f\t%3.2f\n' % (zbot,lwc,ref))
    output.close() 
    if verbose:
        print('..File finished write_cloud_file, closed')

# <codecell>

def write_aerosol_file_explicit(output_file,z_arr,ext,ssa,asy,wvl_arr,verbose=False):
    """
    Purpose:

        Writes the file with profile of aerosol files
        writes one file per layer of aerosol properties, which are wavelength dependent
        outputs 2 columns of layer properties: 
            # z [km]     file_path
            
        ** Translates asymmetry parameter into scattering phase function defined by Henyey-Greenstein
    
    Input: 
  
        output_file: full path to file to be written
        ext: value of aerosol extinction coefficient at each altitude and wavelength [alt,wvl]
            To define the top most layer, ext must be zero at that altitude
        ssa: value of aerosol single scattering albedo at each altitude and wavelength [alt,wvl]
        asy: value of aerosol asymmetry parameter at each altitude and wavelength [alt,wvl]
        z_arr: value of altitudes to use with ext,ssa,asy
        wvl_arr: array of wavelengths in nm
    
    Output:

        output_file
                
    
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:

        numpy
        Run_libradtran (this file)
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-26, NASA Ames, from Santa Cruz, CA
    """
    import numpy as np
    from Run_libradtran import write_aerosol_file_explicit_wvl
    two_z = False
    if not len(z_arr)==np.shape(ext)[0]:
        if len(z_arr)==2:
            #sepcial case for only two points on the z array
            two_z = True
        else:
            raise LookupError('*** Error z_arr not same size as ext ***')
    if two_z:
        if not len(wvl_arr)==len(ext):
            raise LookupError('*** Error wvl_arr not same size as ext ***')
    else:
        if not len(wvl_arr)==np.shape(ext)[1]:
            raise LookupError('*** Error wvl_arr not same size as ext ***')
        
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose:
        print('..printing to file: %s' % output_file)
    
    output.write('# z [km] \t file_path\n')
    izs = np.argsort(z_arr)[::-1]
    zs = np.sort(z_arr)[::-1]
    if verbose:
        print('..printing %i lines onto profile file' % len(zs))
    for iz,z in enumerate(zs):
        file_iz = output_file+'_z%03i' % iz
        if two_z:
            if iz>0:
                write_aerosol_file_explicit_wvl(file_iz,wvl_arr,ext,ssa,asy,verbose=verbose)
                output.write('%4.4f\t%s\n' % (z,file_iz))
            else:
                output.write('%4.4f\t%s\n' % (z,'NULL'))
        else:
            if any(ext[izs[iz],:]):
                write_aerosol_file_explicit_wvl(file_iz,wvl_arr,ext[izs[iz],:],ssa[izs[iz],:],asy[izs[iz],:],verbose=verbose)
                output.write('%4.4f\t%s\n' % (z,file_iz))
            else:
                output.write('%4.4f\t%s\n' % (z,'NULL'))
    output.close() 
    if verbose:
        print('..File finished write_aerosol_file_explicit, closed')

# <codecell>

def write_aerosol_file_explicit_wvl(output_file,wvl_arr,ext,ssa,asy,verbose=False):
    """
    Purpose:

        Writes the file with aerosol properties per wavelength
        outputs many columns of layer properties: 
            # wvl[nm]    ext[km^-1]   ssa[unitless]  legendre_moments
            
            where for Henyey-Greenstein the legendre moments consist of:
                1   asymmetry_parameter, ..., asymmetry_parameter^n
        
        ** Translates asymmetry parameter into scattering phase function defined by Henyey-Greenstein
    
    Input: 
  
        output_file: full path to file to be written
        ext: value of aerosol extinction coefficient at each wavelength [wvl]
        ssa: value of aerosol single scattering albedo at each wavelength [wvl]
        asy: value of aerosol asymmetry parameter at each wavelength [wvl]
        wvl_arr: array of wavelengths in nm
    
    Output:

        output_file
                
    
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:

        numpy
        Run_libradtran (this file)
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-26, NASA Ames, from Santa Cruz, CA
    """
    if not len(wvl_arr)==len(ext):
        raise LookupError("ext and wvl_arr don't have the same size")
    if not len(wvl_arr)==len(ssa):
        raise LookupError("ssa and wvl_arr don't have the same size")  
    if not len(wvl_arr)==len(asy):
        raise LookupError("asy and wvl_arr don't have the same size")  
    
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose:
        print('..printing to explicit aerosol wavelength defined file: %s' % output_file)
    output.write('# wvl[nm]    ext[km^-1]   ssa[unitless]  legendre_moments\n')
    for iw,wvl in enumerate(wvl_arr):
        output.write('%f\t%f\t%1.6f\t%1.6f\t%1.6f\n' % (wvl,ext[iw],ssa[iw],1.0,asy[iw]))
    output.close()
    if verbose:
        print('..File finished write_aerosol_file_explicit_wvl, closed')

# <codecell>

def write_cloud_file_moments(output_file,tau,ref,zbot,ztop,moms_dict=None,verbose=False):
    """
    Purpose:

        Writes the file with profile of cloud moment files
        writes one file per layer of cloud properties, which are wavelength dependent
        outputs 2 columns of layer properties: 
            # z [km]     file_path

    Input: 
  
        output_file: full path to file to be written
        tau: value of cloud optical thickness
        ref: value of cloud particles effective radius [microns]
        zbot: bottom of cloud layer [km]
        ztop: top of cloud layer [km]
        moms_dict: dictionary from saved legendre moments. 
                    Includes: ntheta, pmom, rho, nmom, ssa, nim, nre, ext, wvl, phase, theta, ref
    
    Output:

        output_file              
    
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:

        numpy
        Run_libradtran (this file)
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-29, NASA Ames, from Santa Cruz, CA
    """
    import numpy as np
    from Run_libradtran import write_cloud_file_moments_wvl, get_cloud_ext_ssa_moms
        
    if (zbot >= ztop):
        raise ValueError('*** Error ztop must be larger than zbot ***')
    if (ztop < 1.0):
        print('ztop is smaller than one, check units, ztop should be in km')
        if verbose:
            print('..file preperations continuing')
            
    lwp = 2.0/3.0*tau*ref
    lwc = lwp/(ztop-zbot)*0.001
    ext,ssa,wvl,moments,nmom = get_cloud_ext_ssa_moms(ref,lwc,moms_dict=moms_dict,verbose=False)
    
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose:
        print('..printing to file: %s' % output_file)

    output.write('# z [km] \t file_path\n')
    output.write('%4.4f\t%s\n' % (ztop,'NULL'))
    file_cloud = output_file+'_zbot'
    output.write('%4.4f\t%s\n' % (zbot,file_cloud))
    
    write_cloud_file_moments_wvl(file_cloud,wvl,ext,ssa,moments,nmom,verbose=verbose)
    
    output.close() 
    if verbose:
        print('..File finished write_cloud_file_moments, closed')

# <codecell>

def write_cloud_file_moments_wvl(output_file,wvl,ext,ssa,moments,nmom,verbose=False):
    """
    Purpose:

        Writes the file with cloud properties per wavelength
        outputs many columns of layer properties: 
            # wvl[nm]    ext[km^-1]   ssa[unitless]  legendre_moments
    
    Input: 
  
        output_file: full path to file to be written
        ext: value of aerosol extinction coefficient at each wavelength [wvl]
        ssa: value of aerosol single scattering albedo at each wavelength [wvl]
        moments: array of moments at each wavelength [wvl]
        nmom: number of moments to be written out. 
        wvl: array of wavelengths in nm
    
    Output:

        output_file
                
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:

        numpy
        Run_libradtran (this file)
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-29, NASA Ames, from Santa Cruz, CA
    """
    if not len(wvl)==len(ext):
        raise LookupError("ext and wvl_arr don't have the same size")
    if not len(wvl)==len(ssa):
        raise LookupError("ssa and wvl_arr don't have the same size")  
        
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose:
        print('..printing to cloud moments properties wavelength defined file: %s' % output_file)
    output.write('# wvl[nm]    ext[km^-1]   ssa[unitless]  legendre_moments\n')
    for iw,wv in enumerate(wvl):
        try:
            output.write('%f\t%f\t%1.6f\t%s \n' % 
                         (wv,ext[iw],ssa[iw]," ".join([str(x/(2.0*i+1.0)) for i,x in enumerate(moments[iw,0:nmom[iw]])])))
        except:
            import pdb
            pdb.set_trace()
    output.close()
    if verbose:
        print('..File write_cloud_file_moments_wvl finished, closed')

# <codecell>

def write_albedo_file(output_file,wvl=[],alb=[],verbose=False):
    """
    Purpose:

        Writes the file with surface albedo per wavelength
        outputs 2 column: 
            # wvl[nm]    albedo[unitless]
            
    Input: 
  
        output_file: full path to file to be written
        wvl: wavelength array in nm
        alb: value of albedo at each wavelength, unitless
        
    Output:

        output_file of albedo values 
                    
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:

        none
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-26, NASA Ames, from Santa Cruz, CA
    """
    if not len(wvl)==len(alb):
        raise LookupError("wvl and alb don't have the same size")
    
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose:
        print('..printing to albedo wavelength defined file: %s' % output_file)
    output.write('# wvl[nm]    alb[unitless] \n')
    for iw,w in enumerate(wvl):
        if not w>100.0:
            print '** Possible error with wavelenght, values below 100 nm, check units **'
        output.write('%f\t%f\n' % (w,alb[iw]))
    output.close()
    if verbose:
        print('..File write_albedo_file finished, closed')

# <codecell>

def get_cloud_ext_ssa_moms(ref,lwc,moms_dict=None,verbose=False):
    """
    Purpose:

        Extracts the moments, extinction, ssa 
        
    Input: 
  
        ref: effective radius of cloud particles [microns]
        lwc: liquid/ice water content in cloud [g/m^3]
        moms_dict: dictionary from saved legendre moments. 
                    Includes: ntheta, pmom, rho, nmom, ssa, nim, nre, ext, wvl, phase, theta, ref
        
    Output:

        ext,ssa,wvl,moments,nmom
        ext: extinction coefficient of cloud per wavelength [km^-1]
        ssa: single scattering albedo of cloud per wavelength [unitless]
        wvl: wavelength array of returned properties [nm]
        moments: array of legendre moments at each wavelength
        nmom: number of legendre moments
                
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:

        numpy
        sys
        scipy.io (for loading idl files)
    
    Required files:
   
        none
    
    Example:

        ext,ssa,wvl,moments,nmom = get_cloud_ext_ssa_moms(ref,lwc,moms_dict=moms_dict,verbose=True)
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-30, NASA Ames, from Santa Cruz, CA
    """
    import numpy as np
    if verbose:
        print '..Getting the moments and extinction coefficients from the moms_dict of cloud properties'
    
    if not moms_dict:
        import sys
        if sys.platform=='win32':
            fdict = 'C:/Users/sleblan2/Research/4STAR/rtm_dat/mie_hi.out'
        else:
            fdict = '/u/sleblan2/4STAR/rtm_dat/mie_hi.out'
        print 'No moments dict defined, loading defaults: '+fdict
        import scipy.io as sio
        moms_dict = sio.idl.read(fdict)
    ir = np.argmin(abs(moms_dict['ref']-ref))
    ext = moms_dict['ext'][ir,:]*lwc
    wvl = moms_dict['wvl']
    if wvl[0]<1:
        wvl = wvl*1000.0

    return ext,moms_dict['ssa'][ir,:],wvl,moms_dict['pmom'][ir,:,:],moms_dict['nmom'][ir,:]

# <codecell>

def write_input_aac(output_file,geo={},aero={},cloud={},source={},albedo={},
                    verbose=False,make_base=False,fp_base_file=None,set_quiet=True):
    """
    Name:

        write_input_aac
    
    Purpose:

        Writes the libradtran input file with defined defaults, see below
        outputs libradtran input files in ascii format
        
        ** There is a multitude of other default parameters. These should be changed for any other type of input files to be written
    
    Input: 
  
        output_file: full path to file to be written
        geo: dictionary with geometry details
            sza: solar zenith angle (libradtran default is 0)
            lat: latitude
            lon: longitude
            doy: day of year
            year: year (YYYY format)
            month: month (MM format)
            day: day of the month (DD format)
            hour: hour of the day, 24h format, UTC
            minute: minutes of the hour
            second: seconds of the minute
            zout: at what altitude should be outputted, in km, default is 0 and 100
        aero: dictionary with aerosol properties
            ext: value of aerosol extinction coefficient at each altitude and wavelength [alt,wvl]
              To define the top most layer, ext must be zero at that altitude
            ssa: value of aerosol single scattering albedo at each altitude and wavelength [alt,wvl]
            asy: value of aerosol asymmetry parameter at each altitude and wavelength [alt,wvl]
            z_arr: value of altitudes to use with ext,ssa,asy
            wvl_arr: array of wavelengths in nm
            link_to_mom_file: if True then no moments file is written out, but it is referenced via file_name saved to aero dict
            file_name: file name and paths of explicit file for aerosol defined. 
                        By default it is created in this program if it is not set, if link_to_mom_file is set, this must be defined
        cloud: dictionary with cloud properties
            tau: value of cloud optical thickness
            ref: value of effective radius in microns
            ztop: location in km of cloud top
            zbot: location in km of cloud bottom
            phase: either 'ic' for ice cloud or 'wc' for water cloud
            write_moments_file: (default False) if True, writes out moments into an ascii file instead 
              of reading directly from the netcdf file. Requires the moments to be put in the cloud dict and the nmom.
            moms_dict: is the moment dict structure returned from the mie_hi.out idl.readsav. To be used when building the moments file
            link_to_mom_file: if True then no moments file is written out, but it is referenced via file_name saved to cloud dict
            file_name: file name and paths of moments file for cloud defined. 
                        By default it is created in this program if it is not set, if link_to_mom_file is set, this must be defined
        source: dictionary with source properties
            wvl_range: range of wavelengths to model (default [202,500])
            source: can either be thermal or solar
            dat_path: data path to be used. Defaults to pleaides values (/u/sleblan2/libradtran/libRadtran-2.0-beta/data/)
            integrate_values: if set to True (default), then the resulting output parameters are integrated over the wavelength range
                            if set to False, returns per_nm irradiance values
            wvl_filename: filename and path of wavelength file (second column has wavelengh in nm to be used)
            run_fuliou: if set to True, then runs fu liou instead of sbdart (default is False)
        albedo: dictionary with albedo properties
            create_albedo_file: if true then albedo file is created with the properties defined by alb_wvl and alb (defaults to False)
            albedo_file: path of albedo file to use if already created 
            albedo: value of albedo. Only used if create_albedo_file is false and albedo_file is empty (defaults to 0.29 - Earth's average)
            alb: wavelength dependent value of albedo for use when create_albedo_file is set to True
            alb_wvl: wavelength grid of albedo for use when create_albedo_file is set to True
            sea_surface_albedo: (default False) If True, sets the sea surface to be parameterized by cox and munk, 
                            requires wind_speed to be set.
            wind_speed: [m/s] (default to 10 m/s) Only used if sea_surface_albedo is set to True. 
                            wind speed over water for parameterization of cox_and_munk
        make_base: boolean to set if the base file is to be written out
                    if False, no base file is saved
        fp_base_file: full file path for base file. 
                    If set to a file path and make_base to False, then include path is printed to input file. 
        set_quiet: if True then quiet is set in the input file, if False, quiet is not set. (default True)
    
    Output:

        output_file
    
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:

        numpy
        Run_libradtran (this file)
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-26, NASA Ames, from Santa Cruz, CA
        Modified: Samuel LeBlanc, 2015-06-29, NASA Ames, from Santa Cruz, CA
                - added writing of moments for cloud properties file
        Modified: Samuel LeBlanc, 2015-07-01, NASA Ames, Happy Canada Day!
                - added sea-surface albedo parameterization
        Modified: Samuel LeBlanc, 2015-07-08, Santa Cruz, CA
                - added base file writing
                - added wvl_filename which is to be writen out
        Modified: Samuel LeBlanc, 2015-07-09, NASA Ames, CA
                - possibility of using the fuliou codes instead of sbdart
                - added link_to_mom_file to cloud and aero dict
                - added file_name in cloud and aero dict
    """
    import numpy as np
    from Run_libradtran import write_aerosol_file_explicit,write_cloud_file,write_albedo_file,merge_dicts,write_cloud_file_moments
    
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose: print 'Opening input file %s' % output_file
    
    if make_base:
        if fp_base_file:
            base = file(fp_base_file,'w')
            include_base = False
        else:
            print 'Problem: No fp_base_file set, please set before running'
            return
    else:
        if fp_base_file:
            include_base = True
        else:
            include_base = False
    
    if verbose: print '..setting the dicts to defaults'
    source = merge_dicts({'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/',
                          'source':'solar',
                          'wvl_range':[250.0,500.0],
                          'integrate_values':True,
                          'run_fuliou':False},source)
    albedo = merge_dicts({'create_albedo_file':False,'sea_surface_albedo':False,'wind_speed':10,
                          'albedo':0.29},albedo)
    geo = merge_dicts({'zout':[0,100]},geo)
    cloud = merge_dicts({'write_moments_file':False},cloud)
    
    if source.get('source')=='solar':
        source['source'] = 'solar '+source['dat_path']+'solar_flux/kurudz_1.0nm.dat per_nm'
    
    if verbose: print '..write out general default values'
    if make_base:
        if source['run_fuliou']:
            base.write('mol_abs_param fu\n')
        else:
            base.write('mol_abs_param sbdart\n')
        base.write('rte_solver disort\n')
        if set_quiet:
            base.write('quiet\n')
    elif include_base:
        output.write('include %s\n' % fp_base_file)
    else:
        if set_quiet:
            output.write('quiet\n')
        if source['run_fuliou']:
            output.write('mol_abs_param fu\n')
        else:
            output.write('mol_abs_param sbdart\n')
        output.write('rte_solver disort\n')
    
    if verbose: print '..write out source dict values'
    if make_base:
        if source['integrate_values']:
            if source['run_fuliou']:
                base.write('output_process sum\n')
            else:
                base.write('output_process integrate\n')
        else:
            base.write('output_process per_nm\n')
        base.write('data_files_path \t %s\n' % source['dat_path'])
    elif not include_base:
        if source['integrate_values']:
            if source['run_fuliou']:
                output.write('output_process sum\n')
            else:
                output.write('output_process integrate\n')
        else:
            output.write('output_process per_nm\n')
        output.write('data_files_path \t %s\n' % source['dat_path'])
    output.write('source \t %s \n' % source['source'])
    
    if source.get('wvl_filename'):
        output.write('wavelength\t%s\n' % source['wvl_filename'])
    else:
        if source['wvl_range'][0]>source['wvl_range'][1]:
            print 'wvl_range was set inverse, inversing'
            source['wvl_range'] = list(reversed(source['wvl_range'])) 
        if source['wvl_range'][0]<250:
            print 'wvl_range starting too low, setting to 250 nm'
            source['wvl_range'][0] = 250.0
        output.write('wavelength\t%f\t%f\n' % (source['wvl_range'][0],source['wvl_range'][1]))
    
    if verbose: print '..write out the albedo values'
    if albedo['create_albedo_file']:
        albedo['albedo_file'] = output_file+'_alb'
        write_albedo_file(albedo['albedo_file'],source['alb_wvl'],source['alb'])
        output.write('albedo_file \t%s\n' % albedo['albedo_file'])
    elif albedo.get('albedo_file'):
        output.write('albedo_file \t%s\n' % albedo['albedo_file'])
    elif albedo.get('sea_surface_albedo'):
        output.write('brdf_cam u10\t%i\n' % albedo['wind_speed'])
    else:
        output.write('albedo\t%f\n' % albedo['albedo'])
    
    if verbose: print '..write out the geo values'
    output.write('zout %s \n' % " ".join([str(x) for x in geo['zout']]))
    if geo.get('lat'):
        output.write("latitude\t%s %f\n" % ('S' if geo['lat']<0 else 'N',abs(geo['lat'])))
        output.write("longitude\t%s %f\n" % ('W' if geo['lon']<0 else 'E',abs(geo['lon'])))
    if geo.get('sza'):
        output.write('sza\t%f\n' % geo['sza'])
    if geo.get('doy'):
        output.write('day_of_year\t%i\n' % geo['doy'])
    if geo.get('year'):
        output.write('time\t%04i\t%02i\t%02i\t%02i\t%02i\t%02i\n' 
                     %(geo['year'],geo['month'],geo['day'],geo['hour'],geo['minute'],geo['second']))
        
    if 'ext' in aero:
        if verbose: print '..write out the aerosol parameters'
        if make_base:
            base.write('aerosol_default\n')
            base.write('disort_intcor moments\n') #set to use moments for explicit aerosol file
        elif not include_base:    
            output.write('aerosol_default\n')
            output.write('disort_intcor moments\n') #set to use moments for explicit aerosol file
        if not aero.get('link_to_mom_file'):
            if not aero.get('file_name'):
                aero['file_name'] = output_file+'_aero'
            write_aerosol_file_explicit(aero['file_name'],aero['z_arr'],aero['ext'],aero['ssa'],aero['asy'],aero['wvl_arr'],verbose=verbose)
        output.write('aerosol_file explicit \t%s\n' % aero['file_name'])
    
    if 'tau' in cloud:
        if verbose: print '..write out the cloud properties'
        if not cloud.get('link_to_mom_file'):
            if not cloud.get('file_name'):
                cloud['file_name'] = output_file+'_cloud'
        if cloud['phase']=='ic':
            if verbose: print '..Ice cloud'
            output.write('ic_file %s \t %s\n' % ('moments' if cloud['write_moments_file'] else '1D',cloud['file_name']))
            output.write('ic_properties baum_v36 interpolate\n')
        elif cloud['phase']=='wc':
            if verbose: print '..Liquid water cloud'
            output.write('wc_file %s \t %s\n' % ('moments' if cloud['write_moments_file'] else '1D',cloud['file_name']))
            output.write('wc_properties mie\n')
        else:
            raise ValueError('phase value in cloud dict not recognised')
        if cloud['write_moments_file']:
            if not cloud.get('link_to_mom_file'):
                write_cloud_file_moments(cloud['file_name'],cloud['tau'],cloud['ref'],cloud['zbot'],cloud['ztop'],
                                         verbose=verbose,moms_dict=cloud.get('moms_dict'))
        else:
            if not cloud.get('link_to_mom_file'):
                write_cloud_file(cloud['file_name'],cloud['tau'],cloud['ref'],cloud['zbot'],cloud['ztop'],verbose=verbose)
    output.close()
    if make_base:
        base.close()
    if verbose: print 'Finished printing main input file: Closing file'   

# <codecell>

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

# <codecell>

def make_pmom_inputs(fp_rtm='C:/Users/sleblan2/Research/4STAR/rtm_dat/',source='solar'):
    """
    Purpose:
    
        Create the moments used as input for the write_input_aac function
        
    Input:
    
        fp_rtm: path to the files to load
        source: either solar or thermal (default is solar)
    
    Dependencies:

        numpy
        scipy
    
    Required files:
   
        mie_hi.out
        wc.sol.long.mie.cdf
        wc.sol.short.mie.cdf        
    
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-30, NASA Ames, from Santa Cruz, CA
        Modified: Samuel LeBlanc, 2015-07-02, Santa Cruz, CA
                    - added source keyword
                    - added new file for loading of thermal mie properties
        Modified: Samuel LeBlanc, 2015-07-07, NASA Ames, CA
                    - added the allpmom netcdf file from Claudia Emde for Solar 
        Modified: Samuel LeBlanc, 2015-07-08, Santa Cruz, CA
                    - fixed bugs with Claudia Emde's full pmom netcdf files. Added longer wavelengths to be saved.
    """
    import numpy as np
    import scipy.io as sio
    
    if source=='solar':
        mie = sio.netcdf_file(fp_rtm+'wc_allpmom.sol.mie.cdf','r')
        mie_long = sio.netcdf_file(fp_rtm+'wc.sol.long.mie.cdf','r')
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
    elif source=='solar_sub':
        mie = sio.idl.readsav(fp_rtm+'mie_hi.out')
        mie_long = sio.netcdf_file(fp_rtm+'wc.sol.long.mie.cdf','r')
        pmom = mie
        pmom['wvl'] = np.append(mie['wvl'],mie_long.variables['wavelen'].data)
        pmom['ntheta'] = np.concatenate((mie['ntheta'],np.swapaxes(mie_long.variables['ntheta'].data[:,:,0],0,1)),axis=1)
        pmom['rho'] = np.concatenate((mie['rho'],np.swapaxes(mie_long.variables['rho'].data,0,1)),axis=1)
        pmom['nmom'] = np.concatenate((mie['nmom'],np.swapaxes(mie_long.variables['nmom'].data[:,:,0],0,1)),axis=1)
        pmom['ssa'] = np.concatenate((mie['ssa'],np.swapaxes(mie_long.variables['ssa'].data,0,1)),axis=1)
        pmom['ext'] = np.concatenate((mie['ext'],np.swapaxes(mie_long.variables['ext'].data,0,1)),axis=1)
        pmom['nim'] = np.append(mie['nim'],mie_long.variables['refim'].data)
        pmom['nre'] = np.append(mie['nre'],mie_long.variables['refre'].data)
        pmom['pmom'] = np.concatenate((mie['pmom'],np.concatenate((np.swapaxes(mie_long.variables['pmom'].data[:,:,0,:],0,1),
                                                                   np.zeros((30,79,750))),axis=2)),axis=1)
        pmom['phase'] = np.concatenate((np.concatenate((mie['phase'],np.zeros((30,754,602))),axis=2),
                                        np.swapaxes(mie_long.variables['phase'].data[:,:,0,:],0,1)),axis=1)
        pmom['theta'] = np.concatenate((np.concatenate((mie['theta'],np.zeros((30,754,602))),axis=2),
                                        np.swapaxes(mie_long.variables['theta'].data[:,:,0,:],0,1)),axis=1).shape 
    elif source=='thermal':
        mie_trm = sio.netcdf_file(fp_rtm+'wc_trm_longmie.cdf','r')
        pmom = {'wvl':mie_trm.variables['wavelen'].data*1000.0, 
                'ref':mie_trm.variables['reff'].data,
                'ntheta':np.swapaxes(mie_trm.variables['ntheta'].data[:,:,0],0,1),
                'rho':np.swapaxes(mie_trm.variables['rho'].data,0,1),
                'nmom':np.swapaxes(mie_trm.variables['nmom'].data[:,:,0],0,1),
                'ssa':np.swapaxes(mie_trm.variables['ssa'].data,0,1),
                'ext':np.swapaxes(mie_trm.variables['ext'].data,0,1),
                'nim':mie_trm.variables['refim'].data,
                'nre':mie_trm.variables['refre'].data,
                'pmom':np.swapaxes(mie_trm.variables['pmom'].data[:,:,0,:],0,1),
                'phase':np.swapaxes(mie_trm.variables['phase'].data[:,:,0,:],0,1),
                'theta':np.swapaxes(mie_trm.variables['theta'].data[:,:,0,:],0,1)}
    else:
        print 'Not a correct option for source: select either solar, solar_sub, or thermal'
        return None
    return pmom

# <codecell>

def build_aac_input(fp,fp_alb,fp_out,fp_pmom=None,fp_uvspec='/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec',fp_output=None,
                    wvl_file_sol=None,wvl_file_thm=None,aero_clear=False):
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
        
    Dependencies:
    
        numpy
        scipy
        load_modis
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
        
    """
    import numpy as np
    import scipy.io as sio
    import Run_libradtran as RL
    import load_modis as lm
    import os
    if fp_pmom:
        pmom_solar = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='solar')
        pmom_thermal = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='thermal')
    else:
        pmom_solar = RL.make_pmom_inputs(source='solar')
        pmom_thermal = RL.make_pmom_inputs(source='thermal')
    geo = {'zout':[0,3,100],'year':2007,'day':15,'minute':0,'second':0}
    aero = {'z_arr':[3.0,4.0]}
    cloud = {'ztop':3.0,'zbot':2.0,'phase':'wc','write_moments_file':True}
    source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/','run_fuliou':True}
    albedo = {'create_albedo_file':False}
    if not fp_output:
        change_fp_output = True
    else:
        change_fp_output = False
    
    
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
            
        try:
            file_list = file(fp_out+'AAC_list_file_%s.sh' % mmm,'w')
        except Exception,e:
            print 'Problem with accessing file, return Exception: ',e
            return
        print 'Starting list file'
        fp_out2 = fp_out+mmm+'/'
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
        alb_geo_sub = np.nanmean(np.nanmean(alb_geo['MCD43GF_CMG'].reshape([48,21600/48,75,43200/75]),3),1)
        alb_geo_lat = np.linspace(90,-90,num=48)
        alb_geo_lon = np.linspace(-180,180,num=75)
        
        print 'Running through the files'
        for ilat,lat in enumerate(input_mmm['MODIS_lat'][0,0]):
            for ilon,lon in enumerate(input_mmm['MODIS_lon'][0,0]):
                geo['lat'],geo['lon'] = lat,lon
                # set the aerosol values
                aero['wvl_arr'] = input_mmm['MOC_wavelengths'][0,0][0,:]*1000.0
                aero['ext'] = np.abs(input_mmm['MOC_ext_mean'][0,0][ilat,ilon,:])
                if aero_clear:
                    aero['ext'] = aero['ext']*0.0
                if np.isnan(aero['ext']).all():
                    print 'skipping lat:%i, lon:%i' % (ilat,ilon)
                    continue
                aero['ssa'] = input_mmm['MOC_ssa_mean'][0,0][ilat,ilon,:]
                aero['asy'] = input_mmm['MOC_asym_mean'][0,0][ilat,ilon,:]
                if aero['wvl_arr'].max()<100000.0:
                    aero['wvl_arr'] = np.append(aero['wvl_arr'],100000.0)
                    aero['ext'] = np.append(aero['ext'],aero['ext'][-1])
                    aero['ssa'] = np.append(aero['ssa'],aero['ssa'][-1])
                    aero['asy'] = np.append(aero['asy'],aero['asy'][-1])
                # set the cloud values
                cloud['tau'] = input_mmm['MODIS_COD_mean'][0,0][ilat,ilon]
                cloud['ref'] = input_mmm['MODIS_effrad_mean'][0,0][ilat,ilon]
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
                    RL.write_input_aac(file_out_sol,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                                       make_base=make_base,fp_base_file=fp_base_file,set_quiet=True)
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
                    RL.write_input_aac(file_out_thm,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                                       make_base=False,fp_base_file=fp_base_file,set_quiet=True)
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

# <codecell>

def read_libradtran(fp,zout=[0,3,100]):
    """
    Purpose:
    
        Program to read the output of libradtran irradiance files
        Very simple output of 3 zout, one wavelength
    
    Inputs:
    
        fp: full path of file to read
        zout: array of zout values
        
    Outputs:
    
        out: dictionary with
            wvl:
            zout:
            direct_down: irradiance down from direct beam
            diffuse_down: irradiance down from diffuse soruces
            diffuse_up: irradiance upwwelling
            int_dir_dn: average intensity direct down
            int_dif_dn: average intensity diffuse down
            int_dif_up: average intensity diffuse up
        
    Dependencies:
    
        os
        numpy
        
    Required files:
    
        file output of libradtran
        
    Example:
    
        ...
        
    Modification History:
    
        Written: Samuel LeBlanc, 2015-07-13, Santa Cruz, CA
        
    """
    #import os
    import numpy as np
    #if not os.path.isfile(fp):
    #    raise IOError('File not found')
    #    return
    dat = np.fromfile(fp,sep=' ').reshape((len(zout),7))
    output = {'wvl':dat[:,0],
              'zout':zout,
              'direct_down':dat[:,1],
              'diffuse_down':dat[:,2],
              'diffuse_up':dat[:,3],
              'int_dir_dn':dat[:,4],
              'int_dif_dn':dat[:,5],
              'int_dif_up':dat[:,6]}
    return output

# <codecell>

def read_aac(fp_out,fp_mat,mmm=None):
    """
    Purpose:
    
        Simple program to read all the output of the libradtran runs for AAC
        Program to read the output of libradtran irradiance files
        Very simple output of 3 zout, one wavelength
    
    Inputs:
    
        fp_out: full path of the directory with the files to read
        fp_mat: full path of mat file with lat and lon to use
        mmm: string with the season defined (can be DJF,MAM,JJA, or SON)
        
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
    
        ...
        
    Modification History:
    
        Written: Samuel LeBlanc, 2015-07-13, Santa Cruz, CA
        
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
                    sol = RL.read_libradtran(file_out_sol,zout=output['zout'])
                    thm = RL.read_libradtran(file_out_thm,zout=output['zout'])
                except IOError:
                    print 'File not found skip: lat%02i_lon%02i_%s_HH%02i' %(ilat,ilon,mmm,iutc)
                    continue
                output['SW_irr_dn_utc'][:,ilat,ilon,iutc] = sol['direct_down']+sol['diffuse_down']
                output['SW_irr_up_utc'][:,ilat,ilon,iutc] = sol['diffuse_up']
                output['LW_irr_dn_utc'][:,ilat,ilon,iutc] = thm['direct_down']+thm['diffuse_down']
                output['LW_irr_up_utc'][:,ilat,ilon,iutc] = thm['diffuse_up']
            print mmm,ilat,ilon
            output['SW_irr_dn_avg'][:,ilat,ilon] = np.mean(output['SW_irr_dn_utc'][:,ilat,ilon,:],axis=1)
            output['SW_irr_up_avg'][:,ilat,ilon] = np.mean(output['SW_irr_up_utc'][:,ilat,ilon,:],axis=1)
            output['LW_irr_dn_avg'][:,ilat,ilon] = np.mean(output['LW_irr_dn_utc'][:,ilat,ilon,:],axis=1)
            output['LW_irr_up_avg'][:,ilat,ilon] = np.mean(output['LW_irr_up_utc'][:,ilat,ilon,:],axis=1)
    return output

# <codecell>

if __name__=='__main__':

# <codecell>

    import numpy as np

# <codecell>

    ext = np.array([[0.5,0.4,0.3],[0,0,0]])
    ssa = ext*1.9
    asy = ext*1.7
    ext.shape

# <codecell>

    z_arr = np.array([3,4])
    wvl_arr = np.array([350.,500.,650.])

# <codecell>

    write_aerosol_file_explicit('C:\Users\sleblan2\libradtran/aero.inp',z_arr,ext,ssa,asy,wvl_arr,verbose=True)

# <codecell>

    write_cloud_file('C:\Users\sleblan2\libradtran\cloud.inp',10,10,2,3,verbose=True)

# <codecell>

    write_albedo_file('C:\Users\sleblan2\libradtran/alb.inp',[500.0,600.0],[0.5,0.6],verbose=True)

# <codecell>

    write_input_aac('C:\Users\sleblan2\libradtran/test_input.inp',geo={'zout':[0,3,100],'wvl_range':[202,5600]},aero={},cloud={'tau':10,'ref':10,'phase':'wc','ztop':3,'zbot':2},source={},
                    verbose=True)

# <codecell>

    import Run_libradtran
    reload(Run_libradtran)

# <codecell>

    mie = make_pmom_inputs()

# <codecell>

    geo = {'zout':[0,3,100],
           'lat':-14.0,
           'lon':-85.0,
           'year':2007,'month':2,'day':10,'hour':10,'minute':0,'second':0}
    aero = {'ext':[  4.85364736e-02,   4.66139195e-02,   4.47312609e-02,
         4.30849589e-02,   4.19923201e-02,   4.03355801e-02,
         3.74764159e-02,   3.45595009e-02,   3.19684762e-02,
         2.94772306e-02,   2.74202103e-02,   2.57334360e-02,
         2.39507641e-02,   2.08731768e-02,   1.67933569e-02,
         1.29016393e-02,   9.04034361e-03,   6.65431703e-03,
         4.35758656e-03,   3.47793084e-03,   2.59552084e-03,
         1.92045503e-03,   1.48977972e-03,   1.14460091e-03,
         7.92407241e-04,   5.10274383e-04,   3.17954425e-04,
         1.69683997e-04,   8.10304392e-05,   3.35441191e-05],
            'ssa':[ 0.94027621,  0.94451989,  0.94726128,  0.94923121,  0.95027106,
        0.95172633,  0.95389696,  0.95572622,  0.9572892 ,  0.95868109,
        0.9598034 ,  0.9607617 ,  0.96177931,  0.96326677,  0.95907986,
        0.94878827,  0.92993408,  0.90942679,  0.87588452,  0.85667554,
        0.83186815,  0.80753487,  0.78811889,  0.76859549,  0.741118  ,
        0.70538566,  0.66048003,  0.58568778,  0.46758827,  0.27607095],
            'asy':[ 0.7852096 ,  0.78103743,  0.77684837,  0.77285771,  0.77073574,
        0.767068  ,  0.76026127,  0.75438874,  0.74906518,  0.74499582,
        0.74245588,  0.7404872 ,  0.73816073,  0.73253394,  0.72113882,
        0.70217822,  0.66814234,  0.63567558,  0.58822783,  0.56218055,
        0.52791341,  0.49239392,  0.4624792 ,  0.43165396,  0.38933818,
        0.34021788,  0.28974332,  0.22683931,  0.1585286 ,  0.08400939],
            'z_arr':[3.0,4.0],
            'wvl_arr':np.array([  0.20005,   0.2343 ,   0.2648 ,   0.2921 ,   0.3105 ,   0.34   ,
          0.3975 ,   0.4675 ,   0.54625,   0.6423 ,   0.742  ,   0.8415 ,
          0.9655 ,   1.226  ,   1.6574 ,   2.2024 ,   3.0044 ,   3.7544 ,
          4.9    ,   5.57   ,   6.51   ,   7.57   ,   8.545  ,   9.645  ,
         11.35   ,  13.7    ,  16.7    ,  21.75   ,  30.35   ,  50.     ])*1000.0}
    cloud = {'tau':7.6,'ref':12.47,'ztop':3.0,'zbot':2.0,
             'phase':'wc','write_moments_file':True,'moms_dict':mie}
    source = {'wvl_range':[202,5600],
              'source':'solar',
              'integrate_values':True}
    albedo = {'create_albedo_file':False,
              'albedo':0.2}

# <codecell>

    write_input_aac('C:\Users\sleblan2\libradtran/test_input_aac.inp',geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=True)

# <codecell>

    fp = 'C:\Users\sleblan2/Research/libradtran/testing_new/AAC_input_lat06_lon19_DJF_HH17_sol.out'

# <codecell>

    import pandas as pd

# <codecell>

    d = pd.read_csv(fp,delim_whitespace=True,engine='c')

# <codecell>

    d

# <codecell>

    dat = np.array(d)

# <codecell>

    dat

# <codecell>

    %timeit rr = np.fromfile(fp,sep=' ').reshape((3,7))

# <codecell>

    %timeit dd = np.array(pd.read_csv(fp,delim_whitespace=True,engine='c'))

# <codecell>

    %timeit gg = np.loadtxt(fp)

# <codecell>

    %timeit gh = np.genfromtxt(fp)

# <codecell>

    Run_libradtran.read_libradtran(fp)

# <codecell>


