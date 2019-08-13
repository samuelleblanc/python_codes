
# coding: utf-8

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
#     - load_utils
#   
# Needed Files:
# 
#   - ...
#   
#   
# Modification History:
# 
#     Wrtten: Samuel LeBlanc, NASA Ames, from Santa Cruz, 2015-06-26
#     Modified: Samuel LeBlanc, NASA Ames, 2015-10-06
#             - added function to write out the write_aac dictionaries to a version file

# In[182]:


def write_cloud_file(output_file,tau,ref,zbot,ztop,verbose=False,append_directly_below=False):
    """
    Purpose:

        Program to write out the cloud file profile
        for either ic_file or wc_file in 1D settings
        outputs 3 columns: 
            # z [km]    IWC/LWC [g/m^3]    Reff [um]
            
        ** translates tau and ref to lwc/iwc with the equation LWP = 2/3*tau*ref and LWC=LWP/z_extent*0.001
        ** Writes out file based on layer properties, not level properties. 
           (given value of ref and tau are for a cloud from zbot to ztop)
    
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
        append_directly_below: (default False) if true then opens a already existing file
                               and only writes out the last line as an append
    
    Dependencies:

        numpy
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-26, NASA Ames, from Santa Cruz, CA
        Modified (v1.1): Samuel LeBlanc, DC8 flying above Korea, 2016-05-02
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
    if not append_directly_below:
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
    else:
        try:
            output = file(output_file,'a')
        except Exception,e:
            print 'Problem with accessing file, return Exception: ',e
            return
        output.write('%4.4f\t%4.5f\t%3.2f\n' % (zbot,lwc,ref))
    output.close() 
    if verbose:
        print('..File finished write_cloud_file, closed')


# In[204]:


def write_aerosol_file_explicit(output_file,z_arr,ext,ssa,asy,wvl_arr,verbose=False,expand_hg=False):
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
        expand_hg: (default False) if true expands the Henyey Greenstein legendre moments from the assymetry parameter
    
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


# In[184]:


def write_aerosol_file_explicit_wvl(output_file,wvl_arr,ext,ssa,asy,verbose=False,expand_hg=False):
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
        expand_hg: (default False) if true, then epands the Henyey-Greenstein phase function, not just printing the asymmetry p
    
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
        if expand_hg:
            asys = [asy[iw]**n for n in range(100)]
            st = '%f\t%f\t%1.6f\t' % (wvl,ext[iw],ssa[iw])
            stt = st+'\t'.join(map(str,asys))+'\n'
            output.write(stt)
        else:
            output.write('%f\t%f\t%1.6f\t%1.6f\t%1.6f\n' % (wvl,ext[iw],ssa[iw],1.0,asy[iw]))
    output.close()
    if verbose:
        print('..File finished write_aerosol_file_explicit_wvl, closed')


# In[45]:


def write_cloud_file_moments(output_file,tau,ref,zbot,ztop,moms_dict=None,
                             verbose=False,append_directly_below=False,wvl_range=None):
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
        wvl_range: sets the range of wavelength for the cloud properties to be saved. If not set, prints all values in moms_dict
    
    Output:

        output_file              
    
    Keywords: 

        verbose: (default False) if true prints out info about file writing 
        append_directly_below: (default False) if true appends a new line to an existing file with cloud properties directly 
                               below the last line of the already existing file
    
    Dependencies:

        numpy
        Run_libradtran (this file)
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-29, NASA Ames, from Santa Cruz, CA
        Modified (v1.1): Samuel LeBlanc, DC8 flying above Korea, 2016-05-02
                        - added the append_directly_below keyword
        Modified (v1.2): Samuel LeBlanc, Santa Cruz, CA, 2017-02-17
                        - added the wvl_range keyword
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
    
    if not append_directly_below:
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
    else:
        try:
            output = file(output_file,'a')
        except Exception,e:
            print 'Problem with accessing file, return Exception: ',e
            return
        if verbose:
            print('..printing to file: %s' % output_file)

        file_cloud = output_file+'_zbelow'
        output.write('%4.4f\t%s\n' % (zbot,file_cloud))
    
    write_cloud_file_moments_wvl(file_cloud,wvl,ext,ssa,moments,nmom,verbose=verbose,wvl_range=wvl_range)
    
    output.close() 
    if verbose:
        print('..File finished write_cloud_file_moments, closed')


# In[ ]:


def write_cloud_file_moments_wvl(output_file,wvl,ext,ssa,moments,nmom,verbose=False,wvl_range=None):
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
        wvl: array of wavelengths in nm referring to the properties
        wvl_range: only print the values within this wavelength range
    
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
        
    if not wvl_range:
        wvl_range = [min(wvl),max(wvl)]
    if verbose: print '... wvl_range:',wvl_range[0],wvl_range[1]
        
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose:
        print('..printing to cloud moments properties wavelength defined file: %s' % output_file)
    output.write('# wvl[nm]    ext[km^-1]   ssa[unitless]  legendre_moments\n')
    for iw,wv in enumerate(wvl):
        if (wv<wvl_range[0]-10.0 or wv>wvl_range[1]+10.0) and wvl_range[1]<10000.0:
            continue
        try:
            if len(moments.shape)>1:
                output.write('%f\t%f\t%1.6f\t%s \n' % 
                             (wv,ext[iw],ssa[iw]," ".join([str(x/(2.0*i+1.0)) for i,x in enumerate(moments[iw,0:int(nmom[iw])])])))
            else:
                st = '%f\t%f\t%1.6f\t%s \n' %                             (wv,ext[iw],ssa[iw],
                              " ".join(['{:1.5g}'.format(x/(2.0*i+1.0)) for i,x in enumerate(moments[iw][0,0:int(nmom[iw])])]))
                if len(st)>1000000:
                    st = '%f\t%f\t%1.6f\t%s \n' %                             (wv,ext[iw],ssa[iw],
                              " ".join(['{:1.4g}'.format(x/(2.0*i+1.0)) for i,x in enumerate(moments[iw][0,0:int(nmom[iw])])]))
                    an = len(moments[iw][0,0:int(nmom[iw])])-1
                    print 'trying to reduce size', len(st), an
                    while len(st)>1000000:
                        df = int((len(st)-1000000)/10)
                        an = an-df-10
                        st = '%f\t%f\t%1.6f\t%s \n' %                             (wv,ext[iw],ssa[iw],
                              " ".join(['{:1.4g}'.format(x/(2.0*i+1.0)) for i,x in enumerate(moments[iw][0,0:an])]))
                        print 'in size reduction', len(st),an
                output.write(st)
                #output.write('%f\t%f\t%1.6f\t%s \n' % 
                #             (wv,ext[iw],ssa[iw],
                #" ".join(['{:1.4g}'.format(x/(2.0*i+1.0)) for i,x in enumerate(moments[iw][0,0:nmom[iw]])])))
        except:
            import pdb
            pdb.set_trace()
    output.close()
    if verbose:
        print('..File write_cloud_file_moments_wvl finished, closed')


# In[185]:


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


# In[ ]:


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
        Modified: Samuel LeBlanc, 2017-02-17, Santa Cruz, CA
                  Changed the output commands to easily pass the ice pmoms as well as the wc
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
    nm = moms_dict.get('max_nmoms',-1)
    if moms_dict.get('max_nmoms',False):
        nmom = moms_dict['nmom'][ir,:]*0+nm
    else:
        nmom = moms_dict['nmom'][ir,:]

    return ext,moms_dict['ssa'][ir,:],wvl,moms_dict['pmom'][ir,:],nmom


# In[1]:


def write_input_aac(output_file,geo={},aero={},cloud={},source={},albedo={},
                    verbose=False,make_base=False,fp_base_file=None,set_quiet=True,max_nmom=None,solver='disort'):
    """
    Name:

        write_input_aac
    
    Purpose:

        Writes the libradtran input file with defined defaults, see below
        outputs libradtran input files in ascii format
        
        ** There is a multitude of other default parameters. 
        ** These should be changed for any other type of input files to be written
    
    Input: 
  
        output_file: full path to file to be written
        geo: dictionary with geometry details
            sza: solar zenith angle (libradtran default is 0)
            lat: latitude
            lon: longitude
            doy: day of year  # for calculating sun earth distance
            year: year (YYYY format) # used if you want uvspec to calculate the sza
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
                       By default it is created in this program if it is not set, 
                       if link_to_mom_file is set, this must be defined
            disort_phase: (default False) if true then writes then phase function instead of the moments.
            expand_hg: (default False) if true expands the henyey greenstein phase function into legendre moments.
        cloud: dictionary with cloud properties
            tau: value of cloud optical thickness
            ref: value of effective radius in microns
            ztop: location in km of cloud top
            zbot: location in km of cloud bottom
            phase: either 'ic' for ice cloud or 'wc' for water cloud
            write_moments_file: (default False) if True, writes out moments into an ascii file instead 
                                of reading directly from the netcdf file. 
                                Requires the moments to be put in the cloud dict and the nmom.
            moms_dict: is the moment dict structure returned from the mie_hi.out idl.readsav. or now the mie_hi_delta.mat
                       To be used when building the moments file
            link_to_mom_file: if True then no moments file is written out, but referenced via file_name saved to cloud dict
            file_name: file name and paths of moments file for cloud defined. 
                       By default it is created in this program if not set, if link_to_mom_file is set, this must be defined
            cloud_below: (default False) if set, creates a liquid cloud with constant properties (tau=10,ref=10 microns), 
                         of thickness 1 km, directly below the zbot value. Used for in cloud calculations
        source: dictionary with source properties
            wvl_range: range of wavelengths to model (default [202,500])
            source: can either be thermal or solar
            dat_path: data path to be used. Defaults to pleaides values (/u/sleblan2/libradtran/libRadtran-2.0-beta/data/)
            integrate_values: if True (default), then the resulting output parameters are integrated over wavelength range
                              if set to False, returns per_nm irradiance values
            wvl_filename: filename and path of wavelength file (second column has wavelengh in nm to be used)
            run_fuliou: if set to True, then runs fu liou instead of sbdart (default is False)
            slit_file: (defaults to None) Full file path of slit file to use for calculating the radiative transfer.
            atm_file: (defaults to None) Full file path of the atmosphere file to use, 
                      if not set, libradtran calculates from position and time on Earth
            zenith: (defaults to False) if True, prepares the input file for returning zenith radiances
                    adds a zenith viewing angle (umu=-1), azimuth viewing angle (phi=130) and solar azimuth angle (phi0=130)
        albedo: dictionary with albedo properties
            create_albedo_file: (default False) if true then albedo file is created with properties defined by alb_wvl and alb 
            albedo_file: path of albedo file to use if already created 
            albedo: value of albedo. Only used if create_albedo_file is false and albedo_file is empty 
                    (defaults to 0.29 - Earth's average)
            alb: wavelength dependent value of albedo for use when create_albedo_file is set to True
            alb_wvl: wavelength grid of albedo for use when create_albedo_file is set to True
            sea_surface_albedo: (default False) If True, sets the sea surface to be parameterized by cox and munk, 
                                requires wind_speed to be set.
            wind_speed: [m/s] (default to 10 m/s) Only used if sea_surface_albedo is set to True. 
                        wind speed over water for parameterization of cox_and_munk
        make_base: boolean to set if the base file is to be written out
                   if False, no base file is saved.
                   Base file contains quiet flag, data path, source, rad transfer solver, 
                     output process, wavelength source file, 
        fp_base_file: full file path for base file. 
                      If set to a file path and make_base to False, then include path is printed to input file. 
        set_quiet: if True then quiet is set in the input file, if False, quiet is not set. (default True)
        max_nmom: set the maximum number of moments to write out (default None)
    
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
        
        Example taken from lut creation for NAAMES, 
            uses St John's Newfoundland location, on Nov. 19th for low clouds (from 0.5km to 1.5km) over ocean
            returns zenith radiances and irradiance at 0.2,3.0 and 100 km
            uses internal cloud properties
            for vis spectral range, with slit function
            no aerosol
            cloud tau=2.0, ref=5.0, sza=40.0
        
        >>> geo = {'lat':47.6212167,'lon':52.74245,'doy':322,'zout':[0.2,3.0,100.0]}
        >>> aero = {} # none
        >>> cloud = {'ztop':1.5,'zbot':0.5,'write_moments_file':False}
        >>> source = {'wvl_range':[400.0,981.0],'source':'solar','integrate_values':False,'run_fuliou':False,
                      'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/'}
        >>> albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':14.0}
        >>> cloud['phase'] = 'wc'
        >>> geo['sza'] = 40.0
        >>> cloud['tau'] = 2.0
        >>> cloud['ref'] = 5.0
        
        >>> RL.write_input_aac('/u/sleblan2/NAAMES/runs/NAAMES_v1.dat',geo=geo,aero=aero,cloud=cloud,
                               source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
        
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
        Modified: Samuel LeBlanc, 2015-10-06, NASA Ames, CA
                - modified comments
                - added slit_file to source options
                - added atm_file to source options
                - added zenith radiance keyword
        Modified: Samuel LeBlanc, on the DC8 flying above Korea, 2016-05-02
                - added the cloud_below keyword to the cloud options.
        Modified: Samuel LeBlanc, NASA Ames, CA, 2016-09-22
                - added expansion of henyey greenstein legendre moments for aerosol explict files
        Modified: Samuel LeBlanc, Santa Cruz, CA, 2017-12-06
                - modified use of mie_hi.out for mie_hi_delta.mat with new delta-scaling legendre functions
    """
    import numpy as np
    from Run_libradtran import write_aerosol_file_explicit,write_cloud_file,write_albedo_file,merge_dicts
    from Run_libradtran import write_cloud_file_moments
    
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose: print 'Opening input file %s' % output_file
    
    if make_base:
        if fp_base_file:
            base = file(fp_base_file,'w')
            include_base = True
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
                          'run_fuliou':False,
                          'zenith':False},source)
    albedo = merge_dicts({'create_albedo_file':False,'sea_surface_albedo':False,'wind_speed':10,
                          'albedo':0.29},albedo)
    geo = merge_dicts({'zout':[0,100]},geo)
    cloud = merge_dicts({'write_moments_file':False,'cloud_below':False},cloud)
    
    if source.get('source')=='solar':
        source['source'] = 'solar '+source['dat_path']+'solar_flux/kurudz_1.0nm.dat per_nm'
    
    if verbose: print '..write out general default values'
    if make_base:
        if source['run_fuliou']:
            base.write('mol_abs_param fu\n')
        else:
            base.write('mol_abs_param sbdart\n')
        base.write('rte_solver {}\n'.format(solver))
        if set_quiet:
            base.write('quiet\n')
    else:
        if set_quiet:
            output.write('quiet\n')
        if source['run_fuliou']:
            output.write('mol_abs_param fu\n')
        else:
            output.write('mol_abs_param sbdart\n')
        output.write('rte_solver {}\n'.format(solver))
        
    if include_base:
        output.write('include %s\n' % fp_base_file)
    
    if verbose: print '..write out source dict values'
    if make_base:
        if source['integrate_values']:
            if source['run_fuliou']:
                base.write('output_process sum\n')
            else:
                base.write('output_process integrate\n')
        else:
            base.write('output_process per_nm\n')
        base.write('data_files_path %s\n' % source['dat_path'])
    elif not include_base:
        if source['integrate_values']:
            if source['run_fuliou']:
                output.write('output_process sum\n')
            else:
                output.write('output_process integrate\n')
        else:
            output.write('output_process per_nm\n')
        output.write('data_files_path %s\n' % source['dat_path'])
    output.write('source %s \n' % source['source'])
    
    if source.get('wvl_filename'):
        output.write('wavelength %s\n' % source['wvl_filename'])
    else:
        if source['wvl_range'][0]>source['wvl_range'][1]:
            print 'wvl_range was set inverse, inversing'
            source['wvl_range'] = list(reversed(source['wvl_range'])) 
        if source['wvl_range'][0]<250:
            print 'wvl_range starting too low, setting to 250 nm'
            source['wvl_range'][0] = 250.0
        output.write('wavelength %f %f\n' % (source['wvl_range'][0],source['wvl_range'][1]))
        
    if source.get('slit_file'):
        output.write('slit_function_file %s\n'%source['slit_file'])
        
    if source.get('atm_file'):
        output.write('atmosphere_file %s\n'%source['atm_file'])
        
    if source['zenith']:
        output.write('umu -1.0\n')
        output.write('phi 130.0\n')
        output.write('phi0 130.0\n')
    
    
    
    if verbose: print '..write out the albedo values'
    if albedo['create_albedo_file']:
        albedo['albedo_file'] = output_file+'_alb'
        write_albedo_file(albedo['albedo_file'],albedo['alb_wvl'],albedo['alb'])
        output.write('albedo_file %s\n' % albedo['albedo_file'])
    elif albedo.get('albedo_file'):
        output.write('albedo_file %s\n' % albedo['albedo_file'])
    elif albedo.get('sea_surface_albedo'):
        output.write('brdf_cam u10 %i\n' % albedo['wind_speed'])
    else:
        output.write('albedo %f\n' % albedo['albedo'])
    
    if verbose: print '..write out the geo values'
    output.write('zout %s \n' % " ".join([str(x) for x in geo['zout']]))
    if geo.get('lat'):
        output.write("latitude %s %f\n" % ('S' if geo['lat']<0 else 'N',abs(geo['lat'])))
        output.write("longitude %s %f\n" % ('W' if geo['lon']<0 else 'E',abs(geo['lon'])))
    if geo.get('sza'):
        output.write('sza %f\n' % geo['sza'])
    if geo.get('doy'):
        output.write('day_of_year %i\n' % geo['doy'])
    if geo.get('year'):
        output.write('time %04i %02i %02i %02i %02i %02i\n' 
                     %(geo['year'],geo['month'],geo['day'],geo['hour'],geo['minute'],geo['second']))
        
    if 'ext' in aero:
        if verbose: print '..write out the aerosol parameters'
        if aero.get('disort_phase'):
            dd = 'phase'
        else:
            dd = 'moments'
        aero['expand_hg'] = aero.get('expand_hg',False)
        if make_base:
            base.write('aerosol_default\n')
            base.write('disort_intcor {}\n'.format(dd)) #set to use moments for explicit aerosol file
        elif not include_base:    
            output.write('aerosol_default\n')
            output.write('disort_intcor {}\n'.format(dd)) #set to use moments for explicit aerosol file
        if not aero.get('link_to_mom_file'):
            if not aero.get('file_name'):
                aero['file_name'] = output_file+'_aero'
            write_aerosol_file_explicit(aero['file_name'],aero['z_arr'],aero['ext'],aero['ssa'],
                                        aero['asy'],aero['wvl_arr'],verbose=verbose,expand_hg=aero['expand_hg'])
        output.write('aerosol_file explicit %s\n' % aero['file_name'])
    
    if 'tau' in cloud:
        if verbose: print '..write out the cloud properties'
        if not cloud.get('link_to_mom_file'):
            if not cloud.get('file_name'):
                cloud['file_name'] = output_file+'_cloud'
        if cloud['phase']=='ic':
            if verbose: print '..Ice cloud'
            output.write('ic_file %s %s\n' % ('moments' if cloud['write_moments_file'] else '1D',cloud['file_name']))
            output.write('ic_properties baum_v36 interpolate\n')
        elif cloud['phase']=='wc':
            if verbose: print '..Liquid water cloud'
            output.write('wc_file %s %s\n' % ('moments' if cloud['write_moments_file'] else '1D',cloud['file_name']))
            output.write('wc_properties mie %s\n' %('interpolate' if not source['run_fuliou'] else ' '))
        else:
            raise ValueError('phase value in cloud dict not recognised')
        if cloud['write_moments_file']:
            if not 'ext' in aero:
                output.write('disort_intcor moments\n') 
            if not cloud.get('link_to_mom_file'):
                write_cloud_file_moments(cloud['file_name'],cloud['tau'],cloud['ref'],cloud['zbot'],cloud['ztop'],
                                         verbose=verbose,moms_dict=cloud.get('moms_dict'),wvl_range=source['wvl_range'])
            if cloud['cloud_below']:
                if not cloud.get('link_to_mom_file'):
                    write_cloud_file_moments(cloud['file_name'],10,10,cloud['zbot']-1.0,cloud['zbot'],
                                             verbose=verbose,moms_dict=cloud.get('moms_dict'),
                                             append_directly_below=True,wvl_range=source['wvl_range'])
                
        else:
            if not cloud.get('link_to_mom_file'):
                write_cloud_file(cloud['file_name'],cloud['tau'],cloud['ref'],cloud['zbot'],cloud['ztop'],verbose=verbose)
            if cloud['cloud_below']:
                if not cloud.get('link_to_mom_file'):
                    write_cloud_file(cloud['file_name'],10,10,cloud['zbot']-1.0,cloud['zbot'],
                                     verbose=verbose,append_directly_below=True)
    output.close()
    if make_base:
        base.close()
    if verbose: print 'Finished printing main input file: Closing file'   


# In[178]:


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


# In[1]:


def make_pmom_inputs(fp_rtm='C:/Users/sleblan2/Research/4STAR/rtm_dat/',source='solar',cloudtype='wc',deltascale=True):
    """
    Purpose:
    
        Create the moments used as input for the write_input_aac function
        
    Input:
    
        fp_rtm: path to the files to load
        source: either solar or thermal (default is solar)
        cloudtype: either 'wc' for water cloud or 'ic' for ice cloud 
        deltascale: (default True) if true, uses the mie_hi_delta.mat instead of mie_hi.out
    
    Dependencies:

        numpy
        scipy
    
    Required files:
   
        mie_hi.out or mie_hi_delta.mat
        wc.sol.long.mie.cdf
        wc.sol.short.mie.cdf        
        ic.pmom.ghm.baum.mat
    
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-06-30, NASA Ames, from Santa Cruz, CA
        Modified: Samuel LeBlanc, 2015-07-02, Santa Cruz, CA
                    - added source keyword
                    - added new file for loading of thermal mie properties
        Modified: Samuel LeBlanc, 2015-07-07, NASA Ames, CA
                    - added the allpmom netcdf file from Claudia Emde for Solar 
        Modified: Samuel LeBlanc, 2015-07-08, Santa Cruz, CA
                    - fixed bugs with Claudia Emde's full pmom netcdf files. Added longer wavelengths to be saved.
        Modified: Samuel LeBlanc, 2017-02-15, Santa Cruz, CA
                    - added ic file from calculations of pmom from Create_ic_phase.py file
        Modified: Samuel LeBlanc, 2017-12-06, Santa Cruz, CA
                    - added new mie_hi_delta.mat function, for using delta scaled legendre polynomials
    """
    import numpy as np
    import scipy.io as sio
    
    if cloudtype=='wc':
        if source=='solar':
            mie_long = sio.netcdf_file(fp_rtm+'wc.sol.long.mie.cdf','r')
            if deltascale:
                mie = sio.loadmat(fp_rtm+'mie_hi_delta.mat')
                rho = np.array([1.0])
                pmom = {'wvl':np.append(mie['wvl'],mie_long.variables['wavelen'].data[1:]),
                        'ref':mie['ref'],
                        'ntheta':np.concatenate((mie['ntheta'],np.swapaxes(mie_long.variables['ntheta'].data[1:,:,0],0,1)),axis=1),
                        'rho':mie['rho'],
                        'nmom':np.concatenate((mie['nmom_delta'],np.swapaxes(mie_long.variables['nmom'].data[1:,:,0],0,1)),axis=1),
                        'ssa':np.concatenate((mie['ssa'],np.swapaxes(mie_long.variables['ssa'].data[1:,:],0,1)),axis=1),
                        'ext':np.concatenate((mie['ext'],np.swapaxes(mie_long.variables['ext'].data[1:,:],0,1)),axis=1),
                        'nim':np.append(mie['nim'],mie_long.variables['refim'].data[1:]),
                        'nre':np.append(mie['nre'],mie_long.variables['refre'].data[1:]),
                        'pmom':np.concatenate((np.concatenate((mie['pmom_delta'],np.zeros((30,754,351))),axis=2),
                                               np.swapaxes(mie_long.variables['pmom'].data[1:,:,0,:652],0,1)),axis=1),
                        'phase':np.concatenate((mie['phase'],np.swapaxes(mie_long.variables['phase'].data[1:,:,0,:398],0,1)),axis=1),
                        'theta':np.concatenate((mie['theta'],np.swapaxes(mie_long.variables['theta'].data[1:,:,0,:398],0,1)),axis=1)}
                pmom['file_name'] = [fp_rtm+'mie_hi_delta.mat',fp_rtm+'wc.sol.long.mie.cdf']
            else:
                mie = sio.netcdf_file(fp_rtm+'wc_allpmom.sol.mie.cdf','r')
                try:
                    rho = np.swapaxes(mie.variables['rho'].data,0,1)
                except ValueError:
                    rho = mie.variables['rho'].data
                mie_short = {'wvl':mie.variables['wavelen'].data,
                             'ref':mie.variables['reff'].data,
                             'ntheta':np.swapaxes(mie.variables['ntheta'].data[:,:,0],0,1),
                             'rho':rho,
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
                pmom['file_name'] = [fp_rtm+'wc_allpmom.sol.mie.cdf',fp_rtm+'wc.sol.long.mie.cdf']
        elif source=='solar_sub':
            if deltascale:
                mie = sio.loadmat(fp_rtm+'mie_hi_delta.mat')
                mie['pmom'] = mie['pmom_delta']
                mie['nmom'] = mie['nmom_delta']
            else:
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
            pmom['file_name'] = [fp_rtm+'mie_hi.out',fp_rtm+'wc.sol.long.mie.cdf']
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
            pmom['file_name'] = [fp_rtm+'wc_trm_longmie.cdf']
        else:
            print 'Not a correct option for source: select either solar, solar_sub, or thermal'
            return None
    elif cloudtype =='ic':
        ic = sio.loadmat(fp_rtm+'ic.pmom.ghm.baum.mat')
        pmom = {'wvl':ic['pmom_wvl'][0,:],
                'ref':ic['ref'][0,:],
                'rho':ic['rho'][0,:],
                'ssa':np.swapaxes(ic['ssa'][ic['pmom_iwvl'][0,:],:],0,1),
                'ext':np.swapaxes(ic['ext'][ic['pmom_iwvl'][0,:],:],0,1),
                'pmom':ic['pmom']}#,
                #'phase':ic['phase'],
                #'theta':ic['theta']}
        pmom['nmom'] = np.zeros_like(pmom['pmom'])-1
        pmom['file_name'] = [fp_rtm+'ic.pmom.ghm.baum.mat']
    else:
        print 'Not a correct cloudtype value: either wc or ic'
        return None
    return pmom


# In[3]:


def build_aac_input(fp,fp_alb,fp_out,fp_pmom=None,fp_uvspec='/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec',fp_output=None,
                    wvl_file_sol=None,wvl_file_thm=None,aero_clear=False,version='v1',stdfac_dict={},max_nmom=None,list_only=False,
                    start_lat=None,start_lon=None,mmmlist=['DJF','MAM','JJA','SON'],deltascale=True,zaero=[3.0,4.0],zcld=[2.0,3.0]):
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
        deltascale: (default to True) use the new delta scaled cloud moments if set to True
        
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
        Modofied: Samuel LeBlanc, 2019-02-14, NASA Ames, CA
                - added zaero and zcld keywords in the function call to modify the aerosol and cloud heights.
        
    """
    import numpy as np
    import scipy.io as sio
    import Run_libradtran as RL
    import load_utils as lm
    import os
    if fp_pmom:
        pmom_solar = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='solar',deltascale=deltascale)
        pmom_thermal = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='thermal')
    else:
        pmom_solar = RL.make_pmom_inputs(source='solar',deltascale=deltascale)
        pmom_thermal = RL.make_pmom_inputs(source='thermal')
    if max_nmom:
        pmom_solar['max_nmoms'] = max_nmom
        pmom_thermal['max_nmoms'] = max_nmom
    else:
        pmom_solar['max_nmoms'],pmom_thermal['max_nmoms'] = False,False
        
    geo = {'zout':[0,3,100],'year':2007,'day':15,'minute':0,'second':0}
    aero = {'z_arr':zaero}
    cloud = {'ztop':zcld[1],'zbot':zcld[0],'phase':'wc','write_moments_file':True}
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
        try:
            alb_geo,alb_geo_dict = lm.load_hdf_sd(fpa)
            print 'done loading albedo files'
            alb_geo_sub = np.nanmean(np.nanmean(alb_geo['MCD43GF_CMG'].reshape([48,21600/48,75,43200/75]),3),1)
        except:
            alb_geo = lm.load_hdf_raster1(fpa)
            print 'done loading albedo files using raster 1'
            alb_geo_sub = np.nanmean(np.nanmean(alb_geo.reshape([48,21600/48,75,43200/75]),3),1)        
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
                    #print 'skipping lat:%i, lon:%i' % (ilat,ilon)
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


# In[43]:


def read_libradtran(fp,zout=[0,3,100],num_rad=0):
    """
    Purpose:
    
        Program to read the output of libradtran irradiance and radiance files
        output irradiance for any number of wavelengths, and number of zouts
        in addition to any number of radiance zenith directions, and only one azimuth direction
    
    Inputs:
    
        fp: full path of file to read
        zout: array of zout values
        num_rad: to load radiance set to number of zenith angle directions (defaults to none)
        
    Outputs:
    
        out: dictionary with
            wvl: wavelength array 
            zout: zout values
            direct_down: irradiance down from direct beam
            diffuse_down: irradiance down from diffuse soruces
            diffuse_up: irradiance upwwelling
            int_dir_dn: average intensity direct down
            int_dif_dn: average intensity diffuse down
            int_dif_up: average intensity diffuse up
            umu: (if radiance set to more than one) returns the cosine of the zenith angle direction of the radiance values
            phi: (if rad set) returns a singular azimuth angle of the radiance value
            rad: radiance array at zouts, at a single phi, and at the wavelengths (uu in libradtran notation)
            rad_avg : azimuthally averaged radiance values, not intensity corrected (u0u in libradtran notation)
        
    Dependencies:
    
        os
        numpy
        
    Required files:
    
        file output of libradtran
        
    Example:
    
        ...
        
    Modification History:
    
        Written: Samuel LeBlanc, 2015-07-13, Santa Cruz, CA
        Modified: Samuel LeBlanc, 2015-10-14, NASA Ames, Santa Cruz, CA
                  - added handling of radiance values
        
    """
    #import os
    import numpy as np
    #if not os.path.isfile(fp):
    #    raise IOError('File not found')
    #    return
    zout_len = len(zout)
    if num_rad:
        arr_len = 11
        if num_rad>1:
            arr_len += 3*(num_rad-1)
    else:
        arr_len = 7
    dat = np.fromfile(fp,sep=' ')
    if len(dat)==0 :
        raise IOError('File {} empty'.format(fp))
    try:
        dat = dat.reshape(len(dat)/arr_len/zout_len,zout_len,arr_len)
    except ValueError:
        arr_len = 5
        dat = dat.reshape(len(dat)/arr_len/zout_len,zout_len,arr_len)
    try:
        output = {'wvl':dat[:,0,0].squeeze(),
                  'zout':zout,
                  'direct_down':dat[:,:,1].squeeze(),
                  'diffuse_down':dat[:,:,2].squeeze(),
                  'diffuse_up':dat[:,:,3].squeeze(),
                  'int_dir_dn':dat[:,:,4].squeeze(),
                  'int_dif_dn':dat[:,:,5].squeeze(),
                  'int_dif_up':dat[:,:,6].squeeze()}
    except:
        output = {'wvl':dat[:,0,0].squeeze(),
                  'zout':zout,
                  'direct_down':dat[:,:,1].squeeze(),
                  'diffuse_down':dat[:,:,2].squeeze(),
                  'diffuse_up':dat[:,:,3].squeeze()}
        try:
             output['int_tot'] = dat[:,:,4].squeeze()
        except:
            output['int_tot'] = np.NaN
    if num_rad:
        output['rad'] = dat[:,:,10].squeeze()
        output['rad_avg'] = dat[:,:,9].squeeze()
        output['phi'] = dat[0,0,7]
        output['umu'] = np.array([dat[0,0,8]])
        if num_rad>1:
            output['rad'] = output['rad'][...,np.newaxis]
            output['rad_avg'] = output['rad_avg'][...,np.newaxis]
        for i in range(num_rad-1):
            output['umu'] = np.append(output['umu'],dat[0,0,11+i*3])
            output['rad'][:,:,:,i+1] = dat[:,:,13+i*3]
            output['rad_avg'][:,:,:,i+1] = dat[:,:,12+i*3]
    return output


# In[47]:


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
                    #print 'File not found skip: lat%02i_lon%02i_%s_HH%02i' %(ilat,ilon,mmm,iutc)
                    if iutc==0:
                        print file_out_sol
                    continue
                except ValueError:
                    #print 'Problem with file: lat%02i_lon%02i_%s_HH%02i' %(ilat,ilon,mmm,iutc)
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


# In[ ]:


def combine_wvl(dat1,dat2):
    """
    Program to combine the output of read_libradtran along the wavelength dimension
    """
    import numpy as np
    dat = dat1
    try: 
        if not (dat1['phi']==dat2['phi']):
            print '*** Possible problem, the two arrays do not match in phi'
    except: 
        pass
    if not (dat1['zout']==dat2['zout']):
        print '*** Possible problem, the two arrays do not match in zout'
    try:
        if not (dat1['umu']==dat2['umu']):
            print '*** Possible problem, the two arrays do not match in umu'
    except:
        pass
    try:
        dat['rad_avg'] = np.vstack((dat1['rad_avg'],dat2['rad_avg']))
    except:
        pass
    dat['direct_down'] = np.vstack((dat1['direct_down'],dat2['direct_down']))
    dat['diffuse_up'] = np.vstack((dat1['diffuse_up'],dat2['diffuse_up']))
    dat['diffuse_down'] = np.vstack((dat1['diffuse_down'],dat2['diffuse_down']))
    dat['int_dif_up'] = np.vstack((dat1['int_dif_up'],dat2['int_dif_up']))
    dat['int_dir_dn'] = np.vstack((dat1['int_dir_dn'],dat2['int_dir_dn']))
    dat['int_dif_dn'] = np.vstack((dat1['int_dif_dn'],dat2['int_dif_dn']))
    try:
        dat['rad'] = np.vstack((dat1['rad'],dat2['rad']))
    except:
        pass
    dat['wvl'] = np.hstack((dat1['wvl'],dat2['wvl']))
    return dat


# In[ ]:


def read_lut(fp_out,zout=None,tau=[None],ref=[None],sza=[None],
             phase=['wc','ic'],fmt='lut_sza{sza:02d}_tau{tau:06.2f}_ref{ref:04.1f}_{phase}_w{iwvl:1d}.dat',
             split_wvl=True,numrad=1):
    """
    Purpose:
    
        Simple program to read all the output of the libradtran runs for LUT creation
        Program to read the output of libradtran irradiance and radiance files
    
    Inputs:
    
        fp_out: full path of the directory with the files to read
        zout: array of altitude out values (km)
        tau: array of tau values
        ref: array of ref values
        sza: array of solar zenith angles to read (if not set, will not run through the solar zenith angles)
        phase: array of string values linked to phase
        fmt: format line used to create the files (with new >2.7 format string type, including '{}')
             should include full file name with extension, 
             and named format characters for 'tau', 'ref', 'sza','phase','iwvl' 
             If file name does not require any of these, program doe snot freeze by omission 
        split_wvl: (defaults to True) if set then files are split per wavelength 
                   loads two seperate file per combination of tau, ref, sza, and phase
        numrad: (number of radiance values to read, defaults to 1)
        
    Outputs:
    
        out: dictionary with the saved output 
        
    Dependencies:
    
        Run_libradtran (this file)
        os
        numpy
        
    Required files:
    
        files of libradtran output
        
    Example:
    
        ...
        
    Modification History:
    
        Written: Samuel LeBlanc, 2015-10-14, NASA Ames, Santa Cruz, CA
        Modified: Samuel LeBlanc, 2017-02-17, Santa Cruz, CA
                  Added error checking for testing last file and ice as well as first file
        
    """
    import Run_libradtran as RL
    import os 
    import numpy as np
    
    #test of single load of file
    try:
        dat1 = RL.read_libradtran(os.path.join(fp_out,fmt.format(ref=ref[0],tau=tau[0],sza=sza[0],phase=phase[0],iwvl=0)),
                                  zout=zout,num_rad=numrad)
    except:
        print 'trying backup first run on file:'+fp_out
        dat1 = RL.read_libradtran(os.path.join(fp_out,fmt.format(ref=ref[-1],tau=tau[-1],sza=sza[-1],phase=phase[-1],iwvl=0)),
                                  zout=zout,num_rad=numrad)
    if split_wvl:
        try:
            dat2 = RL.read_libradtran(os.path.join(fp_out,fmt.format(ref=ref[0],tau=tau[0],sza=sza[0],phase=phase[0],iwvl=1)),
                                      zout=zout,num_rad=numrad)
        except:
            dat2 = RL.read_libradtran(os.path.join(fp_out,fmt.format(ref=ref[-1],tau=tau[-1],sza=sza[-1],phase=phase[-1],iwvl=1)),
                                      zout=zout,num_rad=numrad)
        dat = RL.combine_wvl(dat1,dat2)
    else:
        dat = dat1
    output = {'rad':np.zeros((len(phase),len(dat['wvl']),
                              len(zout),len(ref),len(tau),len(sza))),
              'wvl':dat['wvl'],
              'zout':zout,
              'irr_dn':np.zeros((len(phase),len(dat['wvl']),
                                 len(zout),len(ref),len(tau),len(sza))),
              'irr_up':np.zeros((len(phase),len(dat['wvl']),
                                 len(zout),len(ref),len(tau),len(sza))),
              'irr_dn_diff':np.zeros((len(phase),len(dat['wvl']),
                                      len(zout),len(ref),len(tau),len(sza))),
              'tau':tau,
              'ref':ref,
              'sza':sza,
              'phase':phase}
    for iz,s in enumerate(sza):
        for it,t in enumerate(tau):
            for ir,r in enumerate(ref):
                for ip,p in enumerate(phase):
                    try:
                        dat1 = RL.read_libradtran(os.path.join(fp_out,fmt.format(ref=r,tau=t,sza=s,phase=p,iwvl=0)),
                                                  zout=zout,num_rad=numrad)
                        if split_wvl:
                            dat2 = RL.read_libradtran(os.path.join(fp_out,fmt.format(ref=r,tau=t,sza=s,phase=p,iwvl=1)),
                                                      zout=zout,num_rad=numrad)            
                            dat1 = RL.combine_wvl(dat1,dat2)
                    except IOError:
                        continue
                    try:
                        output['rad'][ip,:,:,ir,it,iz] = dat1['rad']
                    except:
                        pass
                    output['irr_up'][ip,:,:,ir,it,iz] = dat1['diffuse_up']
                    output['irr_dn'][ip,:,:,ir,it,iz] = dat1['direct_down']+dat1['diffuse_down']
                    output['irr_dn_diff'][ip,:,:,ir,it,iz] = dat1['diffuse_down']
                    print s,t,r
    return output               
#'lut_sza%02i_tau%06.2f_ref%04.1f'
#'lut_sza%02i_ref%02.1f_tau%03.1f'


# In[4]:


def print_version_details(filename,vv,geo={},aero={},cloud={},source={},albedo={},
                          tau=[None],ref=[None],sza=[None],cloud_pmom_file=None,
                          fmt='lut_sza{sza:02.0f}_tau{tau:06.2f}_ref{ref:04.1f}_{phase}_w{iwvl:1d}.dat',use_json=True):
    'Program to write an ascii file to print out the version and set up info'
    if use_json:
        from load_utils import save_to_json
        lut = {'LUT_input_version':vv,'sza':sza,'tau':tau,'ref':ref,'format':fmt}
        d = {'lut_details':lut,'geo':geo,'aero':aero,'source':source,'albedo':albedo}
        if cloud_pmom_file:
            from copy import deepcopy
            d['cloud'] = deepcopy(cloud)
            d['cloud']['pmom'] = None
            try:
                d['cloud']['moms_dict'] = None
            except:
                pass
            d['cloud']['pmom_file'] = cloud_pmom_file
        else:
            d['cloud'] = cloud
        save_to_json(filename,d)
    else:
        from Run_libradtran import writeDict
        f = open(filename,'w')
        f.write('UVSpec input file version: {} \n'.format(vv))
        f.write('sza = {} \n'.format(sza))
        f.write('tau = {} \n'.format(tau))
        f.write('ref = {} \n'.format(ref))
        f.close()
        writeDict(geo,'geo',filename)
        writeDict(aero,'aero',filename)
        writeDict(cloud,'cloud',filename)
        writeDict(source,'source',filename)
        writeDict(albedo,'albedo',filename)    


# In[3]:


def writeDict(dict_in, dict_name, filename):
    with open(filename, "a") as f:
        f.write("{}\n".format(dict_name))
        for i,v in dict_in.items():            
            f.write("  {}: {}\n".format(i,v))


# In[ ]:


def read_all_aac():
    """
    Simple program to run the read_aac for the files and paths saved
    """
    pass


# In[102]:


def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False


# In[ ]:


def read_24h_aac(d):
    'function to read out the 24hour files (sol and thm) from a list of files'
    import numpy as np
    import Run_libradtran as RL
    d.sort()
    
    zout=[0,3,100]
    nz = len(zout)
    output = {}
    #output = {'zout':zout,'ext':d['aero']['ext'],'asy':d['aero']['asy'],'ssa':d['aero']['ssa']}
    output['SW_irr_dn_utc'] = np.zeros((nz,24))
    output['SW_irr_up_utc'] = np.zeros((nz,24))
    output['LW_irr_dn_utc'] = np.zeros((nz,24))
    output['LW_irr_up_utc'] = np.zeros((nz,24))  
    output['HH'] = np.zeros((24,1))
    output['file_sol'] = []
    output['file_thm'] = []
    
    for HH in xrange(24):
        output['HH'][HH] = HH
        h = 'HH{:02.0f}'.format(HH)
        for fin in d:
            if h in fin:
                if 'sol' in fin:
                    file_out_sol = fin
                if 'thm' in fin:
                    file_out_thm = fin
        
        output['file_sol'].append(file_out_sol)           
        output['file_thm'].append(file_out_thm)
                    
        try:
            sol = RL.read_libradtran(file_out_sol,zout=zout)
        except IOError:
            print 'File {} not found skip'.format(file_out_sol) 
            if HH==0:
                print file_out_sol
            continue
        except ValueError:
            print 'Problem with file: {}'.format(file_out_sol)
            output['SW_irr_dn_utc'][:,HH] = np.nan
            output['SW_irr_up_utc'][:,HH] = np.nan
                    
        try:
            thm = RL.read_libradtran(file_out_thm,zout=zout)
        except IOError:
            print 'File {} not found skip'.format(file_out_thm) 
            if HH==0:
                print file_out_thm
            continue
        except ValueError:
            print 'Problem with file: {}'.format(file_out_thm)
            output['LW_irr_dn_utc'][:,HH] = np.nan
            output['LW_irr_up_utc'][:,HH] = np.nan
            
            
        output['SW_irr_dn_utc'][:,HH] = sol['direct_down']+sol['diffuse_down']
        output['SW_irr_up_utc'][:,HH] = sol['diffuse_up']
        output['LW_irr_dn_utc'][:,HH] = thm['direct_down']+thm['diffuse_down']
        output['LW_irr_up_utc'][:,HH] = thm['diffuse_up']

    output['SW_irr_dn_avg'] = np.mean(output['SW_irr_dn_utc'],axis=1)
    output['SW_irr_up_avg'] = np.mean(output['SW_irr_up_utc'],axis=1)
    output['LW_irr_dn_avg'] = np.mean(output['LW_irr_dn_utc'],axis=1)
    output['LW_irr_up_avg'] = np.mean(output['LW_irr_up_utc'],axis=1)
    return output        


# In[ ]:


if __name__=='__main__':
    import numpy as np
    ext = np.array([[0.5,0.4,0.3],[0,0,0]])
    ssa = ext*1.9
    asy = ext*1.7
    ext.shape

    z_arr = np.array([3,4])
    wvl_arr = np.array([350.,500.,650.])

    write_aerosol_file_explicit('C:\Users\sleblan2\libradtran/aero.inp',z_arr,ext,ssa,asy,wvl_arr,verbose=True)


# In[49]:


if __name__=='__main__':
    write_cloud_file('C:\Users\sleblan2\libradtran\cloud.inp',10,10,2,3,verbose=True)


# In[57]:


if __name__=='__main__':
    write_albedo_file('C:\Users\sleblan2\libradtran/alb.inp',[500.0,600.0],[0.5,0.6],verbose=True)


# In[109]:


if __name__=='__main__':
    write_input_aac('C:\Users\sleblan2\libradtran/test_input.inp',geo={'zout':[0,3,100],'wvl_range':[202,5600]},aero={},cloud={'tau':10,'ref':10,'phase':'wc','ztop':3,'zbot':2},source={},
                    verbose=True)


# In[44]:


if __name__=='__main__':
    import Run_libradtran
    reload(Run_libradtran)


# In[92]:


if __name__=='__main__':
    mie = make_pmom_inputs()


# In[93]:


if __name__=='__main__':
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


# In[94]:


if __name__=='__main__':
    write_input_aac('C:\Users\sleblan2\libradtran/test_input_aac.inp',geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=True)


# In[39]:


if __name__=='__main__':
    fp = 'C:\Users\sleblan2/Research/libradtran/testing_new/AAC_input_lat06_lon19_DJF_HH17_sol.out'

    import pandas as pd

    d = pd.read_csv(fp,delim_whitespace=True,engine='c')


# In[82]:


if __name__=='__main__':
    d


# In[83]:


if __name__=='__main__':
    dat = np.array(d)

    dat

    #if run_from_ipython():
    #    %timeit rr = np.fromfile(fp,sep=' ').reshape((3,7))
    #    %timeit dd = np.array(pd.read_csv(fp,delim_whitespace=True,engine='c'))
    #    %timeit gg = np.loadtxt(fp)
    #    %timeit gh = np.genfromtxt(fp)

    Run_libradtran.read_libradtran(fp)

