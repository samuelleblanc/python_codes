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
        print('..File finished, closed')

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
        print('..File finished, closed')

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
        print('..File finished, closed')

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
        print('..File finished, closed')

# <codecell>

def write_input_aac(output_file,geo={},aero={},cloud={},source={},albedo={},
                    verbose=False):
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
        cloud: dictionary with cloud properties
            tau: value of cloud optical thickness
            ref: value of effective radius in microns
            ztop: location in km of cloud top
            zbot: location in km of cloud bottom
            phase: either 'ic' for ice cloud or 'wc' for water cloud
        source: dictionary with source properties
            wvl_range: range of wavelengths to model (default [202,500])
            source: can either be thermal or solar
            dat_path: data path to be used. Defaults to pleaides values (/u/sleblan2/libradtran/libRadtran-2.0-beta/data/)
            integrate_values: if set to True (default), then the resulting output parameters are integrated over the wavelength range
                            if set to False, returns per_nm irradiance values
        albedo: dictionary with albedo properties
            create_albedo_file: if true then albedo file is created with the properties defined by alb_wvl and alb (defaults to False)
            albedo_file: path of albedo file to use if already created 
            albedo: value of albedo. Only used if create_albedo_file is false and albedo_file is empty (defaults to 0.29 - Earth's average)
            alb: wavelength dependent value of albedo for use when create_albedo_file is set to True
            alb_wvl: wavelength grid of albedo for use when create_albedo_file is set to True
    
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
    from Run_libradtran import write_aerosol_file_explicit,write_cloud_file,write_albedo_file,merge_dicts
    
    try:
        output = file(output_file,'w')
    except Exception,e:
        print 'Problem with accessing file, return Exception: ',e
        return
    if verbose: print 'Opening input file %s' % output_file
    
    if verbose: print '..setting the dicts to defaults'
    source = merge_dicts({'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/',
                          'source':'solar',
                          'wvl_range':[250.0,500.0],
                          'integrate_values':True},source)
    albedo = merge_dicts({'create_albedo_file':False,
                          'albedo':0.29},albedo)
    geo = merge_dicts({'zout':[0,100]},geo)
    
    if source.get('source')=='solar':
        source['source'] = 'solar '+source['dat_path']+'solar_flux/kurudz_1.0nm.dat per_nm'
    
    if verbose: print '..write out general default values'
    output.write('mol_abs_param sbdart\n')
    output.write('rte_solver disort\n')
    
    if verbose: print '..write out source dict values'
    if source['integrate_values']:
        output.write('output_process integrate\n')
    else:
        output.write('output_process per_nm\n')
    output.write('data_files_path \t %s\n' % source['dat_path'])
    output.write('source \t %s \n' % source['source'])
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
        
    if aero.get('ext'):
        if verbose: print '..write out the aerosol parameters'
        output.write('aerosol_default\n')
        output.write('disort_intcor moments\n') #set to use moments for explicit aerosol file
        aero['file_name'] = output_file+'_aero'
        output.write('aerosol_file explicit \t%s\n' % aero['file_name'])
        write_aerosol_file_explicit(aero['file_name'],aero['z_arr'],aero['ext'],aero['ssa'],aero['asy'],aero['wvl_arr'],verbose=verbose)
    
    if cloud.get('tau'):
        if verbose: print '..write out the cloud properties'
        cloud['file_name'] = output_file+'_cloud'
        if cloud['phase']=='ic':
            if verbose: print '..Ice cloud'
            output.write('ic_file 1D \t %s\n' % cloud['file_name'])
            output.write('ic_properties baum_v36 interpolate\n')
        elif cloud['phase']=='wc':
            if verbose: print '..Liquid water cloud'
            output.write('wc_file 1D \t %s\n' % cloud['file_name'])
            output.write('wc_properties mie interpolate\n')
        else:
            raise ValueError('phase value in cloud dict not recognised')
        write_cloud_file(cloud['file_name'],cloud['tau'],cloud['ref'],cloud['zbot'],cloud['ztop'],verbose=verbose)
    
    output.close()
    if verbose: print 'Finished printing: Closing file'   

# <codecell>

def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

# <codecell>

if __name__=='__main__':

# <codecell>

    import numpy as np

# <codecell>

    ext = np.array([[0.5,0.4,0.3],[0,0,0]])
    ssa = ext*1.9
    asy = ext*1.7

# <codecell>

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
    cloud = {'tau':7.6,
             'ref':12.47,
             'ztop':3.0,
             'zbot':2.0,
             'phase':'wc'}
    source = {'wvl_range':[202,5600],
              'source':'solar',
              'integrate_values':True}
    albedo = {'create_albedo_file':False,
              'albedo':0.2}

# <codecell>

    write_input_aac('C:\Users\sleblan2\libradtran/test_input_aac.inp',geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=True)

# <codecell>

    list(reversed(source['wvl_range']))

