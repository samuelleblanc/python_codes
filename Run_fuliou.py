
# coding: utf-8

# # Info
# Name:  
# 
#     Run_fuliou
# 
# Purpose:  
# 
#     Python modules that combine the different modules used for writing and reading fuliou files. 
#     Used with Fuliou version editied for MOC DARF calculations
#     Based on outputs as matlab functions from John Livingston  procedures 'read_FuLiou_input_file_Calipso_revFeb2015.m' 
# 
# Calling Sequence:
# 
#     import Run_fuliou as Rf
#   
# Input:
# 
#     none at command line
#     see methods of module
# 
# Output:
#    
#     input files for fuliou
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
#     Wrtten: Samuel LeBlanc, NASA Ames, from Santa Cruz, 2017-03-10
#     Modified: 

# # Load up the various methods

# In[ ]:

def write_fuliou_input(output_file,geo={},aero={},albedo={},verbose=False):
    """
    Purpose:
        Writes the fuliou input file with defined defaults, see below
        outputs fuliou input files in ascii format
    
    Input: 
        output_file: full path to file to be written
        geo: dictionary with geometry details
            sza: solar zenith angle (at current time and location)
            lat: latitude
            lon: longitude
            doy: day of year  # for calculating sun earth distance
            year: year (YYYY format)
            month: month (MM format)
            day: day of the month (DD format)
            hour: hour of the day, 24h format, UTC
            minute: minutes of the hour
            second: seconds of the minute
            utc: fractional hours since midnight, UTC (can be a substitute to hour, minute, second)
            pressure: sea level pressure (mb, defaults to 1013.25)
            zout: at what altitude should be outputted, in km, default is 0 and 100
        aero: dictionary with aerosol properties
            ext: value of aerosol extinction coefficient at each altitude and wavelength [alt,wvl]
                 To define the top most layer, ext must be zero at that altitude
            ssa: value of aerosol single scattering albedo at each altitude and wavelength [alt,wvl]
            asy: value of aerosol asymmetry parameter at each altitude and wavelength [alt,wvl]
            z_arr: value of altitudes to use with ext,ssa,asy
            wvl_arr: array of wavelengths in nm
        albedo: dictionary with albedo properties
            albedo: value of albedo. Only used if create_albedo_file is false and albedo_file is empty 
                    (defaults to 0.29 - Earth's average)
            alb: wavelength dependent value of albedo for use when create_albedo_file is set to True
            alb_wvl: wavelength grid of albedo for use when create_albedo_file is set to True
            sea_surface_albedo: (default False) If True, sets the albedo as a 
                                sea surface to be parameterized = 0.037/(1.1*cos(sza)^1.4+0.15)
            land_surface_albedo: (default False) If True, sets the albedo as a
                                 land surface to be parameterized = 0.07/(0.3*cos(sza)^1.4+0.15)
            modis_surface_albedo: (default False) If True, needs to be parameterized using Schaaf et al.
                                  requires modis_albedo dict to be completed with MODIS land data
            modis_albedo: dictionary with land albedo values to be calculated from 
        set_quiet: if True then quiet is set in the input file, if False, quiet is not set. (default True)
    
    Output:
        output_file
    
    Keywords: 
        verbose: (default False) if true prints out info about file writing 
    
    Dependencies:
        numpy
        Run_fuliou (this file)
        map_utils
        datetime
    
    Required files:
        if albedo['modis_surface_albedo'] set to true, files from MCD43GF_CMG files must be present, path defined in modis_albedo
    
    Example:
        
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-03-10, Santa Cruz, CA
                        Translated and modified from rad_FuLiou_input_file_Calipso_revFeb2015.m from matlab to python
                        originally written by John Livingston
    
    """
    import map_utils as mu
    from datetime import datetime, timedelta
    from Run_fuliou import calc_sfc_albedo_Schaaf
    import numpy as np

    # Verify inputs
    if not geo.get('utc'):
        geo['utc'] = geo['hour']+geo['minute']/60.0+geo['second']/3600.0
    else:
        geo['hour'] = int(geo['utc'])*1.0
        geo['minute'] = int((geo['utc']-geo['hour'])*60.0)*1.0
        geo['second'] = int((geo['utc']-geo['hour']-geo['minute']/60.0)*3600.0)*1.0
    geo['datetime'] = datetime(geo['year'],geo['month'],geo['day'],geo['hour'],geo['minute'],geo['second'])
    if not geo.get('pressure'):
        geo['pressure'] = 1013.25
    if not geo.get('sza'):
        geo['sza'],geo['azi'] = mu.get_sza_azi(geo['lat'],geo['lon'],geo['datetime'])
        
    if len(aero['wvl_arr'])!=30:
        raise AttributeError('Length of the wavelength array is not 30, check wavelength input')
    if len(aero['ext'][0,:])!=30:
        raise AttributeError('Length of the extinction array is not 30, check extinction wavelength input')
        
    # Calculate SZA for every 5 minutes
    geo['datetime_fine'] = [datetime(geo['year'],geo['month'],geo['day'])+timedelta(minutes=5*i) for i in range(24*60/5)]
    geo['sza_fine'],geo['azi_fine'] = mu.get_sza_azi(geo['lat'],geo['lon'],geo['datetime_fine'])
    
    # Get the surface albedo for every 5 minutes
    if albedo.get('sea_surface_albedo'):
        albedo['albedo_fine'] = 0.037/(1.1*np.cos(np.array(geo['sza_fine'])*np.pi/180.0)**1.4+0.15)
    if albedo.get('land_surface_albedo'):
        albedo['albedo_fine'] = 0.070/(0.3*np.cos(np.array(geo['sza_fine'])*np.pi/180.0)**1.4+0.15)
    if albedo.get('modis_surface_albedo'):
        fiso,fvol,fgeo,frac_diffuse = prep_albedo_calc(albedo['modis_albedo'],geo['sza_fine'],ext=aero['ext'],z=aero['z_arr'])
        albedo['albedo_fine'] = calc_sfc_albedo_Schaaf(fiso,fvol,fgeo,frac_diffuse,geo['sza_fine'])
    if not albedo.get('albedo_fine'):
        raise AttributeError('No albedo defined please review albedo dict input')
    else:
        albedo['albedo_fine'][np.logical_or((np.isfinite(albedo['albedo_fine'])!=True),(albedo['albedo_fine']<0))]=0.0
        
    with open(output_file,'r') as f:
        # write the title setup line
        f.write('%4i %2i %2i %12.5f %11.4f %11.4f %11.4f %11.4f %11.5f %11.5f\n' % (                geo['year'],geo['month'],geo['day'],geo['year'],geo['utc'],geo['sza'],geo['pressure'],
                min(aero['z_arr']),max(aero['z_arr']),geo['lat'],geo['lon']))
        # write the fine sza
        f.write("\n".join(["".join((' %9.4f' % s for s in geo['sza_fine'][i:i+8])) for i in xrange(0,len(geo['sza_fine']),8)]))
        # write the fine surface albedo    
        f.write("\n".join(["".join((' %9.4f' % s for s in albedo['albedo_fine'][i:i+8]))                            for i in xrange(0,len(albedo['albedo_fine']),8)]))
        # write the aerosol properties on one line for each wavelength
        for i,w in enumerate(aero['wvl_arr']):
            f.write('%12.4e %11.4e %11.4e\n'%(aero['ext'][0,i],aero['ssa'][0,i],aero['asy'][0,i]))


# In[ ]:

def prep_albedo_calc(mod,sza,ext=[],z=[],aod=[]):
    """
    Purpose:
        Prepares the indices for MODIS (Schaaf et al.) black-sky and white-sky albedos 
        using assumed and subsequent surface albedo
    
    Input: 
        
        mod: dict of modis_albedo values from files
        sza: array of solar zenith angles to use for calculations
        ext: array of extinction coefficients used for calculating the AOD and the draction of diffuse light
        z: array of altitude values at which the extinction is defined
        aod: (optional instead of ext and z) array of aod values
        
    Output:
        fiso,fvol,fgeo are single values
        frac_diffuse is a vector (with length=length(SZAIN))that is a function of SZA (hence, for a specific AOD)
    
    Keywords: 
        None
    
    Dependencies:
        numpy
    
    Required files:
        None
    
    Example:
        
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-03-15, Santa Cruz, CA
                        Translated and modified from rad_FuLiou_input_file_Calipso_revFeb2015.m from matlab to python
                        originally written by John Livingston
    """


# In[46]:

def calc_sfc_albedo_Schaaf(fiso,fvol,fgeo,frac_diffuse,SZAin):
    """
    Purpose:
        calculates MODIS (Schaaf et al.) black-sky and white-sky albedos using assumed and subsequent surface albedo
    
    Input: 
        fiso,fvol,fgeo are single values
        frac_diffuse is a vector (with length=length(SZAIN))that is a function of SZA (hence, for a specific AOD)
        SZAin is a vector of SZA
        
    Output:
        surface albedo per sza
    
    Keywords: 
        None
    
    Dependencies:
        numpy
    
    Required files:
        None
    
    Example:
        
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-03-15, Santa Cruz, CA
                        Translated calc_sfc_albedo_Schaaf.m from matlab to python
                        originally written by John Livingston
    """
    import numpy as np
    g0bs = np.array([1.0, -0.007574, -1.284909])
    g1bs = np.array([0.0, -0.070987, -0.166314])
    g2bs = np.array([0.0,  0.307588,  0.041840])
    gws  = np.array([1.0,  0.189184, -1.377622])

    albedo_sfc = np.zeros((len(fiso),len(SZAin)))

    SZAinrad = np.array(SZAin)/180.0*np.pi
    SZAsq = SZAinrad*SZAinrad
    SZAcub = SZAinrad*SZAsq

    for i in xrange(len(fiso)):
        alb_bs = fiso[i]*(g0bs[1] + g1bs[1]*SZAsq + g2bs[1]*SZAcub) +                  fvol[i]*(g0bs[2] + g1bs[2]*SZAsq + g2bs[2]*SZAcub) +                  fgeo[i]*(g0bs[3] + g1bs[3]*SZAsq + g2bs[3]*SZAcub)
        alb_ws = fiso[i]*gws[1] + fvol[i]*gws[2] + fgeo[i]*gws[3]

        albedo_sfc[i,:] = alb_ws*frac_diffuse + (1-frac_diffuse)*alb_bs;

    #now for all SZA>=90 set surface albedo=0
    albedo_sfc[:,SZAin>=90]=0.0
    return albedo_sfc


# In[ ]:




# In[58]:

af = 0.070/(0.3*np.cos(np.array(geo['sza_fine'])*np.pi/180.0)**1.4+0.15)


# In[59]:

af[af<0] = 0.0


# In[60]:

af


# In[88]:

af[np.logical_or((np.isfinite(af)!=True),(af<0))]=0.0


# In[89]:

af


# In[ ]:



