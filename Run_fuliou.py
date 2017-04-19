
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

# In[64]:

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
    if verbose: print 'verifying inputs'
    if not geo.get('utc'):
        geo['utc'] = geo['hour']+geo['minute']/60.0+geo['second']/3600.0
    else:
        geo['hour'] = int(geo['utc'])*1.0
        geo['minute'] = int((geo['utc']-geo['hour'])*60.0)*1.0
        geo['second'] = int((geo['utc']-geo['hour']-geo['minute']/60.0)*3600.0)*1.0
    geo['datetime'] = datetime(int(geo['year']),int(geo['month']),int(geo['day']),
                               int(geo['hour']),int(geo['minute']),int(geo['second']))
    if not geo.get('pressure'):
        geo['pressure'] = 1013.25
    if hasattr(geo['lat'],'__len__'):
        geo['lat'], geo['lon'] = geo['lat'][0],geo['lon'][0]
    if not geo.get('sza'):
        geo['sza'],geo['azi'] = mu.get_sza_azi(geo['lat'],geo['lon'],[geo['datetime']])
        geo['sza'] = geo['sza'][0]
        
    if len(aero['wvl_arr'])!=30:
        raise AttributeError('Length of the wavelength array is not 30, check wavelength input')
    if len(aero['ext'][0,:])!=30:
        raise AttributeError('Length of the extinction array is not 30, check extinction wavelength input')
        
    # verify that inputs are physical
    try: aero['ext'][0,aero['ext'][0,:]<0.0] = 0.0
    except: pass
    try: aero['ssa'][0,aero['ssa'][0,:]<0.0] = 0.0
    except: pass
    try: aero['ssa'][0,aero['ssa'][0,:]>1.0] = 1.0
    except: pass
    try: aero['asy'][0,aero['asy'][0,:]<0.0] = 0.0
    except: pass
    try: aero['asy'][0,aero['asy'][0,:]>1.0] = 1.0
    except: pass
        
    # Calculate SZA for every 5 minutes
    if verbose: print 'calculating sza for every 5 minutes'
    geo['datetime_fine'] = [datetime(geo['year'],geo['month'],geo['day'])+timedelta(minutes=5*i) for i in xrange(24*60/5+1)]
    geo['sza_fine'],geo['azi_fine'] = mu.get_sza_azi(geo['lat'],geo['lon'],geo['datetime_fine'])
    
    # Get the surface albedo for every 5 minutes
    if albedo.get('sea_surface_albedo'):
        if verbose: print 'albedo set to sea_surface_albedo'
        albedo['albedo_fine'] = 0.037/(1.1*np.cos(np.array(geo['sza_fine'])*np.pi/180.0)**1.4+0.15)
    if albedo.get('land_surface_albedo'):
        if verbose: print 'albedo set to land_surface_albedo'
        albedo['albedo_fine'] = 0.070/(0.3*np.cos(np.array(geo['sza_fine'])*np.pi/180.0)**1.4+0.15)
    if albedo.get('modis_surface_albedo'):
        if verbose: print 'albedo set to modis_surface_albedo'
        fiso = albedo['modis_albedo']['mean_iso']
        fvol = albedo['modis_albedo']['mean_vol']
        fgeo = albedo['modis_albedo']['mean_geo']
        frac_diffuse = get_frac_diffuse(albedo['modis_albedo'],geo['sza_fine'],ext=aero['ext'],z=aero['z_arr'])
        albedo['albedo_fine'] = calc_sfc_albedo_Schaaf(fiso,fvol,fgeo,frac_diffuse,geo['sza_fine'])
    if not 'albedo_fine' in albedo:
        raise AttributeError('No albedo defined please review albedo dict input')
    else:
        albedo['albedo_fine'][np.logical_or((np.isfinite(albedo['albedo_fine'])!=True),(albedo['albedo_fine']<0))]=0.0
        
    if verbose: print 'writing to file:'+output_file
    with open(output_file,'w') as f:
        # write the title setup line
        f.write('%4i %2i %2i %12.5f %11.4f %11.4f %11.4f %11.4f %11.5f %11.5f\n' % (                geo['year'],geo['month'],geo['day'],geo['utc'],geo['sza'],geo['pressure'],
                min(aero['z_arr']),max(aero['z_arr']),geo['lat'],geo['lon']))
        # write the fine sza
        f.write("\n".join([" ".join(('%9.4f' % s for s in geo['sza_fine'][i:i+8])) for i in xrange(0,len(geo['sza_fine']),8)]))
        f.write('\n')
        # write the fine surface albedo    
        f.write("\n".join([" ".join(('%9.4f' % s for s in albedo['albedo_fine'][i:i+8]))                            for i in xrange(0,len(albedo['albedo_fine']),8)]))
        f.write('\n')
        # write the aerosol properties on one line for each wavelength
        for i,w in enumerate(aero['wvl_arr']):
            f.write('%12.4e %11.4e %11.4e\n'%(aero['ext'][0,i],aero['ssa'][0,i],aero['asy'][0,i]))


# In[438]:

def get_frac_diffuse(mod,sza,ext=[],z=[],aod=None,iwvl_midvis=8):
    """
    Purpose:
        Prepares the indices for MODIS (Schaaf et al.) black-sky and white-sky albedos 
        using assumed and subsequent surface albedo
        Calculates the fractional diffuse protion of the sky from the modis_albedo lut tables.
    
    Input: 
        
        mod: dict of modis_albedo values from files
        sza: array of solar zenith angles to use for calculations
        ext: array of extinction coefficients used for calculating the AOD and the draction of diffuse light
        z: array of altitude values at which the extinction is defined
        aod: (optional instead of ext and z) array of aod values
        iwvl: (defaults to 8) index of the midvisible wavelenght in ext array
        
    Output:
        frac_diffuse is a vector (with length=length(SZAIN))that is a function of SZA (hence, for a specific AOD)
    
    Keywords: 
        None
    
    Dependencies:
        numpy
        scipy.interpolate
    
    Required files:
        None
    
    Example:
        
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-03-15, Santa Cruz, CA
                        Translated and modified from rad_FuLiou_input_file_Calipso_revFeb2015.m from matlab to python
                        originally written by John Livingston
    """
    import numpy as np
    from scipy import interpolate
    
    dz = z[-1]-z[0]
    if not aod:
        aod = dz*ext[0,iwvl_midvis]
    
    #Now interpolate the table at the right AOD
    new_szas = np.zeros_like(mod['table_SZA'])
    for i in xrange(len(mod['table_SZA'])): 
        ftab = interpolate.interp1d(mod['table_AOD'],mod['table_fracdiffuse'][i,:])
        new_szas[i] = ftab(aod)
        
    fftab = interpolate.interp1d(mod['table_SZA'],new_szas)
    frac_diffuse = fftab(sza)
    
    return frac_diffuse


# In[439]:

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
        alb_bs = fiso[i]*(g0bs[0] + g1bs[0]*SZAsq + g2bs[0]*SZAcub) +                  fvol[i]*(g0bs[1] + g1bs[1]*SZAsq + g2bs[2]*SZAcub) +                  fgeo[i]*(g0bs[2] + g1bs[2]*SZAsq + g2bs[3]*SZAcub)
        alb_ws = fiso[i]*gws[0] + fvol[i]*gws[1] + fgeo[i]*gws[2]

        albedo_sfc[i,:] = alb_ws*frac_diffuse + (1-frac_diffuse)*alb_bs;

    #now for all SZA>=90 set surface albedo=0
    albedo_sfc[:,SZAin>=90]=0.0
    return albedo_sfc


# In[440]:

def Prep_DARE_single_sol(fname,f_calipso,fp_rtm,fp_fuliou,fp_alb=None,surface_type='ocean',vv='v1'):
    """
    Purpose:
        Main function to create the files for the DARE calculations for a single solx file defined by fname 
    
    Input: 
        fname: full file path for matlab file solution for single solutions
               (e.g., MOCsolutions20150508T183717_19374_x20080x2D070x2D11120x3A250x2CPoint0x2313387030x2F25645720x2CH0x.mat)
        f_calipso: full file path for Calipso file
        fp_rtm: filepath for baseline rtm folder (subset of it will have input, output) 
        fp_fuliou: full filepath for fuliou executable
        
    Output:
        input files for fuliou
        list file for calling fuliou from NAS
    
    Keywords: 
        fp_alb: File path for MODIS surface albedo, must be defined if surfacetype is land_MODIS
        surface_type: (default to ocean) can be either ocean, land, or land_MODIS
        vv: (default to v1) version number
    
    Dependencies:
        numpy
        Run_fuliou (this file)
        os
        scipy
        datetime
        Sp_parameters
    
    Required files:
        matlab MOC solution for single solx 
        Calipso matlab file
        yohei_MOC_lambda.mat file
    
    Example:
        
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2017-03-17, Santa Cruz, CA
                        migrated to python from matlab based on codes in read_FuLiou_input_file_Calipso_revFeb2015.m
                        originally written by John Livingston
    """
    import numpy as np
    import scipy.io as sio
    import os
    from datetime import datetime
    from Run_fuliou import get_MODIS_surf_albedo
    from Sp_parameters import startprogress, progress, endprogress
    
    # load the materials
    sol = sio.loadmat(fname)
    name = fname.split(os.path.sep)[-1]
    da,num,xt = name.split('_')
    num = int(num)-1 #convert from matlab indexing to python
    da = da.lstrip('MOCsolutions')
    xt = xt.rstrip('.mat')
    
    if verbose: print 'Preparing the input files for pixel: {}'.format(num)
    
    lm = sio.loadmat(fp_rtm+'yohei_MOC_lambda.mat')
    cal = sio.loadmat(f_calipso)
    
    # prep the write out paths
    fp_in = os.path.join(fp_rtm,'input','MOC_1solx_DARE_{vv}_{num}'.format(vv=vv,num=num))
    fp_out = os.path.join(fp_rtm,'output','MOC_1solx_DARE_{vv}_{num}'.format(vv=vv,num=num))
    if not os.path.exists(fp_in):
        os.makedirs(fp_in)
    if not os.path.exists(fp_out):
        os.makedirs(fp_out)
        
    f_list = open(os.path.join(fp_rtm,'run','list_MOC_single_solx_DARE_{vv}_{num}.sh'.format(vv=vv,num=num)),'w')
    print f_list.name
    
    # prep the standard input definitions
    tt = datetime.strptime(cal['calipso_date'][num],'%Y-%m-%d')
    
    geo = {'lat':cal['calipso_lat_oneday'][num], 'lon':cal['calipso_lon_oneday'][num],
           'utc':cal['calipso_time_oneday'][num][0]/100.0,'year':tt.year,'month':tt.month,'day':tt.day}
    aero = {'wvl_arr':sol['solutions']['lambda'][0,0][0]*1000.0,
            'z_arr':[cal['calipso_zmin_oneday'][num][0],cal['calipso_zmax_oneday'][num][0]]}
    if surface_type=='ocean':
        albedo = {'sea_surface_albedo':True,'land_surface_albedo':False,'modis_surface_albedo':False}
    elif surface_type=='land':
        albedo = {'sea_surface_albedo':False,'land_surface_albedo':True,'modis_surface_albedo':False}
    elif surface_type=='land_MODIS':
        albedo = {'sea_surface_albedo':False,'land_surface_albedo':True,'modis_surface_albedo':False}
        if not fp_alb:
            raise ValueError('fp_alb is not set, must be set to find the MCD43GF files')
        albedo['modis_albedo'] = get_MODIS_surf_albedo(fp_alb,tt.timetuple().tm_yday,geo['lat'],geo['lon'],year_of_MODIS=2007)
    else:
        raise ValueError("surface_type can only be 'ocean', 'land', or 'land_MODIS'")

    # Run through the possible solutions to the files and write out the fuliou input files
    startprogress('Writing MOC files')
    for i in xrange(len(sol['ssa'])):
        input_file = os.path.join(fp_in,'MOC_{num}_{i:04d}.datin'.format(num=num,i=i))
        output_file = os.path.join(fp_out,'MOC_{num}_{i:04d}.wrt'.format(num=num,i=i))
        aero['ssa'] = np.array([sol['ssa'][i,:],sol['ssa'][i,:]])
        aero['asy'] = np.array([sol['asym'][i,:],sol['asym'][i,:]])
        aero['ext'] = np.array([sol['ext'][i,:],sol['ext'][i,:]])
        write_fuliou_input(input_file,geo=geo,aero=aero,albedo=albedo,verbose=False)
        
        progress(float(i)/float(len(sol['ssa']))*100.0)
        f_list.write(fp_fuliou+' '+input_file+' '+output_file+'\n')
    
    endprogress()
    sel = sol['solutions']['select'][0,0]
    i_str = ['m1','s0','p1']
    for ie in [-1,0,1]:
        for ia in [-1,0,1]:
            for im in [-1,0,1]:
                form = {'num':num,'e':i_str[ie+1],'a':i_str[ia+1],'s':i_str[im+1]}
                input_file = os.path.join(fp_in,'MOC_{num}_{e}{a}{s}.datin'.format(**form))
                output_file = os.path.join(fp_out,'MOC_{num}_{e}{a}{s}.wrt'.format(**form))
                aero['ssa'] = np.array([sel['ssa'][0,0][0,:]+ia*sel['ssa'][0,0][1,:]]*2)
                aero['asy'] = np.array([sel['asym'][0,0][0,:]+im*sel['asym'][0,0][1,:]]*2)
                aero['ext'] = np.array([sel['ext'][0,0][0,:]+ie*sel['ext'][0,0][1,:]]*2)
                write_fuliou_input(input_file,geo=geo,aero=aero,albedo=albedo,verbose=False)

                print '{e}{a}{s}'.format(**form)
                f_list.write(fp_fuliou+' '+input_file+' '+output_file+'\n')
    
    f_list.close()


# In[441]:

def get_MODIS_surf_albedo(fp,doy,lat,lon,year_of_MODIS=2007):
    """
    Purpose:
        helper function to get the details of the surface albedo from MODIS albedo 
    
    Input: 
        fp: path to get the saved files
        doy: day of year of the requested data
        lat: latitude of the surface albedo to get
        lon: longitude of the surface albedo to get
        
    Output:
        modis_albedo: dict containing the values used to create the albedo from a subset of aod, and sza
                      mean_geo, mean_iso, mean_vol : mean values within the lat and lon space
                      std_geo, std_iso, std_vol: standard deviation of the values within the lat and lon space
                      table_AOD: lut tabel of fractional diffuse proportion, value of the AODs
                      table_SZA: lut table value of sza
                      table_fracdiffuse: lut table of the fractional diffuse proportion linked to AOD and SZA                      
                          
    Keywords: 
        year_of_MODIS: (default to 2007) year as integer of the MODIS albedo files to use
    
    Dependencies:
        numpy
        Run_fuliou (this file)
        map_utils

    
    Required files:
        skyl_lut_bbshortwave.dat file 
        MODIS albedo files for the correct day:
            MCD43GF_geo_shortwave_%03d_2007.hdf
            MCD43GF_iso_shortwave_%03d_2007.hdf
            MCD43GF_vol_shortwave_%03d_2007.hdf
    
    Example:
        
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2017-03-22, Santa Cruz, CA
                        migrated to python from matlab based on codes in read_FuLiou_input_file_Calipso_revFeb2015.m
                        originally written by John Livingston
    """
    import numpy as np
    from map_utils import shoot
    from Run_fuliou import load_hdf_spec
    
    MODIS_start_days = np.append(np.arange(1,367,8),367)
    modis_albedo = {}
    
    # read file of the fractional diffuse for bbshortwave per aod and sza
    f = fp+'skyl_lut_bbshortwave.dat'
    modis_albedo['table_AOD'] = np.arange(0,0.99,0.02)
    r = np.genfromtxt(f,skip_header=2)
    modis_albedo['table_SZA'] = r[:,0]
    modis_albedo['table_fracdiffuse'] = r[:,1:-1]
    
    # get the MODIS file day
    MODIS_select_day = MODIS_start_days[(doy-MODIS_start_days)>=0][-1]
    # get the lat lon range from the center point
    lats,lons = np.arange(4).astype(float),np.arange(4).astype(float)
    for i in range(4): lats[i],lons[i],_ = shoot(lat,lon,45.0+i*90.0,maxdist=20.000*np.sqrt(2.0))
        
    # Select the proper grid points to load from the MCD43GF gapfilled albedo files
    latgrid = np.arange(90,-90,-30.0/3600.0)
    longrid = np.arange(-180.0,180,30.0/3600.0)
    ix = np.where((latgrid>=lats.min())&(latgrid<=lats.max()))[0]
    iy = np.where((longrid>=lons.min())&(longrid<=lons.max()))[0]
    
    # assure that all the grid points are within the lats/lons
    # not used for now
    
    # Load the modis gap filled albedo files
    buffer_geo = load_hdf_spec(fp+'MCD43GF_geo_shortwave_%03d_%d.hdf'%(MODIS_select_day,year_of_MODIS),ix,iy)
    buffer_iso = load_hdf_spec(fp+'MCD43GF_iso_shortwave_%03d_%d.hdf'%(MODIS_select_day,year_of_MODIS),ix,iy)
    buffer_vol = load_hdf_spec(fp+'MCD43GF_vol_shortwave_%03d_%d.hdf'%(MODIS_select_day,year_of_MODIS),ix,iy)
    
    modis_albedo['mean_geo'] = np.nanmean(buffer_geo)
    modis_albedo['mean_iso'] = np.nanmean(buffer_iso)
    modis_albedo['mean_vol'] = np.nanmean(buffer_vol)
    modis_albedo['std_geo'] = np.nanstd(buffer_geo)
    modis_albedo['std_iso'] = np.nanstd(buffer_geo)
    modis_albedo['std_vol'] = np.nanstd(buffer_geo)
    
    return modis_albedo


# In[442]:

def load_hdf_spec(filename,ix,iy,data_name='MCD43GF_CMG'):
    """
     Purpose:
        Simple hdf load file for MCD43GF files. Chose only specific indices to load 
    
    Input: 
        filename: full file path and name of the hdf file to load
        ix: array of indices to be loaded in the x direction
        iy: array of indices to be loaded in the y direction
        
    Output:
        dat: np.array of the requested data at the indices
    
    Keywords: 
        data_name: (defaults to MCD43GF_CMG) the name of the dataset to load within the filename
    
    Dependencies:
        numpy
        pyhdf
    
    Required files:
        filename
        
    Example:
        
        >>b = load_hdf_spec(fp+'MCD43GF_geo_shortwave_193_2007.hdf',[200,201,202],[503,504,505,506])
        >>b
        array([[ nan,  nan,  nan,  nan],
               [ nan,  nan,  nan,  nan],
               [ nan,  nan,  nan,  nan]])
        >>b.shape
        (3L, 4L)
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2017-03-22, Santa Cruz, CA
    """
    import numpy as np
    from pyhdf.SD import SD, SDC
    hdf = SD(filename, SDC.READ)
    if hasattr(ix,'__len__'):
        if (len(ix)-1)<(ix[-1]-ix[0]):
            raise ValueError('ix is not a contiguous array')
        if (len(iy)-1)<(iy[-1]-iy[0]):
            raise ValueError('iy is not a contiguous array')
        dat = hdf.select(data_name).get(start=(ix[0],iy[0]),count=(ix[-1]-ix[0]+1,iy[-1]-iy[0]+1))
    else:
        dat = hdf.select(data_name).get(start=(ix,iy),count=(1,1))
    dat = dat.astype(float)
    dat[dat==32767] = np.nan
    dat = dat/1000.0
    return dat    


# In[5]:

def read_DARE_single_sol(fname,fp_rtm,fp_save,vv='v1',verbose=False):
    """
    Purpose:
        Read the DARE fuliou output files for a single solution
        saves the results as a .mat file
    
    Input: 
        fname: full file path for matlab file solution for single solutions
               (e.g., MOCsolutions20150508T183717_19374_x20080x2D070x2D11120x3A250x2CPoint0x2313387030x2F25645720x2CH0x.mat)
        fp_rtm: base folder path (adds the /output/MOC_...)
        num: number of the single solution
        fp_save: full file path of the folder to contain the saved mat file
        
    Output:
        None
    
    Keywords: 
        vv: (default v1) version of the files 
        verbose: (default False) if true, then prints out to screen messages
    
    Dependencies:
        os
        numpy
        scipy.io
        Run_fuliou (this file)
    
    Required files:
        fuliou file output
        
    Example:
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2017-03-24, Santa Cruz, CA
    """
    import os
    import scipy.io as sio
    from Run_fuliou import read_fuliou_output, analyse_fuliou_output
    
    # load the materials
    if verbose: print 'loading mat file for extra info on source: '+fname
    sol = sio.loadmat(fname)
    name = fname.split(os.path.sep)[-1]
    da,num,xt = name.split('_')
    num = int(num)-1 #convert from matlab indexing to python
    da = da.lstrip('MOCsolutions')
    xt = xt.rstrip('.mat')
    
    if verbose: print 'Opening the file at pixel number: {}'.format(num)
    
    fp_out = os.path.join(fp_rtm,'output','MOC_1solx_DARE_{vv}_{num}'.format(vv=vv,num=num))
    d = []
    
    if verbose: print 'running through all the solutions'
    for i,s in enumerate(sol['ssa']):
        if verbose: print i
        output_file = os.path.join(fp_out,'MOC_{num}_{i:04d}.wrt'.format(num=num,i=i))
        d.append(read_fuliou_output(output_file))
        d[i]['ssa'] = sol['ssa'][i,:]
        d[i]['asy'] = sol['asym'][i,:]
        d[i]['ext'] = sol['ext'][i,:]
        analyse_fuliou_output(d[i])
    
    if verbose: print 'starting on select solutions'
    ds = {}
    sel = sol['solutions']['select'][0,0]
    i_str = ['m1','s0','p1']
    for ie in [-1,0,1]:
        for ia in [-1,0,1]:
            for im in [-1,0,1]:
                form = {'num':num,'e':i_str[ie+1],'a':i_str[ia+1],'s':i_str[im+1]}
                val = '{e}{a}{s}'.format(**form)
                if verbose: print val
                output_file = os.path.join(fp_out,'MOC_{num}_{e}{a}{s}.wrt'.format(**form))
                ds[val] = {}
                ds[val]['ssa'] = [sel['ssa'][0,0][0,:]+ia*sel['ssa'][0,0][1,:]]
                ds[val]['asy'] = [sel['asym'][0,0][0,:]+im*sel['asym'][0,0][1,:]]
                ds[val]['ext'] = [sel['ext'][0,0][0,:]+ie*sel['ext'][0,0][1,:]]
                ds[val] = read_fuliou_output(output_file)
                analyse_fuliou_output(ds[val])
    
    f_out = os.path.join(fp_save,'MOC_1solx_DARE_{vv}_{num}.mat'.format(vv=vv,num=num))
    if verbose: print 'saving to matlab file: '+f_out
    sio.savemat(f_out,{'solutions':d,'select':ds})    


# In[2]:

def read_fuliou_output(fname,verbose=False):
    """
    Purpose:
        Reads a single fuliou output file
    
    Input: 
        fname: full file path to be read (*.wrt)
        
    Output:
        d: dict of returned values in the file
    
    Keywords: 
        verbose: (default False) If True prints out the operation steps 
    
    Dependencies:
        numpy
        sys
    
    Required files:
        fname file 
        
    Example:
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2017-03-24, Santa Cruz, CA
    """
    import numpy as np
    import sys
    d = {}
    if verbose: print 'opening :'+fname
    d['fname'] = fname
    try:
        d['num'] = int(fname.split('.')[0].split('_')[-1])
    except:
        d['num'] = -9999
    with open(fname,'r') as f:
        # read header line with info
        if verbose: print 'opening header lines'
        line = f.readline()
        try:
            s = line.split()
            d['year'],d['month'],d['day'],d['utc'],d['sza'],d['pressure'],d['zmin'],d['zmax'],d['lat'],d['lon'] = map(float,s)
        except:
            w = [0,5,3,3,13,12,12,12,12,12,12]
            s = [ line[ sum(w[0:i]) : ii ].strip() for i,ii in enumerate(np.cumsum(w))][1:]
            d['year'],d['month'],d['day'],d['utc'],d['sza'],d['pressure'],d['zmin'],d['zmax'],d['lat'],d['lon'] = map(float,s)
        d['year'],d['month'],d['day'] = int(d['year']),int(d['month']),int(d['day'])
        # define the variables to be read
        d['cosSZA'],d['AOD550'],d['swdn17lev_aer'],d['swup17lev_aer'] = [],[],[],[]
        d['swdntoa_aer'],d['swuptoa_aer'],d['directsfc_aer'],d['diffusesfc_aer'] = [],[],[],[]
        d['swdn17lev_noaer'],d['swup17lev_noaer'] = [],[]
        d['swdntoa_noaer'],d['swuptoa_noaer'],d['directsfc_noaer'],d['diffusesfc_noaer'] = [],[],[],[]
        i = 0 
        while True:
            line = f.readline()
            if not line: break
            if verbose: sys.stdout.write('{}..'.format(i))
            c,nul = map(float,line.split())
            a = map(float,f.readline().split())
            d['cosSZA'].append(c)
            d['AOD550'].append(a)
            t = map(float,f.readline().split()[1:])
            n = map(float,f.readline().split()[1:])
            d['swdn17lev_aer'].append(t[1:17])
            d['swup17lev_aer'].append(t[18:34])
            d['swdntoa_aer'].append(t[35])
            d['swuptoa_aer'].append(t[36])
            d['directsfc_aer'].append(t[37])
            d['diffusesfc_aer'].append(t[38])
            d['swdn17lev_noaer'].append(n[1:17])
            d['swup17lev_noaer'].append(n[18:34])
            d['swdntoa_noaer'].append(n[35])
            d['swuptoa_noaer'].append(n[36])
            d['directsfc_noaer'].append(n[37])
            d['diffusesfc_noaer'].append(n[38])
            i+=1
            
        # now make sure it returns as numpy arrays
        if verbose: print '...\ntransforming to numpy arrays'
        for m in d.keys():
            if not m in ['year','month','day','utc','sza','pressure','zmin','zmax','lat','lon']:
                d[m] = np.array(d[m])
        return d


# In[3]:

def analyse_fuliou_output(d,smaller=True):
    """
    Purpose:
        Runs typical calculations of the fuliou ouput
    
    Input: 
        d: data dict read from single file
        
    Output:
        None but modifies d containing calculated values (24h integrated, interpolation, surface, toa...)
    
    Keywords: 
        smaller: (default True) removes some variables to make the files smaller
    
    Dependencies:
        numpy
    
    Required files:
        None 
        
    Example:
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2017-03-27, Santa Cruz, CA
    """
    import numpy as np
    from scipy.interpolate import interp1d
    
    dt = np.arange(0,24*60.0+5.0,5.0)/60.0
    isun = d['cosSZA']>=0
    
    # no data, return NaNs
    if len(d['diffusesfc_noaer'])<1:
        d['swdnsfc_aer_118_instant'] = np.nan
        d['swtoaup_aer_118_instant'] = np.nan
        d['swtoaup_aer_118_24hr'] = np.nan
        d['swtoaup_noaer_118_24hr'] = np.nan
        d['dF_toa_24hr'] = np.nan
        d['dF_sfc_24hr'] = np.nan
        d['dF_17lev_24hr'] = np.nan
        
    # get the delta z
    try:
        d['deltaz'] = d['zmax']-d['zmin']
    except:
        pass
    
    #simplify the AOD
    d['AOD550'] = d['AOD550'][0][0]
    
    # get the instant values
    fx_dn = interp1d(dt,d['swdn17lev_aer'][:,0])
    fx_up = interp1d(dt,d['swup17lev_aer'][:,0])
    d['swdnsfc_aer_118_instant'] = float(fx_dn(d['utc']))
    d['swtoaup_aer_118_instant'] = float(fx_up(d['utc']))
    
    d['swnet17lev_aer_118'] = d['swdn17lev_aer']-d['swup17lev_aer']
    d['swnet17lev_noaer_118'] = d['swdn17lev_noaer']-d['swup17lev_noaer']
    
    d['dF_toa_all'] = d['swuptoa_noaer']-d['swuptoa_aer']
    d['dF_17lev_all'] = d['swnet17lev_aer_118']-d['swnet17lev_noaer_118']  
    d['dF_sfc_all'] = d['dF_17lev_all'][:,0]
    
    # handle only the contiguous daytime
    isub = np.split(np.where(isun)[0],np.where(np.diff(np.where(isun)[0])!= 1)[0]+1)
    if len(isub) > 1:
        dt_sub = np.concatenate((dt[isub[1]],dt[isub[0]]+24.0))
        ii = np.concatenate((isub[1],isub[0]))
    else:
        ii = isub[0]
        dt_sub = dt[ii]
        
    # calculate the integrals    
    d['dF_toa_24hr'] = np.trapz(d['dF_toa_all'][ii],x=dt_sub)/24.0
    d['dF_sfc_24hr'] = np.trapz(d['dF_sfc_all'][ii],x=dt_sub)/24.0
    d['dF_17lev_24hr'] = np.trapz(d['dF_17lev_all'][ii,:],x=dt_sub,axis=0)/24.0       
    d['swtoaup_noaer_118_24hr'] = np.trapz(d['swuptoa_noaer'][ii],x=dt_sub)/24.0
    d['swtoaup_aer_118_24hr'] = np.trapz(d['swuptoa_aer'][ii],x=dt_sub)/24.0
    
    # get the instantaneous values
    iutc = np.argmin(abs(dt-d['utc']))
    d['dF_toa_instant'] = d['dF_toa_all'][iutc]
    d['dF_sfc_instant'] = d['dF_sfc_all'][iutc]
    
    
    if smaller:
        nn = ['swup17lev_aer','swdn17lev_aer','swdn17lev_noaer','swup17lev_noaer',
              'swuptoa_noaer','swuptoa_aer','swnet17lev_aer_118','swnet17lev_noaer_118','dF_17lev_all']
        for il in nn:
            d.pop(il,None)


# In[ ]:

def read_analyse(filename):
    'combine the read and analyse statement, return analysed dict'
    import Run_fuliou as rf
    try:
        d = rf.read_fuliou_output(filename)
        rf.analyse_fuliou_output(d)
    except:
        pass
    return d


# In[3]:

def read_fuliou_moc(fp,fp_save):
    """
    Purpose:
        Read the full fuliou solutions for either Ocean, DarkTarget, or DeepBlue
        Uses parallel processing for faster reading times
        Reads the output wrt files
    
    Input: 
        fp: basic folder file path which contains the wrt files
        fp_save: save file path
        
    Output:
        None but saves to file
    
    Keywords: 
        None
    
    Dependencies:
        numpy
        multiprocessing/Pool
        os
        scipy.io
    
    Required files:
        wrt files
        
    Example:
        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2017-04-18, Santa Cruz, CA
    """
    import Run_fuliou as rf
    import numpy as np
    import os
    from multiprocessing import Pool
    import scipy.io as sio
    
    files = [fp+f for f in os.listdir(fp) if f.endswith(".wrt")]
    files.sort()
    n = len(files)
    print 'reading {} files'.format(n)
    
    p = Pool(12)
    results = p.map(rf.read_analyse,files)
    p.close()
    
    
    
    print 'Saving to file:'+fp_save+'intermediate.mat'
    ro = {'results':results}
    sio.savemat(fp_save+'intermediate.mat',ro)
    
    nm_list = {'year_case':'year','month_case':'month','day_case':'day',
               'UT_case':'utc','lon_case':'lon','lat_case':'lat',
               'SZA_case':'sza','deltazkm_case':'deltaz','zkmmin_CALIOP':'zmin',
               'zkmmax_CALIOP':'zmax','AOD550_case':'AOD550',
               'swdnsfc_aer_118_instant':'swdnsfc_aer_118_instant','swtoaup_aer_118_instant':'swtoaup_aer_118_instant',
               'swtoaup_noaer_118_24hr':'swtoaup_noaer_118_24hr','dF_toa_24hr':'dF_toa_24hr',
               'dF_sfc_24hr':'dF_sfc_24hr','swtoaup_aer_118_24hr':'swtoaup_aer_118_24hr'}
    saves = {}
    for a in nm_list.keys():
        saves[a] = np.array([results[i][nm_list[a]] for i in xrange(n)])
    for i in xrange(n):
        if not 'dF_17lev_24hr' in results[i]:
            results[i]['dF_17lev_24hr'] = np.zeros(16)+np.nan
    saves['dF_17lev_24hr'] = np.array([results[i]['dF_17lev_24hr'] for i in xrange(n)])
    
    print 'Saving analysed file: '+fp_save
    sio.savemat(fp_save,saves)
    
    return    


# In[443]:

def run_fuliou_pc():
    from Run_fuliou import Prep_DARE_single_sol
    fname = 'C:\\Users\\sleblan2\\Research\\Calipso\\moc\\MOCsolutions_individual\\'+    'MOCsolutions20150508T183717_19374_x20080x2D070x2D11120x3A250x2CPoint0x2313387030x2F25645720x2CH0x.mat'
    f_calipso = 'C:\\Users\\sleblan2\\Research\\Calipso\\moc\\MOCsolutions_individual\\'+    '2008c_MDQA3p1240nm_OUV388SSAvH_CQACOD0.mat'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\Calipso\\moc\\'
    fp_fuliou = 'C:\\Users\\sleblan2\\Research\\Calipso\\moc\\fuliou.exe'
    Prep_DARE_single_sol(fname,f_calipso,fp_rtm,fp_fuliou,surface_type='ocean',vv='v1')


# In[444]:

def run_fuliou_linux(i=0):
    from Run_fuliou import Prep_DARE_single_sol
    fname = '/nobackup/sleblan2/MOCfolder/moc_single_solution/'
    ff = ['MOCsolutions20150508T183717_19374_x20080x2D070x2D11120x3A250x2CPoint0x2313387030x2F25645720x2CH0x.mat',
          'MOCsolutions20150508T183717_22135_x20080x2D080x2D10170x3A320x2CPoint0x2315421760x2F25645720x2CH0x.mat']
    fname = fname + ff[i]
    f_calipso = '/nobackup/sleblan2/MOCfolder/moc_single_solution/'+    '2008c_MDQA3p1240nm_OUV388SSAvH_CQACOD0.mat'
    fp_rtm = '/nobackup/sleblan2/MOCfolder/'
    fp_fuliou = '/u/sleblan2/fuliou/v20170324/fuliou'
    fp_alb = '/nobackup/sleblan2/AAC_DARF/surface_albedo/'
    
    print 'Starting the prepr DARE single solx for fuliou for file: '+fname
    Prep_DARE_single_sol(fname,f_calipso,fp_rtm,fp_fuliou,fp_alb=fp_alb,surface_type='land_MODIS',vv='v1')


# In[1]:

def read_fuliou_linux(i=0):
    from Run_fuliou import read_DARE_single_sol
    fname = '/nobackup/sleblan2/MOCfolder/moc_single_solution/'
    ff = ['MOCsolutions20150508T183717_19374_x20080x2D070x2D11120x3A250x2CPoint0x2313387030x2F25645720x2CH0x.mat',
          'MOCsolutions20150508T183717_22135_x20080x2D080x2D10170x3A320x2CPoint0x2315421760x2F25645720x2CH0x.mat']
    fname = fname + ff[i]
    fp_rtm = '/nobackup/sleblan2/MOCfolder/'
    
    print 'Starting to read fuliou DARE single solx for file: '+fname
    read_DARE_single_sol(fname,fp_rtm,fp_rtm,vv='v1',verbose=True)


# In[445]:

if __name__ == '__main__':
    import argparse
    long_description = """    Prepare or read the fuliou DARE calculations"""
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

