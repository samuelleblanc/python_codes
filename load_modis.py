# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def __init__():
    """
       Collection of codes to load and analyze various data
       
           - modis 
           - emas 
           - cpl_layers text file
           - apr2 files
           - hdf files
           
        details are in the info of each module
    """
    pass

# <codecell>

def load_modis(geofile,datfile):
    """
    Name:

        load_modis
    
    Purpose:

        to compile functions required to load Modis files
        from within another script
    
    Calling Sequence:

        modis,modis_dict = load_modis(geofile,datfile) 
    
    Input: 
  
        geofile name
        datfile name (hdf files)
    
    Output:

        modis dictionary with tau, ref, etau, eref, phase, qa
        modis_dicts : metadate for each of the variables
    
    Keywords: 

        none
    
    Dependencies:

        gdal
        numpy
        gc: for clearing the garbage
    
    Required files:
   
        geo and dat files
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2014-12-08, NASA Ames
        
    """
    import numpy as np
    from osgeo import gdal
    from Sp_parameters import startprogress, progress, endprogress
    modis_values = (#('cloud_top',57),
                    ('phase',53),
          #          ('cloud_top_temp',58),
                    ('ref',66),
                    ('tau',72),
           #         ('cwp',82),
                    ('eref',90),
                    ('etau',93),
            #        ('ecwp',96),
                    ('multi_layer',105),
                    ('qa',123),
             #       ('cloud_mask',110)
                    )
    geosds = gdal.Open(geofile)
    datsds = gdal.Open(datfile)
    geosub = geosds.GetSubDatasets()
    datsub = datsds.GetSubDatasets()
    print 'Outputting the Geo subdatasets:'
    for i in range(len(geosub)):
        print str(i)+': '+geosub[i][1]
    print 'Outputting the Data subdatasets:'
    for i in range(len(datsub)):
        if any(i in val for val in modis_values):
            print '\x1b[1;36m%i: %s\x1b[0m' %(i,datsub[i][1])
        else:
            print str(i)+': '+datsub[i][1]
    latsds = gdal.Open(geosub[12][0],gdal.GA_ReadOnly)
    lonsds = gdal.Open(geosub[13][0],gdal.GA_ReadOnly)
    szasds = gdal.Open(geosub[21][0],gdal.GA_ReadOnly)
    modis = dict()
    modis['lat'] = latsds.ReadAsArray()
    modis['lon'] = lonsds.ReadAsArray()
    modis['sza'] = szasds.ReadAsArray()
    print modis['lon'].shape
    meta = datsds.GetMetadata() 
    import gc; gc.collect()
    modis_dicts = dict()
    startprogress('Running through modis values')
    for i,j in modis_values:
        sds = gdal.Open(datsub[j][0])
        modis_dicts[i] = sds.GetMetadata()
        modis[i] = np.array(sds.ReadAsArray())
        makenan = True
        bad_points = np.where(modis[i] == float(modis_dicts[i]['_FillValue']))
        scale = float(modis_dicts[i]['scale_factor'])
        offset = float(modis_dicts[i]['add_offset'])
       # print 'MODIS array: %s, type: %s' % (i, modis[i].dtype)
        if scale.is_integer():
            scale = int(scale)
            makenan = False
        if scale != 1 and offset == 0:
            modis[i] = modis[i]*scale+offset
        if makenan:
            modis[i][bad_points] = np.nan
        progress(float(tuple(i[0] for i in modis_values).index(i))/len(modis_values)*100.)
    endprogress()
    print modis.keys()
    del geosds, datsds,sds,lonsds,latsds,geosub,datsub
    return modis,modis_dicts

# <codecell>

def load_ict(fname,header=False):
    """
    Simple ict file loader
    created specifically to load the files from the iwg1 on board the G1 during TCAP, may work with others...
    """
    from datetime import datetime
    import numpy as np
    f = open(fname,'r')
    lines = f.readlines()
    first = lines[0]
    num2skip = int(first.strip().split(',')[0])
    header = lines[0:num2skip]
    factor = map(float,header[10].strip().split(','))
    f.close()
    if any([i!=1 for i in factor]):
        print('Some Scaling factors are not equal to one, Please check the factors:')
        print factor
    def mktime(txt):
        return datetime.strptime(txt,'%Y-%m-%d %H:%M:%S')
    def utctime(seconds_utc):
        return float(seconds_utc)/3600.
    conv = {"Date_Time":mktime, "UTC":utctime, "Start_UTC":utctime, "TIME_UTC":utctime, "UTC_mid":utctime}
    data = np.genfromtxt(fname,names=True,delimiter=',',skip_header=num2skip-1,converters=conv)
    print data.dtype.names
    #scale the values by using the scale factors
    for i,name in enumerate(data.dtype.names):
        if i>0:
            if factor[i-1]!=float(1):
                data[name] = data[name]*factor[i-1]    
    if header:
        return data, header
    else:
        return data

# <codecell>

def modis_qa(qa_array):
    """
    modis qa data parser.
    input of qa numpy array
    output structure of qa arrays
    """
    bin8 = lambda x : ''.join(reversed( [str((x >> i) & 1) for i in range(8)] ) )
    
    
    

# <codecell>

def mat2py_time(matlab_datenum):
    "convert a matlab datenum to a python datetime object. Works on numpy arrays of datenum"
    from datetime import datetime, timedelta
    #matlab_datenum = 731965.04835648148
    m2ptime = lambda tmat: datetime.fromordinal(int(tmat)) + timedelta(days=tmat%1) - timedelta(days = 366)
    try:
        python_datetime = m2ptim2(matlab_datenum)
    except:
        import numpy as np
        python_datetime = np.array([m2ptime(matlab_datenum.flatten()[i]) for i in xrange(matlab_datenum.size)])
    return python_datetime

# <codecell>

def toutc(pydatetime):
    "Convert python datetime to utc fractional hours"
    utc_fx = lambda x: float(x.hour)+float(x.minute)/60.0+float(x.second)/3600.0+float(x.microsecond)/360000000.0
    try: 
        return utc_fx(pydatetime)
    except:
        import numpy as np
        return np.array([utc_fx(pydatetime.flatten()[i])+(pydatetime.flatten()[i].day-pydatetime.flatten()[0].day)*24.0 \
                         for i in xrange(pydatetime.size)]) 

# <codecell>

def load_emas(datfile):
    """
    Name:

        load_emas
    
    Purpose:

        to compile functions required to load emas files
        from within another script.
        Similar to load_modis
    
    Calling Sequence:

        emas,emas_dict = load_emas(datfile) 
    
    Input: 
  
        datfile name (hdf files)
    
    Output:

        emas dictionary with tau, ref, etau, eref, phase, qa
        emas_dicts : metadate for each of the variables
    
    Keywords: 

        none
    
    Dependencies:

        gdal
        numpy
        gc: for clearing the garbage
    
    Required files:
   
        dat files
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2014-12-08, NASA Ames
    """
    import numpy as np
    from osgeo import gdal
    from Sp_parameters import startprogress, progress, endprogress
    emas_values = (#('cloud_top',57),
                    ('phase',53),
          #          ('cloud_top_temp',58),
                    ('ref',66),
                    ('tau',72),
           #         ('cwp',82),
                    ('eref',90),
                    ('etau',93),
            #        ('ecwp',96),
                    ('multi_layer',105),
                    ('qa',123),
             #       ('cloud_mask',110)
                    )
    datsds = gdal.Open(datfile)
    datsub = datsds.GetSubDatasets()
    print 'Outputting the Data subdatasets:'
    for i in range(len(datsub)):
        if any(i in val for val in emas_values):
            print '\x1b[1;36m%i: %s\x1b[0m' %(i,datsub[i][1])
        else:
            print str(i)+': '+datsub[i][1]
    emas = dict()
    meta = datsds.GetMetadata() 
    import gc; gc.collect()
    emas_dicts = dict()
    startprogress('Running through modis values')
    for i,j in emas_values:
        sds = gdal.Open(datsub[j][0])
        emas_dicts[i] = sds.GetMetadata()
        emas[i] = np.array(sds.ReadAsArray())
        makenan = True
        bad_points = np.where(emas[i] == float(emas_dicts[i]['_FillValue']))
        try:
            scale = float(emas_dicts[i]['scale_factor'])
            offset = float(emas_dicts[i]['add_offset'])
            # print 'MODIS array: %s, type: %s' % (i, modis[i].dtype)
            if scale.is_integer():
               scale = int(scale)
               makenan = False
            if scale != 1 and offset == 0:
               emas[i] = emas[i]*scale+offset
        except:
            if issubclass(emas[i].dtype.type, np.integer):
                makenan = False
        if makenan:
            emas[i][bad_points] = np.nan
        progress(float(tuple(i[0] for i in emas_values).index(i))/len(emas_values)*100.)
    endprogress()
    print emas.keys()
    del datsds,sds,datsub
    return emas,emas_dicts

# <codecell>

def load_hdf(datfile,values=None,verbose=True):
    """
    Name:

        load_hdf
    
    Purpose:

        to compile functions required to load emas files
        from within another script.
        Similar to load_modis
    
    Calling Sequence:

        hdf_dat,hdf_dict = load_hdf(datfile,Values=None,verbose=True) 
    
    Input: 
  
        datfile name (hdf files)
    
    Output:

        hdf_dat dictionary with the names of values saved, with associated dictionary values
        hdf_dicts : metadate for each of the variables
    
    Keywords: 

        values: if ommitted, only outputs the names of the variables in file
                needs to be a tuple of 2 element tuples (first element name of variable, second number of record)
                example: modis_values=(('cloud_top',57),('phase',53),('cloud_top_temp',58),('ref',66),('tau',72))
        verbose: if true (default), then everything is printed. if false, nothing is printed
    
    Dependencies:

        gdal
        numpy
        gc: for clearing the garbage
        Sp_parameters for progress issuer
    
    Required files:
   
        dat files
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2014-12-10, NASA Ames
        Modified (v1.1): Samuel LeBlanc, 2015-04-10, NASA Ames
                        - added verbose keyword
    """
    import numpy as np
    from osgeo import gdal
    from Sp_parameters import startprogress, progress, endprogress
    
    datsds = gdal.Open(datfile)
    datsub = datsds.GetSubDatasets()
    if verbose: 
        print 'Outputting the Data subdatasets:'
        for i in range(len(datsub)):
            if values:
                if any(i in val for val in values):
                    print '\x1b[1;36m%i: %s\x1b[0m' %(i,datsub[i][1])
                else:
                    print str(i)+': '+datsub[i][1]
            else:
                print str(i)+': '+datsub[i][1]
    if not values:
        if verbose:
            print 'Done going through file... Please supply pairs of name and index for reading file'
            print " in format values = (('name1',index1),('name2',index2),('name3',index3),...)"
            print " where namei is the nameof the returned variable, and indexi is the index of the variable (from above)"
        return None, None
    hdf = dict()
    meta = datsds.GetMetadata() 
    import gc; gc.collect()
    hdf_dicts = dict()
    if verbose:
        startprogress('Running through data values')
    for i,j in values:
        sds = gdal.Open(datsub[j][0])
        hdf_dicts[i] = sds.GetMetadata()
        hdf[i] = np.array(sds.ReadAsArray())
        if not hdf[i].any():
            import pdb; pdb.set_trace()
        try:
            bad_points = np.where(hdf[i] == float(hdf_dicts[i]['_FillValue']))
            makenan = True
        except KeyError:
            makenan = False
        try:
            scale = float(hdf_dicts[i]['scale_factor'])
            offset = float(hdf_dicts[i]['add_offset'])
            # print 'MODIS array: %s, type: %s' % (i, modis[i].dtype)
            if scale.is_integer():
               scale = int(scale)
               makenan = False
            if scale != 1 and offset == 0:
               hdf[i] = hdf[i]*scale+offset
        except:
            if issubclass(hdf[i].dtype.type, np.integer):
                makenan = False
        if makenan:
            hdf[i][bad_points] = np.nan
        if verbose:
            progress(float(tuple(i[0] for i in values).index(i))/len(values)*100.)
    if verbose:
        endprogress()
        print hdf.keys()
    del datsds,sds,datsub
    return hdf,hdf_dicts

# <codecell>

def load_cpl_layers(datfile,values=None):
    """
    Name:

        load_cpl_layers
    
    Purpose:

        Function to load cpl files of layer properties 'layers_'
        This funciton is called from another script
    
    Calling Sequence:

        cpl_layers = load_cpl_layers(datfile) 
    
    Input: 
  
        datfile name (layers text files)
    
    Output:

        cpl_layers: dictionary with associated properties such as time, A/C altitude, latitude, longitude, and roll
                    number of layers, Ground height (GH) in meters above MSL, top and bottom altitude of each layer, and type of layer
    
    Keywords: 

        none
    
    Dependencies:

        os
        numpy
        re : for regular experessions
    
    Required files:
   
        layers text files
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-03-24, NASA Ames
    """
    import os
    if not(os.path.isfile(datfile)):
        error('File not found!')
    import numpy as np
    import re

    # set variables to read
    num_lines = sum(1 for line in open(datfile))
    head_lines = 14
    d = np.empty(num_lines-head_lines,dtype=[('hh','i4'),
                                             ('mm','i4'),
                                             ('ss','i4'),
                                             ('lat','f8'),
                                             ('lon','f8'),
                                             ('alt','f8'),
                                             ('rol','f8'),
                                             ('num','i4'),
                                             ('gh','f8'),
                                             ('top','f8',(10)),
                                             ('bot','f8',(10)),
                                             ('type','i4',(10)),
                                             ('utc','f8')])
    d[:] = np.NAN
    header = ''
    with open(datfile) as fp:
        for iline, line in enumerate(fp):
            if iline<head_lines:
                header = header + line
            else:
                i = iline-head_lines
                line = line.strip()
                temp = filter(None, re.split(r"[ :()]",line))
                d['hh'][i], d['mm'][i], d['ss'][i] = temp[0], temp[1], temp[2]
                d['lat'][i], d['lon'][i], d['alt'][i] = temp[3], temp[4], temp[5]
                d['rol'][i], d['num'][i], d['gh'][i] = temp[6], temp[7], temp[8]
                for n in range(d['num'][i]):
                    try:
                        d['top'][i][n], d['bot'][i][n], d['type'][i][n] = temp[9+3*n], temp[9+3*n+1], temp[9+3*n+2]
                    except:
                        import pdb; pdb.set_trace()
                d['utc'][i] = float(d['hh'][i])+float(d['mm'][i])/60.0+float(d['ss'][i])/3600.0
                
    return d

# <codecell>

def load_apr(datfiles):
    """
    Name:

        load_apr
    
    Purpose:

        Function to load apr values of zenith dbz from the various files in the datfiles list
    
    Calling Sequence:

        aprout = load_apr(datfiles) 
    
    Input: 
  
        datfiles name (list of .h4 files to combine)
    
    Output:

        aprout: dbz, zenith radar reflectivity
                latz, latitude of the zenith reflectivity
                lonz, longitude of the "
                altflt, actual altitude in the atmosphere of the radar refl.
                utc, time of measurement in utc fractional hours
                
    
    Keywords: 

        none
    
    Dependencies:

        os
        numpy
        load_modis
        datetime
    
    Required files:
   
        hdf APR-2 files from SEAC4RS
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-04-10, NASA Ames
    """
    import os
    import numpy as np
    from load_modis import load_hdf
    import datetime
    
    first = True
    for f in datfiles:
        print 'Running file: ',f
        if not(os.path.isfile(f)):
            print 'Problem with file:', f
            print ' ... Skipping'
            continue
        
        apr_value = (('lat',16),('lon',17),('alt',15),('time',13),('dbz',0),('lat3d',30),('lon3d',31),('alt3d',32),('lat3d_o',24),('lon3d_o',25),('alt3d_o',26),('lat3d_s',27),('lon3d_s',28),('alt3d_s',29))
        apr,aprdicts = load_hdf(f,values=apr_value,verbose=False)
        # transform the 3d latitudes, longitudes, and altitudes to usable values
        apr['latz'] = apr['lat3d']/apr['lat3d_s']+apr['lat3d_o']
        apr['lonz'] = apr['lon3d']/apr['lon3d_s']+apr['lon3d_o']
        apr['altz'] = apr['alt3d']/apr['alt3d_s']+apr['alt3d_o']
        apr['altflt'] = np.copy(apr['altz'])
        try:
            for z in range(apr['altz'].shape[0]):
                apr['altflt'][z,:,:] = apr['altz'][z,:,:]+apr['alt'][z,:]
        except:
            print 'Problem with file: ',f
            print ' ... dimensions do not agree'
            print ' ... Skipping'
            continue
        izen = apr['altz'][:,0,0].argmax() #get the index of zenith
        if first:
            aprout = dict()
            aprout['dbz'] = apr['dbz'][izen,:,:]
            aprout['altflt'] = apr['altz'][izen,:,:]+apr['alt'][izen,:]
            aprout['latz'] = apr['latz'][izen,:,:]
            aprout['lonz'] = apr['lonz'][izen,:,:]
            v = datetime.datetime.utcfromtimestamp(apr['time'][izen,0])
            aprout['utc'] = (apr['time'][izen,:]-(datetime.datetime(v.year,v.month,v.day,0,0,0)-datetime.datetime(1970,1,1)).total_seconds())/3600.
            first = False
        else:
            aprout['dbz'] = np.concatenate((aprout['dbz'].T,apr['dbz'][izen,:,:].T)).T
            aprout['altflt'] = np.concatenate((aprout['altflt'].T,(apr['altz'][izen,:,:]+apr['alt'][izen,:]).T)).T
            aprout['latz'] = np.concatenate((aprout['latz'].T,apr['latz'][izen,:,:].T)).T
            aprout['lonz'] = np.concatenate((aprout['lonz'].T,apr['lonz'][izen,:,:].T)).T
            v = datetime.datetime.utcfromtimestamp(apr['time'][izen,0])
            utc = (apr['time'][izen,:]-(datetime.datetime(v.year,v.month,v.day,0,0,0)-datetime.datetime(1970,1,1)).total_seconds())/3600.
            aprout['utc'] = np.concatenate((aprout['utc'].T,utc.T)).T
            
    print aprout.keys()
    print 'Loaded data points: ', aprout['utc'].shape
    return aprout        

# <codecell>

def load_amsr(datfile,lonlatfile):
    """
    Name:

        load_amsr
    
    Purpose:

        to load amsr data into a sucinct dictionary
    
    Calling Sequence:

        amsr = load_amsr(datfile,lonlatfile) 
    
    Input: 
  
        datfile: path and name of hdf file
        lonlatfile: path and name of hdf files for lat and lon
    
    Output:

        amsr: dictionary with numpy array of values
    
    Keywords: 

       none
    
    Dependencies:

        gdal
        numpy
        gc: for clearing the garbage
        Sp_parameters for progress issuer
        pdb: for debugging
        load_modis: this file
    
    Required files:
   
        dat file
        lonlat file
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-05-04, NASA Ames
        
    """
    import os
    if not(os.path.isfile(datfile)):
        error('Data file not found!')
    if not(os.path.isfile(lonlatfile)):
        error('Lonlat file not found!')
    import numpy as np
    from load_modis import load_hdf
    from osgeo import gdal
    gdat = gdal.Open(datfile)
    dat = dict()
    dat['nfo'] = gdat.GetMetadata()
    dat['ice'] = gdat.GetRasterBand(1).ReadAsArray()
    datll,dicll = load_hdf(lonlatfile,values=(('lon',0),('lat',1)),verbose=False)
    dat['lat'] = datll['lat']
    dat['lon'] = datll['lon']
    return dat

# <codecell>

def load_hdf_sd(FILE_NAME):
    """
    Name:

        load_hdf_sd
    
    Purpose:

        to load everyything in a hdf file using the SD protocol instead of GDAL
    
    Calling Sequence:

        dat,dat_dict = load_hdf_sd(FILE_NAME) 
    
    Input: 
  
        FILE_NAME: path and name of hdf file
    
    Output:

        dat: dictionary with numpy array of values
        dat_dict: dictionary with dictionaries of attributes for each read value
    
    Keywords: 

       none
    
    Dependencies:

        numpy
        pyhdf, SDC, SD
    
    Required files:
   
        dat file
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-05-13, NASA Ames
        
    """
    import numpy as np
    from pyhdf.SD import SD, SDC
    print 'Reading file: '+FILE_NAME
    hdf = SD(FILE_NAME, SDC.READ)
    dat = dict()
    dat_dict = dict()
    for name in hdf.datasets().keys():
        print '  '+name+': %s' % (hdf.datasets()[name],)
        dat[name] = hdf.select(name)[:]
        dat_dict[name] = hdf.select(name).attributes()
        try:
            missing_value = dat_dict[name]['missing_value']
            dat[name][dat[name] == missing_value] = np.nan
        except:
            print '  no missing value for key:'+name
    return dat, dat_dict

# <codecell>

def remove_field_name(a, name):
    names = list(a.dtype.names)
    if name in names:
        names.remove(name)
    b = a[names]
    return b

# <markdowncell>

# Testing of the script:

# <codecell>

if __name__ == "__main__":

# <codecell>

    import os
    import numpy as np
    from osgeo import gdal

# <codecell>

    fp='C:\\Users\\sleblan2\\Research\\TCAP\\'

# <codecell>

    myd06_file = fp+'MODIS\\MYD06_L2.A2013050.1725.006.2014260074007.hdf'
    myd03_file = fp+'MODIS\\MYD03.A2013050.1725.006.2013051163424.hdf'
    print os.path.isfile(myd03_file) #check if it exists
    print os.path.isfile(myd06_file)

# <codecell>

    myd_geo = gdal.Open(myd03_file)
    myd_geo_sub = myd_geo.GetSubDatasets()
    for i in range(len(myd_geo_sub)):
        print str(i)+': '+myd_geo_sub[i][1]

# <codecell>

    latsds = gdal.Open(myd_geo_sub[12][0],gdal.GA_ReadOnly)
    lonsds = gdal.Open(myd_geo_sub[13][0],gdal.GA_ReadOnly)
    szasds = gdal.Open(myd_geo_sub[21][0],gdal.GA_ReadOnly)

# <codecell>

    print latsds.RasterCount # verify that only one raster exists
    lat = latsds.ReadAsArray()
    lon = lonsds.ReadAsArray()
    sza = szasds.ReadAsArray()
    print lon.shape

# <markdowncell>

# Now load the specific data files:

# <codecell>

    myd_dat = gdal.Open(myd06_file)
    myd_dat_sub = myd_dat.GetSubDatasets()
    for i in range(len(myd_dat_sub)):
        print str(i)+': '+myd_dat_sub[i][1]

# <codecell>

    print myd_dat_sub[118]
    retfsds = gdal.Open(myd_dat_sub[118][0])

# <codecell>

    for key,value in myd_dat.GetMetadata_Dict().items():
        print key,value

# <markdowncell>

# Load the different modis values:

# <codecell>

    modis_values = (('cloud_top',57),
                    ('phase',53),
                    ('cloud_top_temp',58),
                    ('ref',66),
                    ('tau',72),
                    ('cwp',82),
                    ('eref',90),
                    ('etau',93),
                    ('ecwp',96),
                    ('multi_layer',105),
                    ('qa',123),
                    ('cloud_mask',110)
                    )

# <markdowncell>

# Testing the metadata dictionary

# <codecell>

    gdal.Open(myd_dat_sub[53][0]).GetMetadata()

# <codecell>

    mm = dict()
    mm['one'] = gdal.Open(myd_dat_sub[72][0]).GetMetadata()
    mm['two'] = gdal.Open(myd_dat_sub[74][0]).GetMetadata()
    mm['two']['_FillValue']

# <codecell>

    from Sp_parameters import startprogress, progress, endprogress
    import gc; gc.collect()

# <codecell>

    tuple(i[0] for i in modis_values).index('etau')

# <codecell>

    modis = dict()
    modis_dicts = dict()
    startprogress('Running through modis values')
    for i,j in modis_values:
        sds = gdal.Open(myd_dat_sub[j][0])
        modis_dicts[i] = sds.GetMetadata()
        modis[i] = np.array(sds.ReadAsArray())*float(modis_dicts[i]['scale_factor'])+float(modis_dicts[i]['add_offset'])
        modis[i][modis[i] == float(modis_dicts[i]['_FillValue'])] = np.nan
        progress(float(tuple(i[0] for i in modis_values).index(i))/len(modis_values)*100.)
    endprogress()

# <codecell>

    print modis.keys()
    print modis_dicts.keys()

