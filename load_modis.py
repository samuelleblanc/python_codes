# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def __init__():
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

def load_ict(fname):
    """
    Simple ict file loader
    created specifically to load the files from the iwg1 on oard the G1 during TCAP, may work with others...
    """
    from datetime import datetime
    import numpy as np
    f = open(fname,'r')
    first = f.readline()
    num2skip = int(first.strip().split(',')[0])
    f.close()
    def mktime(txt):
        return datetime.strptime(txt,'%Y-%m-%d %H:%M:%S')
    def utctime(seconds_utc):
        return float(seconds_utc)/3600.
    conv = {"Date_Time":mktime, "UTC":utctime, "Start_UTC":utctime}
    data = np.genfromtxt(fname,names=True,delimiter=',',skip_header=num2skip-1,dtype=None,converters=conv)
    print data.dtype.names
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

def load_hdf(datfile,values=None):
    """
    Name:

        load_hdf
    
    Purpose:

        to compile functions required to load emas files
        from within another script.
        Similar to load_modis
    
    Calling Sequence:

        hdf_dat,hdf_dict = load_hdf(datfile,Values=None) 
    
    Input: 
  
        datfile name (hdf files)
    
    Output:

        hdf_dat dictionary with the names of values saved, with associated dictionary values
        hdf_dicts : metadate for each of the variables
    
    Keywords: 

        values: if ommitted, only outputs the names of the variables in file
                needs to be a tuple of 2 element tuples (first element name of variable, second number of record)
                example: modis_values=(('cloud_top',57),('phase',53),('cloud_top_temp',58),('ref',66),('tau',72))
    
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
    """
    import numpy as np
    from osgeo import gdal
    from Sp_parameters import startprogress, progress, endprogress
    
    datsds = gdal.Open(datfile)
    datsub = datsds.GetSubDatasets()
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
        print 'Done going through file... Please supply pairs of name and index for reading file'
        print " in format values = (('name1',index1),('name2',index2),('name3',index3),...)"
        print " where namei is the nameof the returned variable, and indexi is the index of the variable (from above)"
        return None, None
    hdf = dict()
    meta = datsds.GetMetadata() 
    import gc; gc.collect()
    hdf_dicts = dict()
    startprogress('Running through modis values')
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
        progress(float(tuple(i[0] for i in values).index(i))/len(values)*100.)
    endprogress()
    print hdf.keys()
    del datsds,sds,datsub
    return hdf,hdf_dicts

# <codecell>

def spherical_dist(pos1, pos2, r=3958.75):
    "Calculate the distance, in km, from one point to another (can use arrays)"
    import numpy as np
    pos1 = pos1 * np.pi / 180
    pos2 = pos2 * np.pi / 180
    cos_lat1 = np.cos(pos1[..., 0])
    cos_lat2 = np.cos(pos2[..., 0])
    cos_lat_d = np.cos(pos1[..., 0] - pos2[..., 0])
    cos_lon_d = np.cos(pos1[..., 1] - pos2[..., 1])
    return r * np.arccos(cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))

# <codecell>

def map_ind(mod_lon,mod_lat,meas_lon,meas_lat,meas_good=None):
    """ Run to get indices in the measurement space of all the closest mod points. Assuming earth geometry."""
    from load_modis import spherical_dist
    from Sp_parameters import startprogress, progress, endprogress
    import numpy as np
    if not any(meas_good):
        meas_good = np.where(meas_lon)
    imodis = np.logical_and(np.logical_and(mod_lon>min(meas_lon[meas_good])-0.02 , mod_lon<max(meas_lon[meas_good])+0.02),
                            np.logical_and(mod_lat>min(meas_lat[meas_good])-0.02 , mod_lat<max(meas_lat[meas_good])+0.02))
    wimodis = np.where(imodis)
    N1 = mod_lon[imodis].size
    modis_grid = np.hstack([mod_lon[imodis].reshape((N1,1)),mod_lat[imodis].reshape((N1,1))])
    N2 = len(meas_good)
    meas_grid = np.hstack([np.array(meas_lon[meas_good]).reshape((N2,1)),np.array(meas_lat[meas_good]).reshape((N2,1))])
    meas_in = meas_grid.astype(int)
    meas_ind = np.array([meas_good.ravel()*0,meas_good.ravel()*0])
    startprogress('Running through flight track')
    for i in xrange(meas_good.size):
        d = spherical_dist(meas_grid[i],modis_grid)
        meas_ind[0,i] = wimodis[0][np.argmin(d)]
        meas_ind[1,i] = wimodis[1][np.argmin(d)]
        progress(float(i)/len(meas_good)*100)
    endprogress()
    return meas_ind

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

