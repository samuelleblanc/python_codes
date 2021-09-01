#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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


# In[1]:


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
                    ('cth',183)
             #       ('cloud_mask',110)
                    )
    geosds = gdal.Open(geofile)
    datsds = gdal.Open(datfile)
    geosub = geosds.GetSubDatasets()
    datsub = datsds.GetSubDatasets()
    print('Outputting the Geo subdatasets:')
    for i in range(len(geosub)):
        print(str(i)+': '+geosub[i][1])
    print('Outputting the Data subdatasets:')
    for i in range(len(datsub)):
        if any(i in val for val in modis_values):
            print('\x1b[1;36m%i: %s\x1b[0m' %(i,datsub[i][1]))
        else:
            print(str(i)+': '+datsub[i][1])
    latsds = gdal.Open(geosub[12][0],gdal.GA_ReadOnly)
    lonsds = gdal.Open(geosub[13][0],gdal.GA_ReadOnly)
    szasds = gdal.Open(geosub[21][0],gdal.GA_ReadOnly)
    modis = dict()
    modis['lat'] = latsds.ReadAsArray()
    modis['lon'] = lonsds.ReadAsArray()
    modis['sza'] = szasds.ReadAsArray()
    print(modis['lon'].shape)
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
    print(list(modis.keys()))
    del geosds, datsds,sds,lonsds,latsds,geosub,datsub
    return modis,modis_dicts


# In[1]:


def load_ict(fname,return_header=False,make_nan=True):
    """
    Simple ict file loader
    created specifically to load the files from the iwg1 on board the G1 during TCAP, may work with others...
    inputs:
       fname: filename with full path
       return_header: (default set to False) if True, returns data, header in that form
       make_nan: (default set to True) if True, the values defined in the header to be missing data, usually -999, is changed to NaNs
    """
    from datetime import datetime
    import numpy as np
    f = open(fname,'r')
    lines = f.readlines()
    first = lines[0]
    sep = ','
    try:
        num2skip = int(first.strip().split(sep)[0])
    except ValueError:
        print('Seperation is set to a space')
        sep = None
        num2skip = int(first.strip().split(sep)[0])
    header = lines[0:num2skip]
    factor = list(map(float,header[10].strip().split(sep)))
    missing = list(map(float,header[11].strip().split(sep)))
    f.close()
    if any([i!=1 for i in factor]):
        print('Some Scaling factors are not equal to one, Please check the factors:')
        print(factor)
    def mktime(txt):
        return datetime.strptime(txt,'%Y-%m-%d %H:%M:%S')
    def utctime(seconds_utc):
        return float(seconds_utc)/3600.
    conv = {"Date_Time":mktime, "UTC":utctime, "Start_UTC":utctime, "TIME_UTC":utctime, "UTC_mid":utctime}
    data = np.genfromtxt(fname,names=True,delimiter=sep,skip_header=num2skip-1,converters=conv)
    print(data.dtype.names)
    #scale the values by using the scale factors
    for i,name in enumerate(data.dtype.names):
        if i>0:
            if factor[i-1]!=float(1):
                data[name] = data[name]*factor[i-1]
            if make_nan:
                data[name][data[name]==missing[i-1]]=np.NaN
    if return_header:
        return data, header
    else:
        return data


# In[3]:


def modis_qa_MOD06(qa_array):
    """
    Purpose:
    
        modis qa data parser for Cloud properties.

            Set for useful COD and REF
            Set for 2nd highest (not highest) confidence level
            
    Input: 
    
        input of qa numpy array from MODIS MOD06/MYD06 hdf files
        qa_array: 3d array of the QA values 'Quality_Assurance_1km' from MOD06 hdf, read by pyhdf.SD
    
    Output:
    
        output qa numpy boolean array (True is high QA)
    
    Dependencies:

        load_utils.bits_stripping
    
    Required files:
   
        none
    
    Example:

        > from pyhdf.SD import SD, SDC
        > fh = SD(str(fp+'data_other/MODIS/MOD06_L2/2016/264/MOD06_L2.A2016264.1030.061.2017328080629.hdf'),SDC.READ)
        > co = fh.select('Cloud_Optical_Thickness')
        > codm = co.get()
        > codm = codm*co.attributes()['scale_factor']
        > codm[codm<0] = np.nan
        > qal = fh.select('Quality_Assurance_1km')
        > qa = modis_qa_MOD06(qal.get())
        > codm[~qa] = np.nan
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2020-11-23, Santa Cruz, CA, ported from earlier implementation in ORACLES_build_DARE
        
    """
    #bin8 = lambda x : ''.join(reversed( [str((x >> i) & 1) for i in range(8)] ) )
    from load_utils import bits_stripping
    
    # for clouds
    qa_cod_useful = bits_stripping(0,1,qa_array[:,:,4])
    qa_cod_conf = bits_stripping(1,2,qa_array[:,:,4])
    qa_ref_useful = bits_stripping(3,1,qa_array[:,:,4])
    qa_ref_conf = bits_stripping(4,2,qa_array[:,:,4])
    qa = (qa_cod_useful>0) & (qa_cod_conf>1) & (qa_ref_useful>0) & (qa_ref_conf>1) 
    
    return qa


# In[4]:


def bits_stripping(bit_start,bit_count,value):
    "Support function to the modis_qa flags (MOD06) to parse out a bit array from hdf files"
    import numpy as np
    bitmask=pow(2,bit_start+bit_count)-1
    return np.right_shift(np.bitwise_and(value,bitmask),bit_start)


# In[ ]:


def read_mod06(fh,set_qa_nan=False):
    """
    Purpose:
    
        Simple modis reader function, with input 'fh' is a SD object. 
        Part of broader analysis which preloads the hdf files for faster parallel analysis.
        
    Input: 
    
        fh: file SD handle for Hdf MOD06/MYD06 files
        set_qa_nan: if set (default False) will nan ou bad cloud QA values
    
    Output:
    
        codm : Cloud optical thickness at 1km  
        ref : cloud effective radius at 1km 
        
    Dependencies:

        load_utils.modis_qa_MOD06
    
    Required files:
   
        MODIS hdf files
    
    Example:

        > from pyhdf.SD import SD, SDC
        > geosy[ad]['fh'] = [SD(str(fp+'data_other/MODIS/MYD06_L2/2016/264/MYD06_L2.A2016264.1315.061.2018062110512.hdf'),SDC.READ)]
        > for igy,gy in enumerate(geosy[ad]['fh']):
                ind = mu.map_ind(geosy[ad]['lat'][igy,:,:],geosy[ad]['lon'][igy,:,:],
                         ar['Latitude'][fla][iaes][idd],ar['Longitude'][fla][iaes][idd])
        >  if np.array(ind).any():
                cody,refy = read_mod06(gy)
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2020-11-23, Santa Cruz, CA, ported from earlier implementation in ORACLES_build_DARE
    
    """
    from load_utils import modis_qa_MOD06
    re = fh.select('Cloud_Effective_Radius')
    ref = re.get()
    ref = ref*re.attributes()['scale_factor']
    ref[ref<0] = np.nan
    
    co = fh.select('Cloud_Optical_Thickness')
    codm = co.get()
    codm = codm*co.attributes()['scale_factor']
    codm[codm<0] = np.nan
        
    qal = fh.select('Quality_Assurance_1km')
    qa = modis_qa_MOD06(qal.get())
    
    if set_qa_nan:
        codm[~qa] = np.nan
        ref[~qa] = np.nan
    
    return codm,ref


# In[ ]:


def mat2py_time(matlab_datenum):
    "convert a matlab datenum to a python datetime object. Works on numpy arrays of datenum"
    from datetime import datetime, timedelta
    #matlab_datenum = 731965.04835648148
    try: #for python 3 compatibility
        xrange
    except NameError:
        xrange = range
    m2ptime = lambda tmat: datetime.fromordinal(int(tmat)) + timedelta(days=tmat%1) - timedelta(days = 366)
    try:
        python_datetime = m2ptim2(matlab_datenum)
    except:
        import numpy as np
        python_datetime = np.array([m2ptime(matlab_datenum.flatten()[i]) for i in xrange(matlab_datenum.size)])
    return python_datetime


# In[1]:


def toutc(pydatetime):
    "Convert python datetime to utc fractional hours"
    utc_fx = lambda x: float(x.hour)+float(x.minute)/60.0+float(x.second)/3600.0+float(x.microsecond)/3600000000.0
    try: #for python 3 compatibility
        xrange
    except NameError:
        xrange = range
    try: 
        return utc_fx(pydatetime)
    except:
        import numpy as np
        try:
            return np.array([utc_fx(pydatetime.flatten()[i])+(pydatetime.flatten()[i].day-pydatetime.flatten()[0].day)*24.0                          for i in xrange(pydatetime.size)]) 
        except AttributeError:
            return np.array([utc_fx(pydatetime[i])+(pydatetime[i].day-pydatetime[0].day)*24.0                          for i in xrange(len(pydatetime))]) 


# In[ ]:


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
    print('Outputting the Data subdatasets:')
    for i in range(len(datsub)):
        if any(i in val for val in emas_values):
            print('\x1b[1;36m%i: %s\x1b[0m' %(i,datsub[i][1]))
        else:
            print(str(i)+': '+datsub[i][1])
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
    print(list(emas.keys()))
    del datsds,sds,datsub
    return emas,emas_dicts


# In[ ]:


def load_hdf(datfile,values=None,verbose=True,all_values=False,i_subdata=1):
    """
    Name:

        load_hdf
    
    Purpose:

        to compile functions required to load emas files
        from within another script.
        Similar to load_modis
    
    Calling Sequence:

        hdf_dat,hdf_dict = load_hdf(datfile,Values=None,verbose=True,all_values=False) 
    
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
        all_values: if True, then outputs all the values with their original names, (defaults to False), overrides values keyword
        i_subdata: value of the subdataset index
    
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
        Modified (v1.2): Samuel LeBlanc, 2016-05-07, Osan AFB, Korea
                        - added error handling for missing fill value
        Modified (v1.3): Samuel LeBlanc, 2016-11-15, NASA Ames
                        - added all_values keyword
        Modified (v1.4): Samuel LeBlanc, 2021-05-28, Santa Cruz, CA
                        - added the i_subdata keyword
    """
    import numpy as np
    from osgeo import gdal
    from Sp_parameters import startprogress, progress, endprogress
    from load_utils import load_hdf_h5py
    
    sds_level = 1
    
    try:
        datsds = gdal.Open(datfile)
    except:
        print('Trying with load_hdf_h5py')
        return load_hdf_h5py(datfile,values=values,verbose=verbose,all_values=all_values)
    datsub = datsds.GetSubDatasets()
    if verbose: 
        print('Outputting the Data subdatasets:')
        for i in range(len(datsub)):
            if values:
                if any(i in val for val in values):
                    print('\x1b[1;36m%i: %s\x1b[0m' %(i,datsub[i][sds_level]))
                else:
                    print(str(i)+': '+datsub[i][sds_level])
            else:
                print(str(i)+': '+datsub[i][sds_level])
    if all_values:
        values = []
        for i in range(len(datsub)):
            values.append((datsub[i][sds_level].split(' ')[1],i))
        values = tuple(values)
    if not values:
        if verbose:
            print('Done going through file... Please supply pairs of name and index for reading file')
            print(" in format values = (('name1',index1),('name2',index2),('name3',index3),...)")
            print(" where namei is the nameof the returned variable, and indexi is the index of the variable (from above)")
        return None, None
    hdf = dict()
    meta = datsds.GetMetadata() 
    import gc; gc.collect()
    hdf_dicts = dict()
    if verbose:
        startprogress('Running through data values')
    for i,j in values:
        sds = gdal.Open(datsub[j][i_subdata])
        hdf_dicts[i] = sds.GetMetadata()
        hdf[i] = np.array(sds.ReadAsArray())
        if not hdf[i].any():
            import pdb; pdb.set_trace()
        try:
            bad_points = np.where(hdf[i] == float(hdf_dicts[i]['_FillValue']))
            makenan = True
        except KeyError:
            makenan = False
        except ValueError:
            makenan = False
            print('*** FillValue not used to replace NANs, will have to do manually ***')
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
            try:
                hdf[i][bad_points] = np.nan
            except ValueError:
                print('*** Can not replace NaN into variable: {}, for the bad points {} ***'.format(i,bad_points))
                
        if verbose:
            progress(float(tuple(i[0] for i in values).index(i))/len(values)*100.)
    if verbose:
        endprogress()
        print(list(hdf.keys()))
    del datsds,sds,datsub
    return hdf,hdf_dicts


# In[ ]:


def load_hdf_h5py(datfile,values=None,verbose=True,all_values=False):
    """
    Name:

        load_hdf_h5py
    
    Purpose:

        to compile functions required to load hdf5 files, through h5py python modules. 
        from within another script.
        Similar to load_modis
    
    Calling Sequence:

        hdf_dat,hdf_dict = load_hdf_h5py(datfile,Values=None,verbose=True,all_values=False) 
    
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
        all_values: if True, then outputs all the values with their original names, (defaults to False), overrides values keyword
    
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
    
        Written (v1.0): Samuel LeBlanc, 2020-06-26, Santa Cruz, CA, COVID-19 Quarantining
        Modified (v1.1): 
    """
    import numpy as np
    import h5py
    from Sp_parameters import startprogress, progress, endprogress
    
    dat = h5py.File(datfile,'r')
    if verbose: print('Reading file: {}'.format(datfile))
    datsub = list(dat.items())
    if verbose: 
        print('Outputting the Data subdatasets:')
        for i in range(len(datsub)):
            if values:
                if any(i in val for val in values):
                    print('\x1b[1;36m%i: %s\x1b[0m' %(i,datsub[i][0]))
                else:
                    print(str(i)+': '+datsub[i][0])
            else:
                print(str(i)+': '+datsub[i][0])
    if all_values:
        values = []
        for i in range(len(datsub)):
            values.append((datsub[i][0],i))
        values = tuple(values)
    if not values:
        if verbose:
            print('Done going through file... Please supply pairs of name and index for reading file')
            print(" in format values = (('name1',index1),('name2',index2),('name3',index3),...)")
            print(" where namei is the nameof the returned variable, and indexi is the index of the variable (from above)")
        return None, None
    hdf = dict()
    meta = dict(dat.attrs)
    import gc; gc.collect()
    hdf_dicts = dict()
    if verbose:
        startprogress('Running through data values')
    for i,j in values:
        sds = dat[i]
        hdf_dicts[i] = dict(sds.attrs)
        try:
            hdf[i] = np.array(dat[i].value)
        except TypeError:
            import pdb; pdb.set_trace()
        try:
            if not hdf[i].any():
                if len(hdf[i])<1:
                    import pdb; pdb.set_trace()
        except TypeError:
            chararray = [''.join(sdf) for sdf in hdf[i]]
            hdf[i] = np.array(chararray)
        try:
            bad_points = np.where(hdf[i] == sds.fillvalue)
            makenan = True
        except KeyError:
            makenan = False
        except ValueError:
            makenan = False
            print('*** FillValue not used to replace NANs, will have to do manually ***')
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
        print(list(hdf.keys()))
    hdf_dicts['Global'] = meta
    del sds,datsub
    return hdf,hdf_dicts


# In[1]:


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
                temp = [_f for _f in re.split(r"[ :()]",line) if _f]
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


# In[1]:


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
        Modified (v1.1): Samuel LeBlanc, 2016-11-15, NASA Ames
                         - added Ka band (35 Ghz) reading from the file
                         - and zenith dbz 
    """
    import os
    import numpy as np
    from load_utils import load_hdf
    import datetime
    
    first = True
    for f in datfiles:
        print('Running file: ',f)
        if not(os.path.isfile(f)):
            print('Problem with file:', f)
            print(' ... Skipping')
            continue
        
        apr_value = (('lat',16),('lon',17),('alt',15),('time',13),('dbz',0),('lat3d',30),
                     ('lon3d',31),('alt3d',32),('lat3d_o',24),('lon3d_o',25),('alt3d_o',26),
                     ('lat3d_s',27),('lon3d_s',28),('alt3d_s',29),('dbz_35',3),('zen_dbz',4),
                     ('alt3dz',41),('alt3dz_o',35),('alt3dz_s',38))
        apr,aprdicts = load_hdf(f,values=apr_value,verbose=False)
        # transform the 3d latitudes, longitudes, and altitudes to usable values
        apr['latz'] = apr['lat3d']/apr['lat3d_s']+apr['lat3d_o']
        apr['lonz'] = apr['lon3d']/apr['lon3d_s']+apr['lon3d_o']
        apr['altz'] = apr['alt3d']/apr['alt3d_s']+apr['alt3d_o']
        apr['altzen'] = apr['alt3dz']/apr['alt3dz_s']+apr['alt3dz_o']
        apr['altflt'] = np.copy(apr['altz'])
        apr['altfltz'] = np.copy(apr['altzen'])
        try:
            for z in range(apr['altz'].shape[0]):
                apr['altflt'][z,:,:] = apr['altz'][z,:,:]+apr['alt'][z,:]
            for z in range(apr['altzen'].shape[0]):
                apr['altfltz'][z,:,:] = apr['altzen'][z,:,:]+apr['alt'][z,:]
        except IndexError:
            try:
                print('swaping axes')
                apr['altflt'] = np.swapaxes(apr['altflt'],0,1)
                apr['altz'] = np.swapaxes(apr['altz'],0,1)
                apr['altfltz'] = np.swapaxes(apr['altfltz'],0,1)
                apr['altzen'] = np.swapaxes(apr['altzen'],0,1)
                apr['latz'] = np.swapaxes(apr['latz'],0,1)
                apr['lonz'] = np.swapaxes(apr['lonz'],0,1)
                apr['dbz'] = np.swapaxes(apr['dbz'],0,1)
                apr['dbz_35'] = np.swapaxes(apr['dbz_35'],0,1)
                apr['zen_dbz'] = np.swapaxes(apr['zen_dbz'],0,1)
                for z in range(apr['altz'].shape[0]):
                    apr['altflt'][z,:,:] = apr['altz'][z,:,:]+apr['alt'][z,:]
                for z in range(apr['altzen'].shape[0]):
                    apr['altfltz'][z,:,:] = apr['altzen'][z,:,:]+apr['alt'][z,:]
            except:
                print('Problem file:',f)
                print('... Skipping')
                continue
        except:
            print('Problem with file: ',f)
            print(' ... Skipping')
            continue
        izen = apr['altz'][:,0,0].argmax() #get the index of zenith
        izenz = apr['altzen'][:,0,0].argmax() #get the index of zenith
        if first:
            aprout = dict()
            aprout['dbz'] = apr['dbz'][izen,:,:]
            aprout['dbz_35'] = apr['dbz_35'][izen,:,:]
            aprout['zen_dbz'] = apr['zen_dbz'][izenz,:,:]
            aprout['altflt'] = apr['altz'][izen,:,:]+apr['alt'][izen,:]
            aprout['altfltz'] = apr['altzen'][izenz,:,:]+apr['alt'][izen,:]
            aprout['latz'] = apr['latz'][izen,:,:]
            aprout['lonz'] = apr['lonz'][izen,:,:]
            v = datetime.datetime.utcfromtimestamp(apr['time'][izen,0])
            aprout['utc'] = (apr['time'][izen,:]-(datetime.datetime(v.year,v.month,v.day,0,0,0)-datetime.datetime(1970,1,1)).total_seconds())/3600.
            first = False
        else:
            aprout['dbz'] = np.concatenate((aprout['dbz'].T,apr['dbz'][izen,:,:].T)).T
            aprout['dbz_35'] = np.concatenate((aprout['dbz_35'].T,apr['dbz'][izen,:,:].T)).T
            aprout['zen_dbz'] = np.concatenate((aprout['zen_dbz'].T,apr['zen_dbz'][izenz,:,:].T)).T
            aprout['altflt'] = np.concatenate((aprout['altflt'].T,(apr['altz'][izen,:,:]+apr['alt'][izen,:]).T)).T
            aprout['altfltz'] = np.concatenate((aprout['altfltz'].T,(apr['altzen'][izenz,:,:]+apr['alt'][izen,:]).T)).T
            aprout['latz'] = np.concatenate((aprout['latz'].T,apr['latz'][izen,:,:].T)).T
            aprout['lonz'] = np.concatenate((aprout['lonz'].T,apr['lonz'][izen,:,:].T)).T
            v = datetime.datetime.utcfromtimestamp(apr['time'][izen,0])
            utc = (apr['time'][izen,:]-(datetime.datetime(v.year,v.month,v.day,0,0,0)-datetime.datetime(1970,1,1)).total_seconds())/3600.
            aprout['utc'] = np.concatenate((aprout['utc'].T,utc.T)).T
            
    print(list(aprout.keys()))
    print('Loaded data points: ', aprout['utc'].shape)
    return aprout


# In[1]:


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
    from load_utils import load_hdf
    from osgeo import gdal
    gdat = gdal.Open(datfile)
    dat = dict()
    dat['nfo'] = gdat.GetMetadata()
    dat['ice'] = gdat.GetRasterBand(1).ReadAsArray()
    datll,dicll = load_hdf(lonlatfile,values=(('lon',0),('lat',1)),verbose=False)
    dat['lat'] = datll['lat']
    dat['lon'] = datll['lon']
    return dat


# In[ ]:


def load_hdf_sd(FILE_NAME,vals=[],verbose=True):
    """
    Name:

        load_hdf_sd
    
    Purpose:

        to load everyything in a hdf file using the SD protocol instead of GDAL
        makes nans out of Missing_value and _FilValue. Scales the values by the scale_factor
    
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
        Modified (v1.1): by Samuel LeBlanc, 2015-07-01, NASA Ames, Happy Canada Day!
                        - added Fill value keyword selection
                        - added scale factor and add offset
        Modified (v1.2): Samuel LeBlanc, 2021-05-28, Santa Cruz, CA
                        - added the vals keyword for only reading certain variables
                        - added verbose keyword
        
    """
    import numpy as np
    from pyhdf.SD import SD, SDC
    print('Reading file: '+FILE_NAME)
    hdf = SD(str(FILE_NAME), SDC.READ)
    dat = dict()
    dat_dict = dict()
    for name in list(hdf.datasets().keys()):
        if not vals or name in vals:
            if verbose: print('  '+name+': %s' % (hdf.datasets()[name],))
            dat[name] = hdf.select(name)[:]
            dat_dict[name] = hdf.select(name).attributes()
            try:
                scale_factor = dat_dict[name].get('scale_factor')
                if not scale_factor:
                    scale_factor = 1.0
                dat[name] = dat[name]*scale_factor
                try:
                    dat[name][dat[name] == dat_dict[name].get('missing_value')*scale_factor] = np.nan
                except TypeError:
                    if verbose: print('No missing_value on '+name)
                try:
                    dat[name][dat[name] == dat_dict[name].get('_FillValue')*scale_factor] = np.nan
                except TypeError:
                    if verbose: print('No FillValue on '+name)
                add_offset = dat_dict[name].get('add_offset')
                if not add_offset:
                    add_offset = 0
                dat[name] = dat[name] + add_offset
            except:
                if verbose: print('Problem in filling with nans and getting the offsets, must do it manually')
    return dat, dat_dict


# In[ ]:


def remove_field_name(a, name):
    names = list(a.dtype.names)
    if name in names:
        names.remove(name)
    b = a[names]
    return b


# In[61]:


def load_netcdf(datfile,values=None,verbose=True,everything=False):
    """
    Name:

        load_netcdf
    
    Purpose:

        To compile the functions required to load a netcdf4 file in a similar manner as hdf files
    
    Calling Sequence:

        cdf_data,cdf_dict = load_netcdf(datfile,Values=None,verbose=True) 
    
    Input: 
  
        datfile name (netcdf files)
    
    Output:

        cdf_data: dictionary with the names of values saved, with associated dictionary values
        cdf_dict: metadate for each of the variables
    
    Keywords: 

        values: if ommitted, only outputs the names of the variables in file
                needs to be a tuple of 2 element tuples (first: name of variable to be outputted,
                second: indes of full name in variables)
                example: modis_values=(('tau',35),('lat',22),('lon',23))
        verbose: if true (default), then everything is printed. if false, nothing is printed
        everything: if true (defaults to false), then everything is saved, without the need for the 'values' variable
    
    Dependencies:

        netcdf4
        
    Required files:
   
        dat files
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-12-07, NASA Ames
        Modified (v1.1): Samuel  LeBlanc, 2018-08-09, NASA Ames Research Center
                        Added the 'everything' keyword for loading all the variables from the netcdf file.
        
    """
    import netCDF4 as nc
    if verbose:
        print('Reading file: '+datfile)
    f = nc.Dataset(datfile,'r')
    varnames = list(f.variables.keys())
    if everything: values = [('nul',-999)]
    if verbose: 
        print('Outputting the Data subdatasets:')
        for i in range(len(varnames)):
            if values:
                if everything: values.append((varnames[i].encode('ascii','ignore'),i))
                if any(i in val for val in values):
                    print('\x1b[1;36m{0}: {1}\x1b[0m'.format(i,varnames[i]))
                else:
                    print('{0}: {1}'.format(i,varnames[i]))
            else:
                print('{0}: {1}'.format(i,varnames[i]))
    if not values:
        if verbose:
            print('Done going through file... Please supply pairs of name and index for reading file')
            print(" in format values = (('name1',index1),('name2',index2),('name3',index3),...)")
            print(" where namei is the name of the returned variable, and indexi is the index of the variable (from above)")
        return None, None
    if everything: values = values[1:]
    cdf_dict = {}
    cdf_data = {}
    for i,j in values:
        cdf_dict[i] = f.variables[varnames[j]]
        cdf_data[i] = f.variables[varnames[j]][:]
    if verbose:
        print(list(cdf_dict.keys()))
    return cdf_data,cdf_dict


# In[ ]:


def load_aeronet(f,verbose=True,version=2):
    """
    Name:

        load_aeronet
    
    Purpose:

        To load the LEV1.0 Aeronet AOD files
    
    Calling Sequence:

        aeronet = load_aeronet(f) 
    
    Input: 
  
        f: path and name of lev10 file
    
    Output:

        aeronet: numpy recarray of values
    
    Keywords: 

       verbose: (default True) if True, then prints out info as data is read
       version: (default 2) version of the AERONET file (version 3 has a different amount of header)
    
    Dependencies:

        numpy
        os
        load_modis: this file
    
    Required files:
   
        LEV10 file
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2016-05-09, Osan AB, Korea
        
    """
    import numpy as np
    from datetime import datetime
    import load_utils as lm
    import os
    if not(os.path.isfile(f)):
        raise IOError('Data file {} not found!'.format(f))
    if f.split('.')[-1].find('lev10')<0 | f.split('.')[-1].find('LEV10')<0:
        raise IOError('Data file {} is not a level 1.0 file - it is not yet available to read'.format(f))
    def makeday(txt):
        return datetime.strptime(txt,'%d:%m:%Y').timetuple().tm_yday
    def maketime(txt):
        return lm.toutc(datetime.strptime(txt,'%H:%M:%S'))
    conv = {'Dateddmmyy':makeday,'Timehhmmss':maketime}
    if verbose:
        print('Opening file: {}'.format(f))
    if version==3:
        header_lines = 6
    else:
        header_lines = 4
    ra = np.genfromtxt(f,skip_header=header_lines,names=True,delimiter=',',converters=conv)
    da = lm.recarray_to_dict(ra)
    ff = open(f,'r')
    lines = ff.readlines()
    da['header'] = lines[0:header_lines]
    if version==2:
        for n in da['header'][2].split(','):
            u = n.split('=')
            try:
                da[u[0]]=float(u[1])
            except:
                da[u[0]] = u[1].strip()
    return da


# In[ ]:


def load_multi_aeronet(dir_path,verbose=True):
    """
    Name:

        load_multi_aeronet
    
    Purpose:

        To load multiple files of the LEV1.0 Aeronet AOD files
    
    Calling Sequence:

        aeronet = load_multi_aeronet(dir_path) 
    
    Input: 
  
        dir_path: path of directory where multiple lev10 file reside
    
    Output:

        aeronet: numpy recarray of combined values from multiple files
    
    Keywords: 

       verbose: (default True) if True, then prints out info as data is read
    
    Dependencies:

        numpy
        os
        load_utils: this file
    
    Required files:
   
        LEV10 file
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2016-05-12, Osan AB, Korea
        
    """
    import os
    import numpy as np
    import load_utils as lm
    f = os.listdir(dir_path)
    aero = []
    for fl in f:
        aero.append(lm.load_aeronet(dir_path+fl,verbose=verbose))
    
    n_max = 0
    for n in range(len(aero)):
        if len(aero[n]['AOT_500'])>n_max: 
            n_max = len(aero[n]['AOT_500'])
            
    anet = {}
    nstations = len(aero)
    for name in list(aero[0].keys()):
        if not isinstance(aero[0][name],np.ndarray):
            anet[name] = []
            for n in range(nstations):
                anet[name].append(aero[n][name])
        else:
            anet[name] = np.zeros((nstations,n_max))+np.nan
            for n in range(nstations):
                anet[name][n,0:len(aero[n][name])] = aero[n][name]
    return anet


# In[ ]:


def aeronet_subset(aero,doy=None,utc=None,julian=None,window=24.0):
    """
    Name:

        aeronet_subset
    
    Purpose:

        Subsets the aero dict created from load_multi_aeronet for returning 
        only a index value linking only one point per aeronet station.
    
    Calling Sequence:

        ii = aeronet_subset(aero,doy=doy,utc=utc,julian=julian,window=24.0) 
    
    Input: 
  
        aero: aeronet dict of compiled aeronet files from many locations
        doy: (default None) put in the day of year value for the time to be selected
        utc: (default None) the utc in fractional hours, for the returned aeronet values
        julian: (default None) the fractional day of year for the returned values
        window: (default 24) the hours of the window, if no value is found within 
                this window, then returned index links to a nan value
        
        if there is no doy, utc julian, or window value, 
        will return the latest value in the dict
    
    Output:

        ii: the index values linking to the searched aero times. 
            it returns a tuple of indexes to be used directly into the aero dict.
    
    Keywords: 

       see above.
    
    Dependencies:

        numpy
    
    Required files:
   
        none
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2016-05-12, Osan AB, Korea
        
    """
    import numpy as np
    
    latest = False
    if julian:
        in_julian = julian
    elif doy:
        if utc:
            in_julian = doy+utc/24.0
        else:
            in_julian = doy+0.0
    else:
        latest = True
               
    if latest:
        ilatest = []
        for i,n in enumerate(aero['Location']):
            ilatest.append(np.nanargmax(aero['Julian_Day'][i,:]))
    else:
        ilatest = []
        for i,n in enumerate(aero['Location']):
            ia = np.nanargmin(abs(aero['Julian_Day'][i,:]-in_julian))
            if np.nanmin(abs(aero['Julian_Day'][i,:]-in_julian))*24.0 < window:
                ilatest.append(ia)
            else:
                ia = len(aero['Julian_Day'][i,:])-1
                ilatest.append(ia)
    iil = []
    for i,n in enumerate(ilatest):
        iil.append(i)
    ii = (iil,ilatest)
    
    return ii


# In[ ]:


def recarray_to_dict(ra):
    'simple function to convert numpy recarray to a dict with numpy arrays. Useful for modifying the output from genfromtxt'
    import numpy as np
    da = {}
    for n in ra.dtype.names:
        da[n]=ra[n]
    return da


# In[ ]:


def save_to_json(f,d):
    'Function to save dictionary with numpy elements (d) to a text file (f) define by the JSON typing'
    from json_tricks.np import dumps
    with open(f,'w') as handle:
        handle.write(dumps(d,indent=2))


# In[ ]:


def load_from_json(f):
    'Function to load JSON file and translate to dictionary with numpy elements from text file(f) define by the JSON typing'
    from json_tricks.np import loads
    with open(f,'r') as handle:
        d = dict(loads(handle.read()))
    return d


# In[ ]:


def deep_convert_dict(layer):
    import collections
    from load_utils import deep_convert_dict
    to_ret = layer
    if isinstance(layer, collections.OrderedDict):
        to_ret = dict(layer)

    try:
        for key, value in to_ret.items():
            to_ret[key] = deep_convert_dict(value)
    except AttributeError:
        pass

    return to_ret


# In[ ]:


def deep_convert_dict(layer):
    """
    Function to transform nested ordereddict to nested dict, 
    
    taken from http://stackoverflow.com/questions/25054003/how-to-convert-a-nested-ordereddict-to-dict
    
    from Patrick Collins
    
    """
    import collections
    to_ret = layer
    if isinstance(layer, collections.OrderedDict):
        to_ret = dict(layer)

    try:
        for key, value in to_ret.items():
            to_ret[key] = deep_convert_dict(value)
    except AttributeError:
        pass

    return to_ret


# In[ ]:


def load_hdf_raster1(datfile):
    """
    Name:

        load_hdf_raster1
    
    Purpose:

        To load the first raster of the hdf file. 
        To be used when load_hdf and load_hdf_sd did not work
        Brute force load first raster using GDAL
    
    Calling Sequence:

        dat = load_hdf_raster1(datfile) 
    
    Input: 
  
        datfile: full path and name of hdf file
    
    Output:

        dat: numpy array of values
    
    Keywords: 

       none
    
    Dependencies:

        numpy
        OSGEO, GDAL
    
    Required files:
   
        dat file
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2019-02-18, Santa Cruz, CA
        
    """
    import numpy as np
    from osgeo import gdal
    import os
    if not(os.path.isfile(datfile)):
        raise IOError('Data file {} not found!'.format(datfile))
    datsds = gdal.Open(datfile)
    rr = datsds.GetRasterBand(1) 
    g = rr.ReadAsArray()
    scale = rr.GetScale()
    offset = rr.GetOffset()
    nan = rr.GetNoDataValue()
    dat = g.astype(float)
    dat[dat==nan] = np.nan
    dat = dat*scale+offset
    del datsds
    return dat    


# Testing of the script:

# In[4]:


if __name__ == "__main__":
    import os
    import numpy as np
    from osgeo import gdal

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


# Testing the metadata dictionary

# In[ ]:


if __name__ == "__main__":

    gdal.Open(myd_dat_sub[53][0]).GetMetadata()

    mm = dict()
    mm['one'] = gdal.Open(myd_dat_sub[72][0]).GetMetadata()
    mm['two'] = gdal.Open(myd_dat_sub[74][0]).GetMetadata()
    mm['two']['_FillValue']

    from Sp_parameters import startprogress, progress, endprogress
    import gc; gc.collect()

    tuple(i[0] for i in modis_values).index('etau')

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

    print(list(modis.keys()))
    print(list(modis_dicts.keys()))

