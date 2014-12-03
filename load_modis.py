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
    """
    import numpy as np
    from osgeo import gdal
    from Sp_parameters import startprogress, progress, endprogress
    geosds = gdal.Open(geofile)
    datsds = gdal.Open(datfile)
    geosub = geosds.GetSubDatasets()
    datsub = datsds.GetSubDatasets()
    print 'Outputting the Geo subdatasets:'
    for i in range(len(geosub)):
        print str(i)+': '+geosub[i][1]
    print 'Outputting the Data subdatasets:'
    for i in range(len(datsub)):
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
    import gc; gc.collect()
    modis_dicts = dict()
    startprogress('Running through modis values')
    for i,j in modis_values:
        sds = gdal.Open(datsub[j][0])
        modis_dicts[i] = sds.GetMetadata()
        modis[i] = nd.array(sds.ReadAsArray())
        modis[i][modis[i] == float(modis_dicts[i]['_FillValue'])] = np.nan
        modis[i] = modis[i]*float(modis_dicts[i]['scale_factor'])+float(modis_dicts[i]['add_offset'])
        progress(float(tuple(i[0] for i in modis_values).index(i))/len(modis_values)*100.)
    endprogress()
    print modis.keys()
    return modis,modis_dicts

# <codecell>


# <markdowncell>

# Testing of the script:

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

