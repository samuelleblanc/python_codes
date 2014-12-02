# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def load_modis(geofile,datfile):
"""
Name:

    load_modis
    
Purpose:

    to compile functions required to load Modis files
    from within another script
    
Calling Sequence:

    modis = load_modis(geofile,datfile) 
    
Input: 
  
    geofile name
    datfile name (hdf files)
    
Output:

    modis class with tau, ref, etau, eref, phase, qa,
    
Keywords: 

    none
    
Dependencies:

    gdal
    numpy
    
Required files:
   
    geo and dat files
    
Example:

    ...
"""
from load_modis import modis_cls
from osgeo import gdal

    

# <codecell>


