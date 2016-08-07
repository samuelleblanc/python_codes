
# coding: utf-8

# Collection of programs to used for reading and analyzing the C-air utils 

# In[50]:

def read_Cair(filein):
    """
    Purpose:  
        read the Cair csv files in the form of URC or UCV
        This reads in the raw data, then parses the times for datetime, creates a utc column

    Input:
        filein: full path to file name to open

    Output:
        dict of values that have been read from the file

    Keywords:
        none

    Dependencies:
        - numpy
        - load_utils
        - dateutil
        - Pysolar v0.6
        - pytz
        - time
        - datetime

    Needed Files:

      - filein 

    Modification History:

        Written: Samuel LeBlanc, Mauna Loa Observatory, Hawaii, 2016-07-05
        Modified: 
    """
    from load_utils import toutc, recarray_to_dict
    import numpy as np
    from dateutil import parser
    from datetime import datetime
    import time
    import pytz
    
    print 'read first'
    c = np.genfromtxt(filein,skip_header=0,names=True,delimiter=',')
    names = c.dtype.names
    co = recarray_to_dict(c)
    
    print 'try for second read'
    try:
        print 'trying datetime parse'
        def dparse(x):
            return datetime.strptime(x,)
            
        def dparse_utc(x):
            return datetime(*time.strptime(str2,'%m/%d/%Y %H:%M:%S.%f %p')[0:6],tzinfo=pytz.timezone('UTC'))
        c1 = np.genfromtxt(filein,skip_header=0,dtype='object',names=True,delimiter=',',                           converters={'DateTime':dparse,'DateTimeUTC':dparse_utc})
    except:
        print 'trying parser'
        def parse_utc(x):
            return parser.parse(x+' UTC')

        c1 = np.genfromtxt(filein,skip_header=0,dtype='object',names=True,delimiter=',',                           converters={'DateTime':parser.parse,'DateTimeUTC':parse_utc})
    print 'toutc'
    #co['utc'] = toutc(c1['DateTimeUTC'])
    co['DateTime'] = c1['DateTime']
    co['DateTimeUTC'] = c1['DateTimeUTC']
    
    # Compile individual Lt arrays into a full table
    co['wvl'] = []
    co['Lt'] = []
    for n in names:
        if n.startswith('Lt'):
            co['wvl'].append(float(n.split('_')[0].strip('Lt')))
            co['Lt'].append(co[n])
    co['Lt'] = np.array(co['Lt'])
    
    return co   


# In[ ]:

def get_darks(c):
    """
    Purpose:  
        Analyze the C-air values and return the ones with darks.

    Input:
        read_Cair dict output, looks for the utc array and DateTimeUTC
        added dict values of lat, lon, and alt unless specified in the call for this function
        
    Output:
        new sza, and airmass calcuations inside the dict

    Keywords:
        lat: latitude of measurement in degrees positive for north
        lon: longitude of measurements in degrees positive for East
        alt: altitude of measurement in meters

    Dependencies:
        - Sun_utils
        - map_utils

    Needed Files:
      - None

    Modification History:

        Written: Samuel LeBlanc, Mauna Loa Observatory, Hawaii, 2016-07-05
        Modified: 
    """


# In[1]:

def calc_rayleigh(c,press=None):
    """
    Purpose:  
        Wrapper to analyze the C-air read data and calculate the rayleigh scattering

    Input:
        read_Cair dict output
        
    Output:
        new sza, and airmass calcuations inside the dict

    Keywords:
        lat: latitude of measurement in degrees positive for north
        lon: longitude of measurements in degrees positive for East
        alt: altitude of measurement in meters

    Dependencies:
        - Sun_utils
        - map_utils

    Needed Files:
      - None

    Modification History:

        Written: Samuel LeBlanc, Mauna Loa Observatory, Hawaii, 2016-07-07
        Modified: 
    """
    import Sun_utils as su
    if not press:
        raise InputError('No pressure indicated please put one in')
    c['tau_ray'] = su.Sun_utils(c['wvl'],press,)
    
    print 'not yet'


# In[ ]:



