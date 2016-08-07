
# coding: utf-8

# # Info

# In[1]:

def __init__():
    """
    Name:  

        Sun_utils

    Purpose:  

        Collection of codes useful for calculating the sunphotometer analysis

    Input:

        none at command line

    Output:

        dict

    Keywords:

        none

    Dependencies:

        - numpy
        - Pysolar v0.6

    Needed Files:

      - ...

    Modification History:

        Written: Samuel LeBlanc, Mauna Loa Observatory, Hawaii, 2016-07-03
        Modified: 
    """
    pass


# # Prepare the airmass

# In[8]:

def airmass(sza,alt,h=21.0):
    """
    Purpose:
        Calculate airmasses, ported from matlab codes 2016-07-03
        Using Kasten and Young [1989], Komhyr [1989], and Kasten [1965]
    
    Input:
        sza: solar zenith angle
        alt: altitude in meter
        h: altitude of ozone layer in km (defaults to 21 km)
        
    Output: 
        m_ray : rayleigh airmass
        m_aero: aerosol airmass
        m_h20: water vapor airmass, same as aerosol
        m_o3: ozone airmass
        m_no2: nitrate airmass
        
    
    """
    import numpy as np

    # calculate airmasses
    if sza<96.07995:
        m_ray = 1./(np.cos(sza*np.pi/180)+0.50572*(96.07995-sza)**(-1.6364))
    else:
        m_ray = 1./(np.cos(sza*np.pi/180)+0.50572*(96.07995-96.0799)**(-1.6364))
    #Kasten and Young (1989)
    if sza<92.65:
        m_h2o = 1./(np.cos(sza*np.pi/180)+0.0548*(92.65-sza)**(-1.452))
    else:
        m_h2o = 1./(np.cos(sza*np.pi/180)+0.0548*(92.65-92.64)**(-1.452))
    #Kasten (1965)
    
    m_aero = m_h2o

    R = 6371.229
    r = alt/1000
    m_o3 = (R+h)/((R+h)**2-(R+r)**2.*(np.sin(sza*np.pi/180))**2)**0.5
    #Komhyr (1989)
    m_no2=m_o3
    return m_ray,m_aero,m_h2o,m_o3,m_no2  


# In[55]:

def tau_rayleigh(wavelength,pressure,date=None,latitude=45.0,atm_id=None):
    """
    Purpose:
        Calculate the rayleigh optical depth by using the time, pressure, latitude
        From Bucholtz A., 1995: Appl. Optics, Vol. 34, No. 15 2765-2773
        
    Input:
        wavelength: wavelength to be calculated, in microns
        pressure: atmospheric pressure at location in hPa
        date: datetime object of the measurement
        latitude: (optional, if atm_id is set, does not read this value, defaults to 45 N)
                  location of the measurement used for determining location
        atm_id: (default none) used to identify the default atmospheric profile
                1=Tropical
                2=Midlat Summer
                3=Midlat Winter
                4=Subarc Summer
                5=Subarc Winter
                6=1962 U.S. Stand Atm
    
    Output: 
        tau: rayleigh optical depth
        tau_err: error in rayleigh optical depth
        
    Dependencies:
        - Pysolar (v0.6)
        - numpy
        
    Needed Files:
        none
        
    Modification History:
        Written: Samuel LeBlanc, Mauna Loa Observatory, Hawaii, 2016-07-04, Happy Fourth of July!
                Ported from rayleighez.m and rayleigh.m malab files written by Yohei Shinozuka and Beat Schmid
    
    """
    from Pysolar import GetDeclination
    import numpy as np
    if not atm_id:
        try: 
            if not date:
                print 'No time set: using US Standard atmosphere for Rayleigh calculations'
                atm_id = 6
            else:
                latlims =[0,30,50,90]
        except:
            latlims =[0,30,50,90]
    elif not atm_id in [1,2,3,4,5,6]:
        print """***  atm_id not valid. Please select from below:
            1=Tropical
            2=Midlat Summer
            3=Midlat Winter
            4=Subarc Summer
            5=Subarc Winter
            6=1962 U.S. Stand Atm
            
        """
    try:
        doy = date.timetuple().tm_yday
        declination = np.array(GetDeclination(doy))
    except AttributeError:
        declination = np.array([GetDeclination(d.timetuple().tm_yday) for d in date])
    
    if abs(latitude) < latlims[1]: 
        atm_id = 1 #tropical
    elif abs(latitude) < latlims[2]: 
        if latitude*declination[0]>0:
            atm_id = 2 # midlatitude summer
        else:
            atm_id = 3 # midlatitude winter
    elif abs(latitude) < latlims[3]:
        if latitude*declination[0]>0:
            atm_id = 4 # subarctic summer
        else:
            atm_id = 5 # subarctic winter
    else: 
        atm_id = 6 # US standard
            
    # following is ported from rayleigh.m code written by Beat S. 1995?
    #Bucholz(95) Table 5 constants
    Coeff_model_atmo = np.array([[0.00652965,0.00868094],                        [0.00651949,0.00866735],                        [0.00653602,0.00868941],                        [0.00648153,0.00861695],                        [0.00649997,0.00864145],                        [0.00650362,0.00864627]])

    A = Coeff_model_atmo[atm_id+1,0]
    B = 3.55212
    C = 1.35579
    D = 0.11563
    wavelength = np.array(wavelength)
    tau = np.array(A/wavelength**(B+C*wavelength+D/wavelength))

    i = wavelength>0.500
    A = Coeff_model_atmo[atm_id+1,1]
    B = 3.99668
    C = 1.10298e-3
    D = 2.71393e-2

    tau[i] = A/wavelength[i]**(B+C*wavelength[i]+D/wavelength[i])
    tau = (pressure/1013.25*tau)
    tau_err = 0.015 # relative error that reflects the estimated accuracy of pressure measurements 
    return tau, tau_err


# In[61]:

def get_season(date,north=True):
    """
    Program that takes input of date to return the season. date is in format of datetime object. 
    If in southern hemisphere set north to False
        
    Does not account for utc
    """
    if north:
        seasons = {'Summer':(datetime(date.year,6,21), datetime(date.year,9,22)),
                   'Autumn':(datetime(date.year,9,23), datetime(date.year,12,20)),
                   'Spring':(datetime(date.year,3,21), datetime(date.year,6,20))}
    else:
        seasons = {'Summer':(datetime(date.year,12,21), datetime(date.year,3,20)),
                   'Spring':(datetime(date.year,9,23), datetime(date.year,12,20)),
                   'Autumn':(datetime(date.year,3,21), datetime(date.year,6,20))}
    for season,(season_start, season_end) in seasons.items():
        if date>=season_start and date<= season_end:
            return season
    return 'Winter'


# In[ ]:

def calc_sza_airmass(datetime_utc,lat,lon,alt,c={}):
    """
    Purpose:  
        Take in an array of datetime in UTC timoezone aware format
        then use it to calculate the solar zenith angle and airmass along the time.
        put in the lattitude and longitude 

    Input:
        datetime_utc: list of datetime objects in utc format
        lat: latitude of measurement (Positive=Northern hemisphere, right now only for ground based)
        lon: longitude of measurements, POsitive=Eastern hemisphere, (only one val, not list yet)
        alt: altitude of measurements, (meters)
        
    Output:
        dict with latitude, longitude, datetime_utc, altitude and 
        new sza, and airmass calcuations inside the dict

    Keywords:
        none

    Dependencies:
        - Sun_utils
        - map_utils

    Needed Files:
      - None

    Modification History:

        Written: Samuel LeBlanc, Mauna Loa Observatory, Hawaii, 2016-07-07
        Modified: 
    """
    import map_utils as mu
    import Sun_utils as su
    c['lat'] = lat
    c['lon'] = lon
    c['alt'] = alt
    c['sza'] = []
    c['azi'] = []
    c['sunearthf'] = []
    c['m_aero'] = []
    c['m_ray'] = []
    c['m_o3'] = []
    if not type(lat) is float:
        lat_list = True
    else:
        lat_list = False
    for i,d in enumerate(datetime_utc):
        if lat_list:
            s,a,sf = mu.get_sza_azi(lat[i],lon[i],d,alt=alt[i],return_sunearthfactor=True)
        else:
            s,a,sf = mu.get_sza_azi(lat,lon,d,alt=alt,return_sunearthfactor=True)
        c['sza'].append(s[0])
        c['azi'].append(a[0]) 
        c['sunearthf'].append(sf[0])
        m_r,m_a,_,m_o,_ = su.airmass(s[0],alt)
        c['m_aero'].append(m_a)
        c['m_ray'].append(m_r)
        c['m_o3'].append(m_o)
    return c


# In[ ]:



