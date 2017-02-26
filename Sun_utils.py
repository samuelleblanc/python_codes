
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
        Modified: Samuel LeBlanc, From SFO to Seoul, 2017-02-25
                  Added tools to calculate the angstrom exponent from a polynomial fit
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

def tau_rayleigh(wavelength,pressure,date=None,latitude=45.0,declination=None,atm_id=None):
    """
    Purpose:
        Calculate the rayleigh optical depth by using the time, pressure, latitude
        From Bucholtz A., 1995: Appl. Optics, Vol. 34, No. 15 2765-2773
        
    Input:
        wavelength: wavelength to be calculated, in microns
        pressure: atmospheric pressure at location in hPa
        date: (optional if declination is not set) datetime object of the measurement
        latitude: (optional, if atm_id is set, does not read this value, defaults to 45 N)
                  location of the measurement used for determining location
        declination: (optional) if set ignores the date keyword, if atm_id is set this is ignored
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
        - Pyephem
        - numpy
        
    Needed Files:
        none
        
    Modification History:
        Written: Samuel LeBlanc, Mauna Loa Observatory, Hawaii, 2016-07-04, Happy Fourth of July!
                Ported from rayleighez.m and rayleigh.m malab files written by Yohei Shinozuka and Beat Schmid
        Modified: Samuel LeBlanc, Santa Cruz, CA, 2016-08-08
                Removed the dependency of Pysolar, and replaced with Pyephem. added the declination keyword
    
    """
    import ephem
    import numpy as np
    if not atm_id:
        try: 
            if not date and not declination:
                print 'No time set: using US Standard atmosphere for Rayleigh calculations'
                atm_id = 6
            else:
                latlims =[0,30,50,90]
        except:
            latlims =[0,30,50,90]
        
        try:
            sun = ephem.Sun()
            sun.compute(date)
            declination = sun.dec*180.0/np.pi
        except AttributeError, ValueError:
            declination = []
            sun = ephem.Sun()
            for d in date:
                sun.compute(d)
                declination.append(sun.dec)
            declination = np.array(declination)

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
    elif not atm_id in [1,2,3,4,5,6]:
        print """***  atm_id not valid. Please select from below:
            1=Tropical
            2=Midlat Summer
            3=Midlat Winter
            4=Subarc Summer
            5=Subarc Winter
            6=1962 U.S. Stand Atm
            
        """
            
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
    c['declination'] = []
    c['m_aero'] = []
    c['m_ray'] = []
    c['m_o3'] = []
    if not type(lat) is float:
        lat_list = True
    else:
        lat_list = False
    for i,d in enumerate(datetime_utc):
        if lat_list:
            s,a,sf,dec = mu.get_sza_azi(lat[i],lon[i],d,alt=alt[i],return_sunf_and_dec=True)
        else:
            s,a,sf,dec = mu.get_sza_azi(lat,lon,d,alt=alt,return_sunf_and_dec=True)
        c['sza'].append(s[0])
        c['azi'].append(a[0]) 
        c['sunearthf'].append(sf[0])
        c['declination'].append(dec[0])
        m_r,m_a,_,m_o,_ = su.airmass(s[0],alt)
        c['m_aero'].append(m_a)
        c['m_ray'].append(m_r)
        c['m_o3'].append(m_o)
    return c


# # Tools for calculating the angstrom exponent

# In[ ]:

def aod_polyfit(wvl,aod,polynum=False):
    """
    Purpose:  
        Take in an array of aod spectra and calculate the polynomials associated with each spectrum
        takes in the wvl

    Input:
         aod: numpy array of time,wvl
         wvl: numpy array of wavelengths in nm
        
    Output:
        array of polynomial coefficients for each time point in the aod

    Keywords:
        polynum: (optional, defaults to N-2 where N is the number of points in the spectrum) 
                 number of polynomials coefficients to return (order of the polynomial), if set to False, uses default N-2
        
    Dependencies:
        - numpy
        - scipy

    Needed Files:
      - None

    Modification History:

        Written: Samuel LeBlanc, Flight between San Francisco and Seoul, 2017-02-25
        Modified: 
    """
    import numpy as np
    from scipy import polyfit
    # sanitize inputs
    shape = aod.shape
    wvl = np.array(wvl)
    if not len(wvl.shape)==1:
        raise ValueError('wvl is not a array with one dimension')
    if not len(wvl) in shape:
        raise ValueError('No size of aod is the same as wvl')
        
    if not polynum:
        polynum = len(wvl)-2
        
    if len(shape)>1:
        ni,n = [(i,j) for (i,j) in enumerate(shape) if not j==len(wvl)][0]
        cc = np.zeros((polynum,n))
        for i in range(n):
            if ni==0: 
                cc[:,i] = polyfit(wvl,aod[i,:],polynum)
            else:
                cc[:,i] = polyfit(wvl,aod[:,i],polynum)
    else:
        cc = polyfit(wvl,aod,polynum)
        
    return cc


# In[86]:

def angstrom_from_poly(c,wvl):
    """
    Purpose:  
        calculate the angstrom at a defined wavelength, can take in arrays

    Input:
        c: polynomial coefficients from output of aod_polyfit (can be an array of coefficients,time)
        wvl: wavelengths (in nm) to calculate the angstrom exponent
        
    Output:
        array of angstrom values, one for each time in c

    Keywords:
        polynum: (optional, defaults to 4) number of polynomials coefficients to return (order of the polynomial) 
        
    Dependencies:
        - numpy
        - scipy

    Needed Files:
      - None

    Modification History:

        Written: Samuel LeBlanc, Flight between San Francisco and Seoul, 2017-02-25
    """
    print 'Not tested yet!!'
    from scipy import polyval
    import numpy as np
    if len(c.shape)>1:
        try:
            ang = np.zeros((len(c[0,:]),len(wvl)))
            for iv,v in enumerate(wvl):
                vu = np.array([v-0.2,v-0.1,v,v+0.1,v+0.2])
                dv = np.gradient(vu)
                for j in xrange(len(c[0,:])):
                    a = np.gradient(np.log(polyval(c[:,j],vu))*(-1.0),dv )
                    ang[j,iv] = a[2]
        except TypeError:
            ang = np.zeros((len(c[0,:])))
            v = wvl
            vu = np.array([v-0.2,v-0.1,v,v+0.1,v+0.2])
            dv = np.gradient(vu)
            for j in xrange(len(c[0,:])):
                a = np.gradient(np.log(polyval(c[:,j],vu))*(-1.0),dv)
                ang[j,iv] = a[2]
    else:
        try: 
            ang = np.zeros((len(wvl)))
            for iv,v in enumerate(wvl):
                vu = np.array([v-0.2,v-0.1,v,v+0.1,v+0.2])
                dv = np.gradient(vu)
                a = np.gradient(np.log(polyval(c,vu))*(-1.0),dv )
                ang[iv] = a[2]
        except TypeError:
            v = wvl
            vu = np.array([v-0.2,v-0.1,v,v+0.1,v+0.2])
            dv = np.gradient(vu)
            a = np.gradient(np.log(polyval(c,vu))*(-1.0),dv)
            ang = a[2]
    return ang


# In[ ]:



