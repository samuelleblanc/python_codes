#!/usr/bin/env python
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


# In[1]:


def logaod_polyfit(wvl,aod,polynum=False):
    """
    Purpose:  
        Take in an array of aod spectra and calculate the polynomials associated with each log of the spectrum
        takes in the wvl

    Input:
         aod: numpy array of time,wvl
         wvl: numpy array of wavelengths in nm
        
    Output:
        array of polynomial coefficients for each time point in the log(aod)

    Keywords:
        polynum: (optional, defaults to N-2 where N is the number of points in the spectrum) 
                 number of polynomials coefficients to return (order of the polynomial), if set to False, uses default N-2
        
    Dependencies:
        - numpy
        - scipy

    Needed Files:
      - None

    Modification History:

        Written: Samuel LeBlanc, NASA Ames Research Center, 2017-11-28
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
                cc[:,i] = polyfit(np.log(wvl),np.log(aod[i,:]),polynum)
            else:
                cc[:,i] = polyfit(np.log(wvl),np.log(aod[:,i]),polynum)
    else:
        cc = polyfit(np.log(wvl),np.log(aod),polynum)
        
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


# In[2]:


def angstrom_from_logpoly(c,wvl,polynum=4):
    """
    Purpose:  
        calculate the angstrom at a defined wavelength, can take in arrays

    Input:
        c: polynomial coefficients from output of logaod_polyfit (can be an array of coefficients,time) 
           ***log(aod) vs. log(wvl) polynomial coefficients only
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

        Written: Samuel LeBlanc, NASA Ames Research Center, 2017-11-28
    """
    from scipy import polyval
    import numpy as np
    co = np.array([c])
    ix = np.where(np.array(co.shape)==polynum+1)[0][0]
    cn = np.rollaxis(co,ix,0)
    nn = np.where(np.array(cn.shape)==1)[0][0]
    cc = np.rollaxis(cn,nn,len(co.shape))
    n = cc.shape[1]
    ang = []
    for i,iv in enumerate(np.array([wvl])): # along wavelength axis
        a = [np.polyval(np.polyder(cc[:,j].flatten()),np.log(iv))*(-1.0) for j in xrange(n)]
        ang.append(a)    
    ang = np.array(ang[0])
    return ang


# In[1]:


def get_angstrom(aod,wvl_in,wvl_out,polynum=4):
    """
    Purpose:  
        calculate the angstrom exponent at a defined wavelength, using a polynomial fit to aerosol optical depth
        and returns the angstrom exponent at specified wavelengths. Combines above function of 
        logaod_polyfit and angstrom_from_logpoly

    Input:
        aod: array of aerosol optical depth, with one of the array length related to 
             the wavelength they are used to be described (wvl_in)
        wvl_in: wavelengths (in nm) relating to the spectra of aod
        wvl_out: wavelengths (in nm) to calculate the angstrom exponent
        
    Output:
        numpy array of angstrom values, one for each time in aod

    Keywords:
        polynum: (optional, defaults to 4) number of polynomials coefficients to use (order of the polynomial) 
        
    Dependencies:
        - numpy
        - Sun_utils (this file)

    Needed Files:
      - None

    Modification History:

        Written: Samuel LeBlanc, Santa Cruz, CA, 2018-03-14
    """
    import numpy as np
    import Sun_utils as su
    from tqdm import tqdm
    
    # reorganise and sanitze the inputs
    ao = np.array([aod])
    ix = np.where(np.array(ao.shape)==len(wvl_in))[0][0]
    an = np.rollaxis(ao,ix,0)
    nn = np.where(np.array(an.shape)==1)[0][0]
    ac = np.rollaxis(an,nn,len(ao.shape))
    n = ac.shape[1]
    
    pbar = tqdm(total=n+1)
    poly = []
    for i in xrange(n):
        pl = su.logaod_polyfit(wvl_in,ac[i,:],polynum=polynum)
        poly.append(pl)
        pbar.update(1)
    poly = np.array(poly)
    
    angs = su.angstrom_from_logpoly(poly,wvl_out,polynum=polynum)
    pbar.update(1)
    return angs


# In[1]:


def calc_angs(time,w,aod):
    'Program to calculate the angstrom exponent by fitting linearly on the aod'
    import numpy as np
    from linfit import linfit
    from tqdm import tqdm_notebook as tqdm
    ang = np.zeros_like(time)
    if not (np.array(w).shape[0] == aod.shape[0]): raise ValueError('Size of AOD and wavelenght arrays not consistent')
    pbar = tqdm(total=len(time)+2)
    pbar.update(1)
    for i,t in enumerate(time):
        try:
            c,cm = linfit(np.log10(w),-np.log10(aod[:,i]))
            p = np.array([c[1],c[0]])
            ang[i] = c[0] 
        except:
            ang[i] = np.nan
        pbar.update(1)
    pbar.update(1)
    if not any(np.isfinite(ang)):
        raise ValueError('No Angstrom Exponent value returned valid')
    return ang


# # The Spectral deconvolution algorithm for fine/coarse mode fraction

# In[3]:


def sda(aod,wvl_in,ref_wvl=0.5,polynum=2,compute_errors=False,physical_forcing=False,bias_corr=True):
    """
    Purpose:  
        Compute the total, fine and coarse mode aerosol optical depth (AOD), the total and fine mode Angstrom exponents (alpha),
        and the total spectral derivative of the angstrom exponent (alpha') at a reference wavelength.
        
        Adapted from tauf_tauc.m matlab function from Norm O'Neill (see below for references and contact info)

    Input:
        aod: array of aerosol optical depth, with one of the array length related to 
             the wavelength they are used to be described (wvl_in)
        wvl_in: wavelengths (in microns) relating to the spectra of aod
        ref_wvl: (defaults to 0.5 micron), the reference wavelength [micron] to use for sda
        polynum: (defaults to 2) number of coefficients for best fit polynomial decomposition of log(aod) vs log(wvl) 
        compute_errors: (defaults to False) if True computes the errors and returns another array
        physical_forcing: (defaults to False) if True determines whether a physical forcing correction is applied
        bias_corr: (defaults to True) if True, sets the bias correction of the calculated values from fine/coarse mode tau
        
    Output:
        dictionary with these values:
            tauf: numpy array of aod fine mode fraction at reference wavelength
            tauc: numpy array of aod coarse mode fraction at reference wavelength
            tau: numpy array of aod best fit at reference wavelength
            etauf: (optional, only returned if compute_errors is True) stochastic error of tauf 
            etauc: (optional, only returned if compute_errors is True) stochastic error of tauc

    Keywords:
        polynum: (optional, defaults to 2) number of polynomials coefficients to use (order of the polynomial) 
        
    Dependencies:
        - numpy
        - Sun_utils (this file)

    Needed Files:
      - None

    Modification History:
        Written: (adpated from original matlab codes) Samuel LeBlanc, Santa Cruz, CA, 2018-03-15
    
    ORIGINAL Comments:
    % HISTORY 
    % 1. VERSION 1 was basically a simplification and merging of CIMEL_process (program employed for ref. 1 and ref. 2) and CIMEL_errors into
    %    the tauf_tauc function.
    % 2. VERSION 2; on Sept. 1, 2004 a corrected option was added corresponding to physical_forcing = 'y'. This version corrects non-physical values
    %    of eta (outside the range [0, 1]). This exclusion is done smoothly by averaging the correction between the extremes of non-physical limits
    %    and allowable values corresponging to the extremes of the estimated stochastic error. This new option also incorporated a new estimate
    %    of the rms AOD error at the reference wavelength (change from 0.01 to 0.006); a lower value was used because the 2nd order polynomial 
    %    AOD fitting error will be less than the raw AOD error as estimated from AERONET Langley plots. Other changes from VERSION 1 included the addition
    %    of a global set of variables and new output parameters in the output variable list. Setting physical_forcing = 'n' and Dtau = 0.01 yields
    %    the outputs obtained in references 1 and 2.
    % 3. VERSION 3; on March 19, 2005 a new version was created. This version incorporated physical forcing "at the upper
    %    alpha_f end" : when AOD input errors were very large then alpha_f could be above the theoretical upper limit of 
    %    alpha_f (= 4 for Rayleigh particles). The new version sets a limit of the theoretical alpha_f value for alpha_f
    %    (in precisely the same way that the lower limit of alpha was set in VERSION 2). Also, wavelength dependent 
    %    coefficients of the alpha'_f versus alpha_f relationship were incorporated (and these relationships, besides
    %    being used in the deconvolution code proper were incorporated in the calculation of the alpha_f dependent 
    %    error of alpha'_f. These relationships are derived in visible.xls. 
    % 4. VERSION 4; on April 26, 2006 a new version was delivered to AERONET.
    %    An overly complicated expression for the stochastic error in alpha' was simplified to Dalphap = k1*Dtau_rel 
    %    where k1 is simply = 10 (the more complicated expression in ref. 1 did not produce signficantly better 
    %    estimates of the empirical errors derived from stochastic simulations). Physical forcing was rendered
    %    less extreme by incorporating a quadratic weighting scheme rather than the MOE (Mean of extrema) averaging
    %    employed in Version 3. This modification (for cases near eta = 1 and eta = 0) produced moderate changes
    %    in tauf, tauc and eta and eliminated a suspicious correlation between derived values of tauf and tauc. These 
    %    changes are discussed in the newest version of the tauf_tauc technical memo (reference 3 below).
    %
    %    VERSION 4.1; on Mar. 19, 2008 and April 1, 2008 some minor code changes were incorporated to respectively (a) correct an inconsistancy 
    %    in the alpha_f_max_theoretical loop (basically it was just moved above the "if alpha_f_min < alpha | alpha_c_max > alpha" loop and a 
    %    few consequentially redundant executable statements were removed) and (b) to correct the error model code (details in ref. 3). Other code 
    %    cleanup was performed Updates were also made to the code documentation.
    %
    %   ** WARNING **
    %   The algorithm has only been physically validated at a reference wavlength of 0.5 um
    %   The input AOD spectrum must have at least 4 optical depths. The wavelengths employed in the 0.5 um reference 
    %   wavelength validation exercise were [1.02 0.870 0.670 0.5 0.44 0.38] um (in other words 0.34 um was 
    %   systematically excluded). In April of 2007 the 1020 nm channel was eliminated from AERONET processing (see ref. 3)
    %
    % REFERENCES
    % 1.    O'Neill, N. T., Eck, T. F., Smirnov, A., Holben B. N., S. Thulasiraman, (2003), Spectral discrimination 
    %       of coarse and fine mode optical depth, JGR, Vol. 108, No. D17, 4559.
    % 2.    O'Neill, N. T., Dubovik, O., Eck, T. F., (2001), A modified Angstrom coefficient for the characterization 
    %       of sub-micron aerosols, App. Opt., Vol. 40, No. 15, pp. 2368-2374.
    % 3.    O'Neill, N. T., Eck, T. F., Smirnov, A., Holben B. N., S. Thulasiraman, (2008), Spectral Deconvolution 
    %       algorithm Technical memo.
    %
    % CONTACT
    % Norm O'Neill
    % Senior Scientist and Professor,
    % CARTEL,Universite de Sherbrooke,Sherbrooke, Quebec, Canada, J1K 2R1
    % tel.; (819)-821-8000, ext 2965, fax; (819)-821-7965
    % Email; norm.oneill@usherbrooke.ca
    % web;  http://pages.usherbrooke.ca/noneill/ 
    
    """
    import numpy as np
    from scipy import polyval
    import Sun_utils as su
    from tqdm import tqdm_notebook as tqdm
    
    # fix coarse mode Angstrom parameters
    alpha_c = -0.15
    alphap_c = 0.0 # generic case as per ref. 2
    
    # Set up fine mode constants for alpha_f estimation
    # a, b, c, b_star, and c_star defined in reference 2 (with the following correction; b* = b + 2 a alpha_c).
    # see visible.xls (sheet "AERONET") for the source of the equations immediately below (available from the author)
    # upper and lower refer to the upper and lower error bounds of the alphap_f vs alpha_f expression
    
    # reorganise and sanitze the inputs
    ao = np.array([aod])
    ix = np.where(np.array(ao.shape)==len(wvl_in))[0][0]
    an = np.rollaxis(ao,ix,0)
    nn = np.where(np.array(an.shape)==1)[0][0]
    ac = np.rollaxis(an,nn,len(ao.shape))
    n = ac.shape[1]
    
    pbar = tqdm(total=n+2)
    pbar.update(1)
    
    #Compute the spectral polynomial (ln-ln fit)
    poly = []
    tau,tauf,tauc,eta = [],[],[],[]
    alphap,alpha = [],[]
    for i in xrange(n):
        #if i>55446:
        #    import pdb; pdb.set_trace()
        t,e,tf,tc,al,alp,pl = calc_orig_tauc_tauf(ac[:,i,0],wvl_in,ref_wvl=ref_wvl,
                                                  polynum=polynum,alpha_c=alpha_c,alphap_c=alphap_c)

        poly.append(pl)
        tau.append(t)
        alpha.append(al)
        alphap.append(alp)
        eta.append(e)
        tauf.append(tf)
        tauc.append(tc)
        
        if bias_corr:
            if polynum == 2: # the bias correction was only developed for the default second order fit 
                tf, tc = bias_corr_tauf_tauc(tau[i],eta[i],alpha[i],alphap[i],ref_wvl=ref_wvl,alpha_c=alpha_c,alphap_c=alphap_c)
                tauf[i] = tf
                tauc[i] = tc
        pbar.update(1)

    tauf = np.array(tauf)
    tauc = np.array(tauc)
    tau = np.array(tau)
    poly = np.array(poly)
    alpha = np.array(alpha)
    alphap = np.array(alphap)
    eta = np.array(eta)
    
    out = {'tauf':tauf,'tauc':tauc,'tau':tau,'poly':poly,'alpha':alpha,'alphap':alphap,'eta':eta}
    pbar.update(1)
    
    return out
    
    


# In[2]:


def calc_orig_tauc_tauf(aod,wvl_in,ref_wvl=0.5,polynum=2,alpha_c=0.15,alphap_c=0.0):
    """
    Purpose:  
        Compute the total, fine and coarse mode aerosol optical depth (AOD), 
        subset function from the original fine mode and coarse mode for single values
        
        fine mode Angstrom exponents (alpha) and the total spectral derivative of the angstrom exponent (alphap)
        at a reference wavelength
        
        Adapted from tauf_tauc.m matlab function from Norm O'Neill (see below for references and contact info)

    Input:
        aod: 1-d array of aerosol optical dept at each the wavelength they are used to be described (wvl_in)
        wvl_in: wavelengths (in nm) relating to the spectra of aod
        ref_wvl: (defaultsto 0.5 micron), the reference wavelength [micron] to use for sda
        polynum: (defaults to 2) number of coefficients for best fit polynomial decomposition of log(aod) vs log(wvl) 
        alpha_c: (defaults to 0.15) angstrom constant for fine/coarse mode seperation
        alphap_c: (defaults to 0.0) derivative of angstrom constant for fine/coarse mode seperation
        
    Output:
        tau: aod best fit at reference wavelength
        eta: ratio of fine mode fraction (of total aod)
        tauf: aod fine mode fraction at reference wavelength
        tauc: aod coarse mode fraction at reference wavelength
        alpha: angstrom exponent
        alphap: angstrome derivative exponent
        poly: polynomial coefficients

    Dependencies:
        - numpy
        - Sun_utils (this file)
        - scipy

    Needed Files:
      - None

    Modification History:
        Written: (adpated from original matlab codes) Samuel LeBlanc, Santa Cruz, CA, 2018-03-22
    
    """
    import numpy as np
    import Sun_utils as su
    from scipy import polyval
    
    # constants defined
    a_upper = -0.22
    b_upper = (10**-0.2388)*(ref_wvl**1.0275)
    c_upper = (10**0.2633)*(ref_wvl**(-0.4683))
    a_lower = -0.3
    b_lower = 0.8
    c_lower = 0.63
    a = (a_lower + a_upper)/2
    b = (b_lower + b_upper)/2
    c = (c_lower + c_upper)/2
    b_star = b + 2*alpha_c*a
    c_star = c - alphap_c + b*alpha_c + a*alpha_c**2.0
    #print 'c_star: {}'.format(c_star)
    
    pl = su.logaod_polyfit(wvl_in,aod,polynum=polynum)
    
    # compute the aod at the reference wavelength
    ar = polyval(pl,np.log(ref_wvl))
    tau = np.exp(ar)

    # compute alpha
    alpha = 0.0
    for i_d in xrange(1,polynum+1):
        expon = polynum-i_d
        alpha = alpha - (expon+1)*pl[i_d-1]*np.log(ref_wvl)**expon

    # compute alphapà
    alphap = 0.0
    for i_d in xrange(2,polynum+1):
        expon = polynum-i_d
        alphap = alphap - (expon+1)*(expon+2)*pl[i_d-2]*np.log(ref_wvl)**expon

    # Compute the derived fine and coarse mode parameters. "t" is the invariant generic parameter defined in reference 2.
    alpha_offset = alpha - alpha_c
    alphap_offset = alphap - alphap_c
    t = alpha_offset - alphap_offset/alpha_offset
    alpha_f = ((t + b_star) + np.sqrt((t + b_star)**2.0 + 4*(1 - a)*c_star))/(2*(1 - a)) + alpha_c
    alpha_f_offset = alpha_f - alpha_c
    eta = alpha_offset/alpha_f_offset
    #print 'alpha_offset: {}, alpha_f_offset: {}'.format(alpha_offset,alpha_f_offset)
    tau_f = eta*tau
    tau_c = tau - tau_f
    
    return tau, eta, tau_f, tau_c, alpha, alphap, pl


# In[ ]:


def bias_corr_tauf_tauc(tau,eta,alpha,alphap,ref_wvl=0.5,alpha_c=0.15,alphap_c=0.0):
    """
    Purpose:  
        Bias correction of the tauf and tauc
        
        Adapted from tauf_tauc.m matlab function from Norm O'Neill (see below for references and contact info)

    Input:
        aod: 1-d array of aerosol optical dept at each the wavelength they are used to be described (wvl_in)
        wvl_in: wavelengths (in nm) relating to the spectra of aod
        ref_wvl: (defaults to 0.5 micron), the reference wavelength [micron] to use for sda
        polynum: (defaults to 2) number of coefficients for best fit polynomial decomposition of log(aod) vs log(wvl) 
        alpha_c: (defaults to 0.15) angstrom constant for fine/coarse mode seperation
        alphap_c: (defaults to 0.0) derivative of angstrom constant for fine/coarse mode seperation
        
    Output:
        tau: aod best fit at reference wavelength
        eta: ratio of fine mode fraction (of total aod)
        tauf: aod fine mode fraction at reference wavelength
        tauc: aod coarse mode fraction at reference wavelength
        alpha: angstrom exponent
        alphap: angstrome derivative exponent
        poly: polynomial coefficients

    Dependencies:
        - numpy
        - Sun_utils (this file)
        - scipy

    Needed Files:
      - None

    Modification History:
        Written: (adpated from original matlab codes) Samuel LeBlanc, Santa Cruz, CA, 2018-03-22
    
    """
    import numpy as np
        
    # constants defined
    a_upper = -0.22
    b_upper = (10**-0.2388)*(ref_wvl**1.0275)
    c_upper = (10**0.2633)*(ref_wvl**(-0.4683))
    a_lower = -0.3
    b_lower = 0.8
    c_lower = 0.63
    a = (a_lower + a_upper)/2
    b = (b_lower + b_upper)/2
    c = (c_lower + c_upper)/2
    b_star = b + 2*alpha_c*a
    c_star = c - alphap_c + b*alpha_c + a*alpha_c**2.0    
    
    # set the offset
    alpha_offset = alpha - alpha_c
        
    # Compute the bias correction
    alphap_bias_correction = 0.65*np.exp(-(eta - 0.78)**2.0/(2*0.18**2.0))

    alphap_bias_corrected = alphap + alphap_bias_correction
    t_bias_corrected =          alpha_offset - (alphap_bias_corrected - alphap_c)/alpha_offset
    alpha_f_bias_corrected =          ((t_bias_corrected + b_star) + np.sqrt((t_bias_corrected + b_star)**2.0 + 4.*(1 - a)*c_star) )/(2.*(1 - a)) + alpha_c
    eta_bias_corrected = alpha_offset/(alpha_f_bias_corrected - alpha_c)
    tau_f_bias_corrected = eta_bias_corrected*tau
    tau_c_bias_corrected = tau - tau_f_bias_corrected
    t = t_bias_corrected
    alphap = alphap_bias_corrected
    alphap_offset = alphap - alphap_c # this corrected offset is used in the error section below
    eta = eta_bias_corrected
    alpha_f = alpha_f_bias_corrected
    alphap_f = a*alpha_f**2.0 + b*alpha_f + c #quadratic relation of ref. 2
    tau_f = tau_f_bias_corrected
    tau_c = tau_c_bias_corrected
    alpha_f_offset = alpha_f - alpha_c
    
    return tau_f, tau_c

