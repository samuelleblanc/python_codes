def get_aeronet(daystr=None,lat_range=[],lon_range=[],lev='LEV10',avg=True,daystr2=None,version='2'):
    """ 
    Purpose:
       Program to go and get the aeronet data on the day defined by daystr
       Returns an numpy named array, with station names, lat, lon, average daily aods per site and wavelength
    Inputs:
       daystr: (optional, defaults to today) day for aeronet data in format yyyy-mm-dd.
               if day is in future, it is changed to today 
       lat_range: (optional) Latitude range for which to load AERONET sites, if none, then returns all AERONET
       lon_range: (optional) longitude range for which to load AEORNET sites, if none, then returns all AERONET
       lev: (defaults to LEV10), determine which level of aeronet to load, choices are LEV10, LEV15, LEV20 for 1.0, 1.5, and 2.0 respectively
            If version is set to 3, LEV is changed to AOD, new options include the Spectral Deconvolution Algorithm (SDA10, SDA15, SDA20), 
            and the Total AOD data (TOT10, TOT15, TOT20)
       avg: (defaults to True), if set to True returns daily averaged values, but if false, returns every measurement point.
       daystr2: (optional) if set, is used to be the last days of the spanned AERONET values
       version: (default 2) version of the direct beam data, can be set to 3
    Outputs:
       numpy structured array with one entry per station
    Dependencies:
       urllib
       StringIO
       BeautifulSoup
       Numpy
       datetime
       warnings
    Example:
       ...
    History:
       Written: Samuel LeBlanc, 2015-10-02, Santa Clara, CA
       Modified: Samuel LeBlanc, 2016-07-12, WFF, CA
       Modified: Samuel LeBlanc, 2018-03-30, Mountain View, CA
                 - added keywords for choosing which aeronet level to read (lev) and if to return the averages or not.
    """
    import numpy as np
    from BeautifulSoup import BeautifulSoup
    from StringIO import StringIO
    from urllib import urlopen
    from datetime import datetime
    from load_utils import recarray_to_dict
    import pandas as pd
    
    # validate input
    if avg:
        avg = '20'
    else:
        avg = '10'
        
    if version=='3':
        lev = lev.replace('LEV','AOD')
        url_name = 'print_web_data_v3'
    else:
        url_name = 'print_web_data_v2_globe'
    
    def safe_list_get (l, idx, default):
        try:
            return l[idx]
        except IndexError:
            return default

    dd_now = datetime.utcnow()
    dd = dd_now.strftime('%Y-%m-%d')
    if not daystr:
        daystr = dd
    else:
        if daystr > dd:
	    daystr = dd
	    import warnings
	    warnings.warn("Date set to future, using today's date")
    if not daystr2: 
        daystr2 = daystr
  
    url = 'http://aeronet.gsfc.nasa.gov/cgi-bin/{urlnm}?year={yyyy}&month={mm:02.0f}&day={dd:02.0f}&year2={yyyy2}'+\
          '&month2={mm2:02.0f}&day2={dd2:02.0f}&{lev}=1&AVG={avg}'
    if lat_range:
        lat_range.sort()
        lon_range.sort()
        url = url+'&lat1={lat1:f}&lat2={lat2:f}&lon1={lon1:f}&lon2={lon2:f}'
    url = url.format(urlnm=url_name,yyyy=daystr[0:4],mm=int(daystr[5:7]),dd=int(daystr[8:10]),lev=lev,avg=avg,lat1=safe_list_get(lat_range,0,None),
                    lat2=safe_list_get(lat_range,1,None),lon1=safe_list_get(lon_range,0,None),lon2=safe_list_get(lon_range,1,None),
                    yyyy2=daystr2[0:4],mm2=int(daystr2[5:7]),dd2=int(daystr2[8:10]))
    print 'Getting file from internet: at aeronet.gsfc.nasa.gov'
    print url
    try:
        htm = urlopen(url)
        html = htm.read()
        soup = BeautifulSoup(html)
    except:
        print 'failed to communicate with AERONET internet site - returning nothing'
        return False
    lines = []
    for ibr,br in enumerate(soup.findAll('br')):
        nt = br.nextSibling
        if (version=='3') & (ibr<3):
            print nt
            continue
        if len(lines)==0:
            if 'Number_of_Wavelengths' in nt:
                nt = nt.strip()+',exact_wvl2,exact_wvl3,exact_wvl4,exact_wvl5'
        lines.append(nt.strip()+'\n')
    s = StringIO(''.join(lines))
    s.seek(0)
    
    try:
        dat_pd = pd.read_csv(s)
        dat = dat_pd.to_records()
    except:
        try:
            dat = np.genfromtxt(s,delimiter=',',names=True,dtype=None)
        except IndexError:
            print 'Failed to read the returned html file'
            #return s
            return False
        except ValueError:
            print 'Failed to format the html file, returning the strings'
            return s
        fields_to_ignore = ['AERONET_Site_Name','Principal_Investigator','PI_Email','Dateddmmyy']
        for label in dat.dtype.names:
            if not label in fields_to_ignore:
                if dat[label].dtype.type is np.str_:
                     dat[label] = np.genfromtxt(dat[label])

    return recarray_to_dict(dat)
    
def plot_aero(m,aero,no_colorbar=True,a_max = 1.5):
    """
    Simple function that takes a basemap plotting object ( m) and plots the points of the aeronet sites with their value in color
    For easy aeronet visualisation
    """
    from matplotlib import cm
    from matplotlib.lines import Line2D
    import numpy as np
    x,y = m(aero['Longitude'].astype('double'),aero['Latitude'].astype('double'))
    
    if no_colorbar:
        colors = np.round(aero['AOT_500'].astype('double')/a_max*7.0)/7.0
        c_ar = np.linspace(0,a_max,7)
        leg_ar = ['{:1.2f} - {:1.2f}'.format(c,c_ar[i+1]) for i,c in enumerate(c_ar[0:-1])]
    else:
        colors = aero['AOT_500'].astype('double')

    cls = cm.gist_ncar(c_ar/a_max)
    
    bb = m.scatter(x,y,c=colors,cmap=cm.gist_ncar,marker='s',
                   vmin=0.0,vmax=a_max,edgecolors='None',s=100)
                   
    if no_colorbar:
        fakepoints = []
        for i,cl in enumerate(cls):
            fakepoints.append(Line2D([0],[0],color=cl,linestyle='None',marker='s',markeredgecolor='None'))
        cbar = m.ax.legend(fakepoints,leg_ar,numpoints=1,frameon=True,loc='lower right',bbox_to_anchor=(0.45,1.04),title='AOD 500 nm',ncol=2)
        try:
            m.ax.add_artist(cbar)
        except:
            pass
    else:
        try:
            cbar = m.colorbar(m.ax,bb)
            cbar.set_label('AOD 500 nm')
        except:
            pass
    u = [bb,cbar]
    return u