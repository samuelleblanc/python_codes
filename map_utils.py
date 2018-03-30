
# coding: utf-8

# In[ ]:


def __init__():
    """
       Collection of codes to run some typical map utilities
       
           - spherical_dist: codes to calculate the distance from a certain lat lon points
           - map_ind: find the indices of the closest point 
           - radius_m2deg: return the radius of a circle defined by meters to lat lon degrees, 
           
        details are in the info of each module
    """
    pass


# In[1]:


def spherical_dist(pos1, pos2, r=6378.1,use_mi=False):
    """
    Calculate the distance, in km if using the default radius, 
    from one point to another (can use arrays)
    pos1 is array of [lat,lon]
    pos2 is array of [lat,lon]
    if use_mi = True, radius set to 3958.75 miles, default to False, with radius of 6378.1 km

    Modified: Samuel LeBlanc, NASA Ames, Santa Cruz, CA, 2015-09-02
    """
    if use_mi:
        r = 3958.75
        print 'using miles'
    import numpy as np
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
    pos1 = pos1 * np.pi / 180
    pos2 = pos2 * np.pi / 180
    cos_lat1 = np.cos(pos1[..., 0])
    cos_lat2 = np.cos(pos2[..., 0])
    cos_lat_d = np.cos(pos1[..., 0] - pos2[..., 0])
    cos_lon_d = np.cos(pos1[..., 1] - pos2[..., 1])
    return r * np.arccos(cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))


# In[1]:


def bearing(pos1,pos2):
    "Calculate the initial bearing, in degrees, to go from one point to another, along a great circle"
    import numpy as np
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
    pos1 = pos1 * np.pi / 180
    pos2 = pos2 * np.pi / 180
    cos_lat1 = np.cos(pos1[..., 0])
    cos_lat2 = np.cos(pos2[..., 0])
    sin_lat1 = np.sin(pos1[...,0])
    sin_lat2 = np.sin(pos2[...,0])
    sin_lon_d = np.sin(pos1[...,1]-pos2[...,1])
    cos_lat_d = np.cos(pos1[..., 0] - pos2[..., 0])
    cos_lon_d = np.cos(pos1[..., 1] - pos2[..., 1])
    return 360.0-((np.arctan2(sin_lon_d*cos_lat2,cos_lat1*sin_lat2-sin_lat1*cos_lat2*cos_lon_d)*180.0/np.pi+360.0) % 360.0)


# In[2]:


def map_ind(mod_lon,mod_lat,meas_lon,meas_lat,meas_good=None):
    """ Run to get indices in the measurement space of all the closest mod points. Assuming earth geometry."""
    from map_utils import spherical_dist
    from Sp_parameters import startprogress, progress, endprogress
    import numpy as np
    try:
        if not meas_good:
            meas_good = np.where(meas_lon)
    except ValueError:
        if not meas_good.any():
            meas_good = np.where(meas_lon)
        
    imodis = np.logical_and(np.logical_and(mod_lon>min(meas_lon[meas_good])-0.02 , mod_lon<max(meas_lon[meas_good])+0.02),
                            np.logical_and(mod_lat>min(meas_lat[meas_good])-0.02 , mod_lat<max(meas_lat[meas_good])+0.02))
    wimodis = np.where(imodis)
    if not wimodis[0].any():
        print '** No points found within range +/- 0.02 in lat and lon, Extending range to +/- 0.2 **'
        imodis = np.logical_and(np.logical_and(mod_lon>min(meas_lon[meas_good])-0.2 , mod_lon<max(meas_lon[meas_good])+0.2),
                                np.logical_and(mod_lat>min(meas_lat[meas_good])-0.2 , mod_lat<max(meas_lat[meas_good])+0.2))
        wimodis = np.where(imodis)
        if not wimodis[0].any():
            print '** No points found in extended range, returning null **'
            return []
    N1 = mod_lon[imodis].size
    modis_grid = np.hstack([mod_lon[imodis].reshape((N1,1)),mod_lat[imodis].reshape((N1,1))])
    try:
        N2 = len(meas_good)
        if N2==1 or N2==2:
            meas_good = meas_good[0]
            N2 = len(meas_good)
        meas_grid = np.hstack([np.array(meas_lon[meas_good]).reshape((N2,1)),np.array(meas_lat[meas_good]).reshape((N2,1))])
    except:
        import pdb; pdb.set_trace()
    meas_in = meas_grid.astype(int)
    meas_ind = np.array([meas_good.ravel()*0,meas_good.ravel()*0])
    startprogress('Running through flight track')
    for i in xrange(meas_good.size):
        d = spherical_dist(meas_grid[i],modis_grid)
        try:
            meas_ind[0,i] = wimodis[0][np.argmin(d)]
        except:
            import pdb; pdb.set_trace()
        meas_ind[1,i] = wimodis[1][np.argmin(d)]
        progress(float(i)/len(meas_good)*100)
    endprogress()
    return meas_ind


# In[ ]:


def radius_m2deg(center_lon,center_lat,radius):
    """ 
    Return the radius in lat lon degrees of a circle centered at the points defined by
      center_lon
      center_lat
    with a radius defined in meters by:
      radius
      
    Dependencies:
        
        - geopy library
    """
    from geopy import Point
    from geopy.distance import VincentyDistance
    origin = Point(center_lat,center_lon)
    destination = VincentyDistance(kilometers=radius/1000.0).destination(origin,0.0)
    radius_degrees = abs(center_lat-destination.latitude)
    return radius_degrees


# In[1]:


def stats_within_radius(lat1,lon1,lat2,lon2,x2,radius,subset=True):
    """
    Run through all points defined by lat1 and lon1 (can be arrays)
    to find the points within defined by lat2 and lon2 that are within a distance in meters defined by radius
    lat2, lon2, x2 can be multidimensional, will be flattened first
    if subset (optional) is set to True, and there are more than 100 points, only every 10th in lat1, lon1 will be used.
    Returns a dicttionary of statistics:
        'index' : array of indices of flattened lat2 and lon2 that are within radius meters of each point of lat1 and lon1
        'std' : array of standard deviation of x2 that are near lat1 and lon1 by radius
        'range' : range of values of x2 near lat1, lon1
        'mean' : mean of values of x2 near lat1, lon1
        'median': median values of x2 near lat1, lon1
    """
    from scipy.spatial import cKDTree
    from map_utils import radius_m2deg
    import numpy as np
    print 'Setting up the lat, lon, localization'
    max_distance = radius_m2deg(lon1[0],lat1[0],radius) #transform to degrees
    if (len(lat1) > 100) & subset:
        points_ref = np.column_stack((lat1[::10],lon1[::10]))
    else:
        points_ref = np.column_stack((lat1,lon1))
    if len(lat2.shape) > 1:
        points = np.column_stack((lat2.reshape(lat2.size),lon2.reshape(lon2.size)))
        xx = x2.reshape(x2.size)
    else:
        points = np.column_stack((lat2,lon2)) 
        xx = x2
    tree = cKDTree(points)
    tree_ref = cKDTree(points_ref)
    out = dict()
    print '... Getting the index points'
    out['index'] = tree_ref.query_ball_tree(tree,max_distance)
    out['std'] = []
    out['range'] = []
    out['mean'] = []
    out['median'] = []
    print '... Running through index points'
    for i in out['index']:
        if not i:
            out['std'].append(np.NaN)
            out['range'].append(np.NaN)
            out['mean'].append(np.NaN)
            out['median'].append(np.NaN)
        else:
            out['std'].append(np.nanstd(xx[i]))
            out['range'].append(np.nanmax(xx[i])-np.nanmin(xx[i]))
            out['mean'].append(np.nanmean(xx[i]))
            out['median'].append(np.median(xx[i]))
    out['std'] = np.array(out['std'])
    out['range'] = np.array(out['range'])
    out['mean'] = np.array(out['mean'])
    out['median'] = np.array(out['median'])
    print out.keys()
    return out


# In[3]:


def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    """
    plot a single circle on a map
    uses the shoot function below
    from: http://www.geophysique.be/2011/02/20/matplotlib-basemap-tutorial-09-drawing-circles/
    by: Thomas Lecocq
    """
    from map_utils import shoot
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])

    #m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X,Y = m(X,Y)
    line = m.ax.plot(X,Y,**kwargs)
    return line


# In[4]:


def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    import numpy as np
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)


# In[ ]:


def great(m, startlon, startlat, azimuth,*args, **kwargs):
    """
    function to draw great circle, takes into account crossing the border
    by: Thomas Lecocq
    """
    glon1 = startlon
    glat1 = startlat
    glon2 = glon1
    glat2 = glat1
 
    step = 50
 
    glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    if azimuth-180 >= 0:
        while glon2 <= startlon:
            line = m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)
 
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    elif azimuth-180 < 0:
        while glon2 >= startlon:
            line = m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)
 
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    return line


# In[1]:


def get_sza_azi(lat,lon,datetimet,alt=None,return_sunearthfactor=False,return_sunf_and_dec=False):
    """
    Program wrapper for pyephem.Sun to get the solar zenith angle and the solar azimuth angle
    can use inputs of list or numpy arrays
    require input of lat,lon,datetimet
    optional input of altitutde (in meters)
    optional output of sun earth distance factor if return_sunearthfactor is set to True
    optional output of sun earth distance factor and sun declination if return_sunf_and_dec is set to True
    """
    import ephem
    from numpy import pi,isscalar
    sun = ephem.Sun()
    obs = ephem.Observer()
    if isscalar(lat):
        if isscalar(datetimet):
            lat = [lat]
            lon = [lon]
            datetime = [datetimet]
        else:
            lati = [lat for i in xrange(len(datetimet))]
            loni = [lon for i in xrange(len(datetimet))]
            lat,lon = lati,loni
    n = len(lat)
    sza = []
    azi = []
    sunf = []
    dec = []
    for i in range(n):
        obs.lat,obs.lon,obs.date = lat[i]/180.0*pi,lon[i]/180.0*pi,datetimet[i]
        if alt:
            obs.elevation = alt
        sun.compute(obs)
        sza.append(90.0-sun.alt*180/pi)
        azi.append(sun.az*180/pi)
        sunf.append(1.0/(sun.earth_distance**2))
        dec.append(sun.dec*180.0/pi)
    if return_sunf_and_dec:
        return sza,azi,sunf,dec
    elif return_sunearthfactor:
        return sza,azi,sunf
    else:
        return sza,azi


# In[1]:


def consecutive(data, stepsize=1):
    'simple program to get consecutive values'
    import numpy as np
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)


# In[3]:


def mplot_spec(m,lon,lat,*args,**kwargs):
    'Program to plot lines on a map, wihtout the extra cross sides of lines because of the dateline problem'
    import numpy as np
    from map_utils import consecutive
    latrange = [m.llcrnrlat,m.urcrnrlat]
    lonrange = [m.llcrnrlon,m.urcrnrlon]
    lon = np.array(lon)
    lat = np.array(lat)
    ii, = np.where((lat<=latrange[1])&(lat>=latrange[0])&(lon>=lonrange[0])&(lon<=lonrange[1]))
    ic = consecutive(ii)
    lines = []
    for c in ic:
        if c[0] != 0 : c = np.insert(c,0,c[0]-1)
        if c[-1] != len(lon)-1: c = np.append(c,c[-1]+1)
        x,y = m(lon[c],lat[c])
        lines.append(m.plot(x,y,*args,**kwargs))
    return lines  


# In[3]:


def PolyArea(x,y):
    'Program to calculate the area within a polygon defined by vertices x,y (both arrays)'
    import numpy as np
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))


# In[2]:


def WithinArea(xpoint,ypoint,xpoly,ypoly):
    'return True if point (at xpoint,ypoint coordinates) within area defined by xpoly,ypoly vertices'
    from map_utils import PolyArea
    import numpy as np
    area = PolyArea(xpoly,ypoly)
    area_sup = 0.0
    for i in xrange(len(xpoly)-1):
        area_sup += PolyArea([xpoint,xpoly[i],xpoly[i+1]],[ypoint,ypoly[i],ypoly[i+1]])
    return area_sup<=area

