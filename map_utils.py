# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def __init__():
    """
       Collection of codes to run some typical map utilities
       
           - spherical_dist: codes to calculate the distance from a certain lat lon points
           - map_ind: find the indices of the closest point 
           - radius_m2deg: return the radius of a circle defined by meters to lat lon degrees, 
           
        details are in the info of each module
    """
    pass

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
    from map_utils import spherical_dist
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

# <codecell>

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
    import geopy
    from geopy.distance import VincentyDistance
    origin = geopy.Point(center_lon,center_lat)
    destination = VincentyDistance(kilometers=radius/1000.0).destination(origin,0.0)
    radius_degrees = abs(center_lat-destination.latitude)
    return radius_degrees

