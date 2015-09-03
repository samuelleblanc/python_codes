# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     calc_sza
# 
# Purpose:  
# 
#     Python script to calculate the solar zenith angle and airmass factor using the pysolar algorithm
#     
# Calling Sequence:
# 
#     python calc_sza.py
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - pysolar
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   
#  Modification History:
#  
#      Written: by Samuel LeBlanc, NASA Ames, 2015-08-20

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
import Pysolar.solar as sol
import datetime
fp='C:/Users/sleblan2/Research/NAAMES/plot/'

# <codecell>

lat = 47.6
lon = np.array([-40.0,-45.0,-50.0,-52.8])
utc = np.arange('2015-11-12T00:00Z','2015-11-13T00:00Z',15,dtype='datetime64[m]').astype(datetime.datetime)

# <codecell>

uu = datetime.datetime(2015,11,12,0,0,10)

# <codecell>

uu.timetuple()

# <codecell>

t = 20.6

# <codecell>

int(floor(t))

# <codecell>

u = utc[0]

# <codecell>

u.timetuple()

# <codecell>

sza = np.zeros((len(lon),len(utc)))

# <codecell>

azi = np.zeros((len(lon),len(utc)))

# <codecell>

utc.size

# <codecell>

utcsub = utc[0::10]

# <codecell>

utcsub

# <codecell>

latsub = [lat,lat,lat,lat,lat,
          lat,lat,lat,lat,lat]

# <codecell>

lonsub = [lon[0],lon[0],lon[0],lon[0],lon[0],
          lon[0],lon[0],lon[0],lon[0],lon[0]]

# <codecell>

szasub = 90.0-sol.GetAltitude(latsub,lonsub,utcsub)

# <codecell>

for i in xrange(len(lon)):
    for j in xrange(len(utc)):
        azi[i,j] = sol.GetAzimuth(lat,lon[i],utc[j])
        sza[i,j] = 90.0-sol.GetAltitude(lat,lon[i],utc[j])

# <codecell>

m = 1.0/np.cos(sza*np.pi/180.0)

# <codecell>

def datetime2utch(dt):
    if len(dt)>1:
        h = np.zeros(len(dt))
        for i,t in enumerate(dt):
            xtup = t.timetuple()
            h[i] = xtup[3]+xtup[4]/60.0+xtup[5]/3600.00+xtup[6]/10**6 
    else:
        xtup = dt.timetuple() 
        h = xtup[3]+xtup[4]/60.0+xtup[5]/3600.00+xtup[6]/10**6 
    return h

# <codecell>

utch = datetime2utch(utc)

# <codecell>

plt.plot(utch,m[0,:])
plt.ylim([1,25])
plt.ylabel('Airmass')
plt.xlabel('UTC [h]')

# <codecell>

plt.plot(utch,sza[0,:])

# <codecell>

plt.plot(utch,azi[0,:])

# <codecell>

sol.GetAltitude(36.6,-122.6,datetime.datetime(2015,8,20,13,0))

# <codecell>

datetime.datetime(2015,8,20,13,0).

# <codecell>


