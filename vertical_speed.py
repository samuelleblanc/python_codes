# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Program to check out previous flights from the P3 to calculate various speeds..

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
import Pysolar.solar as sol
import datetime
fp='C:/Users/sleblan2/Research/flight_planning/p3_flights/'

# <codecell>

import load_modis as lm
reload(lm)

# <codecell>

from Sp_parameters import smooth

# <headingcell level=2>

# P3 during ARCTAS

# <codecell>

arctas,header = lm.load_ict(fp+'pds_p3b_20080419_r2.ict',return_header=True)
header

# <codecell>

vert_speed_ft = np.diff(arctas['GPS_ALT'])
arctas['GPS_ALT']

# <codecell>

vert_speed = vert_speed_ft*0.3084

# <codecell>

plt.plot(smooth(vert_speed,10),smooth(arctas['GPS_ALT'][1:]*0.3084,10),'b.')
plt.xlabel('Vertical speed [m/s]')
plt.ylabel('Altitude [m]')
plt.title('P-3 vertical speed for ARCTAS')

# <headingcell level=2>

# P3 during DISCOVER-AQ Denver

# <codecell>

discover = lm.load_ict(fp+'discoveraq-pds_p3b_20140807_r1.ict')

# <codecell>

d_vert_speed_ft = np.diff(discover['GPS_ALT'])

# <codecell>

d_vert_speed = d_vert_speed_ft*0.3084

# <codecell>

plt.plot(smooth(d_vert_speed,10),smooth(discover['GPS_ALT'][1:]*0.3084,10),'bo',label='DISCOVER-AQ:\n Boulder')
plt.plot(smooth(vert_speed,10),smooth(arctas['GPS_ALT'][1:]*0.3084,10),'r.',alpha=0.3,label='ARCTAS')
plt.xlabel('Vertical speed [m/s]')
plt.ylabel('Altitude [m]')
plt.title('P-3 vertical speed')
plt.legend(frameon=False,numpoints=1)
plt.savefig(fp+'P3_vert_speed.png',dpi=600,transparent=True)

# <codecell>

vs_up = smooth(vert_speed[(vert_speed>1)&(arctas['GPS_ALT'][1:]*0.3084>6000)],10)

# <codecell>

alt_up = smooth(arctas['GPS_ALT'][1:][(vert_speed>1)&(arctas['GPS_ALT'][1:]*0.3084>6000)],10)

# <codecell>

v,xx = linfit.linfit(alt_up,vs_up)

# <codecell>

p3_slope,p3_interp = v

# <codecell>

p3_slope

# <codecell>

alt1

# <codecell>

alt0 = 5000.0

# <codecell>

speed = 5.0

# <codecell>

climb_time = (alt1-alt0)/speed

# <codecell>

climb_time

# <headingcell level=2>

# For ER2 during SEAC4RS

# <codecell>

er2,header = lm.load_ict(fp+'seac4rs-nasdat_er2_20130821_r0.ict',return_header=True)

# <codecell>

er2_2 = lm.load_ict(fp+'seac4rs-nasdat_er2_20130922_r0.ict',return_header=False)

# <codecell>

header

# <codecell>

er2_vs = np.diff(er2['GPS_Altitude'])

# <codecell>

er22_vs = np.diff(er2_2['GPS_Altitude'])

# <codecell>

import plotting_utils as pu

# <codecell>

plt.plot(smooth(er22_vs,20),smooth(er2_2['GPS_Altitude'][1:],20),'.',color='grey')
plt.plot(smooth(er2_vs,20),smooth(er2['GPS_Altitude'][1:],20),'b.')
plt.plot(smooth(er2_vs[er2_vs<0],20),smooth(er2['GPS_Altitude'][1:][er2_vs<0],20),'r.')
pu.plot_lin(smooth(er2_vs[er2_vs>2],20),smooth(er2['GPS_Altitude'][1:][er2_vs>2],20))
pu.plot_lin(smooth(er2_vs[er2_vs<0],20),smooth(er2['GPS_Altitude'][1:][er2_vs<0],20),color='r')
plt.xlabel('Vertical speed [m/s]')
plt.ylabel('Altitude [m]')
plt.title('ER2 vertical speed from SEAC4RS')
plt.legend(frameon=False)
plt.savefig(fp+'ER2_vert_speed.png',dpi=600,transparent=True)

# <markdowncell>

# Get the inverse relationship for alt to vert speed

# <codecell>

import linfit
v = linfit.linfit(smooth(er2['GPS_Altitude'][1:][er2_vs>2],20),smooth(er2_vs[er2_vs>2],20))

# <codecell>

slope = v[0][0]
intercept = v[0][1]

# <codecell>

slope,intercept

# <headingcell level=2>

# For DC8 during SEAC4RS

# <codecell>

dc8,dc8header = lm.load_ict(fp+'nav_dc8_20080320_r1.ict',return_header=True)

# <codecell>

dc8header

# <codecell>

dc8_vs = np.diff(dc8['GPS_ALT'])

# <codecell>

plt.plot(smooth(dc8_vs,10),smooth(dc8['GPS_ALT'][1:],10),'b.')
plt.plot(15.0-0.001*np.linspace(0,12000),np.linspace(0,12000),label='vs = 15-0.001*alt')
plt.title('DC8 during SEAC4RS')
plt.xlabel('Vertical speed [m/s]')
plt.ylabel('Altitude [m/s]')
plt.legend(frameon=False)
plt.savefig(fp+'DC8_vert_speed.png',dpi=600,transparent=True)

# <headingcell level=2>

# C130 during ARISE

# <codecell>

c130,c130header = lm.load_ict(fp+'arise-C130-Hskping_c130_20140911_RA_Preliminary.ict',return_header=True)

# <codecell>

c130header

# <codecell>

plt.plot(smooth(c130['Vertical_Speed'],10),smooth(c130['GPS_Altitude'],10),'b.')
plt.plot(10-0.001*np.linspace(0,7500),np.linspace(0,7500),label='vs = 10-0.001*alt')
plt.title('C130 vertical speed during ARISE')
plt.xlabel('Vertical Speed [m/s]')
plt.ylabel('Altitude [m]')
plt.legend(frameon=False)
plt.savefig(fp+'C130_vert_speed.png',dpi=600,transparent=True)

