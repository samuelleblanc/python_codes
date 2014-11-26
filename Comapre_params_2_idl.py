# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Program that is used to test the results of params in the python vs. idl (original)

# <markdowncell>

# Load the required modules

# <codecell>

%config InlineBackend.rc = {}
import matplotlib
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
from matplotlib import colors
%matplotlib inline

# <codecell>

import numpy as np, h5py
import plotly.plotly as py
import scipy.io as sio
import math
import os
import warnings
warnings.simplefilter('ignore', np.RankWarning)
import Sp_parameters as Sp
py.sign_in("samuelleblanc", "4y3khh7ld4")
#import mpld3
#mpld3.enable_notbeook()

# <codecell>

reload(Sp)

# <codecell>

from Sp_parameters import nanmasked, closestindex, norm2max, find_closest, smooth, deriv
from scipy.interpolate import UnivariateSpline, splev, splrep

# <codecell>

# set the basic directory path
fp='C:\\Users\\sleblan2\\Research\\TCAP\\'

# <codecell>

# load the idl save file containing the modeled radiances and build the Sp modules
s=sio.idl.readsav(fp+'model/sp_v1_20130219_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp, wvl, z, re, ta)'
lut = Sp.Sp(s)
lut.params()

# <codecell>

# load the matlab file containing the measured TCAP radiances
m = sio.loadmat(fp+'4STAR/20130219starzen_rad.mat')
sm = sio.idl.AttrDict(m)
print sm.keys()
print 'Measured radiance Shape: ', sm.rad.shape
meas = Sp.Sp(m)
meas.params()

# <markdowncell>

# Now go through comparing the two methods

# <codecell>

print lut.par.shape
print lut.par[0,5,6,:]
p_0_5_6 = [ 0.169785, -0.672731 ,     7.21242 ,    0.979948 ,    0.194668,     0.109177,     0.338863,      3.11233,   -0.0178171,   -0.0139593,     -1.29179  ,   0.330518
      ,1.17178    ,  1.69885, -0.000469405 ,    0.178697]
print p_0_5_6
print lut.par[1,10,10,:]
p_1_10_10 = [0.132723570, -1.930865765,  0.662532628,  1.008047819,  0.163861081,  0.025557771,
             0.315011382, -0.094412826, -0.016095554, -0.006612694,-1.353369117,  0.306448370, 
             1.179942489,  1.728420854,  0.008580528,  0.051876802,]
print p_1_10_10

# <markdowncell>

# Plot the differences in parameters

# <codecell>

fig4, ax4 = plt.subplots(1,2, figsize=(15,5))
ax4 = ax4.ravel()
ax4[0].plot(np.arange(16),lut.par[0,5,6,:],label='python')
ax4[0].plot(np.arange(16),p_0_5_6,label='idl')
ax4[0].set_title('difference in parameters - liquid example')
ax4[0].set_xlabel('Parameters')
ax4[0].legend()

ax4[1].plot(np.arange(16),lut.par[1,10,10,:],label='python')
ax4[1].plot(np.arange(16),p_1_10_10,label='idl')
ax4[1].set_title('difference in parameters - ice example')
ax4[1].set_xlabel('Parameters')
ax4[1].legend()

# <markdowncell>

# Create the desired arrays which are used in creating the parameters

# <codecell>

s_0_5_6 = lut.sp[0,:,0,5,6]
s,n = nanmasked(s_0_5_6)
snorm = norm2max(s_0_5_6)
[i1000,i1077,i1493,i1600,i1200,i1300,i530,i610,
 i1565,i1634,i1193,i1198,i1236,i1248,i1270,i1644,
 i1050,i1040,i1065,i600,i870,i515] = find_closest(lut.wvl,np.array([1000,1077,1493,1600,1200,1300,530,
                                                     610,1565,1634,1193,1198,1236,1248,
                                                     1270,1644,1050,1040,1065,600,870,515]))
norm2 = s_0_5_6/s_0_5_6[i1000]
dsp = smooth(np.gradient(norm2,lut.wvl/1000.),2)

# <codecell>

norm2_uni = UnivariateSpline(lut.wvl/1000.0,norm2,k=5)
norm2_uni.set_smoothing_factor(1)
dnorm2 = norm2_uni.derivative()

# <codecell>

norm2_bspline = splrep(lut.wvl/1000.0,norm2,k=5)
norm2_b = splev(lut.wvl/1000.0,norm2_bspline,der=0)
dbnorm2 = splev(lut.wvl/1000.0,norm2_bspline,der=1)

# <codecell>

dsp2 = smooth(deriv(norm2,lut.wvl/1000.),2)

# <codecell>

plt.figure()
plt.plot(lut.wvl,norm2)
plt.plot(lut.wvl,norm2_uni(lut.wvl/1000.0))
plt.plot(lut.wvl,norm2_b)
plt.figure()
plt.plot(lut.wvl,dnorm2(lut.wvl/1000.0))
plt.plot(lut.wvl,dbnorm2)

# <codecell>

f1,a1 = plt.subplots(5,sharex=True,figsize=(10,15))
for a in a1:
    a.grid(True)
a1[0].plot(lut.wvl,s_0_5_6)
a1[0].set_xlim([350,1700])
a1[3].set_xlabel('Wavelength [nm]')
a1[1].plot(lut.wvl,snorm)
a1[2].plot(lut.wvl,norm2)
a1[3].plot(lut.wvl,dsp)
a1[4].plot(lut.wvl,dsp2)
a1[0].set_title('original')
a1[1].set_title('normed to max')
a1[2].set_title('normed to 1000 nm')
a1[3].set_title('derivative of normed at 1000 nm using gradient')
a1[4].set_title('another derivative normed at 1000 nm using 3 point lagrangian')

# <markdowncell>

# Compare to IDL plots:
# spectra:
# <img src="./compare_plots/s_0_5_6.png">
# Normed to max:
# <img src="./compare_plots/norm.png">
# Normed to 1000 nm:
# <img src="./compare_plots/norm2.png">
# Derivative of above:
# <img src="./compare_plots/dsp.png">

# <markdowncell>

# Verify distinct pieces of the derivative plot, and the slope fit to them

# <codecell>

f2,a2 = plt.subplots(1)
a2.plot(lut.wvl,dsp2,'b+')
a2.plot(lut.wvl,dsp2,'b')
a2.set_xlim([1000,1077])
a2.set_ylim([-7,2])
a2.grid(True)

# build the linear fit
from linfit import linfit
f,z = linfit(lut.wvl[i1000:i1077],dsp2[i1000:i1077])
print f
print z
print f[1]
a2.plot(lut.wvl[i1000:i1077],f[1]+f[0]*lut.wvl[i1000:i1077],'r')

#linfit output is [m,b], not [b,m] like previously thought

# <codecell>

f2,a2 = plt.subplots(1)
a2.plot(lut.wvl,dsp2,'bx')
a2.plot(lut.wvl,dsp2,'b')
a2.set_xlim([1150,1250])
a2.set_ylim([-5,4])
a2.plot([1193],[-0.67],'r*')
a2.grid(True)

# <codecell>


