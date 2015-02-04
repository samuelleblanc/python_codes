# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     Animate_4STAR_v2
# 
# Purpose:  
# 
#     Python script that is used to view 4STAR full spectral raw counts
#     Animates the raw over time
# 
# Calling Sequence:
# 
#     python Animate_4STAR_v2.py
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     animated figures 
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - load_modis: for mat2py_time and toutc functions
#     - Sp_parameters: for find_closest function
#     - mpl_toolkits
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - General star.mat: or other 4STAR star file containing raw spectra

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
%matplotlib inline
import numpy as np, h5py
#import plotly.plotly as py
import scipy.io as sio
import math
import os
from load_modis import mat2py_time, toutc
from Sp_parameters import find_closest

# <codecell>

# set the basic directory path
fp = 'C:\\Users\\sleblan2\\Research\\ARISE\\c130\\'
datestr = '20141002'
filesep = '\\'
print filesep

# <codecell>

file_name = fp+'20141002_ARISE_Flight_16\\20141002flt16starsun.mat' #'20141002_ARISE_Flight_16\\20141002flt16star.mat'

# <headingcell level=2>

# Load the 4STAR data and select time points

# <codecell>

star = sio.loadmat(file_name)
star.keys()

# <codecell>

if 't' in star.keys():
    star['utc'] = toutc(mat2py_time(star['t']))
else:
    star['utc'] = toutc(mat2py_time(star['vis_sun']['t'][0][0]))
#star['vis_sun']['t'][0][0]

# <codecell>

star['good'] = np.where((star['utc']>25.0) & 
                        (star['utc']<27.6))[0]# & 
                        #(star['Str'].flatten()!=0) & 
                        #(star['sat_time'].flatten()==0) & 
                        #(np.isfinite(star['tau_aero'][:,319])))[0]
len(star['good'])

# <codecell>

if 't' in star.keys():
    rate = star['rateaero'][star['good'],:]
    print rate.shape
else:
    raw_vis = star['vis_sun']['raw'][0][0][star['good'],:]
    print raw_vis.shape

# <headingcell level=2>

# Plot time traces at select wavelengths

# <codecell>

print star['tau_aero'].shape
print star['w'].shape
iw = find_closest(star['w'][0],np.array([380,400,430,500,515,635,875,1000,1140,1200,1600])/1000.0)
print star['w'][0].shape
print iw

# <codecell>

plt.figure()
color.cycle_cmap(len(iw)+1,cmap=plt.cm.gist_ncar)
for i in iw:
    plt.plot(star['utc'][star['good']],star['tau_aero'][star['good'],i], label=str('%.3f $\mu$m' % star['w'][0,i]))
plt.xlabel('UTC [h]')
plt.ylabel('AOD')
plt.legend(bbox_to_anchor=[1,0.1],loc=3)
plt.savefig(fp+datestr+filesep+'AOD_time.png',dpi=600)

# <headingcell level=2>

# Start animation of full spectral aod

# <codecell>

from IPython.display import HTML

def display_animation(anim):
    plt.close(anim._fig)
    return HTML(anim_to_html(anim))

# <codecell>

from tempfile import NamedTemporaryFile

def anim_to_html(anim):
    #mywriter = animation.FFMpegWriter(fps=10)
    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.gif') as f:
            anim.save(f.name, fps=5)
            video = open(f.name, "rb").read()
        anim._encoded_video = video.encode("base64")
    
    return VIDEO_TAG.format(anim._encoded_video)

# <codecell>

from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(0.35, 1.75))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    plt.xlabel('Wavelength [$\mu$m]')
    plt.ylabel('AOD')
    plt.grid()
    return line,

# animation function.  This is called sequentially
def animate(i):
    line.set_data(star['w'][0], star['tau_aero'][star['good'][i]])
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(star['good']), interval=10, blit=True)

# call our new function to display the animation
#display_animation(anim)

# <codecell>

print star['tau_aero'][star['good'][0]]
print star['w'][0]

# <codecell>

def animate(i):
    plt.cla()
    plt.plot(star['w'][0],star['tau_aero'][star['good'][i]])
    plt.ylim(0,1)
    plt.xlim(0.30,1.75)
    plt.grid()
    plt.xlabel('Wavelength [$\mu$m]')
    plt.ylabel('AOD')
    plt.title(str('UTC: %2.4f' % star['utc'][star['good'][i]]))
    plt.tight_layout()

fig = plt.figure(figsize=(5,4))

anim = animation.FuncAnimation(fig, animate, frames=len(star['good']))
anim.save(fp+datestr+filesep+'AOD_sp_anim.gif', writer=animation.ImageMagickWriter(), fps=4);
#anim.save(fp+datestr+filesep+'AOD_sp_anim.mp4', writer=animation.FFMpegWriter(), fps=10,extra_args=['-vcodec', 'libx264']);

# <codecell>

import base64
with open(fp+datestr+filesep+'AOD_sp_anim.m4v', "rb") as image_file:
    encoded_string = base64.b64encode(image_file.read())
VIDEO_TAG = """<video controls>
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""
HTML(VIDEO_TAG.format(encoded_string))

# <codecell>

u = video.encode('base64')

# <codecell>

HTML(data='<video controls><source src="data:video/x-m4v;base64,'+u+'" type="video/mp4"></video>')

# <codecell>

def video(fname, mimetype):
    """Load the video in the file `fname`, with given mimetype, and display as HTML5 video.
    """
    from IPython.display import HTML
    video_encoded = open(fname, "rb").read().encode("base64")
    video_tag = '<video controls alt="test" src="data:video/{0};base64,{1}">'.format(mimetype, video_encoded.decode('ascii'))
    return HTML(data=video_tag)

# <codecell>

u = open(fp+datestr+filesep+'AOD_sp_anim.m4',"rb").read()

# <codecell>

print len(u)

# <codecell>

video(fp+datestr+filesep+'AOD_sp_anim.m4v','x-m4v')

# <codecell>

from JSAnimation.IPython_display import display_animation
display_animation(anim)

# <headingcell level=2>

# Prepare of saving specific variables to netcdf

# <codecell>

wvl = star['w'][0]
aod = star['tau_aero'][star['good']]
utc = star['utc'][star['good']]
lat = star['Lat'][star['good']]
lon = star['Lon'][star['good']]
alt = star['Alt'][star['good']]
sza = star['sza'][star['good']]

# <codecell>

print wvl.shape
print wvl

# <codecell>

print aod.shape
print aod

# <codecell>

print utc.shape

# <headingcell level=2>

# Save to netcdf

# <codecell>

from scipy.io import netcdf
f = netcdf.netcdf_file(fp+datestr+filesep+'20130823_4STAR_AOD.ncf','w')
f.history = 'Created by Samuel LeBlanc on 2015-01-30, 14:54 PST'
f.createDimension('UTC',len(utc))
f.createDimension('Wavelength',len(wvl))
utc_nc = f.createVariable('UTC', 'd',('UTC',))
utc_nc[:] = utc
utc_nc.units = 'Hours since midnight'
aod_nc = f.createVariable('AOD', 'd',('UTC','Wavelength'))
aod_nc[:] = aod
aod_nc.units = 'unitless'
wavelength_nc = f.createVariable('Wavelength', 'd', ('Wavelength',))
wavelength_nc[:] = wvl
wavelength_nc.units = 'microns'
lat_nc = f.createVariable('Latitude', 'd',('UTC',))
lat_nc[:] = lat[:,0]
lat_nc.units = 'Degrees North'
lon_nc = f.createVariable('Longitude', 'd',('UTC',))
lon_nc[:] = lon[:,0]
lon_nc.units = 'Degrees'
sza_nc = f.createVariable('SZA', 'd',('UTC',))
sza_nc[:] = sza[:,0]
sza_nc.units = 'Degrees'
f.close()

# <codecell>

fo = netcdf.netcdf_file(fp+datestr+filesep+'20130823_4STAR_AOD.ncf','r')

# <codecell>

g = fo.variables['AOD']

# <codecell>

print g.units

# <codecell>

fp+datestr+filesep+'20130823_4STAR_AOD.ncf'

# <headingcell level=2>

# Prepare for testing PCA analysis of raw spectra

# <codecell>

from sklearn.decomposition import PCA

# <codecell>

if 'rate' in locals():
    print rate.shape
    rate_fl = where(np.isfinite(rate[:,0]))
    wvl = star['w'][0]
else:
    print raw_vis.shape
    raw_fl = where(np.isfinite(raw_vis[:,0]))
    wvl = np.genfromtxt('C:\\Users\\sleblan2\\Research\\4STAR\\vis_4STAR_wvl.dat')

# <codecell>

wvl.shape

# <markdowncell>

# Get the log of the rate

# <codecell>

ratel = np.log(rate)

# <codecell>

print ratel.shape

# <codecell>

print isfinite(ratel[:,400]).all()

# <codecell>

wflrate = [isfinite(ratel).any(0)][0]
print where(~wflrate)
wflrate = range(4,1530)

# <codecell>

flrate = [isfinite(ratel[:,wflrate]).all(1)][0]
print flrate.shape
print where(~flrate)

# <codecell>

ratel1 = ratel[flrate,:]
print ratel1.shape
ratel2 = ratel1[:,wflrate]
print ratel2.shape
if len(wvl) != len(wflrate):
    wvl = wvl[wflrate]

# <codecell>

print isfinite(ratel2).all()

# <codecell>

if 'rate' in locals():
    ratel1 = ratel[flrate,:]
    pca = PCA(n_components=20)
    pca.fit(ratel2)
    evals = pca.explained_variance_ratio_
else:
    pca = PCA(n_components=20)
    pca.fit(raw_vis)
    evals = pca.explained_variance_ratio_

# <codecell>

plt.plot(evals)
plt.plot(evals.cumsum(),'r')

# <codecell>

print evals.cumsum()

# <codecell>

print evals

# <codecell>

coms = pca.components_

# <codecell>

print coms.shape

# <codecell>

fig3,ax3 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax3 = ax3.ravel()
for i in range(10):
    ax3[i].plot(wvl, coms[i,:])
    ax3[i].set_title('pca %d' % i)
    ax3[i].grid()
    #ax3[i].set_xlim([0.3,1.75])
    if i >7:
        ax3[i].set_xlabel('Wavelength [$\mu$m]')
plt.savefig(fp+datestr+'_langley_pca_log.png',dpi=600)

# <codecell>

print coms[0,:]

# <codecell>


# <codecell>

pca2 = PCA(n_components=10)
pca2.fit(np.log(raw_vis))
evals2 = pca2.explained_variance_ratio_
coms2 = pca2.components_

# <codecell>

fig3,ax3 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax3 = ax3.ravel()
for i in range(10):
    ax3[i].plot(wvl, coms[i,:])
    ax3[i].set_title('pca %d' % i)
    ax3[i].grid()
    ax3[i].set_xlim([0.3,1.05])
    if i >7:
        ax3[i].set_xlabel('Wavelength [$\mu$m]')
plt.savefig(fp+datestr+'_langley_pca.png',dpi=600)

