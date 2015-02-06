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
    raw = star['rawcorr'][star['good'],:]
    m_aero = star['m_aero'][star['good']]
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

video(fp+datestr+filesep+'AOD_sp_anim.m4v','x-m4v')

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
rate_m = ratel/m_aero

# <codecell>

print rate_m.shape
print m_aero.shape
print ratel.shape

# <codecell>

print ratel.shape
print isfinite(ratel[:,400]).all()

# <codecell>

wflrate = [isfinite(ratel).any(0)][0]
print where(~wflrate)
wflrate = range(210,1020)+range(1050,1530)

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
    wvl = star['w'][0][wflrate]

# <codecell>

print isfinite(ratel2).all()

# <codecell>

raw1 = raw[flrate,:]
raw2 = raw1[:,wflrate]
rate_m1 = rate_m[flrate,:]
rate_m2 = rate_m1[:,wflrate]

# <codecell>

if 'rate' in locals():
    ratel1 = ratel[flrate,:]
    pca = PCA(n_components=20)
    pca.fit(ratel2)
    evals = pca.explained_variance_ratio_
    pca_r = PCA(n_components=20)
    pca_r.fit(raw2)
    evals_r = pca_r.explained_variance_ratio_
    pca_m = PCA(n_components=20)
    pca_m.fit(rate_m2)
    evals_m = pca_m.explained_variance_ratio_
else:
    pca = PCA(n_components=20)
    pca.fit(raw_vis)
    evals = pca.explained_variance_ratio_

# <codecell>

plt.plot(evals_r)
plt.plot(evals_r.cumsum(),'r')

# <codecell>

plt.plot(evals)
plt.plot(evals.cumsum(),'r')

# <codecell>

print evals.cumsum()
print evals
coms = pca.components_
print coms.shape
coms_r = pca_r.components_
coms_m = pca_m.components_

# <codecell>

fig2,ax2 = plt.subplots(5,2,figsize=(15,11))
ax2 = ax2.ravel()
for i in range(10):
    if i ==0 :
        ax2[i].plot(wvl, coms_r[i,:]*(-1.0))
    else:
        ax2[i].plot(wvl, coms_r[i,:])
    ax2[i].set_title('pca %d' % i)
    ax2[i].grid()
    #ax3[i].set_xlim([0.3,1.75])
    if i >7:
        ax2[i].set_xlabel('Wavelength [$\mu$m]')
plt.suptitle('Raw counts')
plt.tight_layout()
plt.savefig(fp+datestr+'_langley_pca_raw.png',dpi=600)

# <codecell>

for i,j in list(enumerate(wvl)): 
    print i,j

# <codecell>

fig3,ax3 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax3 = ax3.ravel()
yli = ([-0.2,0.00001],
       [-0.05,0.02],
       [-0.045,0.015],
       [-0.05,0.03],
       [-0.01,0.01],
       [-0.03,0.04],
       [-0.03,0.02],
       [-0.02,0.02],
       [-0.025,0.009],
       [-0.025,0.04])
for i in range(10):
    ax3[i].plot(wvl, coms[i,:])
    ax3[i].set_title('pca %d' % i)
    ax3[i].grid()
    #ax3[i].set_xlim([0.3,1.75])
    if i < 10: 
        ax3[i].set_ylim(yli[i])
    if i >7:
        ax3[i].set_xlabel('Wavelength [$\mu$m]')
plt.savefig(fp+datestr+'_langley_pca_log_zoom.png',dpi=600)

# <codecell>

fig4,ax4 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax4 = ax4.ravel()
for i in range(10):
    ax4[i].plot(wvl, coms[i,:])
    ax4[i].set_title('pca %d' % i)
    ax4[i].grid()
    if i >7:
        ax4[i].set_xlabel('Wavelength [$\mu$m]')
plt.savefig(fp+datestr+'_langley_pca_log.png',dpi=600)

# <codecell>

fig5,ax5 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax5 = ax5.ravel()
yli = ([-0.025,0.0],
       [-0.030,-0.020],
       [-0.010,0.01],
       [0.015,0.03],
       [-0.03,0.03],
       [0.00,0.013],
       [-0.01,0.05],
       [-0.015,0.01],
       [-0.005,0.008],
       [-0.005,0.002])
for i in range(10):
    ax5[i].plot(wvl, coms[i,:])
    ax5[i].set_title('pca %d' % i)
    ax5[i].grid()
    ax5[i].set_xlim([0.3,0.65])
    if i < 10: 
        ax5[i].set_ylim(yli[i])
    if i >7:
        ax5[i].set_xlabel('Wavelength [$\mu$m]')
plt.savefig(fp+datestr+'_langley_pca_log_vis_430nm.png',dpi=600)

# <codecell>

print coms[0,:]

# <codecell>

fig6,ax6 = plt.subplots(5,2,figsize=(15,11))
plt.tight_layout()
ax6 = ax6.ravel()
yli = ([-0.05,0.01],
       [-0.045,0.015],
       [-0.04,0.02],
       [-0.02,0.015],
       [-0.04,0.04],
       [-0.09,0.05],
       [-0.05,0.08],
       [-0.03,0.02],
       [-0.025,0.01],
       [-0.015,0.015])
for i in range(10):
    ax6[i].plot(wvl, coms_m[i,:])
    ax6[i].set_title('pca %d' % i)
    ax6[i].grid()
    #ax3[i].set_xlim([0.3,1.75])
    if i < 10: 
        ax6[i].set_ylim(yli[i])
    if i >7:
        ax6[i].set_xlabel('Wavelength [$\mu$m]')
plt.savefig(fp+datestr+'_langley_pca_log_airmass.png',dpi=600)

# <headingcell level=2>

# Use PCA to decompose and recompose the spectra

# <codecell>

pca_r_low = PCA(n_components=0.999)
pca_r_low.fit(raw2)
evals_r_low = pca_r_low.explained_variance_ratio_
print pca_r_low.n_components

# <codecell>

traw2 = pca_r_low.transform(raw2)
newraw2 = pca_r_low.inverse_transform(traw2)

# <codecell>

print raw2.shape
scores_r = pca_r.transform(raw2)
print scores_r.shape
new_raw2 = pca_r.inverse_transform(scores_r)
print coms_r.shape
print new_raw2.shape

# <markdowncell>

# Run PCA on rate_aero

# <codecell>

rate1 = rate[flrate,:]
rate2 = rate1[:,wflrate[1:]]
print np.isfinite(rate2).shape
print np.isfinite(rate2).all()
print rate2.shape
print np.isfinite(rate2[0,:]).all()

# <codecell>

pp = np.mean(rate2,0)

# <codecell>

print pp.shape

# <codecell>

print where(~np.isfinite(pp))

# <codecell>

print pp

# <codecell>

pca_rate = PCA(n_components=5,whiten=True)
pca_rate.fit(rate2)
print pca_rate.n_components

# <codecell>

trate2 = pca_rate.transform(rate2)
newrate2 = pca_rate.inverse_transform(trate2)

# <codecell>

print rate2.shape
print trate2.shape
print newrate2.shape

# <codecell>

plt.plot(rate2[:,500])
plt.plot(newrate2[:,500],'r')

# <codecell>

plt.plot(rate2[:,900])
plt.plot(newrate2[:,900],'r')

# <codecell>

plt.plot(rate2[500,:])
plt.plot(newrate2[500,:],'r')
coms_r2 = pca_rate.components_
plt.plot()

# <codecell>

print wvl[180:]

# <codecell>


