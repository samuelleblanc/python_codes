# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     Animate_4STAR
# 
# Purpose:  
# 
#     Python script that is used to view 4STAR full spectral AOD
#     Animates the AOD over time
# 
# Calling Sequence:
# 
#     python Animate_4STAR.py
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
#   - 20130823starsunfinal.mat: or other 4STAR starsun file containing aod spectra

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
fp = 'C:\\Users\\sleblan2\\Research\\4STAR\\SEAC4RS\\'
datestr = '20130823'
filesep = '\\'

# <codecell>

print filesep

# <headingcell level=2>

# Load the 4STAR data and select time points

# <codecell>

star = sio.loadmat(fp+datestr+filesep+datestr+'starsunfinal.mat')
star.keys()

# <codecell>

star['utc'] = toutc(mat2py_time(star['t']))

# <codecell>

star['good'] = np.where((star['utc']>20.0667) & (star['utc']<20.4667) & (star['Str'].flatten()!=0) & (star['sat_time'].flatten()==0))[0]
len(star['good'])

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
    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            anim.save(f.name, fps=20, extra_args=['-vcodec', 'libx264'])
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
    return line,

# animation function.  This is called sequentially
def animate(i):
    line.set_data(star['w'][0], star['tau_aero'][star['good']])
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(star['good']), interval=10, blit=True)

# call our new function to display the animation
display_animation(anim)

# <codecell>


