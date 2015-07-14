# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
%matplotlib tk
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap

# <codecell>

def build_basemap(lower_left=[-20,-25],upper_right=[15,0],ax=plt.gca()):
    """
    First try at a building of the basemap with a 'stere' projection
    Must put in the values of the lower left corner and upper right corner (lon and lat)
    
    Defaults to draw 8 meridians and parallels
    """
    m = Basemap(projection='stere',lon_0=(upper_right[0]+lower_left[0]),lat_0=(upper_right[1]+lower_left[1]),
            llcrnrlon=lower_left[0], llcrnrlat=lower_left[1],
            urcrnrlon=upper_right[0], urcrnrlat=upper_right[1],resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(lower_left[0],upper_right[0],8).astype(int),labels=[0,0,0,1])
    m.drawparallels(np.linspace(lower_left[1],upper_right[1],8).astype(int),labels=[1,0,0,0])
    return m

# <codecell>

def format_coord(x, y):
    return 'Lon=%.4f, Lat=%.4f'%(m(x, y, inverse = True))

# <codecell>

m = build_basemap()
ax = plt.gca()
ax.format_coord = format_coord

# <codecell>

fig = plt.gcf()

# <codecell>

def on_resize(event):
    ax.set_title('Resize, button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        event.button, event.x, event.y, event.xdata, event.ydata))

# <codecell>

fig.canvas.mpl_connect('resize_event',on_resize)

# <codecell>

fig.canvas.mpl_connect('button_release_event',define_event)

# <codecell>

fig.canvas.mpl_disconnect('resize_event')

# <codecell>

print 'yes'

# <codecell>

def define_event(event):
    for ax in event.canvas.figure.axes:
        print ax.get_navigate_mode()

# <codecell>

import sys
from PyQt4 import QtGui, QtCore

# <codecell>

plt.figure()
plt.plot([0,1,2,3,4],[1,2,3,4,8])
plt.show()

# <codecell>

plt.show()

# <codecell>

%matplotlib

# <codecell>

%matplotlib qt

# <codecell>


