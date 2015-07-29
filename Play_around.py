# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap

# <codecell>

import os

# <codecell>

os.path.exists

# <codecell>

os.path.exists('C:\\Users\\sleblan2\\Research\\python_codes\\')

# <codecell>

fu = np.array([[550.000 , 1.286162e-06 , 3.552874e+02 , 2.316014e+01 , 1.298847e-07 , 4.876075e+01 , 6.735294e+00], 
     [550.000 , 8.816101e+02 , 1.048253e+02 , 5.573513e+02 , 8.903053e+01 , 2.529635e+01 , 9.411757e+01], 
     [550.000 , 1.114990e+03 , 3.609851e-05 , 5.377008e+02 , 1.125986e+02 , 1.569951e-05 , 8.701024e+01]])

sbdart = np.array([[250.000 , 1.296362e-03 , 3.430699e+05 , 1.931734e+04 , 1.309148e-04 , 4.699350e+04 , 6.008736e+03],
          [250.000 , 8.569704e+05 , 1.053688e+05 , 5.385019e+05 , 8.654226e+04 , 2.511391e+04 , 9.108232e+04],
          [250.000 , 1.110150e+06 , 3.641137e-02 , 5.153732e+05 , 1.121100e+05 , 1.576419e-02 , 8.285177e+04]])

# <codecell>

fu.shape

# <codecell>

plt.plot(fu,'+-')
plt.plot(sbdart/1000.0,'x-')

# <codecell>

import scipy.io as sio

# <codecell>

ffp = 'C:\Users\sleblan2\Dropbox\CALIPSO_Ames\Jens_data_mat\\'

# <codecell>

ss = sio.loadmat(ffp+'DTland_DBland_ocean_right.mat')

# <codecell>

ss.keys()

# <codecell>

import matplotlib.pyplot as plt
import matplotlib.patches as patches
class DraggablePoint:
    lock = None #only one can be animated at a time
    def __init__(self, point):
        self.point = point
        self.press = None
        self.background = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.point.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.point.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.point.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.point.axes: return
        if DraggablePoint.lock is not None: return
        contains, attrd = self.point.contains(event)
        if not contains: return
        self.press = (self.point.center), event.xdata, event.ydata
        DraggablePoint.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.point.figure.canvas
        axes = self.point.axes
        self.point.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.point.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.point)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        if DraggablePoint.lock is not self:
            return
        if event.inaxes != self.point.axes: return
        self.point.center, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.point.center = (self.point.center[0]+dx, self.point.center[1]+dy)

        canvas = self.point.figure.canvas
        axes = self.point.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.point)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if DraggablePoint.lock is not self:
            return

        self.press = None
        DraggablePoint.lock = None

        # turn off the rect animation property and reset the background
        self.point.set_animated(False)
        self.background = None

        # redraw the full figure
        self.point.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.point.figure.canvas.mpl_disconnect(self.cidpress)
        self.point.figure.canvas.mpl_disconnect(self.cidrelease)
        self.point.figure.canvas.mpl_disconnect(self.cidmotion)

fig = plt.figure()
ax = fig.add_subplot(111)
drs = []
circles = [patches.Circle((0.32, 0.3), 0.03, fc='r', alpha=0.5),
               patches.Circle((0.3,0.3), 0.03, fc='g', alpha=0.5)]

for circ in circles:
    ax.add_patch(circ)
    dr = DraggablePoint(circ)
    dr.connect()
    drs.append(dr)

plt.show()

# <codecell>


