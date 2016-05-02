
# coding: utf-8

# In[4]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap


# In[1]:

import os


# In[8]:

os.path.exists


# In[7]:

os.path.exists('C:\\Users\\sleblan2\\Research\\python_codes\\')


# In[7]:

fu = np.array([[550.000 , 1.286162e-06 , 3.552874e+02 , 2.316014e+01 , 1.298847e-07 , 4.876075e+01 , 6.735294e+00], 
     [550.000 , 8.816101e+02 , 1.048253e+02 , 5.573513e+02 , 8.903053e+01 , 2.529635e+01 , 9.411757e+01], 
     [550.000 , 1.114990e+03 , 3.609851e-05 , 5.377008e+02 , 1.125986e+02 , 1.569951e-05 , 8.701024e+01]])

sbdart = np.array([[250.000 , 1.296362e-03 , 3.430699e+05 , 1.931734e+04 , 1.309148e-04 , 4.699350e+04 , 6.008736e+03],
          [250.000 , 8.569704e+05 , 1.053688e+05 , 5.385019e+05 , 8.654226e+04 , 2.511391e+04 , 9.108232e+04],
          [250.000 , 1.110150e+06 , 3.641137e-02 , 5.153732e+05 , 1.121100e+05 , 1.576419e-02 , 8.285177e+04]])


# In[8]:

fu.shape


# In[18]:

plt.plot(fu,'+-')
plt.plot(sbdart/1000.0,'x-')


# In[10]:

import scipy.io as sio


# In[12]:

ffp = 'C:\Users\sleblan2\Dropbox\CALIPSO_Ames\Jens_data_mat\\'


# In[13]:

ss = sio.loadmat(ffp+'DTland_DBland_ocean_right.mat')


# In[14]:

ss.keys()


# In[17]:

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


# In[1]:

import hdf5storage as hs


# In[2]:

fp = 'C:\Users\sleblan2\Research\\4STAR_codes\data_folder\OMIdata\\'


# In[3]:

ff = fp+'OMI-Aura_L2G-OMDOAO3G_2015m1104_v003-2015m1105t054135.he5'


# In[4]:

fp


# In[43]:

help(hs.read)


# In[46]:

import os


# In[48]:

os.path.exists(ff)


# In[5]:

ff


# In[50]:

hs.read(path=fp,filename='OMI-Aura_L2G-OMDOAO3G_2015m1104_v003-2015m1105t054135.he5')


# In[31]:

hs.h5py.is_hdf5(fp+'OMI-Aura_L2G-OMDOAO3G_2015m1112_v003-2015m1113t063226.he5')


# In[6]:

import tables


# In[36]:

ff


# In[26]:

fp


# In[40]:

swathname = 'ColumnAmountO3'


# In[14]:

cd 'C:\\Users\\sleblan2\\Research\\4STAR_codes\\data_folder\\OMIdata\\'


# In[29]:

hdf5ref = tables.openFile(u'C:\\Users\\sleblan2\\Research\\4STAR_codes\\data_folder\\OMIdata\\OMI-Aura_L2G-OMDOAO3G_2015m1118_v003-2015m1119t060112.he5','r')


# In[30]:

fieldname = "ColumnAmountO3"


# In[31]:

fp='C:\\Users\\sleblan2\\Research\\4STAR_codes\\data_folder\\OMIdata\\'
omiO3G = fp+'OMI-Aura_L2G-OMDOAO3G_2015m1118_v003-2015m1119t060112.he5'


# In[32]:

hdf5ref = tables.openFile(omiO3G, mode="r", rootUEP=fieldname)


# In[33]:

import h5py


# In[34]:

fg = h5py.File(ff,'r+')


# In[54]:

hh = fg['HDFEOS']['GRIDS']['ColumnAmountO3']['Data Fields']


# In[57]:

jj = hh['Latitude']


# In[67]:

jj.attrs.keys()


# In[70]:

jj.attrs.get(u'Units')


# In[49]:

g = fg['HDFEOS']['GRIDS']['ColumnAmountO3']['Data Fields'].keys()


# In[71]:

for k in g:
    print fg['HDFEOS']['GRIDS']['ColumnAmountO3']['Data Fields'][k].value


# In[73]:

from datetime import datetime


# In[75]:

datetime.now().timetuple().tm_yday


# In[76]:

fg.close()


# In[77]:

import tables


# In[78]:

tables.File(ff,'r+')


# In[ ]:



