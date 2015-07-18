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

def build_basemap(lower_left=[-20,-30],upper_right=[20,10],ax=plt.gca()):
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

class LineDrawer(object):
    lines = []
    pos = []
    def __init__(self,m):
        self.m = m
    def draw_line(self):
        ax = plt.gca()
        xy = plt.ginput(2)

        x = [p[0] for p in xy]
        y = [p[1] for p in xy]
        line = plt.plot(x,y)
        ax.figure.canvas.draw()
        
        self.lines.append(line)
        self.pos.append(self.m(x,y,inverse=True))

# <codecell>

m = build_basemap()
ax = plt.gca()
ax.format_coord = format_coord

# <codecell>

type(m)

# <codecell>

ld = LineDrawer(m)
ld.draw_line()

# <codecell>

ld.pos

# <codecell>

ld.m()

# <codecell>

ld.draw_line()

# <codecell>

ld.lines

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

def Create_excel():
    """
    Purpose:
    
        Program that creates the link to an excel file
        Starts and populates the first line and titles of the excel workbook
        
    Inputs:
    
        none
        
    Outputs:
    
        wb: workbook instance 
        
    Dependencies:
    
        xlwings
        Excel (win or mac)
        
    Required files:
    
        none
        
    Example:
    
        ...
 
    Modification History:
    
        Written: Samuel LeBlanc, 2015-07-15, Santa Cruz, CA
        
    """
    from xlwings import Workbook, Sheet, Range, Chart
    wb = Workbook()
    Sheet(1).name = 'Flight Path'
    Range('A1').value = ['WP','Lat\n[+-90]','Lon\n[+-180]',
                         'Speed\n[m/s]','delayT\n[min]','Altitude\n[m]',
                         'CumLegT\n[hh:mm]','UTC\n[hh:mm]','LocalT\n[hh:mm]',
                         'LegT\n[hh:mm]','Dist\n[km]','CumDist\n[km]',
                         'Dist\n[nm]','CumDist\n[nm]','Speed\n[kt]',
                         'Altitude\n[kft]','Comments']
    top_line = Range('A1').horizontal
    address = top_line.get_address(False,False)
    import sys
    if sys.platform.startswith('win'):
        from win32com.client import Dispatch
        xl = Dispatch("Excel.Application")
        xl.ActiveWorkbook.Windows(1).SplitColumn = 01
        xl.ActiveWorkbook.Windows(1).SplitRow = 01
        xl.Range(address).Font.Bold = True
    top_line.autofit()
    Range('A2').value = arange(50).reshape((50,1))
    return wb

# <codecell>

from xlwings import Workbook

# <codecell>

wb.Close()

# <codecell>

wb = Create_excel()

# <codecell>


