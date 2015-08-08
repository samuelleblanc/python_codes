import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

class LineBuilder:
    def __init__(self, line,m=None):
        """
        Start the line builder, with line2d object as input, and opitonally the m from basemap object
        """
        self.line = line
        self.m = m
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.connect()
        self.line.axes.format_coord = self.format_position_simple
        self.press = None
        self.contains = False

    def connect(self):
        'Function to connect all events'
        self.cid_onpress = self.line.figure.canvas.mpl_connect('button_press_event', self.onpress)
        self.cid_onrelease = self.line.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.cid_onmotion = self.line.figure.canvas.mpl_connect('motion_notify_event',self.onmotion)
        self.cid_onkeypress = self.line.figure.canvas.mpl_connect('key_press_event',self.onkeypress)
        self.cid_onkeyrelease = self.line.figure.canvas.mpl_connect('key_release_event',self.onkeyrelease)

    def disconnect(self):
        'Function to disconnect all events (except keypress)'
        self.line.figure.canvas.mpl_disconnect(self.cid_onpress)
        self.line.figure.canvas.mpl_disconnect(self.cid_onrelease)
        self.line.figure.canvas.mpl_disconnect(self.cid_onmotion)
        self.line.figure.canvas.mpl_disconnect(self.cid_onkeyrelease)

    def onpress(self,event):
        'Function that enables either selecting a point, or creating a new point when clicked'
        #print 'click', event
        if event.inaxes!=self.line.axes: return
        tb = plt.get_current_fig_manager().toolbar
        if tb.mode!='': return
        self.contains, attrd = line.contains(event)
        if self.contains:
            print 'click is near point:',self.contains,attrd
            self.contains_index = attrd['ind']
            print 'index:', self.contains_index
            if not self.contains_index is 0:
                self.xy = self.xs[self.contains_index-1],self.ys[self.contains_index-1]
                self.line.axes.format_coord = self.format_position_distance
            else:
                self.line.axes.format_coord = self.format_position_simple
            self.line.axes.autoscale(enable=False)
            self.highlight_linepoint, = self.line.axes.plot(self.xs[self.contains_index],
                                                            self.ys[self.contains_index],'bo')
        else:
            self.xy = self.xs[-1],self.ys[-1]
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            self.line.axes.format_coord = self.format_position_distance
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()
        self.press = event.xdata,event.ydata                                    
        sys.stdout.write('moving:')
        sys.stdout.flush()
        
    def onrelease(self,event):
        'Function to set the point location'
        print 'release'#,event
        self.press = None
        self.line.axes.format_coord = self.format_position_simple
        if self.contains:
            hlight = self.highlight_linepoint.findobj()[0]
            while hlight in self.line.axes.lines:
                self.line.axes.lines.remove(hlight)
            self.contains = False
        self.line.figure.canvas.draw()

    def onmotion(self,event):
        'Function that moves the points to desired location'
        if event.inaxes!=self.line.axes: return
        if self.press is None: return
        sys.stdout.write("\r"+" moving: x=%2.5f, y=%2.5f" %(event.xdata,event.ydata))
        sys.stdout.flush()
        if self.contains:
            i = self.contains_index
            self.highlight_linepoint.set_data(event.xdata,event.ydata)
        else:
            i = -1
        self.xs[i] = event.xdata
        self.ys[i] = event.ydata
        self.line.set_data(self.xs,self.ys)
        self.line.figure.canvas.draw()

    def onkeypress(self,event):
        print 'pressed key',event.key,event.xdata,event.ydata
        if event.inaxes!=self.line.axes: return
        if (event.key=='s') | (event.key=='alt+s'):
            print 'Stopping interactive point selection'
            self.disconnect()
        if (event.key=='i') | (event.key=='alt+i'):
            print 'Starting interactive point selection'
            self.connect()
            self.line.axes.format_coord = self.format_position_simple
            self.press = None
            self.contains = False

    def onkeyrelease(self,event):
        #print 'released key',event.key
        if event.inaxes!=self.line.axes: return

    def format_position_simple(self,x,y):
        if self.m:
            return 'Lon=%.7f, Lat=%.7f'%(self.m(x, y, inverse = True))
        else:   
            return 'x=%2.5f, y=%2.5f' % (x,y)

    def format_position_distance(self,x,y):
        if self.m:
            x0,y0 = self.xy
            lon0,lat0 = self.m(x0,y0,inverse=True)
            lon,lat = self.m(x,y,inverse=True)
            r = mu.spherical_dist([lat0,lon0],[lat,lon])
            return 'Lon=%.7f, Lat=%.7f, d=%.2f km'%(lon,lat,r)
        else:
            x0,y0 = self.xy
            self.r = sqrt((x-x0)**2+(y-y0)**2)
            return 'x=%2.5f, y=%2.5f, d=%2.5f' % (x,y,self.r)

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

def pll(string):
    """
    pll for parse_lat_lon
    function that parses a string and converts it to lat lon values
    one space indicates seperation between degree, minutes, or minutes and seconds
    returns decimal degrees
    """
    if type(string) is float:
        return string
    import re
    n = len(string.split())
    str_ls = string.split()
    char_neg = re.findall("[SWsw]+",str_ls[-1])
    char_pos = re.findall("[NEne]+",str_ls[-1])
    if len(char_neg)>0:
        sign = -1
        cr = char_neg[0]
    elif len(char_pos)>0:
        sign = 1
        cr = char_pos[0]
    else:
        sign = 1
        cr = ''
    str_ls[-1] = str_ls[-1].strip(cr)
    deg = float(str_ls[0])*sign
    deg_m = 0.0
    for i in range(n-1,0,-1):
        deg_m = (deg_m + float(str_ls[i])/60.0)/60.0
    return deg+deg_m
