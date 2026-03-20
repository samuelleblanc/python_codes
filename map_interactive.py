# map_interactive software
# to use in combination with moving_lines
# Copyright 2015, Samuel LeBlanc
try:
    from mpl_toolkits.basemap import Basemap
except:
    pass
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    #import matplotlib.pyplot as plt
except:
    pass
import numpy as np
import sys
import re
import copy
from matplotlib.colors import is_color_like

try:
    import map_interactive as mi
    from map_utils import spherical_dist,equi,shoot,bearing
except ModuleNotFoundError:
    from . import map_interactive as mi
    from .map_utils import spherical_dist,equi,shoot,bearing
import time

class LineBuilder:
    """
    Purpose:
        The main interaction between a basemap plot and the functions to create clickable links
        Includes any method to actually plot lines/annotations/figures on the map
    Inputs: (at init)
        line class from a single plot
        m: basemap base class
        ex: excel_interface class
        verbose: (default False) writes out comments along the way
        tb: toolbar instance to use its interactions.
        blit: (optional, default to True) if set, then uses the blit techinque to redraw
    Outputs:
        LineBuilder class 
    Dependencies:
        numpy
        map_utils
        sys
        Basemap
        copy
    Required files:
        kml files for sat tracks
    Example:
        ...
    Modification History:
        Written: Samuel LeBlanc, 2015-08-07, Santa Cruz, CA
        Modified: Samuel LeBlanc, 2015-08-21, Santa Cruz, CA
                 - added new plotting with range circles
        Modified: Samuel LeBlanc, 2015-09-14, NASA Ames, Santa Cruz, CA
                 - added new points, move points from dialog windows
                 - bug fixes
        Modified: Samuel LeBlanc, 2015-09-15, NASA Ames, CA
                 - added handling of blit draw techinque to get gains in speed when drawing
        Modified: Samuel LeBlanc, 2016-06-22, NASA Ames, CA
                 - added flt_modules functionality
        Modified: Samuel LeBlanc, 2016-07-12, WFF, VA
                 - added plotting of aeronet data on map. Taken from internet.
        Modified: Samuel LeBlanc, 2016-07-22, Santa Cruz, CA
                 - added varying circle sizes depending on the size of the plotted map.
                 - added links to toolbar buttons to enable refreshing when zooming or panning, back forward
                 - change plotting of sat tracks (added lines)
                 - fixed legend disappear when showing aeronet and geos
        Modified: Samuel LeBlanc, 2016-08-01, WFF, VA
                 - fixed redraw bug where after a redraw, the corners don't match and cause a problem to WMS loading.
        Modified: Samuel LeBlanc, 2016-08-30, Swakopmund, Namibia
                 - fixed the upside down bocachica image
        Modified: Samuel LeBlanc, 2016-08-31, Swakopmund, Namibia
                 - fixed the bocachica plotting issue
                 - fixed some issue with time indications on wms plotting
        Modified: Samuel LeBlanc, 2016-09-19, Santa Cruz, CA
                 - made the legend draggable
        Modified: Samuel LeBlanc, 2017-09-15, St-John's, NL, Canada
                 - fixed issue with map pickles not having the dp variable.
        Modified: Samuel LeBlanc, 2019-06-03, Santa Cruz, CA
                 - fix linepicker issue with newer matplotlib. 
        Modified: Samuel LeBlanc, 2021-04-09, Santa Cruz, CA
                 - added continuity of meridians and parallels for problematic issues - bug fix in zooming out. 
        Modified: Samuel LeBlanc, 2021-10-21, Santa Cruza, CA
                - Migrated from python 2 to python 3.9, with cartopy and new version of xlwings.
                 
    """
    def __init__(self, line,m=None,ex=None,verbose=False,tb=None, blit=True):
        """
        Start the line builder, with line2d object as input,
        and opitonally the m from basemap object,
        Optionally the ex, dict_position class from the excel_interface,
            for interfacing with Excel spreadsheet
        
        """
        self.line = line
        self.line_arr = []
        self.line_arr.append(line)
        self.iactive = 0
        self.m = m
        self.ex = ex
        self.ex_arr = []
        self.ex_arr.append(ex)
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        if self.m:
            self.lons,self.lats = self.m.convert_latlon(self.xs,self.ys) #self.lons,self.lats = self.m(self.xs,self.ys,inverse=True)
            self.large = self.m.large
            if not self.m.use_cartopy:
                self.par = self.m.par
                self.mer = self.m.mer
        self.connect()
        self.line.axes.format_coord = self.format_position_simple
        self.press = None
        self.contains = False
        self.labelsoff = False
        self.circlesoff = False
        self.moving = False
        self.lbl = None
        self.tb = tb
        self.verbose = verbose
        self.blit = blit
        self.firstrun = True
        try:
            self.get_bg()
        except:
            import pdb; pdb.set_trace()

            

    def connect(self):
        'Function to connect all events'
        self.cid_onpress = self.line.figure.canvas.mpl_connect(
            'button_press_event', self.onpress)
        self.cid_onrelease = self.line.figure.canvas.mpl_connect(
            'button_release_event', self.onrelease)
        self.cid_onmotion = self.line.figure.canvas.mpl_connect(
            'motion_notify_event',self.onmotion)
        self.cid_onkeypress = self.line.figure.canvas.mpl_connect(
            'key_press_event',self.onkeypress)
        self.cid_onkeyrelease = self.line.figure.canvas.mpl_connect(
            'key_release_event',self.onkeyrelease)
        self.cid_onfigureenter = self.line.figure.canvas.mpl_connect(
            'figure_enter_event',self.onfigureenter)
        self.cid_onaxesenter = self.line.figure.canvas.mpl_connect(
            'axes_enter_event',self.onfigureenter)
        self.cid_onzoomscroll = self.line.figure.canvas.mpl_connect('scroll_event',self.zoom_fun)

    def disconnect(self):
        'Function to disconnect all events (except keypress)'
        self.line.figure.canvas.mpl_disconnect(self.cid_onpress)
        self.line.figure.canvas.mpl_disconnect(self.cid_onrelease)
        self.line.figure.canvas.mpl_disconnect(self.cid_onmotion)
        self.line.figure.canvas.mpl_disconnect(self.cid_onkeyrelease)
        self.line.figure.canvas.mpl_disconnect(self.cid_onfigureenter)
        self.line.figure.canvas.mpl_disconnect(self.cid_onaxesenter)
        self.line.figure.canvas.mpl_disconnect(self.cid_onzoomscroll)

    def onpress(self,event):
        'Function that enables either selecting a point, or creating a new point when clicked'
        #print('click', event)
        if event.inaxes!=self.line.axes: return
        if self.tb.mode!='': return
        if self.moving: return
        self.set_alt0 = False
        self.contains_point = False
        self.contains, attrd = self.line.contains(event)
        if self.line.labels_points:
            for lp in self.line.labels_points:
                try:
                    self.contains_point, attrd_point = lp.contains(event)
                except:
                    import pdb; pdb.set_trace()
                if self.contains_point:
                    xs,ys = lp.get_xdata(),lp.get_ydata()
        if self.contains:
            if self.verbose:
                print('click is near point:',self.contains,attrd)
            self.contains_index = attrd['ind']
            if len(self.contains_index)>1:
                self.contains_index = int(self.contains_index[-1])
            else:
                self.contains_index = int(self.contains_index)
            if self.verbose:
                print('index:%i'%self.contains_index)
            if self.contains_index != 0:
                self.xy = self.xs[self.contains_index-1],self.ys[self.contains_index-1]
                self.line.axes.format_coord = self.format_position_distance
                self.line.axes.autoscale(enable=False)
                self.highlight_linepoint, = self.line.axes.plot(self.xs[self.contains_index],
                                                            self.ys[self.contains_index],'bo',zorder=40)
                self.draw_canvas(extra_points=[self.highlight_linepoint])
            else:
                self.line.axes.format_coord = self.format_position_simple
                self.xy = self.xs[-1],self.ys[-1]
                self.xs.append(self.xs[self.contains_index])
                self.ys.append(self.ys[self.contains_index])
                if self.m:
                    lo,la = self.m.convert_latlon(self.xs[self.contains_index],self.ys[self.contains_index]) #self.m(self.xs[self.contains_index],self.ys[self.contains_index],inverse=True)
                    self.lons.append(lo)
                    self.lats.append(la)
                self.set_alt0 = True
                self.contains = False
            ilola = self.contains_index
        elif self.contains_point:
            point_contains_index = attrd_point['ind']
            self.xy = self.xs[-1],self.ys[-1]
            #print(self.xy)
            self.line.axes.format_coord = self.format_position_distance
            self.line.axes.autoscale(enable=False)
            self.highlight_linepoint, = self.line.axes.plot(xs[point_contains_index],
                                                        ys[point_contains_index],'bo',zorder=40)
            self.draw_canvas(extra_points=[self.highlight_linepoint])
            if self.m:
                lo,la = self.m.convert_latlon(xs[point_contains_index],ys[point_contains_index]) #self.m(xs[point_contains_index],ys[point_contains_index],inverse=True)
                self.lons.append(lo)
                self.lats.append(la)
            self.xs.append(xs[point_contains_index])
            self.ys.append(ys[point_contains_index])
            ilola = -2
        else:
            self.xy = self.xs[-1],self.ys[-1]
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            #print(self.xy,event.xdata,event.ydata,self.xs[-1],self.ys[-1])
            if self.m:
                lo,la = self.m.convert_latlon(event.xdata,event.ydata) #self.m(event.xdata,event.ydata,inverse=True)
                self.lons.append(lo)
                self.lats.append(la)
            self.line.axes.format_coord = self.format_position_distance
            ilola = -2
        self.line.set_data(self.xs, self.ys)
        if self.ex:
            try:
                self.azi = self.ex.azi[ilola+1]
            except IndexError:
                self.azi = self.ex.azi[ilola]
        else:
            self.azi = None
        self.line.range_circles,self.line.range_cir_anno = self.plt_range_circles(self.lons[ilola],self.lats[ilola],azi=self.azi)
        self.draw_canvas(extra_points=[self.line.range_circles,self.line.range_cir_anno])
        self.press = event.xdata,event.ydata
        if self.verbose:
            sys.stdout.write('moving:')
            sys.stdout.flush()
        
    def onrelease(self,event):
        'Function to set the point location'
        
        if self.verbose:
            print('release')#,event
        if self.moving: return

        if (self.tb.mode == 'zoom rect') | (self.tb.mode == 'pan/zoom') | (self.tb.mode!=''):
            #self.redraw_pars_mers()
            
            self.m.llcrnrlat,self.m.urcrnrlat = self.m.get_ylim() 
            self.m.llcrnrlon,self.m.urcrnrlon = self.m.get_xlim()
            self.get_bg()
            return
        elif self.tb.mode!='':
            return
        
        if event.inaxes!=self.line.axes: return
        self.press = None
        self.line.axes.format_coord = self.format_position_simple
        
        if self.contains:
            hlight = self.highlight_linepoint.findobj()[0]
            while hlight in self.line.axes.lines:
                self.line.axes.lines.remove(hlight)
            self.contains = False
            if self.ex:
                self.ex.mods(self.contains_index,
                             self.lats[self.contains_index],
                             self.lons[self.contains_index])
                self.ex.calculate()
                self.ex.write_to_excel()
        elif self.contains_point:
            hlight = self.highlight_linepoint.findobj()[0]
            while hlight in self.line.axes.lines:
                self.line.axes.lines.remove(hlight)
            self.contains_point = False
            if self.ex:
                self.ex.appends(self.lats[-1],self.lons[-1])    
                self.ex.calculate()
                self.ex.write_to_excel()           
            self.xy = self.m.invert_lonlat(self.lons[-1],self.lats[-1])
        else:
            if self.ex:
                if self.set_alt0:
                    self.ex.appends(self.lats[-1],self.lons[-1],alt=self.ex.alt[0])
                else:
                    self.ex.appends(self.lats[-1],self.lons[-1])
                self.ex.calculate()
                self.ex.write_to_excel()
        for lrc in self.line.range_circles:
            lrc.remove()
        for arc in self.line.range_cir_anno:
            try:
                arc.remove()
            except:
                continue
        self.update_labels(nodraw=True)
        self.line.figure.canvas.draw()
        self.get_bg()
            

    def onmotion(self,event):
        'Function that moves the points to desired location'
        if event.inaxes!=self.line.axes: return
        if self.press is None: return
        if self.tb.mode!='': return
        if self.moving: return
        if self.verbose:
            sys.stdout.write("\r"+" moving: x=%2.5f, y=%2.5f" %(event.xdata,event.ydata))
            sys.stdout.flush()
        if self.contains:
            i = self.contains_index
            self.highlight_linepoint.set_data(event.xdata,event.ydata)
        else:
            i = -1
        self.xs[i] = event.xdata
        self.ys[i] = event.ydata
        if self.m:
            self.lons[i],self.lats[i] = self.m.convert_latlon(event.xdata,event.ydata)# self.m(event.xdata,event.ydata,inverse=True)
        self.line.set_data(list(self.xs),list(self.ys))
        if self.contains:
            self.draw_canvas(extra_points=[self.highlight_linepoint,self.line.range_circles,self.line.range_cir_anno])
        else:
            self.draw_canvas(extra_points=[self.line.range_circles,self.line.range_cir_anno])       

    def onkeypress(self,event):
        'function to handle keyboard events'
        if self.verbose:
            print('pressed key',event.key,event.xdata,event.ydata)
        if event.inaxes!=self.line.axes: return
        if (event.key=='s') | (event.key=='alt+s'):
            print('Stopping interactive point selection')
            self.disconnect()
        if (event.key=='i') | (event.key=='alt+i'):
            print('Starting interactive point selection')
            self.connect()
            self.line.axes.format_coord = self.format_position_simple
            self.press = None
            self.contains = False

    def onkeyrelease(self,event):
        'function to handle keyboard releases'
        #print('released key',event.key)
        if event.inaxes!=self.line.axes: return

    def onfigureenter(self,event):
        'event handler for updating the figure with excel data'
        if self.firstrun:
            self.get_bg(redraw=True)
            self.firstrun = False
        if self.moving: return
        self.tb.set_message('Recalculating ...')
        if self.verbose:
            print('entered figure')#, event
        if self.ex:
            #self.ex.switchsheet(self.iactive)
            self.ex.check_xl()
            self.lats = list(self.ex.lat)
            self.lons = list(self.ex.lon)
            if self.m:
                x,y = self.m.invert_lonlat(self.ex.lon,self.ex.lat) #self.m(self.ex.lon,self.ex.lat)
                try:
                    self.xs = list(x)
                    self.ys = list(y)
                except TypeError as ei:
                    import pdb; pdb.set_trace()
                self.line.set_data(self.xs,self.ys)
                self.draw_canvas()
        self.update_labels()
        self.tb.set_message('Done Recalculating')
        self.line.axes.format_coord = self.format_position_simple
        
    def zoom_fun(self,event,base_scale=1.25):
        'For zooming with scroll_event'
        # from gist at https://gist.github.com/scott-vsi/522e756d636557ae8f1ef3cdb069cecd
        # make this figure the current figure 
        #plt._figure(ax.get_figure().number)
        # push the current view to define home if stack is empty
        #toolbar = ax.get_figure().canvas.toolbar # only set the home state
        #if self.tb._views.empty():
        #    self.tb.push_current()

        # get the current x and y limits
        cur_xlim = self.m.get_xlim()
        cur_ylim = self.m.get_ylim()
        xdata = event.xdata # get event x location
        ydata = event.ydata # get event y location
        # hat tip to @giadang
        # Get distance from the cursor to the edge of the figure frame
        x_left = xdata - cur_xlim[0]
        x_right = cur_xlim[1] - xdata
        y_top = ydata - cur_ylim[0]
        y_bottom = cur_ylim[1] - ydata
        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1/base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
            #print event.button
        # set new limits
        self.m.set_xlim([xdata - x_left*scale_factor,
                     xdata + x_right*scale_factor])
        self.m.set_ylim([ydata - y_top*scale_factor,
                     ydata + y_bottom*scale_factor])

        # push the current view limits and position onto the stack
        # REVIEW perhaps this should participate in a drag_zoom like
        # matplotlib/backend_bases.py:NavigationToolbar2:drag_zoom()
        self.tb.push_current()
        self.get_bg(redraw=True)
        #self.draw_canvas()#ax.figure.canvas.draw() # force re-draw
                
    def format_position_simple(self,x,y):
        'format the position indicator with only position'
        if self.m:
            if self.m.use_cartopy:
                return 'Lon=%.7f, Lat=%.7f'%self.m.convert_latlon(x,y)
            else:
                return 'Lon=%.7f, Lat=%.7f'%(x,y) #self.m(x, y, inverse = True))
        else:   
            return 'x=%2.5f, y=%2.5f' % (x,y)

    def format_position_distance(self,x,y):
        'format the position indicator with distance from previous point'
        if self.m:
            x0,y0 = self.xy
            lon0,lat0 = self.m.convert_latlon(x0,y0) #self.m(x0,y0,inverse=True)
            lon,lat = self.m.convert_latlon(x,y)
            r = spherical_dist([float(lat0),float(lon0)],[float(lat),float(lon)])
            return 'Lon=%.7f, Lat=%.7f, d=%.2f km'%(lon,lat,r)
        else:
            x0,y0 = self.xy
            self.r = sqrt((x-x0)**2+(y-y0)**2)
            return 'x=%2.5f, y=%2.5f, d=%2.5f' % (x,y,self.r)
        
    def update_labels(self,nodraw=False):
        'method to update the waypoints labels after each recalculations'
        #import matplotlib as mpl
        #if mpl.rcParams['text.usetex']:
        #    s = '\#'
        #else:
        #    s = '#'
        s = '#'
        if self.ex:
            self.wp = self.ex.WP
        else:
            self.n = len(self.xs)
            self.wp = range(1,self.n+1)
        if self.lbl:
           for ll in self.lbl:
                try:
                    ll.remove()
                except:
                    continue
        if self.labelsoff:
            return
        vas = ['bottom', 'baseline', 'center', 'center_baseline', 'top']
        has = ['left', 'right', 'center']
        for i in self.wp:    
            if not self.lbl:
                self.lbl = [self.line.axes.annotate(s+'%i'%i,
                                                    (self.xs[i-1],self.ys[i-1]))]
            else:
                try:
                    if not self.xs[i-1]:
                        continue
                    self.lbl.append(self.line.axes.
                                    annotate(s+'%i'%i,(self.xs[i-1],self.ys[i-1]),ha=has[i%3],va=vas[i%5],zorder=45))
                except IndexError:
                    pass
        if not nodraw:
            self.line.figure.canvas.draw()

    def plt_range_circles(self,lon,lat,azi=None):
        'program to plot range circles starting from the last point selected on the map, with principal plane identified'        
        if self.circlesoff:
            return
        if self.large:
            diam = [50.0,100.0,200.0,500.0,1000.0]
        else:
            diam = [25.0,50.0,100.0,200.0,500.0]
        colors = ['lightgrey','lightgrey','lightgrey','lightsalmon','crimson']
        line = []
        an = []

        for i,d in enumerate(diam):
            ll, = equi(self.m,lon,lat,d,color=colors[i],transform=self.m.merc)
            line.append(ll)
            slon,slat,az = shoot(lon,lat,0.0,d)
            x,y = self.m.invert_lonlat(slon,slat) #self.m(slon,slat)
            ano = self.line.axes.annotate('%i km' %(d),(x,y),color='silver')
            an.append(ano)
        if azi:
            slon,slat,az = shoot(lon,lat,azi,diam[-1])
            mlon,mlat,az = shoot(lon,lat,azi+180.0,diam[-1])
            lazi1, = self.m.plot([slon],[slat],'--*',color='grey',markeredgecolor='#BBBB00',markerfacecolor='#EEEE00',markersize=20,transform=self.m.merc)
            lazi2, = self.m.plot([mlon,lon,slon],[mlat,lat,slat],'--',color='grey',transform=self.m.merc)
            line.append(lazi1)
            line.append(lazi2)
        # plot angle points on the line
        for deg in [0,90,180,270]:
            dlo,dla,az = shoot(lon,lat,deg,diam[-1])
            elo,ela,az = shoot(lon,lat,deg,diam[-1]*0.85)
            ll, = self.m.plot([elo,dlo],[ela,dla],'-',color='grey',transform=self.m.merc)
            line.append(ll)
            for dd in [22.5,45.0,67.5]:
                dlo,dla,az = shoot(lon,lat,deg+dd,diam[-1])
                elo,ela,az = shoot(lon,lat,deg+dd,diam[-1]*0.93)
                ll, = self.m.plot([elo,dlo],[ela,dla],'-',color='grey',transform=self.m.merc)
                line.append(ll)
        return line,an

    def makegrey(self):
        'Program to grey out the entire path'
        self.line.set_color('#AAAAAA')
        self.line.set_zorder(20)
        self.line.figure.canvas.draw()
        self.get_bg()
        
    def colorme(self,c):
        'Program to color the entire path'
        self.line.set_color(c)
        self.line.set_zorder(30)

    def newline(self):
        'Program to do a deep copy of the line object in the LineBuilder class'
        x,y = self.line.get_data()
        line_new, = self.m.plot(x[0],y[0],'o-',linewidth=self.line.get_linewidth()) 
        line_new.labels_points = self.line_arr[-1].labels_points
        self.line_arr.append(line_new)

    def removeline(self,i):
        'Program to remove one line object from the LineBuilder class'
        self.line_arr[i].set_data([],[])
        self.line_arr[i].remove()

    def addfigure_under(self,img,ll_lat,ll_lon,ur_lat,ur_lon,outside=False,text=None,alpha=0.5,name='None',**kwargs):
        'Program to add a figure under the basemap plot'
        try: 
            self.m.figure_under
        except AttributeError:
            self.m.figure_under = {}
        try:
            #import ipdb; ipdb.set_trace()
            if len(img.getbands())>3:
                def invert(image):
                    return image.point(lambda p: 255 - p)
                r0,g0,b0,a0 = img.split()
                img_p = Image.merge(img.mode,(r0,g0,b0,invert(a0)))
                self.m.figure_under[name] = self.m.imshow(img_p,origin='upper',\
                transform=kwargs.get('transform',self.m.proj),extent=[ll_lon,ur_lon,ll_lat,ur_lat])#,\
                #alpha=(1.0-np.rollaxis(np.array(img.getchannel('A')),1)[::-1,:]/255))#1.0-np.array(img.getchannel('A')).reshape(img.size[1],img.size[0])/255)
            else:
                self.m.figure_under[name] = self.m.imshow(img,origin='upper',\
                transform=kwargs.get('transform',self.m.proj),extent=[ll_lon,ur_lon,ll_lat,ur_lat])
            print('...figure {} has been addded. {}'.format(name,text))
        except:
            self.m.figure_under[name] = self.m.imshow(img,origin='upper',\
            transform=kwargs.get('transform',self.m.proj),extent=[ll_lon,ur_lon,ll_lat,ur_lat])
            print('...figure {} has been addded. {}'.format(name,text))
            

        if text:
            try: 
                self.m.figure_under_text
            except AttributeError:
                self.m.figure_under_text = {}
            try:
                self.m.figure_under_text[name] = self.m.ax.text(0.0,-0.15,text,transform=self.m.ax.transAxes,clip_on=False,color='grey')
            except:
                print('Problem adding text on figure, continuning...')
        self.line.figure.canvas.draw()
        self.get_bg()
        
    def addlegend_image_below(self,img):
        'Program to add a image legend to a new axis below the current axis'
        try: 
            self.m.legend_axis = self.line.figure.add_axes([0.1,0.0,0.5,0.1],anchor='NW',zorder=-1)
            self.m.legend_axis.imshow(img)
            self.m.legend_axis.axis('off')
        except:
            return False
        self.line.figure.canvas.draw()
        
    def redraw_pars_mers(self):
        'redraws the parallels and meridians based on the current geometry'
        ylim = self.line.axes.get_ylim()
        xlim = self.line.axes.get_xlim()

        #self.m.llcrnrlon = xlim[0]
        #self.m.llcrnrlat = ylim[0]
        #self.m.urcrnrlon = xlim[1]
        #self.m.urcrnrlat = ylim[1]
        round_to_2 = lambda x:(int(x/2)+1)*2
        round_to_5 = lambda x:(int(x/5)+1)*5
        self.large = True
        if (xlim[1]-xlim[0])<20.0:
            mer = np.arange(round_to_2(xlim[0]),round_to_2(xlim[1])+2,2)
            self.large = False
        else:
            mer = np.arange(round_to_5(xlim[0]),round_to_5(xlim[1])+5,5)
        if (ylim[1]-ylim[0])<20.0:
            par = np.arange(round_to_2(ylim[0]),round_to_2(ylim[1])+2,2)
            self.large = False
        else:
            par = np.arange(round_to_5(ylim[0]),round_to_5(ylim[1])+5,5)
        if len(mer)<2:
            mer = self.mer
            self.line.axes.set_xlim(self.m.orig_xlim)
        if len(par)<2:
            par = self.par
            self.line.axes.set_ylim(self.m.orig_ylim)
        try:
            mi.update_pars_mers(self.m,mer,par,lower_left=(xlim[0],ylim[0]))
        except:
            import traceback
            print('... Problem updating the parallels and meridians')
            traceback.print_exc()
            import pdb; pdb.set_trace()
        self.line.figure.canvas.draw()

    def newpoint(self,Bearing,distance,alt=None,last=True,feet=False,km=True,insert=False,insert_i=-1):
        """
        program to add a new point at the end of the current track with a bearing and distance, optionally an altitude
        if feet is set, altitude is defined in feet not the defautl meters.
        if km is set to True, distance is defined in kilometers (default), if False, it uses nautical miles (nm)
        update to use insert keyword for placing a point in the middle of the line, at the position of insert_i
        """
        if not km:
            dist = distance*2.02 #need to check this value.
        else:
            dist = distance
        newlon,newlat,baz = shoot(self.lons[insert_i],self.lats[insert_i],Bearing,maxdist=distance)
        if self.verbose:
            print('New points at lon: %f, lat: %f' %(newlon,newlat))
        if self.m:
            x,y = self.m.invert_lonlat(newlon,newlat) #self.m(newlon,newlat)
            if not insert:
                self.lons.append(newlon)
                self.lats.append(newlat)
            else:
                self.lons.insert(insert_i+1,newlon)
                self.lats.insert(insert_i+1,newlat)
        else:
            x,y = newlon,newlat
        if not insert:
            self.xs.append(x)
            self.ys.append(y)
        else:
            self.xs.insert(insert_i+1,x)
            self.ys.insert(insert_i+1,y)
        self.line.set_data(self.xs, self.ys)
        
        if self.ex:
            if alt:
                if feet:
                    alt = alt/3.28084
            if not insert:
                self.ex.appends(self.lats[-1],self.lons[-1],alt=alt)
            else:
                self.ex.inserts(insert_i+1,newlat,newlon,alt=alt)
            self.ex.calculate()
            self.ex.write_to_excel()
        self.update_labels()
        self.draw_canvas()

    def movepoint(self,i,Bearing,distance,last=False):
        'Program to move a point a certain distance and bearing'
        newlon,newlat,baz = shoot(self.lons[i],self.lats[i],Bearing,maxdist=distance)
        if self.m: 
            x,y = self.m.invert_lonlat(newlon,newlat) #self.m(newlon,newlat)
            self.lons[i] = newlon
            self.lats[i] = newlat
        else:
            x,y = newlon,newlat
        self.xs[i] = x
        self.ys[i] = y
        self.line.set_data(self.xs, self.ys)
        if self.ex: self.ex.mods(i,self.lats[i],self.lons[i])
        if last:
            if self.ex:
                self.ex.calculate()
                self.ex.write_to_excel()
            self.update_labels()
            self.draw_canvas()
            
    def parse_flt_module_file(self,filename):
        """
        Program that opens a macro file and parses it, runs the move commands on each line of the file
        Added ability to defined either altitude in feet or meters, or the length in km or nm
        
        """
        try:
            from gui import ask
        except ModuleNotFoundError:
            from .gui import ask
        import os
        # set up predefined values
        predef = ['azi','AZI','PP','pp']
        azi = self.ex.azi[-1]
        pp, AZI, PP = azi,azi,azi
        
        # load flt_module
        name = os.path.basename(filename).split('.')[0]
        f = open(filename,'r')
        s = f.readlines()
        # get variables in flt_module
        for l in s:
            if l.startswith('%'):
                names = [u.strip() for u in l.strip().strip('%').split(',')]
                defaults = []
                if any([True for n in names if '=' in n]):
                    defaults = [0.0]*len(names)
                    for j,n in enumerate(names):
                        if n.find('=')>=0:
                            try:
                                n,d = n.split('=')
                                defaults[j] = float(d)
                                names[j] = n
                            except:
                                pass
                            
        choice = ['feet','meters']
        choice2 = ['km','nm']
        # remove values predefined
        for p in predef:
            try:
                names.remove(p)
            except:
                pass
                
        # ask the user supply the variables
        vals = ask(names,choice=choice,choice_title='Altitude:',choice2=choice2,choice2_title='Distance:',
                   title='Enter values for {}'.format(name),defaults=defaults)
        v = vals.names_val
        if vals.choice_val == 'feet': 
            use_feet = True 
        else: 
            use_feet = False
        
        if vals.choice2_val == 'km': 
            use_km = True  #km
            factor = 1000.0
        else: 
            use_km = False #nm
            factor = 1000.0*2.02
            
        # make the values variables
        for i,n in enumerate(names):
            try:
                exec('{}={}'.format(n,vals.names_val[i]))
            except:
                print('problem for {}={}'.format(n,vals.names_val[i]))
        for l in s[1:-1]:
            if not (l.startswith('#') or l.startswith('%')):
                if l.lower().find('time')>=0: #check if there is a time slot instead of distance
                    try:
                        if len(l.split(','))>2:
                            bear,tim,alt = eval(l)
                        else:
                            bear,tim = eval(l)
                            alt = self.ex.alt[-1]
                    except Exception as ie:
                        print('error splitting lines in flt_module: {}'.format(name))
                    try:
                        dist = (tim*self.ex.calcspeed(self.ex.alt[-1],alt)*60.0)/factor # in distance units
                    except:
                        print('problem with distance calculation for {}'.format(l))
                    l = ','.join([l.split(',')[0],'{}'.format(dist),'{}'.format(alt)])
                try:
                    print(eval(l),use_feet,use_km)
                    self.newpoint(*eval(l),last=False,feet=use_feet,km=use_km)
                except:
                    print('problem with {}'.format(l))
        if not s[-1].startswith('#') or not s[-1].startswith('%'):
                try:
                    self.newpoint(*eval(s[-1]),feet=use_feet,km=use_km)
                except:
                    print('problem with last {}'.format(s[-1]))
        f.close()

    def get_bg(self,redraw=False):
        'program to store the canvas background. Used for blit technique'
        if redraw:
            self.line.figure.canvas.draw()
        self.m.llcrnrlat,self.m.urcrnrlat = self.m.get_ylim() 
        self.m.llcrnrlon,self.m.urcrnrlon = self.m.get_xlim()
        self.bg = self.line.figure.canvas.copy_from_bbox(self.line.axes.bbox)

    def draw_canvas(self,extra_points=[]):
        'Program to handle the blit technique or simply a redraw of the canvas'
        if self.blit:
            self.line.figure.canvas.restore_region(self.bg)
            self.line.axes.draw_artist(self.line)
            try:
                for p in extra_points:
                    if type(p) is list:
                        for px in p:
                           self.line.axes.draw_artist(px) 
                    else:
                        self.line.axes.draw_artist(p)
            except Exception as ie:
                print('exception occurred: %s' %ie)
            self.line.figure.canvas.blit(self.line.axes.bbox)
        else:
            self.line.figure.canvas.draw()
            
    def calc_move_from_rot(self,i,angle,lat0,lon0):
        """
        Program to calculate a bearing and angle to rotate the point i based on the lat lon points lat0 and lon0
        Returns a bearing and a distance of the point i to get the correct position
        """
        dist = spherical_dist([lat0,lon0],[self.lats[i],self.lons[i]])
        bearing_start = bearing([lat0,lon0],[self.lats[i],self.lons[i]])
        newlon,newlat,baz = shoot(lon0,lat0,bearing_start+angle,maxdist=dist)
        dist_end = spherical_dist([self.lats[i],self.lons[i]],[newlat,newlon])
        bearing_end = bearing([self.lats[i],self.lons[i]],[newlat,newlon])
        return bearing_end,dist_end        
        
    def calc_dist_from_each_points(self):
        """
        Program to run through the ex_arr waypoints, and calculate the distances to each of the different points
        returns a combined dict of the set of waypoints, utc, cumulative leg times, lat,lon,alt, distances to each other, and indices of the waypoints
        """
        from scipy import interpolate
        combined = {'utc':np.hstack([x.utc for x in self.ex_arr]),
                    'cumlegt':np.hstack([x.cumlegt for x in self.ex_arr]),
                    'index_ex_array':np.hstack([x.utc*0+i for i,x in enumerate(self.ex_arr)]).astype(int)}
        # now sort the combined vals
        i_sort = np.argsort(combined['utc'])
        for k in combined:
            combined[k] = combined[k][i_sort]
        # run through and calculate the positions of each val during these times.
        n_exs = len(self.ex_arr)
        combined['lats'] = np.zeros((len(combined['utc']),n_exs))
        combined['lons'] = np.zeros((len(combined['utc']),n_exs))
        for i,ex in enumerate(self.ex_arr):
            fx_lat = interpolate.interp1d(ex.utc,ex.lat,bounds_error=False)
            fx_lon = interpolate.interp1d(ex.utc,ex.lon,bounds_error=False)
            combined['lats'][:,i] = fx_lat(combined['utc'])
            combined['lons'][:,i] = fx_lon(combined['utc'])
        #now calculate the distance to each lats and lons, based on the index of the value at that moment
        combined['distances'] = np.zeros((len(combined['utc']),n_exs))
        for j,utc in enumerate(combined['utc']):
            for k in list(range(n_exs)):
                combined['distances'][j,k] = spherical_dist([combined['lats'][j,combined['index_ex_array'][j]],combined['lons'][j,combined['index_ex_array'][j]]],[combined['lats'][j,k],combined['lons'][j,k]])
        distances = []
        for i,ex in enumerate(self.ex_arr): 
            indices = np.where(combined['index_ex_array']==i)[0]
            distances.append(combined['distances'][indices,:].tolist())
        
        return combined,distances
        
def build_basemap(lower_left=[-20,-30],upper_right=[20,10],ax=None,fig=None,proj='cyl',profile=None,larger=True, use_cartopy=True):
    """
    First try at a building of the basemap with a 'stere' projection
    Must put in the values of the lower left corner and upper right corner (lon and lat)
    
    Defaults to draw 8 meridians and parallels

    Modified: Samuel LeBlanc, 2015-09-15, NASA Ames
            - added profile keyword that contains the basemap profile dict for plotting the corners
            - added programatic determination of basemap parallels and meridians
    Modified: Samuel LeBlanc, 2016-07-17, Santa Cruz, CA
            - added plotting of larger region, which is then resized, to have space to pan and zoom.
    """
    try:
        from map_interactive import pll
    except ModuleNotFoundError:
        from .map_interactive import pll
    import os
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    if profile:
        upper_right = [pll(profile['Lon_range'][1]),pll(profile['Lat_range'][1])]
        lower_left = [pll(profile['Lon_range'][0]),pll(profile['Lat_range'][0])]
    if profile:
        lat0 = pll(profile['Start_lat'])
        lon0 = pll(profile['Start_lon'])
    else:
        lat0 = (upper_right[1]+lower_left[1])/2.0
        lon0 = (upper_right[0]+lower_left[0])/2.0
    #print(proj,lat0,lon0)
    if larger:
        dp = 30
    else:
        dp = 0
    if not use_cartopy:
        if os.path.isfile('map_{}.pkl'.format(profile['Campaign'])):
            import pickle
            m = pickle.load(open('map_{}.pkl'.format(profile['Campaign']),'rb'))
            #print('doing it pickle style')
            m.ax = ax
                    
    if larger:
        dp = 30
    else:
        dp = 0
    if use_cartopy:
        #proj = ccrs.PlateCarree()
        proj = ccrs.__dict__[profile.get('proj','PlateCarree')]() #central_longitude=profile.get('central_longitude',0),central_latitude=profile.get('central_latitude',30))
        m = fig.add_subplot(111,projection=proj)
        m.proj = proj
        m.proj_name = profile.get('proj','PlateCarree')
        m.use_cartopy = True
    else:
        if proj == 'cyl':
            m = Basemap(projection=proj,lon_0=lon0,lat_0=lat0,
                llcrnrlon=lower_left[0]-dp, llcrnrlat=(lower_left[1]-dp if lower_left[1]-dp>-90 else -90),
                urcrnrlon=upper_right[0]+dp, urcrnrlat=(upper_right[1]+dp if upper_right[1]+dp<90 else 90),resolution='i',ax=ax)
        else:
            dp = 0
            m = Basemap(projection=proj,lon_0=lon0,lat_0=lat0,
                llcrnrlon=lower_left[0]-dp, llcrnrlat=lower_left[1]-dp,
                urcrnrlon=upper_right[0]+dp, urcrnrlat=upper_right[1]+dp,
                resolution='l',ax=ax)
            #m = Basemap(projection='npstere',lat_0=76.5,lon_0=-68.75,
            #    boundinglat=55,resolution='l',ax=ax)
        m.use_cartopy = False
        m.proj_name = proj
                
        m.artists = []
    if m.use_cartopy:
        m.coastlines()
        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',edgecolor='k',
                                                facecolor=cfeature.COLORS['land']+0.0625)
        provinces_50m = cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces_lines','50m',facecolor='none')
        m.gridlines(draw_labels=True,auto_update=True)
        m.add_feature(land_50m,alpha=0.1)
        m.add_feature(cfeature.LAKES, facecolor=[0.69375   , 0.81484375, 0.9828125 ],alpha=0.3)
        m.add_feature(cfeature.RIVERS,alpha=0.2)
        m.add_feature(provinces_50m)
        m.add_feature(cfeature.BORDERS)
        
        m.set_extent([lower_left[0], upper_right[0],lower_left[1],upper_right[1]])
        
        m.orig_xlim = m.get_xlim()
        m.orig_ylim = m.get_ylim()
        m.llcrnrlat,m.urcrnrlat = m.get_ylim() 
        m.llcrnrlon,m.urcrnrlon = m.get_xlim()
        m.large = False
        m.figure.canvas.draw()
        m.ax = m.axes
        
        def converter_latlon(x,y):
            merc = ccrs.PlateCarree()
            return merc.transform_point(x,y,src_crs=m.proj)
        def inverter_lonlat(lon,lat):
            merc = ccrs.PlateCarree()
            
            if hasattr(lon,'__len__'):
                x_tmp,y_tmp = [],[]
                tp = [m.proj.transform_point(lon[ilon],lat[ilon],src_crs=merc) for ilon,llo in enumerate(lon)]
                nul = [(x_tmp.append(t[0]),y_tmp.append(t[1])) for t in tp]
                return x_tmp, y_tmp
            else:
                 return m.proj.transform_point(lon,lat,src_crs=merc)
        m.merc = ccrs.PlateCarree()
        m.convert_latlon = converter_latlon
        m.invert_lonlat = inverter_lonlat
        
    else:
        m.drawcoastlines()
        #m.fillcontinents(color='#AAAAAA')
        m.drawstates()
        m.drawcountries()
        m.large = True
        round_to_5 = lambda x:(int(x/5)+1)*5 
        round_to_2 = lambda x:(int(x/2)+1)*2
        if (upper_right[0]-lower_left[0])<20.0:
            mer = np.arange(round_to_2(lower_left[0]-dp),round_to_2(upper_right[0]+dp)+2,2)
            difx = 0.2
            m.large = False
        else:
            mer = np.arange(round_to_5(lower_left[0]-dp),round_to_5(upper_right[0]+dp)+5,5)
            difx = 1.0
        if (upper_right[1]-lower_left[1])<20.0:
            par = np.arange(round_to_2(lower_left[1]-dp),round_to_2(upper_right[1]+dp)+2,2)
            dify = 0.2
            m.large = False
        else:
            par = np.arange(round_to_5(lower_left[1]-dp),round_to_5(upper_right[1]+dp)+5,5)
            dify = 1.0
        if ax:
            ax.set_xlim(lower_left[0],upper_right[0])
            ax.set_ylim(lower_left[1],upper_right[1])
        m.artists.append(m.drawmeridians(mer,labels=[0,0,0,1]))
        m.artists.append(m.drawparallels(par,labels=[1,0,0,0]))
        m.par = par
        m.mer = mer
        m.orig_xlim = m.ax.get_xlim()
        m.orig_ylim = m.ax.get_ylim()
        def dummy_convert(x,y):
             return(x,y)
        m.convert_latlon = dummy_convert
        m.invert_lonlat = m
        m.merc = None

        # move the meridian labels to a proper position
        for aa in m.artists[0].keys():
            try:
                m.artists[0][aa][1][0].set_position((m.artists[0][aa][1][0].get_position()[0],lower_left[1]-difx))
            except:
                pass
        # move the parallels labels to a proper position
        for aa in m.artists[1].keys():
            try:
                m.artists[1][aa][1][0].set_position((lower_left[0]-dify,m.artists[1][aa][1][0].get_position()[1]))
            except:
                pass
    #import pdb; pdb.set_trace()
    if ax:
        try:
            ax.figure.show()
        except:
            try:
                ax.figure.canvas.draw()
            except:
                pass
    return m

def update_pars_mers(m,meridians,parallels,lower_left=None):
    'Simple program to remove old meridians and parallels and plot new ones'
    r = 0
    for a in m.artists:
        try:
            a.remove()
        except AttributeError:
            for b in a.values():
                try:
                    b.remove()
                except:
                    for c in b:
                        try:
                            c.remove()
                        except:
                            for d in c:
                                try:
                                    d.remove()
                                except:
                                    r = r+1

    m.artists = []
    m.artists.append(m.drawmeridians(meridians,labels=[0,0,0,1]))
    m.artists.append(m.drawparallels(parallels,labels=[1,0,0,0]))
    # move the meridian labels to a proper position
    difx = (meridians[1]-meridians[0])/5.0
    dify = (parallels[1]-parallels[0])/5.0
    if lower_left:
        x0,y0 = lower_left[0],lower_left[1]
    else:
        x0,y0 = parallels[0],meridians[0]
    for aa in m.artists[0].keys():
        m.artists[0][aa][1][0].set_position((m.artists[0][aa][1][0].get_position()[0],y0-dify))
    # move the parallels labels to a proper position
    for aa in m.artists[1].keys():
        m.artists[1][aa][1][0].set_position((x0-difx,m.artists[1][aa][1][0].get_position()[1]))

def pll(string):
    """
    pll for parse_lat_lon
    function that parses a string and converts it to lat lon values
    one space indicates seperation between degree, minutes, or minutes and seconds
    returns decimal degrees
    """
    if type(string) is float:
        return string
    if type(string) is int:
        return float(string)
    if not type(string) is str:
        try:
            return float(string)
        except TypeError:
            print('Error with pll input, trying to return first value')
            return float(string[0])
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
        deg_m = deg_m/60.0
        if str_ls[i]:
            deg_m = deg_m + float(str_ls[i])/60.0
    return deg+(deg_m*sign)

def plot_map_labels(m,filename,marker=None,skip_lines=0,color=None,textcolor='k',alpha=0.3):
    """
    program to plot the map labels on the basemap plot defined by m
    if marker is set, then it will be the default for all points in file
    """
    labels = mi.load_map_labels(filename,skip_lines=skip_lines) 
    
    labels_points = []
    for j,l in enumerate(labels):
        try:
            x,y = m.invert_lonlat(l['lon'],l['lat'])
            xtxt,ytxt = m.invert_lonlat(l['lon']+0.05,l['lat'])
        except:
            x,y = l['lon'],l['lat']
            xtxt,ytxt = l['lon']+0.05,l['lat']
        if marker:
            ma = marker
        else:
            ma = l['marker'] 
        #print('filename: {fa},label:{la},color:{co}'.format(fa=filename,la=l['label'],co=color))
        #if color:
        #    import pdb; pdb.set_trace()
        if color:
            co = color
            tco = textcolor
        else:
            co = l.get('color',color)
            tco = l.get('color',textcolor)
        if not is_color_like(co): 
            print('Color selection error in labels file: {fa}, for label: {la}'.format(fa=filename,la=l['label']))
            co = 'k'
            tco = co
        point = m.plot(x,y,color=co,marker=ma,alpha=alpha)
        m.annotate(l['label'],(xtxt,ytxt),color=tco,alpha=alpha)
        point[0].set_pickradius(10)
        labels_points.append(point[0])
            
        #import pdb; pdb.set_trace()
    return labels_points

def load_map_labels(filename,skip_lines=0):
    """
    Program to load map labels from a text file, csv
    with format: Label, lon, lat, style
    returns list of dictionary with each key as the label
    """
    out = []
    with open(filename,'r') as f:
        for i in range(skip_lines):
            next(f)
        for line in f:
            try:
                sp = line.split(',')
            except Exception as e:
                continue
            if len(sp)<2: continue
            if sp[0].startswith('#'):
                continue
            try:
                marker = sp[3].rstrip('\n').strip()
                if len(marker)>1: 
                    color = marker[1:]
                    marker = marker[0]
                else:
                    color = 'k'
                out.append({'label':sp[0],'lon':mi.pll(sp[1]),'lat':mi.pll(sp[2]),'marker':marker,'color':color})
            except IndexError:
                continue
    return out
    
def load_WMS_file(filename,skip_lines=0):
    """
    Program to load a file with multiple WMS servers
    each line is comma delimited
    first part is the name of the server, second is the website
    returns list of dictionary with each key as a seperate server
    # denotes a comment line
    
    """
    def str2bool(v):
        return v.strip().lower() in ("yes", "true", "t", "1", "on")
    out = []
    with open(filename,'r') as f:
        for i in range(skip_lines):
            next(f)
        for line in f:
            sp = line.split(',')
            if sp[0].startswith('#'):
                continue
            if len(sp)>2:
                out.append({'name':sp[0].strip(),'website':sp[1].strip(),'notime':str2bool(sp[2].rstrip('\n'))})
            else:
                out.append({'name':sp[0].strip(),'website':sp[1].rstrip('\n'),'notime':False})
    return out
    
def load_sat_from_net():
    """
    Program to load the satllite track prediction from the internet
    Checks at the avdc website
    """
    from datetime import datetime,timedelta
    from urllib2 import urlopen
    from pykml import parser
    today = datetime.now().strftime('%Y%m%d')
    site = 'http://avdc.gsfc.nasa.gov/download_2.php?site=98675770&id=25&go=download&path=%2FSubsatellite%2Fkml&file=A-Train_subsatellite_prediction_'+today+'T000000Z.kml'
    print('Satellite tracks url: %s' %site)
    try:
        response = urlopen(site)
        print('Getting the kml prediction file from avdc.gsfc.nasa.gov')
        r = response.read()
        kml = parser.fromstring(r)
    except:
        print('Problem with day, trying previous day...')
        try:
            yesterday = (datetime.now()-timedelta(days=1)).strftime('%Y%m%d')
            site = 'http://avdc.gsfc.nasa.gov/download_2.php?site=98675770&id=25&go=download&path=%2FSubsatellite%2Fkml&file=A-Train_subsatellite_prediction_'+yesterday+'T000000Z.kml'
            response = urlopen(site)
            print('Getting the kml prediction file from avdc.gsfc.nasa.gov')
            r = response.read()
            kml = parser.fromstring(r)
        except:
            import tkMessageBox
            tkMessageBox.showerror('No sat','There was an error communicating with avdc.gsfc.nasa.gov')
            return None
    print('Kml file read...')
    return kml

def load_sat_from_file(filename):
    """
    Program to load the satellite track prediction from a saved file
    """
    from pykml import parser
    f = open(filename,'r')
    r = f.read()
    kml = parser.fromstring(r)
    return kml

def get_sat_tracks(datestr,kml):
    """
    Program that goes and fetches the satellite tracks for the day
    For the day defined with datestr
    kml is the parsed kml structure with pykml
    """
    try:
        from map_interactive import pll
    except ModuleNotFoundError:
        from .map_interactive import pll
    sat = dict()
    # properly format datestr
    day = datestr.replace('-','')
    for i in range(4):
        name = str(kml.Document.Document[i].name).split(':')[1].lstrip(' ')
        for j in range(1,kml.Document.Document[i].countchildren()-1):
            if str(kml.Document.Document[i].Placemark[j].name).find(day) > 0:
                pos_str = str(kml.Document.Document[i].Placemark[j].LineString.coordinates)
                post_fl = pos_str.split(' ')
                lon,lat = [],[]
                for s in post_fl:
                    try:
                        x,y = s.split(',')
                        lon.append(pll(x))
                        lat.append(pll(y))
                    except:
                        pass
        try:
            sat[name] = (lon,lat)
        except UnboundLocalError:
            print('Skipping %s; no points downloaded' %name)
    return sat

def plot_sat_tracks(m,sat,label_every=5,max_num=60): 
    """
    Program that goes through and plots the satellite tracks
    """
    try:
        import map_utils as mu
    except ModuleNotFoundError:
        from . import map_utils as mu
    sat_obj = []
    sat_lines = {}
    sat_text = {}
    sat_keys = []
    for k in sat.keys():
        if type(sat[k]) is dict:
            lon = sat[k]['lon']
            lat = sat[k]['lat']
        else:
            (lon,lat) = sat[k]
        #x,y = m.invert_lonlat(lon,lat) # x,y = m(lon,lat)
        x,y = lon,lat
        tmp_l = m.plot(x,y,marker='+',markersize=1,label=k,linestyle='-',linewidth=0.2,transform=m.merc)
        sat_obj.append(tmp_l)
        sat_lines[k] = tmp_l
        sat_keys.append(k)
        co = sat_obj[-1][-1].get_color()
        #sat_obj.append(mu.mplot_spec(m,lon,lat,'-',linewidth=0.2))

        if type(sat[k]) is dict:
            try:
                latrange = [m.llcrnrlat,m.urcrnrlat]
                lonrange = [m.llcrnrlon,m.urcrnrlon]
            except:
                latrange = m.get_ylim()
                lonrange = m.get_xlim()
            #llcrnr = m.invert_lonlat(lonrange[0],latrange[0])
            #urcrnr = m.invert_lonlat(lonrange[1],latrange[1])
            sat_text[k] = []
            #print('length of sat labels {}, every {} with max of {}'.format(len(sat[k]['d']),label_every,max_num))
            if m.use_cartopy:
                [x0, x1, y0, y1] = m.ax.get_extent()
                transformed_points = m.proj.transform_points(m.merc,np.array(lon[::label_every]), np.array(lat[::label_every]))                   
                j=0
                for i,d in enumerate(sat[k]['d']):
                    if not i%label_every:
                        #transformed_point = m.proj.transform_point(lon[i], lat[i], m.merc)
                        transformed_point = transformed_points[j,:]
                        j = j+1
                        if (transformed_point[0]>x0) & (transformed_point[0]<x1) & (transformed_point[1]>y0) & (transformed_point[1]<y1):
                            sat_obj.append(m.ax.text(x[i],y[i],'%02i:%02i' % (d.tuple()[3],d.tuple()[4]),color=co,transform=m.merc))
                            sat_text[k].append(sat_obj[-1])
            else:
                if (len(sat[k]['d'])/label_every) < (max_num/2):
                    label_every = int(len(sat[k]['d'])/max_num/2)
                for i,d in enumerate(sat[k]['d']):
                    if (not i%label_every) & (i<max_num):
                        if ((lat[i]>=latrange[0])&(lat[i]<=latrange[1])&(lon[i]>=lonrange[0])&(lon[i]<=lonrange[1])):#((lat[i]>=llcrnr[1])&(lat[i]<=urcrnr[1])&(lon[i]>=llcrnr[0])&(lon[i]<=urcrnr[0])):
                            sat_obj.append(m.ax.text(x[i],y[i],'%02i:%02i' % (d.tuple()[3],d.tuple()[4]),color=co,transform=m.merc))
                            sat_text[k].append(sat_obj[-1])
    if len(sat.keys())>4:
        ncol = 2
    else:
        ncol = 1
    sat_legend = m.ax.legend(loc='lower right',bbox_to_anchor=(1.05,1.04),ncol=ncol)
    #try:    
    #    sat_legend.draggable()
    #except:
    #    sat_legend.set_draggable(True)
    #handle the toggle visibility
    sat_leg_lines = sat_legend.get_lines()
    graphs = {}
    texts = {}
    for i,k in enumerate(sat_keys):
        graphs[sat_leg_lines[i]] = sat_lines[k]
        texts[sat_leg_lines[i]] = sat_text[k]
        sat_leg_lines[i].set_picker(5)
        sat_leg_lines[i].set_lw(0.8)
        #sat_leg_lines[i].set_pickradius(5)
    graphs['lastserial'] = 0
    
    def on_pick_legend(event):
        legend0 = event.artist
#        global lastserial
        if event.guiEvent.serial == graphs['lastserial']: return
        graphs['lastserial'] = event.guiEvent.serial
        isVisible = True if graphs[legend0][0].get_alpha() is None else (True if graphs[legend0][0].get_alpha() else False)
        [g.set_visible(not isVisible) for g in graphs[legend0]]
        [g.set_alpha(1.0 if not isVisible else 0) for g in graphs[legend0]]

        for te in texts[legend0]: 
            te.set_visible(not isVisible)
            te.set_alpha(1.0 if not isVisible else 0)
        legend0.set_alpha(1.0 if not isVisible else 0.2) #set_visible(not isVisible)
        #print(legend0,'toggled, for line:',graphs[legend0],isVisible,legend0.get_alpha)
        #import pdb; pdb.set_trace()
        m.figure.canvas.draw()
    
    leg_onpick = m.figure.canvas.mpl_connect('pick_event',on_pick_legend)
    m.ax.add_artist(sat_legend)
    sat_obj.append(sat_legend)
    #sat_obj.append(leg_onpick)
    return sat_obj

def get_sat_tracks_from_tle(datestr,fraction_minute_interval=2,sat_filename='sat.tle'):
    """
    Program to build the satellite tracks from the two line element file
    """
    import ephem
    import numpy as np
    try:
        from map_interactive import get_tle_from_file
    except ModuleNotFoundError:
        from .map_interactive import get_tle_from_file
    import os
    from datetime import datetime, timedelta
    import tkinter.messagebox as tkMessageBox
    try:
        fname = os.path.join('.',sat_filename)
        sat = get_tle_from_file(fname)
    except:
        try:
            try:
                from gui import gui_file_select_fx
            except ModuleNotFoundError:
                from .gui import gui_file_select_fx
            fname = gui_file_select_fx(ext='*.tle',ftype=[('All files','*.*'),('Two Line element','*.tle')])
            sat = get_tle_from_file(fname)
        except:
            tkMessageBox.showerror('No sat','There was an error reading the sat.tle file')
            return None
    dt = datetime.now()-datetime.fromtimestamp(os.path.getmtime(fname))
    if dt>timedelta(days=14):
        tkMessageBox.showerror('Old TLEs','The file {} has been modified more than 14 days ago. Satellite tracks may be innacurate. Please update.'.format(fname))
    else:
        print('... Imported the {} file, now plotting'.format(os.path.abspath(fname)))
    for k in sat.keys():
        try:
            sat[k]['ephem'] = ephem.readtle(k,sat[k]['tle1'],sat[k]['tle2'])
        except ValueError:
            tkMessageBox.showerror('bad TLE line','The Satellite {} line do not conform to TLE standard. Please update.'.format(k))
        sat[k]['d'] = [ephem.Date(datestr+' 00:00')]
        sat[k]['ephem'].compute(sat[k]['d'][0])
        sat[k]['lat'] = [np.rad2deg(sat[k]['ephem'].sublat)]
        sat[k]['lon'] = [np.rad2deg(sat[k]['ephem'].sublong)]
        for t in range(24*60*fraction_minute_interval):
            d = ephem.Date(sat[k]['d'][t]+ephem.minute/float(fraction_minute_interval))
            sat[k]['d'].append(d)
            sat[k]['ephem'].compute(d)
            sat[k]['lat'].append(np.rad2deg(sat[k]['ephem'].sublat))
            sat[k]['lon'].append(np.rad2deg(sat[k]['ephem'].sublong))
    return sat

def get_tle_from_file(filename):
    'Program to load the tle from a file, skips lines with #'
    sat = {}
    name = ''
    first = ''
    second = ''
    for line in open(filename,'r'):
        if line.find('#')>=0:
            continue
        if not name:
            name = line.strip()
            if name.startswith('0'):
                name = name.strip('0').strip()
            continue
        if not first:
            first = line.strip()
            continue
        if not second:
            second = line.strip()
            sat[name] = {'tle1':first,'tle2':second}
            name = ''
            first = ''
            second = ''
    return sat
    
def get_flt_modules():
    'Program to create a list of *.flt files found returns dict of file path, file name, and linked png (is it exists)'
    import os
    fnames_all = os.listdir(os.path.join('.','flt_module'))
    fnames_all.sort()
    fnames = [g for g in fnames_all if g.endswith('flt')] #check correct file ending
    fnames.sort()
    dict = {}
    for f in fnames:
        png = os.path.abspath(os.path.join('.','flt_module',f.split('.')[0]+'.png'))
        if not os.path.isfile(png):
            png = os.path.abspath(os.path.join('.','flt_module',f.split('.')[0]+'.PNG'))
            if not os.path.isfile(png):
                png = None
        dict[f] = {'path':os.path.abspath(os.path.join('.','flt_module',f)),
                   'png':png}
    return dict


def extract_elevation_from_geotiff(latitudes, longitudes,geotiff_path='elevation_10KMmd_GMTEDmd.tif'):
    """
    Load the geotiff with surface elevation. Subset and find the latitudes and longitudes and return the elevation.
    
    Using the DEM by default of 10km, combined and described by:
    
    Amatulli, G., Domisch, S., Tuanmu, M.-N., Parmentier, B., Ranipeta, A., Malczyk, J., and Jetz, W. (2018)
    A suite of global, cross-scale topographic variables for environmental and biodiversity modeling. Scientific Data volume 5, 
    Article number: 180040. DOI: doi:10.1038/sdata.2018.40.
    """
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(geotiff_path) as src:
        # Retrieve the geographic transform (geotransform) information
        transform = src.transform

        # Convert latitude and longitude coordinates to pixel indices
        col_indices, row_indices = ~transform * (longitudes, latitudes)

        # Round the pixel indices to the nearest whole number
        col_indices = col_indices.round().astype(int)
        row_indices = row_indices.round().astype(int)

        # Read the elevation values corresponding to the pixel indices
        elevations = src.read(1)
        #, window=((row_indices.min(), row_indices.max() + 1),
        #                                 (col_indices.min(), col_indices.max() + 1)))

    return elevations[row_indices, col_indices]
        
def get_elev(utc,lat,lon,dt=60,geotiff_path='elevation_10KMmd_GMTEDmd.tif'):
    'Function returning surface elevation and a new latitude, longitude, and time for overplotting at an interval of dt seconds (default 300 seconds)'
    from scipy import interpolate
    import os
    utc = np.array(utc)
    lat = np.array(lat)
    lon = np.array(lon)
    
    utcs = np.arange(utc[0]*3600,utc[-1]*3600,dt)
    fx_lat = interpolate.interp1d(utc*3600,lat,bounds_error=False)
    fx_lon = interpolate.interp1d(utc*3600,lon,bounds_error=False)
    lat_new = fx_lat(utcs)
    lon_new = fx_lon(utcs)
    fname = os.path.join('.',geotiff_path)
    if not os.path.isfile(fname):
        try:
            try:
                from gui import gui_file_select_fx
            except ModuleNotFoundError:
                from .gui import gui_file_select_fx
            fname = gui_file_select_fx(ext='*.tif',ftype=[('All files','*.*'),('GEOTIFF','*.tif')])
        except:
            tkMessageBox.showerror('No GEOTIFF DEM surface Elevation','There was an error reading the DEM surface elevation file: {}'.format(fname))
    
    elev = extract_elevation_from_geotiff(lat_new,lon_new,geotiff_path=fname)
    return elev,lat_new,lon_new,utcs/3600.0,fname
    
def parse_and_plot_kml(kml_content, ax,color='tab:pink'):
    'function to plot the kml content (either kml file or as part of kmz'
    import xml.etree.ElementTree as ET
    import cartopy.crs as ccrs
    # Parse the KML content
    root = ET.fromstring(kml_content)
        # Define the namespaces
    namespaces = {
        'kml': 'http://www.opengis.net/kml/2.2',
        'google_kml': 'http://earth.google.com/kml/2.0',
        'gx': 'http://www.google.com/kml/ext/2.2',
        '': 'http://earth.google.com/kml/2.0'
    }
    # Find the Placemark elements in KML
    placemarks = root.findall('.//kml:Placemark', namespaces) + root.findall('.//google_kml:Placemark', namespaces) + root.findall('.//gx:Placemark', namespaces)
    # Iterate over each Placemark and plot its geometry
    print('... Found {} placemarks in the kml/kmz file'.format(len(placemarks)))
    plcmrks = []
    for placemark in placemarks:
        geom = placemark.findtext('.//coordinates', namespaces=namespaces)
        if geom is None:
            geom = placemark.findtext('.//kml:Point/kml:coordinates', namespaces=namespaces)
        if geom is None:
            geom = placemark.find('.//google_kml:Polygon/google_kml:outerBoundaryIs/google_kml:LinearRing/google_kml:coordinates', namespaces)
        if geom is None:
            geom = placemark.find('.//gx:Polygon/gx:outerBoundaryIs/gx:LinearRing/gx:coordinates', namespaces)
        
        #import ipdb; ipdb.set_trace()
        if geom is not None:
            coords = [tuple(map(float, coord.split(','))) for coord in geom.strip().split()]
            longitudes, latitudes = zip(*coords)
            # Plot the polygon
            pp_tmp = ax.plot(longitudes, latitudes, color=color, linewidth=1,transform=ccrs.Geodetic())
            plcmrks.append(pp_tmp)
        
    return plcmrks

def plot_kml(kml_file, ax,color='tab:pink'):
    'function to plot the insides of the kml file'
    import zipfile
    import os
    print('...opening KML/KMZ: {}'.format(os.path.abspath(kml_file)))
    
    # Extract the contents of KMZ file if provided
    if kml_file.endswith('.kmz'):
        with zipfile.ZipFile(kml_file, 'r') as zfile:
            # Find all KML files within the KMZ archive
            kml_files = [f for f in zfile.namelist() if f.lower().endswith('.kml')]
            if not kml_files:
                raise ValueError('No KML file found in the KMZ archive.')
            # Loop through each KML file
            for kml_file in kml_files:
                with zfile.open(kml_file) as kfile:
                    kml_content = kfile.read()
                    plots = parse_and_plot_kml(kml_content, ax,color=color)
    else:
        with open(kml_file, 'r') as kfile:
            kml_content = kfile.read()
        plots = parse_and_plot_kml(kml_content, ax,color=color)
    print('... Plotted {} lines from KML'.format(len(plots)))
    return plots

def convert_ccrs_to_epsg(ccrs_string):
    'Function to convert the cartoipy projection values to an epsg numrical value'
    projs_dict = {'PlateCarree':4326,
                   'NorthPolarStereo':3995,
                   'AlbersEqualArea':9822,
                   'AzimuthalEquidistant':2163,
                   'LambertCylindrical':9834,
                   'Mercator':3857,
                   'Miller':54003,
                   'Mollweide':54009,
                   'Orthographic':9840,
                   'Robinson':54030,
                   'Stereographic':3995,
                   'SouthPolarStereo':3031,
                   'Geostationary':4121}
    return projs_dict.get(ccrs_string,4326)
    
def alt2pres(altitude):
    '''
    Determine site pressure from altitude.

    Parameters
    ----------
    altitude : numeric
        Altitude above sea level. [m]

    Returns
    -------
    pressure : numeric
        Atmospheric pressure. [Pa]

    Notes
    ------
    The following assumptions are made

    ============================   ================
    Parameter                      Value
    ============================   ================
    Base pressure                  101325 Pa
    Temperature at zero altitude   288.15 K
    Gravitational acceleration     9.80665 m/s^2
    Lapse rate                     -6.5E-3 K/m
    Gas constant for air           287.053 J/(kg K)
    Relative Humidity              0%
    ============================   ================

    References
    -----------
    .. [1] "A Quick Derivation relating altitude to air pressure" from
       Portland State Aerospace Society, Version 1.03, 12/22/2004.
    '''

    press = 100 * ((44331.514 - altitude) / 11880.516) ** (1 / 0.1902632)

    return press

def pres2alt(pressure):
    '''
    Determine altitude from site pressure.

    Parameters
    ----------
    pressure : numeric
        Atmospheric pressure. [Pa]

    Returns
    -------
    altitude : numeric
        Altitude above sea level. [m]

    Notes
    ------
    The following assumptions are made

    ============================   ================
    Parameter                      Value
    ============================   ================
    Base pressure                  101325 Pa
    Temperature at zero altitude   288.15 K
    Gravitational acceleration     9.80665 m/s^2
    Lapse rate                     -6.5E-3 K/m
    Gas constant for air           287.053 J/(kg K)
    Relative Humidity              0%
    ============================   ================

    References
    -----------
    .. [1] "A Quick Derivation relating altitude to air pressure" from
       Portland State Aerospace Society, Version 1.03, 12/22/2004.
    '''

    alt = 44331.5 - 4946.62 * pressure ** (0.190263)

    return alt