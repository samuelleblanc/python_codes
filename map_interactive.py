# map_interactive software
# to use in combination with moving_lines
# Copyright 2015, Samuel LeBlanc
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys
import re
import copy

import map_interactive as mi
from map_utils import spherical_dist,equi,shoot,bearing

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
        Modified: Samuel LeBlanc, 209-06-03, Santa Cruz, CA
                 - fix linepicker issue with newer matplotlib. 
                 
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
            self.lons,self.lats = self.m(self.xs,self.ys,inverse=True)
            self.large = self.m.large
        self.connect()
        self.line.axes.format_coord = self.format_position_simple
        self.press = None
        self.contains = False
        self.labelsoff = False
        self.circlesoff = False
        self.moving = False
        self.lbl = None
        self.verbose = verbose
        self.blit = blit
        self.get_bg()
        if not tb:
         #   import matplotlib.pyplot as plt
         #   self.tb = plt.get_current_fig_manager().toolbar
            print 'No tb set, will not work'
        else:
            self.tb = tb

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

    def disconnect(self):
        'Function to disconnect all events (except keypress)'
        self.line.figure.canvas.mpl_disconnect(self.cid_onpress)
        self.line.figure.canvas.mpl_disconnect(self.cid_onrelease)
        self.line.figure.canvas.mpl_disconnect(self.cid_onmotion)
        self.line.figure.canvas.mpl_disconnect(self.cid_onkeyrelease)
        self.line.figure.canvas.mpl_disconnect(self.cid_onfigureenter)
        self.line.figure.canvas.mpl_disconnect(self.cid_onaxesenter)

    def onpress(self,event):
        'Function that enables either selecting a point, or creating a new point when clicked'
        #print 'click', event
        if event.inaxes!=self.line.axes: return
        if self.tb.mode!='': return
        if self.moving: return
        self.set_alt0 = False
        self.contains, attrd = self.line.contains(event)
        if self.contains:
            if self.verbose:
                print 'click is near point:',self.contains,attrd
            self.contains_index = attrd['ind']
            if len(self.contains_index)>1:
                self.contains_index = int(self.contains_index[-1])
            else:
                self.contains_index = int(self.contains_index)
            if self.verbose:
                print 'index:%i'%self.contains_index
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
                    lo,la = self.m(self.xs[self.contains_index],self.ys[self.contains_index],inverse=True)
                    self.lons.append(lo)
                    self.lats.append(la)
                self.set_alt0 = True
                self.contains = False
            ilola = self.contains_index
        else:
            self.xy = self.xs[-1],self.ys[-1]
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            if self.m:
                lo,la = self.m(event.xdata,event.ydata,inverse=True)
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
            print 'release'#,event
        if self.moving: return

        if (self.tb.mode == 'zoom rect') | (self.tb.mode == 'pan/zoom') | (self.tb.mode!=''):
            self.redraw_pars_mers()
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
        else:
            if self.ex:
                if self.set_alt0:
                    self.ex.appends(self.lats[-1],self.lons[-1],alt=self.ex.alt[0])
                else:
                    self.ex.appends(self.lats[-1],self.lons[-1])
                self.ex.calculate()
                self.ex.write_to_excel()
        for lrc in self.line.range_circles:
            self.m.ax.lines.remove(lrc)
        for arc in self.line.range_cir_anno:
            try:
                arc.remove()
            except:
                continue
        self.update_labels()
        self.line.figure.canvas.draw()
            

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
            self.lons[i],self.lats[i] = self.m(event.xdata,event.ydata,inverse=True)
        self.line.set_data(list(self.xs),list(self.ys))
        if self.contains:
            self.draw_canvas(extra_points=[self.highlight_linepoint,self.line.range_circles,self.line.range_cir_anno])
        else:
            self.draw_canvas(extra_points=[self.line.range_circles,self.line.range_cir_anno])       

    def onkeypress(self,event):
        'function to handle keyboard events'
        if self.verbose:
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
        'function to handle keyboard releases'
        #print 'released key',event.key
        if event.inaxes!=self.line.axes: return

    def onfigureenter(self,event):
        'event handler for updating the figure with excel data'
        if self.moving: return
        self.tb.set_message('Recalculating ...')
        if self.verbose:
            print 'entered figure'#, event
        if self.ex:
            #self.ex.switchsheet(self.iactive)
            self.ex.check_xl()
            self.lats = list(self.ex.lat)
            self.lons = list(self.ex.lon)
            if self.m:
                x,y = self.m(self.ex.lon,self.ex.lat)
                self.xs = list(x)
                self.ys = list(y)
                self.line.set_data(self.xs,self.ys)
                self.draw_canvas()
        self.update_labels()
        self.tb.set_message('Done Recalculating')
        self.line.axes.format_coord = self.format_position_simple
                
    def format_position_simple(self,x,y):
        'format the position indicator with only position'
        if self.m:
            return 'Lon=%.7f, Lat=%.7f'%(self.m(x, y, inverse = True))
        else:   
            return 'x=%2.5f, y=%2.5f' % (x,y)

    def format_position_distance(self,x,y):
        'format the position indicator with distance from previous point'
        if self.m:
            x0,y0 = self.xy
            lon0,lat0 = self.m(x0,y0,inverse=True)
            lon,lat = self.m(x,y,inverse=True)
            r = spherical_dist([lat0,lon0],[lat,lon])
            return 'Lon=%.7f, Lat=%.7f, d=%.2f km'%(lon,lat,r)
        else:
            x0,y0 = self.xy
            self.r = sqrt((x-x0)**2+(y-y0)**2)
            return 'x=%2.5f, y=%2.5f, d=%2.5f' % (x,y,self.r)
        
    def update_labels(self):
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
        for i in self.wp:    
            if not self.lbl:
                self.lbl = [self.line.axes.annotate(s+'%i'%i,
                                                    (self.xs[i-1],self.ys[i-1]))]
            else:
                try:
                    if not self.xs[i-1]:
                        continue
                    self.lbl.append(self.line.axes.
                                    annotate(s+'%i'%i,(self.xs[i-1],self.ys[i-1])))
                except IndexError:
                    pass                    
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
            ll, = equi(self.m,lon,lat,d,color=colors[i])
            line.append(ll)
            slon,slat,az = shoot(lon,lat,0.0,d)
            x,y = self.m(slon,slat)
            ano = self.line.axes.annotate('%i km' %(d),(x,y),color='silver')
            an.append(ano)
        if azi:
            slon,slat,az = shoot(lon,lat,azi,diam[-1])
            mlon,mlat,az = shoot(lon,lat,azi+180.0,diam[-1])
            lazi1, = self.m.plot([slon],[slat],'--*',color='grey',markeredgecolor='#BBBB00',markerfacecolor='#EEEE00',markersize=20)
            lazi2, = self.m.plot([mlon,lon,slon],[mlat,lat,slat],'--',color='grey')
            line.append(lazi1)
            line.append(lazi2)
        # plot angle points on the line
        for deg in [0,90,180,270]:
            dlo,dla,az = shoot(lon,lat,deg,diam[-1])
            elo,ela,az = shoot(lon,lat,deg,diam[-1]*0.85)
            ll, = self.m.plot([elo,dlo],[ela,dla],'-',color='grey')
            line.append(ll)
            for dd in [22.5,45.0,67.5]:
                dlo,dla,az = shoot(lon,lat,deg+dd,diam[-1])
                elo,ela,az = shoot(lon,lat,deg+dd,diam[-1]*0.93)
                ll, = self.m.plot([elo,dlo],[ela,dla],'-',color='grey')
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
	self.line_arr.append(line_new)

    def removeline(self,i):
        'Program to remove one line object from the LineBuilder class'
        self.line_arr[i].set_data([],[])
        self.line_arr[i].remove()

    def addfigure_under(self,img,ll_lat,ll_lon,ur_lat,ur_lon,outside=False,text=None,alpha=0.5,**kwargs):
    	'Program to add a figure under the basemap plot'
     
        left,bottom = self.m(ll_lon,ll_lat)
        right,top = self.m(ur_lon,ur_lat)
        try:
            lons = np.linspace(ll_lon,ur_lon,num=img.shape[1])
            lats = np.linspace(ll_lat,ur_lat,num=img.shape[0])
            print 'in first try addfigure under'
        except AttributeError:
            lons = np.linspace(ll_lon,ur_lon,num=img.size[1])
            lats = np.linspace(ll_lat,ur_lat,num=img.size[0])
        if (ur_lon-ll_lon) < (self.m.urcrnrlon-self.m.llcrnrlon):
            too_big_extent = True
        else:
            too_big_extent = False
        ix = np.where((lats>self.m.llcrnrlat)&(lats<self.m.urcrnrlat))[0]
        iy = np.where((lons>self.m.llcrnrlon)&(lons<self.m.urcrnrlon))[0]
        try:
            si = img.shape[0]
            ix = img.shape[0]-ix
        except AttributeError:
            si = img.size[0]
            ix = img.size[0]-ix
        ix[ix==si] = si-1
        if not outside:
            self.m.figure_under = None
            #self.m.imshow(img,zorder=0,extent=[left,right,top,bottom],**kwargs)
            try:
                if too_big_extent:
                    try:
                        x0,y0 = self.m(ll_lon,ll_lat)
                        x1,y1 = self.m(ur_lon,ur_lat)
                        self.m.figure_under = self.m.ax.imshow(np.flipud(img[ix,:,:][:,iy,:]),zorder=0,alpha=alpha,extent=[x0,x1,y0,y1],**kwargs)
                    except:
                        import pdb
                        pdb.set_trace()
                else:
                    self.m.figure_under = self.m.imshow(img[ix,:,:][:,iy,:],zorder=0,alpha=alpha,**kwargs)
            except Exception as ie:
                print 'problem occurred when placing under figure, trying anyway'
                print ie
                self.m.figure_under = self.m.imshow(img,zorder=0,alpha=alpha,**kwargs)
        else:
            u = self.m.imshow(img,clip_on=False,**kwargs)
        if text:
            try:
                self.m.figure_under_text = self.m.ax.text(0.0,0.0,text,transform=self.m.ax.transAxes)
            except:
                print 'Problem adding text on figure, continuning...'
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
        try:
            mi.update_pars_mers(self.m,mer,par,lower_left=(xlim[0],ylim[0]))
        except:
            print '... Problem updating the parallels and meridians'
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
            print 'New points at lon: %f, lat: %f' %(newlon,newlat)
        if self.m:
            x,y = self.m(newlon,newlat)
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
            x,y = self.m(newlon,newlat)
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
        from gui import ask
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
        choice = ['feet','meters']
        choice2 = ['km','nm']
        # remove values predefined
        for p in predef:
            try:
                names.remove(p)
            except:
                pass
                
        # ask the user supply the variables
        vals = ask(names,choice=choice,choice_title='Altitude:',choice2=choice2,choice2_title='Distance:',title='Enter values for {}'.format(name))
        v = vals.names_val
        if vals.choice_val == 'feet': 
            use_feet = True 
        else: 
            use_feet = False
        
        if vals.choice2_val == 'km': 
            use_km = True 
        else: 
            use_km = False
            
        # make the values variables
        for i,n in enumerate(names):
            try:
                exec('{}={}'.format(n,vals.names_val[i]))
            except:
                print 'problem for {}={}'.format(n,vals.names_val[i])
        for l in s[1:-1]:
            if not (l.startswith('#') or l.startswith('%')):
                try:
                    self.newpoint(*eval(l),last=False,feet=use_feet,km=use_km)
                except:
                    print 'problem with {}'.format(l)
        if not s[-1].startswith('#') or not s[-1].startswith('%'):
                try:
                    self.newpoint(*eval(s[-1]),feet=use_feet,km=use_km)
                except:
                    print 'problem with last {}'.format(s[-1])
        f.close()

    def get_bg(self,redraw=False):
        'program to store the canvas background. Used for blit technique'
        if redraw:
            self.line.figure.canvas.draw()
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
                print 'exception occurred: %s' %ie
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
        
def build_basemap(lower_left=[-20,-30],upper_right=[20,10],ax=None,proj='cyl',profile=None,larger=True):
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
    from map_interactive import pll
    if profile:
        upper_right = [pll(profile['Lon_range'][1]),pll(profile['Lat_range'][1])]
        lower_left = [pll(profile['Lon_range'][0]),pll(profile['Lat_range'][0])]

    if larger:
        dp = 30
    else:
        dp = 0
    try:
        import cPickle as pickle
        m = pickle.load(open('map_{}.pkl'.format(profile['Campaign']),'rb'))
        #print 'doing it pickle style'
        m.ax = ax
    except:
        if larger:
            dp = 30
        else:
            dp = 0
        m = Basemap(projection=proj,lon_0=(upper_right[0]+lower_left[0])/2.0,lat_0=(upper_right[1]+lower_left[1])/2.0,
                llcrnrlon=lower_left[0]-dp, llcrnrlat=lower_left[1]-dp,
                urcrnrlon=upper_right[0]+dp, urcrnrlat=upper_right[1]+dp,resolution='i',ax=ax)
    m.artists = []
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
            print 'Error with pll input, trying to return first value'
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

def plot_map_labels(m,filename,marker=None,skip_lines=0,color='k'):
    """
    program to plot the map labels on the basemap plot defined by m
    if marker is set, then it will be the default for all points in file
    """
    labels = mi.load_map_labels(filename,skip_lines=skip_lines) 
    for l in labels:
        try:
            x,y = m(l['lon'],l['lat'])
            xtxt,ytxt = m(l['lon']+0.05,l['lat'])
        except:
            x,y = l['lon'],l['lat']
            xtxt,ytxt = l['lon']+0.05,l['lat']
        if marker:
            ma = marker
        else:
            ma = l['marker'] 
        m.plot(x,y,color=color,marker=ma)
        m.ax.annotate(l['label'],(xtxt,ytxt))

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
            sp = line.split(',')
            if sp[0].startswith('#'):
                continue
            out.append({'label':sp[0],'lon':mi.pll(sp[1]),'lat':mi.pll(sp[2]),'marker':sp[3].rstrip('\n')})
    return out
    
def load_WMS_file(filename,skip_lines=0):
    """
    Program to load a file with multiple WMS servers
    each line is comma delimited
    first part is the name of the server, second is the website
    returns list of dictionary with each key as a seperate server
    # denotes a comment line
    
    """
    out = []
    with open(filename,'r') as f:
        for i in range(skip_lines):
            next(f)
        for line in f:
            sp = line.split(',')
            if sp[0].startswith('#'):
                continue
            out.append({'name':sp[0].strip(),'website':sp[1].rstrip('\n')})
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
    print 'Satellite tracks url: %s' %site
    try:
        response = urlopen(site)
        print 'Getting the kml prediction file from avdc.gsfc.nasa.gov'
        r = response.read()
        kml = parser.fromstring(r)
    except:
        print 'Problem with day, trying previous day...'
        try:
            yesterday = (datetime.now()-timedelta(days=1)).strftime('%Y%m%d')
            site = 'http://avdc.gsfc.nasa.gov/download_2.php?site=98675770&id=25&go=download&path=%2FSubsatellite%2Fkml&file=A-Train_subsatellite_prediction_'+yesterday+'T000000Z.kml'
            response = urlopen(site)
            print 'Getting the kml prediction file from avdc.gsfc.nasa.gov'
            r = response.read()
            kml = parser.fromstring(r)
        except:
            import tkMessageBox
            tkMessageBox.showerror('No sat','There was an error communicating with avdc.gsfc.nasa.gov')
            return None
    print 'Kml file read...'
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
    from map_interactive import pll
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
            print 'Skipping %s; no points downloaded' %name
    return sat

def plot_sat_tracks(m,sat,label_every=20): 
    """
    Program that goes through and plots the satellite tracks
    """
    import map_utils as mu
    sat_obj = []
    for k in sat.keys():
        if type(sat[k]) is dict:
            lon = sat[k]['lon']
            lat = sat[k]['lat']
        else:
            (lon,lat) = sat[k]
        x,y = m(lon,lat)
        sat_obj.append(m.plot(x,y,marker='.',label=k,linestyle='-',linewidth=0.2))
        co = sat_obj[-1][-1].get_color()
        #sat_obj.append(mu.mplot_spec(m,lon,lat,'-',linewidth=0.2))
        if type(sat[k]) is dict:
            latrange = [m.llcrnrlat,m.urcrnrlat]
            lonrange = [m.llcrnrlon,m.urcrnrlon]
            for i,d in enumerate(sat[k]['d']):
                if not i%label_every:
                    if ((lat[i]>=latrange[0])&(lat[i]<=latrange[1])&(lon[i]>=lonrange[0])&(lon[i]<=lonrange[1])):
                        sat_obj.append(m.ax.text(x[i],y[i],'%02i:%02i' % (d.tuple()[3],d.tuple()[4]),color=co))
    if len(sat.keys())>4:
        ncol = 2
    else:
        ncol = 1
    sat_legend = m.ax.legend(loc='lower right',bbox_to_anchor=(1.05,1.04),ncol=ncol)
    sat_legend.draggable()
    m.ax.add_artist(sat_legend)
    sat_obj.append(sat_legend)
    return sat_obj

def get_sat_tracks_from_tle(datestr):
    """
    Program to build the satellite tracks from the two line element file
    """
    import ephem
    import numpy as np
    from map_interactive import get_tle_from_file
    try:
        sat = get_tle_from_file('.\sat.tle')
    except:
        try:
            from gui import gui_file_select_fx
            fname = gui_file_select_fx(ext='*.tle',ftype=[('All files','*.*'),('Two Line element','*.tle')])
            sat = get_tle_from_file(fname)
        except:
            import tkMessageBox
            tkMessageBox.showerror('No sat','There was an error reading the sat.tle file')
            return None
    for k in sat.keys():
        sat[k]['ephem'] = ephem.readtle(k,sat[k]['tle1'],sat[k]['tle2'])
        sat[k]['d'] = [ephem.Date(datestr+' 00:00')]
        sat[k]['ephem'].compute(sat[k]['d'][0])
        sat[k]['lat'] = [np.rad2deg(sat[k]['ephem'].sublat)]
        sat[k]['lon'] = [np.rad2deg(sat[k]['ephem'].sublong)]
        for t in xrange(24*60*2):
            d = ephem.Date(sat[k]['d'][t]+ephem.minute/2.0)
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
    fnames = [g for g in fnames_all if g.endswith('flt')] #check correct file ending
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
    
        
    
