class dict_position:
    """
    Purpose:
        Class that creates an easy storage for position coordinates.
        Encompasses array values of positions, altitude, time, dist.
        Along with program to update excel spreadsheet with info, read spreadsheet data,
        and update calculations for distances
    Inputs:
        lon0: [degree] initial longitude (optional, defaults to Namibia Walvis bay airport), can be string
        lat0: [degree] initial latitude (optional, defaults to Namibia Walvis bay airport), can be string
        speed: [m/s] speed of aircraft defaults to p3 value of 150 m/s (optional)
        UTC_start: [decimal hours] time of takeoff, defaults to 7.0 UTC (optional)
        UTC_conversion: [decimal hours] conversion (dt) used to change utc to local time (optional), local = utc + dt
        alt0: [m] initial altitude of the plane, airport altitude (optinal) 
    Outputs:
        dict_position class 
    Dependencies:
        xlwings
        Excel (win or mac)
        map_interactive
        map_utils
    Required files:
        none
    Example:
        ...
    Modification History:
        Written: Samuel LeBlanc, 2015-08-07, Santa Cruz, CA
        Modified: Samuel LeBlanc, 2015-08-11, Santa Cruz, CA
                 - update and bug fixes
    """
    def __init__(self,lon0='14 38.717E',lat0='22 58.783S',speed=150.0,UTC_start=7.0,UTC_conversion=+1.0,alt0=0.0):
        import numpy as np
        from xlwings import Range
        import map_interactive as mi
        from map_interactive import pll
        import map_utils as mu
        self.lon = np.array([pll(lon0)])
        self.lat = np.array([pll(lat0)])
        self.speed = np.array([speed])
        self.alt = np.array([alt0])
        self.UTC_conversion = UTC_conversion
        self.utc = np.array([UTC_start])
        self.UTC = self.utc
        self.legt = self.UTC*0.0
        self.dist = self.UTC*0.0
        self.cumdist = self.UTC*0.0
        self.cumlegt = self.legt
        self.delayt = self.legt
        self.bearing = self.lon*0.0
        self.endbearing = self.lon*0.0
        self.turn_deg = self.lon*0.0
        self.turn_time = self.lon*0.0
        self.speed_kts = self.speed*1.94384449246
        self.alt_kft = self.alt*3.28084/1000.0
        self.head = self.legt
        try:
            self.calculate()
        except:
            print 'calculate failed'
        self.wb = self.Create_excel()
        try:
            self.write_to_excel()
        except:
            print 'writing to excel failed'

    def calculate(self):
        """
        Program to fill in all the missing pieces in the dict_position class
        Involves converting from metric to aviation units
        Involves calculating distances
        Involves calculating time of flight local and utc
        Fills in the waypoint numbers

        Assumes that blank spaces/nan are to be filled with new calculations
        """
        import numpy as np
        import map_utils as mu
        default_bank_angle = 15.0
        self.rate_of_turn = 1091.0*np.tan(default_bank_angle*np.pi/180)/self.speed[0] # degree per second
        self.n = len(self.lon)
        self.WP = range(1,self.n+1)
        for i in xrange(self.n-1):
            self.dist[i+1] = mu.spherical_dist([self.lat[i],self.lon[i]],[self.lat[i+1],self.lon[i+1]])
            if np.isfinite(self.speed.astype(float)[i+1]):
                self.speed_kts[i+1] = self.speed[i+1]*1.94384449246
            elif np.isfinite(self.speed_kts.astype(float)[i+1]):
                self.speed[i+1] = seld.speed_kts[i+1]/1.94384449246
            else:
                self.speed[i+1] = self.speed[i]
                self.speed_kts[i+1] = self.speed[i+1]*1.94384449246
            if np.isfinite(self.alt.astype(float)[i+1]):
                self.alt_kft[i+1] = self.alt[i+1]*3.28084/1000.0
            elif np.isfinite(self.alt_kft.astype(float)[i+1]):
                self.alt[i+1] = self.alt_kft[i+1]*1000.0/3.28084
            else:
                self.alt[i+1] = self.alt[i]
                self.alt_kft[i+1] = self.alt[i+1]*3.28084/1000.0
            self.bearing[i] = mu.bearing([self.lat[i],self.lon[i]],[self.lat[i+1],self.lon[i+1]])
            self.endbearing[i] = (mu.bearing([self.lat[i+1],self.lon[i+1]],[self.lat[i],self.lon[i]])+180)%360.0
            try:
                self.bearing[i+1] = mu.bearing([self.lat[i+1],self.lon[i+1]],[self.lat[i+2],self.lon[i+2]])
            except:
                self.bearing[i+1] = self.endbearing[i]
            try:
                self.turn_deg[i+1] = abs(self.endbearing[i]-self.bearing[i+1])
            except:
                self.turn_deg[i+1] = 0.0
            self.turn_time[i+1] = (self.turn_deg[i+1]/self.rate_of_turn)/60.0
            if not np.isfinite(self.delayt.astype(float)[i+1]):
                self.delayt[i+1] = self.turn_time[i+1]
            #else:
            #    self.delayt[i+1] = self.delayt[i+1]+self.turn_time[i+1]
            self.legt[i+1] = (self.dist[i+1]/(self.speed[i+1]/1000.0))/3600.0 + self.delayt[i+1]/60.0
            self.utc[i+1] = self.utc[i]+self.legt[i+1]
            
        self.local = self.utc+self.UTC_conversion
        self.dist_nm = self.dist*0.53996
        self.cumdist = self.dist.cumsum()
        self.cumdist_nm = self.dist_nm.cumsum()
        self.cumlegt = np.nan_to_num(self.legt).cumsum()
        self.time2xl()

    def time2xl(self):
        """
        Convert the UTC fractional hours to hh:mm format for use in excel
        """
        self.cumlegt_xl = self.cumlegt/24.0
        self.utc_xl = self.utc/24.0
        self.local_xl = self.local/24.0
        self.legt_xl = self.legt/24.0

    def write_to_excel(self):
        """
        writes out the dict_position class values to excel spreadsheet
        """
        import numpy as np
        from xlwings import Range
        self.wb.set_current()
        Range('A2').value = np.array([self.WP,
                                      self.lat,
                                      self.lon,
                                      self.speed,
                                      self.delayt,
                                      self.alt,
                                      self.cumlegt_xl,
                                      self.utc_xl,
                                      self.local_xl,
                                      self.legt_xl,
                                      self.dist,
                                      self.cumdist,
                                      self.dist_nm,
                                      self.cumdist_nm,
                                      self.speed_kts,
                                      self.alt_kft
                                      ]).T
        Range('G2:J%i'% (self.n+1)).number_format = 'hh:mm'
        Range('E2:E%i'% (self.n+1)).number_format = '0'
        Range('B:B').autofit('c')
        Range('C:C').autofit('c')

    def check_xl(self):
        """
        wrapper for checking excel updates.
        Reruns check_updates_excel whenever a line is found to be deleted
        """
        while self.check_updates_excel():
            print 'line removed, cutting it out'

    def check_updates_excel(self):
        """
        Check for any change in the excel file
        If there is change, empty out the corresponding calculated areas
        Priority is always given to metric
        """
        from xlwings import Range
        import numpy as np
        self.wb.set_current()
        tmp = Range('A2:P%i'%(self.n+1)).value
        tmp0 = Range('A2:P2').vertical.value
        dim = np.shape(tmp)
        if len(dim)==1:
            tmp = [tmp]
        else:
            dim0 = np.shape(tmp0)
            if len(dim0)==1: dim0 = np.shape([tmp0])
            n0,_ = dim0
            n1,_ = dim
            if n0>n1: tmp = tmp0
            print 'vertical num: %i, range num: %i' %(n0,n1)
        num = 0
        num_del = 0
        for i,t in enumerate(tmp):
            if len(t)<16: continue
            wp,lat,lon,sp,dt,alt,clt,utc,loc,lt,d,cd,dnm,cdnm,spkt,altk = t[0:16]
            if wp > self.n:
                num = num+1
                self.appends(lat,lon,sp,dt,alt,clt,utc,loc,lt,d,cd,dnm,cdnm,spkt,altk)
            elif not wp: # check if empty
                num = num+1
                self.dels(i)
                self.move_xl(i)
                self.n = self.n-1
                return True
            else:
                changed = self.mods(i,lat,lon,sp,spkt,dt,alt,altk)
                if i == 0:
                    if self.utc[i] != utc*24.0:
                        self.utc[i] = utc*24.0
                        changed = True
                if changed: num = num+1
        if self.n>(i+1):
            print 'deleting points'
            for j in range(i+1,self.n-1):
                self.dels(j)
                num = num+1
        if num>0:
            print 'Updated %i lines from Excel, recalculating and printing' % num
            self.calculate()
            self.write_to_excel()
        self.num_changed = num
        return False

    def move_xl(self,i):
        """
        Program that moves up all excel rows by one line overriding the ith line
        """
        from xlwings import Range
        linesbelow = Range('A%i:P%i'%(i+3,self.n+1)).value
        Range('A%i:P%i'%(i+2,i+2)).value = linesbelow
        Range('A%i:P%i'%(self.n+1,self.n+1)).clear_contents()

    def dels(self,i):
        """
        program to remove the ith item in every object
        """
        import numpy as np
        self.lat = np.delete(self.lat,i)
        self.lon = np.delete(self.lon,i)
        self.speed = np.delete(self.speed,i)
        self.delayt = np.delete(self.delayt,i)
        self.alt = np.delete(self.alt,i)
        self.alt_kft = np.delete(self.alt_kft,i)
        self.speed_kts = np.delete(self.speed_kts,i)
        self.cumlegt = np.delete(self.cumlegt,i)
        self.utc = np.delete(self.utc,i)
        self.local = np.delete(self.local,i)
        self.legt = np.delete(self.legt,i)
        self.dist = np.delete(self.dist,i)
        self.cumdist = np.delete(self.cumdist,i)
        self.dist_nm = np.delete(self.dist_nm,i)
        self.cumdist_nm = np.delete(self.cumdist_nm,i)
        self.bearing = np.delete(self.bearing,i)
        self.endbearing = np.delete(self.endbearing,i)
        self.turn_deg = np.delete(self.turn_deg,i)
        self.turn_time = np.delete(self.turn_time,i)

    def appends(self,lat,lon,sp=None,dt=None,alt=None,
                clt=None,utc=None,loc=None,lt=None,d=None,cd=None,
                dnm=None,cdnm=None,spkt=None,altk=None,
                bear=0.0,endbear=0.0,turnd=0.0,turnt=0.0):
        """
        Program that appends to the current class with values supplied, or with defaults from the command line
        """
        import numpy as np
        self.lat = np.append(self.lat,lat)
        self.lon = np.append(self.lon,lon)
        self.speed = np.append(self.speed,sp)
        self.delayt = np.append(self.delayt,dt)
        self.alt = np.append(self.alt,alt)
        if not clt: clt = np.nan
        if not utc: utc = np.nan
        if not loc: loc = np.nan
        if not lt: lt = np.nan
        self.cumlegt = np.append(self.cumlegt,clt*24.0)
        self.utc = np.append(self.utc,utc*24.0)
        self.local = np.append(self.local,loc*24.0)
        self.legt = np.append(self.legt,lt*24.0)
        self.dist = np.append(self.dist,d)
        self.cumdist = np.append(self.cumdist,cd)
        self.dist_nm = np.append(self.dist_nm,dnm)
        self.cumdist_nm = np.append(self.cumdist_nm,cdnm)
        self.speed_kts = np.append(self.speed_kts,spkt)
        self.alt_kft = np.append(self.alt_kft,altk)
        self.bearing = np.append(self.bearing,bear)
        self.endbearing = np.append(self.endbearing,endbear)
        self.turn_deg = np.append(self.turn_deg,turnd)
        self.turn_time = np.append(self.turn_time,turnt)

    def mods(self,i,lat=None,lon=None,sp=None,spkt=None,dt=None,alt=None,altk=None):
        """
        Program to modify the contents of the current class if there is an update on the line, defned by i
        If anything is not input, then the default of NaN is used
        """
        import numpy as np
        changed = False
        self.toempty = {'speed':0,'delayt':0,'alt':0,'speed_kts':0,'alt_kft':0}
        if not lat: lat = np.nan
        if not lon: lon = np.nan
        if not sp: sp = np.nan
        if not spkt: spkt = np.nan
        if not dt: dt = np.nan
        if not alt: alt = np.nan
        if not altk: altk = np.nan
        if self.lat[i] != lat:
            self.lat[i] = lat
            changed = True
        if self.lon[i] != lon:
            self.lon[i] = lon
            changed = True
        if self.speed[i] != sp:
            self.speed[i] = sp
            self.toempty['speed_kts'] = 1
            changed = True
        if self.speed_kts[i] != spkt:
            if np.isfinite(spkt):
                self.toempty['speed'] = 1
                changed = True
        if self.delayt[i] != dt:
            self.delayt[i] = dt
            changed = True
        if self.alt[i] != alt:
            if np.isfinite(alt):
                self.alt[i] = alt
                self.toempty['alt_kft'] = 1
                changed = True
        if self.alt_kft[i] != altk:
            if np.isfinite(altk):
                self.alt_kft[i] = altk
                self.toempty['alt'] = 1
                changed = True
        for s in self.toempty:
            if self.toempty.get(s):
                v = getattr(self,s)
                v[i] = np.nan
                setattr(self,s,v)
        return changed
        
    def Create_excel(self):
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
            Modified: Samuel LeBlanc, 2015-08-07, Santa Cruz, CA
                    - put into the dic_position class, modified slightly
            
        """
        from xlwings import Workbook, Sheet, Range, Chart
        import numpy as np
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
         #   xl.ActiveWorkbook.Windows(1).SplitColumn = 0.4
            xl.ActiveWorkbook.Windows(1).SplitRow = 1.0
            xl.Range(address).Font.Bold = True
        top_line.autofit()
        Range('G2:J2').number_format = 'hh:mm'
        #Range('A2').value = np.arange(50).reshape((50,1))+1
        return wb
  

