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
    """
    def __init__(self,lon0='14 38.717E',lat0='22 58.783S',speed=150.0,UTC_start=7.0,UTC_conversion=+1.0,alt0=0.0):
        import numpy as np
        from xlwings import Range
        import map_interactive as mi
        from map_interactive import pll
        self.lon = np.array([pll(lon0)])
        self.lat = np.array([pll(lat0)])
        self.speed = np.array([speed])
        self.alt = np.array([alt0])
        self.UTC_conversion = UTC_conversion
        self.UTC = np.array(UTC_start)
        self.calculate()
        self.wb = self.Create_excel()
        self.write_to_excel()

    def calculate(self):
        """
        Program to fill in all the missing pieces in the dict_position class
        Involves converting from metric to aviation units
        Involves calculating distances
        Involves calculating time of flight local and utc
        Fills in the waypoint numbers

        Assumes that blank spaces are to be filled with new calculations
        """
        default_bank_angle = 15.0
        self.rate_of_turn = 1091.0*tan(default_bank_angle*np.pi/180)/self.speed[0]

    def write_to_excel(self):
        """
        writes out the dict_position class values to excel spreadsheet
        """
        Range('A2').value = np.array([self.WP,
                                      self.lat,
                                      self.lon,
                                      self.speed,
                                      self.delayt,
                                      self.alt,
                                      self.cumlegt,
                                      self.utc,
                                      self.local,
                                      self.legt,
                                      self.dist,
                                      self.cumdist,
                                      self.dist_nm,
                                      self.cumdist_nm,
                                      self.speed_kts,
                                      self.alt_kft
                                      ]).T

    def read_excel(self):
        """
        raw read and parse of the excel spreadsheet values
        """
        tmp = Range('A2').table.value
        self.tmp

    def check_updates_excel(self):
        """
        Check for any change in the excel file
        If there is change, empty out the corresponding calculated areas
        Priority is always given to metric
        """
        self.read_excel()
        
        print 'not yet'

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
            xl.ActiveWorkbook.Windows(1).SplitColumn = 0.4
            xl.ActiveWorkbook.Windows(1).SplitRow = 1.0
            xl.Range(address).Font.Bold = True
        top_line.autofit()
        Range('A2').value = np.arange(50).reshape((50,1))+1
        return wb
  

