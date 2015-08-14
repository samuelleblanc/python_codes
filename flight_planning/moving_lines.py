import matplotlib
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#if __name__ == "__main__":
    #if run_from_ipython():
        #%matplotlib tk
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap
import sys
import map_utils as mu
import excel_interface as ex
import map_interactive as mi

def Create_interaction():
    fig,ax = plt.subplots()
    m = mi.build_basemap(ax=ax)
    ax.set_title('line segments')
    lat0,lon0 = mi.pll('22 58.783S'), mi.pll('14 38.717E')
    x0,y0 = m(lon0,lat0)
    line, = m.plot([x0],[y0],'ro-')
    text = ('Press s to stop interaction\\n'
            'Press i to restart interaction\\n')
    #plt.text(1.0,0.1,text)
    wb = ex.dict_position()
    lines = mi.LineBuilder(line,m=m,ex=wb)
    plt.show()
    return m

def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

if __name__ == "__main__":
    Create_interaction()
