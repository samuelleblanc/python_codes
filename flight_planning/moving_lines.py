import matplotlib
import os
fp = os.path.dirname(os.path.abspath(__file__))
matplotlib.rc_file(fp+os.path.sep+'file.rc')
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap
import sys
import map_utils as mu
import excel_interface as ex
import map_interactive as mi
import gui

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
    g = gui.gui(lines)
    g.make_gui()
    return lines

if __name__ == "__main__":
    lines = Create_interaction()
