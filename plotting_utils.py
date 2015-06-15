# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def __init__():
    """
       Collection of codes to run some typical plotting utilities
       
           - circles : plot circles with radius in data scale
           
        details are in the info of each module
    """
    pass

# <codecell>

def circles(x, y, s, c='b', ax=None, vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data scale (ie. in data unit)
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or
        RGBA sequence because that is indistinguishable from an array of
        values to be colormapped.  `c` can be a 2-D array in which the
        rows are RGB or RGBA, however.
    ax : Axes object, optional, default: None
        Parent axes of the plot. It uses gca() if not specified.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.  (Note if you pass a `norm` instance, your
        settings for `vmin` and `vmax` will be ignored.)

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Other parameters
    ----------------
    kwargs : `~matplotlib.collections.Collection` properties
        eg. alpha, edgecolors, facecolors, linewidths, linestyles, norm, cmap

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection
    import pylab as plt
    #import matplotlib.colors as colors

    if ax is None:
        ax = plt.gca()    

    if isinstance(c,basestring):
        color = c     # ie. use colors.colorConverter.to_rgba_array(c)
    else:
        color = None  # use cmap, norm after collection is created
    kwargs.update(color=color)

    if isinstance(x, (int, long, float)):
        patches = [Circle((x, y), s),]
    elif isinstance(s, (int, long, float)):
        patches = [Circle((x_,y_), s) for x_,y_ in zip(x,y)]
    else:
        patches = [Circle((x_,y_), s_) for x_,y_,s_ in zip(x,y,s)]
    collection = PatchCollection(patches, **kwargs)

    if color is None:
        collection.set_array(np.asarray(c))
        if vmin is not None or vmax is not None:
            collection.set_clim(vmin, vmax)

    ax.add_collection(collection)
    return collection

# <codecell>

def plot_color_maps(reverse=False):
    """
    Simple plotting function to run through and plot each color map
    Help for choosing which colormap to use
    """
    import pylab as plt
    from numpy import outer
    plt.rc('text', usetex=False)
    a=outer(plt.ones(10,),plt.arange(0,1,0.01))
    plt.figure(figsize=(5,15))
    plt.subplots_adjust(top=0.8,bottom=0.08,left=0.03,right=0.99)
    if reverse:
        maps=[m for m in plt.cm.datad]
        rr = 2
    else:
        maps=[m for m in plt.cm.datad if not m.endswith("_r")]
        rr = 1
    maps.sort()
    l=len(maps)+1
    title_dict = {'fontsize': 10,
                  'verticalalignment': 'center',
                  'horizontalalignment': 'left'}
    for i, m in enumerate(maps):
        plt.subplot(l,rr,i+1)
        plt.axis("off")
        plt.imshow(a,aspect='auto',cmap=plt.get_cmap(m),origin="lower")
        plt.text(1.01,0.5,m,fontdict=title_dict,transform=plt.gca().transAxes)

# <codecell>

def plot_vert_hist(fig,ax1,y,pos,ylim,color='grey',label=None,legend=False,onlyhist=True,loc=2):
    """
    function to plot a 'bean' like vertical histogram
    """
    import Sp_parameters as Sp
    import numpy as np
    from plotting_utils import data2figpoints
    (ymask,iy) = Sp.nanmasked(y)
    ax = fig.add_axes(data2figpoints(pos,0.4,fig=fig,ax1=ax1),frameon=False,ylim=ylim)
    ax.tick_params(axis='both', which='both', labelleft='off', labelright='off',bottom='off',top='off',
               labelbottom='off',labeltop='off',right='off',left='off')
    ax.hist(ymask,orientation='horizontal',normed=True,color=color,edgecolor='None',bins=30,alpha=0.5,label=label)
    if onlyhist:
        label_mean = None
        label_median = None
    else:
        label_mean = 'Mean'
        label_median = 'Median'
    ax.axhline(np.mean(ymask),color='red',linewidth=2,label=label_mean)
    ax.axhline(np.median(ymask),color='k',linewidth=2,linestyle='--',label=label_median)
    if legend:
        ax.legend(frameon=False,loc=loc)
    ax = fig.add_axes(data2figpoints(pos+0.01,-0.4,fig=fig,ax1=ax1),frameon=False,ylim=ylim)
    ax.tick_params(axis='both', which='both', labelleft='off', labelright='off',bottom='off',top='off',
                   labelbottom='off',labeltop='off',right='off',left='off')
    ax.hist(ymask,orientation='horizontal',normed=True,color=color,edgecolor='None',bins=30,alpha=0.5)
    ax.axhline(np.mean(ymask),color='red',linewidth=2)
    ax.axhline(np.median(ymask),color='k',linewidth=2,linestyle='--')

# <codecell>

def data2figpoints(x,dx,fig,ax1):
    "function to tranform data locations to relative figure coordinates (in fractions of total figure"
    flen = fig.transFigure.transform([1,1])
    bot = ax1.transAxes.transform([0,0])/flen
    top = ax1.transAxes.transform([1,1])/flen
    
    start = ax1.transData.transform([x,0])/flen
    end = ax1.transData.transform([x+dx,0])/flen
    left = start[0]
    bottom = bot[1]
    width = end[0]-start[0]
    height = top[1]-bot[1] 
    return left,bottom,width,height

# <codecell>

def plot_lin(x,y,color='b',labels=True):
    """
    function to plot on top of previous a linear fit line, with the line equation in legend.
    """
    from linfit import linfit
    from Sp_parametes import nanmasked
    c = linfit(x,y)
    

