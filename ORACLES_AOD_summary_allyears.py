#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     ORACLES_AOD_summary_allyears
# 
# Purpose:  
# 
#     Prepare an analysis for comparing AOd from all ORACLES years
#   
# Input:
# 
#     none
# 
# Output:
#    
#     plots
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - matplotlib
#     - scipy
# 
#   
# Needed Files:
# 
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc,Santa Cruz, CA, 2018-12-03
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2022-02-28
#              - added the fire counts to the elevation map
#              - updated to use newest HSRL versions
#     

# # Prepare python environment

# In[1]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
import os
matplotlib.rc_file(os.path.join(os.getcwd(),'file.rc'))
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth
from linfit import linfit
from path_utils import getpath
from plotting_utils import make_boxplot
import hdf5storage as hs
from plotting_utils import prelim
from datetime import datetime
from scipy.interpolate import UnivariateSpline
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap


# In[2]:


import scipy.stats as st
from mpl_toolkits.basemap import Basemap
import pandas as pd


# In[3]:


get_ipython().magic(u'matplotlib notebook')


# In[4]:


fp =getpath('ORACLES')#'C:/Userds/sleblan2/Research/ORACLES/'
fp


# In[5]:


vv = 'v2'


# # Load files

# ## Load the 2016 data

# In[6]:


ar6 = hs.loadmat(fp+'/aod_ict/R4/all_aod_ict_R4_2016.mat')


# In[7]:


ar6['flac'] = (ar6['qual_flag']==0)&(ar6['flag_acaod']==1)
ar6['flacr'] = (ar6['qual_flag']==0)&(ar6['flag_acaod']==1)&(ar6['fl_routine'])
ar6['flaco'] = (ar6['qual_flag']==0)&(ar6['flag_acaod']==1)&~(ar6['fl_routine'])


# In[8]:


ar6['flr'] = (ar6['qual_flag']==0) & (ar6['fl_routine'])
ar6['flo'] = (ar6['qual_flag']==0) & ~(ar6['fl_routine'])
ar6['fl'] = (ar6['qual_flag']==0)


# ## Load the 2017 data

# In[9]:


ar7 = hs.loadmat(fp+'/aod_ict_2017/R1/all_aod_ict_R1_2017.mat')


# In[10]:


ar7['flac'] = (ar7['qual_flag']==0)&(ar7['flag_acaod']==1)
ar7['flacr'] = (ar7['qual_flag']==0)&(ar7['flag_acaod']==1)&(ar7['fl_routine'])
ar7['flaco'] = (ar7['qual_flag']==0)&(ar7['flag_acaod']==1)&~(ar7['fl_routine'])


# In[11]:


ar7['flr'] = (ar7['qual_flag']==0) & (ar7['fl_routine'])
ar7['flo'] = (ar7['qual_flag']==0) & ~(ar7['fl_routine'])
ar7['fl'] = (ar7['qual_flag']==0)


# ## Load the 2018 data

# In[12]:


ar8 = hs.loadmat(fp+'/aod_ict_2018/{vv}/all_aod_ict_{vv}_2018.mat'.format(vv='R1'))


# In[13]:


ar8['flac'] = (ar8['qual_flag']==0) & (ar8['flag_acaod']==1)  
ar8['flacr'] = (ar8['qual_flag']==0) & (ar8['flag_acaod']==1)&(ar8['fl_routine']) 
ar8['flaco'] = (ar8['qual_flag']==0) & (ar8['flag_acaod']==1)&~(ar8['fl_routine']) 


# In[14]:


ar8['flr'] = (ar8['qual_flag']==0) & (ar8['fl_routine'])
ar8['flo'] = (ar8['qual_flag']==0) & (ar8['fl_routine']==False)
ar8['fl'] = (ar8['qual_flag']==0)


# ## Load the surface elevation map from SAFARI

# from : https://daac.ornl.gov/daacdata/safari2k/vegetation_wetlands/globe_dem/comp/readme_new.txt
# 
# SAFARI 2000 Digital Elevation Model, 1-km (GLOBE)
# 
# Cite this data set as follows: 
# 
# Hastings, David A., et al. 2002. SAFARI 2000 Digital Elevation Model, 1-km (GLOBE). 
# Available on-line [http://www.daac.ornl.gov/] from Oak Ridge National
# Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, U.S.A. 
# 
# 
# 
# ASCII File Information  
# 
# 
# The data is also available in ASCII Grid (so_africa_dem.dat.gz) format for 
# ArcInfo. The *.dat file contains a single ASCII array with integer values. 
# Coordinates listed below are in decimal degrees.  The ASCII file consists of 
# header information containing a set of keywords, followed by cell values in 
# row-major order. The file format is:
# 
#   NCOLS 4560
#   NROWS 4200
#   XLLCORNER 4.9958333333333
#   YLLCORNER     -34.995833333333
#   CELLSIZE      0.0083333333333333
#   NODATA_value  -9999
#   row 1
#   row 2
#   .
#   .
#   .
#   row n
# 
# Row 1 of the data is at the top of the grid, row 2 is just
# under row 1 and so on. The end of each row of data from the grid is
# terminated with a carriage return in the file. 
# 
# 
# To import this file into ArcInfo use the following command at an ARC prompt: 
# 
#     ASCIIGRID {in_ascii_file} {out_grid} {INT | FLOAT}
# 
# Arguments:
# {in_ascii_file} - the ASCII file to be converted.
# {out_grid} - the name of the grid to be created.
# {INT | FLOAT} - the data type of the output grid.
#     INT - an integer grid will be created.
#     FLOAT - a floating-point grid will be created.
# 

# In[143]:


elev_pd = pd.read_csv(fp+'surface_elevation/so_africa_dem_nohead.dat',delimiter=' ',header=None)


# In[144]:


elev = elev_pd.to_numpy()


# In[145]:


elev.shape


# In[146]:


elev_lat = np.arange(-0.0083333,-35,-0.0083333)


# In[147]:


elev_lat.shape


# In[148]:


elev_lon = np.arange(5.0,43.0,0.00833333)


# In[149]:


elev_lon.shape


# In[150]:


elev_lon[2200]


# In[239]:


plt.figure(figsize=(6,6))
plt.contourf(elev_lon,elev_lat,elev,vmin=0,vmax=3200,extend='max')
plt.colorbar(label='Surface Elevation [m]')
plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('Latitude [$^\\circ$]')
plt.savefig(fp+'plot/SA_Surface_elevation.png',dpi=600,transparent=True)


# In[152]:


elev[elev <0] = np.nan


# In[153]:


strip_elev = np.nanmax(elev[:,0:2200],axis=1)
mean_elev = np.nanmean(elev[:,0:2200],axis=1)


# In[154]:


plt.figure()
plt.plot(strip_elev)
plt.plot(mean_elev)


# In[57]:


plt.figure()
plt.plot(mean_elev,elev_lat)


# ## Load the fire counts

# In[32]:


import load_utils as lu


# In[27]:


fh = os.listdir(fp+'fire/')


# In[30]:


fh = [hh for hh in fh if hh.endswith('hdf')]
fh.sort()


# In[135]:


fi_geo = {}


# In[140]:


nfh = len(fh)
fi_2016,fi_2017,fi_2018 = [],[],[]
fi_doy_2016,fi_doy_2017,fi_doy_2018 = [],[],[]
for jj,ffh in enumerate(fh):
    print('{}/ total {}'.format(jj,nfh))
    grid = ffh.split('.')[2]
    if grid in fi_geo.keys():
        get_geo = False
    else:
        get_geo = True
    fi,fi_dict = lu.load_hdf(fp+'fire/'+ffh,all_values=True,i_subdata=0,get_geo=get_geo,verbose=False)
    if grid in fi_geo.keys():
        fi['FireMask_lat'] = fi_geo[grid]['FireMask_lat']
        fi['FireMask_lon'] = fi_geo[grid]['FireMask_lon']
    else:
        fi_geo[grid] = {}
        fi_geo[grid]['FireMask_lat'] = fi['FireMask_lat']
        fi_geo[grid]['FireMask_lon'] = fi['FireMask_lon']
        nul = fi.pop('QA_lon',None)
        nul = fi.pop('QA_lat',None)
    doy = float(ffh.split('.')[1][-3:])
    if '2016' in ffh.split('.')[1]:
        fi_2016.append(fi)
        fi_doy_2016.append(doy)
    if '2017' in ffh.split('.')[1]:
        fi_2017.append(fi)
        fi_doy_2017.append(doy)
    if '2018' in ffh.split('.')[1]:
        fi_2018.append(fi)
        fi_doy_2018.append(doy)


# In[142]:


print(fi_dict['FireMask']['legend'])


# ### Relate fire counts to surface elevation

# In[183]:


from scipy.interpolate import RectBivariateSpline, interp2d


# In[227]:


elevm = np.ma.array(elev,mask=~np.isfinite(elev))


# In[245]:


elev_fx = interp2d(elev_lon,elev_lat,elevm.filled(0),kind='linear',bounds_error=False)


# In[252]:


for fi in fi_2016:
    fi['elev'] = elev_fx(fi['FireMask_lon'][:,0],fi['FireMask_lat'][0,:])


# In[248]:


for fi in fi_2017:
    fi['elev'] = elev_fx(fi['FireMask_lon'][:,0],fi['FireMask_lat'][0,:])


# In[249]:


for fi in fi_2018:
    fi['elev'] = elev_fx(fi['FireMask_lon'][:,0],fi['FireMask_lat'][0,:])


# ### Calculate the average height of fires per latitude

# In[317]:


lat_bins = np.linspace(0.5,-24.5,26)
center_lat = np.array([(la+lat_bins[i+1])/2.0 for i,la in enumerate(lat_bins[:-1])])


# In[311]:


fires_elev_perlat_2016  = []
fires_num_perlat_2016  = []
for i,la in enumerate(lat_bins[:-1]):
    la1,la2 = la,lat_bins[i+1]
    e2016 = []
    for fi in fi_2016:
        ig = (fi['FireMask_lat']>la2) & (fi['FireMask_lat']<la1) & (fi['FireMask']>6)
        e2016.extend(fi['elev'][ig])
    fires_elev_perlat_2016.append(e2016)
    fires_num_perlat_2016.append(len(e2016))


# In[312]:


fires_elev_perlat_2017  = []
fires_num_perlat_2017  = []
for i,la in enumerate(lat_bins[:-1]):
    la1,la2 = la,lat_bins[i+1]
    e2017 = []
    for fi in fi_2017:
        ig = (fi['FireMask_lat']>la2) & (fi['FireMask_lat']<la1) & (fi['FireMask']>6)
        e2017.extend(fi['elev'][ig])
    fires_elev_perlat_2017.append(e2017)
    fires_num_perlat_2017.append(len(e2017))


# In[313]:


fires_elev_perlat_2018  = []
fires_num_perlat_2018  = []
for i,la in enumerate(lat_bins[:-1]):
    la1,la2 = la,lat_bins[i+1]
    e2018 = []
    for fi in fi_2018:
        ig = (fi['FireMask_lat']>la2) & (fi['FireMask_lat']<la1) & (fi['FireMask']>6)
        e2018.extend(fi['elev'][ig])
    fires_elev_perlat_2018.append(e2018)
    fires_num_perlat_2018.append(len(e2018))


# In[314]:


fires_num_perlat_2016


# In[406]:


plt.figure(figsize=(5,8))
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues

bp6 = plt.boxplot(fires_elev_perlat_2016,positions=center_lat-0.15,vert=False,
                 showfliers=False,widths=0.2,showmeans=True,patch_artist=True)
set_box_whisker_color(gr,bp6,fires_num_perlat_2016)
bp8 = plt.boxplot(fires_elev_perlat_2018,positions=center_lat+0.15,vert=False,
                 showfliers=False,widths=0.2,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,fires_num_perlat_2018)
bp7 = plt.boxplot(fires_elev_perlat_2017,positions=center_lat,vert=False,
                 showfliers=False,widths=0.2,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,fires_num_perlat_2017)

plt.xlabel('Surface elevation of fires [m]')
plt.ylabel('Latitude [$^\\circ$]')
plt.title('Elevation of fires during ORACLES')

plt.legend([bp6['boxes'][5],bp7['boxes'][18],bp8['boxes'][18],bp6['means'][0],bp6['medians'][0],bp6['boxes'][0],
            bp6['whiskers'][0]],
           ['2016','2017','2018','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(fires_num_perlat_2016)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(fires_num_perlat_2016)
scalarmapbr = plt.cm.ScalarMappable(cmap=br)
scalarmapbr.set_array(fires_num_perlat_2016)

cbaxesgr = plt.gcf().add_axes([0.60, 0.43, 0.015, 0.2])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbg.outline.set_visible(False)
cbaxesbl = plt.gcf().add_axes([0.62, 0.43, 0.015, 0.2])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbb.outline.set_visible(False)
cbaxesbr = plt.gcf().add_axes([0.64, 0.43, 0.015, 0.2])
cbr = plt.colorbar(scalarmapbr,cax=cbaxesbr)
cbr.outline.set_visible(False)
cbg.set_ticks([0,10000,20000,30000,40000,50000]),cbg.set_ticklabels(['','','','',''])
cbb.set_ticks([0,10000,20000,30000,40000,50000]),cbb.set_ticklabels(['','','','',''])
cbr.set_ticks([0,10000,20000,30000,40000,50000])
cbaxesgr.yaxis.set_ticks_position('right'),cbaxesbl.yaxis.set_ticks_position('right'),cbaxesbr.yaxis.set_ticks_position('right')
cbaxesgr.spines['right'].set_visible(False), cbaxesbl.spines['left'].set_visible(False)
cbaxesbr.text(8.0,0.5,'Number of fires',rotation=-90,verticalalignment='center')


plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLES_fire_elevation_average.png',
            transparent=True,dpi=500)


# In[407]:


mean_fire_elev6 = [aaa.get_data()[0][0] for aaa in bp6['means']]
mean_fire_elev7 = [aaa.get_data()[0][0] for aaa in bp7['means']]
mean_fire_elev8 = [aaa.get_data()[0][0] for aaa in bp8['means']]


# In[408]:


median_fire_elev6 = [aaa.get_data()[0][0] for aaa in bp6['medians']]
median_fire_elev7 = [aaa.get_data()[0][0] for aaa in bp7['medians']]
median_fire_elev8 = [aaa.get_data()[0][0] for aaa in bp8['medians']]


# In[427]:


plt.figure()
plt.plot(mean_fire_elev6,center_lat,'.-',color='tab:red',label='2016')
plt.plot(mean_fire_elev7,center_lat,'.-',color='tab:green',label='2017')
plt.plot(mean_fire_elev8,center_lat,'.-',color='tab:blue',label='2018')
plt.xlabel('Mean surface elevation of fires [m]')
plt.ylabel('Latitude [$^\\circ$]')
plt.title('Elevation of fires during ORACLES')
plt.legend(frameon=False)
plt.savefig(fp+'plot_all/ORACLES_fire_elevation_average_line.png',
            transparent=True,dpi=500)


# ## Load HSRL bins

# In[15]:


vhsrl = 'v2'


# In[16]:


hsrl = hs.loadmat(fp+'data_other/HSRL/ORACLES_binned_HSRL_allyears_{}.mat'.format(vhsrl))


# In[17]:


hsrl.keys()


# # Calculate the relationship between AOD 550 nm and others

# For Kerry Meyer, 
# 
# From email (Jan 20, 2020):" calculate spectral AOD at the channels I use (0.47, 0.55, 0.67, 0.86, 1.24, 2.1µm). "

# ## Save for 2016

# In[62]:


kma6 = dict()


# In[63]:


kma6[u'AOD0470_rel'] = ar6['AOD0470'][ar6['flac']]/ar6['AOD0550'][ar6['flac']]
kma6[u'AOD0550_abs'] = ar6['AOD0550'][ar6['flac']]
kma6[u'AOD0675_rel'] = ar6['AOD0675'][ar6['flac']]/ar6['AOD0550'][ar6['flac']]
kma6[u'AOD0865_rel'] = ar6['AOD0865'][ar6['flac']]/ar6['AOD0550'][ar6['flac']]
kma6[u'AOD1236_rel'] = ar6['AOD1236'][ar6['flac']]/ar6['AOD0550'][ar6['flac']]


# In[64]:


len(kma6['AOD0470_rel'])


# In[65]:


kma6[u'lat'] = ar6['Latitude'][ar6['flac']]
kma6[u'lon'] = ar6['Longitude'][ar6['flac']]
kma6[u'alt'] = ar6['GPS_Alt'][ar6['flac']]


# In[25]:


days6 = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']
ar6['daysd'] = [days6[i] for i in ar6['days'].astype(int)]
ar6['ndtime'] = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8]),int(ar6['Start_UTC'][i]),
                          int((ar6['Start_UTC'][i]-float(int(ar6['Start_UTC'][i])))*60)) for i,d in enumerate(ar6['daysd'])]
ar6['ndtimes'] = np.array(ar6['ndtime'])


# In[66]:


kma6[u'doy'] = np.array([ad.timetuple().tm_yday+(ad.hour+(ad.minute+ad.second/60.0)/60.0)/24.0 for ad in ar6['ndtimes'][ar6['flac']]])


# In[67]:


kma6.keys()


# In[74]:


kma6[u'details'] = '''
This is a combined product from all 4STAR 2016 ORACLES. The AOD are quality flagged, and only for the above cloud ACAOD. 
Based on the R3 ict files.
AODXXXX_rel - represent the ratio of AOD of the wavelength XXXX (in nm) to the AOD 550 nm
AOD0550_abs - the absolute value of AOD
lon - Longitude
lat - Latitude
alt - GPS Altitude (meters)
doy - fractional day

Produced by Samuel LeBlanc on {}
'''.format(datetime.now())


# In[75]:


hs.savemat(fp+'ORACLES2016_4STAR_ACAOD_relative_to_550.mat',kma6)


# ## Save for 2017

# In[76]:


kma7 = dict()
kma7[u'AOD0470_rel'] = ar7['AOD0470'][ar7['flac']]/ar7['AOD0550'][ar7['flac']]
kma7[u'AOD0550_abs'] = ar7['AOD0550'][ar7['flac']]
kma7[u'AOD0675_rel'] = ar7['AOD0675'][ar7['flac']]/ar7['AOD0550'][ar7['flac']]
kma7[u'AOD0865_rel'] = ar7['AOD0865'][ar7['flac']]/ar7['AOD0550'][ar7['flac']]
kma7[u'AOD1236_rel'] = ar7['AOD1236'][ar7['flac']]/ar7['AOD0550'][ar7['flac']]
kma7[u'lat'] = ar7['Latitude'][ar7['flac']]
kma7[u'lon'] = ar7['Longitude'][ar7['flac']]
kma7[u'alt'] = ar7['GPS_Alt'][ar7['flac']]


# In[26]:


days7 = ['20170801','20170802','20170807','20170809', '20170812','20170813','20170815','20170817','20170818','20170819','20170821',
        '20170824','20170826','20170828','20170830','20170831','20170902','20170903','20170904']
ar7['daysd'] = [days7[i] for i in ar7['days'].astype(int)]
ar7['ndtime'] = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8]),int(ar7['Start_UTC'][i]),
                          int((ar7['Start_UTC'][i]-float(int(ar7['Start_UTC'][i])))*60)) for i,d in enumerate(ar7['daysd'])]
ar7['ndtimes'] = np.array(ar7['ndtime'])


# In[27]:


kma7[u'doy'] = np.array([ad.timetuple().tm_yday+(ad.hour+(ad.minute+ad.second/60.0)/60.0)/24.0 for ad in ar7['ndtimes'][ar7['flac']]])


# In[90]:


kma7[u'details'] = '''
This is a combined product from all 4STAR 2017 ORACLES. The AOD are quality flagged, and only for the above cloud ACAOD. 
Based on the R1 ict files.
AODXXXX_rel - represent the ratio of AOD of the wavelength XXXX (in nm) to the AOD 550 nm
AOD0550_abs - the absolute value of AOD
lon - Longitude
lat - Latitude
alt - GPS Altitude (meters)
doy - fractional day

Produced by Samuel LeBlanc on {}
'''.format(datetime.now())


# In[91]:


kma7.keys()


# In[92]:


hs.savemat(fp+'ORACLES2017_4STAR_ACAOD_relative_to_550.mat',kma7)


# ## Save for 2018

# In[84]:


kma8 = dict()
kma8[u'AOD0470_rel'] = ar8['AOD0470'][ar8['flac']]/ar8['AOD0550'][ar8['flac']]
kma8[u'AOD0550_abs'] = ar8['AOD0550'][ar8['flac']]
kma8[u'AOD0675_rel'] = ar8['AOD0675'][ar8['flac']]/ar8['AOD0550'][ar8['flac']]
kma8[u'AOD0865_rel'] = ar8['AOD0865'][ar8['flac']]/ar8['AOD0550'][ar8['flac']]
kma8[u'AOD1236_rel'] = ar8['AOD1236'][ar8['flac']]/ar8['AOD0550'][ar8['flac']]
kma8[u'lat'] = ar8['Latitude'][ar8['flac']]
kma8[u'lon'] = ar8['Longitude'][ar8['flac']]
kma8[u'alt'] = ar8['GPS_Alt'][ar8['flac']]


# In[28]:


days8 = ['20180921','20180922','20180924','20180927','20180930','20181002','20181003','20181005','20181007','20181010','20181012',
        '20181015','20181017','20181019','20181021','20181023','20181025','20181026','20181027']
ar8['daysd'] = [days8[i] for i in ar8['days'].astype(int)]
ar8['ndtime'] = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8]),int(ar8['Start_UTC'][i]),
                          int((ar8['Start_UTC'][i]-float(int(ar8['Start_UTC'][i])))*60)) for i,d in enumerate(ar8['daysd'])]
ar8['ndtimes'] = np.array(ar8['ndtime'])


# In[86]:


kma8[u'doy'] = np.array([ad.timetuple().tm_yday+(ad.hour+(ad.minute+ad.second/60.0)/60.0)/24.0 for ad in ar8['ndtimes'][ar8['flac']]])


# In[87]:


kma8[u'details'] = '''
This is a combined product from all 4STAR 2018 ORACLES. The AOD are quality flagged, and only for the above cloud ACAOD. 
Based on the R1 ict files.
AODXXXX_rel - represent the ratio of AOD of the wavelength XXXX (in nm) to the AOD 550 nm
AOD0550_abs - the absolute value of AOD
lon - Longitude
lat - Latitude
alt - GPS Altitude (meters)
doy - fractional day

Produced by Samuel LeBlanc on {}
'''.format(datetime.now())


# In[88]:


kma8.keys()


# In[89]:


hs.savemat(fp+'ORACLES2018_4STAR_ACAOD_relative_to_550.mat',kma8)


# # Now plot the data together

# In[25]:


len(ar6['AOD0501'][ar6['fl']]), len(ar7['AOD0501'][ar7['fl']]), len(ar8['AOD0501'][ar8['fl']])


# In[26]:


np.nanmean(ar6['AOD0501'][ar6['fl']]),np.nanmean(ar7['AOD0501'][ar7['fl']]),np.nanmean(ar8['AOD0501'][ar8['fl']])


# ## plot the histograms

# In[17]:


plt.figure()
plt.hist(ar8['AOD0501'][ar8['fl']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='b',label='2018')
plt.hist(ar7['AOD0501'][ar7['fl']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='g',label='2017')
plt.hist(ar6['AOD0501'][ar6['fl']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='r',label='2016')

#plt.hist(ar['AOD0501'][ar['fl1']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='b',label='below 600 m')
#plt.hist(ar['AOD0501'][ar['fl3']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='y',label='800-2200 m')
#plt.hist(ar['AOD0501'][ar['fl2']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='r',label='above 1800 m')
#plt.yscale('log')
plt.axvline(x=np.nanmean(ar8['AOD0501'][ar8['fl']]),ls='-',color='k',lw=2.5,label='Mean')
plt.axvline(x=np.nanmedian(ar8['AOD0501'][ar8['fl']]),ls='--',color='k',label='Median')
plt.axvline(x=np.nanmean(ar8['AOD0501'][ar8['fl']]),ls='-',color='b',lw=2.5)
plt.axvline(x=np.nanmedian(ar8['AOD0501'][ar8['fl']]),ls='--',color='b')

plt.axvline(x=np.nanmean(ar6['AOD0501'][ar6['fl']]),ls='-',color='r',lw=2.5)
plt.axvline(x=np.nanmedian(ar6['AOD0501'][ar6['fl']]),ls='--',color='coral')

plt.axvline(x=np.nanmean(ar7['AOD0501'][ar7['fl']]),ls='-',color='g',lw=2.5)
plt.axvline(x=np.nanmedian(ar7['AOD0501'][ar7['fl']]),ls='--',color='lightgreen')

plt.xlabel('AOD at 501 nm')
plt.ylabel('Normalized Counts')
plt.grid()
plt.title('AOD distribution cloud filtered')
#prelim()
plt.legend(frameon=False)
plt.savefig(fp+'plot_all/AOD_normed_histogram_2018_2017_2016_R1.png',dpi=600,transparent=True)


# In[18]:


plt.figure()
plt.hist(ar8['AOD0501'][ar8['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='b',label='2018')
plt.hist(ar7['AOD0501'][ar7['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='g',label='2017')
plt.hist(ar6['AOD0501'][ar6['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='r',label='2016')

#plt.hist(ar['AOD0501'][ar['fl1']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='b',label='below 600 m')
#plt.hist(ar['AOD0501'][ar['fl3']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='y',label='800-2200 m')
#plt.hist(ar['AOD0501'][ar['fl2']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='r',label='above 1800 m')
#plt.yscale('log')
plt.axvline(x=np.nanmean(ar8['AOD0501'][ar8['flac']]),ls='-',color='k',lw=2.5,label='Mean')
plt.axvline(x=np.nanmedian(ar8['AOD0501'][ar8['flac']]),ls='--',color='k',label='Median')
plt.axvline(x=np.nanmean(ar8['AOD0501'][ar8['flac']]),ls='-',color='b',lw=2.5)
plt.axvline(x=np.nanmedian(ar8['AOD0501'][ar8['flac']]),ls='--',color='b')

plt.axvline(x=np.nanmean(ar6['AOD0501'][ar6['flac']]),ls='-',color='r',lw=2.5)
plt.axvline(x=np.nanmedian(ar6['AOD0501'][ar6['flac']]),ls='--',color='coral')

plt.axvline(x=np.nanmean(ar7['AOD0501'][ar7['flac']]),ls='-',color='g',lw=2.5)
plt.axvline(x=np.nanmedian(ar7['AOD0501'][ar7['flac']]),ls='--',color='lightgreen')

plt.text(np.nanmean(ar6['AOD0501'][ar6['flac']])+0.01,5.5,'{:0.2f}'.format(np.nanmean(ar6['AOD0501'][ar6['flac']])),color='r')
plt.text(np.nanmean(ar7['AOD0501'][ar7['flac']])+0.01,5,'{:0.2f}'.format(np.nanmean(ar7['AOD0501'][ar7['flac']])),color='g')
plt.text(np.nanmean(ar8['AOD0501'][ar8['flac']])+0.01,4.5,'{:0.2f}'.format(np.nanmean(ar8['AOD0501'][ar8['flac']])),color='b')



plt.xlabel('Above Cloud AOD at 501 nm ')
plt.ylabel('Normalized Counts')
plt.grid()
plt.title('Above clouds AOD distribution cloud filtered')
#prelim()
plt.legend(frameon=False)
plt.savefig(fp+'plot_all/AOD_Above_cloud_normed_histogram_2018_2017_2016_R1.png',dpi=600,transparent=True)


# In[19]:


plt.figure()
plt.hist(ar8['AOD0532'][ar8['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='b',label='2018')
plt.hist(ar7['AOD0532'][ar7['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='g',label='2017')
plt.hist(ar6['AOD0532'][ar6['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='r',label='2016')

#plt.hist(ar['AOD0501'][ar['fl1']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='b',label='below 600 m')
#plt.hist(ar['AOD0501'][ar['fl3']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='y',label='800-2200 m')
#plt.hist(ar['AOD0501'][ar['fl2']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='r',label='above 1800 m')
#plt.yscale('log')
plt.axvline(x=np.nanmean(ar8['AOD0532'][ar8['flac']]),ls='-',color='k',lw=2.5,label='Mean')
plt.axvline(x=np.nanmedian(ar8['AOD0532'][ar8['flac']]),ls='--',color='k',label='Median')
plt.axvline(x=np.nanmean(ar8['AOD0532'][ar8['flac']]),ls='-',color='b',lw=2.5)
plt.axvline(x=np.nanmedian(ar8['AOD0532'][ar8['flac']]),ls='--',color='b')

plt.axvline(x=np.nanmean(ar6['AOD0532'][ar6['flac']]),ls='-',color='r',lw=2.5)
plt.axvline(x=np.nanmedian(ar6['AOD0532'][ar6['flac']]),ls='--',color='coral')

plt.axvline(x=np.nanmean(ar7['AOD0532'][ar7['flac']]),ls='-',color='g',lw=2.5)
plt.axvline(x=np.nanmedian(ar7['AOD0532'][ar7['flac']]),ls='--',color='lightgreen')

plt.text(np.nanmean(ar6['AOD0532'][ar6['flac']])+0.01,4.2,'{:0.2f}'.format(np.nanmean(ar6['AOD0532'][ar6['flac']])),color='r')
plt.text(np.nanmean(ar7['AOD0532'][ar7['flac']])+0.01,4,'{:0.2f}'.format(np.nanmean(ar7['AOD0532'][ar7['flac']])),color='g')
plt.text(np.nanmean(ar8['AOD0532'][ar8['flac']])+0.01,3.5,'{:0.2f}'.format(np.nanmean(ar8['AOD0532'][ar8['flac']])),color='b')



plt.xlabel('Above Cloud AOD at 532 nm ')
plt.ylabel('Normalized Counts')
plt.grid()
plt.title('Above clouds AOD distribution cloud filtered')
#prelim()
plt.legend(frameon=False)
plt.savefig(fp+'plot_all/AOD_532_Above_cloud_normed_histogram_2018_2017_2016_R1.png',dpi=600,transparent=True)


# In[19]:


plt.figure()
plt.hist(ar8['AOD0501'][ar8['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='b',label='2018')
plt.hist(ar8['AOD0501'][ar8['flacr']],bins=30,range=(0,1.0),alpha=0.5,
         normed=True,edgecolor='None',color='aqua',label='routine 2018')
plt.hist(ar7['AOD0501'][ar7['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='g',label='2017')
plt.hist(ar7['AOD0501'][ar7['flacr']],bins=30,range=(0,1.0),alpha=0.5,
         normed=True,edgecolor='None',color='lime',label='routine 2017')
plt.hist(ar6['AOD0501'][ar6['flac']],bins=30,range=(0,1.0),alpha=0.5,normed=True,edgecolor='None',color='r',label='2016')
plt.hist(ar6['AOD0501'][ar6['flacr']],bins=30,range=(0,1.0),alpha=0.5,
         normed=True,edgecolor='None',color='orange',label='routine 2016')

plt.axvline(x=np.nanmean(ar8['AOD0501'][ar8['flac']]),ls='-',color='k',lw=2.5,label='Mean')
plt.axvline(x=np.nanmedian(ar8['AOD0501'][ar8['flac']]),ls='--',color='k',label='Median')
plt.axvline(x=np.nanmean(ar8['AOD0501'][ar8['flac']]),ls='-',color='b',lw=2.5)
plt.axvline(x=np.nanmedian(ar8['AOD0501'][ar8['flac']]),ls='--',color='b')
plt.axvline(x=np.nanmean(ar6['AOD0501'][ar6['flac']]),ls='-',color='r',lw=2.5)
plt.axvline(x=np.nanmedian(ar6['AOD0501'][ar6['flac']]),ls='--',color='r')
plt.axvline(x=np.nanmean(ar7['AOD0501'][ar7['flac']]),ls='-',color='g',lw=2.5)
plt.axvline(x=np.nanmedian(ar7['AOD0501'][ar7['flac']]),ls='--',color='g')
plt.axvline(x=np.nanmean(ar8['AOD0501'][ar8['flacr']]),ls='-',color='aqua',lw=2.5)
plt.axvline(x=np.nanmedian(ar8['AOD0501'][ar8['flacr']]),ls='--',color='aqua')
plt.axvline(x=np.nanmean(ar6['AOD0501'][ar6['flacr']]),ls='-',color='orange',lw=2.5)
plt.axvline(x=np.nanmedian(ar6['AOD0501'][ar6['flacr']]),ls='--',color='orange')
plt.axvline(x=np.nanmean(ar7['AOD0501'][ar7['flacr']]),ls='-',color='lime',lw=2.5)
plt.axvline(x=np.nanmedian(ar7['AOD0501'][ar7['flacr']]),ls='--',color='lime')

plt.text(np.nanmean(ar6['AOD0501'][ar6['flac']])+0.01,7.5,'{:0.2f}'.format(np.nanmean(ar6['AOD0501'][ar6['flac']])),color='r')
plt.text(np.nanmean(ar6['AOD0501'][ar6['flacr']])+0.01,7.1,'{:0.2f}'.format(np.nanmean(ar6['AOD0501'][ar6['flacr']])),color='orange')
plt.text(np.nanmean(ar7['AOD0501'][ar7['flac']])+0.01,6,'{:0.2f}'.format(np.nanmean(ar7['AOD0501'][ar7['flac']])),color='g')
plt.text(np.nanmean(ar7['AOD0501'][ar7['flacr']])-0.06,6.5,'{:0.2f}'.format(np.nanmean(ar7['AOD0501'][ar7['flacr']])),color='lime')
plt.text(np.nanmean(ar8['AOD0501'][ar8['flac']])+0.01,5,'{:0.2f}'.format(np.nanmean(ar8['AOD0501'][ar8['flac']])),color='b')
plt.text(np.nanmean(ar8['AOD0501'][ar8['flacr']])+0.01,7.2,'{:0.2f}'.format(np.nanmean(ar8['AOD0501'][ar8['flacr']])),color='aqua')


plt.xlabel('Above Cloud AOD at 501 nm ')
plt.ylabel('Normalized Counts')
plt.grid()
plt.title('Above clouds AOD distribution cloud filtered')
prelim()
plt.legend(frameon=False)
plt.savefig(fp+'plot_all/AOD_Above_cloud_normed_histogram_2018_2017_2016_withroutine_R1.png',dpi=500,transparent=True)


# ## plot against lat and lon

# In[227]:


lim = np.linspace(-15,7,12)
pos = np.array([(l+lim[i+1])/2.0 for i,l in enumerate(lim[0:-1])])


# In[228]:


lima = np.linspace(-15,0,11)
posa = np.array([(l+lima[i+1])/2.0 for i,l in enumerate(lima[0:-1])])


# In[229]:


lim2 = np.linspace(-15,15,16)
pos2 = np.array([(l+lim2[i+1])/2.0 for i,l in enumerate(lim2[0:-1])])


# In[230]:


lim2a = np.linspace(-22,-7,13)
pos2a = np.array([(l+lim2a[i+1])/2.0 for i,l in enumerate(lim2a[0:-1])])


# In[25]:


plt.figure()
plt.plot(ar8['Longitude'][ar8['flac']],ar8['AOD0501'][ar8['flac']],'.',color='lightblue',alpha=0.5,
         markersize=0.4,label='All 2018 data')
plt.xlim([-15,15])
make_boxplot(ar8['AOD0501'][ar8['flac']],ar8['Longitude'][ar8['flac']],lim,pos,color='b',label='2018 Statistics',fliers_off=True)


plt.plot(ar7['Longitude'][ar7['flac']],ar7['AOD0501'][ar7['flac']],'.',color='lightgreen',alpha=0.5,
         markersize=0.4,label='All 2017 data')
plt.xlim([-15,15])
make_boxplot(ar7['AOD0501'][ar7['flac']],ar7['Longitude'][ar7['flac']],lim,pos,color='g',label='2017 Statistics',fliers_off=True)
plt.plot(ar6['Longitude'][ar6['flac']],ar6['AOD0501'][ar6['flac']],'.',color='lightcoral',alpha=0.5,
         markersize=0.4,label='All 2016 data')
make_boxplot(ar6['AOD0501'][ar6['flac']],ar6['Longitude'][ar6['flac']],lim2,pos2,color='red',
             label='2016 Statistics',fliers_off=True)

plt.xlabel('Longitude [$^\\circ$]')
plt.ylabel('4STAR AOD at 501 nm between 600 - 1800 m')
plt.grid()
plt.xlim([-15,15])
plt.ylim([0,1])
plt.legend(frameon=False,numpoints=1,loc=2)
#prelim()

plt.title('ORACLES above cloud AOD')
plt.savefig(fp+'plot_all/AOD_longitude_2016_2017_2018_R1.png',dpi=600,transparent=True)


# In[33]:


plt.figure()
plt.plot(ar7['Latitude'][ar7['flac']],ar7['AOD0501'][ar7['flac']],'.',color='lightsteelblue',alpha=0.5,
         markersize=0.4,label='All 2017 data')
plt.xlim([-22,0])
make_boxplot(ar7['AOD0501'][ar7['flac']],ar7['Latitude'][ar7['flac']],lima,posa,color='blue',label='2017 Statistics',fliers_off=True)
plt.plot(ar7['Latitude'][ar7['flacr']],ar7['AOD0501'][ar7['flacr']],'.',color='lightyellow',alpha=0.5,
         markersize=0.4,label='Routine 2017')
make_boxplot(ar7['AOD0501'][ar7['flacr']],ar7['Latitude'][ar7['flacr']],lima,posa,color='yellow',label='2017 Routine Statistics',
             fliers_off=True)

plt.plot(ar6['Latitude'][ar6['flac']],ar6['AOD0501'][ar6['flac']],'.',color='lightcoral',alpha=0.5,
         markersize=0.4,label='All 2016 data')
make_boxplot(ar6['AOD0501'][ar6['flac']],ar6['Latitude'][ar6['flac']],lim2a,pos2a,color='red',
             label='2016 Statistics',fliers_off=True)

plt.xlabel('Latitude [$^\\circ$]')
plt.ylabel('4STAR AOD at 501 nm between 600 - 1800 m')
plt.grid()
plt.xlim([-22,0])
plt.ylim(0,1)
plt.legend(frameon=False,numpoints=1,loc=2)
#prelim()

plt.title('ORACLES above cloud AOD')
plt.savefig(fp+'plot_all/AOD_latitude_2016_2017.png',dpi=600,transparent=True)


# In[26]:


plt.figure()
plt.plot(ar8['Latitude'][ar8['flac']],ar8['AOD0501'][ar8['flac']],'.',color='lightsteelblue',alpha=0.0,
         markersize=0.4)
plt.xlim([-22,0])
make_boxplot(ar8['AOD0501'][ar8['flac']],ar8['Latitude'][ar8['flac']],lima,posa,color='b',label='2018',fliers_off=True)
#plt.plot(ar8['Latitude'][ar8['flacr']],ar8['AOD0501'][ar8['flacr']],'.',color='aqua',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar8['AOD0501'][ar8['flacr']],ar8['Latitude'][ar8['flacr']],lima,posa,color='aqua',label='2018 Routine',
             fliers_off=True)

#plt.plot(ar7['Latitude'][ar7['flac']],ar7['AOD0501'][ar7['flac']],'.',color='g',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar7['AOD0501'][ar7['flac']],ar7['Latitude'][ar7['flac']],lima,posa,color='g',label='2017',
             fliers_off=True)
#plt.plot(ar7['Latitude'][ar7['flacr']],ar7['AOD0501'][ar7['flacr']],'.',color='lime',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar7['AOD0501'][ar7['flacr']],ar7['Latitude'][ar7['flacr']],lima,posa,color='lime',label='2017 Routine',
             fliers_off=True)

#plt.plot(ar6['Latitude'][ar6['flac']],ar6['AOD0501'][ar6['flac']],'.',color='lightcoral',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar6['AOD0501'][ar6['flac']],ar6['Latitude'][ar6['flac']],lim2a,pos2a,color='red',
             label='2016',fliers_off=True)
#plt.plot(ar6['Latitude'][ar6['flacr']],ar6['AOD0501'][ar6['flacr']],'.',color='orange',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar6['AOD0501'][ar6['flacr']],ar6['Latitude'][ar6['flacr']],lim2a,pos2a,color='orange',
             label='2016 routine',fliers_off=True)

plt.xlabel('Latitude [$^\\circ$]')
plt.ylabel('4STAR AOD at 501 nm')
plt.grid()
plt.xlim([-22,0])
plt.ylim(0,1)
plt.legend(frameon=False,numpoints=1,loc=2)
#prelim()

plt.title('ORACLES above cloud AOD')
plt.savefig(fp+'plot_all/AOD_latitude_2016_2017_2018_R1.png',dpi=600,transparent=True)


# In[27]:


plt.figure(figsize=(8,3))
plt.plot(ar8['Latitude'][ar8['flaco']],ar8['AOD0501'][ar8['flaco']],'.',color='lightsteelblue',alpha=0.0,
         markersize=0.4)
plt.xlim([-22,0])
#make_boxplot(ar8['AOD0501'][ar8['flac']],ar8['Latitude'][ar8['flac']],lima,posa,color='b',label='2018',fliers_off=True)
#plt.plot(ar8['Latitude'][ar8['flacr']],ar8['AOD0501'][ar8['flacr']],'.',color='aqua',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar8['AOD0501'][ar8['flaco']],ar8['Latitude'][ar8['flaco']],lima,posa,color='aqua',label='2018',
             fliers_off=True)

#plt.plot(ar7['Latitude'][ar7['flac']],ar7['AOD0501'][ar7['flac']],'.',color='g',alpha=0.5,
#         markersize=0.4)
#make_boxplot(ar7['AOD0501'][ar7['flac']],ar7['Latitude'][ar7['flac']],lima,posa,color='g',label='2017 coastal',
#             fliers_off=True)
#plt.plot(ar7['Latitude'][ar7['flacr']],ar7['AOD0501'][ar7['flacr']],'.',color='lime',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar7['AOD0501'][ar7['flacr']],ar7['Latitude'][ar7['flacr']],lima,posa,color='lime',label='2017 routine',
             fliers_off=True)

#plt.plot(ar6['Latitude'][ar6['flac']],ar6['AOD0501'][ar6['flac']],'.',color='lightcoral',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar6['AOD0501'][ar6['flaco']],ar6['Latitude'][ar6['flaco']],lim2a,pos2a,color='lightcoral',
             label='2016',fliers_off=True)
#plt.plot(ar6['Latitude'][ar6['flacr']],ar6['AOD0501'][ar6['flacr']],'.',color='orange',alpha=0.5,
#         markersize=0.4)
#make_boxplot(ar6['AOD0501'][ar6['flacr']],ar6['Latitude'][ar6['flacr']],lim2a,pos2a,color='orange',
#             label='2016 routine',fliers_off=True)

plt.xlabel('Latitude [$^\\circ$]')
plt.ylabel('4STAR AOD at 501 nm')
plt.grid()
plt.xlim([-22,0])
plt.ylim(0,1)
plt.legend(frameon=False,numpoints=1,loc=2)
#prelim()

plt.title('ORACLES above cloud AOD - coastal')
plt.savefig(fp+'plot_all/AOD_coastal_latitude_2016_2017_2018_R1.png',dpi=600,transparent=True)


# In[ ]:


plt.figure(figsize=(8,3))
plt.plot(ar8['Latitude'][ar8['flac']],ar8['AOD0501'][ar8['flac']],'.',color='lightsteelblue',alpha=0.0,
         markersize=0.4)
plt.xlim([-22,0])
make_boxplot(ar8['AOD0501'][ar8['flacr']],ar8['Latitude'][ar8['flacr']],lima,posa,color='b',label='2018 routine',fliers_off=True)
#plt.plot(ar8['Latitude'][ar8['flacr']],ar8['AOD0501'][ar8['flacr']],'.',color='aqua',alpha=0.5,
#         markersize=0.4)
#make_boxplot(ar8['AOD0501'][ar8['flacr']],ar8['Latitude'][ar8['flacr']],lima,posa,color='aqua',label='2018',
#             fliers_off=True)

#plt.plot(ar7['Latitude'][ar7['flac']],ar7['AOD0501'][ar7['flac']],'.',color='g',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar7['AOD0501'][ar7['flaco']],ar7['Latitude'][ar7['flaco']],lima,posa,color='g',label='2017',
             fliers_off=True)
#plt.plot(ar7['Latitude'][ar7['flacr']],ar7['AOD0501'][ar7['flacr']],'.',color='lime',alpha=0.5,
#         markersize=0.4)
#make_boxplot(ar7['AOD0501'][ar7['flacr']],ar7['Latitude'][ar7['flacr']],lima,posa,color='lime',label='2017',
#             fliers_off=True)

#plt.plot(ar6['Latitude'][ar6['flac']],ar6['AOD0501'][ar6['flac']],'.',color='lightcoral',alpha=0.5,
#         markersize=0.4)
#make_boxplot(ar6['AOD0501'][ar6['flac']],ar6['Latitude'][ar6['flac']],lim2a,pos2a,color='red',
#             label='2016',fliers_off=True)
#plt.plot(ar6['Latitude'][ar6['flacr']],ar6['AOD0501'][ar6['flacr']],'.',color='orange',alpha=0.5,
#         markersize=0.4)
make_boxplot(ar6['AOD0501'][ar6['flacr']],ar6['Latitude'][ar6['flacr']],lim2a,pos2a,color='red',
             label='2016 routine',fliers_off=True)

plt.xlabel('Latitude [$^\\circ$]')
plt.ylabel('4STAR AOD at 501 nm')
plt.grid()
plt.xlim([-22,0])
plt.ylim(0,1)
plt.legend(frameon=False,numpoints=1,loc=2)
prelim()

plt.title('ORACLES above cloud AOD - oceanic')
plt.savefig(fp+'plot_all/AOD_oceanc_latitude_2016_2017_2018.png',dpi=600,transparent=True)


# # Plot some angstrom exponents

# ## Setup the vertical angtrom profiles

# In[26]:


ar6['fl_QA_angs'] = ar6['fl_QA'] & (ar6['AOD0501']>0.1)
ar6['fl_QA_angs_aca'] = ar6['flac'] & (ar6['AOD0501']>0.1) & (ar6['GPS_Alt']>300.0)


# In[27]:


ar7['fl_QA_angs'] = ar7['fl_QA'] & (ar7['AOD0501']>0.1)
ar7['fl_QA_angs_aca'] = ar7['flac'] & (ar7['AOD0501']>0.1) & (ar7['GPS_Alt']>300.0)


# In[28]:


ar8['fl_QA_angs'] = ar8['fl_QA'] & (ar8['AOD0501']>0.1)
ar8['fl_QA_angs_aca'] = ar8['flac'] & (ar8['AOD0501']>0.1) & (ar8['GPS_Alt']>300.0)


# In[306]:


def make_bined_alt(x,alt,days,fl,n=70):
    binned_ang,binned_alt,binned_num,binned_ndays = [],[],[],[]
    for i in xrange(70):
        flaa = (alt[fl]>=i*100.0) & (alt[fl]<(i+1.0)*100.0)
        binned_ang.append(x[fl][flaa])
        binned_alt.append(np.mean([i*100.0,(i+1.0)*100.0]))
        binned_num.append(len(x[fl][flaa]))
        binned_ndays.append(len(np.unique(days[fl][flaa])))
    return binned_ang,binned_alt,binned_num,binned_ndays


# In[307]:


def set_box_whisker_color(cl,bp,binned_ndays,binmax=None,a=0.6):
    bndm = np.nanmax(binned_ndays)*1.0
    if binmax: bndm = binmax
    for j,b in enumerate(bp['boxes']):
        try:
            b.set_facecolor(cl(binned_ndays[j]*1.0/bndm))
            b.set_edgecolor(cl(binned_ndays[j]*1.0/bndm))
        except:
            b.set_color(cl(binned_ndays[j]*1.0/bndm))
            b.set_fillstyle('full')
        b.set_alpha(a)
    for j,b in enumerate(bp['means']):
        b.set_marker('.')
        b.set_color('None')
        b.set_markerfacecolor('darkgreen')
        b.set_markeredgecolor('darkgreen')
        b.set_alpha(a)
    for j,b in enumerate(bp['whiskers']):
        b.set_linestyle('-')
        b.set_color('pink') #gr(binned_ndays[j]*1.0/bndm))
        b.set_alpha(a+0.1)
    for j,b in enumerate(bp['caps']):
        b.set_alpha(a+0.1)
        b.set_color('pink')#gr(binned_ndays[j]*1.0/bndm))
    for j,b in enumerate( bp['medians']):
        b.set_linewidth(4)
        b.set_color('gold')
        b.set_alpha(a-0.2)
    
    return


# In[31]:


binned_ang6,binned_alt6,binned_num6,binned_ndays6 = make_bined_alt(ar6['AOD_angstrom_470_865'],
                                                                   ar6['GPS_Alt'],ar6['days'],ar6['fl_QA_angs'])
binned_ang7,binned_alt7,binned_num7,binned_ndays7 = make_bined_alt(ar7['AOD_angstrom_470_865'],
                                                                   ar7['GPS_Alt'],ar7['days'],ar7['fl_QA_angs'])
binned_ang8,binned_alt8,binned_num8,binned_ndays8 = make_bined_alt(ar8['AOD_angstrom_470_865'],
                                                                   ar8['GPS_Alt'],ar8['days'],ar8['fl_QA_angs'])
binned_angc6,binned_altc6,binned_numc6,binned_ndaysc6 = make_bined_alt(ar6['AOD_angstrom_470_865'],
                                                                   ar6['GPS_Alt'],ar6['days'],ar6['fl_QA_angs_aca'])
binned_angc7,binned_altc7,binned_numc7,binned_ndaysc7 = make_bined_alt(ar7['AOD_angstrom_470_865'],
                                                                   ar7['GPS_Alt'],ar7['days'],ar7['fl_QA_angs_aca'])
binned_angc8,binned_altc8,binned_numc8,binned_ndaysc8 = make_bined_alt(ar8['AOD_angstrom_470_865'],
                                                                   ar8['GPS_Alt'],ar8['days'],ar8['fl_QA_angs_aca'])


# ## Plot some profiles

# In[34]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(binned_ang6,positions=np.array(binned_alt6)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom of 470 nm / 865 nm')
plt.ylabel('Altitude [m]')
gr = plt.cm.RdPu
bl = plt.cm.Blues
set_box_whisker_color(gr,bp,binned_ndays6)
    
bpc =plt.boxplot(binned_angc6,positions=np.array(binned_altc6)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bpc,binned_ndaysc6)
bpc['boxes'][0].set_color('grey')

ax = plt.gca()
plt.title('2016')
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,2.2)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['All data','Above cloud','Mean','Median','25% - 75%','min-max'],
           frameon=False,loc=2,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(binned_ndays6)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(binned_ndays6)
cbaxesgr = plt.gcf().add_axes([0.33, 0.2, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.35, 0.2, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLES2016_4STAR_Angstrom_2wvl_vertical_cb.png',
            transparent=True,dpi=500)


# In[35]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(binned_ang7,positions=np.array(binned_alt7)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom of 470 nm / 865 nm')
plt.ylabel('Altitude [m]')
gr = plt.cm.RdPu
bl = plt.cm.Blues
set_box_whisker_color(gr,bp,binned_ndays7)
    
bpc =plt.boxplot(binned_angc7,positions=np.array(binned_altc7)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bpc,binned_ndaysc7)
bpc['boxes'][0].set_color('grey')

ax = plt.gca()
plt.title('2017')
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,2.2)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['All data','Above cloud','Mean','Median','25% - 75%','min-max'],
           frameon=False,loc=2,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(binned_ndays7)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(binned_ndays7)
cbaxesgr = plt.gcf().add_axes([0.33, 0.2, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.35, 0.2, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLES2017_4STAR_Angstrom_2wvl_vertical_cb.png',
            transparent=True,dpi=500)


# In[37]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(binned_ang8,positions=np.array(binned_alt8)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom of 470 nm / 865 nm')
plt.ylabel('Altitude [m]')
gr = plt.cm.RdPu
bl = plt.cm.Blues
set_box_whisker_color(gr,bp,binned_ndays8)
    
bpc =plt.boxplot(binned_angc8,positions=np.array(binned_altc8)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bpc,binned_ndaysc8)
bpc['boxes'][0].set_color('grey')

ax = plt.gca()
plt.title('2018')
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,2.2)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][10],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['All data','Above cloud','Mean','Median','25% - 75%','min-max'],
           frameon=False,loc=2,numpoints=1)
#prelim()
scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(binned_ndays8)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(binned_ndays8)
cbaxesgr = plt.gcf().add_axes([0.33, 0.2, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.35, 0.2, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLES2018_4STAR_Angstrom_2wvl_vertical_cb.png',
            transparent=True,dpi=500)


# # Plot some AOD profiles

# ## Prepare the AOD profiles

# In[32]:


binned_aod6,binned_alta6,binned_numa6,binned_ndaysa6 = make_bined_alt(ar6['AOD0501'],
                                                                   ar6['GPS_Alt'],ar6['days'],ar6['flo'])
binned_aod7,binned_alta7,binned_numa7,binned_ndaysa7 = make_bined_alt(ar7['AOD0501'],
                                                                   ar7['GPS_Alt'],ar7['days'],ar7['flo'])
binned_aod8,binned_alta8,binned_numa8,binned_ndaysa8 = make_bined_alt(ar8['AOD0501'],
                                                                   ar8['GPS_Alt'],ar8['days'],ar8['flo'])
binned_aodr6,binned_altar6,binned_numar6,binned_ndaysar6 = make_bined_alt(ar6['AOD0501'],
                                                                   ar6['GPS_Alt'],ar6['days'],ar6['flr'])
binned_aodr7,binned_altar7,binned_numar7,binned_ndaysar7 = make_bined_alt(ar7['AOD0501'],
                                                                   ar7['GPS_Alt'],ar7['days'],ar7['flr'])
binned_aodr8,binned_altar8,binned_numar8,binned_ndaysar8 = make_bined_alt(ar8['AOD0501'],
                                                                   ar8['GPS_Alt'],ar8['days'],ar8['flr'])


# ## Prepare profiles for different wavelengths

# In[20]:


aodnames = np.sort([k for k in ar6.keys() if k.startswith('AOD') & ('_' not in k)])


# In[21]:


aodnames


# In[22]:


aod_wvl = [float(a[3:]) for a in aodnames]


# In[23]:


aod_wvl


# In[27]:


bina6 = {}
binz6 = {}
binn6 = {}
bind6 = {}
bina6r = {}
binz6r = {}
binn6r = {}
bind6r = {}
for i,a in enumerate(aodnames):
    bina6[a],binz6[a],binn6[a],bind6[a] = make_bined_alt(ar6[a],ar6['GPS_Alt'],ar6['days'],ar6['flo'])
    bina6r[a],binz6r[a],binn6r[a],bind6r[a] = make_bined_alt(ar6[a],ar6['GPS_Alt'],ar6['days'],ar6['flr'])


# ## Plot the profiles

# In[39]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(binned_aod6,positions=np.array(binned_alta6)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('AOD_{{501}}')
plt.ylabel('Altitude [m]')
gr = plt.cm.YlGn
bl = plt.cm.YlOrBr
set_box_whisker_color(gr,bp,binned_ndaysa6)
    
bpc =plt.boxplot(binned_aodr6,positions=np.array(binned_altar6)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bpc,binned_ndaysar6)
bpc['boxes'][0].set_color('grey')

ax = plt.gca()
plt.title('2016')
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,0.8)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['Coastal (Other)','Oceanic (Routine)','Mean','Median','25% - 75%','min-max'],
           frameon=False,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(binned_ndaysa6)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(binned_ndaysa6)
cbaxesgr = plt.gcf().add_axes([0.8, 0.3, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.82, 0.3, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

#plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLES2016_4STAR_AOD_vertical_cb.png',
            transparent=True,dpi=500)


# In[40]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(binned_aod7,positions=np.array(binned_alta7)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('AOD_{{501}}')
plt.ylabel('Altitude [m]')
gr = plt.cm.YlGn
bl = plt.cm.YlOrBr
set_box_whisker_color(gr,bp,binned_ndaysa7)
    
bpc =plt.boxplot(binned_aodr7,positions=np.array(binned_altar7)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bpc,binned_ndaysar7)
bpc['boxes'][0].set_color('grey')

ax = plt.gca()
plt.title('2017')
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,0.8)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['Oceanic (Other)','Coastal (Routine)','Mean','Median','25% - 75%','min-max'],
           frameon=False,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(binned_ndaysa7)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(binned_ndaysa7)
cbaxesgr = plt.gcf().add_axes([0.8, 0.3, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.82, 0.3, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

#plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLES2017_4STAR_AOD_vertical_cb.png',
            transparent=True,dpi=500)


# In[98]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(binned_aod8,positions=np.array(binned_alta8)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('AOD_{{501}}')
plt.ylabel('Altitude [m]')
gr = plt.cm.YlGn
bl = plt.cm.YlOrBr
set_box_whisker_color(gr,bp,binned_ndaysa8)
    
bpc =plt.boxplot(binned_aodr8,positions=np.array(binned_altar8)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bpc,binned_ndaysar8)
bpc['boxes'][0].set_color('grey')

ax = plt.gca()
plt.title('2018')
plt.ylim(0,7000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
plt.xlim(0,0.8)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['Coastal (Other)','Oceanic (Routine)','Mean','Median','25% - 75%','min-max'],
           frameon=False,numpoints=1)
#prelim()
scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(binned_ndaysa8)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(binned_ndaysa8)
cbaxesgr = plt.gcf().add_axes([0.8, 0.3, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.82, 0.3, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

#plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLES2018_4STAR_AOD_vertical_cb.png',
            transparent=True,dpi=500)


# ## plot the profiles for 2016 per wavelength

# In[116]:


def extract_from_bp(bp):
    'Function to extract the means, medians, percentile, and range from boxplots'
    d = {}
    d['means'] = np.array([m.get_xdata()[0] for m in bp['means']])
    d['medians'] = np.array([m.get_xdata()[0] for m in bp['medians']])
    d['Q1_percentile'] = np.array([m.get_xdata()[0] for m in bp['boxes']])
    d['Q3_percentile'] = np.array([m.get_xdata()[2] for m in bp['boxes']])
    d['min'] = np.array([m.get_xdata()[0] for i,m in enumerate(bp['whiskers']) if not i%2])
    d['max'] = np.array([m.get_xdata()[1] for i,m in enumerate(bp['whiskers']) if i%2])
    return d


# In[130]:


d = {}


# In[131]:


for i,a in enumerate(aodnames):
    d[a] = {}
    plt.figure(figsize=(4,6))
    bp = plt.boxplot(bina6[a],positions=np.array(binz6[a])-5.0,vert=False,
                    showfliers=False,widths=90,showmeans=True,patch_artist=False)
    d[a]['coastal'] = extract_from_bp(bp)
    [[m.remove() for m in bp[k]] for k in bp]
    bp =plt.boxplot(bina6[a],positions=np.array(binz6[a])-5.0,vert=False,
                    showfliers=False,widths=90,showmeans=True,patch_artist=True)
    plt.xlabel('{}'.format(a))
    plt.ylabel('Altitude [m]')
    gr = plt.cm.YlGn
    bl = plt.cm.YlOrBr
    set_box_whisker_color(gr,bp,binn6[a],binmax=1800)

    bpc =plt.boxplot(bina6r[a],positions=np.array(binz6r[a])+10.0,vert=False,
                     showfliers=False,widths=90,showmeans=True,patch_artist=False)
    d[a]['oceanic'] = extract_from_bp(bpc)
    [[m.remove() for m in bpc[k]] for k in bpc]
    bpc =plt.boxplot(bina6r[a],positions=np.array(binz6r[a])+10.0,vert=False,
                     showfliers=False,widths=90,showmeans=True,patch_artist=True)
    set_box_whisker_color(bl,bpc,binn6r[a],binmax=1800)
    bpc['boxes'][0].set_color('grey')

    ax = plt.gca()
    plt.title('2016')
    plt.ylim(0,7000)
    plt.yticks([0,1000,2000,3000,4000,5000,6000,7000])
    ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000])
    plt.xlim(0,0.8)
    plt.grid()
    plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
               ['Coastal (Other)','Oceanic (Routine)','Mean','Median','25\% - 75\%','min-max'],
               frameon=False,numpoints=1)

    scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
    scalarmapgr.set_array([0,1800])
    scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
    scalarmapbl.set_array([0,1800])
    cbaxesgr = plt.gcf().add_axes([0.8, 0.3, 0.015, 0.3])
    cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr,extend='max')
    cbaxesbl = plt.gcf().add_axes([0.82, 0.3, 0.015, 0.3])
    cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl,extend='max')
    cbg.set_ticks([0,300,600,900,1200,1500])
    cbb.set_ticks([0,300,600,900,1200,1500]),cbb.set_ticklabels(['','','','',''])
    cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
    cbaxesgr.text(-9.0,0.5,'Number of samples',rotation=90,verticalalignment='center')

    #plt.tight_layout()

    plt.savefig(fp+'plot_all/ORACLES2016_4STAR_AOD_vertical_cb_{}.png'.format(a),
            transparent=True,dpi=500)


# ## Save to file the binned vertical AODs

# In[135]:


d['AOD0530']['coastal']['means']


# In[133]:


binz6[a]


# In[139]:


d[a]['coastal'].keys()


# In[140]:


binn6_rg = {'coastal':binn6,'oceanic':binn6r}


# In[149]:


from datetime import datetime


# In[152]:


for i,a in enumerate(aodnames):
    for reg in ['coastal','oceanic']:
        with open(fp+'ORACLES2016_profiles_{}_{}.dat'.format(a,reg),'w') as f:
            f.write('# 4STAR AOD profiles aglomerrated from all valid ORACLES 2016 measurements (R3) for wavelength: {:4.1f} nm in region: {}, by Samuel LeBlanc on {}\n'.format(aod_wvl[i],reg,datetime.now()))
            f.write('#Alt[m],mean,median,25thpercentile,75thpercentile,min,max,number[#]\n')
            for j,z in enumerate(binz6[a]):
                f.write('{:4.1f},{:2.3f},{:2.3f},{:2.3f},{:2.3f},{:2.3f},{:2.3f},{}\n'.format(z,
                        d[a][reg]['means'][j],d[a][reg]['medians'][j],d[a][reg]['Q1_percentile'][j],d[a][reg]['Q3_percentile'][j],
                        d[a][reg]['min'][j],d[a][reg]['max'][j],binn6_rg[reg][a][j]))


# # Plot a 3 year time series for a region

# ## Set the region and times

# In[359]:


lat1,lat2 = -17.0,-10.0
lon1,lon2 = 3.5,6.75


# In[360]:


ar6['flq'] = ar6['flac'] & (ar6['Latitude']>lat1) & (ar6['Latitude']<lat2) & (ar6['Longitude']>lon1) & (ar6['Longitude']<lon2) & (ar6['qual_flag']==0)& (ar6['AOD0501']<1.5)
ar7['flq'] = ar7['flac'] & (ar7['Latitude']>lat1) & (ar7['Latitude']<lat2) & (ar7['Longitude']>lon1) & (ar7['Longitude']<lon2) & (ar7['qual_flag']==0)& (ar7['AOD0501']<1.5)
ar8['flq'] = ar8['flac'] & (ar8['Latitude']>lat1) & (ar8['Latitude']<lat2) & (ar8['Longitude']>lon1) & (ar8['Longitude']<lon2) & (ar8['qual_flag']==0)& (ar8['AOD0501']<1.5)


# In[361]:


days6 = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']
days7 = ['20170801','20170802','20170807','20170809', '20170812','20170813','20170815','20170817','20170818','20170819','20170821',
        '20170824','20170826','20170828','20170830','20170831','20170902','20170903','20170904']
days8 = ['20180921','20180922','20180924','20180927','20180930','20181002','20181003','20181005','20181007','20181010','20181012',
        '20181015','20181017','20181019','20181021','20181023','20181025','20181026','20181027']


# In[362]:


ar6['daysd'] = [days6[i] for i in ar6['days'].astype(int)]
ar7['daysd'] = [days7[i] for i in ar7['days'].astype(int)]
ar8['daysd'] = [days8[i] for i in ar8['days'].astype(int)]


# In[363]:


ar6['ndtime'] = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8]),int(ar6['Start_UTC'][i]),
                          int((ar6['Start_UTC'][i]-float(int(ar6['Start_UTC'][i])))*60)) for i,d in enumerate(ar6['daysd'])]
ar7['ndtime'] = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8]),int(ar7['Start_UTC'][i]),
                          int((ar7['Start_UTC'][i]-float(int(ar7['Start_UTC'][i])))*60)) for i,d in enumerate(ar7['daysd'])]
ar8['ndtime'] = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8]),int(ar8['Start_UTC'][i]),
                          int((ar8['Start_UTC'][i]-float(int(ar8['Start_UTC'][i])))*60)) for i,d in enumerate(ar8['daysd'])]


# In[364]:


ar6['ndtimes'] = np.array(ar6['ndtime'])
ar7['ndtimes'] = np.array(ar7['ndtime'])
ar8['ndtimes'] = np.array(ar8['ndtime'])


# In[365]:


ar6['ndtime2'] = np.array([datetime(2018,int(d[4:6]),int(d[6:8]),int(ar6['Start_UTC'][i]),
                          int((ar6['Start_UTC'][i]-float(int(ar6['Start_UTC'][i])))*60)) for i,d in enumerate(ar6['daysd'])])
ar7['ndtime2'] = np.array([datetime(2018,int(d[4:6]),int(d[6:8]),int(ar7['Start_UTC'][i]),
                          int((ar7['Start_UTC'][i]-float(int(ar7['Start_UTC'][i])))*60)) for i,d in enumerate(ar7['daysd'])])
ar8['ndtime2'] = np.array([datetime(2018,int(d[4:6]),int(d[6:8]),int(ar8['Start_UTC'][i]),
                          int((ar8['Start_UTC'][i]-float(int(ar8['Start_UTC'][i])))*60)) for i,d in enumerate(ar8['daysd'])])


# ## Plot the time series of the region

# In[40]:


bin_aod6,bin_doy6,bin_num6 = [],[],[]
bin_days6 = np.unique(ar6['days'][ar6['flq']])
for d in bin_days6:
    flh = (ar6['days'][ar6['flq']]==d)
    bin_doy6.append(ar6['ndtime2'][ar6['flq']][flh][0].timetuple().tm_yday)
    bin_aod6.append(ar6['AOD0501'][ar6['flq']][flh])
    bin_num6.append(len(ar6['AOD0501'][ar6['flq']][flh]))


# In[41]:


bin_aod7,bin_doy7,bin_num7 = [],[],[]
bin_days7 = np.unique(ar7['days'][ar7['flq']])
for d in bin_days7:
    flh = (ar7['days'][ar7['flq']]==d)
    bin_doy7.append(ar7['ndtime2'][ar7['flq']][flh][0].timetuple().tm_yday)
    bin_aod7.append(ar7['AOD0501'][ar7['flq']][flh])
    bin_num7.append(len(ar7['AOD0501'][ar7['flq']][flh]))


# In[42]:


bin_aod8,bin_doy8,bin_num8 = [],[],[]
bin_days8 = np.unique(ar8['days'][ar8['flq']])
for d in bin_days8:
    flh = (ar8['days'][ar8['flq']]==d)
    bin_doy8.append(ar8['ndtime2'][ar8['flq']][flh][0].timetuple().tm_yday)
    bin_aod8.append(ar8['AOD0501'][ar8['flq']][flh])
    bin_num8.append(len(ar8['AOD0501'][ar8['flq']][flh]))


# In[52]:


plt.figure()
bp = plt.boxplot(bin_aod6,positions=bin_doy6,vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
bl = plt.cm.Blues
set_box_whisker_color(bl,bp,np.array(bin_num6))

bp7 = plt.boxplot(np.array(bin_aod7),positions=np.array(bin_doy7),vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
gr = plt.cm.YlOrBr
set_box_whisker_color(gr,bp7,np.array(bin_num7))

bp8 = plt.boxplot(np.array(bin_aod8),positions=np.array(bin_doy8),vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
by = plt.cm.Greens
set_box_whisker_color(by,bp8,np.array(bin_num8))

xy6 = [[r.get_data()[0],r.get_data()[1]] for r in bp['means']]
xy6 = np.array(xy6)[:,:,0]
xy7 = [[r.get_data()[0],r.get_data()[1]] for r in bp7['means']]
xy7 = np.array(xy7)[:,:,0]
xy8 = [[r.get_data()[0],r.get_data()[1]] for r in bp8['means']]
xy8 = np.array(xy8)[:,:,0]
xx = np.append(np.append(xy7[:,0],xy6[:,0]),xy8[:,0])
yy = np.append(np.append(xy7[:,1],xy6[:,1]),xy8[:,1])
xn = np.linspace(xx.min(),xx.max(),50)
spl = UnivariateSpline(xx,yy, k=5)
lb = plt.plot(xn,spl(xn),'-k',label='fit')

plt.xlim(210,310)
plt.xticks([220,230,240,250,260,270,280,290,300])
plt.gca().set_xticklabels([220,230,240,250,260,270,280,290,300])

plt.ylabel('AOD_{{501}}')
plt.xlabel('Day of Year')
plt.title('ORACLES ACAOD between 17$^{{\circ}}$S - 10$^{{\circ}}$S and 3.5$^{{\circ}}$W - 6.75$^{{\circ}}$W')

#plt.grid()
bp['boxes'][0].set_color('grey')
plt.legend([bp['boxes'][1],bp7['boxes'][0],bp8['boxes'][0],bp['means'][0],bp['medians'][0],bp['boxes'][0],
            bp['whiskers'][0],lb[0]],
           ['2016','2017','2018','Mean','Median','25$\%$ - 75$\%$','min-max','fit'],frameon=False,numpoints=1)
plt.savefig(fp+'plot_all/ORACLESall_4STAR_AOD_monthly_hist.png',
            transparent=True,dpi=500)


# In[200]:


plt.figure()
plt.plot(ar6['ndtime2'][ar6['flq']],ar6['AOD0501'][ar6['flq']],'rs',label='2016',markeredgecolor='None')
plt.plot(ar7['ndtime2'][ar7['flq']],ar7['AOD0501'][ar7['flq']],'g^',label='2017',markeredgecolor='None')
plt.plot(ar8['ndtime2'][ar8['flq']],ar8['AOD0501'][ar8['flq']],'bo',label='2018',markeredgecolor='None')
plt.title('ORACLES ACAOD')
monthyearFmt = mdates.DateFormatter('%B')
monthyearFmt2 = mdates.DateFormatter('%m/%d')
plt.gca().xaxis.set_major_formatter(monthyearFmt2)
#plt.gca().xaxis.set_minor_formatter(monthyearFmt2)
plt.xticks(rotation=90)
ax.set_xlim([datetime(2018, 7, 30), datetime(2018, 10, 30)])
plt.legend(frameon=False,numpoints=1,handletextpad=0.2)
plt.ylabel('AOD_{{501}}')
plt.savefig(fp+'plot_all/ORACLESall_4STAR_AOD_monthly.png',
            transparent=True,dpi=500)


# ## Add a map of the region

# In[59]:


from mpl_toolkits.basemap import Basemap


# In[60]:


def mapfig(ax=plt.gca()):
    m = Basemap(projection='merc',llcrnrlat=-25,urcrnrlat=8,llcrnrlon=-15,urcrnrlon=18,resolution='l',ax=ax)
    m.drawcoastlines()
    #m.drawmeridians(np.linspace(-17,11,8),labels=[0,0,0,1],linewidth=0.1)
    #m.drawparallels(np.linspace(-18,4,12),labels=[1,0,0,0],linewidth=0.1)
    m.drawlsmask(land_color='lightgrey',ocean_color='None',lakes=True)
    #m.shadedrelief(alpha=0.4)
    return m


# In[61]:


fig,ax = plt.subplots(1,1,figsize=(5,5))
ax1 = ax
#ax1 = plt.subplot(1,2,1)
m = mapfig(ax=ax1)
x6,y6 = m(ar6['Longitude'],ar6['Latitude'])
m.plot(x6,y6,'.',color='lightcoral',alpha=0.002,markeredgecolor='None')

x7,y7 = m(ar7['Longitude'],ar7['Latitude'])
m.plot(x7,y7,'.',color='lightgreen',alpha=0.002,markeredgecolor='None')

x8,y8 = m(ar8['Longitude'],ar8['Latitude'])
m.plot(x8,y8,'.',color='lightskyblue',alpha=0.002,markeredgecolor='None')

xss,yss = m([lon1,lon1,lon2,lon2,lon1],[lat1,lat2,lat2,lat1,lat1])
xss2,yss2 = m([lon1,lon1,lon2,lon2,lon1],[lat1,lat1,lat1,lat1,lat1])

plt.fill_between(xss, yss, yss2,color='gold',alpha=0.7)
#plt.Polygon(xss,yss,edgecolor='None',color='gold',alpha=0.3)


# In[63]:


from matplotlib import rc
rc("text", usetex=False)


# In[64]:


plt.figure(figsize=(9,3))
bp = plt.boxplot(bin_aod6,positions=bin_doy6,vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
bl = plt.cm.bwr
set_box_whisker_color(bl,bp,np.array(bin_num6)*0.0+3130)

bp7 = plt.boxplot(np.array(bin_aod7),positions=np.array(bin_doy7),vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
gr = plt.cm.brg
set_box_whisker_color(gr,bp7,np.array(bin_num7)*0.0+870)

bp8 = plt.boxplot(np.array(bin_aod8),positions=np.array(bin_doy8),vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
by = plt.cm.bwr_r
set_box_whisker_color(by,bp8,np.array(bin_num8)*0.0+760)

xy6 = [[r.get_data()[0],r.get_data()[1]] for r in bp['means']]
xy6 = np.array(xy6)[:,:,0]
xy7 = [[r.get_data()[0],r.get_data()[1]] for r in bp7['means']]
xy7 = np.array(xy7)[:,:,0]
xy8 = [[r.get_data()[0],r.get_data()[1]] for r in bp8['means']]
xy8 = np.array(xy8)[:,:,0]
xx = np.append(np.append(xy7[:,0],xy6[:,0]),xy8[:,0])
yy = np.append(np.append(xy7[:,1],xy6[:,1]),xy8[:,1])
xn = np.linspace(xx.min(),xx.max(),50)
spl = UnivariateSpline(xx,yy, k=5)
lb = plt.plot(xn,spl(xn),'-k',label='fit')

plt.xlim(213,320)
plt.ylim(0.1,0.9)
plt.xticks([220,230,240,250,260,270,280,290,300])
plt.gca().set_xticklabels([220,230,240,250,260,270,280,290,300])

plt.plot([213,243],[0.1,0.1],'-',color='lightgreen',lw=6)
plt.plot([244,273],[0.1,0.1],'-',color='lightcoral',lw=6)
plt.plot([274,305],[0.1,0.1],'-',color='lightskyblue',lw=6)
plt.text(215,0.12,'August',color='g')
plt.text(246,0.12,'September',color='r')
plt.text(276,0.12,'October',color='b')

plt.ylabel('AOD_{{501}}')
plt.xlabel('Day of Year')
plt.title('ORACLES ACAOD between 17$^{{\circ}}$S - 10$^{{\circ}}$S and 3.5$^{{\circ}}$W - 6.75$^{{\circ}}$W')

#plt.grid()
bp['boxes'][0].set_color('grey')
plt.legend([bp['boxes'][1],bp7['boxes'][0],bp8['boxes'][0],bp['means'][0],bp['medians'][0],bp['boxes'][0],
            bp['whiskers'][0],lb[0]],
           ['2016','2017','2018','Mean','Median','25$\%$ - 75$\%$','min-max','fit'],frameon=False,numpoints=1)

axb = plt.gcf().add_axes([0.06, 0.6, 0.28, 0.28])
m = mapfig(ax=axb)
x6,y6 = m(ar6['Longitude'],ar6['Latitude'])
m.plot(x6,y6,'.',color='lightcoral',alpha=0.006,markeredgecolor='None')

x7,y7 = m(ar7['Longitude'],ar7['Latitude'])
m.plot(x7,y7,'.',color='lightgreen',alpha=0.006,markeredgecolor='None')

x8,y8 = m(ar8['Longitude'],ar8['Latitude'])
m.plot(x8,y8,'.',color='lightskyblue',alpha=0.05,markeredgecolor='None')

xss,yss = m([lon1,lon1,lon2,lon2,lon1],[lat1,lat2,lat2,lat1,lat1])
xss2,yss2 = m([lon1,lon1,lon2,lon2,lon1],[lat1,lat1,lat1,lat1,lat1])

plt.fill_between(xss, yss, yss2,color='k',alpha=0.6)

plt.savefig(fp+'plot_all/ORACLESall_4STAR_AOD_monthly_hist_map.png',
            transparent=True,dpi=500)


# In[81]:


plt.figure(figsize=(9,3))
bp = plt.boxplot(bin_aod6,positions=bin_doy6,vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
bl = plt.cm.bwr
set_box_whisker_color(bl,bp,np.array(bin_num6)*0.0+3130)

bp7 = plt.boxplot(np.array(bin_aod7),positions=np.array(bin_doy7),vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
gr = plt.cm.brg
set_box_whisker_color(gr,bp7,np.array(bin_num7)*0.0+870)

bp8 = plt.boxplot(np.array(bin_aod8),positions=np.array(bin_doy8),vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
by = plt.cm.bwr_r
set_box_whisker_color(by,bp8,np.array(bin_num8)*0.0+760)


hbp = plt.boxplot(hsrl['bin_aod6'],positions=hsrl['bin_doy6'],vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
hbl = plt.cm.cool
set_box_whisker_color(hbl,hbp,np.array(hsrl['bin_num6'])*0.0+3130,a=0.3)

hbp7 = plt.boxplot(np.array(hsrl['bin_aod7']),positions=np.array(hsrl['bin_doy7']),vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
hgr = plt.cm.summer_r
set_box_whisker_color(hgr,hbp7,np.array(hsrl['bin_num7'])*0.0+870,a=0.3)

hbp8 = plt.boxplot(np.array(hsrl['bin_aod8']),positions=np.array(hsrl['bin_doy8']),vert=True,
                showfliers=False,widths=2,showmeans=True,patch_artist=True)
hby = plt.cm.coolwarm_r
set_box_whisker_color(hby,hbp8,np.array(hsrl['bin_num8'])*0.0+760,a=0.3)


xy6 = [[r.get_data()[0],r.get_data()[1]] for r in bp['means']]
xy6 = np.array(xy6)[:,:,0]
xy7 = [[r.get_data()[0],r.get_data()[1]] for r in bp7['means']]
xy7 = np.array(xy7)[:,:,0]
xy8 = [[r.get_data()[0],r.get_data()[1]] for r in bp8['means']]
xy8 = np.array(xy8)[:,:,0]
xx = np.append(np.append(xy7[:,0],xy6[:,0]),xy8[:,0])
yy = np.append(np.append(xy7[:,1],xy6[:,1]),xy8[:,1])
xn = np.linspace(xx.min(),xx.max(),50)
spl = UnivariateSpline(xx,yy, k=5)
lb = plt.plot(xn,spl(xn),'-k',label='fit')

hxy6 = [[r.get_data()[0],r.get_data()[1]] for r in hbp['means']]
hxy6 = np.array(hxy6)[:,:,0]
hxy7 = [[r.get_data()[0],r.get_data()[1]] for r in hbp7['means']]
hxy7 = np.array(hxy7)[:,:,0]
hxy8 = [[r.get_data()[0],r.get_data()[1]] for r in hbp8['means']]
hxy8 = np.array(hxy8)[:,:,0]
hxx = np.append(np.append(hxy7[:,0],hxy6[:,0]),hxy8[:,0])
hyy = np.append(np.append(hxy7[:,1],hxy6[:,1]),hxy8[:,1])
hxn = np.linspace(hxx.min(),hxx.max(),50)
hspl = UnivariateSpline(hxx,hyy, k=5)
hlb = plt.plot(hxn,hspl(hxn),'--k',alpha=0.3)

plt.xlim(213,320)
plt.ylim(0.1,0.9)
plt.xticks([220,230,240,250,260,270,280,290,300])
plt.gca().set_xticklabels([220,230,240,250,260,270,280,290,300])

plt.plot([213,243],[0.1,0.1],'-',color='lightgreen',lw=6)
plt.plot([244,273],[0.1,0.1],'-',color='lightcoral',lw=6)
plt.plot([274,305],[0.1,0.1],'-',color='lightskyblue',lw=6)
plt.text(215,0.12,'August',color='g')
plt.text(246,0.12,'September',color='r')
plt.text(276,0.12,'October',color='b')

plt.ylabel('ACAOD')
plt.xlabel('Day of Year')
plt.title('ORACLES ACAOD between 17$^{{\circ}}$S - 10$^{{\circ}}$S and 3.5$^{{\circ}}$W - 6.75$^{{\circ}}$W')

#plt.grid()
bp['boxes'][0].set_color('grey')
plt.legend([bp['boxes'][1],bp7['boxes'][0],bp8['boxes'][0],hbp['boxes'][1],hbp7['boxes'][0],hbp8['boxes'][0]],
           ['2016 4STAR','2017 4STAR','2018 4STAR','2016 HSRL','2017 HSRL','2018 HSRL'],frameon=False,numpoints=1)
#plt.legend([bp['means'][0],bp['medians'][0],bp['boxes'][0],
#            bp['whiskers'][0],lb[0]],
#           ['Mean','Median','25$\%$ - 75$\%$','min-max','fit'],frameon=False,numpoints=1,loc='upper center')

axb = plt.gcf().add_axes([0.06, 0.6, 0.28, 0.28])
m = mapfig(ax=axb)
x6,y6 = m(ar6['Longitude'],ar6['Latitude'])
m.plot(x6,y6,'.',color='lightcoral',alpha=0.006,markeredgecolor='None')

x7,y7 = m(ar7['Longitude'],ar7['Latitude'])
m.plot(x7,y7,'.',color='lightgreen',alpha=0.006,markeredgecolor='None')

x8,y8 = m(ar8['Longitude'],ar8['Latitude'])
m.plot(x8,y8,'.',color='lightskyblue',alpha=0.05,markeredgecolor='None')

xss,yss = m([lon1,lon1,lon2,lon2,lon1],[lat1,lat2,lat2,lat1,lat1])
xss2,yss2 = m([lon1,lon1,lon2,lon2,lon1],[lat1,lat1,lat1,lat1,lat1])

plt.fill_between(xss, yss, yss2,color='k',alpha=0.6)

plt.savefig(fp+'plot_all/ORACLESall_4STAR_AOD_monthly_hist_map_with_HSRL.png',
            transparent=True,dpi=500)


# # Plot wavelength AOD spectra

# In[35]:


kar6 = ar6.keys()
kar6.sort()


# In[36]:


nn = [i for i in kar6  if i[0]=='A']


# In[37]:


nn[0:24]


# In[38]:


nn.pop(13)


# In[39]:


ar6['aods'] = np.array([ar6[i] for i in nn[0:23]]).T
ar7['aods'] = np.array([ar7[i] for i in nn[0:23]]).T
ar8['aods'] = np.array([ar8[i] for i in nn[0:23]]).T


# In[40]:


ar6['meanaod'] = np.nanmean(ar6['aods'][ar6['flac'],:],axis=0)
ar6['medianaod'] = np.nanmedian(ar6['aods'][ar6['flac'],:],axis=0)
ar6['stdaod'] = np.nanstd(ar6['aods'][ar6['flac'],:],axis=0)


# In[41]:


ar7['meanaod'] = np.nanmean(ar7['aods'][ar7['flac'],:],axis=0)
ar7['medianaod'] = np.nanmedian(ar7['aods'][ar7['flac'],:],axis=0)
ar7['stdaod'] = np.nanstd(ar7['aods'][ar7['flac'],:],axis=0)


# In[42]:


ar8['meanaod'] = np.nanmean(ar8['aods'][ar8['flac'],:],axis=0)
ar8['medianaod'] = np.nanmedian(ar8['aods'][ar8['flac'],:],axis=0)
ar8['stdaod'] = np.nanstd(ar8['aods'][ar8['flac'],:],axis=0)


# In[43]:


ar6['meanAE'] = np.nanmean(ar6['AOD_angstrom_470_865'][ar6['flac']])
ar7['meanAE'] = np.nanmean(ar7['AOD_angstrom_470_865'][ar7['flac']])
ar8['meanAE'] = np.nanmean(ar8['AOD_angstrom_470_865'][ar8['flac']])
ar6['stdAE'] = np.nanstd(ar6['AOD_angstrom_470_865'][ar6['flac']])
ar7['stdAE'] = np.nanstd(ar7['AOD_angstrom_470_865'][ar7['flac']])
ar8['stdAE'] = np.nanstd(ar8['AOD_angstrom_470_865'][ar8['flac']])


# In[44]:


wvls = [355.0,380.0,452.0,470.0,501.0,520.0,530.0,532.0,550.0,606.0,620.0,660.0,675.0,
        781.0,865.0,1020.0,1040.0,1064.0,1236.0,1250.0,1559.0,1627.0,1650.0,]


# In[66]:


plt.figure()
plt.plot(wvls,ar6['meanaod'],'-xr',label='2016, AE:{:1.2f}$\pm${:1.2f}'.format(ar6['meanAE'],ar6['stdAE']),lw=3)
plt.plot(wvls,ar7['meanaod'],'-xg',label='2017, AE:{:1.2f}$\pm${:1.2f}'.format(ar7['meanAE'],ar7['stdAE']),lw=3)
plt.plot(wvls,ar8['meanaod'],'-xb',label='2018, AE:{:1.2f}$\pm${:1.2f}'.format(ar8['meanAE'],ar8['stdAE']),lw=3)

plt.plot(wvls,ar6['meanaod']+ar6['stdaod'],'--r',alpha=0.5)
plt.plot(wvls,ar7['meanaod']+ar7['stdaod'],'--g',alpha=0.5)
plt.plot(wvls,ar8['meanaod']+ar8['stdaod'],'--b',alpha=0.5)

plt.plot(wvls,ar8['meanaod']-ar8['stdaod'],'--k',alpha=0.5,label='$\pm$ 1 std')

plt.plot(wvls,ar6['meanaod']-ar6['stdaod'],'--r',alpha=0.5)
plt.plot(wvls,ar7['meanaod']-ar7['stdaod'],'--g',alpha=0.5)
plt.plot(wvls,ar8['meanaod']-ar8['stdaod'],'--b',alpha=0.5)

plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.legend(frameon=False)
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.grid()
plt.xticks([350,400,500,600,700,800,900,1000,1200,1500,1650])
plt.gca().set_xticklabels([350,400,500,600,700,800,900,1000,1200,1500,1650])
plt.yticks([0.05,0.1,0.2,0.3,0.4,0.5,0.7])
plt.gca().set_yticklabels([0.05,0.1,0.2,0.3,0.4,0.5,0.7])
plt.xlim(350,1700)
plt.ylim(0.04,0.8)
#prelim()
plt.title('ORACLES ACAOD average spectra')
plt.savefig(fp+'plot_all/ORACLESall_4STAR_AOD_spectr_R1.png',
            transparent=True,dpi=500)


# # Plot mean AOD maps

# ## Setup the statistics

# In[45]:


ar6['bin_mean'],xe,ye,bn = st.binned_statistic_2d(ar6['Latitude'][ar6['flac']],ar6['Longitude'][ar6['flac']],
                                    ar6['AOD0501'][ar6['flac']],
                           bins=36,range=[[-25,2],[-16,16]])
ar6['bin_mean'] = np.ma.masked_array(ar6['bin_mean'],np.isnan(ar6['bin_mean']))


# In[46]:


ar7['bin_mean'],xe,ye,bn = st.binned_statistic_2d(ar7['Latitude'][ar7['flac']],ar7['Longitude'][ar7['flac']],
                                    ar7['AOD0501'][ar7['flac']],
                           bins=36,range=[[-25,2],[-16,16]])
ar7['bin_mean'] = np.ma.masked_array(ar7['bin_mean'],np.isnan(ar7['bin_mean']))


# In[47]:


ar8['bin_mean'],xe,ye,bn = st.binned_statistic_2d(ar8['Latitude'][ar8['flac']],ar8['Longitude'][ar8['flac']],
                                    ar8['AOD0501'][ar8['flac']],
                           bins=36,range=[[-25,2],[-16,16]])
ar8['bin_mean'] = np.ma.masked_array(ar8['bin_mean'],np.isnan(ar8['bin_mean']))


# In[48]:


uniq_cnt = lambda x: len(np.unique(x))


# In[49]:


ar6['bin_nday'],xed,yed,bn = st.binned_statistic_2d(ar6['Latitude'][ar6['flac']],ar6['Longitude'][ar6['flac']],
                                                    ar6['days'][ar6['flac']],
                           bins=36,range=[[-25,2],[-16,16]],statistic=uniq_cnt)
ar6['bin_nday'] = np.ma.masked_array(ar6['bin_nday'],np.isnan(ar6['bin_nday']))


# In[50]:


ar7['bin_nday'],xed,yed,bn = st.binned_statistic_2d(ar7['Latitude'][ar7['flac']],ar7['Longitude'][ar7['flac']],
                                                    ar7['days'][ar7['flac']],
                           bins=36,range=[[-25,2],[-16,16]],statistic=uniq_cnt)
ar7['bin_nday'] = np.ma.masked_array(ar7['bin_nday'],np.isnan(ar7['bin_nday']))


# In[51]:


ar8['bin_nday'],xed,yed,bn = st.binned_statistic_2d(ar8['Latitude'][ar8['flac']],ar8['Longitude'][ar8['flac']],
                                                    ar8['days'][ar8['flac']],
                           bins=36,range=[[-25,2],[-16,16]],statistic=uniq_cnt)
ar8['bin_nday'] = np.ma.masked_array(ar8['bin_nday'],np.isnan(ar8['bin_nday']))


# In[52]:


ia6 = np.where((ar6['bin_nday'].data>0.0))
ia7 = np.where((ar7['bin_nday'].data>0.0))
ia8 = np.where((ar8['bin_nday'].data>0.0))


# ## Plot on a map

# In[53]:


from mpl_toolkits.axes_grid1 import make_axes_locatable


# In[54]:


def mapfig(ax=plt.gca()):
    m = Basemap(projection='merc',llcrnrlat=-25,urcrnrlat=2,llcrnrlon=-16,urcrnrlon=16,resolution='l',ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.linspace(-17,11,8),labels=[0,0,0,1],linewidth=0.1)
    m.drawparallels(np.linspace(-26,4,11),labels=[1,0,0,0],linewidth=0.1)
    m.shadedrelief(alpha=0.4)
    return m


# In[55]:


np.linspace(-26,4,11)


# In[78]:


fig,ax = plt.subplots(1,3,figsize=(15,5))
ax = ax.flatten()

ax1 = ax[0]
m = mapfig(ax=ax1)
mx,my = m(ye,xe)
p = ax1.pcolor(mx,my,ar6['bin_mean'],vmin=0.0,vmax=0.8,cmap='plasma')
ax1.set_title('Sept. 2016')
#cb = m.colorbar(p,extend='both')

ax2 = ax[1]
m = mapfig(ax=ax2)
mx,my = m(ye,xe)
p = ax2.pcolor(mx,my,ar7['bin_mean'],vmin=0.0,vmax=0.8,cmap='plasma')
ax2.set_title('Aug. 2017')
#cb = m.colorbar(p,extend='both')

ax3 = ax[2]
m = mapfig(ax=ax3)
mx,my = m(ye,xe)
p = ax3.pcolor(mx,my,ar8['bin_mean'],vmin=0.0,vmax=0.8,cmap='plasma')
ax3.set_title('Oct. 2018')

#divider = make_axes_locatable(ax3)
#cax = divider.append_axes("right", size="5%", pad=0.05)
go = ax3.get_position()
bot = go.corners()[0,1]
lef = go.corners()[0,0]

cax = fig.add_axes([lef*1.08,bot+go.height*0.1,go.width*0.05,go.height*0.8])
cb = plt.colorbar(p,cax=cax,extend='max')
cb.set_label('Mean AOD$_{{500nm}}$')

plt.tight_layout(rect=[0.02,0,1,1],w_pad=3)


# In[79]:


go.corners()[0,0]


# In[80]:


fig,ax = plt.subplots(1,3,figsize=(15,5))
ax = ax.flatten()

ax1 = ax[0]
m = mapfig(ax=ax1)
mx,my = m(ye,xe)
p = ax1.scatter(mx[ia6[1]],my[ia6[0]],ar6['bin_nday'].data[ia6[0],ia6[1]]*35,
                c=ar6['bin_mean'].data[ia6[0],ia6[1]],
               marker='s',edgecolor='None',cmap='plasma',vmin=0.0,vmax=0.8)
#p = ax1.pcolor(mx,my,ar6['bin_mean'],vmin=0.0,vmax=0.8,cmap='plasma')
ax1.set_title('Sept. 2016')
#cb = m.colorbar(p,extend='both')

sizes = [1,2,3,5] #[10,100,500,1500]
labels = ['{0} day{1}'.format(z,'s' if z>1 else '') for z in sizes]
points = [ax1.scatter([], [], s=z*35, c='grey',marker='s',edgecolor='None') for z in sizes]
ax1.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left',handletextpad=0.1)

ax2 = ax[1]
m = mapfig(ax=ax2)
mx,my = m(ye,xe)
p = ax2.scatter(mx[ia7[1]],my[ia7[0]],ar7['bin_nday'].data[ia7[0],ia7[1]]*35,
                c=ar7['bin_mean'].data[ia7[0],ia7[1]],
               marker='s',edgecolor='None',cmap='plasma',vmin=0.0,vmax=0.8)
#p = ax2.pcolor(mx,my,ar7['bin_mean'],vmin=0.0,vmax=0.8,cmap='plasma')
ax2.set_title('Aug. 2017')

ax3 = ax[2]
m = mapfig(ax=ax3)
mx,my = m(ye,xe)
p = ax3.scatter(mx[ia8[1]],my[ia8[0]],ar8['bin_nday'].data[ia8[0],ia8[1]]*35,
                c=ar8['bin_mean'].data[ia8[0],ia8[1]],
               marker='s',edgecolor='None',cmap='plasma',vmin=0.0,vmax=0.8)
#p = ax3.pcolor(mx,my,ar8['bin_mean'],vmin=0.0,vmax=0.8,cmap='plasma')
ax3.set_title('Oct. 2018')
#prelim()
go = ax3.get_position()
bot = go.corners()[0,1]
lef = go.corners()[0,0]
cax = fig.add_axes([lef*1.08,bot+go.height*0.1,go.width*0.05,go.height*0.8])
cb = plt.colorbar(p,cax=cax,extend='max')
cb.set_label('Mean AOD$_{{500nm}}$')

plt.tight_layout(rect=[0.02,0,1,1],w_pad=3)


plt.savefig(fp+'plot_all/ORACLESall_4STAR_AOD_3map_stats_R1.png',
            transparent=True,dpi=500)


# ## plot a map of differences

# In[81]:


fig,ax = plt.subplots(1,2,figsize=(6,3))
ax = ax.flatten()

ax1 = ax[0]
m = mapfig(ax=ax1)
mx,my = m(ye,xe)
p = ax1.pcolor(mx,my,ar6['bin_mean']-ar7['bin_mean'],vmin=-0.4,vmax=0.4,cmap='seismic_r')
ax1.set_title('Sept. 2016 - Aug. 2017')

ax2 = ax[1]
m = mapfig(ax=ax2)
mx,my = m(ye,xe)
p = ax2.pcolor(mx,my,ar7['bin_mean']-ar8['bin_mean'],vmin=-0.4,vmax=0.4,cmap='seismic_r')
ax2.set_title('Aug. 2017 - Oct. 2018')

go = ax2.get_position()
bot = go.corners()[0,1]
lef = go.corners()[0,0]
cax = fig.add_axes([lef*1.08,bot+go.height*0.1,go.width*0.05,go.height*0.8])
cb = plt.colorbar(p,cax=cax,extend='both')
cb.set_label('Difference in AOD$_{{500nm}}$')

plt.tight_layout(rect=[0.03,0,1,1],w_pad=2.5)
plt.savefig(fp+'plot_all/ORACLESall_4STAR_AOD_3map_stats_difference_R1.png',
            transparent=True,dpi=500)


# # Get the ACAOD gaps

# ## Set up the function and variables

# In[349]:


def get_gap(index,alt,lat,lon,days,aod,ang):
    'Function to get the acaod gap extent (between the top of clouds and bottom of aerosol layer)'
    disc_flacaod_long = np.where(np.diff(index,1)>150)[0]
    
    discontinuity_istart_long =  index[np.append(0,disc_flacaod_long[:-1]+1)]
    discontinuity_iend_long =  index[disc_flacaod_long]
    
    ldelta_alt,ldelta_lon,ldelta_lat,ldelta_lon_days,ldelta_lat_days,ldelta_aod,ldelta_ang = [],[],[],[],[],[],[]
    aero_base_alt,aero_base_ind,meas_ind = [],[],[]
    
    for i,start in enumerate(discontinuity_istart_long):
        try:
            ma = np.nanmax(alt[start:discontinuity_iend_long[i]])
            mi = np.nanmin(alt[start:discontinuity_iend_long[i]])
            mai = np.nanargmax(alt[start:discontinuity_iend_long[i]])
            mii = np.nanargmin(alt[start:discontinuity_iend_long[i]])
            aero_base_ind.append(start+mai)
            meas_ind.append(start+mii)
            ldelta_alt.append(ma-mi)
            ldelta_lon.append(np.nanmean(lon[start:discontinuity_iend_long[i]]))
            ldelta_lat.append(np.nanmean(lat[start:discontinuity_iend_long[i]]))
            ldelta_aod.append(np.nanmean(aod[start:discontinuity_iend_long[i]]))
            ldelta_ang.append(np.nanmean(ang[start:discontinuity_iend_long[i]]))
            ldelta_lat_days.append(np.unique(days[start:discontinuity_iend_long[i]]))
            ldelta_lon_days.append(np.unique(days[start:discontinuity_iend_long[i]]))
            aero_base_alt.append(ma)
        except:
            print('Problem at i:{}'.format(i))
            pass
    ldelta_alt = np.array(ldelta_alt)
    ldelta_lon = np.array(ldelta_lon)
    ldelta_lat = np.array(ldelta_lat)
    ldelta_aod = np.array(ldelta_aod)
    ldelta_ang = np.array(ldelta_ang)
    ldelta_lat_days = np.array(ldelta_lat_days)
    ldelta_lon_days = np.array(ldelta_lon_days)
    aero_base_alt = np.array(aero_base_alt)
    aero_base_ind = np.array(aero_base_ind)
    meas_ind = np.array(meas_ind)
    
    d = {'dalt':ldelta_alt,'dlon':ldelta_lon,'dlat':ldelta_lat,'dlat_ndays':ldelta_lat_days,'dlon_ndays':ldelta_lon_days,
         'daod':ldelta_aod,'dang':ldelta_ang,'aero_base_alt':aero_base_alt,'aero_base_index':aero_base_ind,
         'meas_low_index':meas_ind}
    
    return d


# In[350]:


ar6['fl_acaod_noQA'] = ar6['flag_acaod']==1
ii_flacaod6 = np.where(ar6['fl_acaod_noQA'])[0] 
ii_flacaod6[0]


# In[351]:


ar7['fl_acaod_noQA'] = ar7['flag_acaod']==1
ii_flacaod7 = np.where(ar7['fl_acaod_noQA'])[0] 
ii_flacaod7[0]


# In[352]:


ar8['fl_acaod_noQA'] = ar8['flag_acaod']==1
ii_flacaod8 = np.where(ar8['fl_acaod_noQA'])[0] 
ii_flacaod8[0]


# In[353]:


gap6 = get_gap(ii_flacaod6,ar6['GPS_Alt'],ar6['Latitude'],ar6['Longitude'],ar6['days'],ar6['AOD0501'],ar6['AOD_angstrom_470_865'])
gap7 = get_gap(ii_flacaod7,ar7['GPS_Alt'],ar7['Latitude'],ar7['Longitude'],ar7['days'],ar7['AOD0501'],ar7['AOD_angstrom_470_865'])
gap8 = get_gap(ii_flacaod8,ar8['GPS_Alt'],ar8['Latitude'],ar8['Longitude'],ar8['days'],ar8['AOD0501'],ar8['AOD_angstrom_470_865'])


# In[354]:


gap6.keys()


# In[428]:


ar6['flo'].shape


# In[429]:


ar6['Start_UTC'].shape


# In[355]:


gap6['aero_base_UTC'] = ar6[u'Start_UTC'][gap6['aero_base_index']]
gap6['meas_low_UTC'] = ar6[u'Start_UTC'][gap6['meas_low_index']]


# In[366]:


gap6['aero_base_day'] = ar6['ndtimes'][gap6['aero_base_index']]
gap6['meas_low_day'] = ar6['ndtimes'][gap6['meas_low_index']]


# In[367]:


gap7['aero_base_UTC'] = ar7[u'Start_UTC'][gap7['aero_base_index']]
gap7['meas_low_UTC'] = ar7[u'Start_UTC'][gap7['meas_low_index']]
gap7['aero_base_day'] = ar7['ndtimes'][gap7['aero_base_index']]
gap7['meas_low_day'] = ar7['ndtimes'][gap7['meas_low_index']]


# In[368]:


gap8['aero_base_UTC'] = ar8[u'Start_UTC'][gap8['aero_base_index']]
gap8['meas_low_UTC'] = ar8[u'Start_UTC'][gap8['meas_low_index']]
gap8['aero_base_day'] = ar8['ndtimes'][gap8['aero_base_index']]
gap8['meas_low_day'] = ar8['ndtimes'][gap8['meas_low_index']]


# In[431]:


gap6['flo'] = ar6['flo'][gap6['aero_base_index']]
gap7['flo'] = ar7['flo'][gap7['aero_base_index']]
gap8['flo'] = ar8['flo'][gap8['aero_base_index']]
gap6['flr'] = ar6['flr'][gap6['aero_base_index']]
gap7['flr'] = ar7['flr'][gap7['aero_base_index']]
gap8['flr'] = ar8['flr'][gap8['aero_base_index']]


# ### Save the gap data

# In[37]:


len(gap6['aero_base_UTC']),len(gap6['aero_base_alt'])


# In[38]:


np.save(fp+'ORACLES2016_gap_{}.npy'.format(vv),gap6,allow_pickle=True)
np.save(fp+'ORACLES2017_gap_{}.npy'.format(vv),gap7,allow_pickle=True)
np.save(fp+'ORACLES2018_gap_{}.npy'.format(vv),gap8,allow_pickle=True)


# ### bin out the gaps

# In[433]:


def make_bined_x(x,alt,days,fl,bins=[]):
    binned_ang,binned_alt,binned_num,binned_ndays = [],[],[],[]
    bins.sort()
    for i,y in enumerate(bins[:-1]):
        flaa = (alt[fl]>=y) & (alt[fl]<bins[i+1])
        binned_alt.append(np.mean([y,bins[i+1]]))
        if any(flaa):
            binned_ang.append(x[fl][flaa])
            binned_num.append(len(x[fl][flaa]))
            binned_ndays.append(len(np.unique(np.hstack(days[fl][flaa]))))
        else:
            binned_ang.append(np.array(np.nan))
            binned_num.append(0.0)
            binned_ndays.append(0.0)
    return binned_ang,binned_alt,binned_num,binned_ndays


# In[370]:


i6 = np.isfinite(gap6['dalt'])
i7 = np.isfinite(gap7['dalt'])
i8 = np.isfinite(gap8['dalt'])


# In[371]:


bins = np.linspace(0.5,-24.5,26)


# In[372]:


bgap6_alt,bgap6_lat,bgap6_num,bgap6_ndays = make_bined_x(gap6['dalt'],
                                                           gap6['dlat'],gap6['dlat_ndays'],i6,bins=bins)
bgap7_alt,bgap7_lat,bgap7_num,bgap7_ndays = make_bined_x(gap7['dalt'],
                                                           gap7['dlat'],gap7['dlat_ndays'],i7,bins=bins)
bgap8_alt,bgap8_lat,bgap8_num,bgap8_ndays = make_bined_x(gap8['dalt'],
                                                           gap8['dlat'],gap8['dlat_ndays'],i8,bins=bins)


# In[434]:


io6 = np.isfinite(gap6['dalt']) & gap6['flo']
io7 = np.isfinite(gap7['dalt']) & gap7['flo']
io8 = np.isfinite(gap8['dalt']) & gap8['flo']

ir6 = np.isfinite(gap6['dalt']) & gap6['flr']
ir7 = np.isfinite(gap7['dalt']) & gap7['flr']
ir8 = np.isfinite(gap8['dalt']) & gap8['flr']


# In[435]:


bgapo6_alt,bgapo6_lat,bgapo6_num,bgapo6_ndays = make_bined_x(gap6['dalt'],
                                                           gap6['dlat'],gap6['dlat_ndays'],io6,bins=bins)
bgapo7_alt,bgapo7_lat,bgapo7_num,bgapo7_ndays = make_bined_x(gap7['dalt'],
                                                           gap7['dlat'],gap7['dlat_ndays'],io7,bins=bins)
bgapo8_alt,bgapo8_lat,bgapo8_num,bgapo8_ndays = make_bined_x(gap8['dalt'],
                                                           gap8['dlat'],gap8['dlat_ndays'],io8,bins=bins)
bgapr6_alt,bgapr6_lat,bgapr6_num,bgapr6_ndays = make_bined_x(gap6['dalt'],
                                                           gap6['dlat'],gap6['dlat_ndays'],ir6,bins=bins)
bgapr7_alt,bgapr7_lat,bgapr7_num,bgapr7_ndays = make_bined_x(gap7['dalt'],
                                                           gap7['dlat'],gap7['dlat_ndays'],ir7,bins=bins)
bgapr8_alt,bgapr8_lat,bgapr8_num,bgapr8_ndays = make_bined_x(gap8['dalt'],
                                                           gap8['dlat'],gap8['dlat_ndays'],ir8,bins=bins)


# ### Redo or separating routine and others

# In[436]:


ar6['fl_acaod_noQAr'] = (ar6['flag_acaod']==1) & (ar6['fl_routine'])
ii_flacaodr6 = np.where(ar6['fl_acaod_noQAr'])[0] 

ar7['fl_acaod_noQAr'] = (ar7['flag_acaod']==1) & (ar7['fl_routine'])
ii_flacaodr7 = np.where(ar7['fl_acaod_noQAr'])[0] 

ar8['fl_acaod_noQAr'] = (ar8['flag_acaod']==1) & (ar8['fl_routine'])
ii_flacaodr8 = np.where(ar8['fl_acaod_noQAr'])[0] 

gapr6 = get_gap(ii_flacaodr6,ar6['GPS_Alt'],ar6['Latitude'],ar6['Longitude'],ar6['days'],ar6['AOD0501'],ar6['AOD_angstrom_470_865'])
gapr7 = get_gap(ii_flacaodr7,ar7['GPS_Alt'],ar7['Latitude'],ar7['Longitude'],ar7['days'],ar7['AOD0501'],ar7['AOD_angstrom_470_865'])
gapr8 = get_gap(ii_flacaodr8,ar8['GPS_Alt'],ar8['Latitude'],ar8['Longitude'],ar8['days'],ar8['AOD0501'],ar8['AOD_angstrom_470_865'])


# In[437]:


ir6 = np.isfinite(gapr6['dalt'])
ir7 = np.isfinite(gapr7['dalt'])
ir8 = np.isfinite(gapr8['dalt'])


# In[438]:


bgapr6_alt,bgapr6_lat,bgapr6_num,bgapr6_ndays = make_bined_x(gapr6['dalt'],
                                                           gapr6['dlat'],gapr6['dlat_ndays'],ir6,bins=bins)
bgapr7_alt,bgapr7_lat,bgapr7_num,bgapr7_ndays = make_bined_x(gapr7['dalt'],
                                                           gapr7['dlat'],gapr7['dlat_ndays'],ir7,bins=bins)
bgapr8_alt,bgapr8_lat,bgapr8_num,bgapr8_ndays = make_bined_x(gapr8['dalt'],
                                                           gapr8['dlat'],gapr8['dlat_ndays'],ir8,bins=bins)


# In[439]:


ar6['fl_acaod_noQAo'] = (ar6['flag_acaod']==1) & (~ar6['fl_routine'])
ii_flacaodo6 = np.where(ar6['fl_acaod_noQAo'])[0] 

ar7['fl_acaod_noQAo'] = (ar7['flag_acaod']==1) & (~ar7['fl_routine'])
ii_flacaodo7 = np.where(ar7['fl_acaod_noQAo'])[0] 

ar8['fl_acaod_noQAo'] = (ar8['flag_acaod']==1) & (~ar8['fl_routine'])
ii_flacaodo8 = np.where(ar8['fl_acaod_noQAo'])[0] 

gapo6 = get_gap(ii_flacaodo6,ar6['GPS_Alt'],ar6['Latitude'],ar6['Longitude'],ar6['days'],ar6['AOD0501'],ar6['AOD_angstrom_470_865'])
gapo7 = get_gap(ii_flacaodo7,ar7['GPS_Alt'],ar7['Latitude'],ar7['Longitude'],ar7['days'],ar7['AOD0501'],ar7['AOD_angstrom_470_865'])
gapo8 = get_gap(ii_flacaodo8,ar8['GPS_Alt'],ar8['Latitude'],ar8['Longitude'],ar8['days'],ar8['AOD0501'],ar8['AOD_angstrom_470_865'])


# In[440]:


io6 = np.isfinite(gapo6['dalt'])
io7 = np.isfinite(gapo7['dalt'])
io8 = np.isfinite(gapo8['dalt'])


# In[441]:


bgapo6_alt,bgapo6_lat,bgapo6_num,bgapo6_ndays = make_bined_x(gapo6['dalt'],
                                                           gapo6['dlat'],gapo6['dlat_ndays'],io6,bins=bins)
bgapo7_alt,bgapo7_lat,bgapo7_num,bgapo7_ndays = make_bined_x(gapo7['dalt'],
                                                           gapo7['dlat'],gapo7['dlat_ndays'],io7,bins=bins)
bgapo8_alt,bgapo8_lat,bgapo8_num,bgapo8_ndays = make_bined_x(gapo8['dalt'],
                                                           gapo8['dlat'],gapo8['dlat_ndays'],io8,bins=bins)


# ## Plot the gap distances

# In[71]:


plt.figure()
plt.plot(gap6['dalt'],gap6['dlat'],'r.',label='2016')
plt.plot(gap7['dalt'],gap7['dlat'],'gx',label='2017')
plt.plot(gap8['dalt'],gap8['dlat'],'b+',label='2018')
plt.ylim(-23,0.5)
plt.legend(frameon=False)
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.xlabel('Gap extent [m]')
plt.title('Gap beween cloud top and aerosol bottom')


# ### The histogram distribution of gaps for all years

# In[226]:


plt.figure(figsize=(5,8))
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues

ax0 = plt.gca()

bp6 = plt.boxplot(bgap6_alt,positions=np.array(bgap6_lat)+0.2,vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(gr,bp6,bgap6_ndays)
plt.plot([x.get_data()[0][0] for x in bp6['means']],[x.get_data()[1][0] for x in bp6['means']],'r-')

bp7 = plt.boxplot(bgap7_alt,positions=np.array(bgap7_lat)-0.2,vert=False,widths=0.4,
                  showfliers=False,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,bgap7_ndays)
plt.plot([x.get_data()[0][0] for x in bp7['means']],[x.get_data()[1][0] for x in bp7['means']],'g-')

bp8 = plt.boxplot(bgap8_alt,positions=np.array(bgap8_lat),vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,bgap8_ndays)
plt.plot([x.get_data()[0][0] for x in bp8['means']],[x.get_data()[1][0] for x in bp8['means']],'b-')

plt.xlabel('Gap Extent [m]')
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.legend([bp6['boxes'][5],bp7['boxes'][18],bp8['boxes'][18],bp6['means'][0],bp6['medians'][0],bp6['boxes'][0],
            bp6['whiskers'][0]],
           ['2016','2017','2018','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bgap6_ndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bgap6_ndays)
scalarmapbr = plt.cm.ScalarMappable(cmap=br)
scalarmapbr.set_array(bgap6_ndays)

cbaxesgr = plt.gcf().add_axes([0.63, 0.51, 0.015, 0.2])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbg.outline.set_visible(False)
cbaxesbl = plt.gcf().add_axes([0.65, 0.51, 0.015, 0.2])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbb.outline.set_visible(False)
cbaxesbr = plt.gcf().add_axes([0.67, 0.51, 0.015, 0.2])
cbr = plt.colorbar(scalarmapbr,cax=cbaxesbr)
cbr.outline.set_visible(False)
cbg.set_ticks([0,3,6,9,12,15]),cbg.set_ticklabels(['','','','',''])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbr.set_ticks([0,3,6,9,12,15])
cbaxesgr.yaxis.set_ticks_position('right'),cbaxesbl.yaxis.set_ticks_position('right'),cbaxesbr.yaxis.set_ticks_position('right')
cbaxesgr.spines['right'].set_visible(False), cbaxesbl.spines['left'].set_visible(False)
cbaxesbr.text(4.0,0.5,'Days sampled',rotation=-90,verticalalignment='center')

ax0.set_title('Gap extent through the years')
plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_lat.png',
            transparent=True,dpi=500)


# In[105]:


plt.figure(figsize=(5,8))
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues

ax0 = plt.gca()

bp6 = plt.boxplot(bgap6_alt,positions=np.array(bgap6_lat)+0.2,vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(gr,bp6,bgap6_ndays)
plt.plot([x.get_data()[0][0] for x in bp6['means']],[x.get_data()[1][0] for x in bp6['means']],'r-')

bp7 = plt.boxplot(bgap7_alt,positions=np.array(bgap7_lat)-0.2,vert=False,widths=0.4,
                  showfliers=False,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,bgap7_ndays)
plt.plot([x.get_data()[0][0] for x in bp7['means']],[x.get_data()[1][0] for x in bp7['means']],'g-')

bp8 = plt.boxplot(bgap8_alt,positions=np.array(bgap8_lat),vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,bgap8_ndays)
plt.plot([x.get_data()[0][0] for x in bp8['means']],[x.get_data()[1][0] for x in bp8['means']],'b-')

pse = plt.plot(strip_elev,elev_lat,'-k',label='Surface elevation',alpha=0.7)

plt.xlabel('Gap Extent [m]')
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.legend([bp6['boxes'][5],bp7['boxes'][18],bp8['boxes'][18],bp6['means'][0],bp6['medians'][0],bp6['boxes'][0],
            bp6['whiskers'][0],pse[0]],
           ['2016','2017','2018','Mean','Median','25% - 75%','min-max','Surface\nelevation'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bgap6_ndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bgap6_ndays)
scalarmapbr = plt.cm.ScalarMappable(cmap=br)
scalarmapbr.set_array(bgap6_ndays)

cbaxesgr = plt.gcf().add_axes([0.69, 0.41, 0.015, 0.2])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbg.outline.set_visible(False)
cbaxesbl = plt.gcf().add_axes([0.71, 0.41, 0.015, 0.2])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbb.outline.set_visible(False)
cbaxesbr = plt.gcf().add_axes([0.73, 0.41, 0.015, 0.2])
cbr = plt.colorbar(scalarmapbr,cax=cbaxesbr)
cbr.outline.set_visible(False)
cbg.set_ticks([0,3,6,9,12,15]),cbg.set_ticklabels(['','','','',''])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbr.set_ticks([0,3,6,9,12,15])
cbaxesgr.yaxis.set_ticks_position('right'),cbaxesbl.yaxis.set_ticks_position('right'),cbaxesbr.yaxis.set_ticks_position('right')
cbaxesgr.spines['right'].set_visible(False), cbaxesbl.spines['left'].set_visible(False)
cbaxesbr.text(4.0,0.5,'Days sampled',rotation=-90,verticalalignment='center')

ax0.set_title('Gap extent through the years')
plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_lat_surf_elev.png',
            transparent=True,dpi=500)


# In[386]:


plt.figure(figsize=(5,8))
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues

ax0 = plt.gca()

bp6 = plt.boxplot(bgap6_alt,positions=np.array(bgap6_lat)+0.2,vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(gr,bp6,bgap6_ndays)
mean_gap6 = [x.get_data()[0][0] for x in bp6['means']]
plt.plot([x.get_data()[0][0] for x in bp6['means']],[x.get_data()[1][0] for x in bp6['means']],'r-')

bp7 = plt.boxplot(bgap7_alt,positions=np.array(bgap7_lat)-0.2,vert=False,widths=0.4,
                  showfliers=False,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,bgap7_ndays)
plt.plot([x.get_data()[0][0] for x in bp7['means']],[x.get_data()[1][0] for x in bp7['means']],'g-')
mean_gap7 = [x.get_data()[0][0] for x in bp7['means']]
bp8 = plt.boxplot(bgap8_alt,positions=np.array(bgap8_lat),vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,bgap8_ndays)
plt.plot([x.get_data()[0][0] for x in bp8['means']],[x.get_data()[1][0] for x in bp8['means']],'b-')
mean_gap8 = [x.get_data()[0][0] for x in bp8['means']]
pse = plt.plot(strip_elev,elev_lat,'-k',label='Surface elevation',alpha=0.7)
pse6 = plt.plot(mean_fire_elev6,center_lat,'.--',color='tab:red',label='2016')
pse7 = plt.plot(mean_fire_elev7,center_lat,'.--',color='tab:green',label='2017')
pse8 = plt.plot(mean_fire_elev8,center_lat,'.--',color='tab:blue',label='2018')

plt.xlabel('Gap Extent [m]')
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.legend([bp6['boxes'][5],bp7['boxes'][18],bp8['boxes'][18],bp6['means'][0],bp6['medians'][0],bp6['boxes'][0],
            bp6['whiskers'][0],pse[0],pse6[0],pse7[0],pse8[0]],
           ['2016','2017','2018','Mean','Median','25% - 75%','min-max','Surface\nelevation','Mean fire\nelev 2016',
            'Mean fire\nelev 2017','Mean fire\nelev 2018'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bgap6_ndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bgap6_ndays)
scalarmapbr = plt.cm.ScalarMappable(cmap=br)
scalarmapbr.set_array(bgap6_ndays)

cbaxesgr = plt.gcf().add_axes([0.67, 0.27, 0.015, 0.2])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbg.outline.set_visible(False)
cbaxesbl = plt.gcf().add_axes([0.69, 0.27, 0.015, 0.2])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbb.outline.set_visible(False)
cbaxesbr = plt.gcf().add_axes([0.71, 0.27, 0.015, 0.2])
cbr = plt.colorbar(scalarmapbr,cax=cbaxesbr)
cbr.outline.set_visible(False)
cbg.set_ticks([0,3,6,9,12,15]),cbg.set_ticklabels(['','','','',''])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbr.set_ticks([0,3,6,9,12,15])
cbaxesgr.yaxis.set_ticks_position('right'),cbaxesbl.yaxis.set_ticks_position('right'),cbaxesbr.yaxis.set_ticks_position('right')
cbaxesgr.spines['right'].set_visible(False), cbaxesbl.spines['left'].set_visible(False)
cbaxesbr.text(4.0,0.5,'Days sampled',rotation=-90,verticalalignment='center')

ax0.set_title('Gap extent through the years')
plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_lat_surf_elev_perfire.png',
            transparent=True,dpi=500)


# In[459]:


lat_bins


# In[460]:


binned_elev = []
for i,la1 in enumerate(lat_bins[:-1]):
    ia = (elev_lat<la1) & (elev_lat>lat_bins[i+1])
    binned_elev.append(np.nanmean(strip_elev[ia]))
binned_elev = np.array(binned_elev)


# In[419]:


i6 = np.isfinite(np.flip(mean_gap6)) & np.isfinite(mean_fire_elev6)
r6 = np.corrcoef(np.flip(mean_gap6)[i6],np.array(mean_fire_elev6)[i6])[0,1]
i7 = np.isfinite(np.flip(mean_gap7)) & np.isfinite(mean_fire_elev7)
r7 = np.corrcoef(np.flip(mean_gap7)[i7],np.array(mean_fire_elev7)[i7])[0,1]
i8 = np.isfinite(np.flip(mean_gap8)) & np.isfinite(mean_fire_elev8)
r8 = np.corrcoef(np.flip(mean_gap8)[i8],np.array(mean_fire_elev8)[i8])[0,1]


# In[420]:


r6,r7,r8


# In[421]:


i6 = np.isfinite(np.flip(mean_gap6)) & np.isfinite(median_fire_elev6)
rm6 = np.corrcoef(np.flip(mean_gap6)[i6],np.array(median_fire_elev6)[i6])[0,1]
i7 = np.isfinite(np.flip(mean_gap7)) & np.isfinite(median_fire_elev7)
rm7 = np.corrcoef(np.flip(mean_gap7)[i7],np.array(median_fire_elev7)[i7])[0,1]
i8 = np.isfinite(np.flip(mean_gap8)) & np.isfinite(median_fire_elev8)
rm8 = np.corrcoef(np.flip(mean_gap8)[i8],np.array(median_fire_elev8)[i8])[0,1]


# In[422]:


rm6,rm7,rm8


# In[462]:


i6 = np.isfinite(np.flip(mean_gap6)) & np.isfinite(binned_elev)
rb6 = np.corrcoef(np.flip(mean_gap6)[i6],np.array(binned_elev)[i6])[0,1]
i7 = np.isfinite(np.flip(mean_gap7)) & np.isfinite(binned_elev)
rb7 = np.corrcoef(np.flip(mean_gap7)[i7],np.array(binned_elev)[i7])[0,1]
i8 = np.isfinite(np.flip(mean_gap8)) & np.isfinite(binned_elev)
rb8 = np.corrcoef(np.flip(mean_gap8)[i8],np.array(binned_elev)[i8])[0,1]
rb6,rb7,rb8


# In[423]:


plt.figure()
plt.plot(np.flip(mean_gap6),median_fire_elev6,'.',color='tab:red',label='2016')
plt.plot(np.flip(mean_gap7),median_fire_elev7,'.',color='tab:green',label='2017')
plt.plot(np.flip(mean_gap8),median_fire_elev8,'.',color='tab:blue',label='2018')
plt.xlabel('Mean gap extent [m]')
plt.ylabel('Median fire elevation [m]')
plt.legend(frameon=False)
plt.savefig(fp+'plot_all/ORACLES_gap_vs_fireelev.png',
            transparent=True,dpi=500)


# In[170]:


plt.figure()
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues

bp6 = plt.boxplot(bgap6_alt,positions=np.array(bgap6_lat)+0.2,vert=True,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(gr,bp6,bgap6_ndays)

bp7 = plt.boxplot(bgap7_alt,positions=np.array(bgap7_lat)-0.2,vert=True,widths=0.4,
                  showfliers=False,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,bgap7_ndays)

bp8 = plt.boxplot(bgap8_alt,positions=np.array(bgap8_lat),vert=True,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,bgap8_ndays)


plt.ylabel('Gap Extent [m]')
plt.xlabel('Latitude [$^{{\circ}}$]')
plt.legend([bp6['boxes'][5],bp7['boxes'][18],bp8['boxes'][18],bp6['means'][0],bp6['medians'][0],bp6['boxes'][0],
            bp6['whiskers'][0]],
           ['2016','2017','2018','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)


# ### Subset the gaps for near coast / vs far

# In[444]:


plt.figure(figsize=(5,8))
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues

ax0 = plt.gca()

bp6 = plt.boxplot(bgapo6_alt,positions=np.array(bgapo6_lat)+0.2,vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(gr,bp6,bgapo6_ndays)
plt.plot([x.get_data()[0][0] for x in bp6['means']],[x.get_data()[1][0] for x in bp6['means']],'r-')

bp7 = plt.boxplot(bgapr7_alt,positions=np.array(bgapr7_lat)-0.2,vert=False,widths=0.4,
                  showfliers=False,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,bgapr7_ndays)
plt.plot([x.get_data()[0][0] for x in bp7['means']],[x.get_data()[1][0] for x in bp7['means']],'g-')

bp8 = plt.boxplot(bgapo8_alt,positions=np.array(bgapo8_lat),vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,bgapo8_ndays)
plt.plot([x.get_data()[0][0] for x in bp8['means']],[x.get_data()[1][0] for x in bp8['means']],'b-')


meanc_gap6 = [x.get_data()[0][0] for x in bp6['means']]
meanc_gap7 = [x.get_data()[0][0] for x in bp7['means']]
meanc_gap8 = [x.get_data()[0][0] for x in bp8['means']]

pse = plt.plot(strip_elev,elev_lat,'-k',label='Surface elevation',alpha=0.7)
pse6 = plt.plot(mean_fire_elev6,center_lat,'.--',color='tab:red',label='2016')
pse7 = plt.plot(mean_fire_elev7,center_lat,'.--',color='tab:green',label='2017')
pse8 = plt.plot(mean_fire_elev8,center_lat,'.--',color='tab:blue',label='2018')

plt.xlabel('Gap Extent [m]')
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.legend([bp6['boxes'][5],bp7['boxes'][18],bp8['boxes'][18],bp6['means'][0],bp6['medians'][0],bp6['boxes'][0],
            bp6['whiskers'][0]],
           ['2016','2017','2018','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bgap6_ndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bgap6_ndays)
scalarmapbr = plt.cm.ScalarMappable(cmap=br)
scalarmapbr.set_array(bgap6_ndays)

cbaxesgr = plt.gcf().add_axes([0.63, 0.51, 0.015, 0.2])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbg.outline.set_visible(False)
cbaxesbl = plt.gcf().add_axes([0.65, 0.51, 0.015, 0.2])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbb.outline.set_visible(False)
cbaxesbr = plt.gcf().add_axes([0.67, 0.51, 0.015, 0.2])
cbr = plt.colorbar(scalarmapbr,cax=cbaxesbr)
cbr.outline.set_visible(False)
cbg.set_ticks([0,3,6,9,12,15]),cbg.set_ticklabels(['','','','',''])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbr.set_ticks([0,3,6,9,12,15])
cbaxesgr.yaxis.set_ticks_position('right'),cbaxesbl.yaxis.set_ticks_position('right'),cbaxesbr.yaxis.set_ticks_position('right')
cbaxesgr.spines['right'].set_visible(False), cbaxesbl.spines['left'].set_visible(False)
cbaxesbr.text(4.0,0.5,'Days sampled',rotation=-90,verticalalignment='center')

ax0.set_title('Gap extent Coastal')
plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_lat_coast.png',
            transparent=True,dpi=500)


# In[449]:


i6 = np.isfinite(np.flip(meanc_gap6)) & np.isfinite(median_fire_elev6)
rm6 = np.corrcoef(np.flip(meanc_gap6)[i6],np.array(median_fire_elev6)[i6])[0,1]
i7 = np.isfinite(np.flip(meanc_gap7)) & np.isfinite(median_fire_elev7)
rm7 = np.corrcoef(np.flip(meanc_gap7)[i7],np.array(median_fire_elev7)[i7])[0,1]
i8 = np.isfinite(np.flip(meanc_gap8)) & np.isfinite(median_fire_elev8)
rm8 = np.corrcoef(np.flip(meanc_gap8)[i8],np.array(median_fire_elev8)[i8])[0,1]


# In[450]:


rm6,rm7,rm8


# In[443]:


plt.figure()
plt.plot(np.flip(meanc_gap6),median_fire_elev6,'.',color='tab:red',label='2016')
plt.plot(np.flip(meanc_gap7),median_fire_elev7,'.',color='tab:green',label='2017')
plt.plot(np.flip(meanc_gap8),median_fire_elev8,'.',color='tab:blue',label='2018')
plt.xlabel('Mean gap extent [m]')
plt.ylabel('Median fire elevation [m]')
plt.legend(frameon=False)
plt.savefig(fp+'plot_all/ORACLES_coastalgap_vs_fireelev.png',
            transparent=True,dpi=500)


# In[451]:


from scipy.interpolate import interp1d


# In[454]:


elev_fs = interp1d(elev_lat,strip_elev,fill_value='extrapolate')


# In[455]:


new_elev = elev_fs(center_lat)


# In[456]:


plt.figure()
plt.plot(np.flip(meanc_gap6),new_elev,'.',color='tab:red',label='2016')
plt.plot(np.flip(meanc_gap7),new_elev,'.',color='tab:green',label='2017')
plt.plot(np.flip(meanc_gap8),new_elev,'.',color='tab:blue',label='2018')
plt.xlabel('Mean gap extent [m]')
plt.ylabel('Mean Surface elevation [m]')
plt.legend(frameon=False)
plt.savefig(fp+'plot_all/ORACLES_coastalgap_vs_surfelev.png',
            transparent=True,dpi=500)


# In[457]:


i6 = np.isfinite(np.flip(meanc_gap6)) & np.isfinite(new_elev)
re6 = np.corrcoef(np.flip(meanc_gap6)[i6],np.array(new_elev)[i6])[0,1]
i7 = np.isfinite(np.flip(meanc_gap7)) & np.isfinite(new_elev)
re7 = np.corrcoef(np.flip(meanc_gap7)[i7],np.array(new_elev)[i7])[0,1]
i8 = np.isfinite(np.flip(meanc_gap8)) & np.isfinite(new_elev)
re8 = np.corrcoef(np.flip(meanc_gap8)[i8],np.array(new_elev)[i8])[0,1]


# In[458]:


re6,re7,re8


# In[463]:


i6 = np.isfinite(np.flip(meanc_gap6)) & np.isfinite(binned_elev)
re6 = np.corrcoef(np.flip(meanc_gap6)[i6],np.array(binned_elev)[i6])[0,1]
i7 = np.isfinite(np.flip(meanc_gap7)) & np.isfinite(binned_elev)
re7 = np.corrcoef(np.flip(meanc_gap7)[i7],np.array(binned_elev)[i7])[0,1]
i8 = np.isfinite(np.flip(meanc_gap8)) & np.isfinite(binned_elev)
re8 = np.corrcoef(np.flip(meanc_gap8)[i8],np.array(binned_elev)[i8])[0,1]
re6,re7,re8


# In[464]:


plt.figure(figsize=(5,8))
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues

ax0 = plt.gca()

bp6 = plt.boxplot(bgapr6_alt,positions=np.array(bgapr6_lat)+0.2,vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(gr,bp6,bgapr6_ndays)
plt.plot([x.get_data()[0][0] for x in bp6['means']],[x.get_data()[1][0] for x in bp6['means']],'r-')

bp7 = plt.boxplot(bgapo7_alt,positions=np.array(bgapo7_lat)-0.2,vert=False,widths=0.4,
                  showfliers=False,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,bgapo7_ndays)
plt.plot([x.get_data()[0][0] for x in bp7['means']],[x.get_data()[1][0] for x in bp7['means']],'g-')

bp8 = plt.boxplot(bgapr8_alt,positions=np.array(bgapr8_lat),vert=False,
                 showfliers=False,widths=0.4,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,bgapr8_ndays)
plt.plot([x.get_data()[0][0] for x in bp8['means']],[x.get_data()[1][0] for x in bp8['means']],'b-')

meano_gap6 = [x.get_data()[0][0] for x in bp6['means']]
meano_gap7 = [x.get_data()[0][0] for x in bp7['means']]
meano_gap8 = [x.get_data()[0][0] for x in bp8['means']]

plt.xlabel('Gap Extent [m]')
plt.ylabel('Latitude [$^{{\circ}}$]')
plt.legend([bp6['boxes'][5],bp7['boxes'][18],bp8['boxes'][18],bp6['means'][0],bp6['medians'][0],bp6['boxes'][0],
            bp6['whiskers'][0]],
           ['2016','2017','2018','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bgap6_ndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bgap6_ndays)
scalarmapbr = plt.cm.ScalarMappable(cmap=br)
scalarmapbr.set_array(bgap6_ndays)

cbaxesgr = plt.gcf().add_axes([0.63, 0.51, 0.015, 0.2])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbg.outline.set_visible(False)
cbaxesbl = plt.gcf().add_axes([0.65, 0.51, 0.015, 0.2])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbb.outline.set_visible(False)
cbaxesbr = plt.gcf().add_axes([0.67, 0.51, 0.015, 0.2])
cbr = plt.colorbar(scalarmapbr,cax=cbaxesbr)
cbr.outline.set_visible(False)
cbg.set_ticks([0,3,6,9,12,15]),cbg.set_ticklabels(['','','','',''])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbr.set_ticks([0,3,6,9,12,15])
cbaxesgr.yaxis.set_ticks_position('right'),cbaxesbl.yaxis.set_ticks_position('right'),cbaxesbr.yaxis.set_ticks_position('right')
cbaxesgr.spines['right'].set_visible(False), cbaxesbl.spines['left'].set_visible(False)
cbaxesbr.text(4.0,0.5,'Days sampled',rotation=-90,verticalalignment='center')

ax0.set_title('Gap extent Oceanic')
plt.tight_layout()

plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_lat_oceanic.png',
            transparent=True,dpi=500)


# In[467]:


i6 = np.isfinite(np.flip(meano_gap6)) & np.isfinite(binned_elev)
re6 = np.corrcoef(np.flip(meano_gap6)[i6],np.array(binned_elev)[i6])[0,1]
i7 = np.isfinite(np.flip(meano_gap7)) & np.isfinite(binned_elev)
re7 = np.corrcoef(np.flip(meano_gap7)[i7],np.array(binned_elev)[i7])[0,1]
i8 = np.isfinite(np.flip(meano_gap8)) & np.isfinite(binned_elev)
re8 = np.corrcoef(np.flip(meano_gap8)[i8],np.array(binned_elev)[i8])[0,1]
re6,re7,re8


# ## compare gap distance to AOD

# In[72]:


import plotting_utils as pu


# In[73]:


plt.figure()
plt.plot(gap6['daod'],gap6['dalt'],'r.',label='2016')
plt.plot(gap7['daod'],gap7['dalt'],'gx',label='2017')
plt.plot(gap8['daod'],gap8['dalt'],'b+',label='2018')
plt.xlabel('AOD 501 nm')
plt.ylabel('Gap extent [m]')
plt.title('Gap extent as a function of AOD')
plt.xlim([0,0.8])
plt.legend(frameon=False)


# In[28]:


plt.figure()
plt.plot(gap6['daod'],gap6['dalt'],'r.',label='2016')


# In[29]:


plt.figure()
plt.plot(gap7['daod'],gap7['dalt'],'gx',label='2017')


# In[30]:


plt.figure()
plt.plot(gap8['daod'],gap8['dalt'],'b+',label='2018')


# In[374]:


plt.figure()
plt.plot(1.0/gap6['daod'],gap6['dalt'],'r.',label='2016')
pu.plot_lin(1.0/gap6['daod'],gap6['dalt'],color='r')
plt.plot(1.0/gap7['daod'],gap7['dalt'],'gx',label='2017')
pu.plot_lin(1.0/gap7['daod'],gap7['dalt'],color='g')
plt.plot(1.0/gap8['daod'],gap8['dalt'],'b+',label='2018')
pu.plot_lin(1.0/gap8['daod'],gap8['dalt'],color='b')
plt.xlabel('1/ AOD 501 nm')
plt.ylabel('Gap extent [m]')
plt.title('Gap extent as a function of 1/AOD')
plt.xlim([0,20])
plt.legend(frameon=False)


# In[74]:


plt.figure()
plt.plot(gap6['daod'],gap6['dalt'],'r.',label='2016')
pu.plot_lin(gap6['daod'],gap6['dalt'],color='r')
plt.plot(gap7['daod'],gap7['dalt'],'gx',label='2017')
pu.plot_lin(gap7['daod'],gap7['dalt'],color='g')
plt.plot(gap8['daod'],gap8['dalt'],'b+',label='2018')
pu.plot_lin(gap8['daod'],gap8['dalt'],color='b')
plt.xlabel('AOD 501 nm')
plt.ylabel('Gap extent [m]')
plt.title('Gap extent as a function of AOD')
plt.xlim([0,1.2])
plt.ylim([-5,2000])
plt.legend(frameon=False)

plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_vs_AOD_points.png',
            transparent=True,dpi=500)


# In[85]:


np.corrcoef(gap6['daod'],gap6['dalt'])[0,1]**2,np.corrcoef(gap7['daod'],gap7['dalt'])[0,1]**2,np.corrcoef(gap8['daod'],gap8['dalt'])[0,1]**2


# In[76]:


plt.figure()
plt.plot(gap6['dang'],gap6['dalt'],'r.',label='2016')
pu.plot_lin(gap6['dang'],gap6['dalt'],color='r')
plt.plot(gap7['dang'],gap7['dalt'],'gx',label='2017')
pu.plot_lin(gap7['dang'],gap7['dalt'],color='g')
plt.plot(gap8['dang'],gap8['dalt'],'b+',label='2018')
pu.plot_lin(gap8['dang'],gap8['dalt'],color='b')
plt.xlabel('Angstrom Exponent')
plt.ylabel('Gap extent [m]')
plt.title('Gap extent as a function of Angstrom')
plt.xlim([0,1.2])
plt.ylim([-5,2000])
plt.legend(frameon=False)

plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_vs_Angstrom_points.png',
            transparent=True,dpi=500)


# In[84]:


np.corrcoef(gap6['dang'],gap6['dalt'])[0,1]**2,np.corrcoef(gap7['dang'],gap7['dalt'])[0,1]**2,np.corrcoef(gap8['dang'],gap8['dalt'])[0,1]**2


# ## Build AOD distributions vs. gap extent

# In[318]:


abins = np.logspace(-3,0.30103,30)


# In[319]:


abins


# In[320]:


i6 = np.isfinite(gap6['dalt'])
i7 = np.isfinite(gap7['dalt'])
i8 = np.isfinite(gap8['dalt'])


# In[321]:


bgapa6_alt,bgapa6_aod,bgapa6_num,bgapa6_ndays = make_bined_x(gap6['dalt'],
                                                           gap6['daod'],gap6['dlat_ndays'],i6,bins=abins)
bgapa7_alt,bgapa7_aod,bgapa7_num,bgapa7_ndays = make_bined_x(gap7['dalt'],
                                                           gap7['daod'],gap7['dlat_ndays'],i7,bins=abins)
bgapa8_alt,bgapa8_aod,bgapa8_num,bgapa8_ndays = make_bined_x(gap8['dalt'],
                                                           gap8['daod'],gap8['dlat_ndays'],i8,bins=abins)


# In[349]:


plt.figure(figsize=(12,5))
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues


bp6 = plt.boxplot(bgapa6_alt,positions=np.array(bgapa6_aod)+0.05,vert=True,
                 showfliers=False,widths=0.05,showmeans=True,patch_artist=True)
plt.xlim(0.001,2.0)
set_box_whisker_color(gr,bp6,bgapa6_ndays)
plt.plot([x.get_data()[0][0] for x in bp6['means']],[x.get_data()[1][0] for x in bp6['means']],'r-')
ax = plt.gca()


bp7 = plt.boxplot(bgapa7_alt,positions=np.array(bgapa7_aod)-0.05,vert=True,
                 showfliers=False,widths=0.05,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,bgapa7_ndays)
plt.plot([x.get_data()[0][0] for x in bp7['means']],[x.get_data()[1][0] for x in bp7['means']],'g-')

bp8 = plt.boxplot(bgapa8_alt,positions=np.array(bgapa8_aod),vert=True,
                 showfliers=False,widths=0.05,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,bgapa8_ndays)
plt.plot([x.get_data()[0][0] for x in bp8['means']],[x.get_data()[1][0] for x in bp8['means']],'b-')

plt.xticks([0.001,0.1,0.2,0.3,0.5,0.8,1.0,1.5,2.0])
ax.set_xticklabels([0.001,0.1,0.2,0.3,0.5,0.8,1.0,1.5,2.0])
plt.xlim(0.001,1.0)
plt.ylabel('Gap extent [m]')
plt.xlabel('AOD 501 nm')
plt.title('Gap Extent dependence on AOD')
plt.grid()
plt.legend([bp6['boxes'][17],bp7['boxes'][21],bp8['boxes'][21]],['2016','2017','2018'], frameon=False,loc=1,numpoints=1)


plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_vs_AOD.png',
            transparent=True,dpi=500)


# ## Build Angstrom distributions vs. gap extent
# 

# In[350]:


nbins = np.linspace(0.0,2.1,30)


# In[352]:


nbins


# In[353]:


bgapn6_alt,bgapn6_ang,bgapn6_num,bgapn6_ndays = make_bined_x(gap6['dalt'],
                                                           gap6['dang'],gap6['dlat_ndays'],i6,bins=nbins)
bgapn7_alt,bgapn7_ang,bgapn7_num,bgapn7_ndays = make_bined_x(gap7['dalt'],
                                                           gap7['dang'],gap7['dlat_ndays'],i7,bins=nbins)
bgapn8_alt,bgapn8_ang,bgapn8_num,bgapn8_ndays = make_bined_x(gap8['dalt'],
                                                           gap8['dang'],gap8['dlat_ndays'],i8,bins=nbins)


# In[362]:


plt.figure(figsize=(12,5))
gr = plt.cm.RdPu
bl = plt.cm.YlGn
br = plt.cm.Blues


bp6 = plt.boxplot(bgapn6_alt,positions=np.array(bgapn6_ang)+0.01,vert=True,
                 showfliers=False,widths=0.01,showmeans=True,patch_artist=True)
plt.xlim(0.001,2.0)
set_box_whisker_color(gr,bp6,bgapn6_ndays)
plt.plot([x.get_data()[0][0] for x in bp6['means']],[x.get_data()[1][0] for x in bp6['means']],'r-')
ax = plt.gca()


bp7 = plt.boxplot(bgapn7_alt,positions=np.array(bgapn7_ang)-0.01,vert=True,
                 showfliers=False,widths=0.01,showmeans=True,patch_artist=True)
set_box_whisker_color(bl,bp7,bgapn7_ndays)
plt.plot([x.get_data()[0][0] for x in bp7['means']],[x.get_data()[1][0] for x in bp7['means']],'g-')

bp8 = plt.boxplot(bgapn8_alt,positions=np.array(bgapn8_ang),vert=True,
                 showfliers=False,widths=0.01,showmeans=True,patch_artist=True)
set_box_whisker_color(br,bp8,bgapn8_ndays)
plt.plot([x.get_data()[0][0] for x in bp8['means']],[x.get_data()[1][0] for x in bp8['means']],'b-')

plt.xticks([0.5,0.8,1.2,1.5,1.8,2.0,2.1])
ax.set_xticklabels([0.5,0.8,1.2,1.5,1.8,2.0,2.1])
plt.xlim(0.5,2.2)
plt.ylim(0,2000)
plt.ylabel('Gap extent [m]')
plt.xlabel('Angstrom')
plt.title('Gap Extent dependence on Angstrom')
plt.grid()
plt.legend([bp6['boxes'][17],bp7['boxes'][21],bp8['boxes'][21]],['2016','2017','2018'], frameon=False,loc=1,numpoints=1)


plt.savefig(fp+'plot_all/ORACLESall_Gap_extent_all_vs_Angstrom.png',
            transparent=True,dpi=500)

