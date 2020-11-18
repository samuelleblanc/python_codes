#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     KORUS_AOD_fine_coarse_autocorr
#     
# Purpose:  
# 
#     Analyse some of the AOD values from KORUS AQ
#     Split up between fine mode and coarse mode AOD
#     Subset for level legs only
#         - interpolate within the level legs
#     Run autocorrelation values for the distance/time travelled 
#   
# Input:
# 
#     None at command line
#   
# Output:
# 
#     figures and save files...
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - hdf5 python loader
#     - 
#     - matplotlib
#     - numpy
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - '/aod_ict/all_aod_KORUS_R2_ict.mat'
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-05-18
#     Modified: 

# # Prepare python environment

# In[2]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
import os
matplotlib.rc_file(os.path.join(os.getcwd(),'file.rc'))
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
import load_utils as lu
import plotting_utils as pu
from path_utils import getpath
import hdf5storage as hs
from datetime import datetime
from scipy.interpolate import UnivariateSpline
import matplotlib.dates as mdates
from mpl_toolkits.basemap import Basemap
import scipy.stats as st


# In[3]:


import map_utils as mu
import write_utils as wu
from scipy import interpolate
import math
import sys


# In[4]:


from linfit import linfit
import Sun_utils as su


# In[5]:


get_ipython().magic(u'matplotlib notebook')


# In[6]:


fp =getpath('KORUS')


# # Load files

# Load the KORUS 4STAR AOD ict files for better handling

# ## Load 4STAR AOD ict files

# In[7]:


ar = hs.loadmat(fp+'/aod_ict/all_aod_KORUS_R2_ict.mat')


# In[8]:


ka = ar.keys()


# In[9]:


ka.sort()


# In[10]:


ka


# In[11]:


ar['qual_flag']


# In[11]:


ar['fl_QA']


# In[12]:


ar['qual_flag'].sum()/float(len(ar['qual_flag']))


# In[13]:


fl_ci = (ar['qual_flag']==1.) & (ar['AOD0501']<1.0)


# In[14]:


fl_ci.sum()


# In[15]:


len(fl_ci)*1.0/60.0/60.0


# In[16]:


fl_ci.sum()/float(len(fl_ci))


# In[17]:


ar['fl_QA'].sum()*1.0/60.0/60.0


# In[18]:


ar['qual_flag'].sum()


# In[19]:


fl_ci


# In[14]:


nwl = ka[0:17]


# In[15]:


nwl


# In[16]:


nm = [380.0,452.0,501.0,520.0,532.0,550.0,606.0,620.0,675.0,781.0,865.0,1020.0,1040.0,1064.0,1236.0,1559.0,1627.0]


# In[17]:


days = ['20160501','20160503','20160504','20160506','20160510','20160511',
        '20160512','20160516','20160517','20160519','20160521','20160524',
        '20160526','20160529','20160601','20160602','20160604','20160608',
        '20160609','20160614']


# In[18]:


doys = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in days]


# In[19]:


doys


# In[20]:


fdoys = np.array(doys)


# In[21]:


ar['doys'] = fdoys[ar['days'].astype(int)]+ar['Start_UTC']/24.0


# In[22]:


ar['doys']


# ## Load GOCI AOD files

# In[23]:


gl = os.listdir(fp+'data_other/GOCI')


# In[24]:


gl.sort()


# In[26]:


goci = []


# In[ ]:


for d in days:
    for q in xrange(8):
        print 'loading : GOCI_YAER_V2_AHICLD_AOPs_{d}0{q}.hdf'.format(d=d,q=q)
        g_tmp,g_dict = lu.load_hdf(fp+'data_other/GOCI/GOCI_YAER_V2_AHICLD_AOPs_{d}0{q}.hdf'.format(d=d,q=q),
                                   values=(('lat',1),('lon',0),('aod',2),('fmf',3),('ssa',4),('AE',6),('CF',10),('t',7)),
                                   verbose=False)
        goci.append(g_tmp)  


# In[27]:


for l in gl:
    print 'loading file: {}'.format(l)
    g_tmp,g_dict = lu.load_hdf(fp+'data_other/GOCI/{}'.format(l),
                                   values=(('lat',1),('lon',0),('aod',2),('fmf',3),('ssa',4),('AE',6),('CF',10),('t',7)),
                                   verbose=False)
    goci.append(g_tmp)  


# In[28]:


len(goci)


# In[29]:


np.nanmax(goci[0]['t']),np.nanmean(goci[0]['t']),np.nanmin(goci[0]['t'])


# In[40]:


g_dict


# In[30]:


gl[0][-14:-10],gl[0][-10:-8],gl[0][-8:-6]


# In[31]:


goci_doy = []
for i,l in enumerate(gl):
    goci[i]['doy'] = datetime(int(l[-14:-10]),int(l[-10:-8]),int(l[-8:-6])).timetuple().tm_yday+(int(l[-6:-4])+0.5)/24.0
    print goci[i]['doy']
    goci_doy.append(goci[i]['doy'])
    goci[i]['doys'] = goci[i]['doy']+(int(l[-6:-4])+goci[i]['t']/60.0)/24.0 
goci_doy = np.array(goci_doy)


# ### Build collcation of goci data to 4STAR

# In[32]:


ig4 = Sp.find_closest(goci_doy,ar['doys'])


# In[33]:


ig4


# In[34]:


plt.figure()
plt.plot(goci_doy[ig4],'.',label='goci')
plt.plot(ar['doys'],'x',label='4star')
plt.legend()


# In[35]:


bad_ig4 = abs(goci_doy[ig4]-ar['doys'])>(1.0/24.0)


# In[36]:


sum(bad_ig4)/float(len(ig4))


# In[37]:


len(np.where(bad_ig4)[0])


# In[38]:


len(ig4)


# In[50]:


ar['doys']%1.0


# In[39]:


goci[0].keys()


# In[40]:


goci2ar = {'doys':[],'lat':[],'lon':[],'aod':[],'AE':[],'fmf':[],'CF':[]}
for ii in np.unique(ig4):
    print ii
    imeas = np.where(ig4==ii)[0]
    igoci = mu.map_ind(goci[ii]['lon'],goci[ii]['lat'],ar['Longitude'][imeas],ar['Latitude'][imeas])
    if len(igoci)<1:
        goci2ar['lat'] = np.append(goci2ar['lat'],ar['Latitude'][imeas])
        goci2ar['lon'] = np.append(goci2ar['lon'],ar['Longitude'][imeas])
        goci2ar['aod'] = np.append(goci2ar['aod'],ar['Latitude'][imeas]+np.nan)
        goci2ar['AE'] = np.append(goci2ar['AE'],ar['Latitude'][imeas]+np.nan)
        goci2ar['fmf'] = np.append(goci2ar['fmf'],ar['Latitude'][imeas]+np.nan)
        goci2ar['CF'] = np.append(goci2ar['CF'],ar['Latitude'][imeas]+np.nan)
        goci2ar['doys'] = np.append(goci2ar['doys'],ar['doys'][imeas])
    else:
        goci2ar['lat'] = np.append(goci2ar['lat'],goci[ii]['lat'][igoci[0],igoci[1]])
        goci2ar['lon'] = np.append(goci2ar['lon'],goci[ii]['lon'][igoci[0],igoci[1]])
        goci2ar['aod'] = np.append(goci2ar['aod'],goci[ii]['aod'][igoci[0],igoci[1]])
        goci2ar['AE'] = np.append(goci2ar['AE'],goci[ii]['AE'][igoci[0],igoci[1]])
        goci2ar['fmf'] = np.append(goci2ar['fmf'],goci[ii]['fmf'][igoci[0],igoci[1]])
        goci2ar['CF'] = np.append(goci2ar['CF'],goci[ii]['CF'][igoci[0],igoci[1]])
        goci2ar['doys'] = np.append(goci2ar['doys'],goci[ii]['doys'][igoci[0],igoci[1]])

for k in goci2ar.keys():
    goci2ar[k] = np.array(goci2ar[k])


# In[53]:


igoci


# In[41]:


len(goci2ar['aod'])


# In[42]:


goci2ar['aod'][2000]


# In[57]:


np.unique(ig4)


# In[43]:


imeas = np.where(ig4==8)[0]


# In[44]:


len(imeas)


# In[45]:


np.argmin((goci[8]['lon']-ar['Longitude'][imeas[0]])**2.0+(goci[8]['lat']-ar['Latitude'][imeas[0]])**2.0)


# In[56]:


goci[8]['lon'].flatten()[94639],goci[8]['lat'].flatten()[94639]


# In[61]:


ar['Longitude'][imeas[0]],ar['Latitude'][imeas[0]]


# In[62]:


plt.figure()
plt.plot(goci2ar['aod'])
plt.figure()
plt.plot(goci2ar['AE'])


# In[46]:


goci2ar['aod'][bad_ig4] = np.nan
goci2ar['AE'][bad_ig4] = np.nan
goci2ar['fmf'][bad_ig4] = np.nan
goci2ar['CF'][bad_ig4] = np.nan


# ## Load MERRA2

# ### Old version
# 
# 
# For MERRA2 data
# Downloaded from https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NXGAS.5.12.4/doc/MERRA2.README.pdf
# https://disc.gsfc.nasa.gov/
#     
# on 2020-06-26
# 
# // global attributes:
#                 :History = "Original file generated: Thu May 19 23:35:04 2016 GMT" ;
#                 :Comment = "GMAO filename: d5124_m2_jan10.inst3_2d_gas_Nx.20160501.nc4" ;
#                 :Filename = "MERRA2_400.inst3_2d_gas_Nx.20160501.nc4" ;
#                 :Conventions = "CF-1" ;
#                 :Institution = "NASA Global Modeling and Assimilation Office" ;
#                 :References = "http://gmao.gsfc.nasa.gov" ;
#                 :Format = "NetCDF-4/HDF-5" ;
#                 :SpatialCoverage = "global" ;
#                 :VersionID = "5.12.4" ;
#                 :TemporalRange = "1980-01-01 -> 2016-12-31" ;
#                 :identifier_product_doi_authority = "http://dx.doi.org/" ;
#                 :ShortName = "M2I3NXGAS" ;
#                 :GranuleID = "MERRA2_400.inst3_2d_gas_Nx.20160501.nc4" ;
#                 :ProductionDateTime = "Original file generated: Thu May 19 23:35:04 2016 GMT" ;
#                 :LongName = "MERRA2 inst3_2d_gas_Nx: 2d,3-Hourly,Instantaneous,Single-Level,Assimilation,Aerosol Optical Depth Analysis" ;
#                 :Title = "MERRA2 inst3_2d_gas_Nx: 2d,3-Hourly,Instantaneous,Single-Level,Assimilation,Aerosol Optical Depth Analysis" ;
#                 :SouthernmostLatitude = "-90.0" ;
#                 :NorthernmostLatitude = "90.0" ;
#                 :WesternmostLongitude = "-180.0" ;
#                 :EasternmostLongitude = "179.375" ;
#                 :LatitudeResolution = "0.5" ;
#                 :LongitudeResolution = "0.625" ;
#                 :DataResolution = "0.5 x 0.625" ;
#                 :Source = "CVS tag: GEOSadas-5_12_4_p5" ;
#                 :Contact = "http://gmao.gsfc.nasa.gov" ;
#                 :identifier_product_doi = "10.5067/HNGA0EWW0R09" ;
#                 :RangeBeginningDate = "2016-05-01" ;
#                 :RangeBeginningTime = "00:00:00.000000" ;
#                 :RangeEndingDate = "2016-05-01" ;
#                 :RangeEndingTime = "21:00:00.000000" ;
#                 :DODS_EXTRA.Unlimited_Dimension = "time" ;
#                 :history = "2020-06-26 22:21:20 GMT Hyrax-1.15.1 https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I3NXGAS.5.12.4/2016/05/MERRA2_400.inst3_2d_gas_Nx.20160501.nc4.nc4?AODANA[0:7][243:259][484:499],time,lat[243:259],lon[484:499]" ;

# In[59]:


ml = os.listdir(fp+'data_other/MERRA2')


# In[60]:


ml = [ m for m in ml if m.endswith('nc4')] 
ml.sort()


# In[61]:


merra = []
for m in ml:
    print 'Loading file: {}'.format(m)
    mtmp, mdict = lu.load_hdf_h5py(fp+'data_other/MERRA2/'+m,all_values=True,verbose=False)
    merra.append(mtmp)


# In[62]:


len(merra),len(ml)


# In[201]:


mdict.keys()


# In[63]:


mdict['AODINC']


# In[64]:


merra[0]['time']


# In[65]:


merra[0]['time']/60.0/24.0


# In[62]:


(1-0.875)/2


# ### Collocate to 4STAR

# In[66]:


merra_doy = [datetime(int(l[-16:-12]),int(l[-12:-10]),int(l[-10:-8])).timetuple().tm_yday for l in ml]


# In[67]:


if False:
    for m in merra:
        m['aod'] = m['AODANA']+m['AODINC']
else:
    for m in merra:
        m['aod'] = m['AODANA']


# In[68]:


im4 = Sp.find_closest(np.array(merra_doy),(ar['doys']+0.0625).astype(int))


# In[69]:


np.unique(im4)


# In[70]:


nmerra = len(merra)
nlat = len(merra[0]['lat'])
nlon = len(merra[0]['lon'])
ntime = len(merra[0]['time'])


# In[71]:


merra2ar = {'doys':[],'lat':[],'lon':[],'aod':[],'ind':[]}
for ii in np.unique(im4):
    print ii
    imeas = np.where(im4==ii)[0]
    itime = Sp.find_closest(merra[ii]['time']/60.0/24.0+merra_doy[ii],ar['doys'][imeas])
    ilat = Sp.find_closest(merra[0]['lat'],ar['Latitude'][imeas])
    ilon = Sp.find_closest(merra[0]['lon'],ar['Longitude'][imeas])
    
    if len(imeas)<1:
        merra2ar['lat'] = np.append(merra2ar['lat'],ar['Latitude'][imeas])
        merra2ar['lon'] = np.append(merra2ar['lon'],ar['Longitude'][imeas])
        merra2ar['aod'] = np.append(merra2ar['aod'],ar['Latitude'][imeas]+np.nan)
        merra2ar['doys'] = np.append(merra2ar['doys'],ar['doys'][imeas])
        merra2ar['ind'] = np.append(merra2ar['ind'],ar['doys'][imeas].astype(int)*0) #pixel index
    else:
        merra2ar['lat'] = np.append(merra2ar['lat'],merra[ii]['lat'][ilat])
        merra2ar['lon'] = np.append(merra2ar['lon'],merra[ii]['lon'][ilon])
        merra2ar['aod'] = np.append(merra2ar['aod'],merra[ii]['aod'][itime,ilat,ilon])
        merra2ar['doys'] = np.append(merra2ar['doys'],merra[ii]['time'][itime]/60.0/24.0+merra_doy[ii])
        merra2ar['ind'] = np.append(merra2ar['ind'],
                                    np.ravel_multi_index((imeas*0+ii,itime,ilat,ilon),(nmerra,ntime,nlat,nlon)).astype(int))

for k in merra2ar.keys():
    merra2ar[k] = np.array(merra2ar[k])


# In[72]:


merra2ar['aod'].shape


# In[73]:


merra2ar['ind']


# In[71]:


plt.figure()
plt.plot(merra2ar['doys'],merra2ar['aod'],'.')
plt.plot(ar['doys'],ar['AOD0501'],'.')


# ### New Version with aerosol diagnostic

# From: M2T1NXAER: MERRA-2 tavg1_2d_aer_Nx: 2d,1-Hourly,Time-averaged,Single-Level,Assimilation,Aerosol Diagnostics V5.12.4
# 
# Downloaded at: https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXAER.5.12.4/2016/05/
# 
# DOI:10.5067/KLICLTZ8EM9D
# 
# 
# 
# Global Modeling and Assimilation Office (GMAO) (2015), MERRA-2 tavg1_2d_aer_Nx: 2d,1-Hourly,Time-averaged,Single-Level,Assimilation,Aerosol Diagnostics V5.12.4, Greenbelt, MD, USA, Goddard Earth Sciences Data and Information Services Center (GES DISC), Accessed: 2020-09-17, 10.5067/KLICLTZ8EM9D

# In[47]:


ml = os.listdir(fp+'data_other/MERRA2/aer')


# In[48]:


ml = [ m for m in ml if m.endswith('nc4')] 
ml.sort()


# In[49]:


mtmp, mdict = lu.load_hdf_h5py(fp+'data_other/MERRA2/aer/'+ml[0])


# In[50]:


merra = []
for m in ml:
    print 'Loading file: {}'.format(m)
    mtmp, mdict = lu.load_hdf_h5py(fp+'data_other/MERRA2/aer/'+m,
                                   values=(('lat',50),('lon',51),('time',52),('TOTEXTTAU',48),('TOTSCATAU',49),('TOTANGSTR',47)),
                                   verbose=False)
    merra.append(mtmp)


# In[51]:


len(merra),len(ml)


# In[52]:


mdict.keys()


# In[53]:


mdict['TOTANGSTR']


# In[54]:


mdict['TOTEXTTAU']


# In[72]:


merra[0]['time']/60.0/24.0


# In[73]:


merra[0]['time']


# In[55]:


delta_time = (1.0-merra[0]['time'][-1]/60.0/24.0)/2.0


# ### Collocate new version to 4STAR

# In[56]:


merra_doy = [datetime(int(l[-12:-8]),int(l[-8:-6]),int(l[-6:-4])).timetuple().tm_yday for l in ml]


# In[57]:


for m in merra:
    m['aod'] = m['TOTEXTTAU']
    m['ae'] = m['TOTANGSTR']


# In[58]:


im4 = Sp.find_closest(np.array(merra_doy),(ar['doys']+delta_time).astype(int))


# In[59]:


np.unique(im4)


# In[60]:


nmerra = len(merra)
nlat = len(merra[0]['lat'])
nlon = len(merra[0]['lon'])
ntime = len(merra[0]['time'])


# In[61]:


merra2ar = {'doys':[],'lat':[],'lon':[],'aod':[],'ae':[],'ind':[]}
for ii in np.unique(im4):
    print ii
    imeas = np.where(im4==ii)[0]
    itime = Sp.find_closest(merra[ii]['time']/60.0/24.0+merra_doy[ii],ar['doys'][imeas])
    ilat = Sp.find_closest(merra[0]['lat'],ar['Latitude'][imeas])
    ilon = Sp.find_closest(merra[0]['lon'],ar['Longitude'][imeas])
    
    if len(imeas)<1:
        merra2ar['lat'] = np.append(merra2ar['lat'],ar['Latitude'][imeas])
        merra2ar['lon'] = np.append(merra2ar['lon'],ar['Longitude'][imeas])
        merra2ar['aod'] = np.append(merra2ar['aod'],ar['Latitude'][imeas]+np.nan)
        merra2ar['ae'] = np.append(merra2ar['ae'],ar['Latitude'][imeas]+np.nan)
        merra2ar['doys'] = np.append(merra2ar['doys'],ar['doys'][imeas])
        merra2ar['ind'] = np.append(merra2ar['ind'],ar['doys'][imeas].astype(int)*0) #pixel index
    else:
        merra2ar['lat'] = np.append(merra2ar['lat'],merra[ii]['lat'][ilat])
        merra2ar['lon'] = np.append(merra2ar['lon'],merra[ii]['lon'][ilon])
        merra2ar['aod'] = np.append(merra2ar['aod'],merra[ii]['aod'][itime,ilat,ilon])
        merra2ar['ae'] = np.append(merra2ar['ae'],merra[ii]['ae'][itime,ilat,ilon])
        merra2ar['doys'] = np.append(merra2ar['doys'],merra[ii]['time'][itime]/60.0/24.0+merra_doy[ii])
        merra2ar['ind'] = np.append(merra2ar['ind'],
                                    np.ravel_multi_index((imeas*0+ii,itime,ilat,ilon),(nmerra,ntime,nlat,nlon)).astype(int))

for k in merra2ar.keys():
    merra2ar[k] = np.array(merra2ar[k])


# In[62]:


plt.figure()
plt.plot(merra2ar['doys'],merra2ar['aod'],'.')
plt.plot(ar['doys'],ar['AOD0501'],'.')


# ## Load in situ extinction

# In[63]:


lrgl = os.listdir(fp+'data_other/LARGE')


# In[64]:


lrgl.sort()


# In[65]:


lrgl


# In[66]:


large = []
for g in lrgl:
    print 'Loading file: {}'.format(g)
    gtmp, gdict = lu.load_ict(fp+'data_other/LARGE/'+g,return_header=True)
    large.append(gtmp)


# ### Match to 4STAR

# In[67]:


large[2]['UTC_mid'],large[3]['UTC_mid'] 


# In[68]:


ar['doys']


# In[69]:


large_doys = [datetime(int(l[26:30]),int(l[30:32]),int(l[32:34])).timetuple().tm_yday+np.nanmean(large[i]['UTC_mid'])/24.0 for i,l in enumerate(lrgl)]


# In[70]:


large_doys


# In[71]:


il4 = Sp.find_closest(np.array(large_doys),ar['doys'])


# In[72]:


lg2ar = {'ext_532':[],'SSA_550':[],'scat_450':[],'scat_550':[],'scat_700':[],'AE':[]}
for lll in np.unique(il4):
    imeas = np.where(il4==lll)[0]
    if len(imeas) < 1:
        lg2ar['ext_532'] = np.append(lg2ar['ext_532'],ar['aod'][imeas]+np.nan)
        lg2ar['SSA_550'] = np.append(lg2ar['SSA_550'],ar['aod'][imeas]+np.nan)
        lg2ar['scat_450'] = np.append(lg2ar['scat_450'],ar['aod'][imeas]+np.nan)
        lg2ar['scat_550'] = np.append(lg2ar['scat_550'],ar['aod'][imeas]+np.nan)
        lg2ar['scat_700'] = np.append(lg2ar['scat_700'],ar['aod'][imeas]+np.nan)
        lg2ar['AE'] = np.append(lg2ar['AE'],ar['aod'][imeas]+np.nan)
        
    else: 
        lg2ar['ext_532'] = np.append(lg2ar['ext_532'],
                                     wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambEXT_532_stdPT'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['SSA_550'] = np.append(lg2ar['SSA_550'], 
                                     wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambSSA_550'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['scat_450'] = np.append(lg2ar['scat_450'], 
                                      wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambSC450_stdPT'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['scat_550'] = np.append(lg2ar['scat_550'], 
                                      wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambSC550_stdPT'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['scat_700'] = np.append(lg2ar['scat_700'], 
                                      wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['ambSC700_stdPT'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
        lg2ar['AE'] = np.append(lg2ar['AE'], 
                                wu.nearest_neighbor(large[lll]['UTC_mid'],large[lll]['dryAEscat_550to700'],ar['Start_UTC'][imeas],dist=1.0/60.0/24.0))
  
    


# ## Load the MODIS reflectances

# In[73]:


fp


# In[74]:


fpd = fp+'data_other/MOD021KM.A2016140.0140.061.2017326011557.hdf'


# In[75]:


fpd


# In[9]:


mod,mod_dict = lu.load_hdf(fpd,values=(('refsb',0),('urefsb',1),('emis',2),('uemis',3),('refsb25',4),('refsw',7)))


# In[10]:


mod1,mod1_dict = lu.load_hdf(fp+'data_other/MOD02QKM.A2016140.0140.061.2017326011557.hdf',values=(('lat',2),('lon',3)))


# In[11]:


moda,moda_dict = lu.load_hdf(fp+'data_other/MOD04_L2.A2016140.0140.061.2017326145822.hdf',values=(('aod',135),('lat',71),('lon',70)))


# In[27]:


moda_dict['aod']


# In[12]:


modc,modc_dict = lu.load_hdf(fp+'data_other/MOD06_L2.A2016140.0140.061.2017326153640.hdf',values=(('cod',72),('cir',106)))


# In[11]:


mod1_dict['lat'].keys()


# In[13]:


refl_off = 316.9721985
refl_scal5 = 2.380933438e-06
refl_scal10 = 3.092909765e-06
refl_scal2 = 4.241572969e-05


# In[14]:


ref660 = (mod['refsb'][5,:,:]+refl_off)*refl_scal5
ref860 = (mod['refsb'][10,:,:]+refl_off)*refl_scal10
ref1240 = (mod['refsw'][2,:,:])*refl_scal2


# In[15]:


mod['refsw'].shape


# In[16]:


mod_dict['refsb']['band_names']


# ### Plot some MODIS data

# In[18]:


def make_map1(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=128,lat_0=36.0,
            llcrnrlon=125.0, llcrnrlat=33.0,
            urcrnrlon=130.0, urcrnrlat=38,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(123,133,11),labels=[0,0,0,1])
    m.drawparallels(np.linspace(31,39,17),labels=[1,0,0,0])
    return m


# In[18]:


fig,ax = plt.subplots(1,3)
m1 = make_map1(ax[2])
m1.pcolor(mod1['lon'],mod1['lat'],mod['refsw'][2,:,:],latlon=True)
ax[2].set_title('1240 nm')


m2 = make_map1(ax[0])
m2.pcolor(mod1['lon'],mod1['lat'],mod['refsb'][5,:,:],latlon=True)
ax[0].set_title('660 nm')

m3 = make_map1(ax[1])
m3.pcolor(mod1['lon'],mod1['lat'],mod['refsb'][10,:,:],latlon=True)
ax[1].set_title('860 nm')


# In[33]:


fig,ax = plt.subplots(1,1)
m0 = make_map1(ax)
mta = m0.pcolor(moda['lon'],moda['lat'],moda['aod'],latlon=True,cmap=plt.cm.YlOrRd,vmin=0,vmax=0.8)
mtc = m0.pcolor(mod1['lon'],mod1['lat'],modc['cod'],latlon=True,cmap=plt.cm.gist_earth,vmin=0,vmax=20.0)



cbxr = plt.gcf().add_axes([0.83, 0.2, 0.02, 0.6])
cbxe = plt.gcf().add_axes([0.80, 0.2, 0.02, 0.6])
#cbg.set_ticks([0,6,12,16,18,20])
#cbb.set_ticks([0,6,12,16,18,20]),cbb.set_ticklabels(['','','','',''])
cbxe.xaxis.set_ticks_position('top')#,cbaxesbl.yaxis.set_ticks_position('left')
#cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

bxr = plt.colorbar(mta,label='AOD',cax=cbxr,orientation='vertical',shrink=0.5)
bxe = plt.colorbar(mtc,cax=cbxe,orientation='vertical',shrink=0.5)
cbxe.yaxis.set_ticks_position('left')
#cbxe.set_label(' ')
cbxe.text(-4.0,0.5,'COD',rotation=90,verticalalignment='center')
plt.savefig(fp+'plot/KORUS_MODIS_AOD_COD_20160519.png',dpi=600,transparent=True)


# In[ ]:


fig,ax = plt.subplots(1,4)
m0 = make_map1(ax[0])
mta = m0.pcolor(moda['lon'],moda['lat'],moda['aod'],latlon=True,cmap=plt.cm.plasma)
mtc = m0.pcolor(mod1['lon'],mod1['lat'],modc['cod'],latlon=True,cmap=plt.cm.gist_earth)

m1 = make_map1(ax[3])
mp = m1.pcolor(mod1['lon'],mod1['lat'],ref1240,latlon=True,vmin=0,vmax=0.5)
ax[2].set_title('1240 nm')

m2 = make_map1(ax[1])
m2.pcolor(mod1['lon'],mod1['lat'],ref660,latlon=True,vmin=0,vmax=0.5)
ax[0].set_title('660 nm')

m3 = make_map1(ax[2])
m3.pcolor(mod1['lon'],mod1['lat'],ref860,latlon=True,vmin=0,vmax=0.5)
ax[1].set_title('860 nm')
plt.colorbar(mp,label='Reflectance',ax=ax[1:],shrink=0.6,location='bottom')
#plt.tight_layout()
#plt.savefig(fp+'plot/KORUS_MODIS_660_860_1240_reflect_20160519.png',dpi=600,transparent=True)


# In[38]:


fig,ax = plt.subplots(1,1)
m2 = make_map1(ax)
mp = m2.pcolor(mod1['lon'],mod1['lat'],ref660/ref860,latlon=True,vmin=0,vmax=4)
ax.set_title('660 nm / 860 nm')
plt.colorbar(mp)#label='Radiance ratio 660 nm / 860 nm')


#m3 = make_map1(ax[1])
#m3.pcolor(mod1['lon'],mod1['lat'],mod['refsb'][10,:,:],latlon=True)
#ax[1].set_title('860 nm')


# In[39]:


ref660/ref860


# In[23]:


mod['refsb'][5,:,:]/mod['refsb'][10,:,:]


# ## Load OCO-2 data

# In[59]:


oco, oco_dict = lu.load_hdf(fp+'data_other/acos_L2s_160511_06_B9200_PolB_190817214259.h5')


# In[60]:


lu.load_hdf_h5py(fp+'data_other/acos_L2s_160511_06_B9200_PolB_190817214259.h5')


# In[65]:


import h5py


# In[66]:


foco = h5py.File(fp+'data_other/acos_L2s_160511_06_B9200_PolB_190817214259.h5','r')


# In[81]:


ko = foco.keys()


# In[85]:


for k in ko:
    jk = foco[k].keys()
    print [(k,i) for i in jk if 'lat' in i]
    


# In[68]:


foco['AlbedoResults'].keys()


# In[86]:


foco['Shapes'].keys()


# In[114]:


foco['Shapes']['InputPtr_Array'].attrs.items()


# In[97]:


ra = foco['Shapes']['Retrieval_Array']


# In[102]:


ra.attrs.items()


# In[76]:


foco['AlbedoResults']['albedo_o2_fph'].attrs.items()


# In[120]:


aoo = foco['AlbedoResults']['albedo_o2_fph']


# In[122]:


aoo.value.shape


# In[125]:


foco['SoundingGeometry']['sounding_latitude'].value


# In[132]:


foco['RetrievalResults']['retrieved_o2_column'].value


# In[129]:


plt.figure()
plt.scatter(foco['SoundingGeometry']['sounding_longitude'].value,
          foco['SoundingGeometry']['sounding_latitude'].value,c=foco['AlbedoResults']['albedo_o2_fph'].value)


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# In[76]:


fl1 = ar['days']==ar['days'][0]


# In[77]:


fl1.shape


# In[78]:


fl = (ar['qual_flag']==0) & (np.isfinite(ar['AOD0501'])) 


# In[79]:


fl1.shape


# ## Calculate the Angstrom Exponent

# In[80]:


nwl,nm


# In[81]:


aodrr = np.array([ar[n] for n in nwl])


# In[82]:


aodrr.shape


# In[83]:


angs = su.calc_angs(ar['Start_UTC'],np.array(nm[1:11]),aodrr[1:11,:])


# ## Calculate the fine mode fraction

# In[84]:


fmf = su.sda(aodrr[1:13,:],np.array(nm[1:13])/1000.0)


# In[85]:


fmf.keys()


# In[86]:


fmf['tauc'].shape, ar['GPS_Alt'].shape


# ## Subset the level legs

# In[87]:


def running_std(x,n):
    'Function to do a running standard deviation on array (x) with window size (n)'
    q = x**2
    q = np.convolve(q, np.ones((n, )), mode="same")
    s = np.convolve(x, np.ones((n, )), mode="same")
    o = (q-s**2/n)/float(n-1)
    return o 


# In[88]:


nbox = 20


# In[89]:


std_alt = running_std(ar['GPS_Alt'][fl],nbox)


# In[90]:


std_alt.shape


# In[91]:


ar['GPS_Alt'][fl].shape


# In[92]:


f_level = np.where(std_alt<5.0)[0]


# In[93]:


std_alt1 = running_std(ar['GPS_Alt'][fl1],nbox)


# In[115]:


f_level1 = np.where(std_alt1<5.0)[0]


# In[116]:


ar['Start_UTC'][fl1][f_level1]


# In[117]:


plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(ar['Start_UTC'][fl1],ar['GPS_Alt'][fl1],'.')
plt.plot(ar['Start_UTC'][fl1][f_level1],ar['GPS_Alt'][fl1][f_level1],'r.')


ax2 = plt.subplot(2,1,2,sharex=ax1)
plt.plot(ar['Start_UTC'][fl1],std_alt1,'.')
plt.plot(ar['Start_UTC'][fl1][f_level1],std_alt[f_level1],'r.')
plt.ylim(0,100)


# In[118]:


plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(ar['Start_UTC'][fl],ar['GPS_Alt'][fl],'.')
plt.plot(ar['Start_UTC'][fl][f_level],ar['GPS_Alt'][fl][f_level],'r.')


ax2 = plt.subplot(2,1,2,sharex=ax1)
plt.plot(ar['Start_UTC'][fl],std_alt,'.')
plt.plot(ar['Start_UTC'][fl][f_level],std_alt[[f_level]],'r.')
plt.ylim(0,100)


# ## Seperate each of the level legs into distinct segments

# In[94]:


def get_segments(index,vals_dict,nsep=150,set_nan=True):
    'Function to seperate continuous segments (within a distance of nsep) based on a prior index'
    disc_flacaod_long = np.where(np.diff(index,1)>nsep)[0]
    
    discontinuity_istart_long =  index[np.append(0,disc_flacaod_long[:-1]+1)]
    discontinuity_iend_long =  index[disc_flacaod_long]
    
    kv = vals_dict.keys()
    d = {k:[] for k in kv}
    for i,start in enumerate(discontinuity_istart_long): # loop through discontinuities 
        if discontinuity_iend_long[i]-start < 2: continue
        for k in kv: # loop through keys
            try:
                d[k].append(vals_dict[k][start:discontinuity_iend_long[i]])
            except:
                print start, discontinuity_iend_long[i]
                continue
                #d[k].append([np.nan])
    
    for k in kv:
        d[k] = np.array(d[k])
        
    return d


# In[173]:


np.unique(ar['days'])


# In[174]:


np.where(np.diff(vals['doys'][f_level],1)>0.5)


# In[95]:


def get_segments_by_time(index,doys,vals_dict,tsep=5.0/24.0/60.0/60.0,set_nan=True):
    'Function to seperate continuous segments (within a distance in doys of tsep) based on a prior index (default for 5 seconds in doy)'
    disc_flacaod_long = np.where(np.diff(doys[index],1)>tsep)[0]
    
    discontinuity_istart_long =  index[np.append(0,disc_flacaod_long[:-1]+1)]
    discontinuity_iend_long =  index[disc_flacaod_long]
    
    kv = vals_dict.keys()
    d = {k:[] for k in kv}
    for i,start in enumerate(discontinuity_istart_long): # loop through discontinuities 
        if discontinuity_iend_long[i]-start < 2: continue
        for k in kv: # loop through keys
            try:
                d[k].append(vals_dict[k][start:discontinuity_iend_long[i]])
            except:
                print start, discontinuity_iend_long[i]
                continue
                #d[k].append([np.nan])
    
    for k in kv:
        d[k] = np.array(d[k])
        
    return d


# In[96]:


ar['days']


# In[97]:


angs[angs>5.0] = np.nan


# In[98]:


lg2ar.keys()


# In[99]:


vals = {'utc':ar['Start_UTC'][fl],'alt':ar['GPS_Alt'][fl],'lat':ar['Latitude'][fl],'lon':ar['Longitude'][fl],
        'aod0500':ar['AOD0501'][fl],'aod1040':ar['AOD1040'][fl],'AE':angs[fl],'doys':ar['doys'][fl],
        'aod_fine':fmf['tauf'][fl],'aod_coarse':fmf['tauc'][fl],'fmf':fmf['eta'][fl],
        'GOCI_AOD':goci2ar['aod'][fl],'GOCI_AE':goci2ar['AE'][fl],'GOCI_fmf':goci2ar['fmf'][fl],
        'MERRA_AOD':merra2ar['aod'][fl], 'MERRA_AE':merra2ar['ae'][fl],
        'situ_ext':lg2ar['ext_532'][fl],'situ_ssa':lg2ar['SSA_550'][fl],'situ_ae':lg2ar['AE'][fl]}


# In[100]:


dvals = get_segments(f_level,vals,nsep=100)


# In[101]:


dvalst = get_segments_by_time(f_level,vals['doys'],vals,tsep=200.0/24.0/60.0/60.0)


# In[125]:


dvals.keys()


# In[204]:


len(dvals['AE']),len(dvalst['AE'])


# In[205]:


dvals['doys'][0]


# In[206]:


dvalst['doys'][0]


# In[207]:


dvalst['doys'][1]


# In[186]:


dvals['len_minutes'] = []
dvals['len_doys'] = []
for i,n in enumerate(dvals['utc']):
    try:
        len_min = (n[-1]-n[0])*60.0
        if len_min < 0.0: 
            len_min = (n[-1]-(n[0]-24.0))*60.0
        print i,len_min
        dvals['len_minutes'].append(len_min)
        dvals['len_doys'].append(dvals['doys'][i][-1]-dvals['doys'][i][0])
    except:
        print np.nan
        dvals['len_minutes'].append(np.nan)
        dvals['len_doys'].append(np.nan)


# In[208]:


dvalst['len_minutes'] = []
dvalst['len_doys'] = []
for i,n in enumerate(dvalst['utc']):
    try:
        len_min = (n[-1]-n[0])*60.0
        print i,len_min
        dvalst['len_minutes'].append(len_min)
        dvalst['len_doys'].append(dvalst['doys'][i][-1]-dvalst['doys'][i][0])
    except:
        print np.nan
        dvalst['len_minutes'].append(np.nan)
        dvalst['len_doys'].append(np.nan)


# In[209]:


plt.figure()
plt.hist([dvals['len_doys'],dvalst['len_doys']],bins=50,label=['number based','time based'])
plt.legend()
plt.yscale('log')
plt.ylabel('Number of segments')
plt.xlabel('Fractional day length of segment')


# In[102]:


dvals = dvalst # set the default to use the time seperated segments


# ### discrete colorbar

# In[188]:


def discrete_matshow(data,cmapname='RdBu'):
    ' plotting function for a discrete colormap'
    cmap = plt.get_cmap(cmapname, np.nanmax(data)-np.nanmin(data)+1)
    # set limits .5 outside true range
    scalarmap = plt.cm.ScalarMappable(cmap=cmapname)
    scalarmap.set_array(data)
    #mat = plt.matshow(data,cmap=cmap,vmin = np.min(data)-.5, vmax = np.max(data)+.5)
    #tell the colorbar to tick at integers
    cax = plt.colorbar(scalarmap, ticks=np.arange(np.min(data),np.max(data)+1))
    return cax


# In[ ]:


for q in np.unique(ar['days']):
    flq = ar['days'][fl]==q
    flql = ar['days'][fl][f_level]==q
    plt.figure()
    plt.plot(ar['Start_UTC'][fl][flq],ar['GPS_Alt'][fl][flq],'.')
    plt.plot(ar['Start_UTC'][fl][f_level][flql],ar['GPS_Alt'][fl][f_level][flql],'r.')
    ax = plt.gca()

    ax.set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, len(dvals['utc'])+1)])

    for i,n in enumerate(dvals['utc']):
        plt.plot(n,dvals['alt'][i],'s-',markeredgecolor='None')

    plt.xlabel('UTC [h from midnight]')
    plt.ylabel('Altitude [m]')
    plt.title('Days: {}'.format(q))

    #scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('gist_ncar'))
    #scalarmap.set_array(range(len(dvals['utc'])+1))
    #cb = plt.colorbar(scalarmap)
    cb = discrete_matshow(range(len(dvals['utc'])+1),cmapname='gist_ncar')
    cb.set_label('Level leg number')
#plt.plot(ar['Start_UTC'][fl1][f_level],ar['GPS_Alt'][fl1][f_level],'r.')


# ## Now calculate the distances travelled within each segments

# In[103]:


def get_distances(seg_dict):
    'Function that calculates the cumulative distance and instantaneous change between each point for a set of segments'
    seg_dict['dist'],seg_dict['cumdist'] = [],[]
    for i,l in enumerate(seg_dict['lat']):
        try:
            ckm,km = [],[]
            pos1 = [seg_dict['lat'][i][0],seg_dict['lon'][i][0]] 
            for j,a in enumerate(seg_dict['lat'][i]):
                d = mu.spherical_dist(pos1,[seg_dict['lat'][i][j],seg_dict['lon'][i][j]])
                ckm.append(d)
                km.append(d)
        except:
            cckm,dkm = [np.nan],[np.nan]

        iu = np.where(np.isfinite(ckm))[0]
        try:
            fckm = interpolate.interp1d(seg_dict['utc'][i][iu],np.array(ckm)[iu])
            fkm = interpolate.interp1d(seg_dict['utc'][i][iu],np.array(km)[iu])
            cckm = fckm(seg_dict['utc'][i])
            dkm = fkm(seg_dict['utc'][i])
            seg_dict['cumdist'].append(np.array(cckm))
            seg_dict['dist'].append(np.array(np.diff(dkm)))
        except:
            seg_dict['cumdist'].append(np.array(np.nan))
            seg_dict['dist'].append(np.array(np.nan))

    return seg_dict


# In[104]:


ddv = get_distances(dvals)


# In[105]:


dvals['cumdist'][0:2]


# In[107]:


len(dvals['cumdist'])


# ## Save to file for easier reloading

# In[108]:


fp


# In[110]:


hs.savemat(fp+'KORUS_autocorr_dvals.mat',dvals)


# ## Load from file

# In[8]:


dvals = hs.loadmat(fp+'KORUS_autocorr_dvals.mat')


# # Calculate the autocorrelation of AOD with respect to distance

# **From Shinozuka and Redemann, 2011, Horizontal variability of aerosol optical depth observed during the ARCTAS airborne experiment, ACP**
# 
# Autocorrelation is the correlation coefficient among all
# data pairs xj and xj+k that exist at a separation, or lag, of k. That is,
# 
# ![image.png](attachment:image.png)
# 
# where k indicates the spatial lag (or distance), m+k and std+k denote the mean and standard deviation, respectively, of all data points that are located a distance of +k away from an- other data point, and m−k and std−k are the corresponding quantities for data points located a distance of −k away from another data point (Redemann et al., 2006; Anderson et al., 2003).
# Figure 1c shows pairs of 499nm AOD measured 20km (±0.2 km) away from each other in the Canada phase. The correlation coefficient, r, is 0.37. This is the autocorrelation for 20km.

# Define the different time periods from Met. From: 
#     Peterson, D.A., Hyer, E.J., Han, S.-O., Crawford, J.H., Park, R.J., Holz, R., Kuehn, R.E., Eloranta, E., Knote, C., Jordan, C.E. and Lefer, B.L., 2019. Meteorology influencing springtime air quality, pollution transport, and visibility in Korea. Elem Sci Anth, 7(1), p.57. DOI: http://doi.org/10.1525/elementa.395
# 
#     - Dynamic meteorology and complex aerosol vertical profiles (01–16 May);
#     - Stagnation under a persistent anticyclone (17–22 May);
#     - Dynamic meteorology, low-level transport, and haze development (25–31 May); (Extreme pollution)
#     - Blocking pattern (01–07 June).
# 
# 
# ![image.png](attachment:image.png)

# The distribution and source of pm 2.5 during different met periods
# 
# ![image.png](attachment:image.png)

# ### Build the limits of the autocorrelation

# In[9]:


# for met times
dt1 = ['20160501','20160516']
dt2 = ['20160517','20160522']
dt3 = ['20160523','20160531']
dt4 = ['20160601','20160607']


# In[10]:


t1 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt1]
t2 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt2]
t3 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt3]
t4 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt4]


# In[11]:


# limits of DOY for each of the met times
t1,t2,t3,t4


# In[12]:


#altitude limits in m
z1 = [0.0, 1000.0] 
z2 = [1000.0, 3000.0]
z3 = [3000.0, 15000.0]


# ### Test out Shinozuka & Redemann autocorrelation 

# In[13]:


dvals['cumdist'][2]


# In[14]:


dvals.keys()


# In[121]:


# method for making one autocorrelation per segment. Old not recommended
if False:
    cr = []
    for i,cd in enumerate(dvals['cumdist']):
        #cd = dvals['cumdist'][i]
        corr = {'aod1040':[],'aod0500':[],'AE':[]}
        corr_ks =[0.1,0.25,0.5,0.75,1.0,1.5,2.0,3.0,5.0,7.5,10.0,12.5,15.0,20.0,
                  25.0,30.0,35.0,40.0,50.0,60.0,75.0,100.0,150.0,200.0] 
        for ik, k in enumerate(corr_ks):

        #k = 5.0 # for 5km distance
            if k>np.nanmax(cd):
                [corr[val].append(np.nan) for val in corr.keys()]
                continue
            ipk = np.argmin(abs(cd-k)) #ipk:
            imk = np.argmin(abs(cd-(cd[-1]-k))) #0:imk
            N = len(cd)
            #c = np.sqrt(2.0/(N-1))*math.gamma(N/2.0)/math.gamma((N-1.0)/2.0)

            for val in dvals.keys():
                if val in ['lon','utc','lat','cumdist','cdist_n','dist','alt','autocor','aod1040_r','AE_r','aod_n']: continue
                #print val, len(dvals[val][i])
                mpk = np.nanmean(dvals[val][i][ipk:]) #mean +k
                mmk = np.nanmean(dvals[val][i][0:imk]) #mean -k
                spk = np.nanstd(dvals[val][i][ipk:]) #std +k
                smk = np.nanstd(dvals[val][i][0:imk]) #std -k
                top = [(dvals[val][i][j]-mpk)*(dvals[val][i][j+ipk]-mmk) for j in xrange(N-ipk-1)]
                #dvals[val+'_r'] = []
                corr[val].append(np.sum(top)/((N-1)*spk*smk))
                if (corr[val][-1]>1.0) | (corr[val][-1]<0.0):
                    print '{} has bad corr: {:2.2f} val for key {}: std+k:{:2.2f}, std-k:{:2.2f}, m+k:{:2.2f}, m-k:{:2.2f} '.format(i,
                        corr[val][-1],val,spk,smk,mpk,mmk)

        for val in corr.keys():
            corr[val] = np.array(corr[val])
        cr.append(corr)


# In[122]:


len(cr)


# ### Autocorraltion calculation with all vars

# In[15]:


min_num = 100 # minimum number of points to be considered valid


# In[16]:


types = ['all','t1','t2','t3','t4','z1','z2','z3']


# In[17]:


corr_ks =[0.08,0.1,0.25,0.5,0.75,1.0,1.5,2.0,3.0,5.0,7.5,10.0,12.5,15.0,20.0,
          25.0,30.0,35.0,40.0,50.0,60.0,75.0,100.0,150.0,200.0] 
corr_ks =[0.08,0.1,0.25,0.5,1.0,1.5,3.0,5.0,10.0,15.0,
          35.0,65.0,100.0,150.0,200.0] 
corr_all = [[[{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[],'fmf':[],'GOCI_AOD':[],'GOCI_AE':[],'GOCI_fmf':[],'MERRA_AOD':[],'MERRA_AE':[],'situ_ext':[],'situ_ae':[],'situ_ssa':[]} for i,k in enumerate(corr_ks)],
             [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[],'fmf':[],'GOCI_AOD':[],'GOCI_AE':[],'GOCI_fmf':[],'MERRA_AOD':[],'MERRA_AE':[],'situ_ext':[],'situ_ae':[],'situ_ssa':[]} for i,k in enumerate(corr_ks)]] \
            for j in types] 
corr_alln = [[[{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[],'fmf':[],'GOCI_AOD':[],'GOCI_AE':[],'GOCI_fmf':[],'MERRA_AOD':[],'MERRA_AE':[],'situ_ext':[],'situ_ae':[],'situ_ssa':[]} for i,k in enumerate(corr_ks)],
             [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[],'fmf':[],'GOCI_AOD':[],'GOCI_AE':[],'GOCI_fmf':[],'MERRA_AOD':[],'MERRA_AE':[],'situ_ext':[],'situ_ae':[],'situ_ssa':[]} for i,k in enumerate(corr_ks)]] \
            for j in types] 
#corr_all = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_t1 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_t2 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_t3 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_t4 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_z1 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_z2 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
#corr_z3 = [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]
corr_vals = corr_all[0][0][0].keys()
corr_vals.remove('k')
corr_valsn = corr_alln[0][0][0].keys()
corr_valsn.remove('k')


# In[181]:


corr_all[0][0][0],corr_all[0][0][0]


# In[120]:


np.array(corr_all).shape #type, [minusk,plusk], distance


# In[18]:


len(dvals['cumdist'])


# In[19]:


dv = 0.05
for i,cd in enumerate(dvals['cumdist']):
    mat_dist = np.array([abs(cd-d) for d in cd])
    sys.stdout.write( '*{}*'.format(i))
    for ik,k in enumerate(corr_ks):
        #sys.stdout.write("{} ".format(ik))
        iimg,iipg = np.where((mat_dist>(k*(1.0-dv)))&(mat_dist<(k*(1.0+dv))))
        if len(iimg)==0:
            continue
        N = len(cd)
        for val in corr_vals:
            #all
            corr_alln[0][0][ik][val].append(dvals[val][i][iimg])
            corr_alln[0][1][ik][val].append(dvals[val][i][iipg])
            #type 1
            if (dvals['doys'][i][0]> t1[0]) & (dvals['doys'][i][0]< t1[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[1][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[1][1][ik][val].append(dvals[val][i][iipg])
            #type 2
            if (dvals['doys'][i][0]> t2[0]) & (dvals['doys'][i][0]< t2[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[2][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[2][1][ik][val].append(dvals[val][i][iipg])
            #type 3
            if (dvals['doys'][i][0]> t3[0]) & (dvals['doys'][i][0]< t3[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[3][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[3][1][ik][val].append(dvals[val][i][iipg])
            #type 4
            if (dvals['doys'][i][0]> t4[0]) & (dvals['doys'][i][0]< t4[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[4][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[4][1][ik][val].append(dvals[val][i][iipg])
            #type 5
            if (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z1[1]):
                corr_alln[5][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[5][1][ik][val].append(dvals[val][i][iipg])
            #type 6
            if (dvals['alt'][i][0]> z2[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_alln[6][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[6][1][ik][val].append(dvals[val][i][iipg])
            #type 7
            if (dvals['alt'][i][0]> z3[0]) & (dvals['alt'][i][0]< z3[1]):
                corr_alln[7][0][ik][val].append(dvals[val][i][iimg])
                corr_alln[7][1][ik][val].append(dvals[val][i][iipg])


# In[316]:


if False:
    dv = 0.15 #30% leeway in values
    for ik, k in enumerate(corr_ks):
        for i,cd in enumerate(dvals['cumdist']):
            if k>np.nanmax(cd):
                continue
            if ik==0:
                ipk = 1
                imk = len(cd)-2
                iip = range(1,len(cd)-1)
                iipg = iip
                iimg = range(0,imk)
            else:
                ipk = np.argmin(abs(cd-k)) #ipk:
                imk = np.argmin(abs(cd-(cd[-1]-k))) #0:imk
                if imk<2:
                    continue
                iip = Sp.find_closest(cd,cd[0:imk]+k)
                iip_good = (abs(cd[iip]-cd[0:imk])<k*(1.0+dv)) & (abs(cd[iip]-cd[0:imk])>k*(1.0-dv))
                if not any(iip_good):
                    continue
                iipg = iip[iip_good]
                iimg = np.arange(0,imk)[iip_good]
            N = len(cd)
            for val in corr_vals:
                #all
                corr_all[0][0][ik][val] = np.append(corr_all[0][0][ik][val],dvals[val][i][iimg])
                corr_all[0][1][ik][val] = np.append(corr_all[0][1][ik][val],dvals[val][i][iipg])
                #type 1
                if (dvals['doys'][i][0]> t1[0]) & (dvals['doys'][i][0]< t1[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[1][0][ik][val] = np.append(corr_all[1][0][ik][val],dvals[val][i][iimg])
                    corr_all[1][1][ik][val] = np.append(corr_all[1][1][ik][val],dvals[val][i][iipg])
                #type 2
                if (dvals['doys'][i][0]> t2[0]) & (dvals['doys'][i][0]< t2[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[2][0][ik][val] = np.append(corr_all[2][0][ik][val],dvals[val][i][iimg])
                    corr_all[2][1][ik][val] = np.append(corr_all[2][1][ik][val],dvals[val][i][iipg])
                #type 3
                if (dvals['doys'][i][0]> t3[0]) & (dvals['doys'][i][0]< t3[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[3][0][ik][val] = np.append(corr_all[3][0][ik][val],dvals[val][i][iimg])
                    corr_all[3][1][ik][val] = np.append(corr_all[3][1][ik][val],dvals[val][i][iipg])
                #type 4
                if (dvals['doys'][i][0]> t4[0]) & (dvals['doys'][i][0]< t4[1]) & (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[4][0][ik][val] = np.append(corr_all[4][0][ik][val],dvals[val][i][iimg])
                    corr_all[4][1][ik][val] = np.append(corr_all[4][1][ik][val],dvals[val][i][iipg])
                #type 5
                if (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z1[1]):
                    corr_all[5][0][ik][val] = np.append(corr_all[5][0][ik][val],dvals[val][i][iimg])
                    corr_all[5][1][ik][val] = np.append(corr_all[5][1][ik][val],dvals[val][i][iipg])
                #type 6
                if (dvals['alt'][i][0]> z2[0]) & (dvals['alt'][i][0]< z2[1]):
                    corr_all[6][0][ik][val] = np.append(corr_all[6][0][ik][val],dvals[val][i][iimg])
                    corr_all[6][1][ik][val] = np.append(corr_all[6][1][ik][val],dvals[val][i][iipg])
                #type 7
                if (dvals['alt'][i][0]> z3[0]) & (dvals['alt'][i][0]< z3[1]):
                    corr_all[7][0][ik][val] = np.append(corr_all[7][0][ik][val],dvals[val][i][iimg])
                    corr_all[7][1][ik][val] = np.append(corr_all[7][1][ik][val],dvals[val][i][iipg])


# In[21]:


len(corr_alln[0][0][0]['AE'])


# In[22]:


corr_all[0][0][0]['AE'] = np.hstack(corr_alln[0][0][0]['AE'])


# In[25]:


for j,jt in enumerate(types):
    corr_all[j][1][0]['AE'] = np.hstack(corr_alln[j][1][0]['AE'])
    print corr_all[j][1][0]['AE'].shape


# In[23]:


len(corr_all)


# In[ ]:


autocorr = {}
autocorr_len = {}

for val in corr_vals:
    autocorr[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_len[val] = np.zeros((len(types),len(corr_ks)))
    print val
    for ik,k in enumerate(corr_ks):
        sys.stdout.write("{} ".format(ik))
        for j,jt in enumerate(types):
            corr_all[j][0][ik][val] = np.hstack(corr_alln[j][0][ik][val])
            corr_all[j][1][ik][val] = np.hstack(corr_alln[j][1][ik][val])
            if False: #val is 'AE':
                igd = np.where(corr_all[2][1][3]['AE']<5.0)[0] #np.isfinite(corr_all[j][0][ik][val])
                mmk = np.nanmean(corr_all[j][0][ik][val][igd])
                mpk = np.nanmean(corr_all[j][1][ik][val][igd])
                smk = np.nanstd(corr_all[j][0][ik][val][igd])
                spk = np.nanstd(corr_all[j][1][ik][val][igd])
                top = [(v-mpk)*(corr_all[j][1][ik][val][igd][iv]-mmk) for iv,v in enumerate(corr_all[j][0][ik][val][igd])]
                #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
                autocorr[val][j,ik] = np.nansum(top)/((len(corr_all[j][0][ik][val][igd])-1)*spk*smk)
            else:
                mmk = np.nanmean(corr_all[j][0][ik][val])
                mpk = np.nanmean(corr_all[j][1][ik][val])
                smk = np.nanstd(corr_all[j][0][ik][val])
                spk = np.nanstd(corr_all[j][1][ik][val])
                top = [(v-mpk)*(corr_all[j][1][ik][val][iv]-mmk) for iv,v in enumerate(corr_all[j][0][ik][val])]
                #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
                autocorr[val][j,ik] = np.nansum(top)/((len(corr_all[j][0][ik][val])-1)*spk*smk)
                autocorr_len[val][j,ik] = len(corr_all[j][0][ik][val])
                if autocorr_len[val][j,ik]<min_num:
                    autocorr[val][j,ik] = np.nan


# In[170]:


np.hstack(corr_all[j+1][0][ik][val])


# In[168]:


np.concatenate(corr_all[j+1][0][ik][val])


# In[ ]:





# In[318]:


autocorr['aod0500'].shape


# In[319]:


autocorr_len['aod0500'].shape


# ### Run monte Carlo sub-sampling to calculate variation in autocorrelation

# In[320]:


def fx_autocor(plus_dist,minus_dist,min_num=100):
    'Function to calculate the autocorrelation from 2 populations of the same distribution - plus and minus'
    # minus_dist = corr_all[j][0][ik][val]
    # plus_dist = corr_all[j][1][ik][val]
    mmk = np.nanmean(minus_dist)
    mpk = np.nanmean(plus_dist)
    smk = np.nanstd(minus_dist)
    spk = np.nanstd(plus_dist)
    top = [(v-mpk)*(plus_dist[iv]-mmk) for iv,v in enumerate(minus_dist)]
    #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
    autocorr = np.nansum(top)/((len(minus_dist)-1)*spk*smk)
    autocorr_len = len(minus_dist)
    if autocorr_len<min_num:
        autocorr = np.nan
    return autocorr, autocorr_len


# In[321]:


def rand_autocor(plus_dist,minus_dist,min_num=100,subsamp_ratio=0.5):
    'fx_autocor, but called with random subsampling at a ratio defined by subsamp_ratio [0-1]'
    irand = np.random.randint(len(plus_dist),size=int(subsamp_ratio*len(plus_dist)))
    return fx_autocor(plus_dist[irand],minus_dist[irand],min_num=min_num)


# In[322]:


j,ik,val


# In[323]:


len(corr_all[j][0][ik][val])


# In[324]:


autocorr[val][j,ik],autocorr_len[val][j,ik]


# In[325]:


fx_autocor(corr_all[j][0][ik][val],corr_all[j][1][ik][val])


# In[326]:


rand_autocor(corr_all[j][0][ik][val],corr_all[j][1][ik][val],subsamp_ratio=0.5)


# In[327]:


autocorr_std = {}
autocorrs = {}
num_for_std = 50

for val in corr_vals:
    autocorrs[val] = np.zeros((len(types),len(corr_ks),num_for_std))+np.nan
    autocorr_std[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    print val
    
    for ik,k in enumerate(corr_ks):
        for j,jt in enumerate(types):
            if not len(corr_all[j][0][ik][val]):
                continue
            for u in xrange(num_for_std):
                autocorrs[val][j,ik,u],l = rand_autocor(corr_all[j][0][ik][val],corr_all[j][1][ik][val],subsamp_ratio=0.5)
            autocorr_std[val][j,ik] = np.nanstd(autocorrs[val][j,ik,:])


# In[352]:


autocorr_min = {}
autocorr_max = {}
autocorr_d = {}
autocorr_mean = {}

for val in corr_vals:
    autocorr_min[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_max[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_d[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_mean[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    
    for ik,k in enumerate(corr_ks):
        for j,jt in enumerate(types):
            if not len(corr_all[j][0][ik][val]):
                continue
            autocorr_min[val][j,ik] = np.nanmin(autocorrs[val][j,ik,:])
            autocorr_max[val][j,ik] = np.nanmax(autocorrs[val][j,ik,:])
            autocorr_d[val][j,ik] = autocorr_max[val][j,ik] - autocorr_min[val][j,ik]
            autocorr_mean[val][j,ik] = np.nanmean(autocorrs[val][j,ik,:])


# ### Integrated autocorrelation

# In[193]:


def autocorr(x, t=1):
    return np.corrcoef(np.array([x[:-t], x[t:]]))


# In[194]:


def autocorr2(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]/result.max()


# In[195]:


def autocorr5(x):
    '''numpy.correlate, non partial'''
    n = len(x)
    lags = range(n)
    #mean=x.mean()
    var=np.nanvar(x)
    xp=x-np.nanmean(x)
    corr=np.correlate(xp,xp,'full')[n-1:]/var/n

    return corr[:n]


# In[196]:


authcor = autocorr(dvals['aod0500'][1])
authcor2 = autocorr2(dvals['aod0500'][1])
authcor3 = autocorr5(dvals['aod0500'][1])


# In[197]:


len(authcor2)


# In[237]:


len(dvals['aod0500'][1])


# In[238]:


dvals['dist'][1]


# In[ ]:


[(dvals['dist'][i].mean(),np.size(dvals['dist'][i])) for i in xrange(len(dvals['dist']))]


# ### interpolate AODs to a constant distance grid

# In[204]:


def interp_dist(d,dist=0.12,verbose=False):
    'function to insterpolate the AOD from the dict to an even grid spacing accroding to distance (default 0.12 km)'
    d['cdist_n'],d['aod_n'] = [],[]
    for i,cd in enumerate(d['cumdist']):
        if verbose:
            print i, cd.min(),cd.max(), np.nanmin(cd),np.nanmax(cd)
            if not np.isfinite(cd.min()): print cd
        d['cdist_n'].append(np.arange(cd.min(),cd.max(),dist))
        try:
            fcd = interpolate.interp1d(cd,d['aod0500'][i])
            d['aod_n'].append(fcd(d['cdist_n'][i]))
        except TypeError:
            d['aod_n'].append(np.array(np.nan))


# In[ ]:


dvals['aod0500']


# In[205]:


interp_dist(dvals)


# In[206]:


dvals['autocor'] = [] 
for i,a in enumerate(dvals['aod_n']):
    try:
        dvals['autocor'].append(autocorr5(a))
    except:
        dvals['autocor'].append(np.array(np.nan))


# ### Autocorrelation plots (old)

# In[207]:


plt.figure()
plt.plot(dvals['cdist_n'][1],dvals['autocor'][1],'.')
plt.xlabel('Lag Distance [km]')
plt.ylabel('Correlation')


# In[208]:


plt.figure()
for i,j in enumerate(dvals['cdist_n']):
    try:
        plt.plot(j,dvals['aod_n'][i],'.')
    except:
        pass


# In[209]:


plt.figure()
for i,j in enumerate(dvals['cdist_n']):
    try:
        plt.plot(j,dvals['autocor'][i],'.')
    except:
        pass
plt.ylabel('Correlation Coefficient')
plt.xlabel('Distance [km]')
plt.xscale('log')


# ### Autocorrelation plot with Anderson method

# In[239]:


key_list = ['aod0500','aod1040','AE','aod_fine','aod_coarse']
key_list2 = ['aod0500','aod1040','AE','fmf']
legend_list = ['All','Dynamic','Stagnation','Extreme pollution','Blocking','0-1 km','1-3 km','3+ km']
cl_list = ['k','tab:red','tab:blue','tab:orange','tab:green','tab:olive','tab:cyan','tab:purple']
m_list = ['.','o','s','v','^','*','+','x']
tit = ['AOD$_{{500}}$','AOD$_{{1040}}$','AE','AOD$_{{fine}}$','AOD$_{{coarse}}$']
tit2 = ['AOD$_{{500}}$','AOD$_{{1040}}$','AE','fine-mode fraction']


# In[705]:


fig, ax = plt.subplots(5,3,figsize=(12,9))
for i,k in enumerate(key_list):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks,autocorr[k][j,:],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks,autocorr[k][j,:],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
    if i==1:
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.4])
    if i==4:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_all.png',dpi=600,transparent=True)


# In[240]:


key_list


# In[329]:


fig, ax = plt.subplots(5,3,figsize=(12,9))
for i,k in enumerate(key_list):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
    if i==1:
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.4])
    if i==4:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)


# ## Load the Autocorrelations from Shinozuka & Redemann

# In[241]:


SR_corr_ks = [0.45,1.0,3.0,6.0,10.0,20.0,34.2]
SR_aod_corr_long = [0.998,0.997,0.995,0.985,0.981,0.946,0.917]
SR_aod_corr_loc = [0.975,0.941,0.830,0.712,0.584,0.365,0.328]
SR_AE_corr_long = [1.000,0.977,0.975,0.960,0.944,0.913,0.919]
SR_AE_corr_loc = [0.975,0.956,0.919,0.831,0.747,0.519,0.366]


# In[151]:


#AOD local

0,4376084276355722; 0,9757725145572111
0,9955912496957824; 0,9408298503197231
2,993269282069332; 0,8303668774336445
5,997701705897239; 0,711850105383489
9,981245014250291; 0,5844939806380168
20,005524337736603; 0,36462079805665715
35,166804844012965; 0,32825206301575405


# In[152]:


#AE local
0,43836153494904556; 0,9757721573250459
0,9972623352215253; 0,9556260493694856
2,9925078192756667; 0,9191462151252102
5,985354161428386; 0,8309627406851716
9,993780925398312; 0,7465159146929592
19,996703813167056; 0,5185049833887045
35,22351663341235; 0,36598292430250434


# In[153]:


#AOD long
0,43833272405336876; 0,9987068195620337
0,9971397077970326; 0,9985360625870756
2,9918543848699146; 0,9953484799771375
6,0033122296001995; 0,985586039366985
9,969908496905848; 0,9810416889936774
19,972213345966814; 0,9461254599364131
35,10743016861416; 0,9178948308505701


# In[154]:


#AE long
0,43757781072798063; 1,0001868324223917
0,9989150464154728; 0,9778205265602119
2,992025659700471; 0,9753731289965352
5,9831337971357215; 0,9604326081520385
9,970944322675367; 0,9447901261029547
19,97407661698948; 0,9135730361161722
35,10728129935923; 0,9193744864787629


# In[242]:


note = [['a)','b)'],['c)','d)'],['e)','f)'],['g)','h)'],['i)','j)']]


# In[718]:


fig, ax = plt.subplots(5,3,figsize=(12,9))
for i,k in enumerate(key_list):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks,autocorr[k][j,:],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks,autocorr[k][j,:],label=legend_list[j],color=cl_list[j],marker=m_list[j])
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.1])
    
    if i==3:
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==4:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_all_with_SR2011.png',dpi=600,transparent=True)


# In[244]:


key_list


# In[245]:


autocorr.keys()


# In[1003]:


fig, ax = plt.subplots(4,3,figsize=(12,9))
for i,k in enumerate(key_list2):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'AE':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'AE':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit2[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.1])
    
    if i==3:
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==3:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_rel_all_with_SR2011_with_GOCI.png',dpi=600,transparent=True)


# ### Add MERRA

# In[1237]:


fig, ax = plt.subplots(4,3,figsize=(12,9))
for i,k in enumerate(key_list2):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
        if k is 'AE':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
        if k is 'AE':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit2[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.1])
    
    if i==3:
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==3:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_rel_all_with_SR2011_GOCI_MERRA.png',dpi=600,transparent=True)


# ### Add In situ

# In[227]:


fig, ax = plt.subplots(4,3,figsize=(12,9))
for i,k in enumerate(key_list2):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
                        color=cl_list[j],ls='-.',lw=1)
        if k is 'AE':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'AE':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit2[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].plot([],[],'-.',c='grey',lw=1,label='In situ Ext.')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.1,0.1])
    
    if i==3:
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==3:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_rel_all_with_SR2011_GOCI_MERRA_insitu.png',dpi=600,transparent=True)


# In[228]:


fig, ax = plt.subplots(4,3,figsize=(12,9))
for i,k in enumerate(key_list2):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
                        color=cl_list[j],ls='-.',lw=1)
        if k is 'AE':
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],marker=m_list[j])
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'AE':
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
                         color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr[k][0,1:]/autocorr[k][0,1] - autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(-1,1)
    ax[i,1].set_ylim(-1,1)
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('diff r({})'.format(tit2[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        auk_norm = wu.nearest_neighbor(np.array(corr_ks[1:]),autocorr[k][0,1:]/autocorr[k][0,1],np.array(SR_corr_ks),dist=30.0)
        
        ax[i,0].plot(SR_corr_ks,auk_norm-SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,auk_norm-SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==1:
        ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].plot([],[],'-.',c='grey',lw=1,label='In situ Ext.')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.9])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[1.6,0.1])
    
    if i==3:
        auk_norm = wu.nearest_neighbor(np.array(corr_ks[1:]),autocorr[k][0,1:]/autocorr[k][0,1],np.array(SR_corr_ks),dist=30.0)
        ax[i,0].plot(SR_corr_ks,auk_norm-SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,auk_norm-SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==3:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.savefig(fp+'plot/KORUS_Autocorr_diff_rel_all_with_SR2011_GOCI_MERRA_insitu.png',dpi=600,transparent=True)


# ## Only plot subset of autocorrelation extensive vs intensive

# In[246]:


autocorr_len.keys()


# In[351]:


fig, ax = plt.subplots(2,2,figsize=(6,4))
i = 0
for j in [0,1,2,3,4]:
    ax[i,0].plot(corr_ks[1:],autocorr_len['aod0500'][j,1:],label=legend_list[j],color=cl_list[j],marker=m_list[j])

    ax[i,0].plot(corr_ks[1:],autocorr_len['GOCI_AOD'][j,1:],
                 color=cl_list[j],ls=':',lw=1)
    ax[i,0].plot(corr_ks[1:],autocorr_len['MERRA_AOD'][j,1:],
                color=cl_list[j],ls='--',lw=1)
    ax[i,0].plot(corr_ks[1:],autocorr_len['situ_ext'][j,1:],
                 color=cl_list[j],ls='-.',lw=1)
for j in [0,5,6,7]:    
    ax[i,1].plot(corr_ks[1:],autocorr_len['aod0500'][j,1:],label=legend_list[j],color=cl_list[j],marker=m_list[j])

    ax[i,1].plot(corr_ks[1:],autocorr_len['GOCI_AOD'][j,1:],
                 color=cl_list[j],ls=':',lw=1)
    ax[i,1].plot(corr_ks[1:],autocorr_len['MERRA_AOD'][j,1:],
                 color=cl_list[j],ls='--',lw=1)
    ax[i,1].plot(corr_ks[1:],autocorr_len['situ_ext'][j,1:],
                 color=cl_list[j],ls='-.',lw=1)
#ax[i,0].set_ylim(0,1)
#ax[i,1].set_ylim(0,1)
ax[i,0].set_xscale('log')
ax[i,1].set_xscale('log')
ax[i,0].set_yscale('log')
ax[i,1].set_yscale('log')
ax[i,0].grid()
ax[i,1].grid()

ax[1,0].set_visible(False)
ax[1,1].set_visible(False)


#print 'r({})'.format(k)
ax[i,0].set_ylabel('Number of points'.format(tit2[i]))
pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)

ax[i,0].set_title('Meteorology')
ax[i,1].set_title('Altitude')
ax[i,0].set_xlabel('Distance [km]')
ax[i,1].set_xlabel('Distance [km]')

ax[i,0].legend(frameon=False,bbox_to_anchor=[0.85,-0.25])
ax[i,1].legend(frameon=False,bbox_to_anchor=[0.95,-0.25])
    

plt.savefig(fp+'plot/KORUS_Autocorr_number.png',dpi=600,transparent=True)


# In[248]:


key_list2


# In[249]:


tit2


# In[250]:


key_list3 = ['aod0500','fmf']
tit3 = ['AOD$_{{500}}$','fine-mode fraction']


# In[302]:


fig, ax = plt.subplots(2,3,figsize=(12,3.5))
for i,k in enumerate(key_list3):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
       #     ax[i,0].plot(corr_ks[1:],autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
       #                 color=cl_list[j],ls='-.',lw=1)
       # if k is 'AE':
       #     ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
       #                  color=cl_list[j],ls=':',lw=1)
       #     ax[i,0].plot(corr_ks[1:],autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
       #                  color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
        #    ax[i,1].plot(corr_ks[1:],autocorr['situ_ext'][j,1:]/autocorr['situ_ext'][j,1],
        #                 color=cl_list[j],ls='-.',lw=1)
        #if k is 'AE':
        #    ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AE'][j,1:]/autocorr['GOCI_AE'][j,1],
        #                 color=cl_list[j],ls=':',lw=1)
        #    ax[i,1].plot(corr_ks[1:],autocorr['situ_ae'][j,1:]/autocorr['situ_ae'][j,1],
        #                 color=cl_list[j],ls='-.',lw=1)
        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,1].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit3[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=False,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=False,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==0:
    #    ax[i,0].plot([],[],'d--',c='pink',label='SR 2011 Local')
    #    ax[i,0].plot([],[],'>--',c='yellow',label='SR 2011 Long')
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
    #    ax[i,0].plot([],[],'-.',c='grey',lw=1,label='In situ Ext.')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.1])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[2.3,1.1])

    #if i==1:
    #    ax[i,0].plot(SR_corr_ks,SR_AE_corr_loc,'d--',c='pink',label='SR 2011 Local')
    #    ax[i,0].plot(SR_corr_ks,SR_AE_corr_long,'>--',c='yellow',label='SR 2011 Long')
    
    if i==1:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.subplots_adjust(left=None, bottom=0.15, right=None, top=None,
                wspace=None, hspace=None)
        
plt.savefig(fp+'plot/KORUS_Autocorr_rel_subset_with_SR2011_GOCI_MERRA_insitu.png',dpi=600,transparent=True)


# ### Add the standard deviation subset

# In[251]:


corr_ks = np.array(corr_ks)


# In[406]:


fig, ax = plt.subplots(2,3,figsize=(14,4.5))
for i,k in enumerate(key_list3):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],
                         yerr=autocorr_std[k][j,1:]/autocorr[k][j,1],color=cl_list[j],
                     marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)
            ax[i,0].errorbar(corr_ks[1:]*1.02,autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1], 
                             yerr=autocorr_std['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],capsize=0,ls=':',lw=1)
            ax[i,0].errorbar(corr_ks[1:]*0.98,autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1], 
                             yerr=autocorr_std['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],capsize=0,ls='--',lw=1)

        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].errorbar(corr_ks[1:]*1.02,autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                             yerr=autocorr_std['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1,capsize=0)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,1].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],
                         yerr=autocorr_d[k][j,1:]/autocorr[k][j,1],color=cl_list[j],
                     marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)
            ax[i,1].errorbar(corr_ks[1:]*1.02,autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1], 
                             yerr=autocorr_d['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],capsize=0,ls=':',lw=1)
            ax[i,1].errorbar(corr_ks[1:]*0.98,autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1], 
                             yerr=autocorr_d['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],capsize=0,ls='--',lw=1)

        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].errorbar(corr_ks[1:]*1.02,autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                             yerr=autocorr_d['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                            color=cl_list[j],ls=':',lw=1,capsize=0)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,1].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit3[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=True,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=True,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==0:
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.1])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[2.3,1.1])

    if i==1:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.subplots_adjust(left=None, bottom=0.15, right=None, top=None,
                wspace=None, hspace=None)
        
plt.savefig(fp+'plot/KORUS_Autocorr_std_rel_subset_with_SR2011_GOCI_MERRA.png',dpi=600,transparent=True)


# In[350]:


fig, ax = plt.subplots(2,3,figsize=(12,3.5))
for i,k in enumerate(key_list3):
    for j in [0,1,2,3,4]:
        ax[i,0].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],yerr=autocorr_d[k][j,1:]/autocorr[k][j,1],color=cl_list[j],
                     marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,0].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                        color=cl_list[j],ls='--',lw=1)

        if k is 'fmf':
            ax[i,0].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    for j in [0,5,6,7]:    
        ax[i,1].plot(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],label=legend_list[j],color=cl_list[j],
                     marker=m_list[j],alpha=0.6)
        ax[i,0].errorbar(corr_ks[1:],autocorr[k][j,1:]/autocorr[k][j,1],yerr=autocorr_d[k][j,1:]/autocorr[k][j,1],color=cl_list[j],
                     marker=m_list[j],alpha=0.6,capsize=0)
        if k is 'aod0500':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_AOD'][j,1:]/autocorr['GOCI_AOD'][j,1],
                         color=cl_list[j],ls=':',lw=1)
            ax[i,1].plot(corr_ks[1:],autocorr['MERRA_AOD'][j,1:]/autocorr['MERRA_AOD'][j,1],
                         color=cl_list[j],ls='--',lw=1)

        if k is 'fmf':
            ax[i,1].plot(corr_ks[1:],autocorr['GOCI_fmf'][j,1:]/autocorr['GOCI_fmf'][j,1],
                         color=cl_list[j],ls=':',lw=1)
    ax[i,0].set_ylim(0,1)
    ax[i,1].set_ylim(0,1)
    ax[i,0].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,1].set_yticks([0,0.25,0.5,0.75,1.0])
    ax[i,0].set_xscale('log')
    ax[i,1].set_xscale('log')
    ax[i,0].grid()
    ax[i,1].grid()
    
    ax[i,2].set_visible(False)
    
    #print 'r({})'.format(k)
    ax[i,0].set_ylabel('r({})'.format(tit3[i]))
    plt.setp(ax[i,0].get_xticklabels(), visible=False)
    plt.setp(ax[i,1].get_xticklabels(), visible=False)
    pu.sub_note(note[i][0],ax=ax[i,0],out=False,fontsize=12)
    pu.sub_note(note[i][1],ax=ax[i,1],out=False,fontsize=12)
    
    if i==0:
        ax[i,0].set_title('Meteorology')
        ax[i,1].set_title('Altitude')
        
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_loc,'d--',c='pink',label='SR 2011 Local')
        ax[i,0].plot(SR_corr_ks,SR_aod_corr_long,'>--',c='yellow',label='SR 2011 Long')
        
    if i==0:
        ax[i,0].plot([],[],':',c='k',lw=1,label='GOCI YAER v2')
        ax[i,0].plot([],[],'--',c='grey',lw=1,label='MERRA2 AODANA')
        ax[i,0].legend(frameon=False,bbox_to_anchor=[3.1,1.1])
        ax[i,1].legend(frameon=False,bbox_to_anchor=[2.3,1.1])

    if i==1:
        ax[i,0].set_xlabel('Distance [km]')
        ax[i,1].set_xlabel('Distance [km]')
        plt.setp(ax[i,0].get_xticklabels(), visible=True)
        plt.setp(ax[i,1].get_xticklabels(), visible=True)

plt.subplots_adjust(left=None, bottom=0.15, right=None, top=None,
                wspace=None, hspace=None)
        
#plt.savefig(fp+'plot/KORUS_Autocorr_dif_rel_subset_with_SR2011_GOCI_MERRA.png',dpi=600,transparent=True)


# ## Autocorrelation of extrinsic vs intrisic properties

# In[252]:


autocorr.keys()


# In[349]:


plt.figure(figsize=(6,6))
j = 5
jj = 2
plt.plot(corr_ks[jj:],autocorr['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],label='Point/In situ AE',color='tab:blue',ls='--')
plt.plot(corr_ks[jj:],autocorr['AE'][j,jj:]/autocorr['AE'][j,jj],label='Column/4STAR AE',color='k',ls='--')
plt.plot(corr_ks[jj:],autocorr['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],label='Column/GOCI AE',color='tab:red',ls='--')
plt.plot(corr_ks[jj:],autocorr['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],label='Column/MERRA AE',color='tab:purple',ls='--')

plt.plot(corr_ks[jj:],autocorr['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],label='Point/In situ SSA',color='tab:orange',ls=':')
plt.plot(corr_ks[jj:],autocorr['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],label='Point/In situ Ext.',color='tab:blue',ls='-')

plt.plot(corr_ks[jj:],autocorr['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],label='Column/4STAR AOD',color='k')
plt.plot(corr_ks[jj:],autocorr['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],label='Column/GOCI AOD',color='tab:red')
plt.plot(corr_ks[jj:],autocorr['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],label='Column/MERRA AOD',color='tab:purple')


ax = plt.gca()

ax.errorbar(corr_ks[jj:],autocorr['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],
                 yerr=autocorr_d['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],color='tab:orange',ls=':',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],
                 yerr=autocorr_d['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],color='tab:blue',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],
                 yerr=autocorr_d['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],color='tab:blue',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['AE'][j,jj:]/autocorr['AE'][j,jj],
                 yerr=autocorr_d['AE'][j,jj:]/autocorr['AE'][j,jj],color='k',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],
                 yerr=autocorr_d['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],color='tab:red',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],
                 yerr=autocorr_d['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],color='tab:purple',ls='--',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],
                 yerr=autocorr_d['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],color='k',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],
                 yerr=autocorr_d['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],color='tab:red',
                 alpha=0.6,capsize=0)
ax.errorbar(corr_ks[jj:],autocorr['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],
                 yerr=autocorr_d['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],color='tab:purple',
                 alpha=0.6,capsize=0)

ax.set_ylim(0,1)
ax.set_yticks([0,0.25,0.5,0.75,1.0])
ax.set_xscale('linear')
ax.grid()
plt.legend(frameon=True,ncol=2,loc=(0,1.04))
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Autocorrelation')
plt.subplots_adjust(top=0.7)
#plt.savefig(fp+'plot/KORUS_Autocorr_intrisic_vs_extrinsic.png',dpi=600,transparent=True)


# In[355]:


plt.figure(figsize=(6,6))
j = 0
jj = 1
plt.plot(corr_ks[jj:],autocorr_mean['situ_ae'][j,jj:]/autocorr_mean['situ_ae'][j,jj],label='Point/In situ AE',color='tab:blue',ls='--')
plt.plot(corr_ks[jj:],autocorr_mean['AE'][j,jj:]/autocorr_mean['AE'][j,jj],label='Column/4STAR AE',color='k',ls='--')
plt.plot(corr_ks[jj:],autocorr_mean['GOCI_AE'][j,jj:]/autocorr_mean['GOCI_AE'][j,jj],label='Column/GOCI AE',color='tab:red',ls='--')
plt.plot(corr_ks[jj:],autocorr_mean['MERRA_AE'][j,jj:]/autocorr_mean['MERRA_AE'][j,jj],label='Column/MERRA AE',color='tab:purple',ls='--')

plt.plot(corr_ks[jj:],autocorr_mean['situ_ssa'][j,jj:]/autocorr_mean['situ_ssa'][j,jj],label='Point/In situ SSA',color='tab:orange',ls=':')
plt.plot(corr_ks[jj:],autocorr_mean['situ_ext'][j,jj:]/autocorr_mean['situ_ext'][j,jj],label='Point/In situ Ext.',color='tab:blue',ls='-')

plt.plot(corr_ks[jj:],autocorr_mean['aod0500'][j,jj:]/autocorr_mean['aod0500'][j,jj],label='Column/4STAR AOD',color='k')
plt.plot(corr_ks[jj:],autocorr_mean['GOCI_AOD'][j,jj:]/autocorr_mean['GOCI_AOD'][j,jj],label='Column/GOCI AOD',color='tab:red')
plt.plot(corr_ks[jj:],autocorr_mean['MERRA_AOD'][j,jj:]/autocorr_mean['MERRA_AOD'][j,jj],label='Column/MERRA AOD',color='tab:purple')


ax = plt.gca()
if False:
    ax.errorbar(corr_ks[jj:],autocorr['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],
                     yerr=autocorr_d['situ_ssa'][j,jj:]/autocorr['situ_ssa'][j,jj],color='tab:orange',ls=':',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],
                     yerr=autocorr_d['situ_ae'][j,jj:]/autocorr['situ_ae'][j,jj],color='tab:blue',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],
                     yerr=autocorr_d['situ_ext'][j,jj:]/autocorr['situ_ext'][j,jj],color='tab:blue',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['AE'][j,jj:]/autocorr['AE'][j,jj],
                     yerr=autocorr_d['AE'][j,jj:]/autocorr['AE'][j,jj],color='k',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],
                     yerr=autocorr_d['GOCI_AE'][j,jj:]/autocorr['GOCI_AE'][j,jj],color='tab:red',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],
                     yerr=autocorr_d['MERRA_AE'][j,jj:]/autocorr['MERRA_AE'][j,jj],color='tab:purple',ls='--',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],
                     yerr=autocorr_d['aod0500'][j,jj:]/autocorr['aod0500'][j,jj],color='k',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],
                     yerr=autocorr_d['GOCI_AOD'][j,jj:]/autocorr['GOCI_AOD'][j,jj],color='tab:red',
                     alpha=0.6,capsize=0)
    ax.errorbar(corr_ks[jj:],autocorr['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],
                     yerr=autocorr_d['MERRA_AOD'][j,jj:]/autocorr['MERRA_AOD'][j,jj],color='tab:purple',
                     alpha=0.6,capsize=0)

ax.set_ylim(0,1)
ax.set_yticks([0,0.25,0.5,0.75,1.0])
ax.set_xscale('linear')
ax.grid()
plt.legend(frameon=True,ncol=2,loc=(0,1.04))
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Autocorrelation')
plt.subplots_adjust(top=0.7)
#plt.savefig(fp+'plot/KORUS_Autocorr_intrisic_vs_extrinsic.png',dpi=600,transparent=True)


# ## Plot the comparison of GOCI to 4STAR AOD

# In[258]:


plt.figure()
plt.hist([ar['GPS_Alt'],ar['GPS_Alt'][fl]],label=['All data','QA data'],bins=50)
plt.legend()
plt.xlabel('Altitude [m]')
plt.ylabel('number of samples')
plt.title('KORUS-AQ 4STAR sampling')


# In[259]:


plt.figure()
fla = ar['fl_QA'] & (ar['GPS_Alt']<500.0)
flan = ar['fl_QA'] & (ar['GPS_Alt']<500.0) & np.isfinite(ar['AOD0501']) & np.isfinite(goci2ar['aod'])
r = np.corrcoef(ar['AOD0501'][flan],goci2ar['aod'][flan])[0,1]**2.0
plt.plot(ar['AOD0501'][fla],goci2ar['aod'][fla],'.',label='R$^2$ = {:1.3f}'.format(r))
plt.xlim(0,1.5)
plt.ylim(0,1.5)
plt.xlabel('4STAR AOD$_{{500}}$')
plt.ylabel('GOCI AOD')
plt.plot([0,1.5],[0,1.5],'--k',label='1:1')
pu.plot_lin(ar['AOD0501'][fla],goci2ar['aod'][fla],x_err=ar['UNCAOD0501'][fla],labels=True,shaded_ci=True,ci=95)

plt.legend()
plt.title('All 4STAR samples below 0.5km with nearby GOCI')


# In[175]:


len(np.unique(goci2ar['lon'][flan]))


# In[260]:


nsub = len(np.unique(goci2ar['lon'][flan]))
aod_star = np.zeros((nsub))
aod_goci = np.zeros((nsub))
ae_star = np.zeros((nsub))
ae_goci = np.zeros((nsub))
fmf_star = np.zeros((nsub))
fmf_goci = np.zeros((nsub))
aod_star_std = np.zeros((nsub))
ae_star_std = np.zeros((nsub))
fmf_star_std = np.zeros((nsub))
for j,la in enumerate(np.unique(goci2ar['lat'][flan])):
    ipixel = np.where(goci2ar['lat'][flan]==la)[0]
    aod_goci[j] = np.mean(goci2ar['aod'][flan][ipixel])
    ae_goci[j] = np.mean(goci2ar['AE'][flan][ipixel])
    fmf_goci[j] = np.mean(goci2ar['fmf'][flan][ipixel])
    aod_star[j] = np.mean(ar['AOD0501'][flan][ipixel])
    ae_star[j] = np.mean(angs[flan][ipixel])
    fmf_star[j] = np.mean(fmf['eta'][flan][ipixel])
    aod_star_std[j] = np.std(ar['AOD0501'][flan][ipixel])
    ae_star_std[j] = np.std(angs[flan][ipixel])
    fmf_star_std[j] = np.std(fmf['eta'][flan][ipixel])
    


# In[261]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
rbin = np.corrcoef(aod_star,aod_goci)[0,1]**2.0
ax[0].plot(aod_star,aod_goci,'.',label='R$^2$ = {:1.3f}'.format(rbin))
ax[0].plot([0,1.5],[0,1.5],'--k',label='1:1')
pu.plot_lin(aod_star,aod_goci,labels=True,shaded_ci=True,ci=95,ax=ax[0])
ax[0].legend()
ax[0].set_xlim(0,1.5)
ax[0].set_ylim(0,1.5)
ax[0].set_xlabel('4STAR AOD averaged within GOCI pixel')
ax[0].set_ylabel('GOCI AOD')
ax[0].set_title('KORUS-AQ Average AOD below 0.5 km')

flae = np.isfinite(ae_star) & np.isfinite(ae_goci)
rbinae = np.corrcoef(ae_star[flae],ae_goci[flae])[0,1]**2.0
ax[1].plot(ae_star,ae_goci,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
ax[1].plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(ae_star,ae_goci,labels=True,shaded_ci=True,ci=95,ax=ax[1])
ax[1].legend()
ax[1].set_ylim(0.2,1.8)
ax[1].set_xlim(0.2,1.8)
ax[1].set_xlabel('4STAR AE averaged within GOCI pixel')
ax[1].set_ylabel('GOCI AE')
ax[1].set_title('KORUS-AQ Average AE below 0.5 km')

plt.savefig(fp+'plot/KORUS_GOCI_vs_4STAR_AOD_AE.png',dpi=600,transparent=True)


# In[262]:


plt.figure()
flae = np.isfinite(ae_star) & np.isfinite(ae_goci)
rbinae = np.corrcoef(ae_star[flae],ae_goci[flae])[0,1]**2.0
plt.plot(ae_star,ae_goci,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
plt.plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(ae_star,ae_goci,labels=True,shaded_ci=True,ci=95)
plt.legend()
plt.ylim(0.2,1.8)
plt.xlim(0.2,1.8)
plt.xlabel('4STAR AE averaged within GOCI pixel')
plt.ylabel('GOCI AE')
plt.title('Average AE below 0.5 km')


# In[288]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
ax[0].hist(aod_star-aod_goci,bins=50,normed=True)
ax[0].axvline(0,ls='-',color='k',alpha=1,lw=0.5)
ax[0].axvline(np.nanmean(aod_star-aod_goci),ls='-',color='darkblue',alpha=0.6,lw=3,
            label='mean={:2.2f}'.format(np.nanmean(aod_star-aod_goci)))
ax[0].axvline(np.nanmedian(aod_star-aod_goci),ls='--',color='darkblue',alpha=0.6,lw=3,
            label='median={:2.2f}'.format(np.nanmedian(aod_star-aod_goci)))
ax[0].legend()
flae = np.isfinite(ae_star) & np.isfinite(ae_goci)
ax[1].hist(ae_star[flae]-ae_goci[flae],bins=50,normed=True)
ax[1].axvline(0,ls='-',color='k',alpha=1,lw=0.5)
ax[1].axvline(np.nanmean(ae_star[flae]-ae_goci[flae]),ls='-',color='darkblue',alpha=0.6,lw=3,
            label='mean={:2.2f}'.format(np.nanmean(ae_star[flae]-ae_goci[flae])))
ax[1].axvline(np.nanmedian(ae_star[flae]-ae_goci[flae]),ls='--',color='darkblue',alpha=0.6,lw=3,
            label='median={:2.2f}'.format(np.nanmedian(ae_star[flae]-ae_goci[flae])))
ax[1].legend()
ax[0].set_xlabel('AOD difference (4STAR-GOCI)')
ax[0].set_ylabel('Normalized counts')
ax[1].set_xlabel('AE difference (4STAR-GOCI)')
plt.savefig(fp+'plot/KORUS_hist_diff_GOCI_vs_4STAR_AOD_AE.png',dpi=600,transparent=True)


# ## Compare MERRA2 AOD to 4STAR

# In[264]:


plt.figure()
fla = ar['fl_QA'] & (ar['GPS_Alt']<500.0)
flan = ar['fl_QA'] & (ar['GPS_Alt']<500.0) & np.isfinite(ar['AOD0501']) & np.isfinite(merra2ar['aod'])
r = np.corrcoef(ar['AOD0501'][flan],merra2ar['aod'][flan])[0,1]**2.0
plt.plot(ar['AOD0501'][fla],merra2ar['aod'][fla],'.',label='R$^2$ = {:1.3f}'.format(r))
plt.xlim(0,1.5)
plt.ylim(0,1.5)
plt.xlabel('4STAR AOD$_{{500}}$')
plt.ylabel('MERRA2 AOD')
plt.plot([0,1.5],[0,1.5],'--k',label='1:1')
pu.plot_lin(ar['AOD0501'][fla],merra2ar['aod'][fla],x_err=ar['UNCAOD0501'][fla],labels=True,shaded_ci=True,ci=95)

plt.legend()
plt.title('All 4STAR samples below 0.5km with nearby MERRA2 AOD')


# In[273]:


merra2ar.keys()


# In[274]:


nsub = len(np.unique(merra2ar['ind'][flan]))
aod_star_m = np.zeros((nsub))
aod_merra = np.zeros((nsub))
aod_star_m_std = np.zeros((nsub))
ae_star_m = np.zeros((nsub))
ae_merra = np.zeros((nsub))
ae_star_m_std = np.zeros((nsub))
for j,la in enumerate(np.unique(merra2ar['ind'][flan])):
    ipixel = np.where(merra2ar['ind'][flan]==la)[0]
    aod_merra[j] = np.mean(merra2ar['aod'][flan][ipixel])
    aod_star_m[j] = np.mean(ar['AOD0501'][flan][ipixel])
    aod_star_m_std[j] = np.std(ar['AOD0501'][flan][ipixel])
    ae_merra[j] = np.mean(merra2ar['ae'][flan][ipixel])
    ae_star_m[j] = np.mean(angs[flan][ipixel])
    ae_star_m_std[j] = np.std(angs[flan][ipixel])


# In[267]:


plt.figure(figsize=(7,6))
flae = np.isfinite(aod_star_m) & np.isfinite(aod_merra)
rbinae = np.corrcoef(aod_star_m[flae],aod_merra[flae])[0,1]**2.0
plt.plot(aod_star_m,aod_merra,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
plt.errorbar(aod_star_m,aod_merra,xerr=aod_star_m_std,color='b',marker='.',ls='None',elinewidth=0.4,label='std dev')
plt.plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(aod_star_m,aod_merra,labels=True,shaded_ci=True,ci=95)
plt.legend()
plt.ylim(0.0,1.5)
plt.xlim(0.0,1.5)
plt.xlabel('4STAR AOD averaged within MERRA pixel')
plt.ylabel('MERRA2 AOD')
plt.title('Average AOD below 0.5 km during KORUS-AQ (May-June 2016)')
plt.savefig(fp+'plot/KORUS_MERRA2_vs_4STAR_AOD.png',dpi=600,transparent=True)


# In[277]:


plt.figure(figsize=(7,6))
flae = np.isfinite(ae_star_m) & np.isfinite(ae_merra)
rbinae = np.corrcoef(ae_star_m[flae],ae_merra[flae])[0,1]**2.0
plt.plot(ae_star_m,ae_merra,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
plt.errorbar(ae_star_m,ae_merra,xerr=ae_star_m_std,color='b',marker='.',ls='None',elinewidth=0.4,label='std dev')
plt.plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(ae_star_m,ae_merra,labels=True,shaded_ci=True,ci=95)
plt.legend()
plt.ylim(0.2,2.0)
plt.xlim(0.2,2.0)
plt.xlabel('4STAR AE averaged within MERRA pixel')
plt.ylabel('MERRA2 AE')
plt.title('Average AE below 0.5 km during KORUS-AQ (May-June 2016)')
plt.savefig(fp+'plot/KORUS_MERRA2_vs_4STAR_Ae.png',dpi=600,transparent=True)


# In[287]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
flae = np.isfinite(aod_star_m) & np.isfinite(aod_merra)
rbinae = np.corrcoef(aod_star_m[flae],aod_merra[flae])[0,1]**2.0
ax[0].plot(aod_star_m,aod_merra,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
ax[0].errorbar(aod_star_m,aod_merra,xerr=aod_star_m_std,color='b',marker='.',ls='None',elinewidth=0.4,label='std dev')
ax[0].plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(aod_star_m,aod_merra,labels=True,shaded_ci=True,ci=95,ax=ax[0])
ax[0].legend()
ax[0].set_ylim(0.0,1.5)
ax[0].set_xlim(0.0,1.5)
ax[0].set_xlabel('4STAR AOD averaged within MERRA pixel')
ax[0].set_ylabel('MERRA2 AOD')
fig.suptitle('Average aerosol below 0.5 km during KORUS-AQ (May-June 2016)')

flae = np.isfinite(ae_star_m) & np.isfinite(ae_merra)
rbinae = np.corrcoef(ae_star_m[flae],ae_merra[flae])[0,1]**2.0
ax[1].plot(ae_star_m,ae_merra,'.',label='R$^2$ = {:1.3f}'.format(rbinae))
ax[1].errorbar(ae_star_m,ae_merra,xerr=ae_star_m_std,color='b',marker='.',ls='None',elinewidth=0.4,label='std dev')
ax[1].plot([0,2.0],[0,2.0],'--k',label='1:1')
pu.plot_lin(ae_star_m,ae_merra,labels=True,shaded_ci=True,ci=95,ax=ax[1])
ax[1].legend()
ax[1].set_ylim(0.2,2.0)
ax[1].set_xlim(0.2,2.0)
ax[1].set_xlabel('4STAR AE averaged within MERRA pixel')
ax[1].set_ylabel('MERRA2 AE')

plt.savefig(fp+'plot/KORUS_MERRA2_vs_4STAR_AOD_AE_splitplot.png',dpi=600,transparent=True)


# In[281]:


plt.figure()
flae = np.isfinite(aod_star_m) & np.isfinite(aod_merra)
plt.hist(aod_star_m[flae]-aod_merra[flae],bins=50,normed=True)
plt.axvline(0,ls='-',color='k',alpha=1,lw=0.5)
plt.axvline(np.nanmean(aod_star_m[flae]-aod_merra[flae]),ls='-',color='darkblue',alpha=0.6,
            label='mean={:2.2f}'.format(np.nanmean(aod_star_m[flae]-aod_merra[flae])),lw=3)
plt.axvline(np.nanmedian(aod_star_m[flae]-aod_merra[flae]),ls='--',color='darkblue',alpha=0.6,
            label='median={:2.2f}'.format(np.nanmedian(aod_star_m[flae]-aod_merra[flae])),lw=3)
plt.legend()
plt.xlabel('AOD difference (4STAR-MERRA2)')
plt.ylabel('Normalized counts')
plt.savefig(fp+'plot/KORUS_hist_diff_MERRA2_vs_4STAR_AOD.png',dpi=600,transparent=True)


# In[282]:


plt.figure()
flae = np.isfinite(ae_star_m) & np.isfinite(ae_merra)
plt.hist(ae_star_m[flae]-ae_merra[flae],bins=50,normed=True)
plt.axvline(0,ls='-',color='k',alpha=1,lw=0.5)
plt.axvline(np.nanmean(ae_star_m[flae]-ae_merra[flae]),ls='-',color='darkblue',alpha=0.6,lw=3,
            label='mean={:2.2f}'.format(np.nanmean(ae_star_m[flae]-ae_merra[flae])))
plt.axvline(np.nanmedian(ae_star_m[flae]-ae_merra[flae]),ls='--',color='darkblue',alpha=0.6,lw=3,
            label='median={:2.2f}'.format(np.nanmedian(ae_star_m[flae]-ae_merra[flae])))
plt.legend()
plt.xlabel('AE difference (4STAR-MERRA2)')
plt.ylabel('Normalized counts')
plt.savefig(fp+'plot/KORUS_hist_diff_MERRA2_vs_4STAR_AE.png',dpi=600,transparent=True)


# In[289]:


fig,ax = plt.subplots(1,2,figsize=(10,5))
flae = np.isfinite(aod_star_m) & np.isfinite(aod_merra)
ax[0].hist(aod_star_m[flae]-aod_merra[flae],bins=50,normed=True)
ax[0].axvline(0,ls='-',color='k',alpha=1,lw=0.5)
ax[0].axvline(np.nanmean(aod_star_m[flae]-aod_merra[flae]),ls='-',color='darkblue',alpha=0.6,
            label='mean={:2.2f}'.format(np.nanmean(aod_star_m[flae]-aod_merra[flae])),lw=3)
ax[0].axvline(np.nanmedian(aod_star_m[flae]-aod_merra[flae]),ls='--',color='darkblue',alpha=0.6,
            label='median={:2.2f}'.format(np.nanmedian(aod_star_m[flae]-aod_merra[flae])),lw=3)
ax[0].legend()
ax[0].set_xlabel('AOD difference (4STAR-MERRA2)')
ax[0].set_ylabel('Normalized counts')

flae = np.isfinite(ae_star_m) & np.isfinite(ae_merra)
ax[1].hist(ae_star_m[flae]-ae_merra[flae],bins=50,normed=True)
ax[1].axvline(0,ls='-',color='k',alpha=1,lw=0.5)
ax[1].axvline(np.nanmean(ae_star_m[flae]-ae_merra[flae]),ls='-',color='darkblue',alpha=0.6,lw=3,
            label='mean={:2.2f}'.format(np.nanmean(ae_star_m[flae]-ae_merra[flae])))
ax[1].axvline(np.nanmedian(ae_star_m[flae]-ae_merra[flae]),ls='--',color='darkblue',alpha=0.6,lw=3,
            label='median={:2.2f}'.format(np.nanmedian(ae_star_m[flae]-ae_merra[flae])))
ax[1].legend()
ax[1].set_xlabel('AE difference (4STAR-MERRA2)')
ax[1].set_ylabel('Normalized counts')
plt.savefig(fp+'plot/KORUS_hist_diff_MERRA2_vs_4STAR_AOD_AE_2plt.png',dpi=600,transparent=True)


# ## Now get the angstrom exponent and plot it vertically

# In[212]:


nwl,nm


# In[213]:


aodrr = np.array([ar[n] for n in nwl])


# In[214]:


aodrr.shape


# In[215]:


angs = su.calc_angs(ar['Start_UTC'],np.array(nm[1:11]),aodrr[1:11,:])


# In[282]:


def make_bined_alt(x,alt,days,fl,n=70,rg=None):
    'Function to create binned data for a set range, usually for altitude'
    binned_ang,binned_alt,binned_num,binned_ndays = [],[],[],[]
    if rg:
        dz = (rg[1]-rg[0])/n
    else:
        dz = np.nanmax(alt[fl])/n
        rg = [0.0,np.nanmax(alt[fl])]
    print np.nanmax(alt[fl]),dz
    for i in xrange(n):
        flaa = (alt[fl]>=(i*dz)+rg[0]) & (alt[fl]<((i+1.0)*dz)+rg[0])
        binned_ang.append(x[fl][flaa])
        binned_alt.append(np.mean([(i*dz)+rg[0],((i+1.0)*dz)+rg[0]]))
        binned_num.append(len(x[fl][flaa]))
        binned_ndays.append(len(np.unique(days[fl][flaa])))
    return binned_ang,binned_alt,binned_num,binned_ndays


# ### Plotting of the angstrom vertical dependence over Seoul

# In[217]:


ar['fl_QA_angs'] = ar['fl'] & (ar['AOD0501']>0.05) 


# In[218]:


ar['fl_QA_angs_seoul'] = ar['fl'] & (ar['AOD0501']>0.05) & (ar['Latitude']<37.75) &                        (ar['Latitude']>36.9) & (ar['Longitude']<127.30) & (ar['Longitude']>126.60)


# In[219]:


any(ar['fl_QA_angs_seoul'])


# In[220]:


bang,balt,bnum,bndays = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs'],n=90)


# In[221]:


bangs,balts,bnums,bndayss = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_seoul'],n=90)


# In[222]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(bang,positions=np.array(balt)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom from fit between 452 nm and 865 nm')
plt.ylabel('Altitude [m]')
gr = plt.cm.RdPu
bl = plt.cm.Blues
pu.set_box_whisker_color(gr,bp,bndays)
    
bpc =plt.boxplot(bangs,positions=np.array(balts)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(bl,bpc,bndayss)
bpc['boxes'][0].set_color('grey')

ax = plt.gca()
plt.title('KORUS-AQ Angstrom Exponent')
plt.ylim(0,8000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000,8000])
plt.xlim(-0.1,2.3)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][0],bpc['whiskers'][0]],
           ['All data','Near Seoul','Mean','Median','25% - 75%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bndays)
cbaxesgr = plt.gcf().add_axes([0.83, 0.35, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.85, 0.35, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,3,6,9,12,15])
cbb.set_ticks([0,3,6,9,12,15]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_Angstrom_fit_vertical.png',
            transparent=True,dpi=500)


# ### Plot vertical angstrom exponent for different Met. regimes

# In[226]:


ar['doys']


# In[228]:


ar['fl_QA_angs'] = ar['fl'] & (ar['AOD0501']>0.05) 
ar['fl_QA_angs_met1'] = ar['fl'] & (ar['AOD0501']>0.05) & (ar['doys']> t1[0]) & (ar['doys']< t1[1])
ar['fl_QA_angs_met2'] = ar['fl'] & (ar['AOD0501']>0.05)  & (ar['doys']> t2[0]) & (ar['doys']< t2[1])
ar['fl_QA_angs_met3'] = ar['fl'] & (ar['AOD0501']>0.05)  & (ar['doys']> t3[0]) & (ar['doys']< t3[1])
ar['fl_QA_angs_met4'] = ar['fl'] & (ar['AOD0501']>0.05)  & (ar['doys']> t4[0]) & (ar['doys']< t4[1])

#ar['fl_QA_angs_seoul'] = ar['fl'] & (ar['AOD0501']>0.05) & (ar['Latitude']<37.75) &\
#                        (ar['Latitude']>36.9) & (ar['Longitude']<127.30) & (ar['Longitude']>126.60)


# In[287]:


bang,balt,bnum,bndays = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs'],n=30,rg=[0.0,8000.0])
bangm1,baltm1,bnumm1,bndaysm1 = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_met1'],n=30,rg=[0.0,8000.0])
bangm2,baltm2,bnumm2,bndaysm2 = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_met2'],n=30,rg=[0.0,8000.0])
bangm3,baltm3,bnumm3,bndaysm3 = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_met3'],n=30,rg=[0.0,8000.0])
bangm4,baltm4,bnumm4,bndaysm4 = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_met4'],n=30,rg=[0.0,8000.0])


# In[291]:


plt.figure(figsize=(5.5,7.5))
bp =plt.boxplot(bang,positions=np.array(balt)-25.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('Angstrom from fit between 452 nm and 865 nm')
plt.ylabel('Altitude [m]')

rd = plt.cm.RdPu
bl = plt.cm.Blues
og = plt.cm.YlOrBr
gr = plt.cm.Greens
k = plt.cm.Greys

pu.set_box_whisker_color(k,bp,bndays,mean_color='grey',median_color='grey')
    
bp1 = plt.boxplot(bangm1,positions=np.array(baltm1)+00.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(rd,bp1,bndaysm1,mean_color='tab:red',median_color='r')
bp1['boxes'][0].set_color('grey')
bp2 = plt.boxplot(bangm2,positions=np.array(baltm2)+25.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(bl,bp2,bndaysm2,mean_color='tab:blue',median_color='b')
bp3 = plt.boxplot(bangm3,positions=np.array(baltm3)+50.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(og,bp3,bndaysm3,mean_color='tab:orange',median_color='y')
bp4 = plt.boxplot(bangm4,positions=np.array(baltm4)+75.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(gr,bp4,bndaysm4,mean_color='tab:green',median_color='g')

ax = plt.gca()
plt.title('KORUS-AQ Angstrom Exponent')
plt.ylim(0,8000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000,8000])
plt.xlim(-0.1,2.5)
plt.grid()
plt.legend([bp['boxes'][5],bp1['boxes'][15],bp2['boxes'][16],bp3['boxes'][7],bp4['boxes'][-4],
            bp['means'][0],bp['medians'][0],bp['boxes'][-2],bp['whiskers'][0]],
           ['All data','Dynamic','Stagnation','Extreme\npollution','Blocking','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bndays)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bndays)
scalarmapk = plt.cm.ScalarMappable(cmap=k)
scalarmapk.set_array(bndays)
scalarmaprd = plt.cm.ScalarMappable(cmap=rd)
scalarmaprd.set_array(bndays)
scalarmapog = plt.cm.ScalarMappable(cmap=og)
scalarmapog.set_array(bndays)
cbaxesgr = plt.gcf().add_axes([0.83, 0.15, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.87, 0.15, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbaxesk = plt.gcf().add_axes([0.91, 0.15, 0.015, 0.3])
cbk = plt.colorbar(scalarmapk,cax=cbaxesk)
cbaxesog = plt.gcf().add_axes([0.85, 0.15, 0.015, 0.3])
cbo = plt.colorbar(scalarmapog,cax=cbaxesog)
cbaxesrd = plt.gcf().add_axes([0.89, 0.15, 0.015, 0.3])
cbr = plt.colorbar(scalarmaprd,cax=cbaxesrd)
cbg.set_ticks([0,3,6,9,12,15,18])
cbb.set_ticks([0,3,6,9,12,15,18]),cbb.set_ticklabels(['','','','','',''])
cbk.set_ticks([0,3,6,9,12,15,18]),cbk.set_ticklabels(['','','','','',''])
cbo.set_ticks([0,3,6,9,12,15,18]),cbo.set_ticklabels(['','','','','',''])
cbr.set_ticks([0,3,6,9,12,15,18]),cbr.set_ticklabels(['','','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_Angstrom_fit_vertical_met.png',
            transparent=True,dpi=500)


# ## Analyse the Fine mode fraction

# In[50]:


ar['fl_QA_low'] = ar['fl_QA'] & (ar['GPS_Alt']<500.0)
ar['fl_QA_mid'] = ar['fl_QA'] & (ar['GPS_Alt']>2000.0) & (ar['GPS_Alt']<5000.0) 


# In[51]:


ar['fl_QA_fmf'] = ar['fl_QA'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc']))


# In[52]:


bfaod,baltf,bnumf,bndaysf = make_bined_alt(fmf['tauf'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)
bcaod,baltc,bnumc,bndaysc = make_bined_alt(fmf['tauc'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)
beta,balte,bnume,bndayse = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)


# In[53]:


blat,baltl,bnuml,bndaysl = make_bined_alt(ar['Latitude'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)
blon,baltlo,bnumlo,bndayslo = make_bined_alt(ar['Longitude'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)


# In[54]:


blats = [np.nanmedian(ll) for ll in blat]
blons = [np.nanmedian(ll) for ll in blon]


# In[90]:


blons


# ### Plot the fine mode fraction distribution

# In[87]:


plt.figure()
plt.plot(aodrr[2,:],label='measurement AOD 500 nm')
plt.plot(fmf['tau'],'.k',label='fit AOD 500 nm')
plt.plot(fmf['tauf'], 'ob',label='Fine mode fraction AOD')
plt.plot(fmf['tauc'],'sr',label='Coarse mode fraction AOD')
plt.legend()
plt.ylim(0,1.5)


# In[96]:


plt.figure()
plt.hist(fmf['tauc'][ar['fl_QA']]+fmf['tauf'][ar['fl_QA']],range=[0,1.5],bins=50,label='total')
plt.hist(fmf['tauc'][ar['fl_QA']],range=[0,1.5],bins=50,label='Coarse mode')
plt.legend(frameon=False)


# In[106]:


any(ar['fl_QA_mid'])


# ### Plot the histogram distribution of the fine mode fraction

# In[61]:


plt.figure(figsize=(6,7))
ax1 = plt.subplot(3,1,1)
plt.hist([fmf['tauc'][ar['fl_QA']],fmf['tauf'][ar['fl_QA']]],color=['r','b'],histtype='bar',
            bins=50,range=[0.0,1.5],label=['Coarse','Fine'],edgecolor='None',alpha=0.75,normed=False,stacked=True)
plt.legend(frameon=True,loc=1)
plt.title('KORUS-AQ AOD fine/coarse mode - all data')
#plt.xlabel('AOD 500 nm')
plt.ylabel('Counts')
plt.yscale('log'),plt.xscale('log')
plt.ylim(5,250000),plt.xlim(0.01,1.5)
plt.xticks([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.5])
ax1.set_xticklabels([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,'',1.0,'',1.5])

ax2 = plt.subplot(3,1,2,sharex=ax1)
plt.hist([fmf['tauc'][ar['fl_QA_mid']],fmf['tauf'][ar['fl_QA_mid']]],color=['r','b'],histtype='bar',
            bins=50,range=[0.0,1.5],label=['Coarse','Fine'],edgecolor='None',alpha=0.75,normed=False,stacked=True)
#plt.legend(frameon=False)
plt.title('Between 2 and 5 km')
plt.ylabel('Counts')
plt.yscale('log'),plt.xscale('log')
plt.ylim(5,250000),plt.xlim(0.01,1.5)
plt.xticks([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.5])
ax2.set_xticklabels([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,'',1.0,'',1.5])

ax3 = plt.subplot(3,1,3,sharex=ax2)
plt.hist([fmf['tauc'][ar['fl_QA_low']],fmf['tauf'][ar['fl_QA_low']]],color=['r','b'],histtype='bar',
            bins=50,range=[0.0,1.5],label=['Coarse','Fine'],edgecolor='None',alpha=0.75,normed=False,stacked=True)
#plt.legend(frameon=False)
plt.title('Below 0.5 km')
#plt.xlabel('AOD 500 nm')
plt.ylabel('Counts')
plt.yscale('log'),plt.xscale('log')
plt.ylim(5,250000),plt.xlim(0.01,1.5)
plt.xlabel('AOD 500 nm')
plt.xticks([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.5])
ax3.set_xticklabels([0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,'',1.0,'',1.5])

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_fine_mode_hist.png',
            transparent=True,dpi=500)


# ### Plot the vertical dependence of the fine mode fraction

# In[110]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(bfaod,positions=np.array(baltf)-5.0,vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('AOD 500 nm')
plt.ylabel('Altitude [m]')
bl = plt.cm.YlOrRd
gr = plt.cm.Blues
pu.set_box_whisker_color(gr,bp,bndaysf)
    
bpc =plt.boxplot(bcaod,positions=np.array(baltc)+10.0,vert=False,
                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(bl,bpc,bndaysc)
bpc['boxes'][-1].set_color('grey')

ax = plt.gca()
plt.title('KORUS-AQ Fine/Coarse mode AOD')
plt.ylim(0,8000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000,8000])
plt.xlim(0.0,0.65)
plt.grid()
plt.legend([bp['boxes'][5],bpc['boxes'][18],bpc['means'][0],bpc['medians'][0],bpc['boxes'][-1],bpc['whiskers'][0]],
           ['Fine mode','Coarse mode','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(bndaysf)
scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
scalarmapbl.set_array(bndaysc)
cbaxesgr = plt.gcf().add_axes([0.83, 0.35, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
cbaxesbl = plt.gcf().add_axes([0.85, 0.35, 0.015, 0.3])
cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
cbg.set_ticks([0,6,12,16,18,20])
cbb.set_ticks([0,6,12,16,18,20]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left'),cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-6.0,0.5,'Days sampled',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_fine_mode_AOD_vertical.png',
            transparent=True,dpi=500)


# In[56]:


blats[0]=36.2


# In[57]:


bndm = np.nanmax(blats)*1.0
bndm


# In[100]:


cl = gr
for j,q in enumerate(blats):
    print j, q, cl(blats[j]*1.0/bndm)


# In[111]:


plt.figure(figsize=(4,6))
bp =plt.boxplot(beta,positions=np.array(balte),vert=False,
                showfliers=False,widths=90,showmeans=True,patch_artist=True)
plt.xlabel('fine mode fraction')
plt.ylabel('Altitude [m]')
bl = plt.cm.YlOrRd
gr = plt.cm.Blues
pu.set_box_whisker_color(gr,bp,blats,color_not_start_at_zero=True)
    
#bpc =plt.boxplot(bcaod,positions=np.array(baltc)+10.0,vert=False,
#                 showfliers=False,widths=90,showmeans=True,patch_artist=True)
#pu.set_box_whisker_color(bl,bpc,bndaysc)
bp['boxes'][-1].set_color('grey')

ax = plt.gca()
plt.title('KORUS-AQ fine mode fraction')
plt.ylim(0,8000)
plt.yticks([0,1000,2000,3000,4000,5000,6000,7000,8000])
ax.set_yticklabels([0,1000,2000,3000,4000,5000,6000,7000,8000])
plt.xlim(0.0,1.0)
plt.grid()
plt.legend([bp['boxes'][5],bp['means'][5],bp['medians'][5],bp['boxes'][-1],bp['whiskers'][5]],
           ['AOD fmf','Mean','Median','25\% - 75\%','min-max'],
           frameon=False,loc=1,numpoints=1)

scalarmapgr = plt.cm.ScalarMappable(cmap=gr)
scalarmapgr.set_array(blats)
#scalarmapbl = plt.cm.ScalarMappable(cmap=bl)
#scalarmapbl.set_array(bndays)
cbaxesgr = plt.gcf().add_axes([0.88, 0.35, 0.015, 0.3])
cbg = plt.colorbar(scalarmapgr,cax=cbaxesgr)
#cbaxesbl = plt.gcf().add_axes([0.85, 0.35, 0.015, 0.3])
#cbb = plt.colorbar(scalarmapbl,cax=cbaxesbl)
#cbg.set_ticks([0,6,12,15,18])
#cbb.set_ticks([0,6,12,15,18]),cbb.set_ticklabels(['','','','',''])
cbaxesgr.yaxis.set_ticks_position('left')#,cbaxesbl.yaxis.set_ticks_position('left')
cbaxesgr.text(-9.0,0.5,'Median Latitude',rotation=90,verticalalignment='center')

plt.tight_layout()

plt.savefig(fp+'plot/KORUS_4STAR_fine_mode_AOD_vertical.png',
            transparent=True,dpi=500)


# ### Split fmf for met

# In[293]:


ar['fl_QA_fmf'] = ar['fl_QA'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc']))
ar['fl_QA_fmf_met1'] = ar['fl'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc'])) & (ar['doys']> t1[0]) & (ar['doys']< t1[1])
ar['fl_QA_fmf_met2'] = ar['fl'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc'])) & (ar['doys']> t2[0]) & (ar['doys']< t2[1])
ar['fl_QA_fmf_met3'] = ar['fl'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc'])) & (ar['doys']> t3[0]) & (ar['doys']< t3[1])
ar['fl_QA_fmf_met4'] = ar['fl'] & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc'])) & (ar['doys']> t4[0]) & (ar['doys']< t4[1])


# In[ ]:


beta,balte,bnume,bndayse = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf'],n=90)
beta1,balte1,bnume1,bndayse1 = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf_met1'],n=90)
beta2,balte2,bnume2,bndayse2 = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf_met2'],n=90)
beta3,balte3,bnume3,bndayse3 = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf_met3'],n=90)
beta4,balte4,bnume4,bndayse4 = make_bined_alt(fmf['eta'],ar['GPS_Alt'],ar['days'],ar['fl_QA_fmf_met4'],n=90)


# ## Calculate the autocorrelation of the fine and coarse mode AOD

# In[55]:


fvals = {'utc':ar['Start_UTC'][fl],'alt':ar['GPS_Alt'][fl],'lat':ar['Latitude'][fl],'lon':ar['Longitude'][fl],
        'aod0500':ar['AOD0501'][fl],'aod1040':ar['AOD1040'][fl],'aodf':fmf['tauf'][fl],'aodc':fmf['tauc'][fl],'eta':fmf['eta'][fl]}


# In[56]:


dfvals = get_segments(f_level,fvals,nsep=100)


# In[57]:


ddfv = get_distances(dfvals)


# Now the segments are identified and the cumulative distances are quantified, we must interpolate over the segments, to remove any missing data.

# In[58]:


def interp_dist_fmf(d,dist=0.12):
    'function to insterpolate the AOD from the dict to an even grid spacing accroding to distance (default 0.12 km)'
    d['cdist_n'],d['aod_nf'],d['aod_nc'],d['eta_n'] = [],[],[],[]
    for i,cd in enumerate(d['cumdist']):        
        d['cdist_n'].append(np.arange(cd.min(),cd.max(),dist))
        if np.sum(np.isfinite(d['aodf'][i]))/float(len(d['aodf'][i])) < 0.75: # check if at least 75% of the segment is valid
            af = np.array(np.nan)
            ac = np.array(np.nan)
            et = np.array(np.nan)
        else:
            try:
                fcdf = interpolate.interp1d(cd,d['aodf'][i])
                af = fcdf(d['cdist_n'][i])
                fcdc = interpolate.interp1d(cd,d['aodc'][i])
                ac = fcdc(d['cdist_n'][i])
                fcde = interpolate.interp1d(cd,d['eta'][i])
                et = fcde(d['cdist_n'][i])
            except TypeError:
                af = np.array(np.nan)
                ac = np.array(np.nan)
                et = np.array(np.nan)
        d['aod_nf'].append(af)
        d['aod_nc'].append(ac)
        d['eta_n'].append(et)


# In[59]:


interp_dist_fmf(dfvals)


# In[79]:


dfvals['autocor_f'],dfvals['autocor_c'],dfvals['autocor_e'] = [] ,[],[]
for i,a in enumerate(dfvals['aod_nf']):
    #auf,auc,eut = [],[],[]
    try:
        auf = autocorr5(a)
        auc = autocorr5(dfvals['aod_nc'][i])
        eut = autocorr5(dfvals['eta_n'][i])
    except:
        auf = np.array([np.nan])
        auc = np.array([np.nan])
        eut = np.array([np.nan])
    dfvals['autocor_f'].append(auf[:])
    dfvals['autocor_c'].append(auc[:])
    dfvals['autocor_e'].append(eut[:])


# In[61]:


len(dfvals['aod_nf'])


# In[63]:


dfvals['aod_nf'][i]


# In[64]:


dfvals['cdist_n'][i]


# In[80]:


dfvals['autocor_f'][i]


# In[111]:


mc = np.max([len(m) for m in dfvals['autocor_c']])
imc = np.argmax([len(m) for m in dfvals['autocor_c']])


# In[113]:


cdist = dfvals['cdist_n'][imc]


# In[123]:


autocor_c = np.zeros((len(dfvals['autocor_c']),mc))+np.nan
autocor_f = np.zeros((len(dfvals['autocor_f']),mc))+np.nan
autocor_e = np.zeros((len(dfvals['autocor_e']),mc))+np.nan


# In[124]:


for i,c in enumerate(dfvals['autocor_c']): autocor_c[i,:len(c)]=c
for i,c in enumerate(dfvals['autocor_f']): autocor_f[i,:len(c)]=c
for i,c in enumerate(dfvals['autocor_e']): autocor_e[i,:len(c)]=c


# ### Plot out the autocorrelation 

# In[83]:


plt.figure()
for i,j in enumerate(dfvals['cdist_n']):
    try:
        if len(j)<1: continue
        plt.plot(j,dfvals['autocor_f'][i],'.')
    except:
        continue
plt.ylabel('Correlation Coefficient for fine mode')
plt.xlabel('Distance [km]')
plt.xscale('log')
plt.xlim(0.1,500)


# In[84]:


plt.figure()
for i,j in enumerate(dfvals['cdist_n']):
    try:
        plt.plot(j,dfvals['autocor_c'][i],'.')
    except:
        pass
plt.ylabel('Correlation Coefficient for coarse mode')
plt.xlabel('Distance [km]')
plt.xscale('log')


# In[110]:


dfvals['cdist_n'][1:3]


# In[133]:


autocor_c_ma = np.ma.masked_array(autocor_c,mask=np.isnan(autocor_c))


# In[157]:


def make_binned(x,alt,fl,bn,flb):
    'Function to create binned data for a set range, usually for altitude'
    import numpy as np
    binned_ang,binned_alt,binned_num = [],[],[]
    for i,b in enumerate(bn[:-1]):
        flaa = (alt[flb]>=b) & (alt[flb]<bn[i+1])
        binned_ang.append(x[:,flb][flaa])
        binned_alt.append(np.mean([b,bn[i+1]]))
        binned_num.append(len(x[fl][:,flaa]))
    return binned_ang,binned_alt,binned_num,binned_ndays


# In[154]:


bnc = np.logspace(0.1,3.0)


# In[155]:


bnc


# In[163]:


auc = make_binned(autocor_c,cdist,np.isfinite(autocor_c),bnc)


# In[161]:


b = bnc[0]
i = 0
fl = np.isfinite(autocor_c)


# In[162]:


flaa = (cdist[fl]>=b) & (cdist[fl]<bnc[i+1])


# In[156]:


autocor_c.shape


# In[135]:


plt.figure()
bp = plt.boxplot(autocor_c_ma[:,0:20],positions=cdist[0:20],vert=True,
            showfliers=False,widths=90,showmeans=True,patch_artist=True) 


# In[118]:


autocor_c.shape


# In[134]:


autocor_c_ma[:,0:20]


# In[130]:


np.isfinite(autocor_c)


# In[136]:


autocor_c_ma[:,0]


# ## Map out the level legs and numbers

# In[299]:



#set up a easy plotting function
def make_map(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=128,lat_0=36.0,
            llcrnrlon=123.0, llcrnrlat=32.0,
            urcrnrlon=132.0, urcrnrlat=39,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(123,133,11),labels=[0,0,0,1])
    m.drawparallels(np.linspace(31,39,17),labels=[1,0,0,0])
    return m


# In[303]:


fig,ax = plt.subplots(1,1)
m = make_map(ax)
m.plot(ar['Longitude'],ar['Latitude'],'.',markersize=0.2,color='tab:blue',latlon=True,label='All data')
m.plot(ar['Longitude'][fl][f_level],ar['Latitude'][fl][f_level],'.',markersize=0.5,color='tab:red',latlon=True,label='level legs')
#m.plot(s['Lon'][it],s['Lat'][it],'r+',latlon=True)
plt.legend(markerscale=20)
plt.savefig(fp+'plot/KORUS_map_QA.png',dpi=600,transparent=True)


# ## Plot the AOD histogram by met

# In[294]:


ar['fl_QA_met1'] = ar['fl'] & (ar['doys']> t1[0]) & (ar['doys']< t1[1])
ar['fl_QA_met2'] = ar['fl'] & (ar['doys']> t2[0]) & (ar['doys']< t2[1])
ar['fl_QA_met3'] = ar['fl'] & (ar['doys']> t3[0]) & (ar['doys']< t3[1])
ar['fl_QA_met4'] = ar['fl'] & (ar['doys']> t4[0]) & (ar['doys']< t4[1])


# In[337]:


fig = plt.figure()
n=plt.hist([ar['AOD0501'][ar['fl']],
            ar['AOD0501'][ar['fl_QA_met1']],
          ar['AOD0501'][ar['fl_QA_met2']],
          ar['AOD0501'][ar['fl_QA_met3']],ar['AOD0501'][ar['fl_QA_met4']]],
           bins=20,range=(0,1.4),normed=True,edgecolor='None',alpha=0.7,
         label=['All data','Dynamic','Stagnation','Extreme\npollution','Blocking'])
#y = [(nn+n[1][j+1])/2.0 for j,nn in enumerate(n[1][:-1])]
#for i,p in enumerate(n[-1]):
#    plt.plot(y,n[0][i],'-',color=p[0].get_facecolor(),lw=3)
plt.legend(frameon=False)
plt.grid()
plt.xlim(0,0.5)
plt.xlabel('AOD @ 501 nm')
plt.ylabel('Normalized counts')
#plt.title('Histogram of 4STAR AOD from KORUS-AQ subsetted by altitude')

left, bottom, width, height = [0.63, 0.3, 0.35, 0.2]
ax2 = fig.add_axes([left, bottom, width, height])
n = ax2.hist([ar['AOD0501'][ar['fl']],
            ar['AOD0501'][ar['fl_QA_met1']],
          ar['AOD0501'][ar['fl_QA_met2']],
          ar['AOD0501'][ar['fl_QA_met3']],ar['AOD0501'][ar['fl_QA_met4']]],
             bins=20,range=(0,1.4),normed=True,edgecolor='None',alpha=0.7)
ax2.set_xlim(0.4,1.4)
ax2.set_ylim(0,0.6)
ax2.grid()
ax2.set_xlabel('AOD @ 501 nm')

plt.savefig(fp+'plot/AOD_hist_met_KORUS.png',dpi=600,transparent=True)


# ## Plot the AOD spectra by met

# In[290]:


aod_names = sorted([a for a in ar.keys() if ('AOD' in a) and not ('UNC' in a)])


# In[291]:


aod_names


# In[292]:


wvl = np.array([380,452,501,520,532,550,606,620,675,781,865,1020,1040,1064,1236,1559,1627])
wvl_bins = np.append(wvl[0]-10,wvl+10)


# In[295]:


ars = []
wvs = []
fls = {'fl':[],'fl_QA_met1':[],'fl_QA_met2':[],'fl_QA_met3':[],'fl_QA_met4':[]}
for i,a in enumerate(aod_names):
    ars.append(ar[a])
    wvs.append(ar[a]*0.0+wvl[i])
    for ff in fls.keys():
        fls[ff].append(ar[ff])
ars = np.array(ars)
wvs = np.array(wvs)
for ff in fls.keys():
    fls[ff] = np.array(fls[ff])
    fls[ff] = fls[ff].reshape(fls[ff].size)
arsn = ars.reshape(ars.size)
wvsn = wvs.reshape(wvs.size)


# In[296]:


plt.figure()
plt.plot(wvsn[fls['fl']],arsn[fls['fl']],'.',alpha=0)
plt.yscale('log')
plt.xscale('log')
plt.xticks([350,400,500,600,700,800,900,1000,1200,1400,1700])
plt.xlim(350,1700)
plt.ylim([0.01,2.0])

pu.make_boxplot(arsn[fls['fl']],wvsn[fls['fl']],wvl_bins,wvl,y=1,alpha=0.5,label='All points',
                fliers_off=True,color='k')
pu.make_boxplot(arsn[fls['fl_QA_met1']],wvsn[fls['fl_QA_met1']],wvl_bins,wvl+2,y=1,alpha=0.5,label='Dynamic',
                fliers_off=True,color='tab:orange')
pu.make_boxplot(arsn[fls['fl_QA_met2']],wvsn[fls['fl_QA_met2']],wvl_bins,wvl+4,y=1,alpha=0.5,label='Stagnation',
                fliers_off=True,color='tab:green')
pu.make_boxplot(arsn[fls['fl_QA_met3']],wvsn[fls['fl_QA_met3']],wvl_bins,wvl+6,y=1,alpha=0.5,label='Extreme\npollution',
                fliers_off=True,color='tab:red')
pu.make_boxplot(arsn[fls['fl_QA_met4']],wvsn[fls['fl_QA_met4']],wvl_bins,wvl+8,y=1,alpha=0.5,label='Blocking',
                fliers_off=True,color='tab:purple')

plt.legend(frameon=False)

plt.xlabel('Wavelength [nm]')
plt.ylabel('AOD')
plt.grid()
plt.title('Average AOD spectra')

plt.savefig(fp+'plot/KORUS_AOD_wvl_loglog_bymet.png',dpi=600,transparent=True)


# In[ ]:




