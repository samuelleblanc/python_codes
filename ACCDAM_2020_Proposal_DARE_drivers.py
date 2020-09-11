#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     ACCDAM_2020_Proposal_DARE_drivers
#     
# Purpose:  
# 
#     Make some figures for the DARE ACCDAM proposal, based on SEAC4RS, ORACLES, KORUS-AQ, and ARCTAS 
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2020-09-03
#     Modified: 

# # Prepare python environment

# In[76]:


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
import scipy.io as sio


# In[3]:


import map_utils as mu
import write_utils as wu
from scipy import interpolate
import math


# In[4]:


from linfit import linfit
import Sun_utils as su


# In[77]:


get_ipython().magic(u'matplotlib notebook')


# In[74]:


fpk = getpath('KORUS')
fpo = getpath('ORACLES')
fps = getpath('SEAC4RS')
fpa = getpath('ARCTAS',path='/data/sam/arctas',make_path=True)


# # Load the files

# ## Load the SEAC4RS files

# In[12]:


fls = os.listdir(fps+'merge/')
fls.sort()


# In[13]:


fls


# In[16]:


se,sed = [],[]
for f in fls:
    s,d = lu.load_ict(fps+'merge/'+f,return_header=True)
    se.append(s)
    sed.append(d)


# In[17]:


len(se)


# In[22]:


se[0]['JDAY']


# In[67]:


s = {'AOD380nm_4STAR':[],'AOD452nm_4STAR':[],'AOD501nm_4STAR':[],'AOD865nm_4STAR':[],'LATITUDE':[],'LONGITUDE':[],'UTC':[],'GPS_ALT':[],
     'JDAY':[],'DN500nm_SSFR':[],'UP500nm_SSFR':[],'UP350to2150nm_SSFR':[],'UP350to2150nm_SSFR':[],
     'DNSolarIrrad_BBR':[],'UPSolarIrrad_BBR':[],'Flag_4STAR':[],'ABS660nmTD_PAS_NOAAAeroAbs':[]}
for e in se:
    for k in s.keys():
        s[k] = np.append(s[k],e[k])
        


# In[68]:


fl = (s['Flag_4STAR']==0) & (s['GPS_ALT']<1.0)


# In[69]:


len(s['AOD501nm_4STAR'])


# ### Calculate the Angstrom Exponent

# In[70]:


# calculate the angstrom exponent
s['AE'] = su.calc_angs(s['UTC'],[380.0,452.0,501.0,865.0],
np.array([s['AOD380nm_4STAR'],s['AOD452nm_4STAR'],s['AOD501nm_4STAR'],s['AOD865nm_4STAR']]))


# In[50]:


np.nanmin(s['AE'][fl]),np.nanmax(s['AE'][fl])


# ### Plot out the SEAC4RS ratios

# In[66]:


plt.figure()
plt.scatter(s['UP500nm_SSFR'][fl]/s['DN500nm_SSFR'][fl],s['AOD501nm_4STAR'][fl],c=s['AE'][fl],marker='.',vmin=0,vmax=2.0)
plt.xlim([0,1])
plt.yscale('log')
plt.ylim(0.01,2.0)
plt.ylabel('AOD at 500 nm')
plt.xlabel('Flight-level Albedo at 500 nm')
plt.colorbar(label='Angstrom Exponent')
plt.title('SEAC4RS')

plt.savefig(fps+'plots/SEAC4RS_AOD_vs_Alebdo_AE.png',transparent=True,dpi=600)


# In[60]:


s['GPS_ALT'].min()


# In[72]:


plt.figure()
plt.scatter(s['UP500nm_SSFR'][fl]/s['DN500nm_SSFR'][fl],s['AOD501nm_4STAR'][fl],c=s['ABS660nmTD_PAS_NOAAAeroAbs'][fl],marker='.',vmin=0,vmax=5.0)
plt.xlim([0,1])
plt.yscale('log')
plt.ylim(0.01,2.0)
plt.ylabel('AOD at 500 nm')
plt.xlabel('Flight-level Albedo at 500 nm')
plt.colorbar(label='Aerosol Absorption [1/Mm]')
plt.title('SEAC4RS')

plt.savefig(fps+'plots/SEAC4RS_AOD_vs_Alebdo_abs.png',transparent=True,dpi=600)


# ## Load the ARCTAS retrievals

# In[79]:


a = sio.idl.readsav(fpa+'nasa/20080709/seven/20080709_forcing_abs.out')


# In[80]:


a.keys()


# In[82]:


a['wvl_arr']


# In[83]:


a['asy'].shape


# In[84]:


a['forcing'].shape


# In[86]:


for k in a.keys():
    print k, a[k].shape


# In[85]:


a['eff'].shape


# In[97]:


a['forcing'][1,:,:]


# In[110]:


forc = np.zeros((165,13))+np.nan
for q in xrange(165):
    if any(np.isfinite(a['forcing'][1,q,:])):
        forc[q,:] = Sp.smooth(a['forcing'][1,q,:],1,nan=True,old=True)


# In[112]:


dare = np.trapz(forc,a['wvl_arr'])


# ### Plot the DARE

# In[120]:


plt.figure()
plt.scatter(a['aot'][:,3],dare,c=a['asy'][0:165,3],vmin=0.5,vmax=0.8)
plt.xlabel('AOD at 500nm')
plt.ylabel('DARE [W/m$^2$]')
plt.title('ARCTAS Biomass burning plume')
plt.colorbar(label='Asymmetry Parameter')
plt.savefig(fpa+'plots/ARCTAS_DARE_vs_AOD_asy.png',transparent=True,dpi=600)


# ## Load ARCTAS results and compare to RH

# In[121]:


rh,rhdict = lu.load_ict(fpa+'nasa/pds_p3b_20080709_r2.ict',return_header=True)


# In[122]:


rh['C_RelHumidityWater']


# In[146]:


rh_utc = rh['Time']/3600.0


# In[151]:


il10 = np.argmin(abs(rh_utc-19.613))
il11 = np.argmin(abs(rh_utc-19.686))

il20 = np.argmin(abs(rh_utc-19.245))
il21 = np.argmin(abs(rh_utc-19.343))


# In[152]:


rh['FMS_LAT'][il10],rh['FMS_LAT'][il11],rh['FMS_LAT'][il20],rh['FMS_LAT'][il21]


# In[159]:


fl1 = Sp.find_closest(rh['FMS_LAT'][il10:il11],rtma[:,0])
fl2 = Sp.find_closest(rh['FMS_LAT'][il20:il21],rtma[:,0])


# In[160]:


rh_rt = rh['C_RelHumidityWater'][il10:il11][fl1]
rh_rt2 = rh['C_RelHumidityWater'][il20:il21][fl2]


# In[163]:


rh_rt


# In[136]:


rtma = np.genfromtxt(fpa+'nasa/20080709/seven/rtm_hsrl_20080709_wvl0499.txt',skip_header=1)


# In[140]:


#lat     lon   ssa   asy   asy2    albedo     correction    tau modification    flux divergence   model down      model up       tau    Ft_up   Ft_dn    Fb_up   Fb_dn


# In[137]:


rtma


# In[138]:


rtma.shape


# In[178]:


rtma[:,11]


# In[177]:


plt.figure()
plt.plot(rtma[:,0],rtma[:,2],'v',label='SSA')
plt.plot(rtma[:,0],rtma[:,3],'+',label='ASY')
plt.plot(rtma[:,0],rtma[:,5],'x',label='Albedo')

plt.plot(rtma[:,0],rh_rt/100.0,'.',label='RH/100 [$\%$]')
plt.plot(rtma[:,0],rtma[:,11]*rtma[:,7],'o',label='AOD')


plt.legend(frameon=False)
plt.xlabel('Latitude [$^{{\\circ}}$]')
#plt.xlim(60.9267,61.06)


# ### Load from digitized plots

# In[179]:


eFR_sur = np.array([[60.89086, -19.70582],
[60.89557, -21.68807],
[60.89807, -23.73104],
[60.90017, -22.86450],
[60.90222, -25.25048],
[60.90583, -25.61747],
[60.90887, -24.82233],
[60.91159, -27.02213],
[60.91358, -29.43621],
[60.91358, -31.71722],
[60.91571, -24.87841],
[60.91780, -24.18429],
[60.92007, -28.04244],
[60.91954, -31.13430],
[60.91968, -32.41279],
[60.92217, -25.16222],
[60.92522, -25.66982],
[60.92841, -23.66392],
[60.93045, -25.78450],
[60.93381, -25.43481],
[60.93714, -26.75110],
[60.94066, -26.31281],
[60.94351, -23.40437],
[60.94933, -23.68676],
[60.95266, -23.76326],
[60.95618, -23.80096],
[60.95810, -24.80771],
[60.95986, -27.24801],
[60.96318, -27.15017],
[60.96651, -26.50992],
[60.96964, -26.43704],
[60.97316, -26.06240],
[60.97649, -26.41010],
[60.97920, -27.23452],
[60.98092, -24.66534],
[60.98408, -26.76816],
[60.98547, -28.29201],
[60.98757, -24.98299],
[60.99006, -22.79133],
[61.00336, -23.16706],
[61.00471, -25.94510],
[61.00803, -25.81128],
[61.01144, -26.10915],
[61.01445, -27.57565],
[61.01744, -26.38632],
[61.01896, -29.72664],
[61.02276, -30.49346],
[61.02581, -30.40538],
[61.02926, -28.71374],
[61.03226, -26.85016],
[61.03578, -27.30296],
[61.03911, -26.83706],
[61.04243, -27.70779],
[61.04548, -27.21485],
[61.04881, -30.02081],
[61.05189, -25.65805],
[61.05047, -28.53838],
[61.05463, -27.08428]])


# In[184]:


eFR_toa = np.array([[60.89086, -11.73246],
[60.89507, -12.66105],
[60.89724, -10.08398],
[60.90056, -11.82644],
[60.90402, -11.34209],
[60.90693, -11.96015],
[60.91164, -11.31949],
[60.91469, -12.56224],
[60.91940, -13.78845],
[60.92267, -15.32614],
[60.92384, -16.95301],
[60.92605, -13.22525],
[60.92966, -15.41041],
[60.93160, -13.75439],
[60.93603, -13.42376],
[60.93880, -12.15409],
[60.94268, -13.69298],
[60.94573, -13.84705],
[60.94905, -12.43872],
[60.95252, -12.26625],
[60.95653, -13.71470],
[60.95986, -12.15753],
[60.96318, -11.51728],
[60.96623, -12.63993],
[60.96928, -12.34845],
[60.97427, -12.75377],
[60.97870, -13.68229],
[60.98175, -13.91385],
[60.98535, -14.41644],
[60.98873, -13.95826],
[60.99062, -13.67875],
[61.00475, -15.47612],
[61.00807, -14.51624],
[61.01157, -15.95529],
[61.01417, -14.34008],
[61.01744, -14.43791],
[61.02077, -12.04259],
[61.02387, -12.08816],
[61.02697, -13.99340],
[61.03024, -13.47327],
[61.03274, -10.95422],
[61.03673, -10.36414],
[61.03994, -10.17722],
[61.04327, -11.40827],
[61.04798, -12.82875],
[61.05102, -13.35088],
[61.05352, -9.88262],
[61.05490, -14.93820]])


# In[181]:


eFR_sur.shape


# In[189]:


aod = np.array([[60.89008, 0.40310],
[60.89073, 0.43132],
[60.89260, 0.45282],
[60.89624, 0.45517],
[60.89955, 0.46196],
[60.90128, 0.46810],
[60.90302, 0.46236],
[60.90614, 0.46643],
[60.90884, 0.49038],
[60.91165, 0.51336],
[60.91405, 0.47566],
[60.91773, 0.47631],
[60.91949, 0.49515],
[60.92100, 0.46528],
[60.92117, 0.43938],
[60.92313, 0.41292],
[60.92439, 0.46055],
[60.92397, 0.43737],
[60.92565, 0.49246],
[60.92813, 0.50921],
[60.93013, 0.48339],
[60.93265, 0.51631],
[60.93601, 0.53933],
[60.93938, 0.54168],
[60.94294, 0.52519],
[60.94470, 0.50993],
[60.94878, 0.55701],
[60.95226, 0.56184],
[60.95478, 0.57813],
[60.95646, 0.59963],
[60.95862, 0.65132],
[60.95702, 0.62281],
[60.96166, 0.65578],
[60.96514, 0.65506],
[60.96863, 0.65722],
[60.97187, 0.64398],
[60.97523, 0.65238],
[60.97775, 0.64532],
[60.97921, 0.62019],
[60.97876, 0.58915],
[60.98271, 0.62771],
[60.98580, 0.62050],
[60.98643, 0.64398],
[60.98839, 0.59929],
[60.98951, 0.56805],
[61.00324, 0.55377],
[61.00499, 0.53261],
[61.00792, 0.55960],
[61.00828, 0.54319],
[61.01136, 0.57000],
[61.01336, 0.60399],
[61.01528, 0.58189],
[61.01818, 0.60618],
[61.01948, 0.64095],
[61.02014, 0.66346],
[61.02242, 0.68391],
[61.02529, 0.65434],
[61.02807, 0.66296],
[61.03013, 0.64666],
[61.03041, 0.62684],
[61.03153, 0.68530],
[61.03366, 0.70203],
[61.03629, 0.73401],
[61.03965, 0.73032],
[61.04301, 0.72830],
[61.04609, 0.70566],
[61.04889, 0.68278],
[61.05001, 0.65305],
[61.04983, 0.62483],
[61.05226, 0.74889],
[61.05151, 0.72091],
[61.05113, 0.68127],
[61.05282, 0.77122],
[61.05338, 0.79415]])


# In[190]:


asy = np.array([[60.89036, 0.72145],
[60.89176, 0.74478],
[60.89512, 0.70978],
[60.89736, 0.73322],
[60.90055, 0.70606],
[60.90334, 0.70021],
[60.90636, 0.68252],
[60.90969, 0.68844],
[60.91277, 0.67508],
[60.91407, 0.65501],
[60.91697, 0.67687],
[60.92089, 0.65774],
[60.92061, 0.62371],
[60.92355, 0.61197],
[60.92565, 0.65736],
[60.92845, 0.66342],
[60.92985, 0.63188],
[60.93450, 0.64575],
[60.93713, 0.64758],
[60.94106, 0.64639],
[60.94386, 0.66720],
[60.94862, 0.65988],
[60.95198, 0.66752],
[60.95506, 0.65080],
[60.95822, 0.66017],
[60.96178, 0.65879],
[60.96514, 0.65711],
[60.96851, 0.64786],
[60.97187, 0.65227],
[60.97523, 0.64786],
[60.97871, 0.62701],
[60.97859, 0.60160],
[60.98195, 0.64134],
[60.98527, 0.61836],
[60.98503, 0.60034],
[60.98867, 0.64113],
[60.99063, 0.65585],
[61.00352, 0.61866],
[61.00604, 0.60728],
[61.00940, 0.60328],
[61.01283, 0.60949],
[61.01632, 0.60809],
[61.01948, 0.60665],
[61.02284, 0.62347],
[61.02525, 0.60135],
[61.02845, 0.61085],
[61.03083, 0.62620],
[61.03349, 0.66131],
[61.03685, 0.65879],
[61.04021, 0.66110],
[61.04357, 0.64723],
[61.04637, 0.63030],
[61.04945, 0.61043],
[61.05247, 0.66152],
[61.05450, 0.67098]])


# In[191]:


ssa = np.array([[60.89120, 0.94978],
[60.89540, 0.94045],
[60.89848, 0.92035],
[60.90156, 0.92203],
[60.90464, 0.91067],
[60.90800, 0.92224],
[60.91109, 0.90866],
[60.91361, 0.88954],
[60.91635, 0.93162],
[60.91799, 0.93086],
[60.92005, 0.90992],
[60.92313, 0.93990],
[60.92621, 0.92935],
[60.92929, 0.94537],
[60.93265, 0.93065],
[60.93601, 0.92140],
[60.93938, 0.91530],
[60.94218, 0.93086],
[60.94386, 0.94978],
[60.94890, 0.93086],
[60.95226, 0.93128],
[60.95562, 0.94200],
[60.95898, 0.91698],
[60.96234, 0.90100],
[60.96570, 0.91698],
[60.96907, 0.91467],
[60.97243, 0.92014],
[60.97579, 0.92056],
[60.97915, 0.92497],
[60.98251, 0.93128],
[60.98587, 0.92119],
[60.98895, 0.94272],
[61.00492, 0.93927],
[61.00828, 0.93654],
[61.01164, 0.93675],
[61.01500, 0.92771],
[61.01724, 0.93275],
[61.01892, 0.90079],
[61.02228, 0.88923],
[61.02565, 0.89575],
[61.02901, 0.91152],
[61.03181, 0.91604],
[61.04133, 0.89806],
[61.04665, 0.90134],
[61.04903, 0.89112],
[61.05141, 0.91541],
[61.05394, 0.89175]])


# In[228]:


fig,ax = plt.subplots(2,1,sharex=True,figsize=(4,3))
ax[0].plot(eFR_toa[:,0],eFR_toa[:,1],'-',color='k',label=('TOA'))
ax[0].plot(eFR_sur[:,0],eFR_sur[:,1],'--',color='grey',label=('Surface'))
ax[0].legend(frameon=True,loc=2)
ax[0].set_ylabel('RFari/AOD$_{{500}}$/I$^\downarrow$ \n[\%]')
ax[0].set_title('ARCTAS 20080709 Forest Fire Plume')
ax[1].plot(ssa[:,0],ssa[:,1],'-',marker='.',lw=3,label='SSA')
ax[1].plot(asy[:,0],asy[:,1],'-',lw=3,label='ASY')
ax[1].plot(aod[:,0],aod[:,1],'-',lw=3,label='AOD')
ax[1].legend(frameon=True,loc=2)
ax[1].set_xlabel('Latitude [$^\\circ$]')
ax[1].set_ylabel('Optical property')
plt.tight_layout(h_pad=0.02,pad=0.3)
plt.xlim(60.78,61.06)
plt.savefig(fpa+'ARCTAS_RFari_20080709_results.png',dpi=600,transparent=True)


# In[229]:


fig,ax = plt.subplots(2,1,sharex=True,figsize=(4,3))
ax[0].plot(eFR_toa[:,0],eFR_toa[:,1],'-',color='k',label=('TOA'))
ax[0].plot(eFR_sur[:,0],eFR_sur[:,1],'--',color='grey',label=('Surface'))
ax[0].legend(frameon=True,loc=2)
ax[0].set_ylabel('DARE/AOD$_{{500}}$/I$^\downarrow$ \n[\%]')
ax[0].set_title('ARCTAS 20080709 Forest Fire Plume')
ax[1].plot(ssa[:,0],ssa[:,1],'-',marker='.',lw=3,label='SSA')
ax[1].plot(asy[:,0],asy[:,1],'-',lw=3,label='ASY')
ax[1].plot(aod[:,0],aod[:,1],'-',lw=3,label='AOD')
ax[1].legend(frameon=True,loc=2)
ax[1].set_xlabel('Latitude [$^\\circ$]')
ax[1].set_ylabel('Optical property')
plt.tight_layout(h_pad=0.02,pad=0.3)
plt.xlim(60.78,61.06)
plt.savefig(fpa+'ARCTAS_DARE_20080709_results.png',dpi=600,transparent=True)


# In[ ]:




