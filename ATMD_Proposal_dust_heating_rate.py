
# coding: utf-8

# # Info
# Name:  
# 
#     ATDM Proposal dust heating rate
# 
# Purpose:  
# 
#     Make some figures for the ATMD proposal
#     Focus on the heating rate profiles obtained during ORACLES 2018 transit back, near Cabo Verde, with SSFR
#   
# Input:
# 
#     None
#   
# Output:
# 
#     Figures
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - load_utils.py : for loading OMI HDF5 files
#     - matplotlib
#     - numpy
#     - 
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - for_Sam_20181025.out
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-08-23
#     Modified: 

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


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


fp =getpath('ORACLES')


# # Load files

# ## Load the SSFR files

# In[4]:


ssfr = sio.idl.readsav(fp+'data_other_2018/SSFR/for_Sam_20181025_SSFR.out')


# In[5]:


ssfr.keys()


# Now interpolate the nadir spectra to the zenith wavelengths

# In[6]:


def interp_spline(xold, yold, ynew):
    uv = UnivariateSpline(xold,yold,ext=0,k=1)
    return uv(ynew)


# In[7]:


ssfr['nadspectra'].shape


# In[8]:


ssfr['nadlambda'][0:10]


# In[9]:


ssfr['nspectra'] = np.array([interp_spline(ssfr['nadlambda'],ssfr['nadspectra'][i,:],ssfr['zenlambda']) for i in xrange(len(ssfr['utc']))])


# In[10]:


ssfr['nspectra'].shape


# In[11]:


ssfr['f_net'] = ssfr['zenspectra'] - ssfr['nadspectra']


# In[97]:


plt.figure()
ax = plt.subplot(2,1,1)
ax.plot(ssfr['utc'],ssfr['f_net'][:,200],'.')
ax.plot(ssfr['utc'][i_ssfr],ssfr['f_net'][i_ssfr,200],'go')

ax1 = plt.subplot(2,1,2,sharex=ax)
ax1.plot(ssfr['utc'],ssfr['roll'],'.')
ax1.plot(ssfr['utc'][i_ssfr],ssfr['roll'][i_ssfr],'go')


# ## Load the 4STAR file

# In[12]:


s,sh = lu.load_ict(fp+'aod_ict_2018/4STAR-AOD_P3_20181025_R1.ict',return_header=True)


# In[13]:


sh


# In[14]:


sp = sio.loadmat(fp+'data_2018/4STAR_20181025starsun.mat')


# ## Load the merge file

# In[15]:


mrg,mrg_head = lu.load_netcdf(fp+'data_other_2018/mrg1_P3_20181025_R13.nc',everything=True)


# In[16]:


ssfr['pitch'] = interp_spline(mrg['Start_UTC']/3600.0,mrg['SysPitchAng'],ssfr['utc'])
ssfr['roll'] = interp_spline(mrg['Start_UTC']/3600.0,mrg['Roll_Angle'],ssfr['utc'])


# ## Load the HSRL

# In[17]:


hsrl,hsrl_hed = lu.load_hdf(fp+'data_other_2018/HSRL/HSRL2_P3_20181025_R1.h5',values=(('ext',33),('alt',36),('gps_alt',63),('time',66)))


# In[18]:


hsrl_hed


# In[19]:


hsrl['time']


# In[20]:


np.argmin(abs(hsrl['time']-13.535))


# In[21]:


hsrl['ext'].shape, hsrl['alt'].shape


# In[22]:


plt.figure()
plt.plot(hsrl['ext'][2421,:],hsrl['alt'].T)


# In[22]:


plt.figure()
plt.plot(hsrl['time'],hsrl['ext'][:,200])


# # Build up the plots to select the region
# Do some of the calculations to the data here

# Get the location of the dust profile

# In[14]:


plt.figure()
plt.plot(s['Start_UTC'],s['GPS_Alt'],'.')


# In[10]:


plt.figure()
plt.plot(s['Latitude'],s['GPS_Alt'],'.')


# In[23]:


pfl = [13.565, 13.95]


# In[24]:


i_s = (s['Start_UTC']>pfl[0]) & (s['Start_UTC']<pfl[1]) & (s['qual_flag']==0)


# In[25]:


i_ssfr = (ssfr['utc']>pfl[0]) & (ssfr['utc']<pfl[1]) & (abs(ssfr['roll'])<7.0) & (ssfr['f_net'][:,120]>0.9)


# In[26]:


ni_ssfr = np.where(i_ssfr)[0]
ds = np.where(np.diff(ni_ssfr,1)>3)[0]
disc_ssfr = ni_ssfr[ds]
is_ssfr = np.array([range(i-12,i) for i in disc_ssfr]).flatten()
is_ssfr


# In[104]:


plt.figure()
ax = plt.subplot(5,1,1)
ax.plot(ssfr['utc'],ssfr['zenspectra'][:,200],'.')
ax.plot(ssfr['utc'][i_ssfr],ssfr['zenspectra'][i_ssfr,200],'go')
ax.plot(ssfr['utc'][is_ssfr],ssfr['zenspectra'][is_ssfr,200],'rx')

ax0 = plt.subplot(5,1,2,sharex=ax)
ax0.plot(ssfr['utc'],ssfr['nspectra'][:,200],'.')
ax0.plot(ssfr['utc'][i_ssfr],ssfr['nspectra'][i_ssfr,200],'go')
ax0.plot(ssfr['utc'][is_ssfr],ssfr['nspectra'][is_ssfr,200],'rx')

ax1 = plt.subplot(5,1,3,sharex=ax)
ax1.plot(ssfr['utc'],ssfr['roll'],'.')
ax1.plot(ssfr['utc'][i_ssfr],ssfr['roll'][i_ssfr],'go')

ax2 = plt.subplot(5,1,4,sharex=ax)
ax2.plot(ssfr['utc'],ssfr['pitch'],'.')
ax2.plot(ssfr['utc'][i_ssfr],ssfr['pitch'][i_ssfr],'go')
ax2.set_xlim(13.5787,13.67)

ax3 = plt.subplot(5,1,5,sharex=ax)
ax3.plot(mrg['Start_UTC']/3600.0,mrg['Sun_Azimuth'],'.')
#ax3.plot(ssfr['utc'][i_ssfr],ssfr['pitch'][i_ssfr],'go')
ax3.set_xlim(13.5787,13.67)


# In[27]:


ssfr['alt'] = interp_spline(mrg['Start_UTC']/3600.0,mrg['GPS_Altitude'],ssfr['utc'])


# In[105]:


ssfr['alt'][1:100]


# In[106]:


plt.figure()
pc = plt.pcolor(ssfr['zenlambda'],ssfr['alt'][is_ssfr],ssfr['f_net'][is_ssfr,:],cmap=plt.cm.plasma)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Altitude [m]')
plt.xlim(350,2150)
plt.title('Net Irradiance from SSFR near Cabo Verde 2018-10-25')
cb = plt.colorbar()
cb.set_label('Net Irradiance [W/m$^2$ nm]')


# In[22]:


ssfr['alt'][is_ssfr]


# ## Now calculate the heating rate using the differential method

# In[28]:


from Sp_parameters import deriv,smooth


# In[29]:


ssfr['heat_r'] = np.zeros_like(ssfr['f_net'])


# In[30]:


nb = 75
f_net = np.vstack([ssfr['f_net'][is_ssfr[0:nb],:]*0.0+ssfr['f_net'][is_ssfr[0],:],ssfr['f_net'][is_ssfr,:]])
f_net.shape


# In[31]:


iss_ssfr = np.hstack([range(is_ssfr[0]-nb,is_ssfr[0]),is_ssfr])


# In[32]:


for l,w in enumerate(ssfr['zenlambda']):
    ssfr['heat_r'][iss_ssfr,l] = smooth(deriv(smooth(f_net[:,l],30,nan=False,old=True),
                                  ssfr['alt'][iss_ssfr]),9,nan=False,old=True) 


# In[112]:


plt.figure()
pc = plt.pcolor(ssfr['zenlambda'],ssfr['alt'][is_ssfr],ssfr['heat_r'][is_ssfr,:],cmap=plt.cm.plasma,vmin=-0.0005,vmax=0.002)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Altitude [m]')
plt.xlim(350,2150)
plt.ylim(200,6000)
plt.title('Heating rate profile from SSFR near Cabo Verde 2018-10-25')
cb = plt.colorbar(extend='both')
cb.set_label('Heating rate [W/m$^2$ nm m]')


# In[113]:


plt.figure()
plt.plot(ssfr['utc'][is_ssfr],ssfr['heat_r'][is_ssfr,200],'.')
plt.plot(ssfr['utc'][i_ssfr],ssfr['heat_r'][i_ssfr,200],'go')


# ## Now combine with other data

# In[114]:


plt.figure()
plt.plot(s['AOD0501'][i_s],s['GPS_Alt'][i_s],'.')
plt.plot(ssfr['f_net'][is_ssfr,120],ssfr['alt'][is_ssfr],'x')


# In[33]:


i_m = (mrg['Start_UTC']/3600.0>pfl[0]) & (mrg['Start_UTC']/3600.0<pfl[1])


# In[116]:


mrg.keys()


# In[126]:


plt.figure()
plt.plot(mrg['WINDS_WSPD_ms-1'][i_m],mrg['GPS_Altitude'][i_m],'.')
#plt.plot(mrg['WINDS_w_ms-1'][i_m],mrg['GPS_Altitude'][i_m],'gx')


# ### Make the extinction coefficient

# In[34]:


import scipy.interpolate as si


# In[35]:


alt = mrg['GPS_Altitude'][i_m]


# In[36]:


aod_fx =  si.interp1d(s['GPS_Alt'][i_s],s['AOD0501'][i_s],fill_value="extrapolate")
aod = aod_fx(mrg['GPS_Altitude'][i_m])


# In[37]:


ext0501 = smooth(deriv(smooth(aod,30,nan=False,old=True),
                                  alt),9,nan=False,old=True)*-1.0


# In[154]:


plt.figure()
plt.plot(ext0501,alt,'.')


# In[35]:


mrg['amass_aer'][i_m]


# In[38]:


m_air = 28.9647
R = 287.05
rho = mrg['Static_Pressure']*100.0/(R*(mrg['Static_Air_Temp']+273.15))


# In[37]:


rho[i_m]


# In[39]:


#rho=1.0
cp = 1004.0
inst_to_24h = 1.75
ssfr['heat_r_Kday'] = ssfr['heat_r']*inst_to_24h*3600*24/cp
ssfr['heat_r_Kday'] = ssfr['heat_r_Kday'].T/rho
ssfr['heat_r_Kday'] = ssfr['heat_r_Kday'].T


# In[290]:


plt.figure(figsize=(8,6))
pc = plt.pcolor(ssfr['zenlambda'],ssfr['alt'][is_ssfr],ssfr['heat_r_Kday'][is_ssfr,:],cmap=plt.cm.magma,vmin=0,vmax=0.3)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Altitude [m]')
plt.xlim(350,2150)
plt.ylim(500,6000)
plt.title('Heating rate profile from SSFR near Cabo Verde 2018-10-25')
cb = plt.colorbar(extend='both')
cb.set_label('Heating rate [K/day nm]')
plt.tight_layout()
plt.savefig(fp+'plot/Heating_rate_Caboverde_20181025.png',transparent=True,dpi=600)


# In[46]:


plt.figure(figsize=(8,4.5))
pc = plt.pcolor(ssfr['zenlambda'],ssfr['alt'][is_ssfr],ssfr['heat_r_Kday'][is_ssfr,:],cmap=plt.cm.magma,vmin=0,vmax=0.3)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Altitude [m]')
plt.xlim(350,2150)
plt.ylim(500,6000)
plt.title('Heating rate profile from SSFR near Cabo Verde 2018-10-25')
for rg in shade_rg:
    plt.axvspan(rg[0],rg[1],color='white',alpha=0.25,zorder=200)

cb = plt.colorbar(extend='max')
cb.set_label('Heating rate [K/day nm]')
plt.tight_layout()
plt.savefig(fp+'plot/Heating_rate_Caboverde_20181025_shaded.png',transparent=True,dpi=600)


# In[40]:


fp


# In[182]:


plt.figure()
plt.plot(ssfr['zenlambda'],ssfr['zenspectra'][is_ssfr[-1],:],'.')


# In[41]:


# for water vapor:
shade_rg = [[684.0,705.0],
            [713.0,741.0],
            [756.0,770.0],
            [785.5,805.0],
            [810.0,842.0],
            [890.0,998.0],
            [1077.0,1228.0],
            [1257.0,1277.0],
            [1293.0,1521.0],
            [1760.0,1992.0]]


# In[42]:


irg = [(ssfr['zenlambda']>rg[0]) & (ssfr['zenlambda']<rg[1]) for i,rg in enumerate(shade_rg)]


# In[43]:


len(irg)


# In[44]:


iwv = irg[0] | irg[1] | irg[2] | irg[3] | irg[4] | irg[5] | irg[6] | irg[7] | irg[8] | irg[9]


# In[52]:


iwv


# In[45]:


heating_r_wv = np.zeros_like(ssfr['utc'])
heating_r_aero = np.zeros_like(ssfr['utc'])
for i in is_ssfr:
    heating_r_wv[i] = np.trapz(ssfr['heat_r_Kday'][i,iwv],x=ssfr['zenlambda'][iwv])
    heating_r_aero[i] = np.trapz(ssfr['heat_r_Kday'][i,~iwv],x=ssfr['zenlambda'][~iwv])


# In[54]:


plt.figure()
plt.plot(heating_r_wv[is_ssfr]/1000.0,ssfr['alt'][is_ssfr],'-',label='Water vapor')
plt.plot(heating_r_aero[is_ssfr]/1000.0,ssfr['alt'][is_ssfr],'r-x',label='Aerosol')
plt.ylim(500,6000)
plt.legend(frameon=False)
plt.ylabel('Altitude [m]')
plt.xlabel('Heating Rate [K/day]')


# In[46]:


def interp_spline2(xold, yold, xnew,k=5):
    uv = UnivariateSpline(xold,yold,ext=0,k=k)
    return uv(xnew)


# In[47]:


alt_id = np.arange(0,6200)


# In[48]:


igs =  np.isfinite(heating_r_wv[is_ssfr])


# In[49]:


H_wv = interp_spline2(ssfr['alt'][is_ssfr[igs]][-1:0:-1],heating_r_wv[is_ssfr[igs]][-1:0:-1]/100.0,alt_id,k=4)
H_aero = interp_spline2(ssfr['alt'][is_ssfr[igs]][-1:0:-1],heating_r_aero[is_ssfr[igs]][-1:0:-1]/100.0,alt_id,k=4)


# In[59]:


plt.figure()
plt.plot(heating_r_wv[is_ssfr]/1000.0,ssfr['alt'][is_ssfr],'o',label='Water vapor')
plt.plot(heating_r_aero[is_ssfr]/1000.0,ssfr['alt'][is_ssfr],'rx',label='Aerosol')
plt.plot(H_wv,alt_id,'b-.')
plt.plot(H_aero,alt_id,'r-.')
plt.ylim(500,6000)
plt.xlim(-0.05,0.3)
plt.legend(frameon=False)
plt.ylabel('Altitude [m]')
plt.xlabel('Heating Rate [K/day]')
plt.grid()


# In[80]:


fig = plt.figure(figsize=(12,3))

ax1 = plt.subplot2grid((1,6),(0,0))
ax1.plot(heating_r_wv[is_ssfr]/100.0,ssfr['alt'][is_ssfr],'o',label='Water\nvapor')
ax1.plot(heating_r_aero[is_ssfr]/100.0,ssfr['alt'][is_ssfr],'rx',label='Aerosol')
ax1.plot(H_wv,alt_id,'b-.')
ax1.plot(H_aero,alt_id,'r-.')
ax1.set_ylabel('Altitude [m]')
ax1.set_ylim(500,4500)
ax1.set_xlim(0,2.5)
leg = ax1.legend(frameon=False,loc=2,numpoints=1,markerscale=0,handlelength=0.2,labelspacing=0.1)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
    line.set_linewidth(5.0)
ax1.set_xlabel('Heating Rate [K/day]')
#plt.setp(ax1.get_yticklabels(), visible=False)

ax0 = plt.subplot2grid((1,6),(0,1),sharey=ax1)
ax0.plot(s['AOD0501'][i_s],s['GPS_Alt'][i_s],'.',label='Total')
plt.plot(fmf['tauc'][fmf['iq']],s['GPS_Alt'][fmf['iq']],'r.',label='Coarse')
plt.plot(fmf['tauf'][fmf['iq']],s['GPS_Alt'][fmf['iq']],'g.',label='Fine')
ax0.set_xlabel('AOD at 501 nm')
ax0.set_xlim(0,0.6)
leg = ax0.legend(frameon=False,loc=2,numpoints=1,markerscale=0,handlelength=0.2,labelspacing=0.1)
for line,text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())
    line.set_linewidth(5.0)
plt.setp(ax0.get_yticklabels(), visible=False)

ax2 = plt.subplot2grid((1,6),(0,2),sharey=ax1)
ax2.plot(s['AOD_angstrom_470_865'][i_s],s['GPS_Alt'][i_s],'.',color='slateblue')
ax2.set_xlabel('Angstrom\nexponent (470/865)')
plt.setp(ax2.get_yticklabels(), visible=False)

ax3 = plt.subplot2grid((1,6),(0,3),sharey=ax1)
ax3.plot(ext0501*100000.0,alt,'.',color='lime',label='4STAR')
#ax3.plot(hsrl['ext'][2429,:]*100.0,hsrl['alt'].T,'.',color='m',label='HSRL2')
ax3.set_xlabel('Extinction coef.\nat 501 nm [1/Mm]')
#leg = ax3.legend(frameon=False,loc=5,numpoints=1,markerscale=0,handlelength=0.2,labelspacing=0.1)
#for line,text in zip(leg.get_lines(), leg.get_texts()):
#    text.set_color(line.get_color())
#    line.set_linewidth(5.0)
plt.setp(ax3.get_yticklabels(), visible=False)
ax3.set_xlim(0,100)



ax5 = plt.subplot2grid((1,6),(0,4),sharey=ax1)

ax5.plot(mrg['WINDS_WSPD_ms-1'][i_m][0:-1:15],mrg['GPS_Altitude'][i_m][0:-1:15],'ko',markersize=3)
ax5.barbs(mrg['WINDS_WSPD_ms-1'][i_m][0:-1:15],mrg['GPS_Altitude'][i_m][0:-1:15],
          mrg['WINDS_u_ms-1'][i_m][0:-1:15],mrg['WINDS_v_ms-1'][i_m][0:-1:15],alpha=0.6,color='grey')
ax5.set_xlabel('Wind speed [m/s]')
ax5.set_xlim(0,20)
plt.setp(ax5.get_yticklabels(), visible=False)

ax6 = plt.subplot2grid((1,6),(0,5),sharey=ax1)
ax6.plot(mrg['Static_Air_Temp'][i_m],mrg['GPS_Altitude'][i_m],'.',color='firebrick')
ax6.set_xlabel('Temperature [K]')
ax6.set_xlim(0,30)
plt.setp(ax6.get_yticklabels(), visible=False)

ax6.set_ylim(500,4500)
plt.tight_layout()

plt.savefig(fp+'plot/multipanel_heatingrate_20181025_CaboVerde.png',dpi=600,transparent=True)


# Make some potential temperature figures

# In[139]:


plt.figure()
plt.plot(ext0501*100000.0,alt,'.',color='lime',label='4STAR')
plt.plot(hsrl['ext'][2419,:]*100.0,hsrl['alt'].T,'.',color='m',label='HSRL2')


# In[51]:


pot_temp = mrg['Static_Air_Temp']*(mrg['Static_Pressure']/1013.0)**(0.286)


# In[292]:


plt.figure()
plt.plot(pot_temp[i_m],mrg['GPS_Altitude'][i_m],'.')


# In[52]:


import sun_utils as su


# In[53]:


ka = s.dtype.names


# In[54]:


kk = list(ka)


# In[55]:


kk.sort()


# In[56]:


kk


# In[57]:


nwl = kk[0:23]


# In[58]:


nm = [355.0,380.0,452.0,470.0,501.0,520.0,530.0,532.0,550.0,606.0,620.0,660.0,675.0,700.0,781.0,865.0,1020.0,1040.0,1064.0,1236.0,1250.0,1559.0,1627.0,1650.0]


# In[59]:


aodrr = np.array([s[n] for n in nwl])


# In[60]:


fmf = su.sda(aodrr[1:13,:],np.array(nm[1:13])/1000.0)


# In[61]:


fmf.keys()


# In[62]:


fmf['iq'] = (s['Start_UTC']>pfl[0]) & (s['Start_UTC']<pfl[1]) & (s['qual_flag']==0) & (np.isfinite(fmf['tauf'])) & (np.isfinite(fmf['tauc']))


# In[63]:


s['Start_UTC'].shape


# In[64]:


fmf['tauc'].shape


# In[96]:


plt.figure()
plt.plot(fmf['tau'][i_s],s['GPS_Alt'][i_s],'.')
plt.plot(fmf['tauc'][fmf['iq']],s['GPS_Alt'][fmf['iq']],'r.')
plt.plot(fmf['tauf'][fmf['iq']],s['GPS_Alt'][fmf['iq']],'g.')


# In[102]:


fp

