#!/usr/bin/env python
# coding: utf-8

# # Intro
# Name:  
# 
#     Make_atmos_profile_and_albedo
# 
# Purpose:  
# 
#     Create the input libradtran files for creating the atmospheric profiles and surface albedo for ARCSIX
#   
# Input:
# 
#     none
# 
# Output:
#    
#     input files for libradtran 2.0 (uvspec) 
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - mpl_toolkits for basemap, map plotting
#     - pdb
#     - datetime
# 
#   
# Needed Files:
# 
#   - atmospheric profile file sample
#   - surface albedo file sample
#     
# History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2026-11-13
#              

# # Prepare the python environment

# In[1]:


import numpy as np
import scipy.io as sio
import os
import Run_libradtran as RL
import matplotlib.pyplot as plt


# In[2]:


import load_utils as lu


# In[3]:


fp = '/data2/ARCSIX/data/'
fp_atm = '/home/sam/libradtran/libRadtran-2.0.1/data/atmmod/'
fp_out = '/scratch/rtm/input/ARCSIX_v1/'


# # Load the files

#  Need to focus on the data found at the cloud wall timing, around 15:17 UTC to 15:46 UTC

# In[47]:


os.listdir(fp)


# ## Load the met data

# In[33]:


met,met_header =lu.load_ict(fp+'ARCSIX-MetNav_P3B_20240607_R0.ict',return_header=True)


# In[28]:


met_header


# In[59]:


fig,ax = plt.subplots(3,1,sharex=True)
ax[0].plot(met['Time_Start']/3600.0,met['GPS_Altitude'],'.')
ax[1].plot(met['Time_Start']/3600.0,met['Static_Pressure'],'.')
ax[2].plot(met['Time_Start']/3600.0,met['Mixing_Ratio'],'.')


# ## Load the SSFR data

# In[22]:


ssfr,ssfr_head = lu.load_hdf_h5py(fp+'ARCSIX-SSFR_P3B_20240607_R0.h5',all_values=True)


# In[17]:


import importlib 
importlib.reload(lu)


# In[24]:


ssfr['time']/3600.0


# ## Load the atm uvspec data

# In[34]:


atm = np.loadtxt(fp_atm+'afglss.dat')


# In[37]:


atm.shape


# In[42]:


plt.figure()
plt.plot(atm[:,1],atm[:,0],'x-')


# The first few lines of the subarctic summer atmos profile
# > AFGL atmospheric constituent profile. subarctic summer. ( AFGL-TR-86-0110)       
# >     z(km)      p(mb)        T(K)    air(cm-3)    o3(cm-3)     o2(cm-3)    h2o(cm-3)    co2(cm-3)     no2(cm-3)
# >    120.000     0.00002     380.000 4.307612E+11 2.155500E+02 3.125475E+10 8.622001E+04 1.508850E+07 6.509610E+01

# In[64]:


def gas_law(p,t):
    # p in mb
    # T in Kelvin
    R = 8.314           # Joule/mol/Kelvin
    Na= 6.02297e23 # molecules/mol
    return p/t/R*Na/10000.0


# In[100]:


def h2o_mass_to_number_density(w_g_kg, p_mb, t_k):
    # Convert g/kg to dimensionless volume mixing ratio (mol/mol)
    M_dry = 28.97    # g/mol
    M_h2o = 18.015   # g/mol
    vmr = (w_g_kg / 1000.0) * (M_dry / M_h2o)
    
    # Multiply by total air number density (#/cm^3)
    n_air = gas_law(p_mb, t_k)
    return vmr * n_air


# ## Load the drop sonde data

# In[48]:


sonde,sonde_hdr = lu.load_ict(fp+'ARCSIX-AVAPS_G3_20240607162310_R1.ict',return_header=True)


# In[49]:


sonde_hdr


# In[55]:


fig,ax  = plt.subplots(1,3,sharey=True)
ax[0].plot(sonde['Temperature'],sonde['GPS_Altitude']/1000.0,'.-')
ax[0].set_xlabel('Temp [K]')
ax[0].set_ylabel('Altitude [Km]')
#ax[1].figure()
ax[1].plot(sonde['Pressure'],sonde['GPS_Altitude']/1000.0,'.-')
ax[1].set_xlabel('Pressure [hPa]')
#plt.figure()
ax[2].plot(sonde['MixingRatio'],sonde['GPS_Altitude']/1000.0,'.-')
ax[2].set_xlabel('H20 Mixing ratio [g/kg]')


# # Make ratios to standard atm

# In[82]:


from scipy.interpolate import CubicSpline


# In[97]:


sonde_air = gas_law(sonde['Pressure'],sonde['Temperature']+273.15)


# In[109]:


g = (np.isfinite(sonde['GPS_Altitude']/1000.0)) & (np.isfinite(sonde_air))
fx_air = CubicSpline(np.flip(sonde['GPS_Altitude'][g]/1000.0),np.flip(sonde_air[g]),extrapolate=False)
sonde_air_int = fx_air(atm[:,0])


# In[110]:


sonde_air_int


# In[111]:


sonde_air_ratio = sonde_air_int / atm[:,3]


# In[130]:


ib = ~np.isfinite(sonde_air_ratio)


# In[132]:


sonde_air_ratio[ib] = 1.0


# In[133]:


sonde_air_ratio


# In[113]:


sonde_h20 = h2o_mass_to_number_density(sonde['MixingRatio'],sonde['Pressure'],sonde['Temperature']+273.15)


# In[114]:


sonde_h20


# # plot our the vertical comparisons

# In[116]:


#z(km)      p(mb)        T(K)    air(cm-3)    o3(cm-3)     o2(cm-3)    h2o(cm-3)    co2(cm-3)     no2(cm-3)


# In[129]:


fig,ax = plt.subplots(1,4,sharey=True,figsize=(11,6))
ax[0].set_yscale('log')
ax[0].plot(atm[:,1],atm[:,0],'b.')
ax[0].plot(sonde['Pressure'],sonde['GPS_Altitude']/1000.0,'g+')
ax[0].set_xlabel('Pressure')

ax[1].plot(atm[:,2],atm[:,0],'b.')
ax[1].plot(sonde['Temperature']+273.15,sonde['GPS_Altitude']/1000.0,'g+')
ax[1].set_xlabel('Temperature')

ax[2].plot(atm[:,3],atm[:,0],'b.')
ax[2].plot(sonde_air,sonde['GPS_Altitude']/1000.0,'g+')
ax[2].set_xlabel('Air')

ax[3].plot(atm[:,6],atm[:,0],'b.')
ax[3].plot(sonde_h20,sonde['GPS_Altitude']/1000.0,'g+')
ax[3].set_xlabel('H2O')


# # Build the atm profile from the sonde

# In[ ]:


#z(km)      p(mb)        T(K)    air(cm-3)    o3(cm-3)     o2(cm-3)    h2o(cm-3)    co2(cm-3)     no2(cm-3)


# In[135]:


atm_z = atm[:,0]
atm_air = atm[:,3]*sonde_air_ratio
atm_o3 = atm[:,4]*sonde_air_ratio
atm_o2 = atm[:,5]*sonde_air_ratio
atm_co2 = atm[:,7]*sonde_air_ratio
atm_no2 = atm[:,8]*sonde_air_ratio


# In[136]:


gp = (np.isfinite(sonde['GPS_Altitude']/1000.0)) & (np.isfinite(sonde['Pressure']))
fx_air = CubicSpline(np.flip(sonde['GPS_Altitude'][gp]/1000.0),np.flip(sonde['Pressure'][gp]),extrapolate=False)
atm_p = fx_air(atm[:,0])
p_nan = ~np.isfinite(atm_p)
atm_p[p_nan] = atm[p_nan,1]


# In[146]:


gt = (np.isfinite(sonde['GPS_Altitude']/1000.0)) & (np.isfinite(sonde['Temperature']))
fx_air = CubicSpline(np.flip(sonde['GPS_Altitude'][gt]/1000.0),np.flip(sonde['Temperature'][gt]+273.15),extrapolate=False)
atm_t = fx_air(atm[:,0])
t_nan = ~np.isfinite(atm_t)
atm_t[t_nan] = atm[t_nan,2]


# In[141]:


gh = (np.isfinite(sonde['GPS_Altitude']/1000.0)) & (np.isfinite(sonde_h20))
fx_air = CubicSpline(np.flip(sonde['GPS_Altitude'][gh]/1000.0),np.flip(sonde_h20[gh]),extrapolate=False)
atm_h2o = fx_air(atm[:,0])
h_nan = ~np.isfinite(atm_h2o)
atm_h2o[h_nan] = atm[h_nan,6]


# ## Write the profile to file

# In[144]:


fp


# In[147]:


with open(fp+'atm_20240607.dat','w') as f:
    f.write('# Atmopshere profile based on afgss (subarctic summer) and drop sonde profile from G3 on 20240607\n')
    f.write('#     z(km)      p(mb)        T(K)    air(cm-3)    o3(cm-3)     o2(cm-3)    h2o(cm-3)    co2(cm-3)     no2(cm-3)\n')
    
    for z, p, t, air, o3, o2, h2o, co2, no2 in zip(
        atm_z, atm_p, atm_t, atm_air, atm_o3, atm_o2, atm_h2o, atm_co2, atm_no2
    ):
        line = (
            f"{z:11.3f}"
            f"{p:12.5f}"
            f"{t:12.3f}"
            f"{air:13.6E}"
            f"{o3:13.6E}"
            f"{o2:13.6E}"
            f"{h2o:13.6E}"
            f"{co2:13.6E}"
            f"{no2:13.6E}\n"
        )
        f.write(line)


# # Plot out for the time period in question

# In[32]:


plt.figure()
plt.plot(met['Time_Start']/3600.0,met['Static_Air_Temp'])
plt.plot(met['Time_Start']/3600.0,met['Static_Pressure'])
plt.plot(met['Time_Start']/3600.0,met['Relative_Humidity'])
plt.xlim(14.8,16)


# In[44]:


ival = np.argmin(abs(met['Time_Start']/3600.0 - 15.3))


# In[45]:


ival


# In[46]:


met['Latitude'][ival],met['Longitude'][ival],


# # Now set the albedo

# In[148]:


istart = np.argmin(abs(ssfr['time']/3600.0-15.3)) 
iend = np.argmin(abs(ssfr['time']/3600.0-15.7)) 


# In[149]:


istart,iend


# In[150]:


alb = ssfr['f_up']/ssfr['f_dn']


# In[153]:


ssfr['f_up'][istart,:]


# In[155]:


ssfr['wvl_up']


# In[157]:


from scipy.interpolate import make_interp_spline
interpolator = make_interp_spline(ssfr['wvl_up'], ssfr['f_up'], k=1, axis=1)
f_up_interp = interpolator(ssfr['wvl_dn'])  # Output shape: (19627, 380)


# In[158]:


f_up_interp = np.empty((ssfr['f_up'].shape[0], len(ssfr['wvl_dn'])))

for i in range(ssfr['f_up'].shape[0]):
    row = ssfr['f_up'][i]
    valid = np.isfinite(row)
    
    if np.sum(valid) >= 2:  # Need at least 2 valid points to interpolate
        f_up_interp[i] = np.interp(ssfr['wvl_dn'], ssfr['wvl_up'][valid], row[valid])
    else:
        f_up_interp[i] = np.nan  # Keep row as NaN if completely empty


# In[160]:


albedo = f_up_interp/ssfr['f_dn']


# In[164]:


avg_albedo = np.nanmean(albedo[istart:iend,:],axis=0)


# In[165]:


avg_albedo.shape


# In[173]:


plt.figure()
plt.plot(ssfr['wvl_dn'],albedo[istart:iend,:].T,'.-')
plt.ylim(0,1)


# In[177]:


plt.figure()
plt.plot(ssfr['wvl_dn'],ssfr['f_dn'][istart,:].T,'.-')
plt.ylim(0,0.01)


# In[185]:


ibad = ssfr['f_dn'][istart,:]<0.02


# In[186]:


avg_albedo[ibad] = np.nan


# In[187]:


avg_albedo


# ## Now fill back up with extrapolated values

# In[189]:


# now fill back up with extrapolated values

avg_albedo_good = np.interp(ssfr['wvl_dn'], ssfr['wvl_dn'][~ibad], avg_albedo[~ibad])


# In[191]:


plt.figure()

plt.plot(ssfr['wvl_dn'],avg_albedo_good,'.-')
plt.plot(ssfr['wvl_dn'],avg_albedo,'.-')
plt.ylim(0,1)


# ## Now save the albedo file

# In[194]:


wvl = np.arange(350,2201)


# In[196]:


alb = np.interp(wvl,ssfr['wvl_dn'][~ibad], avg_albedo[~ibad])


# In[199]:


with open(fp+'alb_20240607.dat','w') as f:
    f.write('# wvl[nm]    alb[unitless]\n')
    for i,w in enumerate(wvl):
        aa = alb[i]
        f.write(f'{w:10.6f}\t{aa:10.6f}\n')
    


# In[198]:


alb[i]


# In[201]:


top_time = 54899
bot_time = 55042


# In[202]:


itop = np.argmin(abs(met['Time_Start']-top_time))
ibot = np.argmin(abs(met['Time_Start']-bot_time))


# In[204]:


met['GPS_Altitude'][itop]/1000.0,met['GPS_Altitude'][ibot]/1000.0


# In[207]:


met['Solar_Zenith_Angle'][istart],met['Solar_Zenith_Angle'][iend]


# In[ ]:




