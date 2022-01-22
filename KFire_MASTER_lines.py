#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Exploration of MASTER data for evidence of potassium emission lines in active fire regions
# 
# Input:
# 
#     arguments
# 
# Output:
# 
#     Figure and save files
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - load_utils.py
#     - matplotlib
#     - numpy
#     - Sp_parameters
#     - write_utils
#     - path_utils
#     - hdf5storage
#     - scipy
# 
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-01-19
#     Modified:
# 

# # Prepare python environment

# In[3]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import os


# In[13]:


name = 'KFire'
vv = 'v1'
fp = getpath(name)


# # Load files

# In[37]:


os.listdir(fp+'MASTER/')


# ## 2020 western diversity temp

# In[22]:


hdf,hdf_dict = lu.load_hdf(fp+'MASTER/MASTERL1B_2190500_03_20201013_1804_1829_V01.hdf',all_values=True,i_subdata=0)


# In[25]:


hdf_dict['CalibratedData']


# ## 2016 flights

# In[38]:


hdf2,hdf2_dict = lu.load_hdf(fp+'MASTER/MASTERL1B_1662900_12_20160617_2204_2219_V01.hdf',all_values=True,i_subdata=0)


# In[39]:


hdf2_dict['DataSetHeader']


# # Plot out data

# ## 2020 data

# In[26]:


hdf['CalibratedData'].shape


# In[28]:


plt.figure()
plt.pcolor(hdf['PixelLongitude'],hdf['PixelLatitude'],hdf['CalibratedData'][6,:,:]-hdf['CalibratedData'][5,:,:])


# In[31]:


plt.figure()
plt.hist(hdf['CalibratedData'][6,:,:].flatten(),bins=50)
plt.yscale('log')


# ## 2016 data

# In[34]:


hdf2['CalibratedData'].shape


# In[51]:


fig = plt.figure()
plt.pcolor(hdf2['PixelLongitude'],hdf2['PixelLatitude'],hdf2['CalibratedData'][6,:,:],picker=True,pickradius=5)
#text = plt.gca().text(0,0, "", va="bottom", ha="left")
#def onpick(event):
#    print('onpick:',event.ind)
#    text.set_text('onpick:',event.ind)
#    return True
#fig.canvas.mpl_connect('pick_event', onpick)
#plt.show()


# In[40]:


plt.figure()
plt.pcolor(hdf2['PixelLongitude'],hdf2['PixelLatitude'],hdf2['CalibratedData'][6,:,:]-hdf2['CalibratedData'][5,:,:])


# In[42]:


plt.figure()
plt.pcolor(hdf2['PixelLongitude'],hdf2['PixelLatitude'],hdf2['CalibratedData'][6,:,:]-hdf2['CalibratedData'][7,:,:])
plt.colorbar()


# In[44]:


plt.figure()
plt.hist(hdf2['CalibratedData'][6,:,:].flatten()-hdf2['CalibratedData'][5,:,:].flatten(),bins=50)
plt.yscale('log')


# In[45]:


plt.figure()
plt.hist(hdf2['CalibratedData'][6,:,:].flatten()-hdf2['CalibratedData'][7,:,:].flatten(),bins=50)
plt.yscale('log')


# ## Plot out spectra

# By looking at points on map that seem near the fire. for hdf file #12 
# 
# For 2016 hdf2 point:
#  - lon: -120.025440
#  - lat: 34.51126

# In[52]:


lat_region = [34.507,34.515]
lon_region = [-120.02706,-120.02355]
lat_pnt = 34.51126
lon_pnt = -120.02544


# In[56]:


i_pnt = np.unravel_index(np.argmin(pow(lon_pnt-hdf2['PixelLongitude'],2)+pow(lat_pnt-hdf2['PixelLatitude'],2)),hdf2['PixelLatitude'].shape)


# In[57]:


i_pnt, hdf2['PixelLongitude'][i_pnt], hdf2['PixelLatitude'][i_pnt] 


# In[59]:


i_rg = np.where((hdf2['PixelLongitude']>lon_region[0]) & (hdf2['PixelLongitude']<lon_region[1]) &
         (hdf2['PixelLatitude']>lat_region[0]) & (hdf2['PixelLatitude']<lat_region[1]))


# In[62]:


hdf2.keys()


# In[201]:


wvl = np.array([0.466,0.506,0.546,0.586,0.660,0.716,0.756,0.804,0.872,0.912,0.952,1.6,1.656,1.712,1.766,1.818,1.876,1.924,1.974,
       2.074,2.158,2.206,2.256,2.316,2.384,4.0501,3.27,3.445,3.585,3.74,3.895,4.05,4.231,4.365,4.515,4.67,4.825,4.975,
       5.105,5.255,7.78,8.17,8.61,9.04,9.66,10.08,10.58,11.25,12.11,12.86])


# In[202]:


isort = np.argsort(wvl)


# In[203]:


isort


# In[204]:


scale = [float(r) for r in hdf2_dict['CalibratedData']['scale_factor'].split(', ')]
rad = np.array([hdf2['CalibratedData'][i,:,:]*s for i,s in enumerate(scale)])


# In[205]:


rad.shape


# In[98]:


import matplotlib.ticker as ticker


# In[102]:


def myLogFormat(y,pos):
    # Find the number of decimal places required
    decimalplaces = 2 #int(np.maximum(-np.log10(y),0))     # =0 for numbers >=1
    # Insert that number into a format string
    formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
    # Return the formatted tick label
    return formatstring.format(y)


# In[119]:


fig = plt.figure()

plt.plot(wvl[isort],rad[isort,i_pnt[0],i_pnt[1]],'.-')
plt.xscale('log')
plt.xticks([0.45,0.6,0.766,1.0,1.6,2.1,4.2,12.0])
#plt.set_scientific(False)
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Radiance [{}]'.format(hdf2_dict['CalibratedData']['units']))
plt.grid()
ax = plt.gca()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))
plt.title('MASTER flight over Santa Barbara fire - 2016-06-17')


# In[130]:


fig = plt.figure()

plt.plot(wvl[isort],rad[:,i_rg[0],i_rg[1]][isort,:],'.-')
plt.xscale('log')
plt.xticks([0.45,0.6,0.766,1.0,1.6,2.1,4.2,12.0])
#plt.set_scientific(False)
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Radiance [{}]'.format(hdf2_dict['CalibratedData']['units']))
plt.grid()
ax = plt.gca()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))
plt.title('MASTER flight over Santa Barbara fire region - 2016-06-17')


# In[134]:


SSI_str = '''2002.396
1894.1119
1804.2944
1732.247
1492.3435
1347.9291
1230.9342
1087.0133
954.6004
891.48627
817.17395
241.05524
225.26944
200.20187
175.09557
154.9544
136.57474
128.83186
120.36308
95.52779
79.23871
72.289474
67.94378
58.332928
55.806812
7.965488
17.595697
14.798581
12.533788
10.772283
9.262421
7.965488
6.5911274
5.816839
5.0909443
4.4618354
3.9004765
3.464939
3.0971522
2.7584012
0.6484663
0.0
0.42585844
0.34842965
0.27100083
0.22260782
0.18389341
0.14517902
0.10646461
0.08710741'''
SSI = np.array([float(s) for s in SSI_str.split('\n')])


# In[137]:


refl =  np.array([rad[i,:,:]/s for i,s in enumerate(SSI)])


# In[148]:


fig = plt.figure()
plt.gca().set_prop_cycle(color=[plt.cm.gist_ncar(k) for k in np.linspace(0, 1, len(i_rg[0]))])
plt.plot(wvl[isort],refl[:,i_rg[0],i_rg[1]][isort,:],'.-')

plt.xscale('log')
plt.xticks([0.45,0.6,0.766,1.0,1.6,2.1,4.2,12.0])
plt.ylim([0,1])
#plt.set_scientific(False)
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Reflectance [unit]'.format(hdf2_dict['CalibratedData']['units']))
plt.grid()
ax = plt.gca()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))
plt.title('MASTER flight over Santa Barbara fire region - 2016-06-17')
scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.gist_ncar)
scalarmap.set_array(np.linspace(0, 1, len(i_rg[0])))
plt.colorbar(scalarmap)


# In[150]:


wvl[9]


# In[166]:


fig = plt.figure()
plt.pcolor(hdf2['PixelLongitude'],hdf2['PixelLatitude'],hdf2['CalibratedData'][9,:,:])
plt.colorbar()
xup = np.where(rad[9,:,:]>150.0)
plt.plot(hdf2['PixelLongitude'][xup],hdf2['PixelLatitude'][xup],'xr')
xup2 = np.where(rad[11,:,:]>150.0)
plt.plot(hdf2['PixelLongitude'][xup2],hdf2['PixelLatitude'][xup2],'x',color='orange')


# In[163]:


hdf2['PixelLongitude'][xup]


# ## Make planck balckbody calculations

# From estimates of Wien's displacement
#   
# For fires at 800°C, 1073K, peak wavelength is 2700 nm  
# For a peak wavelength of 2100 nm, temperature should be 1379K -> which is right in the ballpark for forest fire temperatures

# In[152]:


h = 6.626e-34
c = 3.0e+8
k = 1.38e-23

def planck(wav, T):
    a = 2.0*h*c**2
    b = h*c/(wav*k*T)
    intensity = a/ ( (wav**5) * (np.exp(b) - 1.0) )
    return intensity


# In[153]:


wavelengths = np.arange(1e-9, 3e-6, 1e-9) 
intensity2100 = planck(wavelengths, 2100.)
intensity2700 = planck(wavelengths, 2700.)


# In[157]:


plt.figure()
plt.plot(wavelengths,intensity2100,'.b',label='2100 nm')
plt.plot(wavelengths,intensity2700,'.r',label='2700 nm')
plt.legend(frameon=False)
plt.grid()
plt.xlabel('Wavelength [m]')
plt.ylabel('Intensity')
plt.title('Planck distribution for blackbody')


# In[159]:


wvl[isort][9]


# In[160]:


wvl[9]


# In[170]:


from colorpy import ciexyz


# In[168]:


help(colorpy)


# In[171]:


import colour


# In[174]:


sd = colour.SpectralDistribution(rad[isort,i_pnt[0],i_pnt[1]])


# In[206]:


sd.wavelengths = wvl[isort]*1000.0


# In[207]:


sd


# In[212]:


sdi = sd.interpolate(colour.SpectralShape(start=465,end=870,interval=1))


# In[214]:


xyz = colour.sd_to_XYZ(sdi)


# In[246]:


colour.plotting.plot_single_sd(
    sd,
    modulate_colours_with_sd_amplitude=True,
    y_label='Relative Flux / $F_\\lambda$')
plt.plot(wvl[isort],rad[isort,i_pnt[0],i_pnt[1]])


# In[242]:


RGB = colour.XYZ_to_sRGB(xyz)


# In[244]:


colour.plotting.plot_single_colour_swatch(
    colour.plotting.ColourSwatch(
        'The Burn colour', colour.utilities.normalise_maximum(RGB)),
    text_parameters={'size': 'x-large'});


# In[ ]:




