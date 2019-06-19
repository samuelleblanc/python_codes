
# coding: utf-8

# # Info
# Name:  
# 
#     Gas_sepctra_insformation_content
#     
# Purpose:  
# 
#     Loads HiTRAN for NO2 and Ozone, and uses different fwhm to subsample and calculate the information content
#     Adds random noise to sampling
#   
# Input:
# 
#     None
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
#     - load_utils.py : for loading OMI HDF5 files
#     - matplotlib
#     - numpy
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - 
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-06-14
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
from scipy import integrate


# In[2]:


import hapi


# In[3]:


get_ipython().magic(u'matplotlib notebook')


# In[4]:


fp =getpath('4STAR_gas')


# # Load files

# ## Using the HAPI hitran python API

# In[24]:


hapi.getHelp('tutorial')


# In[26]:


hapi.getHelp('data')


# In[27]:


hapi.db_begin('data')


# In[32]:


hapi.getHelp(hapi.fetch)


# In[39]:


hapi.fetch('NO2',0,1,29164,40798)


# ## Load the crosssection files directly

# In[23]:


no2_xsc = []
len(no2_xsc)


# In[ ]:


fname = fp+'no2/NO2_294.0_0.0_15002.0-42002.3_00.xsc'


# In[79]:


def read_xsc(fname):
    'Reading function for xsc files. full file path at fname, returns dict with numpy'
    import numpy as np
    no2_xsc = []
    no2h = []
    with open(fname) as f:
        if len(no2_xsc) < 1:
            no2h = f.readline()
        l = f.readlines()
        l = [tl[:-1] for tl in l]
        ''.join(l).strip().split(' ')
        no2_xsc = np.array([float(g) for g in ''.join(l).strip().split(' ')])
    header = {'molecule':no2h[0:20].strip(),'wvn_min':float(no2h[20:30]),'wvn_max':float(no2h[30:40]),'N':int(no2h[40:47]),
              'T':float(no2h[47:54]),'P':float(no2h[54:60]),'max_xsc':float(no2h[60:70]),'resolution':float(no2h[70:75]),
              'name':no2h[75:90].strip(),'na':no2h[90:94],'broadener':no2h[94:97].strip(),'ref':int(no2h[97:100]),
              'xsc':no2_xsc}
    header['wv'] = np.linspace(header['wvn_min'],header['wvn_max'],header['N'])
    header['nm'] = 1.0/header['wv'] *10.0*1000.0*1000.0
    return header


# In[80]:


no2 = read_xsc(fp+'no2/NO2_294.0_0.0_15002.0-42002.3_00.xsc')


# In[81]:


no2


# In[82]:


plt.figure()
plt.plot(no2['wv'],no2['xsc'])


# In[83]:


plt.figure()
plt.plot(no2['nm'],no2['xsc'])


# ## Make the Shannon information content function

# In[99]:


def sic(p):
    'Shannon Information of probability distribution p'
    n = len(p)
    pd = p / np.sum(p)
    sic = np.nansum(pd*np.log2(pd)) * (-1.0)
    return sic


# In[100]:


no2s = sic(no2['xsc'])


# In[101]:


no2s


# 
# ## Build the gaussian smoothing slit functions

# In[104]:


def gaussian( x, u , s):
    return 1./np.sqrt( 2. * np.pi * s**2 ) * np.exp(-1.0/2.0*((x-u)/s)**2.0)


# In[115]:


plt.figure()
plt.plot(no2['nm'],gaussian(no2['nm'],400.0,30.5)*0.0000000000000001,'.',label='guassian')
plt.plot(no2['nm'],no2['xsc'],'.',label='orig')
plt.plot(no2['nm'],gaussian(no2['nm'],400.0,30.5)*no2['xsc']*100.0,'.',label='filtered')
plt.legend(frameon=False)


# In[122]:


integrate.simps(gaussian(no2['nm'],400.0,30.5),no2['nm'])


# In[120]:


integrate.simps(gaussian(no2['nm'],400.0,30.5)*no2['xsc'],no2['nm'])


# In[125]:


i4 = np.argmin(abs(no2['nm']-400.0))
no2['xsc'][i4],no2['nm'][i4]


# In[130]:


def convolver(xsc,nm,fwhm):
    'function to go through the cross-section (xsc) and convolve a gaussian with full-width-half-max (fwhm), at the sampling grid (nm)'
    from scipy import integrate
    import numpy as np
    return np.array([integrate.simps(gaussian(nm,n,fwhm)*xsc,nm)*-1.0 for n in nm])

        
    


# In[129]:


plt.figure()
plt.plot(no2['nm'],no2['xsc'],'.',label='orig')
plt.plot(no2['nm'],convolver(no2['xsc'],no2['nm'],0.2),label='FWHM 0.2 nm')
plt.plot(no2['nm'],convolver(no2['xsc'],no2['nm'],0.4),label='FWHM 0.4 nm')
plt.plot(no2['nm'],convolver(no2['xsc'],no2['nm'],0.6),label='FWHM 0.6 nm')
plt.plot(no2['nm'],convolver(no2['xsc'],no2['nm'],0.8),label='FWHM 0.8 nm')
plt.plot(no2['nm'],convolver(no2['xsc'],no2['nm'],1.2),label='FWHM 1.2 nm')
plt.plot(no2['nm'],convolver(no2['xsc'],no2['nm'],1.5),label='FWHM 1.5 nm')
plt.plot(no2['nm'],convolver(no2['xsc'],no2['nm'],2.0),label='FWHM 2.0 nm')
plt.legend(frameon=False)


# In[150]:


fwhms = np.around(np.logspace(-1,0.75,20),2)


# In[151]:


fwhms


# In[152]:


subsamp = [convolver(no2['xsc'],no2['nm'],fw) for fw in fwhms]    


# In[153]:


sub_sic = [sic(s) for s in subsamp]


# In[161]:


plt.figure()
plt.plot(no2['nm'],np.array(subsamp).T)


# In[165]:


plt.figure() 
plt.plot(no2['nm'],subsamp[-1])
sic(subsamp[-1]),sic(no2['xsc'])


# In[159]:


plt.figure()
plt.plot(fwhms, sub_sic/sic(no2['xsc']),'+-')


# In[157]:


sic(no2['xsc'])


# In[156]:


sic(gaussian(no2['nm'],400.0,2.0))


# ## Add some white noise

# In[135]:


snr = 1.0/1000.0
np.random.rand(len(no2['nm']))*no2['max_xsc']*snr


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# # Plotting
# Present some fo the early plots here
