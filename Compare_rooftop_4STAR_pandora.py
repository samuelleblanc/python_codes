#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Describe the details ...
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-03-17
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
import glob


# In[361]:


from datetime import datetime,timedelta
import scipy.stats as st


# In[61]:


import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO, BytesIO
else:
    from io import StringIO, BytesIO


# In[9]:


name = 'rooftop'
vv = 'v1'
fp = getpath(name)


# # Load files

# ## Load 4STAR gas summary

# In[12]:


files = glob.glob(fp + '/**/[!4STARB]*gas_summary*.mat', recursive=True)


# In[13]:


files


# In[22]:


gas = {}
for f in files:
    g = sio.loadmat(f)
    daystr = f.replace('/','_').split('_')[-3]
    gas[daystr] = g


# In[25]:


gas['20220208'].keys()


# In[159]:


gas['20220208']['tUTC']


# In[186]:


kg = list(gas.keys())
kg.sort()


# In[211]:


s = {'day':[],'time':[],'tUTC':[],'no2DU':[],'no2resiDU':[],'o3DU':[],'o3resiDU':[]}
for d in kg:
    for kk in ['tUTC','no2DU','no2resiDU','o3DU','o3resiDU']:
        s[kk] = np.append(s[kk],gas[d][kk][:,0])
    s['day'] = np.append(s['day'],[d]*len(gas[d][kk][:,0]))


# In[266]:


s['time'] = pd.to_datetime(s['day']).to_numpy() + [np.timedelta64(timedelta(hours=tt)) for tt in s['tUTC']]


# In[267]:


s['time']


# ## Load Pandora

# ### load files from internet

# In[26]:


import requests


# In[38]:


url_base = 'http://data.pandonia-global-network.org/MountainViewCA/Pandora34s1/L2/'
url_O3 = url_base+'Pandora34s1_MountainViewCA_L2Tot_rout0p1-7.txt'
url_NO2 = url_base+'Pandora34s1_MountainViewCA_L2Tot_rnvs1p1-7.txt'
url_NO2_trop = url_base+'Pandora34s1_MountainViewCA_L2Trop_rnvh1p1-7.txt'


# In[41]:


ro3 = requests.get(url_O3,stream=True)
fo3 = []
for chunk in ro3.iter_content(chunk_size=1024): 
    fo3.append(chunk)


# In[43]:


rno2 = requests.get(url_NO2,stream=True)
fno2 = []
for chunk in rno2.iter_content(chunk_size=1024): 
    fno2.append(chunk)


# ### convert files to something easier to handle

# In[44]:


import pandas as pd


# #### Ozone

# In[93]:


o3data = BytesIO(b''.join(fo3))


# In[132]:


o3data.seek(0)
o3header = []
for i in range(66):
    o3header.append(o3data.readline())
print(*o3header,sep='\n')


# In[143]:


col_name = ['Time','doy_jan2000','measT','sza','saz','lza','laz','o3_du','unc_o3_du','m_o3','diff_corr','qa_o3','sumi2_dq1',
            'sumi_dq2','so2_du','unc_so2_du','m_so2','diff_corr_so2','qa_no2','sumi2_so2_dq1','sumi2_so2_dq2',
            'fit_resi','resi_rms','expmeas_resi_rms','expinst_resi_rms','mean_val','Pres_mbar','dat_proc',
            'cal_version','cal_val_date','L2fit_QA','sumi2_L2fit_dq1','sumi2_L2fit_dq2','L1_QA','sumi2_L1_dq1',
            'sumi2_L1_dq2','effT_wav','resi_straylight','L1_wv_shift','wv_shift_fit','int','darkcount','pos_filter1','pos_filter2']


# In[144]:


for i,h in enumerate(o3header[21:]):
    print(i,'\033[1m'+col_name[i]+'\033[0m',h)


# In[145]:


o3data.seek(0)
pdo3 = pd.read_csv(o3data,encoding='unicode_escape',header=66,delimiter=' ',names=col_name)


# In[146]:


pdo3


# In[147]:


pdo3['datetime'] = pd.to_datetime(pdo3['Time'])


# #### NO2

# In[149]:


no2data = BytesIO(b''.join(fno2))


# In[152]:


no2data.seek(0)
no2header = []
for i in range(59):
    no2header.append(no2data.readline())
print(*no2header,sep='\n')


# In[154]:


col_name_no2 = ['Time','doy_jan2000','measT','sza','saz','lza','laz',
                'no2_du','unc_no2_du','m_no2','diff_corr_no2','qa_no2','sumi2_no2_dq1','sumi_no2_dq2',
                'fit_resi','resi_rms','expmeas_resi_rms','expinst_resi_rms','mean_val','Pres_mbar','dat_proc',
            'cal_version','cal_val_date','L2fit_QA','sumi2_L2fit_dq1','sumi2_L2fit_dq2','L1_QA','sumi2_L1_dq1',
            'sumi2_L1_dq2','effT_wav','resi_straylight','L1_wv_shift','wv_shift_fit','int','darkcount','pos_filter1','pos_filter2']


# In[155]:


for i,h in enumerate(no2header[21:]):
    print('\033[1m'+col_name_no2[i]+'\033[0m',h)


# In[156]:


no2data.seek(0)
pdno2 = pd.read_csv(no2data,encoding='unicode_escape',header=59,delimiter=' ',names=col_name_no2)


# In[157]:


pdno2['datetime'] = pd.to_datetime(pdno2['Time'])


# In[158]:


pdno2


# # Plot out data

# ## Plot time series

# In[270]:


plt.figure()
plt.plot(pdno2['datetime'],pdno2['no2_du'],'.',label='Pandora')
plt.plot(s['time'],s['no2DU'],'.',label='4STAR')
plt.ylim(0,10)
plt.legend()


# In[273]:


plt.figure()
plt.plot(pdo3['datetime'],pdo3['o3_du'],'.',label='Pandora')
plt.plot(s['time'],s['o3DU'],'.',label='4STAR')
plt.ylim(0,500)
plt.legend()


# ## one to one plot, matching measurement times

# ### Ozone

# In[347]:


spd = pd.DataFrame(s)
spd['datetime'] = pd.to_datetime(s['time'],utc=True)
spd2 = spd.sort_values(by='datetime')


# In[320]:


fullpd = pd.merge_asof(spd2,pdo3,direction='nearest')


# In[323]:


fullpd['o3DU'][fullpd['o3DU']<0.0] = np.nan


# In[327]:


fullpd['o3_du'][fullpd['qa_o3']>10] = np.nan


# In[329]:


plt.figure()
plt.plot(fullpd['datetime'],fullpd['o3_du'],'.',label='Pandora')
plt.plot(fullpd['datetime'],fullpd['o3DU'],'.',label='4STAR')
plt.ylim(0,500)
plt.legend()


# In[337]:


import plotting_utils as pu
from sklearn.metrics import mean_squared_error, mean_absolute_error


# In[343]:


fl = np.isfinite(fullpd['o3_du']) & np.isfinite(fullpd['o3DU'])


# In[364]:


r = st.spearmanr(fullpd['o3_du'],fullpd['o3DU'],nan_policy='omit')


# In[367]:


r.correlation


# In[372]:


rmse = mean_squared_error(fullpd['o3_du'][fl],fullpd['o3DU'][fl],squared=True)
rmse


# In[375]:


mae = mean_absolute_error(fullpd['o3_du'][fl],fullpd['o3DU'][fl])
mae


# In[451]:


import importlib
importlib.reload(pu)


# In[387]:


from Sp_parameters import doublenanmask, nanmasked
from scipy import odr
from linfit import linfit


# In[388]:


x = fullpd['o3_du']
y = fullpd['o3DU']
x_err = fullpd['unc_o3_du']
y_err = fullpd['o3resiDU']


# In[389]:


xn,yn,mask = doublenanmask(x,y,return_mask=True)


# In[390]:


model = odr.Model(lin)
dat = odr.RealData(xn,yn,sx=x_err[mask],sy=y_err[mask])


# In[391]:


c,cm = linfit(xn,yn)
p = np.array([c[1],c[0]])


# In[417]:


imask = mask & (y_err>0)


# In[435]:


ri = np.corrcoef(x_err[imask],y_err[imask])[0,1]


# In[442]:


a_bivar, b_bivar, S, cov = pu.bivariate_fit(xn,yn,x_err[imask],y_err[imask],b0=p[1],ri=ri**2)


# In[443]:


a_bivar


# In[444]:


b_bivar


# In[445]:


S


# In[449]:


np.sqrt(cov[1,1])*b_bivar


# In[441]:


plt.figure()
plt.errorbar(fullpd['o3_du'],fullpd['o3DU'],xerr=fullpd['unc_o3_du'],yerr=fullpd['o3resiDU'],marker='.',ls='')
plt.plot(fullpd['o3_du'],fullpd['o3DU'],'.',
         label='all data\nR$_{{spearman}}$={:1.3f}\nRMSE={:1.3f}\nMAE={:1.3f}'.format(r.correlation,rmse,mae))
plt.plot([0,500],[0,500],'--',color='lightgrey',label='1:1')
pu.plot_lin(fullpd['o3_du'],fullpd['o3DU'],x_err=fullpd['unc_o3_du'],y_err=fullpd['o3resiDU'],use_method='york')
plt.legend()
plt.ylabel('4STAR Ozone [DU]')
plt.xlabel('Pandora Ozone [DU]')
plt.ylim(200,400)
plt.xlim(200,400)
plt.savefig(fp+'Winter_2022/plots/4STAR_to_Pandora_O3.png',dpi=600,transparent=True)


# In[401]:


plt.figure()
plt.errorbar(fullpd['o3_du'],fullpd['o3DU'],xerr=fullpd['unc_o3_du'],yerr=fullpd['o3resiDU'],marker='.',ls='')
plt.plot(fullpd['o3_du'],fullpd['o3DU'],'.',
         label='all data\nR$_{{spearman}}$={:1.3f}\nRMSE={:1.3f}\nMAE={:1.3f}'.format(r.correlation,rmse,mae))
plt.plot([0,500],[0,500],'--',color='lightgrey',label='1:1')
pu.plot_lin(fullpd['o3_du'],fullpd['o3DU'],x_err=fullpd['unc_o3_du'],y_err=fullpd['o3resiDU'],use_method='york')
plt.legend()
plt.ylabel('4STAR Ozone [DU]')
plt.xlabel('Pandora Ozone [DU]')
plt.ylim(200,400)
plt.xlim(200,400)
plt.savefig(fp+'Winter_2022/plots/4STAR_to_Pandora_O3.png',dpi=600,transparent=True)


# ### NO2

# In[348]:


spd = pd.DataFrame(s)
spd['datetime'] = pd.to_datetime(s['time'],utc=True)
spd2 = spd.sort_values(by='datetime')


# In[349]:


fullpdn = pd.merge_asof(spd2,pdno2,direction='nearest')


# In[350]:


fullpdn['no2DU'][fullpdn['no2DU']<0.0] = np.nan


# In[352]:


fullpdn['no2_du'][fullpdn['qa_no2']>10] = np.nan


# In[354]:


plt.figure()
plt.plot(fullpdn['datetime'],fullpdn['no2_du'],'.',label='Pandora')
plt.plot(fullpdn['datetime'],fullpdn['no2DU'],'.',label='4STAR')
plt.ylim(0,10)
plt.legend()


# In[355]:


import plotting_utils as pu


# In[356]:


fln = np.isfinite(fullpdn['no2_du']) & np.isfinite(fullpdn['no2DU'])


# In[359]:


plt.figure()
plt.plot(fullpdn['no2_du'],fullpdn['no2DU'],'.',
         label='all data\nR={:1.3f}'.format(np.corrcoef(fullpdn['no2_du'][fln],fullpdn['no2DU'][fln])[0,1]))
plt.plot([0,10],[0,10],'--',color='lightgrey',label='1:1')
pu.plot_lin(fullpdn['no2_du'],fullpdn['no2DU'])
plt.legend()
plt.ylabel('4STAR NO2 [DU]')
plt.xlabel('Pandora NO2 [DU]')
plt.ylim(0,10)
plt.xlim(0,10)

plt.savefig(fp+'Winter_2022/plots/4STAR_to_Pandora_NO2.png',dpi=600,transparent=True)


# In[453]:


plt.figure()
plt.plot(fullpdn['no2_du'],fullpdn['no2DU'],'.',
         label='all data\n'+pu.stats_label(fullpdn['no2_du'],fullpdn['no2DU']))
plt.plot([0,10],[0,10],'--',color='lightgrey',label='1:1')
pu.plot_lin(fullpdn['no2_du'],fullpdn['no2DU'],x_err=fullpdn['unc_no2_du'],y_err=fullpdn['no2resiDU'],use_method='york')
plt.legend()
plt.ylabel('4STAR NO2 [DU]')
plt.xlabel('Pandora NO2 [DU]')
plt.ylim(0,10)
plt.xlim(0,10)

plt.savefig(fp+'Winter_2022/plots/4STAR_to_Pandora_NO2.png',dpi=600,transparent=True)


# In[ ]:




