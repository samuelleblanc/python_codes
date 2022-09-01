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

# In[1]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'notebook')
import os
import glob


# In[2]:


from datetime import datetime,timedelta
import scipy.stats as st
import pandas as pd
import plotting_utils as pu


# In[3]:


import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO, BytesIO
else:
    from io import StringIO, BytesIO


# In[4]:


name = 'rooftop'
vv = 'v1'
fp = getpath(name)


# # Load files

# ## Load 4STAR gas summary

# In[5]:


files = glob.glob(fp + '/gas_summary_v2_post202208/**/[!4STARB]*gas_summary*.mat', recursive=True)


# In[6]:


files.sort()
files


# In[7]:


gas = {}
for f in files:
    g = sio.loadmat(f)
    daystr = f.replace('/','_').split('_')[-3]
    gas[daystr] = g


# In[8]:


gas['20210330'].keys()


# In[9]:


gas['20210330']['tUTC']


# In[10]:


kg = list(gas.keys())
kg.sort()


# In[11]:


s = {'day':[],'time':[],'tUTC':[],'no2DU':[],'no2resiDU':[],'o3DU':[],'o3resiDU':[]}
for d in kg:
    for kk in ['tUTC','no2DU','no2resiDU','o3DU','o3resiDU']:
        s[kk] = np.append(s[kk],gas[d][kk][:,0])
    s['day'] = np.append(s['day'],[d]*len(gas[d][kk][:,0]))


# In[12]:


s['time'] = pd.to_datetime(s['day']).to_numpy() + [np.timedelta64(timedelta(hours=tt)) for tt in s['tUTC']]


# In[13]:


s['time']


# In[14]:


s['no2DU'] = s['no2DU']/10.0
s['no2resiDU'] = s['no2resiDU']/10.0


# ### Filter out bad data and back interpolate for easier handling

# In[161]:


ibad = (s['no2DU']<0.05) | (s['no2DU']>1.0) 


# In[148]:


from Sp_parameters import smooth


# In[162]:


s['no2DU'][ibad] = np.nan


# In[163]:


s_no2DU = smooth(s['no2DU'],1,nan=True,old=True)


# In[164]:


plt.figure()
plt.plot(s_no2DU,'.')


# ### Load previous 4STAR analysis

# In[15]:


files2 = glob.glob(fp + '/gas_summary_v1*/**/[!4STARB]*gas_summary*.mat', recursive=True)
files2.sort()
files2


# In[16]:


gas2 = {}
for f in files2:
    g2 = sio.loadmat(f)
    daystr = f.replace('/','_').split('_')[-3]
    gas2[daystr] = g2


# In[17]:


s2 = {'day':[],'time':[],'tUTC':[],'no2DU':[],'no2resiDU':[],'o3DU':[],'o3resiDU':[]}
kg = list(gas2.keys())
kg.sort()
for d in kg:
    for kk in ['tUTC','no2DU','no2resiDU','o3DU','o3resiDU']:
        s2[kk] = np.append(s2[kk],gas2[d][kk][:,0])
    s2['day'] = np.append(s2['day'],[d]*len(gas2[d][kk][:,0]))


# In[18]:


s2['time'] = pd.to_datetime(s2['day']).to_numpy() + [np.timedelta64(timedelta(hours=tt)) for tt in s2['tUTC']]


# ## Load Pandora

# ### load files from internet

# In[19]:


import requests


# In[20]:


url_base = 'http://data.pandonia-global-network.org/MountainViewCA/Pandora34s1/L2/'
url_O3 = url_base+'Pandora34s1_MountainViewCA_L2Tot_rout0p1-7.txt'
url_NO2 = url_base+'Pandora34s1_MountainViewCA_L2Tot_rnvs1p1-7.txt'
url_NO2_trop = url_base+'Pandora34s1_MountainViewCA_L2Trop_rnvh1p1-7.txt'


# In[21]:


ro3 = requests.get(url_O3,stream=True)
fo3 = []
for chunk in ro3.iter_content(chunk_size=1024): 
    fo3.append(chunk)


# In[22]:


rno2 = requests.get(url_NO2,stream=True)
fno2 = []
for chunk in rno2.iter_content(chunk_size=1024): 
    fno2.append(chunk)


# ### convert files to something easier to handle

# #### Ozone

# In[23]:


o3data = BytesIO(b''.join(fo3))


# In[24]:


o3data.seek(0)
o3header = []
for i in range(66):
    o3header.append(o3data.readline())
print(*o3header,sep='\n')


# In[25]:


col_name = ['Time','doy_jan2000','duration_time','sza','saz','lza','laz','o3_du','unc_o3_du','m_o3','diff_corr','qa_o3','sumi2_dq1',
            'sumi_dq2','so2_du','unc_so2_du','m_so2','diff_corr_so2','qa_no2','sumi2_so2_dq1','sumi2_so2_dq2',
            'fit_resi','resi_rms','expmeas_resi_rms','expinst_resi_rms','mean_val','Pres_mbar','dat_proc',
            'cal_version','cal_val_date','L2fit_QA','sumi2_L2fit_dq1','sumi2_L2fit_dq2','L1_QA','sumi2_L1_dq1',
            'sumi2_L1_dq2','effT_wav','resi_straylight','L1_wv_shift','wv_shift_fit','int','darkcount','pos_filter1','pos_filter2']


# In[26]:


col_name_decription = {}
j = 0
header_start = False
for i,h in enumerate(o3header):
    if h.strip().startswith(b'----'): 
        header_start = ~header_start
        continue
    if header_start:
        print(j,'\033[1m'+col_name[j]+'\033[0m',h)
        col_name_decription[col_name[j]] = h
        j = j+1


# In[27]:


o3data.seek(0)
pdo3 = pd.read_csv(o3data,encoding='unicode_escape',header=66,delimiter=' ',names=col_name)


# In[28]:


pdo3


# In[29]:


pdo3['datetime'] = pd.to_datetime(pdo3['Time'])


# #### NO2

# In[30]:


no2data = BytesIO(b''.join(fno2))


# In[31]:


no2data.seek(0)
no2header = []
for i in range(59):
    no2header.append(no2data.readline())
print(*no2header,sep='\n')


# In[32]:


col_name_no2 = ['Time','doy_jan2000','duration_time','sza','saz','lza','laz',
                'no2_du','unc_no2_du','m_no2','diff_corr_no2','qa_no2','sumi2_no2_dq1','sumi_no2_dq2',
                'fit_resi','resi_rms','expmeas_resi_rms','expinst_resi_rms','mean_val','Pres_mbar','dat_proc',
            'cal_version','cal_val_date','L2fit_QA','sumi2_L2fit_dq1','sumi2_L2fit_dq2','L1_QA','sumi2_L1_dq1',
            'sumi2_L1_dq2','effT_wav','resi_straylight','L1_wv_shift','wv_shift_fit','int','darkcount','pos_filter1','pos_filter2']


# In[33]:


col_name_no2_decription = {}
j = 0
header_start = False
for i,h in enumerate(no2header):
    if h.strip().startswith(b'----'): 
        header_start = ~header_start
        continue
    if header_start:
        print(j,'\033[1m'+col_name_no2[j]+'\033[0m',h)
        col_name_no2_decription[col_name_no2[j]] = h
        j = j+1


# In[34]:


no2data.seek(0)
pdno2 = pd.read_csv(no2data,encoding='unicode_escape',header=59,delimiter=' ',names=col_name_no2)


# In[35]:


pdno2['datetime'] = pd.to_datetime(pdno2['Time'])


# In[36]:


pdno2


# # Plot out data

# ## Plot time series

# In[236]:


plt.figure()
plt.plot(s['time'],s['no2DU'],'.')


# In[132]:


plt.figure()
plt.plot(s['time'],s['no2DU'],'.',label='4STAR')
plt.plot(s2['time'],s2['no2DU'],'+',label='4STAR pre')
plt.legend()


# In[237]:


plt.figure()
plt.plot(pdno2['datetime'],pdno2['no2_du'],'.',label='Pandora')
plt.plot(s['time'],s['no2DU'],'.',label='4STAR')
plt.ylim(0,1)
plt.legend()


# In[43]:


plt.figure()
plt.plot(pdo3['datetime'],pdo3['o3_du'],'.',label='Pandora')
plt.plot(s['time'],s['o3DU'],'.',label='4STAR')
plt.ylim(0,500)
plt.legend()


# ## Match in time to compare each measurement

# ### Time delta fast version - but buggy?

# In[95]:


time_diff = pdno2['datetime'].dt.tz_localize(None).to_numpy()-np.append(s['time'],s['time'][np.zeros((len(pdno2['datetime'])-len(s['time']),1)).astype(int)])


# In[126]:


time_diff = time_diff.astype(float)/1E9 #convert to seconds


# In[130]:


imatch = time_diff<600


# In[165]:


plt.figure()
plt.plot(pdno2['no2_du'][imatch],s_no2DU[imatch[0:len(s['no2DU'])]],'.')


# In[117]:


imatch.sum()


# ### Match in time, slow loop version

# In[178]:


def get_unixtime(dt64):
    return dt64.astype('datetime64[ms]').astype('float')/1000.0


# In[168]:


t = s['time'][10000]


# In[180]:


get_unixtime(t)


# In[170]:


t0 = pdno2['datetime'].dt.tz_localize(None).to_numpy()


# In[179]:


get_unixtime(t0[10000])


# In[184]:


from tqdm import tqdm_notebook as tqdm 


# In[199]:


np.diff(unix_t[190000:191000:100])


# In[205]:


time_span = 600 #in seconds
match_pandora_no2 = []
match_4star_no2 = []

unix_t0 = get_unixtime(t0)
unix_t = get_unixtime(s['time'])

pbar = tqdm(total = len(unix_t))

for t in unix_t[::100]:
    imatch = (unix_t0<=t+time_span) & (unix_t0>=t-time_span)
    if any(imatch):
        match_pandora_no2.append(np.nanmean(pdno2['no2_du'][imatch]))
        match_4star_no2.append(np.nanmean(s_no2DU[(unix_t<=t+time_span) & (unix_t>=t-time_span)]))
    else:
        match_pandora_no2.append(np.nan)
        match_4star_no2.append(np.nan)
    pbar.update(100)
match_pandora_no2 = np.array(match_pandora_no2)
match_4star_no2 = np.array(match_4star_no2)


# In[228]:


days_since_2020 = (unix_t-get_unixtime(np.datetime64('2020-01-01T00:00:00.0')))/60.0/60.0/24.0


# In[234]:


plt.figure()
plt.plot(match_pandora_no2,match_4star_no2,'.',label='averaged 10minute\n'+pu.stats_label(match_pandora_no2,match_4star_no2))
plt.plot([0,1],[0,1],'--',color='lightgrey',label='1:1')
pu.plot_lin(match_pandora_no2,match_4star_no2,x_err=match_pandora_no2*0.2,y_err=match_4star_no2*0.2,use_method='york')
plt.legend()
plt.scatter(match_pandora_no2,match_4star_no2,s=20,c=days_since_2020[::100],marker='o',zorder=10)
plt.colorbar(label='Days since 2020-01-01 [day]')
plt.ylabel('4STAR NO2 [DU]')
plt.xlabel('Pandora NO2 [DU]')
plt.title('NASA Ames rooftop comparison of 4STAR and Pandora NO2\nfrom 2020 to 2022')
plt.ylim(-0.05,1)
plt.xlim(-0.05,1)
plt.savefig(fp+'4STAR_pandora_NO2_comparison.png',dpi=600,transparent=True)


# ## one to one plot, matching measurement times

# ### Ozone

# In[48]:


spd = pd.DataFrame(s)
spd['datetime'] = pd.to_datetime(s['time'],utc=True)
spd2 = spd.sort_values(by='datetime')


# In[49]:


fullpd = pd.merge_asof(spd2,pdo3,direction='nearest')


# In[50]:


fullpd['o3DU'][fullpd['o3DU']<0.0] = np.nan


# In[51]:


fullpd['o3_du'][fullpd['qa_o3']>100] = np.nan


# In[52]:


fullpd['o3_du'] = fullpd['o3_du']/10.0


# In[53]:


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

# In[40]:


spd = pd.DataFrame(s)
spd['datetime'] = pd.to_datetime(s['time'],utc=True)
spd2 = spd.sort_values(by='datetime')


# In[41]:


fullpdn = pd.merge_asof(spd2,pdno2,direction='nearest')


# In[52]:



fullpdn['no2DU'][fullpdn['no2DU']<0.02] = np.nan
fullpdn['no2DU'][fullpdn['no2DU']>0.8] = np.nan


# In[43]:


fullpdn['no2_du'][fullpdn['qa_no2']>10] = np.nan


# In[44]:


plt.figure()
plt.plot(fullpdn['datetime'],fullpdn['no2_du'],'.',label='Pandora')
plt.plot(fullpdn['datetime'],fullpdn['no2DU'],'.',label='4STAR')
plt.ylim(0,10)
plt.legend()


# In[54]:


fln = np.isfinite(fullpdn['no2_du']) & np.isfinite(fullpdn['no2DU'])


# In[57]:


plt.figure()
plt.plot(fullpdn['no2_du'],fullpdn['no2DU'],'.',
         label='all data\nR={:1.3f}'.format(np.corrcoef(fullpdn['no2_du'][fln],fullpdn['no2DU'][fln])[0,1]))
plt.plot([0,10],[0,10],'--',color='lightgrey',label='1:1')
pu.plot_lin(fullpdn['no2_du'],fullpdn['no2DU'])
plt.legend()
plt.ylabel('4STAR NO2 [DU]')
plt.xlabel('Pandora NO2 [DU]')
plt.ylim(0,0.8)
plt.xlim(0,0.8)

#plt.savefig(fp+'Winter_2022/plots/4STAR_to_Pandora_NO2.png',dpi=600,transparent=True)


# In[55]:


plt.figure()
plt.plot(fullpdn['no2_du'],fullpdn['no2DU'],'.',
         label='all data\n'+pu.stats_label(fullpdn['no2_du'],fullpdn['no2DU']))
plt.plot([0,10],[0,10],'--',color='lightgrey',label='1:1')
pu.plot_lin(fullpdn['no2_du'],fullpdn['no2DU'],x_err=fullpdn['unc_no2_du'],y_err=fullpdn['no2resiDU'],use_method='york')
plt.legend()
plt.ylabel('4STAR NO2 [DU]')
plt.xlabel('Pandora NO2 [DU]')
plt.ylim(0,2)
plt.xlim(0,2)

#plt.savefig(fp+'Winter_2022/plots/4STAR_to_Pandora_NO2.png',dpi=600,transparent=True)


# In[ ]:




