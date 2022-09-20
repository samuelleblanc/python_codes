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
import plotting_utils as pu
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
vv = 'v3'
fp = getpath(name)


# # Load files

# ## Load 4STAR gas summary

# In[5]:


files = glob.glob(fp + '/gas_summary_v3/**[!4STARB]*gas_summary*.mat', recursive=True)


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
    gas[d]['time'] = pd.to_datetime([d]*len(gas[d]['tUTC'][:,0])).to_numpy() +                                   [np.timedelta64(timedelta(hours=tt)) for tt in gas[d]['tUTC'][:,0]]


# In[12]:


s['time'] = pd.to_datetime(s['day']).to_numpy() + [np.timedelta64(timedelta(hours=tt)) for tt in s['tUTC']]


# In[13]:


s['time']


# ### Filter out bad data and back interpolate for easier handling

# In[16]:


ibad = (s['no2DU']<0.005) | (s['no2DU']>1.0) | (s['no2resiDU']>0.00005)


# In[17]:


s.keys()


# In[18]:


from Sp_parameters import smooth


# In[19]:


s['no2DU'][ibad] = np.nan


# In[20]:


s_no2DU = smooth(s['no2DU'],1,nan=True,old=True)


# In[21]:


plt.figure()
plt.plot(s['no2DU'],s['no2resiDU'],'.')


# In[22]:


plt.figure()
plt.plot(s['no2resiDU'],'.')


# In[23]:


plt.figure()
plt.plot(s_no2DU,'.')


# ### Load previous 4STAR analysis

# In[22]:


files2 = glob.glob(fp + '/gas_summary_v1*/**/[!4STARB]*gas_summary*.mat', recursive=True)
files2.sort()
files2


# In[23]:


gas2 = {}
for f in files2:
    g2 = sio.loadmat(f)
    daystr = f.replace('/','_').split('_')[-3]
    gas2[daystr] = g2


# In[24]:


s2 = {'day':[],'time':[],'tUTC':[],'no2DU':[],'no2resiDU':[],'o3DU':[],'o3resiDU':[]}
kg = list(gas2.keys())
kg.sort()
for d in kg:
    for kk in ['tUTC','no2DU','no2resiDU','o3DU','o3resiDU']:
        s2[kk] = np.append(s2[kk],gas2[d][kk][:,0])
    s2['day'] = np.append(s2['day'],[d]*len(gas2[d][kk][:,0]))


# In[25]:


s2['time'] = pd.to_datetime(s2['day']).to_numpy() + [np.timedelta64(timedelta(hours=tt)) for tt in s2['tUTC']]


# ## Load Pandora

# ### load files from internet

# In[24]:


import requests


# In[25]:


url_base = 'http://data.pandonia-global-network.org/MountainViewCA/Pandora34s1/L2/'
url_O3 = url_base+'Pandora34s1_MountainViewCA_L2Tot_rout0p1-7.txt'
url_NO2 = url_base+'Pandora34s1_MountainViewCA_L2Tot_rnvs1p1-7.txt'
url_NO2_trop = url_base+'Pandora34s1_MountainViewCA_L2Trop_rnvh1p1-7.txt'


# In[26]:


ro3 = requests.get(url_O3,stream=True)
fo3 = []
for chunk in ro3.iter_content(chunk_size=1024): 
    fo3.append(chunk)


# In[27]:


rno2 = requests.get(url_NO2,stream=True)
fno2 = []
for chunk in rno2.iter_content(chunk_size=1024): 
    fno2.append(chunk)


# ### convert files to something easier to handle

# #### Ozone

# In[28]:


o3data = BytesIO(b''.join(fo3))


# In[29]:


o3data.seek(0)
o3header = []
for i in range(66):
    o3header.append(o3data.readline())
print(*o3header,sep='\n')


# In[30]:


col_name = ['Time','doy_jan2000','duration_time','sza','saz','lza','laz','o3_du','unc_o3_du','m_o3','diff_corr','qa_o3','sumi2_dq1',
            'sumi_dq2','so2_du','unc_so2_du','m_so2','diff_corr_so2','qa_no2','sumi2_so2_dq1','sumi2_so2_dq2',
            'fit_resi','resi_rms','expmeas_resi_rms','expinst_resi_rms','mean_val','Pres_mbar','dat_proc',
            'cal_version','cal_val_date','L2fit_QA','sumi2_L2fit_dq1','sumi2_L2fit_dq2','L1_QA','sumi2_L1_dq1',
            'sumi2_L1_dq2','effT_wav','resi_straylight','L1_wv_shift','wv_shift_fit','int','darkcount','pos_filter1','pos_filter2']


# In[31]:


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


# In[32]:


o3data.seek(0)
pdo3 = pd.read_csv(o3data,encoding='unicode_escape',header=66,delimiter=' ',names=col_name)


# In[33]:


pdo3


# In[34]:


pdo3['datetime'] = pd.to_datetime(pdo3['Time'])


# #### NO2

# In[35]:


no2data = BytesIO(b''.join(fno2))


# In[36]:


no2data.seek(0)
no2header = []
for i in range(59):
    no2header.append(no2data.readline())
print(*no2header,sep='\n')


# In[37]:


col_name_no2 = ['Time','doy_jan2000','duration_time','sza','saz','lza','laz',
                'no2_du','unc_no2_du','m_no2','diff_corr_no2','qa_no2','sumi2_no2_dq1','sumi_no2_dq2',
                'fit_resi','resi_rms','expmeas_resi_rms','expinst_resi_rms','mean_val','Pres_mbar','dat_proc',
            'cal_version','cal_val_date','L2fit_QA','sumi2_L2fit_dq1','sumi2_L2fit_dq2','L1_QA','sumi2_L1_dq1',
            'sumi2_L1_dq2','effT_wav','resi_straylight','L1_wv_shift','wv_shift_fit','int','darkcount','pos_filter1','pos_filter2']


# In[38]:


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


# In[39]:


no2data.seek(0)
pdno2 = pd.read_csv(no2data,encoding='unicode_escape',header=59,delimiter=' ',names=col_name_no2)


# In[40]:


pdno2['datetime'] = pd.to_datetime(pdno2['Time'])


# In[41]:


pdno2


# # Plot out data

# ## Plot time series

# In[42]:


plt.figure()
plt.plot(s['time'],s['no2DU'],'.')


# In[45]:


plt.figure()
plt.plot(s['time'],s['no2DU'],'.',label='4STAR')
plt.plot(s2['time'],s2['no2DU'],'+',label='4STAR pre')
plt.legend()


# In[44]:


plt.figure()
plt.plot(pdno2['datetime'],pdno2['no2_du'],'.',label='Pandora')
plt.plot(s['time'],s['no2DU'],'x',label='4STAR')
plt.ylim(0,1)
plt.legend()


# In[113]:


ddd = s['time'][0]


# In[116]:


ddd -  datetime.datetime(2018,1,1)


# In[51]:


plt.figure()
plt.plot(pdo3['datetime'],pdo3['o3_du'],'.',label='Pandora')
plt.plot(s['time'],s['o3DU'],'.',label='4STAR')
plt.ylim(0,500)
plt.legend()


# ## Plot each day of measurement in new figure

# In[45]:


def get_unixtime(dt64):
    return dt64.astype('datetime64[ms]').astype('float')/1000.0


# In[46]:


unix_t = get_unixtime(s['time'])


# In[47]:


days_since_2018 = (unix_t-get_unixtime(np.datetime64('2018-01-01T00:00:00.0')))/60.0/60.0/24.0


# In[48]:


t0 = pdno2['datetime'].dt.tz_localize(None).to_numpy()


# In[49]:


pandora_days_since_2018 = (get_unixtime(t0)-get_unixtime(np.datetime64('2018-01-01T00:00:00.0')))/60.0/60.0/24.0


# In[50]:


meas_days = np.unique(days_since_2018.astype(int))


# In[51]:


len(meas_days)


# In[52]:


pdd = pdno2['datetime'][0]


# In[53]:


print(pdd.date())


# In[57]:


meas_days


# In[58]:


gas['20200624'].keys()


# In[59]:


day_key = '20200624'
gas[day_key]['time'] = pd.to_datetime([day_key]*len(gas[day_key]['no2DU'][:,0])).to_numpy() +                                      [np.timedelta64(timedelta(hours=tt)) for tt in gas[day_key]['tUTC'][:,0]]


# In[169]:


plt.figure()
plt.plot(gas['20200624']['time'],gas['20200624']['no2DU'],'x')


# In[ ]:





# In[178]:


day = 905
plt.figure()
ip = pandora_days_since_2018.astype(int) == day
iss = days_since_2018.astype(int) == day
plt.plot(pdno2['datetime'][ip],pdno2['no2_du'][ip],'.',label='Pandora')
plt.plot(s['time'][iss],s['no2DU'][iss],'x',label='4STAR')
plt.ylim(0,1)
plt.legend()
plt.title('NO2 roof top on {}'.format(pdno2['datetime'][ip].iloc[0].date()))


# In[159]:


s['no2DU'][iss][~np.isnan(s['no2DU'][iss])]


# In[152]:


for day in meas_days:
    plt.figure()
    ip = pandora_days_since_2018.astype(int) == day
    iss = days_since_2018.astype(int) == day
    if not any(ip): continue
    plt.plot(pdno2['datetime'][ip],pdno2['no2_du'][ip],'.',label='Pandora')
    plt.plot(s['time'][iss],s['no2DU'][iss],'x',label='4STAR')
    plt.ylim(0,1)
    plt.legend()
    plt.title('NO2 roof top on {}'.format(pdno2['datetime'][ip].iloc[0].date()))


# In[68]:


np.max(gas['20200624']['no2DU'])


# In[69]:


np.min(gas['20200624']['no2DU'])


# In[90]:


g = '20200624'
ii = (pdno2['datetime'].dt.tz_localize(None).to_numpy()>gas[g]['time'][0]-8600000000000) &         (pdno2['datetime'].dt.tz_localize(None).to_numpy()<gas[g]['time'][-1]+8600000000000)


# In[74]:


for g in gas:
    plt.figure()
    igood = (gas[g]['no2DU'][:,0] > 0.025) & (gas[g]['no2resiDU'][:,0]<0.05)
    plt.plot(gas[g]['time'],gas[g]['no2DU'],'.',color='lightgrey',label='4STAR all')
    plt.plot(gas[g]['time'][igood],gas[g]['no2DU'][igood],'.',label='4STAR good')
    ii = (pdno2['datetime'].dt.tz_localize(None).to_numpy()>gas[g]['time'][0]-12600000000000) &         (pdno2['datetime'].dt.tz_localize(None).to_numpy()<gas[g]['time'][-1]+12600000000000)
    plt.plot(pdno2['datetime'][ii],pdno2['no2_du'][ii],'x',label='pandora')
    plt.xlim([gas[g]['time'][0]-12600000000000,gas[g]['time'][-1]+12600000000000])
    plt.ylim(0.0,1.6)
    plt.legend()
    
    plt.ylabel('NO2 [DU]')
    plt.title(g)


# In[89]:


no2_star = []
no2_pand = []
titles = []
for g in gas:
    igood = (gas[g]['no2DU'][:,0] > 0.025) & (gas[g]['no2resiDU'][:,0]<0.05)
    ii = (pdno2['datetime'].dt.tz_localize(None).to_numpy()>gas[g]['time'][0]-12600000000000) &         (pdno2['datetime'].dt.tz_localize(None).to_numpy()<gas[g]['time'][-1]+12600000000000)
    no2_star.append(gas[g]['no2DU'][igood])
    no2_pand.append(pdno2['no2_du'][ii].to_numpy())
    titles.append(g)


# In[95]:


no2_star2 = [o.flatten() for o in no2_star]
no2_pand2 = [o.flatten() for o in no2_pand]


# In[126]:


plt.figure()
gr = plt.cm.autumn_r
dummy_n = np.arange(0,len(no2_star2))*0.0+200.0
dummy_n[-1] = 205.0
bp = plt.boxplot(no2_star2,vert=True,positions=np.arange(0,len(no2_star2)),
                   showfliers=False,widths=0.3,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(gr,bp,dummy_n,whisker_color='grey')
grp = plt.cm.winter
bpp = plt.boxplot(no2_pand2,vert=True,positions=np.arange(0,len(no2_star2))+0.3,
                   showfliers=False,widths=0.3,showmeans=True,patch_artist=True)
pu.set_box_whisker_color(grp,bpp,dummy_n,whisker_color='grey')
plt.gca().set_xticks(np.arange(0,len(no2_star2)))
plt.gca().set_xticklabels(titles,rotation=90)
plt.grid()
plt.ylabel('NO2 [DU]')
plt.title('Daily 4STAR and Pandora - Ames roof top measurements')
plt.legend([bp['boxes'][0],bpp['boxes'][0],bp['means'][0],bp['medians'][0],bp['boxes'][0],bp['whiskers'][0]],
               ['4STAR','Pandora','Mean','Median','25% - 75%','min-max'],
               frameon=True,loc=0,numpoints=1)
plt.tight_layout(rect=(0,0.05,1,1))
plt.savefig(fp + '/gas_summary_v3/daily_4STAR_pandora_NO2_{}.png',dpi=600,transparent=True)


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

# In[48]:


def get_unixtime(dt64):
    return dt64.astype('datetime64[ms]').astype('float')/1000.0


# In[49]:


t = s['time'][10000]


# In[52]:


get_unixtime(t)


# In[53]:


t0 = pdno2['datetime'].dt.tz_localize(None).to_numpy()


# In[54]:


get_unixtime(t0[10000])


# In[55]:


from tqdm import tqdm_notebook as tqdm 


# In[71]:


time_span = 300 #in seconds
match_pandora_no2 = []
match_4star_no2 = []

unix_t0 = get_unixtime(t0)
unix_t = get_unixtime(s['time'])

pbar = tqdm(total = len(unix_t))
nstep = 50

for t in unix_t[::nstep]:
    imatch = (unix_t0<=t+time_span) & (unix_t0>=t-time_span)
    if any(imatch):
        match_pandora_no2.append(np.nanmean(pdno2['no2_du'][imatch]))
        match_4star_no2.append(np.nanmean(s['no2DU'][(unix_t<=t+time_span) & (unix_t>=t-time_span)]))
    else:
        match_pandora_no2.append(np.nan)
        match_4star_no2.append(np.nan)
    pbar.update(nstep)
match_pandora_no2 = np.array(match_pandora_no2)
match_4star_no2 = np.array(match_4star_no2)


# In[117]:


days_since_2020 = (unix_t-get_unixtime(np.datetime64('2020-01-01T00:00:00.0')))/60.0/60.0/24.0


# In[79]:


plt.figure()
plt.plot(match_pandora_no2,match_4star_no2,'.',label='averaged 5minute\n'+pu.stats_label(match_pandora_no2,match_4star_no2))
plt.plot([0,1],[0,1],'--',color='lightgrey',label='1:1')
pu.plot_lin(match_pandora_no2,match_4star_no2,x_err=match_pandora_no2*0.2,y_err=match_4star_no2*0.2,use_method='york')
plt.legend()
plt.scatter(match_pandora_no2,match_4star_no2,s=20,c=days_since_2020[::nstep],marker='o',zorder=10)
plt.colorbar(label='Days since 2020-01-01 [day]')
plt.ylabel('4STAR NO2 [DU]')
plt.xlabel('Pandora NO2 [DU]')
plt.title('NASA Ames rooftop comparison of 4STAR and Pandora NO2\nfrom 2020 to 2022 - version {}'.format(vv))
plt.ylim(-0.05,1)
plt.xlim(-0.05,1)
plt.savefig(fp+'4STAR_pandora_NO2_comparison_{}.png'.format(vv),dpi=600,transparent=True)


# ## one to one plot, matching measurement times

# ### Ozone

# In[80]:


spd = pd.DataFrame(s)
spd['datetime'] = pd.to_datetime(s['time'],utc=True)
spd2 = spd.sort_values(by='datetime')


# In[81]:


fullpd = pd.merge_asof(spd2,pdo3,direction='nearest')


# In[82]:


fullpd['o3DU'][fullpd['o3DU']<0.0] = np.nan


# In[83]:


fullpd['o3_du'][fullpd['qa_o3']>100] = np.nan


# In[84]:


fullpd['o3_du'] = fullpd['o3_du']/10.0


# In[85]:


plt.figure()
plt.plot(fullpd['datetime'],fullpd['o3_du'],'.',label='Pandora')
plt.plot(fullpd['datetime'],fullpd['o3DU'],'.',label='4STAR')
plt.ylim(0,500)
plt.legend()


# In[88]:


import plotting_utils as pu
from sklearn.metrics import mean_squared_error, mean_absolute_error


# In[89]:


fl = np.isfinite(fullpd['o3_du']) & np.isfinite(fullpd['o3DU'])


# In[90]:


r = st.spearmanr(fullpd['o3_du'],fullpd['o3DU'],nan_policy='omit')


# In[91]:


r.correlation


# In[93]:


rmse = mean_squared_error(fullpd['o3_du'][fl],fullpd['o3DU'][fl],squared=True)
rmse


# In[94]:


mae = mean_absolute_error(fullpd['o3_du'][fl],fullpd['o3DU'][fl])
mae


# In[95]:


import importlib
importlib.reload(pu)


# In[96]:


from Sp_parameters import doublenanmask, nanmasked
from scipy import odr
from linfit import linfit


# In[97]:


x = fullpd['o3_du']
y = fullpd['o3DU']
x_err = fullpd['unc_o3_du']
y_err = fullpd['o3resiDU']


# In[98]:


xn,yn,mask = doublenanmask(x,y,return_mask=True)


# In[99]:


model = odr.Model(lin)
dat = odr.RealData(xn,yn,sx=x_err[mask],sy=y_err[mask])


# In[100]:


c,cm = linfit(xn,yn)
p = np.array([c[1],c[0]])


# In[101]:


imask = mask & (y_err>0)


# In[102]:


ri = np.corrcoef(x_err[imask],y_err[imask])[0,1]


# In[103]:


a_bivar, b_bivar, S, cov = pu.bivariate_fit(xn,yn,x_err[imask],y_err[imask],b0=p[1],ri=ri**2)


# In[104]:


a_bivar


# In[105]:


b_bivar


# In[106]:


S


# In[107]:


np.sqrt(cov[1,1])*b_bivar


# In[108]:


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
plt.savefig(fp+'Winter_2022/plots/4STAR_to_Pandora_O3_{}.png'.format(vv),dpi=600,transparent=True)


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




