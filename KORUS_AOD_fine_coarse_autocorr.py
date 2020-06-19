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


# In[287]:


import map_utils as mu
from scipy import interpolate
import math


# In[3]:


from linfit import linfit
import Sun_utils as su


# In[4]:


get_ipython().magic(u'matplotlib notebook')


# In[5]:


fp =getpath('KORUS')


# # Load files

# Load the KORUS 4STAR AOD ict files for better handling

# In[6]:


ar = hs.loadmat(fp+'/aod_ict/all_aod_KORUS_R2_ict.mat')


# In[7]:


ka = ar.keys()


# In[8]:


ka.sort()


# In[9]:


ka


# In[10]:


nwl = ka[0:17]


# In[11]:


nwl


# In[12]:


nm = [380.0,452.0,501.0,520.0,532.0,550.0,606.0,620.0,675.0,781.0,865.0,1020.0,1040.0,1064.0,1236.0,1559.0,1627.0]


# In[325]:


days = ['20160501','20160503','20160504','20160506','20160510','20160511',
        '20160512','20160516','20160517','20160519','20160521','20160524',
        '20160526','20160529','20160601','20160602','20160604','20160608',
        '20160609','20160614']


# In[326]:


doys = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in days]


# In[327]:


doys


# In[336]:


fdoys = np.array(doys)


# In[343]:


ar['doys'] = fdoys[ar['days'].astype(int)]+ar['Start_UTC']/24.0


# In[344]:


ar['doys']


# # Run analysis and prepare variables
# Do some of the calculations to the data here

# In[126]:


fl1 = ar['days']==ar['days'][0]


# In[127]:


fl1.shape


# In[128]:


fl = (ar['fl_QA']==0) & (np.isfinite(ar['AOD0501'])) 


# In[129]:


fl1.shape


# ## Calculate the Angstrom Exponent

# In[88]:


nwl,nm


# In[89]:


aodrr = np.array([ar[n] for n in nwl])


# In[90]:


aodrr.shape


# In[91]:


angs = su.calc_angs(ar['Start_UTC'],np.array(nm[1:11]),aodrr[1:11,:])


# ## Calculate the fine mode fraction

# In[345]:


fmf = su.sda(aodrr[1:13,:],np.array(nm[1:13])/1000.0)


# In[346]:


fmf.keys()


# In[348]:


fmf['tauc'].shape, ar['GPS_Alt'].shape


# ## Subset the level legs

# In[140]:


def running_std(x,n):
    'Function to do a running standard deviation on array (x) with window size (n)'
    q = x**2
    q = np.convolve(q, np.ones((n, )), mode="same")
    s = np.convolve(x, np.ones((n, )), mode="same")
    o = (q-s**2/n)/float(n-1)
    return o 


# In[141]:


nbox = 20


# In[142]:


std_alt = running_std(ar['GPS_Alt'][fl],nbox)


# In[143]:


std_alt.shape


# In[144]:


ar['GPS_Alt'][fl].shape


# In[145]:


f_level = np.where(std_alt<5.0)[0]


# In[146]:


std_alt1 = running_std(ar['GPS_Alt'][fl1],nbox)


# In[147]:


f_level1 = np.where(std_alt1<5.0)[0]


# In[149]:


ar['Start_UTC'][fl1][f_level1]


# In[151]:


plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(ar['Start_UTC'][fl1],ar['GPS_Alt'][fl1],'.')
plt.plot(ar['Start_UTC'][fl1][f_level1],ar['GPS_Alt'][fl1][f_level1],'r.')


ax2 = plt.subplot(2,1,2,sharex=ax1)
plt.plot(ar['Start_UTC'][fl1],std_alt1,'.')
plt.plot(ar['Start_UTC'][fl1][f_level1],std_alt[f_level1],'r.')
plt.ylim(0,100)


# In[152]:


plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(ar['Start_UTC'][fl],ar['GPS_Alt'][fl],'.')
plt.plot(ar['Start_UTC'][fl][f_level],ar['GPS_Alt'][fl][f_level],'r.')


ax2 = plt.subplot(2,1,2,sharex=ax1)
plt.plot(ar['Start_UTC'][fl],std_alt,'.')
plt.plot(ar['Start_UTC'][fl][f_level],std_alt[[f_level]],'r.')
plt.ylim(0,100)


# ## Seperate each of the level legs into distinct segments

# In[183]:


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


# In[324]:


ar['days']


# In[349]:


vals = {'utc':ar['Start_UTC'][fl],'alt':ar['GPS_Alt'][fl],'lat':ar['Latitude'][fl],'lon':ar['Longitude'][fl],
        'aod0500':ar['AOD0501'][fl],'aod1040':ar['AOD1040'][fl],'AE':angs[fl],'doys':ar['doys'][fl],
        'aod_fine':fmf['tauf'][fl],'aod_coarse':fmf['tauc'][fl]}


# In[350]:


dvals = get_segments(f_level,vals,nsep=100)


# In[351]:


dvals.keys()


# In[352]:


for n in dvals['utc']:
    try:
        print (n[-1]-n[0])*60.0
    except:
        print np.nan


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

# In[353]:


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


# In[354]:


ddv = get_distances(dvals)


# In[355]:


dvals['cumdist']


# In[356]:


dvals.keys()


# ## Calculate the autocorrelation of AOD with respect to distance

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

# In[357]:


# for met times
dt1 = ['20160501','20160516']
dt2 = ['20160517','20160522']
dt3 = ['20160523','20160531']
dt4 = ['20160601','20160607']


# In[360]:


t1 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt1]
t2 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt2]
t3 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt3]
t4 = [datetime(int(d[0:4]),int(d[4:6]),int(d[6:8])).timetuple().tm_yday for d in dt4]


# In[361]:


# limits of DOY for each of the met times
t1,t2,t3,t4


# In[514]:


#altitude limits in m
z1 = [0.0, 1000.0] 
z2 = [1000.0, 3000.0]
z3 = [3000.0, 15000.0]


# ### Test out Shinozuka & Redemann autocorrelation 

# In[263]:


dvals['cumdist'][2]


# In[363]:


dvals.keys()


# In[ ]:


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


# In[316]:


len(cr)


# In[321]:


plt.figure()
for corr in cr:
   # plt.plot(corr_ks,corr['AE']**2.0,'x-')
    plt.plot(corr_ks,corr['aod1040']**2.0,'s-')
    plt.plot(corr_ks,corr['aod0500']**2.0,'v-')
plt.xscale('log')
plt.ylim(0,1)


# In[574]:


types = ['all','t1','t2','t3','t4','z1','z2','z3']


# In[575]:


corr_ks =[0.1,0.25,0.5,0.75,1.0,1.5,2.0,3.0,5.0,7.5,10.0,12.5,15.0,20.0,
          25.0,30.0,35.0,40.0,50.0,60.0,75.0,100.0,150.0,200.0] 
corr_all = [[[{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)],
             [{'k':k,'aod1040':[],'aod0500':[],'AE':[],'aod_fine':[],'aod_coarse':[]} for i,k in enumerate(corr_ks)]] \
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


# In[576]:


corr_all[0][0][0]


# In[577]:


np.array(corr_all).shape #type, [minusk,plusk], distance


# In[578]:


for ik, k in enumerate(corr_ks):
    for i,cd in enumerate(dvals['cumdist']):
        if k>np.nanmax(cd):
            continue
        ipk = np.argmin(abs(cd-k)) #ipk:
        imk = np.argmin(abs(cd-(cd[-1]-k))) #0:imk
        if imk<2:
            continue
        iip = Sp.find_closest(cd,cd[0:imk]+k)
        N = len(cd)
        for val in corr_vals:
            #all
            corr_all[0][0][ik][val] = np.append(corr_all[0][0][ik][val],dvals[val][i][0:imk])
            corr_all[0][1][ik][val] = np.append(corr_all[0][1][ik][val],dvals[val][i][iip])
            #type 1
            if (dvals['doys'][i][0]> t1[0]) & (dvals['doys'][i][0]< t1[1]):
                corr_all[1][0][ik][val] = np.append(corr_all[1][0][ik][val],dvals[val][i][0:imk])
                corr_all[1][1][ik][val] = np.append(corr_all[1][1][ik][val],dvals[val][i][iip])
            #type 2
            if (dvals['doys'][i][0]> t2[0]) & (dvals['doys'][i][0]< t2[1]):
                corr_all[2][0][ik][val] = np.append(corr_all[2][0][ik][val],dvals[val][i][0:imk])
                corr_all[2][1][ik][val] = np.append(corr_all[2][1][ik][val],dvals[val][i][iip])
            #type 3
            if (dvals['doys'][i][0]> t3[0]) & (dvals['doys'][i][0]< t3[1]):
                corr_all[3][0][ik][val] = np.append(corr_all[3][0][ik][val],dvals[val][i][0:imk])
                corr_all[3][1][ik][val] = np.append(corr_all[3][1][ik][val],dvals[val][i][iip])
            #type 4
            if (dvals['doys'][i][0]> t4[0]) & (dvals['doys'][i][0]< t4[1]):
                corr_all[4][0][ik][val] = np.append(corr_all[4][0][ik][val],dvals[val][i][0:imk])
                corr_all[4][1][ik][val] = np.append(corr_all[4][1][ik][val],dvals[val][i][iip])
            #type 5
            if (dvals['alt'][i][0]> z1[0]) & (dvals['alt'][i][0]< z1[1]):
                corr_all[5][0][ik][val] = np.append(corr_all[5][0][ik][val],dvals[val][i][0:imk])
                corr_all[5][1][ik][val] = np.append(corr_all[5][1][ik][val],dvals[val][i][iip])
            #type 6
            if (dvals['alt'][i][0]> z2[0]) & (dvals['alt'][i][0]< z2[1]):
                corr_all[6][0][ik][val] = np.append(corr_all[6][0][ik][val],dvals[val][i][0:imk])
                corr_all[6][1][ik][val] = np.append(corr_all[6][1][ik][val],dvals[val][i][iip])
            #type 7
            if (dvals['alt'][i][0]> z3[0]) & (dvals['alt'][i][0]< z3[1]):
                corr_all[7][0][ik][val] = np.append(corr_all[7][0][ik][val],dvals[val][i][0:imk])
                corr_all[7][1][ik][val] = np.append(corr_all[7][1][ik][val],dvals[val][i][iip])


# In[591]:


autocorr = {}
for val in corr_vals:
    autocorr[val] = np.zeros((len(types),len(corr_ks)))+np.nan
    
    for ik,k in enumerate(corr_ks):
        for j,jt in enumerate(types):
            if val is 'AE':
                igd = np.where(corr_all[2][1][3]['AE']<5.0)[0] #np.isfinite(corr_all[j][0][ik][val])
                mmk = np.nanmean(corr_all[j][0][ik][val][igd])
                mpk = np.nanmean(corr_all[j][1][ik][val][igd])
                smk = np.nanstd(corr_all[j][0][ik][val][igd])
                spk = np.nanstd(corr_all[j][1][ik][val][igd])
                top = [(v-mpk)*(corr_all[j][1][ik][val][igd][iv]-mmk) for iv,v in enumerate(corr_all[j][0][ik][val][igd])]
                print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
                autocorr[val][j,ik] = np.nansum(top)/((len(corr_all[j][0][ik][val][igd])-1)*spk*smk)
            else:
                mmk = np.nanmean(corr_all[j][0][ik][val])
                mpk = np.nanmean(corr_all[j][1][ik][val])
                smk = np.nanstd(corr_all[j][0][ik][val])
                spk = np.nanstd(corr_all[j][1][ik][val])
                top = [(v-mpk)*(corr_all[j][1][ik][val][iv]-mmk) for iv,v in enumerate(corr_all[j][0][ik][val])]
                print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
                autocorr[val][j,ik] = np.nansum(top)/((len(corr_all[j][0][ik][val])-1)*spk*smk)


# In[582]:


max(corr_all[2][1][3]['AE'])


# In[584]:


len(corr_all[2][1][3]['AE'])


# In[586]:


np.where(corr_all[2][1][3]['AE']>5.0)


# In[590]:


corr_all[2][1][3]['AE'][7091]


# In[585]:


plt.figure()
plt.plot(corr_all[2][1][3]['AE'])


# In[523]:


plt.figure()
plt.plot(corr_ks,autocorr['AE'][0,:],'s-')


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


# In[198]:


len(dvals['aod0500'][1])


# In[199]:


dvals['dist'][1]


# In[200]:


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


# In[202]:


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


# ### Autocorrelation plots

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


# ## Now get the angstrom exponent and plot it vertically

# In[210]:


nwl,nm


# In[89]:


aodrr = np.array([ar[n] for n in nwl])


# In[90]:


aodrr.shape


# In[91]:


angs = su.calc_angs(ar['Start_UTC'],np.array(nm[1:11]),aodrr[1:11,:])


# In[44]:


def make_bined_alt(x,alt,days,fl,n=70):
    'Function to create binned data for a set range, usually for altitude'
    binned_ang,binned_alt,binned_num,binned_ndays = [],[],[],[]
    for i in xrange(70):
        flaa = (alt[fl]>=i*100.0) & (alt[fl]<(i+1.0)*100.0)
        binned_ang.append(x[fl][flaa])
        binned_alt.append(np.mean([i*100.0,(i+1.0)*100.0]))
        binned_num.append(len(x[fl][flaa]))
        binned_ndays.append(len(np.unique(days[fl][flaa])))
    return binned_ang,binned_alt,binned_num,binned_ndays


# In[45]:


ar['fl_QA_angs'] = ar['fl'] & (ar['AOD0501']>0.05) 


# In[46]:


ar['fl_QA_angs_seoul'] = ar['fl'] & (ar['AOD0501']>0.05) & (ar['Latitude']<37.75) &                        (ar['Latitude']>36.9) & (ar['Longitude']<127.30) & (ar['Longitude']>126.60)


# In[47]:


any(ar['fl_QA_angs_seoul'])


# In[70]:


bang,balt,bnum,bndays = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs'],n=90)


# In[71]:


bangs,balts,bnums,bndayss = make_bined_alt(angs,ar['GPS_Alt'],ar['days'],ar['fl_QA_angs_seoul'],n=90)


# ### Plotting of the angstrom vertical dependence

# In[128]:


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


# In[ ]:




