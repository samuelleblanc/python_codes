#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     collection of utilities for handling some common numerical analysis
# 
# Input:
# 
#     see particular functions
# 
# Output:
# 
#     see functions
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2021-02-02
#     Modified:
# 

# In[ ]:


def __init__():
    """
       Collection of codes to do common tasks in numerical analysis of various data
           
        details are in the info of each module
    """
    import numpy as np
    import scipy.stats as st
    from scipy import interpolate
    import map_utils as mu
    pass


# In[ ]:


def running_std(x,n):
    'Function to do a running standard deviation on array (x) with window size (n)'
    q = x**2
    q = np.convolve(q, np.ones((n, )), mode="same")
    s = np.convolve(x, np.ones((n, )), mode="same")
    o = (q-s**2/n)/float(n-1)
    return o 


# In[ ]:


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


# In[ ]:


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


# In[ ]:


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


# In[ ]:


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


# In[ ]:


def rand_autocor(plus_dist,minus_dist,min_num=100,subsamp_ratio=0.5):
    'fx_autocor, but called with random subsampling at a ratio defined by subsamp_ratio [0-1]'
    irand = np.random.randint(len(plus_dist),size=int(subsamp_ratio*len(plus_dist)))
    return fx_autocor(plus_dist[irand],minus_dist[irand],min_num=min_num)


# In[ ]:


def autocorr(x, t=1):
    return np.corrcoef(np.array([x[:-t], x[t:]]))


# In[ ]:


def autocorr2(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]/result.max()


# In[ ]:


def autocorr5(x):
    '''numpy.correlate, non partial'''
    n = len(x)
    lags = range(n)
    #mean=x.mean()
    var=np.nanvar(x)
    xp=x-np.nanmean(x)
    corr=np.correlate(xp,xp,'full')[n-1:]/var/n

    return corr[:n]


# In[ ]:


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


# In[ ]:


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


# In[ ]:


def make_binned(x,alt,fl,bn,flb):
    'Function to create binned data for a set range, usually for altitude'
    
    binned_ang,binned_alt,binned_num = [],[],[]
    for i,b in enumerate(bn[:-1]):
        flaa = (alt[flb]>=b) & (alt[flb]<bn[i+1])
        binned_ang.append(x[:,flb][flaa])
        binned_alt.append(np.mean([b,bn[i+1]]))
        binned_num.append(len(x[fl][:,flaa]))
    return binned_ang,binned_alt,binned_num,binned_ndays


# In[ ]:


def stats_2d(lat,lon,x,fl=[],bins=26,rg=[[-25,-8],[0,16]],days=[],verbose=True):
    "Combined Statistics function to get the mean, median, number, ndays, and std from a 2d dataset"    
    stat = {}
    if not len(fl)>0: fl = np.isfinite(x)
        
    stat['mean'],stat['xm'],stat['ym'],stat['bin'] =           st.binned_statistic_2d(lat[fl],lon[fl],x[fl],bins=bins,range=rg,statistic='mean')
    stat['mean'] = np.ma.masked_array(stat['mean'],np.isnan(stat['mean']))
    
    stat['median'],stat['xe'],stat['ye'],stat['bine'] =           st.binned_statistic_2d(lat[fl],lon[fl],x[fl],bins=bins,range=rg,statistic='median')
    stat['median'] = np.ma.masked_array(stat['median'],np.isnan(stat['median']))

    stat['std'],stat['xs'],stat['ys'],stat['bins'] =           st.binned_statistic_2d(lat[fl],lon[fl],x[fl],bins=bins,range=rg,statistic=np.nanstd)
    stat['std'] = np.ma.masked_array(stat['std'],np.isnan(stat['std']))
    
    stat['cnt'],stat['xn'],stat['yn'],stat['binn'] =           st.binned_statistic_2d(lat[fl],lon[fl],x[fl],bins=bins,range=rg,statistic='count')
    stat['cnt'] = np.ma.masked_array(stat['cnt'],np.isnan(stat['cnt']))

    if len(days)>0:
        uniq_cnt = lambda x: len(np.unique(x))
        stat['dcnt'],stat['xd'],stat['yd'],stat['bind'] =           st.binned_statistic_2d(lat[fl],lon[fl],days[fl],bins=bins,range=rg,statistic=uniq_cnt)
        stat['dcnt'] = np.ma.masked_array(stat['dcnt'],np.isnan(stat['dcnt']))
    else:
        stat['dcnt'] = stat['cnt']*0.0
    
    if verbose:
        print 'Mean values: mean={}, median={}, std={}, num={}, day={}'.format(                    np.nanmean(stat['mean']),np.nanmean(stat['median']),np.nanmean(stat['std']),
                    np.nanmean(stat['cnt']),np.nanmean(stat['dcnt']))
        print 'Median values: mean={}, median={}, std={}, num={}, day={}'.format(                    np.nanmedian(stat['mean']),np.nanmedian(stat['median']),np.nanmedian(stat['std']),
                    np.nanmedian(stat['cnt']),np.nanmedian(stat['dcnt']))
        print 'STD values: mean={}, median={}, std={}, num={}, day={}'.format(                    np.nanstd(stat['mean']),np.nanstd(stat['median']),np.nanstd(stat['std']),
                    np.nanstd(stat['cnt']),np.nanstd(stat['dcnt']))
    return stat


# In[ ]:




