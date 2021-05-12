#!/usr/bin/python/

import sys
import numpy as np

vv = 'v1'
fp = '/data/sam/KORUS-AQ/'
autocorr_dat = np.load(fp+'autocorr_dat_{}.mat.npy'.format(vv),allow_pickle=True)

corr_vals = autocorr_dat['corr_vals']
types = autocorr_dat['types']
corr_ks = autocorr_dat['corr_ks']
corr_all = autocorr_dat['corr_all']
corr_alln = autocorr_dat['corr_alln']
min_num = autocorr_dat['min_num']

autocorr = {}
autocorr_len = {}

def run_autocorr(val,types=types,corr_ks=corr_ks,corr_all=corr_all,corr_alln=corr_alln,min_num=min_num):
    'to run the autocorrelate, and to help with garbage collection'
    print val
    autocorr = np.zeros((len(types),len(corr_ks)))+np.nan
    autocorr_len = np.zeros((len(types),len(corr_ks)))
    for ik,k in enumerate(corr_ks):
        sys.stdout.write("{} ".format(ik))
        for j,jt in enumerate(types):
            corr_all[j][0][ik][val] = np.hstack(corr_alln[j][0][ik][val])
            corr_all[j][1][ik][val] = np.hstack(corr_alln[j][1][ik][val])

            mmk = np.nanmean(corr_all[j][0][ik][val])
            mpk = np.nanmean(corr_all[j][1][ik][val])
            smk = np.nanstd(corr_all[j][0][ik][val])
            spk = np.nanstd(corr_all[j][1][ik][val])
            top = [(v-mpk)*(corr_all[j][1][ik][val][iv]-mmk) for iv,v in enumerate(corr_all[j][0][ik][val])]
            #print val,ik,j,mmk,mpk,smk,spk,np.nansum(top)
            autocorr[j,ik] = np.nansum(top)/((len(corr_all[j][0][ik][val])-1)*spk*smk)
            autocorr_len[j,ik] = len(corr_all[j][0][ik][val])
            if autocorr_len[j,ik]<min_num:
                autocorr[j,ik] = np.nan
    return autocorr, autocorr_len





autocorr = {}
autocorr_len = {}
for val in corr_vals:
    autocorr[val],autocorr_len[val] = run_autocorr(val)

np.save(fp+'autocorr_{}.npy'.format(vv),autocorr)
np.save(fp+'autocorr_len_{}.npy'.format(vv),autocorr_len)

