#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     Check_4STAR_cloud_retr_ORACLES
# 
# Purpose:  
# 
#     To go through the steps to check the zenith cloud retrievals from 4STAR.
#     To address issues with the parameter calculations that don't match modeled values at all. 
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
#     - load_utils.py
#     - matplotlib
#     - numpy
#     - Sp
#   
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2020-01-22
#     Modified: 
# 

# # Prepare python environment

# In[36]:


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


# In[2]:


from datetime import datetime
from matplotlib import animation


# In[3]:


get_ipython().run_cell_magic(u'javascript', u'', u'IPython.notebook.kernel.execute(\'nb_name = "\' + IPython.notebook.notebook_name + \'"\')')


# In[4]:


get_ipython().magic(u'matplotlib notebook')


# In[5]:


name = 'ORACLES'
vv = 'v1'
fp = getpath(name)


# The call for the retrieval results:
# 
# python ..\..\python_codes\Run_zen_cld_retrieval.py -lut ..\model\v5_10aod_ORACLES_lut.mat .\4STAR_20170818starzen.mat -o ..\zen_ict\v5_10aod\ -liq -plotlut -movie -taurange 1 60 -refrange 2 30
# 
# 
# For 20160831 with SZA at 31.0
# ![image.png](attachment:image.png)
#     
# ![20160831_norm_zenrad.png](attachment:20160831_norm_zenrad.png)
# 
# ![20160831_zenrad.png](attachment:20160831_zenrad.png)
# 

# # Load files

# In[38]:


daystr = '20160831'


# ## Load the starzen

# In[39]:


fp_starzen = fp+'/starzen/4STAR_{}starzen.mat'.format(daystr)
mea = sio.loadmat(fp_starzen)


# In[40]:


tt = lu.mat2py_time(mea['t'])
mea['utc'] = lu.toutc(tt)


# In[43]:


print 'Running the parameter calculations on measured spectra'
meas = Sp.Sp(mea,verbose=False)
meas.params(norm_1050=True)


# ### Starzen for different date

# In[44]:


daystr2 = '20160914'


# In[45]:


fp_starzen2 = fp+'/starzen/4STAR_{}starzen.mat'.format(daystr2)
mea2 = sio.loadmat(fp_starzen2)


# In[46]:


tt2 = lu.mat2py_time(mea2['t'])
mea2['utc'] = lu.toutc(tt2)


# In[47]:


print 'Running the parameter calculations on measured spectra'
meas2 = Sp.Sp(mea2,verbose=False)
meas2.params(norm_1050=True)


# ## Load the LUT

# In[48]:


vlut = 'v5_10aod'


# In[49]:


fp_lut_mat = fp+'model/{}_ORACLES_lut.mat'.format(vlut)
luts = hs.loadmat(fp_lut_mat)


# In[50]:


airmass = 1./np.cos(luts['sza']*np.pi/180.0)


# # Start the retrieval pieces

# In[51]:


start_tau = 1.0
end_tau = 60.0
start_ref = 2.0
end_ref = 30.0
forceliq = True
forceice = False
iz = 0


# ## Prep the LUT

# In[52]:


meas.airmass = 1.0/np.cos(meas.sza*np.pi/180.0)
idx = Sp.find_closest(airmass,meas.airmass)


# In[53]:


meas2.airmass = 1.0/np.cos(meas2.sza*np.pi/180.0)
idx2 = Sp.find_closest(airmass,meas2.airmass)


# In[54]:


lut = []


# In[55]:


for s in xrange(len(luts['sza'])):
    sptemp = {}
    sptemp['tau'] = luts['tau']
    sptemp['ref'] = luts['ref']
    sptemp['zout'] = luts['zout']
    sptemp['sza'] = luts['sza']
    sptemp['phase'] = luts['phase']
    sptemp['irr_dn_diff'] = luts['irr_dn_diff'][:,:,:,:,:,s]
    sptemp['irr_dn'] = luts['irr_dn'][:,:,:,:,:,s]
    sptemp['irr_up'] = luts['irr_up'][:,:,:,:,:,s]
    sptemp['wvl'] = [luts['wvl']]
    sptemp['rad'] = luts['rad'][:,:,:,:,:,s]
    ltemp = Sp.Sp(sptemp,verbose=False)
    if s in idx:
        ltemp.params(liq_only=forceliq,ice_only=forceice,iz=iz,norm_1050=True)
        ltemp.param_hires(start_ref=start_ref,end_ref=end_ref,start_tau=start_tau,end_tau=end_tau)
    lut.append(ltemp)


# ## Now run through the retrieval

# In[56]:


(meas.taut,meas.ref,meas.phase,meas.ki) = (np.zeros_like(meas.utc),np.zeros_like(meas.utc),np.zeros_like(meas.utc),np.zeros_like(meas.utc))


# In[57]:


(meas.taut,meas.ref,meas.phase,meas.ki) = (meas.taut*np.nan,meas.ref*np.nan,meas.phase*np.nan,meas.ki*np.nan)


# In[58]:


(meas2.taut,meas2.ref,meas2.phase,meas2.ki) = (np.zeros_like(meas2.utc),np.zeros_like(meas2.utc),np.zeros_like(meas2.utc),np.zeros_like(meas2.utc))
(meas2.taut,meas2.ref,meas2.phase,meas2.ki) = (meas2.taut*np.nan,meas2.ref*np.nan,meas2.phase*np.nan,meas2.ki*np.nan)


# In[59]:


import run_kisq_retrieval as rk


# In[60]:


print 'Running through the airmasses'
for i in np.unique(idx):
    try: 
        print 'airmass: {airmass}, {i}/{i_tot}'.format(airmass=airmass[i],i=i,i_tot=idx.max()-idx.min())
    except:
        if debug:
            import pdb; pdb.set_trace()
        else: 
            print 'Exception occurred at airmass: {} *** skipping'.format(airmass[i])
    meas.good = np.where(idx==i)[0]

    i500 = np.argmin(abs(meas.wvl-500))
    #i980,i995 = np.argmin(abs(meas.wvl-980)),np.argmin(abs(meas.wvl-995))
    #i999,i992 = np.argmin(abs(meas.wvl-999)),np.argmin(abs(meas.wvl-992))
    i981,i982 = 1039,1040
    #ss = np.nanstd(meas.norm[meas.good,i980:i995],axis=1)/np.nanmean(meas.norm[meas.good,i980:i995],axis=1)
    sss = abs(meas.norm[meas.good,i981]-meas.norm[meas.good,i982])
    #flt = (meas.norm[meas.good,i500]>0.4) & (ss<0.05)
    flt = sss<0.1 #(meas.norm[meas.good,i500]>0.4) & (sss<0.1)
    #import pdb; pdb.set_trace()
    meas.good = meas.good[flt]

    print 'meas.good lengh: {},meas.utc length: {}'.format(meas.good.shape,meas.utc.shape)
    tau,ref,phase,ki = rk.run_retrieval(meas,lut[i],force_liq=forceliq,force_ice=forceice)
    meas.taut[meas.good] = tau[meas.good]
    meas.ref[meas.good] = ref[meas.good]
    meas.phase[meas.good] = phase[meas.good]
    meas.ki[meas.good] = ki[meas.good]


# In[61]:


meas.tau = meas.taut


# In[62]:


print 'Running through the airmasses'
for i in np.unique(idx2):
    try: 
        print 'airmass: {airmass}, {i}/{i_tot}'.format(airmass=airmass[i],i=i,i_tot=idx.max()-idx.min())
    except:
        if debug:
            import pdb; pdb.set_trace()
        else: 
            print 'Exception occurred at airmass: {} *** skipping'.format(airmass[i])
    meas2.good = np.where(idx2==i)[0]

    i500 = np.argmin(abs(meas2.wvl-500))
    #i980,i995 = np.argmin(abs(meas.wvl-980)),np.argmin(abs(meas.wvl-995))
    #i999,i992 = np.argmin(abs(meas.wvl-999)),np.argmin(abs(meas.wvl-992))
    i981,i982 = 1039,1040
    #ss = np.nanstd(meas.norm[meas.good,i980:i995],axis=1)/np.nanmean(meas.norm[meas.good,i980:i995],axis=1)
    sss = abs(meas2.norm[meas2.good,i981]-meas2.norm[meas2.good,i982])
    #flt = (meas.norm[meas.good,i500]>0.4) & (ss<0.05)
    flt = sss<0.1 #(meas.norm[meas.good,i500]>0.4) & (sss<0.1)
    #import pdb; pdb.set_trace()
    meas2.good = meas2.good[flt]

    print 'meas2.good lengh: {},meas2.utc length: {}'.format(meas2.good.shape,meas2.utc.shape)
    tau,ref,phase,ki = rk.run_retrieval(meas2,lut[i],force_liq=forceliq,force_ice=forceice)
    meas2.taut[meas2.good] = tau[meas2.good]
    meas2.ref[meas2.good] = ref[meas2.good]
    meas2.phase[meas2.good] = phase[meas2.good]
    meas2.ki[meas2.good] = ki[meas2.good]


# In[63]:


meas2.tau = meas2.taut


# # Now plot out the intermediates and results

# ## Look at overall values

# In[64]:


meas.keys()


# In[65]:


plt.figure()
plt.plot(meas['utc'],meas['par'][:,0],'.')


# In[66]:


lut[0].keys()


# In[67]:


lut[0]['par'].shape


# In[68]:


lut[0]['tau'].shape


# In[69]:


lut[0]['ref'].shape


# In[70]:


lut[0]['sza'].shape


# In[71]:


plt.figure()
plt.plot(lut[0]['tau'],lut[0]['par'][0,10,:,0],'.')


# ## Plot out the lut

# In[72]:


from mpltools import color


# In[73]:


sza = np.arccos(1.0/airmass)*180.0/np.pi


# In[74]:


im = 1


# In[75]:


for im in range(4):
    fig3,ax3 = plt.subplots(5,3,sharex=True,figsize=(15,8))
    ax3 = ax3.ravel()

    for i in range(lut[im].npar-1):
        color.cycle_cmap(len(lut[im].ref[lut[im].ref<30]),cmap=plt.cm.RdBu,ax=ax3[i])
        for j in xrange(len(lut[im].ref)):
            ax3[i].plot(lut[im].tau,lut[im].par[0,j,:,i])
        ax3[i].set_title('Parameter '+str(i+1))
        ax3[i].grid()
        ax3[i].set_xlim([0,60])
        if i > 11: 
            ax3[i].set_xlabel('Tau')

    fig3.tight_layout()
    plt.suptitle('Liquid SZA={:2.1f}$^{{\circ}}$'.format(sza[im]))
    plt.subplots_adjust(top=0.93,right=0.93)

    cbar_ax = fig3.add_axes([0.95,0.10,0.02,0.8])
    scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.RdBu,norm=plt.Normalize(vmin=0,vmax=1))
    scalarmap.set_array(lut[im].ref[lut[im].ref<30])
    cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
    cba.ax.set_ylabel('R$_{ef}$ [$\\mu$m]')
    cba.ax.set_yticklabels(np.linspace(lut[im].ref[0],29,6));

#plt.show()


# In[103]:


for im in range(4):
    fig3,ax3 = plt.subplots(5,3,sharex=True,figsize=(15,8))
    ax3 = ax3.ravel()

    for i in range(lut[im].npar-1):
        color.cycle_cmap(len(lut[im].ref[lut[im].ref<30]),cmap=plt.cm.RdBu,ax=ax3[i])
        for j in xrange(len(lut[im].ref)):
            ax3[i].plot(lut[im].tau,lut[im].par[0,j,:,i])
        ax3[i].set_title('Parameter '+str(i+1))
        ax3[i].grid()
        ax3[i].set_xlim([0,60])
        if i > 11: 
            ax3[i].set_xlabel('Tau')
        ig = np.where(np.isfinite(meas.par[:,0]) & (idx[:,0]==im))[0]
        ax3[i].plot(meas.tau[ig],meas.par[ig,i],'xk')

    fig3.tight_layout()
    plt.suptitle('Liquid SZA={:2.1f}$^{{\circ}}$'.format(sza[im]))
    plt.subplots_adjust(top=0.93,right=0.93)

    cbar_ax = fig3.add_axes([0.95,0.10,0.02,0.8])
    scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.RdBu,norm=plt.Normalize(vmin=0,vmax=1))
    scalarmap.set_array(lut[im].ref[lut[im].ref<30])
    cba = plt.colorbar(scalarmap,ticks=np.linspace(0,1,6),cax=cbar_ax)
    cba.ax.set_ylabel('R$_{ef}$ [$\\mu$m]')
    cba.ax.set_yticklabels(np.linspace(lut[im].ref[0],29,6));
    
#plt.show()


# In[76]:


meas.par.shape


# In[80]:



fig3,ax3 = plt.subplots(5,3,sharex=True,figsize=(15,8))
ax3 = ax3.ravel()

for i in range(lut[im].npar-1):
    #color.cycle_cmap(len(lut[im].ref[lut[im].ref<30]),cmap=plt.cm.RdBu,ax=ax3[i])
    #for j in xrange(len(lut[im].ref)):
    #    ax3[i].plot(lut[im].tau,lut[im].par[0,j,:,i])
    
    ax3[i].plot(meas.utc[igood],meas.par[igood,i],'x-k')
    ax3[i].axhline(lut[1].par[0,0,0,i],label='COD:{},REF:{}'.format(lut[1].tau[0],lut[1].ref[0]),color='r')
    ax3[i].axhline(lut[1].par[0,-1,0,i],label='COD:{},REF:{}'.format(lut[1].tau[0],lut[1].ref[-1]),color='purple')
    ax3[i].axhline(lut[1].par[0,0,-1,i],label='COD:{},REF:{}'.format(lut[1].tau[-1],lut[1].ref[0]),color='pink')
    ax3[i].axhline(lut[1].par[0,-1,-1,i],label='COD:{},REF:{}'.format(lut[1].tau[-1],lut[1].ref[-1]),color='b')
    
    if i<1: 
        ax3[i].legend(frameon=False)
    ax3[i].set_title('Parameter '+str(i+1))
    ax3[i].grid()
    #ax3[i].set_xlim([0,60])
    if i > 11: 
        ax3[i].set_xlabel('UTC [h]')

fig3.tight_layout()
plt.suptitle('Measurements'.format(sza[im]))
plt.subplots_adjust(top=0.93,right=0.93)




# In[78]:


im = 1
pl = Sp.plt_lut_zenrad(lut[im],range_variable='tausp',norm=False,lims=[0,60],cmapname='gist_earth')


# In[80]:


im = 1
pl = Sp.plt_lut_zenrad(lut[im],range_variable='tausp',norm=True,lims=[0,60],cmapname='gnuplot')


# In[85]:


im = 1
pl = Sp.plt_lut_zenrad(lut[im],range_variable='refsp',norm=False,other_index=14,lims=[0,30],cmapname='brg')


# ## Plot out the measurements

# In[86]:


p = Sp.plt_norm_zenrad(meas)


# In[87]:


p = Sp.plt_zenrad(meas)


# In[90]:


p = Sp.plt_norm_zenrad(meas)
pl = Sp.plt_lut_zenrad(lut[im],range_variable='tausp',norm=True,lims=[0,60],cmapname='hsv',fig=p)


# In[133]:


p = Sp.plt_zenrad(meas)
pl = Sp.plt_lut_zenrad(lut[im],range_variable='tausp',norm=False,lims=[0,60],cmapname='spring',fig=p)
pl.set_zorder(-1)
plt.show()


# ### Now for meas2

# In[134]:


p = Sp.plt_zenrad(meas2)
pl = Sp.plt_lut_zenrad(lut[im],range_variable='tausp',norm=False,lims=[0,60],cmapname='spring',fig=p)


# # Got through each parameter

# ## All parameter definition

# In[41]:


from Sp_parameters import nanmasked, norm2max, smooth, deriv, find_closest


# In[31]:


def param(sp,wvlin,iws=None):
    " Calculates the parameters from a spectrum."
    from linfit import linfit
    from Sp_parameters import nanmasked, norm2max, smooth, deriv, find_closest
    npar = 16
    #print 'lenght of sp and wvlin',sp.shape,wvlin.shape
    spc, mask = nanmasked(sp)
    try:
        if len(mask.shape)==3:
            mask = mask[:,0,0]
    except:
        pass
    wvl = wvlin[mask]
    if len(wvl)<600:
        par = np.zeros(npar)*np.nan
        return par
    try:
        norm = norm2max(spc,iws=iws)
    except ValueError:
        par = np.zeros(npar)*np.nan
        return par
    [i1000,i1077,i1493,i1600,i1200,i1300,i530,i610,
     i1565,i1634,i1193,i1198,i1236,i1248,i1270,i1644,
     i1050,i1040,i1065,i600,i870,i515] = find_closest(wvl,np.array([1000,1077,1493,1600,1200,1300,530,
                                                     610,1565,1634,1193,1198,1236,1248,
                                                     1270,1644,1050,1040,1065,600,870,515]))
    if np.isnan(spc[i1000]) or not spc[i1000]:
        par = np.zeros(npar)*np.nan
        return par
    norm2 = spc/spc[i1000]
    try:
        dsp = smooth(deriv(norm2,wvl/1000),2,nan=False,old=True)
    except:
        import pdb; pdb.set_trace()
    imaxwvl = np.argmax(spc)
    maxwvl = wvl[mask[imaxwvl]]
    # now calculate each parameter
    fit0 = np.polyfit(np.array([wvl[i1000],wvl[i1077]]),np.array([norm2[i1000],norm2[i1077]]),1)
    fit0_fn = np.poly1d(fit0)
    fit7 = np.polyfit(np.array([wvl[i1493],wvl[i1600]]),np.array([norm2[i1493],norm2[i1600]]),1)
    fit7_fn = np.poly1d(fit7)
    fit8,z = linfit(wvl[i1000:i1077],dsp[i1000:i1077]) 
    fit9,z = linfit(wvl[i1200:i1300],dsp[i1200:i1300])
    fit10,z = linfit(wvl[i530:i610]/1000,norm[i530:i610])
    fit14,z = linfit(wvl[i1565:i1634],spc[i1565:i1634]/norm[i1565])
    par = [sum(norm2[i1000:i1077]-fit0_fn(wvl[i1000:i1077])),   # 1 curvature of rad normed to 1000 nm for 1000 nm - 1077 nm
           dsp[i1198],                                          # 2 deriv of rad normed to 1000 nm at 1198 nm (!=IDL version)
           dsp[i1493],                                          # 3 deriv of rad normed to 1000 nm at 1493 nm
           norm[i1198]/norm[i1236],                             # 4 ratio of normalized rad of 1198 nm / 1236 nm
           np.nanmean(norm[i1248:i1270]),                       # 5 mean of normalized rad between 1248 nm - 1270 nm
           np.nanmean(norm[i1565:i1644]),                       # 6 mean of normalized rad between 1565 nm - 1644 nm
           np.nanmean(norm[i1000:i1050]),                       # 7 mean of normalized rad between 1000 nm - 1050 nm
           sum(norm2[i1493:i1600]-fit7_fn(wvl[i1493:i1600])),   # 8 curvature of rad normed to 1000 nm for 1493 nm - 1600 nm
           fit8[0],                                             # 9 slope of deriv of rad normed to 1000 nm, 1000 nm - 1077 nm
           fit9[0],                                             # 10 slope of deriv of rad normed to 1000 nm, 1200 nm - 1300 nm
           fit10[0],                                            # 11 slope of normalized radiance between 530 nm - 610 nm
           norm[i1040],                                         # 12 normalized radiance at 1040 nm
           norm[i1000]/norm[i1065],                             # 13 ratio of normalized radiance at 1000 nm / 1065 nm
           norm[i600]/norm[i870],                               # 14 ratio of normalized radiance at 600 nm / 870 nm
           np.nanmin([0.003,fit14[0]]),                         # 15 slope of radiance / rad at 1565 between 1565 nm - 1634 nm
           spc[i515]]                                           # 16 radiance at 515 nm
    # do a check for bad points
    if np.all(np.isnan(par[0:13])): 
        par[14] = np.nan
        par[15] = np.nan
    return param


# In[78]:


igood = np.where((idx[:,0]==1) & np.isfinite(meas.par[:,0]))[0]


# In[79]:


igood2 = np.where((idx2[:,0]==1) & np.isfinite(meas2.par[:,0]))[0]


# ## Parameter 1: Curvature between 1000-1100

# In[33]:


meas.wvl


# In[34]:


i1000 = np.argmin(abs(meas.wvl-1000.0))
i1077 = np.argmin(abs(meas.wvl-1077.0))


# In[35]:


def plt_param1(sp,wvl,fig=None):
    from Sp_parameters import nanmasked
    if not fig:
        fig = plt.figure()
    
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1077 = np.argmin(abs(wvln-1077.0))
    norm = spc / spc[i1000]
    
    line_og = plt.plot(wvln,norm,'.b',label='original')
    line_nd = plt.plot([wvln[i1000],wvln[i1077]],[norm[i1000],norm[i1077]],'sg',label='ends')
    fit0 = np.polyfit(np.array([wvln[i1000],wvln[i1077]]),np.array([norm[i1000],norm[i1077]]),1)
    fit0_fn = np.poly1d(fit0)
    line_ft = plt.plot(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]),'-k',lw=2,label='line fit')
    line_cv = plt.plot(wvln[i1000:i1077],norm[i1000:i1077],'-r',lw=2,label='curve')
    par1 = sum(norm[i1000:i1077]-fit0_fn(wvln[i1000:i1077]))
    fill_pr = plt.fill_between(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]),norm[i1000:i1077],color='grey',alpha=0.4,label='par1= {:2.2f}'.format(par1))
    
    plt.xlim([970,1130])
    plt.ylim([0,1.5])
    plt.legend(frameon=False)
    return fig, line_og[0],line_nd[0],line_ft[0],line_cv[0],fill_pr


# In[36]:


plt_param1(meas.sp[igood][0],meas.wvl)


# ### Make some animations

# In[37]:


from IPython.display import HTML

def display_animation(anim):
    plt.close(anim._fig)
    return HTML(anim_to_html(anim))


# In[277]:


p.vertices[:,1]


# In[282]:


np.hstack([fit0_fn(wvln[i1000:i1077])[0],fit0_fn(wvln[i1000:i1077]),norm[i1000:i1077][-1],norm[i1000:i1077][:0:-1],norm[i1000:i1077][0],norm[i1000:i1077][0]])


# In[255]:


fig, line_og,line_nd,line_ft,line_cv,fill_pr = plt_param1(meas.sp[igood][0],meas.wvl)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Zenith radiance normalized from value at 1000 nm')

def init():
    line_og.set_data([], [])
    line_nd.set_data([], [])
    line_ft.set_data([], [])
    line_cv.set_data([], [])
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = p.vertices[:,1]*0.0
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

def animate(i):
    
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1077 = np.argmin(abs(wvln-1077.0))
    norm = spc / spc[i1000]
    
    line_og.set_data(wvln,norm)
    line_nd.set_data([wvln[i1000],wvln[i1077]],[norm[i1000],norm[i1077]])
    fit0 = np.polyfit(np.array([wvln[i1000],wvln[i1077]]),np.array([norm[i1000],norm[i1077]]),1)
    fit0_fn = np.poly1d(fit0)
    line_ft.set_data(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]))
    line_cv.set_data(wvln[i1000:i1077],norm[i1000:i1077])
    par1 = sum(norm[i1000:i1077]-fit0_fn(wvln[i1000:i1077]))    
    #fill_pr = fig.gca().fill_between(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]),norm[i1000:i1077],color='grey',alpha=0.4,label='par1= {:2.2f}'.format(par1))
    
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = np.hstack([fit0_fn(wvln[i1000:i1077])[0],fit0_fn(wvln[i1000:i1077]),
                                 norm[i1000:i1077][-1],norm[i1000:i1077][:0:-1],norm[i1000:i1077][0],norm[i1000:i1077][0]])
    fill_pr.set_label('par1= {:2.2f}'.format(par1))
    plt.legend(frameon=False)
    
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=20, blit=True)

#HTML(anim.to_jshtml())


# In[138]:


nb_name


# In[139]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 1 animations for {}'.format(daystr),comment='Created on {} with the program {}'.format(datetime.now(),nb_name))


# In[140]:


# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)


# In[141]:


anim.save(fp+'starzen/Parameter1_{}.mp4'.format(daystr),writer=writer)


# ### Plot out measurement day 2

# In[151]:


fig, line_og,line_nd,line_ft,line_cv,fill_pr = plt_param1(meas2.sp[igood2][0],meas2.wvl)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Zenith radiance normalized from value at 1000 nm')

def init():
    line_og.set_data([], [])
    line_nd.set_data([], [])
    line_ft.set_data([], [])
    line_cv.set_data([], [])
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = p.vertices[:,1]*0.0
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

def animate(i):
    
    sp,wvl = meas2.sp[igood2][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr2,meas2.utc[igood2][i],meas2.sza[igood2][i][0],i))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1077 = np.argmin(abs(wvln-1077.0))
    norm = spc / spc[i1000]
    
    line_og.set_data(wvln,norm)
    line_nd.set_data([wvln[i1000],wvln[i1077]],[norm[i1000],norm[i1077]])
    fit0 = np.polyfit(np.array([wvln[i1000],wvln[i1077]]),np.array([norm[i1000],norm[i1077]]),1)
    fit0_fn = np.poly1d(fit0)
    line_ft.set_data(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]))
    line_cv.set_data(wvln[i1000:i1077],norm[i1000:i1077])
    par1 = sum(norm[i1000:i1077]-fit0_fn(wvln[i1000:i1077]))    
    #fill_pr = fig.gca().fill_between(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]),norm[i1000:i1077],color='grey',alpha=0.4,label='par1= {:2.2f}'.format(par1))
    
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = np.hstack([fit0_fn(wvln[i1000:i1077])[0],fit0_fn(wvln[i1000:i1077]),
                                 norm[i1000:i1077][-1],norm[i1000:i1077][:0:-1],norm[i1000:i1077][0],norm[i1000:i1077][0]])
    fill_pr.set_label('par1= {:2.2f}'.format(par1))
    plt.legend(frameon=False)
    
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

anim2 = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas2.utc[igood2]), interval=10, blit=True)


# In[152]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 1 animations for {}'.format(daystr2),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim2.save(fp+'starzen/Parameter1_{}.mp4'.format(daystr2),writer=writer)


# ### Plot out animation for lut

# In[288]:


lut[0].sp.shape


# In[289]:


lut[1].tausp.shape


# In[84]:


lut[1].sza[1]


# In[254]:


fig, line_og,line_nd,line_ft,line_cv,fill_pr = plt_param1(lut[1].sp[0,:,0,10,0],lut[1].wvl)
                                                                    
plt.xlabel('Wavelength [nm]')
plt.ylabel('Zenith radiance normalized from value at 1000 nm')

def init():
    line_og.set_data([], [])
    line_nd.set_data([], [])
    line_ft.set_data([], [])
    line_cv.set_data([], [])
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = p.vertices[:,1]*0.0
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

def animate(i):
    
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1077 = np.argmin(abs(wvln-1077.0))
    norm = spc / spc[i1000]
    
    line_og.set_data(wvln,norm)
    line_nd.set_data([wvln[i1000],wvln[i1077]],[norm[i1000],norm[i1077]])
    fit0 = np.polyfit(np.array([wvln[i1000],wvln[i1077]]),np.array([norm[i1000],norm[i1077]]),1)
    fit0_fn = np.poly1d(fit0)
    line_ft.set_data(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]))
    line_cv.set_data(wvln[i1000:i1077],norm[i1000:i1077])
    par1 = sum(norm[i1000:i1077]-fit0_fn(wvln[i1000:i1077]))    
    #fill_pr = fig.gca().fill_between(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]),norm[i1000:i1077],color='grey',alpha=0.4,label='par1= {:2.2f}'.format(par1))
    
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = np.hstack([fit0_fn(wvln[i1000:i1077])[0],fit0_fn(wvln[i1000:i1077]),
                                 norm[i1000:i1077][-1],norm[i1000:i1077][:0:-1],norm[i1000:i1077][0],norm[i1000:i1077][0]])
    fill_pr.set_label('par1= {:2.2f}'.format(par1))
    plt.legend(frameon=False)
    
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[144]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 1 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),comment='Created on {} with the program {}'.format(datetime.now(),nb_name))


# In[145]:


# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter1_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 2: Derivative at 1200 nm

# ### Define parameter

# In[182]:


def plt_param2(sp,wvl,fig=None):
    from Sp_parameters import nanmasked
    if not fig:
        fig = plt.figure()


    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1198 = np.argmin(abs(wvln-1198.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    par2 = dsp[i1198]
    
    
    ax1 = fig.add_subplot(2,1,1)
    line_og = ax1.plot(wvln,spc,'.b',label='original')
    line_sp = ax1.plot([wvln[i1198],wvln[i1198]],[spc[i1198],spc[i1198]],'ok')
    ax1.set_ylabel('Radiance [W/m$^2$/nm/sr]')
    ax1.set_ylim(0,300)
    
    ax2 = fig.add_subplot(2,1,2,sharex=ax1)
    line_dsp = ax2.plot(wvln,dsp,'.k',label='Differential')
    line_par = ax2.plot([wvln[i1198],wvln[i1198]],[par2,par2],'sg',label='Par2={:2.2f}'.format(par2)) 
    ax2.set_xlim(970,1300)
    ax2.set_ylim(-30,40)
    ax2.set_ylabel('Derivative of\nnormalized radiance')
    ax2.set_xlabel('Wavelength [nm]')
    
    plt.legend(frameon=False)
    return fig, ax1,ax2,line_og[0],line_sp[0],line_dsp[0],line_par[0]


# In[157]:


sp = meas.sp[igood][0]
wvl = meas.wvl


# In[178]:


plt_param2(meas.sp[igood][0],meas.wvl)


# ### Plot the measurements

# In[256]:


fig, ax1,ax2,line_og,line_sp,line_dsp,line_par = plt_param2(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_sp.set_data([], [])
    line_dsp.set_data([], [])
    line_par.set_data([], [])
    return fig,ax1,ax2,line_og,line_sp,line_dsp,line_par

def animate(i):
    sp,wvl = meas.sp[igood][i],meas.wvl
    ax1.set_title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1198 = np.argmin(abs(wvln-1198.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    par2 = dsp[i1198]
    
    line_og.set_data(wvln,spc)
    line_sp.set_data([wvln[i1198],wvln[i1198]],[spc[i1198],spc[i1198]])
    line_dsp.set_data(wvln,dsp)
    line_par.set_data([wvln[i1198],wvln[i1198]],[par2,par2])
    
    line_par.set_label('Par2={:2.2f}'.format(par2))
    ax2.legend(frameon=False)
    
    return fig,ax1,ax2,line_og,line_sp,line_dsp,line_par

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=20, blit=True)

#HTML(anim.to_jshtml())


# In[186]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 2 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))


# In[187]:


# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter2_{}.mp4'.format(daystr),writer=writer)


# ### Plot out the animation for the LUT

# In[257]:


fig, ax1,ax2,line_og,line_sp,line_dsp,line_par = plt_param2(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_sp.set_data([], [])
    line_dsp.set_data([], [])
    line_par.set_data([], [])
    return fig,ax1,ax2,line_og,line_sp,line_dsp,line_par

def animate(i):
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    ax1.set_title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1198 = np.argmin(abs(wvln-1198.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    par2 = dsp[i1198]
    
    line_og.set_data(wvln,spc)
    line_sp.set_data([wvln[i1198],wvln[i1198]],[spc[i1198],spc[i1198]])
    line_dsp.set_data(wvln,dsp)
    line_par.set_data([wvln[i1198],wvln[i1198]],[par2,par2])
    
    line_par.set_label('Par2={:2.2f}'.format(par2))
    ax2.legend(frameon=False)
    
    return fig,ax1,ax2,line_og,line_sp,line_dsp,line_par

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[191]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 2 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))


# In[192]:


# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter2_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 3: Derivative at 1500 nm

# ### Define Parameter 3

# In[205]:


def plt_param3(sp,wvl,fig=None):
    from Sp_parameters import nanmasked
    if not fig:
        fig = plt.figure()

    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1493 = np.argmin(abs(wvln-1493.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    par3 = dsp[i1493]
    
    ax1 = fig.add_subplot(2,1,1)
    line_og = ax1.plot(wvln,spc,'.b',label='original')
    line_sp = ax1.plot([wvln[i1493],wvln[i1493]],[spc[i1493],spc[i1493]],'ok')
    ax1.set_ylabel('Radiance [W/m$^2$/nm/sr]')
    ax1.set_ylim(0,150.0)
    
    ax2 = fig.add_subplot(2,1,2,sharex=ax1)
    line_dsp = ax2.plot(wvln,dsp,'.k',label='Differential')
    line_par = ax2.plot([wvln[i1493],wvln[i1493]],[par3,par3],'sg',label='Par3={:2.2f}'.format(par3)) 
    ax2.set_xlim(1450,1650)
    ax2.set_ylim(-20,20)
    ax2.set_ylabel('Derivative of\nnormalized radiance')
    ax2.set_xlabel('Wavelength [nm]')
    
    plt.legend(frameon=False)
    return fig, ax1,ax2,line_og[0],line_sp[0],line_dsp[0],line_par[0]


# In[206]:


sp = meas.sp[igood][0]
wvl = meas.wvl


# In[208]:


plt_param3(meas.sp[igood][0],meas.wvl)


# ### Plot the measurements

# In[258]:


fig, ax1,ax2,line_og,line_sp,line_dsp,line_par = plt_param3(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_sp.set_data([], [])
    line_dsp.set_data([], [])
    line_par.set_data([], [])
    return fig,ax1,ax2,line_og,line_sp,line_dsp,line_par

def animate(i):
    sp,wvl = meas.sp[igood][i],meas.wvl
    ax1.set_title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1493 = np.argmin(abs(wvln-1493.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    par3 = dsp[i1493]
    
    line_og.set_data(wvln,spc)
    line_sp.set_data([wvln[i1493],wvln[i1493]],[spc[i1493],spc[i1493]])
    line_dsp.set_data(wvln,dsp)
    line_par.set_data([wvln[i1493],wvln[i1493]],[par3,par3])
    
    line_par.set_label('Par3={:2.2f}'.format(par3))
    ax2.legend(frameon=False)
    
    return fig,ax1,ax2,line_og,line_sp,line_dsp,line_par

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[210]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 3 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))


# In[211]:


# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter3_{}.mp4'.format(daystr),writer=writer)


# ### Plot the animations for the LUT

# In[259]:


fig, ax1,ax2,line_og,line_sp,line_dsp,line_par = plt_param3(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_sp.set_data([], [])
    line_dsp.set_data([], [])
    line_par.set_data([], [])
    return fig,ax1,ax2,line_og,line_sp,line_dsp,line_par

def animate(i):
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    ax1.set_title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1493 = np.argmin(abs(wvln-1493.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    par3 = dsp[i1493]
    
    line_og.set_data(wvln,spc)
    line_sp.set_data([wvln[i1493],wvln[i1493]],[spc[i1493],spc[i1493]])
    line_dsp.set_data(wvln,dsp)
    line_par.set_data([wvln[i1493],wvln[i1493]],[par3,par3])
    
    line_par.set_label('Par3={:2.2f}'.format(par3))
    ax2.legend(frameon=False)
    
    return fig,ax1,ax2,line_og,line_sp,line_dsp,line_par

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[213]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 3 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))


# In[214]:


# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter3_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 4: Ratio between 1200 nm / 1237 nm

# ### Define Parameter 4

# In[243]:


def plt_param4(sp,wvl,fig=None):
    import Sp_parameters as Sp
    import numpy as np
    if not fig:
        fig = plt.figure()

    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i1198 = np.argmin(abs(wvln-1198.0))
    i1236 = np.argmin(abs(wvln-1236.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par4 = norm[i1198]/norm[i1236]
    
    
    line_og = plt.plot(wvln,norm,'.b',label='normed to max')
    line_lo = plt.plot([wvln[i1198],wvln[i1198]],[norm[i1198],norm[i1198]],'sg',label='norm at 1200 nm')
    line_loh = plt.axhline(norm[i1198],color='g')
    line_lov = plt.axvline(wvln[i1198],color='g')
    line_hi = plt.plot([wvln[i1236],wvln[i1236]],[norm[i1236],norm[i1236]],'sr',
                       label='norm at 1236 nm\nPar4:{:2.2f}'.format(par4))
    line_hih = plt.axhline(norm[i1236],color='r')
    line_hiv = plt.axvline(wvln[i1236],color='r')
    
    plt.xlim(1100,1300)
    plt.ylim(0,0.5)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Normed Radiance to max')
    plt.grid()
    
    plt.legend(frameon=False)
    return fig, line_og[0],line_lo[0],line_loh,line_hi[0],line_hih


# In[244]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param4(sp,wvl)


# ### Plot the measurements

# In[260]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param4(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_lo.set_data([], [])
    line_loh.set_data([],[])
    line_hi.set_data([], [])
    line_hih.set_data([],[])
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

def animate(i):
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i1198 = np.argmin(abs(wvln-1198.0))
    i1236 = np.argmin(abs(wvln-1236.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par4 = norm[i1198]/norm[i1236]
    
    line_og.set_data(wvln,norm)
    line_lo.set_data([wvln[i1198],wvln[i1198]],[norm[i1198],norm[i1198]])
    line_loh.set_data([0,1],[norm[i1198],norm[i1198]])
    line_hi.set_data([wvln[i1236],wvln[i1236]],[norm[i1236],norm[i1236]])
    line_hih.set_data([0,1],[norm[i1236],norm[i1236]])
    line_hi.set_label('norm at 1236 nm\nPar4:{:2.2f}'.format(par4))
    plt.legend(frameon=False)
    
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[247]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 4 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))


# In[248]:


# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter4_{}.mp4'.format(daystr),writer=writer)


# ### Plot the animations for the LUT

# In[261]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param4(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_lo.set_data([], [])
    line_loh.set_data([],[])
    line_hi.set_data([], [])
    line_hih.set_data([],[])
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

def animate(i):
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i1198 = np.argmin(abs(wvln-1198.0))
    i1236 = np.argmin(abs(wvln-1236.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par4 = norm[i1198]/norm[i1236]
    
    line_og.set_data(wvln,norm)
    line_lo.set_data([wvln[i1198],wvln[i1198]],[norm[i1198],norm[i1198]])
    line_loh.set_data([0,1],[norm[i1198],norm[i1198]])
    line_hi.set_data([wvln[i1236],wvln[i1236]],[norm[i1236],norm[i1236]])
    line_hih.set_data([0,1],[norm[i1236],norm[i1236]])
    line_hi.set_label('norm at 1236 nm\nPar4:{:2.2f}'.format(par4))
    plt.legend(frameon=False)
    
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[252]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 4 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))


# In[253]:


# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter4_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 5, 6, and 7: 
# mean of normalized rad between 1248 nm - 1270 nm; 1565 nm - 1644 nm; 1000 nm - 1050 nm

# ### Definition of parameters

# In[278]:


def plt_param567(sp,wvl,fig=None):
    import Sp_parameters as Sp
    import numpy as np
    if not fig:
        fig = plt.figure()

    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    [i1000,i1565,i1248,i1270,i1644,i1050,] = Sp.find_closest(wvln,np.array([1000,1565,1248,1270,1644,1050]))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par5 = np.nanmean(norm[i1248:i1270])
    par6 = np.nanmean(norm[i1565:i1644])
    par7 = np.nanmean(norm[i1000:i1050])
    
    line_og = plt.plot(wvln,norm,'.b',label='Radiance')
    line_p5 = plt.plot(wvln[i1248:i1270],norm[i1248:i1270],'-',color='chartreuse',lw=15,zorder=-1,label='Par 5 region')
    line_p6 = plt.plot(wvln[i1565:i1644],norm[i1565:i1644],'-',color='lightcoral',lw=15,zorder=-1,label='Par 6 region')
    line_p7 = plt.plot(wvln[i1000:i1050],norm[i1000:i1050],'-',color='aquamarine',lw=15,zorder=-1,label='Par 7 region')
    line_p5h = plt.plot([800.0,np.nanmean(wvln[i1248:i1270]),1800.0],[par5,par5,par5],'-s',color='green',lw=2,zorder=10,label='Par5={:2.2f}'.format(par5))
    line_p6h = plt.plot([800.0,np.nanmean(wvln[i1565:i1644]),1800.0],[par6,par6,par6],'-s',color='red',lw=2,zorder=10,label='Par6={:2.2f}'.format(par6))
    line_p7h = plt.plot([800.0,np.nanmean(wvln[i1000:i1050]),1800.0],[par7,par7,par7],'-s',color='cyan',lw=2,zorder=10,label='Par7={:2.2f}'.format(par7))

    
    plt.xlim(950,1700)
    plt.ylim(0,0.6)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Normed Radiance to max')
    plt.grid()
    
    plt.legend(frameon=False)
    return fig, line_og[0],line_p5[0],line_p6[0],line_p7[0],line_p5h[0],line_p6h[0],line_p7h[0]


# In[279]:


plt_param567(sp,wvl)


# ### Plot the measurements

# In[288]:


fig, line_og,line_p5,line_p6,line_p7,line_p5h,line_p6h,line_p7h = plt_param567(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_p5.set_data([], [])
    line_p6.set_data([],[])
    line_p7.set_data([], [])
    line_p5h.set_data([], [])
    line_p6h.set_data([],[])
    line_p7h.set_data([], [])
    return fig, line_og,line_p5,line_p6,line_p7,line_p5h,line_p6h,line_p7h

def animate(i):
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    [i1000,i1565,i1248,i1270,i1644,i1050,] = Sp.find_closest(wvln,np.array([1000,1565,1248,1270,1644,1050]))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par5 = np.nanmean(norm[i1248:i1270])
    par6 = np.nanmean(norm[i1565:i1644])
    par7 = np.nanmean(norm[i1000:i1050])
    
    line_og.set_data(wvln,norm)
    line_p5.set_data(wvln[i1248:i1270],norm[i1248:i1270])
    line_p6.set_data(wvln[i1565:i1644],norm[i1565:i1644])
    line_p7.set_data(wvln[i1000:i1050],norm[i1000:i1050])
    line_p5h.set_data([800.0,np.nanmean(wvln[i1248:i1270]),1800.0],[par5,par5,par5])
    line_p6h.set_data([800.0,np.nanmean(wvln[i1565:i1644]),1800.0],[par6,par6,par6])
    line_p7h.set_data([800.0,np.nanmean(wvln[i1000:i1050]),1800.0],[par7,par7,par7])
    line_p5h.set_label('Par5={:2.2f}'.format(par5))
    line_p6h.set_label('Par6={:2.2f}'.format(par6))
    line_p7h.set_label('Par7={:2.2f}'.format(par7))
    plt.legend(frameon=False)
    
    return fig, line_og,line_p5,line_p6,line_p7,line_p5h,line_p6h,line_p7h

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[290]:


#HTML(anim.to_jshtml())


# In[284]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 5, 6, and 7 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter567_{}.mp4'.format(daystr),writer=writer)


# ### Plot the LUT animations

# In[292]:


fig, line_og,line_p5,line_p6,line_p7,line_p5h,line_p6h,line_p7h = plt_param567(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_p5.set_data([], [])
    line_p6.set_data([],[])
    line_p7.set_data([], [])
    line_p5h.set_data([], [])
    line_p6h.set_data([],[])
    line_p7h.set_data([], [])
    return fig, line_og,line_p5,line_p6,line_p7,line_p5h,line_p6h,line_p7h

def animate(i):
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    [i1000,i1565,i1248,i1270,i1644,i1050,] = Sp.find_closest(wvln,np.array([1000,1565,1248,1270,1644,1050]))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par5 = np.nanmean(norm[i1248:i1270])
    par6 = np.nanmean(norm[i1565:i1644])
    par7 = np.nanmean(norm[i1000:i1050])
    
    line_og.set_data(wvln,norm)
    line_p5.set_data(wvln[i1248:i1270],norm[i1248:i1270])
    line_p6.set_data(wvln[i1565:i1644],norm[i1565:i1644])
    line_p7.set_data(wvln[i1000:i1050],norm[i1000:i1050])
    line_p5h.set_data([800.0,np.nanmean(wvln[i1248:i1270]),1800.0],[par5,par5,par5])
    line_p6h.set_data([800.0,np.nanmean(wvln[i1565:i1644]),1800.0],[par6,par6,par6])
    line_p7h.set_data([800.0,np.nanmean(wvln[i1000:i1050]),1800.0],[par7,par7,par7])
    line_p5h.set_label('Par5={:2.2f}'.format(par5))
    line_p6h.set_label('Par6={:2.2f}'.format(par6))
    line_p7h.set_label('Par7={:2.2f}'.format(par7))
    plt.legend(frameon=False)
    
    return fig, line_og,line_p5,line_p6,line_p7,line_p5h,line_p6h,line_p7h

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[291]:


#HTML(anim.to_jshtml())


# In[293]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 5, 6, and 7 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter567_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 8: curvature between 1493 nm - 1600 nm
# of rad normed to 1000 nm

# ### Define Parameter 8

# In[297]:


def plt_param8(sp,wvl,fig=None):
    from Sp_parameters import nanmasked
    if not fig:
        fig = plt.figure()
    
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1493 = np.argmin(abs(wvln-1493.0))
    i1600 = np.argmin(abs(wvln-1600.0))
    norm = spc / spc[i1000]
    
    line_og = plt.plot(wvln,norm,'.b',label='original')
    line_nd = plt.plot([wvln[i1493],wvln[i1600]],[norm[i1493],norm[i1600]],'sg',label='ends')
    fit7 = np.polyfit(np.array([wvln[i1493],wvln[i1600]]),np.array([norm[i1493],norm[i1600]]),1)
    fit7_fn = np.poly1d(fit7)
    line_ft = plt.plot(wvln[i1493:i1600],fit7_fn(wvln[i1493:i1600]),'-k',lw=2,label='line fit')
    line_cv = plt.plot(wvln[i1493:i1600],norm[i1493:i1600],'-r',lw=2,label='curve')
        
    par8 = sum(norm[i1493:i1600]-fit7_fn(wvln[i1493:i1600]))
    fill_pr = plt.fill_between(wvln[i1493:i1600],fit7_fn(wvln[i1493:i1600]),norm[i1493:i1600],
                               color='grey',alpha=0.4,label='par8= {:2.2f}'.format(par8))
    
    plt.xlim([1450,1680])
    plt.ylim([0,0.5])
    plt.legend(frameon=False)
    return fig, line_og[0],line_nd[0],line_ft[0],line_cv[0],fill_pr


# In[298]:


plt_param8(sp,wvl)


# ### Plot the Measurement animations

# In[302]:


fig, line_og,line_nd,line_ft,line_cv,fill_pr = plt_param8(meas.sp[igood][0],meas.wvl)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Zenith radiance normalized from value at 1000 nm')

def init():
    line_og.set_data([], [])
    line_nd.set_data([], [])
    line_ft.set_data([], [])
    line_cv.set_data([], [])
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = p.vertices[:,1]*0.0
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

def animate(i):
    
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1493 = np.argmin(abs(wvln-1493.0))
    i1600 = np.argmin(abs(wvln-1600.0))
    norm = spc / spc[i1000]
    
    line_og.set_data(wvln,norm)
    line_nd.set_data([wvln[i1493],wvln[i1600]],[norm[i1493],norm[i1600]])
    fit7 = np.polyfit(np.array([wvln[i1493],wvln[i1600]]),np.array([norm[i1493],norm[i1600]]),1)
    fit7_fn = np.poly1d(fit7)
    line_ft.set_data(wvln[i1493:i1600],fit7_fn(wvln[i1493:i1600]))
    line_cv.set_data(wvln[i1493:i1600],norm[i1493:i1600])
    par8 = sum(norm[i1493:i1600]-fit7_fn(wvln[i1493:i1600]))
    #fill_pr = fig.gca().fill_between(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]),norm[i1000:i1077],color='grey',alpha=0.4,label='par1= {:2.2f}'.format(par1))
    
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = np.hstack([fit7_fn(wvln[i1493:i1600])[0],fit7_fn(wvln[i1493:i1600]),
                                 norm[i1493:i1600][-1],norm[i1493:i1600][:0:-1],norm[i1493:i1600][0],norm[i1493:i1600][0]])
    fill_pr.set_label('par8= {:2.2f}'.format(par8))
    plt.legend(frameon=False)
    
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[305]:


#HTML(anim.to_jshtml())


# In[304]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 8 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter8_{}.mp4'.format(daystr),writer=writer)


# ### Plot the LUT animations

# In[306]:


fig, line_og,line_nd,line_ft,line_cv,fill_pr = plt_param8(lut[1].sp[0,:,0,10,0],lut[1].wvl)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Zenith radiance normalized from value at 1000 nm')

def init():
    line_og.set_data([], [])
    line_nd.set_data([], [])
    line_ft.set_data([], [])
    line_cv.set_data([], [])
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = p.vertices[:,1]*0.0
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

def animate(i):
    
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1493 = np.argmin(abs(wvln-1493.0))
    i1600 = np.argmin(abs(wvln-1600.0))
    norm = spc / spc[i1000]
    
    line_og.set_data(wvln,norm)
    line_nd.set_data([wvln[i1493],wvln[i1600]],[norm[i1493],norm[i1600]])
    fit7 = np.polyfit(np.array([wvln[i1493],wvln[i1600]]),np.array([norm[i1493],norm[i1600]]),1)
    fit7_fn = np.poly1d(fit7)
    line_ft.set_data(wvln[i1493:i1600],fit7_fn(wvln[i1493:i1600]))
    line_cv.set_data(wvln[i1493:i1600],norm[i1493:i1600])
    par8 = sum(norm[i1493:i1600]-fit7_fn(wvln[i1493:i1600]))
    #fill_pr = fig.gca().fill_between(wvln[i1000:i1077],fit0_fn(wvln[i1000:i1077]),norm[i1000:i1077],color='grey',alpha=0.4,label='par1= {:2.2f}'.format(par1))
    
    p = fill_pr.get_paths()[0]
    p.vertices[:,1] = np.hstack([fit7_fn(wvln[i1493:i1600])[0],fit7_fn(wvln[i1493:i1600]),
                                 norm[i1493:i1600][-1],norm[i1493:i1600][:0:-1],norm[i1493:i1600][0],norm[i1493:i1600][0]])
    fill_pr.set_label('par8= {:2.2f}'.format(par8))
    plt.legend(frameon=False)
    
    return fig,line_og,line_nd,line_ft,line_cv,fill_pr

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[309]:


#HTML(anim.to_jshtml())


# In[308]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 8 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter8_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 9 and 10: slope of derivative between 1000 - 1077 and 1200 - 1300 nm

# ### Define parameters plot

# In[348]:


def plt_param910(sp,wvl,fig=None):
    from linfit import linfit
    from Sp_parameters import nanmasked
    if not fig:
        fig = plt.figure()

    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1077 = np.argmin(abs(wvln-1077.0))
    i1200 = np.argmin(abs(wvln-1200.0))
    i1300 = np.argmin(abs(wvln-1300.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    fit8,b8 = linfit(wvln[i1000:i1077],dsp[i1000:i1077]) 
    fit9,b9 = linfit(wvln[i1200:i1300],dsp[i1200:i1300])
    
    par9 = fit8[0]
    par10 = fit9[0]
    
    ax1 = fig.add_subplot(2,1,1)
    line_og = ax1.plot(wvln,spc,'.b',label='original')
    line_sp9 = ax1.plot(wvln[i1000:i1077],spc[i1000:i1077],'-',color='aquamarine',lw=15,zorder=-1,alpha=0.6)
    line_sp10 = ax1.plot(wvln[i1200:i1300],spc[i1200:i1300],'-',color='chartreuse',lw=15,zorder=-1,alpha=0.6)
    ax1.set_ylabel('Radiance [W/m$^2$/nm/sr]')
    ax1.set_ylim(0,200.0)
    
    ax2 = fig.add_subplot(2,1,2,sharex=ax1)
    line_dsp = ax2.plot(wvln,dsp,'.k',label='Differential')
    line_dsp9 = ax2.plot(wvln[i1000:i1077],dsp[i1000:i1077],'-',label='Par 9 range',color='aquamarine',lw=15,zorder=-1,alpha=0.6)
    line_dsp10 = ax2.plot(wvln[i1200:i1300],dsp[i1200:i1300],'-',label='Par 10 range',color='chartreuse',lw=15,zorder=-1,alpha=0.6)
    
    line_par9 = ax2.plot(wvln[i1000:i1077],fit8[0]*wvln[i1000:i1077]+fit8[1],'-',color='b',
                          zorder=10,label='Par9={:2.3f}'.format(par9))
    line_par10 = ax2.plot(wvln[i1200:i1300],fit9[0]*wvln[i1200:i1300]+fit9[1],'-',color='green',
                          zorder=10,label='Par10={:2.3f}'.format(par10))
    ax2.set_xlim(950,1500)
    ax2.set_ylim(-25,25)
    ax2.set_ylabel('Derivative of\nnormalized radiance')
    ax2.set_xlabel('Wavelength [nm]')
    
    plt.legend(frameon=True,loc=1)
    return fig, ax1,ax2,line_og[0],line_sp9[0],line_sp10[0],line_dsp[0],line_dsp9[0],line_dsp10[0],line_par9[0],line_par10[0]


# In[349]:


plt_param910(sp,wvl)


# ### Plot the measurement animation

# In[336]:


from linfit import linfit


# In[350]:


fig, ax1,ax2,line_og,line_sp9,line_sp10,line_dsp,line_dsp9,line_dsp10,line_par9,line_par10 = plt_param910(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_sp9.set_data([], [])
    line_sp10.set_data([], [])
    line_dsp.set_data([], [])
    line_dsp9.set_data([], [])
    line_dsp10.set_data([], [])
    line_par9.set_data([], [])
    line_par10.set_data([], [])
    return fig, ax1,ax2,line_og,line_sp9,line_sp10,line_dsp,line_dsp9,line_dsp10,line_par9,line_par10

def animate(i):
    sp,wvl = meas.sp[igood][i],meas.wvl
    ax1.set_title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = nanmasked(sp)
    i1000 = np.argmin(abs(wvln-1000.0))
    i1077 = np.argmin(abs(wvln-1077.0))
    i1200 = np.argmin(abs(wvln-1200.0))
    i1300 = np.argmin(abs(wvln-1300.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    fit8,b8 = linfit(wvln[i1000:i1077],dsp[i1000:i1077]) 
    fit9,b9 = linfit(wvln[i1200:i1300],dsp[i1200:i1300])
    
    par9 = fit8[0]
    par10 = fit9[0]
    
    line_og.set_data(wvln,spc)
    line_sp9.set_data(wvln[i1000:i1077],spc[i1000:i1077])
    line_sp10.set_data(wvln[i1200:i1300],spc[i1200:i1300])
    line_dsp.set_data(wvln,dsp)
    line_dsp9.set_data(wvln[i1000:i1077],dsp[i1000:i1077])
    line_dsp10.set_data(wvln[i1200:i1300],dsp[i1200:i1300])
    line_par9.set_data(wvln[i1000:i1077],fit8[0]*wvln[i1000:i1077]+fit8[1])
    line_par10.set_data(wvln[i1200:i1300],fit9[0]*wvln[i1200:i1300]+fit9[1])
    
    line_par9.set_label('Par9={:2.3f}'.format(par9))
    line_par10.set_label('Par10={:2.3f}'.format(par10))
    ax2.legend(frameon=True,loc=1)
    
    return fig, ax1,ax2,line_og,line_sp9,line_sp10,line_dsp,line_dsp9,line_dsp10,line_par9,line_par10

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[340]:


#HTML(anim.to_jshtml())


# In[351]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 9, 10 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter910_{}.mp4'.format(daystr),writer=writer)


# ### Plot the LUT animation

# In[352]:


fig, ax1,ax2,line_og,line_sp9,line_sp10,line_dsp,line_dsp9,line_dsp10,line_par9,line_par10 = plt_param910(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_sp9.set_data([], [])
    line_sp10.set_data([], [])
    line_dsp.set_data([], [])
    line_dsp9.set_data([], [])
    line_dsp10.set_data([], [])
    line_par9.set_data([], [])
    line_par10.set_data([], [])
    return fig, ax1,ax2,line_og,line_sp9,line_sp10,line_dsp,line_dsp9,line_dsp10,line_par9,line_par10

def animate(i):
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    ax1.set_title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1077 = np.argmin(abs(wvln-1077.0))
    i1200 = np.argmin(abs(wvln-1200.0))
    i1300 = np.argmin(abs(wvln-1300.0))
    norm2 = spc/spc[i1000]
    dsp = smooth(deriv(norm2,wvln/1000.0),2,nan=False,old=True)

    fit8,b8 = linfit(wvln[i1000:i1077],dsp[i1000:i1077]) 
    fit9,b9 = linfit(wvln[i1200:i1300],dsp[i1200:i1300])
    
    par9 = fit8[0]
    par10 = fit9[0]
    
    line_og.set_data(wvln,spc)
    line_sp9.set_data(wvln[i1000:i1077],spc[i1000:i1077])
    line_sp10.set_data(wvln[i1200:i1300],spc[i1200:i1300])
    line_dsp.set_data(wvln,dsp)
    line_dsp9.set_data(wvln[i1000:i1077],dsp[i1000:i1077])
    line_dsp10.set_data(wvln[i1200:i1300],dsp[i1200:i1300])
    line_par9.set_data(wvln[i1000:i1077],fit8[0]*wvln[i1000:i1077]+fit8[1])
    line_par10.set_data(wvln[i1200:i1300],fit9[0]*wvln[i1200:i1300]+fit9[1])
    
    line_par9.set_label('Par9={:2.3f}'.format(par9))
    line_par10.set_label('Par10={:2.3f}'.format(par10))
    ax2.legend(frameon=True,loc=1)
    
    return fig, ax1,ax2,line_og,line_sp9,line_sp10,line_dsp,line_dsp9,line_dsp10,line_par9,line_par10

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)

#HTML(anim.to_jshtml())


# In[355]:


#HTML(anim.to_jshtml())


# In[354]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 9 and 10 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter910_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 11: slope of normalized radiance between 530 nm - 610 nm

# ### Define parameters plot

# In[363]:


def plt_param11(sp,wvl,fig=None):
    from Sp_parameters import nanmasked
    if not fig:
        fig = plt.figure()
    
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i530 = np.argmin(abs(wvln-530.0))
    i610 = np.argmin(abs(wvln-610.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    fit10,z = linfit(wvln[i530:i610]/1000,norm[i530:i610])
    
    line_og = plt.plot(wvln,norm,'.b',label='original')
    line_sh = plt.plot(wvln[i530:i610],norm[i530:i610],'-',color='violet',alpha=0.6,lw=15,zorder=-1,label='Par 11 range')
    line_ft = plt.plot(wvln[i530:i610],wvln[i530:i610]/1000.0*fit10[0]+fit10[1],'-k',
                       label='line fit\nPar 11={:2.2f}'.format(fit10[0]))
    plt.xlim([430,700])
    plt.ylim([0.5,1.0])
    plt.grid()
    plt.legend(frameon=False)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Normalized radiance')
    return fig, line_og[0],line_sh[0],line_ft[0]


# In[364]:


plt_param11(sp,wvl)


# ### Plot the measurements animation

# In[366]:


fig, line_og,line_sh,line_ft = plt_param11(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_sh.set_data([], [])
    line_ft.set_data([], [])
    return fig,line_og,line_sh,line_ft

def animate(i):
    
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i530 = np.argmin(abs(wvln-530.0))
    i610 = np.argmin(abs(wvln-610.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    fit10,z = linfit(wvln[i530:i610]/1000,norm[i530:i610])
    
    line_og.set_data(wvln,norm)
    line_sh.set_data(wvln[i530:i610],norm[i530:i610])
    line_ft.set_data(wvln[i530:i610],wvln[i530:i610]/1000.0*fit10[0]+fit10[1])
    line_ft.set_label('line fit\nPar 11={:2.2f}'.format(fit10[0]))
    plt.legend(frameon=False)
    
    return fig,line_og,line_sh,line_ft

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)


# In[369]:


#HTML(anim.to_jshtml())


# In[368]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 11 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter11_{}.mp4'.format(daystr),writer=writer)


# ### Plot the LUT animation

# In[371]:


fig, line_og,line_sh,line_ft = plt_param11(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_sh.set_data([], [])
    line_ft.set_data([], [])
    return fig,line_og,line_sh,line_ft

def animate(i):
    
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i530 = np.argmin(abs(wvln-530.0))
    i610 = np.argmin(abs(wvln-610.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    fit10,z = linfit(wvln[i530:i610]/1000,norm[i530:i610])
    
    line_og.set_data(wvln,norm)
    line_sh.set_data(wvln[i530:i610],norm[i530:i610])
    line_ft.set_data(wvln[i530:i610],wvln[i530:i610]/1000.0*fit10[0]+fit10[1])
    line_ft.set_label('line fit\nPar 11={:2.2f}'.format(fit10[0]))
    plt.legend(frameon=False)
    
    return fig,line_og,line_sh,line_ft

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)


# In[376]:


#HTML(anim.to_jshtml())


# In[373]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 11 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter11_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 12: normalized radiance at 1040 nm

# ### Define the parameter plot

# In[381]:


def plt_param12(sp,wvl,fig=None):
    from Sp_parameters import nanmasked
    if not fig:
        fig = plt.figure()
    
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1040 = np.argmin(abs(wvln-1040.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    par12 = norm[i1040]
    
    line_og = plt.plot(wvln,norm,'.k',label='original')
    line_ft = plt.plot([800.0,1040.0,1300.0],[par12,par12,par12],'-sb',
                       label='Par 12={:2.2f}'.format(par12))
    plt.xlim([950,1200])
    plt.ylim([0.0,0.6])
    plt.grid()
    plt.legend(frameon=False)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Normalized radiance')
    return fig, line_og[0],line_ft[0]


# In[382]:


plt_param12(sp,wvl)


# ### Plot the measurements animation

# In[384]:


fig, line_og,line_ft = plt_param12(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_ft.set_data([], [])
    return fig,line_og,line_ft

def animate(i):
    
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1040 = np.argmin(abs(wvln-1040.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    par12 = norm[i1040]
    
    line_og.set_data(wvln,norm)
    line_ft.set_data([800.0,1040.0,1300.0],[par12,par12,par12])
    line_ft.set_label('Par 12={:2.2f}'.format(par12))
    plt.legend(frameon=False)
    
    return fig,line_og,line_ft

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)


# In[387]:


#HTML(anim.to_jshtml())


# In[386]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 12 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter12_{}.mp4'.format(daystr),writer=writer)


# ### Plot the LUT animations

# In[388]:


fig, line_og,line_ft = plt_param12(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_ft.set_data([], [])
    return fig,line_og,line_ft

def animate(i):
    
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1040 = np.argmin(abs(wvln-1040.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    par12 = norm[i1040]
    
    line_og.set_data(wvln,norm)
    line_ft.set_data([800.0,1040.0,1300.0],[par12,par12,par12])
    line_ft.set_label('Par 12={:2.2f}'.format(par12))
    plt.legend(frameon=False)
    
    return fig,line_og,line_ft

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)


# In[391]:


#HTML(anim.to_jshtml())


# In[390]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 12 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter12_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 13: ratio of normalized radiance at 1000 nm / 1065 nm

# ### Define the parameter plot

# In[397]:


def plt_param13(sp,wvl,fig=None):
    import Sp_parameters as Sp
    import numpy as np
    if not fig:
        fig = plt.figure()

    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1065 = np.argmin(abs(wvln-1065.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par13 = norm[i1000]/norm[i1065]
    
    
    line_og = plt.plot(wvln,norm,'.b',label='normed to max')
    line_lo = plt.plot([wvln[i1000],wvln[i1000]],[norm[i1000],norm[i1000]],'sg',label='norm at 1000 nm')
    line_loh = plt.axhline(norm[i1000],color='g')
    line_lov = plt.axvline(wvln[i1000],color='g')
    line_hi = plt.plot([wvln[i1065],wvln[i1065]],[norm[i1065],norm[i1065]],'sr',
                       label='norm at 1065 nm\nPar 13:{:2.3f}'.format(par13))
    line_hih = plt.axhline(norm[i1065],color='r')
    line_hiv = plt.axvline(wvln[i1065],color='r')
    
    plt.xlim(950,1200)
    plt.ylim(0,0.6)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Normed Radiance to max')
    plt.grid()
    
    plt.legend(frameon=False)
    return fig, line_og[0],line_lo[0],line_loh,line_hi[0],line_hih


# In[393]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param13(sp,wvl)


# ### Plot the animations of the measurements 

# In[398]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param13(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_lo.set_data([], [])
    line_loh.set_data([],[])
    line_hi.set_data([], [])
    line_hih.set_data([],[])
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

def animate(i):
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1065 = np.argmin(abs(wvln-1065.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par13 = norm[i1000]/norm[i1065]
    
    line_og.set_data(wvln,norm)
    line_lo.set_data([wvln[i1000],wvln[i1000]],[norm[i1000],norm[i1000]])
    line_loh.set_data([0,1],[norm[i1000],norm[i1000]])
    line_hi.set_data([wvln[i1065],wvln[i1065]],[norm[i1065],norm[i1065]])
    line_hih.set_data([0,1],[norm[i1065],norm[i1065]])
    line_hi.set_label('norm at 1065 nm\nPar 13:{:2.3f}'.format(par13))
    plt.legend(frameon=False)
    
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)


# In[400]:


#HTML(anim.to_jshtml())


# In[399]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 13 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter13_{}.mp4'.format(daystr),writer=writer)


# ### Plot the LUT animations

# In[401]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param13(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_lo.set_data([], [])
    line_loh.set_data([],[])
    line_hi.set_data([], [])
    line_hih.set_data([],[])
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

def animate(i):
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i1000 = np.argmin(abs(wvln-1000.0))
    i1065 = np.argmin(abs(wvln-1065.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par13 = norm[i1000]/norm[i1065]
    
    line_og.set_data(wvln,norm)
    line_lo.set_data([wvln[i1000],wvln[i1000]],[norm[i1000],norm[i1000]])
    line_loh.set_data([0,1],[norm[i1000],norm[i1000]])
    line_hi.set_data([wvln[i1065],wvln[i1065]],[norm[i1065],norm[i1065]])
    line_hih.set_data([0,1],[norm[i1065],norm[i1065]])
    line_hi.set_label('norm at 1065 nm\nPar 13:{:2.3f}'.format(par13))
    plt.legend(frameon=False)
    
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)


# In[404]:


#HTML(anim.to_jshtml())


# In[403]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 13 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter13_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 14: ratio of normalized radiance at 600 nm / 870 nm

# ### Define the parameter plot

# In[411]:


def plt_param14(sp,wvl,fig=None):
    import Sp_parameters as Sp
    import numpy as np
    if not fig:
        fig = plt.figure()

    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i600 = np.argmin(abs(wvln-600.0))
    i870 = np.argmin(abs(wvln-870.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par14 = norm[i600]/norm[i870]
    
    
    line_og = plt.plot(wvln,norm,'.b',label='normed to max')
    line_lo = plt.plot([wvln[i600],wvln[i600]],[norm[i600],norm[i600]],'sg',label='norm at 600 nm')
    line_loh = plt.axhline(norm[i600],color='g')
    line_lov = plt.axvline(wvln[i600],color='g')
    line_hi = plt.plot([wvln[i870],wvln[i870]],[norm[i870],norm[i870]],'sr',
                       label='norm at 870 nm\nPar 14:{:2.3f}'.format(par14))
    line_hih = plt.axhline(norm[i870],color='r')
    line_hiv = plt.axvline(wvln[i870],color='r')
    
    plt.xlim(550,900)
    plt.ylim(0,1.0)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Normed Radiance to max')
    plt.grid()
    
    plt.legend(frameon=False)
    return fig, line_og[0],line_lo[0],line_loh,line_hi[0],line_hih


# In[412]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param14(sp,wvl)


# ### Plot the animations of the measurements 

# In[417]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param14(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_lo.set_data([], [])
    line_loh.set_data([],[])
    line_hi.set_data([], [])
    line_hih.set_data([],[])
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

def animate(i):
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i600 = np.argmin(abs(wvln-600.0))
    i870 = np.argmin(abs(wvln-870.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par14 = norm[i600]/norm[i870]
    
    line_og.set_data(wvln,norm)
    line_lo.set_data([wvln[i600],wvln[i600]],[norm[i600],norm[i600]])
    line_loh.set_data([0,1],[norm[i600],norm[i600]])
    line_hi.set_data([wvln[i870],wvln[i870]],[norm[i870],norm[i870]])
    line_hih.set_data([0,1],[norm[i870],norm[i870]])
    line_hi.set_label('norm at 870 nm\nPar 14:{:2.3f}'.format(par14))
    plt.legend(frameon=False,loc=3)
    
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)


# In[420]:


#HTML(anim.to_jshtml())


# In[419]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 14 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter14_{}.mp4'.format(daystr),writer=writer)


# ### Plot the LUT animations

# In[421]:


fig,line_og,line_lo,line_loh,line_hi,line_hih = plt_param14(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_lo.set_data([], [])
    line_loh.set_data([],[])
    line_hi.set_data([], [])
    line_hih.set_data([],[])
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

def animate(i):
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = Sp.nanmasked(sp)
    wvln = wvl[mask]
    i600 = np.argmin(abs(wvln-600.0))
    i870 = np.argmin(abs(wvln-870.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws]) 
    par14 = norm[i600]/norm[i870]
    
    line_og.set_data(wvln,norm)
    line_lo.set_data([wvln[i600],wvln[i600]],[norm[i600],norm[i600]])
    line_loh.set_data([0,1],[norm[i600],norm[i600]])
    line_hi.set_data([wvln[i870],wvln[i870]],[norm[i870],norm[i870]])
    line_hih.set_data([0,1],[norm[i870],norm[i870]])
    line_hi.set_label('norm at 870 nm\nPar 14:{:2.3f}'.format(par14))
    plt.legend(frameon=False,loc=3)
    
    return fig,line_og,line_lo,line_loh,line_hi,line_hih

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)


# In[424]:


#HTML(anim.to_jshtml())


# In[423]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 14 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter14_LUT_{}.mp4'.format(vlut),writer=writer)


# ## Parameter 15: slope of radiance / rad at 1565 between 1565 nm - 1634 nm

# ### Define parameters plot

# In[447]:


def plt_param15(sp,wvl,fig=None):
    from Sp_parameters import nanmasked
    if not fig:
        fig = plt.figure()
    
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1565 = np.argmin(abs(wvln-1565.0))
    i1634 = np.argmin(abs(wvln-1634.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    fit14,z = linfit(wvl[i1565:i1634],spc[i1565:i1634]/norm[i1565])
    
    line_og = plt.plot(wvln,spc/norm[i1565],'.b',label='original')
    line_sh = plt.plot(wvln[i1565:i1634],spc[i1565:i1634]/norm[i1565],'-',color='lightcoral',alpha=0.6,lw=15,zorder=-1,
                       label='Par 15 range')
    line_ft = plt.plot(wvln[i1565:i1634],wvln[i1565:i1634]*fit14[0]+fit14[1],'-k',
                       label='line fit\nPar 15={:2.2f}'.format(fit14[0]))
    plt.xlim([1520,1650])
    plt.ylim([0,500.0])
    plt.grid()
    plt.legend(frameon=False)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Radiance [W/m$^2$/nm/sr] normalized to norm at 1565 nm')
    return fig, line_og[0],line_sh[0],line_ft[0]


# In[448]:


plt_param15(sp,wvl)


# ### Plot the measurements animation

# In[450]:


fig, line_og,line_sh,line_ft = plt_param15(meas.sp[igood][0],meas.wvl)

def init():
    line_og.set_data([], [])
    line_sh.set_data([], [])
    line_ft.set_data([], [])
    return fig,line_og,line_sh,line_ft

def animate(i):
    
    sp,wvl = meas.sp[igood][i],meas.wvl
    plt.title('For {} at UTC: {:2.4f}, sza:{:2.3f} and i: {}'.format(daystr,meas.utc[igood][i],meas.sza[igood][i][0],i))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1565 = np.argmin(abs(wvln-1565.0))
    i1634 = np.argmin(abs(wvln-1634.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    fit14,z = linfit(wvl[i1565:i1634],spc[i1565:i1634]/norm[i1565])
    
    line_og.set_data(wvln,spc/norm[i1565])
    line_sh.set_data(wvln[i1565:i1634],spc[i1565:i1634]/norm[i1565])
    line_ft.set_data(wvln[i1565:i1634],wvln[i1565:i1634]*fit14[0]+fit14[1])
    line_ft.set_label('line fit\nPar 15={:2.2f}'.format(fit14[0]))
    plt.legend(frameon=False)
    
    return fig,line_og,line_sh,line_ft

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(meas.utc[igood]), interval=10, blit=True)


# In[443]:


#HTML(anim.to_jshtml())


# In[451]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 15 animations for {}'.format(daystr),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter15_{}.mp4'.format(daystr),writer=writer)


# ### Plot the LUT animation

# In[452]:


fig, line_og,line_sh,line_ft = plt_param15(lut[1].sp[0,:,0,10,0],lut[1].wvl)

def init():
    line_og.set_data([], [])
    line_sh.set_data([], [])
    line_ft.set_data([], [])
    return fig,line_og,line_sh,line_ft

def animate(i):
    
    sp,wvl = lut[1].sp[0,:,0,10,i],lut[1].wvl
    plt.title('Modeled For COD: {}, REF: {} $\\mu$m and i: {} at sza:{}'.format(lut[1].tausp[i],lut[1].refsp[10],i,lut[1].sza[1]))
    spc, mask = nanmasked(sp)
    wvln = wvl[mask]
    i1565 = np.argmin(abs(wvln-1565.0))
    i1634 = np.argmin(abs(wvln-1634.0))
    iws = np.where((wvln>=315.0)&(wvln<=940.0))[0]
    norm = spc/np.nanmax(spc[iws])
    
    fit14,z = linfit(wvl[i1565:i1634],spc[i1565:i1634]/norm[i1565])
    
    line_og.set_data(wvln,spc/norm[i1565])
    line_sh.set_data(wvln[i1565:i1634],spc[i1565:i1634]/norm[i1565])
    line_ft.set_data(wvln[i1565:i1634],wvln[i1565:i1634]*fit14[0]+fit14[1])
    line_ft.set_label('line fit\nPar 15={:2.2f}'.format(fit14[0]))
    plt.legend(frameon=False)
    
    return fig,line_og,line_sh,line_ft

anim = animation.FuncAnimation(fig, animate,init_func=init,frames=len(lut[1].tausp), interval=10, blit=True)


# In[454]:


#HTML(anim.to_jshtml())


# In[453]:


metadata = dict(artist='Samuel LeBlanc',title='Parameter 15 animations for LUT file {} for sza: {}'.format(fp_lut_mat.split('/')[-1],lut[1].sza[1]),
                comment='Created on {} with the program {}'.format(datetime.now(),nb_name))
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=metadata, bitrate=1800)
anim.save(fp+'starzen/Parameter15_LUT_{}.mp4'.format(vlut),writer=writer)


# In[ ]:




