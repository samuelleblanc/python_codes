
# coding: utf-8

# # Info
# Name:  
# 
#     Explor_fuliou_DARE_single_sol
# 
# Purpose:  
# 
#     Explore the results of the single solutions calulations from Run_fuliou 
#   
# Input:
# 
#     none at command line
#     see methods of module
# 
# Output:
#    
#     plots
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - math
#     - pdb
#     - datetime
#     - load_utils
#   
# Needed Files:
# 
#     - MOC_1solx_DARE_{vv}_{num}.mat
#   
#   
# Modification History:
# 
#     Wrtten: Samuel LeBlanc, NASA Ames, from Santa Cruz, 2017-03-28
#     Modified: 

# # Load the required modules

# In[1]:

import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')


# In[2]:

import matplotlib.cm as cm


# In[3]:

import numpy as np


# In[4]:

fp = 'C:\\Users\\sleblan2\\Research\\Calipso\\moc\\MOCsolutions_individual\\'


# In[5]:

vv = 'v1'


# In[6]:

num=19373


# # Load the files

# In[7]:

s = sio.loadmat(fp+'MOC_1solx_DARE_{vv}_{num}.mat'.format(vv=vv,num=num))


# In[8]:

s.keys()


# In[9]:

s['solutions'][0,0].dtype.names


# # Format results in useable arrays

# In[10]:

s['solutions'][0,1]['dF_toa_24hr'][0,0][0,0]


# In[11]:

dftoa = np.array([s['solutions'][0,i]['dF_toa_24hr'][0,0][0,0] for i in xrange(len(s['solutions'][0,:]))])
dfsfc = np.array([s['solutions'][0,i]['dF_sfc_24hr'][0,0] for i in xrange(len(s['solutions'][0,:]))])[:,0,0]
dftoai = np.array([s['solutions'][0,i]['dF_toa_instant'][0,0][0,0] for i in xrange(len(s['solutions'][0,:]))])
dfsfci = np.array([s['solutions'][0,i]['dF_sfc_instant'][0,0] for i in xrange(len(s['solutions'][0,:]))])[:,0,0]


# In[12]:

dftoa.shape


# In[13]:

s['select'][0,0][1]['dF_sfc_24hr'][0][0][0,0]


# In[60]:

s['select']['s0s0s0'][0][0]['dF_toa_24hr'][0][0][0][0]


# In[61]:

toa_sel = {}
sfc_sel = {}
toa_seli = {}
sfc_seli = {}
i_str = ['m1','s0','p1']
for ie in [-1,0,1]:
    for ia in [-1,0,1]:
        for im in [-1,0,1]:
            form = {'num':num,'e':i_str[ie+1],'a':i_str[ia+1],'s':i_str[im+1]}
            val = '{e}{a}{s}'.format(**form)
            toa_sel[val] = s['select'][val][0][0]['dF_toa_24hr'][0][0][0][0]
            sfc_sel[val] = s['select'][val][0][0]['dF_sfc_24hr'][0][0][0][0]
            toa_seli[val] = s['select'][val][0][0]['dF_toa_instant'][0][0][0][0]
            sfc_seli[val] = s['select'][val][0][0]['dF_sfc_instant'][0][0][0][0]


# In[62]:

sfc_sel.keys()


# ## set it up for the rouces values

# In[63]:

len(s['solutions'][0,:])


# In[64]:

s['solutions'][0,1]['ssa'][0,0][0]


# In[65]:

ssa = np.array([s['solutions'][0,i]['ssa'][0,0][0,8] for i in xrange(len(s['solutions'][0,:]))])


# In[66]:

ssa[0]


# In[20]:

ext = np.array([s['solutions'][0,i]['ext'][0,0][0,8] for i in xrange(len(s['solutions'][0,:]))])


# In[21]:

asy = np.array([s['solutions'][0,i]['asy'][0,0][0,8] for i in xrange(len(s['solutions'][0,:]))])


# In[33]:

s['solutions'][0,0].dtype.names


# # Start to plot out results

# ## plot the source distribution

# In[22]:

fig,ax = plt.subplots(1,3,figsize=(12,6))

ax[0].hist(ext,bins=40,edgecolor='None')
ax[0].axvline(np.mean(ext),color='k',lw=3,label='Mean')
ax[0].axvline(np.median(ext),color='grey',ls='--',lw=2,label='Median')
ax[0].axvspan(np.mean(ext)-np.std(ext), np.mean(ext)+np.std(ext), alpha=0.2, color='red',label='Std')
ax[0].set_ylabel('Number of points')
ax[0].set_xlabel('Extinction [km$^{{-1}}$]')

ax[1].hist(ssa,bins=40,edgecolor='None')
ax[1].axvline(np.mean(ssa),color='k',lw=3,label='Mean')
ax[1].axvline(np.median(ssa),color='grey',ls='--',lw=2,label='Median')
ax[1].axvspan(np.mean(ssa)-np.std(ssa), np.mean(ssa)+np.std(ssa), alpha=0.2, color='red',label='Std')
ax[1].set_xlabel('Single Scattering Albedo')

ax[2].hist(asy,bins=40,edgecolor='None')
ax[2].axvline(np.mean(asy),color='k',lw=3,label='Mean')
ax[2].axvline(np.median(asy),color='grey',ls='--',lw=2,label='Median')
ax[2].axvspan(np.mean(asy)-np.std(asy), np.mean(asy)+np.std(asy), alpha=0.2, color='red',label='Std')
ax[2].set_xlabel('Asymmetry Parameter')

fig.suptitle('MOC Solutions for single pixel #{num}'.format(num=num))

plt.legend(loc=2,frameon=False)
plt.savefig(fp+'\\plot\\MOC_solutions_hist_{num}.png'.format(num=num),dpi=600,transparent=True)


# ## plot the DARE from theses sources

# In[67]:

dftoa


# In[68]:

plt.figure()
plt.plot(dftoa)


# In[25]:

a = cm.hsv(np.arange(10))


# In[26]:

cs = cm.hsv(np.arange(27)/27.0)


# In[27]:

for i,k in enumerate(toa_sel.keys()):
    print tuple(cs[i,:])


# ### Plot the 24h averages

# In[137]:

plt.figure(figsize=(12,5))
plt.hist(dftoa,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoa),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoa),color='grey',ls='--',lw=2,label='Median')
cs = cm.hsv(np.arange(27)/27.0)
ms = ['s','*','x','d']
for i,k in enumerate(toa_sel.keys()):
    plt.axvline(toa_sel[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dftoa)-np.std(dftoa), np.mean(dftoa)+np.std(dftoa), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE TOA [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE TOA 24hr average for all solutions on single pixel #{}'.format(num))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)

plt.savefig(fp+'plot\\DARE_TOA_24h_{num}.png'.format(num=num),dpi=600,transparent=True)


# In[70]:

toa_sel['s0s0s0']


# In[73]:

toa_sel['s0p1s0']


# In[74]:

toa_sel['s0m1s0']


# In[136]:

plt.figure(figsize=(12,5))
plt.hist(dfsfc,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfc),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfc),color='grey',ls='--',lw=2,label='Median')
cs = cm.hsv(np.arange(27)/27.0)
ms = ['s','*','x','d']
for i,k in enumerate(sfc_sel.keys()):
    plt.axvline(sfc_sel[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dfsfc)-np.std(dfsfc), np.mean(dfsfc)+np.std(dfsfc), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE Surface [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE Surface 24hr average for all solutions on single pixel #{}'.format(num))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)

plt.savefig(fp+'plot\\DARE_SFC_24h_{num}.png'.format(num=num),dpi=600,transparent=True)


# In[76]:

plt.figure(figsize=(12,5))
plt.hist(dfsfc,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfc),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfc),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(27)/27.0)
ms = ['s','*','x','d']
#for i,k in enumerate(sfc_sel.keys()):
#    plt.axvline(sfc_sel[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dfsfc)-np.std(dfsfc), np.mean(dfsfc)+np.std(dfsfc), alpha=0.2, color='red',label='Std')
plt.xlim(-55,-20)
plt.xlabel('DARE Surface [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE Surface 24hr average for all solutions on single pixel #{}'.format(num))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)
plt.savefig(fp+'plot\\DARE_SFC_24h_{num}_nosup.png'.format(num=num),dpi=600,transparent=True)


# ### Plot the instantaneous values

# In[134]:

plt.figure(figsize=(12,5))
plt.hist(dftoai,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoai),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoai),color='grey',ls='--',lw=2,label='Median')
cs = cm.hsv(np.arange(27)/27.0)
ms = ['s','*','x','d']
for i,k in enumerate(toa_sel.keys()):
    plt.axvline(toa_seli[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dftoai)-np.std(dftoai), np.mean(dftoai)+np.std(dftoai), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE TOA [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE TOA instantaneous for all solutions on single pixel #{}'.format(num))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)

plt.savefig(fp+'plot\\DARE_TOA_instantaneous_{num}.png'.format(num=num),dpi=600,transparent=True)


# In[135]:

plt.figure(figsize=(12,5))
plt.hist(dfsfci,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfci),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfci),color='grey',ls='--',lw=2,label='Median')
cs = cm.hsv(np.arange(27)/27.0)
ms = ['s','*','x','d']
for i,k in enumerate(sfc_seli.keys()):
    plt.axvline(sfc_seli[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dfsfci)-np.std(dfsfci), np.mean(dfsfci)+np.std(dfsfci), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE Surface [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE Surface instantaneous for all solutions on single pixel #{}'.format(num))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)

plt.savefig(fp+'plot\\DARE_SFC_instantaneous_{num}.png'.format(num=num),dpi=600,transparent=True)


# In[133]:

fp


# # Redo the analysis for number 22134

# In[138]:

num2 = 22134


# In[139]:

s2 = sio.loadmat(fp+'MOC_1solx_DARE_{vv}_{num}.mat'.format(vv=vv,num=num2))


# In[97]:

dftoa2 = np.array([s2['solutions'][0,i]['dF_toa_24hr'][0,0][0,0] for i in xrange(len(s2['solutions'][0,:]))])
dfsfc2 = np.array([s2['solutions'][0,i]['dF_sfc_24hr'][0,0] for i in xrange(len(s2['solutions'][0,:]))])[:,0,0]
dftoai2 = np.array([s2['solutions'][0,i]['dF_toa_instant'][0,0][0,0] for i in xrange(len(s2['solutions'][0,:]))])
dfsfci2 = np.array([s2['solutions'][0,i]['dF_sfc_instant'][0,0] for i in xrange(len(s2['solutions'][0,:]))])[:,0,0]


# In[98]:

toa_sel2 = {}
sfc_sel2 = {}
toa_seli2 = {}
sfc_seli2 = {}
i_str = ['m1','s0','p1']
for ie in [-1,0,1]:
    for ia in [-1,0,1]:
        for im in [-1,0,1]:
            form = {'num':num,'e':i_str[ie+1],'a':i_str[ia+1],'s':i_str[im+1]}
            val = '{e}{a}{s}'.format(**form)
            toa_sel2[val] = s2['select'][val][0][0]['dF_toa_24hr'][0][0][0][0]
            sfc_sel2[val] = s2['select'][val][0][0]['dF_sfc_24hr'][0][0][0][0]
            toa_seli2[val] = s2['select'][val][0][0]['dF_toa_instant'][0][0][0][0]
            sfc_seli2[val] = s2['select'][val][0][0]['dF_sfc_instant'][0][0][0][0]


# In[99]:

ssa2 = np.array([s2['solutions'][0,i]['ssa'][0,0][0,8] for i in xrange(len(s2['solutions'][0,:]))])
ext2 = np.array([s2['solutions'][0,i]['ext'][0,0][0,8] for i in xrange(len(s2['solutions'][0,:]))])
asy2 = np.array([s2['solutions'][0,i]['asy'][0,0][0,8] for i in xrange(len(s2['solutions'][0,:]))])


# In[100]:

fig,ax = plt.subplots(1,3,figsize=(12,6))

ax[0].hist(ext2,bins=40,edgecolor='None')
ax[0].axvline(np.mean(ext2),color='k',lw=3,label='Mean')
ax[0].axvline(np.median(ext2),color='grey',ls='--',lw=2,label='Median')
ax[0].axvspan(np.mean(ext2)-np.std(ext2), np.mean(ext2)+np.std(ext2), alpha=0.2, color='red',label='Std')
ax[0].set_ylabel('Number of points')
ax[0].set_xlabel('Extinction [km$^{{-1}}$]')

ax[1].hist(ssa2,bins=40,edgecolor='None')
ax[1].axvline(np.mean(ssa2),color='k',lw=3,label='Mean')
ax[1].axvline(np.median(ssa2),color='grey',ls='--',lw=2,label='Median')
ax[1].axvspan(np.mean(ssa2)-np.std(ssa2), np.mean(ssa2)+np.std(ssa2), alpha=0.2, color='red',label='Std')
ax[1].set_xlabel('Single Scattering Albedo')

ax[2].hist(asy2,bins=40,edgecolor='None')
ax[2].axvline(np.mean(asy2),color='k',lw=3,label='Mean')
ax[2].axvline(np.median(asy2),color='grey',ls='--',lw=2,label='Median')
ax[2].axvspan(np.mean(asy2)-np.std(asy2), np.mean(asy2)+np.std(asy2), alpha=0.2, color='red',label='Std')
ax[2].set_xlabel('Asymmetry Parameter')

fig.suptitle('MOC Solutions for single pixel #{num}'.format(num=num2))

plt.legend(loc=2,frameon=False)
plt.savefig(fp+'\\plot\\MOC_solutions_hist_{num}.png'.format(num=num2),dpi=600,transparent=True)


# In[140]:

plt.figure(figsize=(12,5))
plt.hist(dftoa2,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoa2),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoa2),color='grey',ls='--',lw=2,label='Median')
cs = cm.hsv(np.arange(27)/27.0)
ms = ['s','*','x','d']
for i,k in enumerate(toa_sel.keys()):
    plt.axvline(toa_sel2[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dftoa2)-np.std(dftoa2), np.mean(dftoa2)+np.std(dftoa2), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE TOA [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE TOA 24hr average for all solutions on single pixel #{}'.format(num2))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)

plt.savefig(fp+'plot\\DARE_TOA_24h_{num}.png'.format(num=num2),dpi=600,transparent=True)


# In[141]:

plt.figure(figsize=(12,5))
plt.hist(dfsfc2,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfc2),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfc2),color='grey',ls='--',lw=2,label='Median')
cs = cm.hsv(np.arange(27)/27.0)
ms = ['s','*','x','d']
for i,k in enumerate(sfc_sel2.keys()):
    plt.axvline(sfc_sel2[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dfsfc2)-np.std(dfsfc2), np.mean(dfsfc2)+np.std(dfsfc2), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE Surface [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE Surface 24hr average for all solutions on single pixel #{}'.format(num2))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)

plt.savefig(fp+'plot\\DARE_SFC_24h_{num}.png'.format(num=num2),dpi=600,transparent=True)


# In[103]:

len(ext2)


# In[142]:

plt.figure(figsize=(12,5))
plt.hist(dftoai2,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoai2),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoai2),color='grey',ls='--',lw=2,label='Median')
cs = cm.hsv(np.arange(27)/27.0)
ms = ['s','*','x','d']
for i,k in enumerate(toa_seli.keys()):
    plt.axvline(toa_seli2[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dftoai2)-np.std(dftoai2), np.mean(dftoai2)+np.std(dftoai2), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE TOA [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE TOA instantaneous for all solutions on single pixel #{}'.format(num2))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)

plt.savefig(fp+'plot\\DARE_TOA_instant_{num}.png'.format(num=num2),dpi=600,transparent=True)


# In[143]:

plt.figure(figsize=(12,5))
plt.hist(dfsfci2,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfci2),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfci2),color='grey',ls='--',lw=2,label='Median')
cs = cm.hsv(np.arange(27)/27.0)
ms = ['s','*','x','d']
for i,k in enumerate(sfc_seli2.keys()):
    plt.axvline(sfc_seli2[k],lw=1,label=k,ls='-',color=tuple(cs[i,:]),marker=ms[i%4],ms=10)
    
plt.axvspan(np.mean(dfsfci2)-np.std(dfsfci2), np.mean(dfsfci2)+np.std(dfsfci2), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE Surface [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE Surface instantaneous for all solutions on single pixel #{}'.format(num2))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=2)

plt.savefig(fp+'plot\\DARE_SFC_instant_{num}.png'.format(num=num2),dpi=600,transparent=True)


# # Print out all the resulting values

# In[106]:

print 'TOA {}'.format(num)
print 'Mean:',np.mean(dftoa)
print 'Median',np.median(dftoa)
print 'Mean+std:',np.mean(dftoa)+np.std(dftoa)
print 'Mean-std:',np.mean(dftoa)-np.std(dftoa)
for k in np.sort(toa_sel.keys()):
    print k,toa_sel[k]


# In[124]:

print 'TOA SFC {}'.format(num)
print 'Mean:',np.mean(dftoa),np.mean(dfsfc)
print 'Median',np.median(dftoa),np.median(dfsfc)
print 'Mean+std:',np.mean(dftoa)+np.std(dftoa),np.mean(dfsfc)+np.std(dfsfc)
print 'Mean-std:',np.mean(dftoa)-np.std(dftoa),np.mean(dfsfc)-np.std(dfsfc)
for k in np.sort(toa_sel.keys()):
    print k,toa_sel[k],sfc_sel[k]


# In[144]:

import pandas as pd


# In[183]:

pd.DataFrame.from_dict(toa_sel,'index')


# In[184]:

pd.DataFrame.from_dict(sfc_sel,'index')


# In[181]:

for k in np.sort(toa_sel.keys()):
    print k,sfc_sel[k]


# In[126]:

print 'TOA, SFC {}'.format(num2)
print 'Mean:',np.mean(dftoa2),np.mean(dfsfc2)
print 'Median',np.median(dftoa2),np.median(dfsfc2)
print 'Mean+std:',np.mean(dftoa2)+np.std(dftoa2),np.mean(dfsfc2)+np.std(dfsfc2)
print 'Mean-std:',np.mean(dftoa2)-np.std(dftoa2),np.mean(dfsfc2)-np.std(dfsfc2)
for k in np.sort(toa_sel2.keys()):
    print k,toa_sel2[k],sfc_sel2[k]


# In[185]:

pd.DataFrame.from_dict(toa_sel2,'index')


# In[186]:

pd.DataFrame.from_dict(sfc_sel2,'index')


# In[127]:

for k in np.sort(toa_sel2.keys()):
    print k,toa_sel2[k]


# In[128]:

for k in np.sort(toa_sel2.keys()):
    print k,sfc_sel2[k]


# # Select a few representative values and see their DARE

# In[112]:

ii = {'ssa':[],'asy':[],'ext':[],'ssa2':[],'asy2':[],'ext2':[]}
for r in ii.keys():
    v = locals()[r]
    ii[r] = [np.argmin(abs(v-np.nanmean(v))),
             np.argmin(abs(v-np.nanmean(v)+np.std(v))),
             np.argmin(abs(v-np.nanmean(v)-np.std(v)))]


# In[113]:

ii


# ## plot the representative DARES

# In[114]:

str_names = ['mean {nm}','mean {nm}-std','mean {nm}+std']


# In[115]:

plt.figure(figsize=(12,5))
plt.hist(dftoa,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoa),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoa),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(9)/9.0)
ms = ['s','*','x']
io = 0
for i,k in enumerate(['ssa','asy','ext']):
    for jj in [0,1,2]:
        plt.axvline(dftoa[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1
    
plt.axvspan(np.mean(dftoa)-np.std(dftoa), np.mean(dftoa)+np.std(dftoa), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE TOA [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE TOA 24hr average for all solutions on single pixel #{}'.format(num))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=1)

plt.savefig(fp+'plot\\DARE_TOA_24h_selected_{num}.png'.format(num=num),dpi=600,transparent=True)


# In[116]:

plt.figure(figsize=(12,5))
plt.hist(dftoa2,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoa2),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoa2),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(9)/9.0)
ms = ['s','*','x']
io = 0
for i,k in enumerate(['ssa2','asy2','ext2']):
    for jj in [0,1,2]:
        plt.axvline(dftoa2[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1
    
plt.axvspan(np.mean(dftoa2)-np.std(dftoa2), np.mean(dftoa2)+np.std(dftoa2), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE TOA [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE TOA 24hr average for all solutions on single pixel #{}'.format(num2))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=1)

plt.savefig(fp+'plot\\DARE_TOA_24h_selected_{num}.png'.format(num=num2),dpi=600,transparent=True)


# In[117]:

plt.figure(figsize=(12,5))
plt.hist(dfsfc,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfc),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfc),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(9)/9.0)
ms = ['s','*','x']
io = 0
for i,k in enumerate(['ssa','asy','ext']):
    for jj in [0,1,2]:
        plt.axvline(dfsfc[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1
    
plt.axvspan(np.mean(dfsfc)-np.std(dfsfc), np.mean(dfsfc)+np.std(dfsfc), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE Surface [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE Surface 24hr average for all solutions on single pixel #{}'.format(num))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=1)

plt.savefig(fp+'plot\\DARE_SFC_24h_selected_{num}.png'.format(num=num),dpi=600,transparent=True)


# In[118]:

plt.figure(figsize=(12,5))
plt.hist(dfsfc2,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfc2),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfc2),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(9)/9.0)
ms = ['s','*','x']
io = 0
for i,k in enumerate(['ssa2','asy2','ext2']):
    for jj in [0,1,2]:
        plt.axvline(dfsfc2[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1
    
plt.axvspan(np.mean(dfsfc2)-np.std(dfsfc2), np.mean(dfsfc2)+np.std(dfsfc2), alpha=0.2, color='red',label='Std')

plt.xlabel('DARE Surface [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE Surface 24hr average for all solutions on single pixel #{}'.format(num2))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=1)

plt.savefig(fp+'plot\\DARE_SFC_24h_selected_{num}.png'.format(num=num2),dpi=600,transparent=True)


# ## Print out the representative DARE

# In[119]:

print 'DARE TOA {}'.format(num)
print 'DARE Mean: \t{}'.format(np.nanmean(dftoa))
print 'DARE mean+std: \t{}'.format(np.nanmean(dftoa)+np.std(dftoa))
print 'DARE mean-std: \t{}'.format(np.nanmean(dftoa)-np.std(dftoa))

io = 0
for i,k in enumerate(['ssa','asy','ext']):
    for jj in [0,1,2]:
        print str_names[jj].format(nm=k)+': \t{}'.format(dftoa[ii[k][jj]])
        io +=1


# In[120]:

print 'DARE Surface {}'.format(num)
print 'DARE Mean: \t{}'.format(np.nanmean(dfsfc))
print 'DARE mean+std: \t{}'.format(np.nanmean(dfsfc)+np.std(dfsfc))
print 'DARE mean-std: \t{}'.format(np.nanmean(dfsfc)-np.std(dfsfc))

io = 0
for i,k in enumerate(['ssa','asy','ext']):
    for jj in [0,1,2]:
        print str_names[jj].format(nm=k)+': \t{}'.format(dfsfc[ii[k][jj]])
        io +=1


# In[111]:

print 'Properties {}: ssa, asy, ext'.format(num)
io = 0
for i,k in enumerate(['ssa','asy','ext']):
    v = locals()[k]
    for jj in [0,1,2]:
        print str_names[jj].format(nm=k)+        ': \t{s} \t{a} \t{e}'.format(s=ssa[ii[k][jj]],a=asy[ii[k][jj]],e=ext[ii[k][jj]])
        io +=1


# In[121]:

print 'Properties {}: ssa, asy, ext'.format(num2)
io = 0
for i,k in enumerate(['ssa','asy','ext']):
    v = locals()[k]
    for jj in [0,1,2]:
        print str_names[jj].format(nm=k)+        ': \t{s} \t{a} \t{e}'.format(s=ssa2[ii[k][jj]],a=asy2[ii[k][jj]],e=ext2[ii[k][jj]])
        io +=1


# In[122]:

print 'DARE TOA {}'.format(num2)
print 'DARE Mean: \t{}'.format(np.nanmean(dftoa2))
print 'DARE mean+std: \t{}'.format(np.nanmean(dftoa2)+np.std(dftoa2))
print 'DARE mean-std: \t{}'.format(np.nanmean(dftoa2)-np.std(dftoa2))

io = 0
for i,k in enumerate(['ssa2','asy2','ext2']):
    for jj in [0,1,2]:
        print str_names[jj].format(nm=k)+': \t{}'.format(dftoa2[ii[k][jj]])
        io +=1


# In[123]:

print 'DARE Surface {}'.format(num2)
print 'DARE Mean: \t{}'.format(np.nanmean(dfsfc2))
print 'DARE mean+std: \t{}'.format(np.nanmean(dfsfc2)+np.std(dfsfc2))
print 'DARE mean-std: \t{}'.format(np.nanmean(dfsfc2)-np.std(dfsfc2))

io = 0
for i,k in enumerate(['ssa2','asy2','ext2']):
    for jj in [0,1,2]:
        print str_names[jj].format(nm=k)+': \t{}'.format(dfsfc2[ii[k][jj]])
        io +=1


# # Inverse link from mean DARE to input properties

# In[129]:

ii['DARE_TOA'] = [np.argmin(abs(dftoa-np.nanmean(dftoa))),
                 np.argmin(abs(dftoa-np.nanmean(dftoa)+np.std(dftoa))),
                 np.argmin(abs(dftoa-np.nanmean(dftoa)-np.std(dftoa))),
                 np.argmin(abs(dftoa-np.nanmedian(dftoa)))]
ii['DARE_SFC'] = [np.argmin(abs(dfsfc-np.nanmean(dfsfc))),
                 np.argmin(abs(dfsfc-np.nanmean(dfsfc)+np.std(dfsfc))),
                 np.argmin(abs(dfsfc-np.nanmean(dfsfc)-np.std(dfsfc))),
                 np.argmin(abs(dftoa-np.nanmedian(dfsfc)))]
ii['DARE2_TOA'] = [np.argmin(abs(dftoa2-np.nanmean(dftoa2))),
             np.argmin(abs(dftoa2-np.nanmean(dftoa2)+np.std(dftoa2))),
             np.argmin(abs(dftoa2-np.nanmean(dftoa2)-np.std(dftoa2))),
                 np.argmin(abs(dftoa-np.nanmedian(dftoa2)))]
ii['DARE2_SFC'] = [np.argmin(abs(dfsfc2-np.nanmean(dfsfc2))),
             np.argmin(abs(dfsfc2-np.nanmean(dfsfc2)+np.std(dfsfc2))),
             np.argmin(abs(dfsfc2-np.nanmean(dfsfc2)-np.std(dfsfc2))),
                 np.argmin(abs(dftoa-np.nanmedian(dfsfc2)))]


# ## now plot the input with the linked DARE properties

# In[130]:

str_names = ['mean {nm}','mean {nm}-std','mean {nm}+std','median {nm}']


# In[131]:

fig,ax = plt.subplots(1,4,figsize=(15,6))

ax[0].hist(ext,bins=40,edgecolor='None')
ax[0].axvline(np.mean(ext),color='k',lw=3,label='Mean')
ax[0].axvline(np.median(ext),color='grey',ls='--',lw=2,label='Median')
ax[0].axvspan(np.mean(ext)-np.std(ext), np.mean(ext)+np.std(ext), alpha=0.2, color='red',label='Std')
ax[0].set_ylabel('Number of points')
ax[0].set_xlabel('Extinction [km$^{{-1}}$]')

cs = cm.jet(np.arange(9)/9.0)
ms = ['s','*','x','d']
io = 0
for i,k in enumerate(['DARE_TOA','DARE_SFC']):
    for jj in [0,1,2,3]:
        ax[0].axvline(ext[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1

ax[1].hist(ssa,bins=40,edgecolor='None')
ax[1].axvline(np.mean(ssa),color='k',lw=3)
ax[1].axvline(np.median(ssa),color='grey',ls='--',lw=2)
ax[1].axvspan(np.mean(ssa)-np.std(ssa), np.mean(ssa)+np.std(ssa), alpha=0.2, color='red')
ax[1].set_xlabel('Single Scattering Albedo')

io = 0
for i,k in enumerate(['DARE_TOA','DARE_SFC']):
    for jj in [0,1,2,3]:
        ax[1].axvline(ssa[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1
#plt.legend(bbox_to_anchor=(1.03,2.07),loc=2,ncol=1)

ax[2].hist(asy,bins=40,edgecolor='None')
ax[2].axvline(np.mean(asy),color='k',lw=3,label='Mean')
ax[2].axvline(np.median(asy),color='grey',ls='--',lw=2,label='Median')
ax[2].axvspan(np.mean(asy)-np.std(asy), np.mean(asy)+np.std(asy), alpha=0.2, color='red',label='Std')
ax[2].set_xlabel('Asymmetry Parameter')

io = 0
for i,k in enumerate(['DARE_TOA','DARE_SFC']):
    for jj in [0,1,2,3]:
        ax[2].axvline(asy[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1

fig.suptitle('MOC Solutions for single pixel #{num}'.format(num=num))

fig.delaxes(ax[3])

#plt.legend(loc=2,frameon=False)
plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=1)
plt.savefig(fp+'\\plot\\MOC_solutions_hist_linked_{num}.png'.format(num=num),dpi=600,transparent=True)


# In[132]:

fig,ax = plt.subplots(1,4,figsize=(15,6))

ax[0].hist(ext2,bins=40,edgecolor='None')
ax[0].axvline(np.mean(ext2),color='k',lw=3,label='Mean')
ax[0].axvline(np.median(ext2),color='grey',ls='--',lw=2,label='Median')
ax[0].axvspan(np.mean(ext2)-np.std(ext2), np.mean(ext2)+np.std(ext2), alpha=0.2, color='red',label='Std')
ax[0].set_ylabel('Number of points')
ax[0].set_xlabel('Extinction [km$^{{-1}}$]')

cs = cm.jet(np.arange(9)/9.0)
ms = ['s','*','x','d']
io = 0
for i,k in enumerate(['DARE2_TOA','DARE2_SFC']):
    for jj in [0,1,2,3]:
        ax[0].axvline(ext2[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1

ax[1].hist(ssa2,bins=40,edgecolor='None')
ax[1].axvline(np.mean(ssa2),color='k',lw=3)
ax[1].axvline(np.median(ssa2),color='grey',ls='--',lw=2)
ax[1].axvspan(np.mean(ssa2)-np.std(ssa2), np.mean(ssa2)+np.std(ssa2), alpha=0.2, color='red')
ax[1].set_xlabel('Single Scattering Albedo')

io = 0
for i,k in enumerate(['DARE2_TOA','DARE2_SFC']):
    for jj in [0,1,2,3]:
        ax[1].axvline(ssa2[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1
#plt.legend(bbox_to_anchor=(1.03,2.07),loc=2,ncol=1)

ax[2].hist(asy2,bins=40,edgecolor='None')
ax[2].axvline(np.mean(asy2),color='k',lw=3,label='Mean')
ax[2].axvline(np.median(asy2),color='grey',ls='--',lw=2,label='Median')
ax[2].axvspan(np.mean(asy2)-np.std(asy2), np.mean(asy2)+np.std(asy2), alpha=0.2, color='red',label='Std')
ax[2].set_xlabel('Asymmetry Parameter')

io = 0
for i,k in enumerate(['DARE2_TOA','DARE2_SFC']):
    for jj in [0,1,2,3]:
        ax[2].axvline(asy2[ii[k][jj]],lw=1,label=str_names[jj].format(nm=k),
                    ls='-',color=tuple(cs[io,:]),marker=ms[jj],ms=10)
        io +=1

fig.suptitle('MOC Solutions for single pixel #{num}'.format(num=num2))

fig.delaxes(ax[3])

#plt.legend(loc=2,frameon=False)
plt.legend(bbox_to_anchor=(1.03,1.07),loc=2,ncol=1)
plt.savefig(fp+'\\plot\\MOC_solutions_hist_linked_{num}.png'.format(num=num2),dpi=600,transparent=True)


# In[ ]:



