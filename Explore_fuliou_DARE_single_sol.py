
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


# In[106]:

import matplotlib.cm as cm


# In[13]:

import numpy as np


# In[2]:

fp = 'C:\\Users\\sleblan2\\Research\\Calipso\\moc\\MOCsolutions_individual\\'


# In[5]:

vv = 'v1'


# In[77]:

num=19373


# # Load the files

# In[6]:

s = sio.loadmat(fp+'MOC_1solx_DARE_{vv}_{num}.mat'.format(vv=vv,num=num))


# In[7]:

s.keys()


# In[11]:

s['solutions'][0,0].dtype.names


# # Format results in useable arrays

# In[32]:

s['solutions'][0,1]['dF_toa_24hr'][0,0][0,0]


# In[33]:

dftoa = np.array([s['solutions'][0,i]['dF_toa_24hr'][0,0][0,0] for i in xrange(len(s['solutions'][0,:]))])
dfsfc = np.array([s['solutions'][0,i]['dF_sfc_24hr'][0,0] for i in xrange(len(s['solutions'][0,:]))])[:,0,0]
dftoai = np.array([s['solutions'][0,i]['dF_toa_instant'][0,0][0,0] for i in xrange(len(s['solutions'][0,:]))])
dfsfci = np.array([s['solutions'][0,i]['dF_sfc_instant'][0,0] for i in xrange(len(s['solutions'][0,:]))])[:,0,0]


# In[36]:

dftoa.shape


# In[80]:

s['select'][0,0][1]['dF_sfc_24hr'][0][0][0,0]


# In[81]:

toa_sel = {}
sfc_sel = {}
toa_seli = {}
sfc_seli = {}
i_str = ['m1','s0','p1']
ii = 0
for ie in [-1,0,1]:
    for ia in [-1,0,1]:
        for im in [-1,0,1]:
            form = {'num':num,'e':i_str[ie+1],'a':i_str[ia+1],'s':i_str[im+1]}
            val = '{e}{a}{s}'.format(**form)
            toa_sel[val] = s['select'][0,0][ii]['dF_toa_24hr'][0][0][0,0]
            sfc_sel[val] = s['select'][0,0][ii]['dF_sfc_24hr'][0][0][0,0]
            toa_seli[val] = s['select'][0,0][ii]['dF_toa_instant'][0][0][0,0]
            sfc_seli[val] = s['select'][0,0][ii]['dF_sfc_instant'][0][0][0,0]
            ii += 1


# In[83]:

sfc_sel.keys()


# ## set it up for the rouces values

# In[211]:

len(s['solutions'][0,:])


# In[173]:

s['solutions'][0,1]['ssa'][0,0][0]


# In[174]:

ssa = np.array([s['solutions'][0,i]['ssa'][0,0][0,8] for i in xrange(len(s['solutions'][0,:]))])


# In[175]:

ssa[0]


# In[176]:

ext = np.array([s['solutions'][0,i]['ext'][0,0][0,8] for i in xrange(len(s['solutions'][0,:]))])


# In[177]:

asy = np.array([s['solutions'][0,i]['asy'][0,0][0,8] for i in xrange(len(s['solutions'][0,:]))])


# # Start to plot out results

# ## plot the source distribution

# In[208]:

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

# In[37]:

dftoa


# In[38]:

plt.figure()
plt.plot(dftoa)


# In[118]:

a = cm.hsv(np.arange(10))


# In[125]:

cs = cm.hsv(np.arange(27)/27.0)


# In[126]:

for i,k in enumerate(toa_sel.keys()):
    print tuple(cs[i,:])


# ### Plot the 24h averages

# In[209]:

plt.figure(figsize=(12,5))
plt.hist(dftoa,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoa),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoa),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(27)/27.0)
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


# In[210]:

plt.figure(figsize=(12,5))
plt.hist(dfsfc,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfc),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfc),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(27)/27.0)
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


# In[243]:

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

# In[ ]:

plt.figure(figsize=(12,5))
plt.hist(dftoai,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoai),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoai),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(27)/27.0)
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


# In[ ]:

plt.figure(figsize=(12,5))
plt.hist(dfsfci,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfci),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfci),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(27)/27.0)
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


# # Redo the analysis for number 22134

# In[213]:

num2 = 22134


# In[214]:

s2 = sio.loadmat(fp+'MOC_1solx_DARE_{vv}_{num}.mat'.format(vv=vv,num=num2))


# In[215]:

dftoa2 = np.array([s2['solutions'][0,i]['dF_toa_24hr'][0,0][0,0] for i in xrange(len(s2['solutions'][0,:]))])
dfsfc2 = np.array([s2['solutions'][0,i]['dF_sfc_24hr'][0,0] for i in xrange(len(s2['solutions'][0,:]))])[:,0,0]


# In[216]:

toa_sel2 = {}
sfc_sel2 = {}
i_str = ['m1','s0','p1']
ii = 0
for ie in [-1,0,1]:
    for ia in [-1,0,1]:
        for im in [-1,0,1]:
            form = {'num':num,'e':i_str[ie+1],'a':i_str[ia+1],'s':i_str[im+1]}
            val = '{e}{a}{s}'.format(**form)
            toa_sel2[val] = s2['select'][0,0][ii]['dF_toa_24hr'][0][0][0,0]
            sfc_sel2[val] = s2['select'][0,0][ii]['dF_sfc_24hr'][0][0][0,0]
            ii += 1


# In[217]:

ssa2 = np.array([s2['solutions'][0,i]['ssa'][0,0][0,8] for i in xrange(len(s2['solutions'][0,:]))])
ext2 = np.array([s2['solutions'][0,i]['ext'][0,0][0,8] for i in xrange(len(s2['solutions'][0,:]))])
asy2 = np.array([s2['solutions'][0,i]['asy'][0,0][0,8] for i in xrange(len(s2['solutions'][0,:]))])


# In[220]:

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


# In[221]:

plt.figure(figsize=(12,5))
plt.hist(dftoa2,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoa2),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoa2),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(27)/27.0)
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


# In[222]:

plt.figure(figsize=(12,5))
plt.hist(dfsfc2,bins=40,alpha=1.0,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dfsfc2),color='k',lw=3,label='Mean')
plt.axvline(np.median(dfsfc2),color='grey',ls='--',lw=2,label='Median')
cs = cm.jet(np.arange(27)/27.0)
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


# In[223]:

len(ext2)


# # Print out all the resulting values

# In[231]:

print 'TOA {}'.format(num)
print 'Mean:',np.mean(dftoa)
print 'Median',np.median(dftoa)
print 'Mean+std:',np.mean(dftoa)+np.std(dftoa)
print 'Mean-std:',np.mean(dftoa)-np.std(dftoa)
for k in np.sort(toa_sel.keys()):
    print k,toa_sel[k]


# In[235]:

print 'TOA {}'.format(num)
print 'Mean:',np.mean(dfsfc)
print 'Median',np.median(dfsfc)
print 'Mean+std:',np.mean(dfsfc)+np.std(dfsfc)
print 'Mean-std:',np.mean(dfsfc)-np.std(dfsfc)
for k in np.sort(toa_sel.keys()):
    print k,toa_sel[k],sfc_sel[k]


# In[234]:

for k in np.sort(toa_sel.keys()):
    print sfc_sel[k]


# In[237]:

print 'TOA, SFC {}'.format(num2)
print 'Mean:',np.mean(dftoa2),np.mean(dfsfc2)
print 'Median',np.median(dftoa2),np.median(dfsfc2)
print 'Mean+std:',np.mean(dftoa2)+np.std(dftoa2),np.mean(dfsfc2)+np.std(dfsfc2)
print 'Mean-std:',np.mean(dftoa2)-np.std(dftoa2),np.mean(dfsfc2)-np.std(dfsfc2)
for k in np.sort(toa_sel2.keys()):
    print k,toa_sel2[k],sfc_sel2[k]


# In[238]:

for k in np.sort(toa_sel2.keys()):
    print toa_sel2[k]


# In[239]:

for k in np.sort(toa_sel2.keys()):
    print sfc_sel2[k]


# In[ ]:



