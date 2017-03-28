
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


# In[35]:

dfsfc = np.array([s['solutions'][0,i]['dF_sfc_24hr'][0,0] for i in xrange(len(s['solutions'][0,:]))])[:,0,0]


# In[36]:

dftoa.shape


# In[80]:

s['select'][0,0][1]['dF_sfc_24hr'][0][0][0,0]


# In[81]:

toa_sel = {}
sfc_sel = {}
i_str = ['m1','s0','p1']
ii = 0
for ie in [-1,0,1]:
    for ia in [-1,0,1]:
        for im in [-1,0,1]:
            form = {'num':num,'e':i_str[ie+1],'a':i_str[ia+1],'s':i_str[im+1]}
            val = '{e}{a}{s}'.format(**form)
            toa_sel[val] = s['select'][0,0][ii]['dF_toa_24hr'][0][0][0,0]
            sfc_sel[val] = s['select'][0,0][ii]['dF_sfc_24hr'][0][0][0,0]
            ii += 1


# In[83]:

sfc_sel.keys()


# # Start to plot out results

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


# In[127]:

plt.figure(figsize=(12,5))
plt.hist(dftoa,bins=40,alpha=0.7,edgecolor='None',label='all retrievals')
plt.axvline(np.mean(dftoa),color='k',lw=3,label='Mean')
plt.axvline(np.median(dftoa),color='k',lw=2,label='Median')
cs = cm.jet(np.arange(27)/27.0)
for i,k in enumerate(toa_sel.keys()):
    plt.axvline(toa_sel[k],lw=2,label=k,color=tuple(cs[i,:]))

plt.xlabel('DARE TOA [W/m$^2$]')
plt.ylabel('Number of solutions')
plt.title('DARE TOA for all solutions on single pixel #{}'.format(num))

box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])

plt.legend(bbox_to_anchor=(1.03,1.05),loc=2,ncol=2)


# In[ ]:



