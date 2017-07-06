
# coding: utf-8

# # Intro
# Simple Program to load and check the 4STAR archive files.
# 
# For R2 of ORACLES, aod and Gases

# # Load the defaults and imports

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth


# In[2]:

from linfit import linfit


# In[3]:

get_ipython().magic(u'matplotlib notebook')


# In[4]:

fp ='C:/Users/sleblan2/Research/ORACLES/'


# # load the files

# In[5]:

days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[5]:

days = ['20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']#,'20160927','20160929','20160930']#,'20160825']


# In[54]:

days = ['20160902','20160904','20160906','20160908']


# In[5]:

days = ['20160831']


# In[6]:

vv = 'R1'


# In[7]:

vi = 'v7'


# In[38]:

outaod_RA = []
outaod_head_RA = []
outgas_RA = []
outgas_head_RA = []
for i,d in enumerate(days):
    #try:
    #    print 'Doing day: {}'.format(d)
    #    fname_aod = fp+'aod_ict/{vi}/4STAR-AOD_P3_{}_{vv}.ict'.format(d,vv=vv,vi=vi)
    #    tt,th = load_ict(fname_aod,return_header=True)
    #except:
    #    print '*** Problem with day: {} *** Skipping '.format(d)
    #    days.pop(i)
    #    continue
    
    #outaod_RA.append(tt)
    #outaod_head_RA.append(th)
    try:
        print 'Doing day: {}'.format(d)
        fname_gas = fp+'gas_ict/{vi}/4STAR-GAS_P3_{}_{vv}.ict'.format(d,vv=vv,vi=vi)
        print fname_gas
        tt,th = load_ict(fname_gas,return_header=True)
    except:
        print '*** Problem with day: {} *** Skipping '.format(d)
        days.pop(i)
        continue
    
    #outaod_RA.append(tt)
    #outaod_head_RA.append(th)
    #fname_gas = fp+'gas_ict/korusaq-4STAR-GASES_DC8_{}_RA.ict'.format(d)
    #ttr,thr = load_ict(fname_gas,return_header=True)
    outgas_RA.append(tt)
    outgas_head_RA.append(th)


# In[26]:

len(outgas_RA)


# In[27]:

len(days)


# ## Check the files for integrity and header info

# In[123]:

for i,s in enumerate(outaod_head_RA[0]):
    for ig,g in enumerate(outaod_head_RA):
        if not s==g[i]:
            print 'no match on {vv} aod string line {}: {} and {vv} of num {}:{} '.format(i,s,ig,g[i],vv=vv)
#    for ir,r in enumerate(outgas_head_RA):
#        if not s==r[i]:
#            print 'no match on RA gas string line {}: {} and RA of num {}:{} '.format(i,s,ir,r[i])


# In[28]:

print 'day:       AOD {vv}     GAS {vv}'.format(vv=vv)
for i,d in enumerate(days):
    try:
        print '{}: {}  {}'.format(d,len(outaod_RA[i]['Start_UTC']),len(outgas_RA[i]['Start_UTC']))
    except:
        print '{}: missed'.format(d)


# In[29]:

outgas_head_RA[0]


# ## Check the variables in header

# In[32]:

nm = outgas_RA[0].dtype.names


# In[33]:

nm


# # make plots of the gases

# In[34]:

outgas_RA[0].dtype.names


# In[39]:

for i,d in enumerate(days):
    fig,ax = plt.subplots(8,sharex=True,figsize=(10,8))
    ax = ax.ravel()
    ax[0].set_title('Gas {vv} saved file for flight {}'.format(d,vv=vv))
    
    ax[0].plot(outgas_RA[i]['Start_UTC'][outgas_RA[i]['QA_CWV']==0],outgas_RA[i]['CWV'][outgas_RA[i]['QA_CWV']==0],'.')
    ax[0].set_ylabel('CWV')
    #ax[0].set_ylim(0,3)
    ax[0].axhline(0,color='k')
    #box = ax[0].get_position()
    #ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #ax[0].legend(frameon=False,loc='center left',bbox_to_anchor=(1.1,-0.2),numpoints=1)
    ax[1].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['std_CWV'],'xr')
    ax[1].set_ylabel('std_CWV')
    
    ax[2].plot(outgas_RA[i]['Start_UTC'][outgas_RA[i]['QA_O3']==0],outgas_RA[i]['VCD_O3'][outgas_RA[i]['QA_O3']==0],'.')
    ax[2].set_ylabel('VCD_O3')
    
    ax[3].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['resid_O3'],'xr')
    ax[3].set_ylabel('resid_O3')
    
    ax[4].plot(outgas_RA[i]['Start_UTC'][outgas_RA[i]['QA_NO2']==0],outgas_RA[i]['VCD_NO2'][outgas_RA[i]['QA_NO2']==0],'.')
    ax[4].set_ylabel('VCD_NO2')
    
    ax[5].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['resid_NO2'],'xr')
    ax[5].set_ylabel('resid_NO2')
    
    ax[6].plot(outgas_RA[i]['Start_UTC'][outgas_RA[i]['QA_HCOH']==0],outgas_RA[i]['VCD_HCOH'][outgas_RA[i]['QA_HCOH']==0],'.')
    ax[6].set_ylabel('VCD_HCOH')
    
    ax[7].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['resid_HCOH'],'xr')
    ax[7].set_ylabel('resid_HCOH')
    
    ax[7].set_xlabel('UTC [h]')


# In[ ]:



