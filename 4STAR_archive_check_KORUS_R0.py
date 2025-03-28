
# coding: utf-8

# # Intro
# Simple Program to load and check the 4STAR archive files.
# 
# For RA of KORUS-AQ, aod and Gases

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


get_ipython().magic(u'matplotlib notebook')


# In[3]:


fp ='C:/Users/sleblan2/Research/KORUS-AQ/'


# In[4]:


vr='R0'


# # load the files

# In[5]:


days = ['20160501','20160503','20160504','20160506','20160510','20160511',
        '20160512','20160516','20160517','20160519','20160521','20160524',
        '20160526','20160529','20160601','20160602','20160604','20160608',
        '20160609','20160614']#,'20160617','20160618']


# In[23]:


days = ['20160501','20160506','20160512','20160517','20160526','20160614']


# In[6]:


days = ['20160516']


# In[6]:


outaod_RA = []
outaod_head_RA = []
outgas_RA = []
outgas_head_RA = []
for d in days:
    fname_aod = fp+'aod_ict/{vr}/korusaq-4STAR-AOD_DC8_{d}_{vr}.ict'.format(vr=vr,d=d)
    tt,th = load_ict(fname_aod,return_header=True)
    outaod_RA.append(tt)
    outaod_head_RA.append(th)
    
    #fname_gas = fp+'gas_ict/korusaq-4STAR-GASES_DC8_{}_RA.ict'.format(d)
    #ttr,thr = load_ict(fname_gas,return_header=True)
    #outgas_RA.append(ttr)
    #outgas_head_RA.append(thr)


# In[6]:


for i,s in enumerate(outaod_head_RA[0]):
    for ig,g in enumerate(outaod_head_RA):
        if not s==g[i]:
            print 'no match on RA aod string line {}: {} and RA of num {}:{} '.format(i,s,ig,g[i])
    for ir,r in enumerate(outgas_head_RA):
        if not s==r[i]:
            print 'no match on RA gas string line {}: {} and RA of num {}:{} '.format(i,s,ir,r[i])


# In[8]:


print 'day:       AOD RA     GAS RA'
for i,d in enumerate(days):
    try:
        print '{}: {}  {}'.format(d,len(outaod_RA[i]['Start_UTC']),len(outgas_RA[i]['Start_UTC']))
    except:
        print '{}: missed'.format(d)


# In[7]:


outaod_head_RA[0]


# In[9]:


outgas_head_RA[0]


# In[8]:


nm = outaod_RA[0].dtype.names


# In[9]:


nm


# In[10]:


wl = nm[6:-1]


# In[11]:


plt.figure()
plt.plot(out_R2[0][nm[0]],out_R2[0][nm[9]],'.')
for x in out_R2[0][nm[0]][np.where(out_R2[0][nm[4]]==1)[0]]:
    plt.axvline(x,color='#DDDDDD',alpha=0.02)


# In[104]:


import cmaps


# In[105]:


cmaps.cmaps()


# In[12]:


for a in wl:
    print a


# In[13]:


wl = wl[0:16]


# In[14]:


for i,d in enumerate(days):
    fig,ax = plt.subplots(2,sharex=True,figsize=(9,5))
    ax = ax.ravel()
    ax[0].set_title('AOD {vr} saved file for flight {}'.format(d,vr=vr))
    ax[0].set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, len(wl))])
    for aod in wl:
        ax[0].plot(outaod_RA[i][nm[0]],outaod_RA[i][aod],'.',label=aod)
    try:
        for x in outaod_RA[i][nm[0]][np.where(outaod_RA[i][nm[4]]==1)[0]]:
            ax[0].axvline(x,color='#DDDDDD',alpha=0.02)
    except:
        pass
    ax[0].set_ylabel('AOD')
    ax[0].set_ylim(0,2)
    ax[0].axhline(0,color='k')
    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[0].legend(frameon=False,loc='center left',bbox_to_anchor=(1.1,-0.2),numpoints=1)
    ax[1].plot(outaod_RA[i][nm[0]],outaod_RA[i]['GPS_Alt'],'.')
    ax[1].set_ylabel('Alt [m]')
    axy = ax[1].twinx()
    axy.plot(outaod_RA[i][nm[0]],outaod_RA[i]['amass_aer'],'.g')
    axy.set_ylabel('Airmass factor',color='g')
    box = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    box = axy.get_position()
    axy.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[1].set_xlabel('UTC [h]')
    plt.savefig(fp+'aod_ict/{vr}_{}.png'.format(d,vr=vr),dpi=600,transparent=True)


# In[14]:


nm[4]


# In[15]:


fl = np.where(outaod_RA[i][nm[4]]==1)[0]


# In[16]:


fl.shape


# In[18]:


for i,d in enumerate(days):
    if i>0: break
    fig,ax = plt.subplots(2,sharex=True,figsize=(9,5))
    ax = ax.ravel()
    fl = np.where(outaod_RA[i][nm[4]]==1)[0]
    ax[0].set_title('AOD {vr} saved file for flight {}'.format(d,vr=vr))
    ax[0].set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, len(wl))])
    for aod in wl:
        ax[0].plot(outaod_RA[i][nm[0]][fl],outaod_RA[i][aod][fl],'.',label=aod)
    #try:
    #    for x in outaod_RA[i][nm[0]][np.where(outaod_RA[i][nm[4]]==1)[0]]:
    #        ax[0].axvline(x,color='#DDDDDD',alpha=0.02)
    #except:
    #    pass
    ax[0].set_ylabel('AOD')
    ax[0].set_ylim(0,2)
    ax[0].axhline(0,color='k')
    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[0].legend(frameon=False,loc='center left',bbox_to_anchor=(1.1,-0.2),numpoints=1)
    for aod in nm[6+len(wl):-1]:
        ax[1].plot(outaod_RA[i][nm[0]][fl],outaod_RA[i][aod][fl],'.',label=aod)
    ax[1].set_ylabel('AOD Uncertainty')
    #ax[1].plot(outaod_RA[i][nm[0]],outaod_RA[i]['GPS_Alt'],'.')
    #ax[1].set_ylabel('Alt [m]')
    #axy = ax[1].twinx()
    #axy.plot(outaod_RA[i][nm[0]],outaod_RA[i]['amass_aer'],'.g')
    #axy.set_ylabel('Airmass factor',color='g')
    box = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #box = axy.get_position()
    #axy.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[1].set_xlabel('UTC [h]')
    plt.savefig(fp+'aod_ict/{vr}_{}_unc.png'.format(d,vr=vr),dpi=600,transparent=True)


# In[ ]:


for i,d in enumerate(days):
    fig,ax = plt.subplots(2,sharex=True,figsize=(9,5))
    ax = ax.ravel()
    ax[0].set_title('Gas {vr} saved file for flight {}'.format(d,vr=vr))
    ax[0].set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, len(wl))])
    for aod in wl:
        ax[0].plot(outgas_RA[i][nm[0]],outgas_RA[i][aod],'.',label=aod)
    try:
        for x in outgas_RA[i][nm[0]][np.where(outgas_RA[i][nm[4]]==1)[0]]:
            ax[0].axvline(x,color='#DDDDDD',alpha=0.02)
    except:
        pass
    ax[0].set_ylabel('AOD')
    ax[0].set_ylim(0,3)
    ax[0].axhline(0,color='k')
    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[0].legend(frameon=False,loc='center left',bbox_to_anchor=(1.1,-0.2),numpoints=1)
    ax[1].plot(outaod_RA[i][nm[0]],outaod_RA[i]['GPS_Alt'],'.')
    ax[1].set_ylabel('Alt [m]')
    axy = ax[1].twinx()
    axy.plot(outaod_RA[i][nm[0]],outaod_RA[i]['amass_aer'],'.g')
    axy.set_ylabel('Airmass factor',color='g')
    box = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    box = axy.get_position()
    axy.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[1].set_xlabel('UTC [h]')
    plt.savefig(fp+'aod_ict/{vr}_{}.png'.format(d,vr=vr),dpi=600,transparent=True)


# In[17]:


outgas_RA[0].dtype.names


# In[36]:


for i,d in enumerate(days):
    fig,ax = plt.subplots(7,sharex=True,figsize=(9,5))
    ax = ax.ravel()
    ax[0].set_title('Gas RA saved file for flight {}'.format(d))
    
    ax[0].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['CWV'])
    ax[0].set_ylabel('CWV')
    #ax[0].set_ylim(0,3)
    ax[0].axhline(0,color='k')
    #box = ax[0].get_position()
    #ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #ax[0].legend(frameon=False,loc='center left',bbox_to_anchor=(1.1,-0.2),numpoints=1)
    ax[1].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['QA_CWV'],'x')
    ax[1].set_ylabel('QA_CWV')
    
    ax[2].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['VCD_O3'])
    ax[2].set_ylabel('VCD_O3')
    
    ax[3].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['QA_O3'],'x')
    ax[3].set_ylabel('QA_O3')
    
    ax[4].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['VCD_NO2'])
    ax[4].set_ylabel('VCD_NO2')
    
    ax[5].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['QA_NO2'],'x')
    ax[5].set_ylabel('QA_NO2')
    
    ax[6].plot(outgas_RA[i]['Start_UTC'],outgas_RA[i]['qual_flag'],'x')
    ax[6].set_ylabel('qual_flag')
    
    ax[6].set_xlabel('UTC [h]')


# # Combineall the data into a single array

# In[18]:


ar = {}
for n in nm:
    ar[n] = np.array([])


# In[19]:


ar['days'] = np.array([])


# In[20]:


for i,d in enumerate(days):
    ar['days'] = np.append(ar['days'],np.zeros_like(outaod_RA[i]['Start_UTC'])+i)
    for n in nm:
        ar[n] = np.append(ar[n],outaod_RA[i][n])


# In[21]:


ar['GPS_Alt'].shape


# ## Save the combined data

# In[31]:


import hdf5storage as hs


# In[33]:


hs.savemat(fp+'/aod_ict/all_aod_KORUS_ict.mat',ar)


# ## Filterout bad data

# In[22]:


ar['fl_QA'] = ar['qual_flag']==0


# In[23]:


ar['fl'] = ar['fl_QA']


# In[24]:


ar['fl_alt'] = ar['GPS_Alt']<=1500


# In[25]:


ar['fl_alt1'] = ar['GPS_Alt']<=600


# In[26]:


ar['fl1'] = ar['fl_QA']&ar['fl_alt']


# In[27]:


ar['fl2'] = ar['fl_QA']&ar['fl_alt1']


# ## Plot out all the data

# In[28]:


from plotting_utils import prelim


# In[29]:


plt.figure()
plt.hist(ar['AOD0501'][ar['fl']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='g',label='all')
plt.hist(ar['AOD0501'][ar['fl1']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='b',label='below 1500 m')
plt.hist(ar['AOD0501'][ar['fl2']],bins=30,range=(0,1.0),alpha=0.5,normed=False,edgecolor='None',color='r',label='below 600 m')
#plt.yscale('log')
plt.xlabel('AOD at 501 nm')
plt.ylabel('Counts')
plt.grid()
plt.title('AOD distribution')
prelim()
plt.legend(frameon=False)
plt.savefig(fp+'aod_ict/{vv}_AOD_histogram.png'.format(vv=vr),dpi=600,transparent=True)


# In[29]:


np.nanmean(ar['AOD0501'][ar['fl2']])


# In[30]:


np.nanmedian(ar['AOD0501'][ar['fl2']])

