
# coding: utf-8

# Simple Program to load and check the 4STAR archive files.
# 
# For R2 of SEAC4RS

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

fp ='C:/Users/sleblan2/Research/SEAC4RS/'


# In[4]:

days = ['20130806','20130808','20130812','20130814','20130816','20130819',
        '20130821','20130823','20130826','20130827','20130830','20130902',
        '20130904','20130906','20130909','20130911','20130913','20130916',
        '20130918','20130921','20130923']


# In[5]:

help(load_ict)


# In[6]:

out_R2 = []
out_head_R2 = []
out_R1 = []
out_head_R1 = []
for d in days:
    fname = fp+'ict/R2_201601/SEAC4RS-4STAR-AOD-CWV_DC8_{}_R2_YoheiQuick.ict'.format(d)
    tt,th = load_ict(fname,return_header=True)
    out_R2.append(tt)
    out_head_R2.append(th)
    
    fname_R1 = fp+'ict/SEAC4RS-4STAR-AOD-CWV_DC8_{}_R1.ict'.format(d)
    ttr,thr = load_ict(fname_R1,return_header=True)
    out_R1.append(ttr)
    out_head_R1.append(thr)


# In[15]:

for i,s in enumerate(out_head_R2[0]):
    for ig,g in enumerate(out_head_R2):
        if not s==g[i]:
            print 'no match on R2 string line {}: {} and R2 of num {}:{} '.format(i,s,ig,g[i])
    for ir,r in enumerate(out_head_R1):
        if not s==r[i]:
            print 'no match on R1 string line {}: {} and R2 of num {}:{} '.format(i,s,ir,r[i])


# In[26]:

print 'day:       R1     R2'
for i,d in enumerate(days):
    try:
        print '{}: {}  {}'.format(d,len(out_R1[i]['Start_UTC']),len(out_R2[i]['Start_UTC']))
    except:
        print '{}: missed'.format(d)


# In[51]:

out_head_R2[0]


# In[49]:

nm = out_R2[0].dtype.names


# In[50]:

nm


# In[63]:

nm[9:26]


# In[83]:

np.where(out_R2[i][nm[4]]==1)[0]


# In[94]:

plt.figure()
plt.plot(out_R2[0][nm[0]],out_R2[0][nm[9]],'.')
for x in out_R2[0][nm[0]][np.where(out_R2[0][nm[4]]==1)[0]]:
    plt.axvline(x,color='#DDDDDD',alpha=0.02)


# In[104]:

import cmaps


# In[105]:

cmaps.cmaps()


# In[112]:

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, 26-9)])
for i in np.linspace(0,1,26-9):
    plt.plot([i,i*2,i*3,1,2,3,4])


# In[113]:

for i,d in enumerate(days):
    fig,ax = plt.subplots(3,sharex=True,figsize=(12,5))
    ax = ax.ravel()
    ax[0].set_title(d)
    ax[0].set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, 26-9)])
    for aod in nm[9:26]:
        ax[0].plot(out_R2[i][nm[0]],out_R2[i][aod],'.',label=aod)
    try:
        for x in out_R2[i][nm[0]][np.where(out_R2[i][nm[4]]==1)[0]]:
            ax[0].axvline(x,color='#DDDDDD',alpha=0.02)
    except:
        pass
    ax[0].set_ylabel('AOD')
    ax[0].axhline(0,color='k')
    box = ax[0].get_position()
    ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[0].legend(frameon=False,loc='center left',bbox_to_anchor=(1,-0.3),numpoints=1)
    ax[1].set_color_cycle([plt.cm.gist_ncar(k) for k in np.linspace(0, 1, 26-9)])
    for aod in nm[27:]:
        ax[1].plot(out_R2[i][nm[0]],out_R2[i][aod],'.',label=aod)
    ax[1].set_ylabel('AOD unc.')
    ax[1].set_ylim([0,0.2])
    box = ax[1].get_position()
    ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[2].plot(out_R2[i][nm[0]],out_R2[i]['GPS_alt'],'.')
    ax[2].set_ylabel('Alt [m]')
    box = ax[2].get_position()
    ax[2].set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax[2].set_xlabel('UTC [h]')
    plt.savefig(fp+'ict/R2_{}.png'.format(d),dpi=600,transparent=True)


# In[ ]:



