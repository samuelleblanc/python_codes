#!/usr/bin/env python
# coding: utf-8

# # Intro

# Software for reading the .dat files, raw data from 4STAR

# # imports and function

# In[1]:


import numpy as np
import os


# In[2]:


fp = 'C:\Users\Sunsat\Documents\Data\\20160208'


# In[5]:


f = fp+'\\20160208_004_NIR_SUN.dat'


# In[6]:


f


# In[8]:


os.path.isfile(f)


# In[15]:


d = np.genfromtxt(f,names=True,skip_header=9)


# In[39]:


head = []
with open(f,'r') as ff:
    for line in ff:
        if line[0].isdigit():
           break
        else:
            head.append(line)


# In[41]:


def get_header(f):
    head = []
    with open(f,'r') as ff:
        for line in ff:
            if line[0].isdigit():
               break
            else:
                head.append(line)
    numhead = len(head)-1
    return head,numhead


# In[42]:


g = fp+'\\20160208_009_NIR_FORJ.dat'


# In[43]:


h,n = get_header(g)


# In[44]:


h


# In[45]:


n


# In[46]:


h = fp+'\\20160208_008_TRACK.dat'


# In[47]:


hh,nh = get_header(h)


# In[48]:


nh


# In[49]:


hh


# In[67]:


h


# In[ ]:





# In[251]:


typ,fnum =startype(h)


# In[252]:


typ


# In[254]:


int(fnum)


# In[71]:


typ.lower()


# In[193]:


o = datetime.strptime('20:23:09.015','%H:%M:%S.%f')#,'%H:%M:%S.%s')


# In[195]:


utc = o.hour+o.minute/60.0+o.second/3600.0+o.microsecond/3600.0/1000000.0


# In[172]:


def mktime(txt):
    return datetime.strptime(txt,'%H:%M:%S.%f')


# In[197]:


conv1 = {'Time':mkktime}


# In[ ]:





# In[214]:


def mkktime(txt):
    o = datetime.strptime(txt,'%H:%M:%S.%f')
    return o.hour+o.minute/60.0+o.second/3600.0+o.microsecond/3600.0/1000000.0


# In[177]:


conv = {'Time':mktime,1:mktime,'time':mktime}


# In[167]:


convv = {1:mktime}


# In[190]:


nm = ['N','Time','AzAngle','AzError','ElAngle','ElError','Azpos','Elpos','Vb_t','Vl_r','Vtot','Tsetbox',      'Tmonbox','Tbox','Tvnir','Tnir','P1','P2','P3','P4','T1','T2','T3','T4']


# In[222]:


new_nm =  ['N','t','Az_deg','Az_corr','El_deg','El_corr','Az_step','El_step',           'QdV_TB','QdV_LR','QdV_tot','T_box_set','T_box_mon','T_box','T_spec_uvis',           'T_spec_nir','P1','P2','P3','P4','T1','T2','T3','T4']


# In[223]:


def convt(t):
    return float(t)*1000.0-273.15


# In[224]:


conv = {'t':mkktime,'T_box':lambda t:float(t)*100.0-273.15,        'T1':convt,'T2':convt,'T3':convt,'T4':convt}


# In[235]:


dd = np.genfromtxt(h,skip_header=2,names=new_nm,delimiter='\t',converters=conv)


# In[236]:


dd


# In[238]:


dd.dtype.names


# In[244]:


def to_dict(d):
    'Simple function to switch from numpy named array to dict with numpy arrays'
    dd = {}
    for n in d.dtype.names:
        dd[n] = d[n]
    return d


# In[243]:


ddd['El_temp_K'] = 1000.0*ddd['T1']


# In[226]:


import hdf5storage as hs


# In[229]:


import scipy.io as sio


# In[230]:


m = sio.loadmat('C:\Users\Sunsat\Documents\Data\\20151118star.mat')


# In[231]:


m.keys()


# In[233]:


m['track'].dtype.names


# In[280]:


m['track']['fname']


# In[268]:


m['track']['filen'][0][0][:]


# In[264]:


np.arange(len())


# In[261]:


np.zeros_like(m['track']['filen'])


# In[255]:


m['track']['filen']


# In[273]:


s = {}


# In[274]:


s['forj'] = [1,2,3]


# In[275]:


s.keys()


# In[277]:


if 'sun' in s:
    print 'yes'


# In[281]:


fparts


# In[283]:


fparts = os.path.split(h)


# In[284]:


fparts


# In[ ]:




