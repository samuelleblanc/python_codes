
# coding: utf-8

# In[7]:

from urllib2 import urlopen
import pandas as pd
import numpy as np
from StringIO import StringIO


# In[21]:

import matplotlib.pyplot as plt


# In[86]:

site = 'http://10.15.1.1/API/parameter_data/N426NA/4STAR/XML?Start=1502196800&End=1502196805'


# In[87]:

r = urlopen(site)


# In[88]:

rr = r.read()


# In[89]:

rr.find(';')


# In[91]:

rr.splitlines()


# In[96]:




# In[100]:

rn[0]


# In[108]:

[s.split('="')[0] for s in rn.split('" ')]
u = np.array([s.split('="')[1].strip('"') for s in rn[0].split('" ')])



# In[58]:

[s.split('="')[0] for s in rn.split('" ')]
[s.split('=')[0] for s in rr.split() if s.find('=')>-1]
[s.split('=')[1] for s in rr.split() if s.find('=')>-1]


# In[115]:

from datetime import datetime


# In[ ]:

help()


# In[156]:

datetime.strptime(u[16],'%Y-%m-%d %H:%M:%S.%f')


# In[142]:

help(datetime)


# In[119]:

u[16]


# In[117]:

datetime.(u[16])


# In[169]:

def init_packets(name,flight=True):
    """
    Function to load the inital set of packets from the data API
    """
    import numpy as np
    from datetime import datetime
    from urllib2 import urlopen
    
    #Get the names of parameters to init dict
    if flight:
        site = 'http://10.15.1.1/API/parameter_data/N426NA/{}/XML'.format(name)
    else:
        site = 'http://10.15.1.1/API/parameter_data/N426NA/{}/XML'.format(name)
    r = urlopen(site)
    rr = r.read()
    rn = [s[s.rfind('point')+len('point')+1:s.find('/>')] for s in rr.splitlines() if s.find('point')>-1]
    dict_names = [s.split('="')[0] for s in rn[0].split('" ')]
    dict_vals = [s.split('="')[1] for s in rn[0].split('" ')]
    d = {}
    for i,n in enumerate(dict_names):
        try:
            d[n] = np.array([datetime.strptime(dict_vals[i],'%Y-%m-%d %H:%M:%S.%f')])
        except:
            try:
                d[n] = np.array([float(dict_vals[i])])
            except:
                d[n] = np.array([dict_vals[i]])
    return d


# In[ ]:

def get_packets(d,start_time=0):
    """
    Get packets from the start_time (defaults to 0) in seconds from 1970
    appends the packets to the correct values in the dict, and returns the dict
    """
    import numpy as np
    from datetime import datetime
    from urllib2 import urlopen
    
    # get the packets from web using the XML interface
    name = d['head']
    site = 'http://10.15.1.1/API/parameter_data/N426NA/{}/XML?Start={}'.format(name,start_time)
    r = urlopen(site)
    rr = r.read()
    rn = [s[s.rfind('point')+len('point')+1:s.find('/>')] for s in rr.splitlines() if s.find('point')>-1]
    dict_names = [s.split('="')[0] for s in rn[0].split('" ')]
    dict_vals = [s.split('="')[1] for s in rn[0].split('" ')]
    d = {}
    for i,n in enumerate(dict_names):
        try:
            d[n] = np.array([datetime.strptime(dict_vals[i],'%Y-%m-%d %H:%M:%S.%f')])
        except:
            try:
                d[n] = np.array([float(dict_vals[i])])
            except:
                d[n] = np.array([dict_vals[i]])
    return d
    


# In[170]:

d = setup_packets('4STAR')


# In[171]:

d.keys()


# In[173]:

d['timestamp']


# In[ ]:



