
# coding: utf-8

# In[8]:

import matplotlib.pyplot as plt


# In[9]:

get_ipython().magic(u'matplotlib notebook')


# In[10]:

fp = 'C:\Users\NASA_\Documents\Research\ORACLES\plots\\'


# In[11]:

def init_packets(name):
    """
    Function to load the inital set of packets from the data API
    """
    import numpy as np
    from datetime import datetime
    from urllib2 import urlopen
    
    toutc = lambda x: x.hour+x.minute/60.0+x.second/3600.0
    
    #Get the names of parameters to init dict
    site = 'http://10.15.1.1/API/parameter_data/N426NA/{}/XML'.format(name)
    r = urlopen(site)
    rr = r.read()
    rn = [s[s.rfind('<point')+len('<point')+1:s.find('/>')] for s in rr.splitlines() if s.find('point')>-1]
    try:
        dict_names = [s.split('="')[0] for s in rn[0].split('" ')]
        dict_vals = [s.split('="')[1] for s in rn[0].split('" ')]
    except IndexError:
        print 'Problem! returning string'
        return rn
        
    d = {}
    for i,n in enumerate(dict_names):
        try:
            d[n] = np.array([datetime.strptime(dict_vals[i],'%Y-%m-%d %H:%M:%S.%f')])
            if n=='timestamp':
                d['UTC'] = np.array([toutc(tt) for tt in d[n]])
        except:
            try:
                d[n] = np.array([datetime.strptime(dict_vals[i],'%Y-%m-%d %H:%M:%S')])
                if n=='timestamp':
                    d['UTC'] = np.array([toutc(tt) for tt in d[n]])
            except:
                try:
                    d[n] = np.array([float(dict_vals[i])])
                except:
                    d[n] = np.array([dict_vals[i]])
    return d


# In[12]:

def get_packets(d,start_time=None,end_time=None,verbose=True):
    """
    Get packets from the start_time (defaults to 0) in seconds from 1970
    appends the packets to the correct values in the dict, and returns the dict
    """
    import numpy as np
    from datetime import datetime
    from urllib2 import urlopen
    
    # get the packets from web using the XML interface
    name = d['header'][0]
    last = d['timestamp'][-1]
    if len(d['header'])>1:
        epoch_time = datetime(1970,1,1)
        start_time = (last.toordinal()-epoch_time.toordinal())*24*60*60+last.hour*3600+last.minute*60+last.second
    
    site = 'http://10.15.1.1/API/parameter_data/N426NA/{}/XML?Start={}'.format(name,start_time)
    if end_time:
        site = site+'&End={}'.format(end_time)
    r = urlopen(site)
    rr = r.read()
    rn = [s[s.rfind('<point')+len('<point')+1:s.find('/>')] for s in rr.splitlines() if s.find('point')>-1]
    nn = float(len(rn))
    if verbose: print 'number of {} points found: {}'.format(name,nn)
    toutc = lambda x: x.hour+x.minute/60.0+x.second/3600.0
    for j,rnn in enumerate(rn):
        dict_names = [s.split('="')[0] for s in rnn.split('" ')]
        dict_vals = [s.split('="')[1] for s in rnn.split('" ')]
        print u'{}%\r'.format(j/nn*100.0),
        for i,n in enumerate(dict_names):
            #print j,i,n
            try:
                du = np.array([datetime.strptime(dict_vals[i],'%Y-%m-%d %H:%M:%S.%f')])
                d[n] = np.append(d[n],du)
                if n=='timestamp':
                    d['UTC'] = np.append(d['UTC'],np.array([toutc(tt) for tt in du]))
            except:
                try:
                    du = np.array([datetime.strptime(dict_vals[i],'%Y-%m-%d %H:%M:%S')])
                    d[n] = np.append(d[n],du)
                    if n=='timestamp':
                        d['UTC'] = np.append(d['UTC'],np.array([toutc(tt) for tt in du]))
                except:
                    try:
                        d[n] = np.append(d[n],np.array([float(dict_vals[i])]))
                    except:
                        d[n] = np.append(d[n],np.array([dict_vals[i]]))
    return d


# In[13]:

def startup_packet(name):
    'function to run at startup for slowly loading up the data structures, to prohibit overwhelming the data systems'
    from datetime import datetime
    
    d = init_packets(name)
    
    # Try and loop a few times to get the right start time and load slowly
    t = datetime.now()
    epoch_time = datetime(1970,1,1)
    start_time = (t.toordinal()-epoch_time.toordinal())*24*60*60
    now_time = (t.toordinal()-epoch_time.toordinal())*24*60*60+t.hour*3600+t.minute*60+t.second
    for i in xrange(48):
        end_time = start_time+3600/2
        d = get_packets(d,start_time=start_time,end_time=end_time)
        if start_time>now_time:
            break
        else:
            start_time = end_time
    return d


# In[15]:

iwg = startup_packet('IWG1')


# In[16]:

iwg = coma


# In[80]:

d = get_packets(d)


# In[82]:

higear = startup_packet('HIGEAR')


# In[85]:

higear.keys()


# In[17]:

iwg.keys()


# In[86]:

plt.figure()
plt.plot(iwg['r'],higear['cn'],'.')


# In[83]:

iwg = startup_packet('IWG1')


# In[64]:

w = startup_packet('WISPER')


# In[65]:

w.keys()


# In[67]:

w = get_packets(w)


# In[70]:

plt.figure()
plt.plot(w['pic0_qh2o_cal'],w['air_pres'],'.')


# In[75]:

plt.figure()
cb = plt.scatter(w['pic0_qh2o_cal'],w['air_pres'],c=w['UTC'],edgecolor='None')
plt.ylim(1100,500)
b = plt.colorbar(cb)
b.set_label('UTC [h]')


# In[76]:

w = get_packets(w)


# In[ ]:




# In[79]:

plt.figure()
ps = plt.scatter(w['pic0_delo_cal'][w['pic0_qh2o_cal']>200.0],w['pic0_deld_cal'][w['pic0_qh2o_cal']>200.0],
                 c=w['air_pres'][w['pic0_qh2o_cal']>200.0],edgecolor='None')
b = plt.colorbar(ps)
b.set_label('Air pressure')


# In[81]:

plt.figure()
plt.plot(d['UTC'],d['aod_500_nm'],'.')


# In[84]:

plt.figure()
plt.plot(iwg['latitude'],iwg['gps_msl_alt'],'.')
plt.xlabel('Latitude[$^\\circ$]')
plt.ylabel('GPS_MSL_ALT [m]')
plt.grid()
plt.title('2017-08-17 Altitude Profiles PRF04Y17')
plt.savefig(fp+'PRF04_alt_lat.png',dpi=600,transparent=True)


# In[ ]:

instruments = ['IWG1','4STAR','HIGEAR','COMA','UND','OZONE','SSFR','HSRL','WISPER','SP2']


# In[ ]:

data = {}
for na in instruments:
    try:
        data[na] = init_packets(na)
        data[na] = get_packets(na,start_time=)

