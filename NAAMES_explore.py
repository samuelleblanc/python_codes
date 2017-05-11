
# coding: utf-8

# # Intro
# Purpose:  
# 
#     Go through some of the extra data pieces to make a cohesive story
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures and save files...
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - plotting_utils (user defined plotting routines)
#     - hdf5storage
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - 4STAR_cloud retrieval .mat files
#   
#  Modification History:
#  
#      Written: by Samuel LeBlanc, NASA Ames, Moffett Field, CA, 2017-04-25

# # Import of modules

# In[1]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
get_ipython().magic(u'matplotlib notebook')
import numpy as np
import scipy.io as sio
import hdf5storage as hs
import Sp_parameters as Sp


# In[49]:

import load_utils as lu
from write_utils import nearest_neighbor
from load_utils import recarray_to_dict
import plotting_utils as pu


# In[3]:

# set the basic directory path
fp = 'C:/Users/sleblan2/Research/NAAMES/'
fp_plot = 'C:/Users/sleblan2/Research/NAAMES/plot/'


# In[4]:

vv = 'v1'


# # Load some files

# In[5]:

days = ['20151112','20151114','20151117','20151118','20151123']


# ## Load the LARGE AMS files

# In[6]:

ff = fp+'data_other/LARGE/'
rr = 'R0'


# In[7]:

ams = []
ams_header = []
for d in days:
    a,ah = lu.load_ict(ff+'NAAMES-LARGE-AMS_C130_{d}_{r}.ict'.format(d=d,r=rr),return_header=True)
    ams.append(a)
    ams_header.append(ah)


# In[8]:

ah


# ## Load the LARGE CIP (cloud precipitation particles)

# In[52]:

ff = fp+'data_other/LARGE/'
rr = 'R0'


# In[53]:

cip = []
cip_header = []
for d in days:
    ci,cih = lu.load_ict(ff+'NAAMES-LARGE-CIP_C130_{d}_{r}.ict'.format(d=d,r=rr),return_header=True)
    cip.append(ci)
    cip_header.append(cih)


# In[54]:

cih


# ## Load the CCN from LARGE

# In[77]:

ff = fp+'data_other/LARGE/'
rr = 'R0'


# In[78]:

ccn = []
ccn_header = []
for d in days:
    cn,cnh = lu.load_ict(ff+'NAAMES-LARGE-CCN_C130_{d}_{r}.ict'.format(d=d,r=rr),return_header=True)
    ccn.append(cn)
    ccn_header.append(cnh)


# In[79]:

cnh


# ## Load the Winds from LARGE

# In[59]:

ff = fp+'data_other/LARGE/'
rr = 'R0'


# In[61]:

wind = []
wind_header = []
for d in days:
    w,wh = lu.load_ict(ff+'NAAMES-LARGE-WINDS_C130_{d}_{r}.ict'.format(d=d,r=rr),return_header=True)
    wind.append(w)
    wind_header.append(wh)


# In[62]:

wh


# ## Read in the RSP cloud retrievals

# In[9]:

fr = fp+'data_other/RSP/'
rr = 'R1'


# In[10]:

rsp = []
rsp_header = []
for d in days:
    r,rh = lu.load_ict(fr+'NAAMES-RSP-WTRCLD_C130_{d}_{r}.ict'.format(d=d,r=rr),return_header=True)
    rsp.append(r)
    rsp_header.append(rh)


# In[11]:

rh


# ## Read in 4STAR zenith retrievals

# In[12]:

fs = fp+'zen_ict/'
rs = 'R0'


# In[13]:

star = []
star_header = []
for d in days:
    try:
        s,sh = lu.load_ict(fs+'4STAR-CLD_C130_{d}_{r}.ict'.format(d=d,r=rs),return_header=True)
    except:
        s,sh = {'bad':True},['bad=True']
    star.append(s)
    star_header.append(sh)


# In[14]:

sh


# ## Load the C130 Nav data

# In[15]:

fh = fp+'data_other/housekeeping/'
rh = 'R0'


# In[16]:

hsk = []
hsk_header = []
for d in days:
    try:
        h,hh = lu.load_ict(fh+'NAAMES-Hskping_c130_{d}_{r}.ict'.format(d=d,r=rh),return_header=True)
    except:
        h,hh = {'bad':True},['bad=True']
    hsk.append(h)
    hsk_header.append(hh)


# In[17]:

hh


# ## Load the cloud particle insitu data

# In[18]:

fc = fp+'data_other/LARGE/'
rc = 'R0'


# In[19]:

cdp = []
cdp_header = []
for d in days:
    try:
        c,ch = lu.load_ict(fc+'NAAMES-LARGE-CDP_C130_{d}_{r}.ict'.format(d=d,r=rh),return_header=True)
    except:
        c,ch = {'bad':True},['bad=True']
    cdp.append(c)
    cdp_header.append(ch)


# In[20]:

cdp_binsc = [3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,
            34.0,36.0,38.0,40.0,42.0,44.0,46.0,48.0,50.0]


# In[21]:

ch


# In[22]:

cdp_bins = []
for i,d in enumerate(days):
    bins = []
    for j in xrange(30):
        bins.append(cdp[i]['CDP_Bin{:02.0f}'.format(j)])
    cdp_bins.append(np.array(bins))


# In[23]:

cdp_bins[0].shape


# In[24]:

cdp_ref = []
for i,d in enumerate(days):
    


# ## Load cloud flags

# In[25]:

fc = fp+'data_other/LARGE/'
rc = 'R0'


# In[28]:

cld = []
cld_header = []
for d in days:
    try:
        cd,cdh = lu.load_ict(fc+'NAAMES-LARGE-CloudFlag_C130_{d}_{r}.ict'.format(d=d,r=rh),return_header=True)
    except:
        cd,cdh = {'bad':True},['bad=True']
    cld.append(recarray_to_dict(cd))
    cld_header.append(cdh)


# # Now plot out some data

# ## Prepare the data in an appropriate format

# In[89]:

from write_utils import nearest_neighbor
from load_utils import recarray_to_dict


# In[29]:

days


# In[30]:

a.dtype.names


# In[31]:

for i,d in enumerate(days):
    try: 
        print i,len(ams[i]['UTC_mid']),len(hsk[i]['Start_UTC']),len(star[i]['Start_UTC']),len(cld[i]['UTC_mid'])
    except:
        print i,len(ams[i]['UTC_mid']),len(hsk[i]['Start_UTC'])


# In[46]:

for i,d in enumerate(days):
    if i>1:
        star[i] = recarray_to_dict(star[i])


# In[47]:

for i,d in enumerate(days):
    if i>1:
        star[i]['alt'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['GPS_Altitude'],star[i]['Start_UTC'])


# In[33]:

for i,d in enumerate(days):
    cld[i]['alt'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['GPS_Altitude'],cld[i]['UTC_mid'],dist=1.0/3600.0)
    cld[i]['lat'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['Latitude'],cld[i]['UTC_mid'],dist=1.0/3600.0)
    cld[i]['lon'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['Longitude'],cld[i]['UTC_mid'],dist=1.0/3600.0)


# In[34]:

for i,d in enumerate(days):
    cdp[i] = recarray_to_dict(cdp[i])


# In[35]:

for i,d in enumerate(days):
    rsp[i] = recarray_to_dict(rsp[i])


# In[36]:

for i,d in enumerate(days):
    rsp[i]['alt'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['GPS_Altitude'],rsp[i]['UTC_start']/3600.0,dist=1.0/3600.0)


# In[37]:

for i,d in enumerate(days):
    cdp[i]['alt'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['GPS_Altitude'],cdp[i]['UTC_mid'],dist=1.0/3600.0)
    cdp[i]['lat'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['Latitude'],cdp[i]['UTC_mid'],dist=1.0/3600.0)
    cdp[i]['lon'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['Longitude'],cdp[i]['UTC_mid'],dist=1.0/3600.0)


# In[38]:

for i,d in enumerate(days):
    ams[i] = recarray_to_dict(ams[i])


# In[39]:

for i,d in enumerate(days):
    ams[i]['alt'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['GPS_Altitude'],cdp[i]['UTC_mid'],dist=1.0/3600.0)


# In[63]:

for i,d in enumerate(days):
    wind[i] = recarray_to_dict(wind[i])


# In[82]:

for i,d in enumerate(days):
    try:
        ccn[i] = recarray_to_dict(ccn[i])
    except:
        pass
    ccn[i]['alt'] = nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['GPS_Altitude'],ccn[i]['UTC_mid'],dist=1.0/3600.0)


# ## Check out AMS data

# In[185]:

plt.figure()
cs = ['r','g','b','k','c','y']
for i,d in enumerate(days):
    #plt.plot(ams[i]['UTC_mid'],ams[i]['Org_STP'],'x',color=cs[i],label=d,alpha=0.3)
    plt.plot(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],16),'-',color=cs[i],alpha=0.6,label=d)
plt.legend(frameon=False)
plt.xlabel('UTC')
plt.ylabel('Org_STP')


# In[210]:

plt.figure()
plt.plot(ams[2]['UTC_mid'],Sp.smooth(ams[2]['Org_STP'],14))
plt.plot(ams[2]['UTC_mid'],Sp.smooth(ams[2]['SO4_STP'],14),'-g')
plt.plot(ams[2]['UTC_mid'],Sp.smooth(ams[2]['Chl_STP'],14),'-c')

plt.plot(star[2]['Start_UTC'],star[2]['COD']/200.0,'xk')
plt.plot(star[2]['Start_UTC'],star[2]['REF']/100.0,'+r')
plt.xlim(15,15.8)


# In[206]:

plt.figure()
plt.plot(ams[4]['UTC_mid'],Sp.smooth(ams[4]['Org_STP'],14))
plt.plot(ams[4]['UTC_mid'],Sp.smooth(ams[4]['SO4_STP'],14),'-g')
plt.plot(star[4]['Start_UTC'],star[4]['COD']/100.0,'xk')
plt.xlim(12.5,14.1)


# ## Make correlation plots between COD and AMS values

# In[221]:

import plotting_utils as pu


# In[248]:

days


# In[128]:

plt.figure()
for i,d in enumerate(days):

    plt.plot(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],14),label=d)
    
plt.legend()


# In[50]:

fig,ax = plt.subplots(1,2)
ax = ax.ravel()
cs = ['c','y','b','r','g']
dd = ['11/12','11/14','11/17','11/18','11/23']
for i,d in enumerate(days):
    try:
        ic = (np.isfinite(star[i]['COD'])) & (star[i]['alt']<900.0)
    except:
        continue
    org_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    so4_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['SO4_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    
    roo = np.corrcoef(org_utc,star[i]['COD'][ic])
    rss = np.corrcoef(so4_utc,star[i]['COD'][ic])
    ro = roo[0,1]
    rs = rss[0,1]
    #print roo,rss
    
    ax[0].scatter(org_utc,star[i]['COD'][ic],c=cs[i],marker='x',label=dd[i]+' R={:1.2f}'.format(ro))
    ax[1].scatter(so4_utc,star[i]['COD'][ic],c=cs[i],marker='o',label=dd[i]+' R={:1.2f}'.format(rs),edgecolor='None')
    pu.plot_lin(org_utc,star[i]['COD'][ic],color=cs[i],labels=False,shaded_ci=False,ax=ax[0])
    pu.plot_lin(so4_utc,star[i]['COD'][ic],color=cs[i],labels=False,shaded_ci=False,ax=ax[1])
    

ax[0].set_ylabel('4STAR COD')
ax[0].set_xlabel('LARGE AMS Organics [$\\mu$g/m$^3$]')
ax[0].legend(frameon=True)
ax[1].set_xlabel('LARGE AMS SO4 [$\\mu$g/m$^3$]')
plt.legend(frameon=True)
ax[0].set_xlim(0,0.07)
ax[0].set_ylim(0,50)
ax[1].set_xlim(0,0.07)
ax[1].set_ylim(0,50)
#plt.savefig(fp+'plot/COD_vs_AMS.png',transparent=True,dpi=600)


# In[51]:

fig,ax = plt.subplots(1,2)
ax = ax.ravel()
cs = ['c','y','b','r','g']
dd = ['11/12','11/14','11/17','11/18','11/23']
for i,d in enumerate(days):
    try:
        ic = (np.isfinite(star[i]['COD'])) & (star[i]['alt']<900.0)
    except:
        continue
    org_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    so4_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['SO4_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    
    roo = np.corrcoef(org_utc,star[i]['REF'][ic])
    rss = np.corrcoef(so4_utc,star[i]['REF'][ic])
    ro = roo[0,1]
    rs = rss[0,1]
    #print roo,rss
    
    ax[0].scatter(org_utc,star[i]['REF'][ic],c=cs[i],marker='x',label=dd[i]+' R={:1.2f}'.format(ro))
    ax[1].scatter(so4_utc,star[i]['REF'][ic],c=cs[i],marker='o',label=dd[i]+' R={:1.2f}'.format(rs),edgecolor='None')
    pu.plot_lin(org_utc,star[i]['REF'][ic],color=cs[i],labels=False,shaded_ci=False,ax=ax[0])
    pu.plot_lin(so4_utc,star[i]['REF'][ic],color=cs[i],labels=False,shaded_ci=False,ax=ax[1])
    

ax[0].set_ylabel('4STAR REF')
ax[0].set_xlabel('LARGE AMS Organics [$\\mu$g/m$^3$]')
ax[0].legend(frameon=True)
ax[1].set_xlabel('LARGE AMS SO4 [$\\mu$g/m$^3$]')
plt.legend(frameon=True)
ax[0].set_xlim(0,0.07)
ax[0].set_ylim(0,50)
ax[1].set_xlim(0,0.07)
ax[1].set_ylim(0,50)
plt.savefig(fp+'plot/REF_vs_AMS.png',transparent=True,dpi=600)


# In[234]:

np.corrcoef(so4_utc,star[i]['COD'][ic])[0,1]**2


# In[ ]:




# In[213]:

np.isfinite(star[4]['COD'])


# In[ ]:

nearest_neighbor(hsk[i]['Start_UTC'],hsk[i]['GPS_Altitude'],cdp[i]['UTC_mid'],dist=1.0/3600.0)


# ## Plot out drizzling cases

# In[58]:

fig,ax = plt.subplots(1,2)
ax = ax.ravel()
cs = ['c','y','b','r','g']
dd = ['11/12','11/14','11/17','11/18','11/23']

for i,d in enumerate(days):
    try:
        ic = (np.isfinite(star[i]['COD'])) & (star[i]['alt']<900.0)
    except:
        continue
    cipn_utc = nearest_neighbor(cip[i]['UTC_mid'],Sp.smooth(cip[i]['nCIP'],14),star[i]['Start_UTC'][ic],dist=10.0/3600.0)
    print d,np.nanmean(cipn_utc)
    
    org_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    so4_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['SO4_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)

    roo = np.corrcoef(org_utc,cipn_utc)
    rss = np.corrcoef(so4_utc,cipn_utc)
    ro = roo[0,1]
    rs = rss[0,1]
    #print roo,rss
    
    ax[0].scatter(org_utc,cipn_utc,c=cs[i],marker='x',label=dd[i]+' R={:1.2f}'.format(ro))
    ax[1].scatter(so4_utc,cipn_utc,c=cs[i],marker='o',label=dd[i]+' R={:1.2f}'.format(rs),edgecolor='None')
    pu.plot_lin(org_utc,cipn_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[0])
    pu.plot_lin(so4_utc,cipn_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[1])
    

ax[0].set_ylabel('LARGE CIP number density of precipitating particles')
ax[0].set_xlabel('LARGE AMS Organics [$\\mu$g/m$^3$]')
ax[0].legend(frameon=True)
ax[1].set_xlabel('LARGE AMS SO4 [$\\mu$g/m$^3$]')
plt.legend(frameon=True)
ax[0].set_xlim(0,0.07)
ax[0].set_ylim(0,10)
ax[1].set_xlim(0,0.07)
ax[1].set_ylim(0,12)
plt.savefig(fp+'plot/nCIP_vs_AMS.png',transparent=True,dpi=600)


# ## Get the relationship between vertical winds and AMS, COD, REF

# In[66]:

wind[0].keys()


# In[71]:

fig,ax = plt.subplots(1,2)
ax = ax.ravel()
cs = ['c','y','b','r','g']
dd = ['11/12','11/14','11/17','11/18','11/23']

for i,d in enumerate(days):
    try:
        ic = (np.isfinite(star[i]['COD'])) & (star[i]['alt']<900.0)
    except:
        continue
    wind_utc = nearest_neighbor(wind[i]['UTC_mid'],Sp.smooth(wind[i]['w_ms1'],14),star[i]['Start_UTC'][ic],dist=10.0/3600.0)
    print d,np.nanmean(wind_utc)
    
    org_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    so4_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['SO4_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)

    roo = np.corrcoef(org_utc,wind_utc)
    rss = np.corrcoef(so4_utc,wind_utc)
    ro = roo[0,1]
    rs = rss[0,1]
    #print roo,rss
    
    ax[0].scatter(org_utc,wind_utc,c=cs[i],marker='x',label=dd[i]+' R={:1.2f}'.format(ro))
    ax[1].scatter(so4_utc,wind_utc,c=cs[i],marker='o',label=dd[i]+' R={:1.2f}'.format(rs),edgecolor='None')
    pu.plot_lin(org_utc,wind_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[0])
    pu.plot_lin(so4_utc,wind_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[1])
    

ax[0].set_ylabel('LARGE Vertical Winds [m/s]')
ax[0].set_xlabel('LARGE AMS Organics [$\\mu$g/m$^3$]')
ax[0].legend(frameon=True)
ax[1].set_xlabel('LARGE AMS SO4 [$\\mu$g/m$^3$]')
plt.legend(frameon=True)
ax[0].set_xlim(0,0.07)
ax[0].set_ylim(-1,1)
ax[1].set_xlim(0,0.07)
ax[1].set_ylim(-1,1)
plt.savefig(fp+'plot/winds_vs_AMS.png',transparent=True,dpi=600)


# In[76]:

fig,ax = plt.subplots(1,2)
ax = ax.ravel()
cs = ['c','y','b','r','g']
dd = ['11/12','11/14','11/17','11/18','11/23']

for i,d in enumerate(days):
    try:
        ic = (np.isfinite(star[i]['COD'])) & (star[i]['alt']<900.0)
    except:
        continue
    wind_utc = nearest_neighbor(wind[i]['UTC_mid'],Sp.smooth(wind[i]['w_ms1'],14),star[i]['Start_UTC'][ic],dist=10.0/3600.0)
    print d,np.nanmean(wind_utc)
    
    #org_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    #so4_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['SO4_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)

    roo = np.corrcoef(star[i]['COD'][ic],wind_utc)
    rss = np.corrcoef(star[i]['REF'][ic],wind_utc)
    ro = roo[0,1]
    rs = rss[0,1]
    #print roo,rss
    
    ax[0].scatter(star[i]['COD'][ic],wind_utc,c=cs[i],marker='x',label=dd[i]+' R={:1.2f}'.format(ro))
    ax[1].scatter(star[i]['REF'][ic],wind_utc,c=cs[i],marker='o',label=dd[i]+' R={:1.2f}'.format(rs),edgecolor='None')
    pu.plot_lin(star[i]['COD'][ic],wind_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[0])
    pu.plot_lin(star[i]['REF'][ic],wind_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[1])
    

ax[0].set_ylabel('LARGE Vertical Winds [m/s]')
ax[0].set_xlabel('4TAR COD')
ax[0].legend(frameon=True)
ax[1].set_xlabel('4STAR REF [$\\mu$m]')
plt.legend(frameon=True)
ax[0].set_xlim(0,60)
ax[0].set_ylim(-2,2)
ax[1].set_xlim(0,30)
ax[1].set_ylim(-2,2)
plt.savefig(fp+'plot/NAAMES_winds_vs_COD_REF.png',transparent=True,dpi=600)


# In[75]:

plt.figure()
plt.plot(wind[2]['UTC_mid'],wind[2]['w_ms1'],'.',color='grey')
plt.plot(wind[2]['UTC_mid'],Sp.smooth(wind[2]['w_ms1'],30),'-')
plt.xlim(14.8,16.2)


# ## Plot relationship between CCN and COD/REF

# In[85]:

cn_vals = ccn[0].keys()


# In[97]:

try: 
    cn_vals.remove('alt')
except:
    pass
try:
    cn_vals.remove('UTC_mid')
except:
    pass
cn_vals.sort()


# In[98]:

cn_vals


# In[101]:

for i,d in enumerate(days):
    plt.figure()
    plt.title('CCN: '+d)
    for g in cn_vals:
        plt.plot(ccn[i]['UTC_mid'],ccn[i][g],'.',label=g)
    plt.legend(frameon=False,numpoints=1)


# In[124]:

fig,ax = plt.subplots(1,2)
ax = ax.ravel()
cs = ['c','y','b','r','g']
dd = ['11/12','11/14','11/17','11/18','11/23']

for i,d in enumerate(days):
    try:
        ic = (np.isfinite(star[i]['COD'])) & (star[i]['alt']<900.0)
    except:
        continue
    wind_utc = nearest_neighbor(wind[i]['UTC_mid'],Sp.smooth(wind[i]['w_ms1'],14),star[i]['Start_UTC'][ic],dist=10.0/3600.0)
    print d,np.nanmean(wind_utc)
    ccn_utc = nearest_neighbor(ccn[i]['UTC_mid'],Sp.smooth(ccn[i]['CCN40'],10),star[i]['Start_UTC'][ic],dist=10.0/3600.0)
    
    #org_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    #so4_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['SO4_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)

    roo = np.corrcoef(star[i]['COD'][ic],ccn_utc)
    rss = np.corrcoef(star[i]['REF'][ic],ccn_utc)
    ro = roo[0,1]
    rs = rss[0,1]
    #print roo,rss
    
    ax[0].scatter(star[i]['COD'][ic],ccn_utc,c=cs[i],marker='x',label=dd[i]+' R={:1.2f}'.format(ro))
    ax[1].scatter(star[i]['REF'][ic],ccn_utc,c=cs[i],marker='o',label=dd[i]+' R={:1.2f}'.format(rs),edgecolor='None')
    pu.plot_lin(star[i]['COD'][ic],ccn_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[0])
    pu.plot_lin(star[i]['REF'][ic],ccn_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[1])
    

ax[0].set_ylabel('CCN number concentration at 0.375-0.425\% SS [cm$^{{-3}}$]')
ax[0].set_xlabel('4TAR COD')
ax[0].legend(frameon=True)
ax[1].set_xlabel('4STAR REF [$\\mu$m]')
plt.legend(frameon=True)
ax[0].set_xlim(0,60)
ax[0].set_ylim(0,800)
ax[1].set_xlim(0,30)
ax[1].set_ylim(0,800)
plt.savefig(fp+'plot/NAAMES_ccn40_vs_COD_REF.png',transparent=True,dpi=600)


# In[123]:

fig,ax = plt.subplots(1,2)
ax = ax.ravel()
cs = ['c','y','b','r','g']
dd = ['11/12','11/14','11/17','11/18','11/23']

for i,d in enumerate(days):
    try:
        ic = (np.isfinite(star[i]['COD'])) & (star[i]['alt']<900.0)
    except:
        continue
    wind_utc = nearest_neighbor(wind[i]['UTC_mid'],Sp.smooth(wind[i]['w_ms1'],14),star[i]['Start_UTC'][ic],dist=10.0/3600.0)
    print d,np.nanmean(wind_utc)
    ccn_utc = nearest_neighbor(ccn[i]['UTC_mid'],Sp.smooth(ccn[i]['CCN21'],10),star[i]['Start_UTC'][ic],dist=10.0/3600.0)
    
    #org_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['Org_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)
    #so4_utc = nearest_neighbor(ams[i]['UTC_mid'],Sp.smooth(ams[i]['SO4_STP'],14),star[i]['Start_UTC'][ic],dist=40.0/3600.0)

    roo = np.corrcoef(star[i]['COD'][ic],ccn_utc)
    rss = np.corrcoef(star[i]['REF'][ic],ccn_utc)
    ro = roo[0,1]
    rs = rss[0,1]
    #print roo,rss
    
    ax[0].scatter(star[i]['COD'][ic],ccn_utc,c=cs[i],marker='x',label=dd[i]+' R={:1.2f}'.format(ro))
    ax[1].scatter(star[i]['REF'][ic],ccn_utc,c=cs[i],marker='o',label=dd[i]+' R={:1.2f}'.format(rs),edgecolor='None')
    pu.plot_lin(star[i]['COD'][ic],ccn_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[0])
    pu.plot_lin(star[i]['REF'][ic],ccn_utc,color=cs[i],labels=False,shaded_ci=False,ax=ax[1])
    

ax[0].set_ylabel('CCN number concentration at 0.21\% SS [cm$^{{-3}}$]')
ax[0].set_xlabel('4TAR COD')
ax[0].legend(frameon=True)
ax[1].set_xlabel('4STAR REF [$\\mu$m]')
plt.legend(frameon=True)
ax[0].set_xlim(0,60)
ax[0].set_ylim(0,10)
ax[1].set_xlim(0,30)
ax[1].set_ylim(0,10)
plt.savefig(fp+'plot/NAAMES_ccn21_vs_COD_REF.png',transparent=True,dpi=600)


# ## Plot out latitude and altitude dependence

# ### Plot out for 2015-11-17

# In[92]:

plt.figure()
plt.plot(hsk[2]['Start_UTC'],hsk[2]['GPS_Altitude'])
plt.plot(star[2]['Start_UTC'],star[2]['alt'])


# In[103]:

star[2].keys()


# In[192]:

plt.figure()
a = plt.scatter(hsk[2]['Latitude'],hsk[2]['GPS_Altitude'],c=hsk[2]['Start_UTC'],edgecolor='None')
cb = plt.colorbar()
cb.set_label('UTC')

plt.plot(hsk[2]['Latitude'][(hsk[2]['Start_UTC']>14.9)&(hsk[2]['Start_UTC']<15.7)],
         hsk[2]['GPS_Altitude'][(hsk[2]['Start_UTC']>14.9)&(hsk[2]['Start_UTC']<15.7)],'sk')

plt.plot(cld[2]['lat'][cld[2]['CloudFlag']==1],cld[2]['alt'][cld[2]['CloudFlag']==1],'+g')
plt.plot(star[2]['LAT'],star[2]['alt'],'x')


# Lowest legs for 2015-11-17 occur at 14.9 to 15.7 UTCH

# In[265]:

plt.figure()
plt.scatter(hsk[2]['Longitude'],hsk[2]['Latitude'],c=hsk[2]['Start_UTC'],edgecolor='None')
cb = plt.colorbar()
cb.set_label('UTC')

plt.plot(hsk[2]['Longitude'][(hsk[2]['Start_UTC']>14.6)&(hsk[2]['Start_UTC']<16.55)],
         hsk[2]['Latitude'][(hsk[2]['Start_UTC']>14.6)&(hsk[2]['Start_UTC']<16.55)],'sk')


# In[254]:

plt.figure()
a = plt.scatter(hsk[2]['Longitude'],hsk[2]['GPS_Altitude'],c=hsk[2]['Start_UTC'],edgecolor='None')
cb = plt.colorbar()
cb.set_label('UTC')


# In[134]:

icdp = (cdp[2]['UTC_mid']>14.5) & (cdp[2]['UTC_mid']<17.5)


# In[174]:

irsp = (rsp[2]['UTC_start']/3600.0>14.5) & (rsp[2]['UTC_start']/3600.0<17.5)


# In[175]:

rsp[2].keys()


# In[153]:

matplotlib.cm.gist_ncar_r


# In[177]:

plt.figure()
plt.set_cmap('gist_ncar_r')
plt.pcolor(np.array(cdp_binsc)/2.0,cdp[2]['alt'][icdp],
           cdp_bins[2][:,icdp].transpose())
plt.plot(star[2]['REF'],star[2]['REF']*0.0+400.0,'*',label='4STAR')
plt.plot(rsp[2]['Reff'][irsp],rsp[2]['Reff'][irsp]*0.0+800.0,'+',label='RSP')
plt.plot(rsp[2]['Reff_rad_1590'][irsp],rsp[2]['Reff_rad_1590'][irsp]*0.0+800.0,
         'x',label='RSP Radiance')

plt.ylim(0,1500)

plt.xlabel('Particle Radius [$\\mu$m]')
plt.ylabel('Altitude [m]')

cb = plt.colorbar()
cb.set_label('dN/dlogD')


# ### Plot out for 2015-11-23

# In[276]:

plt.figure()
a = plt.scatter(hsk[4]['Latitude'],hsk[4]['GPS_Altitude'],c=hsk[4]['Start_UTC'],edgecolor='None')
cb = plt.colorbar()
cb.set_label('UTC')

plt.plot(hsk[4]['Latitude'][(hsk[4]['Start_UTC']>12.65)&(hsk[4]['Start_UTC']<13.8)],
         hsk[4]['GPS_Altitude'][(hsk[4]['Start_UTC']>12.65)&(hsk[4]['Start_UTC']<13.8)],'sk')

plt.plot(hsk[4]['Latitude'][(hsk[4]['Start_UTC']>12.1)&(hsk[4]['Start_UTC']<12.43)],
         hsk[4]['GPS_Altitude'][(hsk[4]['Start_UTC']>12.1)&(hsk[4]['Start_UTC']<12.43)],'+k')

plt.plot(cld[4]['lat'][cld[4]['CloudFlag']==1],cld[4]['alt'][cld[4]['CloudFlag']==1],'+g')
plt.plot(star[4]['LAT'],star[4]['alt'],'x')


# In[253]:

plt.figure()
a = plt.scatter(hsk[4]['Longitude'],hsk[4]['GPS_Altitude'],c=hsk[4]['Start_UTC'],edgecolor='None')
cb = plt.colorbar()
cb.set_label('UTC')

plt.plot(hsk[4]['Longitude'][(hsk[4]['Start_UTC']>14.9)&(hsk[4]['Start_UTC']<15.7)],
         hsk[4]['GPS_Altitude'][(hsk[4]['Start_UTC']>14.9)&(hsk[4]['Start_UTC']<15.7)],'sk')

plt.plot(cld[4]['lon'][cld[4]['CloudFlag']==1],cld[4]['alt'][cld[4]['CloudFlag']==1],'+g')
plt.plot(star[4]['LON'],star[4]['alt'],'x')


# In[194]:

rsp[4].keys()


# In[198]:

plt.figure()
plt.plot(star[4]['LAT'],star[4]['REF'],'xb',label='4STAR')
plt.plot(rsp[4]['Lat'],rsp[4]['Reff'],'+k',label='RSP')
plt.plot(rsp[4]['Lat'],rsp[4]['Reff_rad_1590'],'+r',label='RSP rad 1590nm')
plt.plot(rsp[4]['Lat'],rsp[4]['Reff_rad_2260'],'sc',label='RSP rad 2260nm')
plt.xlim(42.4,43.2)


# In[201]:

plt.figure()
plt.scatter(rsp[4]['Lat'],rsp[4]['alt'],c=rsp[4]['Reff'],edgecolor='None')
cb = plt.colorbar()
cb.set_label('ref')


# In[204]:

plt.figure()
plt.set_cmap('gist_ncar_r')
plt.pcolor(np.array(cdp_binsc)/2.0,cdp[4]['alt'][icdp],
           cdp_bins[4][:,icdp].transpose())
#plt.plot(star[2]['REF'],star[2]['REF']*0.0+400.0,'*',label='4STAR')
#plt.plot(rsp[2]['Reff'][irsp],rsp[2]['Reff'][irsp]*0.0+800.0,'+',label='RSP')
#plt.plot(rsp[2]['Reff_rad_1590'][irsp],rsp[2]['Reff_rad_1590'][irsp]*0.0+800.0,
#         'x',label='RSP Radiance')

plt.ylim(0,1500)

plt.xlabel('Particle Radius [$\\mu$m]')
plt.ylabel('Altitude [m]')

cb = plt.colorbar()
cb.set_label('dN/dlogD')


# ## Plot out the time dependence

# In[ ]:



