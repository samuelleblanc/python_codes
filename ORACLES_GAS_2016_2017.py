
# coding: utf-8

# # Intro 
# 
# Load the Gas retrievals from ORACLES 2016 and 2017 and make a few figures

# # Load the python modules and set up paths

# In[1]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
import os
matplotlib.rc_file(os.path.join(os.getcwd(),'file.rc'))
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth
from linfit import linfit
from path_utils import getpath
from plotting_utils import make_boxplot


# In[118]:


from mpl_toolkits.basemap import Basemap


# In[38]:


import hdf5storage as hs


# In[68]:


import scipy.stats as st
import Sun_utils as su


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


fp =getpath('ORACLES')#'C:/Userds/sleblan2/Research/ORACLES/'
fp


# # Load the files

# ## Load the 2016 files

# In[29]:


days6 = ['824','825','827','830','831','902','904','906','908','910','912','914','918','920','924','925','927','929','930']


# In[30]:


v6 = 'R1'
outgas6 = []
outgas6_head = []
for i,d in enumerate(days6):
    try:
        print 'Doing day: {}'.format(d)
        fname_aod = fp+'/gas_ict/4STAR-GAS_P3_20160{}_{vv}.ict'.format(d,vv=v6)
        tt,th = load_ict(fname_aod,return_header=True)
    except:
        print '*** Problem with day: {} *** Skipping '.format(d)
        days.pop(i)
        continue
    
    outgas6.append(tt)
    outgas6_head.append(th)


# In[32]:


nm6 = outgas6[0].dtype.names


# ### Save into a single array

# In[93]:


ar6 = {}
for n in nm6:
    ar6[n] = np.array([])


# In[94]:


ar6['idays'] = np.array([])
ar6['days'] = np.array([])


# In[95]:


for i,d in enumerate(days6):
    ar6['idays'] = np.append(ar6['idays'],np.zeros_like(outgas6[i]['Start_UTC'])+i)
    ar6['days'] = np.append(ar6['days'],np.zeros_like(outgas6[i]['Start_UTC'])+float(d))
    for n in nm6:
        ar6[n] = np.append(ar6[n],outgas6[i][n])


# In[96]:


len(ar6['days'])


# In[39]:


hs.savemat(fp+'/gas_ict/all_gas_ict_{vv}_2016.mat'.format(vv=v6),ar6)


# ## Load the 2017 files

# In[16]:


days = ['20170801','20170802','20170807','20170809', '20170812','20170813',
        '20170815','20170817','20170818','20170819','20170821',
        '20170824','20170826','20170828','20170830','20170831','20170902','20170903','20170904']


# In[17]:


vv = 'R0'


# In[18]:


outgas = []
outgas_head = []
for i,d in enumerate(days):
    try:
        print 'Doing day: {}'.format(d)
        fname_aod = fp+'/gas_ict_2017/4STAR-GAS_P3_{}_{vv}.ict'.format(d,vv=vv)
        tt,th = load_ict(fname_aod,return_header=True)
    except:
        print '*** Problem with day: {} *** Skipping '.format(d)
        days.pop(i)
        continue
    
    outgas.append(tt)
    outgas_head.append(th)


# In[20]:


outgas_head[0]


# In[19]:


len(outgas)


# In[23]:


nm = outgas[0].dtype.names


# ### Combine data into single array

# In[40]:


ar = {}
for n in nm:
    ar[n] = np.array([])


# In[41]:


ar['idays'] = np.array([])
ar['days'] = np.array([])


# In[42]:


for i,d in enumerate(days):
    ar['idays'] = np.append(ar['idays'],np.zeros_like(outgas[i]['Start_UTC'])+i)
    ar['days'] = np.append(ar['days'],np.zeros_like(outgas[i]['Start_UTC'])+float(d))
    for n in nm:
        ar[n] = np.append(ar[n],outgas[i][n])


# In[43]:


hs.savemat(fp+'/gas_ict_2017/all_gas_ict_{vv}_2017.mat'.format(vv=vv),ar)


# In[45]:


ar.keys()


# ## Flag by altitude

# In[53]:


ar['fla'] = (ar['GPS_Alt']>600.0) & (ar['GPS_Alt']<1800.0)


# In[99]:


ar6['fla'] = (ar6['GPS_Alt']>600.0) & (ar6['GPS_Alt']<1800.0)


# In[55]:


ar['fl_O3'] = (ar['QA_O3']==0) & ar['fla']
ar['fl_CWV'] = (ar['QA_CWV']==0) & ar['fla']
ar['fl_HCOH'] = (ar['QA_HCOH']==0) & ar['fla']
ar['fl_NO2'] = (ar['QA_NO2']==0) & ar['fla']


# In[100]:


ar6['fl_O3'] = (ar6['QA_O3']==0) & ar6['fla']
ar6['fl_CWV'] = (ar6['QA_CWV']==0) & ar6['fla']
ar6['fl_HCOH'] = (ar6['QA_HCOH']==0) & ar6['fla']
ar6['fl_NO2'] = (ar6['QA_NO2']==0) & ar6['fla']


# # Now plot the statistics

# In[44]:


from plotting_utils import prelim


# ## Plot some histogram

# In[169]:


plt.figure()
plt.hist([ar6['VCD_O3'][ar6['fl_O3']],ar['VCD_O3'][ar['fl_O3']]*1.2],
         bins=15,label=['2016','2017'],edgecolor='None',alpha=0.75,normed=False)
plt.legend(frameon=False)
plt.xlabel('O3 Vertical column density [DU]')
plt.ylabel('O3 Samples')


# In[58]:


plt.figure()
plt.hist([ar6['VCD_NO2'][ar6['fl_NO2']],ar['VCD_NO2'][ar['fl_NO2']]],
         bins=15,label=['2016','2017'],edgecolor='None',alpha=0.75,normed=False)
plt.legend(frameon=False)
plt.xlabel('NO2 Vertical column density [DU]')
plt.ylabel('NO2 Samples')


# In[63]:


np.nanmin(ar['CWV'][ar['fl_CWV']])


# In[64]:


plt.figure()
plt.hist([ar6['CWV'][ar6['fl_CWV']],ar['CWV'][ar['fl_CWV']]],
         bins=15,range=[0,3.7],label=['2016','2017'],edgecolor='None',alpha=0.75,normed=False)
plt.legend(frameon=False)
plt.xlabel('Column Water Vapor [g/cm$^2$]')
plt.ylabel('CWV Samples')


# In[66]:


plt.figure()
plt.hist([ar6['VCD_HCOH'][ar6['fl_HCOH']],ar['VCD_HCOH'][ar['fl_HCOH']]],
         bins=15,label=['2016','2017'],edgecolor='None',alpha=0.75,normed=False)
plt.legend(frameon=False)
plt.xlabel('Formaldehyde Vertical column density [DU]')
plt.ylabel('HCOH Samples')


# ## Prep statistics for 2d histograms on a map

# In[67]:


sta_o3,sta6_o3 = {},{}
sta_no2,sta6_no2 = {},{}
sta_cwv,sta6_cwv = {},{}


# In[184]:


sta_o3['mean'],sta_o3['x'],sta_o3['y'],sta_o3['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_O3']],
                                    ar['Longitude'][ar['fl_O3']],ar['VCD_O3'][ar['fl_O3']]*1.3,bins=26,range=[[-15,2],[-16,8]])
sta_o3['mean'] = np.ma.masked_array(sta_o3['mean'],np.isnan(sta_o3['mean']))
sta_no2['mean'],sta_no2['x'],sta_no2['y'],sta_no2['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_NO2']],
                                    ar['Longitude'][ar['fl_NO2']],ar['VCD_NO2'][ar['fl_NO2']],bins=26,range=[[-15,2],[-16,8]])
sta_no2['mean'] = np.ma.masked_array(sta_no2['mean'],np.isnan(sta_no2['mean']))
sta_cwv['mean'],sta_cwv['x'],sta_cwv['y'],sta_cwv['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_CWV']],
                                    ar['Longitude'][ar['fl_CWV']],ar['CWV'][ar['fl_CWV']],bins=26,range=[[-15,2],[-16,8]])
sta_cwv['mean'] = np.ma.masked_array(sta_cwv['mean'],np.isnan(sta_cwv['mean']))


# In[170]:


sta_o3['std'],sta_o3['x'],sta_o3['y'],sta_o3['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_O3']],
                       ar['Longitude'][ar['fl_O3']],ar['VCD_O3'][ar['fl_O3']],bins=26,range=[[-15,2],[-16,8]],statistic=np.std)
sta_o3['std'] = np.ma.masked_array(sta_o3['std'],np.isnan(sta_o3['std']))
sta_no2['std'],sta_no2['x'],sta_no2['y'],sta_no2['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_NO2']],
                       ar['Longitude'][ar['fl_NO2']],ar['VCD_NO2'][ar['fl_NO2']],bins=26,range=[[-15,2],[-16,8]],statistic=np.std)
sta_no2['std'] = np.ma.masked_array(sta_no2['std'],np.isnan(sta_no2['std']))
sta_cwv['std'],sta_cwv['x'],sta_cwv['y'],sta_cwv['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_CWV']],
                       ar['Longitude'][ar['fl_CWV']],ar['CWV'][ar['fl_CWV']],bins=26,range=[[-15,2],[-16,8]],statistic=np.std)
sta_cwv['std'] = np.ma.masked_array(sta_cwv['std'],np.isnan(sta_cwv['std']))


# In[73]:


uniq_cnt = lambda x: len(np.unique(x))


# In[74]:


sta_o3['dcnt'],sta_o3['x'],sta_o3['y'],sta_o3['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_O3']],
                       ar['Longitude'][ar['fl_O3']],ar['days'][ar['fl_O3']],bins=26,range=[[-15,2],[-16,8]],statistic=uniq_cnt)
sta_o3['dcnt'] = np.ma.masked_array(sta_o3['dcnt'],np.isnan(sta_o3['dcnt']))
sta_no2['dcnt'],sta_no2['x'],sta_no2['y'],sta_no2['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_NO2']],
                       ar['Longitude'][ar['fl_NO2']],ar['days'][ar['fl_NO2']],bins=26,range=[[-15,2],[-16,8]],statistic=uniq_cnt)
sta_no2['dcnt'] = np.ma.masked_array(sta_no2['dcnt'],np.isnan(sta_no2['dcnt']))
sta_cwv['dcnt'],sta_cwv['x'],sta_cwv['y'],sta_cwv['bin'] = st.binned_statistic_2d(ar['Latitude'][ar['fl_CWV']],
                       ar['Longitude'][ar['fl_CWV']],ar['days'][ar['fl_CWV']],bins=26,range=[[-15,2],[-16,8]],statistic=uniq_cnt)
sta_cwv['dcnt'] = np.ma.masked_array(sta_cwv['dcnt'],np.isnan(sta_cwv['dcnt']))


# In[175]:


sta6_o3['mean'],sta6_o3['x'],sta6_o3['y'],sta6_o3['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_O3']],
                                    ar6['Longitude'][ar6['fl_O3']],ar6['VCD_O3'][ar6['fl_O3']]*0.9,bins=26,range=[[-25,-8],[0,16]])
sta6_o3['mean'] = np.ma.masked_array(sta6_o3['mean'],np.isnan(sta6_o3['mean']))
sta6_no2['mean'],sta6_no2['x'],sta6_no2['y'],sta6_no2['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_NO2']],
                                    ar6['Longitude'][ar6['fl_NO2']],ar6['VCD_NO2'][ar6['fl_NO2']],bins=26,range=[[-25,-8],[0,16]])
sta6_no2['mean'] = np.ma.masked_array(sta6_no2['mean'],np.isnan(sta6_no2['mean']))
sta6_cwv['mean'],sta6_cwv['x'],sta6_cwv['y'],sta6_cwv['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_CWV']],
                                    ar6['Longitude'][ar6['fl_CWV']],ar6['CWV'][ar6['fl_CWV']],bins=26,range=[[-25,-8],[0,16]])
sta6_cwv['mean'] = np.ma.masked_array(sta6_cwv['mean'],np.isnan(sta6_cwv['mean']))


# In[77]:


sta6_o3['std'],sta6_o3['x'],sta6_o3['y'],sta6_o3['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_O3']],
                                    ar6['Longitude'][ar6['fl_O3']],ar6['VCD_O3'][ar6['fl_O3']],bins=26,range=[[-25,-8],[0,16]],statistic=np.std)
sta6_o3['std'] = np.ma.masked_array(sta6_o3['std'],np.isnan(sta6_o3['std']))
sta6_no2['std'],sta6_no2['x'],sta6_no2['y'],sta6_no2['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_NO2']],
                                    ar6['Longitude'][ar6['fl_NO2']],ar6['VCD_NO2'][ar6['fl_NO2']],bins=26,range=[[-25,-8],[0,16]],statistic=np.std)
sta6_no2['std'] = np.ma.masked_array(sta6_no2['std'],np.isnan(sta6_no2['std']))
sta6_cwv['std'],sta6_cwv['x'],sta6_cwv['y'],sta6_cwv['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_CWV']],
                                    ar6['Longitude'][ar6['fl_CWV']],ar6['CWV'][ar6['fl_CWV']],bins=26,range=[[-25,-8],[0,16]],statistic=np.std)
sta6_cwv['std'] = np.ma.masked_array(sta6_cwv['std'],np.isnan(sta6_cwv['std']))


# In[104]:


sta6_o3['dcnt'],sta6_o3['x'],sta6_o3['y'],sta6_o3['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_O3']],
                                    ar6['Longitude'][ar6['fl_O3']],ar6['days'][ar6['fl_O3']],bins=26,range=[[-25,-8],[0,16]],statistic=uniq_cnt)
sta6_o3['dcnt'] = np.ma.masked_array(sta6_o3['dcnt'],np.isnan(sta6_o3['dcnt']))
sta6_no2['dcnt'],sta6_no2['x'],sta6_no2['y'],sta6_no2['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_NO2']],
                                    ar6['Longitude'][ar6['fl_NO2']],ar6['days'][ar6['fl_NO2']],bins=26,range=[[-25,-8],[0,16]],statistic=uniq_cnt)
sta6_no2['dcnt'] = np.ma.masked_array(sta6_no2['dcnt'],np.isnan(sta6_no2['dcnt']))
sta6_cwv['dcnt'],sta6_cwv['x'],sta6_cwv['y'],sta6_cwv['bin'] = st.binned_statistic_2d(ar6['Latitude'][ar6['fl_CWV']],
                                    ar6['Longitude'][ar6['fl_CWV']],ar6['days'][ar6['fl_CWV']],bins=26,range=[[-25,-8],[0,16]],statistic=uniq_cnt)
sta6_cwv['dcnt'] = np.ma.masked_array(sta6_cwv['dcnt'],np.isnan(sta6_cwv['dcnt']))


# ## Plot the stats on maps

# In[129]:


def mapfig(ax=plt.gca()):
    m = Basemap(projection='merc',llcrnrlat=-25,urcrnrlat=2,llcrnrlon=-16,urcrnrlon=16,resolution='l',ax=ax)
    m.drawcoastlines()
    m.drawmeridians(np.linspace(-16,16,9),labels=[0,0,0,1],linewidth=0.1)
    m.drawparallels(np.linspace(-25,2,10),labels=[1,0,0,0],linewidth=0.1)
    m.shadedrelief(alpha=0.4)
    return m


# In[160]:


ino2 = np.where(np.isfinite(sta_no2['mean']))
i6no2 = np.where(np.isfinite(sta6_no2['mean']))
io3 = np.where(np.isfinite(sta_o3['mean']))
i6o3 = np.where(np.isfinite(sta6_o3['mean']))


# In[148]:


i6o3


# In[150]:


sta6_o3['mean'].data.shape


# In[153]:


sta6_o3['mean'].data[i6o3,i6o3]


# In[162]:


sta6_o3['x'].data[i6o3[1]]


# In[163]:


plt.figure()
plt.scatter(np.array(sta6_o3['x'])[i6o3[1]],np.array(sta6_o3['y'])[i6o3[0]],50,c=sta6_o3['mean'].data[i6o3[0],i6o3[1]].flatten())


# In[185]:


fig,ax = plt.subplots(1,2,figsize=(12,5))
ax = ax.flatten()
ax1 = ax[0]
#ax1 = plt.subplot(1,2,1)
m = mapfig(ax=ax1)

mx,my = m(sta_no2['y'],sta_no2['x'])
p = ax1.pcolor(mx,my,sta_no2['mean'],vmin=0.0,vmax=2.0,cmap='plasma')
m6x,m6y = m(sta6_no2['y'],sta6_no2['x'])
p2 = ax1.scatter(np.array(m6x)[i6no2[1]],np.array(m6y)[i6no2[0]],50,c=sta6_no2['mean'].data[i6no2[0],i6no2[1]].flatten(),
               marker='o',edgecolor='None',cmap='plasma',vmin=0.0,vmax=2.0)
prelim()
ax1.set_title('NO2 VCD')
cb = m.colorbar(p,extend='max')
cb.set_label('NO2 Column density [DU]')

ax2 = ax[1]
#ax2 = plt.subplot(1,2,2)
m2 = mapfig(ax=ax2)
mxx,myy = m2(sta_o3['y'],sta_o3['x'])
p = ax2.pcolor(mxx,myy,sta_o3['mean'],vmin=200.0,vmax=400,cmap='viridis')
m6xx,m6yy = m(sta6_o3['y'],sta6_o3['x'])
p2 = ax2.scatter(np.array(m6xx)[i6o3[1]],np.array(m6yy)[i6o3[0]],
            50,c=sta6_o3['mean'].data[i6o3[0],i6o3[1]].flatten(),
#p2 = ax2.scatter(m6xx,m6yy,50,c=sta6_o3['mean'].data[i6o3[0],i6o3[1]],
               marker='o',edgecolor='None',cmap='viridis',vmin=200.0,vmax=400)
prelim()
ax2.set_title('O3 VCD')
cb = m2.colorbar(p,extend='both')
cb.set_label('O3 Column density [DU]')

#sizes = [10,100,500,1500]
labels = ['2016','2017']
points = [ax2.scatter([], [], s=50, c='grey',marker='o',edgecolor='None'),
          ax2.scatter([], [], s=50, c='grey',marker='s',edgecolor='None')]
plt.legend(points, labels, scatterpoints=1,frameon=False,loc='lower left')
    
#plt.suptitle('Above cloud Angstrom of 470 - 865 nm\n')
plt.subplots_adjust(left=0.05, right=0.94, top=0.9, bottom=0.06)

plt.savefig(fp+'plot_2017/ORACLES2017_4STAR_stats_O3_NO2_actualmap.png',
            transparent=True,dpi=500)

