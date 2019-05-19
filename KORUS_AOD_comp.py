
# coding: utf-8

# # Info
# Name:  
# 
#     KORUS_AOD_comp
# 
# Purpose:  
# 
#     Comparison of AOD from 4STAR along flight track and GOCI aerosol
#     Additional calculations of the aerosol extinction profile
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
#     - load_utils.py : for loading OMI HDF5 files
#     - matplotlib
#     - numpy
#     - scipy : for saving and reading
#     - pytables
#     - os
#     - datetime
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, OSAN AFB, Korea, 2016-05-06
#     Modified: 

# # Import the required modules and set up base

# In[1]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.io as sio
import Sp_parameters as Sp
import load_utils as lm
import os
from datetime import datetime


# In[2]:


import hdf5storage as hs
import tables


# In[3]:


from mpl_toolkits.basemap import Basemap,cm
get_ipython().magic(u'matplotlib notebook')
fp = 'C:/Users/sleblan2/Research/KORUS-AQ/'


# In[43]:


daystr = '20160511'
daystr_goci = '20160512'


# # Load the various data

# ## Load the 4STAR starsun

# In[44]:


f_star = fp+'data\\{}starsun.mat'.format(daystr)


# In[45]:


try:
    s = sio.loadmat(f_star+'_small.mat')
except IOError:
    s = sio.loadmat(f_star)
except NotImplementedError:
    s = hs.loadmat(f_star+'_small.mat')


# In[46]:


s.keys()


# In[47]:


s['utc'] = lm.toutc(lm.mat2py_time(s['t']))


# In[48]:


s['tau_aero'].shape


# ### Load the starflag for this day

# In[49]:


f_info = 'C:\Users\sleblan2\Research\\4STAR_codes\data_folder\\'


# In[50]:


finf = f_info+'starinfo_{}.m'.format(daystr)


# In[51]:


with open(finf, 'r') as inF:
    for line in inF:
        if 'flagfilename' in line:
            f_flag = line[line.find("'")+1:line.rfind("'")]
            break


# In[52]:


f_flag


# In[53]:


sflag = sio.loadmat(f_info+f_flag)


# In[54]:


flag = sflag['manual_flags']['screen'][0][0][:,0]


# In[55]:


flag.shape


# ### Apply the flags to 4star data

# In[56]:


s['tau_aero'][flag==1,:]=np.nan


# ## Load the GOCI aerosol products

# In[57]:


fp_goci = fp+'sat/GOCI//{}/'.format(daystr_goci)


# In[58]:


fpl = os.listdir(fp_goci)


# In[59]:


fpl


# In[60]:


gg = []
gg_head = []
for f in fpl:
    f_goci = fp_goci+f
    gt,gth = lm.load_hdf(f_goci,values=(('lon',0),('lat',1),('aod550',2),('fmf550',3),('ssa440',4),('type',5),('ang',6),('QA',7),
                                ('obs_time',8),('cf',9),('turbidI',10),('Land_sea_mask',11)),verbose=False)
    gt['year'] = int(f.split('_')[-1].split('.')[0][0:4])
    gt['month'] = int(f.split('_')[-1].split('.')[0][4:6])
    gt['day'] = int(f.split('_')[-1].split('.')[0][6:8])
    gt['hour'] = int(f.split('_')[-1].split('.')[0][8:10])
    gt['minute'] = int(f.split('_')[-1].split('.')[0][10:12])
    gt['seconds'] = int(f.split('_')[-1].split('.')[0][12:14])
    gt['julian'] = float(datetime(gt['year'],gt['month'],gt['day'],gt['hour'],gt['minute'],
                                  gt['seconds']).timetuple().tm_yday)+gt['hour']/24.0+gt['minute']/60.0/24.0+gt['seconds']/3600.0/24.0
    
    gt['aod550_QA'] = gt['aod550']
    gt['aod550_QA'][gt['QA']==3] = np.nan
    gt['aod550_QA'][gt['aod550_QA']<0.0] = np.nan
    gg.append(gt)
    gg_head.append(gth)


# In[61]:


gt['aod550'].shape


# In[62]:


gt['QA'].shape


# In[63]:


gt.keys()


# In[253]:


gth.keys()


# In[254]:


gth['aod550']


# In[64]:


len(gg)


# In[65]:


fpl


# In[66]:


gg[4]['julian']


# ## Get the AERONET data to overlay on plot

# In[67]:


fp


# In[28]:


reload(lm)


# In[29]:


fa = fp+'aeronet/AOT/LEV10/ALL_POINTS/'


# In[30]:


aero = lm.load_multi_aeronet(fa)


# In[31]:


ilatest = lm.aeronet_subset(aero)


# In[32]:


plt.figure()
plt.scatter(aero['long'],aero['lat'],c=aero['AOT_500'][ilatest],
            cmap=plt.cm.rainbow,marker='s',vmin=0.0,vmax=1.5,edgecolors='None',s=30)


# In[33]:


plt.figure()
m = Basemap(projection='stere',lon_0=128,lat_0=36.0,
            llcrnrlon=123.0, llcrnrlat=33.5,
            urcrnrlon=132.0, urcrnrlat=39,resolution='h')
m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
m.drawstates()
m.drawcountries()
m.drawmeridians(np.linspace(123,133,11),labels=[0,0,0,1])
m.drawparallels(np.linspace(33,39,13),labels=[1,0,0,0])
x,y = m(aero['long'],aero['lat'])
bb = m.scatter(x,y,c=aero['AOT_500'][ilatest],
            cmap=plt.cm.rainbow,marker='s',vmin=0.0,vmax=1.5,edgecolors='None',s=30)
cbar = m.colorbar(bb)
cbar.set_label('AOD 500 nm')


# ## Subset the aeronet and 4STAR values to GOCI values

# In[68]:


utcs = []
iaero = []
istar = []
istar_steps = []
dz = 500.0
for i in range(len(gg)):
    utcs.append((gg[i]['julian']-np.floor(gg[i]['julian']))*24.0)
    iaero.append(lm.aeronet_subset(aero,julian=gg[i]['julian'],window=1.0))
    istar.append(((s['utc']-24.0)<utcs[i])&(s['Alt'][:,0]<dz))
    if i >0:
        istar_steps.append(((s['utc']-24.0)<utcs[i])&(s['Alt'][:,0]<dz)&((s['utc']-24.0)>utcs[i-1]))
    else:
        istar_steps.append(((s['utc']-24.0)<utcs[i])&(s['Alt'][:,0]<dz))


# In[69]:


utcs


# In[70]:


for i,u in enumerate(utcs):
    print istar[i].sum()


# # Start making different plots/maps

# In[71]:


#set up a easy plotting function
def make_map(ax=plt.gca()):
    m = Basemap(projection='stere',lon_0=128,lat_0=36.0,
            llcrnrlon=123.0, llcrnrlat=33.5,
            urcrnrlon=132.0, urcrnrlat=39,resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(123,133,11),labels=[0,0,0,1])
    m.drawparallels(np.linspace(33,39,13),labels=[1,0,0,0])
    return m


# ## Start with simple map plot of GOCI

# In[225]:


fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg[0]['lon'],gg[0]['lat'])
clevels = np.linspace(-0.2,0.6,41)

plt.title('GOCI AOD {} {:5.2f}'.format(daystr_goci,utcs[0]))
cs1 = m.contourf(x,y,gg[0]['aod550'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('AOD 550 nm')

#xx,yy = m(star['lon'],star['lat'])
#m.scatter(xx,yy,c=star['tau'],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
#          alpha=0.5,edgecolors='k',linewidth=0.65)
plt.savefig(fp+'plot/{}_GOCI_map_AOD.png'.format(daystr_goci),dpi=600,transparent=True)


# In[263]:


fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg[0]['lon'],gg[0]['lat'])
clevels = np.linspace(0,4,41)

plt.title('GOCI AOD {} {:5.2f}h'.format(daystr_goci,utcs[0]))
cs1 = ax.pcolorfast(x,y,gg[0]['aod550'][:-1,:-1],cmap=plt.cm.rainbow,vmin=-0.1,vmax=0.5)
cbar = m.colorbar(cs1)
cbar.set_label('AOD 550 nm')
plt.savefig(fp+'plot/{}_GOCI_map_AOD.png'.format(daystr_goci),dpi=600,transparent=True)


# In[56]:


fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg['lon'],gg['lat'])
clevels = np.linspace(0,40,41)

plt.title('GOCI AOD 2016-05-05 04:16:44')
cs1 = m.contourf(x,y,gg['obs_time'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('Observation Time')


# In[59]:


fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg['lon'],gg['lat'])
clevels = np.linspace(0,3,4)

plt.title('GOCI AOD 2016-05-05 04:16:44')
cs1 = m.contourf(x,y,gg['QA'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('Quality flag')


# ## Overlay 4STAR values

# In[36]:


plt.figure()
plt.plot(s['utc'],s['tau_aero'][:,450])
plt.ylabel('4STAR AOD at {} nm'.format(s['w'][0][469]*1000.0))
plt.xlabel('UTC [H]')


# In[61]:


ig = gg['QA']==3


# In[62]:


ig.shape


# In[63]:


gg['aod550'].shape


# In[64]:


gg['aod550'][ig]=np.nan


# In[67]:


fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg['lon'],gg['lat'])
clevels = np.linspace(0,4,41)

plt.title('GOCI AOD 2016-05-05 04:16:44')
cs1 = m.contourf(x,y,gg['aod550'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('AOD 550 nm')
m.scatter(x,y,c=gg['aod550'],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],edgecolors='None')


xx,yy = m(s['Lon'],s['Lat'])
m.scatter(xx,yy,c=s['tau_aero'][:,469],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
          alpha=0.5,edgecolors='None')
plt.savefig(fp+'plot/20160505_GOCI_4STAR_map_AOD.png',dpi=600,transparent=True)


# ## Overlay Aeronet AOD

# In[37]:


fp


# In[38]:


fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg['lon'],gg['lat'])
clevels = np.linspace(0,4,41)

plt.title('GOCI AOD 2016-05-05 04:16:44')
cs1 = m.contourf(x,y,gg['aod550'],clevels,cmap=plt.cm.rainbow,extend='max')
cbar = m.colorbar(cs1)
cbar.set_label('AOD 550 nm')
m.scatter(x,y,c=gg['aod550'],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],edgecolors='None')


xx,yy = m(s['Lon'],s['Lat'])
m.scatter(xx,yy,c=s['tau_aero'][:,469],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
          alpha=0.5,edgecolors='None')

xa,ya = m(anet['long'],anet['lat'])
m.scatter(xa,ya,c=anet['AOT_500'][il],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],
          alpha=1.0,edgecolors='m',s=40,linewidth=2)

#plt.savefig(fp+'plot/20160505_GOCI_4STAR_map_AOD.png',dpi=600,transparent=True)


# ## Make a GOCI figure with overlays for every time

# In[151]:


import scipy.ndimage as snim


# In[206]:


istar[i][0]


# In[343]:


plt.figure()
plt.plot(s['utc'][istar[i]],s['tau_aero'][istar[i],469],'x')


# In[361]:


u = utcs[-1]
i = 7
fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg[0]['lon'],gg[0]['lat'])
clevels = np.linspace(0.0,0.6,31)

plt.title('AOD from 4STAR AERONET -- {}'.format(daystr_goci,u))
#cs1 = ax.pcolorfast(x,y,gg[i]['aod550_QA'][:-1,:-1],cmap=plt.cm.rainbow,vmin=clevels[0],vmax=clevels[-1])
#cs2 = m.contour(x,y,gg[i]['cf'],levels=np.linspace(0.0,1.0,11),cmap=plt.cm.gist_gray)
#cs1 = ax.pcolorfast(x,y,gg[i]['aod550_QA'][:-1,:-1],cmap=plt.cm.rainbow,vmin=clevels[0],vmax=clevels[-1])
#cbar = m.colorbar(cs1,pad='10%')
#cbar.set_label('AOD 550 nm')
#cbar2 = m.colorbar(cs2,pad='7%',location='bottom')
#cbar2.set_label('Cloud Fraction ')
#m.scatter(x,y,c=gg[i]['aod550'],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],edgecolors='None')

xa,ya = m(aero['long'],aero['lat'])
m.scatter(xa,ya,c=aero['AOT_500'][iaero[i]],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],
          alpha=1.0,edgecolors='m',s=np.zeros_like(xa)+120,linewidth=3)


xs,ys = m(s['Lon'],s['Lat'])
m.plot(xs,ys,'.b',linewidth=0.1,markersize=0.1)
xx,yy = m(s['Lon'][istar[i]],s['Lat'][istar[i]])
cc = m.scatter(xx,yy,c=s['tau_aero'][istar[i],469],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
          alpha=0.8,edgecolors='k',linewidth=0.0,s=s['tau_aero'][istar[i],469]*30+40)
m.scatter(xx[-1],yy[-1],c='w',edgecolor='k',linewidth=4,marker='o',alpha=1.0,s=100)
cbar = m.colorbar(cc,pad='10%')
cbar.set_label('AOD 550 nm')

plt.savefig(fp+'plot/{}_4STAR_aeronet_map_AOD.png'.format(daystr),dpi=600,transparent=True)


# In[363]:


u = utcs[-1]
i = 7
fig,ax = plt.subplots(1,1,figsize=(11,8))
m = make_map(ax)
x,y = m(gg[0]['lon'],gg[0]['lat'])
clevels = np.linspace(0.0,0.6,31)

plt.title('AOD from 4STAR AERONET -- {}'.format(daystr_goci,u))
#cs1 = ax.pcolorfast(x,y,gg[i]['aod550_QA'][:-1,:-1],cmap=plt.cm.rainbow,vmin=clevels[0],vmax=clevels[-1])
#cs2 = m.contour(x,y,gg[i]['cf'],levels=np.linspace(0.0,1.0,11),cmap=plt.cm.gist_gray)
#cs1 = ax.pcolorfast(x,y,gg[i]['aod550_QA'][:-1,:-1],cmap=plt.cm.rainbow,vmin=clevels[0],vmax=clevels[-1])
#cbar = m.colorbar(cs1,pad='10%')
#cbar.set_label('AOD 550 nm')
#cbar2 = m.colorbar(cs2,pad='7%',location='bottom')
#cbar2.set_label('Cloud Fraction ')
#m.scatter(x,y,c=gg[i]['aod550'],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],edgecolors='None')

cs2 = m.contourf(x,y,gg[i]['cf'],levels=np.linspace(0.0,1.0,11),cmap=plt.cm.gist_gray)
cbar2 = m.colorbar(cs2,pad='7%',location='bottom')
cbar2.set_label('Cloud Fraction ')
xa,ya = m(aero['long'],aero['lat'])
m.scatter(xa,ya,c=aero['AOT_500'][iaero[i]],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],
          alpha=1.0,edgecolors='m',s=np.zeros_like(xa)+120,linewidth=3)


xs,ys = m(s['Lon'],s['Lat'])
m.plot(xs,ys,'.b',linewidth=0.1,markersize=0.1)
xx,yy = m(s['Lon'][istar[i]],s['Lat'][istar[i]])
cc = m.scatter(xx,yy,c=s['tau_aero'][istar[i],469],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
          alpha=0.8,edgecolors='k',linewidth=0.0,s=s['tau_aero'][istar[i],469]*30+40)
m.scatter(xx[-1],yy[-1],c='w',edgecolor='k',linewidth=4,marker='o',alpha=1.0,s=100)
cbar = m.colorbar(cc,pad='10%')
cbar.set_label('AOD 550 nm')

plt.savefig(fp+'plot/{}_4STAR_aeronet_map_AOD_cf.png'.format(daystr),dpi=600,transparent=True)


# In[365]:


for i,u in enumerate(utcs):
    fig,ax = plt.subplots(1,1,figsize=(11,8))
    m = make_map(ax)
    x,y = m(gg[0]['lon'],gg[0]['lat'])
    clevels = np.linspace(0.0,0.6,31)

    plt.title('AOD from 4STAR GOCI AERONET -- {} at {:5.2f}H UTC'.format(daystr_goci,u))
    #cs1 = ax.pcolorfast(x,y,gg[i]['aod550_QA'][:-1,:-1],cmap=plt.cm.rainbow,vmin=clevels[0],vmax=clevels[-1])
    #cs2 = m.contour(x,y,gg[i]['cf'],levels=np.linspace(0.0,1.0,11),cmap=plt.cm.gist_gray)
    cs1 = ax.pcolorfast(x,y,gg[i]['aod550_QA'][:-1,:-1],cmap=plt.cm.rainbow,vmin=clevels[0],vmax=clevels[-1])
    cbar = m.colorbar(cs1,pad='10%')
    cbar.set_label('AOD 550 nm')
    #cbar2 = m.colorbar(cs2,pad='7%',location='bottom')
    #cbar2.set_label('Cloud Fraction ')
    #m.scatter(x,y,c=gg[i]['aod550'],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],edgecolors='None')
    
    xa,ya = m(aero['long'],aero['lat'])
    m.scatter(xa,ya,c=aero['AOT_500'][iaero[i]],cmap=plt.cm.rainbow,marker='s',vmin=clevels[0],vmax=clevels[-1],
              alpha=1.0,edgecolors='m',s=np.zeros_like(xa)+120,linewidth=3)
    
    xx,yy = m(s['Lon'][istar[i]],s['Lat'][istar[i]])
    m.scatter(xx,yy,c=s['tau_aero'][istar[i],469],cmap=plt.cm.rainbow,marker='o',vmin=clevels[0],vmax=clevels[-1],
              alpha=0.5,edgecolors='k',linewidth=0.0,s=s['tau_aero'][istar[i],469]*30+50)
    m.scatter(xx[-1],yy[-1],c='w',edgecolor='k',linewidth=4,marker='o',alpha=1.0,s=100)
    #break
    plt.savefig(fp+'plot/{}_GOCI_4STAR_map_AOD_{:4.2f}h.png'.format(daystr,u),dpi=600,transparent=True)


# In[351]:


ga.shape


# In[237]:


x.shape


# In[238]:


y.shape


# # Get the GOCI values along the flight path

# In[72]:


import map_utils as mu


# In[73]:


goci_ind = mu.map_ind(gg[0]['lon'],gg[0]['lat'],s['Lon'],s['Lat'])


# In[74]:


print goci_ind.shape
print s['Lon'].shape
print gg[0]['lon'].shape


# ## Plot the GOCI vs. 4STAR AOD values

# In[76]:


s_aero = s['tau_aero'][istar[-1],463]
s_utc = s['utc'][istar[-1]]


# In[77]:


plt.figure()
plt.plot(s_utc[s_aero<0.45],s_aero[s_aero<0.45],'xr',label='4STAR')
goci_aod = []
for i,u in enumerate(utcs):
    plt.plot(s['utc'][istar_steps[i]],gg[i]['aod550_QA'][goci_ind[0,istar_steps[i]],goci_ind[1,istar_steps[i]]],'s',
             label='GOCI YAER {:5.2f}h'.format(u))
    goci_aod.append(gg[i]['aod550_QA'][goci_ind[0,istar_steps[i]],goci_ind[1,istar_steps[i]]])
plt.xlabel('UTC [h]')
plt.ylabel('AOD 500 nm')
plt.title('Aerosol optical depth from 4STAR and GOCI on {}'.format(daystr))
box = plt.gca().get_position()
plt.gca().set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.grid()
plt.legend(frameon=True,numpoints=1,bbox_to_anchor=(1.5,0.9))
plt.savefig(fp+'plot/{}_GOCI_4STAR_utc_AOD.png'.format(daystr),dpi=600,transparent=True)


# In[357]:


gaod = np.concatenate(goci_aod).flatten()


# In[358]:


bins = np.linspace(0,0.5,31)


# In[359]:


plt.figure()
plt.hist(s['tau_aero'][istar[-1],463],bins=bins,histtype='stepfilled',normed=True,alpha=0.5,color='r',label='4STAR',edgecolor='None')
plt.hist(gaod,bins=bins,histtype='stepfilled',normed=True,alpha=0.5,color='b',label='GOCI',edgecolor='None')
plt.legend(frameon=False)
plt.xlabel('AOD 500 nm')
plt.ylabel('Normalized distribution')
plt.title('4STAR GOCI normalized distribution for {}'.format(daystr))
plt.savefig(fp+'plot/{}_GOCI_4STAR_AOD_hist.png'.format(daystr),dpi=600,transparent=True)


# # Save the comparison data to a file

# In[136]:


import hdf5storage as hs


# In[137]:


dict_out = {u's':s,u'istar':istar,u'istar_steps':istar_steps,
            u'goci':gg,u'utcs':utcs,u'goci_ind':goci_ind,u'gaod':gaod,
            u'aero':aero,u'iaero':iaero}


# In[ ]:


import cPickle as pickle
pickle.dump(dict_out,open(fp+'{}_GOCI_Aeronet_4STAR.p'.format(daystr),"wb"))


# In[138]:


hs.savemat(fp+'{}_GOCI_Aeronet_4STAR.mat',dict_out)

