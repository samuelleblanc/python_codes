#!/usr/bin/env python
# coding: utf-8

# # Intro
# Name:  
# 
#     ORACLES_Build_DARE
# 
# Purpose:  
# 
#     Build the aerosol radiative effect input files from the SSFR reflectances, 4STAR AOD, and 4STAR skyscan results
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
# Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, CA, 2019-12-02
#      
#      Modified: by Samuel LeBlanc, Santa Cruz, CA, 2020-11-18
#              - updated for R4 ORACLES 2016 and mie 2020-11-12
#      Modified: by Samuel LeBlanc, Santa Cruz, C, 2021-02-04
#              - updated with delta SSA +/-0.03 for v5

# # Import of modules

# In[2]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
import scipy.io as sio
from path_utils import getpath
import matplotlib.pyplot as plt
import load_utils as lu
from write_utils import nearest_neighbor, iterate_dict_unicode


# In[3]:


get_ipython().magic(u'matplotlib notebook')


# In[4]:


from tqdm.notebook import tqdm 
from datetime import datetime, timedelta


# In[5]:


from multiprocessing import Pool, cpu_count
from copy import deepcopy
import signal
import warnings
warnings.simplefilter('ignore')


# In[6]:


from scipy.interpolate import interp1d


# In[7]:


import Run_libradtran as Rl


# In[8]:


name = 'ORACLES'


# In[9]:


vv = 'v5'
vr = 'R4'


# In[10]:


fp = getpath(name)
fp_rtm = getpath('rtm')
fp_uvspec = getpath('uvspecb')+'uvspec'
matfile = fp+'{}_all_cld_ict.mat'.format(vr)
fp_uvspec_dat = getpath('uvspec_dat') 
fp_rtmdat = fp_rtm+'dat/'


# # Load the files

# ## Load the 4STAR AOD

# In[11]:


ar = hs.loadmat(fp+'/aod_ict/{v}/all_aod_ict_{v}_2016.mat'.format(v=vr,vi='v9'))


# In[12]:


ar.keys()


# In[13]:


ar['AOD0501'].shape


# In[12]:


sza = np.arccos(1.0/ar['amass_aer'])*180.0/np.pi


# In[13]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[14]:


len(days)


# In[15]:


ar['days']


# ## Load the retrieved Cloud properties

# In[16]:


cl = hs.loadmat(fp+'data_other/ssfr_2016_retrieved_COD_{}.mat'.format('v4'))
print 'loading: '+fp+'data_other/ssfr_2016_retrieved_COD_{}.mat'.format('v4')


# In[17]:


cl.keys()


# In[18]:


cl['tau'].shape


# In[19]:


dds = ['20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927']


# In[20]:


len(dds)


# In[21]:


cl['days']


# In[22]:


dd = np.unique(cl['days'])


# In[23]:


cod,ref = [],[]
for d in dd:
    print d
    fld = cl['days']==d
    fad = ar['days']==d+3.0
    #nearest neighbor, but not more than a minute away
    cod_tmp = nearest_neighbor(cl['utc'][fld],cl['tau'][fld],ar['Start_UTC'][fad],dist=1.0/60.0) 
    ref_tmp = nearest_neighbor(cl['utc'][fld],cl['ref'][fld],ar['Start_UTC'][fad],dist=1.0/60.0)
    cod = np.append(cod,cod_tmp)
    ref = np.append(ref,ref_tmp)


# In[24]:


cod.shape


# In[25]:


len(dds)


# In[26]:


len(np.unique(ar['days']))


# ## Load the skyscan retrievals

# In[27]:


try:
    ae,ae_dict = lu.load_netcdf(fp+'aeroinv_2016/netcdf4/4STAR-aeroinv_P3_2016_R0.nc',everything=True)
except:
    import h5py as h5
    f5 = h5.File(fp+'aeroinv_2016/netcdf4/4STAR-aeroinv_P3_2016_R0.nc','r')
    ae5 = {}
    ae5_dict = {}
    for ka,kd in f5.iteritems():
        ae5[ka] = kd.value
        ae5_dict[ka] = {}
        for kdict in kd.attrs.iteritems():
            if type(kdict[1])!=type(np.array([])):
                ae5_dict[ka][kdict[0]] = kdict[1]
    ae = ae5
    ae_dict = ae5_dict


# In[28]:


ke = ae.keys()
ke.sort()
ke


# In[29]:


ae['AOD_meas'][0]


# In[30]:


ae_dict['AAOD']


# In[31]:


ae_dict['SSA']


# In[32]:


ae['SSA'].shape


# In[33]:


ae_dict['time']


# In[34]:


ae['time']/3600.0


# In[35]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[36]:


ar['doy'] = np.array([datetime.strptime(days[int(d)],'%Y%m%d').timetuple().tm_yday for d in ar['days']])


# In[37]:


datetime.strptime(days[4],'%Y%m%d').timetuple().tm_yday


# In[38]:


ar['time_ae'] = ar['Start_UTC']+(24.0*(ar['doy']-244))


# In[39]:


ar['time_ae']


# ## Load the MODIS cloud retrievals

# In[40]:


import pytz


# In[41]:


fp


# In[44]:


fmod = os.listdir(fp+'data_other/MODIS/MOD06_L2/2016/268/')
fgeo = os.listdir(fp+'data_other/MODIS/MOD03/2016/268/')
fmod.sort()
fgeo.sort()


# In[45]:


vals = (('QA',123),('CF',38),('ref',66),('COD',72),('mask',110),('t',0))


# In[46]:


geo,geo_dict = lu.load_hdf(fp+'data_other/MODIS/MOD03/2016/268/'+fgeo[0],values=(('lat',12),('lon',13)),verbose=False)


# In[47]:


def bits_stripping(bit_start,bit_count,value):
	bitmask=pow(2,bit_start+bit_count)-1
	return np.right_shift(np.bitwise_and(value,bitmask),bit_start)


# In[48]:


def mod_qa(aq_val):
    qa_cod_useful = bits_stripping(0,1,aq_val[4])
    qa_cod_conf = bits_stripping(1,2,aq_val[4])
    qa_ref_useful = bits_stripping(3,1,aq_val[4])
    qa_ref_conf = bits_stripping(4,2,aq_val[4])
    return (qa_cod_useful>0) & (qa_cod_conf>1) & (qa_ref_useful>0) & (qa_ref_conf>1) 


# In[49]:


def mod_qa_arr(aq_arr):
    qa_cod_useful = bits_stripping(0,1,aq_arr[:,:,4])
    qa_cod_conf = bits_stripping(1,2,aq_arr[:,:,4])
    qa_ref_useful = bits_stripping(3,1,aq_arr[:,:,4])
    qa_ref_conf = bits_stripping(4,2,aq_arr[:,:,4])
    return (qa_cod_useful>0) & (qa_cod_conf>1) & (qa_ref_useful>0) & (qa_ref_conf>1) 


# ### pre-load the geo for easier calcs

# In[50]:


from pyhdf.SD import SD, SDC


# In[51]:


from copy import deepcopy


# In[59]:


if not 'fla' in locals(): 
    # for all 4STAR aerosol arrays
    fla = (ar['flag_acaod']==1) & ar['fl'] & ar['fl_QA'] & (ar['days']>2.0)


# In[60]:


geoso = {}
doys = np.unique(ar['doy'][fla])
for ddy in doys:
    print ddy
    fgeo = os.listdir(fp+'data_other/MODIS/MOD03/2016/{}/'.format(ddy))
    fgeo.sort()
    geo_time = np.array([float(g.split('.')[2][0:2])+float(g.split('.')[2][2:])/60.0 for g in fgeo])
    geoso['{}'.format(ddy)] = {'time':geo_time[:],'fname':fgeo[:],
                               'path':fp+'data_other/MODIS/MOD03/2016/{}/'.format(ddy),
                               'lat':np.zeros((len(geo_time),2030,1354))+np.nan,
                               'lon':np.zeros((len(geo_time),2030,1354))+np.nan}
    for ngeo,fg in enumerate(fgeo):
        geo_tmp,geo_tmp_dict = lu.load_hdf(geoso['{}'.format(ddy)]['path']+fg,values=(('lat',12),('lon',13)),verbose=False)
        geoso['{}'.format(ddy)]['lat'][ngeo,:,:] = geo_tmp['lat'][:2030,:1354]
        geoso['{}'.format(ddy)]['lon'][ngeo,:,:] = geo_tmp['lon'][:2030,:1354]
    # load the MOD06 file pointers
    fmods = os.listdir(fp+'data_other/MODIS/MOD06_L2/2016/{}/'.format(ddy))
    fmods.sort()
    fh = []
    for ff in fmods:
        fh1 = SD(str(fp+'data_other/MODIS/MOD06_L2/2016/{}/'.format(ddy)+ff),SDC.READ)
        fh.append(fh1)
    geoso['{}'.format(ddy)]['fh'] = fh[:]


# In[53]:


geosy = {}
doys = np.unique(ar['doy'][fla])
for ddy in doys:
    print ddy
    fgeo = os.listdir(fp+'data_other/MODIS/MYD03/2016/{}/'.format(ddy))
    fgeo.sort()
    geo_time = np.array([float(g.split('.')[2][0:2])+float(g.split('.')[2][2:])/60.0 for g in fgeo])
    geosy['{}'.format(ddy)] = {'time':deepcopy(geo_time),
                               'fname':deepcopy(fgeo),
                               'path':fp+'data_other/MODIS/MYD03/2016/{}/'.format(ddy),
                   'lat':np.zeros((len(geo_time),2030,1354))+np.nan,'lon':np.zeros((len(geo_time),2030,1354))+np.nan}
    for ngeo,fg in enumerate(fgeo):
        geo_tmp,geo_tmp_dict = lu.load_hdf(geosy['{}'.format(ddy)]['path']+fg,values=(('lat',12),('lon',13)),verbose=False)
        geosy['{}'.format(ddy)]['lat'][ngeo,:,:] = geo_tmp['lat'][:2030,:1354]
        geosy['{}'.format(ddy)]['lon'][ngeo,:,:] = geo_tmp['lon'][:2030,:1354]
    # load the MYD06 file pointers
    fmods = os.listdir(fp+'data_other/MODIS/MYD06_L2/2016/{}/'.format(ddy))
    fmods.sort()
    fh = []
    for ff in fmods:
        fh1 = SD(str(fp+'data_other/MODIS/MYD06_L2/2016/{}/'.format(ddy)+ff),SDC.READ)
        fh.append(fh1)
    geosy['{}'.format(ddy)]['fh'] = fh[:]


# In[79]:


fmod


# In[80]:


geoso['264']['fname']


# In[93]:


def read_mod06_hdf_pixel(fh,ix,iy):
    re = fh.select('Cloud_Effective_Radius')
    ref = re[ix][iy]
    if ref == re.attributes()['_FillValue']: ref = np.nan
    ref = ref*re.attributes()['scale_factor']
    
    co = fh.select('Cloud_Optical_Thickness')
    cod = co[ix][iy]
    if cod == co.attributes()['_FillValue']: cod = np.nan
    cod = cod*co.attributes()['scale_factor']
    
    #ma = fh.select('Cloud_Mask_1km')
    #mask = ma[ix][iy]
    
    qal = fh.select('Quality_Assurance_1km')
    qa = qal[ix][iy]
    
    return cod,ref,qa


# ### Run through the flight data and build colocations

# In[94]:


import map_utils as mu


# In[95]:


def read_mod06(fh):
    re = fh.select('Cloud_Effective_Radius')
    ref = re.get()
    ref = ref*re.attributes()['scale_factor']
    ref[ref<0] = np.nan
    
    co = fh.select('Cloud_Optical_Thickness')
    codm = co.get()
    codm = codm*co.attributes()['scale_factor']
    codm[codm<0] = np.nan
        
    qal = fh.select('Quality_Assurance_1km')
    qa = mod_qa_arr(qal.get())
    
    #cod[~qa] = np.nan
    #ref[~qa] = np.nan
    
    return codm,ref


# In[96]:


doys


# In[97]:


iaes = np.array([it for it,ta in enumerate(ar['time_ae'][fla]) if any(abs(ta - ae['time']/3600)<1.0)])
len(iaes)


# In[100]:


nval = len(ar['Start_UTC'][fla])
mod = {'time':np.zeros((nval))+np.nan,'lat':np.zeros((nval))+np.nan,'lon':np.zeros((nval))+np.nan,
       'cod':np.zeros((nval))+np.nan,
       'ref':np.zeros((nval))+np.nan,'qa':np.zeros((nval))-1,'mask':np.zeros((nval))-1}
myd = {'time':np.zeros((nval))+np.nan,'lat':np.zeros((nval))+np.nan,'lon':np.zeros((nval))+np.nan,
       'cod':np.zeros((nval))+np.nan,
       'ref':np.zeros((nval))+np.nan,'qa':np.zeros((nval))-1,'mask':np.zeros((nval))-1}


# In[101]:


dos = np.unique(ar['doy'][fla][iaes])
for ddd in dos:
    ad = '{}'.format(ddd)
    print ad
    idd = ar['doy'][fla][iaes]==ddd
    for igo,go in enumerate(geoso[ad]['fh']):
        ind = mu.map_ind(geoso[ad]['lat'][igo,:,:],geoso[ad]['lon'][igo,:,:],
                         ar['Latitude'][fla][iaes][idd],ar['Longitude'][fla][iaes][idd])
        if np.array(ind).any():
            codm,ref = read_mod06(go)
            cs = np.array([[np.nanmean(codm[ix-4:ix+4,iy-4:iy+4]),np.nanmean(ref[ix-4:ix+4,iy-4:iy+4])] for ix,iy in ind.T])
            mod['cod'][iaes[idd]] = cs[:,0]
            mod['ref'][iaes[idd]] = cs[:,1]
            mod['time'][iaes[idd]] = geoso[ad]['time'][igo]
            mod['lat'][iaes[idd]] = ar['Latitude'][fla][iaes][idd]
            mod['lon'][iaes[idd]] = ar['Longitude'][fla][iaes][idd]
            
    for igy,gy in enumerate(geosy[ad]['fh']):
        ind = mu.map_ind(geosy[ad]['lat'][igy,:,:],geosy[ad]['lon'][igy,:,:],
                         ar['Latitude'][fla][iaes][idd],ar['Longitude'][fla][iaes][idd])
        if np.array(ind).any():
            cody,refy = read_mod06(gy)
            cs = np.array([[np.nanmean(cody[ix-4:ix+4,iy-4:iy+4]),np.nanmean(refy[ix-4:ix+4,iy-4:iy+4])] for ix,iy in ind.T])
            myd['cod'][iaes[idd]] = cs[:,0]
            myd['ref'][iaes[idd]] = cs[:,1]
            myd['time'][iaes[idd]] = geosy[ad]['time'][igy]
            myd['lat'][iaes[idd]] = ar['Latitude'][fla][iaes][idd]
            myd['lon'][iaes[idd]] = ar['Longitude'][fla][iaes][idd]


# In[102]:


sio.savemat(fp+'data_other/MODIS/MODIS_ORACLES2016_match.mat',{'mod':mod,'myd':myd})


# In[273]:


plt.figure()
plt.plot(mod['cod'])
plt.plot(myd['cod'])


# ### Get the fit to sine function
# Following Min & Zhang,2014 function for cloud fraction diurnal cycle 
# 
# $COD(t) = A sin(\frac{\pi(t+\phi)}{12})+B$
# 
# With local max at 6:50 am (for cloud fraction but applied to COD), giving $\phi=-0.96$
# 
# For two points, we get to solve for A and B:  
# 
# $A = \frac{COD_2(t_2)-COD_1(t_1)}{sin(\frac{\pi(t_2+\phi)}{12})-sin(\frac{\pi(t_1+\phi)}{12})}$
#   
# and 
# 
# $B = COD_1(t_1) - \frac{COD_2(t_2)-COD_1(t_1)}{sin(\frac{\pi(t_2+\phi)}{12})-sin(\frac{\pi(t_1+\phi)}{12})} sin(\frac{\pi(t_1+\phi)}{12})$
#   
# Do the same for ref 
# 

# In[57]:


def fit_sine(f1,f2,t1,t2):
    'get the fitted values for the sine wave f1 anf f2 are either cod or ref, t1 and t2 are the times in hours'
    phi = -0.96
    A = (f2-f1) / (np.sin(np.pi*(t2+phi)/12.0) - np.sin(np.pi*(t1+phi)/12.0))
    B = f1 - A*np.sin(np.pi*(t1+phi)/12.0)
    return A,B


# In[58]:


def get_sines(A,B,f0,t0,utcs=np.arange(0,24.0,0.5),val_range=[0,100]):
    phi = -0.96
    fbs = A*np.sin(np.pi*(utcs+phi)/12.0)+B
    fb0 = A*np.sin(np.pi*(t0+phi)/12.0)+B
    fb_frac = fb0/f0
    new_f = fbs/fb_frac
    if any(new_f<0):
        new_f = new_f+abs(np.nanmin(new_f))
    if np.nanmax(new_f)>val_range[1]:
        frac_fit = val_range[1]/np.nanmax(new_f)
        #print 'fract_fit: ',frac_fit
        new_f = get_sines(A*frac_fit,B+frac_fit,f0,t0,utcs=utcs,val_range=val_range)
    return new_f


# Prep for sine fitting (fill out missing values with interp)

# In[106]:


mod_cod = np.arange(0,len(mod['time'][iaes]))
myd_cod = np.arange(0,len(myd['time'][iaes]))
fmodcod = interp1d(mod_cod[np.isfinite(mod['cod'][iaes])],mod['cod'][iaes[np.isfinite(mod['cod'][iaes])]],bounds_error=False,kind='linear',fill_value='extrapolate')
fmydcod = interp1d(myd_cod[np.isfinite(myd['cod'][iaes])],myd['cod'][iaes[np.isfinite(myd['cod'][iaes])]],bounds_error=False,kind='linear',fill_value="extrapolate")
mod['codb'] = fmodcod(mod_cod)
myd['codb'] = fmydcod(myd_cod)


# In[107]:


mod_ref = np.arange(0,len(mod['time'][iaes]))
myd_ref = np.arange(0,len(myd['time'][iaes]))
fmodref = interp1d(mod_ref[np.isfinite(mod['ref'][iaes])],mod['ref'][iaes[np.isfinite(mod['ref'][iaes])]],bounds_error=False,kind='linear',fill_value='extrapolate')
fmydref = interp1d(myd_ref[np.isfinite(myd['ref'][iaes])],myd['ref'][iaes[np.isfinite(myd['ref'][iaes])]],bounds_error=False,kind='linear',fill_value="extrapolate")
mod['refb'] = fmodref(mod_ref)
myd['refb'] = fmydref(myd_ref)


# In[108]:


A_cod,B_cod = fit_sine(mod['codb'],myd['codb'],mod['time'][iaes],myd['time'][iaes])
A_ref,B_ref = fit_sine(mod['refb'],myd['refb'],mod['time'][iaes],myd['time'][iaes])


# #### Save

# In[113]:


diurn = {'A_cod':A_cod,'B_cod':B_cod,'A_ref':A_ref,'B_ref':B_ref,'utcs':np.arange(0,24,0.5),'iaes':iaes}


# In[114]:


sio.savemat(fp+'data_other/MODIS/MODIS_ORACLES2016_match.mat',{'mod':mod,'myd':myd,'diurn':diurn})


# #### Load

# In[42]:


ddiurn = sio.loadmat(fp+'data_other/MODIS/MODIS_ORACLES2016_match.mat')


# In[43]:


ddiurn.keys()


# In[44]:


def revert_from_mat(m):
    for k in m.dtype.names:
        if len(m[k])==1:
            m[k] = m[k][0,0]
    return m


# In[45]:


mod = ddiurn['mod']
myd = ddiurn['myd']
diurn = ddiurn['diurn']


# In[112]:


mod = revert_from_mat(ddiurn['mod'])
myd = revert_from_mat(ddiurn['myd'])
diurn = revert_from_mat(ddiurn['diurn'])


# In[46]:


A_cod, B_cod, A_ref, B_ref, iaes = diurn['A_cod'][0,0].flatten(),diurn['B_cod'][0,0].flatten(),                                    diurn['A_ref'][0,0].flatten(),diurn['B_ref'][0,0].flatten(),diurn['iaes'][0,0].flatten()


# In[47]:


A_cod


# In[48]:


len(A_cod),len(mod['lat']),len(iaes)


# In[49]:


len(ddiurn['mod'])


# In[50]:


mod.dtype.names


# #### plot

# In[51]:


plt.figure()
plt.plot(A_cod,label='A')
plt.plot(B_cod,label='B')
plt.title('COD')
plt.legend()
plt.figure()
plt.plot(A_ref,label='A')
plt.plot(B_ref,label='B')
plt.legend()
plt.title('ref')


#  Fill out the values to not have NaN

# In[52]:


x_cod = np.arange(0,len(A_cod))
fAcod = interp1d(x_cod[np.isfinite(A_cod)],A_cod[np.isfinite(A_cod)],bounds_error=False,kind='linear',fill_value='extrapolate')
fBcod = interp1d(x_cod[np.isfinite(A_cod)],B_cod[np.isfinite(A_cod)],bounds_error=False,kind='linear',fill_value="extrapolate")
A_codb = fAcod(x_cod)
B_codb = fBcod(x_cod)


# In[53]:


x_ref = np.arange(0,len(A_ref))
fAref = interp1d(x_ref[np.isfinite(A_ref)],A_ref[np.isfinite(A_ref)],bounds_error=False,kind='linear',fill_value='extrapolate')
fBref = interp1d(x_ref[np.isfinite(A_ref)],B_ref[np.isfinite(A_ref)],bounds_error=False,kind='linear',fill_value="extrapolate")
A_refb = fAref(x_ref)
B_refb = fBref(x_ref)


# In[54]:


if not 'fla' in locals():
    # for all 4STAR aerosol arrays
    fla = (ar['flag_acaod']==1) & ar['fl'] & ar['fl_QA'] & (ar['days']>2.0) 


# In[55]:


if not 'flb' in locals():
    # for the cod and ref arrays
    fld = (ar['days']>2.0) & (ar['days']!=17.0) 
    flb = (ar['flag_acaod'][fld]==1) & ar['fl'][fld] & ar['fl_QA'][fld]


# In[59]:


# test
i = 0
print cod[flb][iaes][i],ref[flb][iaes][i],ar['Start_UTC'][fla][iaes][i],A_cod[i],B_cod[i]
get_sines(A_cod[i],B_cod[i],cod[flb][iaes][i],ar['Start_UTC'][fla][iaes][i])


# In[60]:


get_sines(A_ref[i],B_ref[i],ref[flb][iaes][i],ar['Start_UTC'][fla][iaes][i])


# ### Run through the modis match multidimensional (old)

# In[ ]:


if False:
    if isjupyter():
        pbar = tqdm(total=nval)
    for i,u in enumerate(ar['Start_UTC'][fla]):
        if i<22760: 
            pbar.update(1)
            continue
        iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))
        if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 
            ad = '{}'.format(ar['doy'][fla][i])
            for igo,go in enumerate(geoso[ad]['fh']):
                ix,iy = np.unravel_index(np.argmin(abs(geoso[ad]['lat'][igo,:,:]-ar['Latitude'][fla][i])+
                            abs(geoso[ad]['lon'][igo,:,:]-ar['Longitude'][fla][i])),geoso[ad]['lon'][igo,:,:].shape)
                if (abs(geoso[ad]['lat'][igo,ix,iy]-ar['Latitude'][fla][i])<0.1) &                   (abs(geoso[ad]['lon'][igo,ix,iy]-ar['Longitude'][fla][i])<0.1):
                    try:
                        cod,ref,qa = read_mod06_hdf_pixel(go,ix,iy)
                        if mod_qa(qa): 
                            mod['cod'][i] = cod
                            mod['ref'][i] = ref
                            mod['time'][i] = geoso[ad]['time'][igo]
                    except:
                        pass
            for igy,gy in enumerate(geosy[ad]['fh']):
                ix,iy = np.unravel_index(np.argmin(abs(geosy[ad]['lat'][igy,:,:]-ar['Latitude'][fla][i])+
                            abs(geosy[ad]['lon'][igy,:,:]-ar['Longitude'][fla][i])),geosy[ad]['lon'][igy,:,:].shape)
                if (abs(geosy[ad]['lat'][igy,ix,iy]-ar['Latitude'][fla][i])<0.1) &                   (abs(geosy[ad]['lon'][igy,ix,iy]-ar['Longitude'][fla][i])<0.1):
                    try:
                        cody,refy,qay = read_mod06_hdf_pixel(gy,ix,iy)
                        if mod_qa(qay):
                            myd['cod'][i] = cody
                            myd['ref'][i] = refy
                            myd['time'][i] = geosy[ad]['time'][igy]
                    except:
                        pass
        pbar.update(1)


# In[68]:


class KeyboardInterruptError(Exception): pass


# In[574]:


def get_modis_val(i):
    m = {'ocod':np.nan,'oref':np.nan,'omask':np.nan,'otime':np.nan,
         'ycod':np.nan,'yref':np.nan,'ymask':np.nan,'ytime':np.nan}
    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 
        ad = '{}'.format(ar['doy'][fla][i])
        for igo,go in enumerate(geoso[ad]['fh']):
            ix,iy = np.unravel_index(np.argmin(abs(geoso[ad]['lat'][igo,:,:]-ar['Latitude'][fla][i])+
                        abs(geoso[ad]['lon'][igo,:,:]-ar['Longitude'][fla][i])),geoso[ad]['lon'][igo,:,:].shape)
            if (abs(geoso[ad]['lat'][igo,ix,iy]-ar['Latitude'][fla][i])<0.02) &               (abs(geoso[ad]['lon'][igo,ix,iy]-ar['Longitude'][fla][i])<0.02):
                cod,ref,qa = read_mod06_hdf_pixel(go,ix,iy)
                if mod_qa(qa):
                    try: 
                        m['ocod'] = cod
                        m['oref'] = ref
                    except:
                        pass
                    m['otime'] = geoso[ad]['time'][igo]
        for igy,gy in enumerate(geosy[ad]['fh']):
            ix,iy = np.unravel_index(np.argmin(abs(geosy[ad]['lat'][igy,:,:]-ar['Latitude'][fla][i])+
                        abs(geosy[ad]['lon'][igy,:,:]-ar['Longitude'][fla][i])),geosy[ad]['lon'][igy,:,:].shape)
            if (abs(geosy[ad]['lat'][igy,ix,iy]-ar['Latitude'][fla][i])<0.02) &               (abs(geosy[ad]['lon'][igy,ix,iy]-ar['Longitude'][fla][i])<0.02):
                cody,refy,qay = read_mod06_hdf_pixel(gy,ix,iy)
                if mod_qa(qay):
                    try:
                        m['ycod'] = cody
                        m['yref'] = refy
                    except:
                        pass
                    m['ytime'] = geosy[ad]['time'][igy]


# In[575]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[ ]:


p = Pool(cpu_count()-1,worker_init)


# In[ ]:


mod_outputs = []
max_ = len(ar['Start_UTC'][fla])
with tqdm(total=max_) as pbar:
    for i, outs in tqdm(enumerate(p.imap_unordered(get_modis_val, range(0, max_)))):
        pbar.update()
        mod_outputs.append(outs)


# In[318]:


cod[flb]


# # Prepare the base dict and defaults

# In[61]:


from datetime import datetime
datetime(2015,11,17).timetuple().tm_yday


# In[62]:


# for all 4STAR aerosol arrays
fla = (ar['flag_acaod']==1) & ar['fl'] & ar['fl_QA'] & (ar['days']>2.0) 


# In[63]:


# for the cod and ref arrays
fld = (ar['days']>2.0) & (ar['days']!=17.0) 
flb = (ar['flag_acaod'][fld]==1) & ar['fl'][fld] & ar['fl_QA'][fld]


# In[64]:


len(ar['AOD0355'][fla])


# In[65]:


len(cod[flb])


# In[66]:


sum(np.isfinite(cod[~flb])),sum(np.isfinite(cod[flb])),len(cod[flb])


# In[67]:


ka = ar.keys()
ka.sort()
ka


# In[68]:


doy = datetime.strptime(dds[int(ar['days'][fla][0])],'%Y%m%d').timetuple().tm_yday


# In[69]:


doy


# In[70]:


geo = {'lat':ar['Latitude'][0],'lon':ar['Longitude'][0],'doy':doy,'zout':[0,1.5,100.0]}
aero_no = {} # none
cloud = {'ztop':1.0,'zbot':0.5,'write_moments_file':False}
source = {'wvl_range':[201.0,4900.0],'source':'solar','integrate_values':True,'run_fuliou':True,
          'dat_path':fp_uvspec_dat}
albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':5.0}


# In[71]:


cloud['phase'] = 'wc'
geo['sza'] = 40.0
cloud['tau'] = 2.0
cloud['ref'] = 5.0


# In[72]:


pmom = Rl.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',deltascale=True,new=True)
tw = ['ntheta','rho','ssa','nmom','ext']
th = ['phase','pmom','theta']
on = ['wvl','nim','nre']
pmom_new = {}
for k in pmom.keys():
    if k in tw:
        pmom_new[k] = pmom[k][:,::50]
    elif k in th:
        pmom_new[k] = pmom[k][:,::50,:]
    elif k in on:
        pmom_new[k] = pmom[k][::50]
    else:
        pmom_new[k] = pmom[k]
    
    try:
        print k, pmom_new[k].shape        
    except:
        pass
pmom_new['wvl'][0] = 0.250
pmom_new['wvl'][-1] = 4.900
cloud['moms_dict'] = pmom_new


# In[73]:


pmom_new['wvl']


# In[74]:


wvl = np.append(np.append([250.0],ae['wavelength']),4900.0)
wvl


# In[75]:


aero = {'expand_hg':True,'disort_phase':False,'z_arr':[2.0,5.0],
        'wvl_arr':wvl}


# In[76]:


def fx_aero(aprop):
    'Function the aerosol property a 2d matrix for height and spectra, and extend the wavelength from 250 to 4900 nm'
    atmp = np.append([aprop[0]],np.append(aprop,aprop[-1]))
    return np.array([atmp,atmp])


# In[77]:


def fx_ext(a0,a1,a2,wvl=wvl):
    'Function to create the extinction coefficients from 4STAR AODs'
    aod = np.exp(np.polyval([a2,a1,a0],np.log(wvl)))
    aod[-1] = 0.0 # set the last wavelength to zero
    return np.array([aod/3.0,aod*0.0])


# In[78]:


aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][0],ar['AOD_polycoef_a1'][fla][0],ar['AOD_polycoef_a2'][fla][0])


# In[79]:


aero['asy'] = fx_aero(ae['g_total'][0])


# In[80]:


aero['ssa'] = fx_aero(ae['SSA'][0])


# ## Prepare the file list and saving

# In[81]:


from path_utils import isjupyter


# In[82]:


isjupyter()


# ## Write the files

# ### Conventional

# In[80]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[81]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[120]:


# for writing out the files


# In[109]:


ae['time']/3600.0


# In[110]:


ar['time_ae'][fla]


# In[82]:


if isjupyter():
    pbar = tqdm(total=len(ar['Start_UTC'][fla]))
for i,u in enumerate(ar['Start_UTC'][fla]):
    
    f_in = '{name}_{vv}_DARE_{i:03d}_withaero.dat'.format(name=name,vv=vv,i=i)

    geo['lat'],geo['lon'],geo['sza'] = ar['Latitude'][fla][i],ar['Longitude'][fla][i],sza[fla][i]
    day = days[ar['days'][fla][i].astype(int)]
    geo['doy'] = datetime(int(day[0:4]),int(day[4:6]),int(day[6:])).timetuple().tm_yday

    cloud['tau'],cloud['ref'] = cod[flb][i],ref[flb][i]
    cloud['write_moments_file'] = True

    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))

    # Only run for aerosol rertievals within 1 hour
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 

        aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][i],ar['AOD_polycoef_a1'][fla][i],ar['AOD_polycoef_a2'][fla][i])
        aero['ssa'] = fx_aero(ae['SSA'][iae])
        aero['asy'] = fx_aero(ae['g_total'][iae])

        Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in,out=fpp_out+f_in))

        f_in = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vv,i=i)
        Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in,out=fpp_out+f_in))

    if isjupyter(): 
        pbar.update(1)
    else:
        print i

f.close()


# ### Multiprocessing

# In[146]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[60]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[148]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[149]:


if isjupyter():
    pbar = tqdm(total=len(ar['Start_UTC'][fla]))
bb = []
for i,u in enumerate(ar['Start_UTC'][fla]):
    
    f_in = '{name}_{vv}_DARE_{i:03d}_withaero.dat'.format(name=name,vv=vv,i=i)

    geo['lat'],geo['lon'],geo['sza'] = ar['Latitude'][fla][i],ar['Longitude'][fla][i],sza[fla][i]
    day = days[ar['days'][fla][i].astype(int)]
    geo['doy'] = datetime(int(day[0:4]),int(day[4:6]),int(day[6:])).timetuple().tm_yday

    if ~np.isfinite(cod[flb][i]):
        if isjupyter():
            pbar.update(1)
        continue
    cloud['tau'],cloud['ref'] = cod[flb][i],ref[flb][i]
    cloud['write_moments_file'] = True

    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))

    # Only run for aerosol rertievals within 1 hour
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 

        aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][i],ar['AOD_polycoef_a1'][fla][i],ar['AOD_polycoef_a2'][fla][i])
        aero['ssa'] = fx_aero(ae['SSA'][iae])
        aero['asy'] = fx_aero(ae['g_total'][iae])

        #Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
        #                           verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in,out=fpp_out+f_in))

        f_in_noa = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vv,i=i)
        #Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
        #                           verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in_noa,out=fpp_out+f_in_noa))
        
        bb.append({'geo':deepcopy(geo),'cod':cod[flb][i],'ref':ref[flb][i],'aero':deepcopy(aero),
                   'f_in':deepcopy(f_in),'f_in_noa':deepcopy(f_in_noa)})

    if isjupyter(): 
        pbar.update(1)
    else:
        print i

f.close()


# In[150]:


def write_files(d,cloud=cloud,source=source,albedo=albedo,aero_no=aero_no):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    Rl.write_input_aac(fpp_in+d['f_in'],geo=d['geo'],aero=d['aero'],cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
    Rl.write_input_aac(fpp_in+d['f_in_noa'],geo=d['geo'],aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)


# In[151]:


p = Pool(cpu_count()-1,worker_init)


# In[152]:


len(bb)


# In[153]:


results = []
max_ = len(bb)
with tqdm(total=max_) as pbar:
    for i, res in tqdm(enumerate(p.imap_unordered(write_files, bb))):
        pbar.update()
        results.append(res)


# In[257]:


results = p.map(write_files,bb)


# ### Diurnal averaging with multiprocessing

# In[171]:


if vv=='v3':
    vvh = 'v3_24h'
    dssa = 0.0
elif vv =='v4':
    vvh = 'v4_24h'
    dssa = 0.0
elif vv == 'v5':
    vvh = 'v5_aodm'
print vv,vvh


# In[175]:


dssa = 0.0
dasy = 0.0
daod = 0.0


if vvh == 'v5_ssam': dssa = -0.03
elif vvh == 'v5_ssap': dssa = 0.03
elif vvh == 'v5_asym': dasy = -0.03
elif vvh == 'v5_asyp': dasy = 0.03
elif vvh == 'v5_aodm': daod = -0.03
elif vvh == 'v5_aodp': daod = 0.03
daod,dssa,dasy


# In[176]:


vvh = vvh+'_diurnal' # for the version with diurnal cycles from MODIS


# In[177]:


print vvh


# In[178]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[179]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vvh),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vvh)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vvh)


# In[180]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[181]:


if isjupyter():
    pbar = tqdm(total=len(ar['Start_UTC'][fla]))
bb = []
ii = 0 
for i,u in enumerate(ar['Start_UTC'][fla]):
    
    f_in = '{name}_{vv}_DARE_{i:03d}_withaero.dat'.format(name=name,vv=vvh,i=i)

    geo['lat'],geo['lon'],geo['sza'] = ar['Latitude'][fla][i],ar['Longitude'][fla][i],sza[fla][i]
    day = days[ar['days'][fla][i].astype(int)]
    geo['doy'] = datetime(int(day[0:4]),int(day[4:6]),int(day[6:])).timetuple().tm_yday

    if ~np.isfinite(cod[flb][i]):
        if isjupyter():
            pbar.update(1)
        continue
    cloud['tau'],cloud['ref'] = cod[flb][i],ref[flb][i]
    cloud['write_moments_file'] = True

    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))

    # Only run for aerosol rertievals within 1 hour
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 

        aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][i],ar['AOD_polycoef_a1'][fla][i],ar['AOD_polycoef_a2'][fla][i])+daod
        try: aero['ext'][aero['ext']<0.0] = 0.0
        except: pass
        aero['ssa'] = fx_aero(ae['SSA'][iae])+dssa
        try: aero['ssa'][aero['ssa']>1.0] = 1.0
        except: pass
        try: aero['ssa'][aero['ssa']<0.0] = 0.0
        except: pass
        aero['asy'] = fx_aero(ae['g_total'][iae])+dasy
        try: aero['asy'][aero['asy']>1.0] = 1.0
        except: pass
        try: aero['asy'][aero['asy']<0.0] = 0.0
        except: pass

        #Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
        #                           verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in,out=fpp_out+f_in))

        f_in_noa = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vvh,i=i)
        #Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
        #                           verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in_noa,out=fpp_out+f_in_noa))
        
        bb.append({'geo':deepcopy(geo),'cod':cod[flb][i],'ref':ref[flb][i],'aero':deepcopy(aero),
                   'f_in':deepcopy(f_in),'f_in_noa':deepcopy(f_in_noa),
                   'A_cod':A_cod[ii],'B_cod':B_cod[ii],'A_ref':A_ref[ii],'B_ref':B_ref[ii],'UTC':ar['Start_UTC'][fla][i]})
        ii = ii+1

    if isjupyter(): 
        pbar.update(1)
    else:
        print i

f.close()


# In[182]:


# Go through and write out the list file for 24 h
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vvh),'w')
if isjupyter():
    pbar = tqdm(total=len(ar['Start_UTC'][fla]))
errlog = fp_rtm+'{}_DARE_{}_error.log'.format(name,vvh)
for i,u in enumerate(ar['Start_UTC'][fla]):
    if ~np.isfinite(cod[flb][i]):
        if isjupyter():
            pbar.update(1)
        continue
    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 
        for ux in np.arange(0,24,0.5):
            f_in = '{name}_{vv}_DARE_{i:03d}_withaero_{ux:04.1f}.dat'.format(name=name,vv=vvh,i=i,ux=ux)
            f_in_noa = '{name}_{vv}_star_{i:03d}_noaero_{ux:04.1f}.dat'.format(name=name,vv=vvh,i=i,ux=ux)
            f.write('{uv} < {fin} > {out} 2>> {err}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in,out=fpp_out+f_in,err=errlog))
            f.write('{uv} < {fin} > {out} 2>> {err}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in_noa,out=fpp_out+f_in_noa,err=errlog))
    pbar.update(1)  
f.close()


# In[183]:


fp_rtm+'{}_DARE_{}.sh'.format(name,vvh)


# In[184]:


def write_files(d,cloud=cloud,source=source,albedo=albedo,aero_no=aero_no):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    d['f_in'] = d['f_in'].replace('.dat','')+'_{ux:04.1f}.dat'
    d['f_in_noa'] = d['f_in_noa'].replace('.dat','')+'_{ux:04.1f}.dat'
    for ux in np.arange(0,24,0.5):
        d['geo']['utc'] = ux
        d['geo']['hour'] = int(ux)
        d['geo']['minute'] = int((ux-int(ux))*60.0)
        d['geo']['second'] = 0
        d['geo'].pop('sza',None)
        Rl.write_input_aac(fpp_in+d['f_in'].format(ux=ux),geo=d['geo'],aero=d['aero'],cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
        Rl.write_input_aac(fpp_in+d['f_in_noa'].format(ux=ux),geo=d['geo'],aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)


# In[185]:


from datetime import datetime, timedelta
ud = datetime(2016,1,1) + timedelta(bb[0]['geo']['doy'])
bb[0]['geo']['doy']


# In[186]:


timedelta()


# In[187]:


ud.month


# In[188]:


len(bb)


# In[98]:


def write_files(d,cloud=cloud,source=source,albedo=albedo,aero_no=aero_no):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    d['f_in'] = d['f_in'].replace('.dat','')+'_{ux:04.1f}.dat'
    d['f_in_noa'] = d['f_in_noa'].replace('.dat','')+'_{ux:04.1f}.dat'
    ud = datetime.datetime(2016,1,1) + datetime.timedelta(d['geo']['doy'])
    d['geo']['day'] = ud.day
    d['geo']['year'] = ud.year
    d['geo']['month'] = ud.month
    for ux in np.arange(0,24,0.5):
        d['geo']['utc'] = ux
        d['geo']['hour'] = int(ux)
        d['geo']['minute'] = int((ux-int(ux))*60.0)
        d['geo']['second'] = 0
        d['geo'].pop('sza',None)
        Rl.write_input_aac(fpp_in+d['f_in'].format(ux=ux),geo=d['geo'],aero=d['aero'],cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
        Rl.write_input_aac(fpp_in+d['f_in_noa'].format(ux=ux),geo=d['geo'],aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)


# In[99]:


def write_files_sub(d,cloud=cloud,source=source,albedo=albedo,aero_no=aero_no):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    d['f_in'] = d['f_in'].replace('.dat','')+'_{ux:04.1f}.dat'
    d['f_in_noa'] = d['f_in_noa'].replace('.dat','')+'_{ux:04.1f}.dat'
    ud = datetime(2016,1,1) + timedelta(d['geo']['doy'])
    d['geo']['day'] = ud.day
    d['geo']['year'] = ud.year
    d['geo']['month'] = ud.month
    source['run_fuliou'] = True
    d['aero']['link_to_mom_file'] = False
    d['aero']['file_name'] = fpp_in+d['f_in'].format(ux=0.0)+'_aero'
    cloud['link_to_mom_file'] = False
    cloud['file_name'] = fpp_in+d['f_in'].format(ux=0.0)+'_cloud'
    d['geo'].pop('doy',None)
    for ux in np.arange(0,24,0.5):
        d['geo']['utc'] = ux
        d['geo']['hour'] = int(ux)
        d['geo']['minute'] = int((ux-int(ux))*60.0)
        d['geo']['second'] = 0
        d['geo'].pop('sza',None)
        Rl.write_input_aac(fpp_in+d['f_in'].format(ux=ux),geo=d['geo'],aero=d['aero'],cloud=cloud,
                           source=source,albedo=albedo,solver='twostr',
                           verbose=False,make_base=False,set_quiet=True)
        Rl.write_input_aac(fpp_in+d['f_in_noa'].format(ux=ux),geo=d['geo'],aero=aero_no,cloud=cloud,
                           source=source,albedo=albedo,solver='twostr',
                           verbose=False,make_base=False,set_quiet=True)
        cloud['link_to_mom_file'] = True
        d['aero']['link_to_mom_file'] = True


# In[189]:


def write_files_sub_diurn(d,cloud=cloud,source=source,albedo=albedo,aero_no=aero_no):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    d['f_in'] = d['f_in'].replace('.dat','')+'_{ux:04.1f}.dat'
    d['f_in_noa'] = d['f_in_noa'].replace('.dat','')+'_{ux:04.1f}.dat'
    ud = datetime(2016,1,1) + timedelta(d['geo']['doy'])
    d['geo']['day'] = ud.day
    d['geo']['year'] = ud.year
    d['geo']['month'] = ud.month
    source['run_fuliou'] = True
    d['aero']['link_to_mom_file'] = False
    d['aero']['file_name'] = fpp_in+d['f_in'].format(ux=0.0)+'_aero'
    cloud['link_to_mom_file'] = False
    cloud['file_name'] = fpp_in+d['f_in'].format(ux=0.0)+'_cloud'
    d['geo'].pop('doy',None)
    cods = get_sines(d['A_cod'],d['B_cod'],d['cod'],d['UTC'])
    refs = get_sines(d['A_ref'],d['B_ref'],d['ref'],d['UTC'])
    for iux,ux in enumerate(np.arange(0,24,0.5)):
        d['geo']['utc'] = ux
        d['geo']['hour'] = int(ux)
        d['geo']['minute'] = int((ux-int(ux))*60.0)
        d['geo']['second'] = 0
        d['geo'].pop('sza',None)
        cloud['tau'],cloud['ref'] = cods[iux],refs[iux]
        if ~np.isfinite(cloud['tau']): cloud['tau'] = d['cod']
        if ~np.isfinite(cloud['ref']): cloud['ref'] = d['ref']
        Rl.write_input_aac(fpp_in+d['f_in'].format(ux=ux),geo=d['geo'],aero=d['aero'],cloud=cloud,
                           source=source,albedo=albedo,solver='twostr',
                           verbose=False,make_base=False,set_quiet=True)
        Rl.write_input_aac(fpp_in+d['f_in_noa'].format(ux=ux),geo=d['geo'],aero=aero_no,cloud=cloud,
                           source=source,albedo=albedo,solver='twostr',
                           verbose=False,make_base=False,set_quiet=True)
        cloud['link_to_mom_file'] = False
        d['aero']['link_to_mom_file'] = True


# In[190]:


p = Pool(cpu_count()-1,worker_init)


# In[191]:


results = []
max_ = len(bb)
with tqdm(total=max_) as pbar:
    for i, res in tqdm(enumerate(p.imap_unordered(write_files_sub, bb))):
        pbar.update()
        results.append(res)


# In[192]:


results = []
max_ = len(bb)
with tqdm(total=max_) as pbar:
    for i, res in tqdm(enumerate(p.imap_unordered(write_files_sub_diurn, bb))):
        pbar.update()
        results.append(res)


# ### Test out using fifo instead of files

# In[119]:


import os
from multiprocessing import Process, Event
import subprocess as sub


# In[113]:


fp_fifo_in = '/tmp/uvspec_input.fifo'


# In[114]:


os.mkfifo(fp_fifo_in)


# In[118]:


p = open(fp_fifo_in,'w')
p.write('quiet\n')
p.write('mol_abs_param fu\n')
p.write('rte_solver twostr\n')
p.write('output_process sum\n')
p.write('data_files_path /home/sam/libradtran/libRadtran-2.0.2b/data/ \n')
p.write('source solar /home/sam/libradtran/libRadtran-2.0.2b/data/solar_flux/kurudz_1.0nm.dat per_nm\n')
p.write('wavelength 350.000000 4000.000000\n')
p.write('zout 0 3 100\n')
p.write('latitude N 56.938700\n')
p.write('longitude W 111.866900\n')
p.write('time 2018 06 09 15 30 00\n')
p.write('aerosol_default\n')
p.write('aerosol_modify ssa scale 0.85\n')
p.write('disort_intcor moments\n')
p.write('albedo 0.33\n')
p.close()


# In[154]:


def print_fifo():
    p = open(fp_fifo_in,'w')
    p.write('quiet\n')
    p.write('mol_abs_param fu\n')
    p.write('rte_solver twostr\n')
    p.write('output_process sum\n')
    p.write('data_files_path /home/sam/libradtran/libRadtran-2.0.2b/data/ \n')
    p.write('source solar /home/sam/libradtran/libRadtran-2.0.2b/data/solar_flux/kurudz_1.0nm.dat per_nm\n')
    p.write('wavelength 350.000000 4000.000000\n')
    p.write('zout 0 3 100\n')
    p.write('latitude N 56.938700\n')
    p.write('longitude W 111.866900\n')
    p.write('time 2018 06 09 15 30 00\n')
    p.write('aerosol_default\n')
    p.write('aerosol_modify ssa scale 0.85\n')
    p.write('disort_intcor moments\n')
    p.write('albedo 0.33\n')
    p.close()


# In[155]:


def run():
    process = sub.Popen([fp_uvspec],stdin=p, stdout=sub.PIPE,stderr=sub.PIPE)
    stdout,stderr = process.communicate()
    #stderr = process.stderr.read()
    print 'STDOUT:{},{},{}'.format(stdout,stderr,process.poll())


# In[175]:


p = open(fp_fifo_in,'w+')
print 'fifo: ',p
p.flush()
p.write('quiet\n')
p.write('mol_abs_param fu\n')
p.write('rte_solver twostr\n')
p.write('output_process sum\n')
p.write('data_files_path /home/sam/libradtran/libRadtran-2.0.2b/data/ \n')
p.write('source solar /home/sam/libradtran/libRadtran-2.0.2b/data/solar_flux/kurudz_1.0nm.dat per_nm\n')
p.write('wavelength 350.000000 4000.000000\n')
p.write('zout 0 3 100\n')
p.write('latitude N 56.938700\n')
p.write('longitude W 111.866900\n')
p.write('time 2018 06 09 15 30 00\n')
p.write('aerosol_default\n')
p.write('aerosol_modify ssa scale 0.85\n')
p.write('disort_intcor moments\n')
p.write('albedo 0.33\n')
#p.close()
process = sub.Popen([fp_uvspec],stdin=p,stdout=sub.PIPE,stderr=sub.PIPE)
stdout,stderr = process.communicate()
print 'STDOUT:{},{},{}'.format(stdout,stderr,process.poll())
p.close()


# In[216]:


def write_xtrf(fp_fifo_in):
    if not os.path.exists(fp_fifo_in):
        os.mkfifo(fp_fifo_in)
    p = open(fp_fifo_in,'w')
    p.flush()
    g = ['# wvl[nm]    alb[unitless','250.000000      -0.043171','350.000000      -0.010611','400.000000      0.005669','500.000000      0.038229','675.000000      0.058627','870.000000      0.229436','995.000000      0.234727','1200.000000     0.240584','1400.000000     0.246298','1600.000000     0.252013','2100.000000     0.266298','3200.000000     0.297727','4900.000000     0.346298']
    for llj in g:
        p.write('{}\n'.format(llj))
    #p.flush()
    p.close()
    os.unlink(fp_fifo_in)
    #os.remove(fp_fifo_in)


# In[209]:


write_xtrf(fp_fifo_in)


# In[218]:


r, w = os.pipe()
os.write(w,'verbose\n')
os.write(w,'mol_abs_param fu\n')
os.write(w,'rte_solver twostr\n')
os.write(w,'output_process sum\n')
os.write(w,'data_files_path /home/sam/libradtran/libRadtran-2.0.2b/data/ \n')
os.write(w,'source solar /home/sam/libradtran/libRadtran-2.0.2b/data/solar_flux/kurudz_1.0nm.dat per_nm\n')
os.write(w,'wavelength 350.000000 4000.000000\n')
os.write(w,'zout 0 3 100\n')
os.write(w,'latitude N 56.938700\n')
os.write(w,'longitude W 111.866900\n')
os.write(w,'time 2018 06 09 15 30 00\n')
os.write(w,'aerosol_default\n')
os.write(w,'aerosol_modify ssa scale 0.85\n')
os.write(w,'disort_intcor moments\n')
os.write(w,'albedo 0.33\n')
os.write(w,'albedo_file {}'.format(fp_fifo_in))
os.close(w)
#p.close()

p1 = Process(target=write_xtrf,args=(fp_fifo_in,))
p1.start()

process = sub.Popen([fp_uvspec],stdin=r,stdout=sub.PIPE,stderr=sub.PIPE)
stdout,stderr = process.communicate()
print 'STDOUT:{},{},{}'.format(stdout,stderr,process.poll())
#p.close()


# ### Run only the CRE portion in multiprocessing

# In[59]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[60]:


# open the list file
f = open(fp_rtm+'{}_DARE_CRE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_CRE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_CRE_{}/'.format(name,vv)


# In[61]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[63]:


if isjupyter():
    pbar = tqdm(total=len(ar['Start_UTC'][fla]))
bb = []
for i,u in enumerate(ar['Start_UTC'][fla]):
    
    f_in = '{name}_{vv}_DARE_CRE_{i:03d}_withaero_clear.dat'.format(name=name,vv=vv,i=i)

    geo['lat'],geo['lon'],geo['sza'] = ar['Latitude'][fla][i],ar['Longitude'][fla][i],sza[fla][i]
    day = days[ar['days'][fla][i].astype(int)]
    geo['doy'] = datetime(int(day[0:4]),int(day[4:6]),int(day[6:])).timetuple().tm_yday

    if ~np.isfinite(cod[flb][i]):
        if isjupyter():
            pbar.update(1)
        continue
    cloud['tau'],cloud['ref'] = 0.0,ref[flb][i]
    cloud['write_moments_file'] = True

    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))

    # Only run for aerosol rertievals within 1 hour
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 

        aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][i],ar['AOD_polycoef_a1'][fla][i],ar['AOD_polycoef_a2'][fla][i])
        aero['ssa'] = fx_aero(ae['SSA'][iae])
        aero['asy'] = fx_aero(ae['g_total'][iae])

        #Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
        #                           verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in,out=fpp_out+f_in))

     #   f_in_noa = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vv,i=i)
        #Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
        #                           verbose=False,make_base=False,set_quiet=True)
     #   f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in_noa,out=fpp_out+f_in_noa))
        
        bb.append({'geo':deepcopy(geo),'cod':cod[flb][i],'ref':ref[flb][i],'aero':deepcopy(aero),
                   'f_in':deepcopy(f_in)})

    if isjupyter(): 
        pbar.update(1)
    else:
        print i

f.close()


# In[64]:


def write_files_cre(d,cloud=cloud,source=source,albedo=albedo,aero_no=aero_no):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    Rl.write_input_aac(fpp_in+d['f_in'],geo=d['geo'],aero=d['aero'],cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
    #Rl.write_input_aac(fpp_in+d['f_in_noa'],geo=d['geo'],aero=aero_no,cloud=cloud,source=source,albedo=albedo,
    #                               verbose=False,make_base=False,set_quiet=True)


# In[67]:


p = Pool(cpu_count()-1,worker_init)


# In[68]:


results_cre = p.map(write_files_cre,bb)


# In[70]:


len(results_cre)


# ## Run the calculations

# Run the files from command line:
# 
# using command parallel --jobs=20 < ORACLES_DARE_v1.sh

# ### For Single point during each day

# In[172]:


f_list = fp_rtm+'{}_DARE_{}.sh'.format(name,vv)


# In[173]:


get_ipython().system(u' wc -l $f_list')


# In[174]:


f_listout = f_list+'.out'


# In[ ]:


get_ipython().system(u'parallel --jobs=22 --bar < $f_list #2> $f_listout')


# ### For daily averages

# In[193]:


f_list = fp_rtm+'{}_DARE_{}.sh'.format(name,vvh)


# In[194]:


get_ipython().system(u' wc -l $f_list')


# In[195]:


f_listout = f_list+'.out'


# In[ ]:


get_ipython().system(u'parallel --jobs=22 --bar < $f_list #2> $f_listout')


# ### For the CRE

# In[71]:


f_list = fp_rtm+'{}_DARE_CRE_{}.sh'.format(name,vv)


# In[72]:


get_ipython().system(u' wc -l $f_list')


# In[73]:


f_listout = f_list+'.out'


# In[74]:


get_ipython().system(u'parallel --jobs=7 --bar < $f_list 2> $f_listout')


# ## Read the files

# In[144]:


fpp_out,name,vv,geo['zout']


# In[145]:


n = len(ar['Start_UTC'][fla])
nz = len(geo['zout'])
nw = len(aero['wvl_arr'])
nt = 48
dat = {'cod':np.zeros(n)+np.nan,'ref':np.zeros(n)+np.nan,'ext':np.zeros((n,nw))+np.nan,
       'ssa':np.zeros((n,nw))+np.nan,'asy':np.zeros((n,nw))+np.nan,'zout':geo['zout'],
       'wvl':aero['wvl_arr'],'sza':np.zeros(n)+np.nan,
       'dn':np.zeros((n,nz))+np.nan,'up':np.zeros((n,nz))+np.nan,
       'dn_noa':np.zeros((n,nz))+np.nan,'up_noa':np.zeros((n,nz))+np.nan,
       'dn_avg':np.zeros((n,nz))+np.nan,'up_avg':np.zeros((n,nz))+np.nan,
       'dn_noa_avg':np.zeros((n,nz))+np.nan,'up_noa_avg':np.zeros((n,nz))+np.nan}

for i,u in enumerate(ar['Start_UTC'][fla]):
    
    dat['cod'][i] = cod[flb][i]
    dat['ref'][i] = ref[flb][i]
    dat['sza'][i] = sza[fla][i]

    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))
    # Only run for aerosol rertievals within 1 hour
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 

        dat['ext'][i,:] = fx_ext(ar['AOD_polycoef_a0'][fla][i],ar['AOD_polycoef_a1'][fla][i],ar['AOD_polycoef_a2'][fla][i])[0]
        dat['ssa'][i,:] = fx_aero(ae['SSA'][iae])[0]+dssa
        try: dat['ssa'][i,dat['ssa'][i,:]>1.0] = 1.0
        except: pass
        try: dat['ssa'][i,dat['ssa'][i,:]<0.0] = 0.0
        except: pass
        dat['asy'][i,:] = fx_aero(ae['g_total'][iae])[0]


# ### Conventional

# In[133]:



if isjupyter():
    pbar = tqdm(total=n)
    

for i,u in enumerate(ar['Start_UTC'][fla]):
    
    dat['cod'][i] = cod[flb][i]
    dat['ref'][i] = ref[flb][i]
    dat['sza'][i] = sza[fla][i]

    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))
    # Only run for aerosol rertievals within 1 hour
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 

        dat['ext'][i,:] = fx_ext(ar['AOD_polycoef_a0'][fla][i],ar['AOD_polycoef_a1'][fla][i],ar['AOD_polycoef_a2'][fla][i])[0]
        dat['ssa'][i,:] = fx_aero(ae['SSA'][iae])[0]
        dat['asy'][i,:] = fx_aero(ae['g_total'][iae])[0]
        try:
            f_in = '{name}_{vv}_DARE_{i:03d}_withaero.dat'.format(name=name,vv=vv,i=i)
            o = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])
            f_in = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vv,i=i)
            on = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])

            dat['dn'][i,:] = o['diffuse_down']+o['direct_down']
            dat['dn_noa'][i,:] = on['diffuse_down']+on['direct_down']
            dat['up'][i,:] = o['diffuse_up']
            dat['up_noa'][i,:] = on['diffuse_up']
        except:
            pass

    if isjupyter(): 
        pbar.update(1)
    else:
        print i


# ### Multiprocessing

# In[186]:


class KeyboardInterruptError(Exception): pass


# In[187]:


def read_files(i,fpp_out=fpp_out,name=name,vv=vv,zout=geo['zout']):
    'function to feed the pool of workers to read all the files'
    out = {}
    try:
        f_in = '{name}_{vv}_DARE_{i:03d}_withaero.dat'.format(name=name,vv=vv,i=i)
        o = Rl.read_libradtran(fpp_out+f_in,zout=zout)
        f_in = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vv,i=i)
        on = Rl.read_libradtran(fpp_out+f_in,zout=zout)

        #dat['dn'][i,:] = o['diffuse_down']+o['direct_down']
        #dat['dn_noa'][i,:] = on['diffuse_down']+on['direct_down']
        #dat['up'][i,:] = o['diffuse_up']
        #dat['up_noa'][i,:] = on['diffuse_up']
        
        out['dn'] = o['diffuse_down']+o['direct_down']
        out['dn_noa'] = on['diffuse_down']+on['direct_down']
        out['up'] = o['diffuse_up']
        out['up_noa'] = on['diffuse_up']
        out['i'] = i
    except KeyboardInterrupt:
        raise KeyboardInterruptError()
        
    except:
        out['dn'] = np.zeros(len(zout))+np.nan
        out['dn_noa'] = np.zeros(len(zout))+np.nan
        out['up'] = np.zeros(len(zout))+np.nan
        out['up_noa'] = np.zeros(len(zout))+np.nan
        out['i'] = i
    return out


# In[188]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[66]:


p = Pool(cpu_count()-1,worker_init)


# In[67]:


outputs = []
max_ = len(ar['Start_UTC'][fla])
with tqdm(total=max_) as pbar:
    for i, outs in tqdm(enumerate(p.imap_unordered(read_files, range(0, max_)))):
        pbar.update()
        outputs.append(outs)


# In[68]:


outputs[100]


# In[69]:


outputs[0],outputs[2000]


# In[70]:


dat.keys()


# In[71]:


for oo in outputs:
    dat['dn'][oo['i'],:] = oo['dn']
    dat['dn_noa'][oo['i'],:] = oo['dn_noa']
    dat['up'][oo['i'],:] = oo['up']
    dat['up_noa'][oo['i'],:] = oo['up_noa']


# In[72]:


outputs[0]['dn']


# In[73]:


np.nanmean(dat['dn'])


# ### Multiprocessing reading for diurnal averaging

# In[146]:


n = len(ar['Start_UTC'][fla])
nz = len(geo['zout'])
nw = len(aero['wvl_arr'])
nt = 48
dat_24h = {'cod':np.zeros(n)+np.nan,'ref':np.zeros(n)+np.nan,'ext':np.zeros((n,nw))+np.nan,
       'ssa':np.zeros((n,nw))+np.nan,'asy':np.zeros((n,nw))+np.nan,'zout':geo['zout'],
       'wvl':aero['wvl_arr'],'sza':np.zeros(n)+np.nan,
       'dn':np.zeros((n,nt,nz))+np.nan,'up':np.zeros((n,nt,nz))+np.nan,
       'dn_noa':np.zeros((n,nt,nz))+np.nan,'up_noa':np.zeros((n,nt,nz))+np.nan,
       'dn_avg':np.zeros((n,nz))+np.nan,'up_avg':np.zeros((n,nz))+np.nan,
       'dn_noa_avg':np.zeros((n,nz))+np.nan,'up_noa_avg':np.zeros((n,nz))+np.nan,
          'cods':np.zeros((n,nt))+np.nan,'refs':np.zeros((n,nt))+np.nan}
ii = 0 
for i,u in enumerate(ar['Start_UTC'][fla]):

    dat_24h['cod'][i] = cod[flb][i]
    dat_24h['ref'][i] = ref[flb][i]
    dat_24h['sza'][i] = sza[fla][i]


    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))
    # Only run for aerosol rertievals within 1 hour
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 

        dat_24h['ext'][i,:] = fx_ext(ar['AOD_polycoef_a0'][fla][i],ar['AOD_polycoef_a1'][fla][i],ar['AOD_polycoef_a2'][fla][i])[0]+daod
        try: dat_24h['ext'][i,dat_24h['ext'][i,:]<0.0] = 0.0
        except: pass
        #dat_24h['ssa'][i,:] = fx_aero(ae['SSA'][iae])[0]
        dat_24h['ssa'][i,:] = fx_aero(ae['SSA'][iae])[0]+dssa
        try: dat_24h['ssa'][i,dat_24h['ssa'][i,:]>1.0] = 1.0
        except: pass
        try: dat_24h['ssa'][i,dat_24h['ssa'][i,:]<0.0] = 0.0
        except: pass
        dat_24h['asy'][i,:] = fx_aero(ae['g_total'][iae])[0]+dasy
        try: dat_24h['asy'][i,dat_24h['asy'][i,:]>1.0] = 1.0
        except: pass
        try: dat_24h['asy'][i,dat_24h['asy'][i,:]<0.0] = 0.0
        except: pass
        
        dat_24h['cods'][i,:] = get_sines(A_cod[ii],B_cod[ii],cod[flb][i],u)
        dat_24h['refs'][i,:] = get_sines(A_ref[ii],B_ref[ii],ref[flb][i],u)
        ii = ii + 1


# In[147]:


class KeyboardInterruptError(Exception): pass


# In[148]:


v = np.zeros((48,3))+np.nan


# In[149]:


np.nanmean(v,axis=0)


# In[150]:


vvh


# In[151]:


if not 'fpp_out' in locals():
    fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vvh)
    fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vvh)


# In[152]:


def read_files(i,fpp_out=fpp_out,name=name,vv=vvh,zout=geo['zout']):
    'function to feed the pool of workers to read all the files'
    nz = len(zout)
    out = {'dn':np.zeros((48,nz))+np.nan,'dn_noa':np.zeros((48,nz))+np.nan,
           'up':np.zeros((48,nz))+np.nan,'up_noa':np.zeros((48,nz))+np.nan}
    out['i'] = i
    for u,ux in enumerate(np.arange(0,24,0.5)):
        f_in = '{name}_{vv}_DARE_{i:03d}_withaero_{ux:04.1f}.dat'.format(name=name,vv=vv,i=i,ux=ux)
        f_in_noa = '{name}_{vv}_star_{i:03d}_noaero_{ux:04.1f}.dat'.format(name=name,vv=vv,i=i,ux=ux)
        
        try:
            o = Rl.read_libradtran(fpp_out+f_in,zout=zout)
            on = Rl.read_libradtran(fpp_out+f_in_noa,zout=zout)

            #dat['dn'][i,:] = o['diffuse_down']+o['direct_down']
            #dat['dn_noa'][i,:] = on['diffuse_down']+on['direct_down']
            #dat['up'][i,:] = o['diffuse_up']
            #dat['up_noa'][i,:] = on['diffuse_up']

            out['dn'][u,:] = o['diffuse_down']+o['direct_down']
            out['dn_noa'][u,:] = on['diffuse_down']+on['direct_down']
            out['up'][u,:] = o['diffuse_up']
            out['up_noa'][u,:] = on['diffuse_up']
            
        except KeyboardInterrupt:
            raise KeyboardInterruptError()

        except:
            out['dn'][u,:] = np.zeros(len(zout))+np.nan
            out['dn_noa'][u,:] = np.zeros(len(zout))+np.nan
            out['up'][u,:] = np.zeros(len(zout))+np.nan
            out['up_noa'][u,:] = np.zeros(len(zout))+np.nan

    out['dn_avg'] = np.nanmean(out['dn'],axis=0)
    out['dn_noa_avg'] = np.nanmean(out['dn_noa'],axis=0)
    out['up_avg'] = np.nanmean(out['up'],axis=0)
    out['up_noa_avg'] = np.nanmean(out['up_noa'],axis=0)
    return out


# In[153]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[154]:


p = Pool(cpu_count()-1,worker_init)


# In[155]:


outputs = []
max_ = len(ar['Start_UTC'][fla])
with tqdm(total=max_) as pbar:
    for i, outs in tqdm(enumerate(p.imap_unordered(read_files, range(0, max_)))):
        pbar.update()
        outputs.append(outs)


# In[156]:


len(outputs)


# In[157]:


for oo in outputs:
    dat_24h['dn_avg'][oo['i'],:] = oo['dn_avg']
    dat_24h['dn_noa_avg'][oo['i'],:] = oo['dn_noa_avg']
    dat_24h['up_avg'][oo['i'],:] = oo['up_avg']
    dat_24h['up_noa_avg'][oo['i'],:] = oo['up_noa_avg']
    dat_24h['dn'][oo['i'],:,:] = oo['dn']
    dat_24h['dn_noa'][oo['i'],:,:] = oo['dn_noa']
    dat_24h['up'][oo['i'],:,:] = oo['up']
    dat_24h['up_noa'][oo['i'],:,:] = oo['up_noa']


# In[158]:


dat_24h['dare'] = (dat_24h['dn']-dat_24h['up']) - (dat_24h['dn_noa']-dat_24h['up_noa'])


# In[159]:


dat_24h['utc'] = ar['Start_UTC'][fla]
dat_24h['lat'] = ar['Latitude'][fla]
dat_24h['lon'] = ar['Longitude'][fla]
dat_24h['doy'] = ar['doy'][fla]


# In[160]:


dat_24h['dare_avg'] = (dat_24h['dn_avg']-dat_24h['up_avg']) - (dat_24h['dn_noa_avg']-dat_24h['up_noa_avg'])


# ### combine

# In[110]:


dat['dare'] = (dat['dn']-dat['up']) - (dat['dn_noa']-dat['up_noa'])


# In[111]:


dat['utc'] = ar['Start_UTC'][fla]
dat['lat'] = ar['Latitude'][fla]
dat['lon'] = ar['Longitude'][fla]
dat['doy'] = ar['doy'][fla]


# In[206]:


dat['dare_avg'] = (dat['dn_avg']-dat['up_avg']) - (dat['dn_noa_avg']-dat['up_noa_avg'])


# In[207]:


dat['dare_avg'].shape


# ## Save the file

# In[105]:


dat1 = iterate_dict_unicode(dat)
print 'saving file to: '+fp+'{name}_DARE_aero_prop_{vv}.mat'.format(name=name,vv=vv)
hs.savemat(fp+'{name}_DARE_aero_prop_{vv}.mat'.format(name=name,vv=vv),dat1)


# In[78]:


dat1 = iterate_dict_unicode(dat)
print 'saving file to: '+fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vv)
hs.savemat(fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vv),dat1)


# In[161]:


dat1_24h = iterate_dict_unicode(dat_24h)
print 'saving file to: '+fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vvh)
hs.savemat(fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vvh),dat1_24h)


# In[ ]:




