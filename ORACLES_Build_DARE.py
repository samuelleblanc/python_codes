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
#  Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, CA, 2019-12-02

# # Import of modules

# In[1]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
import scipy.io as sio
from path_utils import getpath
import matplotlib.pyplot as plt
import load_utils as lu
from write_utils import nearest_neighbor, iterate_dict_unicode


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


from tqdm.notebook import tqdm 
from datetime import datetime


# In[4]:


from multiprocessing import Pool, cpu_count
from copy import deepcopy
import signal
import warnings
warnings.simplefilter('ignore')


# In[5]:


import Run_libradtran as Rl


# In[6]:


name = 'ORACLES'


# In[7]:


vv = 'v2'
vr = 'R3'


# In[8]:


fp = getpath(name)
fp_rtm = getpath('rtm')
fp_uvspec = getpath('uvspec_bin')+'uvspec'
matfile = fp+'{}_all_cld_ict.mat'.format(vr)
fp_uvspec_dat = getpath('uvspec_dat')
fp_rtmdat = fp_rtm+'dat/'


# # Load the files

# ## Load the 4STAR AOD

# In[9]:


ar = hs.loadmat(fp+'/aod_ict/v8/{v}/all_aod_ict_{v}_2016.mat'.format(v=vr))


# In[10]:


ar.keys()


# In[11]:


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


cl = hs.loadmat(fp+'data_other/ssfr_2016_retrieved_COD_{}.mat'.format(vv))


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
    for ka,kd in f5.iteritems():
        ae5[ka] = kd.value
    ae = ae5


# In[28]:


ke = ae.keys()
ke.sort()
ke


# In[29]:


ae['AOD_meas'][0]


# In[30]:


ae_dict['SSA']


# In[31]:


ae['SSA'].shape


# In[32]:


ae_dict['time']


# In[33]:


ae['time']/3600.0


# In[34]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[35]:


ar['doy'] = np.array([datetime.strptime(days[int(d)],'%Y%m%d').timetuple().tm_yday for d in ar['days']])


# In[36]:


datetime.strptime(days[4],'%Y%m%d').timetuple().tm_yday


# In[37]:


ar['time_ae'] = ar['Start_UTC']+(24.0*(ar['doy']-244))


# In[38]:


ar['time_ae']


# # Prepare the base dict and defaults

# In[39]:


from datetime import datetime
datetime(2015,11,17).timetuple().tm_yday


# In[40]:


# for all 4STAR aerosol arrays
fla = (ar['flag_acaod']==1) & ar['fl'] & ar['fl_QA'] & (ar['days']>2.0) 


# In[41]:


# for the cod and ref arrays
fld = (ar['days']>2.0) & (ar['days']!=17.0) 
flb = (ar['flag_acaod'][fld]==1) & ar['fl'][fld] & ar['fl_QA'][fld]


# In[42]:


len(ar['AOD0355'][fla])


# In[43]:


len(cod[flb])


# In[44]:


sum(np.isfinite(cod[~flb])),sum(np.isfinite(cod[flb])),len(cod[flb])


# In[45]:


ka = ar.keys()
ka.sort()
ka


# In[46]:


doy = datetime.strptime(dds[int(ar['days'][fla][0])],'%Y%m%d').timetuple().tm_yday


# In[47]:


doy


# In[48]:


geo = {'lat':ar['Latitude'][0],'lon':ar['Longitude'][0],'doy':doy,'zout':[0,1.5,100.0]}
aero_no = {} # none
cloud = {'ztop':1.0,'zbot':0.5,'write_moments_file':False}
source = {'wvl_range':[201.0,4900.0],'source':'solar','integrate_values':True,'run_fuliou':True,
          'dat_path':fp_uvspec_dat}
albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':5.0}


# In[49]:


cloud['phase'] = 'wc'
geo['sza'] = 40.0
cloud['tau'] = 2.0
cloud['ref'] = 5.0
pmom = Rl.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',deltascale=False)
cloud['moms_dict'] = pmom


# In[50]:


pmom['wvl'][0] = 0.250


# In[51]:


wvl = np.append(np.append([250.0],ae['wavelength']),4900.0)
wvl


# In[52]:


aero = {'expand_hg':True,'disort_phase':False,'z_arr':[2.0,5.0],
        'wvl_arr':wvl}


# In[53]:


def fx_aero(aprop):
    'Function the aerosol property a 2d matrix for height and spectra, and extend the wavelength from 250 to 4900 nm'
    atmp = np.append([aprop[0]],np.append(aprop,aprop[-1]))
    return np.array([atmp,atmp])


# In[54]:


def fx_ext(a0,a1,a2,wvl=wvl):
    'Function to create the extinction coefficients from 4STAR AODs'
    aod = np.exp(np.polyval([a2,a1,a0],np.log(wvl)))
    aod[-1] = 0.0 # set the last wavelength to zero
    return np.array([aod/3.0,aod*0.0])


# In[55]:


aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][0],ar['AOD_polycoef_a1'][fla][0],ar['AOD_polycoef_a2'][fla][0])


# In[56]:


aero['asy'] = fx_aero(ae['g_total'][0])


# In[57]:


aero['ssa'] = fx_aero(ae['SSA'][0])


# ## Prepare the file list and saving

# In[58]:


def isjupyter():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


# ## Write the files

# ### Conventional

# In[59]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[59]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[120]:


# for writing out the files


# In[89]:


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

# In[64]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[66]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[67]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[168]:


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


# In[169]:


def write_files(d,cloud=cloud,source=source,albedo=albedo,aero_no=aero_no):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    Rl.write_input_aac(fpp_in+d['f_in'],geo=d['geo'],aero=d['aero'],cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
    Rl.write_input_aac(fpp_in+d['f_in_noa'],geo=d['geo'],aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)


# In[65]:


p = Pool(cpu_count()-1,worker_init)


# In[69]:


len(bb)


# In[174]:


results = p.map(write_files,bb)


# ## Run the calculations

# Run the files from command line:
# 
# using command parallel --jobs=20 < ORACLES_DARE_v1.sh

# In[175]:


f_list = fp_rtm+'{}_DARE_{}.sh'.format(name,vv)


# In[176]:


get_ipython().system(u' wc -l $f_list')


# In[179]:


f_listout = f_list+'.out'


# In[191]:


get_ipython().system(u'parallel --jobs=7 --bar < $f_list 2> $f_listout')


# ## Read the files

# In[61]:


fpp_out,name,vv,geo['zout']


# In[62]:


n = len(ar['Start_UTC'][fla])
nz = len(geo['zout'])
nw = len(aero['wvl_arr'])
dat = {'cod':np.zeros(n)+np.nan,'ref':np.zeros(n)+np.nan,'ext':np.zeros((n,nw))+np.nan,
       'ssa':np.zeros((n,nw))+np.nan,'asy':np.zeros((n,nw))+np.nan,'zout':geo['zout'],
       'wvl':aero['wvl_arr'],'sza':np.zeros(n)+np.nan,
       'dn':np.zeros((n,nz))+np.nan,'up':np.zeros((n,nz))+np.nan,
       'dn_noa':np.zeros((n,nz))+np.nan,'up_noa':np.zeros((n,nz))+np.nan}

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

# In[63]:


class KeyboardInterruptError(Exception): pass


# In[64]:


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


# In[65]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[77]:


p = Pool(cpu_count()-1,worker_init)


# In[78]:


outputs = []
max_ = len(ar['Start_UTC'][fla])
with tqdm(total=max_) as pbar:
    for i, outs in tqdm(enumerate(p.imap_unordered(read_files, range(0, max_)))):
        pbar.update()
        outputs.append(outs)


# In[79]:


outputs 


# In[81]:


dat.keys()


# In[83]:


for oo in outputs:
    dat['dn'][oo['i'],:] = oo['dn']
    dat['dn_noa'][oo['i'],:] = oo['dn_noa']
    dat['up'][oo['i'],:] = oo['up']
    dat['up_noa'][oo['i'],:] = oo['up_noa']


# ### combine

# In[87]:


dat['dare'] = (dat['dn']-dat['up']) - (dat['dn_noa']-dat['up_noa'])


# In[88]:


dat['utc'] = ar['Start_UTC'][fla]
dat['lat'] = ar['Latitude'][fla]
dat['lon'] = ar['Longitude'][fla]
dat['doy'] = ar['doy'][fla]


# ## Save the file

# In[89]:


dat1 = iterate_dict_unicode(dat)
print 'saving file to: '+fp+'{name}_DARE_aero_prop_{vv}.mat'.format(name=name,vv=vv)
hs.savemat(fp+'{name}_DARE_aero_prop_{vv}.mat'.format(name=name,vv=vv),dat1)


# In[90]:


dat1 = iterate_dict_unicode(dat)
print 'saving file to: '+fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vv)
hs.savemat(fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vv),dat1)


# In[ ]:




