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


import Run_libradtran as Rl


# In[5]:


name = 'ORACLES'


# In[6]:


vv = 'v2'
vr = 'R3'


# In[7]:


fp = getpath(name)
fp_rtm = getpath('rtm')
fp_uvspec = getpath('uvspec_bin')+'uvspec'
matfile = fp+'{}_all_cld_ict.mat'.format(vr)
fp_uvspec_dat = getpath('uvspec_dat')
fp_rtmdat = fp_rtm+'dat/'


# # Load the files

# ## Load the 4STAR AOD

# In[8]:


ar = hs.loadmat(fp+'/aod_ict/v8/{v}/all_aod_ict_{v}_2016.mat'.format(v=vr))


# In[9]:


ar.keys()


# In[10]:


ar['AOD0501'].shape


# In[11]:


sza = np.arccos(1.0/ar['amass_aer'])*180.0/np.pi


# In[12]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[123]:


len(days)


# In[122]:


ar['days']


# ## Load the retrieved Cloud properties

# In[13]:


cl = hs.loadmat(fp+'data_other/ssfr_2016_retrieved_COD_{}.mat'.format(vv))


# In[14]:


cl.keys()


# In[15]:


cl['tau'].shape


# In[16]:


dds = ['20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927']


# In[124]:


len(dds)


# In[125]:


cl['days']


# In[18]:


dd = np.unique(cl['days'])


# In[128]:


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


# In[129]:


cod.shape


# In[130]:


len(dds)


# In[131]:


len(np.unique(ar['days']))


# ## Load the skyscan retrievals

# In[23]:


try:
    ae,ae_dict = lu.load_netcdf(fp+'aeroinv_2016/netcdf4/4STAR-aeroinv_P3_2016_R0.nc',everything=True)
except:
    import h5py as h5
    f5 = h5.File(fp+'aeroinv_2016/netcdf4/4STAR-aeroinv_P3_2016_R0.nc','r')
    ae5 = {}
    for ka,kd in f5.iteritems():
        ae5[ka] = kd.value
    ae = ae5


# In[24]:


ke = ae.keys()
ke.sort()
ke


# In[192]:


ae['AOD_meas'][0]


# In[26]:


ae_dict['SSA']


# In[27]:


ae['SSA'].shape


# In[28]:


ae_dict['time']


# In[29]:


ae['time']/3600.0


# In[30]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[31]:


ar['doy'] = np.array([datetime.strptime(days[int(d)],'%Y%m%d').timetuple().tm_yday for d in ar['days']])


# In[32]:


datetime.strptime(days[4],'%Y%m%d').timetuple().tm_yday


# In[33]:


ar['time_ae'] = ar['Start_UTC']+(24.0*(ar['doy']-244))


# In[34]:


ar['time_ae']


# # Prepare the base dict and defaults

# In[35]:


from datetime import datetime
datetime(2015,11,17).timetuple().tm_yday


# In[120]:


# for all 4STAR aerosol arrays
fla = (ar['flag_acaod']==1) & ar['fl'] & ar['fl_QA'] & (ar['days']>2.0) 


# In[121]:


# for the cod and ref arrays
fld = (ar['days']>2.0) & (ar['days']!=17.0) 
flb = (ar['flag_acaod'][fld]==1) & ar['fl'][fld] & ar['fl_QA'][fld]


# In[137]:


len(ar['AOD0355'][fla])


# In[138]:


len(cod[flb])


# In[145]:


sum(np.isfinite(cod[~flb])),sum(np.isfinite(cod[flb])),len(cod[flb])


# In[38]:


ka = ar.keys()
ka.sort()
ka


# In[146]:


doy = datetime.strptime(dds[int(ar['days'][fla][0])],'%Y%m%d').timetuple().tm_yday


# In[147]:


doy


# In[148]:


geo = {'lat':ar['Latitude'][0],'lon':ar['Longitude'][0],'doy':doy,'zout':[0,1.5,100.0]}
aero_no = {} # none
cloud = {'ztop':1.0,'zbot':0.5,'write_moments_file':False}
source = {'wvl_range':[201.0,4900.0],'source':'solar','integrate_values':True,'run_fuliou':True,
          'dat_path':fp_uvspec_dat}
albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':5.0}


# In[149]:


cloud['phase'] = 'wc'
geo['sza'] = 40.0
cloud['tau'] = 2.0
cloud['ref'] = 5.0
pmom = Rl.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',deltascale=False)
cloud['moms_dict'] = pmom


# In[150]:


pmom['wvl'][0] = 0.250


# In[151]:


wvl = np.append(np.append([250.0],ae['wavelength']),4900.0)
wvl


# In[152]:


aero = {'expand_hg':True,'disort_phase':False,'z_arr':[2.0,5.0],
        'wvl_arr':wvl}


# In[153]:


def fx_aero(aprop):
    'Function the aerosol property a 2d matrix for height and spectra, and extend the wavelength from 250 to 4900 nm'
    atmp = np.append([aprop[0]],np.append(aprop,aprop[-1]))
    return np.array([atmp,atmp])


# In[154]:


def fx_ext(a0,a1,a2,wvl=wvl):
    'Function to create the extinction coefficients from 4STAR AODs'
    aod = np.exp(np.polyval([a2,a1,a0],np.log(wvl)))
    aod[-1] = 0.0 # set the last wavelength to zero
    return np.array([aod/3.0,aod*0.0])


# In[155]:


aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][0],ar['AOD_polycoef_a1'][fla][0],ar['AOD_polycoef_a2'][fla][0])


# In[156]:


aero['asy'] = fx_aero(ae['g_total'][0])


# In[157]:


aero['ssa'] = fx_aero(ae['SSA'][0])


# ## Prepare the file list and saving

# In[158]:


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

# In[87]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[88]:


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

# In[159]:


from multiprocessing import Pool, cpu_count


# In[160]:


from copy import deepcopy
import signal


# In[161]:


import warnings
warnings.simplefilter('ignore')


# In[162]:


def worker_init():
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[163]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[164]:


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


# In[170]:


p = Pool(cpu_count()-1,worker_init)


# In[171]:


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

# In[195]:


fpp_out,name,vv,geo['zout']


# In[198]:


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


# In[203]:


def read_files(i,dat=dat,fpp_out=fpp_out,name=name,vv=vv,zout=geo['zout']):
    'function to feed the pool of workers to read all the files'
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


# In[ ]:


outputs = p.map(read_files,range(len(ar['Start_UTC'][fla])))


# For the original method (without multiprocessing)

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


# In[134]:


dat['dare'] = (dat['dn']-dat['up']) - (dat['dn_noa']-dat['up_noa'])


# In[52]:


dat['utc'] = ar['Start_UTC'][fla]
dat['lat'] = ar['Latitude'][fla]
dat['lon'] = ar['Longitude'][fla]
dat['doy'] = ar['doy'][fla]


# ## Save the file

# In[53]:


dat1 = iterate_dict_unicode(dat)
print 'saving file to: '+fp+'{name}_DARE_aero_prop_{vv}.mat'.format(name=name,vv=vv)
hs.savemat(fp+'{name}_DARE_aero_prop_{vv}.mat'.format(name=name,vv=vv),dat1)


# In[135]:


dat1 = iterate_dict_unicode(dat)
print 'saving file to: '+fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vv)
hs.savemat(fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vv),dat1)

