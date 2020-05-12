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


# In[16]:


vv = 'v3'
vr = 'R3'


# In[8]:


fp = getpath(name)
fp_rtm = getpath('rtm')
fp_uvspec = getpath('uvspecb')+'uvspec'
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


# In[17]:


sza = np.arccos(1.0/ar['amass_aer'])*180.0/np.pi


# In[18]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[19]:


len(days)


# In[20]:


ar['days']


# ## Load the retrieved Cloud properties

# In[21]:


cl = hs.loadmat(fp+'data_other/ssfr_2016_retrieved_COD_{}.mat'.format(vv))


# In[22]:


cl.keys()


# In[23]:


cl['tau'].shape


# In[24]:


dds = ['20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927']


# In[25]:


len(dds)


# In[26]:


cl['days']


# In[27]:


dd = np.unique(cl['days'])


# In[28]:


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


# In[29]:


cod.shape


# In[30]:


len(dds)


# In[31]:


len(np.unique(ar['days']))


# ## Load the skyscan retrievals

# In[32]:


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


# In[33]:


ke = ae.keys()
ke.sort()
ke


# In[34]:


ae['AOD_meas'][0]


# In[35]:


ae_dict['AAOD']


# In[50]:


ae_dict['SSA']


# In[51]:


ae['SSA'].shape


# In[52]:


ae_dict['time']


# In[36]:


ae['time']/3600.0


# In[37]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[38]:


ar['doy'] = np.array([datetime.strptime(days[int(d)],'%Y%m%d').timetuple().tm_yday for d in ar['days']])


# In[39]:


datetime.strptime(days[4],'%Y%m%d').timetuple().tm_yday


# In[40]:


ar['time_ae'] = ar['Start_UTC']+(24.0*(ar['doy']-244))


# In[41]:


ar['time_ae']


# # Prepare the base dict and defaults

# In[42]:


from datetime import datetime
datetime(2015,11,17).timetuple().tm_yday


# In[43]:


# for all 4STAR aerosol arrays
fla = (ar['flag_acaod']==1) & ar['fl'] & ar['fl_QA'] & (ar['days']>2.0) 


# In[44]:


# for the cod and ref arrays
fld = (ar['days']>2.0) & (ar['days']!=17.0) 
flb = (ar['flag_acaod'][fld]==1) & ar['fl'][fld] & ar['fl_QA'][fld]


# In[45]:


len(ar['AOD0355'][fla])


# In[46]:


len(cod[flb])


# In[47]:


sum(np.isfinite(cod[~flb])),sum(np.isfinite(cod[flb])),len(cod[flb])


# In[236]:


ka = ar.keys()
ka.sort()
ka


# In[48]:


doy = datetime.strptime(dds[int(ar['days'][fla][0])],'%Y%m%d').timetuple().tm_yday


# In[49]:


doy


# In[50]:


geo = {'lat':ar['Latitude'][0],'lon':ar['Longitude'][0],'doy':doy,'zout':[0,1.5,100.0]}
aero_no = {} # none
cloud = {'ztop':1.0,'zbot':0.5,'write_moments_file':False}
source = {'wvl_range':[201.0,4900.0],'source':'solar','integrate_values':True,'run_fuliou':True,
          'dat_path':fp_uvspec_dat}
albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':5.0}


# In[51]:


cloud['phase'] = 'wc'
geo['sza'] = 40.0
cloud['tau'] = 2.0
cloud['ref'] = 5.0
pmom = Rl.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',deltascale=False)
cloud['moms_dict'] = pmom


# In[52]:


pmom['wvl'][0] = 0.250


# In[53]:


wvl = np.append(np.append([250.0],ae['wavelength']),4900.0)
wvl


# In[54]:


aero = {'expand_hg':True,'disort_phase':False,'z_arr':[2.0,5.0],
        'wvl_arr':wvl}


# In[55]:


def fx_aero(aprop):
    'Function the aerosol property a 2d matrix for height and spectra, and extend the wavelength from 250 to 4900 nm'
    atmp = np.append([aprop[0]],np.append(aprop,aprop[-1]))
    return np.array([atmp,atmp])


# In[56]:


def fx_ext(a0,a1,a2,wvl=wvl):
    'Function to create the extinction coefficients from 4STAR AODs'
    aod = np.exp(np.polyval([a2,a1,a0],np.log(wvl)))
    aod[-1] = 0.0 # set the last wavelength to zero
    return np.array([aod/3.0,aod*0.0])


# In[57]:


aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][0],ar['AOD_polycoef_a1'][fla][0],ar['AOD_polycoef_a2'][fla][0])


# In[58]:


aero['asy'] = fx_aero(ae['g_total'][0])


# In[59]:


aero['ssa'] = fx_aero(ae['SSA'][0])


# ## Prepare the file list and saving

# In[60]:


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

# In[250]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[251]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[252]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[253]:


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


# In[254]:


def write_files(d,cloud=cloud,source=source,albedo=albedo,aero_no=aero_no):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    Rl.write_input_aac(fpp_in+d['f_in'],geo=d['geo'],aero=d['aero'],cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
    Rl.write_input_aac(fpp_in+d['f_in_noa'],geo=d['geo'],aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)


# In[255]:


p = Pool(cpu_count()-1,worker_init)


# In[256]:


len(bb)


# In[257]:


results = p.map(write_files,bb)


# ### Diurnal averaging with multiprocessing

# In[61]:


vv = 'v3_24h'


# In[62]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[63]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[64]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[65]:


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


# In[71]:


# Go through and write out the list file for 24 h
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
if isjupyter():
    pbar = tqdm(total=len(ar['Start_UTC'][fla]))
for i,u in enumerate(ar['Start_UTC'][fla]):
    iae = np.argmin(abs(ar['time_ae'][fla][i]-ae['time']/3600.0))
    if abs(ar['time_ae'][fla][i]-ae['time']/3600.0)[iae]<1.0: 
        for ux in np.arange(0,24,0.5):
            f_in = '{name}_{vv}_DARE_{i:03d}_withaero_{ux:04.1f}.dat'.format(name=name,vv=vv,i=i,ux=ux)
            f_in_noa = '{name}_{vv}_star_{i:03d}_noaero_{ux:04.1f}.dat'.format(name=name,vv=vv,i=i,ux=ux)
            f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in,out=fpp_out+f_in))
            f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uvspec,fin=fpp_in+f_in_noa,out=fpp_out+f_in_noa))
        pbar.update(1)  
f.close()


# In[66]:


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


# In[67]:


p = Pool(cpu_count()-1,worker_init)


# In[68]:


results = []
max_ = len(bb)
with tqdm(total=max_) as pbar:
    for i, res in tqdm(enumerate(p.imap_unordered(write_files, bb))):
        pbar.update()
        results.append(res)


# In[ ]:


##results = p.map(write_files,bb)


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

# In[258]:


f_list = fp_rtm+'{}_DARE_{}.sh'.format(name,vv)


# In[259]:


get_ipython().system(u' wc -l $f_list')


# In[260]:


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

# In[263]:


fpp_out,name,vv,geo['zout']


# In[264]:


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

# In[265]:


class KeyboardInterruptError(Exception): pass


# In[266]:


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


# In[267]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[268]:


p = Pool(cpu_count()-1,worker_init)


# In[269]:


outputs = []
max_ = len(ar['Start_UTC'][fla])
with tqdm(total=max_) as pbar:
    for i, outs in tqdm(enumerate(p.imap_unordered(read_files, range(0, max_)))):
        pbar.update()
        outputs.append(outs)


# In[270]:


outputs[0],outputs[2000]


# In[271]:


dat.keys()


# In[272]:


for oo in outputs:
    dat['dn'][oo['i'],:] = oo['dn']
    dat['dn_noa'][oo['i'],:] = oo['dn_noa']
    dat['up'][oo['i'],:] = oo['up']
    dat['up_noa'][oo['i'],:] = oo['up_noa']


# ### combine

# In[273]:


dat['dare'] = (dat['dn']-dat['up']) - (dat['dn_noa']-dat['up_noa'])


# In[274]:


dat['utc'] = ar['Start_UTC'][fla]
dat['lat'] = ar['Latitude'][fla]
dat['lon'] = ar['Longitude'][fla]
dat['doy'] = ar['doy'][fla]


# ## Save the file

# In[275]:


dat1 = iterate_dict_unicode(dat)
print 'saving file to: '+fp+'{name}_DARE_aero_prop_{vv}.mat'.format(name=name,vv=vv)
hs.savemat(fp+'{name}_DARE_aero_prop_{vv}.mat'.format(name=name,vv=vv),dat1)


# In[276]:


dat1 = iterate_dict_unicode(dat)
print 'saving file to: '+fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vv)
hs.savemat(fp+'{name}_DARE_{vv}.mat'.format(name=name,vv=vv),dat1)


# In[ ]:




