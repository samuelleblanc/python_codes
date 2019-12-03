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

# In[126]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
import scipy.io as sio
from path_utils import getpath
import matplotlib.pyplot as plt
import load_utils as lu
from write_utils import nearest_neighbor


# In[127]:


get_ipython().magic(u'matplotlib notebook')


# In[4]:


from tqdm.notebook import tqdm 


# In[22]:


import Run_libradtran as Rl


# In[5]:


name = 'ORACLES'


# In[13]:


vv = 'v1'
vr = 'R3'


# In[40]:


fp = getpath(name)
fp_rtm = getpath('rtm')
fp_uvspec = getpath('uvspec_bin')+'uvspec'
matfile = fp+'{}_all_cld_ict.mat'.format(vr)
fp_uvspec_dat = getpath('uvspec_dat')
fp_rtmdat = fp_rtm+'dat/'


# # Load the files

# ## Load the 4STAR AOD

# In[14]:


ar = hs.loadmat(fp+'/aod_ict/v8/{v}/all_aod_ict_{v}_2016.mat'.format(v=vr))


# In[15]:


ar.keys()


# In[23]:


ar['AOD0501'].shape


# In[189]:


sza = np.arccos(1.0/ar['amass_aer'])*180.0/np.pi


# In[149]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# ## Load the retrieved Cloud properties

# In[17]:


cl = hs.loadmat(fp+'data_other/ssfr_2016_retrieved_COD.mat')


# In[18]:


cl.keys()


# In[24]:


cl['tau'].shape


# In[69]:


dds = ['20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927']


# In[121]:


cl['days']


# In[124]:


dd = np.unique(cl['days'])


# In[150]:


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


# In[151]:


cod.shape


# In[152]:


len(dds)


# In[153]:


len(np.unique(ar['days']))


# ## Load the skyscan retrievals

# In[29]:


ae,ae_dict = lu.load_netcdf(fp+'aeroinv_2016/netcdf4/4STAR-aeroinv_P3_2016_R0.nc',everything=True)


# In[87]:


ke = ae.keys()
ke.sort()
ke


# In[88]:


ae['AOD_meas'][0]


# In[33]:


ae_dict['SSA']


# In[34]:


ae['SSA'].shape


# In[194]:


ae_dict['time']


# In[197]:


ae['time']/3600.0


# In[ ]:


days = ['20160824','20160825','20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927','20160930']


# In[199]:


ar['doy'] = np.array([datetime.strptime(days[int(d)],'%Y%m%d').timetuple().tm_yday for d in ar['days']])


# In[201]:


datetime.strptime(days[4],'%Y%m%d').timetuple().tm_yday


# In[207]:


ar['time_ae'] = ar['Start_UTC']+(24.0*(ar['doy']-244))


# In[204]:


ar['time_ae']


# # Prepare the base dict and defaults

# In[35]:


from datetime import datetime
datetime(2015,11,17).timetuple().tm_yday


# In[164]:


# for all 4STAR aerosol arrays
fla = (ar['flag_acaod']==1) & ar['fl'] & ar['fl_QA'] & (ar['days']>2.0)


# In[165]:


# for the cod and ref arrays
fld = (ar['days']>2.0) & (ar['days']!=17.0) 
flb = (ar['flag_acaod'][fld]==1) & ar['fl'][fld] & ar['fl_QA'][fld]


# In[66]:


ka = ar.keys()
ka.sort()
ka


# In[80]:


doy = datetime.strptime(dds[int(ar['days'][fla][0])],'%Y%m%d').timetuple().tm_yday


# In[82]:


doy


# In[83]:


geo = {'lat':ar['Latitude'][0],'lon':ar['Longitude'][0],'doy':doy,'zout':[0,1.5,100.0]}
aero_no = {} # none
cloud = {'ztop':1.0,'zbot':0.5,'write_moments_file':False}
source = {'wvl_range':[201.0,5600.0],'source':'solar','integrate_values':True,'run_fuliou':True,
          'dat_path':fp_uvspec_dat}
albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':5.0}


# In[209]:


cloud['phase'] = 'wc'
geo['sza'] = 40.0
cloud['tau'] = 2.0
cloud['ref'] = 5.0
pmom = Rl.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',deltascale=False)
cloud['moms_dict'] = pmom


# In[43]:


pmom['wvl'][0] = 0.250


# In[118]:


wvl = np.append(np.append([250.0],ae['wavelength']),4900.0)
wvl


# In[119]:


aero = {'expand_hg':True,'disort_phase':False,'z_arr':[2.0,5.0],
        'wvl_arr':wvl}


# In[112]:


def fx_aero(aprop):
    'Function the aerosol property a 2d matrix for height and spectra, and extend the wavelength from 250 to 4900 nm'
    atmp = np.append([aprop[0]],np.append(aprop,aprop[-1]))
    return np.array([atmp,atmp])


# In[181]:


def fx_ext(a0,a1,a2,wvl=wvl):
    'Function to create the extinction coefficients from 4STAR AODs'
    aod = np.exp(np.polyval([a2,a1,a0],np.log(wvl)))
    aod[-1] = 0.0 # set the last wavelength to zero
    return np.array([aod/3.0,aod*0.0])


# In[182]:


aero['ext'] = fx_ext(ar['AOD_polycoef_a0'][fla][0],ar['AOD_polycoef_a1'][fla][0],ar['AOD_polycoef_a2'][fla][0])


# In[114]:


aero['asy'] = fx_aero(ae['g_total'][0])


# In[117]:


aero['ssa'] = fx_aero(ae['SSA'][0])


# ## Prepare the file list and saving

# In[44]:


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


# In[45]:


# open the list file
f = open(fp_rtm+'{}_DARE_{}.sh'.format(name,vv),'w')
fpp_in = fp_rtm+'input/{}_DARE_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_DARE_{}/'.format(name,vv)


# In[46]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[120]:


# for writing out the files


# In[211]:


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
        aero['ssa'] = fx_aero(ae['SSA'][0])
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


# In[ ]:




