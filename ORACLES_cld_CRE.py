
# coding: utf-8

# # Intro
# Name:  
# 
#     ORACLES_cld_CRE
# 
# Purpose:  
# 
#     Build the cloud radiative effect input files from the cloud retrieval exported from ORACLES_cld_explore file
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
#      Written: by Samuel LeBlanc, Bathurst, NB, 2017-01-06

# # Import of modules

# In[46]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib notebook')
import numpy as np
import hdf5storage as hs
import os


# In[44]:

from load_utils import load_from_json


# In[19]:

name = 'ORACLES'


# In[48]:

vv = 'v1'
vr = 'R0'


# In[49]:

if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\ORACLES\\'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\ORACLES\\rtm\\'
    fp_uvspec = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\bin\\uvspec'
    fp_rtmdat = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\data\\'
    matfile = fp+'..//zen_ict/v3/{}_all_cld_ict.mat'.format(vr)
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/ORACLES/'
    fp_rtm = '/nobackup/sleblan2/rtm/'
    fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
    fp_rtmdat = '/nobackup/sleblan2/AAC_DARF/rtm/' #'/u/sleblan2/4STAR/rtm_dat/'
    matfile = fp+'..//zen_ict/v3/{}_all_cld_ict.mat'.format(vr)
else:
    raise Exception


# # Load the saved files

# In[5]:

ar = hs.loadmat(matfile)


# In[6]:

ar.keys()


# In[12]:

dds = ['20160827','20160830','20160831','20160902','20160904','20160906','20160908',
       '20160910','20160912','20160914','20160918','20160920','20160924','20160925','20160927']


# # Prepare input files for radiative transfer

# In[7]:

import Run_libradtran as Rl


# ## Prepare the defaults

# In[10]:

from datetime import datetime
datetime(2015,11,17).timetuple().tm_yday


# In[14]:

ar['days']


# In[15]:

geo = {'lat':47.6212167,'lon':52.74245,'doy':321,'zout':[0,100.0]}
aero_no = {} # none
cloud = {'ztop':1.0,'zbot':0.5,'write_moments_file':False}
source = {'wvl_range':[201.0,4000.0],'source':'solar','integrate_values':True,'run_fuliou':True,
          'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/'}
albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':5.0}


# In[16]:

cloud['phase'] = 'wc'
geo['sza'] = 40.0
cloud['tau'] = 2.0
cloud['ref'] = 5.0


# In[17]:

phase_star = {0:'wc',1:'ic'}


# In[18]:

phase_modis = {0:'wc',1:'wc',2:'ic',3:'ic',6:'wc'}


# ## Load the aerosol values

# In[50]:

aero = load_from_json(fp+'aero_save.txt')


# ## Prepare the paths and files for input files

# In[71]:

# open the list file
f = open(fp+'rtm/{}_CRE.sh'.format(name),'w')
fpp_in = '/nobackup/sleblan2/rtm/input/{}_CRE/'.format(name)
fpp_out = '/nobackup/sleblan2/rtm/output/{}_CRE/'.format(name)
fp_uv = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
fp_in = fp+'rtm/input/CRE/'


# In[74]:

if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[23]:

ar.keys()


# In[ ]:

for i,l enumerate(ar['lat_fl']):
    
    print i
    
    f_in = '{name}_{vv}_star_{i:03d}_withaero.dat'.format(name=name,vv=vv,i=i)
    geo['lat'],geo['lon'],geo['sza'] = l,ar['lon_fl'][i],ar['sza'][ar['fl'].astype(bool)][i]
    day = dds[ar['days'][ar['fl'].astype(bool)][i].astype(int)]
    geo['doy'] = datetime(int(day[0:4]),int(day[4:6]),int(day[6:])).timetuple().tm_yday
    cloud['tau'],cloud['ref'] = ar['tau_fl'][i],ar['ref_fl'][i]
    Rl.write_input_aac(fp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
    f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
    
    f_in = '{name}_{vv}_star_{i:03d}_withaero_clear.dat'.format(name=name,vv=vv,i=i)
    cloud['tau'] = 0.0
    Rl.write_input_aac(fp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
    f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
    
    f_in = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vv,i=i)
    cloud['tau'] = ar['tau_fl'][i]
    Rl.write_input_aac(fp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
    f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
    
    f_in = '{name}_{vv}_star_{i:03d}_noaero_clear.dat'.format(name=name,vv=vv,i=i)
    cloud['tau'] = 0.0
    Rl.write_input_aac(fp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
    f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
    
f.close()

