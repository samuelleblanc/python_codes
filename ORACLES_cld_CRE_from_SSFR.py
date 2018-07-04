
# coding: utf-8

# # Intro
# Name:  
# 
#     ORACLES_cld_CRE_from_SSFR
# 
# Purpose:  
# 
#     Build the cloud radiative effect input files from the cloud retrieval exported from ORACLES_SSFR_cloud_retrieval file
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
#      Written: by Samuel LeBlanc, Santa Cruz, CA, 2018-07-03

# # Import of modules

# In[4]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu


# In[5]:


from load_utils import load_from_json


# In[6]:


name = 'ORACLES'


# In[7]:


vv = 'v5_irr'


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
    matfile = fp+'ssfr_2016_retrieved_COD.mat'
else:
    raise Exception


# ## Set up for command line arguments

# In[75]:


import argparse


# In[76]:


long_description = """    Prepare the Cloud radiative effect files for calculations and thn save them using the doread argument"""


# In[ ]:


parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-doread','--doread',help='if set, will only read the output, not produce them',
                    action='store_true')


# In[ ]:


in_ = vars(parser.parse_args())
do_read = in_.get('doread',False)


# # Load the saved files

# In[5]:


ar = hs.loadmat(matfile)


# In[6]:


ar.keys()


# In[12]:


dds = ['20160830','20160831','20160902','20160904','20160906','20160908',
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


geo = {'lat':47.6212167,'lon':52.74245,'doy':321,'zout':[0,1.5,100.0]}
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
pmom = Rl.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar')
cloud['moms_dict'] = pmom


# In[17]:


phase_star = {0:'wc',1:'ic'}


# In[18]:


phase_modis = {0:'wc',1:'wc',2:'ic',3:'ic',6:'wc'}


# ## Load the aerosol values

# In[82]:


if os.sys.platform == 'win32':
        fp_aero = fp+'model\\aero_save_v2.txt'
else:
        fp_aero = fp+'aero_save_v2.txt'
aero = load_from_json(fp_aero)


# In[13]:


aero = load_from_json(u'/mnt/c/Users/sleblanc/Research/ORACLES/model/aero_save_v2.txt')


# In[14]:


aero


# In[15]:


wv = np.array(aero['wvl_arr'])


# In[17]:


aero['ext'].shape


# In[18]:


aero['ext'][0,:]


# ## Prepare the paths and files for input files

# In[71]:


# open the list file
f = open(fp+'rtm/{}_CRE_{}.sh'.format(name,vv),'w')
fpp_in = '/nobackup/sleblan2/rtm/input/{}_CRE_ssfr_{}/'.format(name,vv)
fpp_out = '/nobackup/sleblan2/rtm/output/{}_CRE_ssfr_{}/'.format(name,vv)
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


if not do_read:


# In[ ]:


# make input
    for i,l in enumerate(ar['lat']):
        if np.isnan(ar['tau'][i]):
            continue
        print i

        f_in = '{name}_{vv}_ssfr_{i:03d}_withaero.dat'.format(name=name,vv=vv,i=i)
        geo['lat'],geo['lon'],geo['sza'] = l,ar['lon'][i],ar['sza'][i]
        day = dds[ar['days'][i]]
        geo['doy'] = datetime(int(day[0:4]),int(day[4:6]),int(day[6:])).timetuple().tm_yday
        cloud['tau'],cloud['ref'] = ar['tau'][i],ar['ref'][i]
        cloud['write_moments_file'] = True
        ext = np.exp(np.polyval([ar['a2'][i],ar['a1'][i],ar['a0'][i]],np.log(wv)))/3.0
        aero['ext'][0,:] = ext
        
        Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))

        f_in = '{name}_{vv}_ssfr_{i:03d}_withaero_clear.dat'.format(name=name,vv=vv,i=i)
        cloud['tau'] = 0.0
        Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))

       # f_in = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vv,i=i)
       # cloud['tau'] = ar['tau_fl'][i]
       # if cloud['ref']>25.0:
       #     cloud['write_moments_file'] = True
       # else:
       #     cloud['write_moments_file'] = False
       # Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
       #                            verbose=False,make_base=False,set_quiet=True)
       # f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))#

      #  f_in = '{name}_{vv}_star_{i:03d}_noaero_clear.dat'.format(name=name,vv=vv,i=i)
      #  cloud['tau'] = 0.0
      #  Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
      #                             verbose=False,make_base=False,set_quiet=True)
      #  f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))

    f.close()


# In[ ]:


else:


# In[ ]:


# read output
    nstar = len(ar['lat'])
    nz = len(geo['zout'])
    star_aero_CRE = {'dn':np.zeros((nstar,nz))+np.nan,'up':np.zeros((nstar,nz))+np.nan}
    star_aero_CRE_clear = {'dn':np.zeros((nstar,nz))+np.nan,'up':np.zeros((nstar,nz))+np.nan}
    star_aero_C = np.zeros((nstar,nz))+np.nan
 #   star_noaero_CRE = {'dn':np.zeros((nstar,nz))+np.nan,'up':np.zeros((nstar,nz))+np.nan}
 #   star_noaero_CRE_clear = {'dn':np.zeros((nstar,nz))+np.nan,'up':np.zeros((nstar,nz))+np.nan}
 #   star_noaero_C = np.zeros((nstar,nz))+np.nan


# In[ ]:


# run through to read
    print 'SSFR cloud retrievals CRE'
    for i,l in enumerate(ar['lat']):
        if np.isnan(ar['tau'][i]):
            continue
        print '\r{}..'.format(i)
        f_in = '{name}_{vv}_ssfr_{i:03d}_withaero.dat'.format(name=name,vv=vv,i=i)
        s = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])
        f_in = '{name}_{vv}_ssfr_{i:03d}_withaero_clear.dat'.format(name=name,vv=vv,i=i)
        sc = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])

        star_aero_CRE['dn'][i,:] = s['diffuse_down']+s['direct_down']
        star_aero_CRE_clear['dn'][i,:] = sc['diffuse_down']+sc['direct_down']
        star_aero_CRE['up'][i,:] = s['diffuse_up']
        star_aero_CRE_clear['up'][i,:] = sc['diffuse_up']
        star_aero_C[i,:] = (star_aero_CRE['dn'][i,:]-star_aero_CRE['up'][i,:]) -                            (star_aero_CRE_clear['dn'][i,:]-star_aero_CRE_clear['up'][i,:])
        
      #  f_in = '{name}_{vv}_star_{i:03d}_noaero.dat'.format(name=name,vv=vv,i=i)
      #  sn = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])
      #  f_in = '{name}_{vv}_star_{i:03d}_noaero_clear.dat'.format(name=name,vv=vv,i=i)
      #  snc = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])

      #  star_noaero_CRE['dn'][i,:] = sn['diffuse_down']+sn['direct_down']
      #  star_noaero_CRE_clear['dn'][i,:] = snc['diffuse_down']+snc['direct_down']
      #  star_noaero_CRE['up'][i,:] = sn['diffuse_up']
      #  star_noaero_CRE_clear['up'][i,:] = snc['diffuse_up']
      #  star_noaero_C[i,:] = (star_noaero_CRE['dn'][i,:]-star_noaero_CRE['up'][i,:]) - \
      #                       (star_noaero_CRE_clear['dn'][i,:]-star_noaero_CRE_clear['up'][i,:])


# In[ ]:


# save the output
    star1 = {'ssfr_aero_CRE':star_aero_CRE,'ssfr_aero_CRE_clear':star_aero_CRE_clear,'ssfr_aero_C':star_aero_C}
    star = wu.iterate_dict_unicode(star1)
    print 'saving file to: '+fp+'{name}_CRE_{vv}.mat'.format(name=name,vv=vv)
    hs.savemat(fp+'{name}_CRE_{vv}.mat'.format(name=name,vv=vv),star)
    #hs.savemat(fp+'{name}_CRE_{vv}.mat'.format(name=name,vv=vv),star_noaero_CRE,star_noaero_CRE_clear,star_noaero_C,
     #                                                           star_aero_CRE,star_aero_CRE_clear,star_aero_C)

