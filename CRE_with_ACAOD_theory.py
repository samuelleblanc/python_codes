
# coding: utf-8

# # Info
# Name:  
# 
#     CRE_with_ACAOD_theory
# 
# Purpose:  
# 
#     Run CRE calculations with above cloud aerosol
#     Vary in systematic ways to quantify the range of CRE
#   
# Input:
# 
#     command line arguments:
#         - doread : for reading the output files
#   
# Output:
# 
#     libradtran input files
#     save files of the libradtran output
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
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2019-06-19
#     Modified: 

# # Prepare python environment

# In[22]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
from load_utils import load_from_json
import argparse
import Run_libradtran as Rl
from datetime import datetime
from path_utils import getpath


# In[34]:


name = 'ORACLES_theory'


# In[33]:


vv = 'v1'


# In[25]:


fp = getpath('ORACLES')
if os.sys.platform == 'win32':
    fp = getpath('ORACLES')
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\ORACLES\\rtm\\'
    fp_uvspec = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\bin\\uvspec'
    fp_rtmdat = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\data\\'
    #matfile = fp+'..//zen_ict/v3/{}_all_cld_ict.mat'.format(vr)
elif os.sys.platform == 'linux2':
    fp_rtm = '/nobackup/sleblan2/rtm/'
    fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
    fp_rtmdat = '/nobackup/sleblan2/AAC_DARF/rtm/' #'/u/sleblan2/4STAR/rtm_dat/'
    #matfile = fp+'{}_all_cld_ict.mat'.format(vr)
else:
    raise Exception


# ## Setup the command line arguments

# In[8]:


long_description = """    Prepare the Cloud radiative effect files for calculations and thn save them using the doread argument"""


# In[9]:


parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-doread','--doread',help='if set, will only read the output, not produce them',
                    action='store_true')


# In[10]:


in_ = vars(parser.parse_args())
do_read = in_.get('doread',False)


# # Prepare the radiative transfer

# In[13]:


geo = {'lat':47.6212167,'lon':52.74245,'doy':321,'zout':[0,1.5,100.0],'sza':30.0}
aero_no = {} # none
cloud = {'ztop':1.0,'zbot':0.5,'write_moments_file':False}
source = {'wvl_range':[201.0,4000.0],'source':'solar','integrate_values':True,'run_fuliou':True,
          'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/'}
albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':5.0}


# In[18]:


cloud['phase'] = 'wc'
geo['sza'] = 40.0
cloud['tau'] = 2.0
cloud['ref'] = 5.0
pmom = Rl.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar')
cloud['moms_dict'] = pmom


# In[19]:


phase_star = {0:'wc',1:'ic'}


# In[20]:


phase_modis = {0:'wc',1:'wc',2:'ic',3:'ic',6:'wc'}


# ## Load the aerosol values

# In[28]:


try:
    aero = load_from_json(fp+'model/aero_save_v2.txt')
except IOError:
    aero = load_from_json(fp+'aero_save_v2.txt')


# In[29]:


aero


# In[37]:


aero['ext'].shape


# In[30]:


# set the range of ext, ssa, and asy to model at 500 nm
ext_arr = [0.05,0.1,0.15,0.2,0.3]
ssa_arr = [0.75,0.8,0.85,0.875,0.9]
asy_arr = [0.6,0.65,0.7,0.75]


# ## Set up the cloud properties

# In[31]:


cod_arr = [1.0,2.5,5.0,7.5,10.0,12.5,15.0,20.0]
ref_arr = [2.0,5.0,7.5,10.0,12.5,15.0]


# # Prep the libradtran input file writing

# In[35]:


# open the list file
f = open(fp+'rtm/{}_CRE_{}.sh'.format(name,vv),'w')
fpp_in = '/nobackup/sleblan2/rtm/input/{}_CRE_{}/'.format(name,vv)
fpp_out = '/nobackup/sleblan2/rtm/output/{}_CRE_{}/'.format(name,vv)
fp_uv = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
fp_in = fp+'rtm/input/CRE/'


# In[ ]:


if not os.path.isdir(fpp_in):
    os.mkdir(fpp_in)
if not os.path.isdir(fpp_out):
     os.mkdir(fpp_out)


# In[36]:


if not do_read:


# In[ ]:


# make input
    
    for ic,c in enumerate(cod_arr): 
        for ir, r in enumerate(ref_arr):
            for je, e in enumerate(ext_arr):
                for js, s in enumerate(ssa_arr):
                    for ja, a in enumerate(asy_arr):
                        fm = {'ic':ic,'ir':ir,'je':je,'js':js,'ja':ja,'name':name,'vv':vv}
                        f_in = '{name}_{vv}_{ic:02d}{ir:02d}{je:02d}{js:02d}{ja:02d}_withaero.dat'.format(**fm)

                        aero['ext'] = aero['ext']*e/aero['ext'][0,3]
                        aero['ssa'] = aero['ssa']*s/aero['ssa'][0,3]
                        aero['asy'] = aero['asy']*a/aero['asy'][0,3]
                        cloud['tau'],cloud['ref'] = c,r
                        cloud['write_moments_file'] = True
                        
                        Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                                   verbose=False,make_base=False,set_quiet=True)
                        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))

                        f_in = '{name}_{vv}_{ic:02d}{ir:02d}{je:02d}{js:02d}{ja:02d}_withaero_clear.dat'.format(**fm)
                        cloud['tau'] = 0.0
                        Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                                   verbose=False,make_base=False,set_quiet=True)
                        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))

                        f_in = '{name}_{vv}_{ic:02d}{ir:02d}{je:02d}{js:02d}{ja:02d}_noaero.dat'.format(**fm)
                        cloud['tau'] = c
                        if cloud['ref']>25.0:
                            cloud['write_moments_file'] = True
                        else:
                            cloud['write_moments_file'] = False
                        Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                                                   verbose=False,make_base=False,set_quiet=True)
                        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))

                        f_in = '{name}_{vv}_{ic:02d}{ir:02d}{je:02d}{js:02d}{ja:02d}_noaero_clear.dat'.format(**fm)
                        cloud['tau'] = 0.0
                        Rl.write_input_aac(fpp_in+f_in,geo=geo,aero=aero_no,cloud=cloud,source=source,albedo=albedo,
                                                   verbose=False,make_base=False,set_quiet=True)
                        f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
                        
                        print 'cod:{ic:02d}, ref:{ir:02d}, ext:{je:02d}, ssa:{js:02d}, asy:{ja:02d}'.format(**fm)

    f.close()


# # Read the output

# In[ ]:


else:


# In[ ]:


# read output
    nc,nr,ne,ns,na = len(cod_arr),len(ref_arr),len(ext_arr),len(ssa_arr),len(asy_arr)
    nz = len(geo['zout'])
    star_aero_CRE = {'dn':np.zeros((nc,nr,ne,ns,na,nz))+np.nan,'up':np.zeros((nc,nr,ne,ns,na,nz))+np.nan}
    star_aero_CRE_clear = {'dn':np.zeros((nc,nr,ne,ns,na,nz))+np.nan,'up':np.zeros((nc,nr,ne,ns,na,nz))+np.nan}
    star_aero_C = np.zeros((nc,nr,ne,ns,na,nz))+np.nan
    star_noaero_CRE = {'dn':np.zeros((nc,nr,ne,ns,na,nz))+np.nan,'up':np.zeros((nc,nr,ne,ns,na,nz))+np.nan}
    star_noaero_CRE_clear = {'dn':np.zeros((nc,nr,ne,ns,na,nz))+np.nan,'up':np.zeros((nc,nr,ne,ns,na,nz))+np.nan}
    star_noaero_C = np.zeros((nc,nr,ne,ns,na,nz))+np.nan


# In[ ]:


# run through to read
    print '{name}'.format(**fm)
    for ic,c in enumerate(cod_arr): 
        for ir, r in enumerate(ref_arr):
            print 'cod:{ic:02d} - ref:{ir:02d}'.format(ic=ic,ir=ir)
            for je, e in enumerate(ext_arr):
                for js, s in enumerate(ssa_arr):
                    for ja, a in enumerate(asy_arr):
                        fm = {'ic':ic,'ir':ir,'je':je,'js':js,'ja':ja,'name':name,'vv':vv}
                        print '\r{je:02d}{js:02d}{ja:02d}..'.format(**fm)
                        f_in = '{name}_{vv}_{ic:02d}{ir:02d}{je:02d}{js:02d}{ja:02d}_withaero.dat'.format(**fm)
                        s = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])
                        f_in = '{name}_{vv}_{ic:02d}{ir:02d}{je:02d}{js:02d}{ja:02d}_withaero_clear.dat'.format(**fm)
                        sc = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])

                        star_aero_CRE['dn'][ic,ir,je,js,ja,:] = s['diffuse_down']+s['direct_down']
                        star_aero_CRE_clear['dn'][ic,ir,je,js,ja,:] = sc['diffuse_down']+sc['direct_down']
                        star_aero_CRE['up'][ic,ir,je,js,ja,:] = s['diffuse_up']
                        star_aero_CRE_clear['up'][ic,ir,je,js,ja,:] = sc['diffuse_up']
                        star_aero_C[ic,ir,je,js,ja,:] = (star_aero_CRE['dn'][ic,ir,je,js,ja,:]-star_aero_CRE['up'][ic,ir,je,js,ja,:]) -                                            (star_aero_CRE_clear['dn'][ic,ir,je,js,ja,:]-star_aero_CRE_clear['up'][ic,ir,je,js,ja,:])
        
                        f_in = '{name}_{vv}_{ic:02d}{ir:02d}{je:02d}{js:02d}{ja:02d}_noaero.dat'.format(**fm)
                        sn = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])
                        f_in = '{name}_{vv}_{ic:02d}{ir:02d}{je:02d}{js:02d}{ja:02d}_noaero_clear.dat'.format(**fm)
                        snc = Rl.read_libradtran(fpp_out+f_in,zout=geo['zout'])

                        star_noaero_CRE['dn'][ic,ir,je,js,ja,:] = sn['diffuse_down']+sn['direct_down']
                        star_noaero_CRE_clear['dn'][ic,ir,je,js,ja,:] = snc['diffuse_down']+snc['direct_down']
                        star_noaero_CRE['up'][ic,ir,je,js,ja,:] = sn['diffuse_up']
                        star_noaero_CRE_clear['up'][ic,ir,je,js,ja,:] = snc['diffuse_up']
                        star_noaero_C[ic,ir,je,js,ja,:] = (star_noaero_CRE['dn'][ic,ir,je,js,ja,:]-star_noaero_CRE['up'][ic,ir,je,js,ja,:]) -                                              (star_noaero_CRE_clear['dn'][ic,ir,je,js,ja,:]-star_noaero_CRE_clear['up'][ic,ir,je,js,ja,:])


# In[ ]:


# save the output
    star1 = {'star_noaero_CRE':star_noaero_CRE,'star_noaero_CRE_clear':star_noaero_CRE_clear,'star_noaero_C':star_noaero_C,
            'star_aero_CRE':star_aero_CRE,'star_aero_CRE_clear':star_aero_CRE_clear,'star_aero_C':star_aero_C}
    star = wu.iterate_dict_unicode(star1)
    print 'saving file to: '+fp+'{name}_CRE_{vv}.mat'.format(name=name,vv=vv)
    hs.savemat(fp+'{name}_CRE_{vv}.mat'.format(name=name,vv=vv),star)
    #hs.savemat(fp+'{name}_CRE_{vv}.mat'.format(name=name,vv=vv),star_noaero_CRE,star_noaero_CRE_clear,star_noaero_C,
     #                                                           star_aero_CRE,star_aero_CRE_clear,star_aero_C)

