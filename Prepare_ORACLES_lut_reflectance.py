#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     Prepare_ORACLES_lut_reflectance
# 
# Purpose:  
# 
#     Create the input libradtran files for creating a lut of low clouds with aerosol on top to be used in ORACLES operational 
#     cloud retrievals. Using the reflectance values calculated from SSFR measurements. Creating the lut for Irradiance values.
# 
# Calling Sequence:
# 
#     python Prepare_ORACLES_lut_reflectance
#   
# Input:
# 
#     none
# 
# Output:
#    
#     input files for libradtran 2.0 (uvspec) 
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - mplt_toolkits for basemap, map plotting
#     - pdb
#     - datetime
# 
#   
# Needed Files:
# 
#   - aero_file_v4.txt details of the aerosol layer above
#     
# History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2018-06-25
#              based on Prepare_ORACLES_lut
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2020-11-18
#              updates for new Mie scattering calculations: v7_irr
#     

# # Prepare the python environment

# In[1]:


import numpy as np
import scipy.io as sio
import os
import Run_libradtran as RL


# In[2]:


from load_utils import load_from_json
from path_utils import getpath
from tqdm.notebook import tqdm 


# In[3]:


fp = getpath('ORACLES')
fp_rtm = getpath('rtm')
fp_uvspec = getpath('uvspec')+'uvspec'
fp_rtmdat = fp_rtm+'dat/'
fp_uvspec_dat = getpath('uvspec_dat')


# In[6]:


if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\ORACLES\\'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\ORACLES\\rtm\\'
    fp_uvspec = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\bin\\uvspec'
    fp_rtmdat = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\data\\'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/ORACLES/'
    fp_rtm = '/nobackup/sleblan2/rtm/'
    fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
    fp_rtmdat = '/nobackup/sleblan2/AAC_DARF/rtm/' #'/u/sleblan2/4STAR/rtm_dat/'
else:
    raise Exception


# # Setup the variables used to create the lut

# In[4]:


vv = 'v7_irr'
mu = np.arange(1.02,3.4,0.05)
mu.shape


# In[5]:


sza = np.arccos(1.0/mu)*180.0/np.pi
#sza = np.arange(40,91,5)
print(sza)


# In[8]:


tau = np.array([0.1,0.2,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0,50.0,
       60.0,80.0,100.0])
#ref = np.append(np.append(np.arange(1,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,2.5)))
ref = np.arange(0.5,20.5,0.5)


# In[9]:


ref


# In[10]:


print(ref.shape)
print(tau.shape)


# In[169]:


pmom = RL.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',deltascale=True, new=True)


# In[170]:


pmom['file_name']


# In[171]:


aero = load_from_json(fp+'aero_file_v4.txt')


# In[172]:


aero['ext']


# In[173]:


if 'ext' in aero: print 'a'


# In[174]:


#geo = {'lat':-22.979,
#       'lon':14.645,
#       'doy':245,
#       'zout':[0.2,1.5,100.0]}
geo = {'lat':-16.0,
       'lon': 9.0,
       'doy':253, #September 10th
       'zout':[0.2,1.5,100.0]}
#aero = {'z_arr':[2.0,5.0],
#        'ext':np.array([[0.6,0.4,0.10,0.04],[0.0,0.0,0.0,0.0]]),
#        'ssa':np.array([[0.8,0.85,0.9,0.95],[0.9,0.9,0.9,0.9]]),
#        'asy':np.array([[0.8,0.8,0.8,0.8],[0.8,0.8,0.8,0.8]]),
#        'wvl_arr':[400.0,500.0,650.0,940.0],
#        'disort_phase':False,
#        'expand_hg':True}
cloud = {'ztop':1.0,
         'zbot':0.5,
         'write_moments_file':True,
         'moms_dict':pmom}
source = {'wvl_range':[350,1650],
          'source':'solar',
          'integrate_values':False,
          'run_fuliou':False,
          'dat_path':fp_uvspec_dat,
          'atm_file':fp_uvspec_dat + 'atmmod/afglt.dat',
          'zenith':False}
albedo = {'create_albedo_file':False,
          'sea_surface_albedo':True,
          'wind_speed':5.0}


# In[175]:


pmom['file_name']


# In[176]:


RL.print_version_details(fp+'ORACLES_lut_%s.txt'%vv,vv,geo=geo,
                         aero=aero,cloud=cloud,source=source,albedo=albedo,
                         tau=tau,ref=ref,sza=sza,cloud_pmom_file=pmom['file_name'])


# In[178]:


fp_in = os.path.join(fp_rtm,'input','%s_ORACLES'%vv)
fp_out = os.path.join(fp_rtm,'output','%s_ORACLES'%vv)


# In[179]:


f_slit_vis = os.path.join(fp_rtm,'vis_1nm.dat')
f_slit_nir = os.path.join(fp_rtm,'nir_1nm.dat')


# In[180]:


if not os.path.exists(fp_in):
    os.makedirs(fp_in)
if not os.path.exists(fp_out):
    os.makedirs(fp_out)


# In[21]:


f_list = open(os.path.join(fp_rtm,'ORACLES_list_%s.sh'%vv),'w')
print f_list.name


# In[22]:


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


# ## Run through typical method

# In[32]:


if isjupyter():
    pbar = tqdm(total=len(sza)*len(tau)*len(ref))
for s in sza:
    for t in tau:
        for r in ref:
            fname = 'lut_irr_sza%04.1f_tau%06.2f_ref%04.1f' % (s,t,r)
            geo['sza'] = s
            cloud['tau'] = t
            cloud['ref'] = r
            if False: #r>=5.0:
                cloud['phase'] = 'ic'
                fname0 = fname+'_'+cloud['phase']+'_w0.dat'
                source['wvl_range'] = [499.0,501.]
                source['slit_file'] = f_slit_vis
                RL.write_input_aac(os.path.join(fp_in,fname0),geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
                f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname0)+' > '+os.path.join(fp_out,fname0)+'\n')
                fname1 = fname+'_'+cloud['phase']+'_w1.dat'
                source['wvl_range'] = [981.,1700.]
                source['slit_file'] = f_slit_nir
                RL.write_input_aac(os.path.join(fp_in,fname1),geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
                f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname1)+' > '+os.path.join(fp_out,fname1)+'\n')
            if r<=20.0:
                cloud['phase'] = 'wc'
                fname0 = fname+'_'+cloud['phase']+'_w0.dat'
                source['wvl_range'] = [500.0,500.0]
                source['slit_file'] = f_slit_vis
                RL.write_input_aac(os.path.join(fp_in,fname0),geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
                f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname0)+' > '+os.path.join(fp_out,fname0)+'\n')
                fname1 = fname+'_'+cloud['phase']+'_w1.dat'
                source['wvl_range'] = [1650.0,1650.0]
                source['slit_file'] = f_slit_nir
                RL.write_input_aac(os.path.join(fp_in,fname1),geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
                f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname1)+' > '+os.path.join(fp_out,fname1)+'\n')                
            if isjupyter(): 
                pbar.update(1)
            else:
                print s,t,r
        #break
    #break


# In[25]:


f_list.close()


# ## Run for multiprocessing

# In[181]:


from multiprocessing import Pool, cpu_count
from copy import deepcopy
import signal
import warnings
warnings.simplefilter('ignore')


# In[182]:


def worker_init(verbose=True):
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[184]:


f_list = open(fp_rtm+'{}_{}.sh'.format(name,vv),'w')
print f_list.name
fpp_in = fp_rtm+'input/{}_{}/'.format(name,vv)
fpp_out = fp_rtm+'output/{}_{}/'.format(name,vv)


# In[185]:


if not os.path.exists(fp_in):
    os.makedirs(fp_in)
if not os.path.exists(fp_out):
    os.makedirs(fp_out)


# In[186]:


if isjupyter():
    pbar = tqdm(total=len(sza)*len(tau)*len(ref))
bb = []
for s in sza:
    for t in tau:
        for r in ref:
            fname = 'lut_irr_sza%04.1f_tau%06.2f_ref%04.1f' % (s,t,r)
            geo['sza'] = s
            cloud['tau'] = t
            cloud['ref'] = r
            
            cloud['phase'] = 'wc'
            fname0 = fname+'_'+cloud['phase']+'_w0.dat'
            source['wvl_range'] = [500.0,500.0]
            source['slit_file'] = f_slit_vis
            #RL.write_input_aac(os.path.join(fp_in,fname0),geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
            #                   verbose=False,make_base=False,set_quiet=True)
            f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname0)+' > '+os.path.join(fp_out,fname0)+'\n')
            fname1 = fname+'_'+cloud['phase']+'_w1.dat'
            source['wvl_range'] = [1650.0,1650.0]
            source['slit_file'] = f_slit_nir
            #RL.write_input_aac(os.path.join(fp_in,fname1),geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
            #                   verbose=False,make_base=False,set_quiet=True)
            f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname1)+' > '+os.path.join(fp_out,fname1)+'\n')
            
            bb.append({'geo':deepcopy(geo),'cod':t,'ref':r,
                       'f0':deepcopy(fname0),'f1':deepcopy(fname1),'w0':[500.0,500.0],'w1':[1650.0,1650.0],
                       's0':f_slit_vis,'s1':f_slit_nir})
            
            if isjupyter(): 
                pbar.update(1)
            else:
                print s,t,r
f_list.close()


# In[187]:


def write_files(d,cloud=cloud,source=source,albedo=albedo,aero=aero,fp_in=fp_in):
    'function to feed the pool of workers to write out the all the files'
    cloud['tau'],cloud['ref'] = d['cod'],d['ref']
    source['wvl_range'] = d['w0']
    source['slit_file'] = d['s0']
    RL.write_input_aac(fp_in+'/'+d['f0'],geo=d['geo'],aero={},cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
    source['wvl_range'] = d['w1']
    source['slit_file'] = d['s1']
    RL.write_input_aac(fp_in+'/'+d['f1'],geo=d['geo'],aero={},cloud=cloud,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)


# In[188]:


p = Pool(cpu_count()-1,worker_init)


# In[189]:


len(bb)


# In[190]:


results = []
max_ = len(bb)
with tqdm(total=max_) as pbar:
    for i, res in tqdm(enumerate(p.imap_unordered(write_files, bb))):
        pbar.update()
        results.append(res)


# In[191]:


results[0]


# # Run the calculations

# In[193]:


f_list = fp_rtm+'{}_{}.sh'.format(name,vv)


# In[194]:


get_ipython().system(u' wc -l $f_list')


# In[195]:


f_listout = f_list+'.out'


# In[196]:


get_ipython().system(u'parallel --jobs=22 --bar < $f_list #2> $f_listout')


# In[ ]:




