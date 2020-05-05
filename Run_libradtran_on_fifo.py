#!/usr/bin/env python
# coding: utf-8

# # Import modues

# In[18]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
import scipy.io as sio
from path_utils import getpath
import matplotlib.pyplot as plt
import load_utils as lu
from write_utils import nearest_neighbor, iterate_dict_unicode
from tqdm.notebook import tqdm 
from datetime import datetime
from multiprocessing import Pool, cpu_count
from copy import deepcopy
import signal
import warnings
import Run_libradtran as Rl


# In[16]:


get_ipython().magic(u'matplotlib notebook')


# In[17]:


warnings.simplefilter('ignore')


# In[19]:


name = 'ORACLES'
vv = 'v3'
vr = 'R3'


# In[20]:


fp = getpath(name)
fp_rtm = getpath('rtm')
fp_uvspec = getpath('uvspecb')+'uvspec'
matfile = fp+'{}_all_cld_ict.mat'.format(vr)
fp_uvspec_dat = getpath('uvspec_dat') 
fp_rtmdat = fp_rtm+'dat/'


# ### Test out using fifo instead of files

# In[3]:


import os
from multiprocessing import Process, Event
import subprocess as sub


# In[51]:


fp_fifo_in = '/tmp/uvspec_input.fifo'


# In[5]:


os.mkfifo(fp_fifo_in)


# In[8]:


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


# In[9]:


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


# In[10]:


def run():
    process = sub.Popen([fp_uvspec],stdin=p, stdout=sub.PIPE,stderr=sub.PIPE)
    stdout,stderr = process.communicate()
    #stderr = process.stderr.read()
    print 'STDOUT:{},{},{}'.format(stdout,stderr,process.poll())


# In[11]:


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


# In[38]:


def write_xtrf(fp_fifo_in):
    if not os.path.exists(fp_fifo_in):
        os.mkfifo(fp_fifo_in)
    p = open(fp_fifo_in,'w')
    p.flush()
    g = ['# wvl[nm]    alb[unitless]','250.000000      -0.043171','350.000000      -0.010611','400.000000      0.005669','500.000000      0.038229','675.000000      0.058627','870.000000      0.229436','995.000000      0.234727','1200.000000     0.240584','1400.000000     0.246298','1600.000000     0.252013','2100.000000     0.266298','3200.000000     0.297727','4900.000000     0.346298']
    for llj in g:
        p.write('{}\n'.format(llj))
    #p.flush()
    p.close()
    os.unlink(fp_fifo_in)
    #os.remove(fp_fifo_in)


# In[108]:


write_xtrf(fp_fifo_in)


# In[85]:


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


# In[ ]:





# In[ ]:




