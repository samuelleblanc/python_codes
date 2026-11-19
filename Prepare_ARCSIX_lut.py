#!/usr/bin/env python
# coding: utf-8

# # Intro
# Name:  
# 
#     Prepare_ARCSIX_lut
# 
# Purpose:  
# 
#     Create the input libradtran files for creating a lut of low clouds to be used in ARCSIX June 7th
# 
# Calling Sequence:
# 
#     python Prepare_ARCSIX_lut
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
#     - mpl_toolkits for basemap, map plotting
#     - pdb
#     - datetime
# 
#   
# Needed Files:
# 
#   - atmospheric profile file
#   - surface albedo file
#   - cloud mie scattering properties
#     
# History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2026-11-13
#              - from Prepare_ARISE_lut
#     

# # Prepare the python environment

# In[1]:


import numpy as np
import scipy.io as sio
import os
import Run_libradtran as RL


# In[2]:


from load_utils import load_from_json


# In[35]:


fp = '/data2/ARCSIX/data/'
fp_atm = '/home/sam/libradtran/libRadtran-2.0.6/data/atmmod/'
fp_out = '/scratch/rtm/input/ARCSIX_v1/'
fp_rtm = '/scratch/rtm/' #input/ARCSIX_v1/'
fp_uvspec = '/home/sam/libradtran/libRadtran-2.0.6/bin/uvspec'


# # Setup the variables used to create the lut

# In[17]:


vv = 'v1'
mu = np.arange(2.0,2.4,0.1)
mu.shape


# In[18]:


sza = np.arange(60.0,63.0,0.5)


# In[19]:


sza


# In[20]:


tau = np.array([0.1,0.2,0.3,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.3,2.6,3.0,3.5,4.0,4.5,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,22.0,25.0,27.0,30.0,35.0,40.0,50.0,60.0])
ref = np.append(np.append(np.arange(2,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,5.0)))


# In[21]:


ref


# In[22]:


print(ref.shape)
print(tau.shape)


# In[36]:


geo = {'lat':83.48,
       'lon':-42.102,
       'doy':158,
       'zout':[0.1,2.0,100.0]}
aero = {} # none
cloud = {'ztop':0.5,
         'zbot':0.15,
         'write_moments_file':False}
source = {'wvl_range':[400,1750],
          'source':'solar',
          'integrate_values':False,
          'run_fuliou':False,
          'dat_path':'/home/sam/libradtran/libRadtran-2.0.6/data/',
          'atm_file':'/data2/ARCSIX/data/atm_20240607.dat',
          'zenith':True}
albedo = {'create_albedo_file':False,
          'sea_surface_albedo':False,
          'albedo_file':'/data2/ARCSIX/data/alb_20240607.dat'}


# In[37]:


RL.print_version_details(fp+'ARCSIX_lut_%s.txt'%vv,vv,geo=geo,
                         aero=aero,cloud=cloud,source=source,albedo=albedo,tau=tau,ref=ref,sza=sza)


# In[38]:


fp_in = os.path.join(fp_rtm,'input','ARCSIX_%s'%vv)
fp_out = os.path.join(fp_rtm,'output','ARCSIX_%s'%vv)


# In[39]:


fp_in


# In[40]:


f_slit_vis = os.path.join(fp_rtm,'vis_1nm.dat')
f_slit_nir = os.path.join(fp_rtm,'nir_1nm.dat')


# In[41]:


if not os.path.exists(fp_in):
    os.makedirs(fp_in)
if not os.path.exists(fp_out):
    os.makedirs(fp_out)


# In[42]:


f_list = open(os.path.join(fp_rtm,'run','ARCSIX_list_%s.sh'%vv),'w')
print(f_list.name)


# In[43]:


for s in sza:
    for t in tau:
        for r in ref:
            fname = 'lut_sza%04.1f_tau%06.2f_ref%04.1f' % (s,t,r)
            geo['sza'] = s
            cloud['tau'] = t
            cloud['ref'] = r
            if r>=10.0:
                cloud['phase'] = 'ic'
                fname0 = fname+'_'+cloud['phase']+'_w0.dat'
                source['wvl_range'] = [400.,981.]
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
            if r<=30.0:
                cloud['phase'] = 'wc'
                fname0 = fname+'_'+cloud['phase']+'_w0.dat'
                source['wvl_range'] = [400.,981.]
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
            print(s,t,r)


# In[44]:


f_list.close()


# In[ ]:




