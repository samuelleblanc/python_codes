#!/usr/bin/env python
# coding: utf-8

# # Info
# Name:  
# 
#     Prepare_SASZe_lut
# 
# Purpose:  
# 
#     Create the input libradtran files for creating a lut of low clouds with aerosol on top to be used in SaS-Ze operational cloud retrievals
#     This is for test on SGP surface albedo
# 
# Calling Sequence:
# 
#     python Prepare_SASZe_lut
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
#   - ...
#     
# History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-10-17
#              based on Prepare_KORUS_lut from 2017-01-05
#     

# # Prepare the python environment

# In[5]:


import numpy as np
import scipy.io as sio
import os
import Run_libradtran as RL
from path_utils import getpath
reload(RL)


# In[2]:


from load_utils import load_from_json


# In[3]:


name = 'SASZe'
vv = 'v1'


# In[15]:


fp = getpath(name)
fp_rtm = getpath('rtm')
fp_uvspec = getpath('uvspec')+'uvspec'
fp_rtmdat = getpath('rtm')+'dat/'


# # Setup the variables used to create the lut

# In[16]:


mu = np.arange(1.05,4.0,0.2)
mu.shape


# In[17]:


sza = np.round(np.arccos(1.0/mu)*180.0/np.pi)
#sza = np.arange(40,91,5)
print(sza)


# In[18]:


tau = np.array([0.1,0.2,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0,50.0,
       60.0,80.0,100.0])
ref = np.append(np.append(np.arange(1,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,2.5)))


# In[19]:


ref


# In[20]:


print(ref.shape)
print(tau.shape)


# In[22]:


pmom_ic = RL.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',cloudtype='ic')


# In[23]:


pmom_wc = RL.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',cloudtype='wc')


# In[ ]:


#aero = load_from_json(fp+'aero_save.txt')


# In[32]:


#load SGP average spectral albedo
# used to build original LeBlanc et al., 2015 Boulder surtface albedo
# measured by Michalsky et al. (2003), doi:10.1029/2002JD002906

alb_sgp = sio.idl.readsav('/data/sam/SSFR3/surface_albedo/alb.out') 


# In[48]:


geo = {'lat':36.6103,
       'lon':-97.488,
       'doy':120,
       'zout':[0.2,8.0,100.0],
       'name':'SGP'}
aero = {}
cloud_wc = {'ztop':3.0,
         'zbot':1.0,
         'write_moments_file':True,
         'moms_dict':pmom_wc} # for low water cloud mostly
cloud_ic = {'ztop':12.0,
         'zbot':10.0,
         'write_moments_file':True,
         'moms_dict':pmom_ic} # for high ice cloud mostly
source = {'wvl_range':[355,1750],
          'source':'solar',
          'integrate_values':False,
          'run_fuliou':False,
          'dat_path':getpath('uvspec') + '../data/',
          'atm_file':getpath('uvspec') + '../data/atmmod/afglmw.dat',
          'zenith':True}
albedo = {'create_albedo_file':True,
          'sea_surface_albedo':False,
          'alb':alb_sgp['alb'],
          'alb_wvl':alb_sgp['wvl']}


# In[38]:


RL.print_version_details(fp+'{name}_lut_{vv}.txt'.format(name=name,vv=vv),vv,geo=geo,
                         aero=aero,cloud=cloud,source=source,albedo=albedo,tau=tau,ref=ref,sza=sza,
                         cloud_pmom_file=cloud['moms_dict']['file_name'])


# In[39]:


fp_in = os.path.join(fp_rtm,'input','{vv}_{name}'.format(vv=vv,name=name))
fp_out = os.path.join(fp_rtm,'output','{vv}_{name}'.format(vv=vv,name=name))


# In[58]:


f_slit_vis = os.path.join(fp_rtm,'4STAR_vis_slit_1nm.dat')
f_slit_nir = os.path.join(fp_rtm,'nir_1nm.dat')


# In[59]:


if not os.path.exists(fp_in):
    os.makedirs(fp_in)
if not os.path.exists(fp_out):
    os.makedirs(fp_out)


# In[67]:


f_list = open(os.path.join(fp_rtm,'run','{name}_list_{vv}.sh'.format(vv=vv,name=name)),'w')
print f_list.name


# In[63]:


from tqdm import tqdm_notebook as tqdm


# In[68]:


after_first_sza = False
pbar = tqdm(total = len(sza)*len(tau)*len(ref)+1)
for s in sza:
    for t in tau:
        for r in ref:
            fname = 'lut_sza%02i_tau%06.2f_ref%04.1f' % (s,t,r)
            geo['sza'] = s
            cloud_ic['tau'] = t
            cloud_ic['ref'] = r
            cloud_wc['tau'] = t
            cloud_wc['ref'] = r
            if after_first_sza: cloud['link_to_mom_file'] = True
            if r>=25.0:
                cloud_ic['phase'] = 'ic'
                fname0 = fname+'_'+cloud_ic['phase']+'_w0.dat'
                source['wvl_range'] = [400.,981.]
                source['slit_file'] = f_slit_vis
                if after_first_sza: 
                    cloud_ic['file_name'] = os.path.join(fp_in,'lut_sza{s:02.0f}_tau{t:06.2f}_ref{r:04.1f}_{p}_w0.dat_cloud'.format(
                                          s=sza[0],t=t,r=r,p=cloud['phase']))
                RL.write_input_aac(os.path.join(fp_in,fname0),geo=geo,aero=aero,cloud=cloud_ic,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
                f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname0)+' > '+os.path.join(fp_out,fname0)+'\n')
                fname1 = fname+'_'+cloud_ic['phase']+'_w1.dat'
                source['wvl_range'] = [981.,1700.]
                source['slit_file'] = f_slit_nir
                if after_first_sza: cloud_ic['file_name'] = cloud_ic['file_name'].replace('w0','w1')
                RL.write_input_aac(os.path.join(fp_in,fname1),geo=geo,aero=aero,cloud=cloud_ic,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
                f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname1)+' > '+os.path.join(fp_out,fname1)+'\n')
            if r<=30.0:
                cloud_wc['phase'] = 'wc'
                fname0 = fname+'_'+cloud_wc['phase']+'_w0.dat'
                source['wvl_range'] = [400.,981.]
                source['slit_file'] = f_slit_vis
                if after_first_sza: 
                    cloud_wc['file_name'] = os.path.join(fp_in,'lut_sza{s:02.0f}_tau{t:06.2f}_ref{r:04.1f}_{p}_w0.dat_cloud'.format(
                                          s=sza[0],t=t,r=r,p=cloud_wc['phase']))
                RL.write_input_aac(os.path.join(fp_in,fname0),geo=geo,aero=aero,cloud=cloud_wc,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
                f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname0)+' > '+os.path.join(fp_out,fname0)+'\n')
                fname1 = fname+'_'+cloud_wc['phase']+'_w1.dat'
                source['wvl_range'] = [981.,1700.]
                source['slit_file'] = f_slit_nir
                if after_first_sza: cloud_wc['file_name'] = cloud_wc['file_name'].replace('w0','w1')
                RL.write_input_aac(os.path.join(fp_in,fname1),geo=geo,aero=aero,cloud=cloud_wc,source=source,albedo=albedo,
                                   verbose=False,make_base=False,set_quiet=True)
                f_list.write(fp_uvspec+' < '+os.path.join(fp_in,fname1)+' > '+os.path.join(fp_out,fname1)+'\n')                
            #print s,t,r
            pbar.update()
    after_first_sza = True


# In[66]:


f_list.close()

