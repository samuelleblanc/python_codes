
# coding: utf-8

# # Intro
# Name:  
# 
#     Prepare_ARISE_lut
# 
# Purpose:  
# 
#     Create the input libradtran files for creating a lut of low clouds to be used in ARISE cloud retrievals near sea ice edge
# 
# Calling Sequence:
# 
#     python Prepare_ARISE_lut
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2016-10-18
#     

# # Prepare the python environment

# In[3]:


import numpy as np
import scipy.io as sio
import os
import Run_libradtran as RL


# In[12]:


from load_utils import load_from_json


# In[68]:


if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\ARISE\\'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\ARISE\\rtm\\'
    fp_uvspec = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\bin\\uvspec'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/ARISE/'
    fp_rtm = '/nobackup/sleblan2/rtm/'
    fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
else:
    raise Exception


# # Setup the variables used to create the lut

# In[15]:


vv = 'v3_wat1'
mu = np.arange(3.0,4.0,0.1)
mu.shape


# In[17]:


sza = np.arccos(1.0/mu)*180.0/np.pi
#sza = np.arange(40,91,5)
print(sza)


# In[18]:


tau = np.array([0.1,0.2,0.3,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.3,2.6,3.0,3.5,4.0,4.5,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,22.0,25.0,27.0,30.0,35.0,40.0,50.0,60.0])
ref = np.append(np.append(np.arange(2,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,5.0)))


# In[19]:


ref


# In[20]:


print(ref.shape)
print(tau.shape)


# In[22]:


geo = {'lat':72.02,
       'lon':129.3,
       'doy':262,
       'zout':[0.1,2.0,100.0]}
aero = {} # none
cloud = {'ztop':1.0,
         'zbot':0.5,
         'write_moments_file':False}
source = {'wvl_range':[400,1750],
          'source':'solar',
          'integrate_values':False,
          'run_fuliou':False,
          'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/',
          'atm_file':'/nobackup/sleblan2/dat/atmos_20140919.dat',
          'zenith':True}
albedo = {'create_albedo_file':False,
          'sea_surface_albedo':False,
          'albedo_file':'/nobackup/sleblan2/dat/albedo_v3_20140919_ice.dat'}


# In[23]:


if vv is 'v3_wat1':
    cloud['ztop'] = 2.0
    cloud['zbot'] = 0.5
    albedo['albedo_file'] = '/nobackup/sleblan2/dat/20140919_surf_alb_water1_pure.dat'
elif vv is 'v3_wat2':
    cloud['ztop'] = 1.2
    cloud['zbot'] = 0.3
    albedo['albedo_file'] = '/nobackup/sleblan2/dat/20140919_surf_alb_water2_inter.dat'
elif vv is 'v3_ice_top':
    cloud['ztop'] = 0.5
    cloud['zbot'] = 0.1
    albedo['albedo_file'] = '/nobackup/sleblan2/dat/20140919_surf_alb_ice_top.dat'
elif vv is 'v3_ice_mid':
    cloud['ztop'] = 0.8
    cloud['zbot'] = 0.1
    albedo['albedo_file'] = '/nobackup/sleblan2/dat/20140919_surf_alb_ice_mid.dat'
elif vv is 'v3_ice_low':
    cloud['ztop'] = 1.1
    cloud['zbot'] = 0.2
    albedo['albedo_file'] = '/nobackup/sleblan2/dat/20140919_surf_alb_ice_low.dat'
    
    


# In[24]:


RL.print_version_details(fp+'ARISE_lut_%s.txt'%vv,vv,geo=geo,
                         aero=aero,cloud=cloud,source=source,albedo=albedo,tau=tau,ref=ref,sza=sza)


# In[71]:


fp_in = os.path.join(fp_rtm,'input','%s_ARISE'%vv)
fp_out = os.path.join(fp_rtm,'output','%s_ARISE'%vv)


# In[82]:


f_slit_vis = os.path.join(fp_rtm,'4STAR_vis_slit_1nm.dat')
f_slit_nir = os.path.join(fp_rtm,'4STAR_nir_slit_1nm.dat')


# In[72]:


if not os.path.exists(fp_in):
    os.makedirs(fp_in)
if not os.path.exists(fp_out):
    os.makedirs(fp_out)


# In[79]:


f_list = open(os.path.join(fp,'run','ARISE_list_%s.sh'%vv),'w')
print f_list.name


# In[ ]:


for s in sza:
    for t in tau:
        for r in ref:
            fname = 'lut_sza%04.1f_tau%06.2f_ref%04.1f' % (s,t,r)
            geo['sza'] = s
            cloud['tau'] = t
            cloud['ref'] = r
            if False: #r>=5.0:
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
            print s,t,r


# In[ ]:


f_list.close()

