
# coding: utf-8

# Name:  
# 
#     Prepare_NAAMES_lut
# 
# Purpose:  
# 
#     Create the input libradtran files for creating a lut of low clouds to be used in NAAMES operational cloud retrievals
# 
# Calling Sequence:
# 
#     python Prepare_NAAMES_lut
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
#     Written: Samuel LeBlanc, NASA Ames, 2015-10-06
#     

# # Prepare the python environment

# In[66]:

import numpy as np
import scipy.io as sio
import os
import Run_libradtran as RL
reload(RL)


# In[68]:

if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\NAAMES\\'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\NAAMES\\rtm\\'
    fp_uvspec = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\bin\\uvspec'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/NAAMES/'
    fp_rtm = '/nobackup/sleblan2/rtm/'
    fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
else:
    raise Exception


# # Setup the variables used to create the lut

# In[84]:

vv = 'v2'
mu = np.arange(1.05,4.0,0.2)
mu.shape


# In[139]:

sza = np.round(np.arccos(1.0/mu)*180.0/np.pi)
#sza = np.arange(40,91,5)
print(sza)


# In[108]:

tau = np.array([0.1,0.2,0.3,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0,40.0,50.0,
       60.0,70.0,80.0,100.0])
ref = np.append(np.append(np.arange(2,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,2.5)))


# In[109]:

ref


# In[112]:

print(ref.shape)
print(tau.shape)


# In[70]:

geo = {'lat':47.6212167,
       'lon':52.74245,
       'doy':322,
       'zout':[0.2,3.0,100.0]}
aero = {} # none
cloud = {'ztop':1.5,
         'zbot':0.5,
         'write_moments_file':False}
source = {'wvl_range':[350,1750],
          'source':'solar',
          'integrate_values':False,
          'run_fuliou':False,
          'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/',
          'atm_file':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/atmmod/afglmw.dat',
          'zenith':True}
albedo = {'create_albedo_file':False,
          'sea_surface_albedo':True,
          'wind_speed':14.0}


# In[60]:

RL.print_version_details(fp+'NAAMES_lut_%s.txt'%vv,vv,geo=geo,
                         aero=aero,cloud=cloud,source=source,albedo=albedo,tau=tau,ref=ref,sza=sza)


# In[71]:

fp_in = os.path.join(fp_rtm,'input','%s_NAAMES'%vv)
fp_out = os.path.join(fp_rtm,'output','%s_NAAMES'%vv)


# In[82]:

f_slit_vis = os.path.join(fp_rtm,'4STAR_vis_slit_1nm.dat')
f_slit_nir = os.path.join(fp_rtm,'4STAR_nir_slit_1nm.dat')


# In[72]:

if not os.path.exists(fp_in):
    os.makedirs(fp_in)
if not os.path.exists(fp_out):
    os.makedirs(fp_out)


# In[79]:

f_list = open(os.path.join(fp,'run','NAAMES_list_%s.sh'%vv),'w')
print f_list.name


# In[ ]:

for s in sza:
    for t in tau:
        for r in ref:
            fname = 'lut_sza%02i_tau%06.2f_ref%04.1f' % (s,t,r)
            geo['sza'] = s
            cloud['tau'] = t
            cloud['ref'] = r
            if r>=5.0:
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

