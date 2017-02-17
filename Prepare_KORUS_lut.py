
# coding: utf-8

# # Info
# Name:  
# 
#     Prepare_KORUS_lut
# 
# Purpose:  
# 
#     Create the input libradtran files for creating a lut of low clouds with aerosol on top to be used in KORUS operational cloud retrievals
# 
# Calling Sequence:
# 
#     python Prepare_KORUS_lut
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
#     Written: Samuel LeBlanc, Bathurst, NB, 2017-01-05
#              based on Prepare_NAAMES_lut
#     

# # Prepare the python environment

# In[3]:

import numpy as np
import scipy.io as sio
import os
import Run_libradtran as RL
reload(RL)


# In[ ]:

from load_utils import load_from_json


# In[ ]:

name = 'KORUS'


# In[68]:

if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\{}\\'.format(name)
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\{}\\rtm\\'.format(name)
    fp_uvspec = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\bin\\uvspec'
    fp_rtmdat = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\data\\'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/{}/'.format(name)
    fp_rtm = '/nobackup/sleblan2/rtm/'
    fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
    fp_rtmdat = '/nobackup/sleblan2/AAC_DARF/rtm/' #'/u/sleblan2/4STAR/rtm_dat/'
else:
    raise Exception


# # Setup the variables used to create the lut

# In[84]:

vv = 'v3'
mu = np.arange(1.05,4.0,0.2)
mu.shape


# In[139]:

sza = np.round(np.arccos(1.0/mu)*180.0/np.pi)
#sza = np.arange(40,91,5)
print(sza)


# In[108]:

tau = np.array([0.1,0.2,0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,
       6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0,50.0,
       60.0,80.0,100.0])
ref = np.append(np.append(np.arange(1,15),np.arange(15,30,2)),np.ceil(np.arange(30,61,2.5)))


# In[109]:

ref


# In[112]:

print(ref.shape)
print(tau.shape)


# In[ ]:

pmom = RL.make_pmom_inputs(fp_rtm=fp_rtmdat,source='solar',cloudtype='ic')


# In[ ]:

#aero = load_from_json(fp+'aero_save.txt')


# In[4]:

geo = {'lat':36.9921,
       'lon':127.1129,
       'doy':245,
       'zout':[0.2,8.0,100.0],
       'name':'Pyongtaek'}
aero = {'z_arr':[2.0,5.0],
        'ext':np.array([[0.16837,0.12837, 0.08727, 0.0655, 0.0557, 0.0357, 0.0157],[0.0,0.0,0.0,0.0,0.0,0.0,0.0]]), 
        #tau = 0.385100,0.261800,0.196500,0.167100
        'ssa':np.array([[0.939000,0.937000,0.934400,0.930700,0.929600,0.925600,0.920600],[0.9,0.9,0.9,0.9,0.9,0.9,0.9]]),
        'asy':np.array([[0.680998,0.660998,0.629995,0.612484,0.606010,0.596010,0.580010],[0.8,0.8,0.8,0.8,0.8,0.8,0.8]]), 
        'wvl_arr':[400.0,442.0,668.0,870.0,1020.0,1200.0,1700.0], #442,668,870,1020
        'disort_phase':False,
        'expand_hg':True,
        'details':'From AERONET Socheongcho on May 20 2017, 23:14 UTC, extrapolated by hand'}
#aero = {}
cloud = {'ztop':12.0,
         'zbot':10.0,
         'write_moments_file':True,
         'moms_dict':pmom} # for cirrus mostly
source = {'wvl_range':[350,1750],
          'source':'solar',
          'integrate_values':False,
          'run_fuliou':False,
          'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/',
          'atm_file':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/atmmod/afglmw.dat',
          'zenith':True}
albedo = {'create_albedo_file':False,
          'sea_surface_albedo':True,
          'wind_speed':5.0}


# In[60]:

RL.print_version_details(fp+'{name}_lut_{vv}.txt'.format(name=name,vv=vv),vv,geo=geo,
                         aero=aero,cloud=cloud,source=source,albedo=albedo,tau=tau,ref=ref,sza=sza,
                         cloud_pmom_file=cloud['moms_dict']['file_name'])


# In[71]:

fp_in = os.path.join(fp_rtm,'input','{vv}_{name}'.format(vv=vv,name=name))
fp_out = os.path.join(fp_rtm,'output','{vv}_{name}'.format(vv=vv,name=name))


# In[82]:

f_slit_vis = os.path.join(fp_rtm,'4STAR_vis_slit_1nm.dat')
f_slit_nir = os.path.join(fp_rtm,'4STAR_nir_slit_1nm.dat')


# In[72]:

if not os.path.exists(fp_in):
    os.makedirs(fp_in)
if not os.path.exists(fp_out):
    os.makedirs(fp_out)


# In[79]:

f_list = open(os.path.join(fp,'run','{name}_list_{vv}.sh'.format(vv=vv,name=name)),'w')
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
            if False: #r<=30.0:
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

