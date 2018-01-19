
# coding: utf-8

# # Introduction
# Name:  
# 
#     DARE_AAC_spread
# 
# Purpose:  
# 
#     For Meloë's study for global DARE, using calipso, MODIS, and other measurements of Aerosol above clouds
#     Prep and read the radiative transfer files to run the 
#   
# Input:
# 
#     none
# 
# Output:
#    
#     plots
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - scipy : for saving and reading
#     - Run_libradtran
# 
#   
# Needed Files:
# 
#   - matlab input files: Input_to_DARF_mmm.mat
#   - spread of the data information

# # Import the required modules and file paths

# In[4]:


import numpy as np
import scipy.io as sio
import Run_libradtran as RL
import load_utils as lm
import os


# In[5]:


fp = '/u/sleblan2/meloe_AAC/v4_spread/'
fp_alb = '/nobackup/sleblan2/AAC_DARF/surface_albedo/'
fp_out = '/nobackup/sleblan2/AAC_DARF/input/v4_spread/'
fp_pmom = '/nobackup/sleblan2/AAC_DARF/rtm/'
fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
wvl_thm = '/nobackup/sleblan2/AAC_DARF/rtm/wvl.thm.dat'
vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'


# In[6]:


std_label = 'v4_spread'


# ## Check if reading or writing

# In[ ]:


import argparse
long_description = """    Prepare or read the AAC DARE calculations, with the spread for all values.
      - defaults to prepare AAC files """
parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-doread','--doread',help='if set, will only read the output, not produce them',
                    action='store_true')
parser.add_argument('-dowrite','--dowrite',help='if set, will write the input and list files for fuliou',
                    action='store_true')
#parser.add_argument('-i','--index',help='Sets the index of the pixel file to use (19374,22135). Default is 0',type=int)
in_ = vars(parser.parse_args())
doread = in_.get('doread',False)
dowrite = in_.get('dowrite',True)
#i = in_.get('index',0)


# # Load the input files and create the subsets

# In[7]:


mmm = 'JJA'


# In[8]:


fpm = fp+'Input_to_DARF_{mmm}_{vv}.mat'.format(mmm=mmm,vv=version)
print 'in %s months, getting mat file: %s' % (mmm,fpm)


# In[ ]:


input_mmm = sio.loadmat(fpm,mat_dtype=True)['data_input_darf']


# In[10]:


geo = {'zout':[0,3,100],'year':2007,'day':15,'minute':0,'second':0}
geo['month'] = 7
doy = 225


# In[11]:


geo['lat'],geo['lon'] = -7.0,-10.0


# In[19]:


ilat = range(18,26)
ilon = range(37,40)


# In[ ]:


cod = input_mmm['MODIS_COD_mean'][0,0][ilat,ilon]
ext = np.abs(input_mmm['MOC_ext_mean'][0,0][ilat,ilon,:])
ssa = input_mmm['MOC_ssa_mean'][0,0][ilat,ilon,:]
asy = input_mmm['MOC_asym_mean'][0,0][ilat,ilon,:]


# In[ ]:


cloud['ref'] = np.nanmean(input_mmm['MODIS_effrad_mean'][0,0][ilat,ilon])


# # Prepare the inputs

# In[ ]:


pmom_solar = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='solar')
pmom_thermal = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='thermal')
max_nmom=20
pmom_solar['max_nmoms'] = max_nmom
pmom_thermal['max_nmoms'] = max_nmom


# In[ ]:


aero = {'z_arr':[3.0,4.0]}
cloud = {'ztop':3.0,'zbot':2.0,'phase':'wc','write_moments_file':True}
source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/','run_fuliou':True}
albedo = {'create_albedo_file':False}
albedo['sea_surface_albedo'] = True


# # Start the writing out

# In[ ]:


change_fp_output = True
if aero_clear:
    std_label = '_clear'


# In[ ]:


fp_out2 = fp_out+mmm+std_label+'/'
if not os.path.exists(fp_out2):
    os.mkdir(fp_out2)
if change_fp_output:
    fp_output = fp_out2.replace('input','output')
    if not os.path.exists(fp_output):
        os.mkdir(fp_output)
fp_base_file = fp_out2+'base.inp'
make_base = True


# In[ ]:


if dowrite: 
    ff = 'AAC_list_file_{m}_{v}{lbl}.sh'.format(m=mmm,v=version,lbl=std_label)
    file_list = file(fp_out+ff,'w')
    print 'Starting list file: '+fp_out+ff


# In[ ]:


aero['wvl_arr'] = input_mmm['MOC_wavelengths'][0,0][0,:]*1000.0


# In[ ]:


print 'Running through the permutations'
for icod,c in enumerate(cod):
    for iext,e in enumerate(ext):
        for issa,s in enumerate(ssa):
            for iasy,a in enumerate(asy):
                # set the aerosol values
                aero['ext'] = e
                aero['ext'][aero['ext']<0.0] = 0.0
                if np.isnan(aero['ext']).all():
                    print 'skipping cod:%i, ext:%i, ssa:%i, asy:%i' % (icod,iext,issa,iasy)
                    continue
                aero['ssa'] = s
                aero['asy'] = a

                #sanitize inputs after adding subtracting standard deviations
                try: aero['ssa'][aero['ssa']<0.0] = 0.0
                except: pass
                try: aero['ssa'][aero['ssa']>1.0] = 1.0
                except: pass
                try: aero['asy'][aero['asy']<0.0] = 0.0
                except: pass
                try: aero['asy'][aero['asy']>1.0] = 1.0
                except: pass

                if aero['wvl_arr'].max()<100000.0:
                    aero['wvl_arr'] = np.append(aero['wvl_arr'],100000.0)
                    aero['ext'] = np.append(aero['ext'],e[-1])
                    aero['ssa'] = np.append(aero['ssa'],s[-1])
                    aero['asy'] = np.append(aero['asy'],a[-1])
                # set the cloud values
                cloud['tau'] = c
                try: cloud['tau'][cloud['tau']<0.0] = 0.0
                except: pass
                try: cloud['ref'][cloud['ref']<2.0] = 2.0
                except: pass

                cloud['link_to_mom_file'] = False
                aero['link_to_mom_file'] = False
                cloud_file_name_sol = fp_out2+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_sol.inp_cloud' % (icod,iext,issa,iasy,mmm)
                cloud_file_name_thm = fp_out2+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_thm.inp_cloud' % (icod,iext,issa,iasy,mmm)
                aero['file_name'] = fp_out2+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_sol.inp_aero' % (icod,iext,issa,iasy,mmm)

                for HH in xrange(24):
                    geo['hour'] = HH
                    #build the solar input file
                    source['source'] = 'solar'
                    if wvl_file_sol:
                        source['wvl_filename'] = wvl_file_sol
                    else:
                        source['wvl_range'] = [250,5600]
                        source['wvl_filename'] = None
                    cloud['moms_dict'] = pmom_solar
                    cloud['file_name'] = cloud_file_name_sol
                    file_out_sol = fp_out2+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_sol.inp' % (icod,iext,issa,iasy,mmm,HH)
                    if not list_only:
                        RL.write_input_aac(file_out_sol,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                                       make_base=make_base,fp_base_file=fp_base_file,set_quiet=True,solver='rodents')
                    if make_base:
                        make_base = False
                    #build the thermal input file
                    source['source'] = 'thermal'
                    if wvl_file_thm:
                        source['wvl_filename'] = wvl_file_thm
                    else:
                        source['wvl_range'] = [4000,50000-1]
                        source['wvl_filename'] = None
                    cloud['moms_dict'] = pmom_thermal
                    cloud['file_name'] = cloud_file_name_thm
                    file_out_thm = fp_out2+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_thm.inp' % (icod,iext,issa,iasy,mmm,HH)

                    if not list_only:
                        RL.write_input_aac(file_out_thm,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                                       make_base=False,fp_base_file=fp_base_file,set_quiet=True,solver='rodents')
                    file_list.write(fp_uvspec+' < '+file_out_sol+' > '+fp_output
                                    +'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_sol.out\n' % (icod,iext,issa,iasy,mmm,HH))
                    file_list.write(fp_uvspec+' < '+file_out_thm+' > '+fp_output
                                    +'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_thm.out\n' % (icod,iext,issa,iasy,mmm,HH))
                    if not cloud['link_to_mom_file']:
                        cloud['link_to_mom_file'] = True
                    if not aero['link_to_mom_file']:
                        aero['link_to_mom_file'] = True
                    print mmm,icod,iext,issa,iasy,HH
file_list.close()

