#!/usr/bin/env python
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

# In[3]:


import numpy as np
import scipy.io as sio
import Run_libradtran as RL
import load_utils as lm
import os
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import signal
from copy import copy,deepcopy
import os
import warnings


# In[ ]:


warnings.simplefilter('ignore')


# In[20]:


std_label = 'v4_spread'
verbose = False
imax = 36


# In[21]:


fp = '/u/sleblan2/meloe_AAC/' + std_label + '/'
fp_alb = '/nobackup/sleblan2/AAC_DARF/surface_albedo/'
fp_out = '/nobackup/sleblan2/AAC_DARF/input/' + std_label + '/'
fp_pmom = '/nobackup/sleblan2/AAC_DARF/rtm/'
fp_uvspec = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
wvl_thm = '/nobackup/sleblan2/AAC_DARF/rtm/wvl.thm.dat'
vv = 'v3_20170616_CALIOP_AAC_AODxfAAC_With_Without_Sa'


# ## Check if reading or writing and input arguments

# In[6]:


import argparse
long_description = """    Prepare or read the AAC DARE calculations, with the spread for all values.
      - defaults to prepare AAC files """
parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-r','--doread',help='if set, will only read the output, not produce them',
                    action='store_true',default=False)
parser.add_argument('-w','--dowrite',help='if set, will write the input and list files for fuliou',
                    action='store_true',default=False)
parser.add_argument('-s','--shortread',help='if set, will only read the available files outputted in the fp_output file pat',
                    action='store_true',default=False)
parser.add_argument('-i','--index',help='index (from 0 to 36) to split up the work on each node/ COD',type =int)
#parser.add_argument('-i','--index',help='Sets the index of the pixel file to use (19374,22135). Default is 0',type=int)
#parser.add_argument('-t','--tmp_folder',help='Set to use the temporary folder',action='store_true',default=False)
parser.add_argument('-t','--tmp_folder',nargs='?',help='The path to the temp directory to use')
parser.add_argument('-v','--verbose',help='If set, outputs comments about the progress',action='store_true',default=False)
parser.add_argument('-j','--splitindex',
                    help='add extra level of index splitting, on addition to -i, split up the work on each node /COD and per 10/ext',action='store_true',default=False)
parser.add_argument('-m','--match',help='if set, matches the index value to all other ext, ssa, asy, for use in debugging',action='store_true',default=False)
in_ = vars(parser.parse_args())
doread = in_.get('doread',False)
dowrite = in_.get('dowrite',False)
i = in_.get('index',0)
shortread = in_.get('shortread',False)
verbose = in_.get('verbose',False)
j_index = in_.get('splitindex',False)
match = in_.get('match',False)


# In[ ]:


if in_.get('tmp_folder'):
    tmp_folder = in_.get('tmp_folder')
    tmp = True
else:
    tmp = False


# In[22]:


if i>0:
    if j_index:
        std_label = std_label+'_{:02.0f}_e{:02.0f}'.format(i%imax,i/imax)
    else:
        std_label = std_label+'_{:02.0f}'.format(i)


# In[23]:


if tmp:
    fp_out = tmp_folder+'/AAC_DARF/input/' + std_label + '/'


# In[24]:


if verbose: 
    print 'Running on folders: {}'.format(fp_out)


# # Load the input files and create the subsets

# In[22]:


mmm = 'JJA'


# In[31]:


fpm = fp+'Input_to_DARF_{mmm}_{vv}.mat'.format(mmm=mmm,vv=vv)
if verbose:
    print 'in %s months, getting mat file: %s' % (mmm,fpm)


# In[32]:


input_mmm = sio.loadmat(fpm,mat_dtype=True)['data_input_darf']


# In[70]:


geo = {'zout':[0,3,100],'year':2007,'day':15,'minute':0,'second':0}
geo['month'] = 7
doy = 225
aero = {'z_arr':[3.0,4.0]}
cloud = {'ztop':3.0,'zbot':2.0,'phase':'wc','write_moments_file':True}
source = {'integrate_values':True,'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/','run_fuliou':True}
albedo = {'create_albedo_file':False}
albedo['sea_surface_albedo'] = True


# In[34]:


geo['lat'],geo['lon'] = -7.0,-10.0


# In[63]:


ia1,ia2,io1,io2 = 16,25,35,39
ilat = np.arange(18,26)
ilon = np.arange(37,40)


# In[203]:


cod = input_mmm['MODIS_COD_mean'][0,0][ia1:ia2,io1:io2].flatten(order='F')
ext = np.append(np.abs(input_mmm['MOC_ext_mean'][0,0][ia1:ia2,io1:io2,:]).reshape([36,30],order='F'),np.zeros([1,30]),axis=0)
ssa = input_mmm['MOC_ssa_mean'][0,0][ia1:ia2,io1:io2,:].reshape([36,30],order='F')
asy = input_mmm['MOC_asym_mean'][0,0][ia1:ia2,io1:io2,:].reshape([36,30],order='F')


# In[84]:


cloud['ref'] = np.nanmean(input_mmm['MODIS_effrad_mean'][0,0][ia1:ia2,io1:io2])


# # Prepare the inputs

# In[ ]:


pmom_solar = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='solar',deltascale=False)
pmom_thermal = RL.make_pmom_inputs(fp_rtm=fp_pmom,source='thermal')
max_nmom=20
pmom_solar['max_nmoms'] = max_nmom
pmom_thermal['max_nmoms'] = max_nmom


# # Start the writing out

# In[ ]:


change_fp_output = True


# In[ ]:


fp_out2 = fp_out+mmm+std_label+'/'
if not os.path.exists(fp_out2):
    os.makedirs(fp_out2)
if change_fp_output:
    fp_output = fp_out2.replace('input','output')
    if not os.path.exists(fp_output):
        os.makedirs(fp_output)
fp_base_file = fp_out2+'base.inp'
make_base = True


# In[ ]:


if dowrite: 
    ff = 'AAC_list_file_{m}_{v}{lbl}.sh'.format(m=mmm,v=vv,lbl=std_label)
    file_list = file(fp_pmom+ff,'w')
    if verbose: 
        print 'Starting list file: '
    print fp_pmom+ff


# In[217]:


aero['wvl_arr'] = input_mmm['MOC_wavelengths'][0,0][0,:]*1000.0


# In[ ]:


if i>0:
    if j_index:
        codd = [cod[i%imax]]
        iee = np.arange((i/imax)*4,(i/imax+1)*4)
        if (i/imax+1)*4 > 33:
             iee = np.arange((i/imax)*4,37)
    else:
        codd = [cod[i]]
        iee = np.arange(0,37)
else:
    codd = cod
    iee = np.arange(0,37)


# In[ ]:


if match:
    iee = [i]

if verbose:
    print 'Running through the permutations'
b = []
first_time = True
for icod,c in enumerate(codd):
    if i>0: icod = i
    dd = {}
    # set the cloud values
    cloud['tau'] = c
    try: cloud['tau'][cloud['tau']<0.0] = 0.0
    except: pass
    try: cloud['ref'][cloud['ref']<2.0] = 2.0
    except: pass
    
    cloud_file_name_sol = fp_out2+'AAC_input_cod%02i_%s_sol.inp_cloud' % (icod,mmm)
    cloud_file_name_thm = fp_out2+'AAC_input_cod%02i_%s_thm.inp_cloud' % (icod,mmm)
    
    dd['cld_f_sol'] = copy(cloud_file_name_sol)
    dd['cld_f_thm'] = copy(cloud_file_name_thm)

    if dowrite:
        RL.write_cloud_file_moments(cloud_file_name_sol,cloud['tau'],cloud['ref'],cloud['zbot'],cloud['ztop'],
                                    verbose=verbose,moms_dict=pmom_solar,wvl_range=[250,5600])
        RL.write_cloud_file_moments(cloud_file_name_thm,cloud['tau'],cloud['ref'],cloud['zbot'],cloud['ztop'],
                                    verbose=verbose,moms_dict=pmom_thermal,wvl_range=[4000,50000-1])
    dd['cod'],dd['ref'],dd['zbot'],dd['ztop'] = copy(cloud['tau']),copy(cloud['ref']),copy(cloud['zbot']),copy(cloud['ztop'])
    
    if aero['wvl_arr'].max()<100000.0:
        aero['wvl_arr'] = np.append(aero['wvl_arr'],100000.0)
    
    for iext,e in enumerate(ext):
        if not iext in iee:
            continue
        # set the aerosol values
        aero['ext'] = e
        aero['ext'][aero['ext']<0.0] = 0.0
        if np.isnan(aero['ext']).all():
            if verbose:
                print 'skipping cod:%i, ext:%i' % (icod,iext)
            continue
        if not len(aero['ext'])==len(aero['wvl_arr']):
            aero['ext'] = np.append(aero['ext'],e[-1])
        for issa,s in enumerate(ssa):
            if match:
                if not issa in iee:
                    continue
            aero['ssa'] = s
            try: aero['ssa'][aero['ssa']<0.0] = 0.0
            except: pass
            try: aero['ssa'][aero['ssa']>1.0] = 1.0
            except: pass
            if not len(aero['ssa'])==len(aero['wvl_arr']):
                aero['ssa'] = np.append(aero['ssa'],s[-1])
            for iasy,a in enumerate(asy):
                if match:
                    if not iasy in iee:
                        continue
                aero['asy'] = a

                #sanitize inputs after adding subtracting standard deviations
                try: aero['asy'][aero['asy']<0.0] = 0.0
                except: pass
                try: aero['asy'][aero['asy']>1.0] = 1.0
                except: pass
                if not len(aero['asy'])==len(aero['wvl_arr']):
                    aero['asy'] = np.append(aero['asy'],a[-1])
                
                aero['file_name'] = fp_out2+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_sol.inp_aero' % (icod,iext,issa,iasy,mmm)
                dd['aero'] = deepcopy(aero)
                
                fsol = []
                fthm = []
                
                fsol_o = []
                fthm_o = []
                
                for HH in xrange(24):
                    geo['hour'] = HH
                    #build the solar input file
                    source['source'] = 'solar'
                    source['wvl_range'] = [250,5600]
                    source['wvl_filename'] = None
                    file_out_sol = fp_out2+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_sol.inp' % (icod,iext,issa,iasy,mmm,HH)
                    fsol.append(file_out_sol)
                    fsol_o.append(fp_output+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_sol.out' % (icod,iext,issa,iasy,mmm,HH))
                    
                    #build the thermal input file
                    source['source'] = 'thermal'
                    source['wvl_range'] = [4000,50000-1]
                    source['wvl_filename'] = None
                    file_out_thm = fp_out2+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_thm.inp' % (icod,iext,issa,iasy,mmm,HH)
                    fthm.append(file_out_thm)
                    fthm_o.append(fp_output+'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_thm.out' % (icod,iext,issa,iasy,mmm,HH))
                    
                    if dowrite:
                        file_list.write(fp_uvspec+' < "'+file_out_sol+'" > "'+fp_output
                                    +'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_sol.out"\n' % (icod,iext,issa,iasy,mmm,HH))
                        file_list.write(fp_uvspec+' < "'+file_out_thm+'" > "'+fp_output
                                    +'AAC_input_cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i_thm.out"\n' % (icod,iext,issa,iasy,mmm,HH))
                        if first_time:
                            cloud['link_to_mom_file'],cloud['file_name'] = True,cloud_file_name_sol
                            RL.write_input_aac(file_out_sol,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                                               make_base=True,fp_base_file=fp_base_file,set_quiet=True,solver='rodents')
                            first_time = False
        

                #print mmm,icod,iext,issa,iasy
                dd['geo'],dd['source'],dd['albedo'],dd['fp_base_file'] = deepcopy(geo),deepcopy(source),deepcopy(albedo),fp_base_file
                dd['fsol'],dd['fthm'],dd['fsol_o'],dd['fthm_o'] = fsol[:],fthm[:],fsol_o[:],fthm_o[:]
                dd['icod'],dd['iext'],dd['issa'],dd['iasy'],dd['mmm'] = copy(icod),copy(iext),copy(issa),copy(iasy),copy(mmm)
                b.append(deepcopy(dd))
if dowrite: file_list.close()


# In[ ]:


if verbose:
    pbar = tqdm(total=len(b))


# # Functions for use in this program

# In[ ]:


def worker_init():
    # ignore the SIGINI in sub process, just print a log
    def sig_int(signal_num, frame):
        if verbose: 
            print 'signal: %s' % signal_num
        raise IOError
    signal.signal(signal.SIGINT, sig_int)


# In[ ]:


def print_input_aac_24h(d):
    'function to print out the 24hour files (sol and thm) from an input dict'
    cloud = {'tau':d['cod'],'ref':d['ref'],'zbot':d['zbot'],'ztop':d['ztop'],
             'link_to_mom_file':True,'phase':'wc','write_moments_file':True}
    geo = d['geo']
    source = d['source']
    albedo = d['albedo']
    fp_base_file = d['fp_base_file']
    aero = d['aero']
    aero['link_to_mom_file'] = False
    
    #print 'fname: {fsol[0]}, iext: {iext}, issa: {issa}, iasy: {iasy}'.format(**d)
    
    for HH in xrange(24):
        geo['hour'] = HH
        #build the solar input file
        source['source'] = 'solar'
        source['wvl_range'] = [250,5600]
        source['wvl_filename'] = None
        file_out_sol = d['fsol'][HH]
        cloud['file_name'] = d['cld_f_sol']
        RL.write_input_aac(file_out_sol,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                           make_base=False,fp_base_file=fp_base_file,set_quiet=True,solver='rodents')
        
        #build the thermal input file
        source['source'] = 'thermal'
        source['wvl_range'] = [4000,50000-1]
        source['wvl_filename'] = None
        file_out_thm = d['fthm'][HH]
        cloud['file_name'] = d['cld_f_thm']
        
        RL.write_input_aac(file_out_thm,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,verbose=False,
                           make_base=False,fp_base_file=fp_base_file,set_quiet=True,solver='rodents')

        if not aero['link_to_mom_file']:
            aero['link_to_mom_file'] = True
        if verbose:
            pbar.update(1)


# In[ ]:


def read_input_aac_24h(d):
    'function to read out the 24hour files (sol and thm) from an input dict'

    zout=[0,3,100]
    nz = len(zout)
    output = {'zout':zout,'ext':d['aero']['ext'],'asy':d['aero']['asy'],'ssa':d['aero']['ssa']}
    output['SW_irr_dn_utc'] = np.zeros((nz,24))
    output['SW_irr_up_utc'] = np.zeros((nz,24))
    output['LW_irr_dn_utc'] = np.zeros((nz,24))
    output['LW_irr_up_utc'] = np.zeros((nz,24))  
    
    for HH in xrange(24):
        geo['hour'] = HH
        file_out_sol = d['fsol_o'][HH]
        file_out_thm = d['fthm_o'][HH]
        try:
            sol = RL.read_libradtran(file_out_sol,zout=zout)
            thm = RL.read_libradtran(file_out_thm,zout=zout)
        except IOError:
            print 'File not found skip: cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i' %(d['icod'],d['iext'],d['issa'],d['iasy'],d['mmm'],HH)
            if HH==0:
                print file_out_sol
            continue
        except ValueError:
            print 'Problem with file: cod%02i_ext%02i_ssa%02i_asy%02i_%s_HH%02i' %(d['icod'],d['iext'],d['issa'],d['iasy'],d['mmm'],HH)
            output['SW_irr_dn_utc'][:,HH] = np.nan
            output['SW_irr_up_utc'][:,HH] = np.nan
            output['LW_irr_dn_utc'][:,HH] = np.nan
            output['LW_irr_up_utc'][:,HH] = np.nan
            continue
        output['SW_irr_dn_utc'][:,HH] = sol['direct_down']+sol['diffuse_down']
        output['SW_irr_up_utc'][:,HH] = sol['diffuse_up']

        output['LW_irr_dn_utc'][:,HH] = thm['direct_down']+thm['diffuse_down']
        output['LW_irr_up_utc'][:,HH] = thm['diffuse_up']

    output['SW_irr_dn_avg'] = np.mean(output['SW_irr_dn_utc'],axis=1)
    output['SW_irr_up_avg'] = np.mean(output['SW_irr_up_utc'],axis=1)
    output['LW_irr_dn_avg'] = np.mean(output['LW_irr_dn_utc'],axis=1)
    output['LW_irr_up_avg'] = np.mean(output['LW_irr_up_utc'],axis=1)
    if verbose:
        pbar.update(1)

    return output        


# # Core of the reading/writing routines, calling pool of workers

# In[ ]:


if dowrite:
    if verbose:
        print 'Now preparing the pool of workers and making them print the files'
    p = Pool(cpu_count(),worker_init)
    results = p.map(print_input_aac_24h,b)
    p.close()


# In[ ]:


if doread:
    il = len(b)
    if shortread:
        ll = os.listdir(fp_output)
        il = len(ll)
        if verbose:
            print 'Only reading over {} files'.format(il)
    
    if verbose:
        print 'Now preparing the pool of workers and making them read all the files'
    p = Pool(cpu_count(),worker_init)
    results = p.map(read_input_aac_24h,b[0:il])
    p.close()
    
    nm_list = {'SW_irr_dn_avg':'SW_irr_dn_avg','SW_irr_up_avg':'SW_irr_up_avg',
               'LW_irr_dn_avg':'LW_irr_dn_avg','LW_irr_up_avg':'LW_irr_up_avg',
               'ssa':'ssa','asy':'asy','ext':'ext'}
    saves = {}
    for a in nm_list.keys():
        saves[a] = np.array([results[n][nm_list[a]] for n in xrange(il)])

    saves['geo'] = geo
    saves['source'] = source
    saves['source']['wvl_filename'] = []
    saves['albedo'] = albedo
    saves['cloud'] = cloud
    
    fp_save = fp+'AAC_spread_{m}_{v}_{i}.mat'.format(m=mmm,v=vv,i=i)
    try:
        if j_index:
            fp_save = fp+'AAC_spread_{m}_{v}_{i}_e{j}.mat'.format(m=mmm,v=vv,i=i%imax,j=i/imax)
    except:
        pass
    
    if verbose:
        print 'Saving read file: '+fp_save
    sio.savemat(fp_save,saves)

