
# coding: utf-8

# # Info
# Name:  
# 
#     TCAP_cloud
# 
# Purpose:  
# 
#     To compare the various cloud properties retrieved via different methods from TCAP.
#     Looking at the Feb, 19th, 2013 case
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     figures and save files...
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - load_modis.py : for loading modis files
#     - matplotlib
#     - numpy
#     - scipy : for saving and reading
#     - math
#     - os
#     - gc
#     - pdb
#     - datetime
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - 20151117_zen_cld_retrieved.mat: cloud retrieval file
#   - MYD06_L2.A2015321.1540.006.2015322185040.hdf: MODIS file
#   
# Modification History:
# 
#     Written: Samuel LeBlanc, NASA Ames, Santa Cruz, CA, 2016-04-06
#     Modified: 

# # Import initial modules and default paths
# 

# In[2]:

get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.io as sio
import Sp_parameters as Sp
import hdf5storage as hs
import load_modis as lm
import os


# In[3]:

from mpl_toolkits.basemap import Basemap,cm


# In[4]:

import write_utils as wu
import plotting_utils as pu
import map_utils as mu
import Run_libradtran as RL


# In[5]:

get_ipython().magic(u'matplotlib notebook')


# In[6]:

# set the basic directory path
fp='C://Users//sleblan2//Research//TCAP//'


# # Load the various data

# ## Verify the retrievals from 4STAR

# ### Build the proper lut from the idl out file

# In[7]:

vv = 'v1'


# In[8]:

# load the idl save file containing the modeled radiances
s = sio.idl.readsav(fp+'model/sp_v1_20130219_4STAR.out')
print s.keys()
print 'sp', s.sp.shape
print 'sp (wp,   wvl,  z,  re,  ta)'


# ### Compare to lut from NAAMES

# In[9]:

fp_lut_mat = 'C:\\Users\\sleblan2\\Research\\NAAMES\\lut\\v2_NAAMES_lut.mat'
print('Loading the lut file:{}'.format(fp_lut_mat))
if not os.path.isfile(fp_lut_mat):
    print('File {} does not exist'.format(fp_lut_mat))
    raise IOError('LUT File not found: {}'.format(fp_lut_mat))
luts = hs.loadmat(fp_lut_mat)


# In[22]:

luts.keys()


# In[23]:

luts['irr_up'].shape


# In[27]:

luts['rad'].shape, s.sp.shape


# In[28]:

s.sp[:,:,:,:,:,np.newaxis].shape


# In[38]:

s.ref


# In[39]:

luts['ref']


# ### reshape TCAP LUT to fit NAAMES format and save new file

# In[54]:

tcap_lut = {u'tau':s.tau,
            u'rad':s.sp[:,:,:,:,:,np.newaxis],
            u'sza':np.array(s.sza)[np.newaxis],
            u'irr_dn_diff':s.sp_irrdn[:,:,:,:,:,np.newaxis]*0.0,
            u'irr_dn':s.sp_irrdn[:,:,:,:,:,np.newaxis],
            u'irr_up':s.sp_irrup[:,:,:,:,:,np.newaxis],
            u'zout':s.zout,
            u'wvl':s.zenlambda,
            u'phase':['wc','ic'],
            u'ref':s.ref}


# In[55]:

hs.savemat(fp+'model\\{}_TCAP_lut.mat'.format(vv),tcap_lut)


# ### Run the retrieval via the run_zen_cld_retrieval command line
# in the python_codes directory:
# 
# > python Run_zen_cld_retrieval.py -lut C:\Users\sleblan2\Research\TCAP\model\v1_TCAP_lut.mat C:\Users\sleblan2\Research\TCAP\20130219starzen.mat -o C:\Users\sleblan2\Research\TCAP\plots\ 

# ## Check the 4STAR data

# In[8]:

ss = sio.loadmat(fp+'20130219starzen.mat')


# In[13]:

ss.keys()


# In[14]:

ss['t'].shape


# In[15]:

ss['rad'].shape


# In[16]:

plt.figure()
plt.plot(ss['t'],ss['rad'][:,400])


# ## Load the retrieved cloud properties from 4STAR

# In[9]:

star = hs.loadmat(fp+'4STAR//20130219_zen_cld_retrieved.mat')


# In[10]:

star.keys()


# In[11]:

plt.figure()
plt.plot(star['tau'])


# In[11]:

plt.figure()
plt.plot(star['utc'],star['tau'])


# In[13]:

star['ref'] = Sp.smooth(star['ref'],20)


# In[14]:


plt.figure()
plt.plot(star['utc'],star['ref'],'x-')



# ## Load the important MODIS files

# In[12]:

myd06_file = fp+'MODIS\\MYD06_L2.A2013050.1725.006.2014260074007.hdf'
myd03_file = fp+'MODIS\\MYD03.A2013050.1725.006.2013051163424.hdf'
print os.path.isfile(myd03_file) #check if it exists
print os.path.isfile(myd06_file)


# In[13]:

import load_modis as lm
modis,modis_dicts = lm.load_modis(myd03_file,myd06_file)


# ### Subset the MODIS values to match the flight path

# In[14]:

import map_utils as mu


# In[15]:

mod_ind = mu.map_ind(modis['lon'],modis['lat'],star['lon'],star['lat'])


# ### Now compare MODIS vs. 4STAR with bean plots
# 

# In[16]:

import plotting_utils as pu


# In[55]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,60],xlim=[-0.5,1.5])
ax1.set_ylabel('$\\tau$')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['MODIS\n(reflected)','4STAR\n(transmitted)'])
pu.plot_vert_hist(fig,ax1,modis['tau'][mod_ind[0,:],mod_ind[1,:]],0,[0,60],legend=True,onlyhist=False,loc=2,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star['tau'],1,[0,60],legend=True,color='r',bins=50)
plt.savefig(fp+'plots/20130219_COD_bean_modis_4star.png',transparent=True,dpi=600)


# In[56]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[0,60],xlim=[-0.5,1.5])
ax1.set_ylabel('r$_{eff}$ [$\\mu$m]')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['MODIS\n(reflected)','4STAR\n(transmitted)'])
pu.plot_vert_hist(fig,ax1,modis['ref'][mod_ind[0,:],mod_ind[1,:]],0,[0,60],legend=True,onlyhist=False,loc=2,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star['ref'],1,[0,60],legend=True,color='r',bins=50)
plt.savefig(fp+'plots/20130219_ref_bean_modis_4star.png',transparent=True,dpi=600)


# # Save to ict files

# ## Save the ict file for 4STAR cloud data

# In[17]:

d_dict =  {'utc':{'data':star['utc'],'unit':'hours from midnight UTC',
                     'long_description':'Fractional hours starting a midnight, continuous'},
          'COD':{'data':star['tau'],'unit':'None','format':'.1f',
                 'long_description':'Cloud Optical Depth of the cloud above the Aircraft'},
          'REF':{'data':star['ref'],'unit':'microns','format':'.1f',
                 'long_description':'Cloud particle effective radius, pertains to liquid cloud drops and ice crystals'},
          'PHASE':{'data':star['phase'],'unit':'None','format':'.0f',
                   'long_description':'Thermodynamic phase of cloud above, 0: pure liquid cloud, 1: pure ice cloud, mixed phase not retrieved'},
          'LAT':{'data':star['lat'],'unit':'Degrees','format':'.6f',
                 'long_description':'Aircraft position latitude (North positive)'},
          'LON':{'data':star['lon'],'unit':'Degrees','format':'.6f',
                 'long_description':'Aircraft position longitude (East positive)'},
          'ALT':{'data':star['alt'],'unit':'meter','format':'.1f',
                 'long_description':'Aircraft altitude'},
          'SZA':{'data':star['sza'],'unit':'Degrees','format':'.2f',
                 'long_description':'Solar Zenith Angle, angle of the sun between it and zenith'}
          }


# In[18]:

h_dict ={'PI':'Jens Redemann',
         'Institution':'NASA Ames Research Center',
         'Instrument':'4STAR',
         'campaign':'TCAP #1',
         'special_comments':'Preliminary retrieved cloud properties data',
         'PI_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
         'platform':'C130',
         'location':"centered at Cape Cod, Massachussetts, actual location of DOE G-1 described by lat and lon below",
         'instrument_info':'Retrieved products from the 4STAR zenith radiance measurements',
         'data_info':'For references see LeBlanc et al.(2015) AMT, doi:10.5194/amt-8-1361-2015',
         'time_interval':1.0,
         'uncertainty':'Preliminary 7% in REF and 5% in COD',
         'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
         'project_info':'TCAP field mission',
         'stipulations':'Prior OK from PI',
         'rev_comments':"""RA: preliminary retrieved values, may be subject to multiple errors
    including due to clouds influencing presumed surface albedo, non-homogeneous clouds, or mixed phase clouds"""
        }


# In[19]:

order=['LAT','LON','ALT','SZA','COD','REF','PHASE']


# In[20]:

data_dict = wu.prep_data_for_ict(d_dict,time_interval=1.0)


# In[21]:

data_dict.keys()


# In[47]:

wu.make_plots_ict(data_dict,filepath=fp+'plots/',data_id='4STAR-zen-cloud',loc_id='G1',date='20130219',rev='RA',
                  plot_together=['COD','REF','PHASE'],plot_together2=['LAT','LON','SZA'])


# In[22]:

wu.write_ict(h_dict,data_dict,filepath=fp,data_id='4STAR-zen-cloud',loc_id='G1',date='20130219',rev='RA',order=order)


# ## Save the ICT file for the MODIS cloud properties

# In[23]:

md_dict =  {'utc':{'data':star['utc'],'unit':'hours from midnight UTC',
                     'long_description':'Fractional hours starting a midnight, continuous'},
          'COD':{'data':modis['tau'][mod_ind[0,:],mod_ind[1,:]],'unit':'None','format':'.1f',
                 'long_description':'Cloud Optical Depth from MODIS'},
          'REF':{'data':modis['ref'][mod_ind[0,:],mod_ind[1,:]],'unit':'microns','format':'.1f',
                 'long_description':'Cloud particle effective radius, pertains to liquid cloud drops and ice crystals'},
          'PHASE':{'data':modis['phase'][mod_ind[0,:],mod_ind[1,:]],'unit':'None','format':'.0f',
                   'long_description':'Thermodynamic phase of cloud,'+\
                   ' 0 -- cloud free, 1 -- water cloud, 2 -- ice cloud, 3 -- mixed phase cloud, 6 -- undetermined phase'},
          'LAT':{'data':modis['lat'][mod_ind[0,:],mod_ind[1,:]],'unit':'Degrees','format':'.6f',
                 'long_description':'MODIS linked to Aircraft position latitude (North positive)'},
          'LON':{'data':modis['lon'][mod_ind[0,:],mod_ind[1,:]],'unit':'Degrees','format':'.6f',
                 'long_description':'MODIS linked to Aircraft position longitude (East positive)'},
          'SZA':{'data':star['sza'],'unit':'Degrees','format':'.2f',
                 'long_description':'Solar Zenith Angle, angle of the sun between it and zenith'}
          }


# In[24]:

mh_dict ={'PI':'Samuel LeBlanc',
         'Institution':'NASA Ames Research Center',
         'Instrument':'MODIS',
         'campaign':'TCAP #1',
         'special_comments':'MODIS retrieved cloud values linked along the DOE G1 flight path',
         'PI_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
         'platform':'G1',
         'location':"Centered at Cape Cod, Massachussetts, actual location of DOE G1 described by lat and lon below",
         'instrument_info':'Retrieved products from the MODIS, MYD06_L2.A2015321.1540.006.2015322185040.hdf',
         'data_info':'For references see LeBlanc et al.(2015) AMT, doi:10.5194/amt-8-1361-2015',
         'time_interval':1.0,
         'uncertainty':'N\A',
         'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
         'project_info':'TCAP field mission',
         'stipulations':'Prior OK from PI',
         'rev_comments':"""RA: initial go at this, for radiative transfer calculations"""
        }


# In[25]:

order=['LAT','LON','SZA','COD','REF','PHASE']


# In[26]:

mdata_dict = wu.prep_data_for_ict(md_dict,time_interval=1.0)


# In[30]:

wu.write_ict(mh_dict,mdata_dict,filepath=fp,data_id='MODIS-cloud-to-G1',loc_id='G1',date='20130219',rev='RA',order=order)


# # Prepare input files for radiative transfer
# 

# In[27]:

import Run_libradtran as Rl


# ## Prepare the defaults

# In[32]:

from datetime import datetime
datetime(2013,02,19).timetuple().tm_yday


# In[55]:

geo = {'lat':47.6212167,'lon':52.74245,'doy':50,'zout':[0,100.0]}
aero = {} # none
cloud = {'ztop':10.5,'zbot':10.0,'write_moments_file':False}
source = {'wvl_range':[201.0,4000.0],'source':'solar','integrate_values':True,'run_fuliou':True,
          'dat_path':'/u/sleblan2/libradtran/libRadtran-2.0-beta/data/'}
albedo = {'create_albedo_file':False,'sea_surface_albedo':True,'wind_speed':10.0}


# In[56]:

cloud['phase'] = 'wc'
geo['sza'] = 40.0
cloud['tau'] = 2.0
cloud['ref'] = 5.0


# In[57]:

phase_star = {0:'wc',1:'ic'}


# In[58]:

phase_modis = {0:'wc',1:'wc',2:'ic',3:'ic',6:'wc'}


# ## Create input files for 4STAR measurements

# In[59]:

# open the list file
f = open(fp+'rtm/TCAP_20130219_CRE.sh','w')
fpp_in = '/nobackup/sleblan2/rtm/input/TCAP_CRE_20130219/'
fpp_out = '/nobackup/sleblan2/rtm/output/TCAP_CRE_20130219/'
fp_uv = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
fp_in = fp+'rtm/input/CRE/'


# In[60]:

for i,l in enumerate(data_dict['LAT']['data']):
    if l<-100.: # for only valid values
        continue
    if not np.isfinite(data_dict['PHASE']['data'][i]) or not np.isfinite(data_dict['COD']['data'][i]):
        continue
    print i
    
    f_in = 'TCAP_v1_star_{:03d}.dat'.format(i)
    geo['lat'],geo['lon'],geo['sza'] = l,data_dict['LON']['data'][i],data_dict['SZA']['data'][i]
    cloud['tau'],cloud['ref'] = data_dict['COD']['data'][i],data_dict['REF']['data'][i]
    cloud['phase'] = phase_star[data_dict['PHASE']['data'][i]]
    Rl.write_input_aac(fp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
    f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
    
    f_in = 'TCAP_v1_star_{:03d}_clear.dat'.format(i)
    cloud['tau'] = 0.0
    Rl.write_input_aac(fp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
    f.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
f.close()


# ## Create input files for MODIS clouds

# In[61]:

# open the list file
fm = open(fp+'rtm/TCAP_20130219_CRE_modis.sh','w')
fpp_in = '/nobackup/sleblan2/rtm/input/TCAP_CRE_20130219/'
fpp_out = '/nobackup/sleblan2/rtm/output/TCAP_CRE_20130219/'
fp_uv = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/uvspec'
fp_in = fp+'rtm/input/CRE/'


# In[62]:

for i,l in enumerate(mdata_dict['LAT']['data']):
    if l<-100.: # for only valid values
        continue
    if not np.isfinite(mdata_dict['PHASE']['data'][i]) or not np.isfinite(mdata_dict['COD']['data'][i]):
        continue
    print i
    
    f_in = 'TCAP_v1_modis_{:03d}.dat'.format(i)
    geo['lat'],geo['lon'],geo['sza'] = l,mdata_dict['LON']['data'][i],mdata_dict['SZA']['data'][i]
    cloud['tau'],cloud['ref'] = mdata_dict['COD']['data'][i],mdata_dict['REF']['data'][i]
    cloud['phase'] = phase_modis[mdata_dict['PHASE']['data'][i]]
    Rl.write_input_aac(fp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
    fm.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
    
    f_in = 'TCAP_v1_modis_{:03d}_clear.dat'.format(i)
    cloud['tau'] = 0.0
    Rl.write_input_aac(fp_in+f_in,geo=geo,aero=aero,cloud=cloud,source=source,albedo=albedo,
                               verbose=False,make_base=False,set_quiet=True)
    fm.write('{uv} < {fin} > {out}\n'.format(uv=fp_uv,fin=fpp_in+f_in,out=fpp_out+f_in))
fm.close()


# # Run the libradtran runs on the pleiades

# Go to folder on this computer and run the rsync: (make sure the fodler is created on bridge)
# > rsync -auzvb ./ bridge:/nobackup/sleblan2/rtm/input/TCAP_CRE_20130219/
# 
# Then logon to pleaides and run this command in the correct folder
# 
# > sed -i -- 's$C://Users//sleblan2//Research//TCAP//rtm/input/CRE/$/nobackup/sleblan2/rtm/input/TCAP_CRE_20130219$g' *.dat
# 
# Make sure the list files are put in as always
# 
# > scp *.sh bridge:~/TCAP/runs/
# 
# Once copied, run the files
# 
# > ~/TCAP/runs> ./TCAP_20130219_CRE.sh
# > ~/TCAP/runs> ./TCAP_20130219_CRE_modis.sh
# 
# After copy the results to your homw computer

# # Read the output from libradtran

# ## Read the out files

# In[28]:

nstar = len(data_dict['LAT']['data'])
nmodis = len(mdata_dict['LAT']['data'])
star_CRE = {'dn':np.zeros((nstar,2))+np.nan,'up':np.zeros((nstar,2))+np.nan}
star_CRE_clear = {'dn':np.zeros((nstar,2))+np.nan,'up':np.zeros((nstar,2))+np.nan}
modis_CRE = {'dn':np.zeros((nmodis,2))+np.nan,'up':np.zeros((nmodis,2))+np.nan}
modis_CRE_clear = {'dn':np.zeros((nmodis,2))+np.nan,'up':np.zeros((nmodis,2))+np.nan}
star_C = np.zeros((nstar,2))+np.nan
modis_C = np.zeros((nmodis,2))+np.nan


# In[31]:

print 'MODIS'
for i,l in enumerate(mdata_dict['LAT']['data']):
    if l<-100.: # for only valid values
        continue
    if not np.isfinite(mdata_dict['PHASE']['data'][i]) or not np.isfinite(mdata_dict['COD']['data'][i]):
        continue
    print '\r{}..'.format(i)
    try:
        f_in = 'TCAP_v1_modis_{:03d}.dat'.format(i)
        s = Rl.read_libradtran(fp+'rtm/output/CRE/'+f_in,zout=[0,100])
        f_in = 'TCAP_v1_modis_{:03d}_clear.dat'.format(i)
        sc = Rl.read_libradtran(fp+'rtm/output/CRE/'+f_in,zout=[0,100])
    except ValueError:
        print 'Problem with file at: {}'.format(i)
        continue
    
    modis_CRE['dn'][i,:] = s['diffuse_down']+s['direct_down']
    modis_CRE_clear['dn'][i,:] = sc['diffuse_down']+sc['direct_down']
    modis_CRE['up'][i,:] = s['diffuse_up']
    modis_CRE_clear['up'][i,:] = sc['diffuse_up']
    modis_C[i,:] = (modis_CRE['dn'][i,:]-modis_CRE['up'][i,:]) - (modis_CRE_clear['dn'][i,:]-modis_CRE_clear['up'][i,:])


# In[38]:

print '4STAR'
for i,l in enumerate(data_dict['LAT']['data']):
    if l<-100.: # for only valid values
        continue
    if not np.isfinite(data_dict['PHASE']['data'][i]) or not np.isfinite(data_dict['COD']['data'][i]):
        continue
    print '\r{}..'.format(i)
    try:
        f_in = 'TCAP_v1_star_{:03d}.dat'.format(i)
        s = Rl.read_libradtran(fp+'rtm/output/CRE/'+f_in,zout=[0,100])
        f_in = 'TCAP_v1_star_{:03d}_clear.dat'.format(i)
        sc = Rl.read_libradtran(fp+'rtm/output/CRE/'+f_in,zout=[0,100])
    except ValueError:
        print 'Problem with file at: {}'.format(i)
        continue
    except IOError:
        continue
        
    star_CRE['dn'][i,:] = s['diffuse_down']+s['direct_down']
    star_CRE_clear['dn'][i,:] = sc['diffuse_down']+sc['direct_down']
    star_CRE['up'][i,:] = s['diffuse_up']
    star_CRE_clear['up'][i,:] = sc['diffuse_up']
    star_C[i,:] = (star_CRE['dn'][i,:]-star_CRE['up'][i,:]) - (star_CRE_clear['dn'][i,:]-star_CRE_clear['up'][i,:])


# ## Double check the output

# In[39]:

star_CRE


# In[41]:

star_CRE['dn'][34]


# In[36]:

modis_CRE


# In[43]:

modis_C[34]


# ## Plot the CRE

# In[68]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-450,0],xlim=[-0.5,1.5])
ax1.set_ylabel('Cloud Radiative Effect [W/m$^2$]')
ax1.set_title('Cloud Radiative Effect - Surface')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['MODIS\n(reflected)','4STAR\n(transmitted)'])
pu.plot_vert_hist(fig,ax1,modis_C[:,0],0,[-450,0],legend=True,onlyhist=False,loc=2,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star_C[:,0],1,[-450,0],legend=True,color='r',bins=50)
plt.savefig(fp+'plots/20130219_surface_CRE_modis_4STAR.png',transparent=True,dpi=600)


# In[66]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-450,0],xlim=[-0.5,1.5])
ax1.set_ylabel('Cloud Radiative Effect [W/m$^2$]')
ax1.set_title('Cloud Radiative Effect - TOA')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['MODIS\n(reflected)','4STAR\n(transmitted)'])
pu.plot_vert_hist(fig,ax1,modis_C[:,1],0,[-450,0],legend=True,onlyhist=False,loc=2,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star_C[:,1],1,[-450,0],legend=True,color='r',bins=50)
plt.savefig(fp+'plots/20130219_toa_CRE_modis_4STAR.png',transparent=True,dpi=600)


# ## Calculate the relative CRE (rCRE)

# In[77]:

modis_rC = np.zeros_like(modis_C)
star_rC = np.zeros_like(star_C)


# In[79]:

modis_rC[:,0] = modis_C[:,0]/modis_CRE['dn'][:,1]*100.0
star_rC[:,0] = star_C[:,0]/star_CRE['dn'][:,1]*100.0
modis_rC[:,1] = modis_C[:,1]/modis_CRE['dn'][:,1]*100.0
star_rC[:,1] = star_C[:,1]/star_CRE['dn'][:,1]*100.0


# In[80]:

star_rC.shape


# In[85]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-100,0],xlim=[-0.5,1.5])
ax1.set_ylabel('relative Cloud Radiative Effect [%]')
ax1.set_title('relative Cloud Radiative Effect - Surface')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['MODIS\n(reflected)','4STAR\n(transmitted)'])
pu.plot_vert_hist(fig,ax1,modis_rC[:,0],0,[-100,0],legend=True,onlyhist=False,loc=4,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star_rC[:,0],1,[-100,0],legend=True,color='r',bins=50)
plt.savefig(fp+'plots/20130219_surface_rCRE_modis_4STAR.png',transparent=True,dpi=600)


# In[87]:

np.nanmean(modis_rC[:,0]),np.nanmean(star_rC[:,0])


# In[88]:

np.nanstd(modis_rC[:,0]), np.nanstd(star_rC[:,0])


# In[86]:

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8],ylim=[-100,0],xlim=[-0.5,1.5])
ax1.set_ylabel('relative Cloud Radiative Effect [%]')
ax1.set_title('relative Cloud Radiative Effect - TOA')
ax1.set_xticks([0,1])
ax1.set_xticklabels(['MODIS\n(reflected)','4STAR\n(transmitted)'])
pu.plot_vert_hist(fig,ax1,modis_rC[:,1],0,[-100,0],legend=True,onlyhist=False,loc=4,color='g',bins=50)
pu.plot_vert_hist(fig,ax1,star_rC[:,1],1,[-100,0],legend=True,color='r',bins=50)
plt.savefig(fp+'plots/20130219_TOA_rCRE_modis_4STAR.png',transparent=True,dpi=600)


# ## Testing area

# In[51]:

0.5*0.18065/0.001*3/2/7.8


# In[54]:

mdata_dict['COD']['data'][5]


# In[59]:

mdn_cl = 5.950567e+02+5.136170e+01
mup_cl = 4.183247e+01


# In[60]:

mdn =  4.806639e+01+4.617782e+02
mup = 3.633226e+01


# In[61]:

dn_cl = 5.950567e+02+5.136170e+01
up_cl = 4.183247e+01
dn = 1.840833e-02+3.051980e+02
up = 2.012267e+01


# In[62]:

mnet = mdn-mup
mnet_cl = mdn_cl-mup_cl


# In[63]:

mnet,mnet_cl,mnet-mnet_cl


# In[64]:

net = dn-up
net_cl = dn_cl-up_cl
net,net_cl,net-net_cl


# In[ ]:



