#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     To expand the aerosol optical depth for ORACLES skyscans, from Logan's work to Kristina's wavelenght span. 
#     Focus on 2016 skyscans.
#     0.4, 0.47, 0.55, 0.67, 0.86, 1.24, 2.1µm
#     
#     - Using the retrieved size distributions
#     - retrieval results of refractive index (imaginary and real) at wavelengths: 500, 675, 870, 995 nm
# 
# 
# Output:
# 
#     Figure and save files
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - load_utils.py
#     - matplotlib
#     - numpy
#     - Sp_parameters
#     - write_utils
#     - path_utils
#     - libradtran
# 
# Needed Files:
#   - netcdf of aeroinv from ORACLES
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-11-18,  based on TASNPP_Mie_calc_for_aerosol_above_cloud
#              
# 

# # Prepare python environment

# In[1]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'notebook')
import os
import scipy.interpolate as si
import netCDF4 as nc


# In[2]:


name = 'sunsat_ORACLES2016'
vv = 'v3'
fp = getpath(name)


# In[3]:


fp_bin = getpath('uvspec_bin')


# In[4]:


fp_rtm0 = getpath('rtm')
fp_rtm = fp_rtm0 +'ORACLES_mie/'


# In[5]:


if not os.path.exists(fp_rtm): 
    os.mkdir(fp_rtm)


# # Load files

# In[6]:


f = fp + 'data_processed/starskies/skyscans_Logan_20221116/4STAR-aeroinv4wvl_P3_2016_R0.nc'


# In[11]:


ae,ae_dict = lu.load_netcdf(f,everything=True)


# In[12]:


len(ae[b'time'])


# In[13]:


ae_dict[b'scan_tag']


# In[14]:


ae[b'scan_tag']


# ## Extra meta information

# In[15]:


sk_meta = sio.loadmat(fp+'data_processed/starskies/skyscans_ORACLES2016_moredata.mat')


# In[16]:


sk_meta.keys()


# In[17]:


sk_meta['skyscans_meta'].shape


# In[18]:


sk_meta['skyscans_qcflags'].shape


# In[19]:


sk_meta['skyscan_vars'][0]


# In[20]:


sk_meta['skyscans_meta'][0,:]


# In[21]:


sk_meta['skyscans_all'].shape


# In[191]:


sk_meta['skyscans_all'][0,24]


# In[64]:


sk_meta['days'] = [dd[0][0] for dd in sk_meta['skyscans_all'][:,0]]


# In[66]:


sk_meta['filenum'] = [dd[0][0] for dd in sk_meta['skyscans_all'][:,1]]


# ## Map the extra meta information to the skyscans in the file

# In[69]:


ae[b'scan_tag']


# In[87]:


from datetime import datetime


# In[118]:


days = []
for t in ae[b'time']:
    x =  datetime.fromtimestamp(datetime(2016,8,31).timestamp()+t)
    days.append(x.year*10000+x.month*100+x.day)


# In[120]:


days


# In[161]:


qc_flags = [4.0]*len(days)
for i,da in list(enumerate(days)):
    ig = [(da==dk)&(int(ae[b'scan_tag'][i])==sk_meta['filenum'][j]) for j,dk in enumerate(sk_meta['days'])] 
    if any(ig): qc_flags[i] = (sk_meta['skyscans_qcflags'][ig][0][0])


# In[163]:


qc_flags


# # Run through each retrieval and make input files

# ## Prep functions for printing size distribution

# From the libradtran documentation for mie size distribution file:  
# >Specify a two column file, r [micron], dn(r)/dr, which describes a size distribution

# In[22]:


def print_size_dist_file(fname,r,dnr):
    with open(fname,'w') as f:
        for ir,rr in list(enumerate(r)):
            f.write('{:3.10f} {:3.10f}\n'.format(rr,dnr[ir]))


# In[23]:


ae_dict[b'radius']


# In[24]:


ae_dict[b'psd']


# In[25]:


def convert_dvlnr_to_dndr(psd,r):
     # All about that conversion from the volume size distribution of dV(r)/dln(r) to number size distribution dN(r)/dr
    Nr = psd/(4.0*np.pi/3.0)
    for i,rr in list(enumerate(r)):
        Nr[i] = Nr[i]/rr**4.0
    return Nr


# In[26]:


print_size_dist_file(fp_rtm+'mie_tester.psd',ae[b'radius'],convert_dvlnr_to_dndr(ae[b'psd'][0,:],ae[b'radius']))


# ## Prep function for printing index of refraction

# Update from Kerry Meyer's email on Sept 21, 2022.  
# 
# > Thanks, Sam, this is helpful. Do you have the extrapolated phase functions also? And I hate to do this, but can the wavelength range be extended in both directions? I’d like to accommodate the VIIRS 2.25µm channel, and use the spectral response functions to compute band-weighted averages for each MODIS and VIIRS channel. So perhaps adding something like 400nm on the short end (which I guess is the exact 4STAR wavelength there) and 2300nm on the long end – that should encompass most of the response function ranges for MODIS and VIIRS.
# >
# > I’ll also need your input on how best to QA filter these scans. I see a “QA_level” dataset that I assume should be part of that.
# >
# >Kerry

# In[27]:


ae_dict[b'n_real']


# In[28]:


ae_dict[b'n_imag']


# In[29]:


ae_dict[b'wavelength']


# In[30]:


def spline_wavelength_extend(val,wavelen,new_wavelen,su=0.0006,k=2):
    # to calculate a new spline fit to the refractive index
    val_fx = si.CubicSpline(np.append(wavelen,wavelen[-1]+600.0),np.append(val,val[-1]),bc_type='natural',extrapolate=True)
    return val_fx(new_wavelen)


# In[32]:


if vv == 'v1':
    wave_out = np.array([470.0, 550.0, 670.0, 860.0, 1240.0, 2100.0])
elif vv == 'v2':
    wave_out = np.array([400.0, 470.0, 550.0, 670.0, 860.0, 1240.0, 2100.0, 2300.0])
elif vv == 'v3':
    wave_out = np.array([400.0, 470.0, 500.0, 550.0, 670.0, 860.0, 1240.0, 2100.0, 2300.0])


# In[33]:


fig,ax = plt.subplots(1,1)

for it,t in list(enumerate(ae[b'time'])):
    plt.plot(ae[b'wavelength'],ae[b'n_real'][it,:],':',c=plt.cm.viridis(it*3))
    plt.plot(wave_out,spline_wavelength_extend(ae[b'n_real'][it,:],ae[b'wavelength'],wave_out),
             '+-',c=plt.cm.viridis(it*3))
    
plt.title('Real refrac. index')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Refractive index')
    


# In[34]:


plt.figure()
for it,t in list(enumerate(ae[b'time'])):
    plt.plot(ae[b'wavelength'],ae[b'n_imag'][it,:],':',c=plt.cm.viridis(it*3))
    plt.plot(wave_out,spline_wavelength_extend(ae[b'n_imag'][it,:],ae[b'wavelength'],wave_out),
         '+-',c=plt.cm.viridis(it*3))
    
plt.title('Imag refrac. index')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Refractive index')
    


# In[35]:


def print_refrac_file(fname,wavelength,n_real,n_imag):
    #wavelength in nm, n_real and n_imag needs to be positive
    
    with open(fname,'w') as f:
        for iw,w in list(enumerate(wavelength)):
            f.write('{:4.3f} {:3.12f} {:3.12f}\n'.format(w,abs(n_real[iw]),abs(n_imag[iw])))


# In[36]:


itest = 0
print_refrac_file(fp_rtm+'mie_tester.ref',wave_out,
                     spline_wavelength_extend(ae[b'n_real'][itest,:],ae[b'wavelength'],wave_out),
                     spline_wavelength_extend(ae[b'n_imag'][itest,:],ae[b'wavelength'],wave_out))


# ## Mie input file function

# In[37]:


def mie_input(fname,refrac_file_name,size_dist_file_name,program='MIEV0',nmom=None):
    #simple mie input file program
     with open(fname,'w') as f:
            f.write('mie_program {}\n'.format(program))
            f.write('refrac file {}\n'.format(refrac_file_name))
            f.write('size_distribution_file {}\n'.format(size_dist_file_name))
            if nmom:
                f.write('nmom {}\n'.format(nmom))
            
    


# In[38]:


mie_input(fp_rtm+'mie_tester.inp',fp_rtm+'mie_tester.ref',fp_rtm+'mie_tester.psd')


# # Run through and make input files for scans

# In[39]:


base = 'mie_ORACLES2016_4wvl_expansion_{}'.format(vv)
f_list = fp_rtm+base+'_list.sh'
with open(f_list,'w') as f:
    for it,tt in list(enumerate(ae[b'time'])):
        basename = base+'_{:03.0f}'.format(it)
        print_refrac_file(fp_rtm+basename+'.ref',wave_out,
                         spline_wavelength_extend(ae[b'n_real'][it,:],ae[b'wavelength'],wave_out),
                         spline_wavelength_extend(ae[b'n_imag'][it,:],ae[b'wavelength'],wave_out))
        print_size_dist_file(fp_rtm+basename+'.psd',ae[b'radius'],convert_dvlnr_to_dndr(ae[b'psd'][it,:],ae[b'radius']))
        mie_input(fp_rtm+basename+'.inp',fp_rtm+basename+'.ref',fp_rtm+basename+'.psd',nmom=1500)
        f.write('{bin_path}mie < {inp} > {out}\n'.format(bin_path=fp_bin,inp=fp_rtm+basename+'.inp',out=fp_rtm+basename+'.out'))


# In[40]:


get_ipython().run_line_magic('ls', '$fp_rtm')


# In[41]:


get_ipython().system('parallel --jobs=22 --bar < $f_list')


# ## Read the mie output files

# In[42]:


fname = fp_rtm+basename+'.out'


# In[43]:


outs = np.genfromtxt(fname)


# In[45]:


len(outs[0,7:])


# In[46]:


def read_mie_output(fname):
    # read the default ascii output file from mie code
    outs = np.genfromtxt(fname)
    d = {}
    d['wavelength'] = outs[:,0]
    d['n_real'] = outs[:,1]
    d['n_imag'] = outs[:,2]
    d['qext'] = outs[:,3]
    d['ssa'] = outs[:,4]
    d['asym'] = outs[:,5]
    d['spike'] = outs[:,6]
    d['pmom'] = outs[:,7:]
    return d


# In[47]:


dats = []
for it,tt in list(enumerate(ae[b'time'])):
    basename = base+'_{:03.0f}'.format(it)
    dats.append(read_mie_output(fp_rtm+basename+'.out'))
    


# In[48]:


dats[0]['pmom'].shape


# ## Plot out the ssa and asym extrapolations

# In[49]:


fpt = '/data/sam/ORACLES/skyscan/Logan_20221118/'


# In[50]:


plt.figure()
it=0
plt.plot(ae[b'wavelength'],ae[b'SSA'][it,:],':',c=plt.cm.viridis(it*3),label='Original Skyscans')
plt.plot(wave_out,dats[it]['ssa'],'+-',c=plt.cm.viridis(it*3),label='Mie extrapolation')
for it,t in list(enumerate(ae[b'time'])):
    plt.plot(ae[b'wavelength'],ae[b'SSA'][it,:],':',c=plt.cm.viridis(it*3))
    plt.plot(wave_out,dats[it]['ssa'],'+-',c=plt.cm.viridis(it*3))
plt.legend()
sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis,norm=plt.Normalize(vmin=0,vmax=it))
plt.colorbar(sm,label='Skyscan number')
plt.title('SSA extrapolation of ORACLES 2016')
plt.xlabel('Wavelength [nm]')
plt.ylabel('SSA')
plt.savefig(fpt+'SSA_mie_extrapolation_ORACLES2016_{}.png'.format(vv),dpi=600,transparent=True)


# In[51]:


plt.figure()
it=0
plt.plot(ae[b'wavelength'],ae[b'g_total'][it,:],':',c=plt.cm.viridis(it*3),label='Original Skyscans')
plt.plot(wave_out,dats[it]['asym'],'+-',c=plt.cm.viridis(it*3),label='Mie extrapolation')
for it,t in list(enumerate(ae[b'time'])):
    plt.plot(ae[b'wavelength'],ae[b'g_total'][it,:],':',c=plt.cm.viridis(it*3))
    plt.plot(wave_out,dats[it]['asym'],'+-',c=plt.cm.viridis(it*3))
sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis,norm=plt.Normalize(vmin=0,vmax=it))
plt.colorbar(sm,label='Skyscan number')
plt.title('Asymmetry Parameter extrapolation of ORACLES 2016')
plt.xlabel('Wavelength [nm]')
plt.ylabel('Asymmetry Parameter')
plt.legend()
    
plt.savefig(fpt+'ASYM_mie_extrapolation_ORACLES2016_{}.png'.format(vv),dpi=600,transparent=True)


# ## Save as netcdf

# ### Copy file data and attributes from older one

# In[52]:


fpt = '/data/sam/ORACLES/skyscan/Logan_20221118/'
f_out = fpt + '4STAR-aeroinv4wvlexpansion_P3_2016_R1.nc'


# In[53]:


f_out


# In[54]:


f_in = fp + 'data_processed/starskies/skyscans_Logan_20221116/' +  '4STAR-aeroinv4wvl_P3_2016_R0.nc'


# In[55]:


toexclude = ['ExcludeVar1', 'ExcludeVar2']

with nc.Dataset(f_in) as src, nc.Dataset(f_out, "w") as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name not in toexclude:
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)
            dst[name][:] = src[name][:]
            


# ### Update the Dataset with the new calculated ssa/asym

# In[78]:


fn = nc.Dataset(f_out,'a')


# In[79]:


fn


# In[87]:


fn.getncattr('REVISION')


# In[59]:


fn.setncattr('REVISION',"R1")
fn.setncattr('R1',"Same measurement values as R0 4wvl skyscan results from Logan Mitchell. Updated by Samuel LeBlanc to include phase function moments of extrapolated values, and expanded wavelengths. ")
fn.setncattr('History',"Modified By Samuel LeBlanc, 2022-11-18: to add wavelength-extrapolated SSA and Asymmetry Parameter, pahser function using Mie calculations of size distribution and extrapolated index of refraction")


# In[60]:


fn.createDimension('Extrap_wavelength',len(dats[0]['wavelength']))
fn.createVariable('Extrap_wavelength','float64',('Extrap_wavelength'))
fn['Extrap_wavelength'].setncatts({'long_name':'Wavelengths from the extrapolated mie calculations','units':'nm'})
fn['Extrap_wavelength'][:] = dats[0]['wavelength'][:]


# In[61]:


extraps = {}


extraps['n_real'] = {'data':np.array([da['n_real'] for da in dats]),
                     'atts':{'long_name':'Real refractive index, extrapolated in wavelength by Mie calculations',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-11-18'}}
extraps['n_imag'] = {'data':np.array([da['n_imag'] for da in dats]),
                     'atts':{'long_name':'Imaginary refractive index, extrapolated in wavelength by Mie calculations',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-11-18'}}
extraps['qext'] = {'data':np.array([da['qext'] for da in dats]),
                     'atts':{'long_name':'extinction efficiency factor, extrapolated in wavelength by Mie calculations',
                             'units':'cm^3/m^3',
                             'history':'Built by Samuel LeBlanc on 2022-11-18'}}
extraps['ssa'] = {'data':np.array([da['ssa'] for da in dats]),
                     'atts':{'long_name':'Single Scattering Albedo, extrapolated in wavelength by Mie calculations',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-11-18'}}
extraps['asym'] = {'data':np.array([da['asym'] for da in dats]),
                     'atts':{'long_name':'Asymmetry Parameter, extrapolated in wavelength by Mie calculations',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-11-18'}}
extraps['pmom'] = {'data':np.array([da['pmom'] for da in dats]),
                     'atts':{'long_name':'Phase function Legendre moments (km) for the phase function(p(µ)) reconstruction using: p(µ) = Sum((2m + 1) · km · Pm(µ)) for m=0 to infinity',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-11-18'}}


# In[62]:


len(dats[0]['pmom'][0,:])


# In[63]:


fn.createDimension('Extrap_nmom',len(dats[0]['pmom'][0,:]))
fn.createVariable('Extrap_nmom','int',('Extrap_nmom'))
fn['Extrap_nmom'].setncatts({'long_name':'Moment number for phase function moments','units':'None'})
fn['Extrap_nmom'][:] = np.arange(0,len(dats[0]['pmom'][0,:]))


# In[64]:


for k in extraps:
    if k == 'pmom':
        fn.createVariable('Extrap_'+k,'float',('time','Extrap_wavelength','Extrap_nmom'))
        fn['Extrap_'+k].setncatts(extraps[k]['atts'])
        fn['Extrap_'+k][:] = extraps[k]['data']
    else:
        fn.createVariable('Extrap_'+k,'float64',('time','Extrap_wavelength'))
        fn['Extrap_'+k].setncatts(extraps[k]['atts'])
        fn['Extrap_'+k][:] = extraps[k]['data']


# In[199]:


if vv=='v3':
    k = 'QC_flags'
    extraps = {}
    extraps['QC_flags'] = {'data':qc_flags,
                     'atts':{'long_name':'QC of flags, lower is better. From Pistone et al., published on QC values below 1. evaluated on AOD_400>0.2, max scattering angle measured larger than 90, and low mean sky error (below~5)',
                             'units':'None',
                             'history':'Built by Pistone et al., and incorporated into file on 2022-10-04'}}

    fn.createVariable('QC_flags','float',('time'))
    fn['QC_flags'].setncatts(extraps[k]['atts'])
    fn['QC_flags'][:] = extraps[k]['data']


# In[65]:


fn


# In[93]:


fn['Lat'].setncatts({'long_name':'Latitude'})


# In[94]:


fn.close()


# # Calculate the phase functions

# In[67]:


import subprocess


# In[68]:


base = 'mie_ORACLES2016_4wvl_expansion_{}'.format(vv)
f_list = fp_rtm+base+'_list_for_phase.sh'
with open(f_list,'w') as f:
    for it,tt in list(enumerate(ae[b'time'])):
        basename = base+'_{:03.0f}'.format(it)
        # loop through and copy file to delete first line
        for i,w in list(enumerate(wave_out)):
            if i==0:
                subprocess.call("sed '/^$/d' "+fp_rtm+basename+".out > "+fp_rtm+basename+"_w{:1.0f}.out".format(i),shell=True) # remove empty lines
            else:
                subprocess.call("sed '1d' "+fp_rtm+basename+"_w{:1.0f}.out > ".format(i-1)+fp_rtm+basename+"_w{:1.0f}.out".format(i),shell=True)
            f.write('{bin_path}phase {inp} > {out}\n'.format(bin_path=fp_bin,inp=fp_rtm+basename+'_w{:1.0f}.out'.format(i),out=fp_rtm+basename+'_phase_w{:1.0f}.out'.format(i)))


# In[69]:


wave_out


# In[70]:


get_ipython().system('parallel --jobs=22 --bar < $f_list')


# ## Read the output phase function files

# In[71]:


ph1 = np.genfromtxt(fp_rtm+basename+'_phase_w{:1.0f}.out'.format(i))


# In[72]:


phase = np.zeros((len(ae[b'time']),len(wave_out),len(ph1[:,0])))+np.nan
iphase = np.zeros((len(ae[b'time']),len(wave_out),len(ph1[:,0])))+np.nan
mu = ph1[:,0]
scat_angle = np.arccos(mu)*180.0/np.pi


# In[73]:


base = 'mie_ORACLES2016_4wvl_expansion_{}'.format(vv)
for it,tt in list(enumerate(ae[b'time'])):
    basename = base+'_{:03.0f}'.format(it)
    for i,w in list(enumerate(wave_out)):
        ph = np.genfromtxt(fp_rtm+basename+'_phase_w{:1.0f}.out'.format(i))
        phase[it,i,:] = ph[:,1]
        iphase[it,i,:] = ph[:,2]


# In[74]:


phase.shape


# In[75]:


from mpltools import color


# In[76]:


for i,w in list(enumerate(wave_out)):
    plt.figure()
    color.cycle_cmap(len(ae[b'time']),cmap=plt.cm.viridis,ax=plt.gca())
    plt.plot(scat_angle,phase[:,i,:].T)
    plt.xlabel('Scattering Angle [deg]')
    plt.ylabel('Phase function')
    plt.yscale('log')
    plt.title('Wavelength: {:4.1f} nm'.format(w))

    scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    scalarmap.set_array(range(len(ae[b'time'])))
    cba = plt.colorbar(scalarmap,label='skyscan number')


# # Load the calculated reflectances with/withaero for presentation

# In[6]:


import Run_libradtran as rl


# In[10]:


fo = '/scratch/rtm/input/TASNPP_v1/'
zout=[0.2,1.5,100.0]


# In[14]:


aw0 = rl.read_libradtran(fo+'ORACLES_cloud_refl_w0.out',zout=zout)
aw1 = rl.read_libradtran(fo+'ORACLES_cloud_refl_w1.out',zout=zout)
nw0 = rl.read_libradtran(fo+'ORACLES_cloud_refl_w0_noaero.out',zout=zout)
nw1 = rl.read_libradtran(fo+'ORACLES_cloud_refl_w1_noaero.out',zout=zout)


# In[17]:


aw0.keys()


# In[18]:


aw0['direct_down'].shape


# In[20]:


refl_a = np.hstack(aw0['direct_down'],aw1['direct_down'])
#refl_n = 


# In[26]:


refl_a = np.vstack([aw0['diffuse_up'],aw1['diffuse_up']])/(np.vstack([aw0['direct_down'],aw1['direct_down']])+np.vstack([aw0['diffuse_down'],aw1['diffuse_down']]))
refl_n = np.vstack([nw0['diffuse_up'],nw1['diffuse_up']])/(np.vstack([nw0['direct_down'],nw1['direct_down']])+np.vstack([nw0['diffuse_down'],nw1['diffuse_down']]))
wvls = np.hstack([aw0['wvl'], aw1['wvl']])


# In[30]:


up_a = np.vstack([aw0['diffuse_up'],aw1['diffuse_up']])
up_n = np.vstack([nw0['diffuse_up'],nw1['diffuse_up']])


# ## plot the reflectances

# In[27]:


refl_a.shape


# In[49]:


plt.figure(figsize=(6,3))
plt.plot(wvls,refl_a[:,2],label='With Absorbing aerosol')
plt.plot(wvls,refl_n[:,2],label='No aerosol')
plt.plot(wvls,refl_n[:,2]-refl_a[:,2],label='Difference (no aero - aero)')
plt.legend()
plt.grid()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Reflectances at TOA')
plt.tight_layout()
plt.savefig(fo+'ORACLES_cloud_Reflectance_with_without_aero.png',dpi=600,transparent=True)


# In[50]:


plt.figure(figsize=(6,3))
plt.plot(wvls,refl_a[:,2],label='With Absorbing aerosol')
plt.plot(wvls,refl_n[:,2],label='No aerosol')
plt.plot(wvls,refl_n[:,2]-refl_a[:,2],label='Difference (no aero - aero)')
plt.legend()
for x in [470,550,660,860,1240,1630,2100]:
    plt.axvline(x,color='grey',linestyle='--',lw=3)
plt.grid()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Reflectances at TOA')
plt.tight_layout()
plt.savefig(fo+'ORACLES_cloud_Reflectance_with_without_aero_lines.png',dpi=600,transparent=True)


# In[51]:


fo


# # Load and plot the relative ACAOD spectra

# In[7]:


ff = '/data/sam/ORACLES/'


# In[8]:


a6 = hs.loadmat(ff+'ORACLES2016_4STAR_ACAOD_relative_to_550.mat')
a7 = hs.loadmat(ff+'ORACLES2017_4STAR_ACAOD_relative_to_550.mat')
a8 = hs.loadmat(ff+'ORACLES2018_4STAR_ACAOD_relative_to_550.mat')


# In[9]:


a6.keys()


# In[10]:


print(a6['details'])


# In[11]:


wla = [470,550,675,865,1236]
a6['aod'] = np.vstack([a6['AOD0470_rel']*a6['AOD0550_abs'],a6['AOD0550_abs'],a6['AOD0675_rel']*a6['AOD0550_abs'],a6['AOD0865_rel']*a6['AOD0550_abs'],a6['AOD1236_rel']*a6['AOD0550_abs']])
a7['aod'] = np.vstack([a7['AOD0470_rel']*a7['AOD0550_abs'],a7['AOD0550_abs'],a7['AOD0675_rel']*a7['AOD0550_abs'],a7['AOD0865_rel']*a7['AOD0550_abs'],a7['AOD1236_rel']*a7['AOD0550_abs']])
a8['aod'] = np.vstack([a8['AOD0470_rel']*a8['AOD0550_abs'],a8['AOD0550_abs'],a8['AOD0675_rel']*a8['AOD0550_abs'],a8['AOD0865_rel']*a8['AOD0550_abs'],a8['AOD1236_rel']*a8['AOD0550_abs']])


# In[12]:


a6['aod'].shape


# In[71]:


plt.figure()
plt.plot(wla,a6['aod'],'.-',color='grey',lw=0.1)
plt.ylim(0,1.3)


# In[13]:


aod6 = np.ma.masked_array(a6['aod'].T,np.isnan(a6['aod'].T))
aod7 = np.ma.masked_array(a7['aod'].T,np.isnan(a7['aod'].T))
aod8 = np.ma.masked_array(a8['aod'].T,np.isnan(a8['aod'].T))


# In[16]:


aod6.shape


# In[61]:


plt.figure(figsize=(5,3))
plt.plot(wla,np.nanmean(aod6,axis=0),'o-',label='mean')
plt.plot(wla,np.nanmedian(a6['aod'].T,axis=0),'+-',label='median')
plt.plot(wla,np.nanpercentile(a6['aod'].T,75,axis=0),'^:',label='75% percentile')
plt.plot(wla,np.nanpercentile(a6['aod'].T,25,axis=0),'v:',label='25% percentile')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Above Cloud AOD')
plt.title('ORACLES 2016')
plt.tight_layout()
plt.savefig(ff+'ORACLES_2016_acaod_avg_for_TASNPP.png',dpi=600,transparent=True)


# In[67]:


plt.figure(figsize=(5,3))
plt.plot(wla,np.nanmean(a6['aod']/a6['AOD0550_abs'],axis=1),'o-',label='mean')
plt.plot(wla,np.nanmedian(a6['aod']/a6['AOD0550_abs'],axis=1),'+-',label='median')
plt.plot(wla,np.nanpercentile(a6['aod']/a6['AOD0550_abs'],75,axis=1),'^:',label='75% percentile')
plt.plot(wla,np.nanpercentile(a6['aod']/a6['AOD0550_abs'],25,axis=1),'v:',label='25% percentile')
plt.axhline(1,ls='--',color='grey')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Above Cloud AOD Relative to 550 nm')
plt.title('ORACLES 2016')
plt.tight_layout()
plt.savefig(ff+'ORACLES_2016_acaod_avg_for_TASNPP_rel550.png',dpi=600,transparent=True)


# In[80]:


qe_wv = [358.80689111129686,  434.3750056895338, 571.5909028841447, 804.2613372576161, 1046.8750739639393, 1233.8067545624858, 1683.2387765330593, 549.7159654690166]
qe = [0.1507455703364354,  0.12325255219252594, 0.09832245607751045, 0.0640726880950128, 0.047064292802888326, 0.0391425919799973, 0.027027020685372494, 0.10205032909070726]


# In[81]:


qe_rel = np.array(qe)/np.array(qe[-1])
isort = np.argsort(qe_wv)
qe_wv = np.array(qe_wv)[isort]
qe_rel = qe_rel[isort]


# In[82]:


plt.figure(figsize=(5,3))
plt.plot(wla,np.nanmean(a6['aod']/a6['AOD0550_abs'],axis=1),'o-',label='mean')
plt.plot(wla,np.nanmedian(a6['aod']/a6['AOD0550_abs'],axis=1),'+-',label='median')
plt.plot(wla,np.nanpercentile(a6['aod']/a6['AOD0550_abs'],75,axis=1),'^:',label='75% percentile')
plt.plot(wla,np.nanpercentile(a6['aod']/a6['AOD0550_abs'],25,axis=1),'v:',label='25% percentile')
plt.plot(qe_wv,qe_rel,'-k',label='Relative Qe from Haywood et al., 2000')
plt.axhline(1,ls='--',color='grey')
plt.legend()
plt.xlabel('Wavelength [nm]')
plt.ylabel('Above Cloud AOD Relative to 550 nm')
plt.title('ORACLES 2016')
plt.tight_layout()
plt.savefig(ff+'ORACLES_2016_acaod_avg_for_TASNPP_rel550_haywood.png',dpi=600,transparent=True)


# In[ ]:




