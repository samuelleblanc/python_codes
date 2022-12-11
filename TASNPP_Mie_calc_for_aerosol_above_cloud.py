#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     To Calculate aerosol optical properties for aerosol above cloud reitrewvals using MODIS and VIIRS
#     Using the wavelengths: 
#     0.47, 0.55, 0.67, 0.86, 1.24, 2.1µm
#     
#     - Using the retrieved size distributions
#     - retrieval results of refractive index (imaginary and real) at wavelengths: 400, 500, 675, 870, 995 nm
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
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-09-12
#     Modified:Samuel LeBlanc, Santa Cruz, CA, 2022-09-21
#              - added new wavelengths and storing of the phase function
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
fp_rtm = fp_rtm0 +'TASNPP_mie/'


# In[5]:


if not os.path.exists(fp_rtm): 
    os.mkdir(fp_rtm)


# # Load files

# In[6]:


f = fp + 'data_archival/AEROINV_nc_R0/NC_FORMAT_NETCDF4_CLASSIC/4STAR-aeroinv_P3_2016_R0.nc'


# In[7]:


ae,ae_dict = lu.load_netcdf(f,everything=True)


# In[8]:


len(ae[b'time'])


# In[9]:


ae_dict[b'scan_tag']


# In[10]:


ae[b'scan_tag']


# ## Extra meta information

# In[11]:


sk_meta = sio.loadmat(fp+'data_processed/starskies/skyscans_ORACLES2016_moredata.mat')


# In[12]:


sk_meta.keys()


# In[13]:


sk_meta['skyscans_meta'].shape


# In[14]:


sk_meta['skyscans_qcflags'].shape


# In[15]:


sk_meta['skyscan_vars'][0]


# In[16]:


sk_meta['skyscans_meta'][0,:]


# In[17]:


sk_meta['skyscans_all'].shape


# In[18]:


sk_meta['skyscans_all'][0,24]


# In[19]:


sk_meta['days'] = [dd[0][0] for dd in sk_meta['skyscans_all'][:,0]]


# In[20]:


sk_meta['filenum'] = [dd[0][0] for dd in sk_meta['skyscans_all'][:,1]]


# ## Map the extra meta information to the skyscans in the file

# In[21]:


ae[b'scan_tag']


# In[22]:


from datetime import datetime


# In[23]:


days = []
for t in ae[b'time']:
    x =  datetime.fromtimestamp(datetime(2016,8,31).timestamp()+t)
    days.append(x.year*10000+x.month*100+x.day)


# In[24]:


days


# In[25]:


qc_flags = [4.0]*len(days)
for i,da in list(enumerate(days)):
    ig = [(da==dk)&(int(ae[b'scan_tag'][i])==sk_meta['filenum'][j]) for j,dk in enumerate(sk_meta['days'])] 
    if any(ig): qc_flags[i] = (sk_meta['skyscans_qcflags'][ig][0][0])


# In[26]:


qc_flags


# # Run through each retrieval and make input files

# ## Prep functions for printing size distribution

# From the libradtran documentation for mie size distribution file:  
# >Specify a two column file, r [micron], dn(r)/dr, which describes a size distribution

# In[12]:


def print_size_dist_file(fname,r,dnr):
    with open(fname,'w') as f:
        for ir,rr in list(enumerate(r)):
            f.write('{:3.10f} {:3.10f}\n'.format(rr,dnr[ir]))


# In[13]:


ae_dict[b'radius']


# In[14]:


ae_dict[b'psd']


# In[15]:


def convert_dvlnr_to_dndr(psd,r):
     # All about that conversion from the volume size distribution of dV(r)/dln(r) to number size distribution dN(r)/dr
    Nr = psd/(4.0*np.pi/3.0)
    for i,rr in list(enumerate(r)):
        Nr[i] = Nr[i]/rr**4.0
    return Nr


# In[16]:


print_size_dist_file(fp_rtm+'mie_tester.psd',ae[b'radius'],convert_dvlnr_to_dndr(ae[b'psd'][0,:],ae[b'radius']))


# ## Prep function for printing index of refraction

# Update from Kerry Meyer's email on Sept 21, 2022.  
# 
# > Thanks, Sam, this is helpful. Do you have the extrapolated phase functions also? And I hate to do this, but can the wavelength range be extended in both directions? I’d like to accommodate the VIIRS 2.25µm channel, and use the spectral response functions to compute band-weighted averages for each MODIS and VIIRS channel. So perhaps adding something like 400nm on the short end (which I guess is the exact 4STAR wavelength there) and 2300nm on the long end – that should encompass most of the response function ranges for MODIS and VIIRS.
# >
# > I’ll also need your input on how best to QA filter these scans. I see a “QA_level” dataset that I assume should be part of that.
# >
# >Kerry

# In[17]:


ae_dict[b'n_real']


# In[18]:


ae_dict[b'n_imag']


# In[19]:


ae_dict[b'wavelength']


# In[31]:


def spline_wavelength_extend(val,wavelen,new_wavelen,su=0.0006,k=2):
    # to calculate a new spline fit to the refractive index
    val_fx = si.CubicSpline(np.append(wavelen,wavelen[-1]+600.0),np.append(val,val[-1]),bc_type='natural',extrapolate=True)
    return val_fx(new_wavelen)


# In[29]:


if vv == 'v1':
    wave_out = np.array([470.0, 550.0, 670.0, 860.0, 1240.0, 2100.0])
elif vv == 'v2':
    wave_out = np.array([400.0, 470.0, 550.0, 670.0, 860.0, 1240.0, 2100.0, 2300.0])


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

# In[41]:


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

# In[37]:


vv2 = 'v2'
base = 'mie_ORACLES2016_expansion_{}'.format(vv2)


# In[42]:


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


# In[43]:


get_ipython().run_line_magic('ls', '$fp_rtm')


# In[44]:


get_ipython().system('parallel --jobs=22 --bar < $f_list')


# ## Read the mie output files

# In[33]:


fname = fp_rtm+basename+'.out'


# In[46]:


outs = np.genfromtxt(fname)


# In[50]:


len(outs[0,7:]


# In[35]:


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


# In[38]:


dats = []
for it,tt in list(enumerate(ae[b'time'])):
    basename = base+'_{:03.0f}'.format(it)
    dats.append(read_mie_output(fp_rtm+basename+'.out'))
    


# In[39]:


dats[0]['pmom'].shape


# ## Plot out the ssa and asym extrapolations

# In[90]:


fpt = '/data/sam/TASNPP/ORACLES_aerosol_prop/'


# ### SSA

# In[40]:


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


# In[41]:


SSA_mie = np.zeros((len(ae[b'time']),len(wave_out)))
for it,t in list(enumerate(ae[b'time'])):
    SSA_mie[it,:] = dats[it]['ssa']


# In[62]:


import plotting_utils
from plotting_utils import make_boxplot, color_box


# In[71]:


def coloring_box(bo,color='tab:blue',alpha=0.2,label='',mean_marker='s'):
    color_box(bo,color)
    for n in list(bo.keys()):
        nul = [plt.setp(bo[n][idx],alpha=alpha)for idx in range(len(bo[n]))]
    v = [plt.setp(bo['means'][idx],alpha=0.05)for idx in range(len(bo['means']))]
    mean = [a.get_xdata()[0] for a in bo['means']]
    pos = [a.get_ydata()[0] for a in bo['means']]
    plt.plot(mean,pos,mean_marker+'-',zorder=100,color=color,label=label,lw=2.5,alpha=alpha)


# In[93]:


ssa_wv = np.array([307.1024108277745, 2134.6589491879768, 527.8410280538882, 225.56820561077757, 386.6477371001037, 684.9431657840408, 840.0567886996887, 1036.9317716310875, 1150.2840345309835, 1472.443097509636, 460.22733686383543])
ssa_hay = np.array([0.9207865168539326, 0.8999999828552931, 0.9056179775280899, 0.9235954884732708, 0.9168539154395628, 0.8960674071579837, 0.8848314435294505, 0.8780898704957427, 0.8758426880568602, 0.8797752637541696, 0.910112342405855])
issa = np.argsort(ssa_wv)
ssa_wv = ssa_wv[issa]
ssa_hay = ssa_hay[issa]


# In[108]:


wvm = np.array([470,550,660,860,1240,1630,2100])
asy_mod4 = np.array([0.64,0.6,0.56,0.49,0.48,0.56,0.64]) 
asy_hay = np.array([0.65,0.61,0.56,0.49,0.51,0.58,0.63])
ssa_mod4 = np.array([0.88,0.87,0.85,0.80,0.73,0.7,0.7])
ssa_hay = np.array([0.92,0.91,0.90,0.89,0.87,0.88,0.90])
qe_mod4 = np.array([0.156,0.127,0.107,0.077,0.052,0.040,0.032])/(1.0-ssa_mod4)
qe_hay = np.array([0.154,0.123,0.098,0.064,0.038,0.027,0.020])/(1.0-ssa_hay)


# In[118]:


plt.figure(figsize=(5,3))
#bpm = plt.boxplot(SSA_mie,positions=wave_out,showfliers=False,widths=30,showmeans=True,patch_artist=True)
#plt.plot(wave_out,SSA_mie.T,':',alpha=0.1,color='tab:blue')

bpa = plt.boxplot(ae[b'SSA'],positions=ae[b'wavelength'],showfliers=False,widths=30,showmeans=True,patch_artist=True)
plt.plot(ae[b'wavelength'],ae[b'SSA'].T,'--',alpha=0.1,color='tab:orange')

plt.plot(wvm,ssa_hay,'^-',color='grey',label='SAFARI: Haywood et al., 2003')
plt.plot(wvm,ssa_mod4,'v-',color='grey',label='MOD04: Levy et al., 2009')

#coloring_box(bpm,'tab:blue',label='Mie extrapolation',alpha=0.5)
coloring_box(bpa,'tab:orange',label='ORACLES 2016 SSA: Pistone et al. 2019',alpha=0.5)



plt.legend(loc=3)
plt.xlim(350,1300)
plt.xticks([400,500,550,675,860,995,1240],[400,500,550,675,860,995,1240])

plt.ylim(0.65,0.95)
plt.ylabel('SSA')
plt.xlabel('Wavelength [nm]')
plt.savefig(fpt+'SSA_ORACLES2016_{}.png'.format(vv),dpi=600,transparent=True)


# In[119]:


plt.figure(figsize=(5,3))


bpa = plt.boxplot(ae[b'SSA'],positions=ae[b'wavelength'],showfliers=False,widths=30,showmeans=True,patch_artist=True)
plt.plot(ae[b'wavelength'],ae[b'SSA'].T,'--',alpha=0.1,color='tab:orange')
bpm = plt.boxplot(SSA_mie,positions=wave_out,showfliers=False,widths=30,showmeans=True,patch_artist=True)
plt.plot(wave_out,SSA_mie.T,':',alpha=0.1,color='tab:blue')

plt.plot(wvm,ssa_hay,'^-',color='grey',label='SAFARI: Haywood et al., 2003')
plt.plot(wvm,ssa_mod4,'v-',color='grey',label='MOD04: Levy et al., 2009')


coloring_box(bpa,'tab:orange',label='ORACLES 2016 SSA: Pistone et al. 2019',alpha=0.5)
coloring_box(bpm,'tab:blue',label='Mie extrapolation',alpha=0.5)
plt.legend(loc=3)
plt.xlim(350,1300)
plt.xticks([400,500,550,675,860,995,1240],[400,500,550,675,860,995,1240])

plt.ylim(0.65,0.95)
plt.ylabel('SSA')
plt.xlabel('Wavelength [nm]')
plt.savefig(fpt+'SSA_ORACLES2016_with_extrap_{}.png'.format(vv),dpi=600,transparent=True)


# ### Asym

# In[ ]:





# In[58]:


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


# In[99]:


ASY_mie = np.zeros((len(ae[b'time']),len(wave_out)))
for it,t in list(enumerate(ae[b'time'])):
    ASY_mie[it,:] = dats[it]['asym']


# In[ ]:


x = [321.02274278963773, 2246.022697273367, 549.7159654690166, 247.44332509098695, 374.7159199527464, 645.1705026478767, 782.3863998424878, 1046.8750739639393, 1271.5909028841452, 1647.443325090987, 434.3750056895338]
y = [0.658988746900237, 0.6286516682485517, 0.569101089306092, 0.6865168539325843, 0.6325842696629214, 0.5393258255519224, 0.4960674157303371, 0.4853932412822595, 0.5162921348314607, 0.5786516853932584, 0.610112342405855]


# In[105]:


asy_wvl = np.array([321.02274278963773, 2246.022697273367, 549.7159654690166, 247.44332509098695, 374.7159199527464, 645.1705026478767, 782.3863998424878, 1046.8750739639393, 1271.5909028841452, 1647.443325090987, 434.3750056895338])
asy_hay = np.array([0.658988746900237, 0.6286516682485517, 0.569101089306092, 0.6865168539325843, 0.6325842696629214, 0.5393258255519224, 0.4960674157303371, 0.4853932412822595, 0.5162921348314607, 0.5786516853932584, 0.610112342405855])
iasy = np.argsort(asy_wvl)
asy_wvl = asy_wvl[iasy]
asy_hay[iasy] = asy_hay[iasy]


# In[120]:


plt.figure(figsize=(5,3))
#bpm = plt.boxplot(SSA_mie,positions=wave_out,showfliers=False,widths=30,showmeans=True,patch_artist=True)
#plt.plot(wave_out,SSA_mie.T,':',alpha=0.1,color='tab:blue')

bpa = plt.boxplot(ae[b'g_total'],positions=ae[b'wavelength'],showfliers=False,widths=30,showmeans=True,patch_artist=True)
plt.plot(ae[b'wavelength'],ae[b'g_total'].T,'--',alpha=0.1,color='tab:orange')

plt.plot(wvm,asy_hay,'^-',color='grey',label='SAFARI: Haywood et al., 2003')
plt.plot(wvm,asy_mod4,'v-',color='grey',label='MOD04: Levy et al., 2009')

#coloring_box(bpm,'tab:blue',label='Mie extrapolation',alpha=0.5)
coloring_box(bpa,'tab:orange',label='ORACLES 2016 ASY: Pistone et al. 2019',alpha=0.5)



plt.legend()
plt.xlim(350,1300)
plt.xticks([400,500,550,675,860,995,1240],[400,500,550,675,860,995,1240])

plt.ylim(0.45,0.8)
plt.ylabel('Asymmetry Parameter')
plt.xlabel('Wavelength [nm]')
plt.savefig(fpt+'ASY_ORACLES2016_{}.png'.format(vv),dpi=600,transparent=True)


# In[121]:


plt.figure(figsize=(5,3))
bpm = plt.boxplot(ASY_mie,positions=wave_out,showfliers=False,widths=30,showmeans=True,patch_artist=True)
plt.plot(wave_out,ASY_mie.T,':',alpha=0.1,color='tab:blue')

bpa = plt.boxplot(ae[b'g_total'],positions=ae[b'wavelength'],showfliers=False,widths=30,showmeans=True,patch_artist=True)
plt.plot(ae[b'wavelength'],ae[b'g_total'].T,'--',alpha=0.1,color='tab:orange')

plt.plot(wvm,asy_hay,'^-',color='grey',label='SAFARI: Haywood et al., 2003')
plt.plot(wvm,asy_mod4,'v-',color='grey',label='MOD04: Levy et al., 2009')

coloring_box(bpm,'tab:blue',label='Mie extrapolation',alpha=0.5)
coloring_box(bpa,'tab:orange',label='ORACLES 2016 ASY: Pistone et al. 2019',alpha=0.5)



plt.legend()
plt.xlim(350,1300)
plt.xticks([400,500,550,675,860,995,1240],[400,500,550,675,860,995,1240])

plt.ylim(0.45,0.8)
plt.ylabel('Asymmetry Parameter')
plt.xlabel('Wavelength [nm]')
plt.savefig(fpt+'ASY_ORACLES2016_with_extrap_{}.png'.format(vv),dpi=600,transparent=True)


# ## Plot out the size distribution

# In[ ]:


ae[b'radius'],convert_dvlnr_to_dndr(ae[b'psd'][0,:],ae[b'radius'])


# In[124]:


ae_dict[b'radius']


# In[136]:


plt.figure(figsize=(3,2))
plt.plot(ae[b'radius'],ae[b'psd'].T,lw=0.5)
plt.plot(ae[b'radius'],np.nanmean(ae[b'psd'].T,axis=1),'-k',lw=2.5,label='Mean')
plt.plot(ae[b'radius'],np.nanmedian(ae[b'psd'].T,axis=1),'--k',lw=2.5,label='Median')

plt.xscale('log')
plt.xlabel('Radius [micron]')
plt.ylabel('Size distribution dV/dln(r)')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig(fpt+'PSD_ORACLES2016_{}.png'.format(vv),dpi=600,transparent=True)


# ## Save as netcdf

# ### Copy file data and attributes from older one

# In[184]:


fpt = '/data/sam/TASNPP/ORACLES_aerosol_prop/'
f_out = fpt + '4STAR-aeroinv_mie_wavelength_expansion_P3_2016_R2.nc'


# In[185]:


f_out


# In[186]:


f_in = fpt +  '4STAR-aeroinv_mie_wavelength_expansion_P3_2016_R1.nc'


# In[187]:


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

# In[195]:


fn = nc.Dataset(f_out,'a')


# In[196]:


fn


# In[197]:


fn.getncattr('REVISION')


# In[198]:


if vv=='v2':
    fn.setncattr('REVISION',"R1")
    fn.setncattr('R1',"Same measurement values as R0. Update to include phase function moments of extrapolated values, and expanded wavelengths")
    fn.setncattr('History',"Modified By Samuel LeBlanc, 2022-09-21: Expanded wavelengths extrapolation and included phase function moments.\nModified By Samuel LeBlanc, 2022-09-14, to add wavelength-extrapolated SSA and Asymmetry Parameter, using Mie calculations of size distribution and extrapolated index of refraction")
elif vv=='v3':
    fn.setncattr('REVISION',"R2")
    fn.setncattr('R1',"Same measurement values as R0. Update to include phase function moments of extrapolated values, and expanded wavelengths")
    fn.setncattr('R2',"Same measurement and expanded phase moments and extrapolated values from R1. Update to include QC flags from Pistone et al.")
    fn.setncattr('History',"Modified By Samuel LeBlanc, 2022-10-04: Added QC flags from Pistone et al., 2019 .\n Modified By Samuel LeBlanc, 2022-09-21: Expanded wavelengths extrapolation and included phase function moments.\nModified By Samuel LeBlanc, 2022-09-14, to add wavelength-extrapolated SSA and Asymmetry Parameter, using Mie calculations of size distribution and extrapolated index of refraction")
else:
    fn.setncattr('History',"Modified By Samuel LeBlanc, 2022-09-14, to add wavelength-extrapolated SSA and Asymmetry Parameter, using Mie calculations of size distribution and extrapolated index of refraction")
    


# In[69]:


fn.createDimension('Extrap_wavelength',len(dats[0]['wavelength']))
fn.createVariable('Extrap_wavelength','float64',('Extrap_wavelength'))
fn['Extrap_wavelength'].setncatts({'long_name':'Wavelengths from the extrapolated mie calculations','units':'nm'})
fn['Extrap_wavelength'][:] = dats[0]['wavelength'][:]


# In[78]:


extraps = {}


extraps['n_real'] = {'data':np.array([da['n_real'] for da in dats]),
                     'atts':{'long_name':'Real refractive index, extrapolated in wavelength by Mie calculations',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-09-21'}}
extraps['n_imag'] = {'data':np.array([da['n_imag'] for da in dats]),
                     'atts':{'long_name':'Imaginary refractive index, extrapolated in wavelength by Mie calculations',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-09-21'}}
extraps['qext'] = {'data':np.array([da['qext'] for da in dats]),
                     'atts':{'long_name':'extinction efficiency factor, extrapolated in wavelength by Mie calculations',
                             'units':'cm^3/m^3',
                             'history':'Built by Samuel LeBlanc on 2022-09-21'}}
extraps['ssa'] = {'data':np.array([da['ssa'] for da in dats]),
                     'atts':{'long_name':'Single Scattering Albedo, extrapolated in wavelength by Mie calculations',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-09-21'}}
extraps['asym'] = {'data':np.array([da['asym'] for da in dats]),
                     'atts':{'long_name':'Asymmetry Parameter, extrapolated in wavelength by Mie calculations',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-09-21'}}
extraps['pmom'] = {'data':np.array([da['pmom'] for da in dats]),
                     'atts':{'long_name':'Phase function Legendre moments (km) for the phase function(p(µ)) reconstruction using: p(µ) = Sum((2m + 1) · km · Pm(µ)) for m=0 to infinity',
                             'units':'None',
                             'history':'Built by Samuel LeBlanc on 2022-09-21'}}


# In[79]:


len(dats[0]['pmom'][0,:])


# In[80]:


if vv=='v2':
    fn.createDimension('Extrap_nmom',len(dats[0]['pmom'][0,:]))
    fn.createVariable('Extrap_nmom','int',('Extrap_nmom'))
    fn['Extrap_nmom'].setncatts({'long_name':'Moment number for phase function moments','units':'None'})
    fn['Extrap_nmom'][:] = np.arange(0,len(dats[0]['pmom'][0,:]))


# In[81]:


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


# In[200]:


fn


# In[201]:


fn.close()


# # Calculate the phase functions

# In[88]:


import subprocess


# In[101]:


base = 'mie_ORACLES2016_expansion_{}'.format(vv)
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


# In[87]:


wave_out


# In[102]:


get_ipython().system('parallel --jobs=22 --bar < $f_list')


# ## Read the output phase function files

# In[ ]:


ph1 = np.genfromtxt(fp_rtm+basename+'_phase_w{:1.0f}.out'.format(i))


# In[122]:


phase = np.zeros((len(ae[b'time']),len(wave_out),len(ph1[:,0])))+np.nan
iphase = np.zeros((len(ae[b'time']),len(wave_out),len(ph1[:,0])))+np.nan
mu = ph1[:,0]
scat_angle = np.arccos(mu)*180.0/np.pi


# In[123]:


base = 'mie_ORACLES2016_expansion_{}'.format(vv)
for it,tt in list(enumerate(ae[b'time'])):
    basename = base+'_{:03.0f}'.format(it)
    for i,w in list(enumerate(wave_out)):
        ph = np.genfromtxt(fp_rtm+basename+'_phase_w{:1.0f}.out'.format(i))
        phase[it,i,:] = ph[:,1]
        iphase[it,i,:] = ph[:,2]


# In[125]:


phase.shape


# In[129]:


from mpltools import color


# In[131]:


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


# In[ ]:




