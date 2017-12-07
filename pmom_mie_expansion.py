
# coding: utf-8

# # Info
# Name:  
# 
#     pmom_mie_expansion
# 
# Purpose:  
# 
#     Explore the pmom expansion of mie values, with the delta-m scaling factor for the legendre polynomial expansion
#     And testing of the tools within libradtran
#   
# Input:
# 
#     none at command line
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
#     
#   
# Needed Files:
# 
#     - mie_hi.out
#     - phase input and output files from seperate deltam runs and phase runs
#   
#   
# Modification History:
# 
#     Wrtten: Samuel LeBlanc, NASA Ames, from Santa Cruz, 2017-12-06
#     Modified: 

# # Load the required modules and set up the paths

# In[1]:


import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
get_ipython().magic(u'matplotlib notebook')
from path_utils import getpath


# In[4]:


fp = getpath('libradtran',make_path=True,path='/mnt/c/Users/sleblanc/Research/libradtran')
fp


# In[5]:


fs = getpath('4STAR',make_path=True,path='/mnt/c/Users/sleblanc/Research/4STAR')
fs


# # Load the files

# ## Load the original mie_hi.out file with all the pahse functions and non-scaled legendre polynomials

# In[6]:


m = sio.idl.readsav(fs+'rtm_dat/mie_hi.out')


# In[7]:


m.keys()


# In[27]:


with open(fs+'rtm_dat/phase_an.txt','w') as f:
    f.write('{:4.0f}\n'.format(m['ntheta'][25,200]))
    for i in xrange(m['ntheta'][25,200],0,-1):
        f.write('{} {}\n'.format(m['theta'][25,200,i-1],m['phase'][25,200,i-1]))


# ## Load the calculated legendre polynomials with delta-m scaling factor
# After subjecting the phase_an.txt to deltam program from Y. Hu's program compiled to deltam:
# 
# ./deltam3 ../../4STAR/rtm_dat/phase_an.txt pph.out 300
# 
# #sed -i '1d' pph.out

# In[48]:


pmom = np.genfromtxt(fp+'ice/pph.out')


# In[49]:


pmom


# ## Load the calculated phase function
# After subjecting the delta-m scaled legendre polynomials expansion to phase program within libradtran:
# 
# $RES/libradtran/libRadtran-2.0.1/bin/phase -f -c -d ./pph.out > pphase.out

# In[50]:


phase = np.genfromtxt(fp+'ice/pphase.out')


# In[51]:


phase


# In[52]:


phase.shape


# # Now plot the phase functions and legendre polynomials

# In[46]:


nn = m['ntheta'][25,200]


# In[53]:


plt.figure()
plt.plot(m['theta'][25,200,:],m['phase'][25,200,:],label='Original Mie phase function')
plt.plot(phase[:,0],phase[:,1],'-',label='After reconstitution from scaled legendre polynomial')
plt.gca().set_yscale('log')
plt.ylabel('Phase function')
plt.xlabel('Scattering Angle')
plt.xlim(0,180)
plt.legend(frameon=False)


# # Create a loop routine to build the legendre phase functions 

# In[54]:


m['ntheta'].shape


# In[55]:


m['ref'].shape


# In[56]:


m['wvl'].shape


# In[59]:


import subprocess


# In[92]:


fdelta = fp+'ice/deltam3'


# In[62]:


m['pmom'].shape


# In[63]:


m['nmom'].shape


# In[98]:


m['pmom_delta'] = np.zeros((len(m['ref']),len(m['wvl']),301))
m['nmom_delta'] = np.zeros_like(m['nmom'])


# In[101]:


for i,r in enumerate(m['ref']):
    for j,w in enumerate(m['wvl']):
        
        fname_in = fp+'input/phase/mie_r{:02.1f}_w{:04.0f}.in'.format(r,w*1000.0)
        fname_out = fp+'output/phase/mie_r{:02.1f}_w{:04.0f}.out'.format(r,w*1000.0)
        
        with open(fname_in,'w') as f:
            f.write('{:4.0f}\n'.format(m['ntheta'][i,j]))
            for ii in xrange(m['ntheta'][i,j],0,-1):
                f.write('{} {}\n'.format(m['theta'][i,j,ii-1],m['phase'][i,j,ii-1]))
        
        nx = int(2.0*np.pi*r/(w*1000.0)*1000.0)
        if nx<26: nx =26
        if nx>300: nx=300
        
        m['nmom_delta'][i,j] = nx
        
        print fname_in, nx
        subprocess.call([fdelta,fname_in,fname_out,'{}'.format(nx)])
        pmom = np.genfromtxt(fname_out)
        m['pmom_delta'][i,j,:nx+1] = pmom
        
    sio.savemat(fs+'rtm_dat/mie_hi_delta.mat',m)


# In[102]:


sio.savemat(fs+'rtm_dat/mie_hi_delta.mat',m)


# # Now test out the new file

# In[119]:


mie = sio.loadmat(fs+'rtm_dat/mie_hi_delta.mat')


# In[120]:


mie['pmom_delta'].shape


# In[105]:


mie['nmom_delta'].shape


# In[106]:


mie['nmom_delta'][10,200]


# In[109]:


mie['pmom_delta'][10,200,0:mie['nmom_delta'][10,200]]


# In[111]:


fp_rtm = fs+'rtm_dat/'


# In[112]:


fp_rtm


# In[113]:


mie = sio.netcdf_file(fp_rtm+'wc_allpmom.sol.mie.cdf','r')
mie_long = sio.netcdf_file(fp_rtm+'wc.sol.long.mie.cdf','r')
try:
    rho = np.swapaxes(mie.variables['rho'].data,0,1)
except ValueError:
    rho = mie.variables['rho'].data
mie_short = {'wvl':mie.variables['wavelen'].data,
             'ref':mie.variables['reff'].data,
             'ntheta':np.swapaxes(mie.variables['ntheta'].data[:,:,0],0,1),
             'rho':rho,
             'nmom':np.swapaxes(mie.variables['nmom'].data,0,1),
             'ssa':np.swapaxes(mie.variables['ssa'].data,0,1),
             'ext':np.swapaxes(mie.variables['ext'].data,0,1),
             'nim':mie.variables['refim'].data,
             'nre':mie.variables['refre'].data,
             'pmom':np.swapaxes(mie.variables['pmom'].data[:,:,0,:],0,1),
             'phase':np.swapaxes(mie.variables['phase'].data[:,:,0,:],0,1),
             'theta': np.swapaxes(mie.variables['theta'].data[:,:,0,:],0,1)}
pmom = {'wvl':np.append(mie_short['wvl'],mie_long.variables['wavelen'].data[7:]),
        'ref':mie_short['ref'],
        'ntheta':np.concatenate((mie_short['ntheta'],np.swapaxes(mie_long.variables['ntheta'].data[7:,:-5,0],0,1)),axis=1),
        'rho':mie_short['rho'],
        'nmom':np.concatenate((mie_short['nmom'],np.swapaxes(mie_long.variables['nmom'].data[7:,:-5,0],0,1)),axis=1),
        'ssa':np.concatenate((mie_short['ssa'],np.swapaxes(mie_long.variables['ssa'].data[7:,:-5],0,1)),axis=1),
        'ext':np.concatenate((mie_short['ext'],np.swapaxes(mie_long.variables['ext'].data[7:,:-5],0,1)),axis=1),
        'nim':np.append(mie_short['nim'],mie_long.variables['refim'].data[7:]),
        'nre':np.append(mie_short['nre'],mie_long.variables['refre'].data[7:]),
        'pmom':np.concatenate((mie_short['pmom'],np.concatenate((np.swapaxes(mie_long.variables['pmom'].data[7:,:-5,0,:],0,1),
                                                                 np.zeros((25,72,2500))),axis=2)),axis=1),
        'phase':np.concatenate((mie_short['phase'],np.swapaxes(mie_long.variables['phase'].data[7:,:-5,0,:],0,1)),axis=1),
        'theta':np.concatenate((mie_short['theta'],np.swapaxes(mie_long.variables['theta'].data[7:,:-5,0,:],0,1)),axis=1)}
pmom['file_name'] = [fp_rtm+'wc_allpmom.sol.mie.cdf',fp_rtm+'wc.sol.long.mie.cdf']


# In[114]:


pmom['nmom'].shape


# In[116]:


pmom['ref'].shape


# In[121]:


pmom['wvl']


# In[122]:


mie['wvl']


# In[123]:


mie['pmom_delta'].shape


# In[126]:


mie_short['wvl']


# In[127]:


mie.keys()


# In[128]:


rho.shape


# In[129]:


rho


# In[132]:


mie['rho'].max()


# In[136]:


mie_long.variables['wavelen'].data[:]


# In[137]:


np.append(mie['wvl'],mie_long.variables['wavelen'].data[:])


# In[141]:


np.concatenate((mie['ntheta'],np.swapaxes(mie_long.variables['ntheta'].data[:,:,0],0,1)),axis=1).shape


# In[148]:


np.append(mie['nre'],mie_long.variables['refre'].data[:]).shape


# In[156]:


mie['pmom_delta'].shape


# In[149]:


np.concatenate((mie['pmom'],np.concatenate((np.swapaxes(mie_long.variables['pmom'].data[:,:,0,:],0,1),
                                                                             np.zeros((25,72,2500))),axis=2)),axis=1).shape


# In[150]:


np.concatenate((mie_short['phase'],np.swapaxes(mie_long.variables['phase'].data[:,:,0,:],0,1)),axis=1)


# In[151]:


np.concatenate((mie_short['theta'],np.swapaxes(mie_long.variables['theta'].data[:,:,0,:],0,1)),axis=1)


# In[158]:


nmom = np.concatenate((mie['nmom_delta'],np.swapaxes(mie_long.variables['nmom'].data[:,:,0],0,1)),axis=1).max()


# In[161]:


mie['pmom_delta'].shape


# In[166]:


np.swapaxes(mie_long.variables['pmom'].data[:,:,0,:652],0,1).shape


# In[165]:


mie_long.variables['nmom'].data[:,:,:].max()


# In[175]:


np.concatenate((np.concatenate((mie['pmom_delta'],np.zeros((30,754,351))),axis=2),np.swapaxes(mie_long.variables['pmom'].data[:,:,0,:652],0,1)),axis=1).shape


# In[176]:


np.concatenate((mie['phase'],np.swapaxes(mie_long.variables['phase'].data[:,:,0,:],0,1)),axis=1)


# In[177]:


mie['phase'].shape


# In[180]:


np.swapaxes(mie_long.variables['phase'].data[:,:,0,:398],0,1).shape


# In[179]:


mie_long.variables['ntheta'].data[:,:,0].max()


# In[183]:


np.concatenate((mie['phase'],np.swapaxes(mie_long.variables['phase'].data[:,:,0,:398],0,1)),axis=1).shape


# In[184]:


np.concatenate((mie['theta'],np.swapaxes(mie_long.variables['theta'].data[:,:,0,:398],0,1)),axis=1).shape


# In[185]:


np.concatenate((mie['ntheta'],np.swapaxes(mie_long.variables['ntheta'].data[:,:,0],0,1)),axis=1).shape


# In[186]:


pmom = {'wvl':np.append(mie['wvl'],mie_long.variables['wavelen'].data[:]),
                        'ref':mie['ref'],
                        'ntheta':np.concatenate((mie['ntheta'],np.swapaxes(mie_long.variables['ntheta'].data[:,:,0],0,1)),axis=1),
                        'rho':mie['rho'],
                        'nmom':np.concatenate((mie['nmom_delta'],np.swapaxes(mie_long.variables['nmom'].data[:,:,0],0,1)),axis=1),
                        'ssa':np.concatenate((mie['ssa'],np.swapaxes(mie_long.variables['ssa'].data[:,:],0,1)),axis=1),
                        'ext':np.concatenate((mie['ext'],np.swapaxes(mie_long.variables['ext'].data[:,:],0,1)),axis=1),
                        'nim':np.append(mie['nim'],mie_long.variables['refim'].data[:]),
                        'nre':np.append(mie['nre'],mie_long.variables['refre'].data[:]),
                        'pmom':np.concatenate((np.concatenate((mie['pmom_delta'],np.zeros((30,754,351))),axis=2),
                                               np.swapaxes(mie_long.variables['pmom'].data[:,:,0,:652],0,1)),axis=1),
                        'phase':np.concatenate((mie['phase'],np.swapaxes(mie_long.variables['phase'].data[:,:,0,:398],0,1)),axis=1),
                        'theta':np.concatenate((mie['theta'],np.swapaxes(mie_long.variables['theta'].data[:,:,0,:398],0,1)),axis=1)}


# In[188]:


pmom['nmom'].shape

