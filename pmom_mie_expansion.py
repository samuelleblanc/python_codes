#!/usr/bin/env python
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
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2020-11-16
#               - updated for a new set of mie calculations, with nthetamax = 1500, and reff 0.5-20 at 0.5 um resolution

# # Load the required modules and set up the paths

# In[1]:


import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
get_ipython().magic(u'matplotlib notebook')
from path_utils import getpath


# In[66]:


get_ipython().magic(u'matplotlib notebook')


# In[2]:


fp = getpath('libradtran')
fp


# In[6]:


fs = getpath('4STAR')
fs


# In[9]:


frtm = getpath('rtm_home')


# # Load the files

# ## Load the 2017 runs

# ### Load the original mie_hi.out file with all the pahse functions and non-scaled legendre polynomials

# In[345]:


m = sio.idl.readsav(fs+'rtm_dat/mie_hi.out')


# In[346]:


m.keys()


# In[347]:


with open(fs+'rtm_dat/phase_an.txt','w') as f:
    f.write('{:4.0f}\n'.format(m['ntheta'][25,200]))
    for i in xrange(m['ntheta'][25,200],0,-1):
        f.write('{} {}\n'.format(m['theta'][25,200,i-1],m['phase'][25,200,i-1]))


# ### Load the calculated legendre polynomials with delta-m scaling factor
# After subjecting the phase_an.txt to deltam program from Y. Hu's program compiled to deltam:
# 
# ./deltam3 ../../4STAR/rtm_dat/phase_an.txt pph.out 300
# 
# #sed -i '1d' pph.out

# In[48]:


pmom = np.genfromtxt(fp+'ice/pph.out')


# In[49]:


pmom


# ### Load the calculated phase function
# After subjecting the delta-m scaled legendre polynomials expansion to phase program within libradtran:
# 
# $RES/libradtran/libRadtran-2.0.1/bin/phase -f -c -d ./pph.out > pphase.out

# In[50]:


phase = np.genfromtxt(fp+'ice/pphase.out')


# In[51]:


phase


# In[52]:


phase.shape


# ## Load the 2020 runs

# ### Load the netcdf file
# Output from libradtran's mie calculations.  
# 
# Input file:  
# 
# `mie_program MIEV0
# refrac water
# r_eff 0.5 20 0.5
# distribution gamma 7
# wavelength 350 1700
# wavelength_step 1
# nstokes 1
# nthetamax 1500
# nmom 0
# output_user netcdf
# verbose
# basename mie_wc_run2020112_SL`
# 
# 
# And Other input file:  
#   
# `mie_program MIEV0
# refrac water
# r_eff 15 20 0.5
# distribution gamma 7
# wavelength 350 550
# wavelength_step 1
# nstokes 1
# nthetamax 2500
# nmom 0
# output_user netcdf
# verbose
# basename mie_wc_run2020115_largeonly_SL`

# In[10]:


import load_utils as lu


# In[612]:


vv = 'run2020112_SL'


# In[15]:


mi,mi_dict = lu.load_netcdf(frtm+'mie_wc_{}mie.cdf'.format(vv),everything=True)


# In[17]:


mi['ntheta'].shape


# In[20]:


mi['wavelen'][150]


# In[22]:


mi['phase'].shape


# In[23]:


mi['reff'].shape


# In[28]:


mi['theta'][0,0,0,:]


# In[25]:


mi['phase'][0,0,0,:]


# In[38]:


for u,g in enumerate(mi['pmom'][0,0,0,:].compressed()):
    print u,g


# ### Build the inputs to libradtran's pmom

# In[24]:


from tqdm import tqdm_notebook as tqdm


# In[37]:


'r{:04.1f}'.format(7.5)


# In[45]:


theta = mi['theta'][:,:,0,:0:-1]


# In[49]:


for i,p in enumerate(mi['phase'][0,0,0,:0:-1].compressed()):
    print mi['theta'][0,0,0,:0:-1].compressed()[i],i,p


# In[455]:


pbar = tqdm(total=len(mi['wavelen'])*len(mi['reff']))
for nw, w in enumerate(mi['wavelen']):
    for nr,r in enumerate(mi['reff']):
        mitheta = mi['theta'][nw,nr,0,:0:-1].compressed()
        with open(frtm+'pmom_input/phase_{vv}_w{w:01.3f}_r{r:04.1f}.dat'.format(vv=vv,w=w,r=r),'w') as f:
            #f.write('{:4.0f}\n'.format(mi['ntheta'][nw,nr,0]))
            #for i,p in enumerate(mi['phase'][nw,nr,0,:0:-1].compressed()):
            [f.write('{} {}\n'.format(mitheta[i],np.abs(p)))             for i,p in enumerate(mi['phase'][nw,nr,0,:0:-1].compressed())]
        pbar.update()


# ### Now direct run the pmom and save to arrays

# In[54]:


fp_lib = getpath('libradtran')


# In[53]:


import subprocess


# #### Test out a few pmom decomposition and number

# In[464]:


w = 0.354
r = 20.0
pmom = []
for nleg in [25,100,200,300,500,1000,1500,2000,4000,10000]:
    print nleg
    pmoms_tmp = subprocess.check_output([fp_lib+'bin/pmom','-l','{}'.format(nleg),'-n','-c',
                                 frtm+'pmom_input/phase_{vv}_w{w:01.3f}_r{r:04.1f}.dat'.format(vv=vv,w=w,r=r)])
    pmom.append(np.array([float(k) for k in  pmoms_tmp.split()]))                               


# In[465]:


len(mi['phase'][nw,nr,0,:].compressed())


# In[439]:


len(pmom)


# In[466]:


plt.figure()
for pm in pmom:
    plt.plot(pm,'.-',label='N={}'.format(len(pm)))
    print len(pm),np.nanmean(pm[:int(len(pm)*-0.05)])
    plt.xscale('log')
    plt.legend()


# In[467]:


int(len(pm)*-0.1)


# #### Rebuild the phase for testing

# In[468]:


phase = []
for pm in pmom[0:]:
    with open(frtm+'pmom_input/pmom_temp.dat','w') as f:
        [f.write('{}\n'.format(mom/(2*im+1))) for im,mom in enumerate(pm)]
    phase_tmp = subprocess.check_output([fp_lib+'bin/phase','-c','-d','-s','0.5','-n','-f',
                                 frtm+'pmom_input/pmom_temp.dat'])
    phase.append(np.array(phase_tmp.split()))
phase = [ph.reshape(len(phase[0])/3,3).astype(np.float) for ph in phase]


# In[469]:


len(phase)


# In[470]:


nw = 4
nr = 39
mi['reff'][nr], mi['wavelen'][nw]


# In[471]:


plt.figure()
plt.plot(mi['theta'][nw,nr,0,:].compressed(),mi['phase'][nw,nr,0,:].compressed(),'-',label='OG')
for i,ph in enumerate(phase):
    plt.plot(ph[:,0],ph[:,1],'.',label='N={}'.format(len(pmom[i])))
    print len(ph[ph[:,1]>0,1]),len(ph[:,1])
plt.legend()
plt.yscale('log')


# ### Build the functions to quantify the pmoms

# In[ ]:


frtm+'pmom_input/phase_{vv}_w{w:01.3f}_r{r:04.1f}.dat'.format(vv=vv,w=w,r=r)


# In[473]:


def calc_pmom(file_in,nleg):
    'Calculate the pmom, as a wrapper to libradtran pmom, file_in is the phase function file'
    pmoms_tmp = subprocess.check_output([fp_lib+'bin/pmom','-l','{}'.format(nleg),'-n','-c',file_in])
    return np.array([float(k) for k in  pmoms_tmp.split()]) 


# In[474]:


fp_rtm = getpath('rtm')


# In[539]:


def calc_phase(pmom):
    'Calculate the phase function, uses coefficients as input (same needed for libradtran)'
    with open(fp_rtm+'pmom_input/pmom_temp.dat','w') as f:
        [f.write('{}\n'.format(mom/(2*im+1))) for im,mom in enumerate(pmom)]
    phase_tmp = subprocess.check_output([fp_lib+'bin/phase','-c','-d','-s','0.5','-n','-f',
                                         fp_rtm+'pmom_input/pmom_temp.dat'])
    ph = np.array(phase_tmp.split())
    return ph.reshape(len(ph)/3,3).astype(np.float)


# In[493]:


pha = np.array([[-1,-1],[-1,-1]])
i = 0
max_iter = 20
while any(pha[:,1]<0):
    i = i+1
    nleg = 2*nleg
    if i>max_iter: 
        print '** did not converge for'
        break


# In[561]:


mi['reff'][29]


# In[560]:


mi['wavelen'][50]


# In[559]:


max_nmom = 9000
pbar = tqdm(total=len(mi['wavelen'])*len(mi['reff']))
pmoms = np.ma.masked_array(np.zeros((len(mi['wavelen']),len(mi['reff']),max_nmom)),mask=True)
for nw, w in enumerate(mi['wavelen']):
    for nr,r in enumerate(mi['reff']):
        nleg = mi['ntheta'][nw,nr,0]
        if (nleg<1500) & (mi['nmom'][nw,nr,0] <1000):
            pbar.update()
            continue
        pha = np.array([[-1,-1],[-1,-1]])
        while any(pha[:,1]<0):
            pmom_tmp = calc_pmom(frtm+'pmom_input/phase_{vv}_w{w:01.3f}_r{r:04.1f}.dat'.format(vv=vv,w=w,r=r),nleg)
            pha = calc_phase(pmom_tmp)
            nleg = int(2*nleg)
            if nleg>max_nmom: 
                print '** did not converge for nw:{}, nr:{}, nleg:{}'.format(nw,nr,nleg)
                break
        pmoms[nw,nr,:len(pmom_tmp)] = pmom_tmp 
        pbar.update()


# In[565]:


for nw, w in enumerate(mi['wavelen']):
    for nr,r in enumerate(mi['reff']):
        nleg = mi['ntheta'][nw,nr,0]
        if (nleg<1500) & (mi['nmom'][nw,nr,0] <1000):
            pmoms[nw,nr,:mi['nmom'][nw,nr,0]+1] = mi['pmom'][nw,nr,0,:].compressed()


# ### Load from the Hu 
# /data/sam/libradtran/ice  
# 
# ./hu_v3 /home/sam/rtm/pmom_input/phase_run2020112_SL_w0.354_r01.0.dat ./a.dat 2902

#  As test case

# In[281]:


with open('/data/sam/libradtran/ice/a.dat') as f:
    hus = f.read()
hu = np.array(hus.split()).astype(np.float)


# In[282]:


hu


# In[283]:


plt.figure()
plt.plot(hu)


# In[284]:


with open(frtm+'pmom_input/pmom_temp.dat','w') as f:
    [f.write('{}\n'.format(mom)) for mom in hu]
phase_tmp = subprocess.check_output([fp_lib+'bin/phase','-c','-d','-n','-f',
                             frtm+'pmom_input/pmom_temp.dat'])
phasehu = np.array(phase_tmp.split())
phasehu = phasehu.reshape(len(phasehu)/3,3).astype(np.float)


# In[285]:


mi['reff'][1]


# In[342]:


plt.figure()
nw = 4
nr = 1
plt.plot(mi['theta'][nw,nr,0,:].compressed(),mi['phase'][nw,nr,0,:].compressed(),'-',label='OG')
plt.plot(phasehu[:,0],phasehu[:,1],'-')
plt.yscale('log')


# ! Match!

# ### Calculate the pmom using the Hu delta-m

# In[297]:


fp_hu = '/data/sam/libradtran/ice'


# In[312]:


w = 0.354
r = 7.5
pmom = []
for nleg in [25,100,200,500,1000,2500]:
    print nleg
    pmoms_tmp = subprocess.check_output([fp_hu+'/hu_v3',
                                         frtm+'pmom_input/phase_{vv}_w{w:01.3f}_r{r:04.1f}.dat'.format(vv=vv,w=w,r=r),
                                         './dummy.dat','{}'.format(nleg)])
    pmom.append(np.array([float(k) for k in  pmoms_tmp.split()]))      


# In[343]:


frtm+'pmom_input/phase_{vv}_w{w:01.3f}_r{r:04.1f}.dat'.format(vv=vv,w=w,r=r)


# In[316]:


nleg = 3000
pmoms_tmp = subprocess.check_output([fp_hu+'/hu_v3',
                                         frtm+'pmom_input/phase_{vv}_w{w:01.3f}_r{r:04.1f}.dat'.format(vv=vv,w=w,r=r),
                                         './dummy.dat','{}'.format(nleg),'0'],stderr=subprocess.STDOUT)
print pmoms_tmp


# In[303]:


pmom


# In[299]:


phase = []
for pm in pmom:
    with open(frtm+'pmom_input/pmom_temp.dat','w') as f:
        [f.write('{}\n'.format(mom)) for mom in pm]
    phase_tmp = subprocess.check_output([fp_lib+'bin/phase','-c','-d','-n','-f',
                                 frtm+'pmom_input/pmom_temp.dat'])
    phase.append(np.array(phase_tmp.split()))
phase = [ph.reshape(len(phase[0])/3,3).astype(np.float) for ph in phase]


# ### Pull out the pmom from mie

# In[319]:


mi['pmom'][nw,nr,0,:].compressed()


# In[362]:


nw = 1250
nr = 1


# In[462]:


phase = []
with open(frtm+'pmom_input/pmom_temp.dat','w') as f:
    [f.write('{}\n'.format(mom/(2*im+1))) for im,mom in enumerate(mi['pmom'][nw,nr,0,:].compressed())]
phase_tmp = subprocess.check_output([fp_lib+'bin/phase','-c','-d','-s','0.1',
                             frtm+'pmom_input/pmom_temp.dat'])
phase.append(np.array(phase_tmp.split()))
phase = [ph.reshape(len(ph)/3,3).astype(np.float) for ph in phase]


# In[463]:


plt.figure()
nw = 4
nr = 39
plt.plot(mi['theta'][nw,nr,0,:].compressed(),mi['phase'][nw,nr,0,:].compressed(),'-',label='OG')
plt.plot(phase[0][:,0],phase[0][:,1],'.')
plt.yscale('log')


# ## Load longer calcs

# In[567]:


vv = 'run2020115_largeonly_SL'


# In[568]:


mi2,mi2_dict = lu.load_netcdf(frtm+'mie_wc_{}mie.cdf'.format(vv),everything=True)


# In[569]:


mi2['nmom'].shape


# In[570]:


mi2['wavelen']


# In[571]:


mi2['reff']


# In[578]:


mi2['ntheta'][50,:,0]


# In[579]:


mi2['nmom'][50,:,0]


# In[581]:


pmoms.shape


# In[586]:


pmoms[0,39,:]


# In[589]:


for i in xrange(40):
    print i,len(pmoms[0,i,:].compressed()),len(mi['phase'][0,i,:].compressed())


# In[603]:


pha2 = calc_phase(pmoms[0,31,:])


# In[606]:


plt.figure()
plt.plot(pha2[:,0],pha2[:,1],'.',label='pha2')
plt.plot(mi['theta'][0,31,0,:],mi['phase'][0,31,0,:],'.',label='mi')
plt.plot(mi2['theta'][0,2,0,:],mi2['phase'][0,2,0,:],'+',label='mi2')

plt.legend()
plt.yscale('log')


# In[601]:


mi['reff'][31]


# In[602]:


mi2['reff'][2]


# In[607]:


mi['pmom'] = pmoms


# In[608]:


mi.keys()


# In[609]:


mi['pmom'].shape


# In[610]:


frtm


# ## Check and save to file

# In[637]:


for nw,w in enumerate(mi['wavelen']):
    for nr,r in enumerate(mi['reff']):
        mi['nmom'][nw,nr,0] = len(mi['pmom'][nw,nr,:].compressed())


# In[688]:


mii = {}
for k in mi.keys():
    mii[k] = mi[k].data
    print k, mi[k].shape


# In[692]:


sio.savemat(frtm+'dat/mie_{}.mat'.format(vv),mii)


# In[667]:


frtm+'dat/mie_{}.mat'.format(vv)


# In[638]:


for k in mi.keys():
    print k, mi[k].shape


# In[693]:


reload(ru)


# In[694]:


import Run_libradtran as ru
pmom = ru.make_pmom_inputs(fp_rtm=frtm+'dat/',new=True)


# In[695]:


for k in pmom.keys():
    try:
        print k, pmom[k].shape
    except:
        pass


# In[696]:


pmom['pmom'][20,140,:20]


# In[697]:


mi['pmom'].shape


# In[698]:


mi['pmom'][140,20,:20]


# In[699]:


ppm = mi['pmom']


# In[700]:


ppo = np.swapaxes(ppm,0,1)


# In[701]:


ppo[20,140,:20]


# # Now plot the phase functions and legendre polynomials

# In[46]:


nn = m['ntheta'][25,200]


# In[344]:


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


# In[ ]:





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


# In[ ]:




