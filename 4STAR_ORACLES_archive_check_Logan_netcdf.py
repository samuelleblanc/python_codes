#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Describe the details ...
# 
# Input:
# 
#     arguments
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
#     - hdf5storage
#     - scipy
# 
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-12-01
#     Modified:
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
from datetime import datetime 


# In[66]:


name = 'ORACLES_2017_4wvl_nc_daily'
vv = 'v1'
fp = getpath('sunsat_oracles2016')


# In[128]:


fpo = fp+'data_processed/starskies/skyscans_Logan_20221201_nc_2017/'


# In[7]:


fpo = fp+'data_processed/starskies/Logan_skyscans_20221209/'


# In[73]:


name = 'ORACLES_2016_5wvl_nc_daily'
fpo = fp+'data_processed/starskies/skyscans_Logan_20221209_nc_2016/'


# In[36]:


name8 = 'ORACLES_2018_4wvl_nc_daily'
fpo8 = fp+'data_processed/starsky/Logan_starskies_20221209/'


# # Load files

# In[68]:


flist = os.listdir(fpo)
flist.sort()
flist


# In[69]:


nc = []
nc_dict = []
for f in flist:
    nc_tmp,nc_dict_tmp = lu.load_netcdf(fpo+f,everything=True)
    nc.append(nc_tmp)
    nc_dict.append(nc_dict_tmp) 


# In[70]:


for i,n in enumerate(nc):
    if i==0: 
        nc_all = n   
        nc_all[b'time'] = nc_all[b'time']+nc_all[b'base_time']
    else:
        for k in n:
            if not k in [b'wavelength',b'radius',b'PF_angle']:
                if k in [b'base_time']:
                    nc_all[k] = np.append(nc_all[k],n[k])
                elif k in [b'time']:
                    nc_all[k] = np.append(nc_all[k],n[k]+n[b'base_time'])
                else:
                    nc_all[k] = np.append(nc_all[k],n[k],axis=0)


# In[71]:


for k in nc_all:
    print(k,nc_all[k].shape,nc[1][k].shape)


# In[72]:


nc_all[b'datetime'] = np.array([datetime.fromtimestamp(t) for t in nc_all[b'time']]).astype(np.datetime64)


# ## load 2018 files

# In[37]:


flist8 = os.listdir(fpo8)
flist8.sort()
flist8


# In[38]:


nc8 = []
nc8_dict = []
for f in flist8:
    nc8_tmp,nc8_dict_tmp = lu.load_netcdf(fpo8+f,everything=True)
    nc8.append(nc8_tmp)
    nc8_dict.append(nc8_dict_tmp) 


# In[39]:


for i,n in enumerate(nc8):
    if i==0: 
        nc8_all = n   
        nc8_all[b'time'] = nc8_all[b'time']+nc8_all[b'base_time']
    else:
        for k in n:
            if not k in [b'wavelength',b'radius',b'PF_angle']:
                if k in [b'base_time']:
                    nc8_all[k] = np.append(nc8_all[k],n[k])
                elif k in [b'time']:
                    nc8_all[k] = np.append(nc8_all[k],n[k]+n[b'base_time'])
                else:
                    nc8_all[k] = np.append(nc8_all[k],n[k],axis=0)


# In[40]:


nc8_all[b'datetime'] = np.array([datetime.fromtimestamp(t) for t in nc8_all[b'time']]).astype(np.datetime64)


# # Plot out data

# In[74]:


for k in nc_all:
    if nc_all[k].shape == nc_all[b'time'].shape:
        plt.figure()
        plt.plot(nc_all[b'datetime'],nc_all[k],'.')
        try: 
            plt.ylabel(nc_dict_tmp[k].long_name + ' ['+nc_dict_tmp[k].units+']')
        except:
            pass
        plt.xlabel('time')
        plt.title(k.decode())        
        
        plt.savefig(fpo+name+k.decode()+'.png',transparent=True,dpi=500)


# In[75]:


for k in nc_all:
    if nc_all[k].shape == nc_all[b'SSA'].shape:
        plt.figure()
        plt.gca().set_prop_cycle(color=[plt.cm.viridis(k) for k in np.linspace(0, 1, len(nc_all[b'time'])+1)])
        plt.plot(nc_all[b'wavelength'],nc_all[k].T,'.-')
        plt.boxplot(nc_all[k],positions=nc_all[b'wavelength'],widths=20.0,showmeans=True,sym='+',meanprops={'markerfacecolor':'tab:orange','marker':'D'})
        plt.xlabel('Wavelength [nm]')
        plt.ylabel(nc_dict_tmp[k].long_name + ' ['+nc_dict_tmp[k].units+']')
        plt.title(k.decode())
        scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
        scalarmap.set_array(np.linspace(0, 1, len(nc_all[b'time'])+1))
        cb = plt.colorbar(scalarmap)
        indexes = [np.argmin(abs(np.linspace(0, 1, len(nc_all[b'time']))-ti)) for ti in cb.get_ticks()]
        cb.ax.set_yticklabels(nc_all[b'datetime'][indexes])
        
        plt.savefig(fpo+name+k.decode()+'.png',transparent=True,dpi=500)


# In[76]:


for k in nc_all:
    if nc_all[k].shape == nc_all[b'PF_coarse'].shape:
        fig,ax = plt.subplots(4,1)
        for i in range(4):
            ax[i].set_prop_cycle(color=[plt.cm.cividis(k) for k in np.linspace(0, 1, len(nc_all[b'time'])+1)])
            ax[i].plot(nc_all[b'PF_angle'],nc_all[k][:,i,:].T,'.-')
            ax[i].set_ylabel('{} nm'.format(nc_all[b'wavelength'][i]))
            ax[i].set_yscale('log')

        ax[0].set_title(k.decode()+': '+nc_dict_tmp[k].long_name + ' ['+nc_dict_tmp[k].units+']')
        ax[-1].set_xlabel(nc_dict_tmp[b'PF_angle'].long_name + ' ['+nc_dict_tmp[b'PF_angle'].units+']')
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
        scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.cividis)
        scalarmap.set_array(np.linspace(0, 1, len(nc_all[b'time'])+1))
        cb = plt.colorbar(scalarmap,cax=cbar_ax)
        indexes = [np.argmin(abs(np.linspace(0, 1, len(nc_all[b'time']))-ti)) for ti in cb.get_ticks()]
        cb.ax.set_yticklabels(nc_all[b'datetime'][indexes])
        
        
        plt.savefig(fpo+name+k.decode()+'.png',transparent=True,dpi=500)
        


# In[77]:


for k in nc_all:
    if nc_all[k].shape == nc_all[b'normalized_sky_radiance'].shape:
        fig,ax = plt.subplots(4,1)
        for i in range(4):
            ax[i].set_prop_cycle(color=[plt.cm.CMRmap(k) for k in np.linspace(0, 1, len(nc_all[b'time'])+1)])
            ax[i].plot(nc_all[b'sca_angle'].T,nc_all[k][:,i,:].T,'.-',alpha=0.2)
            ax[i].set_ylabel('{} nm'.format(nc_all[b'wavelength'][i]))
            ax[i].set_yscale('log')

        ax[0].set_title(k.decode()+': '+nc_dict_tmp[k].long_name + ' ['+nc_dict_tmp[k].units+']')
        ax[-1].set_xlabel(nc_dict_tmp[b'sca_angle'].long_name + ' ['+nc_dict_tmp[b'sca_angle'].units+']')
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
        scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.CMRmap)
        scalarmap.set_array(np.linspace(0, 1, len(nc_all[b'time'])+1))
        cb = plt.colorbar(scalarmap,cax=cbar_ax)
        indexes = [np.argmin(abs(np.linspace(0, 1, len(nc_all[b'time']))-ti)) for ti in cb.get_ticks()]
        cb.ax.set_yticklabels(nc_all[b'datetime'][indexes])
        
        plt.savefig(fpo+name+k.decode()+'.png',transparent=True,dpi=500)


# In[78]:


plt.figure()
plt.gca().set_prop_cycle(color=[plt.cm.cubehelix(k) for k in np.linspace(0, 1, len(nc_all[b'time'])+1)])
plt.plot(nc_all[b'radius'],nc_all[b'psd'].T)
plt.xscale('log')

plt.xlabel(nc_dict_tmp[b'radius'].long_name + ' ['+nc_dict_tmp[b'radius'].units+']')
plt.ylabel(nc_dict_tmp[b'psd'].long_name + ' ['+nc_dict_tmp[b'psd'].units+']')
plt.title(b'psd'.decode())

scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.cubehelix)
scalarmap.set_array(np.linspace(0, 1, len(nc_all[b'time'])+1))
cb = plt.colorbar(scalarmap)
indexes = [np.argmin(abs(np.linspace(0, 1, len(nc_all[b'time']))-ti)) for ti in cb.get_ticks()]
cb.ax.set_yticklabels(nc_all[b'datetime'][indexes])

plt.savefig(fpo+name+'psd'+'.png',transparent=True,dpi=500)


# In[79]:


ae_fine_mode = [nc_all[b'radius'][np.argmax(p)] for p in nc_all[b'psd'][:,0:9]]
ae_coarse_mode = [nc_all[b'radius'][9:][np.argmax(p)] for p in nc_all[b'psd'][:,9:]]


# In[80]:


fig = plt.figure()
#
bp = plt.boxplot([ae_fine_mode,ae_coarse_mode],positions=[0.15,0.15],vert=False,showmeans=True,widths=0.02)
plt.plot(nc_all[b'radius'],nc_all[b'psd'].T)
plt.xscale('log')

[plt.axvline(m.get_data()[0],color='g',ls=':') for m in bp['means']]
[plt.text(m.get_data()[0]+0.01,0.165,'{:1.3f}'.format(m.get_data()[0][0])) for m in bp['means']]

plt.xlabel('Radius [$\mu$m]')
plt.ylabel('Volume density [dV/dlnR]')
plt.ylim(0,0.18)
plt.yticks([0.0,0.05,0.1,0.15])
plt.gca().set_yticklabels([0.0,0.05,0.1,0.15])

plt.title('ORACLES 2017 4STAR retrievals')

plt.savefig(fpo+name+'psd_avgeraged'+'.png',transparent=True,dpi=500)


# In[81]:


nc_dict_tmp.keys()


# In[82]:


ncf = sio.netcdf_file(fpo+f)


# In[83]:


for k in ncf.__dict__:
    try:
        print(k+':   ',ncf.__dict__[k].decode())
    except: 
        pass


# ## For 2018

# In[41]:


for k in nc8_all:
    if nc8_all[k].shape == nc8_all[b'time'].shape:
        plt.figure()
        plt.plot(nc8_all[b'datetime'],nc8_all[k],'.')
        try: 
            plt.ylabel(nc8_dict_tmp[k].long_name + ' ['+nc8_dict_tmp[k].units+']')
        except:
            pass
        plt.xlabel('time')
        plt.title(k.decode())        
        
        plt.savefig(fpo8+name8+k.decode()+'.png',transparent=True,dpi=500)


# In[50]:


for k in nc8_all:
    if nc8_all[k].shape == nc8_all[b'SSA'].shape:
        plt.figure()
        plt.gca().set_prop_cycle(color=[plt.cm.viridis(k) for k in np.linspace(0, 1, len(nc8_all[b'time'])+1)])
        plt.plot(nc8_all[b'wavelength'],nc8_all[k].T,'.-')
        plt.boxplot(nc8_all[k],positions=nc8_all[b'wavelength'],widths=20.0,showmeans=True,sym='+',meanprops={'markerfacecolor':'tab:orange','marker':'D'})
        plt.xlabel('Wavelength [nm]')
        plt.ylabel(nc8_dict_tmp[k].long_name + ' ['+nc8_dict_tmp[k].units+']')
        plt.title(k.decode())
        scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
        scalarmap.set_array(np.linspace(0, 1, len(nc8_all[b'time'])+1))
        cb = plt.colorbar(scalarmap)
        indexes = [np.argmin(abs(np.linspace(0, 1, len(nc8_all[b'time']))-ti)) for ti in cb.get_ticks()]
        cb.ax.set_yticklabels(nc8_all[b'datetime'][indexes])
        
        plt.savefig(fpo8+name8+k.decode()+'.png',transparent=True,dpi=500)


# In[59]:


for k in nc8_all:
    if nc8_all[k].shape == nc8_all[b'PF_coarse'].shape:
        fig,ax = plt.subplots(4,1)
        for i in range(4):
            ax[i].set_prop_cycle(color=[plt.cm.cividis(k) for k in np.linspace(0, 1, len(nc8_all[b'time'])+1)])
            ax[i].plot(nc8_all[b'PF_angle'],nc8_all[k][:,i,:].T,'.-')
            ax[i].set_ylabel('{} nm'.format(nc8_all[b'wavelength'][i]))
            ax[i].set_yscale('log')

        ax[0].set_title(k.decode()+': '+nc8_dict_tmp[k].long_name + ' ['+nc8_dict_tmp[k].units+']')
        ax[-1].set_xlabel(nc8_dict_tmp[b'PF_angle'].long_name + ' ['+nc8_dict_tmp[b'PF_angle'].units+']')
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
        scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.cividis)
        scalarmap.set_array(np.linspace(0, 1, len(nc8_all[b'time'])+1))
        cb = plt.colorbar(scalarmap,cax=cbar_ax)
        indexes = [np.argmin(abs(np.linspace(0, 1, len(nc8_all[b'time']))-ti)) for ti in cb.get_ticks()]
        cb.ax.set_yticklabels(nc8_all[b'datetime'][indexes])
        
        
        plt.savefig(fpo8+name8+k.decode()+'.png',transparent=True,dpi=500)
        


# In[52]:


for k in nc8_all:
    if nc8_all[k].shape == nc8_all[b'normalized_sky_radiance'].shape:
        fig,ax = plt.subplots(4,1)
        for i in range(4):
            ax[i].set_prop_cycle(color=[plt.cm.CMRmap(k) for k in np.linspace(0, 1, len(nc8_all[b'time'])+1)])
            ax[i].plot(nc8_all[b'sca_angle'].T,nc8_all[k][:,i,:].T,'.-',alpha=0.2)
            ax[i].set_ylabel('{} nm'.format(nc8_all[b'wavelength'][i]))
            ax[i].set_yscale('log')

        ax[0].set_title(k.decode()+': '+nc8_dict_tmp[k].long_name + ' ['+nc8_dict_tmp[k].units+']')
        ax[-1].set_xlabel(nc8_dict_tmp[b'sca_angle'].long_name + ' ['+nc8_dict_tmp[b'sca_angle'].units+']')
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
        scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.CMRmap)
        scalarmap.set_array(np.linspace(0, 1, len(nc8_all[b'time'])+1))
        cb = plt.colorbar(scalarmap,cax=cbar_ax)
        indexes = [np.argmin(abs(np.linspace(0, 1, len(nc8_all[b'time']))-ti)) for ti in cb.get_ticks()]
        cb.ax.set_yticklabels(nc8_all[b'datetime'][indexes])
        
        plt.savefig(fpo8+name8+k.decode()+'.png',transparent=True,dpi=500)


# In[53]:


plt.figure()
plt.gca().set_prop_cycle(color=[plt.cm.cubehelix(k) for k in np.linspace(0, 1, len(nc8_all[b'time'])+1)])
plt.plot(nc8_all[b'radius'],nc8_all[b'psd'].T)
plt.xscale('log')

plt.xlabel(nc8_dict_tmp[b'radius'].long_name + ' ['+nc8_dict_tmp[b'radius'].units+']')
plt.ylabel(nc8_dict_tmp[b'psd'].long_name + ' ['+nc8_dict_tmp[b'psd'].units+']')
plt.title(b'psd'.decode())

scalarmap = plt.cm.ScalarMappable(cmap=plt.cm.cubehelix)
scalarmap.set_array(np.linspace(0, 1, len(nc8_all[b'time'])+1))
cb = plt.colorbar(scalarmap)
indexes = [np.argmin(abs(np.linspace(0, 1, len(nc8_all[b'time']))-ti)) for ti in cb.get_ticks()]
cb.ax.set_yticklabels(nc8_all[b'datetime'][indexes])

plt.savefig(fpo8+name8+'psd'+'.png',transparent=True,dpi=500)


# In[60]:


ae8_fine_mode = [nc8_all[b'radius'][np.argmax(p)] for p in nc8_all[b'psd'][:,0:9]]
ae8_coarse_mode = [nc8_all[b'radius'][9:][np.argmax(p)] for p in nc8_all[b'psd'][:,9:]]


# In[61]:


fig = plt.figure()
#
bp = plt.boxplot([ae8_fine_mode,ae8_coarse_mode],positions=[0.15,0.15],vert=False,showmeans=True,widths=0.02)
plt.plot(nc8_all[b'radius'],nc8_all[b'psd'].T)
plt.xscale('log')

[plt.axvline(m.get_data()[0],color='g',ls=':') for m in bp['means']]
[plt.text(m.get_data()[0]+0.01,0.165,'{:1.3f}'.format(m.get_data()[0][0])) for m in bp['means']]

plt.xlabel('Radius [$\mu$m]')
plt.ylabel('Volume density [dV/dlnR]')
plt.ylim(0,0.18)
plt.yticks([0.0,0.05,0.1,0.15])
plt.gca().set_yticklabels([0.0,0.05,0.1,0.15])

plt.title('ORACLES 2018 4STAR retrievals')

plt.savefig(fpo8+name8+'psd_avgeraged'+'.png',transparent=True,dpi=500)


# In[62]:


nc8_dict_tmp.keys()


# In[64]:


ncf8 = sio.netcdf_file(fpo8+flist8[0])


# In[65]:


for k in ncf8.__dict__:
    try:
        print(k+':   ',ncf8.__dict__[k].decode())
    except: 
        pass


# In[ ]:




