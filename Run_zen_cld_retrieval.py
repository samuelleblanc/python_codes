
# coding: utf-8

# Name:  
# 
#     Run_zen_cld_retrieval
# 
# Purpose:  
# 
#     Run the zenith cloud property retrieval on the supplied 4STAR starzen.mat file
#     Applies the 15 parameters spectral cloud retrieval [LeBlanc et al., 2015, amt]
#     Run with python code versions of the lut saved to .mat and the starzen created with the matlab codes
#     Outputs pngs of the plots of the retrieved values
# 
# Calling Sequence:
# 
#     python Run_zen_cld_retrieval fp_starzen fp_zencld -o fp_zencld_plots -lut fp_lut_mat
#     where:
#         fp_starzen: the full path (with .mat extension) of the starzen file to read
#         fp_zencld: (optional) the full path of the output file (with .mat extension) of the saved retrieval file
#                    if omitted, uses the same directory as starzen.mat, but saves file as *_zencld.mat
#         fp_zencld_plots: (optional) the full path of the directory to save the plots from the retrieval
#                         if omitted uss the same directory as starzen.matà
#         fp_lut_mat: (optional) put in a special look up table file path 
# 
# Input:
# 
#     See Calling Sequence above.
# 
# Output:
#    
#     Save matlab file of retrieved zenith cloud properties
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - numpy
#     - hdf5storage : for saving and reading
#     - Run_Libradtran
#     - os
#     - sys
#     - matplotlib
#     - argparse
#     - Tkinter (for file open dialogs)
#     - Sp_parameter
#     - mpltools for color cycle
#     - load_modis for mat2py_time and toutc
#   
# Needed Files:
# 
#   - lut .mat file created by Read_save_NAAMES_lut
#   - starzen.mat file created by the matlab program : starzen.m
#     
# History:
# 
#     Written: Samuel LeBlanc, NASA Ames, 2015-10-26

# In[ ]:

import argparse


# ## Prepare command line argument parser

# In[39]:

long_description = """    Run the zenith cloud property retrieval on the supplied 4STAR starzen.mat file
    Applies the 15 parameters spectral cloud retrieval [LeBlanc et al., 2015, amt]
    Run with python code versions of the lut saved to .mat and the starzen created with the matlab codes
    Outputs pngs of the plots of the retrieved values"""


# In[69]:

parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('fp_starzen',nargs='?',help='the full path (with .mat extension) of the starzen file to read')
parser.add_argument('fp_zencld',nargs='?',help='''(optional) the full path of the output file (with .mat extension) of the saved retrieval file
                   if omitted, uses the same directory as starzen.mat, but saves file as *_zencld.mat''')
parser.add_argument('-o','--fp_zencld_plot',nargs='?',help="""(optional) the full path of the directory to save the plots from the retrieval
                        if omitted uss the same directory as starzen.mat""")
parser.add_argument('-lut','--fp_lut_mat',nargs='?',help='Put in special look up table file path')
parser.add_argument('-f',nargs='?',help='not used')
parser.add_argument('-n','--noplot',nargs='?',help='if set, will not output plots')


# In[70]:

in_ = vars(parser.parse_args())


# ## Import rest of needed modules

# In[ ]:

import sys
import os
import matplotlib.pyplot as plt
import Run_libradtran as RL
import hdf5storage as hs
import numpy as np
from load_modis import mat2py_time, toutc
import Sp_parameters as Sp
from mpltools import color
import scipy.io as sio
try:
    from mpl_toolkits.basemap import Basemap
    isbasemap = True
except:
    isbasemap = False


# ## Parse command line and get appropriate paths

# In[74]:

if in_['f']:
    print 'Using interactive version hard coded fp_starzen'
    if os.sys.platform == 'win32':
        fp = 'C:\\Users\\sleblan2\\Research\\SEAC4RS\\'
    elif os.sys.platform == 'linux2':
        fp = '/u/sleblan2/SEAC4RS/'
    else:
        raise Exception
    fp_starzen = os.path.abspath(fp+'../4STAR/SEAC4RS/20130913/20130913starzen_3.mat')
    print fp_starzen
else:
    if not in_.get('fp_starzen'):
        from Tkinter import Tk
        from tkFileDialog import asksaveasfilename
        Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
        filename = asksaveasfilename(defaultextension='*.mat',filetypes=[('mat file','*.mat'),('All files','*.*')]) # show an "Open" dialog box and return the path to the selected file
        fp_starzen = os.path.abspath(filename)
    else:
        fp_starzen = in_.get('fp_starzen')


# In[75]:

fp,fname = os.path.split(fp_starzen)
datestr = fname[0:8]


# In[76]:

if in_.get('fp_zencld'):
    fp_zencld = in_.get('fp_zencld')
else:
    fp_zencld = fp+'{datestr}_zencld.mat'.format(datestr=datestr)


# In[77]:

if in_.get('fp_zencld_plot'):
    fp_zencld_plot = in_.get('fp_zencld_plot')
else:
    fp_zencld = fp


# In[93]:

if in_.get('fp_lut_mat'):
    fp_lut_mat = in_.get('fp_lut_mat')
else:
    fp_lut_mat = fp+'/v2_NAAMES_lut.mat'


# In[ ]:

if in_.get('noplot'):
    noplot = True
else:
    noplot = False


# # Load starzen files and the look up table

# In[79]:

if not os.path.isfile(fp_starzen):
    raise IOError('file {} is not found'.format(fp_starzen))


# In[ ]:

print 'loading file {}'.format(fp_starzen)


# In[ ]:

try:
    mea = hs.loadmat(fp_starzen)
except: 
    mea = sio.loadmat(fp_starzen)


# In[ ]:

tt = mat2py_time(mea['t'])
mea['utc'] = toutc(tt)


# In[ ]:

print 'Running the parameter calculations on measured spectra'
meas = Sp.Sp(mea)
meas.params()


# ### plot out the different spectra

# In[ ]:

if not noplot:
    print('Making plots...')
    fig1 = Sp.plt_zenrad(meas)
    fig1.savefig(fp_zencld_plot+'{datestr}_zenrad.png'.format(datestr=datestr),
                 dpi=600,transparent=True)
    print('zenrad...'),

    fig1n = Sp.plt_norm_zenrad(meas)
    fig1n.savefig(fp_zencld_plot+'{datestr}_norm_zenrad.png'.format(datestr=datestr),
                  dpi=600,transparent=True)
    print('norm zenrad...'),

    fig2 = Sp.curtain_zenrad(meas,utc=True)
    fig2.savefig(fp_zencld_plot+'{datestr}_curtain_utc_zenrad.png'.format(datestr=datestr),
                 dpi=600,transparent=True)
    print('utc curtain zenrad...'),

    fig2n = Sp.curtain_norm_zenrad(meas,utc=True)
    fig2n.savefig(fp_zencld_plot+'{datestr}_curtain_utc_norm_zenrad.png'.format(datestr=datestr),
                  dpi=600,transparent=True)
    print('utc curtain norm zenrad...'),

    fig3 = Sp.curtain_zenrad(meas,utc=False)
    fig3.savefig(fp_zencld_plot+'{datestr}_curtain_zenrad.png'.format(datestr=datestr),
                 dpi=600,transparent=True)
    print('curtain zenrad...'),

    fig3n = Sp.curtain_norm_zenrad(meas,utc=False)
    fig3n.savefig(fp_zencld_plot+'{datestr}_curtain_norm_zenrad.png'.format(datestr=datestr),
                  dpi=600,transparent=True)
    print('curtain norm zenrad...'),


# ### Load the lut file

# In[90]:

print('Loading the lut file:{}'.format(fp_lut_mat))
if not os.path.isfile(fp_lut_mat):
    print('File {} does not exist'.format(fp_lut_mat))
    raise IOError('LUT File not found: {}'.format(fp_lut_mat))
try:
    luts = hs.loadmat(fp_lut_mat)
except:
    luts = sio.loadmat(fp_lut_mat)


# ### run through the lut for each sza/airmass and create the lut params in hires

# In[ ]:

lut = []


# In[ ]:

for s in xrange(len(luts['sza'])):
    sptemp = luts
    sptemp['wvl'] = [luts['wvl']]
    sptemp['rad'] = luts['rad'][:,:,:,:,:,i]
    ltemp = Sp.Sp(sptemp)
    ltemp.params()
    ltemp.param_hires()
    lut.append(ltemp)


# In[ ]:

airmass = 1./np.cos(ltemp.sza*np.pi/180.0)


# ## Define the good points and start doing the retrieval

# In[ ]:

meas.airmass = 1.0/np.cos(meas.sza*np.pi/180.0)


# In[ ]:

idx = Sp.find_closest(airmass,meas.airmass)


# In[ ]:

(meas.tau,meas.ref,meas.phase,meas.ki) = (np.zeros_like(utc),np.zeros_like(utc),np.zeros_like(utc),np.zeros_like(utc))


# In[ ]:

print 'Running through the airmasses'
for i in np.unique(idx):
    print 'airmass: {airmass}, {i}/{i_tot}'.format(airmass=airmass[i],i=i,i_tot=idx.max())
    meas.good,= np.where(idx==i)
    tau,ref,phase,ki = k.run_retrieval(meas,lut[i])
    meas.tau[meas.good] = tau
    meas.ref[meas.good] = ref
    meas.phase[meas.good] = phase
    meas.ki[meas.good] = ki


# ## Save the retrieval results

# In[ ]:

fp_out = fp+'{datestr}_zen_cld_retrieved.mat'.format(datestr=datestr)
print 'saving to file: {fp_out}'.format(fp_out=fp_out)
hs.savemat(fp_out,meas)


# ## Make output plots

# In[ ]:

if not noplot:
    print 'making the retrieval plots'
    fig,ax = plt.subplots(4,sharex=True)
    ax[0].set_title('Retrieval results time trace')
    ax[0].plot(meas.utc,meas.tau,'rx')
    ax[0].plot(meas.utc[meas.good],smooth(meas.tau[meas.good],20),'k')
    ax[0].set_ylabel('$\\tau$')
    ax[1].plot(meas.utc,meas.ref,'g+')
    ax[1].set_ylabel('R$_{ef}$ [$\\mu$m]')
    ax[1].plot(meas.utc[meas.good],smooth(meas.ref[meas.good],20),'k')
    ax[2].plot(meas.utc,meas.phase,'k.')
    ax[2].set_ylabel('Phase')
    ax[2].set_ylim([-0.5,1.5])
    ax[2].set_yticks([0,1])
    ax[2].set_yticklabels(['liq','ice'])
    ax[3].plot(meas.utc,meas.ki)
    ax[3].set_ylabel('$\\chi^{2}$')
    ax[3].set_xlabel('UTC [Hours]')
    plt.savefig(fp_zencld_plot+'{datestr}_retrieval_out.png'.format(datestr=datestr),dpi=600,transparent=True)


# In[ ]:

if not noplot:
    fig,ax = plt.subplots(1,2)
    ax[0].plot(meas.lon,meas.lat)
    ss = ax[0].scatter(meas.lon,meas.lat,marker='o',c=meas.tau,cmap=plt.cm.gist_rainbow)
    ax[0].set_ylabel('Latitude')
    ax[0].set_xlabel('Longitude')
    cba = plt.colorbar(ss)
    cba.set_label('Cloud optical thickness')
    
    ax[1].plot(meas.lon,meas.lat)
    sr = ax[1].scatter(meas.lon,meas.lat,marker='o',c=meas.tau,cmap=plt.cm.gist_earth)
    ax[1].set_ylabel('Latitude')
    ax[1].set_xlabel('Longitude')
    cbb = plt.colorbar(sr)
    cbb.set_label('Effective radius [$\mu$m]')
    plt.savefig(fp_zencld_plot+'{datestr}_map_retr_zencld.png'.format(datestr=datestr),dpi=600,transparent=True)

