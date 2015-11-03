
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

# In[81]:

import sys
import os
import matplotlib.pyplot as plt
import Run_libradtran as RL
import hdf5storage as hs
import numpy as np
import argparse
from load_modis import mat2py_time, toutc
import Sp_parameters as Sp
from mpltools import color
import scipy.io as sio


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


# ## Parse command line and get appropriate paths

# In[70]:

in_ = vars(parser.parse_args())


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


# # Load starzen files and the look up table

# In[79]:

if not os.path.isfile(fp_starzen):
    raise IOError()'file {} is not found'.format(fp_starzen))


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

print('Making plots...')


# In[ ]:

fig1 = Sp.plt_zenrad(meas)
fig1.savefig(fp_zencld_plot+'{datestr}_zenrad.png'.format(datestr=datestr),
             dpi=600,transparent=True)
print('zenrad...'),


# In[ ]:

fig1n = Sp.plt_norm_zenrad(meas)
fig1n.savefig(fp_zencld_plot+'{datestr}_norm_zenrad.png'.format(datestr=datestr),
              dpi=600,transparent=True)
print('norm zenrad...'),


# In[ ]:

fig2 = Sp.curtain_zenrad(meas,utc=True)
fig2.savefig(fp_zencld_plot+'{datestr}_curtain_utc_zenrad.png'.format(datestr=datestr),
             dpi=600,transparent=True)
print('utc curtain zenrad...'),


# In[ ]:

fig2n = Sp.curtain_norm_zenrad(meas,utc=True)
fig2n.savefig(fp_zencld_plot+'{datestr}_curtain_utc_norm_zenrad.png'.format(datestr=datestr),
              dpi=600,transparent=True)
print('utc curtain norm zenrad...'),


# In[ ]:

fig3 = Sp.curtain_zenrad(meas,utc=False)
fig3.savefig(fp_zencld_plot+'{datestr}_curtain_zenrad.png'.format(datestr=datestr),
             dpi=600,transparent=True)
print('curtain zenrad...'),


# In[ ]:

fig3n = Sp.curtain_norm_zenrad(meas,utc=False)
fig3n.savefig(fp_zencld_plot+'{datestr}_curtain_norm_zenrad.png'.format(datestr=datestr),
              dpi=600,transparent=True)
print('curtain norm zenrad...'),


# ### Load the lut file

# In[90]:

print('Loading the lut file:{}'.format(fp_lut_mat))
if not os.path.isfile(fp_lu_mat):
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



