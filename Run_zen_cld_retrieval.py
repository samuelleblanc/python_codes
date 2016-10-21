
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
#         forceliq: (optional) if set, then retrieval will only use the liquid cloud
#         plotlut: (optional) if set, then it will create a plot of the lut
#         makemovie: (optional) if set, then it will create a movie of the measured radiances, normalized,
#                               and time series of measured radiances, including retrievals 
#         makeslides: (optional) if set, will create a powerpoint slide of the different figures created
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
#     - load_utils for mat2py_time and toutc
#   
# Needed Files:
# 
#   - lut .mat file created by Read_save_NAAMES_lut
#   - starzen.mat file created by the matlab program : starzen.m
#   
# Example:
# 
#     >> python Run_zen_cld_retrieval.py -lut /nobackup/sleblan2/NAAMES/model/v2_NAAMES_lut.mat /nobackup/sleblan2/NAAMES/c130/20151114starzen.mat -o /u/sleblan2/NAAMES/plot/
# loading file /nobackup/sleblan2/NAAMES/c130/20151114starzen.mat
# Exception occured: unable to create file (File accessibilty: Unable to open file)
# Running the parameter calculations on measured spectra
# Making plots...
# zenrad... norm zenrad... utc curtain zenrad... utc curtain norm zenrad... curtain zenrad... curtain norm zenrad... Loading the lut file: /nobackup/sleblan2/NAAMES/model/v2_NAAMES_lut.mat
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running interpolation on params: [########################################]100% -- Done!
# Running through the airmasses
# airmass: 3.2360679775, 11/3
# ['tau', 'verbose', 'sza', 'stdparn', 'par', 'hiresirrad', 'meansp', 'iwvls', 'stdsp', 'ref', 'norm', 'iset', 'parn', 'stdpar', 'reflect', 'isubwvl', 'wvl', 'npar', 'sphires', 'refsp', 'parhires', 'sp', 'irrad', 'tausp']
# meas.good lengh: (142,),meas.utc length: (850,)
# In Rk, Measurement space length: (850,)
# Retrieval progress over times: [########################################]100% -- Done!
# airmass: 3.42030361983, 12/3
# ['tau', 'verbose', 'sza', 'stdparn', 'par', 'hiresirrad', 'meansp', 'iwvls', 'stdsp', 'ref', 'norm', 'iset', 'parn', 'stdpar', 'reflect', 'isubwvl', 'wvl', 'npar', 'sphires', 'refsp', 'parhires', 'sp', 'irrad', 'tausp']
# meas.good lengh: (48,),meas.utc length: (850,)
# In Rk, Measurement space length: (850,)
# Retrieval progress over times: [########################################]100% -- Done!
# airmass: 3.86370330516, 14/3
# ['tau', 'verbose', 'sza', 'stdparn', 'par', 'hiresirrad', 'meansp', 'iwvls', 'stdsp', 'ref', 'norm', 'iset', 'parn', 'stdpar', 'reflect', 'isubwvl', 'wvl', 'npar', 'sphires', 'refsp', 'parhires', 'sp', 'irrad', 'tausp']
# meas.good lengh: (660,),meas.utc length: (850,)
# In Rk, Measurement space length: (850,)
# Retrieval progress over times: [########################################]100% -- Done!
# saving to file: /nobackup/sleblan2/NAAMES/c13020151114_zen_cld_retrieved.mat
# making the retrieval plots
# making the retrieval histogram plots
# making the map
# 
# History:
# 
#     Written: Samuel LeBlanc, NASA Ames, 2015-10-26
#     Modified: Samuel LeBlanc, NASA Ames, 2016-10-06
#             - added force liquid keyword to call and plot_lut for plotting the lut.
#             - added making a movie with the results
#             - added filtering of data prior to ensure good retrievals
#     Modified: Samuel LeBlanc, NASA Ames, 2016-10-20
#             

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
parser.add_argument('fp_zencld',nargs='?',
                    help='''(optional) the full path of the output file (with .mat extension) of the saved retrieval file
                   if omitted, uses the same directory as starzen.mat, but saves file as *_zencld.mat''')
parser.add_argument('-o','--fp_zencld_plot',nargs='?',
                    help="""(optional) the full path of the directory to save the plots from the retrieval
                        if omitted uss the same directory as starzen.mat""")
parser.add_argument('-lut','--fp_lut_mat',nargs='?',help='Put in special look up table file path')
parser.add_argument('-f',nargs='?',help='not used')
parser.add_argument('-n','--noplot',help='if set, will not output plots',action='store_true')
parser.add_argument('-liq','--forceliq',help='if set, will force to only use liquid lut',action='store_true')
parser.add_argument('-plotlut','--plotlut',help='if set, will plot the look up table',action='store_true')
parser.add_argument('-movie','--makemovie',help='if set, will create a movie of the measurements and retrievals',
                    action='store_true')
parser.add_argument('-slides','--makeslides',help='if set, will create a powerpoint of figures created',
                    action='store_true')
parser.add_argument('-noflt','--nofilter',help='if set, will not filter out bad measurements',
                    action='store_true')
parser.add_argument('-db','--debug',help='if set, turn on debug comments and set_trace() at a few locations',
                    action='store_true')


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
from load_utils import mat2py_time, toutc
import Sp_parameters as Sp
from mpltools import color
import scipy.io as sio
import run_kisq_retrieval as rk
try:
    from mpl_toolkits.basemap import Basemap
    isbasemap = True
except:
    isbasemap = False
try:
    import warnings
    warnings.simplefilter(action = "ignore", category = FutureWarning)
except:
    pass


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
        filename = asksaveasfilename(defaultextension='*.mat',filetypes=[('mat file','*.mat'),('All files','*.*')]) 
        # show an "Open" dialog box and return the path to the selected file
        fp_starzen = os.path.abspath(filename)
    else:
        fp_starzen = in_.get('fp_starzen')


# In[75]:

fp,fname = os.path.split(fp_starzen)
if fname[0:5]=='4STAR':
    datestr = fname[6:14]
else:
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
    fp_lut_mat = fp+os.path.sep+'v2_NAAMES_lut.mat'


# In[ ]:

if in_.get('noplot'):
    noplot = True
else:
    noplot = False


# In[ ]:

if in_.get('forceliq'):
    forceliq = True
    print 'Force liquid cloud'
else:
    forceliq = False


# In[ ]:

plot_lut = in_.get('plotlut',False)
make_movie = in_.get('makemovie',False)
no_filter = in_.get('nofilter',False)
debug = in_.get('debug',False)
make_slides = in_.get('makeslides',False)


# # Load starzen files and the look up table

# In[79]:

if not os.path.isfile(fp_starzen):
    raise IOError('file {} is not found'.format(fp_starzen))


# In[ ]:

print 'loading file {}'.format(fp_starzen)


# In[ ]:

try:
    mea = hs.loadmat(fp_starzen)
except Exception as ei:
    if debug: print 'Exception occured: {}'.format(ei)
    mea = sio.loadmat(fp_starzen)


# In[ ]:

tt = mat2py_time(mea['t'])
mea['utc'] = toutc(tt)


# In[ ]:

print 'Running the parameter calculations on measured spectra'
meas = Sp.Sp(mea,verbose=False)
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
except Exception as ei:
    if debug: print 'Exception occured when loading lut : {} *** Trying second method'.format(ei)
    luts = sio.loadmat(fp_lut_mat)


# ### run through the lut for each sza/airmass and create the lut params in hires

# In[ ]:

airmass = 1./np.cos(luts['sza']*np.pi/180.0)


# In[ ]:

meas.airmass = 1.0/np.cos(meas.sza*np.pi/180.0)
if len(airmass) > 1:
    idx = Sp.find_closest(airmass,meas.airmass)
else: 
    print '*** There is only one airmass in the lut: Using {} ***'.format(airmass)
    idx = np.zeros_like(meas.airmass,dtype=np.int)


# In[ ]:

lut = []


# In[ ]:

for s in xrange(len(luts['sza'])):
    sptemp = {}
    sptemp['tau'] = luts['tau']
    sptemp['ref'] = luts['ref']
    sptemp['zout'] = luts['zout']
    sptemp['sza'] = luts['sza']
    sptemp['phase'] = luts['phase']
    sptemp['irr_dn_diff'] = luts['irr_dn_diff'][:,:,:,:,:,s]
    sptemp['irr_dn'] = luts['irr_dn'][:,:,:,:,:,s]
    sptemp['irr_up'] = luts['irr_up'][:,:,:,:,:,s]
    sptemp['wvl'] = [luts['wvl']]
    sptemp['rad'] = luts['rad'][:,:,:,:,:,s]
    ltemp = Sp.Sp(sptemp,verbose=False)
    if s in idx:
        ltemp.params()
        ltemp.param_hires()
    lut.append(ltemp)


# ## Define the good points and start doing the retrieval

# In[ ]:

meas.airmass = 1.0/np.cos(meas.sza*np.pi/180.0)
if len(airmass) > 1:
    idx = Sp.find_closest(airmass,meas.airmass)
else: 
    print '*** There is only one airmass in the lut: Using {} ***'.format(airmass)
    idx = np.zeros_like(meas.airmass,dtype=np.int)


# In[ ]:

(meas.taut,meas.ref,meas.phase,meas.ki) = (np.zeros_like(meas.utc),np.zeros_like(meas.utc),np.zeros_like(meas.utc),np.zeros_like(meas.utc))


# In[ ]:

(meas.taut,meas.ref,meas.phase,meas.ki) = (meas.taut*np.nan,meas.ref*np.nan,meas.phase*np.nan,meas.ki*np.nan)


# In[ ]:

print 'Running through the airmasses'
for i in np.unique(idx):
    try: 
        print 'airmass: {airmass}, {i}/{i_tot}'.format(airmass=airmass[i],i=i,i_tot=idx.max()-idx.min())
    except:
        if debug:
            import pdb; pdb.set_trace()
        else: 
            print 'Exception occurred at airmass: {} *** skipping'.format(airmass[i])
    meas.good = np.where(idx==i)[0]
    if not no_filter:
        i500 = np.argmin(abs(meas.wvl-500))
        #i980,i995 = np.argmin(abs(meas.wvl-980)),np.argmin(abs(meas.wvl-995))
        #i999,i992 = np.argmin(abs(meas.wvl-999)),np.argmin(abs(meas.wvl-992))
        i981,i982 = 1039,1040
        #ss = np.nanstd(meas.norm[meas.good,i980:i995],axis=1)/np.nanmean(meas.norm[meas.good,i980:i995],axis=1)
        sss = abs(meas.norm[meas.good,i981]-meas.norm[meas.good,i982])
        #flt = (meas.norm[meas.good,i500]>0.4) & (ss<0.05)
        flt = sss<0.1 #(meas.norm[meas.good,i500]>0.4) & (sss<0.1)
        #import pdb; pdb.set_trace()
        meas.good = meas.good[flt]
    if debug:
        try:
            print lut[i].__dict__.keys()
        except Exception as ei:
            print 'exception: {}'.format(ei)
            import pdb; pdb.set_trace()
    print 'meas.good lengh: {},meas.utc length: {}'.format(meas.good.shape,meas.utc.shape)
    tau,ref,phase,ki = rk.run_retrieval(meas,lut[i],force_liq=forceliq)
    meas.taut[meas.good] = tau[meas.good]
    meas.ref[meas.good] = ref[meas.good]
    meas.phase[meas.good] = phase[meas.good]
    meas.ki[meas.good] = ki[meas.good]


# In[ ]:

meas.tau = meas.taut


# ## Save the retrieval results

# In[ ]:

mdict = {'tau':meas.tau,'ref':meas.ref,'phase':meas.phase,'ki':meas.ki,
         'utc':meas.utc,'sza':meas.sza,'lat':meas.lat,'lon':meas.lon,'alt':meas.alt}


# In[ ]:

fp_out = os.path.join(fp,'{datestr}_zen_cld_retrieved.mat'.format(datestr=datestr))
print 'saving to file: {fp_out}'.format(fp_out=fp_out)
hs.savemat(fp_out,mdict)


# ## Make output plots

# In[ ]:

if not noplot:
    print 'making the retrieval plots'
    figz = Sp.plot_zen_cld_retrieval(meas)
    figz.savefig(fp_zencld_plot+'{datestr}_retrieval_out.png'.format(datestr=datestr),dpi=600,transparent=True)


# In[ ]:

if not noplot:
    print 'making the retrieval histogram plots'
    figh = Sp.plot_hist_cld_retrieval(meas)
    figh.savefig(fp_zencld_plot+'{datestr}_retrieval_hist_out.png'.format(datestr=datestr),dpi=600,transparent=True)


# In[ ]:

if not noplot:
    print 'making the map'
    figm = Sp.plot_map_cld_retrieval(meas)
    figm.savefig(fp_zencld_plot+'{datestr}_map_retr_zencld.png'.format(datestr=datestr),dpi=600,transparent=True)


# In[ ]:

if plot_lut:
    print 'making the lut plots'
    for i,ll in enumerate(lut):
        if len(lut[i].par)>4:
            continue
        if not i in np.unique(idx):
            continue
        try: 
            figl,axl = Sp.plot_lut_vs_tau(lut[i],forceliq=forceliq)
            figl.savefig(fp_zencld_plot+'{datestr}_lut_only_sza_{sza:02.1f}.png'.format(datestr=datestr,sza=ll['sza'][i]),
                         dpi=600,transparent=True)
        except:
            print '*** Problem with plotting of sza:{sza:02.1f}'.format(sza=ll['sza'][i])
            if debug:
                import pdb; pdb.set_trace()
        try:
            ik = (idx==np.int32(i))
            npar = len(meas.par[0,:])-1
            if ik.any():
                for ai,ax in enumerate(axl[0:npar]):
                    try:
                        ax.scatter(meas.tau[ik[:,0]],meas.par[ik[:,0],ai],marker='x',zorder=200,color='k')
                    except:
                        print '** problem inside multiple'
                        if debug:
                            import pdb; pdb.set_trace()

                figl.savefig(fp_zencld_plot+'{datestr}_lut_with_data_sza_{sza:02.1f}.png'.format(
                        datestr=datestr,sza=ll['sza'][i]),
                             dpi=600,transparent=True)
        except:
            print '*** Problem with plotting values on top of lut of sza:{sza:02.1f}'.format(sza=ll['sza'][i])
            if debug:
                import pdb; pdb.set_trace()


# In[ ]:

if make_movie:
    print 'Making the movie of measurement and retrievals'
    Sp.plot_sp_movie(meas,fp_zencld_plot,gif=False)


# In[ ]:

if make_slides:
    print 'Making the powerpoint slides of the figures'
    from plotting_utils import make_pptx
    make_pptx(fp_zencld_plot,'Cloud retrievals {}'.format(datestr),title='4STAR Zenith Cloud retrievals for {}'.format(datestr),
              glob_pattern='*{}_????[!only]*?.png'.format(datestr))

