
# coding: utf-8

# # Intro
# Name:  
# 
#     Create_ic_phase_files
# 
# Purpose:  
# 
#     Read the phase functions from libradtran describing the ice phase functions from Baum 2014. Read in the mat file translated     from the netcdf file
#   
# Input:
# 
#     none
# 
# Output:
#    
#     input files for pmom in libradtran 2.0 
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
# Needed Files:
# 
#   - ic.ghm.baum.mat
#     
# History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2017-02-15

# # Prepare the python environment

# In[2]:

import numpy as np
import scipy.io as sio
import os
import Run_libradtran as RL
reload(RL)


# In[3]:

if os.sys.platform == 'win32':
    fp = 'C:\\Users\\sleblan2\\Research\\libradtran\\'
    fp_rtm = 'C:\\Users\\sleblan2\\Research\\libradtran\\ice_pmom\\'
    fp_pmom = 'C:\\Users\\sleblan2\\Research\\libradtran\\libRadtran-2.0-beta\\bin\\pmom'
elif os.sys.platform == 'linux2':
    fp = '/u/sleblan2/libradtran/'
    fp_rtm = '/nobackup/sleblan2/rtm/'
    fp_pmom = '/u/sleblan2/libradtran/libRadtran-2.0-beta/bin/pmom'
else:
    raise Exception


# In[7]:

vv = 'v1'


# In[6]:

name = 'ic_phase'


# ## Setup command line arguments

# In[52]:

import argparse


# In[53]:

long_description = """    Prepare the pmom function using libradtran's pmom function 
      from the Baum 2014 General Habit mixture ice phase functions
      save them using the doread argument"""


# In[54]:

parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-doread','--doread',help='if set, will only read the output, not produce them',
                    action='store_true')


# In[55]:

in_ = vars(parser.parse_args())
doread = in_.get('doread',False)


# # Load the phase functions

# In[4]:

ic = sio.loadmat(fp+'ic.ghm.baum.mat')


# In[5]:

ic.keys()


# In[11]:

ic['phase'].shape


# In[12]:

ic['ref'].shape


# In[13]:

ic['wvl'].shape


# In[15]:

ic['theta'].shape


# # Run through the phase functions and create files

# In[8]:

fp_in = os.path.join(fp_rtm,'input','{vv}_{name}'.format(vv=vv,name=name))
fp_out = os.path.join(fp_rtm,'output','{vv}_{name}'.format(vv=vv,name=name))


# In[9]:

if not os.path.exists(fp_in):
    os.makedirs(fp_in)
if not os.path.exists(fp_out):
    os.makedirs(fp_out)


# In[10]:

if not doread:
    f_list = open(os.path.join(fp,'{name}_list_{vv}.sh'.format(vv=vv,name=name)),'w')
    print f_list.name


# In[29]:

fname = 'ic_ref{r:04.1f}_wvl{w:0.2f}.{e}'


# In[39]:

def write_phase(f,angle,phase):
    'write the phase function file defined by path f, with angle, and phase'
    with open(f,'w') as fo:
        for i in range(len(angle)-1,0,-1):
            fo.write('{} {}\n'.format(angle[i],phase[i]))


# In[43]:

if not doread:
    for ir,r in enumerate(ic['ref'][0]):
        for iw,w in enumerate(ic['wvl'][0]):
            if w > 0.39 and w < 1.71:
                print r,w
                x = 2.0*np.pi*r/w
                write_phase(os.path.join(fp_in,fname.format(r=r,w=w,e='dat')),ic['theta'][iw,ir,0,:],ic['phase'][iw,ir,0,:])
                f_list.write(fp_pmom+' -n -r 3 -c -l {}'.format(int(10*np.pi*r**1.75/w))+os.path.join(fp_in,fname.format(r=r,w=w,e='dat'))+' > '+
                             os.path.join(fp_out,fname.format(r=r,w=w,e='out'))+'\n')
            


# In[ ]:

if not doread:
    f_list.close()


# In[ ]:

if doread:
    pmom = []
    for ir,r in enumerate(ic['ref'][0]):
        pw = []
        wv = []
        iwv = []
        for iw,w in enumerate(ic['wvl'][0]):
            if w > 0.39 and w < 1.71:
                print r,w
                wv.append(w)
                iwv.append(iw)
                with open(os.path.join(fp_out,fname.format(r=r,w=w,e='out')),'r') as f:
                    a = map(float,f.readline().split())
                pw.append(a)
        pmom.append(pw)


# In[ ]:

if doread:
    ic['pmom'] = pmom
    ic['pmom_wvl'] = wv
    ic['pmom_iwvl'] = iwv
    print 'saving to: {}'.format(fp+'ic.pmom.ghm.baum.mat')
    sio.savemat(fp+'ic.pmom.ghm.baum.mat',ic)

