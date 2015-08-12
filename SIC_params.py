# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Name:  
# 
#     SIC_params
# 
# Purpose:  
# 
#     Python script to looking at Shannon Information Content of the zenith cloud retrievals parameters
# 
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
#     - 
#   
# Needed Files:
# 
#   - retr_day_pss_v4.out
#   - retr_day_sub_v4.out
#     
#  Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, CA, 2015-08-12
#             

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpltools import color
%matplotlib inline
import numpy as np, h5py
import scipy.io as sio
import scipy
import math, os, IPython
import Sp_parameters as Sp
import load_modis as lm
IPython.InteractiveShell.cache_size = 0
# set the basic directory path
fp='C:/Users/sleblan2/Research/SSFR3/'

# <codecell>

fppss = fp+'retrieved/cst/retr_day_pss_v4.out'
fpsub = fp+'retrieved/cst/retr_day_sub_v4.out'

# <codecell>

pss = sio.idl.readsav(fppss)
sub = sio.idl.readsav(fpsub)

# <codecell>

pss.keys()

# <codecell>

sub.keys()

# <codecell>


