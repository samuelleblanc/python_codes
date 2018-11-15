
# coding: utf-8

# # Intro
# Loading of files from previous flight to checkout the winds during flight

# In[1]:


get_ipython().magic(u'config InlineBackend.rc = {}')
import matplotlib 
import os
matplotlib.rc_file(os.path.join(os.getcwd(),'file.rc'))
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import Sp_parameters as Sp
from load_utils import mat2py_time, toutc, load_ict
from Sp_parameters import smooth
from linfit import linfit
from path_utils import getpath
from plotting_utils import make_boxplot


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


fp =getpath('ORACLES')#'C:/Userds/sleblan2/Research/ORACLES/'
fp


# # Load files

# In[4]:


d = np.genfromtxt(fp+'flight_planning/2018/RF11/IWG1.17Oct2018-1514.txt')


# In[9]:


d.dtype.names

