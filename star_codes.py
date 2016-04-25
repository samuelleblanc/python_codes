
# coding: utf-8

# # Intro

# In[1]:

def __init__():
    """
Name:  

    star_codes

Purpose:  

    To compile all the 4STAR analysis codes within a single module
    Should be compared to matlab codes created by Yohei
    contains multiple codes including:
        - allstarmat - compile raw ascii files into a matlab binary file
        - starsun = run through typical analysis of direct beam analysis
        - starzen - for zenith radiance analysis
        - others...
    See the help at each individual method

Calling Sequence:

    import star_codes as sc

Input:

    See each individual method within the module

Output:
   
    Save matlab file of raw or processed files
    Some ploting utilities
  
Keywords:

    none
  
Dependencies:

    - numpy
    - hdf5storage : for saving and reading
    - Run_Libradtran
    - os
    - sys
    - matplotlib
    - argparse
    - Tkinter (for file open dialogs)
    - load_modis for mat2py_time and toutc
    - warnings
  
Needed Files:

  - raw ascii files or yyyymmddstar.mat file
  
Example:

   to come...

History:

    Written: Samuel LeBlanc, Flight from Honolulu to SFO over the Pacific, 2016-04-22
    """
    import os
    import sys
    import numpy as np
    impot warnings
    
    pass


# # Create the reading functions

# In[ ]:

def starread(f_in):
    """
Name:  

    starread

Purpose:  

    Reads a single ascii raw data file and creates a dictionary for the type of file

Calling Sequence:

    dict_out = starread(f_in)

Input:

    f_in: full path to filename

Output:
   
    dict with the correct elelements of the file
  
Keywords:

    none
  
Dependencies:

    - numpy
    - os
    - sys
    - load_modis for mat2py_time and toutc
  
Needed Files:

  - raw ascii files
  
Example:

   to come...

History:

    Written: Samuel LeBlanc, Flight from Honolulu to SFO over the Pacific, 2016-04-22
    """
    typ = startype(f_in)
    
    if typ=='nir_sun' or typ=='vis_sun':
        head,numhead = get_header(f_in)
        d = np.genfromtxt(f_in,names=True,skip_header=numhead)


# In[36]:

def startype(f_in,return_daystr=False):
    """
Name:  

    startype

Purpose:  

    returns the type of the file to be read

Calling Sequence:

    typ = starread(f_in)

Input:

    f_in: full path to filename

Output:
   
    typ: ile type as a string - recreates the exact string as in the file name
  
Keywords:

    return_daystr: (default False) if true also returns the daystr from the file name
  
Dependencies:

    - os
    - sys
  
Needed Files:

  - raw ascii files (.dat)
  
Example:

   to come...

History:

    Written: Samuel LeBlanc, Flight from Honolulu to SFO over the Pacific, 2016-04-22
    """
    if not os.path.isfile(f_in):
        raise IOError('File not found! please check the path {}'.format(f_in))
    fparts = os.path.split(f_in)
    fname = fparts[1].split('.')[0].split('_')
    daystr = fname[0]
    typ = fname[-1]
    fnum = fname[1]
    if return_daystr:
        return typ,daystr
    else:
        return typ


# In[ ]:

def get_header(f):
    """
    Simple function to extract the header and number of header lines from a file
    
    input:
     - f: full file path
     
    output:
     - head: list of header lines
     - numhead: number of header lines to be fed into the np.genfromtxt
    """
    head = []
    with open(f,'r') as ff:
        for line in ff:
            if line[0].isdigit():
                break
            else:
                head.append(line)
    numhead = len(head)-1
    return head,numhead


# In[37]:

import numpy as np


# In[66]:

d = np.genfromtxt(f_in,names=True,delimiter=':',skip_header=1)


# In[60]:

d.dtype.names


# In[67]:

help(np.genfromtxt)


# In[ ]:



