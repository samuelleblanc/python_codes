
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
    - laod_utils for mat2py_time and toutc
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
    import warnings
    from datetime import datetime
    
    pass


# # Create the reading functions

# In[4]:

uu = ['a','c','b','t','d','002','001','000']


# In[6]:

uu.sort()


# In[8]:

uu[0]


# In[ ]:

def allstarmat(f_list=None, save_file_path=None):
    """
Name:  

    allstarmat

Purpose:  

    Runs through many raw data file and reads each one, compiles a dict of containing each measurement type
    Output same as matlab

Calling Sequence:

    s = allstarmat(f_list,save_file_path)

Input:

    f_list: (optional) list with full path of files
    save_file_path: (optional) full file path, with filename of the saved combined mat file

Output:
   
    if save_file_path not set, returns the file path
  
Keywords:

    f_list
    save_file_path
  
Dependencies:

    - numpy
    - os
    - sys
  
Needed Files:

  - raw ascii files
  
Example:

   to come...

History:

    Written: Samuel LeBlanc, NASA DC8 on transit from Palmdale to Osan, KORUS-AQ, 2016-04-26
    """
    if not f_list:
        f_list = select_files()
    
    typ,fnum,datestr = startype(f_list.sort()[0],return_daystr=True)  #get the datestr from the first file
    if not save_file_path:
        save_file_path = select_save_path()
    save_file = save_file_path+'{date}star.mat'.format(date=datestr)
    
    s = {}
    for f in f_list:
        typ,fnum = startyp(f)
        if typ in s:
            s[typ] = concat_dict(s[typ],starread(f))
            s[typ]['filename'].append(f)
            fparts = os.path.split(f)
            s[typ]['fname'].append(fparts[0])
            s[typ]['pname'].append(fparts[-1])
        else:
            s[typ] = starread(f)
            s[typ]['filename'] = [f]
            fparts = os.path.split(f)
            s[typ]['fname'] = [fparts[0]]
            s[typ]['pname'] = [fparts[-1]]
    
    


# In[ ]:

def concat_dict(d0,d1):
    """
    Simple program to concatenate the values of two dictionaries
    
    """


# In[ ]:

def select_files():
    """
    Simple gui function to let user select files to read
    
    outputs list of string, with each string a full file path
    """
    
    print 'to come...'


# In[ ]:

def select_save_path():
    """
    Another gui function to let user select the path of the file to be saved
    
    outputs a string with full file path
    """
    print 'to come...'


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
    - load_utils for mat2py_time and toutc
  
Needed Files:

  - raw ascii files
  
Example:

   to come...

History:

    Written: Samuel LeBlanc, Flight from Honolulu to SFO over the Pacific and transit from Palmdale to Osan, 2016-04-26
    """
    typ,fnum = startype(f_in)
    head,numhead = get_header(f_in)
    
    if typ=='track': 
        def mkktime(txt):
            o = datetime.strptime(txt,'%H:%M:%S.%f')
            return o.hour+o.minute/60.0+o.second/3600.0+o.microsecond/3600.0/1000000.0
        new_nm =  ['N','t','Az_deg','Az_corr','El_deg','El_corr','Az_step','El_step',           'QdV_TB','QdV_LR','QdV_tot','T_box_set','T_box_mon','T_box','T_spec_uvis',           'T_spec_nir','P1','P2','P3','P4','T1','T2','T3','T4'] # change the names from the 
        def convt(t):
            return float(t)*1000.0-273.15
        conv = {'t':mkktime,'T_box':lambda t:float(t)*100.0-273.15,'T1':convt,'T2':convt,'T3':convt,'T4':convt}
        dd = np.genfromtxt(h,skip_header=2,naames=new_nm,delimiter='\t',converters=conv)
        d = to_dict(dd)
        
        # special variables calculated from the track file
        d['El_temp_K'] = 1000.0*d['T1']
        d['NearWindow_temp_K'] = 1000.0*d['T2']
        d['Az_temp_K'] = 1000.0*d['T4']
        d['Can_RH'] = (d['P2']/5.0-0.16)/0.0062 # There is a sign error in the offset.  Need confirmation from Roy.
        d['Can_P_mB'] = 1000.0*d['P1']/0.23
        
    else:
        dd = np.genfromtxt(f_in,names=True,skip_header=numhead)
        d = to_dict(dd)
    
    d['filen'] = np.empty(len(d['t']),dtype=int)*0+fnum
    d['head'] = head
    return d


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
    typ = fname[-1].lower()
    fnum = int(fname[1])
    if return_daystr:
        return typ,fnum,daystr
    else:
        return typ,fnum


# In[ ]:

def to_dict(d):
    'Simple function to switch from numpy named array to dict with numpy arrays'
    dd = {}
    for n in d.dtype.names:
        dd[n] = d[n]
    return d


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



