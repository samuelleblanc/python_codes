#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Take incoming folder for sunsat and parse the subfolders and incoming files
# 
# Input:
# 
#     None
# 
# Output:
# 
#     Moved and catagorized files
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - os
#     - dateutil
#     - re
#     - pathlib2
#     - datefinder
# 
# Needed Files:
#   - None
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2020-11-06
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2020-12-02
#              - added support for files without a date in the name, using either the directories' date, or the file's date.
# 

# # Set up the background functions

# In[35]:


from __future__ import print_function 
from parse_and_move_incoming_fx import get_date_and_string,      pull_labels, get_season, filetypes, get_newfilepath, get_filters_from_json,      recurse_through_dir, make_temp_mfile, move_files


# # Prepare the command line argument parser

# In[2]:


import argparse


# In[3]:


long_description = """    Run the incoming file parser and moves the files to the desired sunsat locations
    The File locations and folders are defined by the json file: .filters.json
    Please update the date ranges within that file for any new field mission, 
        if not then assumes rooftop measurements for the season_year
    Can run a call to matlab for any incoming 4STAR raw data"""


# In[4]:


parser = argparse.ArgumentParser(description=long_description)
parser.add_argument('-d','--dry_run',help='if set, turn on dry runs, and not move or delete any file/folder',
                    action='store_true')
parser.add_argument('-q','--quiet',help='if set, quiet the comments',
                    action='store_true')
parser.add_argument('-i','--in_dir',nargs='?',
                    help='Input directory to recurse files, parse, and move',
                    default='/data/sunsat/_incoming_gdrive/')
parser.add_argument('-r','--root_dir',nargs='?',
                    help='full file path of the root directory to save to',
                    default='/data/sunsat/')
parser.add_argument('-m','--run_matlab',help='if set, will run the matlab calls if there is 4STAR/4STARB raw files',
                    action='store_true')


# In[5]:


in_ = vars(parser.parse_known_args()[0])


# # Load the modules and get the defaults

# In[6]:


import os, zipfile
import dateutil.parser
import re
from pathlib2 import Path
from datefinder import find_dates
from datetime import date, datetime
import json
import filecmp
import subprocess
import threading
from aeronet import get_AERONET_file_v2


# In[7]:


in_directory = in_.get('in_dir','/data/sunsat/_incoming_gdrive/')
root_folder = in_.get('root_dir','/data/sunsat/')


# In[8]:


verbose = not in_.get('quiet',False)
dry_run = in_.get('dry_run',True)
run_matlab = in_.get('run_matlab',False)


# In[9]:


if verbose: print( in_)


# In[37]:


# Go through and unzip any folder
prefix = '*DRY RUN*: ' if dry_run else ''
for item in os.listdir(in_directory): # loop through items in dir
    if item.lower().endswith('.zip'): # check for ".zip" extension
        file_name = Path(in_directory+item) # get full path of files
        try: 
            zip_ref = zipfile.ZipFile(str(file_name)) # create zipfile object
            if verbose: 
                print( '{prefix}found zip file: {file_name}, extracting here.'.format(prefix=prefix,file_name=file_name))
            if not dry_run: 
                file_name.parent.joinpath(file_name.stem).mkdir(parents=True,exist_ok=True) # make a dir to extract to
                zip_ref.extractall(str(file_name.parent.joinpath(file_name.stem))) # extract file to dir
            zip_ref.close() # close file
            if not dry_run: os.remove(str(file_name)) # delete zipped file
        except NotImplementedError:
            if verbose: print('{prefix}Error in processing zip file: {file_name}'.format(prefix=prefix,file_name=file_name))


# In[38]:


filters = get_filters_from_json(in_directory)
fl_arr = recurse_through_dir(in_directory,verbose=verbose,filters=filters)


# In[59]:


data_raw_found, data_raw_files, nexact, nmoved, ncreated, ndataraw =                move_files(fl_arr,filters,verbose=verbose,dry_run=dry_run,root_folder=root_folder)


# In[65]:


# check if there are raw data files, and if so, get the aeronet
if data_raw_found:
    daystrss = []
    for k in data_raw_files.keys():
        instname,daystr0 = k.split('_')
        fa_tmp = filetypes('{daystr}_AERONET_NASA_Ames.lev15'.format(instname=instname,daystr=daystr0))
        fmla_tmp = get_newfilepath(fa_tmp,filters=filters,fake_file=True,
                                   root_folder=root_folder,verbose=verbose,dry_run=dry_run)
        if not fa_tmp.newpath.exists() & verbose: print( '{prefix}+Creating new path: {newpath}'.format(**fa_tmp) )
        if not dry_run: fa_tmp.newpath.mkdir(parents=True,exist_ok=True)
        if not daystr0 in daystrss:
            daystrss.append(daystr0)
            if fa_tmp.campaign.find('rooftop') >= 0:
                aeronet_file = get_AERONET_file_v2(date=fa_tmp.fdate,site='NASA_Ames',path=str(fa_tmp.newpath),version=3)
                if verbose: print('Downloaded AERONET file: {}'.format(aeronet_file))
            if fa_tmp.campaign.find('MLO') >= 0:
                aeronet_file = get_AERONET_file_v2(date=fa_tmp.fdate,site='Mauna_Loa',path=str(fa_tmp.newpath),version=3)
                if verbose: print('Downloaded AERONET file: {}'.format(aeronet_file))
            if fa_tmp.campaign.find('SaSa') >= 0:
                aeronet_file = get_AERONET_file_v2(date=fa_tmp.fdate,site='WFF_X-75_Sci_Obs',path=str(fa_tmp.newpath),version=3)
                if verbose: print('Downloaded AERONET file: {}'.format(aeronet_file))

# In[66]:


# clean up folders after move
for dirpath, dirnames, filenames in os.walk(in_directory,topdown=False):
    if not dirpath in in_directory:
        try: 
            if verbose: print( '{pre}-removing :{path}'.format(pre=prefix,path=dirpath))
            if not dry_run:
                os.rmdir(dirpath) 
        except: 
            pass


# In[67]:


nmats = 0
if run_matlab:
    prefix = '*DRY RUN*: ' if dry_run else ''
    for dr,drs in data_raw_files.items():
        # get the position of the new star.mat and starsun.mat files
        f = filetypes('{}star.mat'.format(dr),filters=filters)
        fml = get_newfilepath(f,filters=filters,fake_file=True,root_folder=root_folder)
        if not dry_run: f.newpath.mkdir(parents=True,exist_ok=True)
        fs = filetypes('{}starsun.mat'.format(dr),filters=filters)
        fmls = get_newfilepath(fs,filters=filters,fake_file=True,root_folder=root_folder)
        if not dry_run: fs.newpath.mkdir(parents=True,exist_ok=True)
            
        # make the position of the new quicklook file
        instname,daystr = dr.split('_')
        fq = filetypes('{daystr}_{instname}_Quicklooks.pptx'.format(instname=instname,daystr=daystr))
        fmlq = get_newfilepath(fq,filters=filters,fake_file=True,root_folder=root_folder)
        if not dry_run: fq.newpath.mkdir(parents=True,exist_ok=True)
        
        # make the position of the new figure files
        ff = filetypes('{daystr}_{instname}_plots.png'.format(instname=instname,daystr=daystr))
        fmlf = get_newfilepath(ff,filters=filters,fake_file=True,root_folder=root_folder)
        if not dry_run: ff.newpath.parent.mkdir(parents=True,exist_ok=True)
        
        # make the position of the aeronet files
        fa = filetypes('{daystr}_AERONET_NASA_Ames.lev15'.format(instname=instname,daystr=daystr))
        fmla = get_newfilepath(fa,filters=filters,fake_file=True,root_folder=root_folder)
        if not dry_run: fa.newpath.mkdir(parents=True,exist_ok=True)
        
        # make the position of the gas_summary files
        fg = filetypes('{instname}_{daystr}_gas_summary.mat'.format(instname=instname,daystr=daystr))
        fmlg = get_newfilepath(fg,filters=filters,fake_file=True,root_folder=root_folder)
        if not dry_run: fg.newpath.mkdir(parents=True,exist_ok=True)
        
        # make a string of the raw files    
        filelist = "'"+"';'".join(drs)+"'"
        if not f.instname in ['4STAR','4STARB']: # only for 4STARs for now.
            continue
            
        mfile = make_temp_mfile(in_directory+'temp_{}_{}.m'.format(f.instname,daystr),filelist=filelist,starmat=str(f.newfile),                                starsun=str(fs.newfile),quicklooks=in_directory+str(fq.newfile.name),                                fig_path=str(ff.newpath.parent)+'/', aero_path=str(fa.newpath)+'/',                                gas_path=str(fg.newpath)+'/', sun_path=str(fs.newpath)+'/',incoming_path=in_directory)

        pmfile = Path(mfile)
        def which(pgm):
            path=os.getenv('PATH')
            for p in path.split(os.path.pathsep):
                p=os.path.join(p,pgm)
                if os.path.exists(p):
                    return p

        mat_path = which('matlab')
        if not mat_path: mat_path = '/usr/local/bin/matlab'
        if verbose: 
            print( ' '.join(['{}{}'.format(prefix,mat_path),'-nodisplay','-batch',"cd('{}');{},exit;".format(str(pmfile.parent),pmfile.stem)]))
        if not dry_run:
            os.chdir(str(pmfile.parent))
            process = subprocess.Popen([mat_path,'-nodisplay','-batch',"cd('{}');{},exit;".format(str(pmfile.parent),str(pmfile.stem))],
                                       shell=False, stdout=subprocess.PIPE,stderr=subprocess.PIPE)

            while True:
                # handle output by direct access to stdout and stderr
                output = process.stdout.readline()
                if process.poll() is not None:
                    break
                if output:
                    if verbose: print(output.strip())
            rc = process.poll()
            if verbose: print(rc)
            nmats = nmats + 1
                
            if rc==0:
                os.remove(mfile)
            else:
                print(process.stderr.readline())
                print('ERROR with matlab')
                if not f.newfile.exists():
                     print(' ***  file {} has not been created'.format(str(f.newfile)))
                if not fs.newfile.exists():
                     print(' ***  file {} has not been created'.format(str(fs.newfile)))
                     nmats = nmats-1


# In[68]:


print(datetime.now().strftime("%c")+' :Python moved {nmoved} files, Created {ncreated} folders, found {ndataraw} files, and generated {nmats} starmats/suns'      .format(nmoved=nmoved,ncreated=ncreated,ndataraw=ndataraw,nmats=nmats))


# In[ ]:




