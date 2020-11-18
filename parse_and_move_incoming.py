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
#     Modified:
# 

# # Set up the background functions

# In[1]:


def get_date_and_string(p):
    'p: posix path to extract the date and the string'
    fd = [dtemp for dtemp in find_dates(p.stem,source=True)]
    if fd:
        fdate,f_datestr = fd[0]
    else:
        f_datestr_re = re.search("(\d{4}?\d{2}?\d{2})",p.stem)          
        if f_datestr_re:
            f_datestr = f_datestr_re.group()
            try:
                 fdate = dateutil.parser.parse(f_datestr)
            except:
                if p.is_dir:
                    fdate,f_datestr = get_date_and_string([j for j in p.glob('*')][0])
                else:
                    fdate,f_datestr = get_date_and_string(p.parent)
        else:
            if p.suffix in ['.lev10','.lev15','.lev20']: # special for AERONET
                f_datestr = re.search("(\d{2}?\d{2}?\d{2})",p.stem).group()
                fdate = dateutil.parser.parse(f_datestr,yearfirst=True)
    return fdate,f_datestr
    


# In[2]:


def pull_labels(p,datestr,filters={},dirlabel_in=None):
    'p:posis path to extract any labels, datestr: datestr within the filename'
    f_str = p.stem
    f_strs = f_str.replace(datestr,'').strip('_').split('_')
    f_strs_up = [tstr.upper() for tstr in f_strs]
    instnames_def = {'4STAR':['4STAR','4STARA','4-STAR','4STAR-A'],
                     '4STARB':['4STARB','4STAR2','4STAR-B'],
                     '3STAR':['3STAR','3-STAR'],
                     '2STAR':['2STAR','2-STAR'],
                     'muSSTAR':['MUSSTAR','MUSTAR','MU-STAR','MU-SSTAR'],
                     '5STAR':['5STARG','5STAR','5STAR-G','5STARF','5STAR-F'],
                     'AATS':['AATS','AATS14','AATS-14'],
                     'CAIR':['CAIR','C-AIR']}
    if not filters:
        instnames = instnames_def
    else:
        instnames = filters.get('instrument_names',instnames_def)
    instname,label,dirlabel = None, None, None
    for inm in instnames:
        for fstar in instnames[inm]:
            if fstar in f_strs_up:
                instname = inm
                null = f_strs.pop(f_strs_up.index(fstar))
                f_strs_up.remove(fstar)  
    if p.is_dir():
        if dirlabel_in: f_strs.insert(0,dirlabel_in)
        dirlabel = '_'.join(f_strs)
        dirlabel = dirlabel.strip('_')
    else:
        label = '_'.join(f_strs)
        label = label.strip('_')
    
    return dirlabel,label,instname


# In[3]:


def get_season(now):
    'Function to return a string with the season'
    Y = 2000 # dummy leap year to allow input X-02-29 (leap day)
    seasons = [('Winter', (date(Y,  1,  1),  date(Y,  3, 20))),
               ('Spring', (date(Y,  3, 21),  date(Y,  6, 20))),
               ('Summer', (date(Y,  6, 21),  date(Y,  9, 22))),
               ('Fall', (date(Y,  9, 23),  date(Y, 12, 20))),
               ('Winter', (date(Y, 12, 21),  date(Y, 12, 31)))]
    
    if isinstance(now, datetime):
        now = now.date()
    now = now.replace(year=Y)
    return next(season for season, (start, end) in seasons
                if start <= now <= end)


# In[4]:


class filetypes:
    'Class to identify the filetype, labels, dates, instrumentname, suffixes, input is string indicating full file path'
    def __init__(self,f,dirlabel=None,filters={}):
        p = Path(f)
        self.fdate,self.f_datestr = get_date_and_string(p)
        self.dirlabel,self.label,self.instname = pull_labels(p,self.f_datestr,filters=filters,dirlabel_in=dirlabel)
        self.p = p
        if not self.dirlabel:
            self.dirlabel = dirlabel
        self.daystr = self.fdate.strftime('%Y%m%d')
        
    def _print(self):
        print self.instname, self.dirlabel, self.label, self.fdate, self.f_datestr, self.p.stem
    
    def __getitem__(self,i):
        'Method to call only the variables in the class like a dict'
        return self.__dict__.get(i)

    def keys(self):
        'Method to wrap the dict call to the class object'
        return self.__dict__.keys()


# In[20]:


def get_newfilepath(f,filters={},debug=False):
    'function to build the end file path, input of filetype class f, outputs updated file class f'
    # determine the campaign
    campaign = 'rooftop'
    for campaign_name,date_range in filters['time_filter'].items():
        if (f.fdate >= date_range[0]) & (f.fdate <= date_range[1]):
            campaign = campaign_name
        if f.instname=='4STAR':
            for campaign_name,date_range in filters['time_filter_4STAR'].items():
                if (f.fdate >= date_range[0]) & (f.fdate <= date_range[1]):
                    campaign = campaign_name
        elif f.instname=='4STARB':
            for campaign_name,date_range in filters['time_filter_4STARB'].items():
                if (f.fdate >= date_range[0]) & (f.fdate <= date_range[1]):
                    campaign = campaign_name
        elif f.instname=='AATS':
            for campaign_name,date_range in filters['time_filter_AATS'].items():
                if (f.fdate >= date_range[0]) & (f.fdate <= date_range[1]):
                    campaign = campaign_name
    f.campaign = campaign

    if f.campaign=='rooftop':
        f.season = get_season(f.fdate)
        f.year = f.fdate.year
        f.campaign = os.path.join('rooftop','{season}_{year}'.format(**f))
    
    if debug: print str(f.p)
    if debug: print 'Campaign found to be:', campaign
    
    folders_match_filetype = [ni for ni in filters['directories']                              if (f.p.suffix.lower() in filters['directories'][ni]['filetypes'])]
    folders_match_label = [j for j in folders_match_filetype if                        any([lbl in f.label.lower() for lbl in filters['directories'][j]['label']]) &                        (not any([lbl in f.label.lower() for lbl in filters['directories'][j].get('not_label',[])]))]
    if debug: print 'folders_match_filetype:',folders_match_filetype
    if debug: print 'folders_match_label',folders_match_label
    if len(folders_match_label) == 0:
        folders_match_label = [j for j in folders_match_filetype if                                (not filters['directories'][j]['label'])]
        if len(folders_match_label) == 0:
            if verbose: print '*** Match move directory not found for file {p.stem}, using base path ***'.format(**f)
            folders_match_label = ['']

    f.newpath = Path(root_folder).joinpath('{campaign}'.format(**f),folders_match_label[0],                                 filters['directories'].get(folders_match_label[0],{}).get('folder_name','').format(**f))    
    f.newfile = f.newpath.joinpath(f.p.name)
    
    if debug: print 'newpath:',str(f.newpath),' newfile:',str(f.newfile)
        
    #check if destination file already exists:
    if f.newfile.exists():
        if filecmp.cmp(str(f.newfile),str(f.p),shallow=True):
            if filecmp.cmp(str(f.newfile),str(f.p),shallow=False):
                # they are the same and don't do anything
                if verbose: print '{prefix}Exact same file already exists at: {newfile}, removing incoming file'.format(**f)
                if not dry_run: os.remove(str(f.p))
                return None
        if verbose: print '{prefix}Different file with same name ({p.name}) exists'.format(**f)
        f.newpath = f.newpath.joinpath('Uploaded_on_{}'.format(date.today().strftime('%Y%m%d')))
        f.newfile = f.newpath.joinpath(f.p.name)
        
    return folders_match_label


# In[6]:


def get_filters_from_json(in_directory):
    'function to read in the filters from the json file'
    with open(in_directory+'.filters.json') as fjson: 
        filters = json.load(fjson)

    # sanitize input
    # set the dates
    for nt,lt in filters['time_filter'].items():
        lt[0] = datetime(lt[0][0],lt[0][1],lt[0][2])
        lt[1] = datetime(lt[1][0],lt[1][1],lt[1][2])
    for nt,lt in filters['time_filter_4STAR'].items():
        lt[0] = datetime(lt[0][0],lt[0][1],lt[0][2])
        lt[1] = datetime(lt[1][0],lt[1][1],lt[1][2])
    for nt,lt in filters['time_filter_4STARB'].items():
        lt[0] = datetime(lt[0][0],lt[0][1],lt[0][2])
        lt[1] = datetime(lt[1][0],lt[1][1],lt[1][2])
    for nt,lt in filters['time_filter_AATS'].items():
        lt[0] = datetime(lt[0][0],lt[0][1],lt[0][2])
        lt[1] = datetime(lt[1][0],lt[1][1],lt[1][2])

    # ensure capitalization of the instrument names 
    for ni,li in filters['instrument_names'].items():
        filters['instrument_names'][ni] = [lis.upper() for lis in li]

    # ensure lower of the filetypes and labels
    for ni in filters['directories']:
        filters['directories'][ni]['filetypes'] = [lis.lower() for lis in filters['directories'][ni]['filetypes']]
        filters['directories'][ni]['label'] = [lis.lower() for lis in filters['directories'][ni]['label']]
        if type(filters['directories'][ni]['folder_name']) is list:
            filters['directories'][ni]['folder_name'] = ''
    return filters


# In[7]:


def recurse_through_dir(indir,dirlabel=None,verbose=False,filters={}):
    fl = os.listdir(indir)
    fl_array = []
    dirs = {}
    for f in fl:
        if f.startswith('.'): continue #ignore hidden files
        fl_array.append(filetypes(indir+'/'+f,dirlabel=dirlabel,filters=filters))
        if verbose: fl_array[-1]._print()
        if fl_array[-1].p.is_dir():
            dirs[f] = fl_array[-1].dirlabel
            fla = recurse_through_dir(indir+'/'+f,dirlabel=dirs[f],verbose=verbose,filters=filters)
            fl_array.extend(fla)
    return fl_array


# # Prepare the command line argument parser

# In[8]:


import argparse


# In[9]:


long_description = """    Run the incoming file parser and moves the files to the desired sunsat locations
    The File locations and folders are defined by the json file: .filters.json
    Please update the date ranges within that file for any new field mission, if not then assumes rooftop measurements for the season_year
    Can run a call to matlab for any incoming 4STAR raw data"""


# In[10]:


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


# In[11]:


in_ = vars(parser.parse_known_args()[0])


# # Load the modules and get the defaults

# In[12]:


import os, zipfile
import dateutil.parser
import re
from pathlib2 import Path
from datefinder import find_dates
from datetime import date, datetime
import json
import filecmp
import subprocess


# In[13]:


in_directory = in_.get('in_dir','/data/sunsat/_incoming_gdrive/')
root_folder = in_.get('root_dir','/data/sunsat/')


# In[14]:


verbose = not in_.get('quiet',False)
dry_run = in_.get('dry_run',True)
run_matlab = in_.get('run_matlab',False)


# In[15]:


if verbose: print in_


# In[16]:


# Go through and unzip any folder
prefix = '*DRY RUN*: ' if dry_run else ''
for item in os.listdir(in_directory): # loop through items in dir
    if item.lower().endswith('.zip'): # check for ".zip" extension
        file_name = Path(in_directory+item) # get full path of files
        zip_ref = zipfile.ZipFile(str(file_name)) # create zipfile object
        if verbose: print '{prefix}found zip file: {file_name}, extracting here.'.format(prefix=prefix,file_name=file_name)
        if not dry_run: zip_ref.extractall(in_directory) # extract file to dir
        zip_ref.close() # close file
        if not dry_run: os.remove(str(file_name)) # delete zipped file


# In[23]:


filters = get_filters_from_json(in_directory)
fl_arr = recurse_through_dir(in_directory,verbose=verbose,filters=filters)


# In[24]:


data_raw_found = False
data_raw_files = {}
for f in fl_arr:
    # determine the end path of the file
    if f.p.is_dir(): continue # do nothing for directories
    f.prefix = '*DRY RUN*: ' if dry_run else ''
    
    folders_match_label = get_newfilepath(f,filters=filters,debug=False)
    if not folders_match_label: continue
        
    if 'data_raw' in folders_match_label:
        data_raw_found = True

    # now move the files
    if not f.newpath.exists() & verbose: print '{prefix}+Creating new path: {newpath}'.format(**f) 
    if not dry_run: f.newpath.mkdir(parents=True,exist_ok=True)
    if verbose: print '{prefix}~Moving file from {p}\n   to new path: {newfile}'.format(**f) 
    if not dry_run: f.p.rename(f.newfile)
    if 'data_raw' in folders_match_label: 
        data_raw_files['{instname}_{daystr}'.format(**f)] =                      data_raw_files.get('{instname}_{daystr}'.format(**f),[])
        data_raw_files['{instname}_{daystr}'.format(**f)].append(str(f.newfile))


# In[42]:


# clean up folders after move
for dirpath, dirnames, filenames in os.walk(in_directory,topdown=False):
    if not dirpath in in_directory:
        try: 
            if verbose: print '{pre}-removing :{path}'.format(pre=prefix,path=dirpath)
            if not dry_run:
                os.rmdir(dirpath) 
        except: 
            pass


# In[44]:


if run_matlab:
    prefix = '*DRY RUN*: ' if dry_run else ''
    for dr,drs in data_raw_files.items():
        # get the position of the new star.mat and starsun.mat files
        f = filetypes('{}star.mat'.format(dr),filters=filters)
        fml = get_newfilepath(f,filters=filters)
        if not dry_run: f.newpath.mkdir(parents=True,exist_ok=True)
        fs = filetypes('{}starsun.mat'.format(dr),filters=filters)
        fmls = get_newfilepath(fs,filters=filters)
        if not dry_run: fs.newpath.mkdir(parents=True,exist_ok=True)
            
        # make the position of the new quicklook file
        instname,daystr = dr.split('_')
        fq = filetypes('{daystr}_{instname}_Quicklooks.ppt'.format(instname=instname,daystr=daystr))
        fmlq = get_newfilepath(fq,filters=filters)
        
        # make the position of the new figure files
        ff = filetypes('{daystr}_{instname}_plots.fig'.format(instname=instname,daystr=daystr))
        fmlf = get_newfilepath(ff,filters=filters)
        
        # make the position of the aeronet files
        fa = filetypes('{daystr}_AERONET_NASA_Ames.lev15'.format(instname=instname,daystr=daystr))
        fmla = get_newfilepath(fa,filters=filters)
        
        # make the position of the gas_summary files
        fg = filetypes('{instname}_{daystr}_gas_summary.mat'.format(instname=instname,daystr=daystr))
        fmlg = get_newfilepath(fg,filters=filters)
        
        # make a string of the raw files    
        filelist = "'"+"';'".join(drs)+"'"
        if not f.instname in ['4STAR','4STARB']: # only for 4STARs for now.
            continue
        command_setpath = "setnamedpath('starimg','{fig_path}');setnamedpath('starsun','{aero_path}');setnamedpath('gas_summary','{gas_path}')".                           format(fig_path=str(f.newpath)+'/',aero_path=str(fa.newpath)+'/',gas_path=str(fg.newpath)+'/')
        command = "allstarmat({{{filelist}}},'{starmat}')".format(filelist=filelist,starmat=str(f.newfile))
        command2 = "starsun('{starmat}','{starsun}')".format(starmat=f.newfile,starsun=fs.newfile)
        command3 = "Quicklooks_4STAR('{starsun}','{starmat}','{quicklooks}')".format(starmat=str(f.newfile),
                                       starsun=str(fs.newfile),quicklooks=in_directory+str(fq.newfile.name))
        if verbose: 
            print ' '.join(['{}matlab'.format(prefix),'-nodisplay','-r','"{};{};{};{};exit;"'.                            format(command_setpath,command,command2,command3)])
        if not dry_run:
            process = subprocess.Popen(['matlab','-nodisplay','-r','"{};{};{};{};exit;"'.                                        format(command_setpath,command,command2,command3)],
                                       shell=True,
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)
            while process.returncode is None:
                # handle output by direct access to stdout and stderr
                for line in process.stdout:
                    if verbose: print line
                for liner in process.stderr:
                    if verbose: print liner
                # set returncode if the process has exited
                process.poll()
            #stdout, stderr = process.communicate()
            #if verbose: print stdout, stderr

