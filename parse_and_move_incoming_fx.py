#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Fucntions to support parse_and_move_incoming.py:
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
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2021-01-26
#              - pulled out the functions from parse_and_move_incoming into this standalone
# 

# # Set up the background functions

# In[1]:


from __future__ import print_function
import dateutil.parser
import re, os
from pathlib2 import Path
from datefinder import find_dates
from datetime import date, datetime
import json
import filecmp


# In[2]:


def get_date_and_string(p,dirdate=None):
    'p: posix path to extract the date and the string'
    fdate = dirdate
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
    
    if not fdate: fdate = datetime.fromtimestamp(p.stat().st_ctime)
    if not locals().get('f_datestr'): f_datestr = fdate.strftime('%Y%m%d')
        
    return fdate,f_datestr


# In[3]:


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


# In[4]:


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


# In[5]:


class filetypes:
    'Class to identify the filetype, labels, dates, instrumentname, suffixes, input is string indicating full file path'
    def __init__(self,f,dirlabel=None,dirdate=None,filters={}):
        p = Path(f)
        self.fdate,self.f_datestr = get_date_and_string(p,dirdate=dirdate)
        self.dirlabel,self.label,self.instname = pull_labels(p,self.f_datestr,filters=filters,dirlabel_in=dirlabel)
        self.p = p
        if not self.dirlabel:
            self.dirlabel = dirlabel
        if not dirdate:
            self.dirdate = self.fdate
        else:
            self.dirdate = dirdate
        self.daystr = self.dirdate.strftime('%Y%m%d')
        
    def _print(self):
        print(self.instname, self.dirlabel, self.label, self.fdate, self.f_datestr, self.p.stem)
    
    def __getitem__(self,i):
        'Method to call only the variables in the class like a dict'
        return self.__dict__.get(i)

    def keys(self):
        'Method to wrap the dict call to the class object'
        return self.__dict__.keys()


# In[6]:


def get_newfilepath(f,filters={},debug=False,fake_file=False,nexact=0,root_folder='/nobackupp10/nexprojects/sunsat/'):
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
    
    if debug: print(str(f.p))
    if debug: print('Campaign found to be:', campaign)
    
    folders_match_filetype = [ni for ni in filters['directories']                              if (f.p.suffix.lower() in filters['directories'][ni]['filetypes'])]
    folders_match_label = [j for j in folders_match_filetype if                        any([lbl in f.label.lower() for lbl in filters['directories'][j]['label']]) &                        (not any([lbl in f.label.lower() for lbl in filters['directories'][j].get('not_label',[])]))]
    if debug: print('folders_match_filetype:',folders_match_filetype)
    if debug: print('folders_match_label',folders_match_label)
    if len(folders_match_label) == 0:
        folders_match_label = [j for j in folders_match_filetype if                                (not filters['directories'][j]['label'])]
        if len(folders_match_label) == 0:
            if verbose: print('*** Match move directory not found for file {p.stem}, using base path ***'.format(**f))
            folders_match_label = ['']

    f.newpath = Path(root_folder).joinpath('{campaign}'.format(**f),folders_match_label[0],                                 filters['directories'].get(folders_match_label[0],{}).get('folder_name','').format(**f))   
    f.newfile = f.newpath.joinpath(f.p.name)
    
    if debug: print( 'newpath:',str(f.newpath),' newfile:',str(f.newfile))
        
    #check if destination file already exists:
    if f.newfile.exists() & (not fake_file):
        if filecmp.cmp(str(f.newfile),str(f.p),shallow=True):
            if filecmp.cmp(str(f.newfile),str(f.p),shallow=False):
                # they are the same and don't do anything
                if verbose: 
                    print( '{prefix}Exact same file already exists at: {newfile}, removing incoming file'.format(**f))
                if not dry_run: os.remove(str(f.p))
                nexact = nexact+1
                return None,nexact
        if verbose: print( '{prefix}Different file with same name ({p.name}) exists'.format(**f))
        f.newpath = f.newpath.joinpath('Uploaded_on_{}'.format(date.today().strftime('%Y%m%d')))
        f.newfile = f.newpath.joinpath(f.p.name)
        
    return folders_match_label,nexact


# In[7]:


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


# In[8]:


def recurse_through_dir(indir,dirlabel=None,dirdate=None,verbose=False,filters={}):
    fl = os.listdir(indir)
    fl_array = []
    dirs = {}
    for f in fl:
        if f.startswith('.'): continue #ignore hidden files
        if verbose: print(f+' -> ',end='')
        fl_array.append(filetypes(indir+'/'+f,dirlabel=dirlabel,filters=filters,dirdate=dirdate))
        if verbose: fl_array[-1]._print()
        if fl_array[-1].p.is_dir():
            dirs[f] = fl_array[-1].dirlabel
            fla = recurse_through_dir(indir+'/'+f,dirlabel=dirs[f],verbose=verbose,
                                      filters=filters,dirdate=fl_array[-1].fdate)
            fl_array.extend(fla)
    return fl_array


# In[9]:


def make_temp_mfile(mfilepath,filelist,starmat,starsun,quicklooks,fig_path,aero_path,gas_path,sun_path,incoming_path):
    'Print the matlab commands to a temp m file'
    
    command_setpath = ["setnamedpath('starimg','{fig_path}');".format(fig_path=fig_path),
                       "setnamedpath('starsun','{sun_path}');".format(sun_path=sun_path),
                       "setnamedpath('aeronet','{aero_path}');".format(aero_path=aero_path),
                       "setnamedpath('gas_summary','{gas_path}');".format(gas_path=gas_path)]
    command = "allstarmat({{{filelist}}},'{starmat}')".format(filelist=filelist,starmat=starmat)
    command2 = "starsun('{starmat}','{starsun}')".format(starmat=starmat,starsun=starsun)
    command3 = "Quicklooks_4STAR('{starsun}','{starmat}','{quicklooks}')".format(starmat=starmat,
                                   starsun=starsun,quicklooks=quicklooks)
    
    commands = command_setpath
    commands.append(command)
    commands.append(command2)
    commands.append("pa = getnamedpath('starmat');")
    commands.append("setnamedpath('starmat','{incoming_path}');".format(incoming_path=incoming_path))
    commands.append(command3)
    commands.append("setnamedpath('starmat',pa);")
    
    with open(mfilepath,'w') as f:
        f.write('% Temp matlab script created on : '+str(datetime.now())+' \n')
        for cm_line in commands:
            f.write(cm_line+' \n')
        f.write('exit;\n')
    return mfilepath


# In[10]:


def move_files(fl_arr,filters,verbose=False,dry_run=False):
    'Function to go through and move the files, create folders if necessary, and check if there are any new raw data'
    data_raw_found = False
    data_raw_files = {}
    nexact, nmoved, ncreated, ndataraw = 0,0,0,0
    for f in fl_arr:
        # determine the end path of the file
        if f.p.is_dir(): continue # do nothing for directories
        f.prefix = '*DRY RUN*: ' if dry_run else ''

        folders_match_label,nexact = get_newfilepath(f,filters=filters,debug=False,nexact=nexact,root_folder=root_folder)
        if not folders_match_label: continue

        if 'data_raw' in folders_match_label:
            data_raw_found = True

        # now move the files
        if not f.newpath.exists() & verbose: print( '{prefix}+Creating new path: {newpath}'.format(**f) )
        if not dry_run: 
            if not f.newpath.exists(): ncreated = ncreated+1
            f.newpath.mkdir(parents=True,exist_ok=True)
        if verbose: print( '{prefix}~Moving file from {p}\n   to new path: {newfile}'.format(**f) )
        if not dry_run: 
            nmoved = nmoved+1
            f.p.rename(f.newfile)
        if 'data_raw' in folders_match_label: 
            data_raw_files['{instname}_{daystr}'.format(**f)] =                          data_raw_files.get('{instname}_{daystr}'.format(**f),[])
            data_raw_files['{instname}_{daystr}'.format(**f)].append(str(f.newfile))
            ndataraw = ndataraw+1
    return data_raw_found, data_raw_files, nexact, nmoved, ncreated, ndataraw

