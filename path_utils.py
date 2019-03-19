
# coding: utf-8

# In[63]:


def __init__():
    """
    Purpose:  
        Create a method of handling paths, independent of the machine or os

    Input:
        none here

    Output:
        full path for defined locations

    Keywords:
        none

    Dependencies:
        - os

    Needed Files:
      - python_named_paths.json in home or appdata directory (if not found, will create it) 

    Modification History:
        Written: Samuel LeBlanc, Santa Cruz, 2017-11-15
        Modified: 
    
    """
    f = getpathfile() # run check for existence and if not create the json file
    pass


# In[92]:


def getpathfile():
    """
    Purpose:  
        return the path of the json file to open, platform independent
        creates the file if not existent

    Input:
        None

    Output:
        full path and file name

    Keywords:
        none

    Dependencies:
        - os, sys

    Needed Files:
      - python_named_paths.json in directory (if not found, will create it) 

    Modification History:
        Written: Samuel LeBlanc, Santa Cruz, 2017-11-16
        Modified: 
    
    """
    
    import os, sys
    import json
    name = 'python_named_paths.json'
    if sys.platform.startswith('win'):
        pa = os.path.join(os.environ['APPDATA'],'python')
        f = os.path.join(pa,name)          
    else: #for unixes
        pa = os.path.join(os.environ['HOME'],'.config')
        f = os.path.join(pa,name)
    
    if not os.path.isdir(pa):
        os.makedirs(pa)
    
    if not os.path.isfile(f) or not os.stat(f).st_size:           
        # build the default filepaths
        sep = os.path.sep
        default = {'home':os.environ['HOME']+sep,'~':os.environ['HOME']+sep}
        
        try:
            res = os.environ['RES']+sep
        except:
            res = os.path.join(os.environ['HOME'],'Research')+sep
        
        try:
            nob = os.environ['nob']+sep # for nobackup
        except:
            nob = res
        
        default['res'] = res
        default['research'] = res
        default['nob'] = nob
        default['nobackup'] = nob
        
        with open(f,'w') as outfile:
            json.dump(default,outfile)
    return f


# In[1]:


def listpath():
    """
    Purpose:  
        print out all named paths and their full path

    Input:
        None

    Output:
        None
        
    Keywords:
        None

    Dependencies:
        - os
        - json
        - warnings
        - this file

    Needed Files:
      - python_named_paths.json in .jupyter directory (if not found, will create it) 

    Modification History:
        Written: Samuel LeBlanc, Santa Cruz, 2018-09-05
        Modified: 
    
    """
    import os
    from path_utils import getpathfile, namepath
    import json
    from warnings import warn
    
    f = getpathfile()
    
    with open(f) as json_file:  
        p = json.load(json_file)
    for k in p.keys():
        print "'{}': {}".format(k,p[k])


# In[62]:


def getpath(s='research',make_path=False,path='',verbose=True):
    """
    Purpose:  
        return the path named in the json file, or if not there, an educated guess of the path

    Input:
        s: (string, default research) name of the path to return (string)

    Output:
        f: full path (for directory) for named path
        
    Keywords:
        make_path: (boolean, default False) if True, will name a relative path defined by the research path with s as a subfolder
                    ignored if s is already in path
        path: (string, default empty) if set to a path, and make_path is True, will overwrite default make_path behavior and 
                save a named path defined by path (understands if absolute path or relative to Research directory)
        verbose: (boolean, default True) if True, prints out the the name and the path used

    Dependencies:
        - os
        - json
        - warnings
        - this file

    Needed Files:
      - python_named_paths.json in .jupyter directory (if not found, will create it) 

    Modification History:
        Written: Samuel LeBlanc, Santa Cruz, 2017-11-15
        Modified: 
    
    """
    import os
    from path_utils import getpathfile, namepath
    import json
    from warnings import warn
    
    f = getpathfile()
    
    with open(f) as json_file:  
        p = json.load(json_file)
    
    if len(s)==0:
        return p['research']
    
    if make_path:
        if len(path)==0:
            p = namepath(s,os.path.join(p['research'],s),verbose=verbose)
        else:
            p = namepath(s,path,verbose=verbose)

    if p.get(s.lower()):
        if verbose: print 'Return path named: ',s,p.get(s.lower())
        return p.get(s.lower())
    else:
        if os.path.isdir(os.path.join(p['research'],s)):
            return os.path.isdir(os.path.join(p['research'],s))
        else:
            warn(os.path.join(p['research'],s)+'directory does not exists, returning default research directory')
            return p['research']        
                


# In[90]:


def namepath(name,path,verbose=False):
    """
    Purpose:  
        Name a path to be saved in the json path file

    Input:
        name: Name of the named path, should be similar to actual path
        path: Path (either relative or absolute) of the named path

    Output:
        New path dict

    Keywords:
        verbose: (boolean,default False) if True prints out some info

    Dependencies:
        - os
        - json
        - path_utils (this file)
        - warnings

    Needed Files:
      - python_named_paths.json in home or APPDATA directory (if not found, will create it) 

    Modification History:
        Written: Samuel LeBlanc, Santa Cruz, 2017-11-15
        Modified: 
    
    """
    import os
    from path_utils import getpathfile
    import json
    from warnings import warn
    
    f = getpathfile()
    
    with open(f) as json_file:  
        p = json.load(json_file)
    
    if len(name)==0:
        raise ValueError('name is not set')
    if len(path)==0:
        raise ValueError('path is not set')
    
    if os.path.isabs(path):
        p[name.lower()] = path+os.path.sep
    else:
        p[name.lower()] = os.path.join(p['research'],path)+os.path.sep
    
    if verbose: print 'Created new path: {} under the name: {}'.format(p[name.lower()],name)
    if not os.path.isdir(p[name.lower()]):
        warn('New named path {} does not exist'.format(p[name.lower()]))
    
    with open(f,'w') as json_file:
        json.dump(p,json_file)
    
    return p

