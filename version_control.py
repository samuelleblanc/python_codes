# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#global program_version

# <codecell>

def version_set(v_tag):
    """
    Name:

        version_set
    
    Purpose:

        to set the version of the file currently being called
        Returns a tuple of dictionary with items descripbing 
        the current file name, last time it was modified and by whom, and the v_tag (version tag, eg., v1.3...)
    
    Calling Sequence:

        program_version = version_set(v_tag,program_version)        
    
    Input: 
  
        v_tag: version tag in string (e.g., 'v1.0')
        program_version: tuple of dictionary containing info from all runned programs
    
    Output:

        program_vervion: see above
    
    Keywords: 

        none
    
    Dependencies:

        - os
        - time
        - inspect
        - getpass
    
    Required files:
   
        None
    
    Example:

        ...
        
    Modification History:
    
        Written (v1.0): Samuel LeBlanc, 2015-02-02, NASA Ames
        
    """
    #global program_version
    
    import os, time
    import inspect
    import getpass
    
    def getname(): 
        import IPython
        IPython.display.display(IPython.display.Javascript('IPython.notebook.kernel.execute("theNotebookName = " + "\'"+IPython.notebook.notebook_name+"\'");'))
        IPython.display.display(IPython.display.Javascript('IPython.notebook.kernel.execute("theNotebookPath = " + "\'"+IPython.notebook.notebook_path+"\'");'))
    getname()
  #  thisfilepath = os.getcwd() + os.path.sep+theNotebookPath+theNotebookName
    
    current_file_name = inspect.stack()[1][1]
    #current_file_name = theNotebookName
    print current_file_name
    (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(current_file_name)
    lm = "last modified: %s" % time.ctime(mtime)
    lc = "file created: %s" % time.ctime(ctime)
    user = getpass.getuser()
    
    pro_version = {'Filename':current_file_name,'Version':v_tag,'creation':lc,'modified':lm,'Run by':user}
    #return program_version+pro_version
    return pro_version

