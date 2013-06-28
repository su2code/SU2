__all__ = ['output','folder']

import os, sys, shutil, copy, glob
from .tools import add_suffix, make_link, expand_part

# -------------------------------------------------------------------
#  Output Redirection 
# -------------------------------------------------------------------
# original source: http://stackoverflow.com/questions/6796492/python-temporarily-redirect-stdout-stderr
class output(object):
    ''' Temporarily redirects sys.stdout and sys.stderr when used in
        a 'with' contextmanager
    '''
    def __init__(self, stdout=None, stderr=None):
        
        _newout = False
        _newerr = False
        
        if isinstance(stdout,str):
            stdout = open(stdout,'a')
            _newout = True            
        if isinstance(stderr,str):
            stderr = open(stderr,'a')
            _newerr = True                   
                
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr
        self._newout = _newout
        self._newerr = _newerr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
        
        if self._newout:
            self._stdout.close()
        if self._newerr:
            self._stderr.close()           

#: class output()


# -------------------------------------------------------------------
#  Folder Redirection 
# -------------------------------------------------------------------
class folder(object):
    ''' Temporarily redirects working folder 
        accepts wild cards
        will overwrite pushed files
    '''
    def __init__(self, folder, pull=[], link=[], force=True ):
        ''' folder redirection initialization
            see help( folder ) for more info
        '''
        
        if not isinstance(pull,list) : pull = [pull]
        if not isinstance(link,list) : link = [link]
        
        origin = os.getcwd()
        origin = os.path.abspath(origin).rstrip('/')+'/'
        folder = os.path.abspath(folder).rstrip('/')+'/'
        
        self.origin = origin
        self.folder = folder
        self.pull   = copy.deepcopy(pull)
        self.push   = []
        self.link   = copy.deepcopy(link)
        self.force  = force

    def __enter__(self): 
        
        origin = self.origin  # absolute path
        folder = self.folder  # absolute path
        pull   = self.pull
        push   = self.push
        link   = self.link
        force  = self.force
        
        # check for no folder change
        if folder == origin:
            return []
        
        # relative folder path
        #relative = os.path.relpath(folder,origin)
        
        # check, make folder
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        # copy pull files
        for name in pull:
            new_name = os.path.split(name)[-1]
            new_name = os.path.join(folder,new_name)
            name = os.path.abspath(name)
            if name == new_name: continue
            if os.path.exists( new_name ): 
                if force:
                    os.remove( new_name )
                    shutil.copy(name,new_name)
            else:
                shutil.copy(name,new_name)

        # make links
        for name in link:
            new_name = os.path.split(name)[-1]
            new_name = os.path.join(folder,new_name)
            name = os.path.abspath(name)
            if name == new_name: continue
            if os.path.exists( new_name ): 
                if force:
                    os.remove( new_name )
                    make_link(name,new_name)
            else:
                make_link(name,new_name)
            
                
        # change directory
        os.chdir(folder)
        
        # return empty list to append with files to push to super folder
        return push

    def __exit__(self, exc_type, exc_value, traceback):
        
        origin = self.origin
        folder = self.folder
        push   = self.push
        force  = self.force
        
        # check for no folder change
        if folder == origin:
            return
        
        # move assets
        for name in push:
            new_name = os.path.join(origin,name)
            name = os.path.abspath(name)
            if name == new_name: continue
            if os.path.exists( new_name ): 
                if force:
                    os.remove( new_name )
                    shutil.move(name,new_name)
            else:
                shutil.move(name,new_name)
            
        
        # change directory
        os.chdir(origin)
        
#: class folder()
