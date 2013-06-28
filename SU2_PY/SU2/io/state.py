__all__ = ['State']

import os, sys, shutil, copy
from ..io   import expand_part, get_adjointSuffix, add_suffix, \
                   get_specialCases
from ..util import bunch
from ..util import ordered_bunch


class State(ordered_bunch):
    
    def __init__(self,*args,**kwarg):
        
        super(State,self).__init__(*args,**kwarg)
        
        if not args and not kwarg:
            for key in ['FUNCTIONS','GRADIENTS','VARIABLES','FILES','HISTORY']:
                if not self.has_key(key):
                    self[key] = ordered_bunch()
    
    def update(self,ztate):
        assert isinstance(ztate,State) , 'must update with another State-type'
        for key in self.keys():
            if isinstance(ztate[key],dict):
                self[key].update( ztate[key] )
            elif ztate[key]:
                self[key] = ztate[key]
        
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        output = 'STATE:'
        for k1,v1 in self.iteritems():
            output += '\n    %s:' % k1
            if isinstance(v1,dict):
                for k2,v2 in v1.iteritems():
                    output += '\n        %s: %s' % (k2,v2)
            elif v1:
                output += '\n        %s' % v1
        return output
    
    def pullnlink(self,config):
        
        pull = []; link = []
        
        # choose files to pull and link
        for key,value in self.FILES.iteritems():
            
            # link big files
            if key == 'MESH':
                value = expand_part(value,config)
                link.extend(value)
            elif key == 'DIRECT':
                link.append(value)
            elif 'ADJOINT_' in key:
                link.append(value)
            
            # copy all other files
            else:
                pull.append(value)
        
        #: for each filename
        
        return pull,link
    
    def design_vector(self):
        """ vectorizes State.VARIABLES
        """
        vector = []
        for value in self.VARIABLES.values():
            if isinstance(value,dict):
                for v in value.values():
                    vector.append(v)
            elif not isinstance(value,list):
                value = [value]
            vector.extend(value)
        return vector
    
    def find_files(self,config):
        """ State.find_files(config)
            finds mesh and solution files for a given config
            files already logged in state are not overridden
        """
        
        files = self.FILES
        
        mesh_name     = config.MESH_FILENAME
        direct_name   = config.SOLUTION_FLOW_FILENAME
        adjoint_name  = config.SOLUTION_ADJ_FILENAME
        targetea_name = 'TargetEA.dat'
        
        adj_map = get_adjointSuffix()
        
        restart = config.RESTART_SOL == 'YES'
        special_cases = get_specialCases(config)
        
        def register_file(label,filename):
            if not files.has_key(label):
                if os.path.exists(filename):
                    files[label] = filename
                    print 'found: %s' % filename
            else:
                assert os.path.exists(files[label]) , 'state expected file: %s' % filename
        #: register_file()                

        # mesh
        register_file('MESH',mesh_name)
        
        # direct solution
        if restart:
            register_file('DIRECT',direct_name)
        
        # adjoint solutions
        if restart:
            for obj,suff in adj_map.iteritems():
                ADJ_LABEL = 'ADJOINT_' + obj
                adjoint_name_suffixed = add_suffix(adjoint_name,suff)
                register_file(ADJ_LABEL,adjoint_name_suffixed)
        
        # equivalent area
        if 'EQUIV_AREA' in special_cases:
            register_file('TARGET_EA',targetea_name)
        
        return
    
            
            
    

def default_state(state):
    if not state: state = State()
    assert isinstance(state,State) , 'not a state instance'
    return state