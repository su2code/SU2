## \file state.py
#  \brief python package for state 
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.6
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, time
from ..io   import expand_part, get_adjointSuffix, add_suffix, \
                   get_specialCases
from ..util import bunch
from ..util import ordered_bunch


# ----------------------------------------------------------------------
#  State Factory
# ----------------------------------------------------------------------

def State_Factory(state=None):
    """ state = SU2.io.State()
        
        Starts a state class, an extension of ordered_bunch().
        Stores data generated while traversing SU2 tool chain
        
        Parameters:
            FUNCTIONS - ordered bunch of objective function values
            GRADIENTS - ordered bunch of gradient value lists
            VARIABLES - ordered bunch of variables
            FILES     - ordered bunch of file types
            HISTORY   - ordered bunch of history information
            
        Parameters can be accessed by item or attribute
        ie: state['FUNCTIONS'] or state.FUNCTIONS
        
        Methods:
            update()        - updates self with another state
            pullnlink()     - returns files to pull and link
            design_vector() - vectorizes design variables
            find_files()    - finds existing mesh and solutions
        
        Example of a filled state:
        FUNCTIONS:
            LIFT: 0.2353065809
            DRAG: 0.042149736
            SIDEFORCE: 0.0
            MOMENT_X: 0.0
            MOMENT_Y: 0.0
            MOMENT_Z: 0.046370243
            FORCE_X: 0.0370065195
            FORCE_Y: 0.2361700759
            FORCE_Z: 0.0
            EFFICIENCY: 5.5826347517
        GRADIENTS:
            DRAG: [0.133697, 0.41473, 0.698497, (...)
        VARIABLES:
            DV_VALUE_NEW: [0.002, 0.002, 0.002, (...)
        FILES:
            MESH: mesh.su2
            DIRECT: solution_flow.dat
            ADJOINT_DRAG: solution_adj_cd.dat
        HISTORY:
            DIRECT: {ITERATION=[1.0, 2.0, 3.0, (...)
            ADJOINT_DRAG: {ITERATION=[1.0, 2.0, 3.0, (...)

    """   
    
    if not state is None:
        assert isinstance(state,State) , 'input is must be a state instance'
        return state
    
    NewClass = State()
    
    for key in ['FUNCTIONS','GRADIENTS','VARIABLES','FILES','HISTORY']:
        NewClass[key] = ordered_bunch()
            
    return NewClass


# ----------------------------------------------------------------------
#  State Class
# ----------------------------------------------------------------------

class State(ordered_bunch):
    """ state = SU2.io.state.State()
        
        This is the State class that should be generated with the 
        Factory Function SU2.io.state.State_Factory()
        
        Parameters:
            none, should be loaded with State_Factory()
        
        Methods:
            update()        - updates self with another state
            pullnlink()     - returns files to pull and link
            design_vector() - vectorizes design variables
            find_files()    - finds existing mesh and solutions

    """
    
    _timestamp = 0
    
    def update(self,ztate):
        """ Updates self given another state
        """
        
        if not ztate: return
        assert isinstance(ztate,State) , 'must update with another State-type'
        for key in self.keys():
            if isinstance(ztate[key],dict):
                self[key].update( ztate[key] )
            elif ztate[key]:
                self[key] = ztate[key]
        
        self.set_timestamp()
                
        
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        output = 'STATE:'
        for k1,v1 in self.iteritems():
            output += '\n    %s:' % k1
            if isinstance(v1,dict):
                for k2,v2 in v1.iteritems():
                    output += '\n        %s: %s' % (k2,v2)
            else:
                output += '\n        %s' % v1
        return output
    
    def pullnlink(self,config):
        """ pull,link = SU2.io.State.pullnlink(config)
            returns lists pull and link of files for folder
            redirection, based on a given config
        """
        
        pull = []; link = []
        
        # choose files to pull and link
        for key,value in self.FILES.iteritems():
            
            # link big files
            if key == 'MESH':
                # mesh (merged and partitions)
                if config.DECOMPOSED:
                    value = expand_part(value,config) # hack - twl
                else:
                    value = [value]
                link.extend(value)
            elif key == 'DIRECT':
                #if config.RESTART_SOL == 'YES':
                # direct solution
                link.append(value)
            elif 'ADJOINT_' in key:
                #if config.RESTART_SOL == 'YES':
                # adjoint solution
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
        """ SU2.io.State.find_files(config)
            finds mesh and solution files for a given config.
            updates state.FILES with filenames.
            files already logged in state are not overridden.
            will ignore solutions if config.RESTART_SOL == 'NO'.
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
    
    def __setitem__(self,k,v):
        if self._initialized:
            self.set_timestamp()
        super(State,self).__setitem__(k,v)
    
    def set_timestamp(self):
        self._timestamp = time.time()
    
    def tic(self):
        """ timestamp = State.tic() 
            returns the time that this state was last modified
        """
        return self._timestamp
    
    def toc(self,timestamp):
        """ updated = State.toc(timestamp)
            returns True if state was modified since last timestamp
        """
        return self._timestamp > timestamp
        
    
#: def State
