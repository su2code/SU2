## \file functions.py
#  \brief python package for functions
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

import os, sys, shutil, copy
from .. import run  as su2run
from .. import io   as su2io
from .. import util as su2util
from ..io import redirect_folder, redirect_output


# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------

def function( func_name, config, state=None ):
    """ val = SU2.eval.func(func_name,config,state=None)
    
        Evaluates the aerodynamics and geometry functions.
        
        Wraps:
            SU2.eval.aerodynamics()
            SU2.eval.geometry()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh need not be deformed.
            Updates config and state by reference.
            Redundancy if state.FUNCTIONS is not empty.
        
        Executes in:
            ./DIRECT or ./GEOMETRY
        
        Inputs:
            func_name - SU2 objective function name or 'ALL'
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            If func_name is 'ALL', returns a Bunch() of 
            functions with keys of objective function names
            and values of objective function floats.
            Otherwise returns a float.
    """
    
    # initialize
    state = su2io.State(state)
    
    # redundancy check
    if not state['FUNCTIONS'].has_key(func_name):
        
        # Aerodynamics
        if func_name == "ALL" or func_name in su2io.optnames_aero:
            aerodynamics( config, state )
        
        # Geometry
        elif func_name in su2io.optnames_geo:
            geometry( func_name, config, state )
            
        else:
            raise Exception, 'unknown function name'
        
    #: if not redundant
    
    # prepare output
    if func_name == 'ALL':
        func_out = state['FUNCTIONS']
    else:
        func_out = state['FUNCTIONS'][func_name]
    
    return copy.deepcopy(func_out)

#: def function()


# ----------------------------------------------------------------------
#  Aerodynamic Functions
# ----------------------------------------------------------------------

def aerodynamics( config, state=None ):
    """ vals = SU2.eval.aerodynamics(config,state=None)
    
        Evaluates aerodynamics with the following:
            SU2.run.decompose()
	    SU2.run.deform()
            SU2.run.direct()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Redundancy if state.FUNCTIONS is not empty.
            
        Executes in:
            ./DIRECT
            
        Inputs:
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            Bunch() of functions with keys of objective function names
            and values of objective function floats.
    """
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_direct = 'log_Direct.out'
    else:
        log_direct = None
    
    # ----------------------------------------------------    
    #  Update Mesh
    # ----------------------------------------------------
    
    # does decomposition and deformation
    info = update_mesh(config,state)
                
    # ----------------------------------------------------    
    #  Adaptation (not implemented)
    # ----------------------------------------------------
    
    #if not state.['ADAPTED_FUNC']:
    #    config = su2run.adaptation(config)
    #    state['ADAPTED_FUNC'] = True
    
    # ----------------------------------------------------    
    #  Direct Solution
    # ----------------------------------------------------    
    
    # redundancy check
    direct_done = all( [ state.FUNCTIONS.has_key(key) for key in su2io.optnames_aero ] )
    if direct_done:
        # return aerodynamic function values
        aero = su2util.ordered_bunch()
        for key in su2io.optnames_aero:
            aero[key] = state.FUNCTIONS[key]
        return copy.deepcopy(aero)    
    #: if redundant
    
    # files to pull
    files = state.FILES
    pull = []; link = []
    
    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,config)
    link.extend(name)
    
    # files: direct solution
    if files.has_key('DIRECT'):
        link.append( files['DIRECT'] )
        ##config['RESTART_SOL'] = 'YES' # don't override config file
    else:
        config['RESTART_SOL'] = 'NO'
        
    # files: target equivarea distribution
    if ( 'EQUIV_AREA' in special_cases and 
         'TARGET_EA' in files ) : 
        pull.append( files['TARGET_EA'] )
    
    # output redirection      
    with redirect_folder( 'DIRECT', pull, link ) as push:
        with redirect_output(log_direct):     
            
            # # RUN DIRECT SOLUTION # #
            info = su2run.direct(config)
            su2io.restart2solution(config,info)
            state.update(info)
            
            # files to push
            for this_type in ['DIRECT','WEIGHT_NF']:
                if info['FILES'].has_key(this_type):
                    push.append( info['FILES'][this_type] )
            
    #: with output redirection
    
    # return aerodynamic dictionary
    aero = state['FUNCTIONS'] 
    for key in aero.keys(): 
        if key not in su2io.optnames_aero: del aero[key]
        
    return copy.deepcopy(aero)

#: def aerodynamics()


# ----------------------------------------------------------------------
#  Geometric Functions
# ----------------------------------------------------------------------

def geometry( func_name, config, state=None ):
    """ val = SU2.eval.geometry(config,state=None)
    
        Evaluates geometry with the following:
            SU2.run.decompose()
	    SU2.run.deform()
            SU2.run.geometry()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Redundancy if state.FUNCTIONS does not have func_name.
            
        Executes in:
            ./GEOMETRY
            
        Inputs:
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            Bunch() of functions with keys of objective function names
            and values of objective function floats.
    """
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_geom = 'log_Geometry.out'
    else:
        log_geom = None
    
    # ----------------------------------------------------    
    #  Update Mesh
    # ----------------------------------------------------
    
    # does decomposition and deformation
    info = update_mesh(config,state)    
    
    
    # ----------------------------------------------------    
    #  Geometry Solution
    # ----------------------------------------------------    
    
    # redundancy check
    geometry_done = state.FUNCTIONS.has_key(func_name)
    #geometry_done = all( [ state.FUNCTIONS.has_key(key) for key in su2io.optnames_geo ] )
    if not geometry_done:    
        
        # files to pull
        files = state.FILES
        pull = []; link = []
        
        # files: mesh
        name = files['MESH']
        name = su2io.expand_part(name,config)
        link.extend(name)
        
        # update function name
        ## TODO
        
        # output redirection
        with redirect_folder( 'GEOMETRY', pull, link ) as push:
            with redirect_output(log_geom):     
                
                # setup config
                config.GEO_PARAM = func_name
                
                # # RUN GEOMETRY SOLUTION # #
                info = su2run.geometry(config)
                state.update(info)
                
                # no files to push
                
        #: with output redirection
        
    #: if not redundant 
    
    # return geometry dictionary
    geom = state['FUNCTIONS'] # naive update
    for key in geom.keys(): 
        if key not in su2io.optnames_geo: del geom[key]   
        
    return copy.deepcopy(geom)
    

#: def geometry()



def update_mesh(config,state=None):
    """ SU2.eval.update_mesh(config,state=None)
    
        updates mesh with the following:
            SU2.run.decompose()
	    SU2.run.deform()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            
        Executes in:
            ./DECOMP and ./DEFORM
            
        Inputs:
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            nothing
            
        Modifies:
            config and state by reference
    """
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_decomp = 'log_Decomp.out'
        log_deform = 'log_Deform.out'
    else:
        log_decomp = None
        log_deform = None
        
        
    # ----------------------------------------------------
    #  Decomposition    
    # ----------------------------------------------------
    
    # redundancy check
    if not config.get('DECOMPOSED',False):
        
        # files to pull
        pull = []
        link = config['MESH_FILENAME']
        
        # output redirection
        with redirect_folder('DECOMP',pull,link) as push:    
            with redirect_output(log_decomp):
            
                # # RUN DECOMPOSITION # # 
                info = su2run.decompose(config)
                state.update(info)
                              
                # files to push
                if info.FILES.has_key('MESH'):
                    meshname = info.FILES.MESH
                    names = su2io.expand_part( meshname , config )
                    push.extend( names )
                #: if push
        
        #: with output redirection
        
    #: if not redundant
    
    
    # ----------------------------------------------------
    #  Deformation
    # ----------------------------------------------------
    
    # redundancy check
    deform_set  = config['DV_KIND'] == config['DEFINITION_DV']['KIND']
    deform_todo = not config['DV_VALUE_NEW'] == config['DV_VALUE_OLD']
    if deform_set and deform_todo:
    
        # files to pull
        pull = []
        link = config['MESH_FILENAME']
        link = su2io.expand_part(link,config)
        
        # output redirection
        with redirect_folder('DEFORM',pull,link) as push:
            with redirect_output(log_deform):
                
                # # RUN DEFORMATION # #
                info = su2run.deform(config)
                state.update(info)
                
                # data to push
                meshname = info.FILES.MESH
                names = su2io.expand_part( meshname , config )
                push.extend( names )
        
        #: with redirect output
        
    elif deform_set and not deform_todo:
        state.VARIABLES.DV_VALUE_NEW = config.DV_VALUE_NEW

    #: if not redundant

    return 

