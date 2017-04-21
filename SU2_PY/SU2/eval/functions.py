#!/usr/bin/env python

## \file functions.py
#  \brief python package for functions
#  \author T. Lukaczyk, F. Palacios
#  \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, time
from .. import run  as su2run
from .. import io   as su2io
from .. import util as su2util
from ..io import redirect_folder, redirect_output


# ----------------------------------------------------------------------
#  Main Function Interface
# ----------------------------------------------------------------------

def function( func_name, problem, state=None ):
    """ val = SU2.eval.func(func_name,opt,state=None)
    
        Evaluates the objective function.
        
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
    
    # check for multiple objectives
    multi_objective = (type(func_name)==list)
    # func_name_string is only used to check whether the function has already been evaluated. 
    func_name_string = func_name
    if multi_objective:   func_name_string = func_name[0]  

    # redundancy check
    if not state['FUNCTIONS'].has_key(func_name_string):

        # An objective function is to be evaluated
        problem.kind = 'FUNCTION'

        # Aerodynamics
        if multi_objective or func_name == 'ALL' or func_name in su2io.optnames_aero + su2io.grad_names_directdiff:
            aerodynamics( problem, state )
            
        # Stability
        elif func_name in su2io.optnames_stab:
            stability( problem, state )
        
        # Geometry
        elif func_name in su2io.optnames_geo:
            geometry( func_name, problem, state )
            
        # Structural OF
        elif func_name in su2io.optnames_fea:
            structural( problem, state )
            
        else:
            raise Exception, 'unknown function name, %s' % func_name
        
    #: if not redundant

    # prepare output
    if func_name == 'ALL':
        func_out = state['FUNCTIONS']
    elif (multi_objective):
        # If combine_objective is true, use the 'combo' output.
        objectives=problem.OBJECTIVE_FUNCTION
        func_out = 0.0
        for func in func_name:
            sign = su2io.get_objectiveSign(func)
            func_out+=state['FUNCTIONS'][func]*objectives[func]['SCALE']*sign
        state['FUNCTIONS']['COMBO'] = func_out
    else:
        func_out = state['FUNCTIONS'][func_name]
        
    
    return copy.deepcopy(func_out)

#: def function()


# ----------------------------------------------------------------------
#  Aerodynamic Functions
# ----------------------------------------------------------------------

def aerodynamics( problem, state=None ):
    """ vals = SU2.eval.aerodynamics(config,state=None)
    
        Evaluates aerodynamics with the following:
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

    # Temporary
    config = problem.config

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
    direct_done = all( [ state.FUNCTIONS.has_key(key) for key in su2io.optnames_aero[:9] ] )
    if direct_done:
        # return aerodynamic function values
        aero = su2util.ordered_bunch()
        for key in su2io.optnames_aero:
            if state.FUNCTIONS.has_key(key):
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
        name = files['DIRECT']
        name = su2io.expand_time(name, problem.physics)
        link.extend( name )
        ##config['RESTART_SOL'] = 'YES' # don't override config file
    else:
        config['RESTART_SOL'] = 'NO'
        
    # files: target equivarea distribution
    if ( 'EQUIV_AREA' in special_cases and 
         'TARGET_EA' in files ) : 
        pull.append( files['TARGET_EA'] )

    # files: target pressure distribution
    if ( 'INV_DESIGN_CP' in special_cases and
         'TARGET_CP' in files ) :
        pull.append( files['TARGET_CP'] )

    # files: target heat flux distribution
    if ( 'INV_DESIGN_HEATFLUX' in special_cases and
         'TARGET_HEATFLUX' in files ) :
        pull.append( files['TARGET_HEATFLUX'] )

    # output redirection
    with redirect_folder( 'DIRECT', pull, link ) as push:
        with redirect_output(log_direct):     
            
            # # RUN DIRECT SOLUTION # #
            info = su2run.direct(config)
            su2io.restart2solution(problem,info)
            state.update(info)
            
            # direct files to push
            name = info.FILES['DIRECT']
            name = su2io.expand_time(name, opt.problem)
            push.extend(name)
            
            # equivarea files to push
            if 'WEIGHT_NF' in info.FILES:
                push.append(info.FILES['WEIGHT_NF'])

            # pressure files to push
            if 'TARGET_CP' in info.FILES:
                push.append(info.FILES['TARGET_CP'])

            # heat flux files to push
            if 'TARGET_HEATFLUX' in info.FILES:
                push.append(info.FILES['TARGET_HEATFLUX'])
                
    #: with output redirection
    # return output 
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_aero + su2io.grad_names_directdiff:
        if state['FUNCTIONS'].has_key(key):
            funcs[key] = state['FUNCTIONS'][key]
            
    if 'OUTFLOW_GENERALIZED' in config.OBJECTIVE_FUNCTION:    
        import downstream_function
        state['FUNCTIONS']['OUTFLOW_GENERALIZED']=downstream_function.downstream_function(config,state)

    return funcs

#: def aerodynamics()


# ----------------------------------------------------------------------
#  Stability Functions
# ----------------------------------------------------------------------

def stability( problem, state=None, step=1e-2 ):

    # Temporary
    config = problem.config
    
    folder = 'STABILITY' # os.path.join('STABILITY',func_name) #STABILITY/D_MOMENT_Y_D_ALPHA/
    
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
    #  CENTRAL POINT
    # ----------------------------------------------------    
    
    # will run in DIRECT/
    func_0 = aerodynamics(config,state)      
    
    
    # ----------------------------------------------------    
    #  Run Forward Point
    # ----------------------------------------------------   
    
    # files to pull
    files = state.FILES
    pull = []; link = []
    
    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,config)
    link.extend(name)
    
    # files: direct solution
    if files.has_key('DIRECT'):
        name = files['DIRECT']
        name = su2io.expand_time(name, opt.problem)
        link.extend( name )
        ##config['RESTART_SOL'] = 'YES' # don't override config file
    else:
        config['RESTART_SOL'] = 'NO'
        
    # files: target equivarea distribution
    if ( 'EQUIV_AREA' in special_cases and 
         'TARGET_EA' in files ) : 
        pull.append( files['TARGET_EA'] )

    # files: target pressure distribution
    if ( 'INV_DESIGN_CP' in special_cases and
         'TARGET_CP' in files ) :
        pull.append( files['TARGET_CP'] )

    # files: target heat flux distribution
    if ( 'INV_DESIGN_HEATFLUX' in special_cases and
         'TARGET_HEATFLUX' in files ) :
        pull.append( files['TARGET_HEATFLUX'] )

    # pull needed files, start folder
    with redirect_folder( folder, pull, link ) as push:
        with redirect_output(log_direct):     
            
            konfig = copy.deepcopy(config)
            ztate  = copy.deepcopy(state)
            
            # TODO: GENERALIZE
            konfig.AoA = konfig.AoA + step
            ztate.FUNCTIONS.clear()
            
            func_1 = aerodynamics(konfig,ztate)
                        
            ## direct files to store
            #name = ztate.FILES['DIRECT']
            #if not state.FILES.has_key('STABILITY'):
                #state.FILES.STABILITY = su2io.ordered_bunch()
            #state.FILES.STABILITY['DIRECT'] = name
            
            ## equivarea files to store
            #if 'WEIGHT_NF' in ztate.FILES:
                #state.FILES.STABILITY['WEIGHT_NF'] = ztate.FILES['WEIGHT_NF']
    
    # ----------------------------------------------------    
    #  DIFFERENCING
    # ----------------------------------------------------
        
    for derv_name in su2io.optnames_stab:

        matches = [ k for k in su2io.optnames_aero if k in derv_name ]
        if not len(matches) == 1: continue
        func_name = matches[0]

        obj_func = ( func_1[func_name] - func_0[func_name] ) / step
        
        state.FUNCTIONS[derv_name] = obj_func
    

    # return output 
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_stab:
        if state['FUNCTIONS'].has_key(key):
            funcs[key] = state['FUNCTIONS'][key]    
    
    return funcs
    
    
    
    
# ----------------------------------------------------------------------
#  Geometric Functions
# ----------------------------------------------------------------------

def geometry( func_name, problem, state=None ):
    """ val = SU2.eval.geometry(config,state=None)
    
        Evaluates geometry with the following:
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

    # Temporary
    config = problem.config
    
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
    #  Update Mesh (check with Trent)
    # ----------------------------------------------------
    
    # does decomposition and deformation
    #info = update_mesh(config,state)


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
                config.GEO_MODE  = 'FUNCTION'
                
                # # RUN GEOMETRY SOLUTION # #
                info = su2run.geometry(config)
                state.update(info)
                
                # no files to push
                
        #: with output redirection
        
    #: if not redundant 
    
    # return output 
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_geo:
        if state['FUNCTIONS'].has_key(key):
            funcs[key] = state['FUNCTIONS'][key]
    return funcs
    

#: def geometry()

# ----------------------------------------------------------------------
#  Structural Functions
# ----------------------------------------------------------------------

def structural( problem, state=None ):
    """ vals = SU2.eval.aerodynamics(config,state=None)
    
        Evaluates structural mechanics with the following:
            SU2.run.direct()
        
        Assumptions:
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

    config = problem.config

    # initialize
    state = su2io.State(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = problem.physics.files['MESH']

    # console output
    if problem.config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_direct = 'log_Direct.out'
    else:
        log_direct = None

   
    # ----------------------------------------------------    
    #  Direct Solution
    # ----------------------------------------------------    
    
    # redundancy check
    direct_done = all( [ state.FUNCTIONS.has_key(key) for key in su2io.optnames_fea[:9] ] )
    if direct_done:
        # return structural function values
        struc = su2util.ordered_bunch()
        for key in su2io.optnames_fea:
            if state.FUNCTIONS.has_key(key):
                struc[key] = state.FUNCTIONS[key]
        return copy.deepcopy(struc)
    #: if redundant
    
    # files to pull
    files = state.FILES
    pull = []; link = []

    for key in files:
        if key == 'DIRECT':
            name = files[key]
            name = su2io.expand_time(name, problem.physics)
            link.extend(name)
        else:
            if 'ADJOINT' not in key:
                name = files[key]
                name = su2io.expand_part(name, problem)
                link.extend(name)

    # If direct is not in files means that we are not restarting the solution
    # (This may be removed using CONFIG_DIRECT and CONFIG_ADJOINT)
    if 'DIRECT' not in files:
        config['RESTART_SOL'] = 'NO'

    # output redirection
    with redirect_folder( 'DIRECT', pull, link ) as push:
        with redirect_output(log_direct):     
            
            # # RUN DIRECT SOLUTION # #
            info = su2run.direct(problem)
            su2io.restart2solution(problem, info)
            state.update(info)
            
            # direct files to push
            name = info.FILES['DIRECT']
            name = su2io.expand_time(name, problem.physics)
            push.extend(name)

            # If multizone
            nZone = problem.physics.nZone
            for i in range(1, nZone):
                DIR_LABEL = 'DIRECT_' + str(i)
                name = info.FILES[DIR_LABEL]
                name = su2io.expand_time(name, problem.physics)
                push.extend(name)
                
    #: with output redirection
    # return output 
    funcs = su2util.ordered_bunch()
    for key in su2io.optnames_aero + su2io.grad_names_directdiff:
        if state['FUNCTIONS'].has_key(key):
            funcs[key] = state['FUNCTIONS'][key]
            
    # Update state to have variables
    state.VARIABLES.DV_VALUE_NEW = config.DV_VALUE_NEW

    return funcs

#: def structural()

def update_mesh(opt,state=None):
    """ SU2.eval.update_mesh(config,state=None)
    
        updates mesh with the following:
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

    # Temporary
    config = opt.CONFIG_DIRECT
    
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

