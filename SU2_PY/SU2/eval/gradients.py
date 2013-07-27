## \file gradients.py
#  \brief python package for gradients
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
from .functions import function
from ..io import redirect_folder, redirect_output
import functions

# ----------------------------------------------------------------------
#  Main Gradient Interface
# ----------------------------------------------------------------------

def gradient( func_name, method, config, state=None ):    
    """ val = SU2.eval.grad(func_name,method,config,state=None)
    
        Evaluates the aerodynamic gradients.
        
        Wraps:
            SU2.eval.adjoint()
            SU2.eval.findiff()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh need not be deformed.
            Updates config and state by reference.
            Redundancy if state.GRADIENTS has the key func_name.
            
        Executes in:
            ./ADJOINT_* or ./FINDIFF
        
        Inputs:
            func_name - SU2 objective function name
            method    - 'ADJOINT' or 'FINDIFF'
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            A list of floats of gradient values
    """
    
    # Initialize
    grads = {}
    state = su2io.State(state)
    if func_name == 'ALL':
        raise Exception , "func_name = 'ALL' not yet supported"
    
    # redundancy check
    if not state['GRADIENTS'].has_key(func_name):
        
        # Adjoint Gradients
        if method == 'ADJOINT':
            
            # Aerodynamics
            if func_name in su2io.optnames_aero:
                grads = adjoint( func_name, config, state )
            
            # Geometry (actually a finite difference)
            elif func_name in su2io.optnames_geo:
                grads = geometry( func_name, config, state )
                
            else:
                raise Exception, 'unknown function name'            
            
        # Finite Difference Gradients
        elif method == 'FINDIFF':
            grads = findiff( config, state )
            
        else:
            raise Exception , 'unrecognized gradient method'
        
        # store
        state['GRADIENTS'].update(grads)
    
    # if not redundant
    
    # prepare output
    grads_out = state['GRADIENTS'][func_name]
    
    return copy.deepcopy(grads_out)

#: def gradient()


# ----------------------------------------------------------------------
#  Adjoint Gradients
# ----------------------------------------------------------------------

def adjoint( func_name, config, state=None ):
    """ vals = SU2.eval.adjoint(func_name,config,state=None)
    
        Evaluates the aerodynamics gradients using the 
        adjoint methodology with:
            SU2.eval.func()
               SU2.run.decompose()
	       SU2.run.deform()
               SU2.run.direct()
            SU2.run.adjoint()
            
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Adjoint Redundancy if state.GRADIENTS has key func_name.
            Direct Redundancy if state.FUNCTIONS has key func_name.
        
        Executes in:
            ./ADJOINT_<func_name>
        
        Inputs:
            func_name - SU2 objective function name
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            A Bunch() with keys of objective function names
            and values of list of floats of gradient values
    """
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    special_cases = su2io.get_specialCases(config)
    ADJ_NAME = 'ADJOINT_'+func_name
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_adjoint = 'log_Adjoint.out'
    else:
        log_adjoint = None   
    
    # ----------------------------------------------------
    #  Redundancy Check
    # ----------------------------------------------------    
    
    # master redundancy check
    if state['GRADIENTS'].has_key(func_name):
        grads = state['GRADIENTS']
        return copy.deepcopy(grads)
        
    # ----------------------------------------------------
    #  Direct Solution    
    # ----------------------------------------------------        
        
    # run (includes redundancy checks)
    function( func_name, config, state )   
    
    # ----------------------------------------------------    
    #  Adaptation (not implemented)
    # ----------------------------------------------------
    
    #if not state.['ADAPTED_ADJOINT']:
    #    config = su2run.adaptation(config)
    #    state['ADAPTED_FUNC'] = True
    
    # ----------------------------------------------------    
    #  Adjoint Solution
    # ----------------------------------------------------        
    
    # files to pull
    files = state['FILES']
    pull = []; link = []    
    
    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,config)
    link.extend(name)
    
    # files: direct solution
    link.append( files['DIRECT'] )
    
    # files: adjoint solution
    if files.has_key( ADJ_NAME ):
        link.append( files[ ADJ_NAME ] )
        ##config['RESTART_SOL'] = 'YES' # don't override config file
    else:
        config['RESTART_SOL'] = 'NO'
    
    # files: target equivarea adjoint weights
    if 'EQUIV_AREA' in special_cases:
        pull.append(files['WEIGHT_NF'])    
    
    # output redirection      
    with redirect_folder( ADJ_NAME, pull, link ) as push:
        with redirect_output(log_adjoint):        
            
            # setup config
            config['ADJ_OBJFUNC'] = func_name
            
            # # RUN ADJOINT SOLUTION # #
            info = su2run.adjoint(config)
            su2io.restart2solution(config,info)
            state.update(info)
            
            # Gradient Projection
            info = su2run.projection(config)
            state.update(info)
            
            # files to push
            push.append( state['FILES'][ADJ_NAME] )
            
    #: with output redirection

    # return output 
    grads = state['GRADIENTS']
    return copy.deepcopy(grads)

#: def adjoint()


# ----------------------------------------------------------------------
#  Finite Difference Gradients
# ----------------------------------------------------------------------

def findiff( config, state=None, step=1e-4 ):
    """ vals = SU2.eval.findiff(config,state=None,step=1e-4)
    
        Evaluates the aerodynamics gradients using 
        finite differencing with:
            SU2.eval.func()
               SU2.run.decompose()
	       SU2.run.deform()
               SU2.run.direct()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Gradient Redundancy if state.GRADIENTS has the key func_name.
            Direct Redundancy if state.FUNCTIONS has key func_name.
    
        Executes in:
            ./FINDIFF
            
        Inputs:
            config - an SU2 config
            state  - optional, an SU2 state
            step   - finite difference step size, as a float or
                     list of floats of length n_DV
        
        Outputs:
            A Bunch() with keys of objective function names
            and values of list of floats of gradient values
    """    
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.State(state)
    special_cases = su2io.get_specialCases(config)
    Definition_DV = config['DEFINITION_DV']
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_findiff = 'log_FinDiff.out'
    else:
        log_findiff = None
    
    # ----------------------------------------------------
    #  Redundancy Check
    # ----------------------------------------------------    
        
    # master redundancy check
    opt_names = su2io.optnames_aero + su2io.optnames_geo 
    findiff_todo = all( [ state.GRADIENTS.has_key(key) for key in opt_names ] )
    if findiff_todo:
        grads = state['GRADIENTS']
        return copy.deepcopy(grads)
    
    # ----------------------------------------------------
    #  Zero Step  
    # ----------------------------------------------------   
    
    # run
    func_base = function( 'ALL', config, state )      
    
    # ----------------------------------------------------
    #  Plot Setup
    # ----------------------------------------------------          
    
    grad_filename  = config['GRAD_OBJFUNC_FILENAME']
    grad_filename  = os.path.splitext( grad_filename )[0]
    output_format  = config['OUTPUT_FORMAT']
    plot_extension = su2io.get_extension(output_format)    
    grad_filename  = grad_filename + '_findiff' + plot_extension
    
    # ----------------------------------------------------
    #  Finite Difference Steps
    # ----------------------------------------------------  
    
    # local config
    konfig = copy.deepcopy(config)
    
    # check deformation setup
    n_dv = len(Definition_DV['KIND'])
    deform_set = konfig['DV_KIND'] == Definition_DV['KIND']
    if not deform_set: 
        dvs_base = [0.0] * n_dv
        konfig.unpack_dvs(dvs_base,dvs_base)    
    else:
        dvs_base = konfig['DV_VALUE_NEW']
    
    # initialize gradients
    func_keys = func_base.keys()
    func_keys = ['VARIABLE'] + func_keys + ['FINDIFF_STEP']
    grads = su2util.ordered_bunch.fromkeys(func_keys)
    for key in grads.keys(): grads[key] = []
    
    # step vector
    if isinstance(step,list):
        assert n_dv == len(step) , 'unexpected step vector length'
    else:
        step = [step] * n_dv
    
    # files to pull
    files = state['FILES']
    pull = []; link = []    
    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,konfig)
    link.extend(name)
    # files: direct solution
    link.append( files['DIRECT'] )
    # files: target equivarea distribution
    if 'EQUIV_AREA' in special_cases and 'TARGET_EA' in files:
        pull.append(files['TARGET_EA'])
    
    # output redirection
    with redirect_folder('FINDIFF',pull,link) as push:
        with redirect_output(log_findiff):
            
            # iterate each dv    
            for i_dv in range(n_dv):
                
                this_step = step[i_dv]
                temp_config_name = 'config_FINDIFF_%i.cfg' % i_dv 
                
                this_dvs    = copy.deepcopy(dvs_base)
                this_konfig = copy.deepcopy(konfig)
                this_dvs[i_dv] = this_dvs[i_dv] + this_step
                
                this_state = su2io.State()
                this_state.FILES = copy.deepcopy( state.FILES )
                this_konfig.unpack_dvs(this_dvs,dvs_base)
                
                this_konfig.dump(temp_config_name)
                
                # Direct Solution, findiff step
                func_step = function( 'ALL', this_konfig, this_state )
                
                # remove deform step files
                meshfiles = this_state.FILES.MESH
                meshfiles = su2io.expand_part(meshfiles,this_konfig)
                for name in meshfiles: os.remove(name)
                
                # calc finite difference and store
                for key in grads.keys():
                    if key == 'VARIABLE': 
                        grads[key].append(i_dv)
                    elif key == 'FINDIFF_STEP': 
                        grads[key].append(this_step)
                    else:
                        this_grad = ( func_step[key] - func_base[key] ) / this_step
                        grads[key].append(this_grad)
                #: for each grad name
                
                su2util.write_plot(grad_filename,output_format,grads)
                os.remove(temp_config_name)
            
            #: for each dv
            
    #: with output redirection
    
    # remove plot items
    del grads['VARIABLE']
    del grads['FINDIFF_STEP']
    
    return copy.deepcopy(grads)

#: def findiff()


# ----------------------------------------------------------------------
#  Geometric Gradients
# ----------------------------------------------------------------------

def geometry( func_name, config, state=None ):
    """ val = SU2.eval.gradients.geometry(config,state=None)
    
        Evaluates geometry with the following:
            SU2.eval.functions.geometry()
                SU2.run.decompose()
	        SU2.run.deform()
                SU2.run.geometry()
        
        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Redundancy if state.GRADIENTS does not have func_name.
            
        Executes in:
            ./GEOMETRY
            
        Inputs:
            config    - an SU2 config
            state     - optional, an SU2 state
        
        Outputs:
            A Bunch() with keys of objective function names
            and values of list of floats of gradient values
    """
    
    # # RUN GEOMETRY SOLUTION # #
    if not state.GRADIENTS.has_key(func_name):
        functions.geometry(func_name,config,state)
        
    # return gradient information
    grads = state.GRADIENTS
    return copy.deepcopy(grads)
    
    