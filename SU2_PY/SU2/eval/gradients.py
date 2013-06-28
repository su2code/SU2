
import os, sys, shutil, copy
from .. import run  as su2run
from .. import io   as su2io
from .. import util as su2util
from .functions import function
from ..io import redirect_folder, redirect_output

def gradient( func_name, method, config, state=None ):
    ''' method = "ADJOINT" or "FINDIFF"
        func_name = 'ALL' not yet supported
    '''
    
    # Initialize
    grads = {}
    state = su2io.state.default_state(state)
    
    if func_name == 'ALL':
        raise Exception , "func_name = 'ALL' not yet supported"
    
    # redundancy check
    if not state['GRADIENTS'].has_key(func_name):
        
        # Adjoint Gradients
        if method == 'ADJOINT':
            grads = adjoint(func_name, config, state)
            
        # Finite Difference Gradients
        elif method == 'FINDIFF':
            grads = findiff(config, state)
            
        else:
            raise Exception , 'unrecognized gradient method'
        
        # store
        state['GRADIENTS'].update(grads)
    
    # if not redundant
    
    # prepare output
    grads_out = state['GRADIENTS'][func_name]
    
    return grads_out

def adjoint( func_name, config, state=None ):
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.state.default_state(state)
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
        return grads
        
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
        config['RESTART_SOL'] == 'YES'
    else:
        config['RESTART_SOL'] = 'NO'
    
    # files: target equivarea adjoint weights
    if 'EQUIV_AREA' in special_cases:
        pull.append(files['WEIGHT_NF'])    
    
    # output redirection      
    with redirect_folder( ADJ_NAME, pull, link ) as push:
        with redirect_output(log_adjoint):        
            
            # run
            config['ADJ_OBJFUNC'] = func_name
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
    return grads


def findiff( config, state=None, step=1e-4 ):
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.state.default_state(state)
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
    if len(state['GRADIENTS'].keys()) > 0: # naive check
        grads = state['GRADIENTS']
        return grads
    
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
                
                this_dvs    = copy.deepcopy(dvs_base)
                this_konfig = copy.deepcopy(konfig)
                this_dvs[i_dv] = this_dvs[i_dv] + this_step
                
                this_state = su2io.State()
                this_state.FILES = copy.deepcopy( state.FILES )
                this_konfig.unpack_dvs(this_dvs,dvs_base)
                
                # Direct Solution, findiff step
                func_step = function( 'ALL', this_konfig, this_state )
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
            
            #: for each dv
            
    #: with output redirection
    
    # remove plot items
    del grads['iVar']
    del grads['FinDiff_Step']
    
    return grads