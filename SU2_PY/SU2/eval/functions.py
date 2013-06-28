
import os, sys, shutil, copy
from .. import run as su2run
from .. import io  as su2io
from ..io import redirect_folder, redirect_output

def function( func_name, config, state=None ):
    
    # initialize
    state = su2io.state.default_state(state)
    
    # redundancy check
    if not state['FUNCTIONS'].has_key(func_name):
        
        # Aerodynamics
        aero = aerodynamics( config, state )
        
        # Geometry (todo...)
        # geom = geometry( config, state )
        
    #: if not redundant
    
    # prepare output
    if func_name == 'ALL':
        func_out = state['FUNCTIONS']
    else:
        func_out = state['FUNCTIONS'][func_name]
    
    return copy.deepcopy(func_out)


def aerodynamics( config, state=None ):
    ''' runs with redundancy protection (using state):
            decomposition
	    deformation
            direct solution    
        updates config and state by reference
        deformation must be preset in config
        evaluates each step in its own folder, returning state.FILES to the super folder
    '''
    
    # ----------------------------------------------------
    #  Initialize    
    # ----------------------------------------------------
    
    # initialize
    state = su2io.state.default_state(state)
    if not state.FILES.has_key('MESH'):
        state.FILES.MESH = config['MESH_FILENAME']
    special_cases = su2io.get_specialCases(config)
    
    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_decomp = 'log_Decomp.out'
        log_deform = 'log_Deform.out'
        log_direct = 'log_Direct.out'
    else:
        log_decomp = None
        log_deform = None
        log_direct = None

    # ----------------------------------------------------
    #  Redundancy Check
    # ----------------------------------------------------    
    
    # master redundancy check
    if len( state.FUNCTIONS.keys() ):  # naive check - geoemtry functions?
        aero = state.FUNCTIONS         # naive update
        return aero        
    
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
            
                # run 
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
    if deform_set and deform_todo :
    
        # files to pull
        pull  = []
        point = config['MESH_FILENAME']
        point = su2io.expand_part(point,config)
        
        # output redirection
        with redirect_folder('DEFORM',pull,point) as push:
            with redirect_output(log_deform):
                
                # run
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
                
    # ----------------------------------------------------    
    #  Adaptation (not implemented)
    # ----------------------------------------------------
    
    #if not state.['ADAPTED_FUNC']:
    #    config = su2run.adaptation(config)
    #    state['ADAPTED_FUNC'] = True
    
    # ----------------------------------------------------    
    #  Direct Solution
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
        link.append( files['DIRECT'] )
        config['RESTART_SOL'] = 'YES'
    else:
        config['RESTART_SOL'] = 'NO'
        
    # files: target equivarea distribution
    if ( 'EQUIV_AREA' in special_cases and 
         'TARGET_EA' in files ) : 
        pull.append( files['TARGET_EA'] )
    
    # output redirection      
    with redirect_folder( 'DIRECT', pull, link ) as push:
        with redirect_output(log_direct):     
            
            # run
            info = su2run.direct(config)
            su2io.restart2solution(config,info)
            state.update(info)
            
            # files to push
            for this_type in ['DIRECT','WEIGHT_NF']:
                if info['FILES'].has_key(this_type):
                    push.append( info['FILES'][this_type] )
            
    #: with output redirection
    
    # return aerodynamic dictionary
    aero = state['FUNCTIONS'] # naive update
    return copy.deepcopy(aero)

#: def aerodynamics()

def geometry( func_name, config, state=None ):
    raise NotImplementedError

