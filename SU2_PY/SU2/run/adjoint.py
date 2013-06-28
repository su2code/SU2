
import os, sys, shutil, copy

from .. import io  as su2io
from decompose import decompose as su2decomp
from interface import CFD as SU2_CFD

def adjoint ( config ): # cleanup?
    
    # local copy
    konfig = copy.deepcopy(config)

    # decompose
    su2decomp(konfig)
    
    # setup problem    
    konfig['MATH_PROBLEM']  = 'ADJOINT'
    konfig['CONV_FILENAME'] = konfig['CONV_FILENAME'] + '_adjoint'
    
    # Run Solution
    SU2_CFD(konfig)
    
    # filenames
    plot_format      = konfig['OUTPUT_FORMAT']
    plot_extension   = su2io.get_extension(plot_format)
    history_filename = konfig['CONV_FILENAME'] + plot_extension
    special_cases    = su2io.get_specialCases(konfig)
    
    # get history
    history = su2io.read_history( history_filename )
    
    # update super config
    config.update({ 'DECOMPOSED'   : konfig['DECOMPOSED']   ,
                    'MATH_PROBLEM' : konfig['MATH_PROBLEM'] ,
                    'ADJ_OBJFUNC'  : konfig['ADJ_OBJFUNC']   })
    
    # files out
    objective    = konfig['ADJ_OBJFUNC']
    adj_title    = 'ADJOINT_' + objective
    suffix       = su2io.get_adjointSuffix(objective)
    restart_name = konfig['RESTART_FLOW_FILENAME']
    restart_name = su2io.add_suffix(restart_name,suffix)
    
    # info out
    info = su2io.State()
    info.FILES[adj_title] = restart_name
    info.HISTORY[adj_title] = history
    
    return info
