
import os, sys, shutil, copy

from .. import io  as su2io
from decompose import decompose as su2decomp
from interface import CFD as SU2_CFD

def direct ( config ): # cleanup?
    
    # local copy
    konfig = copy.deepcopy(config)

    # decompose
    su2decomp(konfig)
    
    # setup direct problem
    konfig['MATH_PROBLEM']  = 'DIRECT'
    konfig['CONV_FILENAME'] = konfig['CONV_FILENAME'] + '_direct'    
    
    # Run Solution
    SU2_CFD(konfig)
    
    # filenames
    plot_format      = konfig['OUTPUT_FORMAT']
    plot_extension   = su2io.get_extension(plot_format)
    history_filename = konfig['CONV_FILENAME'] + plot_extension
    special_cases    = su2io.get_specialCases(konfig)

    # get history and objectives
    history      = su2io.read_history( history_filename )
    aerodynamics = su2io.read_aerodynamics( history_filename , special_cases )
    
    # update super config
    config.update({ 'DECOMPOSED'   : konfig['DECOMPOSED']   ,
                    'MATH_PROBLEM' : konfig['MATH_PROBLEM']  })
                    
    # info out
    info = su2io.State()
    info.FUNCTIONS.update( aerodynamics )
    info.FILES.DIRECT = konfig['RESTART_FLOW_FILENAME']
    if 'EQUIV_AREA' in special_cases:
        info.FILES.WEIGHT_NF = 'WeightNF.dat'
    info.HISTORY.DIRECT = history
    
    return info
