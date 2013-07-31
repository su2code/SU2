## \file direct.py
#  \brief python package for running direct solutions
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

from .. import io  as su2io
from decompose import decompose as su2decomp
from merge     import merge     as su2merge
from interface import CFD       as SU2_CFD

# ----------------------------------------------------------------------
#  Direct Simulation
# ----------------------------------------------------------------------

def direct ( config ): 
    """ info = SU2.run.direct(config)
        
        Runs an adjoint analysis with:
            SU2.run.decomp()
            SU2.run.CFD()
            SU2.run.merge()
            
        Assumptions:
            Redundant decomposition if config.DECOMPOSED == True
            Does not rename restart filename to solution filename
            Adds 'direct' suffix to convergence filename
                        
        Outputs:
            info - SU2 State with keys:
                FUNCTIONS
                HISTORY.DIRECT
                FILES.DIRECT
                
        Updates:
            config.DECOMPOSED
            config.MATH_PROBLEM
            
        Executes in:
            ./
    """
    
    # local copy
    konfig = copy.deepcopy(config)

    # decompose
    su2decomp(konfig)
    
    # setup direct problem
    konfig['MATH_PROBLEM']  = 'DIRECT'
    konfig['CONV_FILENAME'] = konfig['CONV_FILENAME'] + '_direct'    
    
    # Run Solution
    SU2_CFD(konfig)
    
    # merge
    konfig['SOLUTION_FLOW_FILENAME'] = konfig['RESTART_FLOW_FILENAME'] 
    su2merge(konfig)
    
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
