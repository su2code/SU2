## \file adjoint.py
#  \brief python package for running adjoint problems 
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
#  Adjoint Simulation
# ----------------------------------------------------------------------

def adjoint( config ): 
    """ info = SU2.run.adjoint(config)
        
        Runs an adjoint analysis with:
            SU2.run.decomp()
            SU2.run.CFD()
            SU2.run.merge()
            
        Assumptions:
            Redundant decomposition if config.DECOMPOSED == True
            Does not run Gradient Projection
            Does not rename restart filename to solution filename
            Adds 'adjoint' suffix to convergence filename
            
        Outputs:
            info - SU2 State with keys:
                HISTORY.ADJOINT_NAME
                FILES.ADJOINT_NAME
                
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
    
    # setup problem    
    konfig['MATH_PROBLEM']  = 'ADJOINT'
    konfig['CONV_FILENAME'] = konfig['CONV_FILENAME'] + '_adjoint'
    
    # Run Solution
    SU2_CFD(konfig)
    
    # merge
    konfig['SOLUTION_ADJ_FILENAME'] = konfig['RESTART_ADJ_FILENAME'] 
    su2merge(konfig)
    
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
