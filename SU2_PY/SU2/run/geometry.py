## \file geometry.py
#  \brief python package for running geometry analyses
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.2
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
from interface import GDC       as SU2_GDC
from ..util import ordered_bunch

# ----------------------------------------------------------------------
#  Direct Simulation
# ----------------------------------------------------------------------

def geometry ( config , step = 1e-3 ): 
    """ info = SU2.run.geometry(config)
        
        Runs an geometry analysis with:
            SU2.run.decomp()
            SU2.run.GDC()
            
        Assumptions:
            Redundant decomposition if config.DECOMPOSED == True
            Performs both function and gradient analysis
                        
        Outputs:
            info - SU2 State with keys:
                FUNCTIONS
                GRADIENTS
                
        Updates:
            config.DECOMPOSED
            
        Executes in:
            ./
    """
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # unpack
    function_name = konfig['GEO_PARAM']
    func_filename = 'of_eval.dat'
    grad_filename = 'of_grad.dat'
    
    # does both direct and gradient, very cheap
    konfig.GEO_MODE  = 'GRADIENT' 
    
    # choose dv values 
    Definition_DV = konfig['DEFINITION_DV']
    n_DV          = len(Definition_DV['KIND'])
    if isinstance(step,list):
        assert len(step) == n_DV , 'unexpected step vector length'
    else:
        step = [step]*n_DV
    dv_old = [0.0]*n_DV # SU2_GPC input requirement, assumes linear superposition of design variables
    dv_new = step
    konfig.unpack_dvs(dv_new,dv_old)    
    
    # decompose
    su2decomp(konfig)
    
    # Run Solution
    SU2_GDC(konfig)
    
    # get function values
    func_file = open(func_filename)
    funcs = float( func_file.readline().strip() )
    func_file.close()
    functions = ordered_bunch({function_name : funcs})
    
    # get gradient_values
    grads = su2io.read_gradients(grad_filename)
    gradients = ordered_bunch({function_name : grads})
    
    
    # update super config
    config.update({ 'DECOMPOSED' : konfig['DECOMPOSED'] })
                    
    # info out
    info = su2io.State()
    info.FUNCTIONS.update( functions )
    info.GRADIENTS.update( gradients )
    
    return info
