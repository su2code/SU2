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

# ----------------------------------------------------------------------
#  Direct Simulation
# ----------------------------------------------------------------------

def geometry ( config ): 
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
    
    # Run Solution
    SU2_GDC(konfig)
    
    # get function values
    functions = {}
    
    # update super config
    config.update({ 'DECOMPOSED' : konfig['DECOMPOSED'] })
                    
    # info out
    info = su2io.State()
    info.FUNCTIONS.update( functions )
    
    return info
