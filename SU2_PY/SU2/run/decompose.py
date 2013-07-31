## \file decompose.py
#  \brief python package for decomposing meshes
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
from interface import DDC as SU2_DDC

# ----------------------------------------------------------------------
#  Decompose Mesh
# ----------------------------------------------------------------------

def decompose( config ):
    """ info = SU2.run.decompose(config)
        
        Decomposes mesh with:
            SU2.run.DDC()
            
        Assumptions:
            config.NUMBER_PART is set
            Redundant decomposition if config.DECOMPOSED == True
            Redundant decomposition if config.NUMBER_PART <= 1
            
        Outputs:
            info - an SU2 State with keys:
                FILES.MESH
                
        Updates:
            config.DECOMPOSED
            
        Executes in:
            ./
    """
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # check if needed
    partitions = konfig['NUMBER_PART']
    decomposed = konfig.get('DECOMPOSED',False)
    if partitions <= 1 or decomposed:
        return su2io.State()
    
    # Run Decomposition
    SU2_DDC(konfig)
    
    # update config super copy
    config['DECOMPOSED'] = True
    
    # info out
    info = su2io.State()
    info.FILES.MESH = config['MESH_FILENAME']
    
    return info
