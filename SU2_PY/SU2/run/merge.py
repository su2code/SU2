## \file merge.py
#  \brief python package for merging meshes
#  \author Tom Economon, Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
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
from interface import SOL as SU2_SOL

# ----------------------------------------------------------------------
#  Merge Mesh
# ----------------------------------------------------------------------

def merge( config ):
    """ info = SU2.run.merge(config)
        
        Merges mesh with:
            SU2.run.SOL()    (volume merging)
            internal scripts (surface merging)
            
        Assumptions:
            config.NUMBER_PART is set 
            Skip if config.NUMBER_PART > 1
            
        Inputs:
            config - an SU2 config
                
        Ouputs:
            info - an empty SU2 State
            
        Executes in:
            ./
    """
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # check if needed
    partitions = konfig['NUMBER_PART']
    if partitions <= 1:
        return su2io.State()
    
    # special cases
    special_cases = su2io.get_specialCases(konfig)
    
    # # MERGING # #
    if 'WRT_UNSTEADY' in special_cases:
        merge_unsteady(konfig)
    else:
        merge_solution(konfig)
        
    # info out (empty)
    info = su2io.State()
    
    return info

#: merge

def merge_unsteady( config, begintime=0, endtime=None ):
    
    if not endtime:
        endtime = config.EXT_ITER
    
    # SU2_SOL handles unsteady volume merge
    merge_solution( config )

    return

#: def merge_unsteady()

def merge_solution( config ):
    """ SU2.io.merge.merge_solution(config)
        general volume surface merging with SU2_SOL
    """
    
    SU2_SOL( config )
    
    return

#: merge_solution( config )
