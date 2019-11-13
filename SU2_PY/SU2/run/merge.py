## \file merge.py
#  \brief python package for merging meshes
#  \author T. Economon, T. Lukaczyk, F. Palacios
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
from .. import io  as su2io
from .interface import SOL as SU2_SOL
from .interface import SOL_FSI as SU2_SOL_FSI

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

    # special cases
    multizone_cases = su2io.get_multizone(konfig)
    
    # # MERGING # #
    if 'FLUID_STRUCTURE_INTERACTION' in multizone_cases:
        merge_multizone(konfig)
    else:
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

def merge_multizone( config, begintime=0, endtime=None ):

    if not endtime:
        endtime = config.TIME_ITER
    
    SU2_SOL_FSI( config )
    
    return

#: merge_solution( config )
