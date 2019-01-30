#!/usr/bin/env python

## \file adjoint.py
#  \brief python package for running adjoint problems
#  \author T. Lukaczyk, F. Palacios
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

import copy

from .. import io   as su2io
from ..io.data import append_nestdict
from .. import mesh as su2mesh

def adaptation ( config , kind='' ):
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # check kind
    if kind: konfig['KIND_ADAPT'] = kind
    kind = konfig.get('KIND_ADAPT','NONE')
    if kind == 'NONE': 
        return {}
    
    # check adapted?
        
    # get adaptation function
    adapt_function = su2mesh.adapt.name_map[kind]
    
    # setup problem
    suffix = 'adapt'
    meshname_orig = konfig['MESH_FILENAME']
    meshname_new  = su2io.add_suffix( konfig['MESH_FILENAME'], suffix )
    konfig['MESH_OUT_FILENAME'] = meshname_new
    
    # Run Adaptation
    info = adapt_function(konfig)
    
    # update super config
    config['MESH_FILENAME'] = meshname_new
    config['KIND_ADAPT']    = kind
    
    # files out
    files = { 'MESH' : meshname_new }
    
    # info out
    append_nestdict( info, { 'FILES' : files } )
    
    return info


 
