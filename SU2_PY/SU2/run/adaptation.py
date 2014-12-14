## \file adjoint.py
#  \brief python package for running adjoint problems
#  \author T. Lukaczyk, F. Palacios
#  \version 3.2.5 "eagle"
#
# Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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

import os, sys, shutil, copy

from .. import io   as su2io
from .. import mesh as su2mesh
from decompose import decompose as su2decomp

def adaptation ( config , kind='' ):
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # check kind
    if kind: konfig['KIND_ADAPT'] = kind
    kind = konfig.get('KIND_ADAPT','NONE')
    if kind == 'NONE': 
        return {}
    
    # check adapted?
    
    # decompose
    su2decomp(konfig)
    
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


 