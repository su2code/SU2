#!/usr/bin/env python

## \file adapt.py
#  \brief mesh functions
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

from .. import io  as su2io
from ..run import CFD as SU2_CFD
from ..run import MSH as SU2_MSH

def full(config):
    config = copy.deepcopy(config)
    
    config.KIND_ADAPT = 'FULL'
    
    raise NotImplementedError

 
def full_flow(config):
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # set config
    konfig.KIND_ADAPT = 'FULL_FLOW'
    
    # run MSH
    SU2_MSH(konfig)
    
    return
    
    
def full_adjoint(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError

def grad_flow(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def grad_adjoint(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def grad_flow_adj(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def robust(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def full_linear(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def computable(config):
    config = copy.deepcopy(config)
    
    pass
 
def computable_robust():
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def remaining(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def wake(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def horizontal_plane(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError

# config name map
name_map = { 'FULL'              : full              ,
             'FULL_FLOW'         : full_adjoint      ,
             'GRAD_FLOW'         : grad_flow         ,
             'FULL_ADJOINT'      : full_adjoint      ,
             'GRAD_ADJOINT'      : grad_adjoint      ,
             'GRAD_FLOW_ADJ'     : grad_flow_adj     ,
             'ROBUST'            : robust            ,
             'FULL_LINEAR'       : full_linear       ,
             'COMPUTABLE'        : computable        ,
             'COMPUTABLE_ROBUST' : computable_robust ,
             'REMAINING'         : remaining         ,
             'WAKE'              : wake              ,
             'HORIZONTAL_PLANE'  : horizontal_plane   }
