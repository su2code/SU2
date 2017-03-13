#!/usr/bin/env python

## \file direct.py
#  \brief python package for running direct solutions
#  \author T. Lukaczyk, F. Palacios
#  \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
from merge     import merge     as su2merge
from interface import CFD       as SU2_CFD

# ----------------------------------------------------------------------
#  Direct Simulation
# ----------------------------------------------------------------------

def direct ( problem ):
    """ info = SU2.run.direct(config)
        
        Runs an adjoint analysis with:
            SU2.run.decomp()
            SU2.run.CFD()
            SU2.run.merge()
            
        Assumptions:
            Does not rename restart filename to solution filename
            Adds 'direct' suffix to convergence filename
                        
        Outputs:
            info - SU2 State with keys:
                FUNCTIONS
                HISTORY.DIRECT
                FILES.DIRECT
                
        Updates:
            config.MATH_PROBLEM
            
        Executes in:
            ./
    """
    
    # local copy
    provlem = copy.deepcopy(problem)
    konfig = provlem.config

    # setup direct problem
    konfig['MATH_PROBLEM']  = 'DIRECT'
    konfig['CONV_FILENAME'] = konfig['CONV_FILENAME'] + '_direct'    
    
    direct_diff = konfig.get('DIRECT_DIFF','NO') == "YES"

    # Run Solution
    SU2_CFD(konfig)
    
    # multizone cases
    multizone_cases = su2io.get_multizone(konfig)

    # merge
    problem.physics.merge_solution(konfig)
    su2merge(konfig)

    # update super config
    problem.config.update({'MATH_PROBLEM': konfig['MATH_PROBLEM']})

    # create state object to send info out
    info = su2io.State()

    # for the kinds of Objective Functions in the problem, read the output (state.history, state.functions)
    for key in problem.ofunction:
        problem.ofunction[key].read_output(provlem, info)
    
    return info
