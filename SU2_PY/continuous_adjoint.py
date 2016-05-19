#!/usr/bin/env python 

## \file continuous_adjoint.py
#  \brief Python script for continuous adjoint computation using the SU2 suite.
#  \author F. Palacios, T. Economon, T. Lukaczyk
#  \version 4.1.2 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

import os, sys, shutil, copy
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    # Command Line Options
    parser=OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-c", "--compute",    dest="compute",    default="True",
                      help="COMPUTE direct and adjoint problem", metavar="COMPUTE")
    parser.add_option("-s", "--step",       dest="step",       default=1E-4,
                      help="DOT finite difference STEP", metavar="STEP")    
    
    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )
    options.step        = float( options.step )
    options.compute     = options.compute.upper() == 'TRUE'
    
    continuous_adjoint( options.filename    ,
                        options.partitions  ,
                        options.compute     ,
                        options.step         )
        
#: def main()


# -------------------------------------------------------------------
#  Continuous Adjoint 
# -------------------------------------------------------------------

def continuous_adjoint( filename           , 
                        partitions  = 0    , 
                        compute     = True ,
                        step        = 1e-4  ):
    
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    
    # State
    state = SU2.io.State()
    
    # Force CSV output in order to compute gradients
    config.WRT_CSV_SOL = 'YES'
    
    # check for existing files
    if not compute:
        config.RESTART_SOL = 'YES'
        state.find_files(config)
    else:
        state.FILES.MESH = config.MESH_FILENAME
    
    # Direct Solution
    if compute:
        info = SU2.run.direct(config) 
        state.update(info)
        SU2.io.restart2solution(config,state)

    # If using chain rule update coefficients using gradients as defined in downstream_function (local file)
    if config.OBJECTIVE_FUNCTION == 'OUTFLOW_GENERALIZED':
        import downstream_function # Must be defined in run folder
        chaingrad = downstream_function.downstream_gradient(config,state,step)
        # Set coefficients for gradients
        config.OBJ_CHAIN_RULE_COEFF = str(chaingrad[0:5])
    
    # Adjoint Solution
    if compute:
        info = SU2.run.adjoint(config)
        state.update(info)
        #SU2.io.restart2solution(config,state)
    
    # Gradient Projection
    info = SU2.run.projection(config,state, step)
    state.update(info)
    
    return state

#: continuous_adjoint()

# -------------------------------------------------------------------
#  Alternate Forumulation
# -------------------------------------------------------------------

def continuous_design( filename           , 
                       partitions  = 0    , 
                       compute     = True ,
                       step        = 1e-4  ):
    
    # TODO: 
    # step
    
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions

    ADJ_NAME = config.OBJECTIVE_FUNCTION
    
    # State
    state = SU2.io.State()
    
    # check for existing files
    if not compute:
        state.find_files(config)
    else:
        state.FILES.MESH = config.MESH_FILENAME    
    
    # Adjoint Gradient
    grads = SU2.eval.grad( ADJ_NAME, 'CONTINUOUS_ADJOINT', config, state )
    
    return state


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

