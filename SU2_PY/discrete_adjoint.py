#!/usr/bin/env python 

## \file discrete_adjoint.py
#  \brief Python script for doing the discrete adjoint computation using the SU2 suite.
#  \author F. Palacios, T. Economon, T. Lukaczyk
#  \version 7.0.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

import os, sys, copy
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
    parser.add_option("-v", "--validate", dest="validate", default="False",
                      help="Validate the gradient using direct diff. mode", metavar="VALIDATION")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")

    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )
    options.step        = float( options.step )
    options.compute     = options.compute.upper() == 'TRUE'
    options.validate    = options.validate.upper() == 'TRUE'
    options.nzones      = int( options.nzones )

    discrete_adjoint( options.filename    ,
                      options.partitions  ,
                      options.compute     ,
                      options.step        ,
                      options.nzones       )

#: def main()


# -------------------------------------------------------------------
#  Discrete Adjoint
# -------------------------------------------------------------------

def discrete_adjoint( filename           ,
                      partitions  = 0    , 
                      compute     = True ,
                      step        = 1e-4 ,
                      nzones      = 1     ):
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.NZONES      = int( nzones )

    # State
    state = SU2.io.State()

    # Force CSV output in order to compute gradients
    config.WRT_CSV_SOL = 'YES'


    config['GRADIENT_METHOD'] = 'DISCRETE_ADJOINT'

    # check for existing files
    if not compute:
        config.RESTART_SOL = 'YES'
        state.find_files(config)
    else:
        state.FILES.MESH = config.MESH_FILENAME

    # Tranfer Convergence Data, if necessary
    konfig = copy.deepcopy(config)

    # Direct Solution
    if compute:
        info = SU2.run.direct(config)
        state.update(info)
        # Update konfig
        konfig = copy.deepcopy(config)

        if konfig.get('WINDOW_CAUCHY_CRIT', 'NO') == 'YES' and konfig.TIME_MARCHING != 'NO':
            konfig['TIME_ITER'] = info.WND_CAUCHY_DATA['TIME_ITER']
            konfig['ITER_AVERAGE_OBJ'] = info.WND_CAUCHY_DATA['ITER_AVERAGE_OBJ']
            konfig['UNST_ADJOINT_ITER'] = info.WND_CAUCHY_DATA['UNST_ADJOINT_ITER']

        SU2.io.restart2solution(konfig,state)

    # Adjoint Solution

    # Run all-at-once
    if compute:
        info = SU2.run.adjoint(konfig)
        state.update(info)
        SU2.io.restart2solution(konfig,state)

    # Gradient Projection
    info = SU2.run.projection(konfig,step)
    state.update(info)

    return state

#: continuous_adjoint()

# -------------------------------------------------------------------
#  Alternate Forumulation
# -------------------------------------------------------------------

def discrete_design( filename           ,
                       partitions  = 0    , 
                       compute     = True ,
                       step        = 1e-4 ,
                       validation  = False):
    
    # TODO: 
    # step
    
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions

    config['GRADIENT_METHOD'] = 'DISCRETE_ADJOINT'

    ADJ_NAME = config.OBJECTIVE_FUNCTION
    
    # State
    state = SU2.io.State()

    state_directdiff = SU2.io.State()

    grads_directdiff = []

#    if validation:
#        state_directdiff.find_files(config)
#        konfig = copy.deepcopy(config)
#        konfig['DIRECT_DIFF'] = "DESIGN_VARIABLES"
#        grad_directdiff = SU2.eval.gradients.directdiff(konfig,state_directdiff)
#        state['FILES']['DIRECT'] = 'DIRECTDIFF/' + state_directdiff['FILES']['DIRECT']
#        state['FUNCTIONS'] = state_directdiff['FUNCTIONS']
    
    # check for existing files
    if any([not compute, validation]) :
        state.find_files(config)
    else:
        state.FILES.MESH = config.MESH_FILENAME

    # Adjoint Gradient
    grads = SU2.eval.grad( ADJ_NAME, config['GRADIENT_METHOD'], config, state )

#    if validation:
#        Definition_DV = config['DEFINITION_DV']
#        n_dv = len(Definition_DV['KIND'])
#        grads_dd  = grad_directdiff[ADJ_NAME]
#        print("Validation Summary")
#        print("--------------------------")
#        print("VARIABLE   " + "DISCRETE ADJOINT"  + "  DIRECT DIFFERENTIATION" + "  ERROR (%)")
#        for idv in range(n_dv):
#            if abs(grads[idv]) > abs(grads_dd[idv]):
#                this_err = abs(grads[idv]/grads_dd[idv])
#            else:
#                this_err = abs(grads_dd[idv]/grads[idv])

#            print(str(idv) + "         " + str(grads[idv]) + "         " + str(grads_dd[idv]) + "        " + str((this_err-1)*100)  + ' %')

    
    return state




# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

