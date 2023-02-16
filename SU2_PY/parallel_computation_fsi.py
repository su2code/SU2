#!/usr/bin/env python

## \file parallel_computation_fsi.py
#  \brief Python script for running FSI simulations using the SU2 suite.
#  \author T. Economon, T. Lukaczyk, F. Palacios, H. Kline, R. Sanchez
#  \version 7.5.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-c", "--compute",    dest="compute",    default="True",
                      help="COMPUTE direct and adjoint problem", metavar="COMPUTE")

    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )
    options.compute     = options.compute.upper() == 'TRUE'

    if options.filename == None:
        raise Exception("No config file provided. Use -f flag")

    parallel_computation( options.filename    ,
                          options.partitions  ,
                          options.compute      )

#: def main()


# -------------------------------------------------------------------
#  CFD Solution
# -------------------------------------------------------------------

def parallel_computation( filename           ,
                          partitions  = 0    ,
                          compute     = True  ):

    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions

    # State
    state = SU2.io.State()

    # check for existing files
    if not compute:
        state.find_files(config)
    else:
        state.FILES.MESH = config.MESH_FILENAME

    # CFD Solution (direct or adjoint)
    info = SU2.run.CFD(config)
    state.update(info)

    # Solution merging
    if config.SOLVER == 'FEM_ELASTICITY':
        config.SOLUTION_FILENAME = config.RESTART_FILENAME
    elif config.SOLVER == 'FLUID_STRUCTURE_INTERACTION':
        config.SOLUTION_FILENAME = config.RESTART_FILENAME
        config.SOLUTION_FILENAME = config.RESTART_FILENAME

    info = SU2.run.merge(config)
    state.update(info)

    return state

#: parallel_computation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

