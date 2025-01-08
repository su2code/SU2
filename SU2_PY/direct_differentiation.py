#!/usr/bin/env python

## \file direct_differentiation.py
#  \brief Python script for doing the direct differentiation computation using the SU2 suite.
#  \author F. Palacios
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

from __future__ import division, print_function, absolute_import
import os, sys, shutil
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

    parser = OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-q", "--quiet",      dest="quiet",      default='False',
                      help="output QUIET to log files", metavar="QUIET")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")

    (options, args)=parser.parse_args()
    options.partitions = int( options.partitions )
    options.quiet      = options.quiet.upper() == 'TRUE'
    options.nzones     = int( options.nzones )

    direct_differentiation( options.filename   ,
                            options.partitions ,
                            options.quiet      ,
                            options.nzones      )
#: def main()


# -------------------------------------------------------------------
#  Direct Differentation Function
# -------------------------------------------------------------------

def direct_differentiation( filename           ,
                            partitions = 0     ,
                            quiet      = False ,
                            nzones     = 1      ):
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.NZONES      = int(nzones)
    config["DIRECT_DIFF"] = 'DESIGN_VARIABLES'

    if quiet:
        config.CONSOLE = 'CONCISE'

    # State
    state = SU2.io.State()
    state.find_files(config)

    foundDerivativeField = False
    for fields in SU2.io.historyOutFields:
        group = SU2.io.historyOutFields[fields]['GROUP']
        if group in config.HISTORY_OUTPUT:
            if SU2.io.historyOutFields[fields]['TYPE'] == 'D_COEFFICIENT':
                foundDerivativeField = True

    if not foundDerivativeField:
        sys.exit('No derivative field found in HISTORY_OUTPUT')

    # link restart files to subfolder DIRECTDIFF, if restart solution is selected
    if config.get('TIME_DOMAIN', 'NO') == 'YES' and config.get('RESTART_SOL', 'NO') == 'YES':
        # check if directory DIRECTDIFF/DIRECT exists, if not, create
        if not os.path.isdir('DIRECTDIFF/DIRECT'):
            if not os.path.isdir('DIRECTDIFF'):
                os.mkdir('DIRECTDIFF')
            os.mkdir('DIRECTDIFF/DIRECT')

        restart_name = config['RESTART_FILENAME'].split('.')[0]
        restart_filename = restart_name + '_' + str(int(config['RESTART_ITER']) - 1).zfill(5) + '.dat'
        if not os.path.isfile('DIRECTDIFF/DIRECT/' + restart_filename):
            #throw, if restart file does not exist
            if not os.path.isfile(restart_filename):
                sys.exit("Error: Restart file <" + restart_filename + "> not found." )
            shutil.copyfile(restart_filename, 'DIRECTDIFF/DIRECT/' + restart_filename)

        # use only, if time integration is second order
        if config.get('TIME_MARCHING', 'NO') == 'DUAL_TIME_STEPPING-2ND_ORDER':
            restart_filename = restart_name + '_' + str(int(config['RESTART_ITER']) - 2).zfill(5) + '.dat'
            if not os.path.isfile('DIRECTDIFF/DIRECT/' + restart_filename):
                # throw, if restart file does not exist
                if not os.path.isfile(restart_filename):
                    sys.exit("Error: Restart file <" + restart_filename + "> not found.")
                shutil.copyfile(restart_filename, 'DIRECTDIFF/DIRECT/' + restart_filename)

    # Direct Differentiation Gradients
    SU2.eval.gradients.directdiff(config,state)

    return state

#: finite_differences()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
