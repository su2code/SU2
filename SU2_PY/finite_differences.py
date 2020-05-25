#!/usr/bin/env python 

## \file finite_differences.py
#  \brief Python script for doing the finite differences computation using the SU2 suite.
#  \author F. Palacios
#  \version 7.0.4 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

import os, sys
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

    finite_differences( options.filename   ,
                        options.partitions ,
                        options.quiet      ,
                        options.nzones      )
#: def main()


# -------------------------------------------------------------------
#  Finite Differences Function 
# -------------------------------------------------------------------

def finite_differences( filename           , 
                        partitions = 0     , 
                        quiet      = False ,
                        nzones     = 1      ):
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.NZONES      = int( nzones )

    if quiet: 
        config.CONSOLE = 'CONCISE'
    
    # State
    state = SU2.io.State()
    state.find_files(config)
    
    # Finite Difference Gradients
    SU2.eval.gradients.findiff(config,state)
    
    return state

#: finite_differences()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
