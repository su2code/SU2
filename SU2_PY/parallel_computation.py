#!/usr/bin/env python 

## \file parallel_computation.py
#  \brief Python script for doing the continuous adjoint computation using the SU2 suite.
#  \author Francisco Palacios, Tom Economon, Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.2
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
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

from optparse import OptionParser
import SU2

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    # Command Line Options
    parser=OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-c", "--compute",    dest="compute",    default="True",
                      help="COMPUTE direct and adjoint problem", metavar="COMPUTE")
    parser.add_option("-d", "--divide_grid",dest="divide_grid",default="True",
                      help="DIVIDE_GRID the numerical grid", metavar="DIVIDE_GRID")
    
    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )
    options.compute     = options.compute.upper() == 'TRUE'
    options.divide_grid = options.divide_grid.upper() == 'TRUE'
    
    parallel_computation( options.filename    ,
                          options.partitions  ,
                          options.compute     ,
                          options.divide_grid  )
        
#: def main()


# -------------------------------------------------------------------
#  Direct Solution
# -------------------------------------------------------------------

def parallel_computation( filename           , 
                          partitions  = 0    , 
                          compute     = True ,
                          divide_grid = True  ):
    
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.DECOMPOSED  = divide_grid
    
    # State
    state = SU2.io.State()
    
    # check for existing files
    if not compute:
        state.find_files(config)
    else:
        state.FILES.MESH = config.MESH_FILENAME
    
    # Direct Solution
    info = SU2.run.direct(config) 
    state.update(info)
    SU2.io.restart2solution(config,state)
    
    return state

#: parallel_computation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

