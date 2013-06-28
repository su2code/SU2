#!/usr/bin/env python 

## \file finite_differences.py
#  \brief Python script for doing the finite differences computation using the SU2 suite.
#  \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
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
        
    parser = OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-s", "--step",       dest="step",       default=1E-4,
                      help="finite difference STEP", metavar="STEP")
    
    (options, args)=parser.parse_args()
    options.partitions = int( options.partitions )
    options.step       = float( options.step )    
    
    finite_differences( options.filename   ,
                        options.partitions ,
                        options.step        )
#: def main()



# -------------------------------------------------------------------
#  Finite Differences Function 
# -------------------------------------------------------------------

def finite_differences( filename          , 
                        partitions = 0    , 
                        step       = 1e-4  ):
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    
    # State
    state = SU2.io.State()
    state.find_files(config)
    
    SU2.eval.gradients.findiff(config,state,step)
    
    return state

#: finite_differences()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
