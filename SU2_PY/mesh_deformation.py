#!/usr/bin/env python 

## \file mesh_deformation.py
#  \brief Python script for doing the parallel deformation using SU2_MDC.
#  \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.6
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

import os, sys, shutil, copy
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-d", "--divide_grid", dest="divide_grid", default="True",
                      help="DIVIDE_GRID the numerical grid", metavar="DIVIDE_GRID")
    parser.add_option("-m", "--merge_grid",     dest="merge_grid",     default="True",
                      help="MERGE_GRID the deformed grid", metavar="MERGE_GRID")

    (options, args)=parser.parse_args()
    options.partitions = int( options.partitions )
    options.divide_grid = options.divide_grid.upper() == 'TRUE'
    options.merge_grid  = options.merge_grid.upper()  == 'TRUE'
    
    # Run Parallel Comutation
    mesh_deformation ( options.filename    ,
                       options.partitions  , 
                       options.divide_grid ,
                       options.merge_grid   )
#: def main()

  
# -------------------------------------------------------------------
#  Parallel Computation Function
# -------------------------------------------------------------------

def mesh_deformation( filename           ,
                      partitions  = 2    , 
                      divide_grid = True ,
                      merge_grid  = True  ):
    
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.DECOMPOSED  = not divide_grid
    
    # State
    state = SU2.io.State()
    state.FILES.MESH = config.MESH_FILENAME
    
    # Deformation
    info = SU2.run.deform(config)
    state.update(info)
    
    return state

#: mesh_deformation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
    
    