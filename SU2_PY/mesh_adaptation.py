#!/usr/bin/env python 

## \file mesh_adaptation.py
#  \brief Python script for doing the grid adaptation using the SU2 suite.
#  \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 3.2.3 "eagle"
#
# SU2, Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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

import SU2
import os, time, sys, shutil
from optparse import OptionParser

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main(): 

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=0,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-p", "--oldpartitions", dest="oldpartitions", default="oldpartitions",
                      help="old number of PARTITIONS (use -n instead)", metavar="OLDPARTITIONS")
    parser.add_option("-c", "--cycle", dest="cycle", default=1,
                      help="number of CYCLE adaptations", metavar="CYCLE")
    parser.add_option("-o", "--overwrite", dest="overwrite", default="False",
                      help="OVERWRITE_MESH the output mesh with the adapted one", metavar="OVERWRITE_MESH")
    parser.add_option("-s", "--save_all", dest="save_all", default="False",
                      help="SAVE_ALL the flow/adjoint/meshes solutions at each adaptation cycle", metavar="SAVE_ALL")

    (options, args)=parser.parse_args()

    options.partitions = int( options.partitions )
    options.cycle      = int( options.cycle      )
    options.overwrite  = options.overwrite == "True"    
    options.save_all   = options.save_all  == "True"

    if options.oldpartitions != "oldpartitions":
        print ("\n IMPORTANT: -p is no longer available in SU2 v3.2.3, use -n flag instead \n")
        sys.exit()
    
    # Run Mesh Adaptation
    mesh_adaptation ( options.filename   ,
                      options.partitions ,
                      options.cycle      ,
                      options.overwrite  ,
                      options.save_all    )

#: def main()


# -------------------------------------------------------------------
#  Mesh Adaptation Function
# -------------------------------------------------------------------

def mesh_adaptation( filename             ,
                     partitions   = 0     , 
                     cycles       = 1     ,
                     overwrite    = False ,
                     save_all     = False  ):

    # Set the name of the configuration file
    config_name = filename
    
    # Read the specified configuration file
    config = SU2.io.Config(config_name)
    
    # Set the number of partitions for parallel computations
    config.NUMBER_PART = partitions
    
    # Run SU2_PRT for parallel computations
    SU2.run.decompose(config)

    # Call CFD to generate a solution
    SU2.run.CFD(config)

    # Rename the output restart to the input solution file
    SU2.io.restart2solution(config)

    # Call MSH
    SU2.run.MSH(config)



#: def mesh_adaptation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
