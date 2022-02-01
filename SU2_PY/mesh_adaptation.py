#!/usr/bin/env python 

## \file mesh_adaptation.py
#  \brief python script for mesh adaptation with SU2-AMG interface
#  \author Brian Mungu\'ia
#  \version 7.3.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
    parser.add_option("-e", "--stderr", dest="stderr", default="False",
                      help="print stderr files", metavar="STDERR")

    (options, args)=parser.parse_args()

    options.partitions = int( options.partitions )
    options.stderr     = options.stderr == "True"
    
    # Run Mesh Adaptation
    mesh_adaptation ( options.filename   ,
                      options.partitions ,
                      options.stderr )

#: def main()


# -------------------------------------------------------------------
#  Mesh Adaptation Function
# -------------------------------------------------------------------

def mesh_adaptation( filename       ,
                     partitions = 0 ,
                     stderr     = False ):
    
    if not filename:
    	sys.stderr.write("  ## ERROR : a .cfg file must be provided.\n");
    	sys.exit(1)
    
    # Set the name of the configuration file
    config_name = filename
    
    # Read the specified configuration file
    config = SU2.io.Config(config_name)
    
    # Set the number of partitions for parallel computations
    config.NUMBER_PART = partitions
    
    # Call CFD to generate a solution
    SU2.run.amg(config, stderr)
    
#: def mesh_adaptation()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
