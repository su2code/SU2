#!/usr/bin/env python 

## \file mesh_adaptation.py
#  \brief Python script for doing the grid adaptation using the SU2 suite.
#  \author F. Palacios
#  \version 6.0.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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
    mesh_adaptation_amg ( options.filename   ,
                          options.partitions ,
                          options.stderr )

#: def main()


# -------------------------------------------------------------------
#  Mesh Adaptation Function
# -------------------------------------------------------------------

def mesh_adaptation_amg( filename       ,
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
