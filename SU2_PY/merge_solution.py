#!/usr/bin/env python 

## \file merge_solution.py
#  \brief Python script for merging of the solution files.
#  \author F. Palacios
#  \version 4.1.2 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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


import os, time, numpy, sys
from optparse import OptionParser
import SU2

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=-1,
                      help="number of PARTITIONS", metavar="PARTITIONS")

    (options, args)=parser.parse_args()
    options.partitions = int(options.partitions)
    
    merge_solution( options.filename   ,
                    options.partitions  )


# -------------------------------------------------------------------
#  MERGE SOLUTION 
# -------------------------------------------------------------------

def merge_solution( filename        ,
                    partitions = -1  ):
    
    config = SU2.io.Config(filename)
    
    if partitions > -1 :
        config.NUMBER_PART = partitions
        
    SU2.run.merge(config)
    
#: def merge_solution()

if __name__ == '__main__':
    main()
