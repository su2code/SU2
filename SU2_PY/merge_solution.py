#!/usr/bin/env python 

## \file merge_solution.py
#  \brief Python script for merging of the solution files.
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
    parser.add_option("-p", "--partitions", dest="partitions", default=-1, 
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
