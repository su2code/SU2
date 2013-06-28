#!/usr/bin/env python 

## \file merge_solution_tecplot.py
#  \brief Python script for merging of the solution files.
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

import os, time
from optparse import OptionParser
from merge_solution import merge_solution


# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-b", "--begintime", dest="begintime", default=-1,
                      help="index of the BEGINTIME", metavar="BEGINTIME")
    parser.add_option("-e", "--endtime", dest="endtime", default=-1,
                      help="index of the ENDTIME", metavar="ENDTIME")
    
    (options, args)=parser.parse_args()

    options.partitions = int(options.partitions)
    options.begintime  = int(options.begintime)
    options.endtime    = int(options.endtime)
    
    merge_unsteady( options.filename   ,
                    options.partitions ,
                    options.begintime  , 
                    options.endtime     )


# -------------------------------------------------------------------
#  Merge Unsteady 
# -------------------------------------------------------------------

def merge_unsteady( filename        ,
                    partitions =  2 , 
                    begintime  = -1 ,
                    endtime    = -1  ):                    
    
    for time in range( endtime-begintime ):
        iteration = begintime + time
        merge_solution( filename, partitions, iteration, "True" )


# -------------------------------------------------------------------
#  Call Main 
# -------------------------------------------------------------------
if __name__=='__main__':
    main()

