#!/usr/bin/python

## \file merge_solution_tecplot.py
#  \brief Python script for merging of the solution files.
#  \author Current Development: Stanford University.
#          Original Structure: CADES 1.0 (2009).
#  \version 1.1.
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

for time in range(int(options.endtime)-int(options.begintime)):
    iteration = int(options.begintime) + time
    if iteration < 10:
        os.system ("merge_solution.py -p %s -f %s -t 0000%s" % (int(options.partitions), options.filename, iteration))
    if iteration >= 10 and iteration < 100:
        os.system ("merge_solution.py -p %s -f %s -t 000%s" % (int(options.partitions), options.filename, iteration))
    if iteration >= 100 and iteration < 1000:
        os.system ("merge_solution.py -p %s -f %s -t 00%s" % (int(options.partitions), options.filename, iteration))
    if iteration >= 1000 and iteration < 10000:
        os.system ("merge_solution.py -p %s -f %s -t 0%s" % (int(options.partitions), options.filename, iteration))
    if iteration >= 10000:
        os.system ("merge_solution.py -p %s -f %s -t %s" % (int(options.partitions), options.filename, iteration))


