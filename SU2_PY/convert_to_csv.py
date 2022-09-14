#!/usr/bin/env python

## \file convert_to_csv.py
#  \brief This script converts SU2 ASCII restart files generated with a version prior v7 to the CSV format
#  \author T. Albring
#  \version 7.3.1 "Blackbird"
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

from optparse import OptionParser
import os

parser = OptionParser(usage = "%prog -i INPUT_FILE",
        description = 'This script converts SU2 ASCII restart files generated with a version prior v7 to the CSV format')
parser.add_option("-i", "--inputfile", dest="infile",
                  help="ASCII restart file (*.dat)", metavar="INPUT_FILE")
(options, args)=parser.parse_args()


infile = open(options.infile, "r")
out_name = options.infile.split('.')[0] + ".csv"
if (os.path.isfile(out_name)):
    print('File ' + out_name + ' already exists.')
    exit(1)
outfile = open(out_name, "w")

while 1:
    line = infile.readline()
    if not line:
        break

    line = line.split() 
    for i, val in enumerate(line):
        outfile.write(val.strip())
        if i != len(line)-1:
            outfile.write(', ')
    outfile.write('\n')

print('Converted ' + options.infile + ' to ' + out_name)
