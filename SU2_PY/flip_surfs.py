#!/usr/bin/env python

## \file flip_surfs.py
#  \brief Python script for flipping surface elements in SU2 mesh.
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

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6

from __future__ import print_function

from optparse import OptionParser
from numpy import *

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read mesh from FILE", metavar="FILE")
parser.add_option("-o", "--output", dest="outfilename",
                  help="write new mesh to FILE", metavar="OUTFILE")

(options, args)=parser.parse_args()

# Process options
filename = str(options.filename)
outfilename = str(options.outfilename)

with open(filename, 'r') as input:
    with open(outfilename, 'w') as output:
        # Output mesh contents up to markers
        for line in input:
            if 'NMARK' in line:
                print(line.rstrip('\n'), file=output)
                nmark = int(line.rstrip('\n').split("=")[1])
                break
            print(line.rstrip('\n'), file=output)

        for marker in range(nmark):
            nelem = 0
            for line in input:
                if 'ELEMS' in line:
                    print(line.rstrip('\n'), file=output)
                    nelem = int(line.rstrip('\n').split("=")[1])
                    break
                print(line.rstrip('\n'), file=output)

            ielem = 0
            elem = []
            for line in input:
                if ielem == nelem: 
                    print(line.rstrip('\n'), file=output)
                    break
                elem = line.split()
                if (int(elem[0]) == 3):
                    tmp = elem[1]
                    elem[1] = elem[2]
                    elem[2] = tmp
                elem_str = ""
                for idx in elem:
                    elem_str = elem_str + f" {idx}"
                print(elem_str, file=output)
                ielem == ielem + 1

        # Output remainder of mesh
        for line in input:
            print(line.rstrip('\n'), file=output)

