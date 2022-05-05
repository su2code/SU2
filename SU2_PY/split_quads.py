#!/usr/bin/env python

## \file split_quads.py
#  \brief Python script for splitting quads into tris in SU2 mesh.
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
        # Output mesh contents up to elements
        for line in input:
            if 'NELEM' in line:
                break
            print(line.rstrip('\n'), file=output)

        # Store quads
        quads = []
        for line in input:
            if '%' in line or 'N' in line:
                end_quads = line
                break
            quads.append(line.split()[1:-1])

        # Split and output quads
        print(f"NELEM= {len(quads)*2}", file=output)
        for i, quad in enumerate(quads):
            print(f" 5 {quad[0]} {quad[1]} {quad[2]} {i*2}", file=output)
            print(f" 5 {quad[2]} {quad[3]} {quad[0]} {i*2+1}", file=output)

        # Output remainder of mesh
        print(end_quads.rstrip('\n'), file=output)
        for line in input:
            print(line.rstrip('\n'), file=output)

