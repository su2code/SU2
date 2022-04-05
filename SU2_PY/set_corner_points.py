#!/usr/bin/env python

## \file set_corner_points.py
#  \brief Python script for adding corners to SU2 mesh.
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

(options, args)=parser.parse_args()

# Process options
options.filename = str(options.filename)

# Open mesh for reading
meshfile = open(options.filename, "r")

# Find number of markers
for num, line in enumerate(meshfile, 1):
    if "NMARK" in line:
        num_mark = num+1
        if not "=" in line:
            raise Exception("\n\n##ERROR : NMARK has no equals sign, line %s.\n\n" % str(num_mark-1))
        # split across equals sign
        line    = line.split("=",1)
        nmarker = int(line[1].strip())
        break

# Build list of points in each marker
marker_points = [[] for _ in range(nmarker)]
for imarker in range(nmarker):
    for num, line in enumerate(meshfile, num_mark):
        if "MARKER_ELEMS" in line:
            num_mark = num+1
            if not "=" in line:
                raise Exception("\n\n##ERROR : MARKER_ELEMS has no equals sign, line %s.\n\n" % str(num_mark-1))
            # split across equals sign
            line  = line.split("=",1)
            nelem = int(line[1].strip())
            break

    ielem = 0
    for num, line in enumerate(meshfile, num_mark):
        if "MARKER_TAG" in line:
            break
        else:
            line = line.split()
            # Determine VTK type
            # Line element
            if int(line[0]) == 3:
                marker_points[imarker].extend([int(line[1]),int(line[2])])
                ielem += 1
            # Tri element
            elif int(line[0]) == 5:
                marker_points[imarker].extend([int(line[1]),int(line[2]),int(line[3])])
                ielem += 1
            # Quad element
            elif int(line[0]) == 9:
                marker_points[imarker].extend([int(line[1]),int(line[2]),int(line[3]),int(line[4])])
                ielem += 1
            # Other element
            else:
                raise ValueError("\n\n##ERROR : Unknown marker type in marker %s, line %s.\n\n" % (str(imarker), str(num)))

    if ielem != nelem:
        raise ValueError("\n\n##ERROR : Incorrect MARKER_ELEMS= %s in marker %s, found %s elements.\n\n" % (str(nelem), str(imarker), str(ielem)))

# Close mesh
meshfile.close()

# Remove duplicate points
for imarker in range(nmarker):
    marker_points[imarker] = set(marker_points[imarker])

# Compare lists for each marker
corners = []
for imarker in range(nmarker):
    for jmarker in range(nmarker):
        if imarker < jmarker:
            corners.extend(list(marker_points[imarker] & marker_points[jmarker]))

# Print corners
tmp_str = "NCORNERS= " + str(len(corners))
print(tmp_str)
for icorner in range(len(corners)):
    tmp_str = "1 " + str(corners[icorner])
    print(tmp_str)
