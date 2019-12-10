#!/usr/bin/env python 

## \file remove_periodic_halos.py
#  \brief Python script to remove halo layers from SU2 native ascii
#         mesh file preprocessed by SU2_MSH for periodic calcs prior to v7.
#  \author T. Economon
#  \version 7.0.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read mesh file to remove halos", metavar="FILE")
(options, args)=parser.parse_args()

# Store the filename for the input mesh

filename  = options.filename

# The first reading of the mesh is to isolate the number of periodic
# points along with their global index values.

periodic_points = []

fid = open(filename,"r"); lines = fid.readlines(); fid.close()
for ii in range(len(lines[:])):
  
  # The periodic nodes are in the receive BC and will have an even index
  if "SEND_RECEIVE" in lines[ii] and int(lines[ii+2].split("=")[1]) == -1:
    nElems = int(lines[ii+1].split("=")[1])
    kk = ii+3
    for jj in range(nElems):
      tokens = lines[kk].replace("\n","").split()
      periodic = int(tokens[2])
      if periodic % 2 == 0:
        periodic_points.append(tokens[1])
      kk += 1

# Store the number of periodic points in an easier way

nPeriodic = len(periodic_points)

# Get the new counts for points, elems, and markers / boundary elems

nDim    = 0
nPoint  = 0
nElems  = 0
nMarker = 0

# Get nDim

for ii in range(len(lines[:])):
  if "NDIME=" in lines[ii]:
    tokens = lines[ii].replace("\n","").split("=")
    nDim   = int(tokens[1])
    break

# To get point count just subtract the number of periodic points from total

for ii in range(len(lines[:])):
  if "NPOIN=" in lines[ii]:
    tokens  = lines[ii].replace("\n","").split()
    nPoints = int(tokens[1]) - nPeriodic
    break

# Count up the interior elements, excluding the added halos

for ii in range(len(lines[:])):
  if "NELEM=" in lines[ii]:
    tokens   = lines[ii].replace("\n","").split("=")
    nElemOld = int(tokens[1])
    nElems   = nElemOld
    kk = ii+1
    for jj in range(nElemOld):
      element = lines[kk].replace("\n","").split()
      conn    = element[1:-1]
      isHalo  = False
      for val in periodic_points:
        if val in conn:
          isHalo = True
      if isHalo:
        nElems -= 1
      kk += 1
    break

# Simply subtract 2 from the marker count to remove SEND_RECEIVE BCs

for ii in range(len(lines[:])):
  if "NMARK=" in lines[ii]:
    tokens  = lines[ii].replace("\n","").split("=")
    nMarker = int(tokens[1]) - 2
    break

# Begin writing the new mesh file

fid = open("mesh_no_halo.su2","w");
fid.write("NDIME= " + str(nDim) + "\n")

# Write the points, excluding the periodic points at the end of the list

for ii in range(len(lines[:])):
  if "NPOIN=" in lines[ii]:
    fid.write("NPOIN= " + str(nPoints) + "\n")
    kk = ii+1
    for jj in range(nPoints):
      fid.write(lines[kk])
      kk += 1
    break

# Write the element connectivity, excluding the periodic halos

for ii in range(len(lines[:])):
  if "NELEM=" in lines[ii]:
    tokens   = lines[ii].replace("\n","").split("=")
    nElemOld = int(tokens[1])
    fid.write("NELEM= " + str(nElems) + "\n")
    kk = ii+1
    for jj in range(nElemOld):
      element = lines[kk].replace("\n","").split()
      conn    = element[1:-1]
      isHalo  = False
      for val in periodic_points:
        if val in conn:
          isHalo = True
      if not isHalo:
        fid.write(lines[kk])
      kk += 1
    break

# Write the markers, excluding an halo elements and the SEND_RECEIVE BCs

for ii in range(len(lines[:])):
  if "NMARK=" in lines[ii]:
    fid.write("NMARK= " + str(nMarker) + "\n")
  if "MARKER_TAG=" in lines[ii] and "SEND_RECEIVE" not in lines[ii]:
    fid.write(lines[ii])
    tokens   = lines[ii+1].replace("\n","").split("=")
    nElemOld = int(tokens[1])
    nElems   = nElemOld
    kk = ii+2

    # First count number of boundary elems to be removed

    for jj in range(nElemOld):
      element = lines[kk].replace("\n","").split()
      conn    = element[1:]
      isHalo  = False
      for val in periodic_points:
        if val in conn:
          isHalo = True
      if isHalo:
        nElems -= 1
      kk += 1

    # Now write the list of boundary elems, excluding halos

    kk = ii+2
    fid.write("MARKER_ELEMS= " + str(nElems) + "\n")
    for jj in range(nElemOld):
      element = lines[kk].replace("\n","").split()
      conn    = element[1:]
      isHalo  = False
      for val in periodic_points:
        if val in conn:
          isHalo = True
      if not isHalo:
        fid.write(lines[kk])
      kk += 1

# Close the file and exit

fid.close()
