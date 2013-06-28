#!/usr/bin/python

## \file merge_grid_su2.py
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

import os, time, numpy
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=2,
                  help="number of PARTITIONS", metavar="PARTITIONS")

(options, args)=parser.parse_args()

# Read the kind of grid adaptation and number of cycles
incompressible_formulation = "NO"

for line in open(options.filename):
    if "MESH_OUT_FILENAME=" in line:
        mesh_out_filename = line.split("=")[1].strip()
    elif "MESH_FILENAME=" in line:
        mesh_filename = line.split("=")[1].strip()
    
        
# Read the total number of points from the different grid partitions, and save the global index.
npoint = 0
nelem = 0
for domain in range(int(options.partitions)):
    input_file = open("%s_%s.%s" % (mesh_out_filename.split(".")[0], domain+1, mesh_out_filename.split(".")[1]))
    while 1:
        line = input_file.readline()
        if "NPOIN" in line:
            local_point = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[1])
            npoint = npoint + local_point
        if not line: break
    input_file.close()

# Create array with the solution
Coordinates = []
for x in range(npoint):
    Coordinates.append("EMPTY")

for domain in range(int(options.partitions)):
    # Read the global coordinates from the input grid file, and create the Local2Global vector.
    input_file = open("%s_%s.%s" % (mesh_out_filename.split(".")[0], domain+1, mesh_out_filename.split(".")[1]))
    while 1:
        line = input_file.readline()
        if "NDIME" in line:
            nDim = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
        if "NPOIN" in line:
            local_point = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
            Local2Global = numpy.zeros(local_point, numpy.int)
            for list_points in range(local_point):
                line_point = input_file.readline()
                LocalIndex = int(line_point.strip().replace("\t"," ").split(" ")[nDim])
                GlobalIndex = int(line_point.strip().replace("\t"," ").split(" ")[nDim+1])
                Local2Global[LocalIndex] = GlobalIndex
        if not line: break
    input_file.close()

    # Read the solution file and save the coordinates solution with the global numbering
    input_file = open("%s_%s.%s" % (mesh_out_filename.split(".")[0], domain+1, mesh_out_filename.split(".")[1]))

    while 1:
        line = input_file.readline()
        if "NPOIN" in line:
            for list_points in range(local_point):
                line = input_file.readline()
                coordinates = line.replace("\t"," ")
                LocalIndex = int(line.strip().replace("\t"," ").split(" ")[nDim])
                GlobalIndex = int(line.strip().replace("\t"," ").split(" ")[nDim+1])
                Coordinates[GlobalIndex] = coordinates
        if not line: break
    input_file.close()

# Write the output file with the coordinates in the right order
output_file =  open("%s.%s" % (mesh_out_filename.split(".")[0], mesh_out_filename.split(".")[1]),"w")

input_file = open("%s.%s" % (mesh_filename.split(".")[0], mesh_filename.split(".")[1]))
while 1:
    line = input_file.readline()
    if "NDIME" in line:
        nDim = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
        output_file.write("NDIM= %s\n" % nDim)
    if "NELEM" in line:
        nElem = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
        output_file.write("NELEM= %s\n" % nElem)
        for list_elem in range(nElem):
            line = input_file.readline()
            output_file.write("%s" % line)
    if "NPOIN" in line:
        nPoint = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
        output_file.write("NPOIN= %s\n" % nPoint)
        for iPoint in range(nPoint):
            line = input_file.readline()
            output_file.write("%s \t %s" %  (iPoint, Coordinates[iPoint]))
    if "NMARK" in line:
        nMark = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
        output_file.write("NMARK= %s\n" % nMark)
        for list_mark in range(nMark):
            line = input_file.readline(); output_file.write("%s" % line)
            line = input_file.readline();
            nElem = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
            output_file.write("MARKER_ELEMS= %s\n" % nElem)
            for list_elem in range(nElem):
                line = input_file.readline()
                output_file.write("%s" % line)  
    if not line: break
    
input_file.close()
output_file.close()

for domain in range(int(options.partitions)):
    os.remove("%s_%s.%s" % (mesh_out_filename.split(".")[0], domain+1, mesh_out_filename.split(".")[1]))
