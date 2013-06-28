#!/usr/bin/env python 

## \file divide_solution_su2.py
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

import os, time, numpy
from optparse import OptionParser
from sys import stdout

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=2,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-t", "--timestep", dest="timestep", default=-1,
                  help="index of the TIMESTEP", metavar="TIMESTEP")

(options, args)=parser.parse_args()

# Read the kind of grid adaptation and number of cycles
incompressible_formulation = "NO"

for line in open(options.filename):
    if "MATH_PROBLEM=" in line:
        math_problem = line.split("=")[1].strip()
    elif "MESH_FILENAME=" in line:
        mesh_filename = line.split("=")[1].strip()
    elif "RESTART_SOL=" in line:
        restart_solution = line.split("=")[1].strip()
    elif "SOLUTION_FLOW_FILENAME=" in line:
        solution_flow_filename = line.split("=")[1].strip()
    elif "ADJ_OBJFUNC=" in line:
        adj_objfunc = line.split("=")[1].strip()
        if adj_objfunc == "DRAG" : adj_prefix = "cd" 
        elif adj_objfunc == "LIFT" : adj_prefix = "cl"
        elif adj_objfunc == "SIDEFORCE" : adj_prefix = "csf"
        elif adj_objfunc == "PRESSURE" : adj_prefix = "cp" 
        elif adj_objfunc == "MOMENT_X" : adj_prefix = "cmx"
        elif adj_objfunc == "MOMENT_Y" : adj_prefix = "cmy"
        elif adj_objfunc == "MOMENT_Z" : adj_prefix = "cmz"
        elif adj_objfunc == "EFFICIENCY" : adj_prefix = "eff"
        elif adj_objfunc == "EQUIVALENT_AREA" : adj_prefix = "ea"
        elif adj_objfunc == "NEARFIELD_PRESSURE" : adj_prefix = "nfp"
        elif adj_objfunc == "FORCE_X" : adj_prefix = "cfx"
        elif adj_objfunc == "FORCE_Y" : adj_prefix = "cfy"
        elif adj_objfunc == "FORCE_Z" : adj_prefix = "cfz"
        elif adj_objfunc == "THRUST" : adj_prefix = "ct"
        elif adj_objfunc == "TORQUE" : adj_prefix = "cq"
        elif adj_objfunc == "FIGURE_OF_MERIT" : adj_prefix = "merit"
        elif adj_objfunc == "FREESURFACE" : adj_prefix = "fs"
    elif "SOLUTION_ADJ_FILENAME=" in line:
        solution_adj_filename = line.split("=")[1].strip()


# Read the total number of points from the different grid partitions.
npoint = 0
for domain in range(int(options.partitions)):

    if restart_solution == "YES" and math_problem == "ADJOINT":
        percentage = int(100*float(domain)/(float(options.partitions)-1.0)*(1.0/3.0))
    else:
        percentage = int(100*float(domain)/(float(options.partitions)-1.0)*0.5)
    stdout.write("\rDividing direct solution files (.dat)... %d%%" % percentage)
    stdout.flush()

    input_file = open("%s_%s.%s" % (mesh_filename.split(".")[0], domain+1, mesh_filename.split(".")[1]))
    while 1:
        line = input_file.readline()
        if "NPOIN" in line:
            local_point = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[1])
            npoint = npoint + local_point
            break
        if not line: break
    input_file.close()

# Read the solution
Solution = []
for x in range(npoint):
    Solution.append("EMPTY")

if int(options.timestep) == -1:
    input_file = open("%s" % (solution_flow_filename))
else:
    input_file = open("%s_%s.%s" % (solution_flow_filename.split(".")[0], options.timestep, solution_flow_filename.split(".")[1]))
for list_points in range(npoint):
    line = input_file.readline()
    solution = line.replace("\t"," ").split(" ",1)[1]
    iPoint = int(line.strip().replace("\t"," ").split(" ")[0])
    Solution[iPoint] = solution
input_file.close()

for domain in range(int(options.partitions)):

    if restart_solution == "YES" and math_problem == "ADJOINT":
        percentage = int(100*(1.0/3.0)+100*float(domain)/(float(options.partitions)-1.0)*(1.0/3.0))
    else:
        percentage = int(50+100*float(domain)/(float(options.partitions)-1.0)*0.5)
    stdout.write("\rDividing direct solution files (.dat)... %d%%" % percentage)
    stdout.flush()

    # Read the global coordinates from the input grid file, and create the Local2Global vector.
    if int(options.timestep) == -1:
        input_file = open("%s_%s.%s" % (mesh_filename.split(".")[0], domain+1, mesh_filename.split(".")[1]))
        output_file = open("%s_%s.%s" % (solution_flow_filename.split(".")[0], domain+1, solution_flow_filename.split(".")[1]),"w")
    else:
        input_file = open("%s_%s.%s" % (mesh_filename.split(".")[0], domain+1, mesh_filename.split(".")[1]))
        output_file = open("%s_%s_%s.%s" % (solution_flow_filename.split(".")[0], domain+1, options.timestep, solution_flow_filename.split(".")[1]),"w")
    while 1:
        line = input_file.readline()
        if "NDIM" in line:
            nDim = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
        if "NPOIN" in line:
            local_point = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[1])
            Local2Global = numpy.zeros(local_point, numpy.int)
            for list_points in range(local_point):
                line_point = input_file.readline()
                LocalIndex = int(line_point.strip().replace("\t"," ").split(" ")[nDim])
                GlobalIndex = int(line_point.strip().replace("\t"," ").split(" ")[nDim+1])
                output_file.write("%s \t %s" %  (LocalIndex, Solution[GlobalIndex]))
        if not line: break
    input_file.close()
    output_file.close()


if restart_solution == "YES" and math_problem == "ADJOINT":
    if int(options.timestep) == -1:
        input_file = open("%s_%s.%s" % (solution_adj_filename.split(".")[0], adj_prefix, solution_adj_filename.split(".")[1]))
    else:
        input_file = open("%s_%s_%s.%s" % (solution_adj_filename.split(".")[0], adj_prefix, options.timestep, solution_adj_filename.split(".")[1]))
    for list_points in range(npoint):
        line = input_file.readline()
        solution = line.replace("\t"," ").partition(" ")[2]
        iPoint = int(line.strip().replace("\t"," ").split(" ")[0])
        Solution[iPoint] = solution
    input_file.close()

    for domain in range(int(options.partitions)):

        percentage = int(100*(2.0/3.0)+100*float(domain)/(float(options.partitions)-1.0)*(1.0/3.0))
        stdout.write("\rDividing adjoint solution files (.dat)... %d%%" % percentage)
        stdout.flush()

        # Read the global coordinates from the input grid file, and create the Local2Global vector.
        if int(options.timestep) == -1:
            input_file = open("%s_%s.%s" % (mesh_filename.split(".")[0], domain+1, mesh_filename.split(".")[1]))
            output_file = open("%s_%s_%s.%s" % (solution_adj_filename.split(".")[0], domain+1, adj_prefix, solution_adj_filename.split(".")[1]),"w")
        else:
            input_file = open("%s_%s.%s" % (mesh_filename.split(".")[0], domain+1, mesh_filename.split(".")[1]))
            output_file = open("%s_%s_%s.%s" % (solution_adj_filename.split(".")[0], domain+1, adj_prefix, options.timestep, solution_adj_filename.split(".")[1]),"w")
        while 1:
            line = input_file.readline()
            if "NDIM" in line:
                nDim = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
            if "NPOIN" in line:
                local_point = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[1])
                Local2Global = numpy.zeros(local_point, numpy.int)
                for list_points in range(local_point):
                    line_point = input_file.readline()
                    LocalIndex = int(line_point.strip().replace("\t"," ").split(" ")[nDim])
                    GlobalIndex = int(line_point.strip().replace("\t"," ").split(" ")[nDim+1])
                    output_file.write("%s \t %s" %  (LocalIndex, Solution[GlobalIndex]))
            if not line: break
        input_file.close()
        output_file.close()

stdout.write("\n")
