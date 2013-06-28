#!/usr/bin/python

## \file parallel_computation.py
#  \brief Python script for doing the parallel computation using SU2_CFD.
#  \author Current Development: Stanford University.
#          Original Structure: CADES 1.0 (2009).
#  \version 1.0.
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
parser.add_option("-d", "--divide", dest="divide", default="True",
                  help="DIVIDE the numerical grid using SU2_DDC", metavar="DIVIDE")
parser.add_option("-o", "--output", dest="output", default="True",
                  help="OUTPUT the domain solution", metavar="OUTPUT")

(options, args)=parser.parse_args()

# General and default parameters
Config_DDC_file = "config_DDC_" + options.filename
Config_CFD_file = "config_CFD_" + options.filename
output_format = "PARAVIEW"

# Change the parameters to do the domain descomposition method
output_file = open(Config_DDC_file,"w")
for line in open(options.filename):
    if "NUMBER_PART" in line:
        output_file.write("NUMBER_PART= %s \n" % options.partitions)
    elif "OUTPUT_FORMAT" in line:
        output_format = line.split("=")[1].strip()
        output_file.write(line)
    else: output_file.write(line)
output_file.close()

# Just in case the domain decomposition is needed
if options.divide != "False":
    os.system ("./SU2_DDC " + Config_DDC_file)

# Do the CFD computation
output_file = open(Config_CFD_file,"w")
for line in open(options.filename):
    if "MATH_PROBLEM" in line:
        math_problem = line.split("=")[1].strip()
        output_file.write(line)
    elif "RESTART_SOL" in line:
        restart_solution = line.split("=")[1].strip()
        output_file.write(line)
    elif "SOLUTION_FLOW_FILENAME=" in line:
        solution_flow_filename = line.split("=")[1].strip()
        output_file.write(line)
    elif "SOLUTION_ADJ_FILENAME=" in line:
        solution_adj_filename = line.split("=")[1].strip()
        output_file.write(line)
    elif "CADJ_OBJFUNC" in line:
        cadj_objfunc = line.split("=")[1].strip()
        if cadj_objfunc == "DRAG" : adj_prefix = "cd" 
	elif cadj_objfunc == "LIFT" : adj_prefix = "cl"
	elif cadj_objfunc == "SIDEFORCE" : adj_prefix = "csf"
	elif cadj_objfunc == "PRESSURE" : adj_prefix = "cp" 
	elif cadj_objfunc == "MOMENT_X" : adj_prefix = "cmx"
	elif cadj_objfunc == "MOMENT_Y" : adj_prefix = "cmy"
	elif cadj_objfunc == "MOMENT_Z" : adj_prefix = "cmz"
	elif cadj_objfunc == "EFFICIENCY" : adj_prefix = "eff"
	elif cadj_objfunc == "EQUIVALENT_AREA" : adj_prefix = "ea"
	elif cadj_objfunc == "NEARFIELD_PRESSURE" : adj_prefix = "nfp"
        output_file.write(line)
    else: output_file.write(line)
output_file.close()

# A restart solution requires the grid splitting 
if restart_solution == "YES" or math_problem == "ADJOINT": 
    os.system ("./divide_solution_su2.py -p %s -f %s" % (int(options.partitions), Config_DDC_file))

# os.system ("mpirun -np %s --prefix /opt/openmpi SU2_CFD %s" % (int(options.partitions)+1, Config_CFD_file))
os.system ("mpirun -np %s SU2_CFD %s" % (int(options.partitions)+1, Config_CFD_file))
# os.system ("mpiexec -n %s SU2_CFD %s" % (int(options.partitions)+1, Config_CFD_file))

if output_format == "PARAVIEW":
    os.system ("./merge_solution_paraview.py -p %s -f %s -o %s" % (int(options.partitions), Config_CFD_file, options.output))
if output_format == "TECPLOT":
    os.system ("./merge_solution_tecplot.py -p %s -f %s -o %s" % (int(options.partitions), Config_CFD_file, options.output))

# All the restart files must be merged in one file
os.system ("./merge_restart_su2.py -p %s -f %s" % (int(options.partitions), Config_CFD_file))

# Remove solution files at each domain
if restart_solution == "YES" or math_problem == "ADJOINT":
    for domain in range(int(options.partitions)):
        os.remove("%s_%s.%s" % (solution_flow_filename.split(".")[0], domain+1, solution_flow_filename.split(".")[1]))

if restart_solution == "YES" and math_problem == "ADJOINT":
    for domain in range(int(options.partitions)):
        os.remove("%s_%s_%s.%s" % (solution_adj_filename.split(".")[0], domain+1, adj_prefix, solution_adj_filename.split(".")[1]))
        
# Remove configuration and mesh files
if options.divide != "False":
    os.remove(Config_DDC_file)
os.remove(Config_CFD_file)
