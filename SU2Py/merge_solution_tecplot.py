#!/usr/bin/python

## \file merge_solution_tecplot.py
#  \brief Python script for merging of the solution files.
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
parser.add_option("-t", "--timestep", dest="timestep", default=0,
                  help="index of the TIMESTEP", metavar="TIMESTEP")
parser.add_option("-o", "--output", dest="output", default="True",
                  help="OUTPUT the domain solution", metavar="OUTPUT")

(options, args)=parser.parse_args()

# Read the kind of grid adaptation and number of cycles
incompressible_formulation = "NO"

for line in open(options.filename):
    if "MATH_PROBLEM" in line:
        math_problem = line.split("=")[1].strip()
    elif "PHYSICAL_PROBLEM" in line:
        physical_problem = line.split("=")[1].strip()
    elif "INCOMPRESSIBLE_FORMULATION" in line:
        incompressible_formulation = line.split("=")[1].strip()
    elif "KIND_TURB_MODEL" in line:
        turb_model = line.split("=")[1].strip()
    elif "VOLUME_FLOW_FILENAME" in line:
        flow_filename = line.split("=")[1].strip()
    elif "VOLUME_ADJ_FILENAME" in line:
        adjoint_filename = line.split("=")[1].strip()
    elif "SURFACE_FLOW_FILENAME" in line:
        surface_flow_filename = line.split("=")[1].strip()
    elif "SURFACE_ADJ_FILENAME" in line:
        surface_adjoint_filename = line.split("=")[1].strip()

# flow_filename = "flow"
# surface_flow_filename = "surface_flow"
# adjoint_filename = "adjoint"
# surface_adjoint_filename = "surface_adjoint"

# Read the total number of points, and elements
if options.output != "False":
    npoint = 0; nelem = 0; nstore = 0; cont = 0
    addpoint = [ 0 ]
    for domain in range(int(options.partitions)):
        if int(options.timestep) == 0:
            if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (flow_filename, domain+1))
            elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (adjoint_filename, domain+1))
        else:
            if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (flow_filename, domain+1, int(options.timestep)))
            elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (adjoint_filename, domain+1, int(options.timestep)))
        for line in input_file:
            if "VARIABLES" in line:
                variables_list = line
            elif "ZONE" in line:
                local_point = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                npoint = npoint + local_point
                if domain == 0: addpoint.append(local_point)
                else: addpoint.append(local_point+addpoint[domain])
                nelem = nelem + int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
                datapacking = line.replace("ZONE"," ").split(",")[2].split("=")[1]
                kindzone = line.replace("ZONE"," ").split(",")[3].split("=")[1]
        input_file.close()
            
    # Open the output file
    if int(options.timestep) == 0:
        if math_problem == "DIRECT": output_file = open("%s.plt" % (flow_filename) ,"w")
        elif math_problem == "ADJOINT": output_file = open("%s.plt" % (adjoint_filename) ,"w")
    else:
        if math_problem == "DIRECT": output_file = open("%s_%s.plt" % (flow_filename, int(options.timestep)) ,"w")
        elif math_problem == "ADJOINT": output_file = open("%s_%s.plt" % (adjoint_filename, int(options.timestep)) ,"w")

    # Write the header
    output_file.write("TITLE = \"Visualization of the volumetric grid\"\n")
        
        # Write the coordinates and solution
    output_file.write("%s" % variables_list)
    output_file.write("ZONE NODES=%s , ELEMENTS=%s, DATAPACKING=%s, ZONETYPE=%s" % (npoint, nelem, datapacking, kindzone))
    for domain in range(int(options.partitions)):
        if int(options.timestep) == 0:
            if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (flow_filename, domain+1))
            elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (adjoint_filename, domain+1))
        else:
            if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (flow_filename, domain+1, int(options.timestep)))
            elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (adjoint_filename, domain+1, int(options.timestep)))
        while 1:
            line = input_file.readline()
            if "ZONE" in line:
                ipoint = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                for list_points in range(int(ipoint)):
                    new_line = input_file.readline()
                    output_file.write("%s" %  new_line)
            if not line: break
        input_file.close()

    # Write the elements
    for domain in range(int(options.partitions)):
        if int(options.timestep) == 0:
            if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (flow_filename, domain+1))
            elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (adjoint_filename, domain+1))
        else:
            if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (flow_filename, domain+1, int(options.timestep)))
            elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (adjoint_filename, domain+1, int(options.timestep)))
        while 1:
            line = input_file.readline()
            if "ZONE" in line:
                ipoint = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                ielem = int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
                for list_points in range(int(ipoint)):
                    new_line = input_file.readline()
                for list_elems in range(int(ielem)):
                    elem_info = input_file.readline().replace("\t"," ").split(" ")
                    cont = 0
                    for iter in elem_info:
                        elem_info[cont] = ("%s" % (int(iter)+addpoint[domain]))
                        cont = cont+1;
                    output_file.write("%s\n" %  ("\t".join(elem_info)))
            if not line: break
        input_file.close()

    output_file.close()

for domain in range(int(options.partitions)):
    if int(options.timestep) == 0:
        if math_problem == "DIRECT": os.remove("%s_%s.plt" % (flow_filename, domain+1))
        elif math_problem == "ADJOINT": os.remove("%s_%s.plt" % (adjoint_filename, domain+1))
    else:
        if math_problem == "DIRECT": os.remove("%s_%s_%s.plt" % (flow_filename, domain+1, int(options.timestep)))
        elif math_problem == "ADJOINT": os.remove("%s_%s_%s.plt" % (adjoint_filename, domain+1, int(options.timestep)))

# Surface file
npoint = 0; nelem = 0; nstore = 0; cont = 0
addpoint = [ 0 ]
for domain in range(int(options.partitions)):
    if int(options.timestep) == 0:
        if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
        elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
    else:
        if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, int(options.timestep)))
        elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, int(options.timestep)))         
    for line in input_file:
        if "VARIABLES" in line:
            variables_list = line
        elif "ZONE" in line:
            local_point = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
            npoint = npoint + local_point
            if domain == 0: addpoint.append(local_point)
            else: addpoint.append(local_point+addpoint[domain])
            nelem = nelem + int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
            datapacking = line.replace("ZONE"," ").split(",")[2].split("=")[1]
            kindzone = line.replace("ZONE"," ").split(",")[3].split("=")[1]
    input_file.close()

# Open the output file
if int(options.timestep) == 0:
    if math_problem == "DIRECT": output_file = open("%s.plt" % (surface_flow_filename) ,"w")
    elif math_problem == "ADJOINT": output_file = open("%s.plt" % (surface_adjoint_filename) ,"w")
else:
    if math_problem == "DIRECT": output_file = open("%s_%s.plt" % (surface_flow_filename, int(options.timestep)) ,"w")
    elif math_problem == "ADJOINT": output_file = open("%s_%s.plt" % (surface_adjoint_filename, int(options.timestep)) ,"w")

# Write the header
output_file.write("TITLE = \"Visualization of the surface grid\"\n")

# Write the coordinates and solution
output_file.write("%s" % variables_list)
output_file.write("ZONE NODES=%s , ELEMENTS=%s, DATAPACKING=%s, ZONETYPE=%s" % (npoint, nelem, datapacking, kindzone))
for domain in range(int(options.partitions)):
    if int(options.timestep) == 0:
        if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
        elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
    else:
        if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, int(options.timestep)))
        elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, int(options.timestep)))  
    while 1:
        line = input_file.readline()
        if "ZONE" in line:
            ipoint = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
            for list_points in range(int(ipoint)):
                new_line = input_file.readline()
                output_file.write("%s" %  new_line)
        if not line: break
    input_file.close()

# Write the elements
for domain in range(int(options.partitions)):
    if int(options.timestep) == 0:
        if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
        elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
    else:
        if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, int(options.timestep)))
        elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, int(options.timestep)))
    while 1:
        line = input_file.readline()
        if "ZONE" in line:
            ipoint = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
            ielem = int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
            for list_points in range(int(ipoint)):
                new_line = input_file.readline()
            for list_elems in range(int(ielem)):
                elem_info = input_file.readline().replace("\t"," ").split(" ")
                cont = 0
                for iter in elem_info:
                    elem_info[cont] = ("%s" % (int(iter)+addpoint[domain]))
                    cont = cont+1;
                output_file.write("%s\n" %  ("\t".join(elem_info)))
        if not line: break
    input_file.close()

for domain in range(int(options.partitions)):
    if int(options.timestep) == 0:
        if math_problem == "DIRECT": os.remove("%s_%s.plt" % (surface_flow_filename, domain+1))
        elif math_problem == "ADJOINT": os.remove("%s_%s.plt" % (surface_adjoint_filename, domain+1))
    else:
        if math_problem == "DIRECT":  os.remove("%s_%s_%s.plt" % (surface_flow_filename, domain+1, int(options.timestep)))
        elif math_problem == "ADJOINT":  os.remove("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, int(options.timestep)))
