#!/usr/bin/env python 

## \file merge_solution.py
#  \brief Python script for merging of the solution files.
#  \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.
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


import os, time, libSU2, numpy
from optparse import OptionParser
from sys import stdout

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=2,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-t", "--timestep", dest="timestep", default=-1,
                  help="index of the TIMESTEP", metavar="TIMESTEP")
parser.add_option("-o", "--output", dest="output", default="True",
                  help="OUTPUT the domain solution", metavar="OUTPUT")

(options, args)=parser.parse_args()

incompressible_formulation = "NO"
output_format              = "TECPLOT"
merge_vol                  = "YES"
merge_srf                  = "YES"

for line in open(options.filename):
    if "MATH_PROBLEM=" in line:
        math_problem = line.split("=")[1].strip()
    elif "PHYSICAL_PROBLEM=" in line:
        physical_problem = line.split("=")[1].strip()
    elif "INCOMPRESSIBLE_FORMULATION=" in line:
        incompressible_formulation = line.split("=")[1].strip()
    elif "KIND_TURB_MODEL=" in line:
        turb_model = line.split("=")[1].strip()
    elif "VOLUME_FLOW_FILENAME=" in line:
        volume_flow_filename = line.split("=")[1].strip()
    elif "VOLUME_ADJ_FILENAME=" in line:
        volume_adjoint_filename = line.split("=")[1].strip()
    elif "SURFACE_FLOW_FILENAME=" in line:
        surface_flow_filename = line.split("=")[1].strip()
    elif "SURFACE_ADJ_FILENAME=" in line:
        surface_adjoint_filename = line.split("=")[1].strip()
    elif "MESH_FILENAME=" in line:
        mesh_filename = line.split("=")[1].strip()
    elif "OUTPUT_FORMAT=" in line:
        output_format = line.split("=")[1].strip()
    elif "WRT_VOL_SOL=" in line:
        merge_vol = line.split("=")[1].strip()
    elif "WRT_SRF_SOL=" in line:
        merge_srf = line.split("=")[1].strip()


# -------------------------------------------------------------------
#  Tecplot Volume File Merging
# -------------------------------------------------------------------

if output_format == "TECPLOT":

    # Read the total number of points, elements, and parsing of the original Tecplot files
    if options.output != "False" and merge_vol != "NO":
        npoint = 0; nelem = 0;
        for domain in range(int(options.partitions)):
            if int(options.timestep) == -1:
                if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (volume_flow_filename, domain+1))
                elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (volume_adjoint_filename, domain+1))
            else:
                if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (volume_flow_filename, domain+1, options.timestep))
                elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (volume_adjoint_filename, domain+1, options.timestep))
            for line in input_file:
                if "TITLE" in line:
                    title = line
                elif "VARIABLES" in line:
                    variables_list = line
                elif "ZONE" in line:
                    if int(options.timestep) == -1:
                        npoint = npoint + int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                        nelem = nelem + int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
                        datapacking = line.replace("ZONE"," ").split(",")[2].split("=")[1]
                        kindzone = line.replace("ZONE"," ").split(",")[3].split("=")[1]
                    else:
                        global_ID = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                        global_time = line.replace("ZONE"," ").split(",")[1].split("=")[1]
                        npoint = npoint + int(line.replace("ZONE"," ").split(",")[2].split("=")[1])
                        nelem = nelem + int(line.replace("ZONE"," ").split(",")[3].split("=")[1])
                        datapacking = line.replace("ZONE"," ").split(",")[4].split("=")[1]
                        kindzone = line.replace("ZONE"," ").split(",")[5].split("=")[1]
            input_file.close()

        # Create an array with the solution, and connectivity
        Solution = []; Connectivity = []
        for ipoint in range(npoint):
            Solution.append("EMPTY")
        for ielem in range(nelem):
            Connectivity.append("EMPTY")

        npoint = 0; nelem = 0
        for domain in range(int(options.partitions)):

            percentage = int(100*float(domain)/(float(options.partitions)-1.0))
            stdout.write("\rMerging solution files (.plt)... %d%%" % percentage)
            stdout.flush()

            # Create the Local2Global vector (using the .su2 files)
            mesh_extsplit = os.path.splitext(mesh_filename)
            input_file = open("%s_%s%s" % (mesh_extsplit[0], domain+1, mesh_extsplit[1]))
            while 1:
                line = input_file.readline()
                if "NDIM" in line:
                    nDim = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
                if "NELEM" in line:
                    local_elem = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
                    for list_elems in range(local_elem):
                        line_elem = input_file.readline()
                if "NPOIN" in line:
                    local_point = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[0])
                    local_point_domain = int(line.replace("\t"," ").split("=")[1].strip().split(" ")[1])
                    npoint = npoint + local_point_domain
                    Local2Global = numpy.zeros(local_point, numpy.int)
                    for list_points in range(local_point):
                        line_point = input_file.readline()
                        LocalIndex = int(line_point.strip().replace("\t"," ").split(" ")[nDim])
                        GlobalIndex = int(line_point.strip().replace("\t"," ").split(" ")[nDim+1])
                        Local2Global[LocalIndex] = GlobalIndex
                if not line: break

            # Close the .su2 input file 
            input_file.close()

            # Open the .plt solution at each domain
            if int(options.timestep) == -1:
                if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (volume_flow_filename, domain+1))
                elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (volume_adjoint_filename, domain+1))
            else:
                if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (volume_flow_filename, domain+1, options.timestep))
                elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (volume_adjoint_filename, domain+1, options.timestep))

            # Three first lines are only the header
            solution = input_file.readline(); solution = input_file.readline(); solution = input_file.readline()

            # Save the solution using the right numbering
            for list_points in range(local_point):                
                Solution[Local2Global[list_points]] = input_file.readline()

            # Save the connectivity using the right numbering
            for list_elems in range(local_elem):
                ConnectivityLocal = input_file.readline().split()
                ConnectivityList = []
                AddElement = False
                for Nodes in ConnectivityLocal:
                    if (int(Nodes) <= local_point_domain): AddElement = True
                    else: AddElement = False
                    ConnectivityList.append(str((Local2Global[int(Nodes)-1])+1))
                if AddElement == True:
                    ConnectivityGlobal = ' '.join(ConnectivityList)
                    Connectivity[nelem] = ConnectivityGlobal
                    nelem = nelem + 1

            # Close the .plt input file 
            input_file.close()

        # Open the final (output) .plt solution
        if int(options.timestep) == -1:
            if math_problem == "DIRECT": output_file = open("%s.plt" % (volume_flow_filename) ,"w")
            elif math_problem == "ADJOINT": output_file = open("%s.plt" % (volume_adjoint_filename) ,"w")
        else:
            if math_problem == "DIRECT": output_file = open("%s_%s.plt" % (volume_flow_filename, options.timestep) ,"w")
            elif math_problem == "ADJOINT": output_file = open("%s_%s.plt" % (volume_adjoint_filename, options.timestep) ,"w")

        # Write the header (title, variables, and zone info)
        output_file.write("%s" % title)
        output_file.write("%s" % variables_list)
        if int(options.timestep) == -1:
            output_file.write("ZONE NODES=%s , ELEMENTS=%s, DATAPACKING=%s, ZONETYPE=%s" % (npoint, nelem, datapacking, kindzone))
        else:
            output_file.write("ZONE STRANDID=%s SOLUTIONTIME=%s NODES=%s , ELEMENTS=%s, DATAPACKING=%s, ZONETYPE=%s" % (global_ID, global_time, npoint, nelem, datapacking, kindzone))

        # Write the solution and connectivity (note that a renumbering is needed)
        for iPoint in range(npoint):
            output_file.write(Solution[iPoint])
        for iElem in range(nelem):
            output_file.write("%s \n" % Connectivity[iElem])

        # Close the .plt output file 
        output_file.close()

        # Remove the .plt input files
        for domain in range(int(options.partitions)):
            if int(options.timestep) == -1:
                if math_problem == "DIRECT": os.remove("%s_%s.plt" % (volume_flow_filename, domain+1))
                elif math_problem == "ADJOINT": os.remove("%s_%s.plt" % (volume_adjoint_filename, domain+1))
            else:
                if math_problem == "DIRECT": os.remove("%s_%s_%s.plt" % (volume_flow_filename, domain+1, options.timestep))
                elif math_problem == "ADJOINT": os.remove("%s_%s_%s.plt" % (volume_adjoint_filename, domain+1, options.timestep))

# -------------------------------------------------------------------
#  Tecplot Surface File Merging
# -------------------------------------------------------------------

    # Write the surface .plt file
    if merge_srf != "NO":
        npoint = 0; nelem = 0; addpoint = [ 0 ]
        for domain in range(int(options.partitions)):
            if int(options.timestep) == -1:
                if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
                elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
            else:
                if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, options.timestep))
                elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, options.timestep))         
            for line in input_file:
                if "VARIABLES" in line:
                    variables_list = line
                elif "ZONE" in line:
                    if int(options.timestep) == -1:
                        local_point = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                        npoint = npoint + local_point
                        if domain == 0: addpoint.append(local_point)
                        else: addpoint.append(local_point+addpoint[domain])
                        nelem = nelem + int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
                        datapacking = line.replace("ZONE"," ").split(",")[2].split("=")[1]
                        kindzone = line.replace("ZONE"," ").split(",")[3].split("=")[1]
                    else:
                        global_ID = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                        global_time = line.replace("ZONE"," ").split(",")[1].split("=")[1]
                        local_point = int(line.replace("ZONE"," ").split(",")[2].split("=")[1])
                        npoint = npoint + local_point
                        if domain == 0: addpoint.append(local_point)
                        else: addpoint.append(local_point+addpoint[domain])
                        nelem = nelem + int(line.replace("ZONE"," ").split(",")[3].split("=")[1])
                        datapacking = line.replace("ZONE"," ").split(",")[4].split("=")[1]
                        kindzone = line.replace("ZONE"," ").split(",")[5].split("=")[1]
            input_file.close()

        # Open the output file
        if int(options.timestep) == -1:
            if math_problem == "DIRECT": output_file = open("%s.plt" % (surface_flow_filename) ,"w")
            elif math_problem == "ADJOINT": output_file = open("%s.plt" % (surface_adjoint_filename) ,"w")
        else:
            if math_problem == "DIRECT": output_file = open("%s_%s.plt" % (surface_flow_filename, options.timestep) ,"w")
            elif math_problem == "ADJOINT": output_file = open("%s_%s.plt" % (surface_adjoint_filename, options.timestep) ,"w")

        # Write the header
        output_file.write("TITLE = \"Visualization of the surface grid\"\n")

        # Write the coordinates and solution
        output_file.write("%s" % variables_list)
        if int(options.timestep) == -1:
            output_file.write("ZONE NODES=%s , ELEMENTS=%s, DATAPACKING=%s, ZONETYPE=%s" % (npoint, nelem, datapacking, kindzone))
        else:
            output_file.write("ZONE STRANDID=%s SOLUTIONTIME=%s NODES=%s , ELEMENTS=%s, DATAPACKING=%s, ZONETYPE=%s" % (global_ID, global_time, npoint, nelem, datapacking, kindzone))
        for domain in range(int(options.partitions)):
            if int(options.timestep) == -1:
                if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
                elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
            else:
                if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, options.timestep))
                elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, options.timestep))  
            while 1:
                line = input_file.readline()
                if "ZONE" in line:
                    if int(options.timestep) == -1:
                        ipoint = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                    else:
                        ipoint = int(line.replace("ZONE"," ").split(",")[2].split("=")[1])
                    for list_points in range(int(ipoint)):
                        new_line = input_file.readline()
                        output_file.write("%s" %  new_line)
                if not line: break
            input_file.close()

        # Write the elements
        for domain in range(int(options.partitions)):
            if int(options.timestep) == -1:
                if math_problem == "DIRECT": input_file = open("%s_%s.plt" % (surface_flow_filename, domain+1))
                elif math_problem == "ADJOINT": input_file = open("%s_%s.plt" % (surface_adjoint_filename, domain+1))
            else:
                if math_problem == "DIRECT": input_file = open("%s_%s_%s.plt" % (surface_flow_filename, domain+1, options.timestep))
                elif math_problem == "ADJOINT": input_file = open("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, options.timestep))
            while 1:
                line = input_file.readline()
                if "ZONE" in line:
                    if int(options.timestep) == -1:
                        ipoint = int(line.replace("ZONE"," ").split(",")[0].split("=")[1])
                        ielem = int(line.replace("ZONE"," ").split(",")[1].split("=")[1])
                    else:
                        ipoint = int(line.replace("ZONE"," ").split(",")[2].split("=")[1])
                        ielem = int(line.replace("ZONE"," ").split(",")[3].split("=")[1])
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
            if int(options.timestep) == -1:
                if math_problem == "DIRECT": os.remove("%s_%s.plt" % (surface_flow_filename, domain+1))
                elif math_problem == "ADJOINT": os.remove("%s_%s.plt" % (surface_adjoint_filename, domain+1))
            else:
                if math_problem == "DIRECT":  os.remove("%s_%s_%s.plt" % (surface_flow_filename, domain+1, options.timestep))
                elif math_problem == "ADJOINT":  os.remove("%s_%s_%s.plt" % (surface_adjoint_filename, domain+1, options.timestep))

        stdout.write("\n")


# -------------------------------------------------------------------
#  Paraview File Merging
# -------------------------------------------------------------------

if output_format == "PARAVIEW":

    # Read the total number of points, and elements
    if math_problem == "DIRECT":
        if options.output != "False" and merge_vol != "NO":
            npoint = 0
            nelem = 0
            nstore = 0
            cont = 0
            addpoint = [ 0 ]
            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    input_file = open("%s_%s.vtk" % (volume_flow_filename, domain+1))
                else:
                    input_file = open("%s_%s_%s.vtk" % (volume_flow_filename, domain+1, int(options.timestep)))
                for line in input_file:
                    if "POINTS" in line:
                        local_point = int(line.replace("\t"," ").split(" ")[1])
                        npoint = npoint + local_point             
                        if domain == 0: addpoint.append(local_point)
                        else: addpoint.append(local_point+addpoint[domain])
                    elif "CELLS" in line:
                        nelem = nelem + int(line.replace("\t"," ").split(" ")[1])
                        nstore = nstore + int(line.replace("\t"," ").split(" ")[2])
                input_file.close()

        # Open the output file
            if int(options.timestep) == -1:
                output_file = open("%s.vtk" % (volume_flow_filename) ,"w")
            else:
                output_file = open("%s_%s.vtk" % (volume_flow_filename, int(options.timestep)) ,"w")

        # Write the header
            output_file.write("# vtk DataFile Version 2.0\n")
            output_file.write("Value of the flow variables\n")
            output_file.write("ASCII\n")
            output_file.write("DATASET UNSTRUCTURED_GRID\n")

        # Write the coordinates
            output_file.write("POINTS %s float\n" % npoint)
            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    input_file = open("%s_%s.vtk" % (volume_flow_filename, domain+1))
                else:
                    input_file = open("%s_%s_%s.vtk" % (volume_flow_filename, domain+1, int(options.timestep)))
                while 1:
                    line = input_file.readline()
                    if "POINTS" in line:
                        ipoint = line.replace("\t"," ").split(" ")[1]
                        for list_points in range(int(ipoint)):
                            new_line = input_file.readline()
                            output_file.write("%s" %  new_line)
                    if not line: break
                input_file.close()

        # Write the elements
            output_file.write("CELLS %s %s\n" % (nelem, nstore))
            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    input_file = open("%s_%s.vtk" % (volume_flow_filename, domain+1))
                else:
                    input_file = open("%s_%s_%s.vtk" % (volume_flow_filename, domain+1, int(options.timestep)))
                while 1:
                    line = input_file.readline()
                    if "CELLS" in line:
                        ielem = line.replace("\t"," ").split(" ")[1]
                        for list_elems in range(int(ielem)):
                            elem_info = input_file.readline().replace("\t"," ").split(" ")
                            cont = 0
                            for iter in elem_info:
                                if cont == 0: elem_info[cont] = ("%s" % int(iter))
                                else: elem_info[cont] = ("%s" % (int(iter)+addpoint[domain]))
                                cont = cont+1;
                            output_file.write("%s\n" %  ("\t".join(elem_info)))
                    if not line: break
                input_file.close()

        # Write the cell types
            output_file.write("CELL_TYPES %s\n" % nelem)
            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    input_file = open("%s_%s.vtk" % (volume_flow_filename, domain+1))
                else:
                    input_file = open("%s_%s_%s.vtk" % (volume_flow_filename, domain+1, int(options.timestep)))
                while 1:
                    line = input_file.readline()
                    if "CELL_TYPES" in line:
                        icell = line.replace("\t"," ").split(" ")[1]
                        for list_cells in range(int(icell)):
                            new_line = input_file.readline()
                            output_file.write("%s" %  new_line)
                    if not line: break
                input_file.close()

        # Write number of points
            output_file.write("POINT_DATA %s\n" % npoint)

        # Write the flow variables
            variables = ["Density", "Density_x_Energy", "Pressure", "Mach_Number"]
            if incompressible_formulation == "YES":
                variables = ["Density", "Pressure"]
            for var in variables:
                output_file.write("SCALARS %s float 1\n" % var)
                output_file.write("LOOKUP_TABLE default\n")
                for domain in range(int(options.partitions)):
                    if int(options.timestep) == -1:
                        input_file = open("%s_%s.vtk" % (volume_flow_filename, domain+1))
                    else:
                        input_file = open("%s_%s_%s.vtk" % (volume_flow_filename, domain+1, int(options.timestep)))
                    while 1:
                        line = input_file.readline()
                        if "POINT_DATA" in line:
                            ipoint = line.replace("\t"," ").split(" ")[1]
                        if ("SCALARS %s float 1" % var) in line:
                            input_file.readline()
                            for list_points in range(int(ipoint)):
                                output_file.write("%s" %  (input_file.readline()))
                        if not line: break
                    input_file.close()

            variables = ["Velocity", "Density_x_Velocity"]
            if incompressible_formulation == "YES":
                variables = ["Velocity"]
            for var in variables:
                output_file.write("VECTORS %s float\n" % var)
                for domain in range(int(options.partitions)):
                    if int(options.timestep) == -1:
                        input_file = open("%s_%s.vtk" % (volume_flow_filename, domain+1))
                    else:
                        input_file = open("%s_%s_%s.vtk" % (volume_flow_filename, domain+1, int(options.timestep)))
                    while 1:
                        line = input_file.readline()
                        if "POINT_DATA" in line:
                            ipoint = line.replace("\t"," ").split(" ")[1]
                        if ("VECTORS %s float" % var) in line:
                            for list_points in range(int(ipoint)):
                                output_file.write("%s" %  (input_file.readline()))
                        if not line: break
                    input_file.close()

            if physical_problem == "NAVIER_STOKES":

                if turb_model == "NONE":
                    variables = ["Temperature", "Laminar_Viscosity"]
                else:
                    variables = ["Temperature", "Laminar_Viscosity", "Eddy_Viscosity", "Nu_Tilde"]

                if incompressible_formulation == "YES":
                    variables = ["Viscosity"]

                for var in variables:
                    output_file.write("SCALARS %s float 1\n" % var)
                    output_file.write("LOOKUP_TABLE default\n")
                    for domain in range(int(options.partitions)):
                        if int(options.timestep) == -1:
                            input_file = open("%s_%s.vtk" % (volume_flow_filename, domain+1))
                        else:
                            input_file = open("%s_%s_%s.vtk" % (volume_flow_filename, domain+1, int(options.timestep)))     
                        while 1:
                            line = input_file.readline()
                            if "POINT_DATA" in line:
                                ipoint = line.replace("\t"," ").split(" ")[1]
                            if ("SCALARS %s float 1" % var) in line:
                                input_file.readline()
                                for list_points in range(int(ipoint)):
                                    output_file.write("%s" %  (input_file.readline()))
                            if not line: break
                        input_file.close()


            output_file.close()

            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    os.remove("%s_%s.vtk" % (volume_flow_filename, domain+1))
                else:
                    os.remove("%s_%s_%s.vtk" % (volume_flow_filename, domain+1, int(options.timestep)))

    # Surface file
        if merge_srf != "NO":
            npoint = 0
            nelem = 0
            nstore = 0
            cont = 0
            addpoint = [ 0 ]
            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    input_file = open("%s_%s.vtk" % (surface_flow_filename, domain+1))
                else:
                    input_file = open("%s_%s_%s.vtk" % (surface_flow_filename, domain+1, int(options.timestep)))         
                for line in input_file:
                    if "POINTS" in line:
                        local_point = int(line.replace("\t"," ").split(" ")[1])
                        npoint = npoint + local_point             
                        if domain == 0: addpoint.append(local_point)
                        else: addpoint.append(local_point+addpoint[domain])
                    elif "CELLS" in line:
                        nelem = nelem + int(line.replace("\t"," ").split(" ")[1])
                        nstore = nstore + int(line.replace("\t"," ").split(" ")[2])
                input_file.close()

        # Open the output file
            if int(options.timestep) == -1:
                output_file = open("%s.vtk" % (surface_flow_filename) ,"w")
            else:
                output_file = open("%s_%s.vtk" % (surface_flow_filename, int(options.timestep)) ,"w")

        # Write the header
            output_file.write("# vtk DataFile Version 2.0\n")
            output_file.write("Value of the surface flow variables\n")
            output_file.write("ASCII\n")
            output_file.write("DATASET UNSTRUCTURED_GRID\n")

        # Write the coordinates
            output_file.write("POINTS %s float\n" % npoint)
            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    input_file = open("%s_%s.vtk" % (surface_flow_filename, domain+1))
                else:
                    input_file = open("%s_%s_%s.vtk" % (surface_flow_filename, domain+1, int(options.timestep)))  
                while 1:
                    line = input_file.readline()
                    if "POINTS" in line:
                        ipoint = line.replace("\t"," ").split(" ")[1]
                        for list_points in range(int(ipoint)):
                            new_line = input_file.readline()
                            output_file.write("%s" %  new_line)
                    if not line: break
                input_file.close()

        # Write the elements
            output_file.write("CELLS %s %s\n" % (nelem, nstore))
            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    input_file = open("%s_%s.vtk" % (surface_flow_filename, domain+1))
                else:
                    input_file = open("%s_%s_%s.vtk" % (surface_flow_filename, domain+1, int(options.timestep)))  
                while 1:
                    line = input_file.readline()
                    if "CELLS" in line:
                        ielem = line.replace("\t"," ").split(" ")[1]
                        for list_elems in range(int(ielem)):
                            elem_info = input_file.readline().rstrip("\n").rstrip().replace("\t",",").split(",")
                            cont = 0
                            for iter in elem_info:
                                if cont == 0: elem_info[cont] = ("%s" % int(iter))
                                else: elem_info[cont] = ("%s" % (int(iter)+addpoint[domain]))
                                cont = cont+1;
                            output_file.write("%s\n" %  ("\t".join(elem_info)))
                    if not line: break
                input_file.close()

        # Write the cell types
            output_file.write("CELL_TYPES %s\n" % nelem)
            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    input_file = open("%s_%s.vtk" % (surface_flow_filename, domain+1))
                else:
                    input_file = open("%s_%s_%s.vtk" % (surface_flow_filename, domain+1, int(options.timestep)))  
                while 1:
                    line = input_file.readline()
                    if "CELL_TYPES" in line:
                        icell = line.replace("\t"," ").split(" ")[1]
                        for list_cells in range(int(icell)):
                            new_line = input_file.readline()
                            output_file.write("%s" %  new_line)
                    if not line: break
                input_file.close()

        # Write number of points
            output_file.write("POINT_DATA %s\n" % npoint)

        # Write the flow variables
            if physical_problem == "EULER":
                variables = ["Pressure_Coefficient", "Mach_Number"]
            if physical_problem == "NAVIER_STOKES":
                variables = ["Pressure_Coefficient", "|Skin_Friction_Coefficient|"]
            for var in variables:
                output_file.write("SCALARS %s float 1\n" % var)
                output_file.write("LOOKUP_TABLE default\n")
                for domain in range(int(options.partitions)):
                    if int(options.timestep) == -1:
                        input_file = open("%s_%s.vtk" % (surface_flow_filename, domain+1))
                    else:
                        input_file = open("%s_%s_%s.vtk" % (surface_flow_filename, domain+1, int(options.timestep)))  
                    while 1:
                        line = input_file.readline()
                        if "POINT_DATA" in line:
                            ipoint = line.replace("\t"," ").split(" ")[1]
                        if ("SCALARS %s float 1" % var) in line:
                            input_file.readline()
                            for list_points in range(int(ipoint)):
                                output_file.write("%s" %  (input_file.readline()))
                        if not line: break
                    input_file.close()

            output_file.close()

            for domain in range(int(options.partitions)):
                if int(options.timestep) == -1:
                    os.remove("%s_%s.vtk" % (surface_flow_filename, domain+1))
                else:
                    os.remove("%s_%s_%s.vtk" % (surface_flow_filename, domain+1, int(options.timestep)))

    # Read the total number of points, and elements
    if math_problem == "ADJOINT":
        if options.output != "False" and merge_vol != "NO":
            npoint = 0
            nelem = 0
            nstore = 0
            cont = 0
            addpoint = [ 0 ]
            for domain in range(int(options.partitions)):
                input_file = open("%s_%s.vtk" % (volume_adjoint_filename, domain+1))
                for line in input_file:
                    if "POINTS" in line:
                        local_point = int(line.replace("\t"," ").split(" ")[1])
                        npoint = npoint + local_point             
                        if domain == 0: addpoint.append(local_point)
                        else: addpoint.append(local_point+addpoint[domain])
                    elif "CELLS" in line:
                        nelem = nelem + int(line.replace("\t"," ").split(" ")[1])
                        nstore = nstore + int(line.replace("\t"," ").split(" ")[2])
                input_file.close()

            output_file = open("adjoint.vtk","w")
        # Write the header
            output_file.write("# vtk DataFile Version 2.0\n")
            output_file.write("Value of the adjoint variables\n")
            output_file.write("ASCII\n")
            output_file.write("DATASET UNSTRUCTURED_GRID\n")

        # Write the coordinates
            output_file.write("POINTS %s float\n" % npoint)
            for domain in range(int(options.partitions)):
                input_file = open("%s_%s.vtk" % (volume_adjoint_filename, domain+1))
                while 1:
                    line = input_file.readline()
                    if "POINTS" in line:
                        ipoint = line.replace("\t"," ").split(" ")[1]
                        for list_points in range(int(ipoint)):
                            new_line = input_file.readline()
                            output_file.write("%s" %  new_line)
                    if not line: break
                input_file.close()

        # Write the elements
            output_file.write("CELLS %s %s\n" % (nelem, nstore))
            for domain in range(int(options.partitions)):
                input_file = open("%s_%s.vtk" % (volume_adjoint_filename, domain+1))
                while 1:
                    line = input_file.readline()
                    if "CELLS" in line:
                        ielem = line.replace("\t"," ").split(" ")[1]
                        for list_elems in range(int(ielem)):
                            elem_info = input_file.readline().replace("\t"," ").split(" ")
                            cont = 0
                            for iter in elem_info:
                                if cont == 0: elem_info[cont] = ("%s" % int(iter))
                                else: elem_info[cont] = ("%s" % (int(iter)+addpoint[domain]))
                                cont = cont+1;
                            output_file.write("%s\n" %  ("\t".join(elem_info)))
                    if not line: break
                input_file.close()

        # Write the cell types
            output_file.write("CELL_TYPES %s\n" % nelem)
            for domain in range(int(options.partitions)):
                input_file = open("%s_%s.vtk" % (volume_adjoint_filename, domain+1))
                while 1:
                    line = input_file.readline()
                    if "CELL_TYPES" in line:
                        icell = line.replace("\t"," ").split(" ")[1]
                        for list_cells in range(int(icell)):
                            new_line = input_file.readline()
                            output_file.write("%s" %  new_line)
                    if not line: break
                input_file.close()

        # Write number of points
            output_file.write("POINT_DATA %s\n" % npoint)

        # Write the adjoint variables
            variables = ["PsiRho", "PsiE"]
            for var in variables:
                output_file.write("SCALARS %s float 1\n" % var)
                output_file.write("LOOKUP_TABLE default\n")
                for domain in range(int(options.partitions)):
                    input_file = open("%s_%s.vtk" % (volume_adjoint_filename, domain+1))
                    while 1:
                        line = input_file.readline()
                        if "POINT_DATA" in line:
                            ipoint = line.replace("\t"," ").split(" ")[1]
                        if ("SCALARS %s float 1" % var) in line:
                            input_file.readline()
                            for list_points in range(int(ipoint)):
                                output_file.write("%s" %  (input_file.readline()))
                        if not line: break
                    input_file.close()

            variables = ["Phi"]
            for var in variables:
                output_file.write("VECTORS %s float\n" % var)
                for domain in range(int(options.partitions)):
                    input_file = open("%s_%s.vtk" % (volume_adjoint_filename, domain+1))
                    while 1:
                        line = input_file.readline()
                        if "POINT_DATA" in line:
                            ipoint = line.replace("\t"," ").split(" ")[1]
                        if ("VECTORS %s float" % var) in line:
                            for list_points in range(int(ipoint)):
                                output_file.write("%s" %  (input_file.readline()))
                        if not line: break
                    input_file.close()

            output_file.close()

            for domain in range(int(options.partitions)):
                os.remove("%s_%s.vtk" % (volume_adjoint_filename, domain+1))

    # Surface file
        if merge_srf != "NO":
            npoint = 0
            nelem = 0
            nstore = 0
            cont = 0
            addpoint = [ 0 ]
            for domain in range(int(options.partitions)):
                input_file = open("%s_%s.vtk" % (surface_adjoint_filename, domain+1))
                for line in input_file:
                    if "POINTS" in line:
                        local_point = int(line.replace("\t"," ").split(" ")[1])
                        npoint = npoint + local_point             
                        if domain == 0: addpoint.append(local_point)
                        else: addpoint.append(local_point+addpoint[domain])
                    elif "CELLS" in line:
                        nelem = nelem + int(line.replace("\t"," ").split(" ")[1])
                        nstore = nstore + int(line.replace("\t"," ").split(" ")[2])
                input_file.close()

            output_file = open("surface_adjoint.vtk","w")
        # Write the header
            output_file.write("# vtk DataFile Version 2.0\n")
            output_file.write("Value of the surface adjoint variables\n")
            output_file.write("ASCII\n")
            output_file.write("DATASET UNSTRUCTURED_GRID\n")

        # Write the coordinates
            output_file.write("POINTS %s float\n" % npoint)
            for domain in range(int(options.partitions)):
                input_file = open("%s_%s.vtk" % (surface_adjoint_filename, domain+1))
                while 1:
                    line = input_file.readline()
                    if "POINTS" in line:
                        ipoint = line.replace("\t"," ").split(" ")[1]
                        for list_points in range(int(ipoint)):
                            new_line = input_file.readline()
                            output_file.write("%s" %  new_line)
                    if not line: break
                input_file.close()

        # Write the elements
            output_file.write("CELLS %s %s\n" % (nelem, nstore))
            for domain in range(int(options.partitions)):
                input_file = open("%s_%s.vtk" % (surface_adjoint_filename, domain+1))
                while 1:
                    line = input_file.readline()
                    if "CELLS" in line:
                        ielem = line.replace("\t"," ").split(" ")[1]
                        for list_elems in range(int(ielem)):
                            elem_info = input_file.readline().rstrip("\n").rstrip().replace("\t",",").split(",")
                            cont = 0
                            for iter in elem_info:
                                if cont == 0: elem_info[cont] = ("%s" % int(iter))
                                else: elem_info[cont] = ("%s" % (int(iter)+addpoint[domain]))
                                cont = cont+1;
                            output_file.write("%s\n" %  ("\t".join(elem_info)))
                    if not line: break
                input_file.close()

        # Write the cell types
            output_file.write("CELL_TYPES %s\n" % nelem)
            for domain in range(int(options.partitions)):
                input_file = open("%s_%s.vtk" % (surface_adjoint_filename, domain+1))
                while 1:
                    line = input_file.readline()
                    if "CELL_TYPES" in line:
                        icell = line.replace("\t"," ").split(" ")[1]
                        for list_cells in range(int(icell)):
                            new_line = input_file.readline()
                            output_file.write("%s" %  new_line)
                    if not line: break
                input_file.close()

        # Write number of points
            output_file.write("POINT_DATA %s\n" % npoint)

        # Write the adjoint variables
            variables = ["Shape_Sensitivity", "PsiRho"]
            for var in variables:
                output_file.write("SCALARS %s float 1\n" % var)
                output_file.write("LOOKUP_TABLE default\n")
                for domain in range(int(options.partitions)):
                    input_file = open("%s_%s.vtk" % (surface_adjoint_filename, domain+1))
                    while 1:
                        line = input_file.readline()
                        if "POINT_DATA" in line:
                            ipoint = line.replace("\t"," ").split(" ")[1]
                        if ("SCALARS %s float 1" % var) in line:
                            input_file.readline()
                            for list_points in range(int(ipoint)):
                                output_file.write("%s" %  (input_file.readline()))
                        if not line: break
                    input_file.close()

            output_file.close()

            for domain in range(int(options.partitions)):
                os.remove("%s_%s.vtk" % (surface_adjoint_filename, domain+1))

        stdout.write("\n")
