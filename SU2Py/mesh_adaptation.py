#!/usr/bin/python

## \file mesh_adaptation.py
#  \brief Python script for doing the grid adaptation using the SU2 suite.
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


import os, time
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-c", "--cycle", dest="cycle", default=1,
                  help="number of CYCLE adaptations", metavar="CYCLE")
parser.add_option("-o", "--overwrite", dest="overwrite", default=False,
                  help="OVERWRITE_MESH the output mesh with the adapted one", metavar="OVERWRITE_MESH")

(options, args)=parser.parse_args()

if int(options.partitions) == 1:
    parallel = False
else : parallel = True
cycles = int(options.cycle)

# Read the kind of grid adaptation and number of cycles
for line in open(options.filename):
    if "KIND_ADAPT=" in line:
        adapt_kind = line.split("=")[1].strip()

# General parameters
Config_MAC_file = "config_MAC_" + options.filename
Config_CFD_file = "config_CFD_" + options.filename
Mesh_MAC_file = "mesh_MAC_" + options.filename.replace(".cfg",".su2")

# Parameters for the grid adaptation
mesh_finest = "mesh_finest.su2"
restart_flow_finest = "restart_flow_finest.dat"
restart_lin_finest = "restart_lin_finest.dat"
restart_adj_finest = "restart_adj_finest.dat"

for iAdaptCycle in range(int(cycles)):

    # Change the parameters to do the first flow simulation
    # In the first iteration we store from the config file: 
    #     -The name of the restart file (restart_flow_file, restart_adj_file, restart_lin_file).
    #     -The kind of objective function (objfunc).
    #     -The name of the original grid (original_mesh_file).
    output_file = open(Config_CFD_file,"w")
    for line in open(options.filename):
        if "MATH_PROBLEM=" in line:
            output_file.write("MATH_PROBLEM= DIRECT \n")
        elif "CADJ_OBJFUNC=" in line:
            if iAdaptCycle == 0: objfunc = line.split("=")[1].strip()
            output_file.write(line)
        elif "RESTART_SOL=" in line:
            if iAdaptCycle == 0: output_file.write(line)
            else: output_file.write("RESTART_SOL= YES \n")
        elif "RESTART_FLOW_FILENAME=" in line:
            if iAdaptCycle == 0: restart_flow_file = line.split("=")[1].strip()
            output_file.write(line)
        elif "RESTART_ADJ_FILENAME=" in line:
            if iAdaptCycle == 0: restart_adj_file = line.split("=")[1].strip()
            output_file.write(line)
        elif "RESTART_LIN_FILENAME=" in line:
            if iAdaptCycle == 0: restart_lin_file = line.split("=")[1].strip()
            output_file.write(line)
        elif "SOLUTION_FLOW_FILENAME=" in line:
            if iAdaptCycle == 0: output_file.write(line)
            else: output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_file)
        elif "MESH_FILENAME=" in line:
            if iAdaptCycle == 0:
                original_mesh_file = line.split("=")[1].strip()
                output_file.write(line)
            else: output_file.write("MESH_FILENAME= %s \n" % Mesh_MAC_file)
        else: output_file.write(line)

    output_file.close()
    if parallel: os.system ("parallel_computation.py -p %s -f %s -o True" % (int(options.partitions), Config_CFD_file))
    else: os.system ("SU2_CFD " + Config_CFD_file)

    # Change the parameters to do the first adjoint simulation
    # In this case we don't store anything, at the first iteration we use the
    # filenames of the original cofiguration file.
    if adapt_kind == "GRAD_ADJOINT" or adapt_kind == "ROBUST" or adapt_kind == "COMPUTABLE_ROBUST" or adapt_kind == "COMPUTABLE" or adapt_kind == "REMAINING":
        output_file = open(Config_CFD_file,"w")
        for line in open(options.filename):
            if "MATH_PROBLEM=" in line:
                output_file.write("MATH_PROBLEM= ADJOINT \n")
            elif "RESTART_SOL=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("RESTART_SOL= YES \n")
            elif "SOLUTION_FLOW_FILENAME=" in line:
                output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_file)
            elif "SOLUTION_ADJ_FILENAME=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("SOLUTION_ADJ_FILENAME= %s \n" % restart_adj_file)
            elif "MESH_FILENAME=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("MESH_FILENAME= %s \n" % Mesh_MAC_file)
            else: output_file.write(line)

        output_file.close()
        if parallel: os.system ("parallel_computation.py -p %s -f %s -o True" % (int(options.partitions), Config_CFD_file))
        else: os.system ("SU2_CFD " + Config_CFD_file)

    # Change the parameters to do the first linear simulation (in case it is necessary)
    # In this case we don't store anything, at the first iteration we use the
    # filenames of the original cofiguration file.
    if adapt_kind == "COMPUTABLE_ROBUST":
        output_file = open(Config_CFD_file,"w")
        for line in open(options.filename):
            if "MATH_PROBLEM=" in line:
                output_file.write("MATH_PROBLEM= LINEARIZED \n")
            elif "RESTART_SOL=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("RESTART_SOL= YES \n")
            elif "SOLUTION_FLOW_FILENAME=" in line:
                output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_file)
            elif "SOLUTION_LIN_FILENAME=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("SOLUTION_LIN_FILENAME= %s \n" % restart_lin_file)
            elif "MESH_FILENAME=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("MESH_FILENAME= %s \n" % Mesh_MAC_file)
            else: output_file.write(line)

        output_file.close()
        if parallel: os.system ("parallel_computation.py -p %s -f %s -o True" % (int(options.partitions), Config_CFD_file))
        else: os.system ("SU2_CFD " + Config_CFD_file)

    # Change the parameters to do a direct and adjoint iteration over a fine grid
    # We will use the following files names:
    #    - mesh_finest
    #    - restart_flow_finest
    #    - restart_adj_finest
    #    - restart_lin_finest

    if ((adapt_kind == "ROBUST" or adapt_kind == "COMPUTABLE" or adapt_kind == "COMPUTABLE_ROBUST" or adapt_kind == "REMAINING") and (iAdaptCycle < int(cycles)-1 or int(cycles) == 1)) :

        # Create the fine grid and extrapolate the flow solution
        # from the coarse to the fine grid and store the extrapolated solution.
        output_file = open(Config_MAC_file,"w")
        for line in open(options.filename):
            if "KIND_ADAPT=" in line:
                output_file.write("KIND_ADAPT= FULL_FLOW \n")
            elif "SOLUTION_FLOW_FILENAME=" in line:
                output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_file)
            elif "RESTART_FLOW_FILENAME=" in line:
                output_file.write("RESTART_FLOW_FILENAME= %s \n" % restart_flow_finest)
            elif "MESH_FILENAME=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("MESH_FILENAME= %s \n" % Mesh_MAC_file)
            elif "MESH_OUT_FILENAME=" in line:
                output_file.write("MESH_OUT_FILENAME= %s \n" % mesh_finest)
            else: output_file.write(line)

        output_file.close()
        os.system ("SU2_MAC " + Config_MAC_file)

        # Create the fine grid and extrapolate the adjoint solution
        # from the coarse to the fine grid and store the extrapolated solution.
        output_file = open(Config_MAC_file,"w")
        for line in open(options.filename):
            if "KIND_ADAPT=" in line:
                output_file.write("KIND_ADAPT= FULL_ADJOINT \n")
            elif "SOLUTION_FLOW_FILENAME=" in line:
                output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_file)
            elif "SOLUTION_ADJ_FILENAME=" in line:
                output_file.write("SOLUTION_ADJ_FILENAME= %s \n" % restart_adj_file)
            elif "RESTART_FLOW_FILENAME=" in line:
                output_file.write("RESTART_FLOW_FILENAME= %s \n" % restart_flow_finest)
            elif "RESTART_ADJ_FILENAME=" in line:
                output_file.write("RESTART_ADJ_FILENAME= %s \n" % restart_adj_finest)
            elif "MESH_FILENAME=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("MESH_FILENAME= %s \n" % Mesh_MAC_file)
            elif "MESH_OUT_FILENAME=" in line:
                output_file.write("MESH_OUT_FILENAME= %s \n" % mesh_finest)
            else: output_file.write(line)

        output_file.close()
        os.system ("SU2_MAC " + Config_MAC_file)

        # Create the fine grid and extrapolate the linear solution
        # from the coarse to the fine grid and store the extrapolated solution.
        if adapt_kind == "COMPUTABLE_ROBUST":
            output_file = open(Config_MAC_file,"w")
            for line in open(options.filename):
                if "KIND_ADAPT=" in line:
                    output_file.write("KIND_ADAPT= FULL_LINEAR \n")
                elif "SOLUTION_FLOW_FILENAME=" in line:
                    output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_file)
                elif "SOLUTION_LIN_FILENAME=" in line:
                    output_file.write("SOLUTION_LIN_FILENAME= %s \n" % restart_lin_file)
                elif "RESTART_LIN_FILENAME=" in line:
                    output_file.write("RESTART_LIN_FILENAME= %s \n" % restart_lin_finest)
                elif "MESH_FILENAME=" in line:
                    if iAdaptCycle == 0: output_file.write(line)
                    else: output_file.write("MESH_FILENAME= %s \n" % Mesh_MAC_file)
                elif "MESH_OUT_FILENAME=" in line:
                    output_file.write("MESH_OUT_FILENAME= %s \n" % mesh_finest)
                else: output_file.write(line)

            output_file.close()
            os.system ("SU2_MAC " + Config_MAC_file)

        # Change the parameters to do one iteration of the flow solver on the finest grid
        # Allways restart with the interpolated solution and store the residual in the 
        # solution file for the fines grid.
        # Never use multigrid or other acceleration technique
        output_file = open(Config_CFD_file,"w")
        for line in open(options.filename):
            if "MATH_PROBLEM=" in line:
                output_file.write("MATH_PROBLEM= DIRECT \n")
            elif "EXT_ITER=" in line:
                output_file.write("EXT_ITER= 2 \n")
            elif "RESTART_SOL=" in line:
                output_file.write("RESTART_SOL= YES \n")
            elif "SOLUTION_FLOW_FILENAME=" in line:
                output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_finest)
            elif "STORE_RESIDUAL=" in line:
                output_file.write("STORE_RESIDUAL= YES \n")
            elif "RESTART_FLOW_FILENAME=" in line:
                output_file.write("RESTART_FLOW_FILENAME= %s \n" % restart_flow_finest)
            elif "MESH_FILENAME=" in line:
                output_file.write("MESH_FILENAME= %s \n" % mesh_finest)
            elif "FULLMG=" in line:
                output_file.write("FULLMG= NO \n")
            elif "MGLEVEL=" in line:
                output_file.write("MGLEVEL= 0 \n")
            elif "MGCYCLE=" in line:
                output_file.write("MGCYCLE= 0 \n")
            elif "MG_PRE_SMOOTH=" in line:
                output_file.write("MG_PRE_SMOOTH= ( 0 )\n")
            elif "MG_POST_SMOOTH=" in line:
                output_file.write("MG_POST_SMOOTH= ( 0 )\n")
            elif "MG_CORRECTION_SMOOTH=" in line:
                output_file.write("MG_CORRECTION_SMOOTH= ( 0 )\n")
            else: output_file.write(line)

        output_file.close()
        if parallel: os.system ("parallel_computation.py -p %s -f %s -o True" % (int(options.partitions), Config_CFD_file))
        else: os.system ("SU2_CFD " + Config_CFD_file)

        # Change the parameters to do one iteration of the adjoint solver on the finest grid
        # Allways restart with the interpolated solution and store the residual in the 
        # solution file for the finest grid. Use the flow solution of the finest grid
        # Never use multigrid or other acceleration technique
        output_file = open(Config_CFD_file,"w")
        for line in open(options.filename):
            if "MATH_PROBLEM=" in line:
                output_file.write("MATH_PROBLEM= ADJOINT \n")
            elif "EXT_ITER=" in line:
                output_file.write("EXT_ITER= 2 \n")
            elif "RESTART_SOL=" in line:
                output_file.write("RESTART_SOL= YES \n")
            elif "SOLUTION_FLOW_FILENAME=" in line:
                output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_finest)
            elif "SOLUTION_ADJ_FILENAME=" in line:
                output_file.write("SOLUTION_ADJ_FILENAME= %s \n" % restart_adj_finest)
            elif "MESH_FILENAME=" in line:
                output_file.write("MESH_FILENAME= %s \n" %  mesh_finest)
            elif "FULLMG=" in line:
                output_file.write("FULLMG= NO \n")
            elif "MGLEVEL=" in line:
                output_file.write("MGLEVEL= 0 \n")
            elif "MGCYCLE=" in line:
                output_file.write("MGCYCLE= 0 \n")
            elif "MG_PRE_SMOOTH=" in line:
                output_file.write("MG_PRE_SMOOTH= ( 0 )\n")
            elif "MG_POST_SMOOTH=" in line:
                output_file.write("MG_POST_SMOOTH= ( 0 )\n")
            elif "MG_CORRECTION_SMOOTH=" in line:
                output_file.write("MG_CORRECTION_SMOOTH= ( 0 )\n")
            else: output_file.write(line)
        output_file.close()
        if parallel: os.system ("parallel_computation.py -p %s -f %s -o True" % (int(options.partitions), Config_CFD_file))
        else: os.system ("SU2_CFD " + Config_CFD_file)

        # Change the parameters to do one iteration of the linear solver on the finest grid
        # Allways restart with the interpolated solution and store the residual in the 
        # solution file for the finest grid. Use the flow solution of the finest grid
        # Never use multigrid or other acceleration technique
        if adapt_kind == "COMPUTABLE_ROBUST":
            output_file = open(Config_CFD_file,"w")
            for line in open(options.filename):
                if "MATH_PROBLEM=" in line:
                    output_file.write("MATH_PROBLEM= LINEARIZED \n")
                elif "EXT_ITER=" in line:
                    output_file.write("EXT_ITER= 2 \n")
                elif "RESTART_SOL=" in line:
                    output_file.write("RESTART_SOL= YES \n")
                elif "SOLUTION_FLOW_FILENAME=" in line:
                    output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_finest)
                elif "SOLUTION_LIN_FILENAME=" in line:
                    output_file.write("SOLUTION_LIN_FILENAME= %s \n" % restart_lin_finest)
                elif "RESTART_LIN_FILENAME=" in line:
                    output_file.write("RESTART_LIN_FILENAME= %s \n" % restart_lin_finest)
                elif "MESH_FILENAME=" in line:
                    output_file.write("MESH_FILENAME= %s \n" %  mesh_finest)
                elif "FULLMG=" in line:
                    output_file.write("FULLMG= NO \n")
                elif "MGLEVEL=" in line:
                    output_file.write("MGLEVEL= 0 \n")
                elif "MGCYCLE=" in line:
                    output_file.write("MGCYCLE= 0 \n")
                else: output_file.write(line)
            output_file.close()
            if parallel: os.system ("parallel_computation.py -p %s -f %s -o True" % (int(options.partitions), Config_CFD_file))
            else: os.system ("SU2_CFD " + Config_CFD_file)

    # Once we have all the solutions files to write
    # the adaptation estimator, we perform the adaptation
    if iAdaptCycle < int(cycles)-1 or int(cycles) == 1: 
        output_file = open(Config_MAC_file,"w")
        for line in open(options.filename):
            if "SOLUTION_FLOW_FILENAME=" in line:
                if adapt_kind == "GRAD_FLOW" or adapt_kind == "GRAD_ADJOINT": 
                    output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_file)
                if adapt_kind == "ROBUST" or adapt_kind == "COMPUTABLE" or adapt_kind == "COMPUTABLE_ROBUST" or adapt_kind == "REMAINING": 
                    output_file.write("SOLUTION_FLOW_FILENAME= %s \n" % restart_flow_finest)
            elif "SOLUTION_ADJ_FILENAME=" in line:
                if adapt_kind == "GRAD_FLOW" or adapt_kind == "GRAD_ADJOINT": 
                    output_file.write("SOLUTION_ADJ_FILENAME= %s \n" % restart_adj_file)
                if adapt_kind == "ROBUST" or adapt_kind == "COMPUTABLE" or adapt_kind == "COMPUTABLE_ROBUST" or adapt_kind == "REMAINING": 
                    output_file.write("SOLUTION_ADJ_FILENAME= %s \n" % restart_adj_finest)
            elif "SOLUTION_LIN_FILENAME=" in line:
                if adapt_kind == "GRAD_FLOW" or adapt_kind == "GRAD_LINEAR": 
                    output_file.write("SOLUTION_LIN_FILENAME= %s \n" % restart_lin_file)
                if adapt_kind == "ROBUST" or adapt_kind == "COMPUTABLE" or adapt_kind == "COMPUTABLE_ROBUST" or adapt_kind == "REMAINING": 
                    output_file.write("SOLUTION_LIN_FILENAME= %s \n" % restart_lin_finest)
            elif "KIND_ADAPT=" in line:
                output_file.write("KIND_ADAPT= " + adapt_kind + "\n")
            elif "MESH_FILENAME=" in line:
                if iAdaptCycle == 0: output_file.write(line)
                else: output_file.write("MESH_FILENAME= %s \n" % Mesh_MAC_file)
            elif "MESH_OUT_FILENAME=" in line:
                output_file.write("MESH_OUT_FILENAME= %s \n" % Mesh_MAC_file)
            else: output_file.write(line)

        output_file.close()
        
        os.system ("SU2_MAC " + Config_MAC_file)

if options.overwrite: os.rename(Mesh_MAC_file, original_mesh_file)

if adapt_kind == "ROBUST" or adapt_kind == "COMPUTABLE" or adapt_kind == "COMPUTABLE_ROBUST" or adapt_kind == "REMAINING":
    os.remove(mesh_finest)
    os.remove(restart_flow_finest)

if adapt_kind == "COMPUTABLE_ROBUST":
    os.remove(restart_lin_finest)

# Remove configuration and mesh files
if int(cycles) != 0:
    os.remove(Config_MAC_file)
    os.remove(Config_CFD_file)
