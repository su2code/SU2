#!/usr/bin/python

## \file parallel_computation.py
#  \brief Python script for doing the parallel computation using SU2_CFD.
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


import os, time, libSU2, shutil
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=2,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-d", "--divide", dest="divide", default="True",
                  help="True/False DIVIDE the numerical grid using SU2_DDC", metavar="DIVIDE")
parser.add_option("-o", "--output", dest="output", default="True",
                  help="OUTPUT the domain solution", metavar="OUTPUT")

(options, args)=parser.parse_args()

# General and default parameters
Config_INP_filename = options.filename
Config_DDC_filename = "config_DDC_" + Config_INP_filename
Config_CFD_filename = "config_CFD_" + Config_INP_filename
unsteady_simulation = "NO"
options.partitions  = int(options.partitions)

# Make working config files
shutil.copy(Config_INP_filename,Config_DDC_filename)
shutil.copy(Config_INP_filename,Config_CFD_filename)

# Get parameters
params_get = libSU2.Get_ConfigParams( Config_CFD_filename )
math_problem           = params_get['MATH_PROBLEM']
restart_solution       = params_get['RESTART_SOL']
solution_flow_filename = params_get['SOLUTION_FLOW_FILENAME']    

# Set parameters
params_set = { 'NUMBER_PART' : options.partitions }
libSU2.Set_ConfigParams( Config_DDC_filename, params_set )

# get adjoint prefix
if math_problem == 'ADJOINT':
    solution_adj_filename  = params_get['SOLUTION_ADJ_FILENAME']
    cadj_objfunc           = params_get['CADJ_OBJFUNC']
    adj_prefix = libSU2.get_AdjointPrefix(cadj_objfunc)

# Case: Parallel Job
if options.partitions > 1:
    
    # Just in case domain decomposition is needed
    if options.divide != "False":
        os.system ("SU2_DDC " + Config_DDC_filename)
    
    # Split restart solutions
    if restart_solution == "YES" or math_problem == "ADJOINT": 
        os.system ( "divide_solution_su2.py -p %s -f %s" 
                    % (options.partitions, Config_DDC_filename) )
    
    # -------------------- #
    # RUN THE PARALLEL JOB #
    # -------------------- #
    os.system ( "mpirun -np %s $SU2_HOME/SU2Py/SU2_CFD %s" 
                % (int(options.partitions), Config_CFD_filename) )

    # Merge solution files
    os.system ( "merge_solution.py -p %s -f %s -o %s" 
                    % (options.partitions, Config_CFD_filename, options.output) )

    # Merge restart files
    os.system ( "merge_restart_su2.py -p %s -f %s" 
                % (options.partitions, Config_CFD_filename) )

    # Remove split solution files
    if restart_solution == "YES" or math_problem == "ADJOINT":
        for domain in range(options.partitions):
            os.remove( "%s_%s.%s" 
                       % (solution_flow_filename.split(".")[0], 
                          domain+1, solution_flow_filename.split(".")[1]) )
    if restart_solution == "YES" and math_problem == "ADJOINT":
        for domain in range(options.partitions):
            os.remove( "%s_%s_%s.%s" 
                       % (solution_adj_filename.split(".")[0], domain+1, 
                          adj_prefix, solution_adj_filename.split(".")[1]) )

    # Remove ddc config file
    os.remove(Config_DDC_filename)

# Otherwise: Serial Job
else:
    os.system ( "mpirun -np 1 $SU2_HOME/SU2Py/SU2_CFD %s"
                % (Config_CFD_filename) )
    
# remove cfd config file
os.remove(Config_CFD_filename)
