#!/usr/bin/python

## \file parallel_deformation.py
#  \brief Python script for doing the parallel deformation using SU2_MDC.
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
Config_MDC_filename = "config_MDC_" + Config_INP_filename
options.partitions  = int(options.partitions)

# Make working config files
shutil.copy(Config_INP_filename,Config_DDC_filename)
shutil.copy(Config_INP_filename,Config_MDC_filename)

# Get parameters
params_get = libSU2.Get_ConfigParams( Config_MDC_filename )
solution_flow_filename = params_get['SOLUTION_FLOW_FILENAME']
visualize_deformation = params_get['VISUALIZE_DEFORMATION']

# Set parameters
params_set = { 'NUMBER_PART' : options.partitions }
libSU2.Set_ConfigParams( Config_DDC_filename, params_set )

# Case: Parallel Job
if options.partitions > 1:
    
    # Just in case domain decomposition is needed
    if options.divide != "False":
        os.system ("SU2_DDC " + Config_DDC_filename)
       
    # -------------------- #
    # RUN THE PARALLEL JOB #
    # -------------------- #
    os.system ( "mpirun -np %s $SU2_HOME/SU2Py/SU2_MDC %s" 
                % (int(options.partitions), Config_MDC_filename) )

    if visualize_deformation == "YES":
        # Merge deformed grid files
        params_set = { 'MATH_PROBLEM'  : 'DIRECT'          ,
                   'VOLUME_FLOW_FILENAME' : "deformed_volumetric_grid",
                   'SURFACE_FLOW_FILENAME' : "deformed_surface_grid" }
        libSU2.Set_ConfigParams( Config_MDC_filename, params_set )
        os.system ( "merge_solution.py -p %s -f %s -o %s" 
                    % (options.partitions, Config_MDC_filename, options.output) )

        # Merge original grid files
        params_set = { 'MATH_PROBLEM'  : 'DIRECT'          ,
                   'VOLUME_FLOW_FILENAME' : "original_volumetric_grid",
                   'SURFACE_FLOW_FILENAME' : "original_surface_grid" }
        libSU2.Set_ConfigParams( Config_MDC_filename, params_set )
        os.system ( "merge_solution.py -p %s -f %s -o %s" 
                    % (options.partitions, Config_MDC_filename, options.output) )

    os.system ( "merge_grid_su2.py -p %s -f %s" 
                    % (options.partitions, Config_MDC_filename) )

    # Remove ddc config file
    os.remove(Config_DDC_filename)

# Otherwise: Serial Job
else:
    os.system ( "mpirun -np 1 $SU2_HOME/SU2Py/SU2_MDC %s"
                % (Config_MDC_filename) )
    
# remove cfd config file
os.remove(Config_MDC_filename)
