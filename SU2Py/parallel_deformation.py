#!/usr/bin/env python 

## \file parallel_deformation.py
#  \brief Python script for doing the parallel deformation using SU2_MDC.
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

import libSU2, libSU2_run
from merge_solution import merge_solution
from merge_grid_su2 import merge_grid_su2
import os, time, sys, shutil
from optparse import OptionParser

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-d", "--divide_grid", dest="divide_grid", default="True",
                      help="DIVIDE_GRID the numerical grid", metavar="DIVIDE_GRID")
    parser.add_option("-m", "--merge_grid",     dest="merge_grid",     default="True",
                      help="MERGE_GRID the deformed grid", metavar="MERGE_GRID")

    (options, args)=parser.parse_args()
    options.partitions = int( options.partitions )  

    # Run Parallel Comutation
    parallel_deformation ( options.filename    ,
                           options.partitions  , 
                           options.divide_grid ,
                           options.merge_grid   )
#: def main()
    
# -------------------------------------------------------------------
#  Parallel Computation Function
# -------------------------------------------------------------------

def parallel_deformation( filename             ,
                          partitions  = 2      , 
                          divide_grid = "True" ,
                          merge_grid  = "True"  ):
    
    # General and default parameters
    Config_INP_filename = filename
    Config_DDC_filename = "config_DDC_" + Config_INP_filename
    Config_MDC_filename = "config_MDC_" + Config_INP_filename

    # Make working config files
    shutil.copy(Config_INP_filename,Config_DDC_filename)
    shutil.copy(Config_INP_filename,Config_MDC_filename)

    # Get parameters
    params_get = libSU2.Get_ConfigParams( Config_MDC_filename )
    visualize_deformation  = params_get.get('VISUALIZE_DEFORMATION','NO')

    # Set parameters
    params_set = { 'NUMBER_PART' : partitions }
    libSU2.Set_ConfigParams( Config_DDC_filename, params_set )
    
    
    # --------------------------------------------------------------
    #  SETUP PARALLEL JOB  
    # --------------------------------------------------------------
    if partitions > 1:
    
        # Just in case domain decomposition is needed
        if divide_grid != "False":
            libSU2_run.SU2_DDC(Config_DDC_filename,partitions)
        
        # ---------------------- #
        #  RUN THE PARALLEL JOB  #
        # ---------------------- #
        libSU2_run.SU2_MDC(Config_MDC_filename,partitions)

        # Merge plt files
        if visualize_deformation == "YES":
            # Merge deformed grid files
            params_set = { 'MATH_PROBLEM'          : 'DIRECT'                   ,
                           'VOLUME_FLOW_FILENAME'  : "deformed_volumetric_grid" ,
                           'SURFACE_FLOW_FILENAME' : "deformed_surface_grid"     }
            libSU2.Set_ConfigParams( Config_MDC_filename, params_set )

            merge_solution( Config_MDC_filename     ,
                            partitions = partitions ,
                            output     = "True"      )  

            # Merge original grid files
            params_set = { 'MATH_PROBLEM'          : 'DIRECT'                   ,
                           'VOLUME_FLOW_FILENAME'  : "original_volumetric_grid" ,
                           'SURFACE_FLOW_FILENAME' : "original_surface_grid"     }
            libSU2.Set_ConfigParams( Config_MDC_filename, params_set )

            merge_solution( Config_MDC_filename     ,
                            partitions = partitions ,
                            output     = "True"      )              

        # Just in case the grid merging is needed
        if merge_grid == "True":
            merge_grid_su2( Config_MDC_filename , partitions )    

    # Otherwise: Serial Job
    else:
        libSU2_run.SU2_MDC(Config_MDC_filename,partitions)
    
    # remove config files
    os.remove(Config_MDC_filename)
    os.remove(Config_DDC_filename)     

#: def parallel_deformation()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
