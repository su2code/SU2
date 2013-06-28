#!/usr/bin/env python 

## \file parallel_adaptation.py
#  \brief Python script for doing the parallel adaptation using SU2_MAC.
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

import os, time, sys, libSU2, shutil
from optparse import OptionParser

SU2_RUN = os.environ['SU2_RUN'] 
sys.path.append( SU2_RUN ) 

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
                      help="MERGE_GRID the adapted grid", metavar="MERGE_GRID")

    (options, args)=parser.parse_args()
    options.partitions = int( options.partitions )  

    # Run Parallel Comutation
    parallel_adaptation ( options.filename        ,
                           options.partitions      , 
                           options.divide_grid     ,
                           options.merge_grid       )
#: def main()
    
# -------------------------------------------------------------------
#  Parallel Computation Function
# -------------------------------------------------------------------

def parallel_adaptation( filename    ,
                          partitions      = 2       , 
                          divide_grid     = "True"  ,
                          merge_grid      = "True"  ,
                          logfile         = None     ):
    
    # General and default parameters
    Config_INP_filename = filename
    Config_DDC_filename = "config_DDC_" + Config_INP_filename
    Config_MAC_filename = "config_MAC_" + Config_INP_filename

    # Make working config files
    shutil.copy(Config_INP_filename,Config_DDC_filename)
    shutil.copy(Config_INP_filename,Config_MAC_filename)

    # Get parameters
    params_get = libSU2.Get_ConfigParams( Config_MAC_filename )
    solution_flow_filename = params_get['SOLUTION_FLOW_FILENAME']
    visualize_deformation  = params_get.get('VISUALIZE_DEFORMATION','NO')

    # Set parameters
    params_set = { 'NUMBER_PART' : partitions }
    libSU2.Set_ConfigParams( Config_DDC_filename, params_set )

    # run commands
    run_SU2_DDC = os.path.join( SU2_RUN , "SU2_DDC " + Config_DDC_filename )
    run_SU2_MAC = os.path.join( SU2_RUN , "SU2_MAC " + Config_MAC_filename )
    run_SU2_DDC = "mpirun -np  1 %s" % (run_SU2_DDC)
    run_SU2_MAC = "mpirun -np %i %s" % (partitions, run_SU2_MAC)
    run_merge_solution = "merge_solution.py -p %s -f %s -o True"  % (partitions, Config_MAC_filename)
    run_merge_grid     = "merge_grid_su2.py -p %s -f %s" % (partitions, Config_MAC_filename)
    
    # log files
    if not logfile is None:
        run_SU2_DDC = run_SU2_DDC + ' >> ' + logfile
        run_SU2_MAC = run_SU2_MAC + ' >> ' + logfile
        run_merge_solution = run_merge_solution + ' >> ' + logfile
        run_merge_grid     = run_merge_grid     + ' >> ' + logfile

    # --------------------------------------------------------------
    #  SETUP PARALLEL JOB  
    # --------------------------------------------------------------
    if partitions > 1:
    
        # Just in case domain decomposition is needed
        if divide_grid != "False":
            os.system ( run_SU2_DDC )
        
        # Remove ddc config file
        os.remove(Config_DDC_filename)   
        
        # ---------------------- #
        #  RUN THE PARALLEL JOB  #
        # ---------------------- #
        os.system ( run_SU2_MAC )

        # Merge plt files
        if visualize_deformation == "YES":
            # Merge adapted grid files
            params_set = { 'MATH_PROBLEM'          : 'DIRECT'                   ,
                           'VOLUME_FLOW_FILENAME'  : "adapted_grid" ,
                           'SURFACE_FLOW_FILENAME' : "adapted_surface"     }
            libSU2.Set_ConfigParams( Config_MAC_filename, params_set )
            os.system ( run_merge_solution )

            # Merge original grid files
            params_set = { 'MATH_PROBLEM'          : 'DIRECT'                   ,
                           'VOLUME_FLOW_FILENAME'  : "original_grid" ,
                           'SURFACE_FLOW_FILENAME' : "original_surface"     }
            libSU2.Set_ConfigParams( Config_MAC_filename, params_set )
            os.system ( run_merge_solution )

        # Just in case the grid merging is needed
        if merge_grid == "True":
            os.system ( run_merge_grid )

    # Otherwise: Serial Job
    else:
        os.system ( run_SU2_MAC )
    
    # remove cfd config file
    os.remove(Config_MAC_filename)

#: def parallel_adaptation()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
