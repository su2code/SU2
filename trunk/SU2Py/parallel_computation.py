#!/usr/bin/env python 

## \file parallel_computation.py
#  \brief Python script for doing the parallel computation using SU2_CFD.
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
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-d", "--divide_grid",     dest="divide_grid",     default="True",
                      help="DIVIDE_GRID the numerical grid", metavar="DIVIDE_GRID")
    parser.add_option("-o", "--output",     dest="output",     default="True",
                      help="OUTPUT the domain solution (.plt)", metavar="OUTPUT")
    
    (options, args)=parser.parse_args()    
    options.partitions = int( options.partitions )  

    # Run Parallel Computation
    parallel_computation ( options.filename            ,
                           options.partitions          , 
                           options.divide_grid         ,
                           options.output               )
#: def main()
    
    

# -------------------------------------------------------------------
#  Parallel Computation Function
# -------------------------------------------------------------------

def parallel_computation( filename                       ,
                          partitions           = 2       , 
                          divide_grid          = "True"  ,
                          output               = "True"  ,
                          logfile              = None     ):

    # General and default parameters
    Config_INP_filename = filename
    Config_DDC_filename = "config_DDC_" + Config_INP_filename
    Config_CFD_filename = "config_CFD_" + Config_INP_filename
    
    # Make working config files
    shutil.copy(Config_INP_filename,Config_DDC_filename)
    shutil.copy(Config_INP_filename,Config_CFD_filename)
    
    # Get parameters
    params_get = libSU2.Get_ConfigParams( Config_CFD_filename )
    math_problem           = params_get['MATH_PROBLEM']
    restart_solution       = params_get['RESTART_SOL']
    solution_flow_filename = params_get['SOLUTION_FLOW_FILENAME']
    solution_adj_filename  = params_get['SOLUTION_ADJ_FILENAME']
    nExtIter               = params_get['EXT_ITER'] 
    
    special_cases = libSU2.get_SpecialCases(params_get)
    unsteady_simulation = 'WRT_UNSTEADY' in special_cases

    # Set parameters
    params_set = { 'NUMBER_PART' : partitions }
    libSU2.Set_ConfigParams( Config_DDC_filename, params_set )
    libSU2.Set_ConfigParams( Config_CFD_filename, params_set )
  
    # get adjoint prefix
    if math_problem == 'ADJOINT':
        adj_objfunc = params_get['ADJ_OBJFUNC']
        adj_prefix   = libSU2.get_AdjointPrefix(adj_objfunc)
        
    # compose function call format and inputs
    run_SU2_DDC = os.path.join( SU2_RUN , "SU2_DDC " + Config_DDC_filename )
    run_SU2_CFD = os.path.join( SU2_RUN , "SU2_CFD " + Config_CFD_filename )
    run_SU2_DDC = "mpirun -np  1 %s" % (run_SU2_DDC)
    run_SU2_CFD = "mpirun -np %i %s" % (partitions, run_SU2_CFD)
    if unsteady_simulation:
        run_merge_solution = "merge_unsteady.py -p %s -f %s -b %s -e %s"  % (partitions, Config_CFD_filename, 0, nExtIter) 
    else:
        run_merge_solution = "merge_solution.py -p %s -f %s -o %s"        % (partitions, Config_CFD_filename, output)
    
    # add log command if desired
    if not logfile is None:
        run_SU2_DDC         = run_SU2_DDC         + ' >> ' + logfile
        run_SU2_CFD         = run_SU2_CFD         + ' >> ' + logfile
        run_merge_solution  = run_merge_solution  + ' >> ' + logfile        
        
    # --------------------------------------------------------------
    #  SETUP PARALLEL JOB  
    # --------------------------------------------------------------
    if partitions > 1:
        
        # Just in case domain decomposition is needed
        if divide_grid != "False":
            os.system ( run_SU2_DDC ) # >> log
            
            
    # --------------------------------------------------------------
    #  RUN PARALLEL JOB  
    # --------------------------------------------------------------    
      
    # MPI Call, Run the Job
    os.system ( run_SU2_CFD ) # >> log      
      
      
    # --------------------------------------------------------------
    #  CLEANUP JOB  
    # --------------------------------------------------------------        
    
    if partitions > 1:
    
        # Merge solution files
        if output == "True":
            os.system ( run_merge_solution )

        # Remove ddc config file
        os.remove(Config_DDC_filename)
        
    #: if partitions>1
    
    # remove cfd config file
    os.remove(Config_CFD_filename)     
    
#: def parallel_computation()



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
