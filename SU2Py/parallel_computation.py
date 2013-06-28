#!/usr/bin/env python 

## \file parallel_computation.py
#  \brief Python script for doing the parallel computation using SU2_CFD.
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

import os, time, sys, shutil, libSU2, libSU2_run
from optparse import OptionParser
from merge_solution import merge_solution
from merge_unsteady import merge_unsteady

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
    parser.add_option("-o", "--output", dest="output", default="True",
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

def parallel_computation( filename              ,
                          partitions  = 2       , 
                          divide_grid = "True"  ,
                          output      = "True"   ):

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
    time_spectral = 'TIME_SPECTRAL' in special_cases

    # Set parameters
    params_set = { 'NUMBER_PART' : partitions }
    libSU2.Set_ConfigParams( Config_DDC_filename, params_set )
    libSU2.Set_ConfigParams( Config_CFD_filename, params_set )

    # get adjoint prefix
    if math_problem == 'ADJOINT':
        adj_objfunc = params_get['ADJ_OBJFUNC']
        adj_prefix   = libSU2.get_AdjointPrefix(adj_objfunc)


    # --------------------------------------------------------------
    #  SETUP PARALLEL JOB  
    # --------------------------------------------------------------
    if partitions > 1:

        # Just in case domain decomposition is needed
        if divide_grid != "False":
            libSU2_run.SU2_DDC(Config_DDC_filename,partitions)


    # --------------------------------------------------------------
    #  RUN PARALLEL JOB  
    # --------------------------------------------------------------    

    # MPI Call, Run the Job
    libSU2_run.SU2_CFD(Config_CFD_filename,partitions)


    # --------------------------------------------------------------
    #  CLEANUP JOB  
    # --------------------------------------------------------------        

    if partitions > 1:

        # Merge solution files
        if output == "True":
            if unsteady_simulation:
                merge_unsteady( Config_CFD_filename, partitions, 0, nExtIter )
            elif time_spectral:	
                merge_unsteady( Config_CFD_filename, partitions, 0, params_get['TIME_INSTANCES'] )
            else:
                merge_solution( Config_CFD_filename, partitions, output=output)

    # Remove ddc config file
    os.remove(Config_DDC_filename)

    #: if partitions>1

    # remove cfd config file
    os.remove(Config_CFD_filename) 

    # new-line character upon exiting merge script    
    print '\n'

#: def parallel_computation()



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()