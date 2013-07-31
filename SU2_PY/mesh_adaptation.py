#!/usr/bin/env python 

## \file mesh_adaptation.py
#  \brief Python script for doing the grid adaptation using the SU2 suite.
#  \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.6
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

import os, time, sys, libSU2, libSU2_run, shutil
from optparse import OptionParser
from parallel_computation import parallel_computation

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main(): 

    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=0,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-c", "--cycle", dest="cycle", default=1,
                      help="number of CYCLE adaptations", metavar="CYCLE")
    parser.add_option("-o", "--overwrite", dest="overwrite", default="False",
                      help="OVERWRITE_MESH the output mesh with the adapted one", metavar="OVERWRITE_MESH")
    parser.add_option("-s", "--save_all", dest="save_all", default="False",
                      help="SAVE_ALL the flow/adjoint/meshes solutions at each adaptation cycle", metavar="SAVE_ALL")

    (options, args)=parser.parse_args()

    options.partitions = int( options.partitions )
    options.cycle      = int( options.cycle      )
    options.overwrite  = options.overwrite == "True"    
    options.save_all   = options.save_all  == "True"

    # Run Mesh Adaptation
    mesh_adaptation ( options.filename   ,
                      options.partitions ,
                      options.cycle      ,
                      options.overwrite  ,
                      options.save_all    )

#: def main()


# -------------------------------------------------------------------
#  Mesh Adaptation Function
# -------------------------------------------------------------------

def mesh_adaptation( filename             ,
                     partitions   = 0     , 
                     cycles       = 1     ,
                     overwrite    = False ,
                     save_all     = False  ):

    # General and default parameters
    Config_INP_filename  = filename
    Config_CFD_filename  = "config_CFD_" + Config_INP_filename
    Config_MAC_filename  = "config_MAC_" + Config_INP_filename
#Mesh_MAC_filename    = "mesh_MAC_" + filename.replace(".cfg",".su2")
    finest_mesh_filename = "mesh_finest.su2"
    finest_flow_filename = "restart_flow_finest.dat"
    finest_lin_filename  = "restart_lin_finest.dat"
    finest_adj_filename  = "restart_adj_finest.dat"

    # assumes serial with partitions = 1
    if partitions == 1: partitions = 0

    # Get parameters
    params_get         = libSU2.Get_ConfigParams( Config_INP_filename )
    kind_adapt         = params_get['KIND_ADAPT']
    objfunc            = params_get['ADJ_OBJFUNC']
    restart_flow_file  = params_get['RESTART_FLOW_FILENAME']
    restart_adj_file   = params_get['RESTART_ADJ_FILENAME']
    original_mesh_file = params_get['MESH_FILENAME']
    #output_mesh_file   = params_get['MESH_OUT_FILENAME']
    Mesh_MAC_filename  = params_get['MESH_OUT_FILENAME']
    cadj_prefix        = libSU2.get_AdjointPrefix(objfunc)

    # Get solution file names
    volume_flow_file  = params_get['VOLUME_FLOW_FILENAME']
    volume_adj_file   = params_get['VOLUME_ADJ_FILENAME']
    surface_flow_file = params_get['SURFACE_FLOW_FILENAME']
    surface_adj_file  = params_get['SURFACE_ADJ_FILENAME']
    history_file      = params_get['CONV_FILENAME']

    # Get mesh filenames and filetypes
    mesh_filetype     = params_get['MESH_FORMAT']
    if mesh_filetype == "CGNS":
        error_str = "Currently cannot support mesh adaptation with CGNS grid files.  Please convert your CGNS mesh to SU2 format using the CGNS_TO_SU2 flag in the configuration file, re-specify the mesh file to the native .su2 file and set the MESH_FORMAT flag to SU2."
        print "\n*****\n" + error_str + "\n*****\n"
        return 1
    elif mesh_filetype == "NETCDF_ASCII":
        error_str ="Currently cannot support mesh adaptation with NETCDF_ASCII grid files.  Please convert your mesh to SU2 format, re-specify the mesh file to the native .su2 file and set the MESH_FORMAT flag to SU2."
        print "\n*****\n" + error_str + "\n*****\n"
        return 1

    # Get output solution filetypes
    output_filetype  = params_get['OUTPUT_FORMAT']
    if output_filetype == "TECPLOT":
        vol_file_ext = ".plt"
    elif output_filetype == "PARAVIEW":
        vol_file_ext = ".vtk"

    if( (kind_adapt == "ROBUST") or (kind_adapt == "COMPUTABLE_ROBUST") ):
        restart_lin_file   = params_get['RESTART_LIN_FILENAME']

    # Loop over number of adaptation cycles
    for iAdaptCycle in range(cycles):

        # Copy original input file to working files
        shutil.copy( Config_INP_filename, Config_MAC_filename )
        shutil.copy( Config_INP_filename, Config_CFD_filename )

        # Run direct flow simulation 
        # For iAdaptCycle == 0, store restart file, objective function and original mesh file
        params_set = { 'MATH_PROBLEM' : 'DIRECT'}

        if iAdaptCycle > 0:
            params_set.update({'RESTART_SOL'            : 'YES' , 
                               'ADJ_OBJFUNC'            : objfunc,
                               'RESTART_FLOW_FILENAME'  : restart_flow_file,
                               'RESTART_ADJ_FILENAME'   : restart_adj_file,
                               'SOLUTION_FLOW_FILENAME' : restart_flow_file,
                               'MESH_FILENAME'          : Mesh_MAC_filename  })
            if( (kind_adapt == "ROBUST") or kind_adapt == ("COMPUTABLE_ROBUST") ):
                params_set.update( {'RESTART_LIN_FILENAME'   : restart_lin_file} )

        # Load the new config file options and run the direct problem
        libSU2.Set_ConfigParams( Config_CFD_filename, params_set )
        if partitions > 1:
            parallel_computation( Config_CFD_filename, partitions )
        else:
            libSU2_run.SU2_CFD( Config_CFD_filename, partitions )

        # Copy flow solution & history file
        if save_all:
            print "Saving cycle " + str(iAdaptCycle) + " flow solution and history files..."
            print "Saving " + volume_flow_file + "_cycle" + str(iAdaptCycle) + vol_file_ext
            print "Saving " + surface_flow_file + "_cycle" + str(iAdaptCycle) + vol_file_ext
            print "Saving " + history_file + "_flow_cycle" + str(iAdaptCycle) + vol_file_ext
            shutil.move( volume_flow_file + vol_file_ext  , volume_flow_file + "_cycle"+str(iAdaptCycle)+vol_file_ext  )
            shutil.move( surface_flow_file + vol_file_ext , surface_flow_file + "_cycle"+str(iAdaptCycle)+vol_file_ext )
            shutil.move( surface_flow_file + ".csv"       , surface_flow_file + "_cycle"+str(iAdaptCycle)+".csv"       )
            shutil.move( history_file + vol_file_ext      , history_file + "_flow_cycle"+str(iAdaptCycle)+vol_file_ext )

        # If needed, run the adjoint simulation
        # For the first adaption cycle, use the filenames of the orignal .cfg file
        if ( kind_adapt == "GRAD_ADJOINT" or kind_adapt == "GRAD_FLOW_ADJ" or kind_adapt == "ROBUST" or kind_adapt == "COMPUTABLE_ROBUST" or kind_adapt == "COMPUTABLE" or kind_adapt == "REMAINING" ):
            params_set = { 'MATH_PROBLEM'           : 'ADJOINT',
                           'SOLUTION_FLOW_FILENAME' : restart_flow_file }

            if iAdaptCycle > 0:
                params_set.update({ 'RESTART_SOL'           : 'YES'            ,
                                    'SOLUTION_ADJ_FILENAME' : restart_adj_file ,
                                    'MESH_FILENAME'         : Mesh_MAC_filename })

            # Load the new config file options and run the adjoint problem
            libSU2.Set_ConfigParams( Config_CFD_filename, params_set )
            if partitions > 1:
                parallel_computation( Config_CFD_filename, partitions )
            else:
                libSU2_run.SU2_CFD( Config_CFD_filename, partitions )

            # Copy adjoint solution & history file
            if save_all:
                print "Saving cycle " + str(iAdaptCycle) + " adjoint solution and history files..."
                print "Saving " + volume_adj_file + "_cycle" + str(iAdaptCycle) + vol_file_ext
                print "Saving " + surface_adj_file + "_adj_cycle" + str(iAdaptCycle) + vol_file_ext
                print "Saving " + history_file + "_flow_cycle" + str(iAdaptCycle) + vol_file_ext
                shutil.move( volume_adj_file + vol_file_ext  , volume_adj_file + "_cycle"+str(iAdaptCycle)+vol_file_ext  )
                shutil.move( surface_adj_file + vol_file_ext , surface_adj_file + "_cycle"+str(iAdaptCycle)+vol_file_ext )
                shutil.move( surface_adj_file + ".csv"       , surface_adj_file + "_cycle"+str(iAdaptCycle)+".csv"       )
                shutil.move( history_file + vol_file_ext     , history_file + "_adj_cycle"+str(iAdaptCycle)+vol_file_ext )

        # If needed, change the parameters to run the first linear simulation
        # For the first adaptation cycle, use the filenames from the original .cfg file
        if kind_adapt == "COMPUTABLE_ROBUST":
            params_set = {'MATH_PROBLEM'           : 'LINEARIZED'     ,
                          'SOLUTION_FLOW_FILENAME' : restart_flow_file }

            if iAdaptCycle > 0:
                params_set.update({'RESTART_SOL'          : 'YES'            ,
                                   'RESTART_LIN_FILENAME' : restart_lin_file ,
                                   'MESH_FILENAME'        : Mesh_MAC_filename })			

            # Load the new config file options and run the linearized problem
            libSU2.Set_ConfigParams(Config_CFD_filename, params_set)
            if partitions > 1:
                parallel_computation( Config_CFD_filename, partitions )
            else:
                libSU2_run.SU2_CFD( Config_CFD_filename, partitions )

        # Change the parameters to do a direct and adjoint iteration over a fine grid
        if (    (kind_adapt == "ROBUST" or kind_adapt == "COMPUTABLE" or 
                 kind_adapt == "COMPUTABLE_ROBUST" or kind_adapt == "REMAINING") 
                and (iAdaptCycle < cycles-1 or cycles == 1) ):

            # Create the fine grid and interpolate the flow solution from coarse to refined grid
            params_set = { 'KIND_ADAPT'             : "FULL_FLOW"          ,
                           'SOLUTION_FLOW_FILENAME' : restart_flow_file    ,
                           'RESTART_FLOW_FILENAME'  : finest_flow_filename ,
                           'MESH_FILENAME'          : original_mesh_file   ,
                           'MESH_OUT_FILENAME'      : finest_mesh_filename  }

            if iAdaptCycle > 0:
                params_set.update( {'MESH_FILENAME' : Mesh_MAC_filename} )

            # Run the mesh adaptation module
            libSU2.Set_ConfigParams( Config_MAC_filename, params_set )
            libSU2_run.SU2_MAC(Config_MAC_filename,partitions)

            # Create the fine grid and interpolate the adjoint solution from coarse to refined grid
            params_set = { 'KIND_ADAPT'             : "FULL_ADJOINT"       ,
                           'SOLUTION_FLOW_FILENAME' : restart_flow_file    ,
                           'SOLUTION_ADJ_FILENAME'  : restart_adj_file     ,
                           'RESTART_FLOW_FILENAME'  : finest_flow_filename ,
                           'RESTART_ADJ_FILENAME'   : finest_adj_filename  ,
                           'MESH_FILENAME'          : original_mesh_file   ,
                           'MESH_OUT_FILENAME'      : finest_mesh_filename  }

            if iAdaptCycle > 0:
                params_set.update( {'MESH_FILENAME' : Mesh_MAC_filename} )

            # Run the mesh adaptation module
            libSU2.Set_ConfigParams( Config_MAC_filename, params_set )
            libSU2_run.SU2_MAC(Config_MAC_filename,partitions)

            # Create the fine grid and interpolate the linear solution from coarse to refined grid
            if kind_adapt == "COMPUTABLE_ROBUST":
                params_set = { 'KIND_ADAPT'             : "FULL_LINEAR"        ,
                               'SOLUTION_FLOW_FILENAME' : restart_flow_file    ,
                               'SOLUTION_LIN_FILENAME'  : restart_lin_file     ,
                               'RESTART_FLOW_FILENAME'  : finest_flow_filename ,
                               'RESTART_ADJ_FILENAME'   : finest_lin_filename  ,
                               'MESH_FILENAME'          : original_mesh_file   ,
                               'MESH_OUT_FILENAME'      : finest_mesh_filename  }

                if iAdaptCycle > 0:
                    params_set.update( {'MESH_FILENAME' : Mesh_MAC_filename} )

                # Run the mesh adaptation module
                libSU2.Set_ConfigParams( Config_MAC_filename, params_set )
                libSU2_run.SU2_MAC( Config_MAC_filename, partitions )

            # Change the parameters to do one iteration of the flow solver on the finest grid
            # Always start from the interpolated solution and store the residual in the solution file for finest grid
            # No multigrid or convergence acceleration
            params_set = { 'MATH_PROBLEM'           : 'DIRECT'             ,
                           'EXT_ITER'               : 2                    ,
                           'RESTART_SOL'            : 'YES'                ,
                           'SOLUTION_FLOW_FILENAME' : finest_flow_filename ,
                           'STORE_RESIDUAL'         : 'YES'                ,
                           'RESTART_FLOW_FILENAME'  : finest_flow_filename ,
                           'MESH_FILENAME'          : finest_mesh_filename ,
                           'FULLMG'                 : 'NO'                 ,
                           'MGLEVEL'                : 0                    ,
                           'MGCYCLE'                : 0                    ,
                           'MG_PRE_SMOOTH'          : '( 0 )'              ,
                           'MG_POST_SMOOTH'         : '( 0 )'              ,
                           'MG_CORRECTION_SMOOTH'   : '( 0 )'               }
            libSU2.Set_ConfigParams( Config_CFD_filename, params_set )
            if partitions > 1:
                parallel_computation( Config_CFD_filename, partitions )
            else:
                libSU2_run.SU2_CFD( Config_CFD_filename, partitions )

            # Change the parameters to do one iteration of the adjoint solver on the finest grid
            # Always start from the interpolated solution and store the residual in the solution file for finest grid
            # No multigrid or convergence acceleration
            params_set = { 'MATH_PROBLEM'           : 'ADJOINT'            ,
                           'EXT_ITER'               : 2                    ,
                           'RESTART_SOL'            : 'YES'                ,
                           'SOLUTION_FLOW_FILENAME' : finest_flow_filename ,
                           'SOLUTION_ADJ_FILENAME'  : finest_adj_filename  , 
                           'MESH_FILENAME'          : finest_mesh_filename ,
                           'FULLMG'                 : 'NO'                 ,
                           'MGLEVEL'                : 0                    ,
                           'MGCYCLE'                : 0                    ,
                           'MG_PRE_SMOOTH'          : '( 0 )'              ,
                           'MG_POST_SMOOTH'         : '( 0 )'              ,
                           'MG_CORRECTION_SMOOTH'   : '( 0 )'               }
            libSU2.Set_ConfigParams( Config_CFD_filename, params_set )
            if partitions > 1:
                parallel_computation( Config_CFD_filename, partitions )
            else:
                libSU2_run.SU2_CFD( Config_CFD_filename, partitions )

            # Change the parameters to do one iteration of the linear solver on the finest grid
            # Always start from the interpolated solution and store the residual in the solution file for finest grid
            # No multigrid or convergence acceleration
            if kind_adapt == "COMPUTABLE_ROBUST":
                params_set = { 'MATH_PROBLEM'           : 'LINEARIZED'         ,
                               'EXT_ITER'               : 2                    ,
                               'RESTART_SOL'            : 'YES'                ,
                               'SOLUTION_FLOW_FILENAME' : finest_flow_filename ,
                               'SOLUTION_LIN_FILENAME'  : finest_lin_filename  , 
                               'RESTART_LIN_FILENAME'   : finest_lin_filename  ,
                               'MESH_FILENAME'          : finest_mesh_filename ,
                               'FULLMG'                 : 'NO'                 ,
                               'MGLEVEL'                : 0                    ,
                               'MGCYCLE'                : 0                    ,
                               'MG_PRE_SMOOTH'          : '( 0 )'              ,
                               'MG_POST_SMOOTH'         : '( 0 )'              ,
                               'MG_CORRECTION_SMOOTH'   : '( 0 )'               }
                libSU2.Set_ConfigParams( Config_CFD_filename, params_set )
                if partitions > 1:
                    parallel_computation( Config_CFD_filename, partitions )
                else:
                    libSU2_run.SU2_CFD( Config_CFD_filename, partitions )

        # Perform adaptation using above solution files			
        if( (kind_adapt == "GRAD_FLOW") or (kind_adapt == "GRAD_ADJOINT") or (kind_adapt == "GRAD_FLOW_ADJ")):
            params_set = { 'SOLUTION_FLOW_FILENAME' : restart_flow_file,
                           'SOLUTION_ADJ_FILENAME'  : restart_adj_file }
        elif( (kind_adapt == "COMPUTABLE") or (kind_adapt == "REMAINING") ):
            params_set = { 'SOLUTION_FLOW_FILENAME' : finest_flow_filename,
                           'SOLUTION_ADJ_FILENAME'  : finest_adj_filename ,
                           'RESTART_FLOW_FILENAME'  : restart_flow_file   ,
                           'RESTART_ADJ_FILENAME'   : restart_adj_file     }
        elif( (kind_adapt == "ROBUST") or (kind_adapt == "COMPUTABLE_ROBUST") ):
            params_set = { 'SOLUTION_FLOW_FILENAME' : finest_flow_filename,
                           'SOLUTION_ADJ_FILENAME'  : finest_adj_filename ,
                           'SOLUTION_LIN_FILENAME'  : finest_lin_filename , 
                           'RESTART_FLOW_FILENAME'  : restart_flow_file   ,
                           'RESTART_ADJ_FILENAME'   : restart_adj_file    ,
                           'RESTART_LIN_FILENAME'   : restart_lin_file     }

        params_set.update({'KIND_ADAPT' : kind_adapt, 'MESH_OUT_FILENAME' : Mesh_MAC_filename})

        if iAdaptCycle > 0:
            params_set.update({'MESH_FILENAME' : Mesh_MAC_filename})

        # Run the mesh adaptation module
        libSU2.Set_ConfigParams( Config_MAC_filename, params_set )
        libSU2_run.SU2_MAC( Config_MAC_filename, partitions )

        # Copy cycle mesh file
        if save_all:
            print "Saving cycle " + str(iAdaptCycle) + " mesh file..."
            shutil.copy( Mesh_MAC_filename, Mesh_MAC_filename.replace(".su2", "_cycle"+str(iAdaptCycle)+".su2"))

    # Clean up
    if overwrite : os.rename( Mesh_MAC_filename, original_mesh_file )
    if( (kind_adapt == "ROBUST") or (kind_adapt == "COMPUTABLE") or
        (kind_adapt == "COMPUTABLE_ROBUST") or (kind_adapt == "REMAINING") ):
        os.remove(finest_mesh_filename)
        os.remove(finest_flow_filename)

    if kind_adapt == "COMPUTABLE_ROBUST":
        os.remove(finest_lin_filename)

    if save_all:
        os.remove(Mesh_MAC_filename)

    os.remove(Config_MAC_filename)
    os.remove(Config_CFD_filename)


#: def mesh_adaptation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()