#!/usr/bin/python

## \file finite_differences.py
#  \brief Python script for doing the finite differences computation using the SU2 suite.
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


import os, time, sys, shutil, numpy, libSU2
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-s", "--step", dest="step", default=1E-4,
                  help="finite difference STEP", metavar="STEP")
parser.add_option("-o", "--output", dest="output", default="True",
                  help="OUTPUT the domain solution", metavar="OUTPUT")

(options, args)=parser.parse_args()


# General and default parameters
Config_INP_filename   = options.filename
Config_CFD_filename   = "config_CFD_"       + Config_INP_filename
Config_MDC_filename   = "config_MDC_"       + Config_INP_filename
Mesh_MDC_filename     = "mesh_MDC_"         + Config_INP_filename.replace(".cfg",".su2")
History_filename      = "history_"          + Config_INP_filename.replace('.cfg','')
objfunc_grad_filename = "objfunc_grad_adj_" + Config_INP_filename.replace('.cfg','.dat')
gradient_filename     = "fin_diff_"         + Config_INP_filename.replace('.cfg','')
options.partitions    = int(options.partitions)

# Finite Difference Step
val_FDStep    = [1E-4] # default value, can add additional steps
val_FDStep[0] = float(options.step)

# running parallel?
if int(options.partitions) == 1:
    parallel = False
else : parallel = True

# Read design variable information
params_get = libSU2.Get_ConfigParams( Config_INP_filename )
output_format         = params_get['OUTPUT_FORMAT']
restart_flow_filename = params_get['RESTART_FLOW_FILENAME']
Definition_DV         = params_get['DEFINITION_DV']
Definition_AoA        = params_get['AoA']
Definition_Mach       = params_get['MACH_NUMBER']
n_DV = len(Definition_DV['Kind'])

# special physics cases
special_cases = []
all_special_cases = ['FREE_SURFACE','ROTATING_FRAME','EQUIV_AREA']
for key in all_special_cases:
    if params_get.has_key(key):
        special_cases.append(key)

# no support for more than one special case
if len(special_cases) > 1:
    error_str = 'Currently cannot support ' + ' and '.join(special_cases) + 'at once'
    raise Exception(error_str)

# plot type extentions
plotfile_ext = libSU2.get_ExtensionName(output_format)
gradient_filename = gradient_filename+plotfile_ext

# copy working config files
shutil.copy(Config_INP_filename,Config_MDC_filename)
shutil.copy(Config_INP_filename,Config_CFD_filename)

# -------------------------------------------------------------------
# INITIAL CFD SOLUTION
# -------------------------------------------------------------------

# write parameters
params_set = { 'MATH_PROBLEM'  : 'DIRECT'         ,
               'CONV_FILENAME' : History_filename  }
libSU2.Set_ConfigParams(Config_CFD_filename,params_set)

# RUN CFD
if parallel: 
    os.system ( "parallel_computation.py -p %s -f %s -o %s. -d True" 
                % (options.partitions, Config_CFD_filename, options.output) )
else: 
    os.system ("SU2_CFD " + Config_CFD_filename)

# Read Objective Function Values
ObjFun_values = libSU2.get_ObjFunVals( History_filename+plotfile_ext, special_cases )

Ini_ObjFun_values = [ ObjFun_values['LIFT']       ,
                      ObjFun_values['DRAG']       ,
                      ObjFun_values['SIDEFORCE']  ,
                      ObjFun_values['MOMENT_X']   ,
                      ObjFun_values['MOMENT_Y']   ,
                      ObjFun_values['MOMENT_Z']   , 
                      ObjFun_values['FORCE_X']    , 
                      ObjFun_values['FORCE_Y']    ,
                      ObjFun_values['FORCE_Z']    ,
                      ObjFun_values['EFFICIENCY']  ]

# Special Physics Cases
for key in special_cases:
    if key == 'FREE_SURFACE':
        Ini_ObjFun_values.append(ObjFun_values['FREE_SURFACE']) # overloaded column
    if key == 'ROTATING_FRAME':
        Ini_ObjFun_values.append( [ ObjFun_values['EQUIVALENT_AREA']    ,
                                    ObjFun_values['NEARFIELD_PRESSURE']  ] )
    if key == 'EQUIV_AREA':
        Ini_ObjFun_values.append( [ ObjFun_values['FIGURE_OF_MERIT'] ,
                                    ObjFun_values['THRUST']          ,
                                    ObjFun_values['TORQUE']           ] )

# make a numpy array of initial function values
Ini_ObjFun_values = numpy.array( Ini_ObjFun_values )   

# Get output header and write format
#   takes the kind of the first design variable
header,write_format =  libSU2.get_GradFileFormat( 'FINITE_DIFFERENCE'      ,
                                                  output_format            ,
                                                  Definition_DV['Kind'][0] ,
                                                  special_cases             )

# Start Output File
Gradient_File = open(gradient_filename,'w')
Gradient_File.write(header)
Gradient_File.close()
    
    
    
# -------------------------------------------------------------------
# ITERATE THROUGH DESIGN VARIABLES
# -------------------------------------------------------------------

# loop over each design variable and finite difference step
for iDesignVar in range(n_DV):
    for this_FDStep in val_FDStep:
        
        this_DV_Kind       = Definition_DV['Kind'][iDesignVar]
        this_DV_Marker     = Definition_DV['Markers'][iDesignVar]
        this_DV_Parameters = Definition_DV['Parameters'][iDesignVar]
        if this_DV_Kind == "AOA" or this_DV_Kind == "MACH_NUMBER":
            this_DV_Kind = "NO_DEFORMATION"
            
        # Get write format for this design variable
        _, this_write_format =  libSU2.get_GradFileFormat( 'FINITE_DIFFERENCE' ,
                                                           output_format       ,
                                                           this_DV_Kind        ,
                                                           special_cases        )            
        
        # -------------------------------------------------------------------
        # MESH DEFORMATION
        # -------------------------------------------------------------------        
        
        # Change the parameters of the design variables
        params_set = { 'DV_KIND'               : [this_DV_Kind]        ,
                       'DV_MARKER'             : this_DV_Marker        ,
                       'DV_PARAM'              : [this_DV_Parameters]  ,
                       'DV_VALUE_OLD'          : [0]                   ,
                       'DV_VALUE_NEW'          : [this_FDStep]         ,       
                       'GRAD_OBJFUNC_FILENAME' : objfunc_grad_filename ,
                       'MESH_OUT_FILENAME'     : Mesh_MDC_filename      }  
        libSU2.Set_ConfigParams(Config_MDC_filename,params_set)

        # Apply grid deformation
        if parallel:
            os.system ("SU2_MDC " + Config_MDC_filename)
#            os.system ( "mpirun -np %s $SU2_HOME/SU2Py/SU2_MDC %s"
#                        % (int(options.partitions), Config_MDC_filename) )
        else: 
            os.system ("SU2_MDC " + Config_MDC_filename)


        # -------------------------------------------------------------------
        # CFD SOLUTION
        # -------------------------------------------------------------------     

        # Change the parameters of the design variables
        params_set = { 'MATH_PROBLEM'  : 'DIRECT'          ,
                       'CONV_FILENAME' : History_filename  ,
                       'MESH_FILENAME' : Mesh_MDC_filename  }  
        
        if this_DV_Kind == 'MACH_NUMBER':
            params_set['MACH_NUMBER'] = Definition_Mach + this_FDStep
        if this_DV_Kind == 'AOA':
            params_set['AoA']         = Definition_AoA  + this_FDStep            
        
        libSU2.Set_ConfigParams(Config_CFD_filename,params_set)
       
        # CFD computation to evaluate nondimensional coefficient
        if parallel: 
            os.system ( "parallel_computation.py -p %s -f %s -o %s -d True" 
                        % (options.partitions, Config_CFD_filename, options.output) )
#            os.system ( "parallel_computation.py -p %s -f %s -o %s -d False" 
#                        % (options.partitions, Config_CFD_filename, options.output) )
        else: 
            os.system ("SU2_CFD " + Config_CFD_filename)
            
        # Read Objective Function Values
        ObjFun_values = libSU2.get_ObjFunVals( History_filename+plotfile_ext, special_cases )        

        This_ObjFun_values = [ ObjFun_values['LIFT']       ,
                               ObjFun_values['DRAG']       ,
                               ObjFun_values['SIDEFORCE']  ,
                               ObjFun_values['MOMENT_X']   ,
                               ObjFun_values['MOMENT_Y']   ,
                               ObjFun_values['MOMENT_Z']   , 
                               ObjFun_values['FORCE_X']    , 
                               ObjFun_values['FORCE_Y']    ,
                               ObjFun_values['FORCE_Z']    ,
                               ObjFun_values['EFFICIENCY']  ]
        
        # Special Physics Cases
        for key in special_cases:
            if key == 'FREE_SURFACE':
                This_ObjFun_values.append( ObjFun_values['FREE_SURFACE'] ) # overloaded column
            if key == 'ROTATING_FRAME':
                This_ObjFun_values.append( [ ObjFun_values['EQUIVALENT_AREA']    ,
                                             ObjFun_values['NEARFIELD_PRESSURE']  ] )
            if key == 'EQUIV_AREA':
                This_ObjFun_values.append( [ ObjFun_values['FIGURE_OF_MERIT'] ,
                                             ObjFun_values['THRUST']          ,
                                             ObjFun_values['TORQUE']           ] )        

        This_ObjFun_values = numpy.array( This_ObjFun_values )
        
        # Finite Difference Gradient
        This_Gradient_values = (This_ObjFun_values - Ini_ObjFun_values) / this_FDStep
        
        # return to list for merging with other outputs
        This_Gradient_values = This_Gradient_values.tolist()
        
        # merge all outputs for writing
        This_Output_Vector = [iDesignVar] + This_Gradient_values + this_DV_Parameters + [this_FDStep]   
        This_Output_Vector = tuple(This_Output_Vector)
        
        # Write to gradient file
        Gradient_File = open(gradient_filename,'a')
        Gradient_File.write( this_write_format % This_Output_Vector )
        Gradient_File.close()        
        
        # Cleanup some files         
        os.remove(Mesh_MDC_filename)
        os.remove(History_filename+plotfile_ext)

    #: for each FDStep   
    
#: for each DesignVariable

# Cleanup some more files
os.remove(Config_CFD_filename)
os.remove(Config_MDC_filename)


        
        
