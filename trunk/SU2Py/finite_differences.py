#!/usr/bin/env python 

## \file finite_differences.py
#  \brief Python script for doing the finite differences computation using the SU2 suite.
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


import os, time, sys, shutil, numpy, libSU2
from optparse import OptionParser
from parallel_computation import parallel_computation
from parallel_deformation import parallel_deformation

SU2_RUN = os.environ['SU2_RUN'] 
sys.path.append( SU2_RUN ) 

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
        
    parser = OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-s", "--step",       dest="step",       default=1E-4,
                      help="finite difference STEP", metavar="STEP")
    
    (options, args)=parser.parse_args()
    options.partitions = int( options.partitions )
    options.step       = float( options.step )    
    
    finite_differences( options.filename   ,
                        options.partitions ,
                        options.step        )
#: def main()



# -------------------------------------------------------------------
#  Finite Differences Function 
# -------------------------------------------------------------------

def finite_differences( filename             , 
                        partitions = 1       , 
                        step       = 1e-4    ,
                        output     = 'False' ,
                        logfile    = None     ):
   
    # General and default parameters
    Config_INP_filename   = filename
    Config_CFD_filename   = "config_CFD_"       + Config_INP_filename
    Config_MDC_filename   = "config_MDC_"       + Config_INP_filename
    Mesh_MDC_filename     = "mesh_MDC_"         + Config_INP_filename.replace(".cfg",".su2")
    History_filename      = "history_"          + Config_INP_filename.replace('.cfg','')
    objfunc_grad_filename = "objfunc_grad_adj_" + Config_INP_filename.replace('.cfg','.dat')
    gradient_filename     = "fin_diff_"         + Config_INP_filename.replace('.cfg','')
    use_logfile           = not logfile is None
    
    # Finite Difference Step
    val_FDStep    = [1E-4] # default value, can add additional steps
    val_FDStep[0] = float(step)
    
    # running parallel?
    if int(partitions) == 1:
        parallel = False
    else : parallel = True
    
    # Read design variable information
    params_get = libSU2.Get_ConfigParams( Config_INP_filename )
    output_format          = params_get['OUTPUT_FORMAT']
    restart_flow_filename  = params_get['RESTART_FLOW_FILENAME']
    solution_flow_filename = params_get['SOLUTION_FLOW_FILENAME']
    Definition_DV          = params_get['DEFINITION_DV']
    Definition_AoA         = params_get['AoA']
    Definition_Mach        = params_get['MACH_NUMBER']
    config_restart         = params_get['RESTART_SOL']
    math_problem           = params_get['MATH_PROBLEM']
    n_DV = len(Definition_DV['KIND'])
    
    # use initial solution?
    restart_frominitial = config_restart == 'NO'
   
    # special physics cases
    special_cases = libSU2.get_SpecialCases(params_get)
    
    # Objective function names for plotting
    ObjFun_names = ['LIFT','DRAG','SIDEFORCE',         
                    'MOMENT_X','MOMENT_Y','MOMENT_Z',
                    'FORCE_X','FORCE_Y','FORCE_Z',
                    'EFFICIENCY']
    
    # special cases
    if 'FREE_SURFACE' in special_cases:
        ObjFun_names.append('FREESURFACE')
    if 'ROTATING_FRAME' in special_cases:
        ObjFun_names.extend( ['FIGURE_OF_MERIT','THRUST','TORQUE'] )    
    if 'EQUIV_AREA' in special_cases:
        ObjFun_names.extend( ['EQUIVALENT_AREA','NEARFIELD_PRESSURE'] )
    if 'AEROACOUSTIC_EULER' in special_cases:
        ObjFun_names.append('NOISE')

    # plot type extentions
    plotfile_ext = libSU2.get_ExtensionName(output_format)
    gradient_filename = gradient_filename+plotfile_ext
    
    # Get output header and write format
    #   takes the kind of the first design variable
    header,write_format =  libSU2.get_GradFileFormat( 'FINITE_DIFFERENCE'      ,
                                                      output_format            ,
                                                      Definition_DV['KIND'][0] ,
                                                      special_cases             )
    
    # Start Output File
    Gradient_File = open(gradient_filename,'w')
    Gradient_File.write(header)
    Gradient_File.close()
        
    # Start Output Dictionaries
    ObjFun_Dict   = dict.fromkeys(ObjFun_names,[])    
    Gradient_Dict = dict.fromkeys(ObjFun_names,[])    
    
    # copy working config files
    shutil.copy(Config_INP_filename,Config_MDC_filename)
    shutil.copy(Config_INP_filename,Config_CFD_filename)
    
    # some os commands
    run_SU2_CFD = "SU2_CFD " + Config_CFD_filename
    run_SU2_MDC = "SU2_MDC " + Config_MDC_filename
    
    if use_logfile:
        run_SU2_CFD = run_SU2_CFD + ' >> ' + logfile  
        run_SU2_MDC = run_SU2_MDC + ' >> ' + logfile  
    
    
    # -------------------------------------------------------------------
    # INITIAL CFD SOLUTION
    # -------------------------------------------------------------------
    
    #if use_logfile:
        #sys.stdout.write( "    Initial Solution  ... ")    
    
    # write parameters
    params_set = { 'MATH_PROBLEM'  : 'DIRECT'         ,
                   'CONV_FILENAME' : History_filename  }
    libSU2.Set_ConfigParams(Config_CFD_filename,params_set)
    
    # RUN CFD
    if parallel: 
        parallel_computation( filename            = Config_CFD_filename ,
                              partitions          = partitions          ,
                              divide_grid         = "True"              ,
                              output              = output              ,
                              logfile             = logfile              )
    else: 
        os.system ( run_SU2_CFD )
          
    # Read Objective Function Values
    ObjFun_values = libSU2.get_ObjFunVals( History_filename+plotfile_ext , special_cases )

    # Keep relevent values
    Ini_ObjFun_values = [ ObjFun_values[this_name] for this_name in ObjFun_names ]
    
    # make a numpy array of initial function values
    Ini_ObjFun_values = numpy.array( Ini_ObjFun_values )   
    
    # store initial objective function values in output dictionary
    for key,value in zip(ObjFun_names,Ini_ObjFun_values):
        ObjFun_Dict[key] = value    
    
    # set cfd config file to restart from initial solution
    if restart_frominitial and not('WRT_UNSTEADY' in special_cases):
        # move restart to solution name
        shutil.move(restart_flow_filename,solution_flow_filename)        
        # set config file
        params_set = { 'RESTART_SOL' : 'YES' }
        libSU2.Set_ConfigParams(Config_CFD_filename,params_set)

    #if use_logfile:
        #sys.stdout.write( "Done \n") 
        
    # -------------------------------------------------------------------
    # ITERATE THROUGH DESIGN VARIABLES
    # -------------------------------------------------------------------
    
    # loop over each design variable and finite difference step
    for iDesignVar in range(n_DV):
        for iFDStep in range(len(val_FDStep)):
            
            this_FDStep        = val_FDStep[iFDStep]
            this_DV_Kind       = Definition_DV['KIND'][iDesignVar]
            this_DV_Marker     = Definition_DV['MARKER'][iDesignVar]
            this_DV_Parameters = Definition_DV['PARAM'][iDesignVar]
            if this_DV_Kind == "AOA" or this_DV_Kind == "MACH_NUMBER":
                this_DV_Kind = "NO_DEFORMATION"

            #if use_logfile:
                #sys.stdout.write( "    Variable %3i      ... " % iDesignVar )   
            
            # -------------------------------------------------------------------
            # MESH DEFORMATION
            # -------------------------------------------------------------------        
            
            # Change the parameters of the design variables
            params_set = { 'DV_KIND'               : [this_DV_Kind]           ,
                           'DV_MARKER'             : this_DV_Marker           ,
                           'DV_PARAM'              : [this_DV_Parameters]     ,
                           'DV_VALUE_OLD'          : [0.0]                    ,
                           'DV_VALUE_NEW'          : [this_FDStep]            ,       
                           'GRAD_OBJFUNC_FILENAME' : objfunc_grad_filename    ,
                           'MESH_OUT_FILENAME'     : Mesh_MDC_filename         }    
            libSU2.Set_ConfigParams(Config_MDC_filename,params_set)
    
            # Apply grid deformation
            if parallel: 
                parallel_deformation( filename            = Config_MDC_filename ,
                                      partitions          = partitions          ,
                                      divide_grid         = "False"             ,
                                      merge_grid          = "False"             ,
                                      logfile             = logfile              )
            else: 
                os.system ( run_SU2_MDC )
    
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
                parallel_computation( filename            = Config_CFD_filename ,
                                      partitions          = partitions          ,
                                      divide_grid         = "False"             ,                                     
                                      output              = "False"             ,
                                      logfile             = logfile              )
            else: 
                os.system ( run_SU2_CFD )                
                
            # Read Objective Function Values
            ObjFun_values = libSU2.get_ObjFunVals( History_filename+plotfile_ext, special_cases )        
    
    
            # Keep relevent values
            This_ObjFun_values = [ ObjFun_values[this_name] for this_name in ObjFun_names ]            

            # Convert to Array
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
            Gradient_File.write( write_format % This_Output_Vector )
            Gradient_File.close()        
             
            # Save to gradient dictionary
            if iFDStep == 0:
                for key,value in zip(ObjFun_names,This_Gradient_values):
                    Gradient_Dict[key] = Gradient_Dict[key] + [value]
                
            # Cleanup some files
            if parallel:
                for domain in range(partitions):
                    Mesh_filename = os.path.splitext(Mesh_MDC_filename)[0] + "_" + str(domain+1) + ".su2"
                    os.remove(Mesh_filename)
            else:
                os.remove(Mesh_MDC_filename)
                
            os.remove(History_filename+plotfile_ext)
            
            #if use_logfile:
                #sys.stdout.write( "Done \n")             
    
        #: for each FDStep   
        
    #: for each DesignVariable
    
    # Cleanup some more files
    os.remove(Config_CFD_filename)
    os.remove(Config_MDC_filename)
    
    # returns a dictionary of gradients, and dictionary of objective function values
    return [ObjFun_Dict,Gradient_Dict]
        
    
#: def finite_differences()



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()


        
        
