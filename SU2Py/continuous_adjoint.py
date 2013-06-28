#!/usr/bin/env python 

## \file continuous_adjoint.py
#  \brief Python script for doing the continuous adjoint computation using the SU2 suite.
#  \author Francisco Palacios, Tom Economon, Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
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

import os, time, sys, shutil, numpy, libSU2, libSU2_run
from optparse import OptionParser
from parallel_computation import parallel_computation
from filter_adjoint import process_surface_adjoint

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():
    
    # Command Line Options
    parser=OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-c", "--compute",    dest="compute",    default="True",
                      help="COMPUTE direct and adjoint problem", metavar="COMPUTE")
    parser.add_option("-o", "--output",     dest="output",     default="True",
                      help="OUTPUT the domain solution", metavar="OUTPUT")    
    parser.add_option("-s", "--step",       dest="step",       default=1E-4,
                      help="finite difference STEP", metavar="STEP")    
    parser.add_option("-d", "--divide_grid",dest="divide_grid",default="True",
                      help="DIVIDE_GRID the numerical grid", metavar="DIVIDE_GRID")
    
    (options, args)=parser.parse_args()
    options.partitions = int( options.partitions )
    options.step       = float( options.step )
    
    continuous_adjoint( options.filename    ,
                        options.partitions  ,
                        options.compute     ,
                        options.output      ,
                        options.step        ,
                        options.divide_grid  )
    
#: def main()


# -------------------------------------------------------------------
#  Continuous Adjoint Function 
# -------------------------------------------------------------------

def continuous_adjoint( filename                , 
                        partitions  = 0         , 
                        compute     = "True"    ,
                        output      = "True"    ,
                        step        = 1e-4      , 
                        divide_grid = "True"     ):
    
    # General and default parameters
    Config_INP_filename   = filename
    Config_CADJ_filename  = "config_CADJ_" + Config_INP_filename
    FinDiff_Step          = float(step)
    partitions            = int(partitions)
    objfunc_grad_filename = "objfunc_grad_adj_" + Config_INP_filename.replace('.cfg','.dat')
    gradient_filename     = "cont_adj_" + Config_INP_filename.replace('.cfg','')  
    
    # do computations?
    if   compute == "True":
        compute_flow = True
        compute_adj  = True
        filter_adj   = False
        compute_grad = True
    elif compute == "Direct":
        compute_flow = True
        compute_adj  = False        
        filter_adj   = False
        compute_grad = False
    elif compute == "Adjoint":
        compute_flow = False
        compute_adj  = True
        filter_adj   = False
        compute_grad = True
    elif compute == "Gradient" or compute == "False":
        compute_flow = False
        compute_adj  = False
        filter_adj   = False
        compute_grad = True    
    elif compute == "Filtered":
        compute_flow = True
        compute_adj  = True
        filter_adj   = True
        compute_grad = True            
    else:
        raise Exception('Unrecognized compute option')
    
    # Make working config files
    shutil.copy(Config_INP_filename,Config_CADJ_filename)
    
    # Get parameters
    params_get = libSU2.Get_ConfigParams( Config_CADJ_filename )
    output_format         = params_get['OUTPUT_FORMAT']
    mesh_filename         = params_get['MESH_FILENAME']
    restart_flow_filename = params_get['RESTART_FLOW_FILENAME']
    surface_adj_filename  = params_get['SURFACE_ADJ_FILENAME']
    objective_function    = params_get['ADJ_OBJFUNC']
    Definition_DV         = params_get['DEFINITION_DV']
    History_filename      = params_get['CONV_FILENAME']
    Conv_Criteria         = params_get['CONV_CRITERIA']
    n_DV                  = len(Definition_DV['KIND'])
    
    # special physics cases
    special_cases = libSU2.get_SpecialCases(params_get)
    
    # special cauchy flag - runs adjoint with same iterations as direct problem
    special_cauchy = Conv_Criteria == 'SPECIAL_CAUCHY'
        
    # history filenames
    History_filename_CFD = History_filename + '_cfd'
    History_filename_ADJ = History_filename + '_adj'
    if 'TIME_SPECTRAL' in special_cases:
	History_filename_CFD = History_filename + '_TS_forces'
	History_filename_ADJ = History_filename + '_TS_forces'

    # plot type extentions
    plotfile_ext = libSU2.get_ExtensionName(output_format)
    
    # objective prefix
    cadj_prefix = libSU2.get_AdjointPrefix(objective_function)
    
    # surface filenames
    surface_adj_filtered = surface_adj_filename+ '_filtered'
    
    # initizlize outputs
    ObjFun_Dict   = []
    Gradient_Dict = []
        
    # -------------------------------------------------------------------
    # DIRECT CFD SOLUTION
    # -------------------------------------------------------------------
    if compute_flow:
        
        # Change the parameters to do a flow simulation
        params_set = { 'MATH_PROBLEM'  : 'DIRECT'                  ,
                       'CONV_FILENAME' : History_filename_CFD       }
	
	# in the case of time-spectral, don't overwite the TS_forces files
        if 'TIME_SPECTRAL' in special_cases:
	    params_set['CONV_FILENAME'] =  History_filename + '_cfd'

	if special_cauchy:
            params_set['CONV_CRITERIA'] = 'CAUCHY'
        libSU2.Set_ConfigParams( Config_CADJ_filename, params_set )

        # RUN THE FLOW SOLUTION        
        parallel_computation( filename    = Config_CADJ_filename ,
                              partitions  = partitions           ,  # handles partitions = 0 as serial
                              divide_grid = divide_grid          ,
                              output      = output                )
        
        # Get objective function values
        ObjFun_Dict = libSU2.get_ObjFunVals( History_filename_CFD+plotfile_ext , special_cases )
        
        # Get number of iterations
        History_Data_CFD = libSU2.read_History( History_filename_CFD+plotfile_ext )
        iterations_direct = int( History_Data_CFD['ITERATION'][-1] )
            
    #: if compute_flow
    
    # -------------------------------------------------------------------
    # ADJOINT CFD SOLUTION
    # -------------------------------------------------------------------
    if compute_adj:
           
        # Change the parameters to do an adjoint simulation
        params_set = { 'MATH_PROBLEM'           : 'ADJOINT'                 ,
                       'SOLUTION_FLOW_FILENAME' : restart_flow_filename     ,
                       'CONV_FILENAME'          : History_filename_ADJ       }
	
	# in the case of time-spectral, don't overwrite the TS_forces file 
	if 'TIME_SPECTRAL' in special_cases:
            params_set['CONV_FILENAME'] = History_filename + '_adj'
        
        # special cauchy => set adjoint iterations to direct iterations
        if special_cauchy:
            params_set['EXT_ITER'] = iterations_direct
            params_set['CONV_CRITERIA'] = 'RESIDUAL'
            
        # update parameters    
        libSU2.Set_ConfigParams( Config_CADJ_filename, params_set )

        # RUN THE ADJOINT SOLUTION
        parallel_computation( filename    = Config_CADJ_filename ,
                              partitions  = partitions           ,
                              divide_grid = 'False'              ,
                              output       = output               )
            
    #: if compute_adj
    
    
    # -------------------------------------------------------------------
    #  ADJOINT SMOOTHING (EXPERIMENTAL)
    # -------------------------------------------------------------------
    if filter_adj:
        # note: hard coded filter type and airfoil marker name
        
        # smooths airfoil surface sensitivities with a laplacian filter
        process_surface_adjoint(Config_CADJ_filename,'LAPLACE','airfoil')  
        
        # update surface filename
        params_set = { 'SURFACE_ADJ_FILENAME' : surface_adj_filtered }
        libSU2.Set_ConfigParams(Config_CADJ_filename,params_set)
        
    #: if filter_adj
    
    
    # -------------------------------------------------------------------
    # GRADIENT PROJECTION
    # -------------------------------------------------------------------
    if compute_grad:
        
        # Set the dv values for the finite difference step
        DV_value_old = numpy.zeros(n_DV)
        DV_value_new = numpy.ones(n_DV)*FinDiff_Step
    
        # Change the parameters of the design variables
        params_set = { 'DV_KIND'               : Definition_DV['KIND']      ,
                       'DV_MARKER'             : Definition_DV['MARKER'][0] ,
                       'DV_PARAM'              : Definition_DV['PARAM']     ,
                       'DV_VALUE_OLD'          : DV_value_old               ,
                       'DV_VALUE_NEW'          : DV_value_new               ,       
                       'GRAD_OBJFUNC_FILENAME' : objfunc_grad_filename       }  
        libSU2.Set_ConfigParams( Config_CADJ_filename, params_set )
        
        # Compute gradient with continuous adjoint
        libSU2_run.SU2_GPC( Config_CADJ_filename , partitions )
    
        # Get raw gradients
        gradients = libSU2.get_GradientVals(objfunc_grad_filename)    
        os.remove(objfunc_grad_filename)
    
    #: if compute_grad
    
    # -------------------------------------------------------------------
    # RETURN DATA
    # -------------------------------------------------------------------
    if compute_grad:
        
        # Get header and write format information
        gradient_filename = gradient_filename+"_"+cadj_prefix+plotfile_ext
        header,_ = libSU2.get_GradFileFormat( 'CONTINUOUS_ADJOINT'      ,
                                               output_format            ,
                                               Definition_DV['KIND'][0]  )
        # Start gradient file
        Grad_File = open(gradient_filename,'w')
        Grad_File.write(header)
        
        # Write output gradients and dv information 
        for i_DV in range(n_DV):
                
            _,write_format = libSU2.get_GradFileFormat( 'CONTINUOUS_ADJOINT',
                                                        output_format,
                                                        Definition_DV['KIND'][i_DV] )        
            
            # Merge information for output
            Output_Vector = [i_DV] + [gradients[i_DV]] + [FinDiff_Step] + Definition_DV['PARAM'][i_DV]
            Output_Vector = tuple(Output_Vector)
        
            # write row
            Grad_File.write(write_format % Output_Vector)
            
        #: for each DV
        
        # Dictionary of Gradients, one objective function
        Gradient_Dict = { objective_function : gradients }
    
    #: if compute_grad
    
    # Remove configuration file
    os.remove(Config_CADJ_filename)
        
    return [ObjFun_Dict,Gradient_Dict]

#: def continuous_adjoint()



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
