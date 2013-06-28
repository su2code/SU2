#!/usr/bin/python

## \file continuous_adjoint.py
#  \brief Python script for doing the continuous adjoint computation using the SU2 suite.
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

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-p", "--partitions", dest="partitions", default=1,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-s", "--step", dest="step", default=1E-4,
                  help="finite difference STEP", metavar="STEP")
parser.add_option("-c", "--compute", dest="compute", default="True",
                  help="COMPUTE direct and adjoint problem", metavar="COMPUTE")
parser.add_option("-d", "--divide", dest="divide", default="True",
                  help="True/False DIVIDE the numerical grid using SU2_DDC", metavar="DIVIDE")

(options, args)=parser.parse_args()

# General and default parameters
Config_INP_filename   = options.filename
Config_CFD_filename   = "config_CFD_" + Config_INP_filename
Config_GPC_filename   = "config_GPC_" + Config_INP_filename
FinDiff_Step          = float(options.step)
options.partitions    = int(options.partitions)
History_filename      = "history_" + Config_INP_filename.replace('.cfg','')
objfunc_grad_filename = "objfunc_grad_adj_" + Config_INP_filename.replace('.cfg','.dat')
gradient_filename     = "cont_adj_" + Config_INP_filename.replace('.cfg','')

# running parallel?
if options.partitions == 1:
    parallel = False
else : parallel = True

# do computations?
if options.compute == "False":
    compute_adj  = False
    compute_flow = False
else:
    compute_adj = True
    compute_flow = True    

# Make working config files
shutil.copy(Config_INP_filename,Config_CFD_filename)
shutil.copy(Config_INP_filename,Config_GPC_filename)

# Get parameters
params_get = libSU2.Get_ConfigParams( Config_CFD_filename )
output_format         = params_get['OUTPUT_FORMAT']
restart_flow_filename = params_get['RESTART_FLOW_FILENAME']
objective_function    = params_get['CADJ_OBJFUNC']
Definition_DV         = params_get['DEFINITION_DV']
n_DV = len(Definition_DV['Kind'])

# plot type extentions
plotfile_ext = libSU2.get_ExtensionName(output_format)
    
    
# -------------------------------------------------------------------
# DIRECT CFD SOLUTION
# -------------------------------------------------------------------

# Change the parameters to do a flow simulation
params_set = {'MATH_PROBLEM' : 'DIRECT'}
libSU2.Set_ConfigParams( Config_CFD_filename, params_set )

if compute_flow:
    if parallel: 
        os.system ( "parallel_computation.py -p %s -f %s -o False -d %s" 
                    % (int(options.partitions), Config_CFD_filename, options.divide) )
    else: 
        os.system ("SU2_CFD " + Config_CFD_filename)


# -------------------------------------------------------------------
# ADJOINT CFD SOLUTION
# -------------------------------------------------------------------

# Change the parameters to do an adjoint simulation
params_set = { 'MATH_PROBLEM'           : 'ADJOINT'             ,
               'SOLUTION_FLOW_FILENAME' : restart_flow_filename ,
               'CONV_FILENAME'          : History_filename       }
libSU2.Set_ConfigParams( Config_CFD_filename, params_set )

if compute_adj:
    if parallel: 
        os.system ( "parallel_computation.py -p %s -f %s -d False -o False" 
                    % (int(options.partitions), Config_CFD_filename) )
    else: 
        os.system ("SU2_CFD " + Config_CFD_filename)


# -------------------------------------------------------------------
# GRADIENT PROJECTION
# -------------------------------------------------------------------

# Change the parameters of the design variables
params_set = { 'DV_KIND'               : Definition_DV['Kind']         ,
               'DV_MARKER'             : Definition_DV['Markers'][0]   ,
               'DV_PARAM'              : Definition_DV['Parameters']   ,
               'DV_VALUE_OLD'          : numpy.zeros(n_DV)             ,
               'DV_VALUE_NEW'          : numpy.ones(n_DV)*FinDiff_Step ,       
               'GRAD_OBJFUNC_FILENAME' : objfunc_grad_filename          }
libSU2.Set_ConfigParams( Config_GPC_filename, params_set )

# Compute gradient with continuous adjoint
if parallel:
    os.system ("SU2_GPC " + Config_GPC_filename)
#    os.system ( "mpirun -np %s $SU2_HOME/SU2Py/SU2_GPC %s"
#                % (int(options.partitions), Config_GPC_filename) )
else: 
    os.system ("SU2_GPC " + Config_GPC_filename)


# -------------------------------------------------------------------
# RETURN DATA
# -------------------------------------------------------------------

# Get raw gradients
gradients = libSU2.get_GradientVals(objfunc_grad_filename)

# Get header and write format information
gradient_filename = gradient_filename+plotfile_ext
header,write_format = libSU2.get_GradFileFormat( 'CONTINUOUS_ADJOINT',
                                                 output_format,
                                                 Definition_DV['Kind'][0] )

# Start gradient file
Grad_File = open(gradient_filename,'w')
Grad_File.write(header)

# Write output gradients and dv information 
for i_DV in range(n_DV):
        
    # Merge information for output
    Output_Vector = [i_DV] + [gradients[i_DV]] + [FinDiff_Step] + Definition_DV['Parameters'][i_DV]
    Output_Vector = tuple(Output_Vector)

    # write row
    Grad_File.write(write_format % Output_Vector)

#: for each DV

# Remove configuration and mesh files            
os.remove(Config_CFD_filename)
os.remove(Config_GPC_filename)
os.remove(objfunc_grad_filename)
if compute_adj or compute_flow:
    os.remove(History_filename+plotfile_ext)
