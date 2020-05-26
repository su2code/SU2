#!/usr/bin/env python3

## \file shape_optimizationFSI.py
#  \ Python script to perform FSI shape optimization.
#  \author Rocco Bombardieri
#  \version 7.0.2 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import os, sys, shutil
from SU2_FSI.FSI_config import OptConfig as OptConfig
from SU2_FSI.FSI_tools import SaveSplineMatrix, readConfig, readDVParam
from SU2_FSI.FSI_project import Project
from optparse import OptionParser  # use a parser for configuration
from scipy.optimize import fmin_slsqp

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

    parser=OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")



    (options, args)=parser.parse_args()
    
    # process inputs
    options.partitions  = int( options.partitions )


    FSIshape_optimization( options.filename    ,
                        options.partitions     )
    
#: main()

def FSIshape_optimization( filename                           ,
                        partitions  = 0                       ):

    # Read opt Config
    config = OptConfig(filename)  # opt configuration file



    config['NUMBER_PART'] = partitions
    config['FOLDER'] = '.'

    its               = int ( config['OPT_ITERATIONS'] )                      # number of opt iterations
    bound_upper       = float ( config['OPT_BOUND_UPPER'] )                   # variable bound to be scaled by the line search
    bound_lower       = float ( config['OPT_BOUND_LOWER'] )                   # variable bound to be scaled by the line search
    relax_factor      = float ( config['OPT_RELAX_FACTOR'] )                  # line search scale
    gradient_factor   = float ( config['OPT_GRADIENT_FACTOR'] )               # objective function and gradient scale
    opt_constraint    = config['OPT_CONSTRAINT']
    folder = config['FOLDER']
    accu              = float ( config['OPT_ACCURACY'] ) * gradient_factor    # optimizer accuracy
    # n_dv
    dv_kind_str = readConfig(config['CONFIG_DEF'], 'DV_KIND')
    n_dv = int(len(dv_kind_str.split(',')))
    x0                = [0.0]*n_dv # initial design
    xb_low            = [float(bound_lower)/float(relax_factor)]*n_dv      # lower dv bound it includes the line search acceleration factor
    xb_up             = [float(bound_upper)/float(relax_factor)]*n_dv      # upper dv bound it includes the line search acceleration fa
    xb                = list(zip(xb_low, xb_up)) # design bounds

    # Instantiate project object
    project = Project(config)
    
    
    # call optimization function  (nned to define some high level elements yet....)
    slsqp(project,x0,xb,its,accu)

    return 

def slsqp(project,x0=None,xb=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_slsqp(project,x0=[],xb=[],its=100,accu=1e-10)
    
        Runs the Scipy implementation of SLSQP with 
        an SU2 project
        
        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10
        
        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # Performing spline
    print('Splining commented')
    #SaveSplineMatrix(project.config)
          
    # function handles
    func           = obj_f
    f_eqcons       = con_ceq
    f_ieqcons      = con_cieq 
    
    # gradient handles
    fprime         = obj_df
    fprime_eqcons  = con_dceq
    fprime_ieqcons = con_dcieq        
    
    # number of design variables  (this is read from the DEF input file)
    n_dv = len(x0)
    project.n_dv = n_dv
        
    # prescale x0 # Not interested at the moment
    #dv_scales = project.config['DEFINITION_DV']['SCALE']
    #k = 0
    #for i, dv_scl in enumerate(dv_scales):
    #    for j in range(dv_size[i]):
    #        x0[k] =x0[k]/dv_scl;
    #        k = k + 1

    # scale accuracy # further scaling of the accuracy (not interested now)
    #obj = project.config['OPT_OBJECTIVE']
    #obj_scale = []
    #for this_obj in obj.keys():
    #    obj_scale = obj_scale + [obj[this_obj]['SCALE']]
    # 
    # Only scale the accuracy for single-objective problems: 
    #if len(obj.keys())==1:
    #    accu = accu*obj_scale[0]
       
    # obj scale
    # I need to look into the Adjoint flow file
    # Adjoint config
    print(project.configFSIAdjoint['SU2_CONFIG'])
    adjointflow_config = project.configFSIAdjoint['SU2_CONFIG']
    line = readConfig(adjointflow_config, 'OBJECTIVE_FUNCTION')
    if line.find('*') != -1:
       line = line.split("*", 1)
       obj_scale = line[1].strip()
    else:
       obj_scale = float(1) 
       
    # scale accuracy: The step size for finite-difference derivative estimates.
    eps = 1.0e-04

    # optimizer summary
    sys.stdout.write('Sequential Least SQuares Programming (SLSQP) parameters:\n')
    sys.stdout.write('Number of design variables: ' + ' ( ' + str(n_dv) + ' ) \n' )
    sys.stdout.write('Objective function scaling factor: ' + str(obj_scale) + '\n')
    sys.stdout.write('Maximum number of iterations: ' + str(its) + '\n')
    sys.stdout.write('Requested accuracy: ' + str(accu) + '\n')
    sys.stdout.write('Initial guess for the independent variable(s): ' + str(x0) + '\n')
    sys.stdout.write('Lower and upper bound for each independent variable: ' + str(xb) + '\n\n')

    # Run Optimizer
    outputs = fmin_slsqp( x0             = x0             ,
                          func           = func           , 
                          f_eqcons       = f_eqcons       , 
                          f_ieqcons      = f_ieqcons      ,
                          fprime         = fprime         ,
                          fprime_eqcons  = fprime_eqcons  , 
                          fprime_ieqcons = fprime_ieqcons , 
                          args           = (project,)     , 
                          bounds         = xb             ,
                          iter           = its            ,
                          iprint         = 2              ,
                          full_output    = True           ,
                          acc            = accu           ,
                          epsilon        = eps            )
    
    # Done
    return outputs


def obj_f(x,project):
    """ obj = obj_f(x,project)
        
        Objective Function
        SU2 Project interface to scipy.fmin_slsqp
        
        scipy_slsqp: minimize f(x), float
    """  
    print('obj_f: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))    
    obj = project.obj_f(x)   
    print('After obj_f: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))   
    return obj

def obj_df(x,project):
    """ dobj = obj_df(x,project)
        
        Objective Function Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        scipy_slsqp: df(x), ndarray[dim]
    """    
    print('obj_df: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))     
    obj_df = project.obj_df(x)   
    print('After obj_df: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))       
    return obj_df

def con_ceq(x,project):
    """ cons = con_ceq(x,project)
        
        Equality Constraint Functions
    """
    if project.design_iter !=-1: print('con_ceq: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))      
    cons = project.con_ceq(x)  
    print('After con_ceq: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))      
    return cons

def con_dceq(x,project):
    """ dcons = con_dceq(x,project)
        
        Equality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        scipy_slsqp: dceq(x), ndarray[nceq x dim]
    """
    print('con_dceq: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))      
    dcons = project.con_dceq(x)
    print('After con_dceq: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))      
    return dcons

def con_cieq(x,project):
    """ cons = con_cieq(x,project)
        
        Inequality Constraints

    """ 
    if project.design_iter !=-1: print('con_cieq: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))      
    cons = project.con_cieq(x)  
    print('After con_cieq: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))          
    return cons
    
    
def con_dcieq(x,project):
    """ dcons = con_dcieq(x,project)
        
        Inequality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        scipy_slsqp: dcieq(x), ndarray[ncieq x dim]
    """
    print('con_dcieq: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))     
    dcons = project.con_dcieq(x)
    print('After con_dcieq: design {} x ={}'.format(project.design_iter,project._design[project.design_iter].getdv()))      
    return dcons












# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

