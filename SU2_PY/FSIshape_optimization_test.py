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
import numpy as np
from SU2_FSI.FSI_config import OptConfig as OptConfig
from SU2_FSI.FSI_tools import SaveSplineMatrix, readConfig, readDVParam
from SU2_FSI.FSI_project import Project
from optparse import OptionParser  # use a parser for configuration
from scipy.optimize import fmin_slsqp

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------


def main(filename= 'config_opt.cfg'     , partitions  = 1                       ):

    # Read opt Config
    config = OptConfig(filename)  # opt configuration file



    config['NUMBER_PART'] = partitions
    config['FOLDER'] = '.'

    # n_dv
    dv_kind_str = readConfig(config['CONFIG_DEF'], 'DV_KIND')
    n_dv = int(len(dv_kind_str.split(',')))
    x0                = np.array([0.0]*n_dv) # initial design

    # Instantiate project object
    project = Project(config)
    
    # number of design variables  (this is read from the DEF input file)
    n_dv = len(x0)
    project.n_dv = n_dv    

    # replichiamo comportamento codice
    con_ceq(x0,project)
    con_cieq(x0,project)
    obj_f(x0,project)
    con_ceq(x0,project)
    con_cieq(x0,project)    
    obj_df(x0,project)    
    x1 = np.array([1.0]*n_dv)
    con_ceq(x1,project)    
    return 


def obj_f(x,project):
    """ obj = obj_f(x,project)
        
        Objective Function
        SU2 Project interface to scipy.fmin_slsqp
        
        scipy_slsqp: minimize f(x), float
    """  
    obj = project.obj_f(x)   

    return obj

def obj_df(x,project):
    """ dobj = obj_df(x,project)
        
        Objective Function Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        scipy_slsqp: df(x), ndarray[dim]
    """    
    
    obj_df = project.obj_df(x)

    
    return obj_df

def con_ceq(x,project):
    """ cons = con_ceq(x,project)
        
        Equality Constraint Functions
    """
    
    cons = project.con_ceq(x)
    
        
    return cons

def con_dceq(x,project):
    """ dcons = con_dceq(x,project)
        
        Equality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        scipy_slsqp: dceq(x), ndarray[nceq x dim]
    """
    
    dcons = project.con_dceq(x)
    
    return dcons

def con_cieq(x,project):
    """ cons = con_cieq(x,project)
        
        Inequality Constraints

    """    
    cons = project.con_cieq(x)
    return cons
    
    
def con_dcieq(x,project):
    """ dcons = con_dcieq(x,project)
        
        Inequality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        scipy_slsqp: dcieq(x), ndarray[ncieq x dim]
    """
    
    dcons = project.con_dcieq(x)
        
    return dcons












# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

