#!/usr/bin/env python 

## \file scipy_tools.py
#  \brief tools for interfacing with scipy
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
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

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, shutil, copy

from .. import eval as su2eval
from numpy import array, zeros
from numpy.linalg import norm


# -------------------------------------------------------------------
#  Scipy SLSQP
# -------------------------------------------------------------------

def scipy_slsqp(project,x0=None,xb=None,its=100,grads=True):
    """ result = scipy_slsqp(project,x0=[],xb=[],its=100)
    
        Runs the Scipy implementation of SLSQP with 
        an SU2 project
        
        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds
            its     - max outer iterations, default 100
        
        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_slsqp

    # handle input cases
    if x0 is None: x0 = []
    if xb is None: xb = []
    
    # function handles
    func           = obj_f
    f_eqcons       = con_ceq
    f_ieqcons      = con_cieq 
    
    # gradient handles
    if project.config.GRADIENT_METHOD == 'NONE': 
        fprime         = None
        fprime_eqcons  = None
        fprime_ieqcons = None
    else:
        fprime         = obj_df
        fprime_eqcons  = con_dceq
        fprime_ieqcons = con_dcieq        
    
    # number of design variables
    n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    project.n_dv = n_dv
    
    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    
    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]    
        
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
                          full_output    = 1              ,
                          acc            = 1e-10          ,
                          epsilon        = 1.0e-06         )
    
    # Done
    return outputs
    
    
def obj_f(x,project):
    """ obj = obj_f(x,project)
        
        Objective Function
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         minimize f(x), list[nobj]
        scipy_slsqp: minimize f(x), float
    """
    
    print ""
    
    obj = project.obj_f(x)
    
    obj = obj[0]
    
    return obj

def obj_df(x,project):
    """ dobj = obj_df(x,project)
        
        Objective Function Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         df(x), list[nobj x dim]
        scipy_slsqp: df(x), ndarray[dim]
    """    
    
    dobj = project.obj_df(x)
    
    dobj = array( dobj[0] )
    
    return dobj

def con_ceq(x,project):
    """ cons = con_ceq(x,project)
        
        Equality Constraint Functions
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         ceq(x) = 0.0, list[nceq]
        scipy_slsqp: ceq(x) = 0.0, ndarray[nceq]
    """
    
    cons = project.con_ceq(x)
    
    if cons: cons = array(cons)
    else:    cons = zeros([0])
        
    return cons

def con_dceq(x,project):
    """ dcons = con_dceq(x,project)
        
        Equality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         dceq(x), list[nceq x dim]
        scipy_slsqp: dceq(x), ndarray[nceq x dim]
    """
    
    dcons = project.con_dceq(x)

    dim = project.n_dv
    if dcons: dcons = array(dcons)
    else:     dcons = zeros([0,dim])
    
    return dcons

def con_cieq(x,project):
    """ cons = con_cieq(x,project)
        
        Inequality Constraints
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         cieq(x) < 0.0, list[ncieq]
        scipy_slsqp: cieq(x) > 0.0, ndarray[ncieq]
    """
    
    cons = project.con_cieq(x)
    
    if cons: cons = array(cons)
    else:    cons = zeros([0])
    
    return -cons
    
def con_dcieq(x,project):
    """ dcons = con_dcieq(x,project)
        
        Inequality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         dcieq(x), list[ncieq x dim]
        scipy_slsqp: dcieq(x), ndarray[ncieq x dim]
    """
    
    dcons = project.con_dcieq(x)
    
    dim = project.n_dv
    if dcons: dcons = array(dcons)
    else:     dcons = zeros([0,dim])
    
    return -dcons
