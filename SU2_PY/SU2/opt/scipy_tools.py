
import os, sys, shutil, copy

from .. import eval as su2eval
from scipy.optimize import fmin_slsqp
from numpy import array, zeros
from numpy.linalg import norm

def scipy_slsqp(project,x0=[],xb=[],its=20):
    
    n_dv = len( project.config['DEFINITION_DV']['KIND'] )
    project.n_dv = n_dv
    
    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    
    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    x0 = [ x0[i]/dv_scl for i,dv_scl in enumerate(dv_scales) ]    
        
    # Run Optimizer
    outputs = fmin_slsqp( x0             = x0         , 
                          func           = obj_f      , 
                          f_eqcons       = con_ceq    , 
                          f_ieqcons      = con_cieq   , 
                          fprime         = obj_df     ,
                          fprime_eqcons  = con_dceq   , 
                          fprime_ieqcons = con_dcieq  , 
                          args           = (project,) , 
                          bounds         = xb         ,
                          iter           = its        ,
                          iprint         = 2          ,
                          full_output    = 1          ,
                          acc            = 1e-10      ,
                          epsilon        = 1.0e-06     )
    
    # Done
    return outputs
    
    
def obj_f(x,project):
    """ su2:         minimize f(x), list[nobj]
        scipy_slsqp: minimize f(x), float
    """
    
    print ""
    
    obj = project.obj_f(x)
    
    obj = obj[0]
    
    return obj

def obj_df(x,project):
    """ su2:         df(x), list[nobj x dim]
        scipy_slsqp: df(x), ndarray[dim]
    """    
    
    dobj = project.obj_df(x)
    
    dobj = array( dobj[0] )
    
    return dobj

def con_ceq(x,project):
    """ su2:         ceq(x) = 0.0, list[nceq]
        scipy_slsqp: ceq(x) = 0.0, ndarray[nceq]
    """
    
    cons = project.con_ceq(x)
    
    if cons: cons = array(cons)
    else:    cons = zeros([0])
        
    return cons

def con_dceq(x,project):
    """ su2:         dceq(x), list[nceq x dim]
        scipy_slsqp: dceq(x), ndarray[nceq x dim]
    """
    
    dcons = project.con_dceq(x)

    dim = project.n_dv
    if dcons: dcons = array(dcons)
    else:     dcons = zeros([0,dim])
    
    return dcons

def con_cieq(x,project):
    """ su2:         cieq(x) < 0.0, list[ncieq]
        scipy_slsqp: cieq(x) > 0.0, ndarray[ncieq]
    """
    
    cons = project.con_cieq(x)
    
    if cons: cons = array(cons)
    else:    cons = zeros([0])
    
    return -cons
    
def con_dcieq(x,project):
    """ su2:         dcieq(x), list[ncieq x dim]
        scipy_slsqp: dcieq(x), ndarray[ncieq x dim]
    """
    
    dcons = project.con_dcieq(x)
    
    dim = project.n_dv
    if dcons: dcons = array(dcons)
    else:     dcons = zeros([0,dim])
    
    return -dcons