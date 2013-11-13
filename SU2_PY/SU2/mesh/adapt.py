
import os, sys, shutil, copy

from .. import io  as su2io
from ..run import decompose as su2decomp
from ..run import CFD as SU2_CFD
from ..run import MAC as SU2_MAC

def full(config):
    config = copy.deepcopy(config)
    
    config.KIND_ADAPT = 'FULL'
    
    raise NotImplementedError

 
def full_flow(config):
    
    # local copy
    konfig = copy.deepcopy(config)

    # decompose
    su2decomp(konfig)    
    
    # set config
    konfig.KIND_ADAPT = 'FULL_FLOW'
    
    # run MAC
    SU2_MAC(konfig)
    
    return
    
    
def full_adjoint(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError

def grad_flow(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def grad_adjoint(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def grad_flow_adj(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def robust(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def full_linear(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def computable(config):
    config = copy.deepcopy(config)
    
    pass
 
def computable_robust():
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def remaining(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def wake(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def horizontal_plane(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError

# config name map
name_map = { 'FULL'              : full              ,
             'FULL_FLOW'         : full_adjoint      ,
             'GRAD_FLOW'         : grad_flow         ,
             'FULL_ADJOINT'      : full_adjoint      ,
             'GRAD_ADJOINT'      : grad_adjoint      ,
             'GRAD_FLOW_ADJ'     : grad_flow_adj     ,
             'ROBUST'            : robust            ,
             'FULL_LINEAR'       : full_linear       ,
             'COMPUTABLE'        : computable        ,
             'COMPUTABLE_ROBUST' : computable_robust ,
             'REMAINING'         : remaining         ,
             'WAKE'              : wake              ,
             'HORIZONTAL_PLANE'  : horizontal_plane   }
