
import os, sys, shutil, copy

from .. import io  as su2io
from ..run import decompose as su2decomp
from ..run import CFD as SU2_CFD
from ..run import MAC as SU2_MAC

def Full(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def FullFlow(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def GradFlow(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def FullAdjoint(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def GradAdjoint(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def GradFlowAdj(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def Robust(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def FullLinear(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def Computable(config):
    config = copy.deepcopy(config)
    
    pass
 
def ComputableRobust():
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def Remaining(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def Wake(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError
 
def HorizontalPlane(config):
    config = copy.deepcopy(config)
    
    raise NotImplementedError

# config name map
name_map = { 'FULL'              : Full             ,
             'FULL_FLOW'         : FullFlow         ,
             'GRAD_FLOW'         : GradFlow         ,
             'FULL_ADJOINT'      : FullAdjoint      ,
             'GRAD_ADJOINT'      : GradAdjoint      ,
             'GRAD_FLOW_ADJ'     : GradFlowAdj      ,
             'ROBUST'            : Robust           ,
             'FULL_LINEAR'       : FullLinear       ,
             'COMPUTABLE'        : Computable       ,
             'COMPUTABLE_ROBUST' : ComputableRobust ,
             'REMAINING'         : Remaining        ,
             'WAKE'              : Wake             ,
             'HORIZONTAL_PLANE'  : HorizontalPlane   }
