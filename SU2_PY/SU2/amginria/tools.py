# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import numpy as np
from itertools import islice
import sys

# --- Prescribed mesh complexities, i.e. desired mesh sizes
def get_mesh_sizes(config):
    return config['ADAP_SIZES'].strip('()').split(",")
    
# --- Use the python interface to amg (YES)? Or the exe (NO)?
def get_python_amg(config):
    
    if 'ADAP_PYTHON' not in config:
        return True
    
    if config['ADAP_PYTHON'] == "YES":
        return True
    elif config['ADAP_PYTHON'] == "NO":
        return False
    else:
        sys.stderr.write("## WARNING : Invalid value for ADAP_PYTHON option. Assuming YES.\n")
        return True

# --- How many sub-iterations per mesh complexity
def get_sub_iterations(config):
    return config['ADAP_SUBITE'].strip('()').split(",")

# --- What residual reduction for each complexity level
def get_residual_reduction(config):
    if 'ADAP_RESIDUAL_REDUCTION' in config:
        return config['ADAP_RESIDUAL_REDUCTION'].strip('()').split(",")
    else:
        nRes = len(config['ADAP_SIZES'].strip('()').split(","))
        res = []
        for i in range(nRes):
            res.append(config['RESIDUAL_REDUCTION'])      
        return res 
        
# --- How many SU2 solver iterations for each complexity level
def get_adj_iter(config):
    if 'ADAP_ADJ_ITER' in config:
        return config['ADAP_ADJ_ITER'].strip('()').split(",")
    else:
        nExt_iter = len(config['ADAP_SIZES'].strip('()').split(","))
        ext_iter = []
        for i in range(nExt_iter):
            ext_iter.append(config['ITER'])        
        return ext_iter

def get_flow_iter(config):
    if 'ADAP_FLOW_ITER' in config:
        return config['ADAP_FLOW_ITER'].strip('()').split(",")
    else:
        nExt_iter = len(config['ADAP_SIZES'].strip('()').split(","))
        flow_iter = []
        for i in range(nExt_iter):
            flow_iter.append(config['ITER'])        
        return flow_iter

def get_cfl(config):
    if 'ADAP_CFL' in config:
        return config['ADAP_CFL'].strip('()').split(",")
    else:
        ncfl = len(config['ADAP_SIZES'].strip('()').split(","))
        cfl = []
        for i in range(ncfl):
            cfl.append(config['CFL_NUMBER'])        
        return cfl
    
def print_adap_options(config, kwds):
    prt = '\nMesh adaptation options:\n'
    for kwd in kwds:
        if kwd in config:
            prt += kwd + ' : ' + config[kwd] + '\n'
    prt += '\n'
    return prt
    
def get_su2_dim(filename):
    
    meshfile = open(filename,'r')
    
    def mesh_readlines(n_lines=1):
        fileslice = islice(meshfile,n_lines)
        return list(fileslice)
    
    dim = -1
    
    keepon = True
    while keepon:
        
        line = mesh_readlines()
        
        if not line: 
            keepon = False
            break
        
        # fix white space
        line = line[0]
        line = line.replace('\t',' ')
        line = line.replace('\n',' ')
    
        # skip comments
        if line[0] == "%":
            pass
    
        # number of dimensions
        elif "NDIME=" in line:
            # save to SU2_MESH data
            dim = int( line.split("=")[1].strip() )
            keepon = False
        
    return dim
        