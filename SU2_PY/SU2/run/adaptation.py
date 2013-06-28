
import os, sys, shutil, copy

from .. import io   as su2io
from .. import mesh as su2mesh
from decompose import decompose as su2decomp

def adaptation ( config , kind='' ):
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # check kind
    if kind: konfig['KIND_ADAPT'] = kind
    kind = konfig.get('KIND_ADAPT','NONE')
    if kind == 'NONE': 
        return {}
    
    # check adapted?
    
    # decompose
    su2decomp(konfig)
    
    # get adaptation function
    adapt_function = su2mesh.adapt.name_map[kind]
    
    # setup problem
    suffix = 'adapt'
    meshname_orig = konfig['MESH_FILENAME']
    meshname_new  = su2io.add_suffix( konfig['MESH_FILENAME'], suffix )
    konfig['MESH_OUT_FILENAME'] = meshname_new
    
    # Run Adaptation
    info = adapt_function(konfig)
    
    # update super config
    config['MESH_FILENAME'] = meshname_new
    config['KIND_ADAPT']    = kind
    
    # files out
    files = { 'MESH' : meshname_new }
    
    # info out
    append_nestdict( info, { 'FILES' : files } )
    
    return info


 