
import os, sys, shutil, copy

from .. import io  as su2io
from decompose import decompose as su2decomp
from interface import MDC as SU2_MDC

def deform ( config, dv_new=[], dv_old=[] ): # cleanup?
    ''' if optional dv_new ommitted, assume config is setup for deformation
        if using dv_old, must provide dv_new
    '''
    
    # error check
    if dv_old and not dv_new: raise Exception, 'must provide dv_old with dv_new'
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # decompose
    su2decomp(konfig)
    
    # unpack design variables
    if dv_new: konfig.unpack_dvs(dv_new,dv_old)
    
    # redundancy check
    if konfig['DV_VALUE_NEW'] == konfig['DV_VALUE_OLD']:
        info = su2io.State()
        info.VARIABLES.DV_VALUE_NEW = konfig.DV_VALUE_NEW        
        return info
    
    # setup mesh name
    suffix = 'deform'
    mesh_name = konfig['MESH_FILENAME']
    meshname_suffixed = su2io.add_suffix( mesh_name , suffix )
    konfig['MESH_OUT_FILENAME'] = meshname_suffixed
    
    # Run Deformation
    SU2_MDC(konfig)
    
    # update super config
    config.update({ 'DECOMPOSED'    : konfig['DECOMPOSED']        ,
                    'MESH_FILENAME' : konfig['MESH_OUT_FILENAME'] , 
                    'DV_KIND'       : konfig['DV_KIND']           ,
                    'DV_MARKER'     : konfig['DV_MARKER']         ,
                    'DV_PARAM'      : konfig['DV_PARAM']          ,
                    'DV_VALUE_OLD'  : konfig['DV_VALUE_NEW']      ,
                    'DV_VALUE_NEW'  : konfig['DV_VALUE_NEW']      })
    # not modified: config['MESH_OUT_FILENAME']
    
    # info out
    info = su2io.State()
    info.FILES.MESH = meshname_suffixed
    info.VARIABLES.DV_VALUE_NEW = konfig.DV_VALUE_NEW
    
    return info

#: def deform()