
import os, sys, shutil, copy
from .. import io  as su2io
from interface import DDC as SU2_DDC

def decompose ( config ):
    
    # local copy
    konfig = copy.deepcopy(config)
    
    # check if needed
    partitions = konfig['NUMBER_PART']
    decomposed = konfig.get('DECOMPOSED',False)
    if partitions <= 1 or decomposed:
        return su2io.State()
    
    # Run Decomposition
    SU2_DDC(konfig)
    
    # update config super copy
    config['DECOMPOSED'] = True
    
    # info out
    info = su2io.State()
    info.FILES.MESH = config['MESH_FILENAME']
    
    return info
