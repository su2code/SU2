## \file deform.py
#  \brief python package for deforming meshes
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy

from .. import io  as su2io
from decompose import decompose as su2decomp
from interface import MDC as SU2_MDC


# ----------------------------------------------------------------------
#  Mesh Deformation
# ----------------------------------------------------------------------

def deform ( config, dv_new=None, dv_old=None ):
    """ info = SU2.run.deform(config,dv_new=[],dv_old=[])
        
        Deforms mesh with:
            SU2.run.decomp()
            SU2.run.MDC()
            
        Assumptions:
            Redundant decomposition if config.DECOMPOSED == True
            If optional dv_new ommitted, config is setup for deformation
            If using dv_old, must provide dv_new
            Adds 'deform' suffix to mesh output name
            
        Outputs:
            info - SU2 State with keys:
                HISTORY.ADJOINT_NAME
                FILES.ADJOINT_NAME
                
        Updates:
            config.DECOMPOSED
            config.MESH_FILENAME
            config.DV_VALUE_OLD = config.DV_VALUE_NEW
            
        Executes in:
            ./
    """    
    
    if dv_new is None: dv_new = []
    if dv_old is None: dv_old = []
    
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