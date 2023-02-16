#!/usr/bin/env python

## \file deform.py
#  \brief python package for deforming meshes
#  \author T. Lukaczyk, F. Palacios
#  \version 7.5.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import copy

from .. import io  as su2io
from .interface import DEF as SU2_DEF


# ----------------------------------------------------------------------
#  Mesh Deformation
# ----------------------------------------------------------------------

def deform ( config, dv_new=None, dv_old=None ):
    """ info = SU2.run.deform(config,dv_new=[],dv_old=[])

        Deforms mesh with:
            SU2.run.decomp()
            SU2.run.DEF()

        Assumptions:
            If optional dv_new ommitted, config is setup for deformation
            If using dv_old, must provide dv_new
            Adds 'deform' suffix to mesh output name

        Outputs:
            info - SU2 State with keys:
                HISTORY.ADJOINT_NAME
                FILES.ADJOINT_NAME

        Updates:
            config.MESH_FILENAME
            config.DV_VALUE_OLD = config.DV_VALUE_NEW

        Executes in:
            ./
    """

    if dv_new is None: dv_new = []
    if dv_old is None: dv_old = []

    # error check
    if dv_old and not dv_new: raise Exception('must provide dv_old with dv_new')

    # local copy
    konfig = copy.deepcopy(config)

    # unpack design variables
    if dv_new: konfig.unpack_dvs(dv_new,dv_old)

    # redundancy check
    if konfig['DV_VALUE_NEW'] == konfig['DV_VALUE_OLD']:
        info = su2io.State()
        info.FILES.MESH = konfig.MESH_FILENAME
        info.VARIABLES.DV_VALUE_NEW = konfig.DV_VALUE_NEW
        return info

    # setup mesh name
    suffix = 'deform'
    mesh_name = konfig['MESH_FILENAME']
    meshname_suffixed = su2io.add_suffix( mesh_name , suffix )
    konfig['MESH_OUT_FILENAME'] = meshname_suffixed

    # Run Deformation
    SU2_DEF(konfig)

    # update super config
    config.update({ 'MESH_FILENAME' : konfig['MESH_OUT_FILENAME'] ,
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
