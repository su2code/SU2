#!/usr/bin/env python

## \file interface.py
#  \brief Wrapper functions for interfacing with the Inria AMG library
#  \author Victorien Menier, Brian Mungu\'ia
#  \version 7.3.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

import numpy as np

import su2gmf
import pyamg

def call_pyamg(mesh, config):
    """Adapt mesh using pyamg module"""
    remesh_options = {}

    if 'hgrad' in config: remesh_options['gradation'] = config['hgrad']

    remesh_options['logfile'] = config['amg_log']
    remesh_options['options'] = config['options']

    Dim = mesh['dimension']

    if Dim == 2 :
        Ver = mesh['xyz']
        mesh['xy'] = np.stack((Ver[:,0],Ver[:,1]), axis=1)

        del mesh['xyz']

    if 'xy' in mesh:  mesh['xy']  = mesh['xy'].tolist()
    if 'xyz' in mesh: mesh['xyz'] = mesh['xyz'].tolist()

    if 'Corners' in mesh:    mesh['Corners']    = mesh['Corners'].tolist()
    if 'Edges' in mesh:      mesh['Edges']      = mesh['Edges'].tolist()
    if 'Triangles' in mesh:  mesh['Triangles']  = mesh['Triangles'].tolist()
    if 'Tetrahedra' in mesh: mesh['Tetrahedra'] = mesh['Tetrahedra'].tolist()

    if 'metric' in mesh: mesh['metric'] = mesh['metric'].tolist()
    if 'sensor' in mesh: mesh['sensor'] = mesh['sensor'].tolist()

    #--- Give pyamg these parameters in case metric intersection violates hmax, hmin, or target
    if 'sensor' in mesh:
        remesh_options['Lp']     = config['Lp']
        remesh_options['hmax']   = config['hmax']
        remesh_options['hmin']   = config['hmin']
        remesh_options['target'] = config['size']

    try:
        mesh_new = pyamg.adapt_mesh(mesh, remesh_options)
    except:
        raise RuntimeError("pyamg failed.")

    return mesh_new

def read_mesh_and_sol(mesh_name, solution_name):
    """Read mesh and solution using su2gmf module"""
    Ver = []
    Tri = []
    Tet = []
    Edg = []
    Cor = []
    Hex = []
    Pyr = []
    Pri = []
    Qua = []
    Sol = []
    SolTag = []

    Markers = []

    su2gmf.ReadMeshAndSol(mesh_name, solution_name, Ver, Cor, Tri, Tet, Edg, Hex, Qua, Pyr, Pri, Sol, SolTag, Markers)

    NbrTet = int(len(Tet)/5)
    Tet = np.reshape(Tet,(NbrTet, 5)).astype(np.int64)

    NbrTri = int(len(Tri)/4)
    Tri = np.reshape(Tri,(NbrTri, 4)).astype(np.int64)

    NbrEdg = int(len(Edg)/3)
    Edg = np.reshape(Edg,(NbrEdg, 3)).astype(np.int64)

    NbrCor = int(len(Cor))
    Cor = np.reshape(Cor,(NbrCor, 1)).astype(np.int64)

    NbrVer = int(len(Ver)/3)
    Ver = np.reshape(Ver,(NbrVer, 3))

    SolSiz = int(len(Sol)/NbrVer)
    Sol = np.reshape(Sol,(NbrVer, SolSiz))

    #--- First row of Markers contains dimension
    Dim = int(Markers[0])

    mesh = dict()

    mesh['dimension']  = Dim

    mesh['xyz']        = Ver

    mesh['Triangles']  = Tri
    mesh['Tetrahedra'] = Tet
    mesh['Edges']      = Edg
    mesh['Corners']    = Cor
    mesh['solution']   = Sol

    mesh['solution_tag'] = SolTag

    mesh['id_solution_tag'] = dict()
    for i, tag in enumerate(SolTag):
        mesh['id_solution_tag'][tag] = i

    mesh['markers'] = Markers

    return mesh

def read_mesh(mesh_name):
    """Read mesh using su2gmf module"""
    Ver = []
    Tri = []
    Tet = []
    Edg = []
    Cor = []
    Hex = []
    Pyr = []
    Pri = []
    Qua = []

    Markers = []

    su2gmf.ReadMesh(mesh_name, Ver, Cor, Tri, Tet, Edg, Hex, Qua, Pyr, Pri, Markers)

    NbrTet = int(len(Tet)/5)
    Tet = np.reshape(Tet,(NbrTet, 5)).astype(np.int64)

    NbrTri = int(len(Tri)/4)
    Tri = np.reshape(Tri,(NbrTri, 4)).astype(np.int64)

    NbrEdg = int(len(Edg)/3)
    Edg = np.reshape(Edg,(NbrEdg, 3)).astype(np.int64)

    NbrCor = int(len(Cor))
    Cor = np.reshape(Cor,(NbrCor, 1)).astype(np.int64)

    NbrVer = int(len(Ver)/3)
    Ver = np.reshape(Ver,(NbrVer, 3))

    #--- First row of Markers contains dimension
    Dim = int(Markers[0])

    mesh = dict()

    mesh['dimension']  = Dim

    mesh['xyz']        = Ver

    mesh['Triangles']  = Tri
    mesh['Tetrahedra'] = Tet
    mesh['Edges']      = Edg
    mesh['Corners']    = Cor

    mesh['markers'] = Markers

    return mesh

def read_sol(solution_name, mesh):
    """Read solution using su2gmf module"""
    NbrVer = len(mesh['xyz'])
    Dim    = mesh['dimension']

    Sol = []
    SolTag = []

    su2gmf.ReadSol(solution_name, Sol, SolTag, NbrVer, Dim)

    SolSiz = int(len(Sol)/NbrVer)
    Sol = np.array(Sol).reshape(NbrVer,SolSiz).tolist()

    sol = dict()

    sol['solution']     = Sol
    sol['solution_tag'] = SolTag

    sol['id_solution_tag'] = dict()
    for i, tag in enumerate(SolTag):
        sol['id_solution_tag'][tag] = i

    return sol

def write_mesh_and_sol(mesh_name, solution_name, mesh):
    """Write mesh and solution using su2gmf module"""
    Tri     = []
    Tet     = []
    Edg     = []
    Cor     = []
    Hex     = []
    Pyr     = []
    Pri     = []
    Qua     = []
    Sol     = []
    Markers = []
    Dim     = 3
    Ver     = []
    SolTag  = []

    if 'Triangles' in mesh:    Tri     = mesh['Triangles']
    if 'Tetrahedra' in mesh:   Tet     = mesh['Tetrahedra']
    if 'Edges' in mesh:        Edg     = mesh['Edges']
    if 'Corners' in mesh:      Cor     = mesh['Corners']
    if 'solution' in mesh:     Sol     = mesh['solution']
    if 'markers' in mesh:      Markers = mesh['markers']
    if 'dimension' in mesh:    Dim     = mesh['dimension']
    if 'solution_tag' in mesh: SolTag  = mesh['solution_tag']
    if 'xyz' in mesh:
        Ver = mesh['xyz']
        Ver = np.array(Ver).reshape(3*len(Ver)).tolist()
    elif 'xy' in mesh:
        Ver = np.array(mesh['xy'])
        z = np.zeros(len(mesh['xy']))
        Ver = np.c_[Ver, z]
        Ver = np.array(Ver).reshape(3*len(Ver)).tolist()

    Tri = np.array(Tri).reshape(4*len(Tri)).tolist()
    Tet = np.array(Tet).reshape(5*len(Tet)).tolist()
    Edg = np.array(Edg).reshape(3*len(Edg)).tolist()
    Cor = np.array(Cor).reshape(len(Cor)).tolist()

    if len(Sol) > 1 :
        SolSiz = len(Sol[1])
        Sol = np.array(Sol).reshape(SolSiz*len(Sol)).tolist()
    else:
        Sol = []

    su2gmf.WriteMeshAndSol(mesh_name, solution_name, Ver, Cor, Tri, Tet, Edg, Hex, Qua, Pyr, Pri, Sol, SolTag, Markers, Dim)

def write_mesh(mesh_name, mesh):
    """Write mesh using su2gmf module"""
    Tri     = []
    Tet     = []
    Edg     = []
    Cor     = []
    Hex     = []
    Pyr     = []
    Pri     = []
    Qua     = []
    Markers = []
    Dim     = 3
    Ver     = []

    if 'Triangles' in mesh:  Tri     = mesh['Triangles']
    if 'Tetrahedra' in mesh: Tet     = mesh['Tetrahedra']
    if 'Edges' in mesh:      Edg     = mesh['Edges']
    if 'Corners' in mesh:    Cor     = mesh['Corners']
    if 'markers' in mesh:    Markers = mesh['markers']
    if 'dimension' in mesh:  Dim     = mesh['dimension']
    if 'xyz' in mesh:
        Ver = mesh['xyz']
        Ver = np.array(Ver).reshape(3*len(Ver)).tolist()
    elif 'xy' in mesh:
        Ver = np.array(mesh['xy'])
        z = np.zeros(len(mesh['xy']))
        Ver = np.c_[Ver, z]
        Ver = np.array(Ver).reshape(3*len(Ver)).tolist()

    Tri = np.array(Tri).reshape(4*len(Tri)).tolist()
    Tet = np.array(Tet).reshape(5*len(Tet)).tolist()
    Edg = np.array(Edg).reshape(3*len(Edg)).tolist()
    Cor = np.array(Cor).reshape(len(Cor)).tolist()

    su2gmf.WriteMesh(mesh_name, Ver, Cor, Tri, Tet, Edg, Hex, Qua, Pyr, Pri, Markers, Dim)

def write_sol(sol_name, sol):
    """Write solution using su2gmf module"""
    Dim = sol['dimension']
    Sol = sol['solution']

    if 'xyz' in sol:
        NbrVer = int(len(sol['xyz']))
        Ver = sol['xyz']
        Ver = np.array(Ver).reshape(3*len(Ver)).tolist()
    elif 'xy' in sol:
        NbrVer = int(len(sol['xy']))
        Ver = np.array(sol['xy'])
        z = np.zeros(len(sol['xy']))
        Ver = np.c_[Ver, z]
        Ver = np.array(Ver).reshape(3*len(Ver)).tolist()

    SolTag = sol['solution_tag']

    if len(Sol) > 1 :
        SolSiz = len(Sol[1])
        Sol = np.array(Sol).reshape(SolSiz*len(Sol)).tolist()
    else:
        raise RuntimeError("Empty solution given to write_sol.")

    su2gmf.WriteSol(sol_name, Ver, Sol, SolTag, NbrVer, Dim)
