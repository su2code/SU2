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

import os, sys
import numpy as np

import _amgio as amgio
import pyamg

def call_pyamg(mesh, config):

    remesh_options                = {}

    if 'hgrad' in config: remesh_options['gradation'] = config['hgrad']
    
    remesh_options['logfile']     = config['amg_log']
    remesh_options['options']     = config['options']
    
    Dim = mesh['dimension']
    
    if Dim == 2 :
        Ver = mesh['xyz']
        mesh['xy'] = np.stack((Ver[:,0],Ver[:,1]), axis=1)
        
        del mesh['xyz']
    
    ''' 
      TO ADD: 
     {'adap_back' 'hmax' 'hmin'
      'sol_in': 'current_sensor.solb', 'sol_itp_in': 'current.solb', 'metric_in': '', 'adap_source': '', 
     'mesh_in': 'current.meshb', 'mesh_out': 'current.new.meshb'}
    ''' 
    
    if 'xy' in mesh:    mesh['xy']  = mesh['xy'].tolist()
    if 'xyz' in mesh:   mesh['xyz'] = mesh['xyz'].tolist()
    
    if 'Corners' in mesh:    mesh['Corners']    = mesh['Corners'].tolist() 
    if 'Edges' in mesh:      mesh['Edges']      = mesh['Edges'].tolist() 
    if 'Triangles' in mesh:  mesh['Triangles']  = mesh['Triangles'].tolist()
    if 'Tetrahedra' in mesh: mesh['Tetrahedra'] = mesh['Tetrahedra'].tolist()   

    if 'metric' in mesh: mesh['metric'] = mesh['metric'].tolist()
    if 'sensor' in mesh:
        mesh['sensor']           = mesh['sensor'].tolist()

    #--- Give pyamg these parameters in case metric intersection violates hmax, hmin, or target
    remesh_options['Lp']     = config['Lp']
    remesh_options['hmax']   = config['hmax']
    remesh_options['hmin']   = config['hmin']
    remesh_options['target'] = config['size']

    try:
        mesh_new = pyamg.adapt_mesh(mesh, remesh_options)        
    except:
        sys.stderr("## ERROR : pyamg failed.\n")
        raise
    
    return mesh_new
    
    
# --- Read mesh and solution using amgio module
def read_mesh_and_sol(mesh_name, solution_name):
    
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
    
    amgio.py_ReadMeshAndSol(mesh_name, solution_name, Ver, Cor, Tri, Tet, Edg, Hex, Qua, Pyr, Pri, Sol, SolTag, Markers)
        
    NbrTet = int(len(Tet)/5)
    Tet = np.reshape(Tet,(NbrTet, 5)).astype(int)
    
    NbrTri = int(len(Tri)/4)
    Tri = np.reshape(Tri,(NbrTri, 4)).astype(int)
    
    NbrEdg = int(len(Edg)/3)
    Edg = np.reshape(Edg,(NbrEdg, 3)).astype(int)

    NbrCor = int(len(Cor))
    Cor = np.reshape(Cor,(NbrCor, 1)).astype(int)

    NbrVer = int(len(Ver)/3)
    Ver = np.reshape(Ver,(NbrVer, 3))
    
    SolSiz = int(len(Sol)/NbrVer)
    Sol = np.array(Sol).reshape(NbrVer,SolSiz).tolist()
    
    # First row of Markers contains dimension
    Dim = int(Markers[0])
    
    mesh = dict()
    
    mesh['dimension']    = Dim
    
    mesh['xyz']          = Ver 
    
    mesh['Triangles']    = Tri
    mesh['Tetrahedra']   = Tet
    mesh['Edges']        = Edg
    mesh['Corners']      = Cor
    mesh['solution']     = Sol
    
    mesh['solution_tag'] = SolTag
    
    mesh['id_solution_tag'] = dict()
    for i in range(len(SolTag)):
        mesh['id_solution_tag'][SolTag[i]] = i
        
    mesh['markers'] = Markers    
    
    return mesh

# --- Read mesh using amgio module
def read_mesh(mesh_name):
    
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
    
    amgio.py_ReadMesh(mesh_name, Ver, Cor, Tri, Tet, Edg, Hex, Qua, Pyr, Pri, Markers)
        
    NbrTet = int(len(Tet)/5)
    Tet = np.reshape(Tet,(NbrTet, 5)).astype(int)
    
    NbrTri = int(len(Tri)/4)
    Tri = np.reshape(Tri,(NbrTri, 4)).astype(int)
    
    NbrEdg = int(len(Edg)/3)
    Edg = np.reshape(Edg,(NbrEdg, 3)).astype(int)

    NbrCor = int(len(Cor))
    Cor = np.reshape(Cor,(NbrCor, 1)).astype(int)

    NbrVer = int(len(Ver)/3)
    Ver = np.reshape(Ver,(NbrVer, 3))
    
    # First row of Markers contains dimension
    Dim = int(Markers[0])
    
    mesh = dict()
    
    mesh['dimension']    = Dim
    
    mesh['xyz']          = Ver 
    
    mesh['Triangles']    = Tri
    mesh['Tetrahedra']   = Tet
    mesh['Edges']        = Edg
    mesh['Corners']      = Cor
        
    mesh['markers'] = Markers    
    
    return mesh

# --- Read mesh using amgio module
def read_sol(solution_name, mesh):
    
    NbrVer = len(mesh['xyz'])
    Dim    = mesh['dimension']

    Sol = []
    SolTag = []
    
    amgio.py_ReadSol(solution_name, Sol, SolTag, NbrVer, Dim)
        
    SolSiz = int(len(Sol)/NbrVer)
    Sol = np.array(Sol).reshape(NbrVer,SolSiz).tolist()
    
    sol = dict()

    sol['solution']     = Sol
    sol['solution_tag'] = SolTag
    
    sol['id_solution_tag'] = dict()
    for i in range(len(SolTag)):
        sol['id_solution_tag'][SolTag[i]] = i
    
    return sol
    
# --- Write mesh and solution using amgio module
def write_mesh_and_sol(mesh_name, solution_name, mesh):
    
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
        
    if 'Triangles' in mesh:     Tri     = mesh['Triangles']
    if 'Tetrahedra' in mesh:    Tet     = mesh['Tetrahedra']
    if 'Edges' in mesh:         Edg     = mesh['Edges']
    if 'Corners' in mesh:       Cor     = mesh['Corners']
    if 'solution' in mesh:      Sol     = mesh['solution']
    if 'markers' in mesh:       Markers = mesh['markers']
    if 'dimension' in mesh:     Dim     = mesh['dimension']
    if 'solution_tag' in mesh:  SolTag  = mesh['solution_tag']
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
    
    amgio.py_WriteMeshAndSol(mesh_name, solution_name, Ver, Cor, Tri, Tet, Edg, Hex, Qua, Pyr, Pri, Sol, SolTag, Markers, Dim)

# --- Write mesh and solution using amgio module
def write_mesh(mesh_name, mesh):
    
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
        
    if 'Triangles' in mesh:     Tri     = mesh['Triangles']
    if 'Tetrahedra' in mesh:    Tet     = mesh['Tetrahedra']
    if 'Edges' in mesh:         Edg     = mesh['Edges']
    if 'Corners' in mesh:       Cor     = mesh['Corners']
    if 'markers' in mesh:       Markers = mesh['markers']
    if 'dimension' in mesh:     Dim     = mesh['dimension']
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
    
    amgio.py_WriteMesh(mesh_name, Ver, Cor, Tri, Tet, Edg, Hex, Qua, Pyr, Pri, Markers, Dim)
    
# --- Write solution using amgio module
def write_sol(sol_name, sol):
    
    Dim = sol['dimension']
    Sol = sol['solution']

    if Dim == 3:
        NbrVer = int(len(sol['xyz']))
    else:
        NbrVer = int(len(sol['xy']))
    
    if 'xyz' in sol:
        Ver = sol['xyz']
        Ver = np.array(Ver).reshape(3*len(Ver)).tolist()
    elif 'xy' in sol:
        Ver = np.array(sol['xy'])
        z = np.zeros(len(sol['xy']))
        Ver = np.c_[Ver, z]
        Ver = np.array(Ver).reshape(3*len(Ver)).tolist()
    
    SolTag = sol['solution_tag']
        
    if len(Sol) > 1 :
        SolSiz = len(Sol[1])
        Sol = np.array(Sol).reshape(SolSiz*len(Sol)).tolist()
    else:
        sys.stderr.write("## ERROR write_sol : No solution.\n")
        sys.exit(1)
        
    amgio.py_WriteSol(sol_name, Ver, Sol, SolTag, NbrVer, Dim)



