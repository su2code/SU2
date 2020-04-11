#!/usr/bin/env python

## \file amg.py
#  \brief python package for interfacing with the AMG Inria library
#  \author Victorien Menier, Brian Mungu\'ia
#  \version 7.0.1 "Blackbird"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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


def return_mesh_size(mesh):
    
    elt_key  = ['xy', 'xyz',  'Triangles', 'Edges', 'Tetrahedra']
    elt_name = {'xy':'vertices', 'xyz':'vertices',  'Triangles':'triangles', 'Edges':'edges', 'Tetrahedra':'tetrahedra'}
    
    tab_out = []
        
    for elt in elt_key:
        
        if elt not in mesh: continue
        
        nbe = len(mesh[elt])
        
        if nbe > 0:
            tab_out.append("%d %s" % (nbe, elt_name[elt]))
    
    return ', '.join(map(str, tab_out))
    
def prepro_back_mesh_2d(config_cfd, config_amg):
    #--- Read initial and background meshes
    sys.stdout.write("Reading initial and background mesh.\n")
    sys.stdout.flush()

    mesh_ini = read_mesh(config_cfd['MESH_FILENAME'])

    Tri_Ini = [mesh_ini['Triangles'][0][i] for i in range(3)]
    Ver_Ini = [mesh_ini['xyz'][i][:2] for i in Tri_Ini]

    V_Ini = [Ver_Ini[1][i] - Ver_Ini[0][i] for i in range(2)]
    W_Ini = [Ver_Ini[2][i] - Ver_Ini[0][i] for i in range(2)]

    N_Ini = V_Ini[0]*W_Ini[1]-V_Ini[1]*W_Ini[0]

    del mesh_ini

    sys.stdout.write("Reading background surface mesh.\n")
    sys.stdout.flush()
    mesh_bak = read_mesh(config_amg['adap_back'])

    Tri_Bak = [mesh_bak['Triangles'][0][i] for i in range(3)]
    Ver_Bak = [mesh_bak['xyz'][i][:2] for i in Tri_Bak]

    V_Bak = [Ver_Bak[1][i] - Ver_Bak[0][i] for i in range(2)]
    W_Bak = [Ver_Bak[2][i] - Ver_Bak[0][i] for i in range(2)]

    N_Bak = V_Bak[0]*W_Bak[1]-V_Bak[1]*W_Bak[0]

    #--- Check orientation
    sys.stdout.write("Checking orientation.\n")
    if (N_Ini*N_Bak < 0.0):
        #--- Flip all triangles if normals are opposite
        sys.stdout.write("Flipping triangles in background surface mesh.\n")
        sys.stdout.flush()

        NbrTri = len(mesh_bak['Triangles'])
        Tri = np.reshape(mesh_bak['Triangles'],(NbrTri, 4)).astype(int)
        for i in range(NbrTri):
            tmp = Tri[i][1]
            Tri[i][1] = Tri[i][2]
            Tri[i][2] = tmp

        mesh_bak['Triangles'] = Tri.tolist()
        Ver = mesh_bak['xyz']
        mesh_bak['xy'] = np.stack((Ver[:,0],Ver[:,1]), axis=1).tolist()
        del mesh_bak['xyz']
        write_mesh(config_amg['adap_back'], mesh_bak)

    del mesh_bak


def amg_call(config):
    
    cmd = ''
    cmd = "amg -in %s -sol %s -p 1 \
         -c %f -hgrad %.2f -hmin %le -hmax %le -out %s \
        -itp  %s  -nordg " \
        % (config['mesh_in'], config['sol_in'],  \
        config['size'],  config['hgrad'], config['hmin'], config['hmax'], \
        config['mesh_out'], config['sol_itp_in'])
        
    if config['adap_source'] != "":
        cmd += ' -source %s ' % config['adap_source']
    
    if config['adap_back'] != "":
        cmd += ' -back %s ' % config['adap_back']
    
    cmd += ' > %s' % config['amg_log']
    os.system(cmd)


def amg_call_met(config):
    
    cmd = ''
    cmd = "amg -in %s -met %s \
        -hgrad %.2f -hmin %le -hmax %le -out %s \
        -itp  %s  -nordg " \
        % (config['mesh_in'], config['sol_in'],  \
        config['hgrad'], config['hmin'], config['hmax'], \
        config['mesh_out'], config['sol_itp_in'])
            
    cmd += ' > %s' % config['amg_log']
    os.system(cmd)


def amg_call_python(mesh, config):
    
    remesh_options                = {}

    remesh_options['gradation']   = config['hgrad']
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
        mesh['sensor']         = mesh['sensor'].tolist()
        remesh_options['Lp']   = config['Lp']
        remesh_options['hmax'] = config['hmax']
        remesh_options['hmin'] = config['hmin']

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
    
    mesh['solution_tag'] = SolTag[:-1]
    
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

    Sol = []
    SolTag = []
    
    amgio.py_ReadSol(solution_name, Sol, SolTag, NbrVer)
        
    SolSiz = int(len(Sol)/NbrVer)
    Sol = np.array(Sol).reshape(NbrVer,SolSiz).tolist()
    
    sol = dict()

    sol['solution']     = Sol
    sol['solution_tag'] = SolTag[:-1]
    
    sol['id_solution_tag'] = dict()
    for i in range(len(SolTag)):
        sol['id_solution_tag'][SolTag[i]] = i
    
    return sol

# --- Merge 2 solutions (e.g. primal and dual)
def merge_sol(mesh0, mesh1):
    mesh0['solution'] = np.hstack((mesh0['solution'],mesh1['solution'])).tolist()
    mesh0['solution_tag'] = np.hstack((mesh0['solution_tag'],mesh1['solution_tag'])).tolist()

# --- Split adjoint solution
def split_adj_sol(mesh):
    nsol = len(mesh['solution_tag'])

    for i in range(nsol):
        if "Adjoint" in mesh['solution_tag'][i]:
            iAdj = i
            break

    adj_sol = dict()

    adj_sol['solution'] = mesh['solution'][:][iAdj:]
    adj_sol['solution_tag'] = mesh['solution_tag'][iAdj:]

    if 'xyz' in mesh:
        adj_sol['xyz'] = mesh['xyz']
    elif 'xy' in mesh:
        adj_sol['xy'] = mesh['xy']

    adj_sol['dimension'] = mesh['dimension']

    np.delete(np.array(mesh['solution']), np.s_[iAdj:nsol], axis=1).tolist()
    np.delete(np.array(mesh['solution_tag']), np.s_[iAdj:nsol]).tolist()

    return adj_sol
    
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


def create_sensor(solution, sensor):
    
    Dim = solution['dimension']
    Sol = np.array(solution['solution'])
    
    if sensor == "MACH":
        
        iMach = solution['id_solution_tag']['Mach']
        sensor = Sol[:,iMach]
        sensor = np.array(sensor).reshape((len(sensor),1))
        sensor_header = ["Mach"]
        
    elif sensor == "PRES":
        
        iPres = solution['id_solution_tag']['Pressure']
        sensor = Sol[:,iPres]
        sensor = np.array(sensor).reshape((len(sensor),1))        
        sensor_header = ["Pres"]
        
    elif sensor == "MACH_PRES":

        iPres  = solution['id_solution_tag']['Pressure']
        iMach  = solution['id_solution_tag']['Mach']
        mach   = np.array(Sol[:,iMach])
        pres   = np.array(Sol[:,iPres])
        sensor = np.stack((mach, pres), axis=1)    
        sensor_header = ["Mach", "Pres"]

    elif sensor == "GOAL":

        if Dim == 2:
            sensor = Sol[:,-3:]
            sensor = np.array(sensor).reshape((len(sensor),3))
        elif Dim == 3:
            sensor = Sol[:,-6:]
            sensor = np.array(sensor).reshape((len(sensor),6))
        sensor_header = ["Goal"]
                
    else :
        sys.stderr.write("## ERROR : Unknown sensor.\n")
        sys.exit(1)
    
    sensor_wrap = dict()
    
    sensor_wrap['solution_tag'] = sensor_header
    sensor_wrap['xyz'] = solution['xyz']
    
    sensor_wrap['dimension']    = solution['dimension']
    sensor_wrap['solution']     = sensor
    
    return sensor_wrap


