#!/usr/bin/env python

## \file back_mesh.py
#  \brief Python script for generating a meshb background mesh from an SU2 mesh.
#  \author Brian Mungu\'ia
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

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6

from __future__ import print_function

from optparse import OptionParser
import numpy as np
from pathlib import Path
import SU2.amginria as su2amg
import pyamg

parser = OptionParser()
parser.add_option("-f", "--file", dest="file",
                  help="read mesh from FILE", metavar="FILE")
parser.add_option("-o", "--output", dest="outfile",
                  help="write new mesh to OUTFILE", metavar="OUTFILE", default="out")
parser.add_option("--hgrad", dest="hgrad",
                  help="gradation", metavar="HGRAD", default=2.5)
parser.add_option("--hmax", dest="hmax",
                  help="max cell size", metavar="HMAX", default=100)
parser.add_option("--hmin", dest="hmin",
                  help="min cell size", metavar="HMIN", default=0.001)
parser.add_option("--iters", dest="iters",
                  help="iterations of coarsening", metavar="ITERS", default=3)

(options, args)=parser.parse_args()

# Process options
file = str(options.file)
outfile = str(options.outfile)
hgrad = float(options.hgrad)
hmax = float(options.hmax)
hmin = float(options.hmin)
iters = int(options.iters)

mesh = su2amg.read_mesh(file)

remesh_options = {}

# Set amg options from input
remesh_options['gradation'] = hgrad
remesh_options['hmax'] = hmax
remesh_options['hmin'] = hmin
remesh_options['logfile'] = 'amg.coarsen.out'
remesh_options['options'] = f'-hmsh-line -noline -nordg'

Dim = mesh['dimension']

if Dim == 2 :
    Ver = mesh['xyz']
    mesh['xy'] = np.stack((Ver[:,0],Ver[:,1]), axis=1)
    
    del mesh['xyz']

# Store initial mesh
if 'xy' in mesh:  mesh['xy']  = mesh['xy'].tolist()
if 'xyz' in mesh: mesh['xyz'] = mesh['xyz'].tolist()

if 'Corners' in mesh:    mesh['Corners']    = mesh['Corners'].tolist() 
if 'Edges' in mesh:      mesh['Edges']      = mesh['Edges'].tolist() 
if 'Triangles' in mesh:  mesh['Triangles']  = mesh['Triangles'].tolist()
if 'Tetrahedra' in mesh: mesh['Tetrahedra'] = mesh['Tetrahedra'].tolist()   

# Remesh
for i in range(iters):
    mesh = pyamg.adapt_mesh(mesh, remesh_options)
    if i == (iters-1): break
    su2amg.write_mesh(f"{outfile}.{i}.meshb", mesh)

# Output final mesh
su2amg.write_mesh(f"{outfile}.meshb", mesh)

# Delete intermediate files
for file in ['back.meshb', 'meshp3_smoo.meshb']:
    Path(file).unlink()
