#!/usr/bin/env python

## \file Compute_multipoint.py
#  \brief Python script for performing a multipoint design.
#  \author Indiana Stokes
#  \version 6.2.0 "Falcon"
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
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

# imports
import numpy as np
from optparse import OptionParser
import os, sys, shutil, copy
import SU2

# Command Line Options
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-n", "--partitions", dest="partitions", default=2,
                  help="number of PARTITIONS", metavar="PARTITIONS")

(options, args)=parser.parse_args()
options.partitions = int( options.partitions )

# load config, start state
config = SU2.io.Config(options.filename)
state  = SU2.io.State()

# prepare config
config.NUMBER_PART = options.partitions

# find solution files if they exist
state.find_files(config)

# run su2
multipoint_drag = SU2.eval.func('MULTIPOINT_DRAG',config,state)
grad_multipoint_drag= SU2.eval.grad('MULTIPOINT_DRAG','CONTINUOUS_ADJOINT',config,state)

print('MULTIPOINT_DRAG     =', multipoint_drag)
print('GRADIENT MULTIPOINT_DRAG =', grad_multipoint_drag)
