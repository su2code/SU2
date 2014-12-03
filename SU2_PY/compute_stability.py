#!/usr/bin/env python

## \file compute_stability.py
#  \brief Python script for performing the shape optimization.
#  \author T. Lukaczyk, F. Palacios
#  \version 3.2.5 "eagle"
#
# Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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
parser.add_option("-p", "--oldpartitions", dest="oldpartitions", default="oldpartitions",
                  help="old number of PARTITIONS (use -n instead)", metavar="OLDPARTITIONS")
parser.add_option("-i", "--iterations", dest="iterations", default=99999,
                  help="number of ITERATIONS", metavar="ITERATIONS")

(options, args)=parser.parse_args()
options.partitions = int( options.partitions )
options.iterations = int( options.iterations )

if options.oldpartitions != "oldpartitions":
  print ("\n IMPORTANT: -p is no longer available in SU2 v3.2.4, use -n flag instead \n")
  sys.exit()

# load config, start state
config = SU2.io.Config(options.filename)
state  = SU2.io.State()

# prepare config
config.NUMBER_PART = options.partitions
config.EXT_ITER    = options.iterations

# find solution files if they exist
state.find_files(config)

# run su2
drag_alpha = SU2.eval.func('D_DRAG_D_ALPHA',config,state)
moment_y_alpha= SU2.eval.func('D_MOMENT_Z_D_ALPHA',config,state)

grad_moment_y_alpha= SU2.eval.grad('D_MOMENT_Z_D_ALPHA','ADJOINT',config,state)

print 'D_DRAG_D_ALPHA     =' , drag_alpha
print 'D_MOMENT_Y_D_ALPHA =' , moment_y_alpha

print 'DD_MOMENT_Y_D_ALPHA_D_X =', grad_moment_y_alpha
