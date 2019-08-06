#!/usr/bin/env python3

## \file meson.py
#  \brief An extended meson script for setting up the environment and running meson
#  \author T. Albring
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
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


import sys, os, subprocess, shutil
from init import init_submodules

def build_ninja():

  ninjapath = sys.path[0] + '/externals/ninja'
    
  try:
    subprocess.run([sys.path[0] + '/ninja', '--version'], stdout=subprocess.PIPE)
  except OSError:   
    print("ninja executable not found. Building ...")
    subprocess.run(['./configure.py', '--bootstrap'], cwd=ninjapath)
    shutil.copy(ninjapath+'/ninja', '.')

if __name__ == '__main__':
  if sys.version_info[0] < 3:
    raise Exception("Script must be run using Python 3")
   
  # Set up the build environment, i.e. clone or download all submodules
  init_submodules('auto')

  # Build ninja if it cannot be found
  build_ninja()

  # Add paths for meson and ninja to environment
  os.environ["NINJA"] = sys.path[0] + "/ninja"
  if os.path.exists(sys.path[0] + '/externals/meson/mesonbuild'):
    sys.path.insert(0, str(sys.path[0] + '/externals/meson'))

  from mesonbuild import mesonmain

  sys.exit(mesonmain.main())
