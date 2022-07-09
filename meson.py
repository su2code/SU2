#!/usr/bin/env python3

## \file meson.py
#  \brief An extended meson script for setting up the environment and running meson
#  \author T. Albring
#  \version 7.3.1 "Blackbird"
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

import sys, os, subprocess, shutil, urllib.request, zipfile
sys.path.append(sys.path[0] + os.path.sep + 'meson_scripts')
from init import init_submodules
from init import remove_file

def build_ninja():

  # If we are on windows, we don't need to compile ninja, we just download the executable
  if os.name == 'nt':
    ninja_exe_url = 'https://github.com/ninja-build/ninja/releases/download/v1.9.0/ninja-win.zip'
    
    # Try to execute ninja, if it fails, download .exe from github
    try:
      subprocess.run([sys.path[0] + os.path.sep + 'ninja.exe', '--version'], stdout=subprocess.PIPE)
    except OSError:
      print ('Downloading ninja ... ')
      try:
        urllib.request.urlretrieve (ninja_exe_url,'ninja-win.zip')
      except:
        print(e)
        print('Download of ninja executable failed.')
        print('Get archive at ' + ninja_exe_url)
        print('extract ninja.exe in the source code root folder.')
        print('Run meson.py again.')
        sys.exit(1)

      zipf = zipfile.ZipFile(sys.path[0] + os.path.sep + 'ninja-win.zip')
      zipf.extractall(sys.path[0])
      remove_file(sys.path[0] + os.path.sep + 'ninja-win.zip')
  else:
    ninjapath = sys.path[0] + os.path.sep + 'externals' + os.path.sep + 'ninja'
    try:
      subprocess.run([sys.path[0] + os.path.sep + 'ninja', '--version'], stdout=subprocess.PIPE)
    except OSError:
      print("ninja executable not found. Building ...")
      subprocess.run(['python3', 'configure.py', '--bootstrap'], cwd=ninjapath)
      shutil.copy(ninjapath+ os.path.sep + 'ninja', '.')

if __name__ == '__main__':
  if sys.version_info[0] < 3:
    raise Exception("Script must be run using Python 3")
   
  # Set up the build environment, i.e. clone or download all submodules
  init_submodules('auto')

  # Build ninja if it cannot be found
  build_ninja()

  # Add paths for meson and ninja to environment
  os.environ["NINJA"] = sys.path[0] + os.path.sep + "ninja"
  if os.name == 'nt':
    os.environ["NINJA"] = os.environ["NINJA"] + '.exe'
  if os.path.exists(sys.path[0] + os.path.sep + 'externals' + os.path.sep + 'meson' + os.path.sep + 'mesonbuild'):
    sys.path.insert(0, str(sys.path[0] + os.path.sep + 'externals' +os.path.sep + 'meson'))

  from mesonbuild import mesonmain

  sys.exit(mesonmain.main())
