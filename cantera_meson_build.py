#!/usr/bin/env python3

## \file cantera_meson_build.py
#  \brief create necessary meson.build file for linking cantera to
#         SU2.
#  \author C.Morales Ubal
#  \version 8.4.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
# License along with SU2. If not, see <http://www.gnu.org/licenses/>
import sys, os

absolute_path = sys.path[0]
relative_path = "subprojects/cantera"
meson_path = os.path.join(absolute_path, relative_path)

# Path to the meson.build under SU2/subprojects/cantera/
meson_build_path = os.path.join(meson_path, "meson.build")

# Create meson.build in the cantera subproject if it does not already exist
if not os.path.exists(meson_build_path):
   print(f"Writing meson.build in {meson_path}")
   meson_build_content = ["project('cantera', 'c', 'cpp', default_options: ['cpp_std=c++17'])\n",
                               "cc = meson.get_compiler('cpp')\n", 
                               "cantera_inc = include_directories('install/include')\n",
                               "cantera_lib = cc.find_library('cantera', dirs: [meson.current_source_dir() / 'install' / 'lib'], required: false, static: true)\n",
                               "cantera_dep = declare_dependency(include_directories: cantera_inc, dependencies: cantera_lib)\n",
                               "meson.override_dependency('cantera', cantera_dep)\n",
                               "cantera_ad_inc = include_directories('install_ad/include')\n",
                               "cantera_ad_lib = cc.find_library('cantera', dirs: [meson.current_source_dir() / 'install_ad' / 'lib'], required: false, static: true)\n",
                               "cantera_ad_dep = declare_dependency(include_directories: cantera_ad_inc, dependencies: cantera_ad_lib)\n",
                               "meson.override_dependency('cantera_ad', cantera_ad_dep)"]
   # Write the meson.build file
   with open(meson_build_path, 'w') as f:
      f.writelines(meson_build_content)
   print(f"meson.build in {meson_path} Created ")
