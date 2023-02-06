#!/usr/bin/env python

## \file which.py
#  \brief looks for where a program is
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


import os

def which(program):
    """ which(program_name)
        finds the location of the program_name if it is on PATH
        returns None if program cannot be found
        does not test for .exe extension on windows

        original source:
        http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            for ext in ['','.exe','.bat']:
                exe_file = os.path.join(path, (program+ext))
                if is_exe(exe_file):
                    return exe_file

    return None

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)



