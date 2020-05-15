#!/usr/bin/env python

## \file FSIOpt_interface.py
#  \python package interfacing with the FSI optimization suite
#  \author Rocco Bombardieri
#  \version 7.0.2 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil
import subprocess

def SaveSplineMatrix(config):
    """
    Spline matrix is computed at the beginning and saved in the main folder to be used at every iteration
    """
    # compose string command
    command = 'mpirun -n ' + str(config['NUMBER_PART']) + ' pyBeamFSI_MLSGen.py -f ' + config['CONFIG_PRIMAL']
    print (command)
    # Compose local output file
    Output_file = config['FOLDER'] + '/Output_Spline.out'

    # Launching shell command
    run_command(command, 'Splining', True,  Output_file)


def run_command(Command, Tool, Output, Output_file = '' ):
    """ runs os command with subprocess
        checks for errors from command
    """

    sys.stdout.flush()
    file = open(Output_file, "w")
    with subprocess.Popen(Command, shell=True, universal_newlines=True,
                            stdout=subprocess.PIPE) as proc:
        if Output == True:
           while True:
               string = proc.stdout.read(1)
               if string:
                   file.write(string)
                   # logfile.flush()
               else:
                   break
    # Waits for code to terminate
    return_code = proc.wait()
    message = string

    if return_code < 0:
        message = Tool + " process was terminated by signal '%s'\n%s" % (-return_code, message)
        raise SystemExit(message)
    elif return_code > 0:
        message = "Path = %s\nCommand = %s\n process returned error '%s'\n%s" % (
        os.path.abspath(','), Tool, return_code, message)
    else:
        sys.stdout.write(message)

    return return_code


def update_mesh():
    '''
            Executes in:
             ./DEFORM
    '''  
    return 


def Primal():
    '''
            Executes in:
             ./Primal
    '''      
    return

def Adjoint():
    '''
            Executes in:
             ./Adjoint
    '''      
    return    

def Geometry():
    '''
            Executes in:
             ./GEO
    '''    
    return


