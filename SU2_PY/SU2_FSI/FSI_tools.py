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

def readConfig(ConfigFileName, voice):
    """
    This function scans an input file looking for a specific voice
    """    
    input_file = open(ConfigFileName)
    while 1:
        line = input_file.readline()
        if not line:
            break
        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
            continue
        # split across equal sign
        line = line.split("=", 1)
        this_param = line[0].strip()
        this_value = line[1].strip()
        if this_param == voice:
            break
    
    if not this_value:
        raise SystemExit(voice + ' is not present in ' + 'ConfigFileName')
    
    return this_value
        

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


def UpdateConfig(ConfigFileName, param, value):   
    """
    This function updates the input param of the given ConfigFileName with the new value
    """ 
    configfile2 = open(ConfigFileName + '_temp',"w")
    with open(ConfigFileName, 'r') as configfile:
      while 1:          
	line = configfile.readline()
        string = line
	if not line:
	  break

        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):  
         configfile2.write(string)
        else: 
         # split across equal sign
         line = line.split("=",1)
         this_param = line[0].strip()
         this_value = line[1].strip()

         #float values
         if this_param == param:
            if this_param == "DV_VALUE":            
                    dv_string = ('%s' % ', '.join(map(str, value)))
                    stringalt = 'DV_VALUE = '+ dv_string + '   \r\n'
                    configfile2.write(stringalt)    
                    
            else:
                    stringalt = this_param + ' = ' + value + '   \r\n'
                    configfile2.write(stringalt) 
                
                    
    configfile.close()    
    configfile2.close()
    # the file is now replaced
    os.remove(ConfigFileName)
    os.rename(ConfigFileName + '_temp', ConfigFileName)                    


def DeformMesh(deform_folder, ConfigFileName):
    '''
            Executes in:
             ./DEFORM
    '''  
    # go to ./DEFORM folder
    os.chdir(deform_folder)
    
    command = 'SU2_DEF ' + ConfigFileName
    Output_file = deform_folder + '/Output_SU2_DEF.out'
    run_command(command, 'SU2_DEF', True, Output_file)
    
    # go back to project folder (two levels up)
    os.chdir(deform_folder + '/../..')

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

def Geometry(geo_folder, ConfigFileName):
    '''
            Executes in:
             ./GEO
    '''    
        # go to ./DEFORM folder
    os.chdir(geo_folder)
    
    command = 'SU2_GEO ' + ConfigFileName
    Output_file = geo_folder + '/Output_SU2_GEO.out'
    run_command(command, 'SU2_DEF', True, Output_file)
    
    # go back to project folder (two levels up)
    os.chdir(geo_folder + '/../..')

    
def ReadGeoConstraints(self, geo_folder,ConsList, sign ):
    """ 
    Fuction that returns the numpy list of geometrical constraints for the given sign
    """
        
    # Open and read output file
    input_file = open(geo_folder + '/' + 'of_func.dat')
    # going directly to line 3
    line = input_file.readline();line = input_file.readline();line = input_file.readline()
     
    # list of constraints given in output
    constr = line.split(',')
    for i in range(len(constr)): constr[i] = constr[i].strip(' "\n')
    # List of values given in output
    # going to line 5
    line = input_file.readline();line = input_file.readline();
    # List of constraint values (float)
    value = line.split(',')
    value = [float(i) for i in value]  
     
    # Initialization of c_eq list
    c_eq_list = []
        
    # quick loop over the constraints
    for i in range(len(ConsList)):
            # Looking for the given constraints
            if ConsList[i][1] == sign:
                # loop over the output constraints 
                for j in range(len(constr)): 
                    if constr[j] == ConsList[i][0]:
                        # value of the constraint
                        a = value[j] - ConsList[i][2]
                        # if there is scale factor it is multiplied
                        if ConsList[i][3] != None:
                            a = a * ConsList[i][3]
                        # adding to list     
                        c_eq_list.append(a)
            
             
    # returning list as numpy list
    return np.append(c_eq_list)  