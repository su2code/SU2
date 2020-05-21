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
import numpy as np

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
    if Output == True:
       file = open(Output_file, "w")
    with subprocess.Popen(Command, shell=True, universal_newlines=True,
                            stdout=subprocess.PIPE) as proc:
           while True:
               string = proc.stdout.read(1)
               if string and Output_file != '':
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
           else:
              configfile2.write(string)
                
                    
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
    # going to ./DEFORM folder
    os.chdir(deform_folder)
    
    command = 'SU2_DEF ' + ConfigFileName
    Output_file =  'Output_SU2_DEF.out'
    run_command(command, 'SU2_DEF', True, Output_file)
    
    # go back to project folder (3 levels up)
    os.chdir( '../../..')

    return 


def FSIPrimal(primal_folder, config):
    '''
            Executes in:
             ./Primal
    '''   
    
    # going to ./GEO folder
    os.chdir(primal_folder)    
    command = 'mpirun -n ' + str(config['NUMBER_PART']) + ' pyBeamFSI_opt.py -f ' + config['CONFIG_PRIMAL']
    print (command)
    # Compose local output file
    Output_file =  'Output_primal.out'

    # Launching shell command
    run_command(command, 'Primal', True,  Output_file)
    
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
    # going to ./GEO folder
    os.chdir(geo_folder)
    
    command = 'SU2_GEO ' + ConfigFileName
    Output_file =  'Output_SU2_GEO.out'
    run_command(command, 'SU2_DEF', True, Output_file)
    
    # go back to project folder (two levels up)
    os.chdir( '../../..')

    
def ReadGeoConstraints( geo_folder,ConsList, sign, iter ):
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
        
    # printing options
    constraint_list = []
    value_list = []
    target_list = []
        
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
                        #appending values for printing options
                        constraint_list.append(ConsList[i][0])
                        value_list.append(value[j])
                        target_list.append(ConsList[i][2])
                        
    # log file printing
    if sign == '=':
        logfile = 'Equality_constr.dat'
    elif sign =='<':
        logfile = 'Inequality_constr_minus.dat'
    else:
        logfile = 'Inequality_constr_plus.dat'
        
        
    if iter ==0:
       log = open(geo_folder + '/../../' + logfile,"w") 
       for i in range(len(constraint_list)):
           log.write('%25s \t' % str(constraint_list[i]) )
       log.write("\n")
       log.write("\n")
       for i in range(len(target_list)):
           log.write('%25s \t' % str(target_list[i]) )
              
       log.write("\n")     
    else:
       log = open(geo_folder + '/../../' + logfile,"a") 
    
    for i in range(len(target_list)):
       log.write('%25s \t' % str(value_list[i]) )
    log.write("\n")        
    
    log.close()         
    # returning list as numpy list

    return np.array(c_eq_list)  