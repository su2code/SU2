#!/usr/bin/env python

## \file FSI_tools.py
#  \python various tools for the FSI optimization suite
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
from math import pow, factorial
import scipy.io

def SaveSplineMatrix(config):
    """
    Spline matrix is computed at the beginning and saved in the main folder to be used at every iteration
    """
    # compose string command
    command = 'mpirun -n ' + str(config['NUMBER_PART']) + ' pyBeamFSI_MLSGen.py -f ' + config['CONFIG_PRIMAL']
    #print (command)
    # Compose local output file
    Output_file = config['FOLDER'] + '/Output_Spline.out'

    # Launching shell command
    run_command(command, 'Splining', True,  Output_file)

def readConfig(ConfigFileName, voice, BreakCode = True):
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
    
    if BreakCode == True:
       if not this_value:
           raise SystemExit(voice + ' is not present in ' + 'ConfigFileName')
    else:
        if not this_value:
           this_value = 'NO'
        
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
    
    # going to ./Primal folder
    os.chdir(primal_folder)    
    command = 'mpirun -n ' + str(config['NUMBER_PART']) + ' pyBeamFSI_opt.py -f ' + config['CONFIG_PRIMAL']
    #print (command)
    # Compose local output file
    Output_file =  'Output_primal.out'

    # Launching shell command
    run_command(command, 'Primal', True,  Output_file)
    
    # go back to project folder (3 levels up)
    os.chdir( '../../..')    
    


def FSIAdjoint(adj_folder, config):
    '''
            Executes in:
             ./Adjoint
    '''      

    # going to ./Primal folder
    os.chdir(adj_folder)    
    command = 'mpirun -n ' + str(config['NUMBER_PART']) + ' pyBeamFSI_AD_opt.py -f ' + config['CONFIG_ADJOINT']
    #print (command)
    # Compose local output file
    Output_file =  'Output_adjoint.out'   
    
    # Launching shell command
    run_command(command, 'Adjoint', True,  Output_file)
    
    # go back to project folder (3 levels up)
    os.chdir( '../../..')        
    
def Geometry(geo_folder, config):
    '''
            Executes in:
             ./GEO
    '''    
    # going to ./GEO folder
    os.chdir(geo_folder)
    
    ConfigFileName = config['CONFIG_GEO']
    
    command = 'mpirun -n ' + str(config['NUMBER_PART']) + ' SU2_GEO ' + ConfigFileName
    Output_file =  'Output_SU2_GEO.out'
    run_command(command, 'SU2_DEF', True, Output_file)
    
    # go back to project folder (two levels up)
    os.chdir( '../../..')

    
def ReadGeoConstraints( geo_folder,ConsList, sign, it ):
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
        
        
    if it ==0:
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

def ReadGeoConstraintGradients( geo_folder,ConsList,n_dv, sign ):
    """ 
    Fuction that returns the numpy list of the gradient of geometric constraints 
    """
        
    # Open and read output file
    input_file = open(geo_folder + '/' + 'geo_gradient.dat')
    # going directly to line 3
    line = input_file.readline();line = input_file.readline();line = input_file.readline()
     
    # list of constraints given in output
    constr = line.split(',')
    for i in range(len(constr)): constr[i] = constr[i].strip(' "\n')  
    
    # going to line 4 
    line = input_file.readline(); 
    
    # Gradients for every control point given in output as a matrix   
    value_arr = []
    for i in range(n_dv):
        line = input_file.readline();
        # List of constraint values (float)
        value = line.split(',')
        value = [float(i) for i in value]  
        value_arr.append(value)
    value_matr = np.array(value_arr)
    
    # Initialization of gradient matrix list
    gradient_array = []  
    
    # quick loop over the constraints
    for i in range(len(ConsList)):
            # Looking for the given constraints
            if ConsList[i][1] == sign:
                # loop over the output constraints 
                for j in range(len(constr)): 
                    if constr[j] == ConsList[i][0]:   
                        a = value_matr[:,j]
                        gradient_array.append(a)
    
    gradient = np.array(gradient_array)
    
    return gradient
    
def PullingPrimalAdjointFiles(configOpt, folder, configFSI, pyBeamMesh, pyBeamProp):
       # pulling primal files
       command = []
       # 1
       config = configOpt['FOLDER'] + '/' + configFSI['SU2_CONFIG']
       command.append('cp ' + config + ' ' + folder + '/')
       # 2
       config = configOpt['FOLDER'] + '/' + configFSI['PYBEAM_CONFIG']
       command.append('cp ' + config + ' ' + folder + '/')
       # 3
       config = configOpt['FOLDER'] + '/' + configFSI['MLS_CONFIG_FILE_NAME']
       command.append('cp ' + config + ' ' + folder + '/')
       # 4
       spline = configOpt['FOLDER'] + '/' + 'Spline.npy'
       command.append('cp ' + spline + ' ' + folder + '/')      
       for i in range(len(command)):
          run_command(command[i], 'Pulling primal config file ' + str(i) , False) 
          
       # pulling files for pyBeam (mesh and properties)
       command = []
       # 1
       config = configOpt['FOLDER'] + '/' + pyBeamMesh
       command.append('cp ' + config + ' ' + folder + '/')
       # 2
       config = configOpt['FOLDER'] + '/' + pyBeamProp
       command.append('cp ' + config + ' ' + folder + '/')       
       for i in range(len(command)):
          run_command(command[i], 'Pulling primal pyBeam file ' + str(i) , False)        
          
def ReadPointInversion(configDef,MeshFile):
    # First it is necessary to look for the DV_MARKER 
    DV_MARKER = readConfig(configDef, 'DV_MARKER')
    DV_MARKER = DV_MARKER.strip('( )')
    # Open mesh file
    PointInv_arr = []
    input_file = open(MeshFile)
    while 1:
        line = input_file.readline()
        if not line:
            break
        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (line[0] == '%'):
            continue
        if line.startswith(DV_MARKER):
            line = line.split()
            a = [int(line[1]),float(line[2]),float(line[3]),float(line[4])]
            PointInv_arr.append(a)
    PointInv = np.array(PointInv_arr) 
    # Sorting according to first column index
    # In case nodes on the boundary are not numbered from 0 to the maximum a argsort is required for PointInv[:,0] 
    PointInv = PointInv[PointInv[:,0].argsort()]
    
    return PointInv

def barnstein(l, mi):
    """ 
    Barnstein polynomials. L is the order of the polinomial. Mi is the parametric point (0,1) around which is evaluated.
    Output is the vector of polynomial coefficients.
    """
    N = np.zeros((l+1, 1))    
    for i in range(l+1):
       N[i] = ( factorial(l)/factorial(i)/factorial(l-i) )* pow((1-mi),(l-i))*pow(mi,i)
       
    return N 

def readDVParam(config_DEF):
    """ 
    This fuction returns the matrix of the i,j,k indices of the DV params. 
    layout: DV_PARAM= ( WING_FFD,    0.0,    0.0,  0.0,   0.0,    0.0,    1.0) ; 
    FFD_CONTROL_POINT ( FFD_BoxTag, i_Ind, j_Ind, k_Ind, x_Disp, y_Disp, z_Disp )    
    For now only one FFD box is considered: oner FFD box name
    For now only z displacements are taken into consideration: ... x_Disp = 0, y_Disp = 0, z_Disp = 1 )
    This is the order in which x_dv are provided by the optimizator and whic hare given to SU2_DEF for deformation (DV_value)
    These indexes are used for the chain rule of surface grid sensitivities
    """    
    
    
    dv_param_str = readConfig(config_DEF, 'DV_PARAM')
    
    dv_param_str = dv_param_str.split(';')
    FFD_indexes = np.empty([len(dv_param_str),3], dtype=int)
    for i in range(len(dv_param_str)):
       dv_param_str[i] = dv_param_str[i].strip('( )')
       a = dv_param_str[i].split(',')
       FFD_indexes[i,:] = [int(float(a[1].strip())),int(float(a[2].strip())),int(float(a[3].strip()))]    
       
    return   FFD_indexes

def readBoundarySensitivities(SensFile):
    # Open mesh file
    GridSensitivities_arr = []
    input_file = open(SensFile)
    while 1:
        line = input_file.readline()
        if not line:
            break
        # remove line returns
        line = line.strip('\r\n')
        line = line.split()
        a =   [float(line[4]), float(line[5]), float(line[6]) ]
        GridSensitivities_arr.append(a)
        
    GridSensitivities = np.array(GridSensitivities_arr)  
    
    return GridSensitivities
          
def ChainRule(adj_folder,FFD_indexes, PointInv,ffd_degree):    
    """ 
    Chain rule to evaluate sensitivity of obj_f with respect to the FFD box control points
    """    
    
    # Reading boundary nodes sensitivities
    SensFile = adj_folder + '/' + 'Boundary_Nodes_Sensitivity.dat'
    GridSensitivities = readBoundarySensitivities(SensFile)
    
    # number of control points with respect to which I evaluate sensitivity
    Nalpha = np.shape(FFD_indexes)[0]
    
    # Number of boundary grid points
    Nmesh = np.shape(PointInv)[0]
    # Initialization of sensitivity matrix this has to be multiplied for each DoF/DeltaX_alpha DoF/DeltaY_alpha DoF/DeltaZ_alpha 
    # where DoF/DeltaX_alpha(i) is the total derivative of the obj_f with respect of the i-th boundary node X displacement
    chain = np.zeros((Nmesh, Nalpha)) 
    
    # Loop over the mesh points
    for i in range(Nmesh):
       Ni =  barnstein(ffd_degree[0], PointInv[i,1] )
       Nj =  barnstein(ffd_degree[1], PointInv[i,2] )
       Nk =  barnstein(ffd_degree[2] , PointInv[i,3] )        
       for j in range(Nalpha):
           chain[i,j] =  Ni[ int(FFD_indexes[j, 0]) ] * Nj[ int(FFD_indexes[j, 1]) ] *Nk[ int(FFD_indexes[j, 2]) ]
    
    
    # chain rule for the sensitivity (for now only Z is allowed)
    obj_df = np.zeros((Nalpha)) 
    
    for i in range(Nalpha):
       obj_df[i] = np.dot(GridSensitivities[:, 2], chain[:,i])
        
    # writing output file with the chain rule
    Output = adj_folder + '/' + 'Control_Points_Sensitivity.dat'
    np.savetxt(Output, obj_df, delimiter =', ')
    
    # DEBUGGING
    '''
    print('adj_folder = '+adj_folder)
    scipy.io.savemat( './PointInvpy.mat', mdict={'PointInvpy': PointInv}) 
    scipy.io.savemat( './chainpy.mat', mdict={'chainpy': chain})    
    scipy.io.savemat( './GridSensitivitiespy.mat', mdict={'GridSensitivitiespy': GridSensitivities}) 
    scipy.io.savemat( './FFD_indexespy.mat', mdict={'FFD_indexespy': FFD_indexes})
    scipy.io.savemat( './obj_dfpy.mat', mdict={'obj_dfpy': obj_df})
    '''
    return obj_df  
    
    
def WriteSolution(folder,x,it):    
    
    logfile = 'Solution.dat'
    if it ==0:
       log = open(folder + '/' + logfile,"w") 
    else:   
       log = open(folder + '/' + logfile,"a") 
       
    log.write('%3s \t' % str(it) ) 
    for i in range(len(x)):
          log.write('%10s \t' % str(x[i]) )
    log.write("\n")      
    log.close() 
    
    return

def SharpEdge(adj_folder,configAdj):
    """ 
    Sharp edge sensitivity option. Removes surface sesntivity in areas with a sharp edge (LE,TE and wing tip).
    It is stated in SU2_CFD adjoint solver config file. Reads cvs surface_adjoint and check where sensitivity is put to 0.
    Then it communicates it to the 'Boundary_Nodes_Sensitivity.dat' which is used for chain rule.
    NB: it requires surface_adjoint.cvs as output of the adjoint solver:
    % Output file surface adjoint coefficient (w/o extension)
    SURFACE_ADJ_FILENAME= surface_adjoint   
    % Files to output 
    % Possible formats : (TECPLOT, TECPLOT_BINARY, SURFACE_TECPLOT,
    %  SURFACE_TECPLOT_BINARY, CSV, SURFACE_CSV, PARAVIEW, PARAVIEW_BINARY, SURFACE_PARAVIEW, 
    %  SURFACE_PARAVIEW_BINARY, MESH, RESTART_BINARY, RESTART_ASCII, CGNS, STL)
    % default : (RESTART, PARAVIEW, SURFACE_PARAVIEW)
    OUTPUT_FILES= (RESTART, PARAVIEW, SURFACE_PARAVIEW,SURFACE_CSV)    
    """      

    # solution file to read
    surface_adj = readConfig(configAdj, 'SURFACE_ADJ_FILENAME')
    ftr = adj_folder + '/' + surface_adj + '.csv'

    CSV = []
    
    # reading cvs file
    try:
       with open(ftr, 'r') as file:
          reader = csv.reader(file)
          for row in reader:

             if row[0] == "PointID":
                continue
             else: 
                ID = int(row[0]);  Sens_surf = float(row[-1]); Sens_z = float(row[-2]); Sens_y = float(row[-3]); Sens_x = float(row[-4]); 
                CSV.append([ID,Sens_x,Sens_y,Sens_z,Sens_surf])
             
    except IOError:
       print("Surface adjoint cvs not available")
    
    
    # Update the  Boundary_Nodes_Sensitivity putting to 0 sharp edge sensitivities
    ConfigFileName = adj_folder + '/' + 'Boundary_Nodes_Sensitivity.dat'
    configfile2 = open(ConfigFileName + '_temp',"w")
    i = 0
    with open(ConfigFileName, 'r') as configfile:
      while 1:          
        line = configfile.readline()
        string = line
        if not line:
           break
        if float(CSV[i][-1]) ==0 and float(CSV[i][-2]) ==0 and float(CSV[i][-3]) ==0 and float(CSV[i][-4]) ==0:
           # remove line returns
           line = line.strip('\r\n')  
           line = line.split()
           string_alt = line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + str(0) + '\t' + str(0) + '\t' + str(0) + '\n'
           configfile2.write(string_alt)
        else:
           configfile2.write(string)
        i = i+1

    configfile.close()    
    configfile2.close()        
    # the file is now replaced
    os.rename(ConfigFileName, 'Boundary_Nodes_Sensitivity_original.dat' )
    os.rename(ConfigFileName + '_temp', ConfigFileName) 
    
    return