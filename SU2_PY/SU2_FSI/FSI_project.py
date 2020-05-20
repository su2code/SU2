#!/usr/bin/env python

## \file FSI_project.py
#  \brief FSI Shape Optimization project orchestrator.
#  \author Rocco Bombardieri based on work of  T. Lukaczyk, F. Palacios
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

import numpy as np
from math import pow, factorial
from SU2_FSI.FSI_config import FSIConfig as FSIConfig
from SU2_FSI import FSI_design
from SU2_FSI.FSI_tools import run_command
# -------------------------------------------------------------------
#  Project Class
# -------------------------------------------------------------------

class Project(object):
    
    """ 
    Project interface class that handles FSI shape optimization
    """

    def __init__(self, config ):
        """
        Class constructor. Declare some variables and do some screen outputs.
        """

        self.config = config  # FSI optimization config object
        
        folder = self.config['FOLDER']  # root folder where optimization is done
        folder = folder.rstrip('/') + '/'
        self.folder  = folder  
        self.deform_folder = ''
        self.geo_folder = ''
        
        self.design_toll = 10**8  # allowable difference into design variable vector to consider the same design
        
        # config objects for primal and adjoint simulations with structural and fluid config files and options
        self.configFSIPrimal = None
        self.configFSIAdjoint = None
        
        # Design container
        self.design = []

        self.design_iter = -1  # optimization iter (design number) [initialization]
        self.magnord_design = 3 # Expected order of magnitude of design number
        
        self.n_dv = 0 # number of design variables
        
        # clean previous designs
        self.clean_previous_designs()

        # Creating design folder
        self.create_design_folder()
        
        # project memorizes config of both adjoin and primal solvers (top level config)
        self.setup_configs()


    def setup_configs(self):
        '''
        Reads config files for primal and Adjoint simulations
        '''
        # Read primal config
        self.configFSIPrimal = FSIConfig(self.config['CONFIG_PRIMAL'])
        # Read Adjoint config
        self.configFSIAdjoint = FSIConfig(self.config['CONFIG_ADJOINT'])

    def obj_f(self,x_new):
        
        # Check if new design is needed
        # In case start new design and deform
        self.CheckNewDesign(x_new)
        
        
        #Primal
        
        # return function
        return
        
    def obj_df(self,x_new):
        
        # Check if new design is needed (it won't as Adjoin is performed after primal)        
        # In case start new design and deform
        self.CheckNewDesign(x_new)
        #Adjoint
        
        # return d function
        return 

    def con_ceq(self,x_new):
        
        # Check if new design is needed        
        # In case start new design and deform
        self.CheckNewDesign(x_new)
        
        #Check if Geo has been executed, if it hasn't execute Geo
        self.CheckGeo()
        
        # pulls constraint equality
        c_eq = self.design[self.design_iter].pull_c_eq(self.geo_folder)
        # return ceq
        
        return c_eq
    
    def con_dceq(self,x_new):
        
        # Check if new design is needed (it won't as geo gradient is calculated after ge       
        # In case start new design and deform
        self.CheckNewDesign(x_new)
        #If Geo hasn't been executed execute Geo
        
        # pull gradient of constraint equality
        
        # return dceq
        return
    
    def con_cieq(self,x_new):
        
        # Check if new design is needed        
        # In case start new design and deform
        self.CheckNewDesign(x_new)
        
        #Check if Geo has been executed, if it hasn't execute Geo
        self.CheckGeo()
        
        # pull constraint inequality
        c_ieq = self.design[self.design_iter].pull_c_ieq(self.geo_folder)
        
        # return cieq
        return c_ieq
    
    def con_dcieq(self,x_new):
        
        # Check if new design is needed (it won't as geo gradient is calculated after ge       
        # In case start new design and deform
        self.CheckNewDesign(x_new)
        #If Geo hasn't been executed execute Geo
        
        # pull gradient of constraint inequality
        
        # return dcieq    
        return
        
        
    def clean_previous_designs(self):

        # I remove old designs
        print('Removing old designs...')
        command = 'rm -r ' + self.folder + '/DESIGNS'
        
        # Executes shell command
        run_command(command, 'Remove old designs', False)     
    
    
    def create_design_folder(self):
              
        command = 'mkdir ' + self.folder + '/DESIGNS'

        # Executes shell command
        run_command(command, 'Create design folder', False)  
            
            
    def CheckNewDesign(self, x_in):
        
       if self.design_iter == -1:
           print('Starting new design')
           self.design_iter += 1
           # starting new design
           self.InitializeNewDesign(x_in)
       else:    
          delta = self.design[self.design_iter].x - x_in
          module = np.linalg.norm(delta)
          if module > self.design_toll:
             print('Starting new design')
             self.design_iter += 1
             # starting new design
             self.InitializeNewDesign(x_in)
             # performing mesh deform
             self.DeformMesh()
          else:
             continue
            
    
    def InitializeNewDesign(self,x_in):  
        
        # old design
        x_old = self.design[self.design_iter-1].x
        
        # create design folder
        design_folder = self.folder + '/DESIGNS' + 'DSN_'+ str(int(self.design_iter)).zfill(self.magnord_design)
        command = 'mkdir ' + design_folder

        # Executes shell command
        run_command(command, 'Creating design ' + str(int(self.design_iter)).zfill(self.magnord_design) + ' directory', False)  
            
        # initialize and append new design object    
        self.design.append(Design(self.config,self.configFSIPrimal,self.configFSIAdjoint, design_folder, self.design_iter ,x_in, x_old ))    
        

    def DeformMesh(self):    
        
        # old design
        #x = self.design[self.design_iter].x
        
        # Check if there is the need to deform the mesh
        # It is FALSE if any x of the current design is different than 0
        #all_zeros = not np.any( x )
        
        #if all_zeros == False:
            
        # create folder for analysis
        self.deform_folder = self.folder + '/DESIGNS' + 'DSN_'+ str(int(self.design_iter)).zfill(self.magnord_design) + '/DEFORM'
        command = 'mkdir ' + self.deform_folder        
        # Executes shell command
        run_command(command, 'Creating deform directory for design ' + str(int(self.design_iter)).zfill(self.magnord_design), False)   
           
        # pull config deformation file
        config_deform = self.config['FOLDER'] + '/' + self.config['CONFIG_DEF'] 
        command = 'cp ' + config_deform + ' ' + self.deform_folder + '/'
        run_command(command, 'Pulling deformation config', False)
        
        # pull mesh file       
        mesh_filename = readConfig(self.config['CONFIG_DEF'], 'MESH_FILENAME')
        command = 'cp ' + mesh_filename + ' ' + self.deform_folder + '/'
        run_command(command, 'Pulling mesh config for deformation', False)
        
        # Performing mesh deformation
        self.design[self.design_iter].SU2_DEF(self.deform_folder)
        
        
    def CheckGeo(self):    
       """
       If Geo sensitivities (constraints and gradients). are required, checks if GEO run has been done. If not sets up the folder and performs it.
       """ 
       
       if self.design[self.design_iter].geo == False:
            
           # create folder for analysis
           self.geo_folder = self.folder + '/DESIGNS' + 'DSN_'+ str(int(self.design_iter)).zfill(self.magnord_design) + '/GEO'
           command = 'mkdir ' + self.geo_folder        
           # Executes shell command
           run_command(command, 'Creating GEO directory for design ' + str(int(self.design_iter)).zfill(self.magnord_design), False)      
           
           # pull geo deformation file
           config_geo = self.config['FOLDER'] + '/' + self.config['CONFIG_GEO'] 
           command = 'cp ' + config_geo + ' ' + self.geo_folder + '/'
           run_command(command, 'Pulling geo config', False)     
           # pull mesh file 
           self.SetMesh(self.geo_folder)
           
           # Runs SU2_GEO
           self.design[self.design_iter].SU2_GEO(self.geo_folder)
           
           
    def SetMesh(self, destination_folder):
       """
       Pulls mesh file. If optimization iter is 1, pulling is from project folder. If a deformation occurred, pulling is done from DEFORM folder
       """ 
       
       if self.design[self.design_iter].deformation == False:
          # In case deformation hasn't occurred (first iteration) we need the original mesh file
          mesh_filename = readConfig(self.config['CONFIG_DEF'], 'MESH_FILENAME')
          command = 'cp ' + mesh_filename + ' ' + destination_folder + '/'
          run_command(command, 'Pulling mesh config for deformation', False)
          
       else:
          # in case deform has occurred mesh file is named as output of SU2_DEF and needs to be pulled from the dedicated folder
          mesh_filename = readConfig(self.config['CONFIG_DEF'], 'MESH_OUT_FILENAME')
          command = 'cp ' + self.deform_folder + '/' + mesh_filename + ' ' + destination_folder + '/'
          
           