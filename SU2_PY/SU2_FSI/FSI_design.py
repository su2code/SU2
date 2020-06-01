#!/usr/bin/env python

## \file FSI_design.py
#  \brief FSI Shape Optimization design.
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
'''
NB:  In the Armijo-GOldstein method, some optimization iteration only forsee the evaluation of the primal (and constraint equation evaluation) 
several times (at different scaled dv_values) to set up a correct step size. 
In those iterations neither adjoint nor graident evaluation of the constraint is evaluated. Only the primal with new intermediate dv values is required + constraint evaluation. 
This operation (mesh_deform + primal) can be performed in the same design folder. 
Can we do that?
Each design is the characterized by:
dv_value_new
dv_value_old
Deform folder
Primal folder
Geo_folder (in principle only for geo constraints, for SOME iterations we can read also gradients)
Adjoint folder

The first between geo and primal which is called every iteration calls the deform
Deform needs only to be done every design (every interation iter)
Usually Primal is done first but at iteration 0 there is a geo extra run (this complicates things)
'''

'''
Design for every iteration needs to memorize
   config         - optimization config file
   config_primal  - Primal config file
   config_adjoint - adjoint config file
   x              - DV_value
   x_old          - DV_Value old (previous design)
   obj_func       - Value of objective function
   gradient       - value of the gradient to pass to potimizator
   constraints    - value of constraints
   constr_grad    - value of constraints gradient
   deform         - booleian for deformation
'''


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys
import numpy as np
from math import pow, factorial
from SU2_FSI.FSI_config import FSIConfig as FSIConfig
from SU2_FSI.FSI_tools import readConfig, run_command, UpdateConfig, DeformMesh, Geometry, ReadGeoConstraints, FSIPrimal, FSIAdjoint, ChainRule, ReadGeoConstraintGradients
# -------------------------------------------------------------------
#  Project Class
# -------------------------------------------------------------------

class Design:
    """

        Starts a design class to manage multiple designs

        Attributes:
             config  - base config for optimization
             folder  - project working folder
             iter    - Optimization iteration number
             designs - list of designs
             folder  - project working folder

    """

    def __init__( self, config, configFSIPrimal,configFSIAdjoint, folder, design_folder, nbr ,dv, x_old ):
        # Attributes:
        self.config  = config      # base config
        self.design_folder  = design_folder      # design folder
        self.folder =  folder            # run folder 
        self.design_nbr = nbr      # current design number
        self.__dv = None
        self._setdv(dv)              # collection of dv_variables of the current design (private)
        self.x_old = x_old          # collection of dv_variable of the previous design (private)
        self.n_dv = len(self.__dv)
        self.configFSIPrimal = configFSIPrimal
        self.configFSIAdjoint = configFSIAdjoint   
        # booleians for every analysis
        self.deformation = False
        self.geo = False
        self.primal = False
        self.adjoint = False
        
        
        self.c_eq = None
        self.c_ieq_plus = None
        self.c_ieq_minus = None
        self.obj_f = None
        self.obj_df = None
        self.c_dieq_plus = None
        self.c_dieq_minus = None

    def _setdv(self,dv):
        self.__dv = dv

    def getdv(self):
        
        return self.__dv
    
    #def getxold(self):
        
    #    return self._x_old   

    def SU2_DEF(self,deform_folder):    
        
        # booleian for deformation
        self.deformation = True
        # before running the command the new DV_VALUE needs to be specified inside the deformation config
        ConfigFileName = deform_folder + '/' + self.config['CONFIG_DEF']
        # It includes the relaxation prameter as done in SU2
        UpdateConfig(ConfigFileName, 'DV_VALUE', self.__dv*float ( self.config['OPT_RELAX_FACTOR'] ))
        
        # performing SU2_DEF
        DeformMesh(deform_folder,self.config['CONFIG_DEF'])


    def FSIPrimal(self,primal_folder):
        
        # booleian for primal
        self.primal = True        

        # Set up current mesh file name
        self.UpdateMeshFilename( primal_folder, self.configFSIPrimal['SU2_CONFIG'])
        
        # perform primal
        FSIPrimal(primal_folder, self.config)
        
    def FSIAdjoint(self,adj_folder):   
        
        self.adjoint = True
        
        # Set up current mesh file name
        self.UpdateMeshFilename( adj_folder, self.configFSIAdjoint['SU2_CONFIG'])
        
        # perform adjoint
        FSIAdjoint(adj_folder, self.config)        

    def SU2_GEO(self, geo_folder):

        # booleian for deformation
        self.geo = True
        # Set up current mesh file name
        self.UpdateMeshFilename( geo_folder, self.config['CONFIG_GEO'])
                
        # Performing SU2_GEO
        Geometry(geo_folder, self.config)
        
    
    def pull_c_eq(self, geo_folder):
        """ 
        Fuction that returns the numpy list of c_eq
        """
        self.c_eq = ReadGeoConstraints(geo_folder, self.config['OPT_CONSTRAINT'], '=', self.design_nbr)

        return self.c_eq
    
    def pull_c_deq(self, geo_folder):
        """ 
        Fuction that returns the numpy matrix of c_deq
        """
        if self.geo == True:
        
           self.dc_eq = ReadGeoConstraintGradients( geo_folder,self.config['OPT_CONSTRAINT'],self.n_dv, '=' )

           # scaling obj_function with gradient factor
           global_factor = float(self.config['OPT_GRADIENT_FACTOR']) 
        
        else:
           print('Cant evaluate c_deq. Geo hasn t run.[pull_c_deq in FSI_design.py]')
           sys.exit()            
            
        return self.dc_eq, global_factor    
                           
    def pull_c_ieq(self, geo_folder):
        """ 
        Fuction that returns the numpy list of c_ieq
        """
        # First > constraints
        self.c_ieq_plus = ReadGeoConstraints(geo_folder, self.config['OPT_CONSTRAINT'], '>', self.design_nbr)

        # After < constraints
        self.c_ieq_minus = ReadGeoConstraints(geo_folder, self.config['OPT_CONSTRAINT'], '<', self.design_nbr)     


        return np.concatenate((self.c_ieq_plus, - self.c_ieq_minus), axis=0) 
    
    def pull_c_dieq(self, geo_folder):
        """ 
        Fuction that returns the numpy matrix of c_dieq
        """
        
        if self.geo == True:
           # First > constraints
           self.c_dieq_plus = ReadGeoConstraintGradients( geo_folder,self.config['OPT_CONSTRAINT'],self.n_dv, '>' )
           # After < constraints   
           self.c_dieq_minus = ReadGeoConstraintGradients( geo_folder,self.config['OPT_CONSTRAINT'],self.n_dv, '<' )

           # scaling obj_function with gradient factor
           global_factor = float(self.config['OPT_GRADIENT_FACTOR']) 
                
           if len(self.c_dieq_plus) !=0 and len(self.c_dieq_minus) !=0: 
              return np.concatenate((self.c_dieq_plus, - self.c_dieq_minus), axis=0), global_factor 
           elif len(self.c_dieq_plus) == 0:
              return -self.c_dieq_minus, global_factor
           elif len(self.c_dieq_minus) == 0:
              return self.c_dieq_plus, global_factor
           else:
               print('Empty constraint gradient.[pull_c_dieq in FSI_design.py]')
               sys.exit()           
           
        else:
           print('Cant evaluate c_dieq. Geo hasn t run.[pull_c_dieq in FSI_design.py]')
           sys.exit()   
           
           
    def pull_obj_df(self,adj_folder,FFD_indexes, PointInv,ffd_degree):
    
        if self.adjoint == True:
            
            obj_df = ChainRule(adj_folder,FFD_indexes, PointInv,ffd_degree)
            
        else:
            print('Can t evaluate obj_df as Adjoint hasn t run')
            sys.exit()    
        
        # scaling obj_function with gradient factor
        global_factor = float(self.config['OPT_GRADIENT_FACTOR']) 
                  
        return obj_df, global_factor

    def UpdateMeshFilename(self, folder, SU2configFile):
        """ 
        Fuction that updates the name of the mesh file input from the config file of GEO, SU2_primal and SU2_Adjoint
        according to possible deformations
        """
        # Set up current mesh file name
        if self.design_nbr == 0:
           mesh_filename = readConfig(self.folder + '/' + self.config['CONFIG_DEF'], 'MESH_FILENAME')
        else:   
           mesh_filename = readConfig(self.folder + '/' + self.config['CONFIG_DEF'], 'MESH_OUT_FILENAME')
           
        ConfigFileName =  folder + '/' + SU2configFile
        UpdateConfig(ConfigFileName, 'MESH_FILENAME', mesh_filename)
        
    def pull_obj_f(self, primal_folder):    
        """ 
        Fuction that returns the numpy list of obj_func multiplied by the scale term
        """           
        
        # check which is the objective function (and scale)
        ConfigFileName =  self.folder + '/' + self.configFSIPrimal['SU2_CONFIG']
        obj_scale = readConfig(ConfigFileName, 'OBJECTIVE_FUNCTION')
        if obj_scale.find('*') != -1:
           string = obj_scale.split('*')
           obj = string[0]
           scale = float(string[1])
        else:   
           obj = obj_scale
           scale = float(1)
          
        # also considering OPT_GRADIENT_FACTOR
        global_factor = float(self.config['OPT_GRADIENT_FACTOR'])
           
        # Now read file
        input_file = open(primal_folder + '/' + 'Objectives.dat')
        # go to line 2
        line = input_file.readline();line = input_file.readline();
        line = line.split()
        cd = float(line[0])
        cl = float(line[1])
        input_file.close()
        # removing this temporary file so that it cannot be accessed in case in the future
        os.remove(primal_folder + '/' + 'Objectives.dat')
        # writing history
        logfileCD = 'HistoryCD.dat'
        logfileCL = 'HistoryCL.dat'
        if iter ==0:
          logCD = open(primal_folder + '/../../' + logfileCD,"w") 
          logCL = open(primal_folder + '/../../' + logfileCL,"w") 
        else:
          logCD = open(primal_folder + '/../../' + logfileCD,"a")  
          logCL = open(primal_folder + '/../../' + logfileCL,"a")
          
        logCD.write( str(cd) + '\n' )
        logCL.write( str(cl) + '\n' )        
          
        logCD.close()
        logCL.close()
        
        
        if obj =='DRAG':
           self.obj_f = cd*scale*global_factor 
           return cd, scale, global_factor
        elif obj =='LIFT':
           self.obj_f = cl*scale* global_factor
           return cl, scale, global_factor
        else:
            print('Not implemented objective function.[pull_obj_f in FSI_design.py]')
            sys.exit()
            
    