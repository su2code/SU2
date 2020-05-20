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

import numpy as np
from math import pow, factorial
from SU2_FSI.FSI_config import FSIConfig as FSIConfig
from SU2_FSI.FSI_tools import run_command, UpdateConfig, DeformMesh, Geometry
# -------------------------------------------------------------------
#  Project Class
# -------------------------------------------------------------------

class Design(object):
    """

        Starts a design class to manage multiple designs

        Attributes:
             config  - base config for optimization
             folder  - project working folder
             iter    - Optimization iteration number
             designs - list of designs
             folder  - project working folder

    """

    def __init__( self, config, configFSIPrimal,configFSIAdjoint,  folder, nbr ,x, x_old ):

        # Attributes:
        self.config  = config      # base config
        self.folder  = folder      # design folder
        self.design_nbr = nbr      # current design number
        self.x = x              # collection of dv_variables of the current design
        self.x_old = x_old          # collection of dv_variable of the previous design
        self.configFSIPrimal = configFSIPrimal
        self.configFSIAdjoint = configFSIAdjoint   
        # booleians for every analysis
        self.deformation = False
        self.geo = False
        self.primal = False
        self.Adjoint = False
        self.c_eq
        self.c_ieq

    def SU2_DEF(self,deform_folder):    
        
        # booleian for deformation
        self.deformation = True
        # before running the command the new DV_VALUE needs to be specified inside the deformation config
        ConfigFileName = deform_folder + '/' + self.config['CONFIG_DEF']
        UpdateConfig(ConfigFileName, 'DV_VALUE', self.x)
        
        # performing SU2_DEF
        DeformMesh(deform_folder,ConfigFileName)





    def SU2_GEO(self,geo_folder):

        # booleian for deformation
        self.geo = True
        # Set up current mesh file name
        mesh_filename = readConfig(self.config['CONFIG_DEF'], 'MESH_OUT_FILENAME')
        ConfigFileName = geo_folder + '/' + self.config['CONFIG_GEO']
        UpdateConfig(ConfigFileName, 'DV_VALUE', mesh_filename)
        
        # Performing SU2_GEO
        Geometry(geo_folder,ConfigFileName)
        
    
    def pull_c_eq(self,geo_folder):
     """ 
     Fuction that returns the numpy list of c_eq
     """
     c_eq = ReadGeoConstraints(geo_folder,self.config['OPT_CONSTRAINT'], '=')
     
     return c_eq
 
    def pull_c_ieq(self,geo_folder):
     """ 
     Fuction that returns the numpy list of c_ieq
     """
     # First > constraints
     c_ieq_plus = ReadGeoConstraints(geo_folder,self.config['OPT_CONSTRAINT'], '>')
     
     # After < constraints
     c_ieq_minus = ReadGeoConstraints(geo_folder,self.config['OPT_CONSTRAINT'], '<')     
     
     return np.concatenate((c_ieq_plus, - c_ieq_minus), axis = 0)
     
        