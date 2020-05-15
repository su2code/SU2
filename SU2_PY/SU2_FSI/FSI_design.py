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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from math import pow, factorial
from SU2_FSI.FSI_config import FSIConfig as FSIConfig
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

    def __init__( self, config,configFSIPrimal,configFSIAdjoint,  folder, nbr ,x, x_old ):

        
        self.config  = config      # base config
        self.folder  = folder      # design folder
        self.design_nbr = nbr      # current design number
        self.x = x              # collection of dv_variables of the current design
        self.x_old = x_old          # collection of dv_variable of the previous design
        self.configFSIPrimal = configFSIPrimal
        self.configFSIAdjoint = configFSIAdjoint        

    def create_design_folders(self):

        # I create design folders
        command_list = []
        command_list.append('mkdir ' + self.folder + '/DESIGNS/DSN_' + str(self.iter) )
        command_list.append('mkdir ' + self.folder + '/DESIGNS/DSN_' + str(self.iter) + '/Primal'  )
        command_list.append('mkdir ' + self.folder + '/DESIGNS/DSN_' + str(self.iter) + '/Adjoint' )
        command_list.append('mkdir ' + self.folder + '/DESIGNS/DSN_' + str(self.iter) + '/GEO'     )
        command_list.append('mkdir ' + self.folder + '/DESIGNS/DSN_' + str(self.iter) + '/DEFORM')

        for command in command_list:
           try:
               os.system(command)
           except TypeWarning as exception:
               print('Creating new design folder failed with the following warning: {}'.format(TypeWarning))

        return

    def update_iter(self):

        self.iter = self.iter + 1

        return


    def copy_files(self ,analysis):
        '''
        Copy in the dedicated folder input files to perform primal or Adjoint (everything but
        mesh file which is updated in deform)
        analysis = Primal/Adjoint
        '''
        if (analysis == 'Primal'):
            config = self.configFSIPrimal
        elif (analysis == 'Adjoint'):
            config = self.configFSIAdjoint

        command_list = []
        #config
        command_list.append('cp ' + self.config['CONFIG_PRIMAL'] + ' ' + folder + '/DESIGNS/DSN_' + str(self.iter) + '/' + analysis + '/' )
        # Flow file
        command_list.append('cp ' + config['SU2_CONFIG'] + ' ' + folder + '/DESIGNS/DSN_' + str(self.iter) + '/' + analysis + '/' )
        # structural
        command_list.append('cp ' + config['PYBEAM_CONFIG'] + ' ' + folder + '/DESIGNS/DSN_' + str(self.iter) + '/' + analysis + '/' )
        # Spline
        command_list.append('cp ' + config['MLS_CONFIG_FILE_NAME'] + ' ' + folder + '/DESIGNS/DSN_' + str(self.iter) + '/' + analysis + '/' )
        # structural mesh (let's assume this name doesn't change for now)
        command_list.append('cp ' + 'mesh.pyBeam' + ' ' + folder + '/DESIGNS/DSN_' + str(self.iter) + '/' + analysis + '/' )
        # structural properties (let's assume this name doesn't change for now)
        command_list.append('cp ' + 'property.pyBeam' + ' ' + folder + '/DESIGNS/DSN_' + str(self.iter) + '/' + analysis + '/' )


    def SU2_GEO(self):
         '''
         Copy in the dedicated folder input files to perform primal or Adjoint (everything but
         mesh file which is updated in deform)
         analysis = Primal/Adjoint
         '''

         return