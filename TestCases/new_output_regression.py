#!/usr/bin/env python

## \file serial_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 6.1.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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

from __future__ import print_function, division, absolute_import
import sys
from TestCase import TestCase

def main():
    '''This program runs SU2 and ensures that the output matches specified values. 
       This will be used to do checks when code is pushed to github 
       to make sure nothing is broken. '''

    test_list = []

    ##########################
    ### FEA - FSI          ###
    ##########################

    # Static beam, 3d
    statbeam3d           = TestCase('statbeam3d')
    statbeam3d.cfg_dir   = "fea_fsi/StatBeam_3d"
    statbeam3d.cfg_file  = "configBeam_3d.cfg"
    statbeam3d.new_output= True
    statbeam3d.test_iter = 0
    statbeam3d.test_vals = [-8.498274, -8.230638, -8.123824, 6.4095e+04] #last 4 columns
    statbeam3d.su2_exec  = "SU2_RAD"
    statbeam3d.timeout   = 1600
    statbeam3d.tol       = 0.00001
    test_list.append(statbeam3d)

    # # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.unsteady  = True
    dynbeam2d.new_output= True  
    dynbeam2d.test_iter = 6 
    dynbeam2d.test_vals = [-9.420640, -5.365872, -12.430382, 6.5210e+04] #last 4 columns
    dynbeam2d.su2_exec  = "SU2_RAD"
    dynbeam2d.timeout   = 1600
    dynbeam2d.tol       = 0.00001
    test_list.append(dynbeam2d)

    # # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI.cfg"
    fsi2d.unsteady  = True
    fsi2d.new_output= True     
    fsi2d.test_iter = 4
    fsi2d.test_vals = [2.000000, 4.000000, -4.017814, -9.089686] #last 4 columns
    fsi2d.su2_exec  = "SU2_RAD"
    fsi2d.timeout   = 1600
    fsi2d.tol       = 0.00001
    test_list.append(fsi2d)    

    # # FSI, 2D airfoil with RBF interpolation
    # airfoilRBF           = TestCase('airfoil_fsi_rbf')
    # airfoilRBF.cfg_dir   = "fea_fsi/Airfoil_RBF"
    # airfoilRBF.cfg_file  = "config.cfg"
    # airfoilRBF.test_iter = 50
    # airfoilRBF.test_vals = [-8.000964, -2.600088, 0.276433, 0.000824] #last 4 columns
    # airfoilRBF.su2_exec  = "SU2_CFD"
    # airfoilRBF.timeout   = 1600
    # airfoilRBF.tol       = 0.00001
    # test_list.append(airfoilRBF)
   
    # ##########################
    # ### Zonal multiphysics ###
    # ##########################

    # # CHT incompressible
    # cht_incompressible           = TestCase('cht_incompressible')
    # cht_incompressible.cfg_dir   = "coupled_cht/incompressible"
    # cht_incompressible.cfg_file  = "config.cfg"
    # cht_incompressible.test_iter = 10
    # cht_incompressible.test_vals = [0.000000, 0.000000, -7.685301, -12947.783696] #last 4 columns
    # cht_incompressible.su2_exec  = "SU2_CFD"
    # cht_incompressible.timeout   = 1600
    # cht_incompressible.tol       = 0.00001
    # test_list.append(cht_incompressible)

    # ###############################
    # ### Radiative Heat Transfer ###
    # ###############################    

    # # FSI, 2d
    p1rad           = TestCase('p1rad')
    p1rad.cfg_dir   = "radiation/p1model"
    p1rad.cfg_file  = "config.cfg"
    p1rad.new_output= True     
    p1rad.test_iter = 100
    p1rad.test_vals = [-7.750053, -7.917666, -2.118093, 0.092236] #last 4 columns
    p1rad.su2_exec  = "SU2_RAD"
    p1rad.timeout   = 1600
    p1rad.tol       = 0.00001
    test_list.append(p1rad)       

    ######################################
    ### RUN TESTS                      ###
    ######################################  

    pass_list = [ test.run_test() for test in test_list ]

    # Tests summary
    print('==================================================================')
    print('Summary of the serial tests')
    print('python version:', sys.version)
    for i, test in enumerate(test_list):
        if (pass_list[i]):
            print('  passed - %s'%test.tag)
        else:
            print('* FAILED - %s'%test.tag)
    
    if all(pass_list):
        sys.exit(0)
    else:
        sys.exit(1)
    # done

if __name__ == '__main__':
    main()
