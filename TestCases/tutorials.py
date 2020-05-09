#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.0.4 "Blackbird"
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

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6
from __future__ import print_function

import sys
from TestCase import TestCase    

def main():
    '''This program runs SU2 and ensures that the output matches specified values. 
       This will be used to do checks when code is pushed to github 
       to make sure nothing is broken. '''

    test_list = []
    
    ######################################
    ### RUN TUTORIAL CASES             ###
    ######################################
     
    ### Compressible Flow
     
    # Inviscid Bump
    tutorial_inv_bump            = TestCase('inviscid_bump_tutorial')
    tutorial_inv_bump.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_Bump"
    tutorial_inv_bump.cfg_file   = "inv_channel.cfg"
    tutorial_inv_bump.test_iter  = 0
    tutorial_inv_bump.test_vals  = [-1.437425, 4.075857, 0.003000, 0.012720] #last 4 columns
    tutorial_inv_bump.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_inv_bump.timeout    = 1600
    tutorial_inv_bump.tol        = 0.00001
    tutorial_inv_bump.no_restart = True
    test_list.append(tutorial_inv_bump)
    
    # Inviscid Wedge
    tutorial_inv_wedge            = TestCase('inviscid_wedge_tutorial')
    tutorial_inv_wedge.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_Wedge"
    tutorial_inv_wedge.cfg_file   = "inv_wedge_HLLC.cfg"
    tutorial_inv_wedge.test_iter  = 0
    tutorial_inv_wedge.test_vals  = [-0.481460, 5.253008, -0.292159, 0.052922] #last 4 columns
    tutorial_inv_wedge.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_inv_wedge.timeout    = 1600
    tutorial_inv_wedge.tol        = 0.00001
    tutorial_inv_wedge.no_restart = True
    test_list.append(tutorial_inv_wedge)
    
    # Inviscid ONERA M6
    tutorial_inv_onera            = TestCase('inviscid_onera_tutorial')
    tutorial_inv_onera.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_ONERAM6"
    tutorial_inv_onera.cfg_file   = "inv_ONERAM6.cfg"
    tutorial_inv_onera.test_iter  = 0
    tutorial_inv_onera.test_vals  = [-5.204928, -4.597762, 0.247124, 0.085734] #last 4 columns
    tutorial_inv_onera.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_inv_onera.timeout    = 1600
    tutorial_inv_onera.tol        = 0.00001
    tutorial_inv_onera.no_restart = True
    test_list.append(tutorial_inv_onera)
    
    # Laminar Cylinder
    tutorial_lam_cylinder            = TestCase('laminar_cylinder_tutorial')
    tutorial_lam_cylinder.cfg_dir    = "../Tutorials/compressible_flow/Laminar_Cylinder"
    tutorial_lam_cylinder.cfg_file   = "lam_cylinder.cfg"
    tutorial_lam_cylinder.test_iter  = 0
    tutorial_lam_cylinder.test_vals  = [-6.162141, -0.699617, 0.186570, 69.267308] #last 4 columns
    tutorial_lam_cylinder.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_lam_cylinder.timeout    = 1600
    tutorial_lam_cylinder.tol        = 0.00001
    tutorial_lam_cylinder.no_restart = True
    test_list.append(tutorial_lam_cylinder)

    # Laminar Flat Plate
    tutorial_lam_flatplate            = TestCase('laminar_flatplate_tutorial')
    tutorial_lam_flatplate.cfg_dir    = "../Tutorials/compressible_flow/Laminar_Flat_Plate"
    tutorial_lam_flatplate.cfg_file   = "lam_flatplate.cfg"
    tutorial_lam_flatplate.test_iter  = 0
    tutorial_lam_flatplate.test_vals  = [-2.821818, 2.657591, -0.400044, 0.029413] #last 4 columns
    tutorial_lam_flatplate.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_lam_flatplate.timeout    = 1600
    tutorial_lam_flatplate.tol        = 0.00001
    tutorial_lam_flatplate.no_restart = True
    test_list.append(tutorial_lam_flatplate)
    
    # Turbulent Flat Plate
    tutorial_turb_flatplate            = TestCase('turbulent_flatplate_tutorial')
    tutorial_turb_flatplate.cfg_dir    = "../Tutorials/compressible_flow/Turbulent_Flat_Plate"
    tutorial_turb_flatplate.cfg_file   = "turb_SA_flatplate.cfg"
    tutorial_turb_flatplate.test_iter  = 0
    tutorial_turb_flatplate.test_vals  = [-2.258584, -4.899502, -0.429387, 0.201236] #last 4 columns
    tutorial_turb_flatplate.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_turb_flatplate.timeout    = 1600
    tutorial_turb_flatplate.tol        = 0.00001
    tutorial_turb_flatplate.no_restart = True
    test_list.append(tutorial_turb_flatplate)
    
    # Transitional FlatPlate
    tutorial_trans_flatplate            = TestCase('transitional_flatplate_tutorial')
    tutorial_trans_flatplate.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Flat_Plate"
    tutorial_trans_flatplate.cfg_file   = "transitional_BC_model_ConfigFile.cfg"
    tutorial_trans_flatplate.test_iter  = 0
    tutorial_trans_flatplate.test_vals  = [-22.021786, -15.330906, 0.000000, 0.023952] #last 4 columns
    tutorial_trans_flatplate.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_trans_flatplate.timeout    = 1600
    tutorial_trans_flatplate.tol        = 0.00001
    tutorial_trans_flatplate.no_restart = True
    test_list.append(tutorial_trans_flatplate)

    # Turbulent ONERA M6
    tutorial_turb_oneram6            = TestCase('turbulent_oneram6_tutorial')
    tutorial_turb_oneram6.cfg_dir    = "../Tutorials/compressible_flow/Turbulent_ONERAM6"
    tutorial_turb_oneram6.cfg_file   = "turb_ONERAM6.cfg"
    tutorial_turb_oneram6.test_iter  = 0
    tutorial_turb_oneram6.test_vals  = [-4.499497, -11.473637, 0.332666, 0.098280] #last 4 columns
    tutorial_turb_oneram6.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_turb_oneram6.timeout    = 1600
    tutorial_turb_oneram6.tol        = 0.00001
    test_list.append(tutorial_turb_oneram6)

    # NICD Nozzle
    tutorial_nicfd_nozzle           = TestCase('nicfd_nozzle')
    tutorial_nicfd_nozzle.cfg_dir   = "../Tutorials/compressible_flow/NICFD_nozzle"
    tutorial_nicfd_nozzle.cfg_file  = "NICFD_nozzle.cfg"
    tutorial_nicfd_nozzle.test_iter = 20
    tutorial_nicfd_nozzle.test_vals = [-1.909272, -6.407306, 4.284191, 0.000000, 0.000000]
    tutorial_nicfd_nozzle.su2_exec  = "mpirun -np 2 SU2_CFD"
    tutorial_nicfd_nozzle.timeout   = 1600
    tutorial_nicfd_nozzle.tol       = 0.00001
    tutorial_nicfd_nozzle.no_restart = True
    test_list.append(tutorial_nicfd_nozzle)
    

    # Unsteady NACA0012
    tutorial_unst_naca0012               = TestCase('unsteady_naca0012')
    tutorial_unst_naca0012.cfg_dir       = "../Tutorials/compressible_flow/Unsteady_NACA0012"
    tutorial_unst_naca0012.cfg_file      = "unsteady_naca0012.cfg"
    tutorial_unst_naca0012.test_iter     = 500
    tutorial_unst_naca0012.test_vals     = [500, 0, 0.302003, 0.665069, -5.300141, 0.000000,  0.0000e+00,  0.0000e+00]
    tutorial_unst_naca0012.su2_exec      = "mpirun -np 2 SU2_CFD"
    tutorial_unst_naca0012.timeout       = 1600
    tutorial_unst_naca0012.tol           = 0.00001
    tutorial_unst_naca0012.unsteady      = True
    test_list.append(tutorial_unst_naca0012)

    ### Design
    
    # Inviscid NACA 0012 Design
    tutorial_design_inv_naca0012            = TestCase('design_inv_naca0012')
    tutorial_design_inv_naca0012.cfg_dir    = "../Tutorials/design/Inviscid_2D_Unconstrained_NACA0012"
    tutorial_design_inv_naca0012.cfg_file   = "inv_NACA0012_basic.cfg"
    tutorial_design_inv_naca0012.test_iter  = 0
    tutorial_design_inv_naca0012.test_vals  = [-3.585391, -2.989014, 0.134515, 0.208523] #last 4 columns
    tutorial_design_inv_naca0012.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_design_inv_naca0012.timeout    = 1600
    tutorial_design_inv_naca0012.tol        = 0.00001
    tutorial_design_inv_naca0012.no_restart = True
    test_list.append(tutorial_design_inv_naca0012)

    # Turbulent RAE 2822 Design
    tutorial_design_turb_rae2822            = TestCase('design_turb_rae2822')
    tutorial_design_turb_rae2822.cfg_dir    = "../Tutorials/design/Turbulent_2D_Constrained_RAE2822"
    tutorial_design_turb_rae2822.cfg_file   = "turb_SA_RAE2822.cfg"
    tutorial_design_turb_rae2822.test_iter  = 0
    tutorial_design_turb_rae2822.test_vals  = [-1.700114, -4.941261, 0.218432, 0.190639] #last 4 columns
    tutorial_design_turb_rae2822.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_design_turb_rae2822.timeout    = 1600
    tutorial_design_turb_rae2822.tol        = 0.00001
    tutorial_design_turb_rae2822.no_restart = True
    test_list.append(tutorial_design_turb_rae2822)

    # Multi Objective Design
    tutorial_design_multiobj            = TestCase('design_multiobj')
    tutorial_design_multiobj.cfg_dir    = "../Tutorials/design/Multi_Objective_Shape_Design"
    tutorial_design_multiobj.cfg_file   = "inv_wedge_ROE_multiobj_combo.cfg"
    tutorial_design_multiobj.test_iter  = 0
    tutorial_design_multiobj.test_vals  = [2.657333, -3.020635, 324840.000000, 0.000000] #last 4 columns
    tutorial_design_multiobj.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_design_multiobj.timeout    = 1600
    tutorial_design_multiobj.tol        = 0.00001
    tutorial_design_multiobj.no_restart = True
    test_list.append(tutorial_design_multiobj)

    ######################################
    ### RUN TESTS                      ###
    ######################################
    
    pass_list = [ test.run_test() for test in test_list ]

    # Tests summary
    print('==================================================================')
    print('Summary of the parallel tests')
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
