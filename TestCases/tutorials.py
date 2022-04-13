#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.3.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

    ### CHT

    # CHT incompressible unsteady
    cht_incompressible_unsteady           = TestCase('cht_incompressible_unsteady')
    cht_incompressible_unsteady.cfg_dir   = "../Tutorials/multiphysics/unsteady_cht/"
    cht_incompressible_unsteady.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible_unsteady.test_iter = 2
    cht_incompressible_unsteady.test_vals = [-2.659390, -2.533160, -0.080399, -0.080399, -0.080399, -12.421450,  0.0000e+00, 0.0, 0.0, 0.0, 0.0000e+00, 2.3824e+02] #last columns
    cht_incompressible_unsteady.su2_exec  = "mpirun -n 2 SU2_CFD"
    cht_incompressible_unsteady.timeout   = 1600
    cht_incompressible_unsteady.multizone = True
    cht_incompressible_unsteady.unsteady  = True
    cht_incompressible_unsteady.tol       = 0.00001
    test_list.append(cht_incompressible_unsteady)

    # CHT incompressible, 2D, 3 pins in crossflow
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "../Tutorials/multiphysics/steady_cht"
    cht_incompressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [-2.128826, -0.588813, -0.588813, -0.588813] #last 4 columns
    cht_incompressible.su2_exec  = "SU2_CFD"
    cht_incompressible.timeout   = 1600
    cht_incompressible.multizone = True
    cht_incompressible.tol       = 0.00001
    test_list.append(cht_incompressible)

    ### Incompressible Flow

    # 2D pin case massflow periodic with heatflux BC and prescribed extracted outlet heat
    sp_pinArray_2d_mf_hf           = TestCase('sp_pinArray_2d_mf_hf')
    sp_pinArray_2d_mf_hf.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Streamwise_Periodic"
    sp_pinArray_2d_mf_hf.cfg_file  = "sp_pinArray_2d_mf_hf.cfg"
    sp_pinArray_2d_mf_hf.test_iter = 25
    sp_pinArray_2d_mf_hf.test_vals = [-4.626384, 1.444465, -0.750978, 241.757337] #last 4 lines
    sp_pinArray_2d_mf_hf.su2_exec  = "mpirun -n 2 SU2_CFD"
    sp_pinArray_2d_mf_hf.timeout   = 1600
    sp_pinArray_2d_mf_hf.tol       = 0.00001
    test_list.append(sp_pinArray_2d_mf_hf)

    # 2D pin case pressure drop periodic with heatflux BC and temperature periodicity
    sp_pinArray_2d_dp_hf_tp           = TestCase('sp_pinArray_2d_dp_hf_tp')
    sp_pinArray_2d_dp_hf_tp.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Streamwise_Periodic"
    sp_pinArray_2d_dp_hf_tp.cfg_file  = "sp_pinArray_2d_dp_hf_tp.cfg"
    sp_pinArray_2d_dp_hf_tp.test_iter = 25
    sp_pinArray_2d_dp_hf_tp.test_vals = [-4.667133, 1.395801, -0.709306, 208.023676] #last 4 lines
    sp_pinArray_2d_dp_hf_tp.su2_exec  = "mpirun -n 2 SU2_CFD"
    sp_pinArray_2d_dp_hf_tp.timeout   = 1600
    sp_pinArray_2d_dp_hf_tp.tol       = 0.00001
    test_list.append(sp_pinArray_2d_dp_hf_tp)

    ### Species Transport

    # 3 species (2 eq) primitive venturi mixing
    species3_primitiveVenturi           = TestCase('species3_primitiveVenturi')
    species3_primitiveVenturi.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Species_Transport"
    species3_primitiveVenturi.cfg_file  = "species3_primitiveVenturi.cfg"
    species3_primitiveVenturi.test_iter = 50
    species3_primitiveVenturi.test_vals = [-6.028145, -5.258104, -5.107927, -5.922051, -1.582604, -6.314220, -6.431771, 5, -0.808615, 5, -2.351160, 5, -0.288300, 1.645644, 0.499064, 0.601230, 0.545351]
    species3_primitiveVenturi.su2_exec  = "mpirun -n 2 SU2_CFD"
    species3_primitiveVenturi.timeout   = 1600
    species3_primitiveVenturi.new_output = True
    species3_primitiveVenturi.tol       = 0.00001
    test_list.append(species3_primitiveVenturi)


    # 3 species (2 eq) primitive venturi mixing
    DAspecies3_primitiveVenturi           = TestCase('DAspecies3_primitiveVenturi')
    DAspecies3_primitiveVenturi.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Species_Transport"
    DAspecies3_primitiveVenturi.cfg_file  = "DAspecies3_primitiveVenturi.cfg"
    DAspecies3_primitiveVenturi.test_iter = 50
    DAspecies3_primitiveVenturi.test_vals = [-8.519150, -7.786969, -7.774848, -7.474167, -12.127149, -12.262476, -11.456643]
    DAspecies3_primitiveVenturi.su2_exec  = "mpirun -n 2 SU2_CFD_AD"
    DAspecies3_primitiveVenturi.timeout   = 1600
    DAspecies3_primitiveVenturi.new_output = True
    DAspecies3_primitiveVenturi.tol       = 0.00001
    test_list.append(DAspecies3_primitiveVenturi)

    ### Compressible Flow

    # Inviscid Bump
    tutorial_inv_bump            = TestCase('inviscid_bump_tutorial')
    tutorial_inv_bump.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_Bump"
    tutorial_inv_bump.cfg_file   = "inv_channel.cfg"
    tutorial_inv_bump.test_iter  = 0
    tutorial_inv_bump.test_vals  = [-1.437425, 4.075857, 0.005439, 0.012998]
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
    tutorial_inv_wedge.test_vals  = [-0.481460, 5.253008, -0.291747, 0.052515]
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
    tutorial_inv_onera.test_vals  = [-5.204928, -4.597762, 0.247451, 0.085770]
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
    tutorial_lam_cylinder.test_vals  = [-6.162141, -0.699617, 0.125776, 69.613563]
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
    tutorial_turb_flatplate.test_vals  = [-2.258584, -4.899502, -0.429375, 0.201236]
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
    tutorial_trans_flatplate.test_vals  = [-22.021786, -15.330766, 0.000000, 0.023952] #last 4 columns
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
    tutorial_turb_oneram6.test_vals  = [-4.564441, -11.524277, 0.327954, 0.097349]
    tutorial_turb_oneram6.su2_exec   = "mpirun -np 2 SU2_CFD"
    tutorial_turb_oneram6.timeout    = 1600
    tutorial_turb_oneram6.tol        = 0.00001
    test_list.append(tutorial_turb_oneram6)

    # NICD Nozzle
    tutorial_nicfd_nozzle           = TestCase('nicfd_nozzle')
    tutorial_nicfd_nozzle.cfg_dir   = "../Tutorials/compressible_flow/NICFD_nozzle"
    tutorial_nicfd_nozzle.cfg_file  = "NICFD_nozzle.cfg"
    tutorial_nicfd_nozzle.test_iter = 20
    tutorial_nicfd_nozzle.test_vals = [-2.187400, -9.409241, 3.477513, 0.000000, 0.000000]
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

    # PROPELLER VARIBLE LOAD
    propeller_var_load           = TestCase('propeller_variable_load')
    propeller_var_load.cfg_dir   = "../Tutorials/compressible_flow/ActuatorDisk_VariableLoad"
    propeller_var_load.cfg_file  = "propeller_variable_load.cfg"
    propeller_var_load.test_iter = 20
    propeller_var_load.test_vals = [-1.830252, -4.535038, -0.000323, 0.171648]
    propeller_var_load.su2_exec  = "mpirun -np 2 SU2_CFD"
    propeller_var_load.timeout   = 3200
    propeller_var_load.tol       = 0.00001
    test_list.append(propeller_var_load)

    ### Design

    # Inviscid NACA 0012 Design
    tutorial_design_inv_naca0012            = TestCase('design_inv_naca0012')
    tutorial_design_inv_naca0012.cfg_dir    = "../Tutorials/design/Inviscid_2D_Unconstrained_NACA0012"
    tutorial_design_inv_naca0012.cfg_file   = "inv_NACA0012_basic.cfg"
    tutorial_design_inv_naca0012.test_iter  = 0
    tutorial_design_inv_naca0012.test_vals  = [-3.585391, -2.989014, 0.135070, 0.208565]
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
    tutorial_design_turb_rae2822.test_vals  = [-1.700114, -4.941305, 0.218348, 0.190357]
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
