#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 8.1.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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
    cht_incompressible_unsteady.test_vals = [-2.659390, -2.533160, -0.080399, -0.080399, -0.080399, -12.421450, 0.000000, 0, 0, 0, 0, 2.3824e+02] #last columns
    cht_incompressible_unsteady.multizone = True
    cht_incompressible_unsteady.unsteady  = True
    test_list.append(cht_incompressible_unsteady)

    # CHT incompressible, 2D, 3 pins in crossflow
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "../Tutorials/multiphysics/steady_cht"
    cht_incompressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [-2.128826, -0.588813, -0.588813, -0.588813] #last 4 columns
    cht_incompressible.command   = TestCase.Command(exec = "SU2_CFD")
    cht_incompressible.multizone = True
    test_list.append(cht_incompressible)

    # Solid-to-solid and solid-to-fluid CHT with contact resistance
    cht_CR           = TestCase('cht_solid_solid')
    cht_CR.cfg_dir   = "../Tutorials/multiphysics/contact_resistance_cht"
    cht_CR.cfg_file  = "master.cfg"
    cht_CR.test_iter = 80
    cht_CR.test_vals = [ -8.857438, -9.377593, -10.097769, -2.122358]
    cht_CR.multizone = True
    test_list.append(cht_CR)

    ### Incompressible Flow

    # 2D pin case massflow periodic with heatflux BC and prescribed extracted outlet heat
    sp_pinArray_2d_mf_hf           = TestCase('sp_pinArray_2d_mf_hf')
    sp_pinArray_2d_mf_hf.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Streamwise_Periodic"
    sp_pinArray_2d_mf_hf.cfg_file  = "sp_pinArray_2d_mf_hf.cfg"
    sp_pinArray_2d_mf_hf.test_iter = 25
    sp_pinArray_2d_mf_hf.test_vals = [-4.613682, 1.471278, -0.748987, 241.667177]
    test_list.append(sp_pinArray_2d_mf_hf)

    # 2D pin case pressure drop periodic with heatflux BC and temperature periodicity
    sp_pinArray_2d_dp_hf_tp           = TestCase('sp_pinArray_2d_dp_hf_tp')
    sp_pinArray_2d_dp_hf_tp.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Streamwise_Periodic"
    sp_pinArray_2d_dp_hf_tp.cfg_file  = "sp_pinArray_2d_dp_hf_tp.cfg"
    sp_pinArray_2d_dp_hf_tp.test_iter = 25
    sp_pinArray_2d_dp_hf_tp.test_vals = [-4.640621, 1.436697, -0.707302, 208.023676]
    test_list.append(sp_pinArray_2d_dp_hf_tp)

    ### Species Transport

    # 3 species (2 eq) primitive venturi mixing
    species3_primitiveVenturi           = TestCase('species3_primitiveVenturi')
    species3_primitiveVenturi.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Species_Transport"
    species3_primitiveVenturi.cfg_file  = "species3_primitiveVenturi.cfg"
    species3_primitiveVenturi.test_iter = 50
    species3_primitiveVenturi.test_vals = [-5.869509, -5.252493, -5.127926, -5.912790, -1.767067, -6.152558, -6.304196, 5.000000, -0.933280, 5.000000, -2.314730, 5.000000, -0.680255, 1.649865, 0.500678, 0.596475, 0.552712]
    test_list.append(species3_primitiveVenturi)

    # 3 species (2 eq) primitive venturi mixing
    DAspecies3_primitiveVenturi           = TestCase('DAspecies3_primitiveVenturi')
    DAspecies3_primitiveVenturi.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Species_Transport"
    DAspecies3_primitiveVenturi.cfg_file  = "DAspecies3_primitiveVenturi.cfg"
    DAspecies3_primitiveVenturi.test_iter = 50
    DAspecies3_primitiveVenturi.test_vals         = [-7.584508, -7.211527, -6.740742, -6.896386, -11.472089, -10.865347, -10.096770]
    DAspecies3_primitiveVenturi.test_vals_aarch64 = [-7.865411, -7.548131, -7.347978, -7.217536, -11.822422, -10.968444, -10.193225]
    DAspecies3_primitiveVenturi.command   = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    test_list.append(DAspecies3_primitiveVenturi)

    # 2 species (1 eq) kenics static mixer for composition-dependent model
    kenics_mixer_tutorial           = TestCase('kenics_mixer_tutorial')
    kenics_mixer_tutorial.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Species_Transport_Composition_Dependent_Model"
    kenics_mixer_tutorial.cfg_file  = "kenics_mixer_tutorial.cfg"
    kenics_mixer_tutorial.test_iter = 10
    kenics_mixer_tutorial.test_vals = [-7.489796, -6.823473, -6.838032, -5.157793, -7.911512, -3.166212, -7.447901, 5.000000, -1.862212, 4.000000, -5.131048, 3.000000, -6.380193, 0.025472, 0.000000, 0.025472, 0.000000, 64.105000, 8.479200, 48.099000, 7.526700]
    kenics_mixer_tutorial.command   = TestCase.Command("mpirun -n 2", "SU2_CFD")
    test_list.append(kenics_mixer_tutorial)

    # 90 degree pipe bend with wall functions from the experiments of Sudo et al.
    sudo_tutorial           = TestCase('sudo_bend')
    sudo_tutorial.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Turbulent_Bend_Wallfunctions"
    sudo_tutorial.cfg_file  = "sudo.cfg"
    sudo_tutorial.test_iter = 10
    sudo_tutorial.test_vals = [-14.579462, -13.203791, -13.601782, -12.616876, -14.005299, -10.817605, 15.000000, -2.296083]
    sudo_tutorial.command   = TestCase.Command("mpirun -n 2", "SU2_CFD")
    test_list.append(sudo_tutorial)

    ### Incompressible Combustion

    # Pre-mixed, laminar hydrogen flame with heat loss
    premixed_hydrogen = TestCase('premixed_hydrogen')
    premixed_hydrogen.cfg_dir = "../Tutorials/incompressible_flow/Inc_Combustion/1__premixed_hydrogen"
    premixed_hydrogen.cfg_file = "H2_burner.cfg"
    premixed_hydrogen.test_iter = 10
    premixed_hydrogen.test_vals = [-9.809794, -10.369804, -11.044267, -4.332945, -11.883789]
    test_list.append(premixed_hydrogen)

    ### Compressible Flow

    # Inviscid Bump
    tutorial_inv_bump            = TestCase('inviscid_bump_tutorial')
    tutorial_inv_bump.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_Bump"
    tutorial_inv_bump.cfg_file   = "inv_channel.cfg"
    tutorial_inv_bump.test_iter  = 0
    tutorial_inv_bump.test_vals  = [-1.548003, 3.983585, 0.020973, 0.071064]
    test_list.append(tutorial_inv_bump)

    # Inviscid Wedge
    tutorial_inv_wedge            = TestCase('inviscid_wedge_tutorial')
    tutorial_inv_wedge.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_Wedge"
    tutorial_inv_wedge.cfg_file   = "inv_wedge_HLLC.cfg"
    tutorial_inv_wedge.test_iter  = 0
    tutorial_inv_wedge.test_vals  = [-0.864206, 4.850246, -0.245674, 0.043209]
    tutorial_inv_wedge.no_restart = True
    test_list.append(tutorial_inv_wedge)

    # Inviscid ONERA M6
    tutorial_inv_onera            = TestCase('inviscid_onera_tutorial')
    tutorial_inv_onera.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_ONERAM6"
    tutorial_inv_onera.cfg_file   = "inv_ONERAM6.cfg"
    tutorial_inv_onera.test_iter  = 0
    tutorial_inv_onera.test_vals  = [-5.504789, -4.895776, 0.249157, 0.118834]
    tutorial_inv_onera.no_restart = True
    test_list.append(tutorial_inv_onera)

    # Laminar Cylinder
    tutorial_lam_cylinder            = TestCase('laminar_cylinder_tutorial')
    tutorial_lam_cylinder.cfg_dir    = "../Tutorials/compressible_flow/Laminar_Cylinder"
    tutorial_lam_cylinder.cfg_file   = "lam_cylinder.cfg"
    tutorial_lam_cylinder.test_iter  = 0
    tutorial_lam_cylinder.test_vals  = [-6.162141, -0.699617, 0.126007, 69.619462]
    tutorial_lam_cylinder.no_restart = True
    test_list.append(tutorial_lam_cylinder)

    # Laminar Flat Plate
    tutorial_lam_flatplate            = TestCase('laminar_flatplate_tutorial')
    tutorial_lam_flatplate.cfg_dir    = "../Tutorials/compressible_flow/Laminar_Flat_Plate"
    tutorial_lam_flatplate.cfg_file   = "lam_flatplate.cfg"
    tutorial_lam_flatplate.test_iter  = 0
    tutorial_lam_flatplate.test_vals  = [-2.821818, 2.657591, -0.400044, 0.029365] #last 4 columns
    tutorial_lam_flatplate.no_restart = True
    test_list.append(tutorial_lam_flatplate)

    # Turbulent Flat Plate
    tutorial_turb_flatplate            = TestCase('turbulent_flatplate_tutorial')
    tutorial_turb_flatplate.cfg_dir    = "../Tutorials/compressible_flow/Turbulent_Flat_Plate"
    tutorial_turb_flatplate.cfg_file   = "turb_SA_flatplate.cfg"
    tutorial_turb_flatplate.test_iter  = 0
    tutorial_turb_flatplate.test_vals  = [-2.258584, -4.901015, -0.429373, 0.201034]
    tutorial_turb_flatplate.no_restart = True
    test_list.append(tutorial_turb_flatplate)

    # Transitional FlatPlate
    tutorial_trans_flatplate            = TestCase('transitional_flatplate_tutorial')
    tutorial_trans_flatplate.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Flat_Plate"
    tutorial_trans_flatplate.cfg_file   = "transitional_BC_model_ConfigFile.cfg"
    tutorial_trans_flatplate.test_iter  = 0
    tutorial_trans_flatplate.test_vals  = [-22.021786, -15.330766, 0.000000, 0.023952] #last 4 columns
    tutorial_trans_flatplate.no_restart = True
    test_list.append(tutorial_trans_flatplate)

    # Transitional FlatPlate T3A
    tutorial_trans_flatplate_T3A            = TestCase('transitional_flatplate_tutorial_T3A')
    tutorial_trans_flatplate_T3A.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Flat_Plate/Langtry_and_Menter/T3A"
    tutorial_trans_flatplate_T3A.cfg_file   = "transitional_LM_model_ConfigFile.cfg"
    tutorial_trans_flatplate_T3A.test_iter  = 20
    tutorial_trans_flatplate_T3A.test_vals  = [-5.837399, -2.092246, -3.983493, -0.302381, -1.920868, 1.667180, -3.496278, 0.391608]
    tutorial_trans_flatplate_T3A.test_vals_aarch64 = [-5.837368, -2.092246, -3.984172, -0.302357, -1.928108, 1.667157, -3.496279, 0.391610]
    tutorial_trans_flatplate_T3A.no_restart = True
    test_list.append(tutorial_trans_flatplate_T3A)

    # Transitional FlatPlate T3Am
    tutorial_trans_flatplate_T3Am            = TestCase('transitional_flatplate_tutorial_T3Am')
    tutorial_trans_flatplate_T3Am.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Flat_Plate/Langtry_and_Menter/T3A-"
    tutorial_trans_flatplate_T3Am.cfg_file   = "transitional_LM_model_ConfigFile.cfg"
    tutorial_trans_flatplate_T3Am.test_iter  = 20
    tutorial_trans_flatplate_T3Am.test_vals  = [-6.063700, -1.945073, -3.946836, -0.549147, -3.863792, 2.664440, -2.517610, 1.112977]
    tutorial_trans_flatplate_T3Am.test_vals_aarch64 = [-6.063726, -1.945088, -3.946923, -0.549166, -3.863794, 2.664439, -2.517601, 1.112978]
    tutorial_trans_flatplate_T3Am.no_restart = True
    test_list.append(tutorial_trans_flatplate_T3Am)

    # Transitional E387 SA
    tutorial_trans_e387_sa            = TestCase('tutorial_trans_e387_sa')
    tutorial_trans_e387_sa.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Airfoil/Langtry_and_Menter/E387"
    tutorial_trans_e387_sa.cfg_file   = "transitional_SA_LM_model_ConfigFile.cfg"
    tutorial_trans_e387_sa.test_iter  = 20
    tutorial_trans_e387_sa.test_vals  = [-6.527027, -5.081541, -0.795261, 1.022607, 0.150175, 2.000000, -9.580660]
    tutorial_trans_e387_sa.no_restart = True
    test_list.append(tutorial_trans_e387_sa)

    # Transitional E387 SST
    tutorial_trans_e387_sst            = TestCase('tutorial_trans_e387_sst')
    tutorial_trans_e387_sst.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Airfoil/Langtry_and_Menter/E387"
    tutorial_trans_e387_sst.cfg_file   = "transitional_SST_LM_model_ConfigFile.cfg"
    tutorial_trans_e387_sst.test_iter  = 20
    tutorial_trans_e387_sst.test_vals  = [-6.532424, -2.932985, 0.401483, 1.078077, 0.188212, 2.000000, -10.005652]
    tutorial_trans_e387_sst.no_restart = True
    test_list.append(tutorial_trans_e387_sst)

    # Turbulent ONERA M6
    tutorial_turb_oneram6            = TestCase('turbulent_oneram6_tutorial')
    tutorial_turb_oneram6.cfg_dir    = "../Tutorials/compressible_flow/Turbulent_ONERAM6"
    tutorial_turb_oneram6.cfg_file   = "turb_ONERAM6.cfg"
    tutorial_turb_oneram6.test_iter  = 0
    tutorial_turb_oneram6.test_vals  = [-4.564441, -11.524295, 0.327905, 0.097340]
    test_list.append(tutorial_turb_oneram6)

    # NICD Nozzle
    tutorial_nicfd_nozzle           = TestCase('nicfd_nozzle')
    tutorial_nicfd_nozzle.cfg_dir   = "../Tutorials/compressible_flow/NICFD_nozzle"
    tutorial_nicfd_nozzle.cfg_file  = "NICFD_nozzle.cfg"
    tutorial_nicfd_nozzle.test_iter = 20
    tutorial_nicfd_nozzle.test_vals = [-2.056675, -2.124123, 3.687027, 0.000000, 0.000000]
    tutorial_nicfd_nozzle.no_restart = True
    test_list.append(tutorial_nicfd_nozzle)

    # Unsteady NACA0012
    tutorial_unst_naca0012               = TestCase('unsteady_naca0012')
    tutorial_unst_naca0012.cfg_dir       = "../Tutorials/compressible_flow/Unsteady_NACA0012"
    tutorial_unst_naca0012.cfg_file      = "unsteady_naca0012.cfg"
    tutorial_unst_naca0012.test_iter     = 520
    tutorial_unst_naca0012.test_vals         = [520.000000, 0.000000, -5.290694, 0.000000, 0.317272, 0.820972, 0.002144, 0.012805]
    tutorial_unst_naca0012.test_vals_aarch64 = [520.000000, 0.000000, -5.298777, 0.000000, 0.288956, 0.736706, 0.002419, 0.007134]
    tutorial_unst_naca0012.unsteady      = True
    test_list.append(tutorial_unst_naca0012)

    # PROPELLER VARIBLE LOAD
    propeller_var_load           = TestCase('propeller_variable_load')
    propeller_var_load.cfg_dir   = "../Tutorials/compressible_flow/ActuatorDisk_VariableLoad"
    propeller_var_load.cfg_file  = "propeller_variable_load.cfg"
    propeller_var_load.test_iter = 20
    propeller_var_load.test_vals = [-1.830257, -4.535041, -0.000323, 0.171647]
    propeller_var_load.timeout   = 3200
    test_list.append(propeller_var_load)

    ### Design

    # Inviscid NACA 0012 Design
    tutorial_design_inv_naca0012            = TestCase('design_inv_naca0012')
    tutorial_design_inv_naca0012.cfg_dir    = "../Tutorials/design/Inviscid_2D_Unconstrained_NACA0012"
    tutorial_design_inv_naca0012.cfg_file   = "inv_NACA0012_basic.cfg"
    tutorial_design_inv_naca0012.test_iter  = 0
    tutorial_design_inv_naca0012.test_vals  = [-3.918503, -3.332494, 0.134359, 0.218097]
    tutorial_design_inv_naca0012.no_restart = True
    test_list.append(tutorial_design_inv_naca0012)

    # Turbulent RAE 2822 Design
    tutorial_design_turb_rae2822            = TestCase('design_turb_rae2822')
    tutorial_design_turb_rae2822.cfg_dir    = "../Tutorials/design/Turbulent_2D_Constrained_RAE2822"
    tutorial_design_turb_rae2822.cfg_file   = "turb_SA_RAE2822.cfg"
    tutorial_design_turb_rae2822.test_iter  = 0
    tutorial_design_turb_rae2822.test_vals  = [-1.700114, -4.941834, 0.218348, 0.190357]
    tutorial_design_turb_rae2822.no_restart = True
    test_list.append(tutorial_design_turb_rae2822)

    # Multi Objective Design
    tutorial_design_multiobj            = TestCase('design_multiobj')
    tutorial_design_multiobj.cfg_dir    = "../Tutorials/design/Multi_Objective_Shape_Design"
    tutorial_design_multiobj.cfg_file   = "inv_wedge_ROE_multiobj_combo.cfg"
    tutorial_design_multiobj.test_iter  = 0
    tutorial_design_multiobj.test_vals  = [2.657333, -3.020635, 370220.000000, 0.000000]
    tutorial_design_multiobj.no_restart = True
    test_list.append(tutorial_design_multiobj)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    # set suitable defaults unless something else has been specified
    # command: "mpirun -n 2 SU2_CFD"
    # timeout: 1600
    # tol:     0.00001
    for test in test_list:
        if test.command.empty():
            test.command = TestCase.Command("mpirun -n 2", "SU2_CFD")
        if test.timeout == 0:
            test.timeout = 1600
        if test.tol == 0.0:
            test.tol = 0.00001

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
