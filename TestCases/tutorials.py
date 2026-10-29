#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 8.5.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
    cht_incompressible_unsteady.test_vals = [-3.204738, -0.080399, -0.080399, -0.080399, -11.531676, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 238.240000]
    cht_incompressible_unsteady.multizone = True
    cht_incompressible_unsteady.unsteady  = True
    test_list.append(cht_incompressible_unsteady)

    # CHT incompressible, 2D, 3 pins in crossflow
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "../Tutorials/multiphysics/steady_cht"
    cht_incompressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [-1.376345, -0.591212, -0.591212, -0.591212]
    cht_incompressible.command   = TestCase.Command(exec = "SU2_CFD")
    cht_incompressible.multizone = True
    test_list.append(cht_incompressible)

    # Solid-to-solid and solid-to-fluid CHT with contact resistance
    cht_CR           = TestCase('cht_solid_solid')
    cht_CR.cfg_dir   = "../Tutorials/multiphysics/contact_resistance_cht"
    cht_CR.cfg_file  = "master.cfg"
    cht_CR.test_iter = 80
    cht_CR.test_vals = [-8.606938, -9.227901, -10.410107, -2.129377]
    cht_CR.multizone = True
    test_list.append(cht_CR)

    ### Incompressible Flow

    # 2D pin case massflow periodic with heatflux BC and prescribed extracted outlet heat
    sp_pinArray_2d_mf_hf           = TestCase('sp_pinArray_2d_mf_hf')
    sp_pinArray_2d_mf_hf.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Streamwise_Periodic"
    sp_pinArray_2d_mf_hf.cfg_file  = "sp_pinArray_2d_mf_hf.cfg"
    sp_pinArray_2d_mf_hf.test_iter = 25
    sp_pinArray_2d_mf_hf.test_vals = [-4.674156, 0.512176, -0.756308, 241.767693]
    sp_pinArray_2d_mf_hf.test_vals_aarch64 = [-4.683630, -0.755570, 241.872160]
    test_list.append(sp_pinArray_2d_mf_hf)

    # 2D pin case pressure drop periodic with heatflux BC and temperature periodicity
    sp_pinArray_2d_dp_hf_tp           = TestCase('sp_pinArray_2d_dp_hf_tp')
    sp_pinArray_2d_dp_hf_tp.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Streamwise_Periodic"
    sp_pinArray_2d_dp_hf_tp.cfg_file  = "sp_pinArray_2d_dp_hf_tp.cfg"
    sp_pinArray_2d_dp_hf_tp.test_iter = 25
    sp_pinArray_2d_dp_hf_tp.test_vals = [-4.727583, 0.453505, -0.714476, 208.023676]
    sp_pinArray_2d_dp_hf_tp.test_vals_aarch64 = [-4.739709, 0.520033, -0.713547, 208.023676]
    test_list.append(sp_pinArray_2d_dp_hf_tp)

    # 90 degree pipe bend with wall functions from the experiments of Sudo et al.
    sudo_tutorial = TestCase('sudo_bend')
    sudo_tutorial.cfg_dir = "../Tutorials/incompressible_flow/Inc_Turbulent_Bend_Wallfunctions"
    sudo_tutorial.cfg_file = "sudo.cfg"
    sudo_tutorial.test_iter = 10
    sudo_tutorial.test_vals = [-14.007950, -12.597368, -13.189199, -12.791612, -12.928621, -9.436284, 15.000000, -1.574072]
    sudo_tutorial.command = TestCase.Command("mpirun -n 2", "SU2_CFD")
    test_list.append(sudo_tutorial)

    # design-primal: 90 degree pipe bend with wall functions from the experiments of Sudo et al.
    sudo_design_primal = TestCase('sudo_bend_design_primal')
    sudo_design_primal.cfg_dir = "../Tutorials/design/Inc_Turbulent_Bend_Wallfunctions"
    sudo_design_primal.cfg_file = "sudo_primal.cfg"
    sudo_design_primal.test_iter = 10
    sudo_design_primal.test_vals = [-12.284057, -10.882408, -11.358456, -10.866265, -11.009547, -7.753843, 64.578000]
    sudo_design_primal.command  = TestCase.Command("mpirun -n 2", "SU2_CFD")
    test_list.append(sudo_design_primal)

    # design-adjoint: 90 degree pipe bend with wall functions from the experiments of Sudo et al.
    sudo_design_adjoint = TestCase('sudo_bend_design_adjoint')
    sudo_design_adjoint.cfg_dir = "../Tutorials/design/Inc_Turbulent_Bend_Wallfunctions"
    sudo_design_adjoint.cfg_file = "sudo_adjoint.cfg"
    sudo_design_adjoint.test_iter = 10
    sudo_design_adjoint.test_vals = [-4.352938, -3.312788, -3.140720, -3.237621, -3.789619, -7.366722]
    sudo_design_adjoint.command  = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    test_list.append(sudo_design_adjoint)

    # Laminar vortex shedding behind a cylinder (Re=120)
    von_karman_cylinder = TestCase('von_karman_cylinder')
    von_karman_cylinder.cfg_dir = "../Tutorials/incompressible_flow/Inc_Von_Karman_Cylinder"
    von_karman_cylinder.cfg_file  = "unsteady_incomp_cylinder.cfg"
    von_karman_cylinder.test_iter = 10
    von_karman_cylinder.test_vals = [-7.845765, -7.681042, -8.736704, -0.002581, 1.423652]
    test_list.append(von_karman_cylinder)

    # Lid Driven Cavity Flow (Re=400), pressure-based. Disabled: the config and mesh live in
    # su2code/Tutorials#86, which had not landed as of this commit - ../Tutorials/incompressible_flow/
    # Inc_Lid_Driven_Cavity/incomp_pb_liddrivencavity.cfg does not exist yet, and enabling this
    # entry would fail CI outright rather than fail the regression check. Re-enable once that PR
    # merges (test_vals recorded below were measured against the config in that PR).
    # lid_driven_cavity = TestCase('lid_driven_cavity')
    # lid_driven_cavity.cfg_dir = "../Tutorials/incompressible_flow/Inc_Lid_Driven_Cavity"
    # lid_driven_cavity.cfg_file  = "incomp_pb_liddrivencavity.cfg"
    # lid_driven_cavity.test_iter = 2
    # lid_driven_cavity.test_vals = [-5.798956, -4.393539, -5.070494]
    # test_list.append(lid_driven_cavity)


    ### Species Transport

    # 3 species (2 eq) primitive venturi mixing
    species3_primitiveVenturi           = TestCase('species3_primitiveVenturi')
    species3_primitiveVenturi.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Species_Transport"
    species3_primitiveVenturi.cfg_file  = "species3_primitiveVenturi.cfg"
    species3_primitiveVenturi.test_iter = 50
    species3_primitiveVenturi.test_vals = [-5.679419, -4.675598, -4.699798, -5.463760, -1.072282, -5.965552, -6.075164, 5.000000, -0.558583, 5.000000, -2.514945, 5.000000, -0.514439, 1.660063, 0.502195, 0.603290, 0.554578]
    test_list.append(species3_primitiveVenturi)

    # 3 species (2 eq) primitive venturi mixing
    DAspecies3_primitiveVenturi           = TestCase('DAspecies3_primitiveVenturi')
    DAspecies3_primitiveVenturi.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Species_Transport"
    DAspecies3_primitiveVenturi.cfg_file  = "DAspecies3_primitiveVenturi.cfg"
    DAspecies3_primitiveVenturi.test_iter = 50
    DAspecies3_primitiveVenturi.test_vals = [-9.819107, -8.643459, -8.676921, -8.347339, -12.926241, -9.739487, -8.947991]
    DAspecies3_primitiveVenturi.command   = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    test_list.append(DAspecies3_primitiveVenturi)

    # 2 species (1 eq) kenics static mixer for composition-dependent model
    kenics_mixer_tutorial           = TestCase('kenics_mixer_tutorial')
    kenics_mixer_tutorial.cfg_dir   = "../Tutorials/incompressible_flow/Inc_Species_Transport_Composition_Dependent_Model"
    kenics_mixer_tutorial.cfg_file  = "kenics_mixer_tutorial.cfg"
    kenics_mixer_tutorial.test_iter = 10
    kenics_mixer_tutorial.test_vals = [-7.490448, -6.823978, -6.838620, -6.383873, -7.878123, -3.003654, -7.452198, 5.000000, -1.857226, 4.000000, -5.336949, 3.000000, -6.363373, 0.025668, 0.000000, 0.025668, 0.000000, 62.717000, 8.462600, 46.725000, 7.529400]
    kenics_mixer_tutorial.command   = TestCase.Command("mpirun -n 2", "SU2_CFD")
    test_list.append(kenics_mixer_tutorial)


    ### Incompressible Combustion

    # Pre-mixed, laminar hydrogen flame with heat loss
    premixed_hydrogen = TestCase('premixed_hydrogen')
    premixed_hydrogen.cfg_dir = "../Tutorials/incompressible_flow/Inc_Combustion/1__premixed_hydrogen"
    premixed_hydrogen.cfg_file = "H2_burner.cfg"
    premixed_hydrogen.test_iter = 10
    premixed_hydrogen.test_vals = [-8.771332, -9.721253, -11.354950, -4.381656, -12.832815]
    test_list.append(premixed_hydrogen)

    ### Compressible Flow

    # Inviscid Bump
    tutorial_inv_bump            = TestCase('inviscid_bump_tutorial')
    tutorial_inv_bump.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_Bump"
    tutorial_inv_bump.cfg_file   = "inv_channel.cfg"
    tutorial_inv_bump.test_iter  = 0
    tutorial_inv_bump.test_vals = [-1.437425, 4.075857, 0.085674, 0.123798]
    test_list.append(tutorial_inv_bump)

    # Inviscid Wedge
    tutorial_inv_wedge            = TestCase('inviscid_wedge_tutorial')
    tutorial_inv_wedge.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_Wedge"
    tutorial_inv_wedge.cfg_file   = "inv_wedge_HLLC.cfg"
    tutorial_inv_wedge.test_iter  = 0
    tutorial_inv_wedge.test_vals = [-0.481460, 5.253008, -0.191841, 0.033699]
    tutorial_inv_wedge.no_restart = True
    test_list.append(tutorial_inv_wedge)

    # Inviscid ONERA M6
    tutorial_inv_onera            = TestCase('inviscid_onera_tutorial')
    tutorial_inv_onera.cfg_dir    = "../Tutorials/compressible_flow/Inviscid_ONERAM6"
    tutorial_inv_onera.cfg_file   = "inv_ONERAM6.cfg"
    tutorial_inv_onera.test_iter  = 0
    tutorial_inv_onera.test_vals = [-5.204928, -4.597762, 0.243530, 0.111160]
    tutorial_inv_onera.no_restart = True
    test_list.append(tutorial_inv_onera)

    # Laminar Cylinder
    tutorial_lam_cylinder            = TestCase('laminar_cylinder_tutorial')
    tutorial_lam_cylinder.cfg_dir    = "../Tutorials/compressible_flow/Laminar_Cylinder"
    tutorial_lam_cylinder.cfg_file   = "lam_cylinder.cfg"
    tutorial_lam_cylinder.test_iter  = 0
    tutorial_lam_cylinder.test_vals = [-6.162141, -0.699617, -0.080602, 42.990503]
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
    tutorial_turb_flatplate.test_vals  = [-2.258584, -4.901015, -0.429401, 0.201034]
    tutorial_turb_flatplate.no_restart = True
    test_list.append(tutorial_turb_flatplate)

    # Transitional FlatPlate
    tutorial_trans_flatplate            = TestCase('transitional_flatplate_tutorial')
    tutorial_trans_flatplate.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Flat_Plate"
    tutorial_trans_flatplate.cfg_file   = "transitional_BC_model_ConfigFile.cfg"
    tutorial_trans_flatplate.test_iter  = 0
    tutorial_trans_flatplate.test_vals  = [-22.025101, -15.330766, 0.000000, 0.023944]
    tutorial_trans_flatplate.no_restart = True
    test_list.append(tutorial_trans_flatplate)

    # Transitional FlatPlate T3A
    tutorial_trans_flatplate_T3A            = TestCase('transitional_flatplate_tutorial_T3A')
    tutorial_trans_flatplate_T3A.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Flat_Plate/Langtry_and_Menter/T3A"
    tutorial_trans_flatplate_T3A.cfg_file   = "transitional_LM_model_ConfigFile.cfg"
    tutorial_trans_flatplate_T3A.test_iter  = 20
    tutorial_trans_flatplate_T3A.test_vals  = [-5.790325, -2.054834, -3.896287, -0.255197, -1.746976, 5.119339, -3.493229, 0.393180]
    tutorial_trans_flatplate_T3A.test_vals_aarch64 = [-5.808996, -2.070606, -3.969765, -0.277943, -1.953289, 1.708472, -3.514943, 0.357411]
    tutorial_trans_flatplate_T3A.no_restart = True
    test_list.append(tutorial_trans_flatplate_T3A)

    # Transitional FlatPlate T3Am
    tutorial_trans_flatplate_T3Am            = TestCase('transitional_flatplate_tutorial_T3Am')
    tutorial_trans_flatplate_T3Am.cfg_dir    = "../Tutorials/compressible_flow/Transitional_Flat_Plate/Langtry_and_Menter/T3A-"
    tutorial_trans_flatplate_T3Am.cfg_file   = "transitional_LM_model_ConfigFile.cfg"
    tutorial_trans_flatplate_T3Am.test_iter  = 20
    tutorial_trans_flatplate_T3Am.test_vals  = [-5.589599, -1.700867, -3.097326, -0.104799, -3.750523, 3.287643, -2.394576, 1.119623]
