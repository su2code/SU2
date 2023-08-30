#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

    #####################################
    ### Disc. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    discadj_naca0012           = TestCase('discadj_naca0012')
    discadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    discadj_naca0012.cfg_file  = "inv_NACA0012_discadj.cfg"
    discadj_naca0012.test_iter = 100
    discadj_naca0012.test_vals = [-3.561506, -8.926634, -0.000000, 0.005587]
    test_list.append(discadj_naca0012)

    # Inviscid Cylinder 3D (multiple markers)
    discadj_cylinder3D           = TestCase('discadj_cylinder3D')
    discadj_cylinder3D.cfg_dir   = "disc_adj_euler/cylinder3D"
    discadj_cylinder3D.cfg_file  = "inv_cylinder3D.cfg"
    discadj_cylinder3D.test_iter = 5
    discadj_cylinder3D.test_vals = [-3.734502, -3.839637, 0.000000, 0.000000]
    test_list.append(discadj_cylinder3D)

    # Arina nozzle 2D
    discadj_arina2k              = TestCase('discadj_arina2k')
    discadj_arina2k.cfg_dir      = "disc_adj_euler/arina2k"
    discadj_arina2k.cfg_file     = "Arina2KRS.cfg"
    discadj_arina2k.test_iter    = 20
    discadj_arina2k.test_vals    = [-3.111181, -3.501516, 6.8705e-02, 0]
    test_list.append(discadj_arina2k)

    # Equivalent area NACA64-206
    ea_naca64206              = TestCase('ea_naca64206')
    ea_naca64206.cfg_dir      = "optimization_euler/equivalentarea_naca64206"
    ea_naca64206.cfg_file     = "NACA64206.cfg"
    ea_naca64206.test_iter    = 10
    ea_naca64206.test_vals    = [3.181093, 2.471539, -5487700.0, 8.3604]
    test_list.append(ea_naca64206)

    ####################################
    ### Disc. adj. compressible RANS ###
    ####################################

    # Adjoint turbulent NACA0012 SA
    discadj_rans_naca0012_sa           = TestCase('discadj_rans_naca0012_sa')
    discadj_rans_naca0012_sa.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    discadj_rans_naca0012_sa.test_iter = 10
    discadj_rans_naca0012_sa.test_vals = [-2.230568, 0.644202, 0.181590, -0.000018, 5.000000, -3.421717, 5.000000, -6.769530]
    test_list.append(discadj_rans_naca0012_sa)

    # Adjoint turbulent NACA0012 SST
    discadj_rans_naca0012_sst           = TestCase('discadj_rans_naca0012_sst')
    discadj_rans_naca0012_sst.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    discadj_rans_naca0012_sst.test_iter = 10
    discadj_rans_naca0012_sst.test_vals         = [-2.221848, -0.502430, 0.182000, -0.000018]
    discadj_rans_naca0012_sst.test_vals_aarch64 = [-2.221848, -0.502430, 0.182000, -0.000018]
    test_list.append(discadj_rans_naca0012_sst)

    #######################################
    ### Disc. adj. incompressible Euler ###
    #######################################

    # Adjoint Incompressible Inviscid NACA0012
    discadj_incomp_NACA0012           = TestCase('discadj_incomp_NACA0012')
    discadj_incomp_NACA0012.cfg_dir   = "disc_adj_incomp_euler/naca0012"
    discadj_incomp_NACA0012.cfg_file  = "incomp_NACA0012_disc.cfg"
    discadj_incomp_NACA0012.test_iter = 20
    discadj_incomp_NACA0012.test_vals = [20.000000, -4.095412, -2.690483, 0.000000]
    test_list.append(discadj_incomp_NACA0012)

    #####################################
    ### Disc. adj. incompressible N-S ###
    #####################################

    # Adjoint Incompressible Viscous Cylinder (Heated)
    discadj_incomp_cylinder           = TestCase('discadj_incomp_cylinder')
    discadj_incomp_cylinder.cfg_dir   = "disc_adj_incomp_navierstokes/cylinder"
    discadj_incomp_cylinder.cfg_file  = "heated_cylinder.cfg"
    discadj_incomp_cylinder.test_iter = 20
    discadj_incomp_cylinder.test_vals = [20.000000, -2.195581, -2.162081, 0.000000]
    test_list.append(discadj_incomp_cylinder)

    ######################################
    ### Disc. adj. incompressible RANS ###
    ######################################

    # Adjoint Incompressible Turbulent NACA 0012 SA
    discadj_incomp_turb_NACA0012_sa           = TestCase('discadj_incomp_turb_NACA0012_sa')
    discadj_incomp_turb_NACA0012_sa.cfg_dir   = "disc_adj_incomp_rans/naca0012"
    discadj_incomp_turb_NACA0012_sa.cfg_file  = "turb_naca0012_sa.cfg"
    discadj_incomp_turb_NACA0012_sa.test_iter = 10
    discadj_incomp_turb_NACA0012_sa.test_vals         = [10.000000, -3.846018, -1.031079, 0.000000]
    discadj_incomp_turb_NACA0012_sa.test_vals_aarch64 = [10.000000, -3.846020, -1.031078, 0.000000]
    test_list.append(discadj_incomp_turb_NACA0012_sa)

    # Adjoint Incompressible Turbulent NACA 0012 SST
    discadj_incomp_turb_NACA0012_sst           = TestCase('discadj_incomp_turb_NACA0012_sst')
    discadj_incomp_turb_NACA0012_sst.cfg_dir   = "disc_adj_incomp_rans/naca0012"
    discadj_incomp_turb_NACA0012_sst.cfg_file  = "turb_naca0012_sst.cfg"
    discadj_incomp_turb_NACA0012_sst.test_iter = 10
    discadj_incomp_turb_NACA0012_sst.test_vals         = [-4.183391, -2.287284, -7.892044, 0, -1.230919]
    discadj_incomp_turb_NACA0012_sst.test_vals_aarch64 = [-4.183391, -2.287284, -7.892044, 0, -1.230919]
    test_list.append(discadj_incomp_turb_NACA0012_sst)

    ####################################################################
    ###  Disc. Adj. Axisymmetric RANS                                ###
    ####################################################################

    # Adjoint Axisymmetric RANS
    discadj_axisymmetric_rans_nozzle            = TestCase('discadj_axisymmetric_rans')
    discadj_axisymmetric_rans_nozzle.cfg_dir    = "axisymmetric_rans/air_nozzle"
    discadj_axisymmetric_rans_nozzle.cfg_file   = "air_nozzle_restart.cfg"
    discadj_axisymmetric_rans_nozzle.test_iter  = 10
    discadj_axisymmetric_rans_nozzle.test_vals  = [9.516464, 5.022130, 7.370705, 2.744551]
    discadj_axisymmetric_rans_nozzle.no_restart = True
    test_list.append(discadj_axisymmetric_rans_nozzle)

    #######################################################
    ### Unsteady Disc. adj. compressible RANS           ###
    #######################################################

    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder.cfg"
    discadj_cylinder.test_iter = 9
    discadj_cylinder.test_vals = [3.746909, -1.544883, -0.008321, 0.000014] #last 4 columns
    discadj_cylinder.unsteady  = True
    test_list.append(discadj_cylinder)

    ##############################################################
    ### Unsteady Disc. adj. compressible RANS Windowed Average ###
    ##############################################################

    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder_windowed_average_AD')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder_Windowing_AD.cfg"
    discadj_cylinder.test_iter = 9
    discadj_cylinder.test_vals = [3.004406] #last column
    discadj_cylinder.unsteady  = True
    test_list.append(discadj_cylinder)

    ##############################################################
    ### Unsteady Disc. adj. compressible RANS Windowed Average ###
    ##############################################################

    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder_windowed_average')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder_Windowing.cfg"
    discadj_cylinder.test_iter = 6
    discadj_cylinder.test_vals = [0.202349, -0.000119, 1.899933, -0.000050, 1.067900]
    discadj_cylinder.tol       = 0.0001
    discadj_cylinder.command   = TestCase.Command("mpirun -n 2", "SU2_CFD_DIRECTDIFF")
    discadj_cylinder.unsteady  = True
    test_list.append(discadj_cylinder)

    ##########################################################################
    ### Unsteady Disc. adj. compressible RANS DualTimeStepping 1st order   ###
    ##########################################################################

    # Turbulent Cylinder
    discadj_DT_1ST_cylinder           = TestCase('unsteady_cylinder_DT_1ST')
    discadj_DT_1ST_cylinder.cfg_dir   = "disc_adj_rans/cylinder_DT_1ST"
    discadj_DT_1ST_cylinder.cfg_file  = "cylinder.cfg"
    discadj_DT_1ST_cylinder.test_iter = 9
    discadj_DT_1ST_cylinder.test_vals = [3.698168, -1.607050, -0.002159, 0.000028] #last 4 columns
    discadj_DT_1ST_cylinder.unsteady  = True
    test_list.append(discadj_DT_1ST_cylinder)

    ######################################################
    ### Unsteady Disc. adj. compressible pitching NACA ###
    ######################################################

    # compressible pitching NACA0012
    discadj_pitchingNACA0012           = TestCase('pitchingNACA0012')
    discadj_pitchingNACA0012.cfg_dir   = "disc_adj_euler/naca0012_pitching"
    discadj_pitchingNACA0012.cfg_file  = "inv_NACA0012_pitching.cfg"
    discadj_pitchingNACA0012.test_iter = 4
    discadj_pitchingNACA0012.test_vals = [-1.223480, -1.639387, -0.007591, 0.000013]
    discadj_pitchingNACA0012.unsteady  = True
    test_list.append(discadj_pitchingNACA0012)

    #######################################################
    ### Disc. adj. turbomachinery                       ###
    #######################################################

    # Transonic Stator 2D
    discadj_trans_stator           = TestCase('transonic_stator')
    discadj_trans_stator.cfg_dir   = "disc_adj_turbomachinery/transonic_stator_2D"
    discadj_trans_stator.cfg_file  = "transonic_stator.cfg"
    discadj_trans_stator.test_iter = 79
    discadj_trans_stator.test_vals = [79, 0.769972, 0.374328, 0.474458, -0.996522, 2.154079, -4.446062]
    discadj_trans_stator.test_vals_aarch64 = [79.000000, 0.769908, 0.374328, 0.474687, -0.996522, 2.154079, -4.446062]
    test_list.append(discadj_trans_stator)

    ###################################
    ### Structural Adjoint          ###
    ###################################

    # Structural model
    discadj_fea           = TestCase('discadj_fea')
    discadj_fea.cfg_dir   = "disc_adj_fea"
    discadj_fea.cfg_file  = "configAD_fem.cfg"
    discadj_fea.test_iter = 4
    discadj_fea.test_vals         = [-2.849453, -3.238429, -0.000364, -8.708700] #last 4 columns
    discadj_fea.test_vals_aarch64 = [-2.849573, -3.238516, -0.000364, -8.708700] #last 4 columns
    test_list.append(discadj_fea)

    ###################################
    ### Disc. adj. heat             ###
    ###################################

    # Discrete adjoint for heated cylinder
    discadj_heat           = TestCase('discadj_heat')
    discadj_heat.cfg_dir   = "disc_adj_heat"
    discadj_heat.cfg_file  = "disc_adj_heat.cfg"
    discadj_heat.test_iter = 10
    discadj_heat.test_vals         = [-2.226539, 0.605868, 0.000000, -6.256400]
    discadj_heat.test_vals_aarch64 = [-2.226539, 0.605868, 0.000000, -6.256400]
    test_list.append(discadj_heat)

    ###################################
    ### Coupled FSI Adjoint         ###
    ###################################

    # Legacy driver
    discadj_fsi           = TestCase('discadj_fsi')
    discadj_fsi.cfg_dir   = "disc_adj_fsi"
    discadj_fsi.cfg_file  = "config.cfg"
    discadj_fsi.test_iter = 6
    discadj_fsi.test_vals = [6.000000, -1.949946, -3.080711, 0.000440, -1.063100]
    test_list.append(discadj_fsi)

    # Multi physics framework
    discadj_fsi2           = TestCase('discadj_fsi_airfoil')
    discadj_fsi2.cfg_dir   = "disc_adj_fsi/Airfoil_2d"
    discadj_fsi2.cfg_file  = "config.cfg"
    discadj_fsi2.test_iter = 8
    discadj_fsi2.test_vals         = [-4.349377, 0.192713, -1.303589, 0.75407, 2.3244]
    discadj_fsi2.test_vals_aarch64 = [-4.349372, 0.190601, -1.303589, 0.754070, 2.324400]
    discadj_fsi2.tol       = 0.00001
    test_list.append(discadj_fsi2)

    ###################################
    ### Coupled CHT Adjoint         ###
    ###################################

    # Coupled discrete adjoint for heatflux in heated cylinder array
    discadj_cht           = TestCase('discadj_cht')
    discadj_cht.cfg_dir   = "coupled_cht/disc_adj_incomp_2d"
    discadj_cht.cfg_file  = "cht_2d_3cylinders.cfg"
    discadj_cht.test_iter = 10
    discadj_cht.test_vals = [-2.955510, -3.085549, -3.085516, -3.085511]
    test_list.append(discadj_cht)

    # 2D DA cht streamwise periodic case, 2 zones, avg temp objective
    da_sp_pinArray_cht_2d_dp_hf           = TestCase('da_sp_pinArray_cht_2d_dp_hf')
    da_sp_pinArray_cht_2d_dp_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    da_sp_pinArray_cht_2d_dp_hf.cfg_file  = "DA_configMaster.cfg"
    da_sp_pinArray_cht_2d_dp_hf.test_iter = 100
    da_sp_pinArray_cht_2d_dp_hf.test_vals = [-4.743218, -4.059800, -4.138231]
    da_sp_pinArray_cht_2d_dp_hf.multizone = True
    test_list.append(da_sp_pinArray_cht_2d_dp_hf)

    # 2D DA cht streamwise periodic case, 2 zones, PressureDrop objective, additional pressure drop adjoint equation
    da_sp_pinArray_cht_2d_mf           = TestCase('da_sp_pinArray_cht_2d_mf')
    da_sp_pinArray_cht_2d_mf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/dp-adjoint_chtPinArray_2d"
    da_sp_pinArray_cht_2d_mf.cfg_file  = "configMaster.cfg"
    da_sp_pinArray_cht_2d_mf.test_iter = 100
    da_sp_pinArray_cht_2d_mf.test_vals = [-4.611237, -1.290668, -1.453678, -18.507826, -0.859013, -5.755971, -19.067984, -47.993931]
    da_sp_pinArray_cht_2d_mf.multizone = True
    test_list.append(da_sp_pinArray_cht_2d_mf)

    # 2D unsteady CHT vortex shedding at RE=200. TAVG_Temperature OF
    da_unsteadyCHT_cylinder           = TestCase('da_unsteadyCHT_cylinder')
    da_unsteadyCHT_cylinder.cfg_dir   = "coupled_cht/disc_adj_unsteadyCHT_cylinder"
    da_unsteadyCHT_cylinder.cfg_file  = "chtMaster.cfg"
    da_unsteadyCHT_cylinder.test_iter = 2
    da_unsteadyCHT_cylinder.test_vals         = [-3.508906, -4.317739, -4.241558, -11.836892, -12.862650, 0.000000, 3.688000, 0.295190]
    da_unsteadyCHT_cylinder.test_vals_aarch64 = [-3.508906, -4.317739, -4.241558, -11.836892, -12.862650, 0.000000, 3.688000, 0.295190]
    da_unsteadyCHT_cylinder.unsteady  = True
    da_unsteadyCHT_cylinder.multizone = True
    test_list.append(da_unsteadyCHT_cylinder)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    # set suitable defaults unless something else has been specified
    # command: "mpirun -n 2 SU2_CFD_AD"
    # timeout: 1600
    # tol:     0.00001
    for test in test_list:
        if test.command.empty():
            test.command = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
        if test.timeout == 0:
            test.timeout = 1600
        if test.tol == 0.0:
            test.tol = 0.00001

    pass_list = [ test.run_test() for test in test_list ]

    ##################################
    ### Disc. adj. flamelet solver ###
    ##################################

    # 2D planar laminar premixed flame on isothermal burner (restart)
    discadj_flamelet_ch4_hx                  = TestCase('discadj_flamelet_ch4_hx')
    discadj_flamelet_ch4_hx.command          = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    discadj_flamelet_ch4_hx.cfg_dir          = "flamelet/02_laminar_premixed_ch4_flame_hx_ad"
    discadj_flamelet_ch4_hx.cfg_file         = "lam_prem_ch4_hx_ad.cfg"
    discadj_flamelet_ch4_hx.multizone        = False
    discadj_flamelet_ch4_hx.test_iter        = 10
    discadj_flamelet_ch4_hx.timeout          = 20000
    discadj_flamelet_ch4_hx.reference_file   = "restart_adj_custom.csv.ref"
    discadj_flamelet_ch4_hx.test_file        = "restart_adj_custom.csv"
    discadj_flamelet_ch4_hx.comp_threshold   = 1e-6
    discadj_flamelet_ch4_hx.tol_file_percent = 0.1
    pass_list.append(discadj_flamelet_ch4_hx.run_filediff())
    test_list.append(discadj_flamelet_ch4_hx)

    # 2D planar laminar premixed flame on isothermal burner with conjugate heat transfer (restart)
    discadj_flamelet_ch4_cht                  = TestCase('discadj_flamelet_ch4_cht')
    discadj_flamelet_ch4_cht.command          = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    discadj_flamelet_ch4_cht.cfg_dir          = "flamelet/04_laminar_premixed_ch4_flame_cht_ad"
    discadj_flamelet_ch4_cht.cfg_file         = "lam_prem_ch4_cht_ad_master.cfg"
    discadj_flamelet_ch4_cht.multizone        = True
    discadj_flamelet_ch4_cht.test_iter        = 5
    discadj_flamelet_ch4_cht.reference_file   = "restart_adj_T_0.csv.ref"
    discadj_flamelet_ch4_cht.test_file        = "restart_adj_T_0.csv"
    discadj_flamelet_ch4_cht.comp_threshold   = 1e-6
    discadj_flamelet_ch4_cht.tol_file_percent = 0.1
    discadj_flamelet_ch4_cht.timeout          = 20000
    pass_list.append(discadj_flamelet_ch4_cht.run_filediff())
    test_list.append(discadj_flamelet_ch4_cht)

    ################################################
    ### Gradient check (dot) for flamelet solver ###
    ################################################

    # 2D planar laminar premixed flame on isothermal burner (restart)
    # This test restarts on the output of test discadj_flamelet_ch4_hx and
    # will only pass if test discadj_flamelet_ch4_hx passes.
    dot_flamelet_ch4_hx                  = TestCase('dot_flamelet_ch4_hx')
    dot_flamelet_ch4_hx.cfg_dir          = "flamelet/02_laminar_premixed_ch4_flame_hx_ad"
    dot_flamelet_ch4_hx.cfg_file         = "lam_prem_ch4_hx_dot.cfg"
    dot_flamelet_ch4_hx.test_iter        = 10
    dot_flamelet_ch4_hx.command          = TestCase.Command("mpirun -n 2", "SU2_DOT_AD")
    dot_flamelet_ch4_hx.timeout          = 20000
    dot_flamelet_ch4_hx.reference_file   = "of_grad.csv.ref"
    dot_flamelet_ch4_hx.test_file        = "of_grad.csv"
    dot_flamelet_ch4_hx.comp_threshold   = 1e-6
    dot_flamelet_ch4_hx.tol_file_percent = 0.1
    pass_list.append(dot_flamelet_ch4_hx.run_filediff())
    test_list.append(dot_flamelet_ch4_hx)

    # 2D planar laminar premixed flame on isothermal burner with conjugate heat transfer (restart)
    # This test restarts on the output of test discadj_flamelet_ch4_cht and
    # will only pass if test discadj_flamelet_ch4_cht passes.
    dot_flamelet_ch4_cht                  = TestCase('dot_flamelet_ch4_cht')
    dot_flamelet_ch4_cht.cfg_dir          = "flamelet/04_laminar_premixed_ch4_flame_cht_ad"
    dot_flamelet_ch4_cht.cfg_file         = "lam_prem_ch4_cht_dot_master.cfg"
    dot_flamelet_ch4_cht.multizone        = True
    dot_flamelet_ch4_cht.test_iter        = 10
    dot_flamelet_ch4_cht.command          = TestCase.Command("mpirun -n 2", "SU2_DOT_AD")
    dot_flamelet_ch4_cht.timeout          = 20000
    dot_flamelet_ch4_cht.reference_file   = "of_grad.csv.ref"
    dot_flamelet_ch4_cht.test_file        = "of_grad.csv"
    dot_flamelet_ch4_cht.comp_threshold   = 1e-6
    dot_flamelet_ch4_cht.tol_file_percent = 0.1
    pass_list.append(dot_flamelet_ch4_cht.run_filediff())
    test_list.append(dot_flamelet_ch4_cht)

    ##################################################
    ### Structural Adjoint - Topology Optimization ###
    ##################################################

    # test discrete_adjoint.py
    discadj_topol_optim = TestCase('discadj_topol_optim')
    discadj_topol_optim.cfg_dir = "fea_topology"
    discadj_topol_optim.cfg_file  = "config.cfg"
    discadj_topol_optim.test_iter = 0
    discadj_topol_optim.command = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    discadj_topol_optim.timeout   = 1600
    discadj_topol_optim.reference_file         = "grad_ref_node.dat.ref"
    discadj_topol_optim.reference_file_aarch64 = "grad_ref_node_aarch64.dat.ref"
    discadj_topol_optim.test_file = "grad_ref_node.dat"
    pass_list.append(discadj_topol_optim.run_filediff())
    test_list.append(discadj_topol_optim)

    ####################################################################################
    ### Unsteady Disc. adj. compressible RANS Windowed Average with restart solution ###
    ####################################################################################

    # NACA0012 Airfoil
    unsteady_naca0012           = TestCase('unsteady_NACA0012_restart_adjoint')
    unsteady_naca0012.cfg_dir   = "disc_adj_rans/naca0012"
    unsteady_naca0012.cfg_file  = "naca0012.cfg"
    unsteady_naca0012.test_iter = 14
    unsteady_naca0012.command   = TestCase.Command(exec = "discrete_adjoint.py", param = "-f")
    unsteady_naca0012.timeout   = 1600
    unsteady_naca0012.reference_file = "of_grad_cd.csv.ref"
    unsteady_naca0012.test_file = "of_grad_cd.csv"
    unsteady_naca0012.unsteady  = True
    pass_list.append(unsteady_naca0012.run_filediff())
    test_list.append(unsteady_naca0012)

    ####################################################################################
    ### Unsteady Disc. adj. compressible RANS Windowed Average  only adjoint 		 ###
    ####################################################################################

    # NACA0012 Airfoil (Test depends on results of "unsteady_NACA0012_restart_adjoint")
    unsteady_naca0012           = TestCase('unsteady_NACA0012_adjoint_only')
    unsteady_naca0012.cfg_dir   = "disc_adj_rans/naca0012"
    unsteady_naca0012.cfg_file  = "naca0012.cfg"
    unsteady_naca0012.test_iter = 14
    unsteady_naca0012.command   = TestCase.Command(exec = "discrete_adjoint.py", param = "-m adj -f")
    unsteady_naca0012.timeout   = 1600
    unsteady_naca0012.reference_file = "of_grad_cd.csv.ref"
    unsteady_naca0012.test_file = "of_grad_cd.csv"
    unsteady_naca0012.unsteady  = True
    pass_list.append(unsteady_naca0012.run_filediff())
    test_list.append(unsteady_naca0012)

    ####################################################################
    ###  Python Wrapper                                              ###
    ####################################################################

    # FEA AD Flow Load Sensitivity
    pywrapper_FEA_AD_FlowLoad               = TestCase('pywrapper_FEA_AD_FlowLoad')
    pywrapper_FEA_AD_FlowLoad.cfg_dir       = "py_wrapper/disc_adj_fea/flow_load_sens"
    pywrapper_FEA_AD_FlowLoad.cfg_file      = "configAD_fem.cfg"
    pywrapper_FEA_AD_FlowLoad.test_iter     = 100
    pywrapper_FEA_AD_FlowLoad.test_vals     = [-0.13945587401785386, -0.5859858866132448, -0.00036377840086080694, -0.0031005670174756366] #last 4 columns
    pywrapper_FEA_AD_FlowLoad.command       = TestCase.Command("mpirun -n 2", "python", "run_adjoint.py --parallel -f")
    pywrapper_FEA_AD_FlowLoad.timeout       = 1600
    pywrapper_FEA_AD_FlowLoad.tol           = 0.000001
    pywrapper_FEA_AD_FlowLoad.new_output    = False
    test_list.append(pywrapper_FEA_AD_FlowLoad)
    pass_list.append(pywrapper_FEA_AD_FlowLoad.run_test())

    # FEA unsteady AD Load Sensitivity
    pywrapper_Unst_FEA_AD = TestCase('pywrapper_Unst_FEA_AD')
    pywrapper_Unst_FEA_AD.cfg_dir = "py_wrapper/custom_load_fea"
    pywrapper_Unst_FEA_AD.cfg_file = "config.cfg"
    pywrapper_Unst_FEA_AD.test_iter = 100
    pywrapper_Unst_FEA_AD.test_vals = [0.256684, 0.256684, 0.319877, 0.320149, -0.184491, -0.184509]
    pywrapper_Unst_FEA_AD.command = TestCase.Command("mpirun -n 2", "python", "run_ad.py")
    pywrapper_Unst_FEA_AD.timeout = 1600
    pywrapper_Unst_FEA_AD.tol = 0.00001
    pywrapper_Unst_FEA_AD.new_output = False
    test_list.append(pywrapper_Unst_FEA_AD)
    pass_list.append(pywrapper_Unst_FEA_AD.run_test())

    # Heat solver unsteady AD
    pywrapper_Unst_Heat_AD = TestCase('pywrapper_Unst_Heat_AD')
    pywrapper_Unst_Heat_AD.cfg_dir = "py_wrapper/custom_heat_flux"
    pywrapper_Unst_Heat_AD.cfg_file = "run_ad.py"
    pywrapper_Unst_Heat_AD.test_iter = 100
    pywrapper_Unst_Heat_AD.test_vals = [0.776365, 0.776430, 1.000003]
    pywrapper_Unst_Heat_AD.command = TestCase.Command("mpirun -n 2", "python", "run_ad.py")
    pywrapper_Unst_Heat_AD.timeout = 1600
    pywrapper_Unst_Heat_AD.tol = 0.00001
    pywrapper_Unst_Heat_AD.new_output = False
    test_list.append(pywrapper_Unst_Heat_AD)
    pass_list.append(pywrapper_Unst_Heat_AD.run_test())

    # Flow AD Mesh Displacement Sensitivity
    pywrapper_CFD_AD_MeshDisp               = TestCase('pywrapper_CFD_AD_MeshDisp')
    pywrapper_CFD_AD_MeshDisp.cfg_dir       = "py_wrapper/disc_adj_flow/mesh_disp_sens"
    pywrapper_CFD_AD_MeshDisp.cfg_file      = "configAD_flow.cfg"
    pywrapper_CFD_AD_MeshDisp.test_iter     = 1000
    pywrapper_CFD_AD_MeshDisp.test_vals     = [30.000000, -2.520972485907894, 1.3848377455328362, 0.000000] #last 4 columns
    pywrapper_CFD_AD_MeshDisp.command       = TestCase.Command("mpirun -n 2", "python", "run_adjoint.py --parallel -f")
    pywrapper_CFD_AD_MeshDisp.timeout       = 1600
    pywrapper_CFD_AD_MeshDisp.tol           = 0.000001
    pywrapper_CFD_AD_MeshDisp.new_output    = False
    test_list.append(pywrapper_CFD_AD_MeshDisp)
    pass_list.append(pywrapper_CFD_AD_MeshDisp.run_test())

    # Flow AD Mesh Displacement and Initial Coordinates Sensitivity
    pywrapper_wavy_wall_steady = TestCase('pywrapper_wavy_wall_steady')
    pywrapper_wavy_wall_steady.cfg_dir = "py_wrapper/wavy_wall"
    pywrapper_wavy_wall_steady.cfg_file = "run_steady.py"
    pywrapper_wavy_wall_steady.test_iter = 100
    pywrapper_wavy_wall_steady.test_vals = [-1.360044, 2.580709, -2.892473]
    pywrapper_wavy_wall_steady.command = TestCase.Command("mpirun -n 2", "python", "run_steady.py")
    pywrapper_wavy_wall_steady.timeout = 1600
    pywrapper_wavy_wall_steady.tol = 0.00001
    pywrapper_wavy_wall_steady.new_output = False
    test_list.append(pywrapper_wavy_wall_steady)
    pass_list.append(pywrapper_wavy_wall_steady.run_test())

    ####################################################################
    ###  Unsteady Disc. adj. compressible RANS restart optimization  ###
    ####################################################################

    # test shape_optimization.py
    naca_restart_shape_opt      = TestCase('restart_shape_optimization')
    naca_restart_shape_opt.cfg_dir    = "optimization_rans/naca0012"
    naca_restart_shape_opt.cfg_file   = "naca0012.cfg"
    naca_restart_shape_opt.test_iter  = 1
    naca_restart_shape_opt.test_vals  = [1.000000, 1.000000, 0.007046, 0.196671]
    naca_restart_shape_opt.command    = TestCase.Command(exec = "shape_optimization.py", param = "-f")
    naca_restart_shape_opt.timeout    = 1600
    naca_restart_shape_opt.tol       = 0.00001
    pass_list.append(naca_restart_shape_opt.run_opt())
    test_list.append(naca_restart_shape_opt)

    ####################################################################
    ###  Unsteady Disc. Adj. Coupled FSI                             ###
    ####################################################################

    # Unsteady multi physics framework
    dyn_discadj_fsi           = TestCase('dyn_discadj_fsi')
    dyn_discadj_fsi.cfg_dir   = "disc_adj_fsi/dyn_fsi"
    dyn_discadj_fsi.cfg_file  = "config.cfg"
    dyn_discadj_fsi.test_iter = 2
    dyn_discadj_fsi.command = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    dyn_discadj_fsi.timeout   = 1600
    dyn_discadj_fsi.reference_file = "grad_dv.opt.ref"
    dyn_discadj_fsi.reference_file_aarch64 = "grad_dv_aarch64.opt.ref"
    dyn_discadj_fsi.test_file = "grad_young.opt"
    dyn_discadj_fsi.unsteady  = True
    pass_list.append(dyn_discadj_fsi.run_filediff())
    test_list.append(dyn_discadj_fsi)

    ####################################################################
    ###  Sobolev Gradient Smoothing                                  ###
    ####################################################################

    # Sobolev gradient smoothing solver
    grad_smooth_oneram6           = TestCase('grad_smooth_sob')
    grad_smooth_oneram6.cfg_dir   = "grad_smooth/oneram6"
    grad_smooth_oneram6.cfg_file  = "ONERAM6_gradsmooth.cfg"
    grad_smooth_oneram6.test_iter = 2
    grad_smooth_oneram6.command   = TestCase.Command("mpirun -n 2", "SU2_DOT_AD")
    grad_smooth_oneram6.timeout   = 1600
    grad_smooth_oneram6.reference_file = "of_hess.dat.ref"
    grad_smooth_oneram6.reference_file_aarch64 = "of_hess_aarch64.dat.ref"
    grad_smooth_oneram6.test_file = "of_hess.dat"
    pass_list.append(grad_smooth_oneram6.run_filediff())
    test_list.append(grad_smooth_oneram6)


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
