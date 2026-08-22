#!/usr/bin/env python

## \file hybrid_regression.py
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
from TestCase import parse_args

def main():
    '''This program runs SU2 and ensures that the output matches specified values.
       This will be used to do checks when code is pushed to github
       to make sure nothing is broken. '''

    args = parse_args('Hybrid Regression Tests')

    test_list = []
    file_diff_list = []

    ##########################
    ### Compressible Euler ###
    ##########################

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 20
    channel.test_vals = [-2.000214, 3.536937, 0.034200, 0.193725]
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.387446, -3.845358, 0.296291, 0.025143]
    test_list.append(naca0012)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-3.690258, 2.034169, -0.249531, 0.043953]
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [0.280803, 0.008625]
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-4.010006, 1.521019, 0.300952, 0.019472]
    test_list.append(fixedCL_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.666182, 7.055173, -0.000631, 3.917770]
    test_list.append(bluntbody)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 100
    flatplate.test_vals = [-6.543281, -1.065152, 0.001196, 0.029390, 2.361500, -2.332100, 0.000000, 0.000000]
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-8.500691, -3.024888, 0.037975, 1.664042, 0.000000]
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals         = [-6.476618, -1.014896, 0.594402, 70.856953, 0.000000]
    cylinder_lowmach.test_vals_aarch64 = [-6.830996, -1.368850, -0.143956, 73.963354, 0]
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-7.595535, -6.993332, 0.007670, 13.591507, 0]
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals         = [-12.005126, -7.582352, -0.000000, 2.089953]
    poiseuille_profile.test_vals_aarch64 = [-12.009012, -7.262530, -0.000000, 2.089953]
    test_list.append(poiseuille_profile)

    # 2D Rotational Periodic
    periodic2d           = TestCase('periodic2d')
    periodic2d.cfg_dir   = "navierstokes/periodic2D"
    periodic2d.cfg_file  = "config.cfg"
    periodic2d.test_iter = 1400
    periodic2d.test_vals = [-10.817608, -8.363541, -8.287458, -5.334101, -1.088411, -2945.200000]
    test_list.append(periodic2d)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.190528, -5.335496, 0.383385, 0.077605, 0.000000]
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-1.028081, 5.870628, 0.370417, 0.074982, 0.000000]
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.487683, 5.870614, 0.370417, 0.074982]
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.958328, -7.438303, -0.187472, 0.015052]
    test_list.append(turb_flatplate)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.418706, -6.631575, 0.238586, 0.159599, 0.000000]
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 5
    turb_naca0012_sa.test_vals = [-12.038067, -16.332088, 1.080346, 0.018385, 20.000000, -2.873679, 0.000000, -14.250269, 0.000000]
    turb_naca0012_sa.test_vals_aarch64 = [-12.038091, -16.332090, 1.080346, 0.018385, 20.000000, -2.873236, 0.000000, -14.250271, 0.000000]
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-12.093958, -15.250719, -5.906323, 1.070413, 0.015775, -2.855776, 0.000000]
    turb_naca0012_sst.test_vals_aarch64 = [-12.075928, -15.246732, -5.861249, 1.070036, 0.015841, -2.835263, 0]
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals = [-12.080818, -14.837175, -5.732906, 1.000893, 0.019109, -2.119717]
    turb_naca0012_sst_sust.test_vals_aarch64 = [-12.073210, -14.836724, -5.732627, 1.000050, 0.019144, -2.629689]
    test_list.append(turb_naca0012_sst_sust)

    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [-5.192390, -10.448218, 0.773965, 1.022535, 0.040529, -2.383435]
    test_list.append(turb_naca0012_sst_fixedvalues)

    # NACA0012 (SST, explicit Euler for flow and turbulence equations)
    turb_naca0012_sst_expliciteuler           = TestCase('turb_naca0012_sst_expliciteuler')
    turb_naca0012_sst_expliciteuler.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_expliciteuler.cfg_file  = "turb_NACA0012_sst_expliciteuler.cfg"
    turb_naca0012_sst_expliciteuler.test_iter = 10
    turb_naca0012_sst_expliciteuler.test_vals = [-3.532365, -3.157224, 3.743381, 1.124798, 0.501715, -16.000000]
    test_list.append(turb_naca0012_sst_expliciteuler)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389718, -8.410477, 0.000048, 0.056343]
    test_list.append(propeller)

    #######################################
    ### Axisymmetric Compressible RANS  ###
    #######################################

    # Axisymmetric air nozzle (transonic)
    axi_rans_air_nozzle_restart           = TestCase('axi_rans_air_nozzle_restart')
    axi_rans_air_nozzle_restart.cfg_dir   = "axisymmetric_rans/air_nozzle"
    axi_rans_air_nozzle_restart.cfg_file  = "air_nozzle_restart.cfg"
    axi_rans_air_nozzle_restart.test_iter = 10
    axi_rans_air_nozzle_restart.test_vals = [-11.083070, -5.374685, -8.880088, -4.073510, 0.000000]
    axi_rans_air_nozzle_restart.test_vals_aarch64 = [-14.140441, -9.154674, -10.886121, -5.806594, 0.000000]
    test_list.append(axi_rans_air_nozzle_restart)

    #################################
    ## Compressible RANS Restart  ###
    #################################

    # NACA0012 SST Multigrid restart
    turb_naca0012_sst_restart_mg           = TestCase('turb_naca0012_sst_restart_mg')
    turb_naca0012_sst_restart_mg.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_restart_mg.cfg_file  = "turb_NACA0012_sst_multigrid_restart.cfg"
    turb_naca0012_sst_restart_mg.test_iter = 20
    turb_naca0012_sst_restart_mg.ntest_vals = 5
    turb_naca0012_sst_restart_mg.test_vals = [-6.589092, -5.057151, 0.830239, -0.008809, 0.078148]
    test_list.append(turb_naca0012_sst_restart_mg)

    #############################
    ### Compressibele RANS UQ ###
    #############################

    # NACA0012 1c
    turb_naca0012_1c           = TestCase('turb_naca0012_1c')
    turb_naca0012_1c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_1c.cfg_file  = "turb_NACA0012_uq_1c.cfg"
    turb_naca0012_1c.test_iter = 10
    turb_naca0012_1c.test_vals = [-4.978722, 1.343530, 0.443909, -0.029125]
    turb_naca0012_1c.test_vals_aarch64 = [-4.976620, 1.345983, 0.433171, -0.033685]
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals = [-5.482844, 1.260870, 0.404377, -0.040361]
    turb_naca0012_2c.test_vals_aarch64 = [-5.485484, 1.263406, 0.411442, -0.040859]
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals         = [-5.583730, 1.228732, 0.381968, -0.046233]
    turb_naca0012_3c.test_vals_aarch64 = [-5.583737, 1.232005, 0.390258, -0.046305]
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals = [-5.134179, 1.283465, 0.553971, 0.011916]
    turb_naca0012_p1c1.test_vals_aarch64 = [-5.114189, 1.285037, 0.406851, -0.043003]
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals         = [-5.553990, 1.234030, 0.424168, -0.033501]
    turb_naca0012_p1c2.test_vals_aarch64 = [-5.548245, 1.236384, 0.381821, -0.050337]
    test_list.append(turb_naca0012_p1c2)

    ######################################
    ### Harmonic Balance               ###
    ######################################

    # Description of the regression test
    harmonic_balance           = TestCase('harmonic_balance')
    harmonic_balance.cfg_dir   = "harmonic_balance"
    harmonic_balance.cfg_file  = "HB.cfg"
    harmonic_balance.test_iter = 25
    harmonic_balance.test_vals = [-1.559187, 0.829574, 0.931511, 3.954440]
    test_list.append(harmonic_balance)

    # Turbulent pitching NACA 64a010 airfoil
    hb_rans_preconditioning           = TestCase('hb_rans_preconditioning')
    hb_rans_preconditioning.cfg_dir   = "harmonic_balance/hb_rans_preconditioning"
    hb_rans_preconditioning.cfg_file  = "davis.cfg"
    hb_rans_preconditioning.test_iter = 25
    hb_rans_preconditioning.test_vals = [-1.905220, 0.481910, 0.598990, 3.605349, -5.945852]
    test_list.append(hb_rans_preconditioning)

    #############################
    ### Incompressible Euler  ###
    #############################

    # NACA0012 Hydrofoil
    inc_euler_naca0012           = TestCase('inc_euler_naca0012')
    inc_euler_naca0012.cfg_dir   = "incomp_euler/naca0012"
    inc_euler_naca0012.cfg_file  = "incomp_NACA0012.cfg"
    inc_euler_naca0012.test_iter = 20
    inc_euler_naca0012.test_vals = [-5.858104, -4.937295, 0.519817, 0.008958]
    test_list.append(inc_euler_naca0012)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals = [-6.171116, -5.413059, 0.008862, 0.127559]
    test_list.append(inc_nozzle)

    #############################
    ### Incompressible N-S    ###
    #############################

    # Laminar cylinder
    inc_lam_cylinder          = TestCase('inc_lam_cylinder')
    inc_lam_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_lam_cylinder.cfg_file  = "incomp_cylinder.cfg"
    inc_lam_cylinder.test_iter = 10
    inc_lam_cylinder.test_vals = [-4.161196, -3.573053, 0.025533, 4.944646]
    test_list.append(inc_lam_cylinder)

    # Buoyancy-driven cavity
    inc_buoyancy          = TestCase('inc_buoyancy')
    inc_buoyancy.cfg_dir   = "incomp_navierstokes/buoyancy_cavity"
    inc_buoyancy.cfg_file  = "lam_buoyancy_cavity.cfg"
    inc_buoyancy.test_iter = 20
    inc_buoyancy.test_vals = [-3.860465, -3.122402, 1.322722, -23.413000]
    test_list.append(inc_buoyancy)

    # Laminar heated cylinder with polynomial fluid model
    inc_poly_cylinder          = TestCase('inc_poly_cylinder')
    inc_poly_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_poly_cylinder.cfg_file  = "poly_cylinder.cfg"
    inc_poly_cylinder.test_iter = 20
    inc_poly_cylinder.test_vals         = [-8.230355, -2.426276, 0.018487, 1.902089, -173.000000]
    inc_poly_cylinder.test_vals_aarch64 = [-8.260165, -2.445453, 0.027209, 1.915447, -171.620000]
    test_list.append(inc_poly_cylinder)

    # X-coarse laminar bend as a mixed element CGNS test
    inc_lam_bend          = TestCase('inc_lam_bend')
    inc_lam_bend.cfg_dir   = "incomp_navierstokes/bend"
    inc_lam_bend.cfg_file  = "lam_bend.cfg"
    inc_lam_bend.test_iter = 10
    inc_lam_bend.test_vals = [-3.558795, -3.051040, -0.014118, 1.101254]
    test_list.append(inc_lam_bend)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012, SA
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.758114, -10.974548, -0.000004, -0.028637, 5.000000, -4.080155, 2.000000, -4.491724]
    test_list.append(inc_turb_naca0012)

    # NACA0012, SST_SUST
    inc_turb_naca0012_sst_sust           = TestCase('inc_turb_naca0012_sst_sust')
    inc_turb_naca0012_sst_sust.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_sst_sust.cfg_file  = "naca0012_SST_SUST.cfg"
    inc_turb_naca0012_sst_sust.test_iter = 20
    inc_turb_naca0012_sst_sust.test_vals = [-7.170018, 0.332638, 0.000002, 0.312117]
    test_list.append(inc_turb_naca0012_sst_sust)

    # Weakly coupled heat equation
    inc_weakly_coupled = TestCase('inc_weakly_coupled')
    inc_weakly_coupled.cfg_dir = "disc_adj_heat"
    inc_weakly_coupled.cfg_file = "primal.cfg"
    inc_weakly_coupled.test_iter = 10
    inc_weakly_coupled.test_vals = [-18.106270, -16.302974, -16.484882, -15.006575, -17.858050, -14.024836, 5.609100]
    test_list.append(inc_weakly_coupled)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-7.927912, -2.479033, 0.012837, 0.005130]
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals         = [-7.533969, -2.066689, 1.832252, 1.843016]
    spinning_cylinder.test_vals_aarch64 = [-8.008023, -2.611064, 1.497308, 1.487483]
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-2.560678, -1.175981, 0.062204, 1.399351, 2.219223, 1.399297, 2.217470, 0.000000]
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977498, 3.481817, -0.010188, -0.004302]
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Cosine gust in z-direction
    cosine_gust           = TestCase('cosine_gust_zdir')
    cosine_gust.cfg_dir   = "gust"
    cosine_gust.cfg_file  = "cosine_gust_zdir.cfg"
    cosine_gust.test_iter = 79
    cosine_gust.test_vals = [-2.418805, 0.002212, -0.001249, 0.000439, -0.000592]
    cosine_gust.unsteady  = True
    cosine_gust.enabled_with_tsan = False
    test_list.append(cosine_gust)

    # Gust with mesh deformation
    gust_mesh_defo           = TestCase('gust_with_mesh_deformation')
    gust_mesh_defo.cfg_dir   = "gust"
    gust_mesh_defo.cfg_file  = "gust_with_mesh_deformation.cfg"
    gust_mesh_defo.test_iter = 6
    gust_mesh_defo.test_vals = [-1.844761, 0.001100, -0.000310]
    gust_mesh_defo.unsteady  = True
    gust_mesh_defo.enabled_with_tsan = False
    test_list.append(gust_mesh_defo)

    # Aeroelastic
    aeroelastic           = TestCase('aeroelastic')
    aeroelastic.cfg_dir   = "aeroelastic"
    aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    aeroelastic.test_iter = 2
    aeroelastic.test_vals = [-1.876631, 4.021073, 0.081270, 0.027580, -0.001642, -0.000127, -1.052233]
    aeroelastic.unsteady  = True
    aeroelastic.enabled_on_cpu_arch = ["x86_64"] # Requires AVX-capable architecture
    test_list.append(aeroelastic)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714713, -5.763300, -0.214960, 0.023758, 0.000000]
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.008640, -6.889022, 1.435181, 0.433567]
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # unsteady pitching NACA0012, Euler, Deforming
    unst_deforming_naca0012           = TestCase('unst_deforming_naca0012')
    unst_deforming_naca0012.cfg_dir   = "disc_adj_euler/naca0012_pitching_def"
    unst_deforming_naca0012.cfg_file  = "inv_NACA0012_pitching_deform.cfg"
    unst_deforming_naca0012.test_iter = 5
    unst_deforming_naca0012.test_vals = [-3.665299, -3.794240, -3.717018, -3.148592]
    unst_deforming_naca0012.unsteady  = True
    unst_deforming_naca0012.enabled_with_tsan = False
    test_list.append(unst_deforming_naca0012)

    ######################################
    ### NICFD                          ###
    ######################################

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 30
    edge_VW.test_vals = [-7.124287, -0.922783, -0.000009, 0.000000]
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 20
    edge_PPR.test_vals         = [-12.028835, -5.863587, -0.000034, 0.000000]
    edge_PPR.test_vals_aarch64 = [ -7.139211, -0.980821, -0.000034, 0.000000]
    test_list.append(edge_PPR)

    ######################################
    ### Turbomachinery                 ###
    ######################################

	# Jones APU Turbocharger restart
    Jones_tc_restart           = TestCase('jones_turbocharger_restart')
    Jones_tc_restart.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_restart.cfg_file  = "Jones_restart.cfg"
    Jones_tc_restart.test_iter = 5
    Jones_tc_restart.test_vals = [-11.907788, -12.215989, -19.151836, -13.452316, -19.083454, -13.445506, 73286.000000, 73286.000000, 0.020056, 82.286000]
    test_list.append(Jones_tc_restart)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals = [1.167512, 1.598494, -2.928579, 2.573642, -2.527390, 3.016171, 106370.000000, 106370.000000, 5.726800, 64.383000]
    test_list.append(axial_stage2D)

    # 2D transonic stator restart
    transonic_stator_restart           = TestCase('transonic_stator_restart')
    transonic_stator_restart.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_restart.cfg_file  = "transonic_stator_restart.cfg"
    transonic_stator_restart.test_iter = 20
    transonic_stator_restart.test_vals         = [-4.367141, -2.487739, -2.079071, 1.728176, -1.464951, 3.225028, -471620.000000, 94.839000, -0.051125]
    transonic_stator_restart.test_vals_aarch64 = [-4.442510, -2.561369, -2.165778, 1.652750, -1.355494, 3.172712, -471620.000000, 94.843000, -0.043825]
    test_list.append(transonic_stator_restart)

    # Multiple turbomachinery interface restart
    multi_interface                    = TestCase('multi_interface')
    multi_interface.cfg_dir            = "turbomachinery/multi_interface"
    multi_interface.cfg_file           = "multi_interface_rst.cfg"
    multi_interface.test_iter          = 5
    multi_interface.test_vals          = [-8.632242, -8.894741, -9.348730]
    multi_interface.test_vals_aarch64  = [-8.632229, -8.894737, -9.348730]
    test_list.append(multi_interface)

    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 5
    uniform_flow.test_vals = [5.000000, 0.000000, -0.195002, -10.624446]
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow)

    # Channel_2D, native SU2 binary mesh format (.su2b)
    # channel_2D_WA.cfg loads channel_2D_su2bin.su2b directly, so SU2_DEF must
    # first convert channel_2D.su2 (3 zones) into that binary mesh.
    channel_2D_su2bin_convert           = TestCase('channel_2D_su2bin_convert')
    channel_2D_su2bin_convert.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D_su2bin_convert.cfg_file  = "mesh_su2_to_su2bin.cfg"
    channel_2D_su2bin_convert.command   = TestCase.Command(exec = "SU2_DEF")
    channel_2D_su2bin_convert.timeout   = 600
    test_list.append(channel_2D_su2bin_convert)

    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 2
    channel_2D.test_vals = [2.000000, 0.000000, 0.466194, 0.350086, 0.398979]
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 2
    channel_3D.test_vals = [2.000000, 0.000000, 0.632251, 0.534200, 0.431953]
    channel_3D.test_vals_aarch64 = [2.000000, 0.000000, 0.629112, 0.524948, 0.422396]
    channel_3D.unsteady  = True
    channel_3D.multizone = True
    channel_3D.enabled_with_tsan = False
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [0.092415, 0.568971, 0.692864, 0.989451, 1.048245]
    pipe.unsteady  = True
    pipe.multizone = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [3.000000, 0.000000, 0.664817, 1.125803, 1.117608]
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone  = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 0.899645, 1.076181]
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone  = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [13.000000, -0.391571, -1.460876]
    bars_SST_2D.multizone = True
    test_list.append(bars_SST_2D)

    # Sliding mesh with incompressible flows (steady)
    slinc_steady           = TestCase('slinc_steady')
    slinc_steady.cfg_dir   = "sliding_interface/incompressible_steady"
    slinc_steady.cfg_file  = "config.cfg"
    slinc_steady.test_iter = 19
    slinc_steady.test_vals         = [19.000000, -1.144102, -1.424986]
    slinc_steady.test_vals_aarch64 = [19.000000, -1.154874, -1.378120]
    slinc_steady.multizone = True
    test_list.append(slinc_steady)

    ##########################
    ### FEA - FSI          ###
    ##########################

    # Static beam, 3d
    statbeam3d           = TestCase('statbeam3d')
    statbeam3d.cfg_dir   = "fea_fsi/StatBeam_3d"
    statbeam3d.cfg_file  = "configBeam_3d.cfg"
    statbeam3d.test_iter = 0
    statbeam3d.test_vals = [-2.862081, -1.795143, -2.380136, 110350.000000]
    statbeam3d.test_vals_aarch64 = [-2.777602, -1.710976, -2.445072, 110350]
    test_list.append(statbeam3d)

    # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.test_iter = 6
    dynbeam2d.test_vals = [-3.240013, 2.895060, -0.353141, 76220]
    dynbeam2d.unsteady  = True
    test_list.append(dynbeam2d)

    # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [4.000000, 0.000000, -3.726029, -4.277532]
    fsi2d.multizone= True
    fsi2d.unsteady = True
    fsi2d.enabled_with_tsan = False
    test_list.append(fsi2d)

    # FSI, Dynamic, 2D, new mesh solver
    dyn_fsi           = TestCase('dyn_fsi')
    dyn_fsi.cfg_dir   = "fea_fsi/dyn_fsi"
    dyn_fsi.cfg_file  = "config.cfg"
    dyn_fsi.test_iter = 4
    dyn_fsi.test_vals = [-4.330727, -4.152808, 0.000000, 102.000000]
    dyn_fsi.multizone = True
    dyn_fsi.unsteady  = True
    test_list.append(dyn_fsi)

    # FSI+CHT, Static, 2D, new mesh solver, restart
    fsi_cht_restart           = TestCase('fsi_cht_restart')
    fsi_cht_restart.cfg_dir   = "fea_fsi/stat_fsi"
    fsi_cht_restart.cfg_file  = "config_restart.cfg"
    fsi_cht_restart.test_iter = 0
    fsi_cht_restart.test_vals = [5.000000, 0.006352, -1.960362, -9.327033, -9.643180, -9.319159, 608.380000, -0.012974, 0.000000, 20.000000]
    fsi_cht_restart.multizone = True
    test_list.append(fsi_cht_restart)

    ##############################################
    ### Method of Manufactured Solutions (MMS) ###
    ##############################################

    # FVM, compressible, laminar N-S
    mms_fvm_ns           = TestCase('mms_fvm_ns')
    mms_fvm_ns.cfg_dir   = "mms/fvm_navierstokes"
    mms_fvm_ns.cfg_file  = "lam_mms_roe.cfg"
    mms_fvm_ns.test_iter = 20
    mms_fvm_ns.test_vals = [-2.808514, 2.152654, 0.000000, 0.000000]
    test_list.append(mms_fvm_ns)

    # FVM, incompressible, euler
    mms_fvm_inc_euler           = TestCase('mms_fvm_inc_euler')
    mms_fvm_inc_euler.cfg_dir   = "mms/fvm_incomp_euler"
    mms_fvm_inc_euler.cfg_file  = "inv_mms_jst.cfg"
    mms_fvm_inc_euler.test_iter = 20
    mms_fvm_inc_euler.test_vals         = [-9.128034, -9.441406, 0.000000, 0.000000]
    mms_fvm_inc_euler.test_vals_aarch64 = [-9.128034, -9.441406, 0.000000, 0.000000]
    test_list.append(mms_fvm_inc_euler)

    # FVM, incompressible, laminar N-S
    mms_fvm_inc_ns           = TestCase('mms_fvm_inc_ns')
    mms_fvm_inc_ns.cfg_dir   = "mms/fvm_incomp_navierstokes"
    mms_fvm_inc_ns.cfg_file  = "lam_mms_fds.cfg"
    mms_fvm_inc_ns.test_iter = 20
    mms_fvm_inc_ns.test_vals = [-7.414944, -7.631546, 0.000000, 0.000000]
    test_list.append(mms_fvm_inc_ns)

    ##########################
    ###   Python wrapper   ###
    ##########################

    # NACA0012
    pywrapper_translating_naca0012 = TestCase('pywrapper_translating_naca0012')
    pywrapper_translating_naca0012.cfg_dir = "py_wrapper/translating_NACA0012"
    pywrapper_translating_naca0012.cfg_file = "config.cfg"
    pywrapper_translating_naca0012.command = TestCase.Command(exec = "python", param = "run_su2.py")
    pywrapper_translating_naca0012.timeout = 60
    pywrapper_translating_naca0012.reference_file = "forces_0.csv.ref"
    pywrapper_translating_naca0012.reference_file_aarch64 = "forces_0_aarch64.csv.ref"
    pywrapper_translating_naca0012.test_file = "forces_0.csv"
    pywrapper_translating_naca0012.tol_file_percent = 5.0
    pywrapper_translating_naca0012.enabled_on_cpu_arch = ["x86_64"]
    pywrapper_translating_naca0012.enabled_with_tsan = False
    file_diff_list.append(pywrapper_translating_naca0012)

    # NACA0012 with updated moving frame
    pywrapper_updated_moving_frame_naca0012 = TestCase('pywrapper_updated_moving_frame_naca0012')
    pywrapper_updated_moving_frame_naca0012.cfg_dir = "py_wrapper/updated_moving_frame_NACA12"
    pywrapper_updated_moving_frame_naca0012.cfg_file = "config.cfg"
    pywrapper_updated_moving_frame_naca0012.command = TestCase.Command(exec = "python", param = "run_su2.py")
    pywrapper_updated_moving_frame_naca0012.timeout = 60
    pywrapper_updated_moving_frame_naca0012.reference_file = "forces_0.csv.ref"
    pywrapper_updated_moving_frame_naca0012.reference_file_aarch64 = "forces_0_aarch64.csv.ref"
    pywrapper_updated_moving_frame_naca0012.test_file = "forces_0.csv"
    pywrapper_updated_moving_frame_naca0012.tol_file_percent = 5.0
    pywrapper_updated_moving_frame_naca0012.enabled_on_cpu_arch = ["x86_64"]
    pywrapper_updated_moving_frame_naca0012.enabled_with_tsan = False
    file_diff_list.append(pywrapper_updated_moving_frame_naca0012)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    for test in test_list:
        if test.command.empty():
            test.command = TestCase.Command(exec = "SU2_CFD", param = "-t 2")
        if test.timeout == 0:
            test.timeout = 600
        test.tol = 1e-4
    #end

    pass_list = [ test.run_test(args.tsan) for test in test_list ]
    pass_list += [ test.run_filediff(args.tsan) for test in file_diff_list ]

    # Tests summary
    print('==================================================================')
    print('Summary of the hybrid parallel tests')
    print('python version:', sys.version)
    for i, test in enumerate(test_list+file_diff_list):
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
