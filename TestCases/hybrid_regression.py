#!/usr/bin/env python

## \file hybrid_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 8.2.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
    channel.test_vals = [-2.965642, 2.459171, 0.016012, 0.042270]
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.766168, -4.287699, 0.326688, 0.022661]
    test_list.append(naca0012)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-1.379426, 4.288828, -0.245341, 0.043244]
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [0.280800, 0.008623]
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-3.896832, 1.637749, 0.301084, 0.019485]
    test_list.append(fixedCL_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.475463, 6.835018, 0.000226, 1.784354]
    test_list.append(bluntbody)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 100
    flatplate.test_vals = [-7.680034, -2.207912, 0.001084, 0.036233, 2.361500, -2.325300, 0.000000, 0.000000]
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-8.266513, -2.783904, -0.019899, 1.615668, 0]
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals         = [-6.830996, -1.368850, -0.143956, 73.963354, 0]
    cylinder_lowmach.test_vals_aarch64 = [-6.830996, -1.368850, -0.143956, 73.963354, 0]
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-5.046131, 0.652984, 0.008355, 13.735818, 0]
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals         = [-12.009016, -7.262446, -0.000000, 2.089953]
    poiseuille_profile.test_vals_aarch64 = [-12.494717, -7.711274, -0.000000, 2.085796]
    test_list.append(poiseuille_profile)

    # 2D Rotational Periodic
    periodic2d           = TestCase('periodic2d')
    periodic2d.cfg_dir   = "navierstokes/periodic2D"
    periodic2d.cfg_file  = "config.cfg"
    periodic2d.test_iter = 1400
    periodic2d.test_vals = [-10.817611, -8.363544, -8.287460, -5.334104, -1.088411, -2945.2]
    test_list.append(periodic2d)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.020123, -5.269339, 0.807147, 0.060499, 0]
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510349, 4.950289, 0.811983, 0.061600, 0.000000]
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.678140, 4.950289, 0.811983, 0.061600]
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.316134, -6.737979, -0.187461, 0.057468]
    test_list.append(turb_flatplate)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.408675, -6.662904, 0.238578, 0.158968, 0.000000]
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 5
    turb_naca0012_sa.test_vals         = [-12.098325, -14.149988, 1.057665, 0.022971, 20.000000, -2.292707, 0.000000, -12.068169, 0]
    turb_naca0012_sa.test_vals_aarch64 = [-12.098325, -14.149988, 1.057665, 0.022971, 20.000000, -2.292707, 0.000000, -12.068169, 0]
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-12.232646, -14.434809, -6.037672, 1.047444, 0.019214, -1.652532, 0.000000]
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals = [-12.152294, -14.781143, -6.357944, 1.000270, 0.019123, -1.417707]
    test_list.append(turb_naca0012_sst_sust)

    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [-5.192404, -10.250818, -1.560208, 1.022562, 0.040530, -2.382731]
    test_list.append(turb_naca0012_sst_fixedvalues)

    # NACA0012 (SST, explicit Euler for flow and turbulence equations)
    turb_naca0012_sst_expliciteuler           = TestCase('turb_naca0012_sst_expliciteuler')
    turb_naca0012_sst_expliciteuler.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_expliciteuler.cfg_file  = "turb_NACA0012_sst_expliciteuler.cfg"
    turb_naca0012_sst_expliciteuler.test_iter = 10
    turb_naca0012_sst_expliciteuler.test_vals = [-3.533765, -3.157766, 3.364026, 1.124757, 0.501700, -float("inf")]
    test_list.append(turb_naca0012_sst_expliciteuler)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389724, -8.409502, 0.000048, 0.056344]
    test_list.append(propeller)

    #######################################
    ### Axisymmetric Compressible RANS  ###
    #######################################

    # Axisymmetric air nozzle (transonic)
    axi_rans_air_nozzle_restart           = TestCase('axi_rans_air_nozzle_restart')
    axi_rans_air_nozzle_restart.cfg_dir   = "axisymmetric_rans/air_nozzle"
    axi_rans_air_nozzle_restart.cfg_file  = "air_nozzle_restart.cfg"
    axi_rans_air_nozzle_restart.test_iter = 10
    axi_rans_air_nozzle_restart.test_vals = [-12.065954, -6.836372, -8.889803, -3.832665, 0.000000]
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
    turb_naca0012_sst_restart_mg.test_vals = [-7.633091, -7.181942, -0.627082, -0.000020, 0.078737]
    test_list.append(turb_naca0012_sst_restart_mg)

    #############################
    ### Compressibele RANS UQ ###
    #############################

    # NACA0012 1c
    turb_naca0012_1c           = TestCase('turb_naca0012_1c')
    turb_naca0012_1c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_1c.cfg_file  = "turb_NACA0012_uq_1c.cfg"
    turb_naca0012_1c.test_iter = 10
    turb_naca0012_1c.test_vals         = [-4.976639, 1.141468, 0.243064, -0.112166]
    turb_naca0012_1c.test_vals_aarch64 = [-4.981105, 1.138873, 0.248013, -0.117248]
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals         = [-5.486005, 0.968502, 0.233564, -0.114514]
    turb_naca0012_2c.test_vals_aarch64 = [-5.483345, 0.968720, 0.214914, -0.124932]
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals         = [-5.584366, 0.931883, 0.223356, -0.116142]
    turb_naca0012_3c.test_vals_aarch64 = [-5.584300, 0.931293, 0.207447, -0.125691]
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals         = [-5.114311, 1.076961, 0.224322, -0.118874]
    turb_naca0012_p1c1.test_vals_aarch64 = [-5.132358, 1.075658, 0.337268, -0.082827]
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals         = [-5.548834, 0.946383, 0.211109, -0.121304]
    turb_naca0012_p1c2.test_vals_aarch64 = [-5.548775, 0.945962, 0.211150, -0.121291]
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
    hb_rans_preconditioning.test_vals = [-1.902111, 0.484080, 0.601469, 3.608991, -5.949369]
    test_list.append(hb_rans_preconditioning)

    #############################
    ### Incompressible Euler  ###
    #############################

    # NACA0012 Hydrofoil
    inc_euler_naca0012           = TestCase('inc_euler_naca0012')
    inc_euler_naca0012.cfg_dir   = "incomp_euler/naca0012"
    inc_euler_naca0012.cfg_file  = "incomp_NACA0012.cfg"
    inc_euler_naca0012.test_iter = 20
    inc_euler_naca0012.test_vals = [-7.127256, -6.466554, 0.531991, 0.008466]
    test_list.append(inc_euler_naca0012)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals         = [-6.595923, -5.820518, -0.018365, 0.126411]
    inc_nozzle.test_vals_aarch64 = [-5.624385, -4.988472, -0.000096, 0.137032]
    test_list.append(inc_nozzle)

    #############################
    ### Incompressible N-S    ###
    #############################

    # Laminar cylinder
    inc_lam_cylinder          = TestCase('inc_lam_cylinder')
    inc_lam_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_lam_cylinder.cfg_file  = "incomp_cylinder.cfg"
    inc_lam_cylinder.test_iter = 10
    inc_lam_cylinder.test_vals = [-4.004277, -3.227956, 0.003851, 7.626583]
    test_list.append(inc_lam_cylinder)

    # Buoyancy-driven cavity
    inc_buoyancy          = TestCase('inc_buoyancy')
    inc_buoyancy.cfg_dir   = "incomp_navierstokes/buoyancy_cavity"
    inc_buoyancy.cfg_file  = "lam_buoyancy_cavity.cfg"
    inc_buoyancy.test_iter = 20
    inc_buoyancy.test_vals = [-4.432484, 0.507522, 0.000000, 0.000000]
    test_list.append(inc_buoyancy)

    # Laminar heated cylinder with polynomial fluid model
    inc_poly_cylinder          = TestCase('inc_poly_cylinder')
    inc_poly_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_poly_cylinder.cfg_file  = "poly_cylinder.cfg"
    inc_poly_cylinder.test_iter = 20
    inc_poly_cylinder.test_vals         = [-7.827942, -2.061513, 0.029525, 1.953498, -174.780000]
    inc_poly_cylinder.test_vals_aarch64 = [-7.827942, -2.061513, 0.029525, 1.953498, -174.780000]
    test_list.append(inc_poly_cylinder)

    # X-coarse laminar bend as a mixed element CGNS test
    inc_lam_bend          = TestCase('inc_lam_bend')
    inc_lam_bend.cfg_dir   = "incomp_navierstokes/bend"
    inc_lam_bend.cfg_file  = "lam_bend.cfg"
    inc_lam_bend.test_iter = 10
    inc_lam_bend.test_vals         = [-3.560185, -3.051988, -0.013972, 1.102842]
    inc_lam_bend.test_vals_aarch64 = [-3.437996, -3.086188, -0.015600, 1.142213]
    test_list.append(inc_lam_bend)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012, SA
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.788405, -11.040560, 0.000008, 0.309505]
    test_list.append(inc_turb_naca0012)

    # NACA0012, SST_SUST
    inc_turb_naca0012_sst_sust           = TestCase('inc_turb_naca0012_sst_sust')
    inc_turb_naca0012_sst_sust.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_sst_sust.cfg_file  = "naca0012_SST_SUST.cfg"
    inc_turb_naca0012_sst_sust.test_iter = 20
    inc_turb_naca0012_sst_sust.test_vals = [-7.291693, 0.132607, 0.000002, 0.312092]
    test_list.append(inc_turb_naca0012_sst_sust)

    # Weakly coupled heat equation
    inc_weakly_coupled = TestCase('inc_weakly_coupled')
    inc_weakly_coupled.cfg_dir = "disc_adj_heat"
    inc_weakly_coupled.cfg_file = "primal.cfg"
    inc_weakly_coupled.test_iter = 10
    inc_weakly_coupled.test_vals = [-18.894811, -17.879327, -18.412938, -17.855948, -18.343462, -15.659612, 5.545700]
    test_list.append(inc_weakly_coupled)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.627869, -0.164403, 0.054734, 2.545856]
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals         = [-8.008048, -2.611074, 1.497289, 1.487468]
    spinning_cylinder.test_vals_aarch64 = [-8.006541, -2.609759, 1.495662, 1.486341]
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-2.560838, -1.175929, 0.062081, 1.399401, 2.220361, 1.399349, 2.218600, 0.000000]
    square_cylinder.test_vals_aarch64 = [-2.557902, -1.173574, 0.058050, 1.399794, 2.220402, 1.399748, 2.218604, 0]
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977498, 3.481818, -0.010484, -0.008178]
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Cosine gust in z-direction
    cosine_gust           = TestCase('cosine_gust_zdir')
    cosine_gust.cfg_dir   = "gust"
    cosine_gust.cfg_file  = "cosine_gust_zdir.cfg"
    cosine_gust.test_iter = 79
    cosine_gust.test_vals = [-2.418805, 0.001949, -0.001254, 0.000425, -0.000593]
    cosine_gust.unsteady  = True
    cosine_gust.enabled_with_tsan = False
    test_list.append(cosine_gust)

    # Gust with mesh deformation
    gust_mesh_defo           = TestCase('gust_with_mesh_deformation')
    gust_mesh_defo.cfg_dir   = "gust"
    gust_mesh_defo.cfg_file  = "gust_with_mesh_deformation.cfg"
    gust_mesh_defo.test_iter = 6
    gust_mesh_defo.test_vals = [-1.844761, 0.001095, -0.000273]
    gust_mesh_defo.unsteady  = True
    gust_mesh_defo.enabled_with_tsan = False
    test_list.append(gust_mesh_defo)

    # Aeroelastic
    aeroelastic           = TestCase('aeroelastic')
    aeroelastic.cfg_dir   = "aeroelastic"
    aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    aeroelastic.test_iter = 2
    aeroelastic.test_vals         = [0.074052, 0.027623, -0.001641, -0.000128]
    aeroelastic.test_vals_aarch64 = [0.074170, 0.027590, -0.001579, -0.000160]
    aeroelastic.unsteady  = True
    aeroelastic.enabled_on_cpu_arch = ["x86_64"] # Requires AVX-capable architecture
    aeroelastic.enabled_with_tsan = False
    test_list.append(aeroelastic)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714713, -5.788302, -0.214960, 0.023758, 0.000000]
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.008629, -6.888996, 1.435193, 0.433537]
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # unsteady pitching NACA0012, Euler, Deforming
    unst_deforming_naca0012           = TestCase('unst_deforming_naca0012')
    unst_deforming_naca0012.cfg_dir   = "disc_adj_euler/naca0012_pitching_def"
    unst_deforming_naca0012.cfg_file  = "inv_NACA0012_pitching_deform.cfg"
    unst_deforming_naca0012.test_iter = 5
    unst_deforming_naca0012.test_vals = [-3.665168, -3.793307, -3.716526, -3.148348]
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
    edge_VW.test_iter = 40
    edge_VW.test_vals = [-5.681149, 0.463233, -0.000009, 0.000000]
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 40
    edge_PPR.test_vals         = [-7.139177, -0.980792, -0.000034, 0.000000]
    edge_PPR.test_vals_aarch64 = [-8.573595, -2.391849, -0.000034, 0.000000]
    test_list.append(edge_PPR)

    ######################################
    ### Turbomachinery                 ###
    ######################################

	# Jones APU Turbocharger restart
    Jones_tc_restart           = TestCase('jones_turbocharger_restart')
    Jones_tc_restart.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_restart.cfg_file  = "Jones_restart.cfg"
    Jones_tc_restart.test_iter = 5
    Jones_tc_restart.test_vals = [-7.308010, -5.332065, -14.895822, -9.330703, -12.071733, -6.548623, 73291.000000, 73291.000000, 0.020111, 82.896000]
    test_list.append(Jones_tc_restart)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals         = [1.090053, 1.550934, -2.895064, 2.607596, -2.479704, 3.063740, 106380.000000, 106380.000000, 5.733600, 64.747000]
    axial_stage2D.test_vals_aarch64 = [0.983739, 1.534333, -2.888521, 2.606770, -2.418339, 3.087275, 106380, 106380, 5.7325,  64.711]
    test_list.append(axial_stage2D)

    # 2D transonic stator restart
    transonic_stator_restart           = TestCase('transonic_stator_restart')
    transonic_stator_restart.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_restart.cfg_file  = "transonic_stator_restart.cfg"
    transonic_stator_restart.test_iter = 20
    transonic_stator_restart.test_vals         = [-4.357591, -2.480153, -2.075008, 1.737627, -1.441579, 3.240658, -471620.000000, 94.840000, -0.054589]
    transonic_stator_restart.test_vals_aarch64 = [-4.357748, -2.480402, -2.075152, 1.737469, -1.440919, 2.727299, -471620.000000, 94.840000, -0.054603]
    test_list.append(transonic_stator_restart)

    # Multiple turbomachinery interface restart
    multi_interface                    = TestCase('multi_interface')
    multi_interface.cfg_dir            = "turbomachinery/multi_interface"
    multi_interface.cfg_file           = "multi_interface_rst.cfg"
    multi_interface.test_iter          = 5
    multi_interface.test_vals          = [-8.632374, -8.895124, -9.350417]
    multi_interface.test_vals_aarch64  = [-8.632374, -8.895124, -9.350417]
    test_list.append(multi_interface)

    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 5
    uniform_flow.test_vals = [5.000000, 0.000000, -0.195002, -10.624444]
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow)

    # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 2
    channel_2D.test_vals = [2.000000, 0.000000, 0.464907, 0.348052, 0.397452]
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 2
    channel_3D.test_vals         = [2.000000, 0.000000, 0.629106, 0.524901, 0.422380]
    channel_3D.test_vals_aarch64 = [2.000000, 0.000000, 0.620558, 0.504323, 0.412729]
    channel_3D.unsteady  = True
    channel_3D.multizone = True
    channel_3D.enabled_with_tsan = False
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [0.080826, 0.547321, 0.655098, 0.968229, 1.049129]
    pipe.unsteady  = True
    pipe.multizone = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [3.000000, 0.000000, 0.717065, 1.119817, 1.160326]
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone  = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 1.207118, 1.065254]
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone  = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [13.000000, 1.167903, -1.660901]
    bars_SST_2D.multizone = True
    test_list.append(bars_SST_2D)

    # Sliding mesh with incompressible flows (steady)
    slinc_steady           = TestCase('slinc_steady')
    slinc_steady.cfg_dir   = "sliding_interface/incompressible_steady"
    slinc_steady.cfg_file  = "config.cfg"
    slinc_steady.test_iter = 19
    slinc_steady.test_vals         = [19.000000, -1.799803, -2.128200]
    slinc_steady.test_vals_aarch64 = [19.000000, -1.799803, -2.128199]
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
    statbeam3d.test_vals         = [-2.787802, -1.721974, -2.438436, 110350]
    statbeam3d.test_vals_aarch64 = [-2.787802, -1.721974, -2.438436, 110350]
    test_list.append(statbeam3d)

    # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.test_iter = 6
    dynbeam2d.test_vals = [-3.240016, 2.895057, -0.353147, 66127.000000]
    dynbeam2d.unsteady  = True
    test_list.append(dynbeam2d)

    # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [4.000000, 0.000000, -3.726029, -4.277769]
    fsi2d.multizone= True
    fsi2d.unsteady = True
    fsi2d.enabled_with_tsan = False
    test_list.append(fsi2d)

    # FSI, Dynamic, 2D, new mesh solver
    dyn_fsi           = TestCase('dyn_fsi')
    dyn_fsi.cfg_dir   = "fea_fsi/dyn_fsi"
    dyn_fsi.cfg_file  = "config.cfg"
    dyn_fsi.test_iter = 4
    dyn_fsi.test_vals         = [-4.330725, -4.152983, 0.000000, 103.000000]
    dyn_fsi.test_vals_aarch64 = [-4.332167, -4.057742, 0.000000, 102.000000]
    dyn_fsi.multizone = True
    dyn_fsi.unsteady  = True
    test_list.append(dyn_fsi)

    # FSI, Static, 2D, new mesh solver, restart
    stat_fsi_restart           = TestCase('stat_fsi_restart')
    stat_fsi_restart.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi_restart.cfg_file  = "config_restart.cfg"
    stat_fsi_restart.test_iter = 1
    stat_fsi_restart.test_vals         = [-3.549586, -4.460650, 0.000000, 35.000000]
    stat_fsi_restart.test_vals_aarch64 = [-3.442878, -4.228058, 0.000000, 37.000000]
    stat_fsi_restart.multizone = True
    test_list.append(stat_fsi_restart)

    ##############################################
    ### Method of Manufactured Solutions (MMS) ###
    ##############################################

    # FVM, compressible, laminar N-S
    mms_fvm_ns           = TestCase('mms_fvm_ns')
    mms_fvm_ns.cfg_dir   = "mms/fvm_navierstokes"
    mms_fvm_ns.cfg_file  = "lam_mms_roe.cfg"
    mms_fvm_ns.test_iter = 20
    mms_fvm_ns.test_vals = [-2.808514, 2.152655, 0.000000, 0.000000]
    test_list.append(mms_fvm_ns)

    # FVM, incompressible, euler
    mms_fvm_inc_euler           = TestCase('mms_fvm_inc_euler')
    mms_fvm_inc_euler.cfg_dir   = "mms/fvm_incomp_euler"
    mms_fvm_inc_euler.cfg_file  = "inv_mms_jst.cfg"
    mms_fvm_inc_euler.test_iter = 20
    mms_fvm_inc_euler.test_vals         = [-9.128033, -9.441406, 0.000000, 0.000000]
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
    pywrapper_translating_naca0012.tol_file_percent = 0.1
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
    pywrapper_updated_moving_frame_naca0012.tol_file_percent = 0.1
    pywrapper_updated_moving_frame_naca0012.enabled_on_cpu_arch = ["x86_64"]
    pywrapper_updated_moving_frame_naca0012.enabled_with_tsan = False
    file_diff_list.append(pywrapper_updated_moving_frame_naca0012)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    for test in test_list:
        test.command = TestCase.Command(exec = "SU2_CFD", param = "-t 2")
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
