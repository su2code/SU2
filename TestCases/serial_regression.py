#!/usr/bin/env python

## \file serial_regression.py
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

from __future__ import print_function, division, absolute_import
import sys
from TestCase import TestCase
from TestCase import parse_args

def main():
    '''This program runs SU2 and ensures that the output matches specified values.
       This will be used to do checks when code is pushed to github
       to make sure nothing is broken. '''

    args = parse_args('Serial Regression Tests')

    test_list = []

    #########################
    ## NEMO solver ###
    #########################

    # Adiabatic thermal bath
    thermalbath           = TestCase('thermalbath')
    thermalbath.cfg_dir   = "nonequilibrium/thermalbath/finitechemistry"
    thermalbath.cfg_file  = "thermalbath.cfg"
    thermalbath.test_iter = 10
    thermalbath.test_vals = [0.945997, 0.945997, -12.018025, -12.217291, -32.000000, 10.013239]
    test_list.append(thermalbath)

    # Adiabatic frozen thermal bath
    thermalbath_frozen           = TestCase('thermalbath_frozen')
    thermalbath_frozen.cfg_dir   = "nonequilibrium/thermalbath/frozen"
    thermalbath_frozen.cfg_file  = "thermalbath_frozen.cfg"
    thermalbath_frozen.test_iter = 10
    thermalbath_frozen.test_vals = [-32.000000, -32.000000, -12.018022, -11.978730, -32.000000, 10.013545]
    test_list.append(thermalbath_frozen)

    # Inviscid single wedge, implicit
    invwedge = TestCase('invwedge')
    invwedge.cfg_dir = "nonequilibrium/invwedge"
    invwedge.cfg_file = "invwedge_ausm.cfg"
    invwedge.test_iter = 10
    invwedge.test_vals = [-1.073699, -1.598462, -18.299911, -18.627322, -18.573334, 2.241760, 1.868575, 5.286072, 0.843741]
    invwedge.test_vals_aarch64 = [-1.046323, -1.571086, -18.301361, -18.628744, -18.574788, 2.271778, 1.875687, 5.315769, 0.870008]
    test_list.append(invwedge)

    # Viscous single cone - axisymmetric
    visc_cone = TestCase('visc_cone')
    visc_cone.cfg_dir = "nonequilibrium/visc_wedge"
    visc_cone.cfg_file = "axi_visccone.cfg"
    visc_cone.test_iter = 10
    visc_cone.test_vals = [-5.215239, -5.739373, -20.560910, -20.517094, -20.406632, 1.262779, -3.205484, -0.015695, 0.093205, 32641.000000]
    visc_cone.test_vals_aarch64 = [-5.215229, -5.739368, -20.556662, -20.517022, -20.437459, 1.262784, -3.205455, -0.015696, 0.093207, 32656.000000]
    test_list.append(visc_cone)

    #########################
    ## Compressible Euler ###
    #########################

    # Dry run test Euler
    channel_d           = TestCase('dry run Euler')
    channel_d.cfg_dir   = "euler/channel"
    channel_d.cfg_file  = "inv_channel_RK.cfg"
    channel_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(channel_d)

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 10
    channel.test_vals = [-2.691660, 2.781413, -0.009401, 0.011862]
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.766184, -4.287722, 0.326688, 0.022661]
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
    oneram6.test_vals = [-11.498143, -10.969216, 0.280800, 0.008623]
    oneram6.timeout   = 9600
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-3.837516, 1.700577, 0.301169, 0.019490]
    test_list.append(fixedCL_naca0012)

    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals         = [-1.077848, 4.386916, -0.000333, 0.029678]
    polar_naca0012.test_vals_aarch64 = [-1.063447, 4.401847, 0.000291, 0.031696]
    polar_naca0012.command   = TestCase.Command(exec = "compute_polar.py", param = "-n 1 -i 11")
    # flaky test on arm64
    polar_naca0012.enabled_on_cpu_arch = ["x86_64"]
    test_list.append(polar_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.475378, 6.834898, 0.000000, 1.783956]
    test_list.append(bluntbody)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Dry run test NS
    flatplate_d           = TestCase('dry run NS')
    flatplate_d.cfg_dir   = "navierstokes/flatplate"
    flatplate_d.cfg_file  = "lam_flatplate.cfg"
    flatplate_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(flatplate_d)

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 20
    flatplate.test_vals = [-5.097199, 0.382306, 0.001326, 0.027904, 2.361500, -2.333600, 0.000000, 0.000000]
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-8.363995, -2.882536, -0.017780, 1.607979, 0]
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.830989, -1.368842, -0.143838, 73.962440, 0]
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-5.050753, 0.648333, 0.012273, 13.643141, 0]
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals         = [-12.009021, -7.262796, -0.000000, 2.089953]
    poiseuille_profile.test_vals_aarch64 = [-12.494684, -7.711379, -0.000000, 2.085796] #last 4 columns
    test_list.append(poiseuille_profile)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # Dry run RANS
    rae2822_sa_d           = TestCase('dry run RANS')
    rae2822_sa_d .cfg_dir   = "rans/rae2822"
    rae2822_sa_d .cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa_d .command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(rae2822_sa_d)

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
    rae2822_sst.test_vals = [-0.510326, 4.944645, 0.812737, 0.061081, 0.000000]
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.680659, 4.944645, 0.812737, 0.061080]
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.316135, -6.737979, -0.187461, 0.057468]
    test_list.append(turb_flatplate)

    # FLAT PLATE, WALL FUNCTIONS, COMPRESSIBLE SST
    turb_wallfunction_flatplate_sst           = TestCase('turb_sst_wallfunction_flatplate')
    turb_wallfunction_flatplate_sst.cfg_dir   = "wallfunctions/flatplate/compressible_SST"
    turb_wallfunction_flatplate_sst.cfg_file  = "turb_SST_flatplate.cfg"
    turb_wallfunction_flatplate_sst.test_iter = 10
    turb_wallfunction_flatplate_sst.test_vals = [-4.587729, -1.908761, -1.968077, 0.903297, -1.608475, 1.497162, 10.000000, -1.733571, 0.033374, 0.002430]
    test_list.append(turb_wallfunction_flatplate_sst)

    # FLAT PLATE, WALL FUNCTIONS, COMPRESSIBLE SA
    turb_wallfunction_flatplate_sa           = TestCase('turb_sa_wallfunction_flatplate')
    turb_wallfunction_flatplate_sa.cfg_dir   = "wallfunctions/flatplate/compressible_SA"
    turb_wallfunction_flatplate_sa.cfg_file  = "turb_SA_flatplate.cfg"
    turb_wallfunction_flatplate_sa.test_iter = 10
    turb_wallfunction_flatplate_sa.test_vals = [-4.487562, -2.016144, -2.169439, 0.834724, -5.382532, 10.000000, -1.508186, 0.034484, 0.002639]
    test_list.append(turb_wallfunction_flatplate_sa)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.408685, -6.662908, 0.238580, 0.158968, 0.000000]
    turb_oneram6.timeout   = 3200
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D results for finest grid: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 5
    turb_naca0012_sa.test_vals = [-12.091778, -14.685322, 1.057665, 0.022971, 20.000000, -2.686203, 20.000000, -4.459916, 0]
    turb_naca0012_sa.timeout   = 3200
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals         = [-12.229889, -14.434883, -6.037177, 1.047444, 0.019214, -2.108366, 0.000000]
    turb_naca0012_sst.test_vals_aarch64 = [-12.229890, -14.434837, -6.410709, 1.047444, 0.019214, -2.107944, 0.000000]
    turb_naca0012_sst.timeout   = 3200
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST V2003m, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_2003m           = TestCase('turb_naca0012_sst_2003m')
    turb_naca0012_sst_2003m.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_2003m.cfg_file  = "turb_NACA0012_sst_2003m.cfg"
    turb_naca0012_sst_2003m.test_iter = 10
    turb_naca0012_sst_2003m.test_vals = [-8.617103, -10.718898, -4.021920, 1.046418, 0.019227, -1.708431, 0.000000]
    turb_naca0012_sst_2003m.timeout   = 3200
    test_list.append(turb_naca0012_sst_2003m)

    # NACA0012 (SST_SUST, FUN3D results for finest grid: CL=1.0840, CD=0.01253) restart
    turb_naca0012_sst_sust_restart           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust_restart.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust_restart.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust_restart.test_iter = 10
    turb_naca0012_sst_sust_restart.test_vals = [-12.157353, -14.781748, -6.359242, 1.000270, 0.019123, -1.780616]
    turb_naca0012_sst_sust_restart.test_vals_aarch64 = [-12.157374, -14.782027, -6.726462, 1.000270, 0.019123, -1.780624]
    turb_naca0012_sst_sust_restart.timeout   = 3200
    test_list.append(turb_naca0012_sst_sust_restart)

    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [-5.206634, -10.233383, -1.557783, 1.022024, 0.040554, -3.477756]
    turb_naca0012_sst_fixedvalues.timeout   = 3200
    test_list.append(turb_naca0012_sst_fixedvalues)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389724, -8.409502, 0.000048, 0.056344]
    propeller.timeout   = 3200
    test_list.append(propeller)

    #######################################
    ### Axisymmetric Compressible RANS  ###
    #######################################

    # Axisymmetric air nozzle (transonic) restart
    axi_rans_air_nozzle_restart           = TestCase('axi_rans_air_nozzle_restart')
    axi_rans_air_nozzle_restart.cfg_dir   = "axisymmetric_rans/air_nozzle"
    axi_rans_air_nozzle_restart.cfg_file  = "air_nozzle_restart.cfg"
    axi_rans_air_nozzle_restart.test_iter = 10
    axi_rans_air_nozzle_restart.test_vals = [-12.067274, -6.841605, -8.843340, -3.783551, 0.000000]
    axi_rans_air_nozzle_restart.test_vals_aarch64 = [-12.067037, -6.840810, -8.843191, -3.783401, 0.000000]
    axi_rans_air_nozzle_restart.tol       = 0.0001
    test_list.append(axi_rans_air_nozzle_restart)

    #################################
    ## Compressible RANS Restart  ###
    #################################

    # NACA0012 SST Multigrid restart
    turb_naca0012_sst_restart_mg           = TestCase('turb_naca0012_sst_restart_mg')
    turb_naca0012_sst_restart_mg.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_restart_mg.cfg_file  = "turb_NACA0012_sst_multigrid_restart.cfg"
    turb_naca0012_sst_restart_mg.test_iter = 50
    turb_naca0012_sst_restart_mg.ntest_vals = 5
    turb_naca0012_sst_restart_mg.test_vals = [-7.645513, -7.268192, -1.665656, -0.000027, 0.078643]
    turb_naca0012_sst_restart_mg.timeout   = 3200
    turb_naca0012_sst_restart_mg.tol       = 0.000001
    test_list.append(turb_naca0012_sst_restart_mg)

    #############################
    ### Incompressible Euler  ###
    #############################

    # Dry run Inc Euler
    inc_euler_naca0012_d           = TestCase('dry run Inc Euler')
    inc_euler_naca0012_d.cfg_dir   = "incomp_euler/naca0012"
    inc_euler_naca0012_d.cfg_file  = "incomp_NACA0012.cfg"
    inc_euler_naca0012_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(inc_euler_naca0012_d)

    # NACA0012 Hydrofoil
    inc_euler_naca0012           = TestCase('inc_euler_naca0012')
    inc_euler_naca0012.cfg_dir   = "incomp_euler/naca0012"
    inc_euler_naca0012.cfg_file  = "incomp_NACA0012.cfg"
    inc_euler_naca0012.test_iter = 20
    inc_euler_naca0012.test_vals = [-7.140809, -6.485990, 0.531993, 0.008466]
    test_list.append(inc_euler_naca0012)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals = [-6.623301, -5.844127, -0.023181, 0.126370]
    test_list.append(inc_nozzle)

    #############################
    ### Incompressible N-S    ###
    #############################

    # Dry Run Inc. NS
    inc_lam_cylinder_d          = TestCase('dry run Inc. NS')
    inc_lam_cylinder_d.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_lam_cylinder_d.cfg_file  = "incomp_cylinder.cfg"
    inc_lam_cylinder_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(inc_lam_cylinder_d)

    # Laminar cylinder
    inc_lam_cylinder          = TestCase('inc_lam_cylinder')
    inc_lam_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_lam_cylinder.cfg_file  = "incomp_cylinder.cfg"
    inc_lam_cylinder.test_iter = 10
    inc_lam_cylinder.test_vals = [-4.004277, -3.227956, 0.003852, 7.626578]
    test_list.append(inc_lam_cylinder)

    # Buoyancy-driven cavity
    inc_buoyancy          = TestCase('inc_buoyancy')
    inc_buoyancy.cfg_dir   = "incomp_navierstokes/buoyancy_cavity"
    inc_buoyancy.cfg_file  = "lam_buoyancy_cavity.cfg"
    inc_buoyancy.test_iter = 20
    inc_buoyancy.test_vals = [-4.436657, 0.507847, 0.000000, 0.000000]
    test_list.append(inc_buoyancy)

    # Laminar heated cylinder with polynomial fluid model
    inc_poly_cylinder          = TestCase('inc_poly_cylinder')
    inc_poly_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_poly_cylinder.cfg_file  = "poly_cylinder.cfg"
    inc_poly_cylinder.test_iter = 20
    inc_poly_cylinder.test_vals = [-8.083556, -2.134369, 0.018999, 1.932938, -173.730000]
    test_list.append(inc_poly_cylinder)

    # X-coarse laminar bend as a mixed element CGNS test
    inc_lam_bend          = TestCase('inc_lam_bend')
    inc_lam_bend.cfg_dir   = "incomp_navierstokes/bend"
    inc_lam_bend.cfg_file  = "lam_bend.cfg"
    inc_lam_bend.test_iter = 10
    inc_lam_bend.test_vals = [-3.647474, -3.230291, -0.016108, 1.085750]
    test_list.append(inc_lam_bend)

    ############################
    ### Incompressible RANS  ###
    ############################

    # Dry run Inc. RANS
    inc_turb_naca0012_d           = TestCase('dry run Inc. RANS')
    inc_turb_naca0012_d.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_d.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(inc_turb_naca0012_d)

    # NACA0012, SA
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.788495, -11.040578, 0.000023, 0.309503]
    test_list.append(inc_turb_naca0012)

    # NACA0012, SST_SUST
    inc_turb_naca0012_sst_sust           = TestCase('inc_turb_naca0012_sst_sust')
    inc_turb_naca0012_sst_sust.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_sst_sust.cfg_file  = "naca0012_SST_SUST.cfg"
    inc_turb_naca0012_sst_sust.test_iter = 20
    inc_turb_naca0012_sst_sust.test_vals = [-7.291334, 0.132662, 0.000021, 0.312093]
    test_list.append(inc_turb_naca0012_sst_sust)

    # FLAT PLATE, WALL FUNCTIONS, INCOMPRESSIBLE SST
    inc_turb_wallfunction_flatplate_sst           = TestCase('inc_turb_sst_wallfunction_flatplate')
    inc_turb_wallfunction_flatplate_sst.cfg_dir   = "wallfunctions/flatplate/incompressible_SST"
    inc_turb_wallfunction_flatplate_sst.cfg_file  = "turb_SST_flatplate.cfg"
    inc_turb_wallfunction_flatplate_sst.test_iter = 10
    inc_turb_wallfunction_flatplate_sst.test_vals = [-6.881306, -5.729111, -6.724803, -4.242636, -7.162846, -2.044959, 10.000000, -2.877924, 0.001151, 0.003161, 0.000000]
    test_list.append(inc_turb_wallfunction_flatplate_sst)

    # FLAT PLATE, WALL FUNCTIONS, INCOMPRESSIBLE SA
    inc_turb_wallfunction_flatplate_sa           = TestCase('inc_turb_sa_wallfunction_flatplate')
    inc_turb_wallfunction_flatplate_sa.cfg_dir   = "wallfunctions/flatplate/incompressible_SA"
    inc_turb_wallfunction_flatplate_sa.cfg_file  = "turb_SA_flatplate.cfg"
    inc_turb_wallfunction_flatplate_sa.test_iter = 10
    inc_turb_wallfunction_flatplate_sa.test_vals = [-6.894429, -5.717193, -6.743475, -4.242550, -9.587276, 10.000000, -2.879802, 0.001021, 0.003759]
    test_list.append(inc_turb_wallfunction_flatplate_sa)

    ####################
    ### DG-FEM Euler ###
    ####################

    # Dry run DG Euler
    fem_euler_naca0012_d           = TestCase('dry run DG Euler')
    fem_euler_naca0012_d.cfg_dir   = "hom_euler/NACA0012_5thOrder"
    fem_euler_naca0012_d.cfg_file  = "fem_NACA0012_reg.cfg"
    fem_euler_naca0012_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(fem_euler_naca0012_d)

    # NACA0012
    fem_euler_naca0012           = TestCase('fem_euler_naca0012')
    fem_euler_naca0012.cfg_dir   = "hom_euler/NACA0012_5thOrder"
    fem_euler_naca0012.cfg_file  = "fem_NACA0012_reg.cfg"
    fem_euler_naca0012.test_iter = 10
    fem_euler_naca0012.test_vals = [-6.519946, -5.976944, 0.255551, 0.000028]
    test_list.append(fem_euler_naca0012)

    ############################
    ### DG-FEM Navier-Stokes ###
    ############################

    # Dry run DG NS
    fem_ns_flatplate_d           = TestCase('dry run DG NS')
    fem_ns_flatplate_d.cfg_dir   = "hom_navierstokes/FlatPlate/nPoly4"
    fem_ns_flatplate_d.cfg_file  = "lam_flatplate_reg.cfg"
    fem_ns_flatplate_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(fem_ns_flatplate_d)

    # Flat plate
    fem_ns_flatplate           = TestCase('fem_ns_flatplate')
    fem_ns_flatplate.cfg_dir   = "hom_navierstokes/FlatPlate/nPoly4"
    fem_ns_flatplate.cfg_file  = "lam_flatplate_reg.cfg"
    fem_ns_flatplate.test_iter = 25
    fem_ns_flatplate.test_vals = [1.383727, 3.175247, 0.058387, 0.257951]
    test_list.append(fem_ns_flatplate)

    # Steady cylinder
    fem_ns_cylinder           = TestCase('fem_ns_cylinder')
    fem_ns_cylinder.cfg_dir   = "hom_navierstokes/CylinderViscous/nPoly3"
    fem_ns_cylinder.cfg_file  = "fem_Cylinder_reg.cfg"
    fem_ns_cylinder.test_iter = 10
    fem_ns_cylinder.test_vals = [0.454960, 0.979123, -0.000028, 79.984799]
    test_list.append(fem_ns_cylinder)

    # Steady sphere
    fem_ns_sphere           = TestCase('fem_ns_sphere')
    fem_ns_sphere.cfg_dir   = "hom_navierstokes/SphereViscous/nPoly3_QuadDominant"
    fem_ns_sphere.cfg_file  = "fem_Sphere_reg.cfg"
    fem_ns_sphere.test_iter = 10
    fem_ns_sphere.test_vals = [-0.288121, 0.240324, 0.000258, 21.797363]
    test_list.append(fem_ns_sphere)

    # Unsteady sphere ADER
    fem_ns_sphere_ader           = TestCase('fem_ns_sphere_ader')
    fem_ns_sphere_ader.cfg_dir   = "hom_navierstokes/SphereViscous/nPoly3_QuadDominant"
    fem_ns_sphere_ader.cfg_file  = "fem_Sphere_reg_ADER.cfg"
    fem_ns_sphere_ader.test_iter = 10
    fem_ns_sphere_ader.test_vals = [-35.000000, -35.000000, 0.000047, 31.110911]
    fem_ns_sphere_ader.unsteady   = True
    test_list.append(fem_ns_sphere_ader)

    # Unsteady cylinder
    fem_ns_unsteady_cylinder           = TestCase('fem_ns_unsteady_cylinder')
    fem_ns_unsteady_cylinder.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder.cfg_file  = "fem_unst_cylinder.cfg"
    fem_ns_unsteady_cylinder.test_iter = 11
    fem_ns_unsteady_cylinder.test_vals = [-3.558582, -3.014464, -0.038927, 1.383983]
    fem_ns_unsteady_cylinder.unsteady   = True
    test_list.append(fem_ns_unsteady_cylinder)

    # Unsteady cylinder ADER
    fem_ns_unsteady_cylinder_ader           = TestCase('fem_ns_unsteady_cylinder_ader')
    fem_ns_unsteady_cylinder_ader.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder_ader.cfg_file  = "fem_unst_cylinder_ADER.cfg"
    fem_ns_unsteady_cylinder_ader.test_iter = 11
    fem_ns_unsteady_cylinder_ader.test_vals = [-35.000000, -35.000000, -0.041003, 1.391339]
    fem_ns_unsteady_cylinder_ader.unsteady   = True
    test_list.append(fem_ns_unsteady_cylinder_ader)

    #########################
    ###    Transition     ###
    #########################

    # Schubauer-Klebanoff Natural Transition
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.test_vals    = [-8.284308, -13.240273, 0.000057, 0.007982]
    test_list.append(schubauer_klebanoff_transition)

    #####################################
    ### Cont. adj. compressible Euler ###
    #####################################

    # Dry run Cont. Adj. Euler
    contadj_naca0012_d           = TestCase('dry run Cont. Adj. Euler')
    contadj_naca0012_d.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012_d.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(contadj_naca0012_d)

    # Inviscid NACA0012
    contadj_naca0012           = TestCase('contadj_naca0012')
    contadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012.test_iter = 5
    contadj_naca0012.test_vals = [-9.748339, -15.067997, -0.726250, 0.020280]
    contadj_naca0012.tol       = 0.001
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.034680, -12.592674, -1.086100, 0.007556]
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10
    contadj_wedge.test_vals = [2.872064, -2.756210, 1010800.000000, 0.000000]
    test_list.append(contadj_wedge)

    # Inviscid fixed CL NACA0012
    contadj_fixedCL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixedCL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixedCL_naca0012.test_iter = 100
    contadj_fixedCL_naca0012.test_vals = [0.754936, -4.793625, -0.524550, -0.000227]
    test_list.append(contadj_fixedCL_naca0012)

    ###################################
    ### Cont. adj. compressible N-S ###
    ###################################

    # Dry run Cont. Adj. NS
    contadj_ns_cylinder_d           = TestCase('dry run Cont. Adj. NS')
    contadj_ns_cylinder_d.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder_d.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(contadj_ns_cylinder_d)

    # Adjoint laminar cylinder
    contadj_ns_cylinder           = TestCase('contadj_ns_cylinder')
    contadj_ns_cylinder.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder.test_iter = 20
    contadj_ns_cylinder.test_vals = [-3.665842, -9.132048, 2.056700, -0.000000]
    test_list.append(contadj_ns_cylinder)

    # Adjoint laminar naca0012 subsonic
    contadj_ns_naca0012_sub           = TestCase('contadj_ns_naca0012_sub')
    contadj_ns_naca0012_sub.cfg_dir   = "cont_adj_navierstokes/naca0012_sub"
    contadj_ns_naca0012_sub.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_sub.test_iter = 20
    contadj_ns_naca0012_sub.test_vals = [-2.743268, -8.215193, 0.518810, 0.001210]
    test_list.append(contadj_ns_naca0012_sub)

    # Adjoint laminar naca0012 transonic
    contadj_ns_naca0012_trans           = TestCase('contadj_ns_naca0012_trans')
    contadj_ns_naca0012_trans.cfg_dir   = "cont_adj_navierstokes/naca0012_trans"
    contadj_ns_naca0012_trans.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_trans.test_iter = 20
    contadj_ns_naca0012_trans.test_vals = [-1.039664, -6.575019, 1.772300, 0.012495]
    test_list.append(contadj_ns_naca0012_trans)

    #######################################################
    ### Cont. adj. compressible RANS (frozen viscosity) ###
    #######################################################

    # Adjoint turbulent NACA0012
    contadj_rans_naca0012           = TestCase('contadj_rans_naca0012')
    contadj_rans_naca0012.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012.cfg_file  = "turb_nasa.cfg"
    contadj_rans_naca0012.test_iter = 20
    contadj_rans_naca0012.test_vals = [-0.794162, -5.761722, 19.214000, -0.000000]
    test_list.append(contadj_rans_naca0012)

    # Adjoint turbulent NACA0012 with binary restarts
    contadj_rans_naca0012_bin           = TestCase('contadj_rans_naca0012_bin')
    contadj_rans_naca0012_bin.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012_bin.cfg_file  = "turb_nasa_binary.cfg"
    contadj_rans_naca0012_bin.test_iter = 18
    contadj_rans_naca0012_bin.test_vals = [-0.794169, -5.761671, 19.214000, -0.000000]
    test_list.append(contadj_rans_naca0012_bin)

    # Adjoint turbulent RAE2822
    contadj_rans_rae2822           = TestCase('contadj_rans_rae2822')
    contadj_rans_rae2822.cfg_dir   = "cont_adj_rans/rae2822"
    contadj_rans_rae2822.cfg_file  = "turb_SA_RAE2822.cfg"
    contadj_rans_rae2822.test_iter = 20
    contadj_rans_rae2822.test_vals = [-5.369688, -10.872211, -0.212470, 0.005448]
    test_list.append(contadj_rans_rae2822)

    #############################
    ### Compressibele RANS UQ ###
    #############################

    # NACA0012 1c
    turb_naca0012_1c           = TestCase('turb_naca0012_1c')
    turb_naca0012_1c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_1c.cfg_file  = "turb_NACA0012_uq_1c.cfg"
    turb_naca0012_1c.test_iter = 10
    turb_naca0012_1c.test_vals = [-4.992341, 1.135213, 0.352815, -0.084443]
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals = [-5.485201, 0.968804, 0.270412, -0.105544]
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals = [-5.584365, 0.931915, 0.247973, -0.108789]
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals = [-5.120106, 1.075360, 0.291660, -0.104760]
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals = [-5.549184, 0.946308, 0.249527, -0.112742]
    test_list.append(turb_naca0012_p1c2)

    ######################################
    ### Harmonic Balance               ###
    ######################################

    # Description of the regression test
    harmonic_balance           = TestCase('harmonic_balance')
    harmonic_balance.cfg_dir   = "harmonic_balance"
    harmonic_balance.cfg_file  = "HB.cfg"
    harmonic_balance.test_iter = 25
    harmonic_balance.test_vals = [-1.559187, 0.829575, 0.931512, 3.954440]
    test_list.append(harmonic_balance)

    # Turbulent pitching NACA 64a010 airfoil
    hb_rans_preconditioning           = TestCase('hb_rans_preconditioning')
    hb_rans_preconditioning.cfg_dir   = "harmonic_balance/hb_rans_preconditioning"
    hb_rans_preconditioning.cfg_file  = "davis.cfg"
    hb_rans_preconditioning.test_iter = 25
    hb_rans_preconditioning.test_vals = [-1.902097, 0.484069, 0.601483, 3.609005, -5.949355]
    test_list.append(hb_rans_preconditioning)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Rotating NACA 0012
    rot_naca0012           = TestCase('rot_naca0012')
    rot_naca0012.cfg_dir   = "rotating/naca0012"
    rot_naca0012.cfg_file  = "rot_NACA0012.cfg"
    rot_naca0012.test_iter = 25
    rot_naca0012.test_vals = [-2.610126, 2.922809, -0.080488, 0.002170]
    test_list.append(rot_naca0012)

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.627870, -0.164404, 0.054706, 2.545833]
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.894516, -2.469078, 1.703821, 1.670411]
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-2.560839, -1.175927, 0.062080, 1.399401, 2.220371, 1.399349, 2.218609, 0.000000]
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977498, 3.481817, -0.010657, -0.007976]
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Aeroelastic
    aeroelastic         = TestCase('aeroelastic')
    aeroelastic.cfg_dir   = "aeroelastic"
    aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    aeroelastic.test_iter = 2
    aeroelastic.test_vals = [0.074208, 0.027599, -0.001641, -0.000128]
    aeroelastic.unsteady  = True
    test_list.append(aeroelastic)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714713, -5.788301, -0.214960, 0.023758, 0.000000]
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.007635, -6.879809, 1.445300, 0.419281]
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # unsteady pitching NACA0012, Euler, Deforming
    unst_deforming_naca0012           = TestCase('unst_deforming_naca0012')
    unst_deforming_naca0012.cfg_dir   = "disc_adj_euler/naca0012_pitching_def"
    unst_deforming_naca0012.cfg_file  = "inv_NACA0012_pitching_deform.cfg"
    unst_deforming_naca0012.test_iter = 5
    unst_deforming_naca0012.test_vals = [-3.665152, -3.793306, -3.716483, -3.148336]
    unst_deforming_naca0012.unsteady  = True
    test_list.append(unst_deforming_naca0012)

    ######################################
    ### NICFD                          ###
    ######################################

    # ls89_sa
    ls89_sa           = TestCase('ls89_sa')
    ls89_sa.cfg_dir   = "nicf/LS89"
    ls89_sa.cfg_file  = "turb_SA_PR.cfg"
    ls89_sa.test_iter = 20
    ls89_sa.test_vals = [-5.050483, -13.371790, 0.174939, 0.430757]
    test_list.append(ls89_sa)

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 20
    edge_VW.test_vals = [-0.772250, 5.429879, -0.000470, 0.000000]
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 20
    edge_PPR.test_vals = [-2.126694, 4.066051, -0.000013, 0.000000]
    test_list.append(edge_PPR)


    ######################################
    ### Turbomachinery                 ###
    ######################################

    # Aachen Turbine restart
    Aachen_3D_restart           = TestCase('aachen_turbine_restart')
    Aachen_3D_restart.cfg_dir   = "turbomachinery/Aachen_turbine"
    Aachen_3D_restart.cfg_file  = "aachen_3D_MP_restart.cfg"
    Aachen_3D_restart.test_iter = 5
    Aachen_3D_restart.test_vals = [-7.701448, -8.512359, -6.014939, -6.468741, -5.801762, -4.607173, -5.551041, -5.300777, -3.804188, -5.256055, -5.765225, -3.609601, -2.229277, -2.883896, -0.563470]
    test_list.append(Aachen_3D_restart)

    # Jones APU Turbocharger restart
    Jones_tc_restart           = TestCase('jones_turbocharger_restart')
    Jones_tc_restart.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_restart.cfg_file  = "Jones_restart.cfg"
    Jones_tc_restart.test_iter = 5
    Jones_tc_restart.test_vals = [-7.308013, -5.332076, -14.895815, -9.330694, -12.071730, -6.548620, 73291.000000, 73291.000000, 0.020111, 82.896000]
    test_list.append(Jones_tc_restart)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals = [1.090052, 1.551005, -2.895066, 2.607594, -2.479704, 3.063740, 106380.000000, 106380.000000, 5.733600, 64.747000]
    test_list.append(axial_stage2D)

    # 2D transonic stator restart
    transonic_stator_restart           = TestCase('transonic_stator_restart')
    transonic_stator_restart.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_restart.cfg_file  = "transonic_stator_restart.cfg"
    transonic_stator_restart.test_iter = 20
    transonic_stator_restart.test_vals         = [-4.359175, -2.486698, -2.079382, 1.736046, -1.443172, 3.240049, -471620.000000, 94.840000, -0.055450]
    transonic_stator_restart.test_vals_aarch64 = [-4.359175, -2.486698, -2.079382, 1.736046, -1.443172, 3.240049, -471620.000000, 94.840000, -0.055450]
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

    # Dry run Multizone
    uniform_flow_d         = TestCase('dry run Multizone')
    uniform_flow_d.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow_d.cfg_file  = "uniform_NN.cfg"
    uniform_flow_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(uniform_flow_d)

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 2
    uniform_flow.test_vals = [2.000000, 0.000000, -0.230641, -13.250022]
    uniform_flow.tol       = 0.000001
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow)

   # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 2
    channel_2D.test_vals = [2.000000, 0.000000, 0.464920, 0.348054, 0.397553]
    channel_2D.timeout   = 100
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 1
    channel_3D.test_vals = [1.000000, 0.000000, 0.611989, 0.798999, 0.702438]
    channel_3D.unsteady  = True
    channel_3D.multizone = True
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [0.547322, 0.655095, 0.968232, 1.049121]
    pipe.unsteady  = True
    pipe.multizone = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [3.000000, 0.000000, 0.717067, 1.119816, 1.160327]
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 1.207119, 1.065262]
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [13.000000, 1.167903, -1.660900]
    bars_SST_2D.multizone = True
    test_list.append(bars_SST_2D)

    # Sliding mesh with incompressible flows (steady)
    slinc_steady           = TestCase('slinc_steady')
    slinc_steady.cfg_dir   = "sliding_interface/incompressible_steady"
    slinc_steady.cfg_file  = "config.cfg"
    slinc_steady.test_iter = 19
    slinc_steady.test_vals = [19.000000, -1.803732, -2.108492]
    slinc_steady.timeout   = 100
    slinc_steady.multizone = True
    test_list.append(slinc_steady)

    # Sliding mesh with incompressible flows (unsteady)
    # slinc_unsteady           = TestCase('slinc_unsteady')
    # slinc_unsteady.cfg_dir   = "sliding_interface/incompressible_unsteady"
    # slinc_unsteady.cfg_file  = "config.cfg"
    # slinc_unsteady.test_iter = 19
    # slinc_unsteady.test_vals = [-3.515218,1.930028,0.000000,0.000000] #last 4 columns
    # slinc_unsteady.timeout   = 100
    # slinc_unsteady.unsteady  = True
    # test_list.append(slinc_unsteady)

    ##########################
    ### FEA - FSI          ###
    ##########################

    # Dry run FEA
    statbeam3d_d           = TestCase('dry run FEA')
    statbeam3d_d.cfg_dir   = "fea_fsi/StatBeam_3d"
    statbeam3d_d.cfg_file  = "configBeam_3d.cfg"
    statbeam3d_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(statbeam3d_d)

    # Static beam, 3d
    statbeam3d           = TestCase('statbeam3d')
    statbeam3d.cfg_dir   = "fea_fsi/StatBeam_3d"
    statbeam3d.cfg_file  = "configBeam_3d.cfg"
    statbeam3d.test_iter = 0
    statbeam3d.test_vals         = [-6.168640, -5.939035, -6.071159, 110190]
    statbeam3d.test_vals_aarch64 = [-6.168640, -5.939035, -6.071159, 110190] #last 4 columns
    test_list.append(statbeam3d)

    # Mix elem, 3d beam, Knowles
    knowlesbeam           = TestCase('mixelemknowles')
    knowlesbeam.cfg_dir   = "fea_fsi/MixElemsKnowles"
    knowlesbeam.cfg_file  = "config.cfg"
    knowlesbeam.test_iter = 0
    knowlesbeam.test_vals         = [-14.513598, -13.577350, -28.126416, 0.445030, 9.730600]
    knowlesbeam.test_vals_aarch64 = [-14.475326, -13.54641, -28.057487, 0.44503, 9.7306] #last 5 columns
    knowlesbeam.tol       = 0.0001
    test_list.append(knowlesbeam)

    # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.unsteady  = True
    dynbeam2d.test_iter = 6
    dynbeam2d.test_vals = [-3.240015, 2.895057, -0.353146, 66127.000000]
    test_list.append(dynbeam2d)

    # # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [4.000000, 0.000000, -3.726014, -4.277767]
    fsi2d.multizone = True
    fsi2d.unsteady  = True
    test_list.append(fsi2d)

    # FSI, Static, 2D, new mesh solver
    stat_fsi           = TestCase('stat_fsi')
    stat_fsi.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi.cfg_file  = "config.cfg"
    stat_fsi.test_iter = 7
    stat_fsi.test_vals = [-3.336320, -4.991964, 0.000000, 7.000000]
    stat_fsi.multizone = True
    test_list.append(stat_fsi)

    # FSI, Static, 2D, new mesh solver, restart
    stat_fsi_restart           = TestCase('stat_fsi_restart')
    stat_fsi_restart.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi_restart.cfg_file  = "config_restart.cfg"
    stat_fsi_restart.test_iter = 1
    stat_fsi_restart.test_vals = [-3.401553, -4.672932, 0.000000, 26.000000]
    stat_fsi_restart.multizone = True
    test_list.append(stat_fsi_restart)

    # FSI, Dynamic, 2D, new mesh solver
    dyn_fsi           = TestCase('dyn_fsi')
    dyn_fsi.cfg_dir   = "fea_fsi/dyn_fsi"
    dyn_fsi.cfg_file  = "config.cfg"
    dyn_fsi.test_iter = 4
    dyn_fsi.test_vals = [-4.330728, -4.152995, 0.000000, 85.000000]
    dyn_fsi.multizone = True
    dyn_fsi.unsteady  = True
    test_list.append(dyn_fsi)

    # FSI, 2D airfoil with RBF interpolation
    airfoilRBF           = TestCase('airfoil_fsi_rbf')
    airfoilRBF.cfg_dir   = "fea_fsi/Airfoil_RBF"
    airfoilRBF.cfg_file  = "config.cfg"
    airfoilRBF.test_iter = 1

    airfoilRBF.test_vals = [1.000000, 0.086004, -3.365759]
    airfoilRBF.tol       = 0.0001
    airfoilRBF.multizone = True
    test_list.append(airfoilRBF)

    # ###############################
    # ### Radiative Heat Transfer ###
    # ###############################

    # Radiative heat transfer
    p1rad           = TestCase('p1rad')
    p1rad.cfg_dir   = "radiation/p1model"
    p1rad.cfg_file  = "configp1.cfg"
    p1rad.test_iter = 100
    p1rad.test_vals = [-7.751309, -7.923059, -2.119084, 0.091733, -47.387000]
    test_list.append(p1rad)

    # ###############################
    # ### Conjugate heat transfer ###
    # ###############################

    # Dry run CHT
    cht_incompressible_d           = TestCase('dry run CHT')
    cht_incompressible_d.cfg_dir   = "coupled_cht/incomp_2d"
    cht_incompressible_d.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible_d.command   = TestCase.Command(exec = "SU2_CFD", param = "-d")
    test_list.append(cht_incompressible_d)

    # CHT incompressible
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "coupled_cht/incomp_2d"
    cht_incompressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [-2.128826, -0.588813, -0.588813, -0.588813]
    cht_incompressible.multizone = True
    test_list.append(cht_incompressible)

     # CHT compressible
    cht_incompressible           = TestCase('cht_compressible')
    cht_incompressible.cfg_dir   = "coupled_cht/comp_2d"
    cht_incompressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [-4.256032, -0.532728, -0.532729, -0.532728]
    cht_incompressible.multizone = True
    test_list.append(cht_incompressible)

    ##############################################
    ### Method of Manufactured Solutions (MMS) ###
    ##############################################

    # FVM, compressible, laminar N-S
    mms_fvm_ns           = TestCase('mms_fvm_ns')
    mms_fvm_ns.cfg_dir   = "mms/fvm_navierstokes"
    mms_fvm_ns.cfg_file  = "lam_mms_roe.cfg"
    mms_fvm_ns.test_iter = 20
    mms_fvm_ns.test_vals = [-2.808514, 2.152654, 0.000000, 0.000000]
    mms_fvm_ns.tol       = 0.0001
    test_list.append(mms_fvm_ns)

    # FVM, incompressible, euler
    mms_fvm_inc_euler           = TestCase('mms_fvm_inc_euler')
    mms_fvm_inc_euler.cfg_dir   = "mms/fvm_incomp_euler"
    mms_fvm_inc_euler.cfg_file  = "inv_mms_jst.cfg"
    mms_fvm_inc_euler.test_iter = 20
    mms_fvm_inc_euler.test_vals = [-9.128345, -9.441741, 0.000000, 0.000000]
    mms_fvm_inc_euler.tol       = 0.0001
    test_list.append(mms_fvm_inc_euler)

    # FVM, incompressible, laminar N-S
    mms_fvm_inc_ns           = TestCase('mms_fvm_inc_ns')
    mms_fvm_inc_ns.cfg_dir   = "mms/fvm_incomp_navierstokes"
    mms_fvm_inc_ns.cfg_file  = "lam_mms_fds.cfg"
    mms_fvm_inc_ns.test_iter = 20
    mms_fvm_inc_ns.test_vals = [-7.414944, -7.631546, 0.000000, 0.000000]
    mms_fvm_inc_ns.tol       = 0.0001
    test_list.append(mms_fvm_inc_ns)

    # DG, compressible, euler
    ringleb_dg_euler           = TestCase('ringleb_dg_euler')
    ringleb_dg_euler.cfg_dir   = "mms/dg_ringleb"
    ringleb_dg_euler.cfg_file  = "ringleb_dg.cfg"
    ringleb_dg_euler.test_iter = 100
    ringleb_dg_euler.test_vals = [-5.136652, -4.724941, 0.000000, 0.000000]
    ringleb_dg_euler.tol       = 0.0001
    test_list.append(ringleb_dg_euler)

    # DG, compressible, laminar N-S
    mms_dg_ns           = TestCase('mms_dg_ns')
    mms_dg_ns.cfg_dir   = "mms/dg_navierstokes"
    mms_dg_ns.cfg_file  = "lam_mms_dg.cfg"
    mms_dg_ns.test_iter = 100
    mms_dg_ns.test_vals = [-1.845393, 3.520699, 0.000000, 0.000000]
    mms_dg_ns.tol       = 0.0001
    test_list.append(mms_dg_ns)

    # DG, compressible, laminar N-S 3D
    mms_dg_ns_3d           = TestCase('mms_dg_ns_3d')
    mms_dg_ns_3d.cfg_dir   = "mms/dg_navierstokes_3d"
    mms_dg_ns_3d.cfg_file  = "lam_mms_dg_3d.cfg"
    mms_dg_ns_3d.test_iter = 100
    mms_dg_ns_3d.test_vals = [-0.146826, 5.356413, 0.000000, 0.000000]
    mms_dg_ns_3d.tol       = 0.0001
    test_list.append(mms_dg_ns_3d)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    # set suitable defaults unless something else has been specified
    # command: "SU2_CFD"
    # timeout: 1600
    # tol:     0.00001
    for test in test_list:
        if test.command.empty():
            test.command = TestCase.Command(exec = "SU2_CFD")
        if test.timeout == 0:
            test.timeout = 1600
        if test.tol == 0.0:
            test.tol = 0.00001

    pass_list = [ test.run_test(args.tsan, args.asan) for test in test_list ]


    ######################################
    ### RUN SU2_GEO TESTS              ###
    ######################################

    # NACA0012
    naca0012_geo           = TestCase('naca0012_geo')
    naca0012_geo.cfg_dir   = "optimization_euler/steady_naca0012"
    naca0012_geo.cfg_file  = "inv_NACA0012_adv.cfg"
    naca0012_geo.test_vals = [1.000000, 62.045500, 0.120011, 0.000000]
    naca0012_geo.command   =  TestCase.Command(exec = "SU2_GEO")
    naca0012_geo.timeout   = 1600
    naca0012_geo.tol       = 0.00001
    pass_list.append(naca0012_geo.run_geo(args.tsan, args.asan))
    test_list.append(naca0012_geo)

    ######################################
    ### RUN SU2_DEF TESTS              ###
    ######################################

    # intersection prevention
    intersect_def            = TestCase('intersectionprevention')
    intersect_def.cfg_dir   = "deformation/intersection_prevention"
    intersect_def.cfg_file  = "def_intersect.cfg"
    intersect_def.test_iter = 10
    intersect_def.test_vals = [0.000112]
    intersect_def.command   =  TestCase.Command(exec = "SU2_DEF")
    intersect_def.timeout   = 1600
    intersect_def.tol       = 1e-04

    pass_list.append(intersect_def.run_def(args.tsan, args.asan))
    test_list.append(intersect_def)

    # Inviscid NACA0012 (triangles)
    naca0012_def            = TestCase('naca0012_def')
    naca0012_def.cfg_dir   = "deformation/naca0012"
    naca0012_def.cfg_file  = "def_NACA0012.cfg"
    naca0012_def.test_iter = 10
    naca0012_def.test_vals = [0.0034470]
    naca0012_def.command   =  TestCase.Command(exec = "SU2_DEF")
    naca0012_def.timeout   = 1600
    naca0012_def.tol       = 1e-06

    pass_list.append(naca0012_def.run_def(args.tsan, args.asan))
    test_list.append(naca0012_def)

    # Inviscid NACA0012 based on SURFACE_FILE input (surface_bump.dat)
    naca0012_def_file            = TestCase('naca0012_def_file')
    naca0012_def_file.cfg_dir   = "deformation/naca0012"
    naca0012_def_file.cfg_file  = "surface_file_NACA0012.cfg"
    naca0012_def_file.test_iter = 10
    naca0012_def_file.test_vals = [0.0034470]
    naca0012_def_file.command   =  TestCase.Command(exec = "SU2_DEF")
    naca0012_def_file.timeout   = 1600
    naca0012_def_file.tol       = 1e-6

    pass_list.append(naca0012_def_file.run_def(args.tsan, args.asan))
    test_list.append(naca0012_def_file)

    # RAE2822 (mixed tris + quads)
    rae2822_def            = TestCase('rae2822_def')
    rae2822_def.cfg_dir   = "deformation/rae2822"
    rae2822_def.cfg_file  = "def_RAE2822.cfg"
    rae2822_def.test_iter = 10
    rae2822_def.test_vals = [0.0000000]
    rae2822_def.command   =  TestCase.Command(exec = "SU2_DEF")
    rae2822_def.timeout   = 1600
    rae2822_def.tol       = 1e-06

    pass_list.append(rae2822_def.run_def(args.tsan, args.asan))
    test_list.append(rae2822_def)

    # Turb NACA4412 (quads, wall distance)
    naca4412_def            = TestCase('naca4412_def')
    naca4412_def.cfg_dir   = "deformation/naca4412"
    naca4412_def.cfg_file  = "def_NACA4412.cfg"
    naca4412_def.test_iter = 10
    naca4412_def.test_vals = [0.0000000]
    naca4412_def.command   =  TestCase.Command(exec = "SU2_DEF")
    naca4412_def.timeout   = 1600
    naca4412_def.tol       = 1e-06

    pass_list.append(naca4412_def.run_def(args.tsan, args.asan))
    test_list.append(naca4412_def)

    # Brick of tets (inverse volume)
    brick_tets_def            = TestCase('brick_tets_def')
    brick_tets_def.cfg_dir   = "deformation/brick_tets"
    brick_tets_def.cfg_file  = "def_brick_tets.cfg"
    brick_tets_def.test_iter = 10
    brick_tets_def.test_vals = [0.0008970]
    brick_tets_def.command   =  TestCase.Command(exec = "SU2_DEF")
    brick_tets_def.timeout   = 1600
    brick_tets_def.tol       = 1e-06

    pass_list.append(brick_tets_def.run_def(args.tsan, args.asan))
    test_list.append(brick_tets_def)

    # Brick of isotropic hexas (inverse volume)
    brick_hex_def           = TestCase('brick_hex_def')
    brick_hex_def.cfg_dir   = "deformation/brick_hex"
    brick_hex_def.cfg_file  = "def_brick_hex.cfg"
    brick_hex_def.test_iter = 10
    brick_hex_def.test_vals = [0.0002080]
    brick_hex_def.command   =  TestCase.Command(exec = "SU2_DEF")
    brick_hex_def.timeout   = 1600
    brick_hex_def.tol       = 1e-06

    pass_list.append(brick_hex_def.run_def(args.tsan, args.asan))
    test_list.append(brick_hex_def)

    # Brick with a pyramid layer (inverse volume)
    brick_pyra_def           = TestCase('brick_pyra_def')
    brick_pyra_def.cfg_dir   = "deformation/brick_pyra"
    brick_pyra_def.cfg_file  = "def_brick_pyra.cfg"
    brick_pyra_def.test_iter = 10
    brick_pyra_def.test_vals = [0.0015010]
    brick_pyra_def.command   =  TestCase.Command(exec = "SU2_DEF")
    brick_pyra_def.timeout   = 1600
    brick_pyra_def.tol       = 1e-06

    pass_list.append(brick_pyra_def.run_def(args.tsan, args.asan))
    test_list.append(brick_pyra_def)

    # Brick of isotropic prisms (inverse volume)
    brick_prism_def           = TestCase('brick_prism_def')
    brick_prism_def.cfg_dir   = "deformation/brick_prism"
    brick_prism_def.cfg_file  = "def_brick_prism.cfg"
    brick_prism_def.test_iter = 10
    brick_prism_def.test_vals = [0.0021210]
    brick_prism_def.command   =  TestCase.Command(exec = "SU2_DEF")
    brick_prism_def.timeout   = 1600
    brick_prism_def.tol       = 1e-06

    pass_list.append(brick_prism_def.run_def(args.tsan, args.asan))
    test_list.append(brick_prism_def)

    # Brick of prisms with high aspect ratio cells near the wall (wall distance)
    brick_prism_rans_def           = TestCase('brick_prism_rans_def')
    brick_prism_rans_def.cfg_dir   = "deformation/brick_prism_rans"
    brick_prism_rans_def.cfg_file  = "def_brick_prism_rans.cfg"
    brick_prism_rans_def.test_iter = 10
    brick_prism_rans_def.test_vals = [0.0000000]
    brick_prism_rans_def.command   =  TestCase.Command(exec = "SU2_DEF")
    brick_prism_rans_def.timeout   = 1600
    brick_prism_rans_def.tol       = 1e-06

    pass_list.append(brick_prism_rans_def.run_def(args.tsan, args.asan))
    test_list.append(brick_prism_rans_def)

    # Brick of hexas with high aspect ratio cells near the wall (inverse volume)
    brick_hex_rans_def           = TestCase('brick_hex_rans_def')
    brick_hex_rans_def.cfg_dir   = "deformation/brick_hex_rans"
    brick_hex_rans_def.cfg_file  = "def_brick_hex_rans.cfg"
    brick_hex_rans_def.test_iter = 10
    brick_hex_rans_def.test_vals = [0.0000000]
    brick_hex_rans_def.command   =  TestCase.Command(exec = "SU2_DEF")
    brick_hex_rans_def.timeout   = 1600
    brick_hex_rans_def.tol       = 1e-06

    pass_list.append(brick_hex_rans_def.run_def(args.tsan, args.asan))
    test_list.append(brick_hex_rans_def)

    # Cylindrical FFD test
    cylinder_ffd_def           = TestCase('cylinder_ffd_def')
    cylinder_ffd_def.cfg_dir   = "deformation/cylindrical_ffd"
    cylinder_ffd_def.cfg_file  = "def_cylindrical.cfg"
    cylinder_ffd_def.test_iter = 10
    cylinder_ffd_def.test_vals = [0.0004700]
    cylinder_ffd_def.command   =  TestCase.Command(exec = "SU2_DEF")
    cylinder_ffd_def.timeout   = 1600
    cylinder_ffd_def.tol       = 1e-06

    pass_list.append(cylinder_ffd_def.run_def(args.tsan, args.asan))
    test_list.append(cylinder_ffd_def)

    # Spherical FFD test
    sphere_ffd_def           = TestCase('sphere_ffd_def')
    sphere_ffd_def.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def.cfg_file  = "def_spherical.cfg"
    sphere_ffd_def.test_iter = 10
    sphere_ffd_def.test_vals = [0.0035670]
    sphere_ffd_def.command   =  TestCase.Command(exec = "SU2_DEF")
    sphere_ffd_def.timeout   = 1600
    sphere_ffd_def.tol       = 1e-06

    pass_list.append(sphere_ffd_def.run_def(args.tsan, args.asan))
    test_list.append(sphere_ffd_def)

    # Spherical FFD test using BSplines
    sphere_ffd_def_bspline           = TestCase('sphere_ffd_def_bspline')
    sphere_ffd_def_bspline.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def_bspline.cfg_file  = "def_spherical_bspline.cfg"
    sphere_ffd_def_bspline.test_iter = 10
    sphere_ffd_def_bspline.test_vals = [0.0020680]
    sphere_ffd_def_bspline.command   =  TestCase.Command(exec = "SU2_DEF")
    sphere_ffd_def_bspline.timeout   = 1600
    sphere_ffd_def_bspline.tol       = 1e-06

    pass_list.append(sphere_ffd_def_bspline.run_def(args.tsan, args.asan))
    test_list.append(sphere_ffd_def_bspline)

    ######################################
    ### RUN PYTHON TESTS               ###
    ######################################

    # test continuous_adjoint.py
    contadj_euler_py = TestCase('contadj_euler_py')
    contadj_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    contadj_euler_py.cfg_file  = "inv_NACA0012.cfg"
    contadj_euler_py.test_iter = 10
    contadj_euler_py.command   =  TestCase.Command(exec = "continuous_adjoint.py", param = "-f")
    contadj_euler_py.timeout   = 1600
    contadj_euler_py.reference_file = "of_grad_cd.dat.ref"
    contadj_euler_py.test_file = "of_grad_cd.dat"
    contadj_euler_py.enabled_with_asan = False
    pass_list.append(contadj_euler_py.run_filediff(args.tsan, args.asan))
    test_list.append(contadj_euler_py)

    # test shape_optimization.py
    shape_opt_euler_py           = TestCase('shape_opt_euler_py')
    shape_opt_euler_py.cfg_dir   = "optimization_euler/steady_naca0012"
    shape_opt_euler_py.cfg_file  = "inv_NACA0012_adv.cfg"
    shape_opt_euler_py.test_iter = 1
    shape_opt_euler_py.test_vals = [1.000000, 1.000000, 0.000021, 0.003643]
    shape_opt_euler_py.command   =  TestCase.Command(exec = "shape_optimization.py", param = "-g CONTINUOUS_ADJOINT -f")
    shape_opt_euler_py.timeout   = 1600
    shape_opt_euler_py.tol       = 0.00001
    shape_opt_euler_py.enabled_with_asan = False
    pass_list.append(shape_opt_euler_py.run_opt(args.tsan, args.asan))
    test_list.append(shape_opt_euler_py)

    # Multiple functionals with the continuous adjoint
    contadj_multi_py            = TestCase('contadj_multi_py')
    contadj_multi_py.cfg_dir    = "cont_adj_euler/wedge"
    contadj_multi_py.cfg_file   = "inv_wedge_ROE_multiobj.cfg"
    contadj_multi_py.test_iter  = 10
    contadj_multi_py.command    =  TestCase.Command(exec = "continuous_adjoint.py", param = "-f")
    contadj_multi_py.timeout    = 1600
    contadj_multi_py.reference_file = "of_grad_combo.dat.ref"
    contadj_multi_py.test_file  = "of_grad_combo.dat"
    contadj_multi_py.enabled_with_asan = False
    pass_list.append(contadj_multi_py.run_filediff(args.tsan, args.asan))
    test_list.append(contadj_multi_py)

    # Optimization with multiple objectives, with gradients evaluated individually
    # the difference in gradient value relative to combined case
    # is due to lack of solution file for the adjoint and small number of iterations
#    opt_multiobj_py            = TestCase('opt_multiobj_py')
#    opt_multiobj_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
#    opt_multiobj_py.cfg_file   = "inv_wedge_ROE_multiobj.cfg"
#    opt_multiobj_py.test_iter  = 1
#    opt_multiobj_py.test_vals = [1, 1, 1.084701E+02, 3.799222E+00] #last 4 columns
#    opt_multiobj_py.command    =  TestCase.Command(exec = "shape_optimization.py", param = "-g CONTINUOUS_ADJOINT -f")
#    opt_multiobj_py.timeout    = 1600
#    opt_multiobj_py.tol       = 0.00001
#    pass_list.append(opt_multiobj_py.run_opt())
#    test_list.append(opt_multiobj_py)
#
#    # test optimization, with multiple objectives and gradient evaluated as 'combo'
#    opt_multiobjcombo_py            = TestCase('opt_multiobjcombo_py')
#    opt_multiobjcombo_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
#    opt_multiobjcombo_py.cfg_file   = "inv_wedge_ROE_multiobj_combo.cfg"
#    opt_multiobjcombo_py.test_iter  = 1
#    opt_multiobjcombo_py.test_vals = [1, 1, 1.084701E+02, 3.789322E+00] #last 4 columns
#    opt_multiobjcombo_py.command    =  TestCase.Command(exec = "shape_optimization.py", param = "-g CONTINUOUS_ADJOINT -f")
#    opt_multiobjcombo_py.timeout    = 1600
#    opt_multiobjcombo_py.tol       = 0.00001
#    pass_list.append(opt_multiobjcombo_py.run_opt())
#    test_list.append(opt_multiobjcombo_py)

    # test optimization, with multiple objectives evaluated on a single surface
    opt_multiobj1surf_py            = TestCase('opt_multiobj1surf_py')
    opt_multiobj1surf_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
    opt_multiobj1surf_py.cfg_file   = "inv_wedge_ROE_multiobj_1surf.cfg"
    opt_multiobj1surf_py.test_iter  = 1
    opt_multiobj1surf_py.test_vals = [1.000000, 1.000000, 36.122920, 4.611743]
    opt_multiobj1surf_py.command    =  TestCase.Command(exec = "shape_optimization.py", param = "-g CONTINUOUS_ADJOINT -f")
    opt_multiobj1surf_py.timeout    = 1600
    opt_multiobj1surf_py.tol       = 0.00001
    opt_multiobj1surf_py.enabled_with_asan = False
    pass_list.append(opt_multiobj1surf_py.run_opt(args.tsan, args.asan))
    test_list.append(opt_multiobj1surf_py)

    # test optimization, with a single objective evaluated on multiple surfaces
    opt_2surf1obj_py            = TestCase('opt_2surf1obj_py')
    opt_2surf1obj_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
    opt_2surf1obj_py.cfg_file   = "inv_wedge_ROE_2surf_1obj.cfg"
    opt_2surf1obj_py.test_iter  = 1
    opt_2surf1obj_py.test_vals = [1.000000, 1.000000, 2.005099, 0.000383]
    opt_2surf1obj_py.command    =  TestCase.Command(exec = "shape_optimization.py", param = "-g CONTINUOUS_ADJOINT -f")
    opt_2surf1obj_py.timeout    = 1600
    opt_2surf1obj_py.tol       = 0.00001
    opt_2surf1obj_py.enabled_with_asan = False
    pass_list.append(opt_2surf1obj_py.run_opt(args.tsan, args.asan))
    test_list.append(opt_2surf1obj_py)

    ##########################
    ###   Python wrapper   ###
    ##########################

    # NACA0012
    pywrapper_naca0012           = TestCase('pywrapper_naca0012')
    pywrapper_naca0012.cfg_dir   = "euler/naca0012"
    pywrapper_naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    pywrapper_naca0012.test_iter = 20
    pywrapper_naca0012.test_vals = [-4.766184, -4.287722, 0.326688, 0.022661]
    pywrapper_naca0012.command   =  TestCase.Command(exec = "SU2_CFD.py", param = "-f")
    pywrapper_naca0012.timeout   = 1600
    pywrapper_naca0012.tol       = 0.00001
    pywrapper_naca0012.enabled_with_asan = False
    test_list.append(pywrapper_naca0012)
    pass_list.append(pywrapper_naca0012.run_test(args.tsan, args.asan))

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    pywrapper_turb_naca0012_sst           = TestCase('pywrapper_turb_naca0012_sst')
    pywrapper_turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    pywrapper_turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    pywrapper_turb_naca0012_sst.test_iter = 10
    pywrapper_turb_naca0012_sst.test_vals         = [-12.229889, -14.434883, -6.037177, 1.047444, 0.019214, -2.108366, 0.000000]
    pywrapper_turb_naca0012_sst.test_vals_aarch64 = [-12.229889, -14.434883, -6.037177, 1.047444, 0.019214, -2.108366, 0.000000]
    pywrapper_turb_naca0012_sst.command   =  TestCase.Command(exec = "SU2_CFD.py", param = "-f")
    pywrapper_turb_naca0012_sst.timeout   = 3200
    pywrapper_turb_naca0012_sst.tol       = 0.00001
    pywrapper_turb_naca0012_sst.enabled_with_asan = False
    test_list.append(pywrapper_turb_naca0012_sst)
    pass_list.append(pywrapper_turb_naca0012_sst.run_test(args.tsan, args.asan))

    # Square cylinder
    pywrapper_square_cylinder           = TestCase('pywrapper_square_cylinder')
    pywrapper_square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    pywrapper_square_cylinder.cfg_file  = "turb_square.cfg"
    pywrapper_square_cylinder.test_iter = 3
    pywrapper_square_cylinder.test_vals = [-2.560839, -1.175927, 0.062080, 1.399401, 2.220371, 1.399349, 2.218609, 0.000000]
    pywrapper_square_cylinder.command   =  TestCase.Command(exec = "SU2_CFD.py", param = "-f")
    pywrapper_square_cylinder.timeout   = 1600
    pywrapper_square_cylinder.tol       = 0.00001
    pywrapper_square_cylinder.unsteady  = True
    pywrapper_square_cylinder.enabled_with_asan = False
    test_list.append(pywrapper_square_cylinder)
    pass_list.append(pywrapper_square_cylinder.run_test(args.tsan, args.asan))

    # Aeroelastic
    pywrapper_aeroelastic         = TestCase('pywrapper_aeroelastic')
    pywrapper_aeroelastic.cfg_dir   = "aeroelastic"
    pywrapper_aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    pywrapper_aeroelastic.test_iter = 2
    pywrapper_aeroelastic.test_vals = [0.074208, 0.027599, -0.001641, -0.000128]
    pywrapper_aeroelastic.command   =  TestCase.Command(exec = "SU2_CFD.py", param = "-f")
    pywrapper_aeroelastic.timeout   = 1600
    pywrapper_aeroelastic.tol       = 0.00001
    pywrapper_aeroelastic.unsteady  = True
    pywrapper_aeroelastic.enabled_with_asan = False
    test_list.append(pywrapper_aeroelastic)
    pass_list.append(pywrapper_aeroelastic.run_test(args.tsan, args.asan))

    # FSI, 2d
    pywrapper_fsi2d           = TestCase('pywrapper_fsi2d')
    pywrapper_fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    pywrapper_fsi2d.cfg_file  = "configFSI.cfg"
    pywrapper_fsi2d.test_iter = 4
    pywrapper_fsi2d.test_vals = [4.000000, 0.000000, -3.726014, -4.277767]
    pywrapper_fsi2d.command   =  TestCase.Command(exec = "SU2_CFD.py", param = "--nZone 2 --fsi True -f")
    pywrapper_fsi2d.unsteady  = True
    pywrapper_fsi2d.multizone   = True
    pywrapper_fsi2d.timeout   = 1600
    pywrapper_fsi2d.tol       = 0.00001
    pywrapper_fsi2d.enabled_with_asan = False
    test_list.append(pywrapper_fsi2d)
    pass_list.append(pywrapper_fsi2d.run_test(args.tsan, args.asan))

    # Unsteady CHT
    pywrapper_unsteadyCHT               = TestCase('pywrapper_unsteadyCHT')
    pywrapper_unsteadyCHT.cfg_dir       = "py_wrapper/flatPlate_unsteady_CHT"
    pywrapper_unsteadyCHT.cfg_file      = "unsteady_CHT_FlatPlate_Conf.cfg"
    pywrapper_unsteadyCHT.test_iter     = 5
    pywrapper_unsteadyCHT.test_vals     = [-1.614167, 2.247359, 0.000771, 0.172997]
    pywrapper_unsteadyCHT.command       =  TestCase.Command(exec = "python", param = "launch_unsteady_CHT_FlatPlate.py -f")
    pywrapper_unsteadyCHT.timeout       = 1600
    pywrapper_unsteadyCHT.tol           = 0.00001
    pywrapper_unsteadyCHT.unsteady      = True
    pywrapper_unsteadyCHT.enabled_with_asan = False
    test_list.append(pywrapper_unsteadyCHT)
    pass_list.append(pywrapper_unsteadyCHT.run_test(args.tsan, args.asan))

    # Rigid motion
    pywrapper_rigidMotion               = TestCase('pywrapper_rigidMotion')
    pywrapper_rigidMotion.cfg_dir       = "py_wrapper/flatPlate_rigidMotion"
    pywrapper_rigidMotion.cfg_file      = "flatPlate_rigidMotion_Conf.cfg"
    pywrapper_rigidMotion.test_iter     = 5
    pywrapper_rigidMotion.test_vals     = [-1.614166, 2.243100, 0.350208, 0.089498]
    pywrapper_rigidMotion.command       = TestCase.Command(exec = "python", param = "launch_flatPlate_rigidMotion.py -f")
    pywrapper_rigidMotion.timeout       = 1600
    pywrapper_rigidMotion.tol           = 0.00001
    pywrapper_rigidMotion.unsteady      = True
    pywrapper_rigidMotion.enabled_with_asan = False
    test_list.append(pywrapper_rigidMotion)
    pass_list.append(pywrapper_rigidMotion.run_test(args.tsan, args.asan))

    # Custom inlet
    pywrapper_custom_inlet = TestCase('pywrapper_custom_inlet')
    pywrapper_custom_inlet.cfg_dir = "py_wrapper/custom_inlet"
    pywrapper_custom_inlet.cfg_file = "lam_flatplate.cfg"
    pywrapper_custom_inlet.test_iter = 20
    pywrapper_custom_inlet.test_vals = [-4.123206, -1.543215, -3.735006, 1.339481, -0.793478, 0.161210, -0.007009, 0.513560, -0.520570]
    pywrapper_custom_inlet.command = TestCase.Command(exec = "python", param = "run.py")
    pywrapper_custom_inlet.timeout = 1600
    pywrapper_custom_inlet.tol = 0.0001
    pywrapper_custom_inlet.unsteady = True
    pywrapper_custom_inlet.enabled_with_asan = False
    test_list.append(pywrapper_custom_inlet)
    pass_list.append(pywrapper_custom_inlet.run_test(args.tsan, args.asan))

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
