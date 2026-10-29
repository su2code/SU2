#!/usr/bin/env python

## \file serial_regression.py
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

from __future__ import print_function, division, absolute_import
import os
import subprocess
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
    thermalbath.test_vals = [0.945997, 0.945997, -11.897761, -11.867806, -32.000000, 10.013239]
    test_list.append(thermalbath)

    # Adiabatic frozen thermal bath
    thermalbath_frozen           = TestCase('thermalbath_frozen')
    thermalbath_frozen.cfg_dir   = "nonequilibrium/thermalbath/frozen"
    thermalbath_frozen.cfg_file  = "thermalbath_frozen.cfg"
    thermalbath_frozen.test_iter = 10
    thermalbath_frozen.test_vals = [-32.000000, -32.000000, -11.810719, -11.867806, -32.000000, 10.013545]
    test_list.append(thermalbath_frozen)

    # Inviscid single wedge, implicit
    invwedge = TestCase('invwedge')
    invwedge.cfg_dir = "nonequilibrium/invwedge"
    invwedge.cfg_file = "invwedge_ausm.cfg"
    invwedge.test_iter = 10
    invwedge.test_vals = [-1.085516, -1.610279, -18.299901, -18.627313, -18.573325, 2.230829, 1.850866, 5.274272, 0.833227]
    invwedge.test_vals_aarch64 = [-1.085516, -1.610279, -18.299901, -18.627313, -18.573325, 2.230829, 1.850866, 5.274272, 0.833227]
    test_list.append(invwedge)

    # Viscous single cone - axisymmetric
    visc_cone = TestCase('visc_cone')
    visc_cone.cfg_dir = "nonequilibrium/visc_wedge"
    visc_cone.cfg_file = "axi_visccone.cfg"
    visc_cone.test_iter = 10
    visc_cone.test_vals = [-5.275039, -5.799179, -20.180464, -20.157030, -20.112497, 1.156763, -3.217449, -0.013710, 0.092250, 27878.000000]
    visc_cone.test_vals_aarch64 = [-5.275039, -5.799179, -20.180464, -20.157030, -20.112496, 1.156764, -3.217451, -0.013710, 0.092250, 27878.000000]
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
    channel.test_vals = [-1.906687, 3.625653, -0.011054, 0.205254, 1.000000]
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.031446, -3.524530, 0.273056, 0.009703]
    test_list.append(naca0012)

    # NACA0012 - FMG test
    naca0012_FMG           = TestCase('naca0012_FMG')
    naca0012_FMG.cfg_dir   = "euler/naca0012"
    naca0012_FMG.cfg_file  = "inv_NACA0012.cfg"
    naca0012_FMG.test_iter = 20
    naca0012_FMG.test_vals = [-4.160287, -3.553906, 0.296154, 0.016890]
    test_list.append(naca0012_FMG)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-4.350003, 1.379848, -0.249526, 0.043953]
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-4.825509, -4.205271, 0.277086, 0.006619]
    oneram6.timeout   = 9600
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-3.886939, 1.644825, 0.300989, 0.019478]
    test_list.append(fixedCL_naca0012)

    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals = [-1.311796, 4.165508, 0.000092, 0.113306]
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
    bluntbody.test_vals = [0.565185, 6.945875, 0.000007, 1.970309]
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
    flatplate.test_vals = [-5.623336, -0.144687, 0.001202, 0.028141, 2.361500, -2.333400, 0.000000, 0.000000]
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-8.725171, -3.234867, -0.028017, 1.605526, 0.000000]
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-8.698687, -3.236566, 0.401347, 1.754864, 0.000000]
    test_list.append(cylinder_lowmach)

    # Laminar NACA0012 with linelet preconditioner (see #2877)
    lam_naca0012_linelet = TestCase('lam_naca0012_linelet')
    lam_naca0012_linelet.cfg_dir = "navierstokes/naca0012"
    lam_naca0012_linelet.cfg_file = "lam_NACA0012_LINELET.cfg"
    lam_naca0012_linelet.test_iter = 20
    lam_naca0012_linelet.test_vals = [-5.345007, -2.936878, -2.943454, 0.096356]
    test_list.append(lam_naca0012_linelet)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-7.597002, -6.994814, 0.012154, 13.596296, 0]
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals         = [-12.003115, -7.626023, -0.000000, 2.089953]
    poiseuille_profile.test_vals_aarch64 = [-12.009012, -7.262299, -0.000000, 2.089953] #last 4 columns
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
    rae2822_sa.test_vals = [-3.117674, -5.153287, 0.777070, 0.019412, 1000.000000]
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-1.916371, -1.461911, 5.894621, 0.611793, 0.016873, 100.000000]
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.464309, 5.852571, 0.490309, 0.042519]
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-5.286997, -6.483481, -0.187358, 0.003712, 10.000000, -2.181754]
    test_list.append(turb_flatplate)

    # FLAT PLATE, WALL FUNCTIONS, COMPRESSIBLE SST
    turb_wallfunction_flatplate_sst           = TestCase('turb_sst_wallfunction_flatplate')
    turb_wallfunction_flatplate_sst.cfg_dir   = "wallfunctions/flatplate/compressible_SST"
    turb_wallfunction_flatplate_sst.cfg_file  = "turb_SST_flatplate.cfg"
    turb_wallfunction_flatplate_sst.test_iter = 10
    turb_wallfunction_flatplate_sst.test_vals = [-4.036300, -1.916264, -1.821790, 1.443878, -1.579554, 1.521535, 10.000000, -2.351460, 0.030011, 0.002409]
    test_list.append(turb_wallfunction_flatplate_sst)

    # FLAT PLATE, ROUGHNESS BC WILCOX2006 SST
    turb_flatplate_sst_roughBCWilcox2006           = TestCase('turb_sst_flatplate_roughBCWilcox2006')
    turb_flatplate_sst_roughBCWilcox2006.cfg_dir   = "rans/flatplate/roughness/bc_wilcox2006"
    turb_flatplate_sst_roughBCWilcox2006.cfg_file  = "turb_SST_flatplate_roughBCWilcox2006.cfg"
    turb_flatplate_sst_roughBCWilcox2006.test_iter = 10
    turb_flatplate_sst_roughBCWilcox2006.test_vals = [-5.117900, -2.534224, -2.904279, 0.381710, -3.100346, 1.180161, -0.188797, 0.004029]
    test_list.append(turb_flatplate_sst_roughBCWilcox2006)

    # FLAT PLATE, WALL FUNCTIONS, COMPRESSIBLE SA
    turb_wallfunction_flatplate_sa           = TestCase('turb_sa_wallfunction_flatplate')
    turb_wallfunction_flatplate_sa.cfg_dir   = "wallfunctions/flatplate/compressible_SA"
    turb_wallfunction_flatplate_sa.cfg_file  = "turb_SA_flatplate.cfg"
    turb_wallfunction_flatplate_sa.test_iter = 10
    turb_wallfunction_flatplate_sa.test_vals = [-3.580017, -1.860982, -1.506109, 1.898317, -5.258602, 10.000000, -2.557123, -0.295793, 0.002533]
    test_list.append(turb_wallfunction_flatplate_sa)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.419977, -6.631993, 0.238577, 0.159694, 0.000000]
    turb_oneram6.timeout   = 3200
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D results for finest grid: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 5
    turb_naca0012_sa.test_vals = [-12.037319, -16.384159, 1.080346, 0.018385, 20.000000, -3.456939, 20.000000, -4.641254, 0.000000]
    turb_naca0012_sa.test_vals_aarch64 = [-12.037297, -16.384158, 1.080346, 0.018385, 20.000000, -3.455886, 20.000000, -4.641247, 0.000000]
    turb_naca0012_sa.timeout   = 3200
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-12.094366, -15.251082, -5.906366, 1.070413, 0.015775, -3.178304, 0.000000]
    turb_naca0012_sst.test_vals_aarch64 = [-12.076068, -15.246740, -5.861280, 1.070036, 0.015841, -3.297854, 0.000000]
    turb_naca0012_sst.timeout   = 3200
    test_list.append(turb_naca0012_sst)

    # E387 transitional SST+LM tutorial config, re-run here as a sanitizer-only probe.
    # Covers the density gradient not being available for MUSCL_TURB=YES with a flow scheme
    # that does not store that gradient.
    tutorial_trans_e387_sst_asan                  = TestCase('tutorial_trans_e387_sst_asan')
    tutorial_trans_e387_sst_asan.cfg_dir          = "../Tutorials/compressible_flow/Transitional_Airfoil/Langtry_and_Menter/E387"
    tutorial_trans_e387_sst_asan.cfg_file         = "transitional_SST_LM_model_ConfigFile.cfg"
    tutorial_trans_e387_sst_asan.test_iter        = 2
    tutorial_trans_e387_sst_asan.test_vals        = [-6.418119, -4.827573, -2.220229, 3.029787, 3.123846, 5.000000, -5.604338]
    tutorial_trans_e387_sst_asan.timeout          = 1600
    tutorial_trans_e387_sst_asan.no_restart       = True
    tutorial_trans_e387_sst_asan.enabled_with_regular = False
    test_list.append(tutorial_trans_e387_sst_asan)

    # NACA0012 (SST V2003m, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_2003m           = TestCase('turb_naca0012_sst_2003m')
    turb_naca0012_sst_2003m.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_2003m.cfg_file  = "turb_NACA0012_sst_2003m.cfg"
    turb_naca0012_sst_2003m.test_iter = 10
    turb_naca0012_sst_2003m.test_vals = [-7.129622, -10.410728, -3.699440, 1.061671, 0.019047, -2.435494, 0.000000]
    turb_naca0012_sst_2003m.timeout   = 3200
    test_list.append(turb_naca0012_sst_2003m)

    # NACA0012 (SST_SUST, FUN3D results for finest grid: CL=1.0840, CD=0.01253) restart
    turb_naca0012_sst_sust_restart           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust_restart.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust_restart.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust_restart.test_iter = 10
    turb_naca0012_sst_sust_restart.test_vals = [-12.080536, -14.837169, -5.733461, 1.000893, 0.019109, -2.633948]
    turb_naca0012_sst_sust_restart.test_vals_aarch64 = [-12.074189, -14.836725, -5.732398, 1.000050, 0.019144, -3.315560]
    turb_naca0012_sst_sust_restart.timeout   = 3200
    test_list.append(turb_naca0012_sst_sust_restart)

    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [-5.206618, -10.436758, 0.774095, 1.021995, 0.040553, -3.477597]
    turb_naca0012_sst_fixedvalues.timeout   = 3200
    test_list.append(turb_naca0012_sst_fixedvalues)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389724, -8.410479, 0.000048, 0.056344]
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
    axi_rans_air_nozzle_restart.test_vals = [-5.121517, 0.452804, -2.910659, 1.662143, 0.000000]
