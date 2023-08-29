#!/usr/bin/env python

## \file hybrid_regression.py
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
    channel.test_vals = [-2.667328, 2.797437, 0.018714, 0.006906]
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.023999, -3.515034, 0.339426, 0.022217]
    test_list.append(naca0012)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-0.942862, 4.784581, -0.208106, 0.036665]
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [0.281703, 0.011821]
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-7.374806, -1.872330, 0.300000, 0.019471]
    test_list.append(fixedCL_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.540010, 6.916656, 0.000027, 1.869004]
    test_list.append(bluntbody)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 100
    flatplate.test_vals         = [-9.154123, -3.663185, 0.001112, 0.036277, 2.361500, -2.325300, -2.278800, -2.278800]
    flatplate.test_vals_aarch64 = [-9.154130, -3.663197, 0.001112, 0.036277, 2.361500, -2.325300, -2.278800, -2.278800]
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals         = [-6.765429, -1.297425, 0.019571, 0.310233, 0.123270]
    cylinder.test_vals_aarch64 = [-6.765429, -1.297425, 0.019571, 0.310231, 0.123270]
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals         = [-6.850130, -1.388096, -0.056036, 108.140809, 0.007988]
    cylinder_lowmach.test_vals_aarch64 = [-6.850130, -1.388096, -0.056036, 108.140813, 0.007988]
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals         = [-5.048283, 0.650813, 0.008713, 13.677671, -2.054800]
    poiseuille.test_vals_aarch64 = [-5.048282, 0.650814, 0.008713, 13.677691, -2.054800]
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals         = [-12.494728, -7.712546, -0.000000, 2.085796]
    poiseuille_profile.test_vals_aarch64 = [-12.494717, -7.711274, -0.000000, 2.085796]
    test_list.append(poiseuille_profile)

    # 2D Rotational Periodic
    periodic2d           = TestCase('periodic2d')
    periodic2d.cfg_dir   = "navierstokes/periodic2D"
    periodic2d.cfg_file  = "config.cfg"
    periodic2d.test_iter = 1400
    periodic2d.test_vals         = [-10.818509, -8.363386, -8.287481, -5.334812, -1.087925, -2945.200000]
    periodic2d.test_vals_aarch64 = [-10.818510, -8.363388, -8.287480, -5.334814, -1.087922, -2945.2]
    test_list.append(periodic2d)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.020123, -5.269324, 0.807147, 0.060499, -80602.000000]
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510363, 4.872736, 0.815617, 0.060920, -73391.000000]
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.535141, 4.872736, 0.815617, 0.060920]
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.157358, -6.736289, -0.176258, 0.057431]
    test_list.append(turb_flatplate)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.388836, -6.689426, 0.230320, 0.157640, -32539.000000]
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 5
    turb_naca0012_sa.test_vals         = [-10.451625, -13.859808, 1.057622, 0.022916, 20.000000, -1.358306, 20.000000, -2.512316, -44.540000]
    turb_naca0012_sa.test_vals_aarch64 = [-10.451625, -13.859809, 1.057622, 0.022916, 20.000000, -1.358307, 20.000000, -2.512316, -44.540000]
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals         = [-12.216989, -14.439370, -7.108639, 1.050109, 0.019148, -1.483537, -38.510000]
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals         = [-12.148204, -14.755890, -6.342986, 1.001189, 0.019384, -1.432333]
    test_list.append(turb_naca0012_sst_sust)

    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals         = [-5.192598, -10.042489, -1.617721, 1.022073, 0.040198, -2.381138]
    turb_naca0012_sst_fixedvalues.test_vals_aarch64 = [-5.192599, -10.042490, -1.617720, 1.022073, 0.040198, -2.381138]
    test_list.append(turb_naca0012_sst_fixedvalues)

    # NACA0012 (SST, explicit Euler for flow and turbulence equations)
    turb_naca0012_sst_expliciteuler           = TestCase('turb_naca0012_sst_expliciteuler')
    turb_naca0012_sst_expliciteuler.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_expliciteuler.cfg_file  = "turb_NACA0012_sst_expliciteuler.cfg"
    turb_naca0012_sst_expliciteuler.test_iter = 10
    turb_naca0012_sst_expliciteuler.test_vals = [-3.532289, -3.157766, 3.364024, 1.122901, 0.500798, -float("inf")]
    test_list.append(turb_naca0012_sst_expliciteuler)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389575, -8.409251, 0.000048, 0.056329]
    test_list.append(propeller)

    #######################################
    ### Axisymmetric Compressible RANS  ###
    #######################################

    # Axisymmetric air nozzle (transonic)
    axi_rans_air_nozzle_restart           = TestCase('axi_rans_air_nozzle_restart')
    axi_rans_air_nozzle_restart.cfg_dir   = "axisymmetric_rans/air_nozzle"
    axi_rans_air_nozzle_restart.cfg_file  = "air_nozzle_restart.cfg"
    axi_rans_air_nozzle_restart.test_iter = 10
    axi_rans_air_nozzle_restart.test_vals         = [-12.060634, -6.992273, -8.808757, -4.097753, -2019.800000]
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
    turb_naca0012_sst_restart_mg.test_vals = [-7.644702, -7.433663, -1.918163, -0.000003, 0.079111]
    test_list.append(turb_naca0012_sst_restart_mg)

    #############################
    ### Compressibele RANS UQ ###
    #############################

    # NACA0012 1c
    turb_naca0012_1c           = TestCase('turb_naca0012_1c')
    turb_naca0012_1c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_1c.cfg_file  = "turb_NACA0012_uq_1c.cfg"
    turb_naca0012_1c.test_iter = 10
    turb_naca0012_1c.test_vals         = [-4.980878, 1.138865, 0.247731, -0.117335]
    turb_naca0012_1c.test_vals_aarch64 = [-4.981105, 1.138873, 0.248013, -0.117248]
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals         = [-5.483313, 0.968731, 0.215420, -0.124767]
    turb_naca0012_2c.test_vals_aarch64 = [-5.483345, 0.968720, 0.214914, -0.124932]
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals         = [-5.584300, 0.931293, 0.207446, -0.125692]
    turb_naca0012_3c.test_vals_aarch64 = [-5.584300, 0.931293, 0.207447, -0.125691]
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals         = [-5.132363, 1.075634, 0.337251, -0.082829]
    turb_naca0012_p1c1.test_vals_aarch64 = [-5.132358, 1.075658, 0.337268, -0.082827]
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals         = [-5.554392, 0.943703, 0.229483, -0.121058]
    turb_naca0012_p1c2.test_vals_aarch64 = [-5.554425, 0.943683, 0.229293, -0.121125]
    test_list.append(turb_naca0012_p1c2)

    ######################################
    ### Harmonic Balance               ###
    ######################################

    # Description of the regression test
    harmonic_balance           = TestCase('harmonic_balance')
    harmonic_balance.cfg_dir   = "harmonic_balance"
    harmonic_balance.cfg_file  = "HB.cfg"
    harmonic_balance.test_iter = 25
    harmonic_balance.test_vals = [-1.589740, 0.790168, 0.937052, 3.922578]
    test_list.append(harmonic_balance)

    # Turbulent pitching NACA 64a010 airfoil
    hb_rans_preconditioning           = TestCase('hb_rans_preconditioning')
    hb_rans_preconditioning.cfg_dir   = "harmonic_balance/hb_rans_preconditioning"
    hb_rans_preconditioning.cfg_file  = "davis.cfg"
    hb_rans_preconditioning.test_iter = 25
    hb_rans_preconditioning.test_vals = [-1.902111, 0.484080, 0.601469, 3.608991, -5.949373]
    test_list.append(hb_rans_preconditioning)

    #############################
    ### Incompressible Euler  ###
    #############################

    # NACA0012 Hydrofoil
    inc_euler_naca0012           = TestCase('inc_euler_naca0012')
    inc_euler_naca0012.cfg_dir   = "incomp_euler/naca0012"
    inc_euler_naca0012.cfg_file  = "incomp_NACA0012.cfg"
    inc_euler_naca0012.test_iter = 20
    inc_euler_naca0012.test_vals = [-4.858287, -3.810487, 0.491850, 0.007002]
    test_list.append(inc_euler_naca0012)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals         = [-5.971249, -4.910844, -0.000196, 0.121635]
    inc_nozzle.test_vals_aarch64 = [-5.971248, -4.910844, -0.000196, 0.121635]
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
    inc_poly_cylinder.test_vals         = [-7.851512, -2.093420, 0.029974, 1.921595, -175.300000]
    inc_poly_cylinder.test_vals_aarch64 = [-7.851510, -2.093419, 0.029974, 1.921595, -175.300000]
    test_list.append(inc_poly_cylinder)

    # X-coarse laminar bend as a mixed element CGNS test
    inc_lam_bend          = TestCase('inc_lam_bend')
    inc_lam_bend.cfg_dir   = "incomp_navierstokes/bend"
    inc_lam_bend.cfg_file  = "lam_bend.cfg"
    inc_lam_bend.test_iter = 10
    inc_lam_bend.test_vals         = [-3.437996, -3.086189, -0.015600, 1.142212]
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
    inc_turb_naca0012.test_vals = [-4.788405, -11.039465, 0.000008, 0.309490]
    test_list.append(inc_turb_naca0012)

    # NACA0012, SST_SUST
    inc_turb_naca0012_sst_sust           = TestCase('inc_turb_naca0012_sst_sust')
    inc_turb_naca0012_sst_sust.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_sst_sust.cfg_file  = "naca0012_SST_SUST.cfg"
    inc_turb_naca0012_sst_sust.test_iter = 20
    inc_turb_naca0012_sst_sust.test_vals = [-7.270637, 0.018416, 0.000004, 0.307678]
    test_list.append(inc_turb_naca0012_sst_sust)

    # Weakly coupled heat equation
    inc_weakly_coupled = TestCase('inc_weakly_coupled')
    inc_weakly_coupled.cfg_dir = "disc_adj_heat"
    inc_weakly_coupled.cfg_file = "primal.cfg"
    inc_weakly_coupled.test_iter = 10
    inc_weakly_coupled.test_vals = [-16.498562, -15.335587, -15.570565, -12.867787, -18.185331, -13.071436, 5.545800]
    test_list.append(inc_weakly_coupled)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.627934, -0.164469, 0.052000, 2.547063]
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals         = [-8.001291, -2.607959, 1.501321, 1.488559]
    spinning_cylinder.test_vals_aarch64 = [-8.001291, -2.607959, 1.501321, 1.488559]
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-2.557949, -1.173575, 0.058030, 1.399794, 2.220402, 1.399748, 2.218603, -0.453270]
    square_cylinder.test_vals_aarch64 = [-2.557902, -1.173574, 0.058050, 1.399794, 2.220402, 1.399748, 2.218604, -0.453270]
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977520, 3.481804, -0.012402, -0.007454]
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Cosine gust in z-direction
    cosine_gust           = TestCase('cosine_gust_zdir')
    cosine_gust.cfg_dir   = "gust"
    cosine_gust.cfg_file  = "cosine_gust_zdir.cfg"
    cosine_gust.test_iter = 79
    cosine_gust.test_vals = [-2.418813, 0.004650, -0.001878, -0.000637, -0.000271]
    cosine_gust.unsteady  = True
    cosine_gust.enabled_with_tsan = False
    test_list.append(cosine_gust)

    # Gust with mesh deformation
    gust_mesh_defo           = TestCase('gust_with_mesh_deformation')
    gust_mesh_defo.cfg_dir   = "gust"
    gust_mesh_defo.cfg_file  = "gust_with_mesh_deformation.cfg"
    gust_mesh_defo.test_iter = 6
    gust_mesh_defo.test_vals = [-1.844778, 0.000846, -0.000408]
    gust_mesh_defo.unsteady  = True
    gust_mesh_defo.enabled_with_tsan = False
    test_list.append(gust_mesh_defo)

    # Aeroelastic
    aeroelastic           = TestCase('aeroelastic')
    aeroelastic.cfg_dir   = "aeroelastic"
    aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    aeroelastic.test_iter = 2
    aeroelastic.test_vals         = [0.074432, 0.033108, -0.001650, -0.000127]
    aeroelastic.test_vals_aarch64 = [0.074836, 0.033102, -0.001650, -0.000127]
    aeroelastic.unsteady  = True
    aeroelastic.enabled_on_cpu_arch = ["x86_64"] # Requires AVX-capable architecture
    aeroelastic.enabled_with_tsan = False
    test_list.append(aeroelastic)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714757, -5.882779, -0.215005, 0.023783, -618.130000]
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.008629, -6.888963, 1.435186, 0.433529]
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # unsteady pitching NACA0012, Euler, Deforming
    unst_deforming_naca0012           = TestCase('unst_deforming_naca0012')
    unst_deforming_naca0012.cfg_dir   = "disc_adj_euler/naca0012_pitching_def"
    unst_deforming_naca0012.cfg_file  = "inv_NACA0012_pitching_deform.cfg"
    unst_deforming_naca0012.test_iter = 5
    unst_deforming_naca0012.test_vals = [-3.665126, -3.793497, -3.716496, -3.148290]
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
    edge_VW.test_iter = 100
    edge_VW.test_vals = [-5.040246, 1.124521, -0.000009, 0.000000]
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 100
    edge_PPR.test_vals         = [-5.401604, 0.738202, -0.000035, 0.000000]
    edge_PPR.test_vals_aarch64 = [-5.401642, 0.738164, -0.000035, 0.000000]
    test_list.append(edge_PPR)

    ######################################
    ### Turbomachinery                 ###
    ######################################

	# Jones APU Turbocharger restart
    Jones_tc_restart           = TestCase('jones_turbocharger_restart')
    Jones_tc_restart.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_restart.cfg_file  = "Jones_restart.cfg"
    Jones_tc_restart.test_iter = 5
    Jones_tc_restart.test_vals = [-6.604542, -2.792279, -14.328530, -8.769313, -11.371439, -5.845633, 73273.000000, 73273.000000]
    test_list.append(Jones_tc_restart)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals         = [0.974788, 1.534326, -2.897693, 2.599374, -2.418314, 3.087291, 106380.000000, 106380.000000]
    test_list.append(axial_stage2D)

    # 2D transonic stator restart
    transonic_stator_restart           = TestCase('transonic_stator_restart')
    transonic_stator_restart.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_restart.cfg_file  = "transonic_stator_restart.cfg"
    transonic_stator_restart.test_iter = 20
    transonic_stator_restart.test_vals         = [-5.346503, -3.510597, -3.156830, 0.752511, -3.834382, 1.816610, -471690]
    test_list.append(transonic_stator_restart)

    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 5
    uniform_flow.test_vals = [5.000000, 0.000000, -0.188748, -10.631533]
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow)

    # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 2
    channel_2D.test_vals = [2.000000, 0.000000, 0.397975, 0.352765, 0.405420]
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 2
    channel_3D.test_vals         = [2.000000, 0.000000, 0.620170, 0.505179, 0.415317]
    channel_3D.test_vals_aarch64 = [2.000000, 0.000000, 0.620189, 0.505311, 0.415246]
    channel_3D.unsteady  = True
    channel_3D.multizone = True
    channel_3D.enabled_with_tsan = False
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [0.150024, 0.491949, 0.677759, 0.963991, 1.006947]
    pipe.unsteady  = True
    pipe.multizone = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [3.000000, 0.000000, 0.777568, 1.134807, 1.224137]
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone  = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 1.214344, 1.663912]
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone  = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [13.000000, -0.773513, -1.700012]
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
    statbeam3d.test_vals         = [-2.378370, -1.585252, -2.028505, 64359.000000]
    statbeam3d.test_vals_aarch64 = [-2.382650, -1.561882, -2.045083, 64433.000000]
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
    fsi2d.test_vals = [4.000000, 0.000000, -3.743227, -4.133479]
    fsi2d.multizone= True
    fsi2d.unsteady = True
    fsi2d.enabled_with_tsan = False
    test_list.append(fsi2d)

    # FSI, Static, 2D, new mesh solver
    stat_fsi           = TestCase('stat_fsi')
    stat_fsi.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi.cfg_file  = "config.cfg"
    stat_fsi.test_iter = 7
    stat_fsi.test_vals         = [-5.397954, -5.719688, 0.000000, 10.000000]
    stat_fsi.test_vals_aarch64 = [-5.423016, -5.753459, 0.000000, 10.000000]
    stat_fsi.multizone = True
    test_list.append(stat_fsi)

    # FSI, Dynamic, 2D, new mesh solver
    dyn_fsi           = TestCase('dyn_fsi')
    dyn_fsi.cfg_dir   = "fea_fsi/dyn_fsi"
    dyn_fsi.cfg_file  = "config.cfg"
    dyn_fsi.test_iter = 4
    dyn_fsi.test_vals         = [-4.355806, -4.060582, 0.000000, 103.000000]
    dyn_fsi.test_vals_aarch64 = [-4.355806, -4.060582, 0.000000, 103.000000]
    dyn_fsi.multizone = True
    dyn_fsi.unsteady  = True
    test_list.append(dyn_fsi)

    # FSI, Static, 2D, new mesh solver, restart
    stat_fsi_restart           = TestCase('stat_fsi_restart')
    stat_fsi_restart.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi_restart.cfg_file  = "config_restart.cfg"
    stat_fsi_restart.test_iter = 1
    stat_fsi_restart.test_vals         = [-3.474078, -4.242240, 0.000000, 36.000000]
    stat_fsi_restart.test_vals_aarch64 = [-3.474081, -4.242372, 0.000000, 37.000000]
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
    mms_fvm_ns.test_vals = [-2.851428, 2.192348, 0.000000, 0.000000]
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
