#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.5.1 "Blackbird"
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

    #########################
    ## FLAMELET solver ###
    #########################

    # 2D planar laminar premixed flame on isothermal burner restar
    fgm_planar_restart = TestCase('fgm_planar')
    fgm_planar_restart.cfg_dir = "flamelet/laminar_premixed_flame"
    fgm_planar_restart.cfg_file = "fgm_planar_restart.cfg"
    fgm_planar_restart.test_iter = 10
    fgm_planar_restart.test_vals = [-15.229677, -15.060158, -15.304078, -8.446904, -15.011948, -15.920890]
    fgm_planar_restart.new_output = True
    test_list.append(fgm_planar_restart)

    #########################
    ## NEMO solver ###
    #########################

    # Adiabatic thermal bath
    thermalbath = TestCase('thermalbath')
    thermalbath.cfg_dir = "nonequilibrium/thermalbath/finitechemistry"
    thermalbath.cfg_file = "thermalbath.cfg"
    thermalbath.test_iter = 10
    thermalbath.test_vals = [0.945997, 0.945997, -12.039262, -12.171767, -32.000000, 10.013239]
    thermalbath.new_output = True
    test_list.append(thermalbath)

    # Adiabatic thermal bath
    ionized = TestCase('ionized')
    ionized.cfg_dir = "nonequilibrium/thermalbath/finitechemistry"
    ionized.cfg_file = "weakly_ionized.cfg"
    ionized.test_iter = 10
    ionized.test_vals = [-29.806157, -11.130797, -11.337264, -17.235059, -17.578729, -15.190274, -25.013626, -32.000000, -5.174887, 0.000000, 0.000000]
    ionized.test_vals_aarch64 = [-29.816386, -10.729986, -11.720016, -17.484469, -18.237891, -15.241605, -24.956918, -32.000000, -5.727244, 0.000000, 0.000000]
    ionized.new_output = True
    test_list.append(ionized)

    # Adiabatic frozen thermal bath
    thermalbath_frozen = TestCase('thermalbath_frozen')
    thermalbath_frozen.cfg_dir = "nonequilibrium/thermalbath/frozen"
    thermalbath_frozen.cfg_file = "thermalbath_frozen.cfg"
    thermalbath_frozen.test_iter = 10
    thermalbath_frozen.test_vals = [-32.000000, -32.000000, -11.962477, -11.962477, -32.000000, 10.013545]
    thermalbath_frozen.new_output = True
    test_list.append(thermalbath_frozen)

    # Inviscid single wedge, ausm, implicit
    invwedge_a = TestCase('invwedge_ausm')
    invwedge_a.cfg_dir = "nonequilibrium/invwedge"
    invwedge_a.cfg_file = "invwedge_ausm.cfg"
    invwedge_a.test_iter = 10
    invwedge_a.test_vals = [-1.042842, -1.567605, -18.301264, -18.628631, -18.574668, 2.275192, 1.879772, 5.319421, 0.873699]
    invwedge_a.test_vals_aarch64 = [-1.042842, -1.567605, -18.301264, -18.628631, -18.574668, 2.275192, 1.879772, 5.319421, 0.873699]
    invwedge_a.new_output = True
    test_list.append(invwedge_a)

    # Inviscid single wedge, ausm+-up2, implicit
    invwedge_ap2 = TestCase('invwedge_ap2')
    invwedge_ap2.cfg_dir = "nonequilibrium/invwedge"
    invwedge_ap2.cfg_file = "invwedge_ausmplusup2.cfg"
    invwedge_ap2.test_iter = 10
    invwedge_ap2.test_vals = [-0.952589, -1.477352, -16.736014, -17.064021, -17.009120, 2.387086, 1.287286, 5.403046, 0.956402]
    invwedge_ap2.test_vals_aarch64 = [-0.952589, -1.477352, -16.736014, -17.064021, -17.009120, 2.387086, 1.287286, 5.403046, 0.956402]
    invwedge_ap2.new_output = True
    test_list.append(invwedge_ap2)

    # Inviscid single wedge, msw, implicit
    invwedge_msw = TestCase('invwedge_msw')
    invwedge_msw.cfg_dir = "nonequilibrium/invwedge"
    invwedge_msw.cfg_file = "invwedge_msw.cfg"
    invwedge_msw.test_iter = 10
    invwedge_msw.test_vals = [-1.165957, -1.690720, -18.298756, -18.626164, -18.572159, 2.151638, 1.721236, 5.193813, 0.751584]
    invwedge_msw.test_vals_aarch64 = [-1.165957, -1.690720, -18.298756, -18.626164, -18.572159, 2.151638, 1.721236, 5.193813, 0.751584]
    invwedge_msw.new_output = True
    test_list.append(invwedge_msw)

    # Inviscid single wedge, roe, implicit
    invwedge_roe = TestCase('invwedge_roe')
    invwedge_roe.cfg_dir = "nonequilibrium/invwedge"
    invwedge_roe.cfg_file = "invwedge_roe.cfg"
    invwedge_roe.test_iter = 10
    invwedge_roe.test_vals = [-1.038582, -1.563344, -18.300307, -18.627706, -18.573706, 2.278987, 1.861307, 5.323753, 0.874900]
    invwedge_roe.test_vals_aarch64 = [-1.038582, -1.563344, -18.300307, -18.627706, -18.573706, 2.278987, 1.861307, 5.323753, 0.874900]
    invwedge_roe.new_output = True
    test_list.append(invwedge_roe)

    # Inviscid single wedge, lax, implicit
    invwedge_lax = TestCase('invwedge_lax')
    invwedge_lax.cfg_dir = "nonequilibrium/invwedge"
    invwedge_lax.cfg_file = "invwedge_lax.cfg"
    invwedge_lax.test_iter = 10
    invwedge_lax.test_vals = [-1.075662, -1.600425, -32.000000, -32.000000, -24.972431, 2.252952, 1.725158, 5.282140, 0.848823]
    invwedge_lax.test_vals_aarch64 = [-1.075662, -1.600425, -32.000000, -32.000000, -24.972431, 2.252952, 1.725158, 5.282140, 0.848823]
    invwedge_lax.new_output = True
    test_list.append(invwedge_lax)

    # Inviscid single wedge, implicit, AUSM+M scheme
    invwedge_ausm_m = TestCase('invwedge_ausm_m')
    invwedge_ausm_m.cfg_dir = "nonequilibrium/invwedge"
    invwedge_ausm_m.cfg_file = "invwedge_am.cfg"
    invwedge_ausm_m.test_iter = 10
    invwedge_ausm_m.test_vals = [-1.055083, -1.579845, -16.739725, -17.063618, -17.012831, 2.265430, 1.797602, 5.302740, 0.85654]
    invwedge_ausm_m.test_vals_aarch64 = [-1.055083, -1.579845, -16.739725, -17.063618, -17.012831, 2.265430, 1.797602, 5.302740, 0.85654]
    invwedge_ausm_m.new_output = True
    test_list.append(invwedge_ausm_m)

    # Inviscid single wedge, implicit, NEMO supersonic inlet
    invwedge_ss_inlet = TestCase('invwedge_ss_inlet')
    invwedge_ss_inlet.cfg_dir = "nonequilibrium/invwedge"
    invwedge_ss_inlet.cfg_file = "invwedge_ss_inlet.cfg"
    invwedge_ss_inlet.test_iter = 10
    invwedge_ss_inlet.test_vals = [-1.042718, -1.567481, -18.250175, -18.579516, -18.523248, 2.275305, 1.880068, 5.319548, 0.873821]
    invwedge_ss_inlet.test_vals_aarch64 = [-1.042718, -1.567481, -18.250175, -18.579516, -18.523248, 2.275305, 1.880068, 5.319548, 0.873821]
    invwedge_ss_inlet.new_output = True
    test_list.append(invwedge_ss_inlet)

    # Viscous single cone - axisymmetric
    visc_cone = TestCase('visc_cone')
    visc_cone.cfg_dir = "nonequilibrium/viscous"
    visc_cone.cfg_file = "axi_visccone.cfg"
    visc_cone.test_iter = 10
    visc_cone.test_vals = [-5.222212, -5.746462, -20.569425, -20.633786, -20.547642, 1.255865, -3.208363, -0.016006, 0.093455, 32633.000000]
    visc_cone.test_vals_aarch64 = [-5.222212, -5.746462, -20.569425, -20.633786, -20.547642, 1.255865, -3.208363, -0.016006, 0.093455, 32633.000000]
    visc_cone.new_output = True
    test_list.append(visc_cone)

    # Viscous single wedge with Mutation++
    #viscwedge_mpp = TestCase('viscwedge_mpp')
    #viscwedge_mpp.cfg_dir = "nonequilibrium/viscwedge_mpp"
    #viscwedge_mpp.cfg_file = "viscwedge_mpp.cfg"
    #viscwedge_mpp.test_iter = 10
    #viscwedge_mpp.test_vals = [-20.608474, -20.586446,-20.707524, -5.171304,-5.696067,-1.548350,-2.071211,2.231054,-2.545494]
    #viscwedge_mpp.new_output = True
    #test_list.append(viscwedge_mpp)

    # Viscous single wedge - super catalytic walls
    super_cat = TestCase('super_cat')
    super_cat.cfg_dir = "nonequilibrium/viscous"
    super_cat.cfg_file = "super_cat.cfg"
    super_cat.test_iter = 10
    super_cat.test_vals = [-5.232590, -5.757884, -20.727046, -20.748136, -20.564044, 1.246889, -3.205235, -0.028406, 0.250857, 3.2459e+04]
    super_cat.su2_exec = "mpirun -n 2 SU2_CFD"
    super_cat.timeout = 1600
    super_cat.new_output = True
    super_cat.tol = 0.00001
    test_list.append(super_cat)

    # Viscous single wedge - partially catalytic walls
    partial_cat = TestCase('partial_cat')
    partial_cat.cfg_dir = "nonequilibrium/viscous"
    partial_cat.cfg_file = "partial_cat.cfg"
    partial_cat.test_iter = 10
    partial_cat.test_vals = [-5.210300, -5.735063, -20.880374, -20.825890, -23.475263, 1.806281, -2.813924, -0.078469, 0.496017, 2.9021e+04]
    partial_cat.su2_exec = "mpirun -n 2 SU2_CFD"
    partial_cat.timeout = 1600
    partial_cat.new_output = True
    partial_cat.tol = 0.00001
    test_list.append(partial_cat)

    ##########################
    ### Compressible Euler ###
    ##########################

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 20
    channel.test_vals = [-2.647975, 2.818090, 0.022280, 0.004644]
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.014140, -3.537888, 0.333403, 0.021227]
    test_list.append(naca0012)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-0.944740, 4.782451, -0.208522, 0.036742]
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-9.277150, -8.694005, 0.281703, 0.011821]
    oneram6.timeout   = 3200
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-7.379831, -1.886302, 0.300000, 0.019471]
    test_list.append(fixedCL_naca0012)

    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals         = [-1.217981, 4.256386, 0.009084, 0.016823]
    polar_naca0012.test_vals_aarch64 = [-2.936433, 2.478784, 0.005113, 0.008684]
    polar_naca0012.command   = TestCase.Command(exec = "compute_polar.py", param = "-i 11")
    test_list.append(polar_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.540009, 6.916653, 0.000000, 1.868976] #last 4 columns
    test_list.append(bluntbody)

    # Equivalent area NACA64-206
    ea_naca64206           = TestCase('ea_naca64206')
    ea_naca64206.cfg_dir   = "optimization_euler/equivalentarea_naca64206"
    ea_naca64206.cfg_file  = "NACA64206.cfg"
    ea_naca64206.test_iter = 10
    ea_naca64206.test_vals = [-1.076215, -0.391987, -0.000701, 67775.0]
    test_list.append(ea_naca64206)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 100
    flatplate.test_vals = [-9.336395, -3.849426, 0.001112, 0.036276, 2.361500, -2.325300, -2.279700, -2.279700]
    test_list.append(flatplate)

    # Custom objective function
    flatplate_udobj           = TestCase('flatplate_udobj')
    flatplate_udobj.cfg_dir   = "user_defined_functions"
    flatplate_udobj.cfg_file  = "lam_flatplate.cfg"
    flatplate_udobj.test_iter = 20
    flatplate_udobj.test_vals = [-6.653802, -1.181430, -0.794887, 0.000611, -0.000369, 0.000736, -0.001104, 596.690000, 299.800000, 296.890000, 21.492000, 0.563990, 37.148, 2.278700]
    test_list.append(flatplate_udobj)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-6.754517, -1.286785, -0.213640, 0.706519, 0.158870]
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.858484, -1.396528, -1.854558, 110.033249, 0.001951]
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-5.050847, 0.648238, 0.000200, 13.639839, -2.047000]
    poiseuille.tol       = 0.001
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals         = [-12.492939, -7.672950, -0.000000, 2.085796]
    poiseuille_profile.test_vals_aarch64 = [-12.492842, -7.672800, -0.000000, 2.085796]
    test_list.append(poiseuille_profile)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.004689, -5.265793, 0.809463, 0.062016, -80577.000000]
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510640, 4.868654, 0.813724, 0.062439, -80115.000000]
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.425725, 4.868653, 0.813724, 0.062438]
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.147548, -6.729213, -0.176227, 0.057731]
    test_list.append(turb_flatplate)

    # Flat plate (compressible) with species inlet
    turb_flatplate_species           = TestCase('turb_flatplate_species')
    turb_flatplate_species.cfg_dir   = "rans/flatplate"
    turb_flatplate_species.cfg_file  = "turb_SA_flatplate_species.cfg"
    turb_flatplate_species.test_iter = 20
    turb_flatplate_species.test_vals = [-4.147548, -0.634735, -1.770801, 1.335176, -3.250308, 9, -6.700992, 5, -6.999234, 10, -6.033847, 0.996033, 0.996033]
    test_list.append(turb_flatplate_species)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.388839, -6.689413, 0.230321, 0.157640, -32539.000000]
    turb_oneram6.timeout   = 3200
    test_list.append(turb_oneram6)

    # ONERA M6 Wing - vorticity confinement
    turb_oneram6_vc = TestCase('turb_oneram6_vc')
    turb_oneram6_vc.cfg_dir = "rans/oneram6"
    turb_oneram6_vc.cfg_file = "turb_ONERAM6_vc.cfg"
    turb_oneram6_vc.test_iter = 15
    turb_oneram6_vc.test_vals = [-2.262387, -6.626467, 0.228393, 0.140799, -2.7107e+04]
    turb_oneram6_vc.timeout = 3200
    test_list.append(turb_oneram6_vc)

    # ONERA M6 Wing - Newton-Krylov
    turb_oneram6_nk           = TestCase('turb_oneram6_nk')
    turb_oneram6_nk.cfg_dir   = "rans/oneram6"
    turb_oneram6_nk.cfg_file  = "turb_ONERAM6_nk.cfg"
    turb_oneram6_nk.test_iter = 20
    turb_oneram6_nk.test_vals = [-4.892253, -4.514006, -11.432312, 0.221026, 0.045570, 2.000000, -0.899460, 31.384000]
    turb_oneram6_nk.timeout   = 600
    turb_oneram6_nk.tol       = 0.0001
    test_list.append(turb_oneram6_nk)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 10
    turb_naca0012_sa.test_vals = [-8.621456, -10.378269, 1.064502, 0.019710, 20.000000, -1.811700, 20.000000, -5.171326, -46.506000]
    turb_naca0012_sa.timeout   = 3200
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-11.422619, -12.803419, -5.867375, 1.049989, 0.019163, -1.827695, -38.695000]
    turb_naca0012_sst.timeout   = 3200
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals = [-11.366100, -12.643576, -5.749377, 1.005234, 0.019017, -1.818746]
    turb_naca0012_sst_sust.timeout   = 3200
    test_list.append(turb_naca0012_sst_sust)

    # NACA0012 (SST, 2003m, Vorticity)
    turb_naca0012_sst_2003_Vm           = TestCase('turb_naca0012_sst_2003_Vm')
    turb_naca0012_sst_2003_Vm.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_2003_Vm.cfg_file  = "turb_NACA0012_sst_2003-Vm.cfg"
    turb_naca0012_sst_2003_Vm.test_iter = 10
    turb_naca0012_sst_2003_Vm.test_vals = [-7.672926, -10.025010, -3.365892, 1.048735, 0.019723, -2.052543]
    turb_naca0012_sst_2003_Vm.timeout   = 3200
    test_list.append(turb_naca0012_sst_2003_Vm)

    # NACA0012 (SST, 1994m Kato-Launder)
    turb_naca0012_sst_1994_KLm           = TestCase('turb_naca0012_sst_1994_KLm')
    turb_naca0012_sst_1994_KLm.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_1994_KLm.cfg_file  = "turb_NACA0012_sst_1994-KLm.cfg"
    turb_naca0012_sst_1994_KLm.test_iter = 10
    turb_naca0012_sst_1994_KLm.test_vals = [-8.567222, -10.798741, -3.990574, 1.049274, 0.019199, -1.809143]
    turb_naca0012_sst_1994_KLm.timeout   = 3200
    test_list.append(turb_naca0012_sst_1994_KLm)


    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [-5.216685, -9.561914, -1.565777, 1.022393, 0.040542, -3.729635]
    turb_naca0012_sst_fixedvalues.timeout   = 3200
    test_list.append(turb_naca0012_sst_fixedvalues)

    # NACA0012 (SST, explicit Euler for flow and turbulence equations)
    turb_naca0012_sst_expliciteuler           = TestCase('turb_naca0012_sst_expliciteuler')
    turb_naca0012_sst_expliciteuler.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_expliciteuler.cfg_file  = "turb_NACA0012_sst_expliciteuler.cfg"
    turb_naca0012_sst_expliciteuler.test_iter = 10
    turb_naca0012_sst_expliciteuler.test_vals = [-3.532228, -3.157766, 3.364025, 1.124824, 0.501717, -float("inf")]
    turb_naca0012_sst_expliciteuler.timeout   = 3200
    test_list.append(turb_naca0012_sst_expliciteuler)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389575, -8.409529, 0.000048, 0.056329] #last 4 columns
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
    axi_rans_air_nozzle_restart.test_vals = [-12.089268, -7.493381, -8.716391, -4.021218, -1924.800000]
    axi_rans_air_nozzle_restart.tol       = 0.0001
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
    turb_naca0012_sst_restart_mg.test_vals = [-7.619889, -7.729499, -1.981039, -0.000016, 0.079062]
    turb_naca0012_sst_restart_mg.timeout   = 3200
    turb_naca0012_sst_restart_mg.tol       = 0.000001
    test_list.append(turb_naca0012_sst_restart_mg)

    #############################
    ### Incompressible Euler  ###
    #############################

    # NACA0012 Hydrofoil
    inc_euler_naca0012           = TestCase('inc_euler_naca0012')
    inc_euler_naca0012.cfg_dir   = "incomp_euler/naca0012"
    inc_euler_naca0012.cfg_file  = "incomp_NACA0012.cfg"
    inc_euler_naca0012.test_iter = 20
    inc_euler_naca0012.test_vals = [-4.801273, -3.773079, 0.495236, 0.007346]
    test_list.append(inc_euler_naca0012)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals = [-5.982321, -4.953536, 0.000454, 0.121390]
    test_list.append(inc_nozzle)

    #############################
    ### Incompressible N-S    ###
    #############################

    # Laminar cylinder
    inc_lam_cylinder          = TestCase('inc_lam_cylinder')
    inc_lam_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_lam_cylinder.cfg_file  = "incomp_cylinder.cfg"
    inc_lam_cylinder.test_iter = 10
    inc_lam_cylinder.test_vals = [-4.004072, -3.194881, -0.076553, 7.780048]
    test_list.append(inc_lam_cylinder)

    # Buoyancy-driven cavity
    inc_buoyancy          = TestCase('inc_buoyancy')
    inc_buoyancy.cfg_dir   = "incomp_navierstokes/buoyancy_cavity"
    inc_buoyancy.cfg_file  = "lam_buoyancy_cavity.cfg"
    inc_buoyancy.test_iter = 20
    inc_buoyancy.test_vals = [-4.435827, 0.508037, 0.000000, 0.000000]
    test_list.append(inc_buoyancy)

    # Laminar heated cylinder with polynomial fluid model
    inc_poly_cylinder          = TestCase('inc_poly_cylinder')
    inc_poly_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_poly_cylinder.cfg_file  = "poly_cylinder.cfg"
    inc_poly_cylinder.test_iter = 20
    inc_poly_cylinder.test_vals = [-7.791831, -2.062292, 0.013040, 1.913997, -171.120000]
    test_list.append(inc_poly_cylinder)

    # X-coarse laminar bend as a mixed element CGNS test
    inc_lam_bend          = TestCase('inc_lam_bend')
    inc_lam_bend.cfg_dir   = "incomp_navierstokes/bend"
    inc_lam_bend.cfg_file  = "lam_bend.cfg"
    inc_lam_bend.test_iter = 10
    inc_lam_bend.test_vals = [-3.447746, -3.085237, -0.020816, 1.147373]
    test_list.append(inc_lam_bend)

    # 3D laminar channnel with 1 cell in flow direction, streamwise periodic
    sp_pipeSlice_3d_dp_hf_tp           = TestCase('sp_pipeSlice_3d_dp_hf_tp')
    sp_pipeSlice_3d_dp_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/pipeSlice_3d"
    sp_pipeSlice_3d_dp_hf_tp.cfg_file  = "sp_pipeSlice_3d_dp_hf_tp.cfg"
    sp_pipeSlice_3d_dp_hf_tp.test_iter = 10
    sp_pipeSlice_3d_dp_hf_tp.test_vals = [-11.119796, -11.234737, -8.694310, -0.000023] #last 4 lines
    test_list.append(sp_pipeSlice_3d_dp_hf_tp)

    # 2D pin array with heat transfer BC on pin surfaces
    inc_heatTransfer_BC           = TestCase('inc_heatTransfer_BC')
    inc_heatTransfer_BC.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    inc_heatTransfer_BC.cfg_file  = "BC_HeatTransfer.cfg"
    inc_heatTransfer_BC.test_iter = 50
    inc_heatTransfer_BC.test_vals = [-8.242651, -7.341179, -7.407346, -0.152603, -1667.300000] #last 5 lines
    test_list.append(inc_heatTransfer_BC)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.788595, -11.040557, -0.000002, 0.309519]
    test_list.append(inc_turb_naca0012)

    # NACA0012, SST_SUST
    inc_turb_naca0012_sst_sust           = TestCase('inc_turb_naca0012_sst_sust')
    inc_turb_naca0012_sst_sust.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_sst_sust.cfg_file  = "naca0012_SST_SUST.cfg"
    inc_turb_naca0012_sst_sust.test_iter = 20
    inc_turb_naca0012_sst_sust.test_vals = [-7.274050, 0.145887, -0.000001, 0.312023]
    test_list.append(inc_turb_naca0012_sst_sust)

    ####################
    ### DG-FEM Euler ###
    ####################

    # NACA0012
    fem_euler_naca0012           = TestCase('fem_euler_naca0012')
    fem_euler_naca0012.cfg_dir   = "hom_euler/NACA0012_5thOrder"
    fem_euler_naca0012.cfg_file  = "fem_NACA0012_reg.cfg"
    fem_euler_naca0012.test_iter = 10
    fem_euler_naca0012.test_vals = [-6.519946,-5.976944,0.255551,0.000028] #last 4 columns
    test_list.append(fem_euler_naca0012)

    ############################
    ### DG-FEM Navier-Stokes ###
    ############################

    # Flat plate
    fem_ns_flatplate           = TestCase('fem_ns_flatplate')
    fem_ns_flatplate.cfg_dir   = "hom_navierstokes/FlatPlate/nPoly4"
    fem_ns_flatplate.cfg_file  = "lam_flatplate_reg.cfg"
    fem_ns_flatplate.test_iter = 25
    fem_ns_flatplate.test_vals = [1.383727,3.175247,0.058387,0.257951] #last 4 columns
    test_list.append(fem_ns_flatplate)

    # Steady cylinder
    fem_ns_cylinder           = TestCase('fem_ns_cylinder')
    fem_ns_cylinder.cfg_dir   = "hom_navierstokes/CylinderViscous/nPoly3"
    fem_ns_cylinder.cfg_file  = "fem_Cylinder_reg.cfg"
    fem_ns_cylinder.test_iter = 10
    fem_ns_cylinder.test_vals = [0.454960,0.979123,-0.000028,79.984799] #last 4 columns
    test_list.append(fem_ns_cylinder)

    # Steady sphere
    fem_ns_sphere           = TestCase('fem_ns_sphere')
    fem_ns_sphere.cfg_dir   = "hom_navierstokes/SphereViscous/nPoly3_QuadDominant"
    fem_ns_sphere.cfg_file  = "fem_Sphere_reg.cfg"
    fem_ns_sphere.test_iter = 10
    fem_ns_sphere.test_vals = [-0.288121,0.240324,0.000258,21.797363] #last 4 columns
    fem_ns_sphere.command   = TestCase.Command(exec = "SU2_CFD")
    test_list.append(fem_ns_sphere)

    # Unsteady sphere ADER
    fem_ns_sphere_ader           = TestCase('fem_ns_sphere_ader')
    fem_ns_sphere_ader.cfg_dir   = "hom_navierstokes/SphereViscous/nPoly3_QuadDominant"
    fem_ns_sphere_ader.cfg_file  = "fem_Sphere_reg_ADER.cfg"
    fem_ns_sphere_ader.test_iter = 10
    fem_ns_sphere_ader.test_vals = [-35.000000,-35.000000,0.000047,31.110911] #last 4 columns
    fem_ns_sphere_ader.command   = TestCase.Command(exec = "SU2_CFD")
    test_list.append(fem_ns_sphere_ader)

    # Unsteady cylinder
    fem_ns_unsteady_cylinder           = TestCase('fem_ns_unsteady_cylinder')
    fem_ns_unsteady_cylinder.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder.cfg_file  = "fem_unst_cylinder.cfg"
    fem_ns_unsteady_cylinder.test_iter = 11
    fem_ns_unsteady_cylinder.test_vals = [-3.558582,-3.014464,-0.038927,1.383983] #last 4 columns
    fem_ns_unsteady_cylinder.command   = TestCase.Command(exec = "SU2_CFD")
    fem_ns_unsteady_cylinder.unsteady  = True
    test_list.append(fem_ns_unsteady_cylinder)

    # Unsteady cylinder ADER
    fem_ns_unsteady_cylinder_ader           = TestCase('fem_ns_unsteady_cylinder_ader')
    fem_ns_unsteady_cylinder_ader.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder_ader.cfg_file  = "fem_unst_cylinder_ADER.cfg"
    fem_ns_unsteady_cylinder_ader.test_iter = 11
    fem_ns_unsteady_cylinder_ader.test_vals = [-35.000000,-35.000000,-0.041003,1.391339] #last 4 columns
    fem_ns_unsteady_cylinder_ader.command   = TestCase.Command(exec = "SU2_CFD")
    fem_ns_unsteady_cylinder_ader.unsteady  = True
    test_list.append(fem_ns_unsteady_cylinder_ader)

    ###########################
    ### Turbulence modeling ###
    ###########################

    # SA Baseline (Identical to RANS SA RAE2822)
    turbmod_sa_bsl_rae2822           = TestCase('turbmod_sa_bsl_rae2822')
    turbmod_sa_bsl_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_bsl_rae2822.cfg_file  = "turb_SA_BSL_RAE2822.cfg"
    turbmod_sa_bsl_rae2822.test_iter = 20
    turbmod_sa_bsl_rae2822.test_vals = [-2.004689, 0.742306, 0.497308, -5.265793, 0.809463, 0.062016]
    turbmod_sa_bsl_rae2822.new_output = True
    test_list.append(turbmod_sa_bsl_rae2822)

    # SA Negative
    turbmod_sa_neg_rae2822           = TestCase('turbmod_sa_neg_rae2822')
    turbmod_sa_neg_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_neg_rae2822.cfg_file  = "turb_SA_NEG_RAE2822.cfg"
    turbmod_sa_neg_rae2822.test_iter = 10
    turbmod_sa_neg_rae2822.test_vals         = [-1.374695, 1.976506, 1.898195, 4.831133, 1.187310, 0.426019, -86764]
    turbmod_sa_neg_rae2822.test_vals_aarch64 = [-1.347530, 1.439078, 1.306846, -1.928774, 1.480543, 0.571601, -91503]
    turbmod_sa_neg_rae2822.new_output = True
    test_list.append(turbmod_sa_neg_rae2822)

    # SA Compressibility Correction
    turbmod_sa_comp_rae2822           = TestCase('turbmod_sa_comp_rae2822')
    turbmod_sa_comp_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_comp_rae2822.cfg_file  = "turb_SA_COMP_RAE2822.cfg"
    turbmod_sa_comp_rae2822.test_iter = 20
    turbmod_sa_comp_rae2822.test_vals = [-2.004687, 0.742304, 0.497309, -5.266081, 0.809467, 0.062029]
    turbmod_sa_comp_rae2822.new_output = True
    test_list.append(turbmod_sa_comp_rae2822)

    # SA Edwards
    turbmod_sa_edw_rae2822           = TestCase('turbmod_sa_edw_rae2822')
    turbmod_sa_edw_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_edw_rae2822.cfg_file  = "turb_SA_EDW_RAE2822.cfg"
    turbmod_sa_edw_rae2822.test_iter = 20
    turbmod_sa_edw_rae2822.test_vals = [-2.004687, 0.742306, 0.497310, -5.290769, 0.809485, 0.062036]
    turbmod_sa_edw_rae2822.new_output = True
    test_list.append(turbmod_sa_edw_rae2822)

    # SA Compressibility and Edwards
    turbmod_sa_comp_edw_rae2822           = TestCase('turbmod_sa_comp_edw_rae2822')
    turbmod_sa_comp_edw_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_comp_edw_rae2822.cfg_file  = "turb_SA_COMP_EDW_RAE2822.cfg"
    turbmod_sa_comp_edw_rae2822.test_iter = 20
    turbmod_sa_comp_edw_rae2822.test_vals = [-2.004685, 0.742307, 0.497311, -5.290750, 0.809487, 0.062045]
    turbmod_sa_comp_edw_rae2822.new_output = True
    test_list.append(turbmod_sa_comp_edw_rae2822)

    # SA QCR
    turbmod_sa_qcr_rae2822           = TestCase('turbmod_sa_qcr_rae2822')
    turbmod_sa_qcr_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_qcr_rae2822.cfg_file  = "turb_SA_QCR_RAE2822.cfg"
    turbmod_sa_qcr_rae2822.test_iter = 20
    turbmod_sa_qcr_rae2822.test_vals = [-2.004793, 0.742353, 0.497315, -5.265974, 0.807841, 0.062027]
    turbmod_sa_qcr_rae2822.new_output = True
    test_list.append(turbmod_sa_qcr_rae2822)

    ############################
    ###      Transition      ###
    ############################

    # Schubauer-Klebanoff Natural Transition Case
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.test_vals    = [-7.994740, -13.240225, 0.000046, 0.007987]
    test_list.append(schubauer_klebanoff_transition)

    #####################################
    ### Cont. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    contadj_naca0012           = TestCase('contadj_naca0012')
    contadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012.test_iter = 5
    contadj_naca0012.test_vals = [-9.300815, -14.587362, 0.300920, 0.019552]
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.130993, -12.702085, 0.685900, 0.007594]
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10
    contadj_wedge.test_vals = [2.872691, -2.755572, 853010.000000, 0.000000] #last 4 columns
    test_list.append(contadj_wedge)

    # Inviscid fixed CL NACA0012
    contadj_fixed_CL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixed_CL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixed_CL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixed_CL_naca0012.test_iter = 100
    contadj_fixed_CL_naca0012.test_vals = [0.275856, -5.200511, 0.342710, 0.000105]
    test_list.append(contadj_fixed_CL_naca0012)

    ###################################
    ### Cont. adj. compressible N-S ###
    ###################################

    # Adjoint laminar cylinder
    contadj_ns_cylinder           = TestCase('contadj_ns_cylinder')
    contadj_ns_cylinder.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder.test_iter = 20
    contadj_ns_cylinder.test_vals = [-3.651430, -9.113079, 2.056700, -0.000000]
    test_list.append(contadj_ns_cylinder)

    # Adjoint laminar naca0012 subsonic
    contadj_ns_naca0012_sub           = TestCase('contadj_ns_naca0012_sub')
    contadj_ns_naca0012_sub.cfg_dir   = "cont_adj_navierstokes/naca0012_sub"
    contadj_ns_naca0012_sub.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_sub.test_iter = 20
    contadj_ns_naca0012_sub.test_vals = [-2.743268, -8.215193, 0.518810, 0.001210] #last 4 columns
    test_list.append(contadj_ns_naca0012_sub)

    # Adjoint laminar naca0012 transonic
    contadj_ns_naca0012_trans           = TestCase('contadj_ns_naca0012_trans')
    contadj_ns_naca0012_trans.cfg_dir   = "cont_adj_navierstokes/naca0012_trans"
    contadj_ns_naca0012_trans.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_trans.test_iter = 20
    contadj_ns_naca0012_trans.test_vals = [ -1.039664, -6.575019, 1.772300, 0.012495] #last 4 columns
    test_list.append(contadj_ns_naca0012_trans)

    #######################################################
    ### Cont. adj. compressible RANS (frozen viscosity) ###
    #######################################################

    # Adjoint turbulent NACA0012
    contadj_rans_naca0012           = TestCase('contadj_rans_naca0012')
    contadj_rans_naca0012.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012.cfg_file  = "turb_nasa.cfg"
    contadj_rans_naca0012.test_iter = 20
    contadj_rans_naca0012.test_vals = [-0.794162, -5.761722, 19.214000, -0.000000] #last 4 columns
    test_list.append(contadj_rans_naca0012)

    # Adjoint turbulent NACA0012 with binary restarts
    contadj_rans_naca0012_bin           = TestCase('contadj_rans_naca0012_bin')
    contadj_rans_naca0012_bin.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012_bin.cfg_file  = "turb_nasa_binary.cfg"
    contadj_rans_naca0012_bin.test_iter = 18
    contadj_rans_naca0012_bin.test_vals = [-0.794169, -5.761671, 19.214000, -0.000000] #last 4 columns
    test_list.append(contadj_rans_naca0012_bin)

    # Adjoint turbulent RAE2822
    contadj_rans_rae2822           = TestCase('contadj_rans_rae822')
    contadj_rans_rae2822.cfg_dir   = "cont_adj_rans/rae2822"
    contadj_rans_rae2822.cfg_file  = "turb_SA_RAE2822.cfg"
    contadj_rans_rae2822.test_iter = 20
    contadj_rans_rae2822.test_vals = [-5.372407, -10.874841, -0.212470, 0.005448]
    test_list.append(contadj_rans_rae2822)

    #############################
    ### Compressibele RANS UQ ###
    #############################

    # NACA0012 1c
    turb_naca0012_1c           = TestCase('turb_naca0012_1c')
    turb_naca0012_1c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_1c.cfg_file  = "turb_NACA0012_uq_1c.cfg"
    turb_naca0012_1c.test_iter = 10
    turb_naca0012_1c.test_vals         = [-4.976246, 1.141479, 0.459999, -0.078853]
    turb_naca0012_1c.test_vals_aarch64 = [-4.983548, 1.138789, 0.456805, -0.079667]
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals         = [-5.483323, 0.968820, 0.304757, -0.113468]
    turb_naca0012_2c.test_vals_aarch64 = [-5.483345, 0.968812, 0.305213, -0.113317]
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals = [-5.584310, 0.931348, 0.279056, -0.113209]
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals         = [-5.128867, 1.077141, 0.586532, -0.047632]
    turb_naca0012_p1c1.test_vals_aarch64 = [-5.130279, 1.076643, 0.587076, -0.047445]
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals         = [-5.554659, 0.943705, 0.399234, -0.095799]
    turb_naca0012_p1c2.test_vals_aarch64 = [-5.554645, 0.943709, 0.398620, -0.096021]
    test_list.append(turb_naca0012_p1c2)

    ######################################
    ### Harmonic Balance               ###
    ######################################

    # Description of the regression test
    harmonic_balance           = TestCase('harmonic_balance')
    harmonic_balance.cfg_dir   = "harmonic_balance"
    harmonic_balance.cfg_file  = "HB.cfg"
    harmonic_balance.test_iter = 25
    harmonic_balance.test_vals = [-1.589739, 3.922579, 0.006702, 0.099632]
    harmonic_balance.new_output = False
    test_list.append(harmonic_balance)

    # Turbulent pitching NACA 64a010 airfoil
    hb_rans_preconditioning           = TestCase('hb_rans_preconditioning')
    hb_rans_preconditioning.cfg_dir   = "harmonic_balance/hb_rans_preconditioning"
    hb_rans_preconditioning.cfg_file  = "davis.cfg"
    hb_rans_preconditioning.test_iter = 25
    hb_rans_preconditioning.test_vals = [-1.902098, -5.949275, 0.007768, 0.128061]
    hb_rans_preconditioning.new_output = False
    test_list.append(hb_rans_preconditioning)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Rotating NACA 0012
    rot_naca0012           = TestCase('rot_naca0012')
    rot_naca0012.cfg_dir   = "rotating/naca0012"
    rot_naca0012.cfg_file  = "rot_NACA0012.cfg"
    rot_naca0012.test_iter = 25
    rot_naca0012.test_vals = [-2.698005, 2.845328, -0.079439, 0.002128]
    test_list.append(rot_naca0012)

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.611007, -0.146826, 1.113206, 1.491678]
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.802803, -2.362844, 1.687705, 1.519676]
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-1.158117, 0.067945, 1.399789, 2.220404, 1.399743, 2.218605, -0.453170]
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977520, 3.481804, -0.012377, -0.007389]
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Aeroelastic
    aeroelastic           = TestCase('aeroelastic')
    aeroelastic.cfg_dir   = "aeroelastic"
    aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    aeroelastic.test_iter = 2
    aeroelastic.test_vals = [0.076550, 0.033042, -0.001650, -0.000127]
    aeroelastic.unsteady  = True
    test_list.append(aeroelastic)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714758, -5.883004, -0.215005, 0.023783, -618.160000]
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.004011, -6.876230, 1.487888, 0.421869]
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # Flat plate
    flatplate_unsteady           = TestCase('flatplate_unsteady')
    flatplate_unsteady.cfg_dir   = "navierstokes/flatplate"
    flatplate_unsteady.cfg_file  = "lam_flatplate_unst.cfg"
    flatplate_unsteady.test_iter = 3
    flatplate_unsteady.test_vals = [7.9509e-06, -8.868859, -8.231652, -6.283262, -5.466675, -3.391163, 0.002078, -0.343642]
    flatplate_unsteady.unsteady  = True
    test_list.append(flatplate_unsteady)

    ######################################
    ### NICFD                          ###
    ######################################

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 100
    edge_VW.test_vals = [-5.048044, 1.115667, -0.000009, 0.000000]
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 100
    edge_PPR.test_vals = [-5.400790, 0.739723, -0.000035, 0.000000]
    test_list.append(edge_PPR)

    # Rarefaction Q1D nozzle, include CoolProp fluid model
    coolprop_fluidModel           = TestCase('coolprop_fluidModel')
    coolprop_fluidModel.cfg_dir   = "nicf/coolprop"
    coolprop_fluidModel.cfg_file  = "fluidModel.cfg"
    coolprop_fluidModel.test_iter = 5
    coolprop_fluidModel.test_vals = [-4.525458, -1.578441, 3.439057, 0.000000, 0.000000]
    test_list.append(coolprop_fluidModel)

    # Rarefaction Q1D nozzle, include CoolProp fluid model
    datadriven_fluidModel           = TestCase('datadriven_fluidModel')
    datadriven_fluidModel.cfg_dir   = "nicf/datadriven"
    datadriven_fluidModel.cfg_file  = "datadriven_nozzle.cfg"
    datadriven_fluidModel.test_iter = 20
    datadriven_fluidModel.test_vals = [-1.873575, 0.964932, 5.467181, 0.000000, 0.000000]
    test_list.append(datadriven_fluidModel)

    # Rarefaction Q1D nozzle, include CoolProp transport model
    coolprop_transportModel           = TestCase('coolprop_transportModel')
    coolprop_transportModel.cfg_dir   = "nicf/coolprop"
    coolprop_transportModel.cfg_file  = "transportModel.cfg"
    coolprop_transportModel.test_iter = 5
    coolprop_transportModel.test_vals = [-4.527922, -1.308741, 4.630469, 0.000000, 0.000000]
    test_list.append(coolprop_transportModel)

    ######################################
    ### Turbomachinery                 ###
    ######################################

    # Jones APU Turbocharger restart
    Jones_tc_restart           = TestCase('jones_turbocharger_restart')
    Jones_tc_restart.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_restart.cfg_file  = "Jones_restart.cfg"
    Jones_tc_restart.test_iter = 5
    Jones_tc_restart.test_vals = [-10.691504, -7.643703, 85.827890, 2.277151]
    Jones_tc_restart.new_output = False
    test_list.append(Jones_tc_restart)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals = [-1.937007, 5.338943, 73.357200, 0.915725]
    axial_stage2D.new_output = False
    test_list.append(axial_stage2D)

    # 2D transonic stator restart
    transonic_stator_restart           = TestCase('transonic_stator_restart')
    transonic_stator_restart.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_restart.cfg_file  = "transonic_stator_restart.cfg"
    transonic_stator_restart.test_iter = 20
    transonic_stator_restart.test_vals = [-6.801911, -0.746307, 5.003453, 0.002946]
    transonic_stator_restart.new_output = False
    test_list.append(transonic_stator_restart)

    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 5
    uniform_flow.test_vals = [5.000000, 0.000000, -0.188747, -10.631538]
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow)

    # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 2
    channel_2D.test_vals         = [2.000000, 0.000000, 0.398011, 0.352778, 0.405461]
    channel_2D.test_vals_aarch64 = [2.000000, 0.000000, 0.398036, 0.352783, 0.405462]
    channel_2D.timeout   = 100
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 2
    channel_3D.test_vals         = [2.000000, 0.000000, 0.620176, 0.505161, 0.415248]
    channel_3D.test_vals_aarch64 = [2.000000, 0.000000, 0.620182, 0.505302, 0.415257]
    channel_3D.unsteady  = True
    channel_3D.multizone = True
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [0.150024, 0.491953, 0.677755, 0.963980, 1.006936]
    pipe.unsteady  = True
    pipe.multizone = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [3.000000, 0.000000, 0.777572, 1.134804, 1.224137]
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone  = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 1.214350, 1.663914]
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone  = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [13.000000, -0.604409, -1.523885]
    bars_SST_2D.command   = TestCase.Command(exec = "SU2_CFD")
    bars_SST_2D.multizone = True
    test_list.append(bars_SST_2D)

    # Sliding mesh with incompressible flows (steady)
    slinc_steady           = TestCase('slinc_steady')
    slinc_steady.cfg_dir   = "sliding_interface/incompressible_steady"
    slinc_steady.cfg_file  = "config.cfg"
    slinc_steady.test_iter = 19
    slinc_steady.test_vals = [19.000000, -1.803732, -2.108492] #last 3 columns
    slinc_steady.command   = TestCase.Command(exec = "SU2_CFD")
    slinc_steady.timeout   = 100
    slinc_steady.tol       = 0.00002
    slinc_steady.multizone = True
    test_list.append(slinc_steady)

    # Sliding mesh with incompressible flows (unsteady)
    # slinc_unsteady           = TestCase('slinc_unsteady')
    # slinc_unsteady.cfg_dir   = "sliding_interface/incompressible_unsteady"
    # slinc_unsteady.cfg_file  = "config.cfg"
    # slinc_unsteady.test_iter = 19
    # slinc_unsteady.test_vals = [-3.513701,1.931626,0.000000,0.000000] #last 4 columns
    # slinc_unsteady.command   = TestCase.Command(exec = "SU2_CFD")
    # slinc_unsteady.timeout   = 100
    # slinc_unsteady.unsteady  = True
    # test_list.append(slinc_unsteady)

    ##########################
    ### FEA - FSI          ###
    ##########################

    # Static beam, 3d
    statbeam3d           = TestCase('statbeam3d')
    statbeam3d.cfg_dir   = "fea_fsi/StatBeam_3d"
    statbeam3d.cfg_file  = "configBeam_3d.cfg"
    statbeam3d.test_iter = 0
    statbeam3d.test_vals         = [-8.396797, -8.162206, -8.156102, 64095.0] #last 4 columns
    statbeam3d.test_vals_aarch64 = [-8.396793, -8.162255, -8.156118, 64095.0] #last 4 columns
    statbeam3d.command   = TestCase.Command(exec = "parallel_computation_fsi.py", param = "-f")
    test_list.append(statbeam3d)

    # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.test_iter = 6
    dynbeam2d.test_vals = [-3.240015, 2.895057, -0.353146, 6.6127e+04] #last 4 columns
    dynbeam2d.unsteady  = True
    test_list.append(dynbeam2d)

    # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [4, 0, -3.743210, -4.133483] #last 4 columns
    fsi2d.command   = TestCase.Command(exec = "parallel_computation_fsi.py", param = "-f")
    fsi2d.multizone= True
    fsi2d.unsteady = True
    test_list.append(fsi2d)

    # FSI, Static, 2D, new mesh solver
    stat_fsi           = TestCase('stat_fsi')
    stat_fsi.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi.cfg_file  = "config.cfg"
    stat_fsi.test_iter = 7
    stat_fsi.test_vals = [-3.296605, -4.934646, 0.000000, 7.000000] #last 4 columns
    stat_fsi.multizone = True
    test_list.append(stat_fsi)

    # FSI, Dynamic, 2D, new mesh solver
    dyn_fsi           = TestCase('dyn_fsi')
    dyn_fsi.cfg_dir   = "fea_fsi/dyn_fsi"
    dyn_fsi.cfg_file  = "config.cfg"
    dyn_fsi.test_iter = 4
    dyn_fsi.test_vals = [-4.355829, -4.060587, 5.3837e-08, 98]
    dyn_fsi.multizone = True
    dyn_fsi.unsteady  = True
    test_list.append(dyn_fsi)

    # FSI, Static, 2D, new mesh solver, restart
    stat_fsi_restart           = TestCase('stat_fsi_restart')
    stat_fsi_restart.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi_restart.cfg_file  = "config_restart.cfg"
    stat_fsi_restart.test_iter = 1
    stat_fsi_restart.test_vals = [-3.435926, -4.264912, 0.000000, 28.000000] #last 4 columns
    stat_fsi_restart.multizone = True
    test_list.append(stat_fsi_restart)

    # ###############################
    # ### Radiative Heat Transfer ###
    # ###############################

    # Radiative heat transfer
    p1rad           = TestCase('p1rad')
    p1rad.cfg_dir   = "radiation/p1model"
    p1rad.cfg_file  = "configp1.cfg"
    p1rad.new_output= True
    p1rad.test_iter = 100
    p1rad.test_vals = [-7.743666, -7.921411, -2.111848, 0.098302, -45.023000]
    test_list.append(p1rad)


    # #############################
    # ### Solid Heat Conduction ###
    # #############################

    # 2D pins, periodically connected
    solid_periodic_pins           = TestCase('solid_periodic_pins')
    solid_periodic_pins.cfg_dir   = "solid_heat_conduction/periodic_pins"
    solid_periodic_pins.cfg_file  = "configSolid.cfg"
    solid_periodic_pins.test_iter = 750
    solid_periodic_pins.test_vals = [-15.878977, -14.569206, 300.900000, 425.320000, 0.000000, 5.000000, -1.672737] #last 7 lines
    solid_periodic_pins.test_vals_aarch64 = [-15.879010, -14.569206, 300.900000, 425.320000, 0.000000, 5.000000, -1.672630] #last 7 lines
    test_list.append(solid_periodic_pins)

    # ###############################
    # ### Conjugate heat transfer ###
    # ###############################

    # CHT incompressible
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "coupled_cht/incomp_2d"
    cht_incompressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [-2.128826, -0.588813, -0.588813, -0.588813] #last 4 columns
    cht_incompressible.command   = TestCase.Command(exec = "SU2_CFD")
    cht_incompressible.multizone = True
    test_list.append(cht_incompressible)

    # CHT compressible
    cht_compressible           = TestCase('cht_compressible')
    cht_compressible.cfg_dir   = "coupled_cht/comp_2d"
    cht_compressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_compressible.test_iter = 10
    cht_compressible.test_vals = [-4.256032, -0.532728, -0.532729, -0.532728]
    cht_compressible.command   = TestCase.Command(exec = "SU2_CFD")
    cht_compressible.multizone = True
    test_list.append(cht_compressible)

    # 2D CHT case streamwise periodicity. Also test Multizone PerSurface screen output.
    sp_pinArray_cht_2d_dp_hf           = TestCase('sp_pinArray_cht_2d_dp_hf')
    sp_pinArray_cht_2d_dp_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    sp_pinArray_cht_2d_dp_hf.cfg_file  = "configMaster.cfg"
    sp_pinArray_cht_2d_dp_hf.test_iter = 100
    sp_pinArray_cht_2d_dp_hf.test_vals = [0.255546, -0.778969, -0.964351, -0.752492, 208.023676, 351.010000, -0.000000, -0.752490, 0.752490]
    sp_pinArray_cht_2d_dp_hf.multizone = True
    test_list.append(sp_pinArray_cht_2d_dp_hf)

    # simple small 3D pin case massflow periodic with heatflux BC
    sp_pinArray_3d_cht_mf_hf_tp           = TestCase('sp_pinArray_3d_cht_mf_hf_tp')
    sp_pinArray_3d_cht_mf_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_3d"
    sp_pinArray_3d_cht_mf_hf_tp.cfg_file  = "configMaster.cfg"
    sp_pinArray_3d_cht_mf_hf_tp.test_iter = 30
    sp_pinArray_3d_cht_mf_hf_tp.test_vals         = [-3.773085, -4.220555, -4.811282, -0.009675, 99.879858, 419.200000, 0.000000]
    sp_pinArray_3d_cht_mf_hf_tp.test_vals_aarch64 = [-13.400623, -7.476945, -7.025285, -0.009675, 99.879812, 419.200000, 0.0]
    sp_pinArray_3d_cht_mf_hf_tp.multizone = True
    test_list.append(sp_pinArray_3d_cht_mf_hf_tp)

    ##########################
    ###   Python wrapper   ###
    ##########################

    # NACA0012
    pywrapper_naca0012           = TestCase('pywrapper_naca0012')
    pywrapper_naca0012.cfg_dir   = "euler/naca0012"
    pywrapper_naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    pywrapper_naca0012.test_iter = 100
    pywrapper_naca0012.test_vals = [-6.747210, -6.149915, 0.333445, 0.021241] #last 4 columns
    pywrapper_naca0012.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--parallel -f")
    test_list.append(pywrapper_naca0012)

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    pywrapper_turb_naca0012_sst           = TestCase('pywrapper_turb_naca0012_sst')
    pywrapper_turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    pywrapper_turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    pywrapper_turb_naca0012_sst.test_iter = 10
    pywrapper_turb_naca0012_sst.test_vals = [-11.422619, -12.803419, -5.867375, 1.049989, 0.019163, -1.827695, -38.695000]
    pywrapper_turb_naca0012_sst.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--parallel -f")
    pywrapper_turb_naca0012_sst.timeout   = 3200
    test_list.append(pywrapper_turb_naca0012_sst)

    # Square cylinder
    pywrapper_square_cylinder           = TestCase('pywrapper_square_cylinder')
    pywrapper_square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    pywrapper_square_cylinder.cfg_file  = "turb_square.cfg"
    pywrapper_square_cylinder.test_iter = 10
    pywrapper_square_cylinder.test_vals = [-1.136553, -0.347305, 1.407915, 2.358881, 1.404192, 2.301559, -0.348120]
    pywrapper_square_cylinder.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--parallel -f")
    pywrapper_square_cylinder.unsteady  = True
    test_list.append(pywrapper_square_cylinder)

    # Aeroelastic
    pywrapper_aeroelastic         = TestCase('pywrapper_aeroelastic')
    pywrapper_aeroelastic.cfg_dir   = "aeroelastic"
    pywrapper_aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    pywrapper_aeroelastic.test_iter = 2
    pywrapper_aeroelastic.test_vals = [0.076550, 0.033042, -0.001650, -0.000127]
    pywrapper_aeroelastic.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--parallel -f")
    pywrapper_aeroelastic.unsteady  = True
    test_list.append(pywrapper_aeroelastic)

    # FSI, 2d
    pywrapper_fsi2d           = TestCase('pywrapper_fsi2d')
    pywrapper_fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    pywrapper_fsi2d.cfg_file  = "configFSI.cfg"
    pywrapper_fsi2d.test_iter = 4
    pywrapper_fsi2d.test_vals = [4, 0, -3.743210, -4.133483] #last 4 columns
    pywrapper_fsi2d.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--nZone 2 --fsi True --parallel -f")
    pywrapper_fsi2d.unsteady  = True
    pywrapper_fsi2d.multizone = True
    test_list.append(pywrapper_fsi2d)

    # Unsteady CHT
    pywrapper_unsteadyCHT               = TestCase('pywrapper_unsteadyCHT')
    pywrapper_unsteadyCHT.cfg_dir       = "py_wrapper/flatPlate_unsteady_CHT"
    pywrapper_unsteadyCHT.cfg_file      = "unsteady_CHT_FlatPlate_Conf.cfg"
    pywrapper_unsteadyCHT.test_iter     = 5
    pywrapper_unsteadyCHT.test_vals     = [-1.614167, 2.245725, -0.001241, 0.175713]
    pywrapper_unsteadyCHT.command       = TestCase.Command("mpirun -np 2", "python", "launch_unsteady_CHT_FlatPlate.py --parallel -f")
    pywrapper_unsteadyCHT.unsteady      = True
    pywrapper_unsteadyCHT.new_output    = True
    test_list.append(pywrapper_unsteadyCHT)

    # Rigid motion
    pywrapper_rigidMotion               = TestCase('pywrapper_rigidMotion')
    pywrapper_rigidMotion.cfg_dir       = "py_wrapper/flatPlate_rigidMotion"
    pywrapper_rigidMotion.cfg_file      = "flatPlate_rigidMotion_Conf.cfg"
    pywrapper_rigidMotion.test_iter     = 5
    pywrapper_rigidMotion.test_vals     = [-1.614170, 2.242953, 0.350037, 0.093116]
    pywrapper_rigidMotion.command       = TestCase.Command("mpirun -np 2", "python", "launch_flatPlate_rigidMotion.py --parallel -f")
    pywrapper_rigidMotion.unsteady      = True
    test_list.append(pywrapper_rigidMotion)

    ##############################################
    ### Method of Manufactured Solutions (MMS) ###
    ##############################################

    # FVM, compressible, laminar N-S
    mms_fvm_ns           = TestCase('mms_fvm_ns')
    mms_fvm_ns.cfg_dir   = "mms/fvm_navierstokes"
    mms_fvm_ns.cfg_file  = "lam_mms_roe.cfg"
    mms_fvm_ns.test_iter = 20
    mms_fvm_ns.test_vals = [-2.851428, 2.192348, 0.000000, 0.000000] #last 4 columns
    mms_fvm_ns.tol       = 0.0001
    test_list.append(mms_fvm_ns)

    # FVM, incompressible, euler
    mms_fvm_inc_euler           = TestCase('mms_fvm_inc_euler')
    mms_fvm_inc_euler.cfg_dir   = "mms/fvm_incomp_euler"
    mms_fvm_inc_euler.cfg_file  = "inv_mms_jst.cfg"
    mms_fvm_inc_euler.test_iter = 20
    mms_fvm_inc_euler.test_vals = [-9.128660, -9.441806, 0.000000, 0.000000] #last 4 columns
    mms_fvm_inc_euler.tol       = 0.0001
    test_list.append(mms_fvm_inc_euler)

    # FVM, incompressible, laminar N-S
    mms_fvm_inc_ns           = TestCase('mms_fvm_inc_ns')
    mms_fvm_inc_ns.cfg_dir   = "mms/fvm_incomp_navierstokes"
    mms_fvm_inc_ns.cfg_file  = "lam_mms_fds.cfg"
    mms_fvm_inc_ns.test_iter = 20
    mms_fvm_inc_ns.test_vals = [-7.414944, -7.631546, 0.000000, 0.000000] #last 4 columns
    mms_fvm_inc_ns.tol       = 0.0001
    test_list.append(mms_fvm_inc_ns)

    # DG, compressible, euler
    ringleb_dg_euler           = TestCase('ringleb_dg_euler')
    ringleb_dg_euler.cfg_dir   = "mms/dg_ringleb"
    ringleb_dg_euler.cfg_file  = "ringleb_dg.cfg"
    ringleb_dg_euler.test_iter = 100
    ringleb_dg_euler.test_vals = [-5.136652, -4.724941, 0.000000, 0.000000] #last 4 columns
    ringleb_dg_euler.command   = TestCase.Command(exec = "SU2_CFD")
    ringleb_dg_euler.tol       = 0.0001
    test_list.append(ringleb_dg_euler)

    # DG, compressible, laminar N-S
    mms_dg_ns           = TestCase('mms_dg_ns')
    mms_dg_ns.cfg_dir   = "mms/dg_navierstokes"
    mms_dg_ns.cfg_file  = "lam_mms_dg.cfg"
    mms_dg_ns.test_iter = 100
    mms_dg_ns.test_vals = [-1.845393, 3.520699, 0.000000, 0.000000] #last 4 columns
    mms_dg_ns.command   = TestCase.Command(exec = "SU2_CFD")
    mms_dg_ns.tol       = 0.0001
    test_list.append(mms_dg_ns)

    # DG, compressible, laminar N-S 3D
    mms_dg_ns_3d           = TestCase('mms_dg_ns_3d')
    mms_dg_ns_3d.cfg_dir   = "mms/dg_navierstokes_3d"
    mms_dg_ns_3d.cfg_file  = "lam_mms_dg_3d.cfg"
    mms_dg_ns_3d.test_iter = 100
    mms_dg_ns_3d.test_vals = [-0.146826, 5.356413, 0.000000, 0.000000] #last 4 columns
    mms_dg_ns_3d.command   = TestCase.Command(exec = "SU2_CFD")
    mms_dg_ns_3d.tol       = 0.0001
    test_list.append(mms_dg_ns_3d)

    #####################
    ## Species solver ###
    #####################

    # 2 species (1 eq) primitive venturi mixing using mixing model
    species2_primitiveVenturi_mixingmodel           = TestCase('species2_primitiveVenturi_mixingmodel')
    species2_primitiveVenturi_mixingmodel.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel.cfg_file  = "species2_primitiveVenturi_mixingmodel.cfg"
    species2_primitiveVenturi_mixingmodel.test_iter = 50
    species2_primitiveVenturi_mixingmodel.test_vals = [-5.477236, -4.589782, -4.582531, -5.690416, -0.068943, -5.631254, 5.000000, -1.668303, 5.000000, -4.866978, 5.000000, -1.168368, 0.000425, 0.000404, 0.000021, 0.000000]
    species2_primitiveVenturi_mixingmodel.new_output = True
    test_list.append(species2_primitiveVenturi_mixingmodel)

    # 2 species (1 eq) primitive venturi mixing using mixing model and bounded scalar transport
    species2_primitiveVenturi_mixingmodel_boundedscalar           = TestCase('species2_primitiveVenturi_mixingmodel_boundedscalar')
    species2_primitiveVenturi_mixingmodel_boundedscalar.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_boundedscalar.cfg_file  = "species2_primitiveVenturi_mixingmodel_boundedscalar.cfg"
    species2_primitiveVenturi_mixingmodel_boundedscalar.test_iter = 50
    species2_primitiveVenturi_mixingmodel_boundedscalar.test_vals = [-5.419823, -4.490921, -4.496340, -5.724950, -0.120331, -5.693224, 5.000000, -1.766882, 5.000000, -4.958487, 5.000000, -2.150247, 0.000300, 0.000300, 0.000000, 0.000000]
    species2_primitiveVenturi_mixingmodel_boundedscalar.new_output = True
    test_list.append(species2_primitiveVenturi_mixingmodel_boundedscalar)

    # 2 species (1 eq) primitive venturi mixing using mixing model including viscosity and thermal conductivity
    species2_primitiveVenturi_mixingmodel_viscosity           = TestCase('species2_primitiveVenturi_mixingmodel_viscosity')
    species2_primitiveVenturi_mixingmodel_viscosity.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_viscosity.cfg_file  = "species2_primitiveVenturi_mixingmodel_viscosity.cfg"
    species2_primitiveVenturi_mixingmodel_viscosity.test_iter = 50
    species2_primitiveVenturi_mixingmodel_viscosity.test_vals = [-4.738778, -4.325766, -4.610914, -5.834431, 0.521730, -4.934375, 5.000000, -1.887191, 5.000000, -5.499917, 5.000000, -1.770845, 2.292904, 0.971941, 0.608500, 0.712464]
    species2_primitiveVenturi_mixingmodel_viscosity.new_output = True
    test_list.append(species2_primitiveVenturi_mixingmodel_viscosity)

    # 2 species (1 eq) primitive venturi mixing using mixing model including heat capacity and mass diffusivity
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2           = TestCase('species2_primitiveVenturi_mixingmodel_heatcapacity_H2.cfg')
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.cfg_file  = "species2_primitiveVenturi_mixingmodel_heatcapacity_H2.cfg"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.test_iter = 50
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.test_vals = [-6.105905, -4.979417, -4.870704, -7.338344, 2.435662, -5.621107, 30.000000, -5.737172, 12.000000, -8.148164, 9.000000, -8.040415, 2.084589, 1.000000, 0.600000, 0.484589]
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.su2_exec  = "mpirun -n 2 SU2_CFD"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.timeout   = 1600
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.new_output = True
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.tol       = 0.00001
    test_list.append(species2_primitiveVenturi_mixingmodel_heatcapacity_H2)

    # 2 species (1 eq) primitive venturi mixing using mixing model including heat capacity and mass diffusivity NonDimensional case
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND           = TestCase('species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.cfg')
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.cfg_file  = "species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.cfg"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.test_iter = 50
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.test_vals = [-5.711760, -5.284237, -5.175532, -8.342082, 2.130839, -5.226944, 30.000000, -5.737157, 12.000000, -8.148020, 9.000000, -8.040460, 2.084591, 1.000000, 0.600000, 0.484591]
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.su2_exec  = "mpirun -n 2 SU2_CFD"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.timeout   = 1600
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.new_output = True
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.tol       = 0.00001
    test_list.append(species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND)

    # 2 species (1 eq) primitive venturi mixing
    species2_primitiveVenturi           = TestCase('species2_primitiveVenturi')
    species2_primitiveVenturi.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi.cfg_file  = "species2_primitiveVenturi.cfg"
    species2_primitiveVenturi.test_iter = 50
    species2_primitiveVenturi.test_vals = [-6.004340, -5.236017, -5.080329, -5.888785, -1.554476, -6.076144, 5.000000, -0.808416, 5.000000, -2.325030, 5.000000, -0.238222, 0.000089, 0.000088, 0.000001, 0.000000]
    species2_primitiveVenturi.new_output = True
    test_list.append(species2_primitiveVenturi)

    # 2 species (1 eq) primitive venturi mixing with bounded scalar transport
    species_primitiveVenturi_boundedscalar             = TestCase('species2_primitiveVenturi_bounded_scalar')
    species_primitiveVenturi_boundedscalar.cfg_dir     = "species_transport/venturi_primitive_3species"
    species_primitiveVenturi_boundedscalar.cfg_file    = "species2_primitiveVenturi_boundedscalar.cfg"
    species_primitiveVenturi_boundedscalar.test_iter   = 50
    species_primitiveVenturi_boundedscalar.test_vals   = [-5.297585, -4.397797, -4.377086, -5.593131, -1.011782, -5.623540, 5.000000, -1.775123, 5.000000, -4.086339, 5.000000, -2.080187, 0.000424, 0.000424, 0.000000, 0.000000]
    species_primitiveVenturi_boundedscalar.new_output  = True
    test_list.append(species_primitiveVenturi_boundedscalar)

    # 2 species (1 eq) primitive venturi mixing using mixing model including inlet markers for turbulent intensity and viscosity ratios
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS           = TestCase('species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.cfg')
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.cfg_file  = "species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.cfg"
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.test_iter = 50
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.test_vals = [-4.019621, -1.652733, -1.414138, -0.971744, 1.590456, -3.762843, 23.000000, -5.066221, 12.000000, -5.359931, 4.000000, -6.078888, 2.000000, 1.000000, 0.000000, 1.000000]
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.su2_exec  = "mpirun -n 2 SU2_CFD"
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.timeout   = 1600
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.new_output = True
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.tol       = 0.00001
    test_list.append(species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS)

    # 3 species (2 eq) primitive venturi mixing with inlet files.
    # Note that the residuals are exactly the same as for the non-inlet case which should be the case for a fresh inlet file.
    species3_primitiveVenturi_inletFile           = TestCase('species3_primitiveVenturi_inletFile')
    species3_primitiveVenturi_inletFile.cfg_dir   = "species_transport/venturi_primitive_3species"
    species3_primitiveVenturi_inletFile.cfg_file  = "species3_primitiveVenturi_inletFile.cfg"
    species3_primitiveVenturi_inletFile.test_iter = 50
    species3_primitiveVenturi_inletFile.test_vals = [-6.074971, -5.306648, -5.150960, -5.959416, -1.625107, -6.343704, -6.460033, 5.000000, -0.808413, 5.000000, -2.325029, 5.000000, -0.274923]
    species3_primitiveVenturi_inletFile.new_output = True
    test_list.append(species3_primitiveVenturi_inletFile)

    # rectangle passive transport validation
    species_passive_val           = TestCase('species_passive_val')
    species_passive_val.cfg_dir   = "species_transport/passive_transport_validation"
    species_passive_val.cfg_file  = "passive_transport.cfg"
    species_passive_val.test_iter = 50
    species_passive_val.test_vals =         [-16.559189, -16.315116, -16.908670, -4.257599, 10, -4.523292, 8, -5.19335, 0.18661, 0]
    species_passive_val.test_vals_aarch64 = [-16.538551, -16.312552, -16.882823, -4.257599, 10, -4.585464, 8, -5.19335, 0.18661, 0]
    species_passive_val.new_output = True
    test_list.append(species_passive_val)

    # species transport, 3 species with multizone (2 fluid regions)
    species3_multizone_restart           = TestCase('species3_multizone_restart')
    species3_multizone_restart.cfg_dir   = "species_transport/multizone"
    species3_multizone_restart.cfg_file  = "configMaster.cfg"
    species3_multizone_restart.test_iter = 5
    species3_multizone_restart.test_vals = [-7.169637, -6.381576]
    species3_multizone_restart.new_output = True
    species3_multizone_restart.multizone = True
    test_list.append(species3_multizone_restart)

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

    ######################################
    ### RUN SU2_SOL TESTS              ###
    ######################################

    # parallel STL output using
    stl_writer_test                = TestCase('stl_writer_test')
    stl_writer_test.cfg_dir        = "rans/oneram6"
    stl_writer_test.cfg_file       = "turb_ONERAM6.cfg"
    stl_writer_test.test_iter      = 1
    stl_writer_test.command        = TestCase.Command("mpirun -n 2", "SU2_SOL")
    stl_writer_test.timeout        = 1600
    stl_writer_test.reference_file = "surface_flow.stl.ref"
    stl_writer_test.test_file      = "surface_flow.stl"
    pass_list.append(stl_writer_test.run_filediff())
    test_list.append(stl_writer_test)

    ######################################
    ### RUN SU2_DEF TESTS              ###
    ######################################

    # Inviscid NACA0012 (triangles)
    naca0012_def            = TestCase('naca0012_def')
    naca0012_def.cfg_dir   = "deformation/naca0012"
    naca0012_def.cfg_file  = "def_NACA0012.cfg"
    naca0012_def.test_iter = 10
    naca0012_def.test_vals = [0.00352488] #residual
    naca0012_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    naca0012_def.timeout   = 1600
    naca0012_def.tol       = 1e-8

    pass_list.append(naca0012_def.run_def())
    test_list.append(naca0012_def)

    # Inviscid NACA0012 based on SURFACE_FILE input (surface_bump.dat)
    naca0012_def_file            = TestCase('naca0012_def_file')
    naca0012_def_file.cfg_dir   = "deformation/naca0012"
    naca0012_def_file.cfg_file  = "surface_file_NACA0012.cfg"
    naca0012_def_file.test_iter = 10
    naca0012_def_file.test_vals = [0.00352488] #residual
    naca0012_def_file.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    naca0012_def_file.timeout   = 1600
    naca0012_def_file.tol       = 1e-8

    pass_list.append(naca0012_def_file.run_def())
    test_list.append(naca0012_def_file)

    # RAE2822 (mixed tris + quads)
    rae2822_def            = TestCase('rae2822_def')
    rae2822_def.cfg_dir   = "deformation/rae2822"
    rae2822_def.cfg_file  = "def_RAE2822.cfg"
    rae2822_def.test_iter = 10
    rae2822_def.test_vals = [8.24002e-09] #residual
    rae2822_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    rae2822_def.timeout   = 1600
    rae2822_def.tol       = 1e-13

    pass_list.append(rae2822_def.run_def())
    test_list.append(rae2822_def)

    # Turb NACA4412 (quads, wall distance)
    naca4412_def            = TestCase('naca4412_def')
    naca4412_def.cfg_dir   = "deformation/naca4412"
    naca4412_def.cfg_file  = "def_NACA4412.cfg"
    naca4412_def.test_iter = 10
    naca4412_def.test_vals = [2.210380e-12] #residual
    naca4412_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    naca4412_def.timeout   = 1600
    naca4412_def.tol       = 1e-12

    pass_list.append(naca4412_def.run_def())
    test_list.append(naca4412_def)

    # Brick of tets (inverse volume)
    brick_tets_def            = TestCase('brick_tets_def')
    brick_tets_def.cfg_dir   = "deformation/brick_tets"
    brick_tets_def.cfg_file  = "def_brick_tets.cfg"
    brick_tets_def.test_iter = 10
    brick_tets_def.test_vals = [0.000955394] #residual
    brick_tets_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    brick_tets_def.timeout   = 1600
    brick_tets_def.tol       = 1e-9

    pass_list.append(brick_tets_def.run_def())
    test_list.append(brick_tets_def)

    # Brick of isotropic hexas (inverse volume)
    brick_hex_def           = TestCase('brick_hex_def')
    brick_hex_def.cfg_dir   = "deformation/brick_hex"
    brick_hex_def.cfg_file  = "def_brick_hex.cfg"
    brick_hex_def.test_iter = 10
    brick_hex_def.test_vals = [0.000166575] #residual
    brick_hex_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    brick_hex_def.timeout   = 1600
    brick_hex_def.tol       = 1e-9

    pass_list.append(brick_hex_def.run_def())
    test_list.append(brick_hex_def)

    # Brick with a pyramid layer (inverse volume)
    brick_pyra_def           = TestCase('brick_pyra_def')
    brick_pyra_def.cfg_dir   = "deformation/brick_pyra"
    brick_pyra_def.cfg_file  = "def_brick_pyra.cfg"
    brick_pyra_def.test_iter = 10
    brick_pyra_def.test_vals = [0.00161454] #residual
    brick_pyra_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    brick_pyra_def.timeout   = 1600
    brick_pyra_def.tol       = 1e-8

    pass_list.append(brick_pyra_def.run_def())
    test_list.append(brick_pyra_def)

    # Brick of isotropic prisms (inverse volume)
    brick_prism_def           = TestCase('brick_prism_def')
    brick_prism_def.cfg_dir   = "deformation/brick_prism"
    brick_prism_def.cfg_file  = "def_brick_prism.cfg"
    brick_prism_def.test_iter = 10
    brick_prism_def.test_vals = [0.00254732] #residual
    brick_prism_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    brick_prism_def.timeout   = 1600
    brick_prism_def.tol       = 1e-8

    pass_list.append(brick_prism_def.run_def())
    test_list.append(brick_prism_def)

    # Brick of prisms with high aspect ratio cells near the wall (wall distance)
    brick_prism_rans_def           = TestCase('brick_prism_rans_def')
    brick_prism_rans_def.cfg_dir   = "deformation/brick_prism_rans"
    brick_prism_rans_def.cfg_file  = "def_brick_prism_rans.cfg"
    brick_prism_rans_def.test_iter = 10
    brick_prism_rans_def.test_vals = [2.99461e-07] #residual
    brick_prism_rans_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    brick_prism_rans_def.timeout   = 1600
    brick_prism_rans_def.tol       = 1e-12

    pass_list.append(brick_prism_rans_def.run_def())
    test_list.append(brick_prism_rans_def)

    # Brick of hexas with high aspect ratio cells near the wall (inverse volume)
    brick_hex_rans_def           = TestCase('brick_hex_rans_def')
    brick_hex_rans_def.cfg_dir   = "deformation/brick_hex_rans"
    brick_hex_rans_def.cfg_file  = "def_brick_hex_rans.cfg"
    brick_hex_rans_def.test_iter = 10
    brick_hex_rans_def.test_vals = [3.54213e-06] #residual
    brick_hex_rans_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    brick_hex_rans_def.timeout   = 1600
    brick_hex_rans_def.tol       = 1e-11

    pass_list.append(brick_hex_rans_def.run_def())
    test_list.append(brick_hex_rans_def)

    # Cylindrical FFD test
    cylinder_ffd_def           = TestCase('cylinder_ffd_def')
    cylinder_ffd_def.cfg_dir   = "deformation/cylindrical_ffd"
    cylinder_ffd_def.cfg_file  = "def_cylindrical.cfg"
    cylinder_ffd_def.test_iter = 10
    cylinder_ffd_def.test_vals = [0.000902348] #residual
    cylinder_ffd_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    cylinder_ffd_def.timeout   = 1600
    cylinder_ffd_def.tol       = 1e-9

    pass_list.append(cylinder_ffd_def.run_def())
    test_list.append(cylinder_ffd_def)

    # Spherical FFD test
    sphere_ffd_def           = TestCase('sphere_ffd_def')
    sphere_ffd_def.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def.cfg_file  = "def_spherical.cfg"
    sphere_ffd_def.test_iter = 10
    sphere_ffd_def.test_vals = [0.00360367] #residual
    sphere_ffd_def.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    sphere_ffd_def.timeout   = 1600
    sphere_ffd_def.tol       = 1e-8

    pass_list.append(sphere_ffd_def.run_def())
    test_list.append(sphere_ffd_def)

    # Spherical FFD test using BSplines
    sphere_ffd_def_bspline           = TestCase('sphere_ffd_def_bspline')
    sphere_ffd_def_bspline.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def_bspline.cfg_file  = "def_spherical_bspline.cfg"
    sphere_ffd_def_bspline.test_iter = 10
    sphere_ffd_def_bspline.test_vals = [0.00208393] #residual
    sphere_ffd_def_bspline.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    sphere_ffd_def_bspline.timeout   = 1600
    sphere_ffd_def_bspline.tol       = 1e-8

    pass_list.append(sphere_ffd_def_bspline.run_def())
    test_list.append(sphere_ffd_def_bspline)

    # 2D FD streamwise periodic cht, avg temp obj func
    fd_sp_pinArray_cht_2d_dp_hf                = TestCase('fd_sp_pinArray_cht_2d_dp_hf')
    fd_sp_pinArray_cht_2d_dp_hf.cfg_dir        = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    fd_sp_pinArray_cht_2d_dp_hf.cfg_file       = "FD_configMaster.cfg"
    fd_sp_pinArray_cht_2d_dp_hf.test_iter      = 100
    fd_sp_pinArray_cht_2d_dp_hf.command        = TestCase.Command(exec = "finite_differences.py", param = "-z 2 -n 2 -f")
    fd_sp_pinArray_cht_2d_dp_hf.timeout        = 1600
    fd_sp_pinArray_cht_2d_dp_hf.reference_file         = "of_grad_findiff.csv.ref"
    fd_sp_pinArray_cht_2d_dp_hf.reference_file_aarch64 = "of_grad_findiff_aarch64.csv.ref"
    fd_sp_pinArray_cht_2d_dp_hf.test_file      = "FINDIFF/of_grad_findiff.csv"
    fd_sp_pinArray_cht_2d_dp_hf.multizone      = True

    pass_list.append(fd_sp_pinArray_cht_2d_dp_hf.run_filediff())
    test_list.append(fd_sp_pinArray_cht_2d_dp_hf)


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
