#!/usr/bin/env python

## \file parallel_regression.py
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

def main():
    '''This program runs SU2 and ensures that the output matches specified values.
       This will be used to do checks when code is pushed to github
       to make sure nothing is broken. '''

    test_list = []

    #######################
    ### Flamelet solver ###
    #######################

    # 2D planar laminar premixed methane flame on isothermal burner (restart)
    cfd_flamelet_ch4 = TestCase('cfd_flamelet_ch4')
    cfd_flamelet_ch4.cfg_dir = "flamelet/01_laminar_premixed_ch4_flame_cfd"
    cfd_flamelet_ch4.cfg_file = "lam_prem_ch4_cfd.cfg"
    cfd_flamelet_ch4.test_iter = 10
    cfd_flamelet_ch4.test_vals = [-14.297190, -12.593356, -14.249342, -6.069757, -14.816955, -17.057821]
    cfd_flamelet_ch4.new_output = True
    test_list.append(cfd_flamelet_ch4)

   # axisymmetric 2D planar laminar premixed methane flame on isothermal burner (restart)
    cfd_flamelet_ch4_axi = TestCase('cfd_flamelet_ch4_axi')
    cfd_flamelet_ch4_axi.cfg_dir = "flamelet/05_laminar_premixed_ch4_flame_cfd_axi"
    cfd_flamelet_ch4_axi.cfg_file = "lam_prem_ch4_cfd_axi.cfg"
    cfd_flamelet_ch4_axi.test_iter = 10
    cfd_flamelet_ch4_axi.test_vals = [-10.935482, -11.477701, -11.302709, -12.734435, -6.110542]
    cfd_flamelet_ch4_axi.new_output = True
    test_list.append(cfd_flamelet_ch4_axi)

    # 2D planar laminar partially premixed flame on isothermal burner and heat exchanger (restart)
    cfd_flamelet_ch4_partial_premix = TestCase('cfd_flamelet_ch4_partial_premix')
    cfd_flamelet_ch4_partial_premix.cfg_dir = "flamelet/06_laminar_partial_premixed_ch4_flame_cfd"
    cfd_flamelet_ch4_partial_premix.cfg_file = "lam_partial_prem_ch4_cfd.cfg"
    cfd_flamelet_ch4_partial_premix.test_iter = 10
    cfd_flamelet_ch4_partial_premix.test_vals = [-9.645814, -11.304194, -3.675580, -13.158500, -11.087904]
    cfd_flamelet_ch4_partial_premix.new_output = True
    test_list.append(cfd_flamelet_ch4_partial_premix)

    # 2D planar laminar premixed hydrogen flame on isothermal burner with heat exchanger emulator (restart)
    cfd_flamelet_h2 = TestCase('cfd_flamelet_h2')
    cfd_flamelet_h2.cfg_dir = "flamelet/07_laminar_premixed_h2_flame_cfd"
    cfd_flamelet_h2.cfg_file = "laminar_premixed_h2_flame_cfd.cfg"
    cfd_flamelet_h2.test_iter = 5
    cfd_flamelet_h2.test_vals = [-10.006096, -9.844001, -3.290043, -11.338465]
    cfd_flamelet_h2.new_output = True
    test_list.append(cfd_flamelet_h2)

    #########################
    ## NEMO solver ###
    #########################

    # Adiabatic thermal bath
    thermalbath = TestCase('thermalbath')
    thermalbath.cfg_dir = "nonequilibrium/thermalbath/finitechemistry"
    thermalbath.cfg_file = "thermalbath.cfg"
    thermalbath.test_iter = 10
    thermalbath.test_vals = [0.945997, 0.945997, -12.039262, -12.171767, -32.000000, 10.013239]
    test_list.append(thermalbath)

    # Adiabatic thermal bath
    ionized = TestCase('ionized')
    ionized.cfg_dir = "nonequilibrium/thermalbath/finitechemistry"
    ionized.cfg_file = "weakly_ionized.cfg"
    ionized.test_iter = 10
    ionized.test_vals = [-29.806157, -11.130797, -11.337264, -17.235059, -17.578729, -15.190274, -25.013626, -32.000000, -5.174887, 0.000000, 0.000000]
    ionized.test_vals_aarch64 = [-29.816386, -10.729986, -11.720016, -17.484469, -18.237891, -15.241605, -24.956918, -32.000000, -5.727244, 0.000000, 0.000000]
    test_list.append(ionized)

    # Adiabatic frozen thermal bath
    thermalbath_frozen = TestCase('thermalbath_frozen')
    thermalbath_frozen.cfg_dir = "nonequilibrium/thermalbath/frozen"
    thermalbath_frozen.cfg_file = "thermalbath_frozen.cfg"
    thermalbath_frozen.test_iter = 10
    thermalbath_frozen.test_vals = [-32.000000, -32.000000, -11.962477, -11.962477, -32.000000, 10.013545]
    test_list.append(thermalbath_frozen)

    # Inviscid single wedge, ausm, implicit
    invwedge_a = TestCase('invwedge_ausm')
    invwedge_a.cfg_dir = "nonequilibrium/invwedge"
    invwedge_a.cfg_file = "invwedge_ausm.cfg"
    invwedge_a.test_iter = 10
    invwedge_a.test_vals = [-1.069675, -1.594438, -18.299923, -18.627316, -18.573325, 2.245721, 1.874105, 5.290285, 0.847729]
    invwedge_a.test_vals_aarch64 = [-1.070904, -1.595667, -18.299980, -18.627372, -18.573382, 2.244654, 1.871030, 5.289134, 0.846502]
    test_list.append(invwedge_a)

    # Inviscid single wedge, ausm+-up2, implicit
    invwedge_ap2 = TestCase('invwedge_ap2')
    invwedge_ap2.cfg_dir = "nonequilibrium/invwedge"
    invwedge_ap2.cfg_file = "invwedge_ausmplusup2.cfg"
    invwedge_ap2.test_iter = 10
    invwedge_ap2.test_vals = [-0.982059, -1.506822, -16.735977, -17.063993, -17.009083, 2.354326, 1.482256, 5.373931, 0.927155]
    invwedge_ap2.test_vals_aarch64 = [-0.982059, -1.506822, -16.735977, -17.063993, -17.009083, 2.354326, 1.482256, 5.373931, 0.927155]
    test_list.append(invwedge_ap2)

    # Inviscid single wedge, msw, implicit
    invwedge_msw = TestCase('invwedge_msw')
    invwedge_msw.cfg_dir = "nonequilibrium/invwedge"
    invwedge_msw.cfg_file = "invwedge_msw.cfg"
    invwedge_msw.test_iter = 10
    invwedge_msw.test_vals = [-1.212335, -1.737098, -18.299220, -18.626618, -18.572623, 2.106171, 1.651949, 5.143958, 0.704444]
    invwedge_msw.test_vals_aarch64 = [-1.224649, -1.749412, -18.299151, -18.626550, -18.572552, 2.094106, 1.635779, 5.131012, 0.692821]
    test_list.append(invwedge_msw)

    # Inviscid single wedge, roe, implicit
    invwedge_roe = TestCase('invwedge_roe')
    invwedge_roe.cfg_dir = "nonequilibrium/invwedge"
    invwedge_roe.cfg_file = "invwedge_roe.cfg"
    invwedge_roe.test_iter = 10
    invwedge_roe.test_vals = [-1.062496, -1.587259, -17.208314, -17.538065, -17.481444, 2.255933, 1.853450, 5.293725, 0.890784]
    invwedge_roe.test_vals_aarch64 = [-1.069128, -1.593891, -17.208222, -17.537969, -17.481352, 2.249020, 1.852904, 5.287143, 0.879852]
    test_list.append(invwedge_roe)

    # Inviscid single wedge, lax, implicit
    invwedge_lax = TestCase('invwedge_lax')
    invwedge_lax.cfg_dir = "nonequilibrium/invwedge"
    invwedge_lax.cfg_file = "invwedge_lax.cfg"
    invwedge_lax.test_iter = 10
    invwedge_lax.test_vals = [-0.877280, -1.402043, -32.000000, -32.000000, -24.952631, 2.451869, 1.857084, 5.486158, 1.051580]
    invwedge_lax.test_vals_aarch64 = [-0.877280, -1.402043, -32.000000, -32.000000, -24.952631, 2.451869, 1.857084, 5.486158, 1.051580]
    test_list.append(invwedge_lax)

    # Inviscid single wedge, implicit, AUSM+M scheme
    invwedge_ausm_m = TestCase('invwedge_ausm_m')
    invwedge_ausm_m.cfg_dir = "nonequilibrium/invwedge"
    invwedge_ausm_m.cfg_file = "invwedge_am.cfg"
    invwedge_ausm_m.test_iter = 10
    invwedge_ausm_m.test_vals = [-1.173033, -1.697796, -16.739586, -17.063491, -17.012692, 2.124519, 1.963804, 5.182881, 0.747539]
    invwedge_ausm_m.test_vals_aarch64 = [-1.171654, -1.696417, -16.739585, -17.063491, -17.012691, 2.125633, 1.966511, 5.184281, 0.749068]
    test_list.append(invwedge_ausm_m)

    # Inviscid single wedge, implicit, NEMO supersonic inlet
    invwedge_ss_inlet = TestCase('invwedge_ss_inlet')
    invwedge_ss_inlet.cfg_dir = "nonequilibrium/invwedge"
    invwedge_ss_inlet.cfg_file = "invwedge_ss_inlet.cfg"
    invwedge_ss_inlet.test_iter = 10
    invwedge_ss_inlet.test_vals = [-1.068592, -1.593355, -18.250183, -18.579524, -18.523255, 2.246972, 1.874197, 5.291273, 0.848771]
    invwedge_ss_inlet.test_vals_aarch64 = [-1.069892, -1.594654, -18.250175, -18.579516, -18.523248, 2.245827, 1.871123, 5.290054, 0.847476]
    test_list.append(invwedge_ss_inlet)

    # Viscous single cone - axisymmetric
    visc_cone = TestCase('visc_cone')
    visc_cone.cfg_dir = "nonequilibrium/visc_wedge"
    visc_cone.cfg_file = "axi_visccone.cfg"
    visc_cone.test_iter = 10
    visc_cone.test_vals = [-5.222270, -5.746525, -20.560278, -20.510152, -20.409102, 1.255758, -3.208382, -0.016014, 0.093462, 32619.000000]
    visc_cone.test_vals_aarch64 = [-5.222267, -5.746523, -20.560279, -20.510152, -20.409102, 1.255758, -3.208380, -0.016014, 0.093462, 32633.000000]
    test_list.append(visc_cone)

    # Viscous single wedge with Mutation++
    #viscwedge_mpp = TestCase('viscwedge_mpp')
    #viscwedge_mpp.cfg_dir = "nonequilibrium/viscwedge_mpp"
    #viscwedge_mpp.cfg_file = "viscwedge_mpp.cfg"
    #viscwedge_mpp.test_iter = 10
    #viscwedge_mpp.test_vals = [-20.608474, -20.586446,-20.707524, -5.171304,-5.696067,-1.548350,-2.071211,2.231054,-2.545494]
    #test_list.append(viscwedge_mpp)

    # Viscous single wedge - super catalytic walls
    super_cat = TestCase('super_cat')
    super_cat.cfg_dir = "nonequilibrium/visc_wedge"
    super_cat.cfg_file = "super_cat.cfg"
    super_cat.test_iter = 10
    super_cat.test_vals = [-5.232595, -5.757889, -20.641415, -20.640623, -20.541670, 1.246866, -3.205258, -0.028372, 0.250647, 32440.000000]
    test_list.append(super_cat)

    # Viscous single wedge - partially catalytic walls
    partial_cat = TestCase('partial_cat')
    partial_cat.cfg_dir = "nonequilibrium/visc_wedge"
    partial_cat.cfg_file = "partial_cat.cfg"
    partial_cat.test_iter = 10
    partial_cat.test_vals = [-5.210302, -5.735065, -20.880448, -20.825971, -23.475263, 1.806201, -2.813952, -0.078400, 0.495606, 29020.000000]
    test_list.append(partial_cat)

    # Viscous cylinder, ionization, Gupta-Yos
    ion_gy = TestCase('ion_gy')
    ion_gy.cfg_dir = "nonequilibrium/visc_cylinder"
    ion_gy.cfg_file = "cyl_ion_gy.cfg"
    ion_gy.test_iter = 10
    ion_gy.test_vals = [-11.629873, -4.165563, -4.702662, -4.950351, -5.146155, -4.993878, -6.893332, 5.990109, 5.990004, -0.014849, 0.000000, 90090.000000]
    test_list.append(ion_gy)

    ##########################
    ### Compressible Euler ###
    ##########################

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 20
    channel.test_vals = [-2.904401, 2.536139, 0.020924, 0.042340]
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.607257, -4.138750, 0.327682, 0.022685]
    test_list.append(naca0012)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-1.387215, 4.281574, -0.245443, 0.043264]
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-11.500303, -10.971884, 0.280800, 0.008623]
    oneram6.timeout   = 3200
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-3.840915, 1.693979, 0.301144, 0.019489]
    test_list.append(fixedCL_naca0012)

    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals         = [-1.096894, 4.372054, 0.002074, 0.031515]
    polar_naca0012.test_vals_aarch64 = [-1.083394, 4.386134, 0.001588, 0.033513]
    polar_naca0012.command   = TestCase.Command(exec = "compute_polar.py", param = "-i 11")
    # flaky test on arm64
    polar_naca0012.enabled_on_cpu_arch = ["x86_64"]
    test_list.append(polar_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.475144, 6.834602, -0.000007, 1.783980]
    test_list.append(bluntbody)

    # Equivalent area NACA64-206
    ea_naca64206           = TestCase('ea_naca64206')
    ea_naca64206.cfg_dir   = "optimization_euler/equivalentarea_naca64206"
    ea_naca64206.cfg_file  = "NACA64206.cfg"
    ea_naca64206.test_iter = 10
    ea_naca64206.test_vals = [-1.125893, -0.474056, -0.002414, 67775.000000]
    test_list.append(ea_naca64206)

    # SUPERSONIC FLOW PAST A RAMP IN A CHANNEL
    ramp = TestCase('ramp')
    ramp.cfg_dir = "euler/ramp"
    ramp.cfg_file = "inv_ramp.cfg"
    ramp.test_iter = 10
    ramp.test_vals = [-13.650727, -8.014514, -0.076277, 0.054839]
    ramp.test_vals_aarch64 = [-13.398422, -7.786461, -0.081064, 0.056474]
    test_list.append(ramp)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 100
    flatplate.test_vals = [-7.613292, -2.141207, 0.001084, 0.036230, 2.361500, -2.325300, 0.000000, 0.000000]
    test_list.append(flatplate)

    # Custom objective function
    flatplate_udobj           = TestCase('flatplate_udobj')
    flatplate_udobj.cfg_dir   = "user_defined_functions"
    flatplate_udobj.cfg_file  = "lam_flatplate.cfg"
    flatplate_udobj.test_iter = 20
    flatplate_udobj.test_vals = [-6.760101, -1.283906, -0.745653, 0.000587, -0.000038, 0.000977, -0.001015, 596.450000, 299.550000, 296.900000, 21.318000, 0.586640, 36.553000, 2.188800]
    test_list.append(flatplate_udobj)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-8.422091, -2.930561, -0.003396, 1.608418, 0.000000]
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.841604, -1.379532, -1.266739, 76.118218, 0.000000]
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-5.050889, 0.648196, 0.000199, 13.639173, 0.000000]
    poiseuille.tol       = 0.001
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals = [-12.007512, -7.227061, -0.000000, 2.089953]
    poiseuille_profile.test_vals_aarch64 = [-12.492864, -7.671632, -0.000000, 2.085796]
    test_list.append(poiseuille_profile)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.004689, -5.265797, 0.809463, 0.062016, 0]
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510313, 4.941942, 0.813793, 0.062425, 0.000000]
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.677866, 4.941942, 0.813793, 0.062425]
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.297198, -6.730442, -0.187632, 0.057700]
    test_list.append(turb_flatplate)

    # Flat plate (compressible) with species inlet
    turb_flatplate_species           = TestCase('turb_flatplate_species')
    turb_flatplate_species.cfg_dir   = "rans/flatplate"
    turb_flatplate_species.cfg_file  = "turb_SA_flatplate_species.cfg"
    turb_flatplate_species.test_iter = 20
    turb_flatplate_species.test_vals = [-4.249479, -0.634915, -1.716291, 1.223196, -3.307930, 9.000000, -6.634088, 5.000000, -6.985954, 10.000000, -6.255640, 0.999903, 0.999903]
    test_list.append(turb_flatplate_species)

    # Flat plate SST compressibility correction Wilcox
    turb_flatplate_CC_Wilcox = TestCase('turb_flatplate_CC_Wilcox')
    turb_flatplate_CC_Wilcox.cfg_dir   = "rans/flatplate"
    turb_flatplate_CC_Wilcox.cfg_file  = "turb_SST_flatplate_compressibility_Wilcox.cfg"
    turb_flatplate_CC_Wilcox.test_iter = 20
    turb_flatplate_CC_Wilcox.test_vals = [-1.280934, 1.974155, 1.440363, 5.038339, -3.791118, 8.297543]
    test_list.append(turb_flatplate_CC_Wilcox)

    # Flat plate SST compressibility correction Sarkar
    turb_flatplate_CC_Sarkar = TestCase('turb_flatplate_CC_Sarkar')
    turb_flatplate_CC_Sarkar.cfg_dir   = "rans/flatplate"
    turb_flatplate_CC_Sarkar.cfg_file  = "turb_SST_flatplate_compressibility_Sarkar.cfg"
    turb_flatplate_CC_Sarkar.test_iter = 20
    turb_flatplate_CC_Sarkar.test_vals = [-1.280934, 1.974155, 1.440363, 5.038339, -3.791121, 8.297543]
    test_list.append(turb_flatplate_CC_Sarkar)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.408685, -6.662907, 0.238579, 0.158968, 0.000000]
    turb_oneram6.timeout   = 3200
    test_list.append(turb_oneram6)

    # ONERA M6 Wing - vorticity confinement
    turb_oneram6_vc = TestCase('turb_oneram6_vc')
    turb_oneram6_vc.cfg_dir = "rans/oneram6"
    turb_oneram6_vc.cfg_file = "turb_ONERAM6_vc.cfg"
    turb_oneram6_vc.test_iter = 15
    turb_oneram6_vc.test_vals = [-2.282318, -6.614780, 0.234330, 0.143024, 0.000000]
    turb_oneram6_vc.timeout = 3200
    test_list.append(turb_oneram6_vc)

    # ONERA M6 Wing - Newton-Krylov
    turb_oneram6_nk           = TestCase('turb_oneram6_nk')
    turb_oneram6_nk.cfg_dir   = "rans/oneram6"
    turb_oneram6_nk.cfg_file  = "turb_ONERAM6_nk.cfg"
    turb_oneram6_nk.test_iter = 20
    turb_oneram6_nk.test_vals = [-4.828066, -4.426250, -11.417591, 0.224679, 0.044309, 1.000000, -0.642977, 31.384000]
    turb_oneram6_nk.timeout   = 600
    turb_oneram6_nk.tol       = 0.0001
    test_list.append(turb_oneram6_nk)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 5
    turb_naca0012_sa.test_vals = [-12.094695, -14.685268, 1.057665, 0.022971, 20.000000, -1.692967, 20.000000, -4.037673, 0]
    turb_naca0012_sa.timeout   = 3200
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-12.232530, -14.435294, -6.037147, 1.047444, 0.019214, -1.704457, 0.000000]
    turb_naca0012_sst.test_vals_aarch64 = [-12.232530, -14.435294, -6.037147, 1.047444, 0.019214, -1.704457, 0.000000]
    turb_naca0012_sst.timeout   = 3200
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals = [-12.146171, -14.781592, -6.358755, 1.000270, 0.019123, -1.593320]
    turb_naca0012_sst_sust.test_vals_aarch64 = [-12.146171, -14.781592, -6.358755, 1.000270, 0.019123, -1.593320]
    turb_naca0012_sst_sust.timeout   = 3200
    test_list.append(turb_naca0012_sst_sust)

    # NACA0012 (SST, 2003m, Vorticity)
    turb_naca0012_sst_2003_Vm           = TestCase('turb_naca0012_sst_2003_Vm')
    turb_naca0012_sst_2003_Vm.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_2003_Vm.cfg_file  = "turb_NACA0012_sst_2003-Vm.cfg"
    turb_naca0012_sst_2003_Vm.test_iter = 10
    turb_naca0012_sst_2003_Vm.test_vals = [-8.263591, -10.306174, -3.626544, 1.045289, 0.019325, -1.551890]
    turb_naca0012_sst_2003_Vm.timeout   = 3200
    test_list.append(turb_naca0012_sst_2003_Vm)

    # NACA0012 (SST, 1994m Kato-Launder)
    turb_naca0012_sst_1994_KLm           = TestCase('turb_naca0012_sst_1994_KLm')
    turb_naca0012_sst_1994_KLm.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_1994_KLm.cfg_file  = "turb_NACA0012_sst_1994-KLm.cfg"
    turb_naca0012_sst_1994_KLm.test_iter = 10
    turb_naca0012_sst_1994_KLm.test_vals = [-8.563675, -10.803009, -4.036981, 1.046738, 0.019251, -1.806475]
    turb_naca0012_sst_1994_KLm.timeout   = 3200
    test_list.append(turb_naca0012_sst_1994_KLm)


    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [-5.216568, -10.236279, -1.556755, 1.022392, 0.040547, -3.738510]
    turb_naca0012_sst_fixedvalues.timeout   = 3200
    test_list.append(turb_naca0012_sst_fixedvalues)

    # NACA0012 (SST, explicit Euler for flow and turbulence equations)
    turb_naca0012_sst_expliciteuler           = TestCase('turb_naca0012_sst_expliciteuler')
    turb_naca0012_sst_expliciteuler.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_expliciteuler.cfg_file  = "turb_NACA0012_sst_expliciteuler.cfg"
    turb_naca0012_sst_expliciteuler.test_iter = 10
    turb_naca0012_sst_expliciteuler.test_vals = [-3.533765, -3.157766, 3.364026, 1.124757, 0.501700, -float("inf")]
    turb_naca0012_sst_expliciteuler.timeout   = 3200
    test_list.append(turb_naca0012_sst_expliciteuler)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389724, -8.409502, 0.000048, 0.056344]
    propeller.timeout   = 3200
    test_list.append(propeller)

    # Actuator disk BEM method for propeller
    actuatordisk_bem = TestCase('actuatordisk_bem')
    actuatordisk_bem.cfg_dir = "rans/actuatordisk_bem"
    actuatordisk_bem.cfg_file = "actuatordisk_bem.cfg"
    actuatordisk_bem.test_iter = 15
    actuatordisk_bem.test_vals = [-5.388943, -10.318617, 0.001362, -0.376520]
    actuatordisk_bem.timeout = 3200
    actuatordisk_bem.tol = 0.001
    test_list.append(actuatordisk_bem)

    #######################################
    ### Axisymmetric Compressible RANS  ###
    #######################################

    # Axisymmetric air nozzle (transonic) restart
    axi_rans_air_nozzle_restart           = TestCase('axi_rans_air_nozzle_restart')
    axi_rans_air_nozzle_restart.cfg_dir   = "axisymmetric_rans/air_nozzle"
    axi_rans_air_nozzle_restart.cfg_file  = "air_nozzle_restart.cfg"
    axi_rans_air_nozzle_restart.test_iter = 10
    axi_rans_air_nozzle_restart.test_vals = [-12.067375, -6.838713, -8.845034, -3.786225, 0.000000]
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
    turb_naca0012_sst_restart_mg.test_vals = [-7.600526, -7.181924, -0.627078, -0.000016, 0.078729]
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
    inc_euler_naca0012.test_vals = [-7.154146, -6.507095, 0.532001, 0.008466]
    test_list.append(inc_euler_naca0012)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals = [-6.576818, -5.796867, -0.003141, 0.126481]
    test_list.append(inc_nozzle)

    # Laminar wall mounted cylinder, Euler walls, cylinder wall diagonally split
    inc_cylinder_split = TestCase('inc_cylinder_split')
    inc_cylinder_split.cfg_dir = "incomp_navierstokes/cylinder_split"
    inc_cylinder_split.cfg_file = "cylinder_split.cfg"
    inc_cylinder_split.test_iter = 10
    inc_cylinder_split.test_vals = [-10.207060, -11.197502, -11.296847, 10.000000]
    test_list.append(inc_cylinder_split)

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

    # Laminar sphere, Re=1. Last column: Cd=24/Re
    inc_lam_sphere          = TestCase('inc_lam_sphere')
    inc_lam_sphere.cfg_dir   = "incomp_navierstokes/sphere"
    inc_lam_sphere.cfg_file  = "sphere.cfg"
    inc_lam_sphere.test_iter = 5
    inc_lam_sphere.test_vals = [-8.342048, -9.328063, 0.121003, 25.782687]
    test_list.append(inc_lam_sphere)

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
    inc_poly_cylinder.test_vals = [-7.786564, -2.036735, 0.012980, 1.944887, -170.930000]
    test_list.append(inc_poly_cylinder)

    # X-coarse laminar bend as a mixed element CGNS test
    inc_lam_bend            = TestCase('inc_lam_bend')
    inc_lam_bend.cfg_dir    = "incomp_navierstokes/bend"
    inc_lam_bend.cfg_file   = "lam_bend.cfg"
    inc_lam_bend.test_iter = 10
    inc_lam_bend.test_vals = [-3.588863, -3.101606, -0.022302, 1.062971]
    test_list.append(inc_lam_bend)

    # 3D laminar channnel with 1 cell in flow direction, streamwise periodic
    sp_pipeSlice_3d_dp_hf_tp           = TestCase('sp_pipeSlice_3d_dp_hf_tp')
    sp_pipeSlice_3d_dp_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/pipeSlice_3d"
    sp_pipeSlice_3d_dp_hf_tp.cfg_file  = "sp_pipeSlice_3d_dp_hf_tp.cfg"
    sp_pipeSlice_3d_dp_hf_tp.test_iter = 10
    sp_pipeSlice_3d_dp_hf_tp.test_vals = [-11.119796, -11.234737, -9.429479, -0.000023]
    test_list.append(sp_pipeSlice_3d_dp_hf_tp)

    # 2D pin array with heat transfer BC on pin surfaces
    inc_heatTransfer_BC           = TestCase('inc_heatTransfer_BC')
    inc_heatTransfer_BC.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    inc_heatTransfer_BC.cfg_file  = "BC_HeatTransfer.cfg"
    inc_heatTransfer_BC.test_iter = 50
    inc_heatTransfer_BC.test_vals = [-8.747445, -7.653823, -8.079512, -0.654920, -1670.100000]
    test_list.append(inc_heatTransfer_BC)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.788595, -11.040625, -0.000002, 0.309519]
    test_list.append(inc_turb_naca0012)

    # NACA0012, SST_SUST
    inc_turb_naca0012_sst_sust           = TestCase('inc_turb_naca0012_sst_sust')
    inc_turb_naca0012_sst_sust.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_sst_sust.cfg_file  = "naca0012_SST_SUST.cfg"
    inc_turb_naca0012_sst_sust.test_iter = 20
    inc_turb_naca0012_sst_sust.test_vals = [-7.291511, 0.132644, -0.000001, 0.312109]
    test_list.append(inc_turb_naca0012_sst_sust)

    ####################
    ### DG-FEM Euler ###
    ####################

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
    test_list.append(fem_ns_sphere_ader)

    # Unsteady cylinder
    fem_ns_unsteady_cylinder           = TestCase('fem_ns_unsteady_cylinder')
    fem_ns_unsteady_cylinder.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder.cfg_file  = "fem_unst_cylinder.cfg"
    fem_ns_unsteady_cylinder.test_iter = 11
    fem_ns_unsteady_cylinder.test_vals = [-3.558582, -3.014464, -0.038927, 1.383983]
    fem_ns_unsteady_cylinder.unsteady  = True
    test_list.append(fem_ns_unsteady_cylinder)

    # Unsteady cylinder ADER
    fem_ns_unsteady_cylinder_ader           = TestCase('fem_ns_unsteady_cylinder_ader')
    fem_ns_unsteady_cylinder_ader.cfg_dir   = "hom_navierstokes/UnsteadyCylinder/nPoly4"
    fem_ns_unsteady_cylinder_ader.cfg_file  = "fem_unst_cylinder_ADER.cfg"
    fem_ns_unsteady_cylinder_ader.test_iter = 11
    fem_ns_unsteady_cylinder_ader.test_vals = [-35.000000, -35.000000, -0.041003, 1.391339]
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
    turbmod_sa_bsl_rae2822.test_vals = [-2.004689, 0.742306, 0.497308, -5.265797, 0.809463, 0.062016]
    test_list.append(turbmod_sa_bsl_rae2822)

    # SA Negative
    turbmod_sa_neg_rae2822           = TestCase('turbmod_sa_neg_rae2822')
    turbmod_sa_neg_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_neg_rae2822.cfg_file  = "turb_SA_NEG_RAE2822.cfg"
    turbmod_sa_neg_rae2822.test_iter = 10
    turbmod_sa_neg_rae2822.test_vals         = [-1.204800, 1.611685, 1.349330, 1.489602, 1.263603, 0.466487, 0]
    turbmod_sa_neg_rae2822.test_vals_aarch64 = [-1.359612, 1.493629, 1.218367, -1.441703, 1.248499, 0.457987, 0]
    test_list.append(turbmod_sa_neg_rae2822)

    # SA Compressibility Correction
    turbmod_sa_comp_rae2822           = TestCase('turbmod_sa_comp_rae2822')
    turbmod_sa_comp_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_comp_rae2822.cfg_file  = "turb_SA_COMP_RAE2822.cfg"
    turbmod_sa_comp_rae2822.test_iter = 20
    turbmod_sa_comp_rae2822.test_vals = [-2.004687, 0.742304, 0.497309, -5.266084, 0.809467, 0.062029]
    test_list.append(turbmod_sa_comp_rae2822)

    # SA Edwards
    turbmod_sa_edw_rae2822           = TestCase('turbmod_sa_edw_rae2822')
    turbmod_sa_edw_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_edw_rae2822.cfg_file  = "turb_SA_EDW_RAE2822.cfg"
    turbmod_sa_edw_rae2822.test_iter = 20
    turbmod_sa_edw_rae2822.test_vals = [-2.004687, 0.742306, 0.497310, -5.290769, 0.809485, 0.062036]
    test_list.append(turbmod_sa_edw_rae2822)

    # SA Compressibility and Edwards
    turbmod_sa_comp_edw_rae2822           = TestCase('turbmod_sa_comp_edw_rae2822')
    turbmod_sa_comp_edw_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_comp_edw_rae2822.cfg_file  = "turb_SA_COMP_EDW_RAE2822.cfg"
    turbmod_sa_comp_edw_rae2822.test_iter = 20
    turbmod_sa_comp_edw_rae2822.test_vals = [-2.004685, 0.742307, 0.497311, -5.290750, 0.809487, 0.062045]
    test_list.append(turbmod_sa_comp_edw_rae2822)

    # SA QCR
    turbmod_sa_qcr_rae2822           = TestCase('turbmod_sa_qcr_rae2822')
    turbmod_sa_qcr_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_qcr_rae2822.cfg_file  = "turb_SA_QCR_RAE2822.cfg"
    turbmod_sa_qcr_rae2822.test_iter = 20
    turbmod_sa_qcr_rae2822.test_vals = [-2.004793, 0.742353, 0.497315, -5.265977, 0.807841, 0.062027]
    test_list.append(turbmod_sa_qcr_rae2822)

    ############################
    ###      Transition      ###
    ############################

    # Schubauer-Klebanoff Natural Transition Case
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.test_vals    = [-8.215651, -13.240283, 0.000048, 0.007983]
    test_list.append(schubauer_klebanoff_transition)

    #####################################
    ### Cont. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    contadj_naca0012           = TestCase('contadj_naca0012')
    contadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012.test_iter = 5
    contadj_naca0012.test_vals = [-9.662585, -14.998832, -0.726250, 0.020280]
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.032190, -12.587083, -1.086100, 0.007556]
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10
    contadj_wedge.test_vals = [2.872064, -2.756210, 1010800.000000, 0.000000]
    test_list.append(contadj_wedge)

    # Inviscid fixed CL NACA0012
    contadj_fixed_CL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixed_CL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixed_CL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixed_CL_naca0012.test_iter = 100
    contadj_fixed_CL_naca0012.test_vals = [0.748407, -4.810872, -0.520110, -0.000291]
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
    turb_naca0012_1c.test_vals = [-4.981445, 1.140044, 0.466241, -0.076886]
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals = [-5.485045, 0.968951, 0.317727, -0.109873]
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals = [-5.584375, 0.931875, 0.274184, -0.113514]
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals = [-5.122421, 1.074501, 0.412268, -0.091538]
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals = [-5.549646, 0.946077, 0.284451, -0.112299]
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
    hb_rans_preconditioning.tol       = 0.00001
    hb_rans_preconditioning.test_vals = [-1.902098, 0.484070, 0.601481, 3.609002, -5.949356]
    test_list.append(hb_rans_preconditioning)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Rotating NACA 0012
    rot_naca0012           = TestCase('rot_naca0012')
    rot_naca0012.cfg_dir   = "rotating/naca0012"
    rot_naca0012.cfg_file  = "rot_NACA0012.cfg"
    rot_naca0012.test_iter = 25
    rot_naca0012.test_vals = [-2.603551, 2.924633, -0.081272, 0.002162]
    test_list.append(rot_naca0012)

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.610923, -0.146741, 1.115860, 1.490430]
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.806056, -2.364884, 1.685228, 1.518276]
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-1.175927, 0.062080, 1.399402, 2.220364, 1.399350, 2.218602, 0.000000]
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977498, 3.481817, -0.010773, -0.008068]
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Aeroelastic
    aeroelastic           = TestCase('aeroelastic')
    aeroelastic.cfg_dir   = "aeroelastic"
    aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    aeroelastic.test_iter = 2
    aeroelastic.test_vals = [0.075023, 0.027483, -0.001643, -0.000126]
    aeroelastic.unsteady  = True
    test_list.append(aeroelastic)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714713, -5.788290, -0.214960, 0.023758, 0.000000]
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.004011, -6.876250, 1.487888, 0.421869]
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # Flat plate
    flatplate_unsteady           = TestCase('flatplate_unsteady')
    flatplate_unsteady.cfg_dir   = "navierstokes/flatplate"
    flatplate_unsteady.cfg_file  = "lam_flatplate_unst.cfg"
    flatplate_unsteady.test_iter = 3
    flatplate_unsteady.test_vals = [-8.875128, -8.250204, -6.305788, -5.469452, -3.398230, 0.002075, -0.325535]
    flatplate_unsteady.unsteady  = True
    test_list.append(flatplate_unsteady)

    ######################################
    ### NICFD                          ###
    ######################################

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 50
    edge_VW.test_vals = [-9.057409, -2.833203, -0.000009, 0.000000]
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 50
    edge_PPR.test_vals = [-9.781896, -3.630892, -0.000034, 0.000000]
    test_list.append(edge_PPR)

    # Rarefaction Q1D nozzle, include CoolProp fluid model
    coolprop_fluidModel           = TestCase('coolprop_fluidModel')
    coolprop_fluidModel.cfg_dir   = "nicf/coolprop"
    coolprop_fluidModel.cfg_file  = "fluidModel.cfg"
    coolprop_fluidModel.test_iter = 5
    coolprop_fluidModel.test_vals = [-4.684483, -1.583073, 3.724768, 0.000000, 0.000000]
    coolprop_fluidModel.enabled_on_cpu_arch = ["x86_64"]
    test_list.append(coolprop_fluidModel)

    # Rarefaction Q1D nozzle, include CoolProp transport model
    coolprop_transportModel           = TestCase('coolprop_transportModel')
    coolprop_transportModel.cfg_dir   = "nicf/coolprop"
    coolprop_transportModel.cfg_file  = "transportModel.cfg"
    coolprop_transportModel.test_iter = 5
    coolprop_transportModel.test_vals = [-4.684988, -1.314608, 4.668105, 0.000000, 0.000000]
    coolprop_transportModel.enabled_on_cpu_arch = ["x86_64"]
    test_list.append(coolprop_transportModel)

    # Rarefaction Q1D nozzle, include data-driven fluid model
    datadriven_fluidModel           = TestCase('datadriven_fluidModel')
    datadriven_fluidModel.cfg_dir   = "nicf/datadriven"
    datadriven_fluidModel.cfg_file  = "datadriven_nozzle.cfg"
    datadriven_fluidModel.test_iter = 50
    datadriven_fluidModel.test_vals = [-6.338898, -3.837472, -4.351292, -1.860262, -2.700991, 0.691801]
    test_list.append(datadriven_fluidModel)

    ######################################
    ### Turbomachinery                 ###
    ######################################

    # Aachen Turbine restart
    Aachen_3D_restart = TestCase('aachen_turbine_restart')
    Aachen_3D_restart.cfg_dir = "turbomachinery/Aachen_turbine"
    Aachen_3D_restart.cfg_file = "aachen_3D_MP_restart.cfg"
    Aachen_3D_restart.test_iter = 5
    Aachen_3D_restart.test_vals = [-7.701448, -8.512241, -6.014939, -6.468738, -5.801759, -4.607179, -5.551037, -5.300771, -3.804188, -5.256055, -5.765160, -3.609605, -2.229276, -2.883962, -0.563469]
    test_list.append(Aachen_3D_restart)

    # Jones APU Turbocharger restart
    Jones_tc_restart           = TestCase('jones_turbocharger_restart')
    Jones_tc_restart.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_restart.cfg_file  = "Jones_restart.cfg"
    Jones_tc_restart.test_iter = 5
    Jones_tc_restart.test_vals = [-7.308012, -5.332077, -14.895814, -9.330692, -12.071730, -6.548620, 73291.000000, 73291.000000, 0.020111, 82.896000]
    test_list.append(Jones_tc_restart)

    # 2D axial stage
    axial_stage2D           = TestCase('axial_stage2D')
    axial_stage2D.cfg_dir   = "turbomachinery/axial_stage_2D"
    axial_stage2D.cfg_file  = "Axial_stage2D.cfg"
    axial_stage2D.test_iter = 20
    axial_stage2D.test_vals = [1.090054, 1.550990, -2.895064, 2.607596, -2.479704, 3.063740, 106380.000000, 106380.000000, 5.733600, 64.747000]
    test_list.append(axial_stage2D)

    # 2D transonic stator restart
    transonic_stator_restart           = TestCase('transonic_stator_restart')
    transonic_stator_restart.cfg_dir   = "turbomachinery/transonic_stator_2D"
    transonic_stator_restart.cfg_file  = "transonic_stator_restart.cfg"
    transonic_stator_restart.test_iter = 20
    transonic_stator_restart.test_vals = [-4.354756, -2.473198, -2.076102, 1.740556, -1.441690, 3.246152, -471620.000000, 94.839000, -0.050276]
    transonic_stator_restart.test_vals_aarch64 = [-5.011834, -3.091110, -2.757795, 1.087934, -3.544707, 2.166101, -471630, 94.868, -0.035888]
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
    uniform_flow.test_vals = [5.000000, 0.000000, -0.195002, -10.624448]
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow)

    # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 2
    channel_2D.test_vals = [2.000000, 0.000000, 0.464931, 0.348057, 0.397535]
    channel_2D.timeout   = 100
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 2
    channel_3D.test_vals = [2.000000, 0.000000, 0.629113, 0.524906, 0.422425]
    channel_3D.unsteady  = True
    channel_3D.multizone = True
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [0.080827, 0.547324, 0.655095, 0.968235, 1.049121]
    pipe.unsteady  = True
    pipe.multizone = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [3.000000, 0.000000, 0.717065, 1.119815, 1.160330]
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone  = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 1.207118, 1.065260]
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone  = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [13.000000, 1.167840, -1.660900]
    bars_SST_2D.multizone = True
    test_list.append(bars_SST_2D)

    # Sliding mesh with incompressible flows (steady)
    slinc_steady           = TestCase('slinc_steady')
    slinc_steady.cfg_dir   = "sliding_interface/incompressible_steady"
    slinc_steady.cfg_file  = "config.cfg"
    slinc_steady.test_iter = 19
    slinc_steady.test_vals = [19.000000, -1.778863, -2.040209]
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
    statbeam3d.test_vals         = [-6.058758, -5.750933, -5.892188, 110190]
    statbeam3d.test_vals_aarch64 = [-6.058758, -5.750933, -5.892188, 110190]
    statbeam3d.command   = TestCase.Command(exec = "parallel_computation_fsi.py", param = "-f")
    test_list.append(statbeam3d)

    # Static beam, 3d with coupled temperature
    thermal_beam_3d = TestCase('thermal_beam_3d')
    thermal_beam_3d.cfg_dir = "fea_fsi/ThermalBeam_3d"
    thermal_beam_3d.cfg_file = "configBeam_3d.cfg"
    thermal_beam_3d.test_iter = 0
    thermal_beam_3d.test_vals = [-6.140220, -5.842734, -5.972391, -8.091358, 262, -8.246755, 81, -8.298569, 135620, 144.65]
    thermal_beam_3d.command = TestCase.Command(exec = "parallel_computation_fsi.py", param = "-f")
    test_list.append(thermal_beam_3d)

    # Rotating cylinder, 3d
    rotating_cylinder_fea           = TestCase('rotating_cylinder_fea')
    rotating_cylinder_fea.cfg_dir   = "fea_fsi/rotating_cylinder"
    rotating_cylinder_fea.cfg_file  = "config.cfg"
    rotating_cylinder_fea.test_iter = 0
    # For a thin disk with the inner and outer radius of this geometry, from
    # "Formulas for Stress, Strain, and Structural Matrices", 2nd Edition, figure 19-4,
    # the maximum stress is 165.6MPa, we get a Von Misses stress very close to that.
    rotating_cylinder_fea.test_vals = [-6.861940, -6.835550, -6.895497, 22, -8.313847, 1.6502e+08]
    test_list.append(rotating_cylinder_fea)

    # Dynamic beam, 2d
    dynbeam2d           = TestCase('dynbeam2d')
    dynbeam2d.cfg_dir   = "fea_fsi/DynBeam_2d"
    dynbeam2d.cfg_file  = "configBeam_2d.cfg"
    dynbeam2d.test_iter = 6
    dynbeam2d.test_vals = [-3.240015, 2.895057, -0.353146, 66127.000000]
    dynbeam2d.unsteady  = True
    test_list.append(dynbeam2d)

    # FSI, 2d
    fsi2d           = TestCase('fsi2d')
    fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    fsi2d.cfg_file  = "configFSI.cfg"
    fsi2d.test_iter = 4
    fsi2d.test_vals = [4.000000, 0.000000, -3.726013, -4.277768]
    fsi2d.command   = TestCase.Command(exec = "parallel_computation_fsi.py", param = "-f")
    fsi2d.multizone= True
    fsi2d.unsteady = True
    test_list.append(fsi2d)

    # FSI, Static, 2D, new mesh solver
    stat_fsi           = TestCase('stat_fsi')
    stat_fsi.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi.cfg_file  = "config.cfg"
    stat_fsi.test_iter = 7
    stat_fsi.test_vals = [-3.301938, -4.971986, 0.000000, 11.000000]
    stat_fsi.multizone = True
    test_list.append(stat_fsi)

    # FSI, Dynamic, 2D, new mesh solver
    dyn_fsi           = TestCase('dyn_fsi')
    dyn_fsi.cfg_dir   = "fea_fsi/dyn_fsi"
    dyn_fsi.cfg_file  = "config.cfg"
    dyn_fsi.test_iter = 4
    dyn_fsi.test_vals = [-4.330741, -4.153001, 0.000000, 97.000000]
    dyn_fsi.multizone = True
    dyn_fsi.unsteady  = True
    test_list.append(dyn_fsi)

    # FSI, Static, 2D, new mesh solver, restart
    stat_fsi_restart           = TestCase('stat_fsi_restart')
    stat_fsi_restart.cfg_dir   = "fea_fsi/stat_fsi"
    stat_fsi_restart.cfg_file  = "config_restart.cfg"
    stat_fsi_restart.test_iter = 1
    stat_fsi_restart.test_vals = [-3.486655, -4.425104, 0.000000, 27.000000]
    stat_fsi_restart.multizone = True
    test_list.append(stat_fsi_restart)

    # ###############################
    # ### Radiative Heat Transfer ###
    # ###############################

    # Radiative heat transfer
    p1rad           = TestCase('p1rad')
    p1rad.cfg_dir   = "radiation/p1model"
    p1rad.cfg_file  = "configp1.cfg"
    p1rad.test_iter = 100
    p1rad.test_vals = [-7.743666, -7.921411, -2.111848, 0.098302, -47.897000]
    test_list.append(p1rad)


    # #############################
    # ### Solid Heat Conduction ###
    # #############################

    # 2D pins, periodically connected
    solid_periodic_pins           = TestCase('solid_periodic_pins')
    solid_periodic_pins.cfg_dir   = "solid_heat_conduction/periodic_pins"
    solid_periodic_pins.cfg_file  = "configSolid.cfg"
    solid_periodic_pins.test_iter = 750
    solid_periodic_pins.test_vals = [-15.878977, -14.569206, 300.900000, 425.320000, 5.000000, -1.672737]
    solid_periodic_pins.test_vals_aarch64 = [-15.879016, -14.569206, 300.900000, 425.320000, 5.000000, -1.672666]
    test_list.append(solid_periodic_pins)

    # ###############################
    # ### Conjugate heat transfer ###
    # ###############################

    # CHT incompressible
    cht_incompressible           = TestCase('cht_incompressible')
    cht_incompressible.cfg_dir   = "coupled_cht/incomp_2d"
    cht_incompressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_incompressible.test_iter = 10
    cht_incompressible.test_vals = [-2.128827, -0.588812, -0.588812, -0.588812]
    cht_incompressible.multizone = True
    test_list.append(cht_incompressible)

    # CHT compressible
    cht_compressible           = TestCase('cht_compressible')
    cht_compressible.cfg_dir   = "coupled_cht/comp_2d"
    cht_compressible.cfg_file  = "cht_2d_3cylinders.cfg"
    cht_compressible.test_iter = 10
    cht_compressible.test_vals = [-4.256053, -0.532725, -0.532725, -0.532726]
    cht_compressible.multizone = True
    test_list.append(cht_compressible)

    # 2D CHT case streamwise periodicity. Also test Multizone PerSurface screen output.
    sp_pinArray_cht_2d_dp_hf           = TestCase('sp_pinArray_cht_2d_dp_hf')
    sp_pinArray_cht_2d_dp_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    sp_pinArray_cht_2d_dp_hf.cfg_file  = "configMaster.cfg"
    sp_pinArray_cht_2d_dp_hf.test_iter = 100
    sp_pinArray_cht_2d_dp_hf.test_vals = [0.500399, -0.667466, -0.984103, -0.726712, 208.023676, 350.140000, -0.000000, -0.726710, 0.726710]
    sp_pinArray_cht_2d_dp_hf.multizone = True
    test_list.append(sp_pinArray_cht_2d_dp_hf)

    # simple small 3D pin case massflow periodic with heatflux BC
    sp_pinArray_3d_cht_mf_hf_tp           = TestCase('sp_pinArray_3d_cht_mf_hf_tp')
    sp_pinArray_3d_cht_mf_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_3d"
    sp_pinArray_3d_cht_mf_hf_tp.cfg_file  = "configMaster.cfg"
    sp_pinArray_3d_cht_mf_hf_tp.test_iter = 30
    sp_pinArray_3d_cht_mf_hf_tp.test_vals         = [-1.116189, -1.881686, -2.561500, -0.009805, 104.600449, 418.360000, 0.000000]
    sp_pinArray_3d_cht_mf_hf_tp.test_vals_aarch64 = [-1.117102, -1.880628, -2.561816, -0.009804, 104.600540, 418.360000, 0.000000]
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
    pywrapper_naca0012.test_vals = [-9.249009, -8.546597, 0.335769, 0.023275]
    pywrapper_naca0012.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--parallel -f")
    test_list.append(pywrapper_naca0012)

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    pywrapper_turb_naca0012_sst           = TestCase('pywrapper_turb_naca0012_sst')
    pywrapper_turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    pywrapper_turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    pywrapper_turb_naca0012_sst.test_iter = 10
    pywrapper_turb_naca0012_sst.test_vals = [-12.232530, -14.435294, -6.037147, 1.047444, 0.019214, -1.704457, 0.000000]
    pywrapper_turb_naca0012_sst.test_vals_aarch64 = [-12.232530, -14.435294, -6.037147, 1.047444, 0.019214, -1.704457, 0.000000]
    pywrapper_turb_naca0012_sst.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--parallel -f")
    pywrapper_turb_naca0012_sst.timeout   = 3200
    test_list.append(pywrapper_turb_naca0012_sst)

    # Square cylinder
    pywrapper_square_cylinder           = TestCase('pywrapper_square_cylinder')
    pywrapper_square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    pywrapper_square_cylinder.cfg_file  = "turb_square.cfg"
    pywrapper_square_cylinder.test_iter = 10
    pywrapper_square_cylinder.test_vals = [-1.178764, -0.349891, 1.401055, 2.358075, 1.401419, 2.301174, 0.000000]
    pywrapper_square_cylinder.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--parallel -f")
    pywrapper_square_cylinder.unsteady  = True
    test_list.append(pywrapper_square_cylinder)

    # Aeroelastic
    pywrapper_aeroelastic         = TestCase('pywrapper_aeroelastic')
    pywrapper_aeroelastic.cfg_dir   = "aeroelastic"
    pywrapper_aeroelastic.cfg_file  = "aeroelastic_NACA64A010.cfg"
    pywrapper_aeroelastic.test_iter = 2
    pywrapper_aeroelastic.test_vals = [0.075023, 0.027483, -0.001643, -0.000126]
    pywrapper_aeroelastic.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--parallel -f")
    pywrapper_aeroelastic.unsteady  = True
    test_list.append(pywrapper_aeroelastic)

    # Custom FEA load
    pywrapper_custom_fea_load = TestCase('pywrapper_custom_fea_load')
    pywrapper_custom_fea_load.cfg_dir = "py_wrapper/custom_load_fea"
    pywrapper_custom_fea_load.cfg_file = "config.cfg"
    pywrapper_custom_fea_load.test_iter = 13
    pywrapper_custom_fea_load.test_vals = [-7.263559, -4.946814, -14.165142, 34.000000, -6.380144, 320.580000]
    pywrapper_custom_fea_load.test_vals_aarch64 = [-7.263558, -4.946814, -14.165142, 35.000000, -6.802790, 320.580000]
    pywrapper_custom_fea_load.command = TestCase.Command("mpirun -np 2", "python", "run.py")
    test_list.append(pywrapper_custom_fea_load)

    # FSI, 2d
    pywrapper_fsi2d           = TestCase('pywrapper_fsi2d')
    pywrapper_fsi2d.cfg_dir   = "fea_fsi/WallChannel_2d"
    pywrapper_fsi2d.cfg_file  = "configFSI.cfg"
    pywrapper_fsi2d.test_iter = 4
    pywrapper_fsi2d.test_vals = [4.000000, 0.000000, -3.726013, -4.277768]
    pywrapper_fsi2d.command   = TestCase.Command("mpirun -np 2", "SU2_CFD.py", "--nZone 2 --fsi True --parallel -f")
    pywrapper_fsi2d.unsteady  = True
    pywrapper_fsi2d.multizone = True
    test_list.append(pywrapper_fsi2d)

    # Unsteady FSI with custom load
    pywrapper_unsteadyFSI = TestCase('pywrapper_unsteadyFSI')
    pywrapper_unsteadyFSI.cfg_dir = "py_wrapper/dyn_fsi"
    pywrapper_unsteadyFSI.cfg_file = "config.cfg"
    pywrapper_unsteadyFSI.test_iter = 4
    pywrapper_unsteadyFSI.test_vals = [0.000000, 31.000000, 5.000000, 58.000000, -1.756677, -2.828286, -7.638545, -6.863959, 0.000156]
    pywrapper_unsteadyFSI.command = TestCase.Command("mpirun -np 2", "python", "run.py")
    pywrapper_unsteadyFSI.unsteady = True
    pywrapper_unsteadyFSI.multizone = True
    test_list.append(pywrapper_unsteadyFSI)

    # Unsteady CHT
    pywrapper_unsteadyCHT = TestCase('pywrapper_unsteadyCHT')
    pywrapper_unsteadyCHT.cfg_dir = "py_wrapper/flatPlate_unsteady_CHT"
    pywrapper_unsteadyCHT.cfg_file = "unsteady_CHT_FlatPlate_Conf.cfg"
    pywrapper_unsteadyCHT.test_iter = 5
    pywrapper_unsteadyCHT.test_vals = [-1.614167, 2.247364, -0.001388, 0.172995]
    pywrapper_unsteadyCHT.command = TestCase.Command("mpirun -np 2", "python", "launch_unsteady_CHT_FlatPlate.py --parallel -f")
    pywrapper_unsteadyCHT.unsteady = True
    test_list.append(pywrapper_unsteadyCHT)

    # Rigid motion
    pywrapper_rigidMotion = TestCase('pywrapper_rigidMotion')
    pywrapper_rigidMotion.cfg_dir = "py_wrapper/flatPlate_rigidMotion"
    pywrapper_rigidMotion.cfg_file = "flatPlate_rigidMotion_Conf.cfg"
    pywrapper_rigidMotion.test_iter = 5
    pywrapper_rigidMotion.test_vals = [-1.614166, 2.243100, 0.350194, 0.089497]
    pywrapper_rigidMotion.command = TestCase.Command("mpirun -np 2", "python", "launch_flatPlate_rigidMotion.py --parallel -f")
    pywrapper_rigidMotion.unsteady = True
    test_list.append(pywrapper_rigidMotion)

    # Deforming Bump in Channel
    pywrapper_deformingBump = TestCase('pywrapper_deformingBump')
    pywrapper_deformingBump.cfg_dir = "py_wrapper/deforming_bump_in_channel"
    pywrapper_deformingBump.cfg_file = "config.cfg"
    pywrapper_deformingBump.test_iter = 1
    pywrapper_deformingBump.test_vals = [0.500000, 0.000000, -2.811520, -1.603562, -2.074259, 2.424289, 7.616891, -0.205655]
    pywrapper_deformingBump.command = TestCase.Command("mpirun -np 2", "python", "run.py")
    pywrapper_deformingBump.unsteady = True
    test_list.append(pywrapper_deformingBump)

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
    mms_fvm_inc_euler.test_vals = [-9.128660, -9.441806, 0.000000, 0.000000]
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

    #####################
    ## Species solver ###
    #####################

    # 2 species (1 eq) primitive venturi mixing using mixing model
    species2_primitiveVenturi_mixingmodel           = TestCase('species2_primitiveVenturi_mixingmodel')
    species2_primitiveVenturi_mixingmodel.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel.cfg_file  = "species2_primitiveVenturi_mixingmodel.cfg"
    species2_primitiveVenturi_mixingmodel.test_iter = 50
    species2_primitiveVenturi_mixingmodel.test_vals = [-5.736233, -4.561289, -4.666830, -5.863760, -0.071055, -5.584743, 5.000000, -1.376449, 5.000000, -4.869012, 5.000000, -1.452251, 0.000372, 0.000356, 0.000016, 0.000000]
    test_list.append(species2_primitiveVenturi_mixingmodel)

    # 2 species (1 eq) primitive venturi mixing using mixing model and bounded scalar transport
    species2_primitiveVenturi_mixingmodel_boundedscalar           = TestCase('species2_primitiveVenturi_mixingmodel_boundedscalar')
    species2_primitiveVenturi_mixingmodel_boundedscalar.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_boundedscalar.cfg_file  = "species2_primitiveVenturi_mixingmodel_boundedscalar.cfg"
    species2_primitiveVenturi_mixingmodel_boundedscalar.test_iter = 50
    species2_primitiveVenturi_mixingmodel_boundedscalar.test_vals = [-5.689876, -4.507214, -4.611632, -6.120372, -0.118388, -5.705973, 5.000000, -1.437234, 5.000000, -4.924953, 5.000000, -1.768691, 0.000313, 0.000313, 0.000000, 0.000000]
    test_list.append(species2_primitiveVenturi_mixingmodel_boundedscalar)

    # 2 species (1 eq) primitive venturi mixing using mixing model including viscosity, thermal conductivity and inlet markers for SA turbulence model
    species2_primitiveVenturi_mixingmodel_viscosity           = TestCase('species2_primitiveVenturi_mixingmodel_viscosity')
    species2_primitiveVenturi_mixingmodel_viscosity.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_viscosity.cfg_file  = "species2_primitiveVenturi_mixingmodel_viscosity.cfg"
    species2_primitiveVenturi_mixingmodel_viscosity.test_iter = 50
    species2_primitiveVenturi_mixingmodel_viscosity.test_vals = [-5.232339, -3.617118, -3.857221, -7.533847, -5.126646, 5.000000, -1.682962, 5.000000, -3.474078, 5.000000, -2.086859, 2.495548, 0.985490, 0.600234, 0.909824]
    test_list.append(species2_primitiveVenturi_mixingmodel_viscosity)

    # 2 species (1 eq) primitive venturi mixing using mixing model including heat capacity and mass diffusivity
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2           = TestCase('species2_primitiveVenturi_mixingmodel_heatcapacity_H2.cfg')
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.cfg_file  = "species2_primitiveVenturi_mixingmodel_heatcapacity_H2.cfg"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.test_iter = 50
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2.test_vals = [-5.830112, -4.504576, -4.654774, -7.014975, 2.317603, -5.509555, 30.000000, -6.906408, 11.000000, -8.227485, 8.000000, -9.294144, 2.077533, 1.000000, 0.600000, 0.477533]
    test_list.append(species2_primitiveVenturi_mixingmodel_heatcapacity_H2)

    # 2 species (1 eq) primitive venturi mixing using mixing model including heat capacity and mass diffusivity NonDimensional case
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND           = TestCase('species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.cfg')
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.cfg_file  = "species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.cfg"
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.test_iter = 50
    species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND.test_vals = [-5.429311, -4.805894, -4.943400, -7.978127, 2.012275, -5.115722, 10.000000, -2.264471, 2.000000, -5.061843, 4.000000, -5.055035, 2.077660, 1.000000, 0.600000, 0.477660]
    test_list.append(species2_primitiveVenturi_mixingmodel_heatcapacity_H2_ND)

    # 2 species (1 eq) primitive venturi mixing
    species2_primitiveVenturi           = TestCase('species2_primitiveVenturi')
    species2_primitiveVenturi.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi.cfg_file  = "species2_primitiveVenturi.cfg"
    species2_primitiveVenturi.test_iter = 50
    species2_primitiveVenturi.test_vals = [-5.919132, -4.941085, -4.919763, -5.717215, -1.405791, -6.087254, 5.000000, -0.524653, 5.000000, -2.608809, 5.000000, -0.351723, 0.000113, 0.000112, 0.000000, 0.000000]
    test_list.append(species2_primitiveVenturi)

    # 2 species (1 eq) primitive venturi mixing with bounded scalar transport
    species_primitiveVenturi_boundedscalar             = TestCase('species2_primitiveVenturi_bounded_scalar')
    species_primitiveVenturi_boundedscalar.cfg_dir     = "species_transport/venturi_primitive_3species"
    species_primitiveVenturi_boundedscalar.cfg_file    = "species2_primitiveVenturi_boundedscalar.cfg"
    species_primitiveVenturi_boundedscalar.test_iter   = 50
    species_primitiveVenturi_boundedscalar.test_vals   = [-5.534441, -4.368794, -4.465911, -5.935827, -0.867763, -5.633723, 5.000000, -1.470749, 5.000000, -4.161492, 5.000000, -1.724769, 0.000434, 0.000434, 0.000000, 0.000000]
    test_list.append(species_primitiveVenturi_boundedscalar)

    # 2 species (1 eq) primitive venturi mixing using mixing model including inlet markers for turbulent intensity and viscosity ratios
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS           = TestCase('species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.cfg')
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.cfg_file  = "species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.cfg"
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.test_iter = 50
    species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS.test_vals = [-4.565556, -1.773987, -1.569972, -0.978893, 1.693213, -3.843844, 23.000000, -5.040721, 10.000000, -5.559387, 3.000000, -5.370252, 2.000000, 1.000000, 0.000000, 1.000000]
    test_list.append(species2_primitiveVenturi_mixingmodel_TURBULENT_MARKERS)

    # 3 species (2 eq) primitive venturi mixing with inlet files.
    # Note that the residuals are exactly the same as for the non-inlet case which should be the case for a fresh inlet file.
    species3_primitiveVenturi_inletFile           = TestCase('species3_primitiveVenturi_inletFile')
    species3_primitiveVenturi_inletFile.cfg_dir   = "species_transport/venturi_primitive_3species"
    species3_primitiveVenturi_inletFile.cfg_file  = "species3_primitiveVenturi_inletFile.cfg"
    species3_primitiveVenturi_inletFile.test_iter = 50
    species3_primitiveVenturi_inletFile.test_vals = [-5.989944, -5.011655, -4.990373, -5.787315, -1.475219, -6.316068, -6.446277, 5.000000, -0.525858, 5.000000, -2.609037, 5.000000, -0.362429]
    test_list.append(species3_primitiveVenturi_inletFile)

    # rectangle passive transport validation
    species_passive_val           = TestCase('species_passive_val')
    species_passive_val.cfg_dir   = "species_transport/passive_transport_validation"
    species_passive_val.cfg_file  = "passive_transport.cfg"
    species_passive_val.test_iter = 50
    species_passive_val.test_vals = [-16.561013, -16.280437, -16.910003, -4.257599, 10.000000, -4.457215, 8.000000, -5.193350, 0.186610, 0.000000]
    species_passive_val.test_vals_aarch64 = [-16.538551, -16.312552, -16.882823, -4.257599, 10, -4.585464, 8, -5.19335, 0.18661, 0]
    test_list.append(species_passive_val)

    # species transport, 3 species with multizone (2 fluid regions)
    species3_multizone_restart           = TestCase('species3_multizone_restart')
    species3_multizone_restart.cfg_dir   = "species_transport/multizone"
    species3_multizone_restart.cfg_file  = "configMaster.cfg"
    species3_multizone_restart.test_iter = 5
    species3_multizone_restart.test_vals = [-3.715625, -3.070578]
    species3_multizone_restart.multizone = True
    test_list.append(species3_multizone_restart)

    #####################
    ## CGNS writer ###
    #####################

    # CGNS writer
    cgns_writer = TestCase('cgns_writer')
    cgns_writer.cfg_dir = "cgns_writer"
    cgns_writer.cfg_file = "config.cfg"
    cgns_writer.test_iter = 1
    cgns_writer.test_vals = [-2.974473, 0.640256, 5.371028, -6.732060]
    cgns_writer.new_output = True
    test_list.append(cgns_writer)

    ######################################
    ### RUN CHT TEST WITH FILEDIFF     ###
    ######################################

    # 2D planar laminar premixed methane flame on isothermal burner with conjugate heat transfer in cooling fin (restart)
    cfd_flamelet_ch4_cht = TestCase('cfd_flamelet_ch4_cht')
    cfd_flamelet_ch4_cht.cfg_dir = "flamelet/03_laminar_premixed_ch4_flame_cht_cfd"
    cfd_flamelet_ch4_cht.cfg_file = "lam_prem_ch4_cht_cfd_master.cfg"
    cfd_flamelet_ch4_cht.test_iter = 5
    cfd_flamelet_ch4_cht.test_vals = [-10.166224, -9.299576, -10.022482, -11.907558, -3.447626, -12.290421, -13.719415, -6.864920]
    cfd_flamelet_ch4_cht.timeout = 1600
    cfd_flamelet_ch4_cht.multizone = True
    test_list.append(cfd_flamelet_ch4_cht)

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

    # Inviscid NACA0012 (triangles)
    naca0012_cst = TestCase('naca0012_cst')
    naca0012_cst.cfg_dir = "deformation/cst"
    naca0012_cst.cfg_file = "naca0012.cfg"
    naca0012_cst.test_iter = 10
    naca0012_cst.test_vals = [0.000385514] #residual
    naca0012_cst.command = TestCase.Command("mpirun -n 2", "SU2_DEF")
    naca0012_cst.timeout = 1600
    naca0012_cst.tol = 1e-8

    pass_list.append(naca0012_cst.run_def())
    test_list.append(naca0012_cst)

    # 2D FD streamwise periodic cht, avg temp obj func
    fd_sp_pinArray_cht_2d_dp_hf = TestCase('fd_sp_pinArray_cht_2d_dp_hf')
    fd_sp_pinArray_cht_2d_dp_hf.cfg_dir = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    fd_sp_pinArray_cht_2d_dp_hf.cfg_file = "FD_configMaster.cfg"
    fd_sp_pinArray_cht_2d_dp_hf.test_iter = 100
    fd_sp_pinArray_cht_2d_dp_hf.command = TestCase.Command(exec = "finite_differences.py", param = "-z 2 -n 2 -f")
    fd_sp_pinArray_cht_2d_dp_hf.timeout = 1600
    fd_sp_pinArray_cht_2d_dp_hf.comp_threshold  = 1e-6
    fd_sp_pinArray_cht_2d_dp_hf.tol_file_percent = 0.2
    fd_sp_pinArray_cht_2d_dp_hf.reference_file = "of_grad_findiff.csv.ref"
    fd_sp_pinArray_cht_2d_dp_hf.reference_file_aarch64 = "of_grad_findiff_aarch64.csv.ref"
    fd_sp_pinArray_cht_2d_dp_hf.test_file = "FINDIFF/of_grad_findiff.csv"
    fd_sp_pinArray_cht_2d_dp_hf.multizone = True

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
