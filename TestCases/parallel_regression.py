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

# make print(*args) function available in PY2.6+, doesn't work on PY < 2.6
from __future__ import print_function

import sys
from TestCase import TestCase

def main():
    '''This program runs SU2 and ensures that the output matches specified values.
       This will be used to do checks when code is pushed to github
       to make sure nothing is broken. '''

    test_list = []

    ##########################################
    ### Combustion with Detailed Chemistry ###
    ##########################################

    # 2D planar laminar premixed methane flame on isothermal burner
    cfd_cantera_dc_ch4 = TestCase('cfd_cantera_dc_ch4')
    cfd_cantera_dc_ch4.cfg_dir = "combustion/1step_methane_air"
    cfd_cantera_dc_ch4.cfg_file = "2d_slit_burner_fluid.cfg"
    cfd_cantera_dc_ch4.test_iter = 10
    cfd_cantera_dc_ch4.test_vals = [-12.897453, -15.198159, -12.810575, 10.000000, -1.881353, 1000.000000]
    cfd_cantera_dc_ch4.new_output = True
    test_list.append(cfd_cantera_dc_ch4)

    # 2D planar laminar premixed hydrogen flame on isothermal burner (restart)
    cfd_cantera_dc_h2 = TestCase('cfd_cantera_dc_h2')
    cfd_cantera_dc_h2.cfg_dir = "combustion/h2o2_hydrogen_air"
    cfd_cantera_dc_h2.cfg_file = "2d_slit_burner_h2.cfg"
    cfd_cantera_dc_h2.test_iter = 10
    cfd_cantera_dc_h2.test_vals = [-15.993819, -17.134091, -32.000000, 5.000000, -2.029177, 100.000000]
    cfd_cantera_dc_h2.new_output = True
    test_list.append(cfd_cantera_dc_h2)

    #######################
    ### Flamelet solver ###
    #######################

    # 2D planar laminar premixed methane flame on isothermal burner (restart)
    cfd_flamelet_ch4 = TestCase('cfd_flamelet_ch4')
    cfd_flamelet_ch4.cfg_dir = "flamelet/01_laminar_premixed_ch4_flame_cfd"
    cfd_flamelet_ch4.cfg_file = "lam_prem_ch4_cfd.cfg"
    cfd_flamelet_ch4.test_iter = 10
    cfd_flamelet_ch4.test_vals = [-10.539588, -12.650786, -5.649832, -13.470491, -17.249262, 94.854000]
    cfd_flamelet_ch4.new_output = True
    test_list.append(cfd_flamelet_ch4)

   # axisymmetric 2D planar laminar premixed methane flame on isothermal burner (restart)
    cfd_flamelet_ch4_axi = TestCase('cfd_flamelet_ch4_axi')
    cfd_flamelet_ch4_axi.cfg_dir = "flamelet/05_laminar_premixed_ch4_flame_cfd_axi"
    cfd_flamelet_ch4_axi.cfg_file = "lam_prem_ch4_cfd_axi.cfg"
    cfd_flamelet_ch4_axi.test_iter = 10
    cfd_flamelet_ch4_axi.test_vals = [-9.774148, -11.255837, -10.017962, -11.850828, -4.916261, 72.669000]
    cfd_flamelet_ch4_axi.new_output = True
    test_list.append(cfd_flamelet_ch4_axi)

    # 2D planar laminar partially premixed flame on isothermal burner and heat exchanger (restart)
    cfd_flamelet_ch4_partial_premix = TestCase('cfd_flamelet_ch4_partial_premix')
    cfd_flamelet_ch4_partial_premix.cfg_dir = "flamelet/06_laminar_partial_premixed_ch4_flame_cfd"
    cfd_flamelet_ch4_partial_premix.cfg_file = "lam_partial_prem_ch4_cfd.cfg"
    cfd_flamelet_ch4_partial_premix.test_iter = 10
    cfd_flamelet_ch4_partial_premix.test_vals = [-10.482722, -3.621055, -12.932224, -10.154391, 10.000000]
    cfd_flamelet_ch4_partial_premix.new_output = True
    test_list.append(cfd_flamelet_ch4_partial_premix)

    # 2D planar laminar premixed hydrogen flame on isothermal burner with heat exchanger emulator (restart)
    cfd_flamelet_h2 = TestCase('cfd_flamelet_h2')
    cfd_flamelet_h2.cfg_dir = "flamelet/07_laminar_premixed_h2_flame_cfd"
    cfd_flamelet_h2.cfg_file = "laminar_premixed_h2_flame_cfd.cfg"
    cfd_flamelet_h2.test_iter = 5
    cfd_flamelet_h2.test_vals = [-8.036794, -8.372668, -1.842800, -9.388447]
    test_list.append(cfd_flamelet_h2)

    # Flame ignition methods
    flame_init_methods = TestCase("flame_init_methods")
    flame_init_methods.cfg_dir = "flamelet/08_flame_init_method"
    flame_init_methods.cfg_file ="flamelet_initialization.py"
    flame_init_methods.command = TestCase.Command("mpirun -np 2", "python", "flamelet_initialization.py")
    flame_init_methods.test_vals=[1600, 1600]
    flame_init_methods.new_output = True
    test_list.append(flame_init_methods)

    # 2D laminar premixed ch4-air flame, transient flame propagation
    cfd_flamelet_ch4_unsteady = TestCase('cfd_flamelet_ch4_unsteady')
    cfd_flamelet_ch4_unsteady.cfg_dir = "flamelet/09_laminar_premixed_ch4_flame_unsteady"
    cfd_flamelet_ch4_unsteady.cfg_file = "lam_prem_ch4_unsteady.cfg"
    cfd_flamelet_ch4_unsteady.test_iter = 5
    cfd_flamelet_ch4_unsteady.test_vals = [-8.856417, -8.095249, -9.153745, -9.321678]
    cfd_flamelet_ch4_unsteady.test_vals_aarch64 = [-8.855500, -8.095195, -9.153704, -9.321686]
    test_list.append(cfd_flamelet_ch4_unsteady)

    #########################
    ## NEMO solver ###
    #########################

    # Adiabatic thermal bath
    thermalbath = TestCase('thermalbath')
    thermalbath.cfg_dir = "nonequilibrium/thermalbath/finitechemistry"
    thermalbath.cfg_file = "thermalbath.cfg"
    thermalbath.test_iter = 10
    thermalbath.test_vals = [0.945997, 0.945997, -11.860991, -11.890652, -32.000000, 10.013239]
    test_list.append(thermalbath)

    # Adiabatic thermal bath
    ionized = TestCase('ionized')
    ionized.cfg_dir = "nonequilibrium/thermalbath/finitechemistry"
    ionized.cfg_file = "weakly_ionized.cfg"
    ionized.test_iter = 10
    ionized.test_vals = [-28.173855, -10.350260, -11.719875, -16.100105, -16.850783, -13.485322, -23.431375, -32.000000, -4.558277, 0.000000, 0.000000]
    ionized.test_vals_aarch64 = [-29.816386, -10.729986, -11.720016, -17.484469, -18.237891, -15.241605, -24.956918, -32.000000, -5.727244, 0.000000, 0.000000]
    test_list.append(ionized)

    # Adiabatic frozen thermal bath
    thermalbath_frozen = TestCase('thermalbath_frozen')
    thermalbath_frozen.cfg_dir = "nonequilibrium/thermalbath/frozen"
    thermalbath_frozen.cfg_file = "thermalbath_frozen.cfg"
    thermalbath_frozen.test_iter = 10
    thermalbath_frozen.test_vals = [-32.000000, -32.000000, -11.818647, -11.848553, -32.000000, 10.013545]
    test_list.append(thermalbath_frozen)

    # Inviscid single wedge, ausm, implicit
    invwedge_a = TestCase('invwedge_ausm')
    invwedge_a.cfg_dir = "nonequilibrium/invwedge"
    invwedge_a.cfg_file = "invwedge_ausm.cfg"
    invwedge_a.test_iter = 10
    invwedge_a.test_vals = [-1.081764, -1.606527, -18.299923, -18.627313, -18.573325, 2.234595, 1.854912, 5.278328, 0.837083]
    invwedge_a.test_vals_aarch64 = [-1.081764, -1.606527, -18.299923, -18.627313, -18.573325, 2.234595, 1.854912, 5.278328, 0.837083]
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
    invwedge_msw.test_vals = [-1.206414, -1.731177, -18.301825, -18.629206, -18.575226, 2.111923, 1.660408, 5.150501, 0.710186]
    invwedge_msw.test_vals_aarch64 = [-1.206414, -1.731177, -18.299771, -18.627181, -18.573171, 2.111923, 1.660408, 5.150501, 0.710186]
    test_list.append(invwedge_msw)

    # Inviscid single wedge, roe, implicit
    invwedge_roe = TestCase('invwedge_roe')
    invwedge_roe.cfg_dir = "nonequilibrium/invwedge"
    invwedge_roe.cfg_file = "invwedge_roe.cfg"
    invwedge_roe.test_iter = 10
    invwedge_roe.test_vals = [-0.856717, -1.381479, -16.733080, -17.063919, -17.006196, 2.492274, 1.625077, 5.504576, 1.052309]
    invwedge_roe.test_vals_aarch64 = [-1.156325, -1.681088, -17.726749, -18.056513, -17.999860, 2.159488, 1.793498, 5.198638, 0.769831]
    test_list.append(invwedge_roe)

    # Inviscid single wedge, lax, implicit
    invwedge_lax = TestCase('invwedge_lax')
    invwedge_lax.cfg_dir = "nonequilibrium/invwedge"
    invwedge_lax.cfg_file = "invwedge_lax.cfg"
    invwedge_lax.test_iter = 10
    invwedge_lax.test_vals = [-0.882120, -1.406883, -32.000000, -32.000000, -24.953606, 2.447071, 1.853423, 5.480164, 1.047015]
    invwedge_lax.test_vals_aarch64 = [-0.882120, -1.406883, -32.000000, -32.000000, -24.953606, 2.447071, 1.853423, 5.480164, 1.047015]
    test_list.append(invwedge_lax)

    # Inviscid single wedge, implicit, AUSM+M scheme
    invwedge_ausm_m = TestCase('invwedge_ausm_m')
    invwedge_ausm_m.cfg_dir = "nonequilibrium/invwedge"
    invwedge_ausm_m.cfg_file = "invwedge_am.cfg"
    invwedge_ausm_m.test_iter = 10
    invwedge_ausm_m.test_vals = [-1.173033, -1.697796, -16.739586, -17.063491, -17.012692, 2.124519, 1.963804, 5.182881, 0.747539]
    invwedge_ausm_m.test_vals_aarch64 = [-1.173033, -1.697796, -16.739586, -17.063491, -17.012692, 2.124519, 1.963804, 5.182881, 0.747539]
    test_list.append(invwedge_ausm_m)

    # Inviscid single wedge, implicit, NEMO supersonic inlet
    invwedge_ss_inlet = TestCase('invwedge_ss_inlet')
    invwedge_ss_inlet.cfg_dir = "nonequilibrium/invwedge"
    invwedge_ss_inlet.cfg_file = "invwedge_ss_inlet.cfg"
    invwedge_ss_inlet.test_iter = 10
    invwedge_ss_inlet.test_vals = [-1.081061, -1.605824, -18.246258, -18.575522, -18.519330, 2.235398, 1.855022, 5.278994, 0.837755]
    invwedge_ss_inlet.test_vals_aarch64 = [-1.081061, -1.605824, -18.246260, -18.575547, -18.519333, 2.235398, 1.855022, 5.278994, 0.837755]
    test_list.append(invwedge_ss_inlet)

    # Viscous single cone - axisymmetric
    visc_cone = TestCase('visc_cone')
    visc_cone.cfg_dir = "nonequilibrium/visc_wedge"
    visc_cone.cfg_file = "axi_visccone.cfg"
    visc_cone.test_iter = 10
    visc_cone.test_vals = [-5.275015, -5.799160, -20.292303, -20.436143, -20.888373, 1.156776, -3.217423, -0.013712, 0.092252, 27878.000000]
    visc_cone.test_vals_aarch64 = [-5.275015, -5.799159, -20.292303, -20.436143, -20.888374, 1.156775, -3.217422, -0.013713, 0.092252, 27878.000000]
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
    super_cat.test_vals = [-5.296256, -5.820629, -21.102797, -21.154470, -21.132197, 1.137393, -3.220148, -0.025962, 0.247096, 27884.000000]
    test_list.append(super_cat)

    # Viscous single wedge - partially catalytic walls
    partial_cat = TestCase('partial_cat')
    partial_cat.cfg_dir = "nonequilibrium/visc_wedge"
    partial_cat.cfg_file = "partial_cat.cfg"
    partial_cat.test_iter = 10
    partial_cat.test_vals = [-5.210302, -5.735065, -20.880448, -20.825971, -23.475263, 1.806201, -2.813952, -0.078400, 0.495606, 29020.000000]
    test_list.append(partial_cat)

    # Viscous cylinder, ionization, Gupta-Yos, marched 100 iterations from its restart.
    ion_gy = TestCase('ion_gy')
    ion_gy.cfg_dir = "nonequilibrium/visc_cylinder"
    ion_gy.cfg_file = "cyl_ion_gy.cfg"
    ion_gy.test_iter = 99
    ion_gy.test_vals = [-12.191649, -4.245292, -4.904190, -5.585803, -5.472455, -5.057619, -7.442352, 3.429183, 3.433581, -0.014861, 0.000001, 90357.000000]
    ion_gy.tol = 0.01
    test_list.append(ion_gy)

    ##########################
    ### Compressible Euler ###
    ##########################

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 20
    channel.test_vals = [3.612169, 0.104888, 0.172110, 1.000000]
    test_list.append(channel)

    # NACA0012, native SU2 binary mesh format (.su2b)
    # First, SU2_DEF converts mesh_NACA0012_inv.su2 into the native SU2 binary
    # mesh format. The regression case below then loads that freshly generated
    # .su2b file, using a copy of inv_NACA0012_Roe.cfg with MESH_FILENAME and
    # MESH_FORMAT swapped to point at it.
    naca0012_su2bin_convert           = TestCase('naca0012_su2bin_convert')
    naca0012_su2bin_convert.cfg_dir   = "euler/naca0012"
    naca0012_su2bin_convert.cfg_file  = "mesh_su2_to_su2bin.cfg"
    naca0012_su2bin_convert.command   = TestCase.Command("mpirun -n 2", "SU2_DEF")
    test_list.append(naca0012_su2bin_convert)

    naca0012_su2bin_cfg_path = "euler/naca0012/inv_NACA0012_Roe_su2bin.cfg"
    with open("euler/naca0012/inv_NACA0012_Roe.cfg", 'r') as f:
        naca0012_su2bin_cfg = f.read()
    naca0012_su2bin_cfg = naca0012_su2bin_cfg.replace(
        "MESH_FILENAME= mesh_NACA0012_inv.su2", "MESH_FILENAME= mesh_NACA0012_inv_su2bin")
    naca0012_su2bin_cfg = naca0012_su2bin_cfg.replace(
        "MESH_FORMAT= SU2\n", "MESH_FORMAT= SU2B\n")
    with open(naca0012_su2bin_cfg_path, 'w') as f:
        f.write(naca0012_su2bin_cfg)

    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe_su2bin.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.024996, -3.520626, 0.275033, 0.009889]
    test_list.append(naca0012)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-4.188807, 1.541867, -0.249526, 0.043953]
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-4.787340, -4.166981, 0.277687, 0.006533]
    oneram6.timeout   = 3200
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-3.887564, 1.642912, 0.300963, 0.019476]
    test_list.append(fixedCL_naca0012)

    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals = [-1.315882, 4.158313, 0.002698, 0.114644]
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
    bluntbody.test_vals = [0.565207, 6.945909, -0.000001, 1.970304]
    test_list.append(bluntbody)

    # Equivalent area NACA64-206
    ea_naca64206           = TestCase('ea_naca64206')
    ea_naca64206.cfg_dir   = "optimization_euler/equivalentarea_naca64206"
    ea_naca64206.cfg_file  = "NACA64206.cfg"
    ea_naca64206.test_iter = 10
    ea_naca64206.test_vals = [-1.126282, -0.474499, -0.002405, 67775.000000]
    test_list.append(ea_naca64206)

    # SUPERSONIC FLOW PAST A RAMP IN A CHANNEL
    ramp = TestCase('ramp')
    ramp.cfg_dir = "euler/ramp"
    ramp.cfg_file = "inv_ramp.cfg"
    ramp.test_iter = 10
    ramp.test_vals = [-13.648309, -8.010607, -0.076277, 0.054839]
    ramp.test_vals_aarch64 = [-13.648406, -8.014579, -0.076277, 0.054839]
    test_list.append(ramp)

    ramp_msw = TestCase('ramp_msw')
    ramp_msw.cfg_dir = "euler/ramp"
    ramp_msw.cfg_file = "inv_ramp_msw.cfg"
    ramp_msw.test_iter = 100
    ramp_msw.test_vals = [-7.242600, -1.471141, -0.077507, 0.054419]
    ramp_msw.tol = [0.2, 0.2, 0.00001, 0.00001]
    test_list.append(ramp_msw)

    # MFR_coupling
    MFR_coupling           = TestCase('MFR_coupling')
    MFR_coupling.cfg_dir   = "euler/turbofan_MFR_coupling"
    MFR_coupling.cfg_file  = "MFR_coupling.cfg"
    MFR_coupling.test_iter = 100
    MFR_coupling.test_vals = [-211.250000, 150.030000, 20.151000]
    test_list.append(MFR_coupling)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 100
    flatplate.test_vals = [-8.160319, -2.686084, 0.001084, 0.036277, 2.361500, -2.325300, 0.000000, 0.000000]
    test_list.append(flatplate)

    # Supersonic laminar flat plate
    flatplate_supersonic           = TestCase('flatplate_supersonic')
    flatplate_supersonic.cfg_dir   = "navierstokes/flatplate"
    flatplate_supersonic.cfg_file  = "flatplate_supersonic.cfg"
    flatplate_supersonic.test_iter = 100
    flatplate_supersonic.test_vals = [100.000000, -9.034577, -8.611322, -9.029955, -3.700944, 0.002511]
    test_list.append(flatplate_supersonic)

    # Custom objective function
    flatplate_udobj           = TestCase('flatplate_udobj')
    flatplate_udobj.cfg_dir   = "user_defined_functions"
    flatplate_udobj.cfg_file  = "lam_flatplate.cfg"
    flatplate_udobj.test_iter = 20
    flatplate_udobj.test_vals = [-0.001015, 596.450000, 299.560000, 296.900000, 21.317000, 0.586750, 36.553000, 2.188800]
    test_list.append(flatplate_udobj)

    # Probe performance test (11 probes, ADT path)
    probe_performance_11 = TestCase('probe_performance_11')
    probe_performance_11.cfg_dir = "user_defined_functions"
    probe_performance_11.cfg_file = "test_11_probes.cfg"
    probe_performance_11.test_iter = 4
    probe_performance_11.test_vals = [-6.290996, 101020.000000, 101050.000000, 99125.000000] # RMS_DENSITY, probe1, probe6, probe11
    # Tolerances are typically 0.001 in TestCase.py
    test_list.append(probe_performance_11)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-8.557718, -3.072847, -0.018912, 1.610317, 0.000000]
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-8.737675, -3.275569, 0.050297, 1.839111, 0.000000]

    test_list.append(cylinder_lowmach)

    # Laminar NACA0012 with linelet preconditioner (see #2877), 4 ranks so that some partitions have no linelets
    lam_naca0012_linelet = TestCase('lam_naca0012_linelet')
    lam_naca0012_linelet.cfg_dir = "navierstokes/naca0012"
    lam_naca0012_linelet.cfg_file = "lam_NACA0012_LINELET.cfg"
    lam_naca0012_linelet.test_iter = 20
    lam_naca0012_linelet.test_vals = [-5.335764, -2.936176, -2.943918, 0.102456]
    lam_naca0012_linelet.command = TestCase.Command("mpirun --oversubscribe -n 4", "SU2_CFD")
    test_list.append(lam_naca0012_linelet)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-6.994786, 0.000197, 13.596249, 0]
    poiseuille.tol       = 0.001
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals = [-12.003186, -7.623957, -0.000000, 2.089953]
    poiseuille_profile.test_vals_aarch64 = [-12.007498, -7.226926, -0.000000, 2.089953]
    poiseuille_profile.tol       = [0.001, 0.001, 1e-5, 1e-5, 1e-5]
    test_list.append(poiseuille_profile)

    # 2D rotational periodic sector without limiter (Jacobian of the periodic points)
    periodic2d_no_limiter = TestCase('periodic2d_no_limiter')
    periodic2d_no_limiter.cfg_dir = "navierstokes/periodic2D"
    periodic2d_no_limiter.cfg_file = "no_limiter.cfg"
    periodic2d_no_limiter.test_iter = 100
    periodic2d_no_limiter.test_vals = [-3.217331, -0.583709, -0.629299, 2.265684, -1.057572, -810.150000]
    test_list.append(periodic2d_no_limiter)

    # The fully coupled nonlinear boundary derivative used to diverge before iteration 100 at CFL 100.
    periodic2d_no_limiter_riemann = TestCase('periodic2d_no_limiter_riemann')
    periodic2d_no_limiter_riemann.cfg_dir = "navierstokes/periodic2D"
    periodic2d_no_limiter_riemann.cfg_file = "no_limiter_riemann.cfg"
    periodic2d_no_limiter_riemann.test_iter = 100
    periodic2d_no_limiter_riemann.test_vals = [-3.216878, -0.582314, -0.628636, 2.266530, -1.057965, -811.500000]
    test_list.append(periodic2d_no_limiter_riemann)

    # 2D rotational periodic sector with multigrid
    periodic2d_multigrid = TestCase('periodic2d_multigrid')
    periodic2d_multigrid.cfg_dir = "navierstokes/periodic2D"
    periodic2d_multigrid.cfg_file = "multigrid.cfg"
    periodic2d_multigrid.test_iter = 300
    periodic2d_multigrid.test_vals = [-4.422875, -1.615366, -1.453862, 1.047230, -1.410284, -2060.900000]
    test_list.append(periodic2d_multigrid)

    # 2D rotational periodic sector with limiter
    periodic2d = TestCase('periodic2d')
    periodic2d.cfg_dir = "navierstokes/periodic2D"
    periodic2d.cfg_file = "config.cfg"
    periodic2d.test_iter = 100
    periodic2d.test_vals = [-3.266432, -0.617956, -0.620044, 2.216317, -1.007813, -1035.000000]
    test_list.append(periodic2d)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-3.030366, -5.142544, 0.790275, 0.019666, 1000.000000]
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-1.485957, 5.891491, 0.675200, 0.023427, 100.000000]
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.465087, 5.851605, 0.496009, 0.041100]
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-0.187388, 0.003749, 10.000000, -1.406182]
    test_list.append(turb_flatplate)

    # Flat plate (compressible) with species inlet
    turb_flatplate_species           = TestCase('turb_flatplate_species')
    turb_flatplate_species.cfg_dir   = "rans/flatplate"
    turb_flatplate_species.cfg_file  = "turb_SA_flatplate_species.cfg"
    turb_flatplate_species.test_iter = 20
    turb_flatplate_species.test_vals = [-5.359058, -1.881677, -3.075574, 0.124434, -3.602252, 2.000000, -1.239891, 1.000000, -1.158675, 3.000000, -0.849018, 0.999876, 0.999876]
    test_list.append(turb_flatplate_species)

    # Flat plate SST compressibility correction Wilcox
    turb_flatplate_CC_Wilcox = TestCase('turb_flatplate_CC_Wilcox')
    turb_flatplate_CC_Wilcox.cfg_dir   = "rans/flatplate"
    turb_flatplate_CC_Wilcox.cfg_file  = "turb_SST_flatplate_compressibility_Wilcox.cfg"
    turb_flatplate_CC_Wilcox.test_iter = 20
    turb_flatplate_CC_Wilcox.test_vals = [-1.204997, 2.149390, 1.396996, 5.216302, -3.821981, 9.707535]
    test_list.append(turb_flatplate_CC_Wilcox)

    # Flat plate SST compressibility correction Sarkar
    turb_flatplate_CC_Sarkar = TestCase('turb_flatplate_CC_Sarkar')
    turb_flatplate_CC_Sarkar.cfg_dir   = "rans/flatplate"
    turb_flatplate_CC_Sarkar.cfg_file  = "turb_SST_flatplate_compressibility_Sarkar.cfg"
    turb_flatplate_CC_Sarkar.test_iter = 20
    turb_flatplate_CC_Sarkar.test_vals = [-1.204997, 2.149390, 1.396996, 5.216302, -3.821978, 9.707535]
    test_list.append(turb_flatplate_CC_Sarkar)

    # FLAT PLATE, ROUGHNESS BC KNOPP SST
    turb_flatplate_sst_roughBCKnopp           = TestCase('turb_sst_flatplate_roughBCKnopp')
    turb_flatplate_sst_roughBCKnopp.cfg_dir   = "rans/flatplate/roughness/bc_knopp"
    turb_flatplate_sst_roughBCKnopp.cfg_file  = "turb_SST_flatplate_roughBCKnopp.cfg"
    turb_flatplate_sst_roughBCKnopp.test_iter = 10
    turb_flatplate_sst_roughBCKnopp.test_vals = [-5.066141, -2.470493, -2.851854, 0.439238, -2.595926, 1.496827, -0.188074, 0.004571]
    test_list.append(turb_flatplate_sst_roughBCKnopp)

    # FLAT PLATE, ROUGHNESS BC AUPOIX SST
    turb_flatplate_sst_roughBCAupoix           = TestCase('turb_sst_flatplate_roughBCAupoix')
    turb_flatplate_sst_roughBCAupoix.cfg_dir   = "rans/flatplate/roughness/bc_aupoix"
    turb_flatplate_sst_roughBCAupoix.cfg_file  = "turb_SST_flatplate_roughBCAupoix.cfg"
    turb_flatplate_sst_roughBCAupoix.test_iter = 10
    turb_flatplate_sst_roughBCAupoix.test_vals = [-5.290263, -2.305835, -2.892120, 0.214802, -1.393793, 3.192602, -0.188752, 0.007182]
    test_list.append(turb_flatplate_sst_roughBCAupoix)

    # FLAT PLATE, WALL FUNCTIONS, COMPRESSIBLE SA, restart from a solution without wall functions
    turb_wallfunction_flatplate_sa_restart = TestCase('turb_sa_wallfunction_flatplate_restart')
    turb_wallfunction_flatplate_sa_restart.cfg_dir = "wallfunctions/flatplate/compressible_SA"
    turb_wallfunction_flatplate_sa_restart.cfg_file = "turb_SA_flatplate_restart.cfg"
    turb_wallfunction_flatplate_sa_restart.test_iter = 10
    turb_wallfunction_flatplate_sa_restart.test_vals = [-3.603551, -1.991100, -1.848182, 1.875994, -5.284153, 10.000000, -2.255956, -0.064942, 0.002054]
    test_list.append(turb_wallfunction_flatplate_sa_restart)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.419976, -6.631992, 0.238577, 0.159694, 0.000000]
    turb_oneram6.timeout   = 3200
    test_list.append(turb_oneram6)

    # ONERA M6 Wing - vorticity confinement
    turb_oneram6_vc = TestCase('turb_oneram6_vc')
    turb_oneram6_vc.cfg_dir = "rans/oneram6"
    turb_oneram6_vc.cfg_file = "turb_ONERAM6_vc.cfg"
    turb_oneram6_vc.test_iter = 15
    turb_oneram6_vc.test_vals = [-2.295077, -6.572734, 0.234428, 0.145063, 0.000000]
    turb_oneram6_vc.timeout = 3200
    test_list.append(turb_oneram6_vc)

    # ONERA M6 Wing - Newton-Krylov
    turb_oneram6_nk           = TestCase('turb_oneram6_nk')
    turb_oneram6_nk.cfg_dir   = "rans/oneram6"
    turb_oneram6_nk.cfg_file  = "turb_ONERAM6_nk.cfg"
    turb_oneram6_nk.test_iter = 20
    turb_oneram6_nk.test_vals = [-4.848121, -4.450658, -11.426212, 0.221522, 0.049358, 2.000000, -0.879182, 10.000000]
    turb_oneram6_nk.timeout   = 600
    turb_oneram6_nk.tol       = 0.0001
    test_list.append(turb_oneram6_nk)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 5
    turb_naca0012_sa.test_vals = [-12.037552, -16.376949, 1.080346, 0.018385, 20.000000, -1.564082, 20.000000, -4.180917, 0.000000]
    turb_naca0012_sa.test_vals_aarch64 = [-12.037489, -16.376949, 1.080346, 0.018385, 20.000000, -1.564143, 20.000000, -4.180945, 0.000000]
    turb_naca0012_sa.timeout   = 3200
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-12.094736, -15.251092, -5.906365, 1.070413, 0.015775, -2.376001, 0.000000]
    turb_naca0012_sst.test_vals_aarch64 = [-12.075620, -15.246688, -5.861276, 1.070036, 0.015841, -1.991001, 0.000000]
    turb_naca0012_sst.timeout   = 3200
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals = [-12.081966, -14.837177, -5.733436, 1.000893, 0.019109, -2.241123]
    turb_naca0012_sst_sust.test_vals_aarch64 = [-12.073964, -14.836726, -5.732390, 1.000050, 0.019144, -2.229074]
    turb_naca0012_sst_sust.timeout   = 3200
    test_list.append(turb_naca0012_sst_sust)

    # NACA0012 (SST, 2003m, Vorticity)
    turb_naca0012_sst_2003_Vm           = TestCase('turb_naca0012_sst_2003_Vm')
    turb_naca0012_sst_2003_Vm.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_2003_Vm.cfg_file  = "turb_NACA0012_sst_2003-Vm.cfg"
    turb_naca0012_sst_2003_Vm.test_iter = 10
    turb_naca0012_sst_2003_Vm.test_vals = [-10.168739, -3.668709, 1.060563, 0.019147, -2.282800]
    turb_naca0012_sst_2003_Vm.timeout   = 3200
    test_list.append(turb_naca0012_sst_2003_Vm)

    # NACA0012 (SST, 1994m Kato-Launder)
    turb_naca0012_sst_1994_KLm           = TestCase('turb_naca0012_sst_1994_KLm')
    turb_naca0012_sst_1994_KLm.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_1994_KLm.cfg_file  = "turb_NACA0012_sst_1994-KLm.cfg"
    turb_naca0012_sst_1994_KLm.test_iter = 10
    turb_naca0012_sst_1994_KLm.test_vals = [-10.390871, -3.701584, 1.062164, 0.019076, -2.426123]
    turb_naca0012_sst_1994_KLm.timeout   = 3200
    test_list.append(turb_naca0012_sst_1994_KLm)

    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [-10.440016, 0.774146, 1.022363, 0.040546, -3.736437]
    turb_naca0012_sst_fixedvalues.timeout   = 3200
    test_list.append(turb_naca0012_sst_fixedvalues)

    # NACA0012 (SST, explicit Euler for flow and turbulence equations)
    turb_naca0012_sst_expliciteuler           = TestCase('turb_naca0012_sst_expliciteuler')
    turb_naca0012_sst_expliciteuler.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_expliciteuler.cfg_file  = "turb_NACA0012_sst_expliciteuler.cfg"
    turb_naca0012_sst_expliciteuler.test_iter = 10
    turb_naca0012_sst_expliciteuler.test_vals = [-3.532365, -3.157224, 3.743381, 1.124798, 0.501715, -16.000000]
    turb_naca0012_sst_expliciteuler.timeout   = 3200
    test_list.append(turb_naca0012_sst_expliciteuler)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389724, -8.410479, 0.000048, 0.056344]
    propeller.timeout   = 3200
    test_list.append(propeller)

    # Actuator disk BEM method for propeller
    actuatordisk_bem = TestCase('actuatordisk_bem')
    actuatordisk_bem.cfg_dir = "rans/actuatordisk_bem"
    actuatordisk_bem.cfg_file = "actuatordisk_bem.cfg"
    actuatordisk_bem.test_iter = 15
    actuatordisk_bem.test_vals = [-5.389236, -10.319122, 0.001362, -0.376528]
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
    axi_rans_air_nozzle_restart.test_vals = [-2.662833, 2.912056, -2.707804, 1.889981, 0.000000]
