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
    wedge.test_vals = [-4.017956, 1.709065, -0.249527, 0.043953]
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
    ramp.test_vals = [-13.648343, -8.010490, -0.076277, 0.054839]
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
    flatplate.test_vals = [-8.165135, -2.690978, 0.001084, 0.036277, 2.361500, -2.325300, 0.000000, 0.000000]
    test_list.append(flatplate)

    # Supersonic laminar flat plate
    flatplate_supersonic           = TestCase('flatplate_supersonic')
    flatplate_supersonic.cfg_dir   = "navierstokes/flatplate"
    flatplate_supersonic.cfg_file  = "flatplate_supersonic.cfg"
    flatplate_supersonic.test_iter = 100
    flatplate_supersonic.test_vals = [100.000000, -2.940836, -0.677210, -0.570631, 2.495392, 0.001469]
    test_list.append(flatplate_supersonic)

    # Custom objective function
    flatplate_udobj           = TestCase('flatplate_udobj')
    flatplate_udobj.cfg_dir   = "user_defined_functions"
    flatplate_udobj.cfg_file  = "lam_flatplate.cfg"
    flatplate_udobj.test_iter = 20
    flatplate_udobj.test_vals = [-0.001015, 596.450000, 299.550000, 296.900000, 21.318000, 0.586640, 36.553000, 2.188800]
    test_list.append(flatplate_udobj)

    # Probe performance test (11 probes, ADT path)
    probe_performance_11 = TestCase('probe_performance_11')
    probe_performance_11.cfg_dir = "user_defined_functions"
    probe_performance_11.cfg_file = "test_11_probes.cfg"
    probe_performance_11.test_iter = 4
    probe_performance_11.test_vals = [-6.290748, 101020, 101050, 99123] # RMS_DENSITY, probe1, probe6, probe11
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
    poiseuille_profile.test_vals = [-12.003181, -7.626400, -0.000000, 2.089953]
    poiseuille_profile.test_vals_aarch64 = [-12.007498, -7.226926, -0.000000, 2.089953]
    poiseuille_profile.tol       = [0.001, 0.001, 1e-5, 1e-5, 1e-5]
    test_list.append(poiseuille_profile)

    # 2D rotational periodic sector without limiter (Jacobian of the periodic points)
    periodic2d_no_limiter = TestCase('periodic2d_no_limiter')
    periodic2d_no_limiter.cfg_dir = "navierstokes/periodic2D"
    periodic2d_no_limiter.cfg_file = "no_limiter.cfg"
    periodic2d_no_limiter.test_iter = 100
    periodic2d_no_limiter.test_vals = [-3.216881, -0.582813, -0.629049, 2.266141, -1.056877, -811.470000]
    test_list.append(periodic2d_no_limiter)

    # 2D rotational periodic sector with multigrid
    periodic2d_multigrid = TestCase('periodic2d_multigrid')
    periodic2d_multigrid.cfg_dir = "navierstokes/periodic2D"
    periodic2d_multigrid.cfg_file = "multigrid.cfg"
    periodic2d_multigrid.test_iter = 300
    periodic2d_multigrid.test_vals = [-4.422778, -1.615251, -1.453736, 1.047350, -1.410225, -2060.600000]
    test_list.append(periodic2d_multigrid)

    # 2D rotational periodic sector with limiter
    periodic2d = TestCase('periodic2d')
    periodic2d.cfg_dir = "navierstokes/periodic2D"
    periodic2d.cfg_file = "config.cfg"
    periodic2d.test_iter = 100
    periodic2d.test_vals = [-3.266248, -0.617784, -0.620252, 2.216476, -1.007761, -1035.000000]
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
    turb_flatplate_species.test_vals = [-5.359043, -1.881674, -3.075839, 0.124416, -3.602253, 2.000000, -1.239286, 1.000000, -1.158672, 3.000000, -0.849003, 0.999876, 0.999876]
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
    turb_flatplate_sst_roughBCKnopp.test_vals = [-5.058634, -2.460850, -2.847064, 0.447200, -2.595042, 1.497149, -0.188079, 0.004571]
    test_list.append(turb_flatplate_sst_roughBCKnopp)

    # FLAT PLATE, ROUGHNESS BC AUPOIX SST
    turb_flatplate_sst_roughBCAupoix           = TestCase('turb_sst_flatplate_roughBCAupoix')
    turb_flatplate_sst_roughBCAupoix.cfg_dir   = "rans/flatplate/roughness/bc_aupoix"
    turb_flatplate_sst_roughBCAupoix.cfg_file  = "turb_SST_flatplate_roughBCAupoix.cfg"
    turb_flatplate_sst_roughBCAupoix.test_iter = 10
    turb_flatplate_sst_roughBCAupoix.test_vals = [-5.278728, -2.302982, -2.890419, 0.227585, -1.393786, 3.192601, -0.188729, 0.0071821]
    test_list.append(turb_flatplate_sst_roughBCAupoix)

    # FLAT PLATE, WALL FUNCTIONS, COMPRESSIBLE SA, restart from a solution without wall functions
    turb_wallfunction_flatplate_sa_restart = TestCase('turb_sa_wallfunction_flatplate_restart')
    turb_wallfunction_flatplate_sa_restart.cfg_dir = "wallfunctions/flatplate/compressible_SA"
    turb_wallfunction_flatplate_sa_restart.cfg_file = "turb_SA_flatplate_restart.cfg"
    turb_wallfunction_flatplate_sa_restart.test_iter = 10
    turb_wallfunction_flatplate_sa_restart.test_vals = [-3.603550, -1.991111, -1.848187, 1.875995, -5.284151, 10.000000, -2.255938, -0.064941, 0.002054]
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
    axi_rans_air_nozzle_restart.test_vals = [-5.144831, 0.426691, -2.894552, 1.607260, 0.000000]
    axi_rans_air_nozzle_restart.test_vals_aarch64 = [-14.143310, -9.163287, -10.858232, -5.787715, 0.000000]
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
    turb_naca0012_sst_restart_mg.test_vals = [-6.560190, -5.057149, 0.830238, -0.008812, 0.078167]
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
    inc_euler_naca0012.test_vals = [-6.920909, -6.189139, 0.531865, 0.008472]
    test_list.append(inc_euler_naca0012)

    # NACA0012 Hydrofoil
    inc_euler_naca0012_pb = TestCase('inc_euler_naca0012_pb')
    inc_euler_naca0012_pb.cfg_dir = "incomp_euler/naca0012"
    inc_euler_naca0012_pb.cfg_file = "incomp_pb_NACA0012.cfg"
    inc_euler_naca0012_pb.test_iter = 20
    inc_euler_naca0012_pb.test_vals = [-4.454834, -4.784269, 0.427449, 0.012084]
    test_list.append(inc_euler_naca0012_pb)

    # NACA0012 Hydrofoil, pressure-based with SIMPLEC
    inc_euler_naca0012_pb_simplec = TestCase('inc_euler_naca0012_pb_simplec')
    inc_euler_naca0012_pb_simplec.cfg_dir = "incomp_euler/naca0012"
    inc_euler_naca0012_pb_simplec.cfg_file = "incomp_pb_NACA0012_simplec.cfg"
    inc_euler_naca0012_pb_simplec.test_iter = 20
    inc_euler_naca0012_pb_simplec.test_vals = [-4.400610, -4.580884, 0.480955, 0.021030]
    test_list.append(inc_euler_naca0012_pb_simplec)

    # C-D nozzle with pressure inlet and mass flow outlet
    inc_nozzle           = TestCase('inc_nozzle')
    inc_nozzle.cfg_dir   = "incomp_euler/nozzle"
    inc_nozzle.cfg_file  = "inv_nozzle.cfg"
    inc_nozzle.test_iter = 20
    inc_nozzle.test_vals = [-6.051005, -5.361765, 0.006999, 0.127295]
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
    inc_lam_cylinder.test_vals = [-4.816099, -5.029676, -0.000155, 3.113777]
    test_list.append(inc_lam_cylinder)

    # Laminar cylinder, pressure-based
    inc_lam_cylinder_pb = TestCase('inc_lam_cylinder_pb')
    inc_lam_cylinder_pb.cfg_dir = "incomp_navierstokes/cylinder"
    inc_lam_cylinder_pb.cfg_file = "incomp_pb_cylinder.cfg"
    inc_lam_cylinder_pb.test_iter = 10
    inc_lam_cylinder_pb.test_vals = [-3.486601, -3.776871, 0.012481, 6.205488]
    test_list.append(inc_lam_cylinder_pb)

    # Laminar heated cylinder with polynomial fluid model, pressure-based, coupled energy
    # equation and variable density. Convergence is genuine but slow (needs ~40k iterations
    # for rms[h] to reach a low residual); this only checks a short trajectory guard, matching
    # the pattern used for other hard-to-converge cases in this suite.
    inc_poly_cylinder_pb = TestCase('inc_poly_cylinder_pb')
    inc_poly_cylinder_pb.cfg_dir = "incomp_navierstokes/cylinder"
    inc_poly_cylinder_pb.cfg_file = "pb_poly_cylinder.cfg"
    inc_poly_cylinder_pb.test_iter = 20
    inc_poly_cylinder_pb.test_vals = [-9.820675, -3.000423, 0.002745, 1.778349, -172.280000]
    test_list.append(inc_poly_cylinder_pb)

    # Unsteady heated cylinder with a polynomial fluid model, pressure-based with SIMPLEC
    inc_poly_cylinder_pb_unsteady = TestCase('inc_poly_cylinder_pb_unsteady')
    inc_poly_cylinder_pb_unsteady.cfg_dir = "incomp_navierstokes/cylinder"
    inc_poly_cylinder_pb_unsteady.cfg_file = "pb_poly_cylinder_unsteady.cfg"
    inc_poly_cylinder_pb_unsteady.test_iter = 1
    inc_poly_cylinder_pb_unsteady.test_vals = [-6.211300, -0.188735, 50.390122, -3363.100000]
    inc_poly_cylinder_pb_unsteady.unsteady = True
    test_list.append(inc_poly_cylinder_pb_unsteady)

    # Laminar sphere, Re=1. Last column: Cd=24/Re
    inc_lam_sphere          = TestCase('inc_lam_sphere')
    inc_lam_sphere.cfg_dir   = "incomp_navierstokes/sphere"
    inc_lam_sphere.cfg_file  = "sphere.cfg"
    inc_lam_sphere.test_iter = 5
    inc_lam_sphere.test_vals = [-7.599894, -8.023038, -8.135420, -8.783151, 0.121003, 25.782689, -1.574889]
    test_list.append(inc_lam_sphere)

    # Laminar sphere, Re=1, pressure-based. Only 3D pressure-based case in the regression suite.
    inc_lam_sphere_pb = TestCase('inc_lam_sphere_pb')
    inc_lam_sphere_pb.cfg_dir = "incomp_navierstokes/sphere"
    inc_lam_sphere_pb.cfg_file = "pb_sphere.cfg"
    inc_lam_sphere_pb.test_iter = 9
    inc_lam_sphere_pb.test_vals = [-6.104968, -2.307896, -2.481829, -2.570422, 0.189097, 170.764256, -6.848045]
    test_list.append(inc_lam_sphere_pb)

    # Laminar sphere, Re=1, pressure-based, automatic relaxation factors. The only case
    # in the suite that exercises USE_AUTOMATIC_RELAXATION_FACTORS, and 3D since the alpha_p
    # bug this guards against is invisible in 2D (it sums the wrong set of Jacobian diagonals).
    inc_lam_sphere_pb_urf = TestCase('inc_lam_sphere_pb_urf')
    inc_lam_sphere_pb_urf.cfg_dir = "incomp_navierstokes/sphere"
    inc_lam_sphere_pb_urf.cfg_file = "pb_sphere_urf.cfg"
    inc_lam_sphere_pb_urf.test_iter = 9
    inc_lam_sphere_pb_urf.test_vals = [-4.745234, -2.241022, -2.374688, -1.689495, 0.216842, 206.109402, -6.281761]
    test_list.append(inc_lam_sphere_pb_urf)

    # Buoyancy-driven cavity
    inc_buoyancy          = TestCase('inc_buoyancy')
    inc_buoyancy.cfg_dir   = "incomp_navierstokes/buoyancy_cavity"
    inc_buoyancy.cfg_file  = "lam_buoyancy_cavity.cfg"
    inc_buoyancy.test_iter = 20
    inc_buoyancy.test_vals = [-3.860462, -3.122401, 1.322722, -23.414000]
    test_list.append(inc_buoyancy)

    # Laminar heated cylinder with polynomial fluid model
    inc_poly_cylinder          = TestCase('inc_poly_cylinder')
    inc_poly_cylinder.cfg_dir   = "incomp_navierstokes/cylinder"
    inc_poly_cylinder.cfg_file  = "poly_cylinder.cfg"
    inc_poly_cylinder.test_iter = 20
    inc_poly_cylinder.test_vals = [-2.363389, 0.005552, 1.922971, -172.380000]
    test_list.append(inc_poly_cylinder)

    # X-coarse laminar bend as a mixed element CGNS test
    inc_lam_bend            = TestCase('inc_lam_bend')
    inc_lam_bend.cfg_dir    = "incomp_navierstokes/bend"
    inc_lam_bend.cfg_file   = "lam_bend.cfg"
    inc_lam_bend.test_iter = 10
    inc_lam_bend.test_vals = [-3.585943, -3.096592, -0.022111, 1.064110]
    test_list.append(inc_lam_bend)

    # X-coarse laminar bend as a mixed element CGNS test, pressure-based
    inc_lam_bend_pb = TestCase('inc_lam_bend_pb')
    inc_lam_bend_pb.cfg_dir = "incomp_navierstokes/bend"
    inc_lam_bend_pb.cfg_file = "pb_lam_bend.cfg"
    inc_lam_bend_pb.test_iter = 10
    inc_lam_bend_pb.test_vals = [-4.137886, -3.477106, -0.014099, 0.306589]
    test_list.append(inc_lam_bend_pb)

    # X-coarse laminar bend, pressure-based, variable density with a heated inlet
    inc_lam_bend_pb_heated = TestCase('inc_lam_bend_pb_heated')
    inc_lam_bend_pb_heated.cfg_dir = "incomp_navierstokes/bend"
    inc_lam_bend_pb_heated.cfg_file = "pb_lam_bend_heated.cfg"
    inc_lam_bend_pb_heated.test_iter = 10
    inc_lam_bend_pb_heated.test_vals = [-3.938675, -3.367579, 4.059184, 400.000000, 400.000000]
    test_list.append(inc_lam_bend_pb_heated)

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
    inc_heatTransfer_BC.test_vals = [-8.904266, -7.745636, -8.064003, -1.079197, -1671.000000]
    test_list.append(inc_heatTransfer_BC)

    # 2D pin array, periodic with a body force, flow and weakly coupled heat equation
    inc_periodic_weak_heat = TestCase('inc_periodic_weak_heat')
    inc_periodic_weak_heat.cfg_dir = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    inc_periodic_weak_heat.cfg_file = "periodic_weak_heat.cfg"
    inc_periodic_weak_heat.test_iter = 10
    inc_periodic_weak_heat.test_vals = [-4.712559, -5.747243, -5.751727, 873.702714]
    test_list.append(inc_periodic_weak_heat)

    ############################
    ### Incompressible RANS  ###
    ############################

    # NACA0012
    inc_turb_naca0012           = TestCase('inc_turb_naca0012')
    inc_turb_naca0012.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012.cfg_file  = "naca0012.cfg"
    inc_turb_naca0012.test_iter = 20
    inc_turb_naca0012.test_vals = [-4.758064, -10.974498, -0.000005, -0.028653, 3, -5.740545, 2, -5.036110]
    test_list.append(inc_turb_naca0012)

    # NACA0012, SST_SUST
    inc_turb_naca0012_sst_sust           = TestCase('inc_turb_naca0012_sst_sust')
    inc_turb_naca0012_sst_sust.cfg_dir   = "incomp_rans/naca0012"
    inc_turb_naca0012_sst_sust.cfg_file  = "naca0012_SST_SUST.cfg"
    inc_turb_naca0012_sst_sust.test_iter = 20
    inc_turb_naca0012_sst_sust.test_vals = [-7.169837, 0.332730, -0.000001, 0.312131]
    test_list.append(inc_turb_naca0012_sst_sust)

    # Flat plate, pressure-based
    inc_flatplate_pb = TestCase('inc_flatplate_pb')
    inc_flatplate_pb.cfg_dir = "incomp_rans/rough_flatplate"
    inc_flatplate_pb.cfg_file = "pb_rough_flatplate_incomp.cfg"
    inc_flatplate_pb.test_iter = 10
    inc_flatplate_pb.test_vals = [-4.063342, -9.884401, 0.000010, 0.228472]
    test_list.append(inc_flatplate_pb)

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
    turbmod_sa_bsl_rae2822.test_vals = [-3.062213, -0.409580, -0.834363, -5.133434, 0.790438, 0.018107]
    test_list.append(turbmod_sa_bsl_rae2822)

    # SA Negative
    turbmod_sa_neg_rae2822           = TestCase('turbmod_sa_neg_rae2822')
    turbmod_sa_neg_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_neg_rae2822.cfg_file  = "turb_SA_NEG_RAE2822.cfg"
    turbmod_sa_neg_rae2822.test_iter = 10
    turbmod_sa_neg_rae2822.test_vals = [1.275906, 1.062762, -2.067335, 1.301619, 0.481422, 0.000000]
    turbmod_sa_neg_rae2822.test_vals_aarch64 = [-1.345593, 1.448310, 1.208721, -0.846597, 1.248410, 0.489117, 0.000000]
    test_list.append(turbmod_sa_neg_rae2822)

    # SA Compressibility Correction
    turbmod_sa_comp_rae2822           = TestCase('turbmod_sa_comp_rae2822')
    turbmod_sa_comp_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_comp_rae2822.cfg_file  = "turb_SA_COMP_RAE2822.cfg"
    turbmod_sa_comp_rae2822.test_iter = 20
    turbmod_sa_comp_rae2822.test_vals = [-3.063286, -0.408577, -0.834597, -5.159754, 0.790991, 0.018085]
    test_list.append(turbmod_sa_comp_rae2822)

    # SA Edwards
    turbmod_sa_edw_rae2822           = TestCase('turbmod_sa_edw_rae2822')
    turbmod_sa_edw_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_edw_rae2822.cfg_file  = "turb_SA_EDW_RAE2822.cfg"
    turbmod_sa_edw_rae2822.test_iter = 20
    turbmod_sa_edw_rae2822.test_vals = [-3.069523, -0.405788, -0.835321, -5.185709, 0.794436, 0.017856]
    test_list.append(turbmod_sa_edw_rae2822)

    # SA Compressibility and Edwards
    turbmod_sa_comp_edw_rae2822           = TestCase('turbmod_sa_comp_edw_rae2822')
    turbmod_sa_comp_edw_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_comp_edw_rae2822.cfg_file  = "turb_SA_COMP_EDW_RAE2822.cfg"
    turbmod_sa_comp_edw_rae2822.test_iter = 20
    turbmod_sa_comp_edw_rae2822.test_vals = [-3.067904, -0.338707, -0.811157, -5.279475, 0.809906, 0.019797]
    test_list.append(turbmod_sa_comp_edw_rae2822)

    # SA QCR
    turbmod_sa_qcr_rae2822           = TestCase('turbmod_sa_qcr_rae2822')
    turbmod_sa_qcr_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_qcr_rae2822.cfg_file  = "turb_SA_QCR_RAE2822.cfg"
    turbmod_sa_qcr_rae2822.test_iter = 20
    turbmod_sa_qcr_rae2822.test_vals = [-3.028293, -0.418748, -0.797999, -5.088756, 0.761887, 0.016830]
    test_list.append(turbmod_sa_qcr_rae2822)

    ############################
    ###      Transition      ###
    ############################

    # Schubauer-Klebanoff Natural Transition Case
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.test_vals    = [-8.215651, -13.239431, 0.000048, 0.007983]
    test_list.append(schubauer_klebanoff_transition)

    #####################################
    ### Cont. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    contadj_naca0012           = TestCase('contadj_naca0012')
    contadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012.test_iter = 5
    contadj_naca0012.test_vals = [-9.547544, -15.069350, -0.726250, 0.020280]
    contadj_naca0012.test_vals_aarch64 = [-9.662546, -14.998818, -0.726250, 0.020280]
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.085026, -12.651321, -1.086100, 0.007556]
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10
    contadj_wedge.test_vals = [2.872065, -2.756214, 1010800.000000, 0.000000]
    test_list.append(contadj_wedge)

    # Inviscid fixed CL NACA0012
    contadj_fixed_CL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixed_CL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixed_CL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixed_CL_naca0012.test_iter = 100
    contadj_fixed_CL_naca0012.test_vals = [1.408903, -4.010709, 0.005808, 0.004028]
    test_list.append(contadj_fixed_CL_naca0012)

    ###################################
    ### Cont. adj. compressible N-S ###
    ###################################

    # Adjoint laminar cylinder
    contadj_ns_cylinder           = TestCase('contadj_ns_cylinder')
    contadj_ns_cylinder.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder.test_iter = 20
    contadj_ns_cylinder.test_vals = [-3.626865, -9.079504, 2.056700, -0.000000]
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
    contadj_rans_rae2822.test_vals = [-5.399270, -10.904277, -0.212470, 0.005448]
    test_list.append(contadj_rans_rae2822)

    #############################
    ### Compressibele RANS UQ ###
    #############################

    # NACA0012 1c
    turb_naca0012_1c           = TestCase('turb_naca0012_1c')
    turb_naca0012_1c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_1c.cfg_file  = "turb_NACA0012_uq_1c.cfg"
    turb_naca0012_1c.test_iter = 10
    turb_naca0012_1c.test_vals = [-4.979554, 1.344291, 0.666805, 0.010247]
    turb_naca0012_1c.test_vals_aarch64 = [-4.981036, 1.345868, 0.673232, 0.010091]
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals = [-5.482691, 1.262918, 0.529442, -0.023584]
    turb_naca0012_2c.test_vals_aarch64 = [-5.484365, 1.264701, 0.501741, -0.033109]
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals = [-5.583617, 1.229594, 0.463384, -0.035444]
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals = [-5.129555, 1.283892, 0.808074, 0.047784]
    turb_naca0012_p1c1.test_vals_aarch64 = [-5.122100, 1.284478, 0.608744, -0.008593]
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals = [-5.554045, 1.234446, 0.603731, -0.008717]
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
    hb_rans_preconditioning.test_vals = [-1.905509, 0.482119, 0.598486, 3.605120, -5.946067]
    test_list.append(hb_rans_preconditioning)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Rotating NACA 0012
    rot_naca0012           = TestCase('rot_naca0012')
    rot_naca0012.cfg_dir   = "rotating/naca0012"
    rot_naca0012.cfg_file  = "rot_NACA0012.cfg"
    rot_naca0012.test_iter = 25
    rot_naca0012.test_vals = [-0.023360, 0.082942, 2.000000, -1.182142]
    test_list.append(rot_naca0012)

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-7.770113, -2.308185, 0.011072, 0.008947]
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.471322, -2.007082, 2.183565, 1.847848]
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-1.175979, 0.062203, 1.399352, 2.219226, 1.399298, 2.217472, 0.000000]
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Aeroacoustics postprocessing of a converged square-cylinder pressure history
    aeroacoustics_square_cylinder = TestCase('aeroacoustics_square_cylinder')
    aeroacoustics_square_cylinder.cfg_dir = "unsteady/square_cylinder"
    aeroacoustics_square_cylinder.cfg_file = "aeroacoustics.cfg"
    aeroacoustics_square_cylinder.test_iter = 0
    aeroacoustics_square_cylinder.test_vals = [103.1512634765, 1.9531250000]
    aeroacoustics_square_cylinder.command = TestCase.Command(
        exec = "python", param = "run_aeroacoustics.py"
    )
    aeroacoustics_square_cylinder.timeout = 60
    aeroacoustics_square_cylinder.tol = 1e-6
    test_list.append(aeroacoustics_square_cylinder)

    # Gust
    sine_gust           = TestCase('sine_gust')
    sine_gust.cfg_dir   = "gust"
    sine_gust.cfg_file  = "inv_gust_NACA0012.cfg"
    sine_gust.test_iter = 5
    sine_gust.test_vals = [-1.977498, 3.481817, -0.011360, -0.005241]
    sine_gust.unsteady  = True
    test_list.append(sine_gust)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714713, -5.763284, -0.214960, 0.023758, 0.000000]
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)

    # Stochastic BackScatter (SBS) model
    sbs_backward_step            = TestCase('sbs_backward_step')
    sbs_backward_step.cfg_dir    = "backscatter/backward_step"
    sbs_backward_step.cfg_file   = "backwardStep.cfg"
    sbs_backward_step.test_iter  = 3
    sbs_backward_step.test_vals = [-6.352885, -3.465372, -5.507901, -3.906544, -9.506305, -6.365236, -6.331021, -6.331028]
    sbs_backward_step.unsteady   = True
    sbs_backward_step.decompress = True
    sbs_backward_step.grid_file  = "backward_step.su2"
    test_list.append(sbs_backward_step)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.004011, -6.876257, 1.487888, 0.421869]
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # Flat plate
    flatplate_unsteady           = TestCase('flatplate_unsteady')
    flatplate_unsteady.cfg_dir   = "navierstokes/flatplate"
    flatplate_unsteady.cfg_file  = "lam_flatplate_unst.cfg"
    flatplate_unsteady.test_iter = 3
    flatplate_unsteady.test_vals = [0.000008, -8.875128, -8.250204, -6.305788, -5.469452, -3.398230, 0.002075, -0.325535]
    flatplate_unsteady.unsteady  = True
    test_list.append(flatplate_unsteady)

    ######################################
    ### NICFD                          ###
    ######################################

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 20
    edge_VW.test_vals = [-4.478306, 1.723006, -0.000009, 0.000000]
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 20
    edge_PPR.test_vals = [-8.395110, -2.236894, -0.000034, 0.000000]
    test_list.append(edge_PPR)

    # Rarefaction Q1D nozzle, include CoolProp fluid model
    coolprop_fluidModel           = TestCase('coolprop_fluidModel')
    coolprop_fluidModel.cfg_dir   = "nicf/coolprop"
    coolprop_fluidModel.cfg_file  = "fluidModel.cfg"
    coolprop_fluidModel.test_iter = 5
    coolprop_fluidModel.test_vals = [-4.666146, -0.172864, 4.492817, 0.000000, 0.000000]
    coolprop_fluidModel.enabled_on_cpu_arch = ["x86_64"]
    test_list.append(coolprop_fluidModel)

    # Rarefaction Q1D nozzle, include CoolProp transport model
    coolprop_transportModel           = TestCase('coolprop_transportModel')
    coolprop_transportModel.cfg_dir   = "nicf/coolprop"
    coolprop_transportModel.cfg_file  = "transportModel.cfg"
    coolprop_transportModel.test_iter = 5
    coolprop_transportModel.test_vals = [-4.667122, -0.172865, 4.687652, 0.000000, 0.000000]
    coolprop_transportModel.enabled_on_cpu_arch = ["x86_64"]
    test_list.append(coolprop_transportModel)

    # Rarefaction Q1D nozzle, include data-driven fluid model
    datadriven_fluidModel           = TestCase('datadriven_fluidModel')
    datadriven_fluidModel.cfg_dir   = "nicf/datadriven"
    datadriven_fluidModel.cfg_file  = "datadriven_nozzle.cfg"
    datadriven_fluidModel.test_iter = 50
    datadriven_fluidModel.test_vals = [-6.283594, -3.253961, -4.172935, -0.949415, -2.046039, 1.546933]
    test_list.append(datadriven_fluidModel)

    ######################################
    ### Turbomachinery                 ###
    ######################################

    # Aachen Turbine restart
    Aachen_3D_restart = TestCase('aachen_turbine_restart')
    Aachen_3D_restart.cfg_dir = "turbomachinery/Aachen_turbine"
    Aachen_3D_restart.cfg_file = "aachen_3D_MP_restart.cfg"
    Aachen_3D_restart.test_iter = 5
    Aachen_3D_restart.tol = 0.00001
    Aachen_3D_restart.test_vals = [-7.688401, -8.466746, -6.034652, -6.427491, -5.818347, -4.597687, -5.530255, -5.303027, -3.819060, -5.242491, -5.755143, -3.629921, -2.213103, -2.859285, -0.576337]
    test_list.append(Aachen_3D_restart)

    # Jones APU Turbocharger restart
    Jones_tc_restart           = TestCase('jones_turbocharger_restart')
    Jones_tc_restart.cfg_dir   = "turbomachinery/APU_turbocharger"
    Jones_tc_restart.cfg_file  = "Jones_restart.cfg"
    Jones_tc_restart.test_iter = 5
    Jones_tc_restart.test_vals = [-11.903628, -12.211875, -19.212536, -13.508088, -19.079442, -13.440735, 73286.000000, 73286.000000, 0.020056, 82.286000]
