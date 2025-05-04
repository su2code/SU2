#!/usr/bin/env python

## \file hybrid_regression_AD.py
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

    args = parse_args('Hybrid Regression AD Tests')

    test_list = []

    #####################################
    ### Disc. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    discadj_naca0012           = TestCase('discadj_naca0012')
    discadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    discadj_naca0012.cfg_file  = "inv_NACA0012_discadj.cfg"
    discadj_naca0012.test_iter = 100
    discadj_naca0012.test_vals = [-3.562611, -8.932638, -0.000000, 0.005608]
    test_list.append(discadj_naca0012)

    # Inviscid Cylinder 3D (multiple markers)
    discadj_cylinder3D           = TestCase('discadj_cylinder3D')
    discadj_cylinder3D.cfg_dir   = "disc_adj_euler/cylinder3D"
    discadj_cylinder3D.cfg_file  = "inv_cylinder3D.cfg"
    discadj_cylinder3D.test_iter = 5
    discadj_cylinder3D.test_vals = [-3.689810, -3.883743, -0.000000, 0.000000]
    test_list.append(discadj_cylinder3D)

    # Arina nozzle 2D
    discadj_arina2k              = TestCase('discadj_arina2k')
    discadj_arina2k.cfg_dir      = "disc_adj_euler/arina2k"
    discadj_arina2k.cfg_file     = "Arina2KRS.cfg"
    discadj_arina2k.test_iter    = 20
    discadj_arina2k.test_vals    = [-3.254503, -3.495599, 0.052373, 0.000000]
    test_list.append(discadj_arina2k)

    ####################################
    ### Disc. adj. compressible RANS ###
    ####################################

    # Adjoint turbulent NACA0012 SA
    discadj_rans_naca0012_sa           = TestCase('discadj_rans_naca0012_sa')
    discadj_rans_naca0012_sa.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    discadj_rans_naca0012_sa.test_iter = 10
    discadj_rans_naca0012_sa.test_vals = [-2.997064, -0.196172, 0.000003, -0.000000, 5.000000, -2.919675, 5.000000, -7.323218]
    test_list.append(discadj_rans_naca0012_sa)

    # Adjoint turbulent NACA0012 SST
    discadj_rans_naca0012_sst           = TestCase('discadj_rans_naca0012_sst')
    discadj_rans_naca0012_sst.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    discadj_rans_naca0012_sst.test_iter = 10
    discadj_rans_naca0012_sst.test_vals = [-2.237413, -0.200125, 2.763200, -0.039612]
    test_list.append(discadj_rans_naca0012_sst)

    #######################################
    ### Disc. adj. incompressible Euler ###
    #######################################

    # Adjoint Incompressible Inviscid NACA0012
    discadj_incomp_NACA0012           = TestCase('discadj_incomp_NACA0012')
    discadj_incomp_NACA0012.cfg_dir   = "disc_adj_incomp_euler/naca0012"
    discadj_incomp_NACA0012.cfg_file  = "incomp_NACA0012_disc.cfg"
    discadj_incomp_NACA0012.test_iter = 20
    discadj_incomp_NACA0012.test_vals = [20.000000, -4.091640, -2.655563, 0.000000]
    test_list.append(discadj_incomp_NACA0012)

    #####################################
    ### Disc. adj. incompressible N-S ###
    #####################################

    # Adjoint Incompressible Viscous Cylinder (Heated)
    discadj_incomp_cylinder           = TestCase('discadj_incomp_cylinder')
    discadj_incomp_cylinder.cfg_dir   = "disc_adj_incomp_navierstokes/cylinder"
    discadj_incomp_cylinder.cfg_file  = "heated_cylinder.cfg"
    discadj_incomp_cylinder.test_iter = 20
    discadj_incomp_cylinder.test_vals         = [20.000000, -2.746353, -2.934792, 0.000000]
    discadj_incomp_cylinder.test_vals_aarch64 = [20.000000, -2.746353, -2.934792, 0.000000]
    test_list.append(discadj_incomp_cylinder)

    ######################################
    ### Disc. adj. incompressible RANS ###
    ######################################

    # Adjoint Incompressible Turbulent NACA 0012 SA
    discadj_incomp_turb_NACA0012_sa           = TestCase('discadj_incomp_turb_NACA0012_sa')
    discadj_incomp_turb_NACA0012_sa.cfg_dir   = "disc_adj_incomp_rans/naca0012"
    discadj_incomp_turb_NACA0012_sa.cfg_file  = "turb_naca0012_sa.cfg"
    discadj_incomp_turb_NACA0012_sa.test_iter = 10
    discadj_incomp_turb_NACA0012_sa.test_vals = [10.000000, -3.845995, -1.031096, 0.000000]
    test_list.append(discadj_incomp_turb_NACA0012_sa)

    # Adjoint Incompressible Turbulent NACA 0012 SST
    discadj_incomp_turb_NACA0012_sst           = TestCase('discadj_incomp_turb_NACA0012_sst')
    discadj_incomp_turb_NACA0012_sst.cfg_dir   = "disc_adj_incomp_rans/naca0012"
    discadj_incomp_turb_NACA0012_sst.cfg_file  = "turb_naca0012_sst.cfg"
    discadj_incomp_turb_NACA0012_sst.test_iter = 10
    discadj_incomp_turb_NACA0012_sst.test_vals = [-3.597708, -2.983823, -8.354806, 0.000000, -0.916525]
    test_list.append(discadj_incomp_turb_NACA0012_sst)

    #######################################################
    ### Unsteady Disc. adj. compressible RANS           ###
    #######################################################

    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder.cfg"
    discadj_cylinder.test_iter = 9
    discadj_cylinder.test_vals = [1.639345, -2.834279, -0.009538, 0.000020]
    discadj_cylinder.unsteady  = True
    discadj_cylinder.enabled_with_tsan = False
    test_list.append(discadj_cylinder)

    ##############################################################
    ### Unsteady Disc. adj. compressible RANS Windowed Average ###
    ##############################################################

    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder_windowed_average_AD')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder_Windowing_AD.cfg"
    discadj_cylinder.test_iter = 9
    discadj_cylinder.test_vals = [2.183376]
    discadj_cylinder.unsteady  = True
    discadj_cylinder.enabled_with_tsan = False
    test_list.append(discadj_cylinder)

    ##########################################################################
    ### Unsteady Disc. adj. compressible RANS DualTimeStepping 1st order   ###
    ##########################################################################

    # Turbulent Cylinder
    discadj_DT_1ST_cylinder           = TestCase('unsteady_cylinder_DT_1ST')
    discadj_DT_1ST_cylinder.cfg_dir   = "disc_adj_rans/cylinder_DT_1ST"
    discadj_DT_1ST_cylinder.cfg_file  = "cylinder.cfg"
    discadj_DT_1ST_cylinder.test_iter = 9
    discadj_DT_1ST_cylinder.test_vals = [1.196346, -3.339016, -0.006212, 0.000020]
    discadj_DT_1ST_cylinder.unsteady  = True
    discadj_DT_1ST_cylinder.enabled_with_tsan = False
    test_list.append(discadj_DT_1ST_cylinder)

    ######################################################
    ### Unsteady Disc. adj. compressible pitching NACA ###
    ######################################################

    # compressible pitching NACA0012
    discadj_pitchingNACA0012           = TestCase('pitchingNACA0012')
    discadj_pitchingNACA0012.cfg_dir   = "disc_adj_euler/naca0012_pitching"
    discadj_pitchingNACA0012.cfg_file  = "inv_NACA0012_pitching.cfg"
    discadj_pitchingNACA0012.test_iter = 4
    discadj_pitchingNACA0012.test_vals = [-1.220333, -1.646832, -0.007539, 0.000013]
    discadj_pitchingNACA0012.unsteady  = True
    discadj_pitchingNACA0012.enabled_with_tsan = False
    test_list.append(discadj_pitchingNACA0012)

    #######################################################
    ### Disc. adj. turbomachinery                       ###
    #######################################################

    # Transonic Stator 2D
    discadj_trans_stator           = TestCase('transonic_stator')
    discadj_trans_stator.cfg_dir   = "disc_adj_turbomachinery/transonic_stator_2D"
    discadj_trans_stator.cfg_file  = "transonic_stator.cfg"
    discadj_trans_stator.test_iter = 79
    discadj_trans_stator.test_vals         = [79.000000, 0.770295, 0.383672, 0.472433, -0.996122, 2.153513, -4.444080]
    discadj_trans_stator.test_vals_aarch64 = [79, 0.769987, 0.383135, 0.472391, -0.996504, 2.153296, -4.444301]
    discadj_trans_stator.enabled_with_tsan = False
    test_list.append(discadj_trans_stator)

    ###################################
    ### Structural Adjoint          ###
    ###################################

    # Structural model
    discadj_fea           = TestCase('discadj_fea')
    discadj_fea.cfg_dir   = "disc_adj_fea"
    discadj_fea.cfg_file  = "configAD_fem.cfg"
    discadj_fea.test_iter = 4
    discadj_fea.test_vals         = [1.774569, 1.928023, -0.000364, -8.690300]
    discadj_fea.test_vals_aarch64 = [1.939275, 1.989717, -0.000364, -8.708200]
    test_list.append(discadj_fea)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    for test in test_list:
        test.command = TestCase.Command(exec = "SU2_CFD_AD", param = "-t 2")
        test.timeout = 600
        test.tol = 1e-4
    #end

    pass_list = [ test.run_test(args.tsan) for test in test_list ]

    ###################################
    ### Python Wrapper              ###
    ###################################

    # FEA AD Flow Load Sensitivity
    pywrapper_FEA_AD_FlowLoad               = TestCase('pywrapper_FEA_AD_FlowLoad')
    pywrapper_FEA_AD_FlowLoad.cfg_dir       = "py_wrapper/disc_adj_fea/flow_load_sens"
    pywrapper_FEA_AD_FlowLoad.cfg_file      = "configAD_fem.cfg"
    pywrapper_FEA_AD_FlowLoad.test_iter     = 100
    pywrapper_FEA_AD_FlowLoad.test_vals     = [-0.132037, -0.554472, -0.000364, -0.003101]
    pywrapper_FEA_AD_FlowLoad.test_vals_aarch64 = [-0.131745, -0.553214, -0.000364, -0.003101]
    pywrapper_FEA_AD_FlowLoad.command       = TestCase.Command(exec = "python", param = "run_adjoint.py --parallel -f")
    pywrapper_FEA_AD_FlowLoad.timeout       = 1600
    pywrapper_FEA_AD_FlowLoad.tol           = 1e-2
    pywrapper_FEA_AD_FlowLoad.new_output    = False
    pywrapper_FEA_AD_FlowLoad.enabled_with_tsan = False
    test_list.append(pywrapper_FEA_AD_FlowLoad)
    pass_list.append(pywrapper_FEA_AD_FlowLoad.run_test(args.tsan))

    # Flow AD Mesh Displacement Sensitivity
    pywrapper_CFD_AD_MeshDisp               = TestCase('pywrapper_CFD_AD_MeshDisp')
    pywrapper_CFD_AD_MeshDisp.cfg_dir       = "py_wrapper/disc_adj_flow/mesh_disp_sens"
    pywrapper_CFD_AD_MeshDisp.cfg_file      = "configAD_flow.cfg"
    pywrapper_CFD_AD_MeshDisp.test_iter     = 1000
    pywrapper_CFD_AD_MeshDisp.test_vals     = [30.000000, -2.496218, 1.441643, 0.000000]
    pywrapper_CFD_AD_MeshDisp.test_vals_aarch64 = [30.000000, -2.499079, 1.440068, 0.000000]
    pywrapper_CFD_AD_MeshDisp.command       = TestCase.Command(exec = "python", param = "run_adjoint.py --parallel -f")
    pywrapper_CFD_AD_MeshDisp.timeout       = 1600
    pywrapper_CFD_AD_MeshDisp.tol           = 1e-2
    pywrapper_CFD_AD_MeshDisp.new_output    = False
    pywrapper_CFD_AD_MeshDisp.enabled_with_tsan = False
    test_list.append(pywrapper_CFD_AD_MeshDisp)
    pass_list.append(pywrapper_CFD_AD_MeshDisp.run_test(args.tsan))


    # Tests summary
    print('==================================================================')
    print('Summary of the hybrid parallel AD tests')
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
