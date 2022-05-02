#!/usr/bin/env python

## \file serial_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.3.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

    #####################################
    ### Disc. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    discadj_naca0012           = TestCase('discadj_naca0012')
    discadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    discadj_naca0012.cfg_file  = "inv_NACA0012_discadj.cfg"
    discadj_naca0012.test_iter = 100
    discadj_naca0012.test_vals = [-3.561506, -8.926634, -0.000000, 0.005587]
    discadj_naca0012.su2_exec  = "SU2_CFD_AD"
    discadj_naca0012.timeout   = 1600
    discadj_naca0012.tol       = 0.00001
    test_list.append(discadj_naca0012)

    # Inviscid Cylinder 3D (multiple markers)
    discadj_cylinder3D           = TestCase('discadj_cylinder3D')
    discadj_cylinder3D.cfg_dir   = "disc_adj_euler/cylinder3D"
    discadj_cylinder3D.cfg_file  = "inv_cylinder3D.cfg"
    discadj_cylinder3D.test_iter = 5
    discadj_cylinder3D.test_vals = [-3.737675, -3.842311, -0.000000, 0.000000]
    discadj_cylinder3D.su2_exec  = "SU2_CFD_AD"
    discadj_cylinder3D.timeout   = 1600
    discadj_cylinder3D.tol       = 0.00001
    test_list.append(discadj_cylinder3D)

    # Arina nozzle 2D  
    discadj_arina2k              = TestCase('discadj_arina2k')
    discadj_arina2k.cfg_dir      = "disc_adj_euler/arina2k"
    discadj_arina2k.cfg_file     = "Arina2KRS.cfg"
    discadj_arina2k.test_iter    = 20
    discadj_arina2k.test_vals    = [-3.087863, -3.481496, 6.8879e-02, 0]
    discadj_arina2k.su2_exec     = "SU2_CFD_AD"
    discadj_arina2k.timeout      = 1600
    discadj_arina2k.tol          = 0.00001
    test_list.append(discadj_arina2k)

    #######################################################
    ### Disc. adj. compressible RANS                    ###
    #######################################################
    
    # Adjoint turbulent NACA0012 SA
    discadj_rans_naca0012_sa           = TestCase('discadj_rans_naca0012_sa')
    discadj_rans_naca0012_sa.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    discadj_rans_naca0012_sa.test_iter = 10
    discadj_rans_naca0012_sa.test_vals = [-2.230555, 0.645023, 0.180740, -0.000018, 5.000000, -4.275184, 5.000000, -8.892454] #last 8 columns
    discadj_rans_naca0012_sa.su2_exec  = "SU2_CFD_AD"
    discadj_rans_naca0012_sa.timeout   = 1600
    discadj_rans_naca0012_sa.tol       = 0.00001
    test_list.append(discadj_rans_naca0012_sa)

    # Adjoint turbulent NACA0012 SST
    discadj_rans_naca0012_sst           = TestCase('discadj_rans_naca0012_sst')
    discadj_rans_naca0012_sst.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    discadj_rans_naca0012_sst.test_iter = 10
    discadj_rans_naca0012_sst.test_vals = [-2.221792, -0.491367, 0.182010, -0.000018]
    discadj_rans_naca0012_sst.su2_exec  = "SU2_CFD_AD"
    discadj_rans_naca0012_sst.timeout   = 1600
    discadj_rans_naca0012_sst.tol       = 0.00001
    test_list.append(discadj_rans_naca0012_sst)

    #######################################
    ### Disc. adj. incompressible Euler ###
    #######################################

    # Adjoint Incompressible Inviscid NACA0012
    discadj_incomp_NACA0012           = TestCase('discadj_incomp_NACA0012')
    discadj_incomp_NACA0012.cfg_dir   = "disc_adj_incomp_euler/naca0012"
    discadj_incomp_NACA0012.cfg_file  = "incomp_NACA0012_disc.cfg"
    discadj_incomp_NACA0012.test_iter = 20
    discadj_incomp_NACA0012.test_vals = [20.000000, -4.092007, -2.652750, 0.000000]
    discadj_incomp_NACA0012.su2_exec  = "SU2_CFD_AD"
    discadj_incomp_NACA0012.timeout   = 1600
    discadj_incomp_NACA0012.tol       = 0.00001
    test_list.append(discadj_incomp_NACA0012)

    #####################################
    ### Disc. adj. incompressible N-S ###
    #####################################

    # Adjoint Incompressible Viscous Cylinder (Heated)
    discadj_incomp_cylinder           = TestCase('discadj_incomp_cylinder')
    discadj_incomp_cylinder.cfg_dir   = "disc_adj_incomp_navierstokes/cylinder"
    discadj_incomp_cylinder.cfg_file  = "heated_cylinder.cfg"
    discadj_incomp_cylinder.test_iter = 20
    discadj_incomp_cylinder.test_vals = [20.000000, -2.373367, -2.368305, 0.000000] #last 4 columns
    discadj_incomp_cylinder.su2_exec  = "SU2_CFD_AD"
    discadj_incomp_cylinder.timeout   = 1600
    discadj_incomp_cylinder.tol       = 0.00001
    test_list.append(discadj_incomp_cylinder)

    ######################################
    ### Disc. adj. incompressible RANS ###
    ######################################

    # Adjoint Incompressible Turbulent NACA 0012 SA
    discadj_incomp_turb_NACA0012_sa           = TestCase('discadj_incomp_turb_NACA0012_sa')
    discadj_incomp_turb_NACA0012_sa.cfg_dir   = "disc_adj_incomp_rans/naca0012"
    discadj_incomp_turb_NACA0012_sa.cfg_file  = "turb_naca0012_sa.cfg"
    discadj_incomp_turb_NACA0012_sa.test_iter = 10
    discadj_incomp_turb_NACA0012_sa.test_vals = [10.000000, -3.845979, -1.031095, 0.000000] #last 4 columns
    discadj_incomp_turb_NACA0012_sa.su2_exec  = "SU2_CFD_AD"
    discadj_incomp_turb_NACA0012_sa.timeout   = 1600
    discadj_incomp_turb_NACA0012_sa.tol       = 0.00001
    test_list.append(discadj_incomp_turb_NACA0012_sa)

    # Adjoint Incompressible Turbulent NACA 0012 SST
    discadj_incomp_turb_NACA0012_sst           = TestCase('discadj_incomp_turb_NACA0012_sst')
    discadj_incomp_turb_NACA0012_sst.cfg_dir   = "disc_adj_incomp_rans/naca0012"
    discadj_incomp_turb_NACA0012_sst.cfg_file  = "turb_naca0012_sst.cfg"
    discadj_incomp_turb_NACA0012_sst.test_iter = 10
    discadj_incomp_turb_NACA0012_sst.test_vals = [-3.845577, -2.414021, -8.420183, 0.000000] #last 4 columns
    discadj_incomp_turb_NACA0012_sst.su2_exec  = "SU2_CFD_AD"
    discadj_incomp_turb_NACA0012_sst.timeout   = 1600
    discadj_incomp_turb_NACA0012_sst.tol       = 0.00001
    test_list.append(discadj_incomp_turb_NACA0012_sst)

    #######################################################
    ### Unsteady Disc. adj. compressible RANS           ###
    #######################################################
   
    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder.cfg" 
    discadj_cylinder.test_iter = 9
    discadj_cylinder.test_vals = [3.746909, -1.544883, -0.008321, 0.000014] #last 4 columns
    discadj_cylinder.su2_exec  = "SU2_CFD_AD"
    discadj_cylinder.timeout   = 1600
    discadj_cylinder.tol       = 0.00001
    discadj_cylinder.unsteady  = True
    test_list.append(discadj_cylinder)
    
    ##########################################################################
    ### Unsteady Disc. adj. compressible RANS DualTimeStepping 1st order   ###
    ##########################################################################

    # Turbulent Cylinder
    discadj_DT_1ST_cylinder           = TestCase('unsteady_cylinder_DT_1ST')
    discadj_DT_1ST_cylinder.cfg_dir   = "disc_adj_rans/cylinder_DT_1ST"
    discadj_DT_1ST_cylinder.cfg_file  = "cylinder.cfg"
    discadj_DT_1ST_cylinder.test_iter = 9
    discadj_DT_1ST_cylinder.test_vals = [3.698168, -1.607050, -0.002159, 0.000028] #last 4 columns
    discadj_DT_1ST_cylinder.su2_exec  = "SU2_CFD_AD"
    discadj_DT_1ST_cylinder.timeout   = 1600
    discadj_DT_1ST_cylinder.tol       = 0.00001
    discadj_DT_1ST_cylinder.unsteady  = True
    test_list.append(discadj_DT_1ST_cylinder)

    ######################################################
    ### Unsteady Disc. adj. compressible pitching NACA ###
    ######################################################

    # compressible pitching NACA0012
    discadj_pitchingNACA0012           = TestCase('pitchingNACA0012')
    discadj_pitchingNACA0012.cfg_dir   = "disc_adj_euler/naca0012_pitching"
    discadj_pitchingNACA0012.cfg_file  = "inv_NACA0012_pitching.cfg"
    discadj_pitchingNACA0012.test_iter = 4
    discadj_pitchingNACA0012.test_vals = [-1.218846, -1.645199, -0.007645, 0.000013]
    discadj_pitchingNACA0012.su2_exec  = "SU2_CFD_AD"
    discadj_pitchingNACA0012.timeout   = 1600
    discadj_pitchingNACA0012.tol       = 0.00001
    discadj_pitchingNACA0012.unsteady  = True
    test_list.append(discadj_pitchingNACA0012)

    # deforming pitching NACA0012
    unst_deforming_naca0012           = TestCase('unst_deforming_naca0012')
    unst_deforming_naca0012.cfg_dir   = "disc_adj_euler/naca0012_pitching_def"
    unst_deforming_naca0012.cfg_file  = "inv_NACA0012_pitching_deform_ad.cfg"
    unst_deforming_naca0012.test_iter = 4
    unst_deforming_naca0012.test_vals = [-1.958006, -1.841808, 1081.700000, 0.000004]
    unst_deforming_naca0012.su2_exec  = "SU2_CFD_AD"
    unst_deforming_naca0012.timeout   = 1600
    unst_deforming_naca0012.tol       = 0.00001
    unst_deforming_naca0012.unsteady  = True
    test_list.append(unst_deforming_naca0012)

    ###################################
    ### Structural Adjoint          ###
    ###################################
   
    # Structural model
    discadj_fea           = TestCase('discadj_fea')
    discadj_fea.cfg_dir   = "disc_adj_fea"
    discadj_fea.cfg_file  = "configAD_fem.cfg" 
    discadj_fea.test_iter = 4
    discadj_fea.test_vals = [-2.849531, -3.238474, -3.6413e-04, -8.7087] #last 4 columns
    discadj_fea.su2_exec  = "SU2_CFD_AD"
    discadj_fea.timeout   = 1600
    discadj_fea.tol       = 0.00001
    test_list.append(discadj_fea)    

    ###################################
    ### Disc. adj. heat             ###
    ###################################

    # Discrete adjoint for heated cylinder
    discadj_heat           = TestCase('discadj_heat')
    discadj_heat.cfg_dir   = "disc_adj_heat"
    discadj_heat.cfg_file  = "disc_adj_heat.cfg"
    discadj_heat.test_iter = 10
    discadj_heat.test_vals = [-2.271569, 0.671288, -3.172000, -8.231500] #last 4 columns
    discadj_heat.su2_exec  = "SU2_CFD_AD"
    discadj_heat.timeout   = 1600
    discadj_heat.tol       = 0.00001
    test_list.append(discadj_heat)

    ###################################
    ### Coupled FSI Adjoint         ###
    ###################################
   
    # Structural model
    discadj_fsi           = TestCase('discadj_fsi')
    discadj_fsi.cfg_dir   = "disc_adj_fsi"
    discadj_fsi.cfg_file  = "config.cfg"
    discadj_fsi.test_iter = 6
    discadj_fsi.test_vals = [6, -1.572702, -3.084381, 4.3990e-04, -1.0631] #last 5 columns
    discadj_fsi.su2_exec  = "SU2_CFD_AD"
    discadj_fsi.timeout   = 1600
    discadj_fsi.tol       = 0.00001
    test_list.append(discadj_fsi)

    ###################################
    ### Coupled CHT Adjoint         ###
    ###################################

    # Coupled discrete adjoint for heatflux in heated cylinder array
    discadj_cht           = TestCase('discadj_cht')
    discadj_cht.cfg_dir   = "coupled_cht/disc_adj_incomp_2d"
    discadj_cht.cfg_file  = "cht_2d_3cylinders.cfg"
    discadj_cht.test_iter = 10
    discadj_cht.test_vals = [-2.364405, -3.085551, -3.085518, -3.085513] #last 4 columns
    discadj_cht.su2_exec  = "SU2_CFD_AD"
    discadj_cht.timeout   = 1600
    discadj_cht.tol       = 0.00001
    test_list.append(discadj_cht)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    pass_list = [ test.run_test() for test in test_list ]
    
    ###################################
    ### Coupled RHT-CFD Adjoint     ###
    ###################################

    # Coupled discrete adjoint for radiative heat transfer in heated cylinder
    discadj_rht                = TestCase('discadj_rht')
    discadj_rht.cfg_dir        = "radiation/p1adjoint"
    discadj_rht.cfg_file       = "configp1adjoint.cfg"
    discadj_rht.test_iter      = 10
    discadj_rht.su2_exec       = "discrete_adjoint.py -f"
    discadj_rht.timeout        = 1600
    discadj_rht.reference_file = "of_grad_cd.csv.ref"
    discadj_rht.test_file      = "of_grad_cd.csv"
    pass_list.append(discadj_rht.run_filediff())
    test_list.append(discadj_rht)

    ######################################
    ### RUN PYTHON TESTS               ###
    ######################################
    
    # test discrete_adjoint.py
    discadj_euler_py = TestCase('discadj_euler_py')
    discadj_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    discadj_euler_py.cfg_file  = "inv_NACA0012.cfg"
    discadj_euler_py.test_iter = 10
    discadj_euler_py.su2_exec  = "discrete_adjoint.py -f"
    discadj_euler_py.timeout   = 1600
    discadj_euler_py.reference_file = "of_grad_cd_disc.dat.ref"
    discadj_euler_py.test_file = "of_grad_cd.dat"
    pass_list.append(discadj_euler_py.run_filediff())
    test_list.append(discadj_euler_py)
    
    # test discrete_adjoint with multiple ffd boxes
    discadj_multiple_ffd_py = TestCase('discadj_multiple_ffd_py')
    discadj_multiple_ffd_py.cfg_dir = "multiple_ffd/naca0012"
    discadj_multiple_ffd_py.cfg_file  = "inv_NACA0012_ffd.cfg"
    discadj_multiple_ffd_py.test_iter = 9
    discadj_multiple_ffd_py.su2_exec  = "discrete_adjoint.py -f"
    discadj_multiple_ffd_py.timeout   = 1600
    discadj_multiple_ffd_py.reference_file = "of_grad_cd.dat.ref"
    discadj_multiple_ffd_py.test_file = "of_grad_cd.dat"
    pass_list.append(discadj_multiple_ffd_py.run_filediff())
    test_list.append(discadj_multiple_ffd_py)

    # test direct_differentiation.py
    directdiff_euler_py = TestCase('directdiff_euler_py')
    directdiff_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    directdiff_euler_py.cfg_file  = "inv_NACA0012_FD.cfg"
    directdiff_euler_py.test_iter = 10
    directdiff_euler_py.su2_exec  = "direct_differentiation.py -f"
    directdiff_euler_py.timeout   = 1600
    directdiff_euler_py.reference_file = "of_grad_directdiff.dat.ref"
    directdiff_euler_py.test_file = "DIRECTDIFF/of_grad_directdiff.dat"
    pass_list.append(directdiff_euler_py.run_filediff())
    test_list.append(directdiff_euler_py)

    # test direct_differentiation.py with multiple ffd boxes
    directdiff_multiple_ffd_py = TestCase('directdiff_multiple_ffd_py')
    directdiff_multiple_ffd_py.cfg_dir = "multiple_ffd/naca0012"
    directdiff_multiple_ffd_py.cfg_file  = "inv_NACA0012_ffd.cfg"
    directdiff_multiple_ffd_py.test_iter = 9
    directdiff_multiple_ffd_py.su2_exec  = "direct_differentiation.py -f"
    directdiff_multiple_ffd_py.timeout   = 1600
    directdiff_multiple_ffd_py.reference_file = "of_grad_directdiff.dat.ref"
    directdiff_multiple_ffd_py.test_file = "DIRECTDIFF/of_grad_directdiff.dat"
    pass_list.append(directdiff_multiple_ffd_py.run_filediff())
    test_list.append(directdiff_multiple_ffd_py)

    # test continuous_adjoint.py, with multiple objectives
#    discadj_multi_py            = TestCase('discadj_multi_py')
#    discadj_multi_py.cfg_dir    = "cont_adj_euler/wedge"
#    discadj_multi_py.cfg_file   = "inv_wedge_ROE_multiobj.cfg"
#    discadj_multi_py.test_iter  = 10
#    discadj_multi_py.su2_exec   = "discrete_adjoint.py"
#    discadj_multi_py.timeout    = 1600
#    discadj_multi_py.reference_file = "of_grad_combo.dat.refdiscrete"
#    discadj_multi_py.test_file  = "of_grad_combo.dat"
#    pass_list.append(discadj_multi_py.run_filediff())
#    test_list.append(discadj_multi_py)

    # FEA AD Flow Load Sensitivity
    pywrapper_FEA_AD_FlowLoad               = TestCase('pywrapper_FEA_AD_FlowLoad')
    pywrapper_FEA_AD_FlowLoad.cfg_dir       = "py_wrapper/disc_adj_fea/flow_load_sens"
    pywrapper_FEA_AD_FlowLoad.cfg_file      = "configAD_fem.cfg"
    pywrapper_FEA_AD_FlowLoad.test_iter     = 100
    pywrapper_FEA_AD_FlowLoad.test_vals     = [-0.13945587401579657, -0.585985886606256, -0.00036377840086080753, -0.0031005670174756375] #last 4 columns
    pywrapper_FEA_AD_FlowLoad.su2_exec      = "python run_adjoint.py -f"
    pywrapper_FEA_AD_FlowLoad.timeout       = 1600
    pywrapper_FEA_AD_FlowLoad.tol           = 0.000001
    pywrapper_FEA_AD_FlowLoad.new_output    = False
    test_list.append(pywrapper_FEA_AD_FlowLoad)
    pass_list.append(pywrapper_FEA_AD_FlowLoad.run_test())

    # Flow AD Mesh Displacement Sensitivity
    pywrapper_FEA_AD_FlowLoad               = TestCase('pywrapper_CFD_AD_MeshDisp')
    pywrapper_FEA_AD_FlowLoad.cfg_dir       = "py_wrapper/disc_adj_flow/mesh_disp_sens"
    pywrapper_FEA_AD_FlowLoad.cfg_file      = "configAD_flow.cfg"
    pywrapper_FEA_AD_FlowLoad.test_iter     = 1000
    pywrapper_FEA_AD_FlowLoad.test_vals     = [30.000000, -2.518695, 1.390150, 0.000000] #last 4 columns
    pywrapper_FEA_AD_FlowLoad.su2_exec      = "python run_adjoint.py -f"
    pywrapper_FEA_AD_FlowLoad.timeout       = 1600
    pywrapper_FEA_AD_FlowLoad.tol           = 0.000001
    pywrapper_FEA_AD_FlowLoad.new_output    = False
    test_list.append(pywrapper_FEA_AD_FlowLoad)
    pass_list.append(pywrapper_FEA_AD_FlowLoad.run_test())


    ###################################
    ### Sobolev Gradient Smoothing  ###
    ###################################

    grad_smooth_naca0012           = TestCase('grad_smooth_naca0012')
    grad_smooth_naca0012.cfg_dir   = "grad_smooth/naca0012"
    grad_smooth_naca0012.cfg_file  = "inv_NACA0012_gradsmooth.cfg"
    grad_smooth_naca0012.test_iter = 1
    grad_smooth_naca0012.su2_exec  = "SU2_DOT_AD"
    grad_smooth_naca0012.timeout   = 1600
    grad_smooth_naca0012.tol       = 0.00001
    grad_smooth_naca0012.reference_file = "of_hess.dat.ref"
    grad_smooth_naca0012.test_file = "of_hess.dat"
    pass_list.append(grad_smooth_naca0012.run_filediff())
    test_list.append(grad_smooth_naca0012)

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
