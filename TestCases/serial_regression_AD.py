#!/usr/bin/env python

## \file serial_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 5.0.0 "Raven"
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
    discadj_naca0012.test_vals = [-3.606839, -9.035212, -0.000000, 0.005688] #last 4 columns
    discadj_naca0012.su2_exec  = "SU2_CFD_AD"
    discadj_naca0012.timeout   = 1600
    discadj_naca0012.tol       = 0.00001
    test_list.append(discadj_naca0012)

    #######################################################
    ### Disc. adj. compressible RANS                    ###
    #######################################################
    
    # Adjoint turbulent NACA0012 SA
    discadj_rans_naca0012_sa           = TestCase('discadj_rans_naca0012_sa')
    discadj_rans_naca0012_sa.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    discadj_rans_naca0012_sa.test_iter = 10
    discadj_rans_naca0012_sa.test_vals = [-1.751962, 0.485775, 0.182122, -0.000018] #last 4 columns
    discadj_rans_naca0012_sa.su2_exec  = "SU2_CFD_AD"
    discadj_rans_naca0012_sa.timeout   = 1600
    discadj_rans_naca0012_sa.tol       = 0.00001
    test_list.append(discadj_rans_naca0012_sa)

    # Adjoint turbulent NACA0012 SST
    discadj_rans_naca0012_sst           = TestCase('discadj_rans_naca0012_sst')
    discadj_rans_naca0012_sst.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    discadj_rans_naca0012_sst.test_iter = 10
    discadj_rans_naca0012_sst.test_vals = [-1.654903, -0.491416, 0.109157, 0.000011] #last 4 columns
    discadj_rans_naca0012_sst.su2_exec  = "SU2_CFD_AD"
    discadj_rans_naca0012_sst.timeout   = 1600
    discadj_rans_naca0012_sst.tol       = 0.00001
    test_list.append(discadj_rans_naca0012_sst)

    #######################################
    ### Disc. adj. incompressible Euler ###
    #######################################

    # Adjoint Incompressible Inviscid NACA0012
    discadj_incomp_NACA0012           = TestCase('discadj_incomp_NACA0012')
    discadj_incomp_NACA0012.cfg_dir   = "cont_adj_incomp_euler/naca0012"
    discadj_incomp_NACA0012.cfg_file  = "incomp_NACA0012_disc.cfg"
    discadj_incomp_NACA0012.test_iter = 20
    discadj_incomp_NACA0012.test_vals = [-2.911277, -2.705448, 0.000000, 0.000000] #last 4 columns
    discadj_incomp_NACA0012.su2_exec  = "SU2_CFD_AD"
    discadj_incomp_NACA0012.timeout   = 1600
    discadj_incomp_NACA0012.tol       = 0.00001
    test_list.append(discadj_incomp_NACA0012)

    #####################################
    ### Disc. adj. incompressible N-S ###
    #####################################

    # Adjoint Incompressible Viscous Cylinder
    discadj_incomp_cylinder           = TestCase('discadj_incomp_cylinder')
    discadj_incomp_cylinder.cfg_dir   = "cont_adj_incomp_navierstokes/cylinder"
    discadj_incomp_cylinder.cfg_file  = "lam_incomp_cylinder_disc.cfg"
    discadj_incomp_cylinder.test_iter = 20
    discadj_incomp_cylinder.test_vals = [-2.713325, -1.751175, 0.000000, 0.000000] #last 4 columns
    discadj_incomp_cylinder.su2_exec  = "SU2_CFD_AD"
    discadj_incomp_cylinder.timeout   = 1600
    discadj_incomp_cylinder.tol       = 0.00001
    test_list.append(discadj_incomp_cylinder)

    ######################################
    ### Disc. adj. incompressible RANS ###
    ######################################

    # Adjoint Incompressible Turbulent NACA 0012
    discadj_incomp_turb_NACA0012           = TestCase('discadj_incomp_turb_NACA0012')
    discadj_incomp_turb_NACA0012.cfg_dir   = "incomp_rans/naca0012"
    discadj_incomp_turb_NACA0012.cfg_file  = "naca0012_disc.cfg"
    discadj_incomp_turb_NACA0012.test_iter = 100
    discadj_incomp_turb_NACA0012.test_vals = [-3.627673, -1.624120, 0.000000, 0.000000] #last 4 columns
    discadj_incomp_turb_NACA0012.su2_exec  = "SU2_CFD_AD"
    discadj_incomp_turb_NACA0012.timeout   = 1600
    discadj_incomp_turb_NACA0012.tol       = 0.00001
    test_list.append(discadj_incomp_turb_NACA0012)

    #######################################################
    ### Unsteady Disc. adj. compressible RANS           ###
    #######################################################
   
    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder.cfg" 
    discadj_cylinder.test_iter = 9
    discadj_cylinder.test_vals = [3.746904, -1.544886, -0.008345, 0.000014] #last 4 columns
    discadj_cylinder.su2_exec  = "SU2_CFD_AD"
    discadj_cylinder.timeout   = 1600
    discadj_cylinder.tol       = 0.00001
    discadj_cylinder.unsteady  = True
    test_list.append(discadj_cylinder)
    
    ###################################
    ### Structural Adjoint          ###
    ###################################
   
    # Structural model
    discadj_fea           = TestCase('discadj_fea')
    discadj_fea.cfg_dir   = "disc_adj_fea"
    discadj_fea.cfg_file  = "configAD_fem.cfg" 
    discadj_fea.test_iter = 9
    discadj_fea.test_vals = [-6.403771, -6.404196, -3.6413e-04, -8.7087e+00] #last 4 columns
    discadj_fea.su2_exec  = "SU2_CFD_AD"
    discadj_fea.timeout   = 1600
    discadj_fea.tol       = 0.00001
    test_list.append(discadj_fea)    
    
    ###################################
    ### Coupled FSI Adjoint         ###
    ###################################
   
    # Structural model
    discadj_fsi           = TestCase('discadj_fsi')
    discadj_fsi.cfg_dir   = "disc_adj_fsi"
    discadj_fsi.cfg_file  = "configAD_fsi.cfg" 
    discadj_fsi.test_iter = 3000
    discadj_fsi.test_vals = [0.958848,-0.157183,0.658415,1.302076] #last 4 columns
    discadj_fsi.su2_exec  = "SU2_CFD_AD"
    discadj_fsi.timeout   = 1600
    discadj_fsi.tol       = 0.00001
    test_list.append(discadj_fsi)      

    ######################################
    ### RUN TESTS                      ###
    ######################################  

    pass_list = [ test.run_test() for test in test_list ]
    
    ######################################
    ### RUN PYTHON TESTS               ###
    ######################################
    
    # test discrete_adjoint.py
    discadj_euler_py = TestCase('discadj_euler_py')
    discadj_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    discadj_euler_py.cfg_file  = "inv_NACA0012.cfg"
    discadj_euler_py.test_iter = 10
    discadj_euler_py.su2_exec  = "discrete_adjoint.py"
    discadj_euler_py.timeout   = 1600
    discadj_euler_py.reference_file = "of_grad_cd_disc.dat.ref"
    discadj_euler_py.test_file = "of_grad_cd.dat"
    pass_list.append(discadj_euler_py.run_filediff())
    test_list.append(discadj_euler_py)
    
    # test direct_differentiation.py
    directdiff_euler_py = TestCase('directdiff_euler_py')
    directdiff_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    directdiff_euler_py.cfg_file  = "inv_NACA0012_FD.cfg"
    directdiff_euler_py.test_iter = 10
    directdiff_euler_py.su2_exec  = "direct_differentiation.py"
    directdiff_euler_py.timeout   = 1600
    directdiff_euler_py.reference_file = "of_grad_directdiff.dat.ref"
    directdiff_euler_py.test_file = "DIRECTDIFF/of_grad_directdiff.dat"
    pass_list.append(directdiff_euler_py.run_filediff())
    test_list.append(directdiff_euler_py)

    # test continuous_adjoint.py, with multiple objectives
    discadj_multi_py            = TestCase('discadj_multi_py')
    discadj_multi_py.cfg_dir    = "cont_adj_euler/wedge"
    discadj_multi_py.cfg_file   = "inv_wedge_ROE_multiobj.cfg"
    discadj_multi_py.test_iter  = 10
    discadj_multi_py.su2_exec   = "discrete_adjoint.py"
    discadj_multi_py.timeout    = 1600
    discadj_multi_py.reference_file = "of_grad_combo.dat.refdiscrete"
    discadj_multi_py.test_file  = "of_grad_combo.dat"
    pass_list.append(discadj_multi_py.run_filediff())
    test_list.append(discadj_multi_py)

    # Tests summary
    print '=================================================================='
    print 'Summary of the serial tests'
    for i, test in enumerate(test_list):
        if (pass_list[i]):
            print '  passed - %s'%test.tag
        else:
            print '* FAILED - %s'%test.tag
    
    if all(pass_list):
        sys.exit(0)
    else:
        sys.exit(1)
    # done

if __name__ == '__main__':
    main()
