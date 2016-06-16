#!/usr/bin/env python

## \file serial_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 4.2.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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
    discadj_naca0012.test_vals = [-3.606839,-9.035212,-3.3386e-07,1.8777e-01] #last 4 columns
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
    discadj_rans_naca0012_sa.test_vals = [-1.751947, 0.485758, 0.181440, -0.385110] #last 4 columns
    discadj_rans_naca0012_sa.su2_exec  = "SU2_CFD_AD"
    discadj_rans_naca0012_sa.timeout   = 1600
    discadj_rans_naca0012_sa.tol       = 0.00001
    test_list.append(discadj_rans_naca0012_sa)

    # Adjoint turbulent NACA0012 SST
    discadj_rans_naca0012_sst           = TestCase('discadj_rans_naca0012_sst')
    discadj_rans_naca0012_sst.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    discadj_rans_naca0012_sst.test_iter = 10
    discadj_rans_naca0012_sst.test_vals = [-1.658566, -0.487694, 0.087556, -0.537010] #last 4 columns
    discadj_rans_naca0012_sst.su2_exec  = "SU2_CFD_AD"
    discadj_rans_naca0012_sst.timeout   = 1600
    discadj_rans_naca0012_sst.tol       = 0.00001
    test_list.append(discadj_rans_naca0012_sst)

    #######################################################
    ### Unsteady Disc. adj. compressible RANS           ###
    #######################################################
   
    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder.cfg" 
    discadj_cylinder.test_iter = 10
    discadj_cylinder.test_vals = [3.522068,-1.787841,-1.2030e-02,1.1156e-03] #last 4 columns
    discadj_cylinder.su2_exec  = "SU2_CFD_AD"
    discadj_cylinder.timeout   = 1600
    discadj_cylinder.tol       = 0.00001
    discadj_cylinder.unsteady  = True
    test_list.append(discadj_cylinder)

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
