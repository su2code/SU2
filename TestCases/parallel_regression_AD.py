#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.5.0 "Blackbird"
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

    # Arina nozzle 2D
    discadj_arina2k              = TestCase('discadj_arina2k')
    discadj_arina2k.cfg_dir      = "disc_adj_euler/arina2k"
    discadj_arina2k.cfg_file     = "Arina2KRS.cfg"
    discadj_arina2k.test_iter    = 20
    discadj_arina2k.test_vals    = [-3.111181, -3.501516, 6.8705e-02, 0]
    test_list.append(discadj_arina2k)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    # set suitable defaults unless something else has been specified
    # command: "mpirun -n 2 SU2_CFD_AD"
    # timeout: 1600
    # tol:     0.00001
    for test in test_list:
        if test.command.empty():
            test.command = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
        if test.timeout == 0:
            test.timeout = 1600
        if test.tol == 0.0:
            test.tol = 0.00001

    pass_list = [ test.run_test() for test in test_list ]

    ##################################################
    ### Structural Adjoint - Topology Optimization ###
    ##################################################

    # test discrete_adjoint.py
    discadj_topol_optim = TestCase('discadj_topol_optim')
    discadj_topol_optim.cfg_dir = "fea_topology"
    discadj_topol_optim.cfg_file  = "config.cfg"
    discadj_topol_optim.test_iter = 0
    discadj_topol_optim.command = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    discadj_topol_optim.timeout   = 1600
    discadj_topol_optim.reference_file         = "grad_ref_node.dat.ref"
    discadj_topol_optim.reference_file_aarch64 = "grad_ref_node_aarch64.dat.ref"
    discadj_topol_optim.test_file = "grad_ref_node.dat"
    pass_list.append(discadj_topol_optim.run_filediff())
    test_list.append(discadj_topol_optim)

    ####################################################################################
    ### Unsteady Disc. adj. compressible RANS Windowed Average with restart solution ###
    ####################################################################################

    # NACA0012 Airfoil
    unsteady_naca0012           = TestCase('unsteady_NACA0012_restart_adjoint')
    unsteady_naca0012.cfg_dir   = "disc_adj_rans/naca0012"
    unsteady_naca0012.cfg_file  = "naca0012.cfg"
    unsteady_naca0012.test_iter = 14
    unsteady_naca0012.command   = TestCase.Command(exec = "discrete_adjoint.py", param = "-f")
    unsteady_naca0012.timeout   = 1600
    unsteady_naca0012.reference_file = "of_grad_cd.csv.ref"
    unsteady_naca0012.test_file = "of_grad_cd.csv"
    unsteady_naca0012.unsteady  = True
    pass_list.append(unsteady_naca0012.run_filediff())
    test_list.append(unsteady_naca0012)

    ####################################################################################
    ### Unsteady Disc. adj. compressible RANS Windowed Average  only adjoint 		 ###
    ####################################################################################

    # NACA0012 Airfoil (Test depends on results of "unsteady_NACA0012_restart_adjoint")
    unsteady_naca0012           = TestCase('unsteady_NACA0012_adjoint_only')
    unsteady_naca0012.cfg_dir   = "disc_adj_rans/naca0012"
    unsteady_naca0012.cfg_file  = "naca0012.cfg"
    unsteady_naca0012.test_iter = 14
    unsteady_naca0012.command   = TestCase.Command(exec = "discrete_adjoint.py", param = "-m adj -f")
    unsteady_naca0012.timeout   = 1600
    unsteady_naca0012.reference_file = "of_grad_cd.csv.ref"
    unsteady_naca0012.test_file = "of_grad_cd.csv"
    unsteady_naca0012.unsteady  = True
    pass_list.append(unsteady_naca0012.run_filediff())
    test_list.append(unsteady_naca0012)

    ####################################################################
    ###  Unsteady Disc. adj. compressible RANS restart optimization  ###
    ####################################################################

    # test shape_optimization.py
    naca_restart_shape_opt      = TestCase('restart_shape_optimization')
    naca_restart_shape_opt.cfg_dir    = "optimization_rans/naca0012"
    naca_restart_shape_opt.cfg_file   = "naca0012.cfg"
    naca_restart_shape_opt.test_iter  = 1
    naca_restart_shape_opt.test_vals  = [1.000000, 1.000000, 0.007046, 0.196671]
    naca_restart_shape_opt.command    = TestCase.Command(exec = "shape_optimization.py", param = "-f")
    naca_restart_shape_opt.timeout    = 1600
    naca_restart_shape_opt.tol       = 0.00001
    pass_list.append(naca_restart_shape_opt.run_opt())
    test_list.append(naca_restart_shape_opt)

    ####################################################################
    ###  Unsteady Disc. Adj. Coupled FSI                             ###
    ####################################################################

    # Unsteady multi physics framework
    dyn_discadj_fsi           = TestCase('dyn_discadj_fsi')
    dyn_discadj_fsi.cfg_dir   = "disc_adj_fsi/dyn_fsi"
    dyn_discadj_fsi.cfg_file  = "config.cfg"
    dyn_discadj_fsi.test_iter = 2
    dyn_discadj_fsi.command = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    dyn_discadj_fsi.timeout   = 1600
    dyn_discadj_fsi.reference_file = "grad_dv.opt.ref"
    dyn_discadj_fsi.reference_file_aarch64 = "grad_dv_aarch64.opt.ref"
    dyn_discadj_fsi.test_file = "grad_young.opt"
    dyn_discadj_fsi.unsteady  = True
    pass_list.append(dyn_discadj_fsi.run_filediff())
    test_list.append(dyn_discadj_fsi)

    ####################################################################
    ###  Sobolev Gradient Smoothing                                  ###
    ####################################################################

    # Sobolev gradient smoothing solver
    grad_smooth_oneram6           = TestCase('grad_smooth_sob')
    grad_smooth_oneram6.cfg_dir   = "grad_smooth/oneram6"
    grad_smooth_oneram6.cfg_file  = "ONERAM6_gradsmooth.cfg"
    grad_smooth_oneram6.test_iter = 2
    grad_smooth_oneram6.command   = TestCase.Command("mpirun -n 2", "SU2_DOT_AD")
    grad_smooth_oneram6.timeout   = 1600
    grad_smooth_oneram6.reference_file = "of_hess.dat.ref"
    grad_smooth_oneram6.reference_file_aarch64 = "of_hess_aarch64.dat.ref"
    grad_smooth_oneram6.test_file = "of_hess.dat"
    pass_list.append(grad_smooth_oneram6.run_filediff())
    test_list.append(grad_smooth_oneram6)


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
