#!/usr/bin/env python

## \file serial_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.0.4 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

from __future__ import print_function, division, absolute_import
import sys
from TestCase import TestCase

def main():
    '''This program runs SU2 and ensures that the output matches specified values.
       This will be used to do checks when code is pushed to github
       to make sure nothing is broken. '''

    test_list = []

    #################################
    ## Streamwise Periodic primal ###
    #################################

    # Laminar cylinder in channel, streamwise periodic
    streamwise_periodic_cylinder           = TestCase('streamwise_periodic_cylinder')
    streamwise_periodic_cylinder.cfg_dir   = "incomp_navierstokes/streamwise_periodic/half_cylinder_2D"
    streamwise_periodic_cylinder.cfg_file  = "half_cylinder_2D.cfg"
    streamwise_periodic_cylinder.test_iter = 30
    streamwise_periodic_cylinder.test_vals = [30, -7.818388, -6.797497, -6.968131] #last 4 lines
    streamwise_periodic_cylinder.su2_exec  = "parallel_computation.py -f"
    streamwise_periodic_cylinder.timeout   = 1600
    streamwise_periodic_cylinder.tol       = 0.00001
    test_list.append(streamwise_periodic_cylinder)

    # 3D laminar channnel with 1 cell in flow direction, streamwise periodic
    sp_pipeSlice_3d_dp_hf_tp           = TestCase('sp_pipeSlice_3d_dp_hf_tp')
    sp_pipeSlice_3d_dp_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/pipeSlice_3d"
    sp_pipeSlice_3d_dp_hf_tp.cfg_file  = "sp_pipeSlice_3d_dp_hf_tp.cfg"
    sp_pipeSlice_3d_dp_hf_tp.test_iter = 10
    sp_pipeSlice_3d_dp_hf_tp.test_vals = [10, -10.352122, -10.185237, -10.185237] #last 4 lines
    sp_pipeSlice_3d_dp_hf_tp.su2_exec  = "parallel_computation.py -f"
    sp_pipeSlice_3d_dp_hf_tp.timeout   = 1600
    sp_pipeSlice_3d_dp_hf_tp.tol       = 0.00001
    test_list.append(sp_pipeSlice_3d_dp_hf_tp)

    # 2D pin case pressure drop periodic with heatflux BC and temperature periodicity
    sp_pinArray_2d_dp_hf_tp           = TestCase('sp_pinArray_2d_dp_hf_tp')
    sp_pinArray_2d_dp_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/pinArray_2d"
    sp_pinArray_2d_dp_hf_tp.cfg_file  = "sp_pinArray_2d_dp_hf_tp.cfg"
    sp_pinArray_2d_dp_hf_tp.test_iter = 25
    sp_pinArray_2d_dp_hf_tp.test_vals = [-4.667133, 1.395801, -0.709306, 208.023676] #last 4 lines
    sp_pinArray_2d_dp_hf_tp.su2_exec  = "parallel_computation.py -f"
    sp_pinArray_2d_dp_hf_tp.timeout   = 1600
    sp_pinArray_2d_dp_hf_tp.tol       = 0.00001
    test_list.append(sp_pinArray_2d_dp_hf_tp)

    # 2D pin case massflow periodic with heatflux BC and prescribed heat
    sp_pinArray_2d_mf_hf           = TestCase('sp_pinArray_2d_mf_hf')
    sp_pinArray_2d_mf_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/pinArray_2d"
    sp_pinArray_2d_mf_hf.cfg_file  = "sp_pinArray_2d_mf_hf.cfg"
    sp_pinArray_2d_mf_hf.test_iter = 25
    sp_pinArray_2d_mf_hf.test_vals = [-4.666406, 1.398210, -0.710070, 208.677550] #last 4 lines
    sp_pinArray_2d_mf_hf.su2_exec  = "parallel_computation.py -f"
    sp_pinArray_2d_mf_hf.timeout   = 1600
    sp_pinArray_2d_mf_hf.tol       = 0.00001
    test_list.append(sp_pinArray_2d_mf_hf)

    # 2D CHT case with HF BC and
    sp_pinArray_cht_2d_mf_hf           = TestCase('sp_pinArray_cht_2d_mf_hf')
    sp_pinArray_cht_2d_mf_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    sp_pinArray_cht_2d_mf_hf.cfg_file  = "configMaster.cfg"
    sp_pinArray_cht_2d_mf_hf.test_iter = 100
    sp_pinArray_cht_2d_mf_hf.test_vals = [0.250241, -0.743036, -1.049060, -0.753332, 208.023676, 355.360000] #last 7 lines
    sp_pinArray_cht_2d_mf_hf.su2_exec  = "mpirun -n 2 SU2_CFD"
    sp_pinArray_cht_2d_mf_hf.timeout   = 1600
    sp_pinArray_cht_2d_mf_hf.tol       = 0.00001
    sp_pinArray_cht_2d_mf_hf.multizone = True
    test_list.append(sp_pinArray_cht_2d_mf_hf)

    # simple small 3D pin case massflow periodic with heatflux BC
    sp_pinArray_3d_cht_mf_hf_tp           = TestCase('sp_pinArray_3d_cht_mf_hf_tp')
    sp_pinArray_3d_cht_mf_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_3d"
    sp_pinArray_3d_cht_mf_hf_tp.cfg_file  = "configMaster.cfg"
    sp_pinArray_3d_cht_mf_hf_tp.test_iter = 30
    sp_pinArray_3d_cht_mf_hf_tp.test_vals = [0.511984, -3.063453, -0.462699, -0.008477, 214.707868, 429.350000, 368.310000] #last 7 lines
    sp_pinArray_3d_cht_mf_hf_tp.su2_exec  = "mpirun -n 2 SU2_CFD"
    sp_pinArray_3d_cht_mf_hf_tp.timeout   = 1600
    sp_pinArray_3d_cht_mf_hf_tp.tol       = 0.00001
    sp_pinArray_3d_cht_mf_hf_tp.multizone = True
    test_list.append(sp_pinArray_3d_cht_mf_hf_tp)

    ##################################
    ## Streamwise Periodic adjoint ###
    ##################################

    # 2D DA case single zone pressure drop
    da_sp_pinArray_cht_2d_dp_hf           = TestCase('da_sp_pinArray_cht_2d_dp_hf')
    da_sp_pinArray_cht_2d_dp_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    da_sp_pinArray_cht_2d_dp_hf.cfg_file  = "DA_configMaster.cfg"
    da_sp_pinArray_cht_2d_dp_hf.test_iter = 100
    da_sp_pinArray_cht_2d_dp_hf.test_vals = [-4.709021, -3.993726, -3.804347, -3.993726] #last 4 lines
    da_sp_pinArray_cht_2d_dp_hf.su2_exec  = "mpirun -n 2 SU2_CFD_AD"
    da_sp_pinArray_cht_2d_dp_hf.timeout   = 1600
    da_sp_pinArray_cht_2d_dp_hf.tol       = 0.00001
    da_sp_pinArray_cht_2d_dp_hf.multizone = True
    test_list.append(da_sp_pinArray_cht_2d_dp_hf)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    pass_list = [ test.run_test() for test in test_list ]

    # 2D DA case cht pressure drop, heat obj function
    fd_sp_pinArray_cht_2d_dp_hf                = TestCase('fd_sp_pinArray_cht_2d_dp_hf')
    fd_sp_pinArray_cht_2d_dp_hf.cfg_dir        = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    fd_sp_pinArray_cht_2d_dp_hf.cfg_file       = "FD_configMaster.cfg"
    fd_sp_pinArray_cht_2d_dp_hf.test_iter      = 100
    fd_sp_pinArray_cht_2d_dp_hf.su2_exec       = "finite_differences.py -z 2 -n 2 -f"
    fd_sp_pinArray_cht_2d_dp_hf.timeout        = 1600
    fd_sp_pinArray_cht_2d_dp_hf.reference_file = "of_grad_findiff.csv.ref"
    fd_sp_pinArray_cht_2d_dp_hf.test_file      = "FINDIFF/of_grad_findiff.csv"
    fd_sp_pinArray_cht_2d_dp_hf.multizone      = True
    pass_list.append(fd_sp_pinArray_cht_2d_dp_hf.run_filediff())
    test_list.append(fd_sp_pinArray_cht_2d_dp_hf)

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
