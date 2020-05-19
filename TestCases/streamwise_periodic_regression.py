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
    streamwise_periodic_cylinder.test_vals = [30, -7.841567, -6.794739, -6.997455] #last 4 lines
    streamwise_periodic_cylinder.su2_exec  = "parallel_computation.py -f"
    streamwise_periodic_cylinder.timeout   = 1600
    streamwise_periodic_cylinder.tol       = 0.00001
    test_list.append(streamwise_periodic_cylinder)

    # 3D laminar channnel with 1 cell in flow direction, streamwise periodic
    sp_pipeSlice_3d_dp_hf_tp     = TestCase('sp_pipeSlice_3d_dp_hf_tp')
    sp_pipeSlice_3d_dp_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/sp_pipeSlice_3d_dp_hf_tp"
    sp_pipeSlice_3d_dp_hf_tp.cfg_file  = "sp_pipeSlice_3d_dp_hf_tp.cfg"
    sp_pipeSlice_3d_dp_hf_tp.test_iter = 10
    sp_pipeSlice_3d_dp_hf_tp.test_vals = [10, -10.352122, -10.185237, -10.185237] #last 4 lines
    sp_pipeSlice_3d_dp_hf_tp.su2_exec  = "parallel_computation.py -f"
    sp_pipeSlice_3d_dp_hf_tp.timeout   = 1600
    sp_pipeSlice_3d_dp_hf_tp.tol       = 0.00001
    test_list.append(sp_pipeSlice_3d_dp_hf_tp)

    # create 2D pin case pressure drop periodic with heatflux BC and temperature periodicity (without turbulence model for now)
    sp_pinArray_2d_dp_hf_tp           = TestCase('sp_pinArray_2d_dp_hf_tp')
    sp_pinArray_2d_dp_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/sp_pinArray_pinArray_2d_dp_hf_tp"
    sp_pinArray_2d_dp_hf_tp.cfg_file  = "sp_pinArray_2d_dp_hf_tp.cfg"
    sp_pinArray_2d_dp_hf_tp.test_iter = 10
    sp_pinArray_2d_dp_hf_tp.test_vals = [10, -10.352122, -10.185237, -10.185237] #last 4 lines
    sp_pinArray_2d_dp_hf_tp.su2_exec  = "parallel_computation.py -f"
    sp_pinArray_2d_dp_hf_tp.timeout   = 1600
    sp_pinArray_2d_dp_hf_tp.tol       = 0.00001
    test_list.append(sp_pinArray_2d_dp_hf_tp)

    # create 2D pin case massflow periodic with heatflux BC and prescribed heat (without turbulence model for now)
    sp_pinArray_2d_mf_hf           = TestCase('sp_pinArray_2d_mf_hf')
    sp_pinArray_2d_mf_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/sp_pinArray_2d_mf_hf"
    sp_pinArray_2d_mf_hf.cfg_file  = "sp_pinArray_2d_mf_hf.cfg"
    sp_pinArray_2d_mf_hf.test_iter = 10
    sp_pinArray_2d_mf_hf.test_vals = [10, -10.352122, -10.185237, -10.185237] #last 4 lines
    sp_pinArray_2d_mf_hf.su2_exec  = "parallel_computation.py -f"
    sp_pinArray_2d_mf_hf.timeout   = 1600
    sp_pinArray_2d_mf_hf.tol       = 0.00001
    test_list.append(sp_pinArray_2d_mf_hf)

    # create simple small 3D pin case massflow periodic with heatflux BC and temperature periodicity (without turbulence model for now)
    sp_pinArray_3d_mf_hf_tp           = TestCase('sp_pinArray_3d_mf_hf_tp')
    sp_pinArray_3d_mf_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/sp_pinArray_3d_mf_hf_tp"
    sp_pinArray_3d_mf_hf_tp.cfg_file  = "sp_pinArray_3d_mf_hf_tp.cfg"
    sp_pinArray_3d_mf_hf_tp.test_iter = 10
    sp_pinArray_3d_mf_hf_tp.test_vals = [10, -10.352122, -10.185237, -10.185237] #last 4 lines
    sp_pinArray_3d_mf_hf_tp.su2_exec  = "parallel_computation.py -f"
    sp_pinArray_3d_mf_hf_tp.timeout   = 1600
    sp_pinArray_3d_mf_hf_tp.tol       = 0.00001
    test_list.append(sp_pinArray_3d_mf_hf_tp)

    # create 2D CHT case with HF BC and  
    sp_pinArray_cht_2d_mf_hf           = TestCase('sp_pinArray_cht_2d_mf_hf')
    sp_pinArray_cht_2d_mf_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/sp_pinArray_cht_2d_mf_hf"
    sp_pinArray_cht_2d_mf_hf.cfg_file  = "sp_pinArray_cht_2d_mf_hf.cfg"
    sp_pinArray_cht_2d_mf_hf.test_iter = 10
    sp_pinArray_cht_2d_mf_hf.test_vals = [10, -10.352122, -10.185237, -10.185237] #last 4 lines
    sp_pinArray_cht_2d_mf_hf.su2_exec  = "parallel_computation.py -f"
    sp_pinArray_cht_2d_mf_hf.timeout   = 1600
    sp_pinArray_cht_2d_mf_hf.tol       = 0.00001
    test_list.append(sp_pinArray_cht_2d_mf_hf)

    ##################################
    ## Streamwise Periodic adjoint ###
    ##################################

    # 2D DA case single zone pressure drop
    sp_da_pinArray_2d_dp_hf_tp           = TestCase('sp_pinArray_2d_dp_hf_tp')
    sp_da_pinArray_2d_dp_hf_tp.cfg_dir   = "incomp_navierstokes/streamwise_periodic/sp_pinArray_pinArray_2d_dp_hf_tp"
    sp_da_pinArray_2d_dp_hf_tp.cfg_file  = "sp_pinArray_2d_dp_hf_tp.cfg"
    sp_da_pinArray_2d_dp_hf_tp.test_iter = 10
    sp_da_pinArray_2d_dp_hf_tp.test_vals = [10, -10.352122, -10.185237, -10.185237] #last 4 lines
    sp_da_pinArray_2d_dp_hf_tp.su2_exec  = "parallel_computation.py -f"
    sp_da_pinArray_2d_dp_hf_tp.timeout   = 1600
    sp_da_pinArray_2d_dp_hf_tp.tol       = 0.00001
    test_list.append(sp_pinArray_2d_dp_hf_tp)
    
    # 2D DA case cht pressure drop, heat obj function
    sp_da_pinArray_cht_2d_mf_hf           = TestCase('sp_pinArray_cht_2d_mf_hf')
    sp_da_pinArray_cht_2d_mf_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/sp_pinArray_cht_2d_mf_hf"
    sp_da_pinArray_cht_2d_mf_hf.cfg_file  = "sp_pinArray_cht_2d_mf_hf.cfg"
    sp_da_pinArray_cht_2d_mf_hf.test_iter = 10
    sp_da_pinArray_cht_2d_mf_hf.test_vals = [10, -10.352122, -10.185237, -10.185237] #last 4 lines
    sp_da_pinArray_cht_2d_mf_hf.su2_exec  = "parallel_computation.py -f"
    sp_da_pinArray_cht_2d_mf_hf.timeout   = 1600
    sp_da_pinArray_cht_2d_mf_hf.tol       = 0.00001
    test_list.append(sp_pinArray_cht_2d_mf_hf)

    pass_list = [ test.run_test() for test in test_list ]

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
