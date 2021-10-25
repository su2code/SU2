#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.1.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

    ######################################
    ### FLAMELET BRANCH                ###
    ######################################

    # insert tests here

    ######################################
    ### NOT FLAMELET                   ###
    ######################################

    # 2D laminar premixed flame with FGM table 
    flamelet_laminar_premixed           = TestCase('flamelet_laminar_premixed')
    flamelet_laminar_premixed.cfg_dir   = "flamelet/laminar_premixed_flame"
    flamelet_laminar_premixed.cfg_file  = "premixed_burner.cfg"
    flamelet_laminar_premixed.test_iter = 10
    flamelet_laminar_premixed.test_vals = [-0.716925, -8.784123, -12.251280, 2317.000000]
    flamelet_laminar_premixed.su2_exec  = "parallel_computation.py -f"
    flamelet_laminar_premixed.timeout   = 1600
    flamelet_laminar_premixed.tol       = 0.00001
    test_list.append(flamelet_laminar_premixed)

    # 2D CHT case streamwise periodicity
    sp_pinArray_cht_2d_dp_hf           = TestCase('sp_pinArray_cht_2d_dp_hf')
    sp_pinArray_cht_2d_dp_hf.cfg_dir   = "incomp_navierstokes/streamwise_periodic/chtPinArray_2d"
    sp_pinArray_cht_2d_dp_hf.cfg_file  = "configMaster.cfg"
    sp_pinArray_cht_2d_dp_hf.test_iter = 100
    sp_pinArray_cht_2d_dp_hf.test_vals = [0.247026, -0.811632, -0.982066, -0.753312, 208.023676, 350.180000] #last 7 lines
    sp_pinArray_cht_2d_dp_hf.su2_exec  = "mpirun -n 2 SU2_CFD"
    sp_pinArray_cht_2d_dp_hf.timeout   = 1600
    sp_pinArray_cht_2d_dp_hf.tol       = 0.00001
    sp_pinArray_cht_2d_dp_hf.multizone = True
    test_list.append(sp_pinArray_cht_2d_dp_hf)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    pass_list = [ test.run_test() for test in test_list ]

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
