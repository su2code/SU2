#!/usr/bin/env python

## \file vandv.py
#  \brief Regression tests for the V&V repository.
#  \note Rules for adding cases here:
#   - Use the SU2 --dry_run mode for configs of large tests.
#   - Restart from converged results for medium problems.
#   - Run small cases (<20s) to convergence.
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

import sys
from TestCase import TestCase

def main():
    '''This program runs a subset of the V&V cases.'''

    test_list = []

    ##########################
    ### Compressible RANS ###
    ##########################

    # 30P30N
    p30n30           = TestCase('30P30N')
    p30n30.cfg_dir   = "vandv/rans/30p30n"
    p30n30.cfg_file  = "config.cfg"
    p30n30.test_iter = 20
    p30n30.test_vals = [-10.642687, -10.303442, -10.496806, -10.253229, -13.517216, 0.050962, 2.828563, 1.317849, -0.215792]
    test_list.append(p30n30)

    #################
    ### RUN TESTS ###
    #################

    for test in test_list:
        test.su2_exec = "mpirun -n 2 SU2_CFD"
        test.timeout = 300
        test.tol = 1e-5
    #end

    pass_list = [ test.run_test() for test in test_list ]

    # Tests summary
    print('==================================================================')
    print('Summary of the V&V tests')
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
