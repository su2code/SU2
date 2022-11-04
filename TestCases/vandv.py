#!/usr/bin/env python

## \file vandv.py
#  \brief Regression tests for the V&V repository.
#  \note Rules for adding cases here:
#   - Use the SU2 --dry_run mode for configs of large tests.
#   - Restart from converged results for medium problems.
#   - Run small cases (<20s) to convergence.
#  \version 7.4.0 "Blackbird"
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
    p30n30.test_vals         = [-10.639125, -10.302345, -10.493880, -10.249452, -13.517221, 0.050962, 2.828563, 1.317849, -0.217586]
    p30n30.test_vals_aarch64 = [-10.636663, -10.298256, -10.479415, -10.246700, -13.517161, 0.050962, 2.828563, 1.317849, -0.209200]
    test_list.append(p30n30)

    # flat plate - sst-v1994m
    flatplate_sst1994m_restart           = TestCase('flatplate_sst1994m')
    flatplate_sst1994m_restart.cfg_dir   = "vandv/rans/flatplate"
    flatplate_sst1994m_restart.cfg_file  = "turb_flatplate_sst.cfg"
    flatplate_sst1994m_restart.test_iter = 5
    flatplate_sst1994m_restart.test_vals         = [-10.451866, -6.714739, -3.397794, -0.002535, 0.002809]
    flatplate_sst1994m_restart.test_vals_aarch64 = [-13.022949, -10.035600, -5.142348, -0.002535, 0.002809]
    test_list.append(flatplate_sst1994m_restart)

    # bump in channel - sst-v1994m 
    bump_sst1994m           = TestCase('bump_sst1994m')
    bump_sst1994m.cfg_dir   = "vandv/rans/bump_in_channel"
    bump_sst1994m.cfg_file  = "turb_bump_sst.cfg"
    bump_sst1994m.test_iter = 50
    bump_sst1994m.test_vals         = [-9.357032, -5.944060, -2.100896, 0.024576, 0.004967]
    bump_sst1994m.test_vals_aarch64 = [-13.085614, -10.325319, -5.559183, 0.024576, 0.004967]
    test_list.append(bump_sst1994m)

    # bump in channel - sst-v1994m restart
    bump_sst1994m_restart           = TestCase('bump_sst1994m_restart')
    bump_sst1994m_restart.cfg_dir   = "vandv/rans/bump_in_channel"
    bump_sst1994m_restart.cfg_file  = "turb_bump_sst_restart.cfg"
    bump_sst1994m_restart.test_iter = 5
    bump_sst1994m_restart.test_vals         = [-9.357032, -5.944060, -2.100896, 0.024576, 0.004967]
    bump_sst1994m_restart.test_vals_aarch64 = [-13.085614, -10.325319, -5.559183, 0.024576, 0.004967]
    test_list.append(bump_sst1994m_restart)

    # SWBLI SA
    swbli_sa           = TestCase('swbli_sa')
    swbli_sa.cfg_dir   = "vandv/rans/swbli"
    swbli_sa.cfg_file  = "config_sa.cfg"
    swbli_sa.test_iter = 20
    swbli_sa.test_vals         = [-11.029255, -10.511982, -11.400926, -10.128471, -14.536798, 0.002233, -2.608466, 2.786]
    test_list.append(swbli_sa)

    # SWBLI - sst-v2003m restart
    swbli_sst_restart           = TestCase('swbli_sst')
    swbli_sst_restart.cfg_dir   = "vandv/rans/swbli"
    swbli_sst_restart.cfg_file  = "config_sst.cfg"
    swbli_sst_restart.test_iter = 20
    swbli_sst_restart.test_vals = [-8.059128, -7.393263, -7.956248, -6.897129, -7.268286, -2.491667, 0.002315, -2.435645, -3.182638, 1.000000]
    test_list.append(swbli_sst_restart)

    #################
    ### RUN TESTS ###
    #################

    for test in test_list:
        test.command = TestCase.Command("mpirun -n 2", "SU2_CFD")
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
