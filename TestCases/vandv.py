#!/usr/bin/env python

## \file vandv.py
#  \brief Regression tests for the V&V repository.
#  \note Rules for adding cases here:
#   - Use the SU2 --dry_run mode for configs of large tests.
#   - Restart from converged results for medium problems.
#   - Run small cases (<20s) to convergence.
#  \version 8.1.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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
    p30n30.test_vals         = [-10.806343, -10.326374, -10.559106, -10.432754, -13.517105, 0.050962, 2.828563, 1.317849, -0.228843]
    p30n30.test_vals_aarch64 = [-10.801521, -10.325747, -10.557163, -10.427274, -13.517118, 0.050962, 2.828563, 1.317849, -0.207763]
    test_list.append(p30n30)

    # flat plate - sst-v1994m
    flatplate_sst1994m           = TestCase('flatplate_sst1994m')
    flatplate_sst1994m.cfg_dir   = "vandv/rans/flatplate"
    flatplate_sst1994m.cfg_file  = "turb_flatplate_sst.cfg"
    flatplate_sst1994m.test_iter = 5
    flatplate_sst1994m.test_vals         = [-13.027926, -10.276119, -11.311717, -8.137517, -10.520065, -5.127385, 0.002775]
    flatplate_sst1994m.test_vals_aarch64 = [-13.028095, -11.271115, -11.532461, -8.387610, -11.417974, -5.116988, 0.002808]
    test_list.append(flatplate_sst1994m)

    # bump in channel - sst-v1994m
    bump_sst1994m           = TestCase('bump_sst1994m')
    bump_sst1994m.cfg_dir   = "vandv/rans/bump_in_channel"
    bump_sst1994m.cfg_file  = "turb_bump_sst.cfg"
    bump_sst1994m.test_iter = 5
    bump_sst1994m.test_vals         = [-13.022054, -9.882710, -10.557148, -7.605034, -10.172437, -5.549948, 0.004904]
    bump_sst1994m.test_vals_aarch64 = [-13.034665, -10.510699, -10.627802, -7.661320, -10.680337, -5.749566, 0.004972]
    test_list.append(bump_sst1994m)

    # SWBLI SA
    swbli_sa           = TestCase('swbli_sa')
    swbli_sa.cfg_dir   = "vandv/rans/swbli"
    swbli_sa.cfg_file  = "config_sa.cfg"
    swbli_sa.test_iter = 5
    swbli_sa.test_vals         = [-11.564511, -10.836187, -11.792765, -10.383947, -15.718717, 0.002212, -2.993991, 1.340100]
    swbli_sa.test_vals_aarch64 = [-11.564511, -10.836187, -11.792765, -10.383947, -15.718717, 0.002212, -2.993991, 1.340100]
    test_list.append(swbli_sa)


    # SWBLI - sst-v2003m
    swbli_sst           = TestCase('swbli_sst')
    swbli_sst.cfg_dir   = "vandv/rans/swbli"
    swbli_sst.cfg_file  = "config_sst.cfg"
    swbli_sst.test_iter = 5
    swbli_sst.test_vals = [-11.528112, -10.961624, -11.903226, -10.630539, -11.117619, -4.573066, 0.002318, -2.905628, -4.037947, 1.340100]
    test_list.append(swbli_sst)

    ##########################
    ### Incompressible RANS ###
    ##########################

    # Sandia jet - sst-v2003m
    sandiajet_sst           = TestCase('sandiajet_sst')
    sandiajet_sst.cfg_dir   = "vandv/species_transport/sandia_jet"
    sandiajet_sst.cfg_file  = "validation.cfg"
    sandiajet_sst.test_iter = 5
    sandiajet_sst.test_vals = [-17.169907, -13.518707, -15.442566, -12.021165, -9.660040, -15.289842, 5.000000, -2.746249, 5.000000, -4.836800, 5.000000, -3.966350, 0.000259, 0.000000, 0.000000, 0.000259, 4047.400000, 3946.800000, 49.161000, 51.433000]
    sandiajet_sst.test_vals_aarch64 = [-17.069026, -13.156800, -15.290567, -11.689831, -9.349978, -14.907311, 5.000000, -2.738947, 5.000000, -4.813747, 5.000000, -3.981740, 0.000259, 0.000000, 0.000000, 0.000259, 4047.400000, 3946.800000, 49.161000, 51.433000]
    test_list.append(sandiajet_sst)

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
