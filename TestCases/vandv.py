#!/usr/bin/env python

## \file vandv.py
#  \brief Regression tests for the V&V repository.
#  \note Rules for adding cases here:
#   - Use the SU2 --dry_run mode for configs of large tests.
#   - Restart from converged results for medium problems.
#   - Run small cases (<20s) to convergence.
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
    flatplate_sst1994m.test_vals         = [-13.023358, -9.956752, -11.099910, -7.933220, -10.206577, -5.132343, 0.002808]
    flatplate_sst1994m.test_vals_aarch64 = [-13.022835, -9.956652, -11.102384, -7.928197, -10.206580, -5.132317, 0.002808]
    test_list.append(flatplate_sst1994m)

    # bump in channel - sst-v1994m
    bump_sst1994m           = TestCase('bump_sst1994m')
    bump_sst1994m.cfg_dir   = "vandv/rans/bump_in_channel"
    bump_sst1994m.cfg_file  = "turb_bump_sst.cfg"
    bump_sst1994m.test_iter = 5
    bump_sst1994m.test_vals         = [-12.986182, -10.719941, -10.556276, -7.606531, -10.774915, -5.605156, 0.004972]
    bump_sst1994m.test_vals_aarch64 = [-13.039365, -10.729085, -10.609923, -7.682911, -10.774915, -5.605087, 0.004972]
    test_list.append(bump_sst1994m)

    # SWBLI SA
    swbli_sa           = TestCase('swbli_sa')
    swbli_sa.cfg_dir   = "vandv/rans/swbli"
    swbli_sa.cfg_file  = "config_sa.cfg"
    swbli_sa.test_iter = 5
    swbli_sa.test_vals         = [-11.530796, -10.915564, -12.034495, -10.538719, -15.922522, 0.002233, -3.359164, 1.340100]
    swbli_sa.test_vals_aarch64 = [-11.530796, -10.915564, -12.034495, -10.538719, -15.922522, 0.002233, -3.359164, 1.340100]
    test_list.append(swbli_sa)


    # SWBLI - sst-v2003m
    swbli_sst           = TestCase('swbli_sst')
    swbli_sst.cfg_dir   = "vandv/rans/swbli"
    swbli_sst.cfg_file  = "config_sst.cfg"
    swbli_sst.test_iter = 5
    swbli_sst.test_vals = [-11.527743, -11.150388, -11.944923, -10.750834, -11.116769, -4.030059, 0.002339, -2.730391, -4.067274, 1.276300]
    test_list.append(swbli_sst)

    ##########################
    ### Incompressible RANS ###
    ##########################

    # Sandia jet - sst-v2003m
    sandiajet_sst           = TestCase('sandiajet_sst')
    sandiajet_sst.cfg_dir   = "vandv/species_transport/sandia_jet"
    sandiajet_sst.cfg_file  = "validation.cfg"
    sandiajet_sst.test_iter = 5
    sandiajet_sst.test_vals = [-16.249917, -13.835991, -14.303372, -13.276035, -10.074262, -14.027223, 5, -1.672359, 5, -4.938477, 5, -3.462217, 2.5859e-04, 2.8215e-32, 4.5010e-68, 2.5859e-04, 4.0474e+03, 3.9468e+03, 4.9170e+01, 5.1441e+01]
    sandiajet_sst.test_vals_aarch64 = [-16.249289, -13.833785, -14.303058, -13.276559, -10.267928, -14.027240, 5, -1.676412, 5, -4.815216, 5, -3.462247, 0.000259, 0, 0, 0.000259, 4047.4, 3946.8, 49.17, 51.441]
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
