#!/usr/bin/env python

## \file vandv.py
#  \brief Regression tests for the V&V repository.
#  \note Rules for adding cases here:
#   - Use the SU2 --dry_run mode for configs of large tests.
#   - Restart from converged results for medium problems.
#   - Run small cases (<20s) to convergence.
#  \version 8.2.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
    p30n30.test_iter = 5
    p30n30.test_vals         = [-11.502310, -11.511459, -11.981995, -11.704990, -14.235571, 0.052235, 2.830394, 1.318894, -0.291726]
    p30n30.test_vals_aarch64 = [-11.502310, -11.511459, -11.981995, -11.704990, -14.235571, 0.052235, 2.830394, 1.318894, -0.291726]
    test_list.append(p30n30)

    # flat plate - sst-v1994m
    flatplate_sst1994m           = TestCase('flatplate_sst1994m')
    flatplate_sst1994m.cfg_dir   = "vandv/rans/flatplate"
    flatplate_sst1994m.cfg_file  = "turb_flatplate_sst.cfg"
    flatplate_sst1994m.test_iter = 5
    flatplate_sst1994m.test_vals         = [-13.027324, -10.054934, -11.134650, -7.976571, -10.260300, -4.778563, 0.002796]
    flatplate_sst1994m.test_vals_aarch64 = [-13.024930, -9.634457, -10.707600, -7.558080, -9.926634, -4.910704, 0.002786]
    test_list.append(flatplate_sst1994m)

    # bump in channel - sst-v1994m
    bump_sst1994m           = TestCase('bump_sst1994m')
    bump_sst1994m.cfg_dir   = "vandv/rans/bump_in_channel"
    bump_sst1994m.cfg_file  = "turb_bump_sst.cfg"
    bump_sst1994m.test_iter = 5
    bump_sst1994m.test_vals         = [-13.058028, -10.648326, -10.605014, -7.606233, -10.707705, -5.453705, 0.004903]
    bump_sst1994m.test_vals_aarch64 = [-13.025265, -10.669816, -10.615338, -7.577125, -10.709448, -5.453868, 0.004903]
    test_list.append(bump_sst1994m)

    # SWBLI SA
    swbli_sa           = TestCase('swbli_sa')
    swbli_sa.cfg_dir   = "vandv/rans/swbli"
    swbli_sa.cfg_file  = "config_sa.cfg"
    swbli_sa.test_iter = 5
    swbli_sa.test_vals         = [-11.504424, -10.941741, -12.049925, -10.586263, -16.090385, 0.002242, -1.614365, 1.340100]
    swbli_sa.test_vals_aarch64 = [-11.504424, -10.941741, -12.049925, -10.586263, -16.090385, 0.002242, -1.614365, 1.340100]
    test_list.append(swbli_sa)


    # SWBLI - sst-v2003m
    swbli_sst           = TestCase('swbli_sst')
    swbli_sst.cfg_dir   = "vandv/rans/swbli"
    swbli_sst.cfg_file  = "config_sst.cfg"
    swbli_sst.test_iter = 5
    swbli_sst.test_vals = [-11.501567, -10.850906, -11.566872, -10.370341, -11.409505, -2.181479, 0.001796, -1.452456, -2.944206, 10.000000]
    test_list.append(swbli_sst)

    # DSMA661 - SA   
    dsma661_sa            = TestCase('dsma661_sa')
    dsma661_sa.cfg_dir    = "vandv/rans/dsma661"
    dsma661_sa.cfg_file   = "dsma661_sa_config.cfg"
    dsma661_sa.test_iter  = 5
    dsma661_sa.test_vals  = [-11.270155, -8.240208, -9.000574, -5.954878, -10.737828, 0.155687, 0.024232]
    test_list.append(dsma661_sa)

    # DSMA661 - SST-V2003m
    dsma661_sst           = TestCase('dsma661_sst')
    dsma661_sst.cfg_dir   = "vandv/rans/dsma661"
    dsma661_sst.cfg_file  = "dsma661_sst_config.cfg"
    dsma661_sst.test_iter = 5
    dsma661_sst.test_vals = [-10.994972, -8.484014, -9.144459, -6.025780, -10.590255, -8.031694, 0.155875, 0.023353]
    test_list.append(dsma661_sst)

    ##########################
    ### Incompressible RANS ###
    ##########################

    # Sandia jet - sst-v2003m
    sandiajet_sst           = TestCase('sandiajet_sst')
    sandiajet_sst.cfg_dir   = "vandv/species_transport/sandia_jet"
    sandiajet_sst.cfg_file  = "validation.cfg"
    sandiajet_sst.test_iter = 5
    sandiajet_sst.test_vals = [-17.167460, -14.133874, -15.538854, -14.038830, -10.311748, -15.739547, 5.000000, -2.916316, 5.000000, -5.380194, 5.000000, -4.153689, 0.000258, 0.000000, 0.000000, 0.000258, 4019.500000, 3918.900000, 49.151000, 51.436000]
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
