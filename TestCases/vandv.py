#!/usr/bin/env python

## \file vandv.py
#  \brief Regression tests for the V&V repository.
#  \note Rules for adding cases here:
#   - Use the SU2 --dry_run mode for configs of large tests.
#   - Restart from converged results for medium problems.
#   - Run small cases (<20s) to convergence.
#  \version 8.5.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

import sys, shutil
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
    p30n30.test_vals = [-11.267106, -11.168215, -11.182822, -10.949673, -14.233489, 0.052235, 2.830394, 1.318894, -1.210646, 1.000000, 12.763000]
    test_list.append(p30n30)

    # This is not part of the V&V cases yet, its tested in this script because it is a relatively long test (~1 min).
    shutil.copy("vandv/rans/30p30n/solution.dat", "vandv/rans/30p30n/solution_0.dat")
    p30n30_ad = TestCase('30P30N_ad')
    p30n30_ad.cfg_dir = "vandv/rans/30p30n"
    p30n30_ad.cfg_file = "config_ad.cfg"
    p30n30_ad.test_iter = 9
    p30n30_ad.test_vals = [-7.283709, -6.072615, -5.995304, -7.197048, -4.568373, -1.167146, -2.316777, 1.1791e-01, 3.4123e+01]
    p30n30_ad.command = TestCase.Command("mpirun -n 2", "SU2_CFD_AD")
    test_list.append(p30n30_ad)

    # flat plate - sst-v1994m
    flatplate_sst1994m           = TestCase('flatplate_sst1994m')
    flatplate_sst1994m.cfg_dir   = "vandv/rans/flatplate"
    flatplate_sst1994m.cfg_file  = "turb_flatplate_sst.cfg"
    flatplate_sst1994m.test_iter = 5
    flatplate_sst1994m.test_vals = [-6.757337, -3.881283, -4.244894, -1.293607, -4.474870, -0.678315, 0.002801]
