#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.2.1 "Blackbird"
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

    #####################
    ## Species solver ###
    #####################

    # 2 species (1 eq) primitive venturi mixing
    species2_primitiveVenturi           = TestCase('species2_primitiveVenturi')
    species2_primitiveVenturi.cfg_dir   = "species_transport/venturi_primitive_3species"
    species2_primitiveVenturi.cfg_file  = "species2_primitiveVenturi.cfg"
    species2_primitiveVenturi.test_iter = 50
    species2_primitiveVenturi.test_vals = [-5.957517, -5.187476, -5.037298, -5.851420, -1.511976, -6.046002, 5, -0.808614, 5 , -2.351161, 5, -0.247992, 0.000092, 0.000090, 0.000002, 0.000000]
    species2_primitiveVenturi.su2_exec  = "mpirun -n 2 SU2_CFD"
    species2_primitiveVenturi.timeout   = 1600
    species2_primitiveVenturi.new_output = True
    species2_primitiveVenturi.tol       = 0.00001
    test_list.append(species2_primitiveVenturi)

    # 3 species (2 eq) primitive venturi mixing
    species3_primitiveVenturi           = TestCase('species3_primitiveVenturi')
    species3_primitiveVenturi.cfg_dir   = "species_transport/venturi_primitive_3species"
    species3_primitiveVenturi.cfg_file  = "species3_primitiveVenturi.cfg"
    species3_primitiveVenturi.test_iter = 50
    species3_primitiveVenturi.test_vals = [-6.028145, -5.258104, -5.107927, -5.922051, -1.582604, -6.314220, -6.431771, 5, -0.808615, 5, -2.351160, 5, -0.288300, 1.645644, 0.499064, 0.601230, 0.545351]
    species3_primitiveVenturi.su2_exec  = "mpirun -n 2 SU2_CFD"
    species3_primitiveVenturi.timeout   = 1600
    species3_primitiveVenturi.new_output = True
    species3_primitiveVenturi.tol       = 0.00001
    test_list.append(species3_primitiveVenturi)

    # 3 species (2 eq) primitive venturi mixing with inlet files.
    # Note that the residuals are exactly the same as for the non-inlet case which should be the case for a fresh inlet file.
    species3_primitiveVenturi_inletFile           = TestCase('species3_primitiveVenturi_inletFile')
    species3_primitiveVenturi_inletFile.cfg_dir   = "species_transport/venturi_primitive_3species"
    species3_primitiveVenturi_inletFile.cfg_file  = "species3_primitiveVenturi_inletFile.cfg"
    species3_primitiveVenturi_inletFile.test_iter = 50
    species3_primitiveVenturi_inletFile.test_vals = [-6.028145, -5.258104, -5.107927, -5.922051, -1.582604, -6.314220, -6.431771, 5, -0.808615, 5, -2.351160, 5, -0.288300]
    species3_primitiveVenturi_inletFile.su2_exec  = "mpirun -n 2 SU2_CFD"
    species3_primitiveVenturi_inletFile.timeout   = 1600
    species3_primitiveVenturi_inletFile.new_output = True
    species3_primitiveVenturi_inletFile.tol       = 0.00001
    test_list.append(species3_primitiveVenturi_inletFile)

    # rectangle passive transport validation
    species_passive_val           = TestCase('species_passive_val')
    species_passive_val.cfg_dir   = "species_transport/passive_transport_validation"
    species_passive_val.cfg_file  = "passive_transport.cfg"
    species_passive_val.test_iter = 50
    species_passive_val.test_vals = [-16.559189, -16.315116, -16.908670, -4.316833, 10.000000, -4.523292, 8.000000, -5.173152]
    species_passive_val.su2_exec  = "mpirun -n 2 SU2_CFD"
    species_passive_val.timeout   = 1600
    species_passive_val.new_output = True
    species_passive_val.tol       = 0.00001
    test_list.append(species_passive_val)

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
