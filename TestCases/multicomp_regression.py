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

    #########################
    ## NEMO solver ###
    #########################

    # Multicomponent flow with variable fluid properties 
    multicompflow_variableprop           = TestCase('multicompflow_variableprop')
    multicompflow_variableprop.cfg_dir   = "incomp_rans/multicomponentflow_variablefluidproperties"
    multicompflow_variableprop.cfg_file  = "c6_cms_venturi_conference.cfg"
    multicompflow_variableprop.test_iter = 5
    multicompflow_variableprop.test_vals = [-5.560560,   -4.412340,   -4.090246,   -5.335476,          50,   -6.729831,          14,   -8.156767]
    multicompflow_variableprop.su2_exec  = "mpirun -n 2 SU2_CFD"
    multicompflow_variableprop.timeout   = 1600
    multicompflow_variableprop.new_output = True
    multicompflow_variableprop.tol       = 0.00001
    test_list.append(multicompflow_variableprop)

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
