#!/usr/bin/env python 

## \file convert_routines_c2c.py
#  \brief ________________________.
#  \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.6
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# import modules
import sys
import string
import os

# Set up constants:
# SU2_HOME
global su2_home
su2_home = os.environ.get('SU2_HOME')

# Source directory
global source_directory
source_directory = su2_home + 'SU2_CFD/src/c_routines_d'

# Target directory
global target_directory
target_directory = su2_home + 'SU2_CFD/src/c_routines_d'

# define conversion subroutine
# ----
def convert_c2c(source_routine):

    source_routine += '_d'
    source_location = source_directory + '/' + source_routine + '.c'

    # check to see if converted file has already been created (if it exists, skip conversion)
    target_routine = source_routine

    target_location = target_directory + '/' + target_routine + '_debug.c'

    if os.path.isfile(target_location):
        print target_routine + ' exists!'

    else:

        # open and read in source file
        source_file = open(source_location, 'r')

        source_file_data = source_file.readlines()

        # convert code and save to new_data
        new_data = []
        temp_data = []
        redefines = []
        
        # add new function definition
        for data_line in source_file_data:
            if data_line.strip().startswith('//SU2_C2CPP'):
                data_words = data_line.split()

                # redefine defines...
                if data_words[1] == 'REDEFINE':
                    
                    redefines.append(data_words[2])
            #include code
            else:
                temp_data.append(data_line)
        
        # redefine defines...
        for data_line in temp_data:
            for redef in redefines:
                pos = redef.find('=')
                old = redef[:pos]
                new = redef[pos+1:]
                new = new[0] + new[2:]

                data_line = data_line.replace(old, new)

            new_data.append(data_line)

        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)

# ----

# Convert Euler Roe Set Residual
#source_routine = 'CUpwRoe_Flow__SetResidual'
#convert_c2c(source_routine)

source_routine = 'CTurbSolution__CalcGradient_GG'
convert_c2c(source_routine)

source_routine = 'CTurbSolution__CalcGradient_LS'
convert_c2c(source_routine)

source_routine = 'CTurbSolution__CalcPrimVar_Compressible'
convert_c2c(source_routine)

source_routine = 'CTurbSolution__CalcLaminarViscosity'
convert_c2c(source_routine)

source_routine = 'CTurbSASolution__CalcEddyViscosity'
convert_c2c(source_routine)

source_routine = 'CUpwSca_TurbSA__SetResidual'
convert_c2c(source_routine)

source_routine = 'CAvgGrad_TurbSA__SetResidual'
convert_c2c(source_routine)

source_routine = 'CSourcePieceWise_TurbSA__SetResidual'
convert_c2c(source_routine)
