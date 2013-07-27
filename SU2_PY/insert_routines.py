#!/usr/bin/env python 

## \file insert_routines.py
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
source_directory = su2_home + 'SU2_CFD/src/cpp_routines_d'

# Target directory
global target_directory
target_directory = su2_home + 'SU2_CFD/src'

# define conversion subroutine
# ----
def insert_routine(source_routine, target_routine, target_file):

    source_routine += '_d'
    source_location = source_directory + '/' + source_routine + '.cpp'

    target_location = target_directory + '/' + target_file


    # open and read in source and target files
    source_file = open(source_location, 'r')

    source_file_data = source_file.readlines()

    target_file = open(target_location, 'r')

    target_file_data = target_file.readlines()

    new_data = []
    temp_data = []
    routine_open = 0
        
    # get function to add in
    for data_line in source_file_data:
        if data_line.strip().startswith('//SU2_INSERT'):
            data_words = data_line.split()

            if data_words[1] == 'START' and routine_open == 0:
                routine_open = 1
            elif data_words[1] == 'END' and routine_open == 1:
                routine_open = 0
        elif routine_open == 1:      
            #include code
            temp_data.append(data_line)

    routine_open = 0

    


    # add in function
    for data_line in target_file_data:
        if data_line.startswith('//SU2_DIFF'):

            new_data.append(data_line)

            data_words = data_line.split()
            
            if data_words[1] == 'START' and routine_open == 0:
                if data_words[2] == target_routine:
                    routine_open = 1

                    new_data.append('\n')

                    for new_data_line in temp_data:
                        new_data.append(new_data_line)
 
                    new_data.append('\n')

            elif data_words[1] == 'END' and routine_open == 1:
                if data_words[2] == target_routine:
                    routine_open = 0

        elif routine_open == 0:
             new_data.append(data_line)

        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)

# ----

# Insert Euler Roe Set Residual
#source_routine = 'CSourcePieceWise_Plasma__SetResidual_Axisymmetric'
#target_routine = 'CSourcePieceWise_Plasma__SetResidual_Axisymmetric'
#target_file = 'numerics_source_ad.cpp'

#insert_routine(source_routine, target_routine, target_file)

# Insert GG Gradient
#source_routine = 'CTurbSolution__CalcGradient_GG'
#target_routine = 'CTurbSolution__CalcGradient_GG'
#target_file = 'solution_direct_turbulent_ad.cpp'

#insert_routine(source_routine, target_routine, target_file)

# Insert LS Gradient
#source_routine = 'CTurbSolution__CalcGradient_LS'
#target_routine = 'CTurbSolution__CalcGradient_LS'
#target_file = 'solution_direct_turbulent_ad.cpp'

#insert_routine(source_routine, target_routine, target_file)

# Insert Primitive Variables
#source_routine = 'CTurbSolution__CalcPrimVar_Compressible'
#target_routine = 'CTurbSolution__CalcPrimVar_Compressible'
#target_file = 'solution_direct_turbulent_ad.cpp'

#insert_routine(source_routine, target_routine, target_file)

# Insert Laminar Viscosity
#source_routine = 'CTurbSolution__CalcLaminarViscosity'
#target_routine = 'CTurbSolution__CalcLaminarViscosity'
#target_file = 'solution_direct_turbulent_ad.cpp'

#insert_routine(source_routine, target_routine, target_file)

# Insert Eddy Viscosity
#source_routine = 'CTurbSASolution__CalcEddyViscosity'
#target_routine = 'CTurbSASolution__CalcEddyViscosity'
#target_file = 'solution_direct_turbulent_ad.cpp'

#insert_routine(source_routine, target_routine, target_file)

# Insert SA Upwind
source_routine = 'CUpwSca_TurbSA__SetResidual'
target_routine = 'CUpwSca_TurbSA__SetResidual'
target_file = 'numerics_convective_ad.cpp'

insert_routine(source_routine, target_routine, target_file)

# Insert SA Viscous
source_routine = 'CAvgGrad_TurbSA__SetResidual'
target_routine = 'CAvgGrad_TurbSA__SetResidual'
target_file = 'numerics_viscous_ad.cpp'

insert_routine(source_routine, target_routine, target_file)

# Insert SA Source
source_routine = 'CSourcePieceWise_TurbSA__SetResidual'
target_routine = 'CSourcePieceWise_TurbSA__SetResidual'
target_file = 'numerics_source_ad.cpp'

insert_routine(source_routine, target_routine, target_file)
