#!/usr/bin/env python 

## \file insert_math.py
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

# Target directory(ies)
global target_directory, header_directory
target_directory = su2_home + 'SU2_CFD/src'
header_directory = su2_home + 'SU2_CFD/include'

# Target file(s)
global target_location, header_location
target_location = target_directory + '/math_ad.cpp'
header_location = header_directory + '/math_ad.hpp'


# define insertion subroutine
# ----
def insert_math(source_routine):

    source_routine += '_d'
    source_location = source_directory + '/' + source_routine + '.c'

    # open and read in source and target files
    source_file = open(source_location, 'r')

    source_file_data = source_file.readlines()

    header_file = open(header_location, 'r')

    header_file_data = header_file.readlines()

    target_file = open(target_location, 'r')

    target_file_data = target_file.readlines()

    new_data = []
    temp_data = []
    temp_data2 = []
    routine_open = 0
    comment_open = 0


    # get function to add in
    for data_line in source_file_data:
        if data_line.strip().startswith('/*'):
            comment_open = 1
        elif data_line.strip().endswith('*/'):
            comment_open = 0
        elif comment_open == 0:
            #include code
            temp_data.append(data_line)

    for data_line in temp_data:
        if not data_line.strip().startswith('//'):
            temp_data2 = data_line
            break

    # add in header
    for data_line in header_file_data:
        if data_line.strip().startswith('//SU2_DIFF'):

            new_data.append(data_line)

            data_words = data_line.split()

            if data_words[1] == 'START' and routine_open == 0:
                if data_words[2] == source_routine:
                    routine_open = 1

                    new_data.append('\n')

                    new_data_line = temp_data2.replace('{',';') 

                    new_data.append(new_data_line)

                    new_data.append('\n')

            elif data_words[1] == 'END' and routine_open == 1:
                if data_words[2] == source_routine:
                    routine_open = 0

        elif routine_open == 0:
            new_data.append(data_line)

        # open and write header file
        header_file = open(header_location, 'w')

        header_file.writelines(new_data)

    new_data = []

    # add in function
    for data_line in target_file_data:
        if data_line.strip().startswith('//SU2_DIFF'):

            new_data.append(data_line)

            data_words = data_line.split()

            if data_words[1] == 'START' and routine_open == 0:
                if data_words[2] == source_routine:
                    routine_open = 1

                    new_data.append('\n')

                    for new_data_line in temp_data:
                        new_data.append(new_data_line)

                    new_data.append('\n')

            elif data_words[1] == 'END' and routine_open == 1:
                if data_words[2] == source_routine:
                    routine_open = 0

        elif routine_open == 0:
            new_data.append(data_line)

        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)

# ----

# Insert differentiated fabs
source_routine = 'fabs'
source_directory = su2_home + 'SU2_CFD/src/c_routines_d'
insert_math(source_routine)

# Insert exact differentiated pow
source_routine = 'pow'
source_directory = su2_home + 'SU2_CFD/src/c_routines_d/exact'
insert_math(source_routine)