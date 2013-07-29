#!/usr/bin/env python 

## \file convert_routines_cpp2c.py
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
import sys, string, os

# Set up constants:
# SU2_HOME
global su2_home
su2_home = os.environ.get('SU2_HOME')

# Source directory
global source_directory
source_directory = su2_home + 'SU2_CFD/src'

# Target directory
global target_directory
target_directory = source_directory + '/c_routines'

# define conversion subroutine
# ----
def convert_cpp2c(source_routine, source_location):

    # check to see if converted file has already been created (if it exists, skip conversion)
    target_routine = source_routine.replace('::', '__')

    target_location = target_directory + '/' + target_routine + '.c'

    if os.path.isfile(target_location):
        print target_routine + ' exists!'

    else:

        # open and read in source file
        source_file = open(source_location, 'r')

        source_file_data = source_file.readlines()

        # convert code and save to new_data
        new_data = []
        routine_open = 0
        comment_open = 0
        inputs_open = 0
        declaration_open = 0
        subroutine_open = 0
        sub_open = 0
#        subroutine_list = []
        var_type = ''
        define_count = 0
        replace_line_open = 0

        for data_line in source_file_data:
            # check to open routine
            # check for cpp2c comment
            if data_line.strip().startswith('//SU2_CPP2C') and routine_open == 0 :
                # split into words
                data_words = data_line.split()

                # check for start of routine
                if data_words[1] == 'START':
                    if data_words[2] == source_routine:
                        # start declaration:
                        routine_open = 1
                        new_line = 'void ' + target_routine
                        #new_data.append(new_line)

            # if routine is open
            elif routine_open == 1:
                # check for cpp2c comment
                if data_line.strip().startswith('//SU2_CPP2C'):
                    # split into words
                    data_words = data_line.split()

                    # check to start/end commented section
                    if data_words[1] == 'COMMENT':
                        if data_words[2] == 'START':
                            comment_open = 1
                        elif data_words[2] == 'END':
                            comment_open = 0

                    # check for end of routine
                    elif data_words[1] == 'END':
                        if data_words[2] == source_routine:
                            # end routine
                            new_line = '}'
                            new_data.append(new_line)
                            routine_open = 0
                            break

                    # open/close input variable declaration
                    elif data_words[1] == 'CALL_LIST':
                        if (data_words[2] == 'START' and inputs_open == 0):
                            new_line += '('
                            extra_line = ''
                            inputs_open = 1
                        elif (data_words[2] == 'END' and inputs_open == 1):
                            new_line = new_line[:-2]
                            new_line += ')\r\n'
                            new_data.append(new_line)
                            inputs_open = 0

                            new_line = '{\r\n'
                            new_data.append(new_line)
                            divider_line = '\r\n//SU2_C2CPP DIVIDER\r\n' + extra_line
#                            new_data.append(new_line)
#                            new_data.append(extra_line)

#                            divider_line = extra_line

                    # get list of independent variables
                    elif (data_words[1] == 'INVARS' and inputs_open == 1):

                        extra_line += '//SU2_C2CPP INVARS'

                        var_type = 'double'

                        # declare inputs
                        for word in data_words[2:]:
                            new_line += var_type + ' ' + word + ', '
                            extra_line += ' ' + word

                        extra_line += '\r\n' 

                    # get list of dependent variables
                    elif (data_words[1] == 'OUTVARS' and inputs_open == 1):

                        extra_line += '//SU2_C2CPP OUTVARS'

                        var_type = 'double'

                        # declare outputs
                        for word in data_words[2:]:
                            new_line += var_type + ' ' + word + ', '
                            extra_line += ' ' + word

                        extra_line += '\r\n'

                    # get list of doubles/ints that need declaration (either input or internally declare)
                    elif data_words[1] == 'DECL_LIST':
                        if (data_words[2] == 'START' and declaration_open == 0):
                            declaration_open = 1
                            new_data.append(data_line)
                        elif (data_words[2] == 'END' and declaration_open == 1):
                            declaration_open = 0
                            new_data.append(data_line)
                            new_data.append(divider_line)

                    elif data_words[1] == 'VARS':
                        if data_words[2] == 'INT':
                            var_type = 'int'
                        elif data_words[2] == 'DOUBLE':
                            var_type = 'double'

                        if inputs_open == 1:
                            for word in data_words[3:]:
                                new_line += var_type + ' ' + word + ', '
                        elif declaration_open == 1:
                            if data_words[3] == 'SCALAR':
                                new_line = var_type + ' '
                                for word in data_words[4:]:
                                    new_line += word + ', '

                                new_line = new_line[:-2]
                                new_line += ';\r\n'
                                new_data.append(new_line)
                                new_data.append(data_line)
                            elif data_words[3] == 'MATRIX':
                                dimensions = []
                                for word in data_words[4:]:
                                    if word.startswith('SIZE='):
                                        dimensions.append(word[5:])

                                new_line = ''

                                for word in data_words[4:]:
                                    if not word.startswith('SIZE='):
                                        new_line += var_type + ' '
                                        new_line += word

                                        for dim in dimensions:

                                            new_line += '[' + dim + ']'

                                        new_line +=  ';\r\n'

                                new_data.append(new_line)
                                new_data.append(data_line)
            
                    # get #define constants
                    elif data_words[1] == 'DEFINE':
                        extra_line = ''
                        for word in data_words[2:]:
                            define_count = define_count + 1
                            extra_line += '#define ' + word + ' ' + str(1234000 + define_count) + '\r\n'
                            divider_line += '//SU2_C2CPP REDEFINE ' + str(1234000 + define_count) + '=' + word[0] + '_' + word[1:] + '\r\n' # inserts '_' in name to stop Tapenade replacing it with #define
                        new_data.insert(0, extra_line)

                    # get sub routine
                    elif data_words[1] == 'SUBROUTINE':
                        if data_words[2] == 'START' and subroutine_open == 0:
#                            new_line = data_words[3] + '('
                            subroutine_open = 1
                            sub_routine = data_words[3]
                        elif subroutine_open == 1:
                            if data_words[2] == 'END':
#                                new_line = new_line[:-2]
#                                new_line += ');\r\n'
#                                new_data.append(new_line)
                                new_sub_data = []
                                for sub_data_line in sub_data:
                                    var_count = 0
                                    for old_var in sub_vars_old:
                                        sub_data_line = sub_data_line.replace(old_var, sub_vars_new[var_count])

                                        var_count = var_count +1

                                    new_sub_data.append(sub_data_line)
                                for sub_data_line in new_sub_data:
                                    new_data.append(sub_data_line)
                                subroutine_open = 0
                            elif data_words[2] == 'LOCATION':
                                sub_location = source_directory + '/' + data_words[3]
                                sub_file = open(sub_location, 'r')
                                sub_file_data = sub_file.readlines()
                                sub_data = []
                                sub_vars_old = []

                                # now loop through to get sub routine code
                                for sub_data_line in sub_file_data:
                                    if sub_data_line.strip().startswith('//SU2_CPP2C') and sub_open == 0 :
                                        # split into words
                                        sub_data_words = sub_data_line.split()

                                        # check for start of subroutine
                                        if sub_data_words[1] == 'SUB':
                                            if sub_data_words[2] == 'START':
                                                if sub_data_words[3] == sub_routine:
                                                # start declaration:
                                                    sub_open = 1
                                    elif sub_open == 1:
                                        if sub_data_line.strip().startswith('//SU2_CPP2C'):
                                            # split into words
                                            sub_data_words = sub_data_line.split()

                                            if sub_data_words[1] == 'SUB':
                                                if sub_data_words[2] == 'END' and sub_data_words[3] == sub_routine:

                                                    sub_open = 0
                                                elif sub_data_words[2] == 'VARS':
                                                    for sub_var in sub_data_words[3:]:
                                                        sub_vars_old.append(sub_var)

                                        else:
                                            sub_data.append(sub_data_line)

                            elif data_words[2] == 'VARS':
                                sub_vars_new = []

                                for word in data_words[3:]:
                                    sub_vars_new.append(word)

                    # get new line
                    if data_words[1] == 'REPLACE_LINE':
                        if data_words[2] == 'START' and replace_line_open == 0:
                            replace_line_open = 1
                        elif data_words[2] == 'END' and replace_line_open == 1:
                            replace_line_open = 0
                        elif data_words[2] == 'NEW_LINE':
                            new_data_line = ''
                            for data_word in data_words[3:]:
                                new_data_line += data_word + ' '
                            new_data.append(new_data_line)

                    # get any required #includes (differentiable forms of non-differentiable standard functions)
#                    elif data_words[1] == 'INCLUDE':
#                        subroutine_known = 0
#                        for subroutine in subroutine_list:
#                            if (subroutine == data_words[2]):
#                                 subroutine_known = 1
#                        if (subroutine_known == 0):
#                            subroutine_list.append(data_words[2])

                # deal with non comment line
                elif comment_open == 0 and subroutine_open == 0 and replace_line_open == 0:
                    # cut and paste
                    new_data.append(data_line)

        # add in includes:
#        for subroutine in subroutine_list:
#            extra_line = '#include "' + subroutine + '.c"\r\n'
#            new_data.insert(0, extra_line)



        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)

##

# Convert Euler Roe Set Residual
# routine to convert
source_routine = 'CUpwRoe_Flow::SetResidual'
source_location = source_directory + '/numerics_convective.cpp'
convert_cpp2c(source_routine, source_location)

#source_routine = 'CSourcePieceWise_Plasma::SetResidual_Axisymmetric'
#source_location = source_directory + '/numerics_source.cpp'
#convert_cpp2c(source_routine, source_location)

# Convert Green Gauss routine
source_routine = 'CTurbSolution::CalcGradient_GG'
source_location = source_directory + '/solution_direct_turbulent.cpp'
convert_cpp2c(source_routine, source_location)

# Convert Least Square routine
source_routine = 'CTurbSolution::CalcGradient_LS'
source_location = source_directory + '/solution_direct_turbulent.cpp'
convert_cpp2c(source_routine, source_location)

# Convert conservative to primitive routine
source_routine = 'CTurbSolution::CalcPrimVar_Compressible'
source_location = source_directory + '/solution_direct_turbulent.cpp'
convert_cpp2c(source_routine, source_location)

# Convert Laminar Viscosity routine
source_routine = 'CTurbSolution::CalcLaminarViscosity'
source_location = source_directory + '/solution_direct_turbulent.cpp'
convert_cpp2c(source_routine, source_location)

# Convert Eddy Viscosity routine
source_routine = 'CTurbSASolution::CalcEddyViscosity'
source_location = source_directory + '/solution_direct_turbulent.cpp'
convert_cpp2c(source_routine, source_location)

# Convert SA upwind routine
source_routine = 'CUpwSca_TurbSA::SetResidual'
source_location = source_directory + '/numerics_convective.cpp'
convert_cpp2c(source_routine, source_location)

# Convert  routine
source_routine = 'CAvgGrad_TurbSA::SetResidual'
source_location = source_directory + '/numerics_viscous.cpp'
convert_cpp2c(source_routine, source_location)

# Convert  routine
source_routine = 'CSourcePieceWise_TurbSA::SetResidual'
source_location = source_directory + '/numerics_source.cpp'
convert_cpp2c(source_routine, source_location)