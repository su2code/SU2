#!/usr/bin/env python 

## \file convert_routines_c2cpp.py
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
target_directory = su2_home + 'SU2_CFD/src/cpp_routines_d'

# define conversion subroutine
# ----
def convert_c2cpp(source_routine):

    source_routine += '_d'
    source_location = source_directory + '/' + source_routine + '.c'

    # check to see if converted file has already been created (if it exists, skip conversion)
    target_routine = source_routine

    target_location = target_directory + '/' + target_routine + '.cpp'

    if os.path.isfile(target_location):
        print target_routine + ' exists!'

    else:

        # open and read in source file
        source_file = open(source_location, 'r')

        source_file_data = source_file.readlines()

        # convert code and save to new_data
        new_data = []
        temp_data = []
        temp2_data = []
        redefines = []
        ok_vars = []
        
        start = 0
        start_decl = 0
        counter = 0
        end_decl = 0
        declaration_open = 0
        
        # first need to check for variables in the DECL_LIST
        for data_line in source_file_data:
        	if data_line.strip().startswith('//SU2_CPP2C'):
        		data_words = data_line.split()
        		if data_words[1] == 'DECL_LIST':
        			if (data_words[2] == 'START' and declaration_open == 0):
        				declaration_open = 1
        			elif (data_words[2] == 'END' and declaration_open == 1):
        				declaration_open = 0
        		elif data_words[1] == 'VARS' and declaration_open == 1:
        			if data_words[3] == 'SCALAR':
        				for word in data_words[4:]:
        					ok_vars.append(word)
        			elif data_words[3] == 'MATRIX':
        				for word in data_words[4:]:
        					if not word.startswith('SIZE='):
        						ok_vars.append(word)

        for data_line in source_file_data:
            # check to open routine
            # check for open bracket
            if data_line.strip().endswith('{'):
                start = counter

            # check for cpp2c comment
            if data_line.strip().startswith('//SU2_C2CPP DIVIDER'):
                start_decl = start
                end_decl = counter

            counter = counter + 1


        # keep new vars from tapenade
 #       for data_line in source_file_data[start_decl+1:end_decl]:
 #           if data_line.strip().startswith('int') or data_line.strip().startswith('double'):
 #               data_words = data_line.split()
 #               if data_words[1].startswith('arg') or data_words[1].startswith('tmp') or data_words[1].startswith('ii') or data_words[1].startswith('result'):
 #                   temp_data.append(data_line)
 
        for data_line in source_file_data[start_decl+1:end_decl]:
			if data_line.strip().startswith('int') or data_line.strip().startswith('double'):
				data_words = data_line.split()
				if data_words[1].endswith(';'):
					var_match = 0
					if data_words[1].find('['):
						test_word = data_words[1][0:data_words[1].find('[')]
					else:
						test_word = data_words[1][0:-1]
					for var_test in ok_vars:
						if test_word == var_test or test_word == (var_test + 'd'):
							var_match = 1
					if var_match == 0:
						temp_data.append(data_line)
						

        # add new function definition
        for data_line in source_file_data[end_decl:]:
            if data_line.strip().startswith('//SU2_C2CPP'):
                data_words = data_line.split()

                if data_words[1] == 'INVARS': 
                    new_line = 'void ' + target_routine + '('
                    for word in data_words[2:]:
                        new_line += 'double ' + word + ', double ' + word + 'd, '

                elif data_words[1] == 'OUTVARS':
                    for word in data_words[2:]:
                        new_line += 'double ' + word + ', double ' + word + 'd, '
                    new_line = new_line[:-2]

                    new_line += ')\n{\n'
                    new_line += '//SU2_INSERT START\n'
                    temp_data.insert(0, new_line)

                # redefine defines...
                elif data_words[1] == 'REDEFINE':
                    
                    redefines.append(data_words[2])
            #include rest of code
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

            temp2_data.append(data_line)

        # add ending comment for insert

        end = 0
        counter = 0

        for data_line in temp2_data:
            if data_line.strip().endswith('}'):
                end = counter
            counter = counter + 1

        counter = 0
        for data_line in temp2_data:
            if counter == end:
                data_line = data_line.strip()
                data_line = data_line[:-1]
                data_line += '\n//SU2_INSERT END\n}'

            new_data.append(data_line)
            counter = counter + 1

        # open and write target file
        target_file = open(target_location, 'w')

        target_file.writelines(new_data)

# ----

# Convert Euler Roe Set Residual
#source_routine = 'CUpwRoe_Flow__SetResidual'
#convert_c2cpp(source_routine)

#source_routine = 'CSourcePieceWise_Plasma__SetResidual_Chemistry'
#convert_c2cpp(source_routine)

#source_routine = 'CSourcePieceWise_Plasma__SetResidual_Axisymmetric'
#convert_c2cpp(source_routine)

# Convert Green Gauss Gradient
source_routine = 'CTurbSolution__CalcGradient_GG'
convert_c2cpp(source_routine)

# Convert Least Squares Gradient
source_routine = 'CTurbSolution__CalcGradient_LS'
convert_c2cpp(source_routine)

# Convert Primitive Vars
source_routine = 'CTurbSolution__CalcPrimVar_Compressible'
convert_c2cpp(source_routine)

# Convert Laminar Viscosity
source_routine = 'CTurbSolution__CalcLaminarViscosity'
convert_c2cpp(source_routine)

# Convert Eddy Viscosity
source_routine = 'CTurbSASolution__CalcEddyViscosity'
convert_c2cpp(source_routine)

# Convert SA Upwind Residual
source_routine = 'CUpwSca_TurbSA__SetResidual'
convert_c2cpp(source_routine)

# Convert SA Viscous Residual
source_routine = 'CAvgGrad_TurbSA__SetResidual'
convert_c2cpp(source_routine)

# Convert SA Source Residual
source_routine = 'CSourcePieceWise_TurbSA__SetResidual'
convert_c2cpp(source_routine)
