#!/usr/bin/env python 

## \file differentiate_routines.py
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

import os

# Note: Where included, spaces (' ') are important to ensure the call to the command line makes sense)

# Set up constants:
# SU2_HOME
su2_home = os.environ.get('SU2_HOME')

# Tapenade command
global tapenade
tapenade = 'tapenade' + ' '

# Tapenade mode
global mode
mode = '-tangent' + ' '

# Source directory
source_directory = su2_home + 'SU2_CFD/src/c_routines'

# Output directory
global output_directory
output_directory = '-outputdirectory' + ' ' + su2_home + 'SU2_CFD/src/c_routines_d' + ' '


# Move to SU2_HOME
print 'Moving to SU2_HOME\n'
os.chdir(su2_home)

# define differentiate subroutine:
# ----
def diff_routine(routine_name, invars, outvars, file_list):

# Routine
    root = '-root ' + routine_name + ' '

# Variables
    independents = '-vars "' + invars + '"' + ' '
    dependents = '-outvars "' + outvars + '"' + ' '

    # Send command
    tapenade_call = tapenade + mode + root + independents + dependents + \
                    output_directory + file_list
    print 'Differentiating ' + routine_name
    print 'Command: $ ' + tapenade_call
    os.system(tapenade_call)

# ----


# Differentiate CUpwRoe_Flow::SetResidual:
#routine_name = 'CUpwRoe_Flow__SetResidual'
#invars = 'U_i U_j '
#outvars = 'val_residual '
#file_list = source_directory + '/' + 'CUpwRoe_Flow__SetResidual.c '
#file_list += source_directory + '/' + 'fabs.c '
#diff_routine(routine_name, invars, outvars, file_list)

#routine_name = 'CSourcePieceWise_Plasma__SetResidual_Axisymmetric'
#file_list = ''
#file_list += source_directory + '/' + 'CSourcePieceWise_Plasma__SetResidual_Axisymmetric.c '
#invars = ''
#invars += 'U_i' + ' '
#outvars = ''
#outvars += 'val_residual' + ' '

# Differentiate CTurbSolution::CalcGradient_GG:
routine_name = 'CTurbSolution__CalcGradient_GG'
invars = 'val_U_i val_U_js '
outvars = 'grad_U_i '
file_list = source_directory + '/' + 'CTurbSolution__CalcGradient_GG.c '
diff_routine(routine_name, invars, outvars, file_list)

# Differentiate CTurbSolution::CalcGradient_LS:
routine_name = 'CTurbSolution__CalcGradient_LS'
invars = 'val_U_i val_U_js '
outvars = 'grad_U_i '
file_list = source_directory + '/' + 'CTurbSolution__CalcGradient_LS.c '
diff_routine(routine_name, invars, outvars, file_list)

# Differentiate CTurbSolution::CalcPrimVar_Compressible:
routine_name = 'CTurbSolution__CalcPrimVar_Compressible'
invars = 'Solution '
outvars = 'Primitive '
file_list = source_directory + '/' + 'CTurbSolution__CalcPrimVar_Compressible.c '
diff_routine(routine_name, invars, outvars, file_list)

# Differentiate CTurbSolution::CalcLaminarViscosity:
routine_name = 'CTurbSolution__CalcLaminarViscosity'
invars = 'val_U_i '
outvars = 'val_laminar_viscosity_i '
file_list = source_directory + '/' + 'CTurbSolution__CalcLaminarViscosity.c '
file_list += source_directory + '/' + 'fabs.c '
file_list += source_directory + '/' + 'pow.c '
diff_routine(routine_name, invars, outvars, file_list)

# Differentiate CTurbSolution::CalcEddyViscosity:
routine_name = 'CTurbSASolution__CalcEddyViscosity'
invars = 'val_FlowVars val_TurbVar val_laminar_viscosity '
outvars = 'val_eddy_viscosity '
file_list = source_directory + '/' + 'CTurbSASolution__CalcEddyViscosity.c '
diff_routine(routine_name, invars, outvars, file_list)

# Differentiate CUpwSca_TurbSA::SetResidual:
routine_name = 'CUpwSca_TurbSA__SetResidual'
invars = 'U_i U_j TurbVar_i TurbVar_j '
outvars = 'val_residual '
file_list = source_directory + '/' + 'CUpwSca_TurbSA__SetResidual.c '
file_list += source_directory + '/' + 'fabs.c '
diff_routine(routine_name, invars, outvars, file_list)

# Differentiate CAvgGrad_TurbSA::SetResidual:
routine_name = 'CAvgGrad_TurbSA__SetResidual'
invars = 'U_i U_j Laminar_Viscosity_i Laminar_Viscosity_j TurbVar_i TurbVar_j TurbVar_Grad_i TurbVar_Grad_j '
outvars = 'val_residual '
file_list = source_directory + '/' + 'CAvgGrad_TurbSA__SetResidual.c '
diff_routine(routine_name, invars, outvars, file_list)

# Differentiate CSourcePieceWise_TurbSA::SetResidual:
routine_name = 'CSourcePieceWise_TurbSA__SetResidual'
invars = 'U_i PrimVar_Grad_i Laminar_Viscosity_i TurbVar_i TurbVar_Grad_i '
outvars = 'val_residual '
file_list = source_directory + '/' + 'CSourcePieceWise_TurbSA__SetResidual.c '
file_list += source_directory + '/' + 'fabs.c '
file_list += source_directory + '/' + 'pow.c '
diff_routine(routine_name, invars, outvars, file_list)
