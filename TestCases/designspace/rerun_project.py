#!/usr/bin/env python 

## \file rerun_project.py
#  \brief Example project re-run file for SU2
#  \author Trent Lukaczyk.
#  \version 2.0 (beta).
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   PROJECT: NACA 0012 1-DV Sweep Example
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, numpy
import libSU2
from tasks_project import Project

try:
    import scipy.io
    scipy_loaded = True
except ImportError:
    scipy_loaded = False
    
    
# -------------------------------------------------------------------
#  Setup
# -------------------------------------------------------------------    

# filenames
config_filename  = 'config_NACA0012.cfg'
design_filename  = 'design_NACA0012.pkl'
project_filename = 'project_NACA0012.pkl'

# load project
The_Project = libSU2.load_data(project_filename)

# update project folder
The_Project.folder_self = os.getcwd()

# read config
config_data = libSU2.Get_ConfigParams(config_filename)
n_DV = len(config_data['DEFINITION_DV']['KIND'])

# design variable to change
i_DV = 9 # lower surface, half-chord

# design variable values 
DV_vals = numpy.linspace(-0.02,0.02, 11 )

# setup config changes
config_delta = []
for X in DV_vals:
    DV_X = numpy.zeros(n_DV)
    DV_X[i_DV] = X
    config_delta.append( {'VARIABLES':DV_X} )


# -------------------------------------------------------------------
#  Run Project
# -------------------------------------------------------------------    

# evaluate project (saves all design data in process)
design_new,_,_ = The_Project.evaluate(config_delta)

# save project
libSU2.save_data(project_filename,The_Project)

# save matlab data
if scipy_loaded:
    design_name_mat = os.path.splitext( design_filename )[0] + '.mat'
    libSU2.save_data(design_name_mat,design_new)

# done
print '\nDONE !'





