#!/usr/bin/env python 

## \file optimize_project.py
#  \brief Example project optimization file for SU2
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
#   PROJECT: NACA 0012 Optimization Example
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, numpy
from scipy.optimize import fmin_slsqp 
import libSU2, tasks_su2
from tasks_project import Project
    
try:
    import scipy.io
    scipy_loaded = True
except ImportError:
    scipy_loaded = False
    

# -------------------------------------------------------------------
#  Setup Project
# -------------------------------------------------------------------    

# filenames
config_filename  = 'config_NACA0012.cfg'  # SU2 config file
design_filename  = 'design_NACA0012.pkl'  # design data (objectives,gradients,variables)
project_filename = 'project_NACA0012.pkl' # project data (a class structure)

# load old project
if os.path.exists(project_filename):
    The_Project = libSU2.load_data(project_filename)
    The_Project.folder_self = os.getcwd()        

# start new project
else:
    # new design data
    design_init = { 'VARIABLES'  : [] ,
                    'OBJECTIVES' : {} ,
                    'GRADIENTS'  : {}  }
    libSU2.save_data(design_filename,design_init)
    # start project
    The_Project = Project( config_name = config_filename ,
                           design_name = design_filename  )        
#: if load/start project


# -------------------------------------------------------------------
#  Setup Optimizer
# -------------------------------------------------------------------     

n_DV = len( The_Project.config_current['DEFINITION_DV']['KIND'] )

# Initial guess
x0    = numpy.zeros(n_DV)    

# Bounds
xb    = [] #numpy.array([-0.02,0.02]) * numpy.ones([n_DV,2])

# Functions
f     = tasks_su2.eval_f
df    = tasks_su2.eval_df
ceq   = tasks_su2.eval_ceq
dceq  = tasks_su2.eval_dceq
cieq  = tasks_su2.eval_cieq
dcieq = tasks_su2.eval_dcieq
args  = (The_Project,)

# Max Iterations
its   = 20

# -------------------------------------------------------------------
#  Run Optimizer
# -------------------------------------------------------------------     

fmin_slsqp( x0             = x0      , 
            func           = f       , 
            f_eqcons       = ceq     , 
            f_ieqcons      = cieq    , 
            fprime         = df      ,
            fprime_eqcons  = dceq    , 
            fprime_ieqcons = dcieq   , 
            args           = args    , 
            bounds         = xb      ,
            iter           = its     ,
            iprint         = 2       ,
            full_output    = 1       ,
            acc            = 1e-10   ,
            epsilon        = 1.0e-06  )


# -------------------------------------------------------------------
#  Save Data
# -------------------------------------------------------------------     

# pull all design data
design_data = The_Project.design_current

# save project data
libSU2.save_data( project_filename,The_Project)

# save python data
libSU2.save_data( design_filename,design_data)

# save matlab data
if scipy_loaded:
    design_matname = os.path.splitext( design_filename )[0] + '.mat'
    libSU2.save_data( design_matname,design_data)
