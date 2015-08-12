#!/usr/bin/env python 

## \file plot_project.py
#  \brief Example project plotting for SU2
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
import pylab as plt


# -------------------------------------------------------------------
#  Setup
# -------------------------------------------------------------------

# filenames
config_filename  = 'config_NACA0012.cfg'
design_filename  = 'design_NACA0012.pkl'
project_filename = 'project_NACA0012.pkl'

# load design data
design_data = libSU2.load_data(design_filename)

# objectives to plot
obj_plot = ['LIFT','DRAG','MOMENT_Z']

# -------------------------------------------------------------------
#  Plot
# -------------------------------------------------------------------

idv_plot = 9 
x_plot = numpy.array( design_data['VARIABLES']  )
i_plot = numpy.argsort( x_plot[:,idv_plot] , 0).squeeze()
n_obj  = len(obj_plot)

# start plot
plt.figure('OBJECTIVES',(8,10)); plt.clf()
plt.figure('GRADIENTS',(8,10)); plt.clf()    
    
# plot each objective 
for i_obj,objective in enumerate(obj_plot):

    y_plot  = numpy.array( design_data['OBJECTIVES'][objective] )
    
    # plot objectives
    plt.figure('OBJECTIVES');
    plt.subplot(n_obj,1,i_obj+1)
    plt.plot( x_plot[i_plot,idv_plot] , y_plot[i_plot] , 'o-' )
    plt.xlabel('VARIABLE %i (1/CHORD)' % idv_plot)
    plt.ylabel(objective)
    
    if i_obj == 0 : plt.title('NACA0012 OBJECTIVES')
    
    if design_data['GRADIENTS'].has_key(objective):
        
        dy_plot = numpy.array( design_data['GRADIENTS'][objective] )
        
        # plot gradients
        plt.figure('GRADIENTS');
        plt.subplot(n_obj,1,i_obj+1)
        plt.plot( x_plot[i_plot,idv_plot] , dy_plot[i_plot,idv_plot] , 'o-' )
        plt.xlabel('VARIABLE %i (1/CHORD)' % idv_plot)
        plt.ylabel(objective+' SENSITIVITY')            
        
        if i_obj == 0 : plt.title('NACA0012 SENSITIVITIES')

    #: if gradients

#: for objective

# -------------------------------------------------------------------
#  Show
# -------------------------------------------------------------------

# show plot
#plt.show()

# save plot
plt.figure('OBJECTIVES');
plt.savefig('NACA0012_OBJECTIVES.png',dpi=300)
plt.figure('GRADIENTS');
plt.savefig('NACA0012_GRADIENTS.png',dpi=300)

# done
print 'DONE!'
