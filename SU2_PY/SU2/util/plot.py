#!/usr/bin/env python 

## \file plot.py
#  \brief python package for plotting
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
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


def write_plot(filename,plot_format,data_plot,keys_plot=None):
    """ write_plot(filename,plot_format,data_plot,keys_plot=[])
        writes a tecplot or paraview plot of dictionary data 
        data_plot is a dictionary of lists with equal length
        if data_plot is an ordered dictionary, will output in order
        otherwise use keys_plot to specify the order of output
    """
    
    if keys_plot is None: keys_plot = []
    
    if not keys_plot:
        keys_plot = data_plot.keys()
    
    header = ('"') + ('","').join(keys_plot) + ('"') + (' \n')
    
    if plot_format == 'TECPLOT':
        header = 'VARIABLES=' + header
    
    n_lines = 0
    
    for i,key in enumerate(keys_plot):
        value = data_plot[key]
        if i == 0:
            n_lines = len(value)
        else:
            assert n_lines == len(value) , 'unequal plot vector lengths'
        
    plotfile = open(filename,'w')
    plotfile.write(header)
        
    for i_line in range(n_lines):
        for j,key in enumerate(keys_plot):
            value = data_plot[key]
            if j > 0: plotfile.write(", ")
            plotfile.write('%s' % value[i_line])
        plotfile.write('\n')
    
    plotfile.close()
    
    return
    
def tecplot(filename,data_plot,keys_plot=[]):
    write_plot(filename,'TECPLOT',data_plot,keys_plot)

def paraview(filename,data_plot,keys_plot=[]):
    write_plot(filename,'PARAVIEW',data_plot,keys_plot)
        
