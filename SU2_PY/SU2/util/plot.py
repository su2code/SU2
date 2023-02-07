#!/usr/bin/env python

## \file plot.py
#  \brief python package for plotting
#  \author T. Lukaczyk, F. Palacios
#  \version 7.5.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.


def write_plot(filename,plot_format,data_plot,keys_plot=None):
    """ write_plot(filename,plot_format,data_plot,keys_plot=[])
        writes a tecplot or paraview plot of dictionary data
        data_plot is a dictionary of lists with equal length
        if data_plot is an ordered dictionary, will output in order
        otherwise use keys_plot to specify the order of output
    """

    default_spacing = 16
    indent_spacing  = 0

    if keys_plot is None: keys_plot = []

    if not keys_plot:
        keys_plot = data_plot.keys()

    keys_print  = [ '"'+key+'"' for key in keys_plot ]
    keys_space = [default_spacing] * len(keys_plot)

    header = ''
    if (plot_format == 'TECPLOT'):
        header = 'VARIABLES='
        indent_spacing += 10
    indent_spacing = ' '*indent_spacing

    n_lines = 0
    for i,key in enumerate(keys_plot):
        # check vector lengths
        value = data_plot[key]
        if i == 0:
            n_lines = len(value)
        else:
            assert n_lines == len(value) , 'unequal plot vector lengths'

        # check spacing
        if len(key) > keys_space[i]:
            keys_space[i] = len(key)
        keys_space[i] = "%-" + str(keys_space[i]) + "s"

    plotfile = open(filename,'w')
    plotfile.write(header)
    for i,key in enumerate(keys_print):
        if i > 0: plotfile.write(", ")
        plotfile.write(keys_space[i] % key)
    plotfile.write('\n')

    for i_line in range(n_lines):
        plotfile.write(indent_spacing)
        for j,key in enumerate(keys_plot):
            value = data_plot[key]
            if j > 0: plotfile.write(", ")
            plotfile.write(keys_space[j] % value[i_line])
        plotfile.write('\n')

    plotfile.close()

    return

def tecplot(filename,data_plot,keys_plot=[]):
    write_plot(filename,'TECPLOT',data_plot,keys_plot)

def paraview(filename,data_plot,keys_plot=[]):
    write_plot(filename,'CSV',data_plot,keys_plot)

