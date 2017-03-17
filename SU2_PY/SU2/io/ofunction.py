#!/usr/bin/env python

## \file ofunction.py
#  \brief python package for properties of objective functions
#  \author R. Sanchez
#  \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import numpy as np
from ..util import bunch, ordered_bunch, switch
from .tools import *
from config_options import *

try:
    from collections import OrderedDict
except ImportError:
    from ..util.ordered_dict import OrderedDict

inf = 1.0e20


# ----------------------------------------------------------------------
#  Configuration Class
# ----------------------------------------------------------------------

class ofunction(object):
    """ config = SU2.io.phys_problem()

        Starts a generic problem class, an object that hosts
        filenames and other properties that are problem-dependent

        use 1: initialize by reading config object
            dv = SU2.io.problem(config_direct, config_adjoint)

        Parameters can be accessed by item
        ie: config['MESH_FILENAME']

        Methods:
            read()       - read from a dv file
            write()      - write to a dv file (requires existing file)
    """

    def __init__(self, list_of):

        list_of_aero = []

        for key in list_of:
            if key in optnames_aero:
                list_of_aero.append(key)

        self.list_of = list_of_aero



class of_aero(ofunction):
    """ config = SU2.io.physics()

        Starts a generic problem class, an object that hosts
        filenames and other properties that are problem-dependent

        use 1: initialize by reading config object
            dv = SU2.io.problem(config_direct, config_adjoint)

        Parameters can be accessed by item
        ie: config['MESH_FILENAME']

        Methods:
            read()       - read from a dv file
            write()      - write to a dv file (requires existing file)
    """

    def __init__(self, list_of):

        list_of_aero = []

        for key in list_of:
            if key in optnames_aero:
                list_of_aero.append(key)

        self.list_of = list_of_aero

    def read_output(self, problem, state):

        # Read config
        config = problem.config

        plot_format = config['OUTPUT_FORMAT']
        plot_extension = get_extension(plot_format)
        history_filename = config['CONV_FILENAME'] + plot_extension

        # averaging final iterations
        final_avg = problem.config.get('ITER_AVERAGE_OBJ', 0)

        # get history and objectives

        history = read_history(history_filename)
        # this should be removed (the whole idea is to delete "special cases"
        special_cases = []
        aerodynamics = read_aerodynamics(history_filename, special_cases, final_avg)

        # info out
        state.FUNCTIONS.update(aerodynamics)
        state.FILES.DIRECT = config['RESTART_FLOW_FILENAME']
        if 'EQUIV_AREA' in problem.physics.files.keys():
            state.FILES.WEIGHT_NF = 'WeightNF.dat'
        if 'INV_DESIGN_CP' in problem.physics.files.keys():
            state.FILES.TARGET_CP = 'TargetCp.dat'
        if 'INV_DESIGN_HEATFLUX' in problem.physics.files.keys():
            state.FILES.TARGET_HEATFLUX = 'TargetHeatFlux.dat'
            state.HISTORY.DIRECT = history


class of_fea(ofunction):
    """ config = SU2.io.ofunction()

        Starts a FEA problem class, an object that hosts
        filenames and other properties that are problem-dependent

        use 1: initialize by reading config object
            dv = SU2.io.fea_problem(config_direct, config_adjoint)

        Parameters can be accessed by item
        ie: fea_problem.files['MESH']

        Methods:
            read()       - read from a dv file
            write()      - write to a dv file (requires existing file)
    """

    def __init__(self, list_of):

        list_of_fea = []

        for key in list_of:
            if key in optnames_fea:
                list_of_fea.append(key)

        self.list_of = list_of_fea

    def read_output(self, problem, state):

        # Read config
        config = problem.config

#        plot_format = config['OUTPUT_FORMAT']
#        plot_extension = get_extension(plot_format)
#        history_filename = config['CONV_FILENAME'] + plot_extension
# TODO: fix history to read from multicore FEA
#        history = read_history(history_filename)
#        state.HISTORY.DIRECT = history

        state.FILES.DIRECT = problem.physics.files['DIRECT']
        # TODO: add loop for multiphysics

        if 'REFERENCE_GEOMETRY' in self.list_of:
            ofunction_filename = 'of_refgeom.dat'
            plot_file = open(ofunction_filename)
            line = plot_file.readline()
            state.FUNCTIONS.REFERENCE_GEOMETRY = float(line)

        if 'REFERENCE_NODE' in self.list_of:
            ofunction_filename = 'of_refnode.dat'
            plot_file = open(ofunction_filename)
            line = plot_file.readline()
            state.FUNCTIONS.REFERENCE_NODE = float(line)
