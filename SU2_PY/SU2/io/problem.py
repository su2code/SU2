#!/usr/bin/env python

## \file opt.py
#  \brief python package for optimization files
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
from .config import *
from .dv import *
from .physics import *
from .ofunction import *
from config_options import *

try:
    from collections import OrderedDict
except ImportError:
    from ..util.ordered_dict import OrderedDict

inf = 1.0e20


# ----------------------------------------------------------------------
#  Configuration Class
# ----------------------------------------------------------------------

class Problem(ordered_bunch):
    """ opt = SU2.io.opt(filename="")

        Starts an optimization class, an extension of
        ordered_bunch()

        use 1: initialize by reading opt file
            config = SU2.io.Opt('filename')
        use 2: initialize from dictionary or bunch
            config = SU2.io.Opt(param_dict)
        use 3: initialize empty
            config = SU2.io.Opt()

        Parameters can be accessed by item or attribute
        ie: opt['DESIGN_VARIABLES'] or opt.DESIGN_VARIABLES

        Methods:
            read()       - read from a config file
            write()      - write to a config file (requires existing file)
    """

    _filename = 'config.opt'

    def __init__(self, *args, **kwarg):

        # look for filename in inputs
        if args and isinstance(args[0], str):
            filename = args[0]
            args = args[1:]
        elif kwarg.has_key('filename'):
            filename = kwarg['filename']
            del kwarg['filename']
        else:
            filename = ''

        # initialize ordered bunch
        super(Problem, self).__init__(*args, **kwarg)

        # read config if it exists
        # determined by the extension: .cfg or .opt
        if filename:
            try:
                self.read(filename)
            except IOError:
                print 'Could not find config file: %s' % filename
            except:
                print 'Unexpected error: ', sys.exc_info()[0]
                raise

    def read(self, filename):
        """ reads from a config or an opt file """

        file_split = filename.split('.')
        extension = file_split[1]

        if extension == 'opt':
            provlem = self.read_opt(filename)
        else:
            provlem = self.read_cfg(filename)

        self.update(provlem)

    def read_opt(self, filename):
        """ reads from an optimization file """

        # Initialize an optimization problem
        opt = Problem()

        # Read properties from .opt file
        opti = read_opt_file(filename)
        opt.update(opti)

        # Add references to config, config_adj and filename
        opt.config = Config(opt.DIRECT_FILENAME)
        opt.config_adj = Config(opt.ADJOINT_FILENAME)
        opt._filename = filename

        # Read physics and store in physics object
        opt.physics = opt.read_physics(opt.config, opt.OBJECTIVE_FUNCTION)

        # Read objective function and store in OF object
        opt.ofunction = opt.set_ofunction()

        # Set properties of the problem
        opt.nKind_DV = len(opt.DESIGN_VARIABLES)
        opt.nOF = len(opt.OBJECTIVE_FUNCTION)

        # Read Design Variables
        opt.DV_KIND = opt.read_dv()

        # Copy structure for DV_NEW and DV_OLD
        opt.DV_NEW  = copy.deepcopy(opt.DV_KIND)
        opt.DV_OLD  = copy.deepcopy(opt.DV_KIND)

        # Set up initial guess and boundaries
        (opt.x0, opt.xb) = opt.pack_dvs()

        # Total number of design variables
        opt.nDV = len(opt.x0)

        # Stores the kind of problem (function/gradient)
        opt.kind = 'FUNCTION'

        return opt

    def read_cfg(self, filename):
        """ reads from a config or an opt file """

        # Initialize a regular, direct or adjoint problem
        problem = Problem()

        # Read in config: access will be done through problem.config
        problem.config = Config(filename)

        # Read physics and store in physics object
        problem.physics = problem.read_physics(problem.config, problem.OBJECTIVE_FUNCTION)

        # Stores the kind of problem (function/gradient)
        problem.kind = 'FUNCTION'

        # Room for more stuff

        return problem

    def read_dv(self):
        """ reads from a dv file """
        # initialize output dictionary
        data_dict = OrderedDict()

        # For the kind of design variables
        for i in range(0, self.nKind_DV):
            if self.DESIGN_VARIABLES[i] == 'YOUNG_MODULUS':
                data_dict[self.DESIGN_VARIABLES[i]] = DV_FEA('dv_young.opt')
            if self.DESIGN_VARIABLES[i] == 'ELECTRIC_FIELD':
                data_dict[self.DESIGN_VARIABLES[i]] = DV_FEA('dv_efield.opt')
            if self.DESIGN_VARIABLES[i] == 'STRUCTURAL_DENSITY':
                data_dict[self.DESIGN_VARIABLES[i]] = DV_FEA('dv_density.opt')

        return data_dict

    def dump_dv(self, folder):
        """ writes a group of dv files """

        # For the kind of design variables
        for i in range(0, self.nKind_DV):
            self.DV_KIND[self.DESIGN_VARIABLES[i]].write(folder)

    def dump_dv_new(self, folder):
        """ writes a group of dv files """

        # For the kind of design variables
        for i in range(0, self.nKind_DV):
            self.DV_NEW[self.DESIGN_VARIABLES[i]].write(folder)

    def read_gradients(self, objective):
        """ reads from a dv file """

        # initialize output vector
        df_dx = []

        # For the kind of design variables
        for i in range(0, self.nKind_DV):

            # Choose appropriate file name
            if self.DESIGN_VARIABLES[i] == 'YOUNG_MODULUS': grad_file = 'grad_young.opt'
            elif self.DESIGN_VARIABLES[i] == 'ELECTRIC_FIELD': grad_file = 'grad_efield.opt'
            elif self.DESIGN_VARIABLES[i] == 'STRUCTURAL_DENSITY': grad_file = 'grad_density.opt'

            # Read variables from file (each DV kind has a child class associated, with a method "read")
            nDV, data_list = self.DV_KIND[self.DESIGN_VARIABLES[i]].read_gradient(grad_file)

            # For the number of variables in each case
            for j in range(0, nDV):
                # For the size of each variable
                for k in range(0, data_list[j]['SIZE']):
                    df_dx.append(float(data_list[j]['GRADIENT'][k]))

        gradients = {objective: df_dx}

        return gradients

    def pack_dvs(self):

        x0 = []
        xb_low = []
        xb_up = []

        # For the kind of design variables
        for i in range(0, self.nKind_DV):
            # Shallow copy should be enough here
            dvKind = copy.copy(self.DV_KIND[self.DESIGN_VARIABLES[i]])
            # For the number of design variables in the current DV kind
            for j in range(0,dvKind.nDV):
                # For the number of values per design variable
                for k in range(0,dvKind.DV[j]['SIZE']):
                    x0.append(float(dvKind.DV[j]['VALUE'][k]))
                    xb_low.append(float(dvKind.DV[j]['LOWER_BOUND'][k]))
                    xb_up.append(float(dvKind.DV[j]['UPPER_BOUND'][k]))
                    # SCALE HERE
        xb = zip(xb_low, xb_up)  # design bounds

        return x0, xb

    def pack_new_dvs(self):

        x = []

        # For the kind of design variables
        for i in range(0, self.nKind_DV):
            # Shallow copy should be enough here
            dvKind = copy.copy(self.DV_NEW[self.DESIGN_VARIABLES[i]])
            # For the number of design variables in the current DV kind
            for j in range(0,dvKind.nDV):
                # For the number of values per design variable
                for k in range(0,dvKind.DV[j]['SIZE']):
                    x.append(float(dvKind.DV[j]['VALUE'][k]))
                    # SCALE HERE

        return x

    def unpack_dvs(self, dv_new, dv_old=None):
        """ updates config with design variable vectors
            will scale according to each DEFINITION_DV scale parameter

            Modifies:
                DV_KIND
                DV_MARKER
                DV_PARAM
                DV_VALUE_OLD
                DV_VALUE_NEW

            Inputs:
                dv_new - list or array of new dv values
                dv_old - optional, list or array of old dv values, defaults to zeros

        """

        dv_new = copy.deepcopy(dv_new)
        dv_old = copy.deepcopy(dv_old)

        # handle unpacking cases
        dv_names = self.DESIGN_VARIABLES
        dv_kinds = self.DV_KIND
        n_dv_total = self.nDV

        if not dv_old: dv_old = [0.0] * n_dv_total
        assert len(dv_new) == len(dv_old), 'unexpected design vector length'

        k = 0
        # For the kind of design variables
        for i in range(0, self.nKind_DV):
            # Deep copy here -> Then need to update
            new_dv = copy.deepcopy(self.DV_NEW[self.DESIGN_VARIABLES[i]])
            old_dv = copy.deepcopy(self.DV_OLD[self.DESIGN_VARIABLES[i]])
            # For the number of design variables in the current DV kind
            for j in range(0,new_dv.nDV):
                for iDV in range(0,new_dv.DV[j]['SIZE']):
                    # Read scale
                    dv_scl = float(dv_kinds[dv_names[i]].DV[j]['SCALE'])
                    # Store new and old values in structure
                    new_dv.DV[j]['VALUE'][iDV] = str(dv_new[k] * dv_scl)
                    old_dv.DV[j]['VALUE'][iDV] = str(dv_old[k] * dv_scl)
                    k += 1

            # Update old and new DVs
            self.DV_NEW[self.DESIGN_VARIABLES[i]].update(new_dv)
            self.DV_OLD[self.DESIGN_VARIABLES[i]].update(old_dv)

    def dist(self, provlem):
        """ calculates a distance to another problem

            Inputs:
                provlem    - the problem we want to test the difference from

            Outputs:
                distance   - a float

            Currently only tests the value of DV_NEW

        """

        current_dv = self.pack_new_dvs()
        this_dv    = provlem.pack_new_dvs()

        # Pack as a numpy array
        currdv = np.array(current_dv)
        thisdv = np.array(this_dv)

        # Compute distance
        distance = np.sqrt(np.sum((currdv - thisdv) ** 2))

        return distance

    def set_ofunction(self):
        """ builds an objective function object """
        # initialize output dictionary
        ofunction = OrderedDict()
        # retrieve list of objective functions
        list_of = self.OBJECTIVE_FUNCTION.keys()

        if any(i in optnames_aero for i in list_of):
            ofunction['AERO'] = of_aero(list_of)

        if any(i in optnames_fea for i in list_of):
            ofunction['FEA'] = of_fea(list_of)

        return ofunction

    def read_physics(self, config, ofunction):
        """ reads the properties of a problem from config files """

        # Reads the kind of problem that we want to solve
        kind_problem = config.PHYSICAL_PROBLEM

        for case in switch(kind_problem):
            if case("FLUID_STRUCTURE_INTERACTION"):
                phys = fsi(config, ofunction)
                break
            if case("FEM_ELASTICITY"):
                phys = fea(config, ofunction)
                break
            if case():
                phys = flow(config, ofunction)

        return phys


    def write(self, folder):
        """ writes a new design variable file """
        write_dv_opt(self, folder)

    def __getattr__(self, k):
        try:
            return super(Problem, self).__getattr__(k)
        except AttributeError:
            raise AttributeError, 'Config parameter not found'

    def __getitem__(self, k):
        try:
            return super(Problem, self).__getitem__(k)
        except KeyError:
            raise KeyError, 'Config parameter not found: %s' % k

    def __eq__(self, konfig):
        return super(Problem, self).__eq__(konfig)

    def __ne__(self, konfig):
        return super(Problem, self).__ne__(konfig)

    def local_files(self):
        """ removes path prefix from all *_FILENAME params
        """
        for key, value in self.iteritems():
            if key.split('_')[-1] == 'FILENAME':
                self[key] = os.path.basename(value)

    def __repr__(self):
        # return '<Opt> %s' % self._filename
        return self.__str__()

    def __str__(self):
        output = 'Problem: %s' % self._filename
        for k, v in self.iteritems():
            output += '\n    %s= %s' % (k, v)
        return output


#: class Problem

# -------------------------------------------------------------------
#  Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def read_opt_file(filename):
    """ reads a config file """

    # initialize output dictionary
    data_dict = OrderedDict()

    input_file = open(filename)

    # process each line
    while 1:
        # read the line
        line = input_file.readline()
        if not line:
            break

        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
            continue
        # split across equals sign
        line = line.split("=", 1)
        this_param = line[0].strip()
        this_value = line[1].strip()

        assert not data_dict.has_key(this_param), ('Config file has multiple specifications of %s' % this_param)
        for case in switch(this_param):

            # comma delimited lists of strings with or without paren's
            if case("DESIGN_VARIABLES"):
                # remove white space
                this_value = ''.join(this_value.split())
                # split by comma
                data_dict[this_param] = this_value.split(",")
                break

                # unitary objective definition
            if case('OBJECTIVE_FUNCTION'):
                # remove white space
                this_value = ''.join(this_value.split())
                # split by +
                this_def = {}
                this_value = this_value.split(",")

                for this_obj in this_value:
                    # split by scale
                    this_obj = this_obj.split("*")
                    this_name = this_obj[0]
                    this_scale = 1.0
                    if len(this_obj) > 1:
                        this_scale = float(this_obj[1])
                    this_def.update({this_name: {'SCALE': this_scale}})
                # save to output dictionary
                data_dict[this_param] = this_def
                break

            # unitary constraint definition
            if case('CONSTRAINTS'):
                # remove white space
                this_value = ''.join(this_value.split())
                # check for none case
                if this_value == 'NONE':
                    data_dict[this_param] = {'EQUALITY': OrderedDict(), 'INEQUALITY': OrderedDict()}
                    break
                    # split definitions
                this_value = this_value.split(';')
                this_def = OrderedDict()
                for this_con in this_value:
                    if not this_con: continue  # if no definition
                    # defaults
                    this_obj = 'NONE'
                    this_sgn = '='
                    this_scl = 1.0
                    this_val = 0.0
                    # split scale if present
                    this_con = this_con.split('*')
                    if len(this_con) > 1:
                        this_scl = float(this_con[1])
                    this_con = this_con[0]
                    # find sign
                    for this_sgn in ['<', '>', '=']:
                        if this_sgn in this_con: break
                    # split sign, store objective and value
                    this_con = this_con.strip('()').split(this_sgn)
                    assert len(this_con) == 2, 'incorrect constraint definition'
                    this_obj = this_con[0]
                    this_val = float(this_con[1])
                    # store in dictionary
                    this_def[this_obj] = {'SIGN': this_sgn,
                                          'VALUE': this_val,
                                          'SCALE': this_scl}
                #: for each constraint definition
                # sort constraints by type
                this_sort = {'EQUALITY': OrderedDict(),
                             'INEQUALITY': OrderedDict()}
                for key, value in this_def.iteritems():
                    if value['SIGN'] == '=':
                        this_sort['EQUALITY'][key] = value
                    else:
                        this_sort['INEQUALITY'][key] = value
                #: for each definition
                # save to output dictionary
                data_dict[this_param] = this_sort
                break

                # otherwise
            # string parameters
            if case():
                data_dict[this_param] = this_value
                break

    #: for line

    # hack - twl
    if not data_dict.has_key('ITERATIONS'):
        data_dict['ITERATIONS'] = 100
    if not data_dict.has_key('ACCURACY'):
        data_dict['ACCURACY'] = 1e-10

    return data_dict


#: def read_config()

