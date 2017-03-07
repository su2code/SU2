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
from config_options import *

try:
    from collections import OrderedDict
except ImportError:
    from ..util.ordered_dict import OrderedDict

inf = 1.0e20


# ----------------------------------------------------------------------
#  Configuration Class
# ----------------------------------------------------------------------

class Opt(ordered_bunch):
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
        super(Opt, self).__init__(*args, **kwarg)

        # read config if it exists
        if filename:
            try:
                self.read(filename)
            except IOError:
                print 'Could not find config file: %s' % filename
            except:
                print 'Unexpected error: ', sys.exc_info()[0]
                raise

        self.CONFIG_DIRECT = Config(self.DIRECT_FILENAME)
        self.CONFIG_ADJOINT = Config(self.ADJOINT_FILENAME)
        self._filename = filename

        # Set properties of the problem
        self.nKind_DV = len(self.DESIGN_VARIABLES)
        self.nOF = len(self.OBJECTIVE_FUNCTION)

        # Read Design Variables
        self.DV_KIND = self.read_dv()

        # Set up initial guess and boundaries
        (self.x0, self.xb) = self.setx0()

        # Total number of design variables
        self.nDV = len(self.x0)

    def read(self, filename):
        """ reads from a config file """
        opti = read_opt(filename)
        self.update(opti)

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

        return data_dict

    def reread_dv(self,folder):
        """ reads from a dv file """
        # initialize output vector
        x_dv = []

        # For the kind of design variables
        for i in range(0, self.nKind_DV):

            # Choose appropriate file name
            if self.DESIGN_VARIABLES[i] == 'YOUNG_MODULUS': dv_file = folder + 'dv_young.opt'
            elif self.DESIGN_VARIABLES[i] == 'ELECTRIC_FIELD': dv_file = folder + 'dv_efield.opt'

            # Read variables from file (each DV kind has a child class associated, with a method "read")
            nDV, data_list = self.DV_KIND[self.DESIGN_VARIABLES[i]].read(dv_file)

            # For the number of variables in each case
            for j in range(0, nDV):
                x_dv.append(float(data_list[j]['VALUE']))

        return x_dv

    def unpack_dvs(self,dv_new,dv_old=None):
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
        dv_names   = self.DESIGN_VARIABLES
        dv_kinds   = self.DV_KIND
        n_dv_total = self.nDV
        n_dv_kind  = self.nKind_DV

        if not dv_old: dv_old = [0.0] * n_dv_total
        assert len(dv_new) == len(dv_old), 'unexpected design vector length'

        # handle param
 #       param_dv = self['DV_PARAM']

        print dv_names

        k = 0
        # apply scale
        for i in range(0, n_dv_kind):
            n_dv = dv_kinds[dv_names[i]].nDV
            for j in range(0, n_dv):
                dv_scl = float(dv_kinds[dv_names[i]].DV[j]['SCALE'])
                dv_new[k] = dv_new[k] * dv_scl
                dv_old[k] = dv_new[k] * dv_scl

        # Change the parameters of the design variables

#        self['DV_KIND'] = def_dv['KIND']
#        param_dv['PARAM'] = def_dv['PARAM']
#        param_dv['FFDTAG'] = def_dv['FFDTAG']
#        param_dv['SIZE'] = def_dv['SIZE']

        # self.update({'DV_MARKER': def_dv['MARKER'][0],
        #              'DV_VALUE_OLD': dv_old,
        #              'DV_VALUE_NEW': dv_new})

        self.CONFIG_DIRECT.update({'DV_VALUE_OLD': dv_old,
                                   'DV_VALUE_NEW': dv_new})

    def setx0(self):

        x0 = []
        xb_low = []
        xb_up = []

        # For the kind of design variables
        for i in range(0, self.nKind_DV):
            # Shallow copy should be enough here
            dvKind = copy.copy(self.DV_KIND[self.DESIGN_VARIABLES[i]])
            # For the number of design variables in the current DV kind
            for j in range(0,dvKind.nDV):
                x0.append(float(dvKind.DV[j]['VALUE']))
                xb_low.append(float(dvKind.DV[j]['LOWER_BOUND']))
                xb_up.append(float(dvKind.DV[j]['UPPER_BOUND']))

        xb = zip(xb_low, xb_up)  # design bounds

        return x0, xb


    def read_of(self, filename):
        """ reads from a config file """
        opti = read_opt(filename)
        self.update(opti)


    def write(self, filename=''):
        """ updates an existing config file """
        if not filename: filename = self._filename
        assert os.path.exists(filename), 'must write over an existing config file'
        write_opt(filename, self)

    def __getattr__(self, k):
        try:
            return super(Opt, self).__getattr__(k)
        except AttributeError:
            raise AttributeError, 'Config parameter not found'

    def __getitem__(self, k):
        try:
            return super(Opt, self).__getitem__(k)
        except KeyError:
            raise KeyError, 'Config parameter not found: %s' % k

    def __eq__(self, konfig):
        return super(Opt, self).__eq__(konfig)

    def __ne__(self, konfig):
        return super(Opt, self).__ne__(konfig)

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
        output = 'Config: %s' % self._filename
        for k, v in self.iteritems():
            output += '\n    %s= %s' % (k, v)
        return output


#: class Config


# -------------------------------------------------------------------
#  Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def read_opt(filename):
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
                this_value = this_value.split(";")

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


def read_dv_opt(opt):
    """ reads from a dv file """
    # initialize output dictionary
    data_dict = OrderedDict()

    # For the kind of design variables
    for i in range(0, opt.nKind_DV):
        if opt.DESIGN_VARIABLES[i] == 'YOUNG_MODULUS':
            data_dict[opt.DESIGN_VARIABLES[i]] = DV_FEA('dv_young.opt')
        if opt.DESIGN_VARIABLES[i] == 'ELECTRIC_FIELD':
            data_dict[opt.DESIGN_VARIABLES[i]] = DV_FEA('dv_efield.opt')

    return data_dict




# -------------------------------------------------------------------
#  Set SU2 Configuration Parameters
# -------------------------------------------------------------------

def write_opt(filename, param_dict):
    """ updates an existing config file """

    temp_filename = "temp.cfg"
    shutil.copy(filename, temp_filename)
    output_file = open(filename, "w")

    # break pointers
    param_dict = copy.deepcopy(param_dict)

    for raw_line in open(temp_filename):
        # remove line returns
        line = raw_line.strip('\r\n')

        # make sure it has useful data
        if not "=" in line:
            output_file.write(raw_line)
            continue

        # split across equals sign
        line = line.split("=")
        this_param = line[0].strip()
        old_value = line[1].strip()

        # skip if parameter unwanted
        if not param_dict.has_key(this_param):
            output_file.write(raw_line)
            continue

        # start writing parameter
        new_value = param_dict[this_param]
        output_file.write(this_param + "= ")

        # handle parameter types
        for case in switch(this_param):

            # comma delimited list of floats
            if case("DV_VALUE_NEW"): pass
            if case("DV_VALUE_OLD"): pass
            if case("DV_VALUE"):
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write("%s" % new_value[i_value])
                    if i_value + 1 < n_lists:
                        output_file.write(", ")
                break

            # comma delimited list of strings no paren's
            if case("DV_KIND"): pass
            if case("TASKS"): pass
            if case("GRADIENTS"):
                if not isinstance(new_value, list):
                    new_value = [new_value]
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(new_value[i_value])
                    if i_value + 1 < n_lists:
                        output_file.write(", ")
                break

                # comma delimited list of strings inside paren's
            if case("MARKER_EULER"): pass
            if case("MARKER_FAR"): pass
            if case("MARKER_PLOTTING"): pass
            if case("MARKER_MONITORING"): pass
            if case("MARKER_SYM"): pass
            if case("DV_MARKER"):
                if not isinstance(new_value, list):
                    new_value = [new_value]
                output_file.write("( ")
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(new_value[i_value])
                    if i_value + 1 < n_lists:
                        output_file.write(", ")
                output_file.write(" )")
                break

                # semicolon delimited lists of comma delimited lists
            if case("DV_PARAM"):

                assert isinstance(new_value['PARAM'], list), 'incorrect specification of DV_PARAM'
                if not isinstance(new_value['PARAM'][0], list): new_value = [new_value]

                for i_value in range(len(new_value['PARAM'])):

                    output_file.write("( ")
                    this_param_list = new_value['PARAM'][i_value]
                    this_ffd_list = new_value['FFDTAG'][i_value]
                    n_lists = len(this_param_list)

                    if this_ffd_list != []:
                        output_file.write("%s, " % this_ffd_list)
                        for j_value in range(1, n_lists):
                            output_file.write("%s" % this_param_list[j_value])
                            if j_value + 1 < n_lists:
                                output_file.write(", ")
                    else:
                        for j_value in range(n_lists):
                            output_file.write("%s" % this_param_list[j_value])
                            if j_value + 1 < n_lists:
                                output_file.write(", ")

                    output_file.write(") ")
                    if i_value + 1 < len(new_value['PARAM']):
                        output_file.write("; ")
                break

            # int parameters
            if case("NUMBER_PART"): pass
            if case("ADAPT_CYCLES"): pass
            if case("TIME_INSTANCES"): pass
            if case("AVAILABLE_PROC"): pass
            if case("UNST_ADJOINT_ITER"): pass
            if case("EXT_ITER"):
                output_file.write("%i" % new_value)
                break

            if case("DEFINITION_DV"):
                n_dv = len(new_value['KIND'])
                if not n_dv:
                    output_file.write("NONE")
                for i_dv in range(n_dv):
                    this_kind = new_value['KIND'][i_dv]
                    output_file.write("( ")
                    output_file.write("%i , " % get_dvID(this_kind))
                    output_file.write("%s " % new_value['SCALE'][i_dv])
                    output_file.write("| ")
                    # markers
                    n_mark = len(new_value['MARKER'][i_dv])
                    for i_mark in range(n_mark):
                        output_file.write("%s " % new_value['MARKER'][i_dv][i_mark])
                        if i_mark + 1 < n_mark:
                            output_file.write(", ")
                    #: for each marker
                    if not this_kind in ['AOA', 'MACH_NUMBER']:
                        output_file.write(" | ")
                        # params
                        if this_kind in ['FFD_SETTING', 'FFD_ANGLE_OF_ATTACK', 'FFD_CONTROL_POINT', 'FFD_NACELLE',
                                         'FFD_GULL', 'FFD_TWIST_ANGLE', 'FFD_TWIST', 'FFD_TWIST_2D', 'FFD_ROTATION',
                                         'FFD_CAMBER', 'FFD_THICKNESS', 'FFD_CONTROL_POINT_2D', 'FFD_CAMBER_2D',
                                         'FFD_THICKNESS_2D']:
                            n_param = len(new_value['PARAM'][i_dv])
                            output_file.write("%s , " % new_value['FFDTAG'][i_dv])
                            for i_param in range(1, n_param):
                                output_file.write("%s " % new_value['PARAM'][i_dv][i_param])
                                if i_param + 1 < n_param:
                                    output_file.write(", ")
                        else:
                            n_param = len(new_value['PARAM'][i_dv])
                            for i_param in range(n_param):
                                output_file.write("%s " % new_value['PARAM'][i_dv][i_param])
                                if i_param + 1 < n_param:
                                    output_file.write(", ")

                                    #: for each param
                    output_file.write(" )")
                    if i_dv + 1 < n_dv:
                        output_file.write("; ")
                #: for each dv
                break

            if case("OPT_OBJECTIVE"):
                i_name = 0
                for name, value in new_value.iteritems():
                    if i_name > 0: output_file.write("; ")
                    output_file.write("%s * %s" % (name, value['SCALE']))
                    i_name += 1
                break

            if case("OPT_CONSTRAINT"):
                i_con = 0
                for con_type in ['EQUALITY', 'INEQUALITY']:
                    this_con = new_value[con_type]
                    for name, value in this_con.iteritems():
                        if i_con > 0: output_file.write("; ")
                        output_file.write("( %s %s %s ) * %s"
                                          % (name, value['SIGN'], value['VALUE'], value['SCALE']))
                        i_con += 1
                        #: for each constraint
                #: for each constraint type
                if not i_con: output_file.write("NONE")
                break

            if case("ELECTRIC_FIELD_MOD"):
                if not isinstance(new_value, list):
                    new_value = [new_value]
                output_file.write("( ")
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(str(new_value[i_value]))
                    if i_value + 1 < n_lists:
                        output_file.write(", ")
                output_file.write(" )")
                break

            if case("ELECTRIC_FIELD_MIN"):
                if not isinstance(new_value, list):
                    new_value = [new_value]
                output_file.write("( ")
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(str(new_value[i_value]))
                    if i_value + 1 < n_lists:
                        output_file.write(", ")
                output_file.write(" )")
                break

            if case("ELECTRIC_FIELD_MAX"):
                if not isinstance(new_value, list):
                    new_value = [new_value]
                output_file.write("( ")
                n_lists = len(new_value)
                for i_value in range(n_lists):
                    output_file.write(str(new_value[i_value]))
                    if i_value + 1 < n_lists:
                        output_file.write(", ")
                output_file.write(" )")
                break

                # default, assume string, integer or unformatted float
            if case():
                output_file.write('%s' % new_value)
                break

                #: for case

        # remove from param dictionary
        del param_dict[this_param]

        # next line
        output_file.write("\n")

        #: for each line

    # check that all params were used
    for this_param in param_dict.keys():
        if not this_param in ['JOB_NUMBER']:
            print ('Warning: Parameter %s not found in config file and was not written' % (this_param))

    output_file.close()
    os.remove(temp_filename)


#: def write_config()


