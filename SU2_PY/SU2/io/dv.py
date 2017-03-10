#!/usr/bin/env python

## \file dv.py
#  \brief python package for reading in design variables
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

class DV(ordered_bunch):
    """ config = SU2.io.DV(filename="")

        Starts a generic DV class, an extension of
        ordered_bunch()

        use 1: initialize by reading config file
            dv = SU2.io.DV('filename')

        Parameters can be accessed by item
        ie: config['MESH_FILENAME']

        Methods:
            read()       - read from a dv file
            write()      - write to a dv file (requires existing file)
    """

    _filename = 'dv_kind.opt'

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
        super(DV, self).__init__(*args, **kwarg)

        # read config if it exists
        if filename:
            try:
                (self.nDV, self.DV) = self.read(filename)
            except IOError:
                print 'Could not find DV file: %s' % filename
            except:
                print 'Unexpected error: ', sys.exc_info()[0]
                raise

        self._filename = filename

    def read(self, filename):

        """ reads from a dv file """
        # initialize output list and dictionary
        data_list = []

        input_file = open(filename)

        nDV = 0

        # process each line
        while 1:
            # read the line
            line = input_file.readline()
            if not line:
                break

            # remove line returns
            line = line.strip('\r\n')
            # skip first line (starts with ID, INDEX)
            if line[0] == 'I':
                continue
            # split across equals sign
            line = line.split("\t")

            # store in a dictionary
            dv_definitions = {'INDEX': line[0].strip(),
                              'VALUE': [line[1].strip()],
                              'SCALE': line[2].strip(),
                              'PROP1': line[3].strip(),
                              'PROP2': line[4].strip(),
                              'SIZE': 1}

            # save to output dictionary
            data_dict = dv_definitions

            data_list.append(data_dict.copy())
            nDV += 1

        return nDV, data_list

    def write(self, filename=''):
        """ updates an existing config file """
        if not filename: filename = self._filename
        assert os.path.exists(filename), 'must write over an existing config file'
        write_config(filename, self)

    def __getattr__(self, k):
        try:
            return super(DV, self).__getattr__(k)
        except AttributeError:
            raise AttributeError, 'Config parameter not found'

    def __getitem__(self, k):
        try:
            return super(DV, self).__getitem__(k)
        except KeyError:
            raise KeyError, 'Config parameter not found: %s' % k

    def __eq__(self, konfig):
        return super(DV, self).__eq__(konfig)

    def __ne__(self, konfig):
        return super(DV, self).__ne__(konfig)

    def local_files(self):
        """ removes path prefix from all *_FILENAME params
        """
        for key, value in self.iteritems():
            if key.split('_')[-1] == 'FILENAME':
                self[key] = os.path.basename(value)

    def __repr__(self):
        # return '<Config> %s' % self._filename
        return self.__str__()

    def __str__(self):
        output = 'File: %s' % self._filename
        for k, v in self.iteritems():
            output += '\n    %s= %s' % (k, v)
        return output


#: class Config

# -------------------------------------------------------------------
#  Get SU2 Configuration Parameters
# -------------------------------------------------------------------

def read_dv(filename):
    """ reads a config file """




# -------------------------------------------------------------------
#  Set SU2 Configuration Parameters
# -------------------------------------------------------------------

def write_config(filename, param_dict):
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

class DV_FEA(DV):

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
        super(DV, self).__init__(*args, **kwarg)

        # read config if it exists
        if filename:
            try:
                (self.nDV, self.DV) = self.read(filename)
            except IOError:
                print 'Could not find DV file: %s' % filename
            except:
                print 'Unexpected error: ', sys.exc_info()[0]
                raise

        self._filename = filename

    def read(self, filename):

        """ reads from a dv file """
        # initialize output list and dictionary
        data_list = []
        data_dict = OrderedDict()

        input_file = open(filename)

        nDV = 0

        # process each line
        while 1:
            # read the line
            line = input_file.readline()
            if not line:
                break

            # remove line returns
            line = line.strip('\r\n')
            # skip first line (starts with ID, INDEX)
            if line[0] == 'I':
                continue
            # split across equals sign
            line = line.split("\t")

            # store in a dictionary
            dv_definitions = {'INDEX': line[0].strip(),
                              'VALUE': [line[1].strip()],
                              'SCALE': line[2].strip(),
                              'LOWER_BOUND': [line[3].strip()],
                              'UPPER_BOUND': [line[4].strip()],
                              'SIZE': 1}

            # save to output dictionary
            data_dict = dv_definitions

            data_list.append(data_dict.copy())
            nDV += 1

        return nDV, data_list

    def write(self, folder):

        """ writes to a dv file """
        # initialize output list and dictionary
        data_list = []
        data_dict = OrderedDict()

        output_filename = folder + self._filename
        output_file = open(output_filename, "w")

        nDV = self.nDV

        for i in range(0, nDV):
            output_file.write("%s\t" % self.DV[i]['INDEX'])
            output_file.write("%s\t" % self.DV[i]['VALUE'][0])
            output_file.write("%s\t" % self.DV[i]['SCALE'])
            output_file.write("%s\t" % self.DV[i]['LOWER_BOUND'][0])
            output_file.write("%s\t" % self.DV[i]['UPPER_BOUND'][0])
            output_file.write("\n")

        output_file.close()