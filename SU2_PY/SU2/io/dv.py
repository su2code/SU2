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

    def read_gradient(self, filename):

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
                              'GRADIENT': [line[1].strip()],
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