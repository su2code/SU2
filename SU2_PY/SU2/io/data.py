#!/usr/bin/env python

## \file data.py
#  \brief python package for data utility functions
#  \author T. Lukaczyk, F. Palacios
#  \version 8.0.0 "Harrier"
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy

if sys.version_info[0] > 2:
    # Py3 pickle now manage both accelerated cPickle and pure python pickle
    # See https://docs.python.org/3/whatsnew/3.0.html#library-changes, 4th item.
    import pickle
else:
    import cPickle as pickle


from .filelock import filelock

# -------------------------------------------------------------------
#  Load a Dictionary of Data
# -------------------------------------------------------------------


def load_data(file_name, var_names=None, file_format="infer", core_name="python_data"):
    """data = load_data( file_name, var_names=None   ,
                      file_format = 'infer'       ,
                      core_name   = 'python_data'  )

    loads dictionary of data from python pickle or matlab struct

    Inputs:
        file_name   - data file name
        var_names   - variable names to read
        file_format - 'infer', 'pickle', or 'matlab'
        core_name   - data is stored under a dictionary with this name

    default looks for variable 'python_data' in file_name
    file_format = pickle, will return any python object
    file_format = matlab, will return strings or float lists and
    requires scipy.io.loadmat
    file_format = infer (default), will infer format from extention
    ('.mat','.pkl')
    """

    try:
        import scipy.io

        scipy_loaded = True
    except ImportError:
        scipy_loaded = False

    if not os.path.exists(file_name):
        raise Exception("File does not exist: %s" % file_name)

    # process file format
    if file_format == "infer":
        if os.path.splitext(file_name)[1] == ".mat":
            file_format = "matlab"
        elif os.path.splitext(file_name)[1] == ".pkl":
            file_format = "pickle"
    assert file_format in ["matlab", "pickle"], "unsupported file format"

    # get filelock
    with filelock(file_name):

        # LOAD MATLAB
        if file_format == "matlab" and scipy_loaded:
            input_data = scipy.io.loadmat(
                file_name=file_name,
                squeeze_me=False,
                chars_as_strings=True,
                struct_as_record=True,
            )
            # pull core variable
            assert core_name in input_data, "core data not found"
            input_data = input_data[core_name]

            # convert recarray to dictionary
            input_data = rec2dict(input_data)

        # LOAD PICKLE
        elif file_format == "pickle":
            input_data = load_pickle(file_name)
            # pull core variable
            assert core_name in input_data, "core data not found"
            input_data = input_data[core_name]

        #: if file_format

    #: with filelock

    # load specified varname into dictionary
    if var_names != None:
        # check for one item name array
        if isinstance(var_names, str):
            var_names = [
                var_names,
            ]
        for key in input_data.keys():
            if not key in var_names:
                del input_data[key]
        #: for key
    #: if var_names

    return input_data


#: def load()


# -------------------------------------------------------------------
#  Save a Dictionary of Data
# -------------------------------------------------------------------


def save_data(
    file_name, data_dict, append=False, file_format="infer", core_name="python_data"
):
    """save_data( file_name, data_dict, append=False ,
               file_format = 'infer'              ,
               core_name='python_data'             ):

    Inputs:
        file_name   - data file name
        data_dict   - a dictionary or bunch to write
        append      - True/False to append existing data
        file_format - 'infer', 'pickle', or 'matlab'
        core_name   - data is stored under a dictionary with this name

    file_format = pickle, will save any pickleable python object
    file_format = matlab, will save strings or float lists and
    requires scipy.io.loadmat
    file_format = infer (default), will infer format from extention
    ('.mat','.pkl')

    matlab format saves data file from matlab 5 and later
    will save nested dictionaries into nested matlab structures
    cannot save classes and modules
    uses scipy.io.loadmat
    """

    try:
        import scipy.io

        scipy_loaded = True
    except ImportError:
        scipy_loaded = False

    # process file format
    if file_format == "infer":
        if os.path.splitext(file_name)[1] == ".mat":
            file_format = "matlab"
        elif os.path.splitext(file_name)[1] == ".pkl":
            file_format = "pickle"
    assert file_format in ["matlab", "pickle"], "unsupported file format"

    # get filelock
    with filelock(file_name):

        # if appending needed
        # TODO: don't overwrite other core_names
        if append == True and os.path.exists(file_name):
            # check file exists
            if not os.path.exists(file_name):
                raise Exception("Cannot append, file does not exist: %s" % file_name)
            # load old data
            data_dict_old = load(
                file_name=file_name,
                var_names=None,
                file_format=file_format,
                core_name=core_name,
            )
            # check for keys not in new data
            for key, value in data_dict_old.iteritems():
                if not (key in data_dict):
                    data_dict[key] = value
            #: for each dict item
        #: if append

        # save to core name
        data_dict = {core_name: data_dict}

        # SAVE MATLAB
        if file_format == "matlab":
            # bunch it
            data_dict = mat_bunch(data_dict)
            # save it
            scipy.io.savemat(
                file_name=file_name,
                mdict=data_dict,
                format="5",  # matlab 5 .mat format
                oned_as="column",
            )
        elif file_format == "pickle":
            # save it
            save_pickle(file_name, data_dict)

        #: if file_format

    #: with filelock

    return


#: def save()


# -------------------------------------------------------------------
#  Load Pickle
# -------------------------------------------------------------------


def load_pickle(file_name):
    """data = load_pickle(file_name)
    loads a pickle with core_data dictionaries
    assumes first entry is a list of all following data names
    returns dictionary of data
    """
    pkl_file = open(file_name, "rb")
    # names = safe_unpickle.loadf(pkl_file)
    names = pickle.load(pkl_file)
    data_dict = dict.fromkeys(names, [])
    for key in names:
        # data_dict[key] = safe_unpickle.loadf(pkl_file)
        data_dict[key] = pickle.load(pkl_file)
    pkl_file.close()
    return data_dict


# -------------------------------------------------------------------
#  Save Pickle
# -------------------------------------------------------------------


def save_pickle(file_name, data_dict):
    """save_pickle(file_name, data_dict)
    saves a core data dictionary
    first pickle entry is a list of all following data names
    """
    pkl_file = open(file_name, "wb")
    names = list(data_dict.keys())
    pickle.dump(names, pkl_file)
    for key in names:
        pickle.dump(data_dict[key], pkl_file)
    pkl_file.close()


# -------------------------------------------------------------------
#  Safe UnPickle
# -------------------------------------------------------------------

# class safe_unpickle(pickle.Unpickler):
#''' adds some safety to unpickling
# checks that only supported classes are loaded
# original source from http://nadiana.com/python-pickle-insecure#comment-144
#'''

## modules : classes considered safe
# PICKLE_SAFE = {
#'copy_reg'        : ['_reconstructor']  ,
#'__builtin__'     : ['object']          ,
#'numpy'           : ['dtype','ndarray'] ,
#'numpy.core.multiarray' : ['scalar','_reconstruct'] ,
#'collections'     : ['OrderedDict']     ,
#'SU2.io.state'    : ['State']           , # SU2 Specific
#'SU2.io.config'   : ['Config']          ,
#'SU2.eval.design' : ['Design']          ,
#'SU2.opt.project' : ['Project']         ,
#'SU2.util.ordered_bunch' : ['OrderedBunch'] ,
#'SU2.util.bunch'  : ['Bunch']           ,
#'tasks_general'   : ['General_Task']    ,
#'tasks_project'   : ['Project','Job']   ,
#'tasks_su2'       : ['Decomp','Deform','Direct','Cont_Adjoint',
#'Multiple_Cont_Adjoint','Finite_Diff','Adapt'] ,
# }

## make sets
# for key in PICKLE_SAFE.keys():
# PICKLE_SAFE[key] = set(PICKLE_SAFE[key])

## check for save module/class
# def find_class(self, module, name):
# if not module in self.PICKLE_SAFE:
# raise pickle.UnpicklingError(
#'Attempting to unpickle unsafe module %s' % module
# )
# __import__(module)
# mod = sys.modules[module]
# if not name in self.PICKLE_SAFE[module]:
# raise pickle.UnpicklingError(
#'Attempting to unpickle unsafe class %s' % name
# )
# klass = getattr(mod, name)
# return klass

## extend the load() and loads() methods
# @classmethod
# def loadf(self, pickle_file): # loads a file like pickle.load()
# return self(pickle_file).load()
# @classmethod
# def loads(self, pickle_string): #loads a string like pickle.loads()
# return self(StringIO.StringIO(pickle_string)).load()


# -------------------------------------------------------------------
#  Convert Record Array to Dictionary
# -------------------------------------------------------------------


def rec2dict(array_in):
    """converts numpy record array to dictionary of lists
    needed for loading matlab data
    assumes array comes from scipy.io.loadmat, with
    squeeze_me = False and struct_as_record = True
    """

    import numpy

    assert isinstance(array_in, numpy.ndarray), "input must be a numpy record array"

    # make sure it's not an object array
    if array_in.dtype == numpy.dtype("object"):
        array_in = array_in.tolist()

    # get record keys/names
    keys = array_in.dtype.names

    # start output dictionary
    dataout = dict.fromkeys(keys, [])

    for key in keys:

        # squeeze_me option puts all items in a two-dim array
        value = array_in[key].tolist()[0][0]

        # convert string
        if isinstance(value[0], unicode):
            value = str(value[0])

        # convert array
        elif isinstance(value, numpy.ndarray):
            # check for another struct level
            if value.dtype.names == None:
                value = value.tolist()
            # telescoping
            else:
                value = rec2dict(value)

        # store value
        dataout[key] = value

    return dataout


#: def rec2dict()


# -------------------------------------------------------------------
#  Flatten a List
# -------------------------------------------------------------------


def flatten_list(input_list):
    """flatten an irregular list of lists of any depth"""
    output_list = []
    for value in input_list:
        if isinstance(value, list):
            output_list.extend(flatten_list(value))  # telescope
        else:
            output_list.append(value)
    return output_list


#: def flatten_list()


# -------------------------------------------------------------------
#  Append Lists in a Nested Dictionary
# -------------------------------------------------------------------


def append_nestdict(base_dict, add_dict):
    """append_nestdict(base_dict,add_dict)
    appends base_dict with add_dict, allowing for
    updating nested dictionaries
    will update base_dict in place
    """

    # break pointer
    add_dict = copy.deepcopy(add_dict)

    # append add_dict keys
    for key in add_dict.keys():

        # ensure base_dict key exists and is a list
        if not base_dict.has_key(key):
            if isinstance(add_dict[key], dict):
                base_dict[key] = {}
            else:
                base_dict[key] = []
        elif not (isinstance(base_dict[key], list) or isinstance(base_dict[key], dict)):
            assert not isinstance(
                add_dict[key], dict
            ), "base[key] is not a dictionary while add[key] is"
            base_dict[key] = [base_dict[key]]

        # append list or telescope
        if isinstance(base_dict[key], dict):
            append_nestdict(base_dict[key], add_dict[key])  # telescope
        else:
            base_dict[key].append(add_dict[key])

    #: for add_dict[key]

    # base_dict will be updated through its pointer
    return


#: def append_nestdict()


# -------------------------------------------------------------------
#  Matlab Bunch Class
# -------------------------------------------------------------------


class mat_bunch:
    """replicates dictionary functionality with class dot structure
    for output of dictionaries to matlab
    """

    def __init__(self, d):
        for k, v in d.items():
            if isinstance(v, dict):
                if len(v):
                    v = mat_bunch(v)
                else:
                    v = []
            self.__dict__[k] = v

    def __dict__(self):
        return self.__dict__

    # items
    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def items(self):
        return self.__dict__.items()

    # dictionary get/set/etc
    def __getitem__(self, k):
        return self.__dict__[k]

    def __setitem__(self, k, v):
        self.__dict__[k] = v

    def __delitem__(self, k):
        del self.__dict__[k]

    def __str__(self):
        print_format = "%s: %s"
        state = []
        for k, v in self.__dict__.items():
            if isinstance(v, mat_bunch):
                v = "%i-item mat_bunch" % len(v.items())
            state.append(print_format % (k, v))
        return "\n".join(state)


#: class mat_bunch
