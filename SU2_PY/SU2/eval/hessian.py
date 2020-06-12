#!/usr/bin/env python

## \file hessian.py
#  \brief python package for hessians
#  \author T. Dick
#  \version 7.0.5 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

import os, sys, shutil, copy, subprocess
from .. import run  as su2run
from .. import io   as su2io
from .. import util as su2util
from .functions import function, update_mesh
from .gradients import gradient, adjoint
from ..io import redirect_folder, redirect_output
from SU2.eval import functions

# ----------------------------------------------------------------------
#  Main Hessian Interface
# ----------------------------------------------------------------------

def hessian( func_name, method, config, state=None ):
    """ val = SU2.eval.hess(func_name,method,config,state=None)

        Evaluates the Hessian approximation.

        TODO: Additional writes the gradient into state, therefore one call to this is sufficient for both.

        Assumptions:
            Config is already setup for deformation.
            Mesh need not be deformed.
            Updates config and state by reference.
            Redundancy if state.HESSIAN has the key func_name.

        Executes in:
            ./HESSIAN

        Inputs:
            func_name - SU2 objective function name
            method    - needs to be 'DISCRETE_ADJOINT'
            config    - an SU2 config
            state     - optional, an SU2 state

        Outputs:
            A list of floats of hessian values
    """

    # Initialize
    hess = {}
    state = su2io.State(state)
    if func_name == 'ALL':
        raise Exception("func_name = 'ALL' not yet supported")
    func_output = func_name
    if (type(func_name)==list):
        raise Exception("No combo objectives.")
    else:
        config.OPT_COMBINE_OBJECTIVE="NO"
        config.OBJECTIVE_WEIGHT = "1.0"

    # redundancy check
    if not func_output in state['HESSIAN']:

        # Adjoint Hessian
        # ensure Discrete Adjoint is used
        if any([method == 'DISCRETE_ADJOINT']):

            # enable the smoothing
            # activate the corresponding config option
            config.SMOOTH_GRADIENT="YES"

            # Aerodynamics
            if func_output in su2io.historyOutFields:
                if su2io.historyOutFields[func_output]['TYPE'] == 'COEFFICIENT':
                    hess = adjoint_hessian( func_name, config, state )

            elif func_name in su2io.historyOutFields:
                if su2io.historyOutFields[func_name]['TYPE'] == 'COEFFICIENT':
                    hess = adjoint_hessian( func_name, config, state )

            else:
                raise Exception('unknown function name: %s' % func_name)

        else:
            raise Exception('unrecognized hessian method')

        # store
        state['HESSIAN'].update(hess)

    # if not redundant

    # prepare output
    hessian_out = state['HESSIAN'][func_output]

    return copy.deepcopy(hessian_out)

#: def hessian()


# ----------------------------------------------------------------------
#  Adjoint Hessian
# ----------------------------------------------------------------------

def adjoint_hessian( func_name, config, state=None ):
    """ vals = SU2.eval.hessian(func_name,config,state=None)

        Evaluates the hessian for smoothing using the
        adjoint methodology with:
            SU2.eval.func()
            SU2.run.deform()
            SU2.run.direct()
            SU2.run.adjoint()

        Assumptions:
            Config is already setup for deformation.
            Mesh may or may not be deformed.
            Updates config and state by reference.
            Adjoint Redundancy if state.HESSIAN has key func_name.
            Direct Redundancy if state.FUNCTIONS has key func_name.

        Executes in:
            ./HESSIAN_<func_name>

        Inputs:
            func_name - SU2 objective function name
            config    - an SU2 config
            state     - optional, an SU2 state

        Outputs:
            A Bunch() with keys of objective function names
            and values of list of floats of gradient values
    """

    # ----------------------------------------------------
    #  Initialize
    # ----------------------------------------------------

    # initialize
    state = su2io.State(state)
    special_cases = su2io.get_specialCases(config)

    # list of objectives is not allowed at the moment
    func_output = func_name

    HESS_NAME = 'HESSIAN_'+func_output

    # console output
    if config.get('CONSOLE','VERBOSE') in ['QUIET','CONCISE']:
        log_hessian = 'log_Hessian.out'
    else:
        log_hessian = None

    # ----------------------------------------------------
    #  Redundancy Check
    # ----------------------------------------------------

    # master redundancy check
    if func_output in state['HESSIAN']:
        hess = state['HESSIAN']
        return copy.deepcopy(hess)

    # ----------------------------------------------------
    #  Direct Solution
    # ----------------------------------------------------

    # run (includes redundancy checks)
    function( func_name, config, state )

    # ----------------------------------------------------
    #  Adjoint Run for Solution
    # ----------------------------------------------------

    # files to pull
    files = state['FILES']
    pull = []; link = []

    # files: mesh
    name = files['MESH']
    name = su2io.expand_part(name,config)
    link.extend(name)

    # files: direct solution
    name = files['DIRECT']
    name = su2io.expand_zones(name,config)
    name = su2io.expand_time(name,config)
    link.extend(name)

    if 'FLOW_META' in files:
        pull.append(files['FLOW_META'])

    # files: adjoint solution
    if HESS_NAME in files:
        name = files[HESS_NAME]
        name = su2io.expand_zones(name,config)
        name = su2io.expand_time(name,config)
        link.extend(name)
    else:
        config['RESTART_SOL'] = 'NO'

    if not 'OUTPUT_FILES' in config:
        config['OUTPUT_FILES'] = ['RESTART']

    if not 'SURFACE_CSV' in config['OUTPUT_FILES']:
      config['OUTPUT_FILES'].append('SURFACE_CSV')


    # output redirection
    with redirect_folder( HESS_NAME, pull, link ) as push:
        with redirect_output(log_hessian):

            # Format objective in config
            config['OBJECTIVE_FUNCTION'] = func_name

            # # RUN ADJOINT SOLUTION # #

            # local copy
            konfig = copy.deepcopy(config)

            konfig['CONV_FILENAME'] = konfig['CONV_FILENAME'] + '_hessian'

            # Run Solution
            su2run.interface.CFD_SERIAL(konfig)

            # merge
            konfig['SOLUTION_ADJ_FILENAME'] = konfig['RESTART_ADJ_FILENAME']
            su2run.merge(konfig)

            # filenames
            plot_format      = konfig.get('TABULAR_FORMAT', 'CSV')
            plot_extension   = su2io.get_extension(plot_format)
            history_filename = konfig['CONV_FILENAME'] + plot_extension
            special_cases    = su2io.get_specialCases(konfig)

            # get history
            history = su2io.read_history( history_filename, config.NZONES )

            # update super config
            config.update({ 'MATH_PROBLEM' : konfig['MATH_PROBLEM'] ,
                            'OBJECTIVE_FUNCTION'  : konfig['OBJECTIVE_FUNCTION']   })

            # files out
            objective    = konfig['OBJECTIVE_FUNCTION']
            adj_title    = 'ADJOINT_' + objective
            suffix       = su2io.get_adjointSuffix(objective)
            restart_name = konfig['RESTART_FILENAME']
            restart_name = su2io.add_suffix(restart_name,suffix)

            # info out
            info = su2io.State()
            info.FILES[adj_title] = restart_name
            info.HISTORY[adj_title] = history

            # end of run adjoint #

            su2io.restart2solution(config,info)
            state.update(info)

            # Read in the Hessian from the file!

            # get filenames
            objective      = config['OBJECTIVE_FUNCTION']
            hess_filename  = config['HESS_OBJFUNC_FILENAME']

            # read raw hessian
            raw_hessian = su2io.read_hessian(hess_filename)

            objective = objective.split(',')
            hessian = { objective[0] : raw_hessian }

            # update state
            state.HESSIAN.update( hessian )

    #: with output redirection

    # return output
    hess = su2util.ordered_bunch()
    hess[func_output] = state['HESSIAN'][func_output]
    return hess

#: def adjoint_hessian()

