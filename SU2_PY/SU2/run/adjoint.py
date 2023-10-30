#!/usr/bin/env python

## \file adjoint.py
#  \brief python package for running adjoint problems
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

import copy

from .. import io as su2io
from .merge import merge as su2merge
from .interface import CFD as SU2_CFD

# ----------------------------------------------------------------------
#  Adjoint Simulation
# ----------------------------------------------------------------------


def adjoint(config):
    """info = SU2.run.adjoint(config)

    Runs an adjoint analysis with:
        SU2.run.decomp()
        SU2.run.CFD()
        SU2.run.merge()

    Assumptions:
        Does not run Gradient Projection
        Does not rename restart filename to solution filename
        Adds 'adjoint' suffix to convergence filename

    Outputs:
        info - SU2 State with keys:
            HISTORY.ADJOINT_NAME
            FILES.ADJOINT_NAME

    Updates:
        config.MATH_PROBLEM

    Executes in:
        ./
    """

    # local copy
    konfig = copy.deepcopy(config)

    # setup problem
    if konfig.get("GRADIENT_METHOD", "CONTINUOUS_ADJOINT") == "DISCRETE_ADJOINT":
        konfig["MATH_PROBLEM"] = "DISCRETE_ADJOINT"
    else:
        konfig["MATH_PROBLEM"] = "CONTINUOUS_ADJOINT"

    konfig["CONV_FILENAME"] = konfig["CONV_FILENAME"] + "_adjoint"

    # Run Solution
    SU2_CFD(konfig)

    # merge
    konfig["SOLUTION_ADJ_FILENAME"] = konfig["RESTART_ADJ_FILENAME"]
    su2merge(konfig)

    # filenames
    plot_format = konfig.get("TABULAR_FORMAT", "CSV")
    plot_extension = su2io.get_extension(plot_format)
    history_filename = konfig["CONV_FILENAME"] + plot_extension
    special_cases = su2io.get_specialCases(konfig)

    # get history
    history = su2io.read_history(history_filename, config.NZONES)

    # update super config
    config.update(
        {
            "MATH_PROBLEM": konfig["MATH_PROBLEM"],
            "OBJECTIVE_FUNCTION": konfig["OBJECTIVE_FUNCTION"],
        }
    )

    # files out
    objective = konfig["OBJECTIVE_FUNCTION"]
    if "," in objective:
        objective = "COMBO"
    adj_title = "ADJOINT_" + objective
    suffix = su2io.get_adjointSuffix(objective)
    restart_name = konfig["RESTART_FILENAME"]
    restart_name = su2io.add_suffix(restart_name, suffix)

    # info out
    info = su2io.State()
    info.FILES[adj_title] = restart_name
    info.HISTORY[adj_title] = history

    return info
