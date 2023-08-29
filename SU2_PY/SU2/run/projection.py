#!/usr/bin/env python

## \file projection.py
#  \brief python package for running gradient projection
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

from .. import io as su2io
from .. import util as su2util
from .interface import DOT as SU2_DOT


# ----------------------------------------------------------------------
#  Gradient Projection
# ----------------------------------------------------------------------


def projection(config, state={}, step=1e-3):
    """info = SU2.run.projection(config,state,step=1e-3)

    Runs an gradient projection with:
        SU2.run.decomp()
        SU2.run.DOT()

    Assumptions:
        Writes tecplot file of gradients
        Adds objective suffix to gradient plot filename

    Inputs:
        config - an SU2 config
        state  - only required when using external custom DV
        step   - a float or list of floats for geometry sensitivity
                 finite difference step

    Outputs:
        info - SU2 State with keys:
            GRADIENTS.<config.OBJECTIVE_FUNCTION>

    Updates:
        config.MATH_PROBLEM

    Executes in:
        ./
    """
    # local copy
    konfig = copy.deepcopy(config)

    # choose dv values
    Definition_DV = konfig["DEFINITION_DV"]
    n_DV = sum(Definition_DV["SIZE"])
    if isinstance(step, list):
        assert len(step) == n_DV, "unexpected step vector length"
    else:
        step = [step] * n_DV
    dv_old = [
        0.0
    ] * n_DV  # SU2_DOT input requirement, assumes linear superposition of design variables
    dv_new = step
    konfig.unpack_dvs(dv_new, dv_old)

    # filenames
    objective = konfig["OBJECTIVE_FUNCTION"]
    grad_filename = konfig["GRAD_OBJFUNC_FILENAME"]
    output_format = konfig.get("TABULAR_FORMAT", "CSV")
    plot_extension = su2io.get_extension(output_format)
    adj_suffix = su2io.get_adjointSuffix(objective)
    grad_plotname = (
        os.path.splitext(grad_filename)[0] + "_" + adj_suffix + plot_extension
    )

    # Run Projection
    SU2_DOT(konfig)

    # read raw gradients
    raw_gradients = su2io.read_gradients(grad_filename)
    os.remove(grad_filename)

    info = su2io.State()

    # Write Gradients
    data_plot = su2util.ordered_bunch()
    data_plot["VARIABLE"] = range(len(raw_gradients))
    data_plot["GRADIENT"] = raw_gradients
    data_plot["FINDIFF_STEP"] = step
    su2util.write_plot(grad_plotname, output_format, data_plot)

    # gradient output dictionary
    objective = objective.split(",")
    if len(objective) > 1:
        objective = ["COMBO"]

    gradients = {objective[0]: raw_gradients}

    # info out
    info.GRADIENTS.update(gradients)

    return info
