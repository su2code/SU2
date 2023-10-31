#!/usr/bin/env python

## \file direct.py
#  \brief python package for running direct solutions
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
#  Direct Simulation
# ----------------------------------------------------------------------


def direct(config):
    """info = SU2.run.direct(config)

    Runs an adjoint analysis with:
        SU2.run.decomp()
        SU2.run.CFD()
        SU2.run.merge()

    Assumptions:
        Does not rename restart filename to solution filename
        Adds 'direct' suffix to convergence filename

    Outputs:
        info - SU2 State with keys:
            FUNCTIONS
            HISTORY.DIRECT
            FILES.DIRECT

    Updates:
        config.MATH_PROBLEM

    Executes in:
        ./
    """

    # local copy
    konfig = copy.deepcopy(config)

    # setup direct problem
    konfig["MATH_PROBLEM"] = "DIRECT"
    konfig["CONV_FILENAME"] = konfig["CONV_FILENAME"] + "_direct"

    direct_diff = konfig.get("DIRECT_DIFF", "NO") == "YES"

    # Run Solution
    SU2_CFD(konfig)

    # multizone cases
    multizone_cases = su2io.get_multizone(konfig)

    # merge
    konfig["SOLUTION_FILENAME"] = konfig["RESTART_FILENAME"]
    if "FLUID_STRUCTURE_INTERACTION" in multizone_cases:
        konfig["SOLUTION_FILENAME"] = konfig["RESTART_FILENAME"]

    # filenames
    plot_format = konfig.get("TABULAR_FORMAT", "CSV")
    plot_extension = su2io.get_extension(plot_format)

    # adapt the history_filename, if a restart solution is chosen
    # check for 'RESTART_ITER' is to avoid forced restart situation in "compute_polar.py"...
    if konfig.get("RESTART_SOL", "NO") == "YES" and konfig.get("RESTART_ITER", 1) != 1:
        if konfig.get("CONFIG_LIST", []) != []:
            konfig[
                "CONV_FILENAME"
            ] = "config_CFD"  # master cfg is always config_CFD. Hardcoded names are prob nt ideal.
        restart_iter = "_" + str(konfig["RESTART_ITER"]).zfill(5)
        history_filename = konfig["CONV_FILENAME"] + restart_iter + plot_extension
    else:
        if konfig.get("CONFIG_LIST", []) != []:
            konfig["CONV_FILENAME"] = "config_CFD"
        history_filename = konfig["CONV_FILENAME"] + plot_extension

    special_cases = su2io.get_specialCases(konfig)

    # averaging final iterations
    final_avg = config.get("ITER_AVERAGE_OBJ", 0)
    # get chosen windowing function, default is square
    wnd_fct = config.get("WINDOW_FUNCTION", "SQUARE")

    # get history and objectives
    history = su2io.read_history(history_filename, config.NZONES)
    aerodynamics = su2io.read_aerodynamics(
        history_filename, config.NZONES, special_cases, final_avg, wnd_fct
    )

    # update super config
    config.update({"MATH_PROBLEM": konfig["MATH_PROBLEM"]})

    # info out
    info = su2io.State()
    info.FUNCTIONS.update(aerodynamics)
    info.FILES.DIRECT = konfig["RESTART_FILENAME"]
    if "INV_DESIGN_CP" in special_cases:
        info.FILES.TARGET_CP = "TargetCp.dat"
    if "INV_DESIGN_HEATFLUX" in special_cases:
        info.FILES.TARGET_HEATFLUX = "TargetHeatFlux.dat"
    info.HISTORY.DIRECT = history

    """If WINDOW_CAUCHY_CRIT is activated and the time marching converged before the final time has been reached,
       store the information for the adjoint run"""
    if config.get("WINDOW_CAUCHY_CRIT", "NO") == "YES" and config.TIME_MARCHING != "NO":
        konfig["TIME_ITER"] = int(
            info.HISTORY.DIRECT.Time_Iter[-1] + 1
        )  # update the last iteration
        if konfig["UNST_ADJOINT_ITER"] > konfig["TIME_ITER"]:
            konfig["ITER_AVERAGE_OBJ"] = max(
                0,
                konfig["ITER_AVERAGE_OBJ"]
                - (konfig["UNST_ADJOINT_ITER"] - konfig["TIME_ITER"]),
            )
            konfig["UNST_ADJOINT_ITER"] = konfig["TIME_ITER"]

        info["WND_CAUCHY_DATA"] = {
            "TIME_ITER": konfig["TIME_ITER"],
            "UNST_ADJOINT_ITER": konfig["UNST_ADJOINT_ITER"],
            "ITER_AVERAGE_OBJ": konfig["ITER_AVERAGE_OBJ"],
        }

    su2merge(konfig)

    return info
