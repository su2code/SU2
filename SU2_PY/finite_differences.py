#!/usr/bin/env python

## \file finite_differences.py
#  \brief Python script for doing the finite differences computation using the SU2 suite.
#  \author F. Palacios
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

import os, sys
from optparse import OptionParser

sys.path.append(os.environ["SU2_RUN"])
import SU2

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------


def main():

    parser = OptionParser()
    parser.add_option(
        "-f", "--file", dest="filename", help="read config from FILE", metavar="FILE"
    )
    parser.add_option(
        "-n",
        "--partitions",
        dest="partitions",
        default=1,
        help="number of PARTITIONS",
        metavar="PARTITIONS",
    )
    parser.add_option(
        "-q",
        "--quiet",
        dest="quiet",
        default="False",
        help="output QUIET to log files",
        metavar="QUIET",
    )
    parser.add_option(
        "-z",
        "--zones",
        dest="nzones",
        default="1",
        help="Number of Zones",
        metavar="ZONES",
    )

    (options, args) = parser.parse_args()
    options.partitions = int(options.partitions)
    options.quiet = options.quiet.upper() == "TRUE"
    options.nzones = int(options.nzones)

    finite_differences(
        options.filename, options.partitions, options.quiet, options.nzones
    )


#: def main()


# -------------------------------------------------------------------
#  Finite Differences Function
# -------------------------------------------------------------------


def finite_differences(filename, partitions=0, quiet=False, nzones=1):
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.NZONES = int(nzones)

    if quiet:
        config.CONSOLE = "CONCISE"

    # State
    state = SU2.io.State()
    state.find_files(config)

    # add restart files to state.FILES
    if (
        config.get("TIME_DOMAIN", "NO") == "YES"
        and config.get("RESTART_SOL", "NO") == "YES"
    ):
        restart_name = config["RESTART_FILENAME"].split(".")[0]
        restart_filename = (
            restart_name + "_" + str(int(config["RESTART_ITER"]) - 1).zfill(5) + ".dat"
        )
        if not os.path.isfile(
            restart_filename
        ):  # throw, if restart files does not exist
            sys.exit("Error: Restart file <" + restart_filename + "> not found.")
        state["FILES"]["RESTART_FILE_1"] = restart_filename

        # use only, if time integration is second order
        if config.get("TIME_MARCHING", "NO") == "DUAL_TIME_STEPPING-2ND_ORDER":
            restart_filename = (
                restart_name
                + "_"
                + str(int(config["RESTART_ITER"]) - 2).zfill(5)
                + ".dat"
            )
            if not os.path.isfile(
                restart_filename
            ):  # throw, if restart files does not exist
                sys.exit("Error: Restart file <" + restart_filename + "> not found.")
            state["FILES"]["RESTART_FILE_2"] = restart_filename

    # Finite Difference Gradients
    SU2.eval.gradients.findiff(config, state)

    return state


#: finite_differences()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == "__main__":
    main()
