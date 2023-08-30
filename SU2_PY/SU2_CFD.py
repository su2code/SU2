#!/usr/bin/env python

## \file SU2_CFD.py
#  \brief Python script to launch SU2_CFD through the Python Wrapper.
#  \author David Thomas
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

from __future__ import division, print_function, absolute_import
from optparse import OptionParser  # use a parser for configuration
import SU2  # imports SU2 python tools
import pysu2  # imports the SU2 wrapped module

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------


def main():

    # Command line options
    parser = OptionParser()
    parser.add_option(
        "-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE"
    )
    parser.add_option(
        "--nDim",
        dest="nDim",
        default=2,
        help="Define the number of DIMENSIONS",
        metavar="DIMENSIONS",
    )
    parser.add_option(
        "--nZone",
        dest="nZone",
        default=1,
        help="Define the number of ZONES",
        metavar="NZONE",
    )
    parser.add_option(
        "--parallel",
        action="store_true",
        help="Specify if we need to initialize MPI",
        dest="with_MPI",
        default=False,
    )
    parser.add_option(
        "--fsi",
        dest="fsi",
        default="False",
        help="Launch the FSI driver",
        metavar="FSI",
    )
    parser.add_option(
        "--fem",
        dest="fem",
        default="False",
        help="Launch the FEM driver (General driver)",
        metavar="FEM",
    )
    parser.add_option(
        "--harmonic_balance",
        dest="harmonic_balance",
        default="False",
        help="Launch the Harmonic Balance (HB) driver",
        metavar="HB",
    )
    parser.add_option(
        "--poisson_equation",
        dest="poisson_equation",
        default="False",
        help="Launch the poisson equation driver (General driver)",
        metavar="POIS_EQ",
    )
    parser.add_option(
        "--wave_equation",
        dest="wave_equation",
        default="False",
        help="Launch the wave equation driver (General driver)",
        metavar="WAVE_EQ",
    )
    parser.add_option(
        "--heat_equation",
        dest="heat_equation",
        default="False",
        help="Launch the heat equation driver (General driver)",
        metavar="HEAT_EQ",
    )

    (options, args) = parser.parse_args()
    options.nDim = int(options.nDim)
    options.nZone = int(options.nZone)
    options.fsi = options.fsi.upper() == "TRUE"
    options.fem = options.fem.upper() == "TRUE"
    options.harmonic_balance = options.harmonic_balance.upper() == "TRUE"
    options.poisson_equation = options.poisson_equation.upper() == "TRUE"
    options.wave_equation = options.wave_equation.upper() == "TRUE"
    options.heat_equation = options.heat_equation.upper() == "TRUE"

    if options.filename == None:
        raise Exception("No config file provided. Use -f flag")

    if options.with_MPI == True:
        from mpi4py import MPI  # use mpi4py for parallel run (also valid for serial)

        comm = MPI.COMM_WORLD
    else:
        comm = 0

    # Initialize the corresponding driver of SU2, this includes solver preprocessing
    try:
        if (options.nZone == 1) and (
            options.fem
            or options.poisson_equation
            or options.wave_equation
            or options.heat_equation
        ):
            SU2Driver = pysu2.CSinglezoneDriver(options.filename, options.nZone, comm)
        elif options.harmonic_balance:
            SU2Driver = pysu2.CHBDriver(options.filename, options.nZone, comm)
        elif options.nZone >= 2:
            SU2Driver = pysu2.CMultizoneDriver(options.filename, options.nZone, comm)
        else:
            SU2Driver = pysu2.CSinglezoneDriver(options.filename, options.nZone, comm)
    except TypeError as exception:
        print("A TypeError occured in pysu2.CDriver : ", exception)
        if options.with_MPI == True:
            print(
                "ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build."
            )
        else:
            print(
                "ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper."
            )
        return

    # Launch the solver for the entire computation
    SU2Driver.StartSolver()

    # Finalize the solver and exit cleanly
    SU2Driver.Finalize()

    if SU2Driver != None:
        del SU2Driver


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == "__main__":
    main()
