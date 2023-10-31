#!/usr/bin/env python

## \file fsi_computation.py
#  \brief Python wrapper code for FSI computation by coupling a third-party structural solver to SU2.
#  \authors Nicola Fonzi, Vittorio Cavalieri based on the work of David Thomas
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

import os
import sys
import shutil
import copy
import time as timer
from math import *  # use mathematical expressions
from optparse import OptionParser  # use a parser for configuration

# imports the CFD (SU2) module for FSI computation
import pysu2
import FSI_tools as FSI  # imports FSI python tools

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------


def main():

    # --- Get the FSI conig file name form the command line options --- #
    parser = OptionParser()
    parser.add_option(
        "-f", "--file", dest="filename", help="read config from FILE", metavar="FILE"
    )
    parser.add_option(
        "--parallel",
        action="store_true",
        help="Specify if we need to initialize MPI",
        dest="with_MPI",
        default=False,
    )

    (options, args) = parser.parse_args()

    if options.with_MPI:
        from mpi4py import (
            MPI,
        )  # MPI is initialized from now by python and can be continued in C++

        comm = MPI.COMM_WORLD
        myid = comm.Get_rank()
        numberPart = comm.Get_size()
        have_MPI = True
    else:
        comm = 0
        myid = 0
        numberPart = 1
        have_MPI = False

    rootProcess = 0

    # --- Set the working directory --- #
    if myid == rootProcess:
        if os.getcwd() not in sys.path:
            sys.path.append(os.getcwd())
            print("Setting working directory : {}".format(os.getcwd()))
        else:
            print("Working directory is set to {}".format(os.getcwd()))

    if have_MPI:
        comm.barrier()

    # starts timer
    start = timer.time()

    confFile = str(options.filename)

    FSI_config = FSI.FSIConfig(confFile, comm)  # FSI configuration file
    CFD_ConFile = FSI_config["CFD_CONFIG_FILE_NAME"]  # CFD configuration file
    CSD_ConFile = FSI_config["CSD_CONFIG_FILE_NAME"]  # CSD configuration file

    CSD_Solver = FSI_config["CSD_SOLVER"]  # CSD solver

    if have_MPI:
        comm.barrier()

    # --- Initialize the fluid solver --- #
    if myid == rootProcess:
        print("\n")
        print(" Initializing fluid solver ".center(80, "*"))
    try:
        FluidSolver = pysu2.CSinglezoneDriver(CFD_ConFile, 1, comm)
    except TypeError as exception:
        print("A TypeError occured in pysu2.CSinglezoneDriver : ", exception)
        if have_MPI:
            print(
                "ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build."
            )
        else:
            print(
                "ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper."
            )
        return

    if have_MPI:
        comm.barrier()

    # --- Initialize the solid solver --- #
    # Serial solvers
    if CSD_Solver in ["NATIVE"]:
        if myid == rootProcess:
            print("\n")
            print(" Initializing solid solver ".center(80, "*"))
            if CSD_Solver == "NATIVE":
                from SU2_Nastran import pysu2_nastran

                if FSI_config["IMPOSED_MOTION"] == "NO":
                    SolidSolver = pysu2_nastran.Solver(CSD_ConFile, False)
                else:
                    SolidSolver = pysu2_nastran.Solver(CSD_ConFile, True)
        else:
            SolidSolver = None
    # Parallel solvers
    # For now we are only using serial solvers
    else:
        raise Exception("\n Invalid solid solver option")

    if have_MPI:
        comm.barrier()

    # --- Initialize and set the FSI interface (coupling environement) --- #
    if myid == rootProcess:
        print("\n")
        print(" Initializing FSI interface ".center(80, "*"))
    if have_MPI:
        comm.barrier()
    FSIInterface = FSI.Interface(FSI_config, FluidSolver, SolidSolver, have_MPI)

    if myid == rootProcess:
        print("\n")
        print(" Connect fluid and solid solvers ".center(80, "*"))
    if have_MPI:
        comm.barrier()
    FSIInterface.connect(FSI_config, FluidSolver, SolidSolver)

    if myid == rootProcess:
        print("\n")
        print(" Mapping fluid-solid interfaces ".center(80, "*"))
    if have_MPI:
        comm.barrier()
    FSIInterface.interfaceMapping(FluidSolver, SolidSolver, FSI_config)

    if have_MPI:
        comm.barrier()

    if FSI_config["MAPPING_MODES"] == "NO":
        # --- Launch a steady or unsteady FSI computation --- #
        if FSI_config["TIME_MARCHING"] == "YES":
            try:
                FSIInterface.UnsteadyFSI(FSI_config, FluidSolver, SolidSolver)
            except NameError as exception:
                if myid == rootProcess:
                    print(
                        "An NameError occured in FSIInterface.UnsteadyFSI : ", exception
                    )
            except TypeError as exception:
                if myid == rootProcess:
                    print(
                        "A TypeError occured in FSIInterface.UnsteadyFSI : ", exception
                    )
            except KeyboardInterrupt as exception:
                if myid == rootProcess:
                    print(
                        "A KeyboardInterrupt occured in FSIInterface.UnsteadyFSI : ",
                        exception,
                    )
        else:
            try:
                FSIInterface.SteadyFSI(FSI_config, FluidSolver, SolidSolver)
            except NameError as exception:
                if myid == rootProcess:
                    print(
                        "An NameError occured in FSIInterface.SteadyFSI : ", exception
                    )
            except TypeError as exception:
                if myid == rootProcess:
                    print("A TypeError occured in FSIInterface.SteadyFSI : ", exception)
            except KeyboardInterrupt as exception:
                if myid == rootProcess:
                    print(
                        "A KeyboardInterrupt occured in FSIInterface.SteadyFSI : ",
                        exception,
                    )
    else:
        try:
            FSIInterface.MapModes(FSI_config, FluidSolver, SolidSolver)
        except NameError as exception:
            if myid == rootProcess:
                print("An NameError occured in FSIInterface.MapModes : ", exception)
        except TypeError as exception:
            if myid == rootProcess:
                print("A TypeError occured in FSIInterface.MapModes : ", exception)
        except KeyboardInterrupt as exception:
            if myid == rootProcess:
                print(
                    "A KeyboardInterrupt occured in FSIInterface.MapModes : ", exception
                )

    if have_MPI:
        comm.barrier()

    # --- Exit cleanly the fluid and solid solvers --- #
    FluidSolver.Finalize()
    if myid == rootProcess:
        SolidSolver.exit()

    if have_MPI:
        comm.barrier()

    # stops timer
    stop = timer.time()
    elapsedTime = stop - start

    if myid == rootProcess:
        print(
            "\n Computation successfully performed in {} seconds.".format(elapsedTime)
        )

    return


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == "__main__":
    main()
