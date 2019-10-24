#!/usr/bin/env python

## \file fsi_computation.py
#  \brief Python wrapper code for adjoint computation by coupling pyBeam and SU2.
#  \author Ruben Sanchez based on work by David Thomas and Rocco Bombardieri
#  \version 7.0.0
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
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
import time as timer

from optparse import OptionParser  # use a parser for configuration

from SU2_FSI import FSI_config as io        # imports FSI config tools
from SU2_FSI import AdjointInterface as FSI # imports FSI python tools
import pyBeamInterface as pyBeamInterface
import pyMLSInterface as Spline_Module

# imports the CFD (SU2) module for FSI computation
import pysu2ad as pysu2
import pyBeam


# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():
    # --- Get the FSI config file name form the command line options --- #
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("--serial", action="store_true",
                      help="Specify if we need to initialize MPI", dest="serial", default=False)

    (options, args) = parser.parse_args()

    if options.serial:
        comm = 0
        myid = 0
        numberPart = 1
        have_MPI = False
    else:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        myid = comm.Get_rank()
        numberPart = comm.Get_size()
        have_MPI = True

    rootProcess = 0

    # --- Set the working directory --- #
    if myid == rootProcess:
        if os.getcwd() not in sys.path:
            sys.path.append(os.getcwd())
            print("Setting working directory : {}".format(os.getcwd()))
        else:
            print("Working directory is set to {}".format(os.getcwd()))

    # starts timer
    start = timer.time()

    confFile = str(options.filename)

    FSI_config = io.FSIConfig(confFile)  # FSI configuration file
    CFD_ConFile = FSI_config['SU2_CONFIG']  # CFD configuration file
    CSD_ConFile = FSI_config['PYBEAM_CONFIG']  # CSD configuration file
    MLS_confFile = FSI_config['MLS_CONFIG_FILE_NAME']  # MLS configuration file

    if have_MPI:
        comm.barrier()

    # --- Initialize the fluid solver: SU2 --- #
    if myid == rootProcess:
        print('\n***************************** Initializing SU2 **************************************')
    try:
        FluidSolver = pysu2.CDiscAdjSinglezoneDriver(CFD_ConFile, 1, comm)
    except TypeError as exception:
        print('A TypeError occured in pysu2.CSingleZoneDriver : ', exception)
        if serial:
            print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please remove the --serial option that is incompatible with a parallel build.')
        else:
            print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, add --serial to launch your simulation.')
        return

    if have_MPI:
        comm.barrier()

    # --- Initialize the solid solver: pyBeam --- #
    if myid == rootProcess:
        print('\n***************************** Initializing pyBeam ************************************')
        try:
            SolidSolver = pyBeamInterface.pyBeamADSolver(CSD_ConFile)
        except TypeError as exception:
            print('ERROR building the Solid Solver: ', exception)
    else:
        SolidSolver = None

    if have_MPI:
        comm.barrier()

    # --- Initialize and set the coupling environment --- #
    if myid == rootProcess:
        print('\n***************************** Initializing FSI interface *****************************')
    try:
        FSIInterface = FSI.AdjointInterface(FSI_config, FluidSolver, SolidSolver, None, have_MPI)
    except TypeError as exception:
        print('ERROR building the FSI Interface: ', exception)

    if have_MPI:
        comm.barrier()


    if myid == rootProcess:
        print('\n***************************** Connect fluid and solid solvers *****************************')
    try:
        FSIInterface.connect(FSI_config, FluidSolver, SolidSolver)
    except TypeError as exception:
        print('ERROR building the Interpolation Interface: ', exception)

    if have_MPI:
        comm.barrier()

    if myid == rootProcess:  # we perform this calculation on the root core
        print('\n***************************** Initializing MLS Interpolation *************************')
        try:
            MLS = Spline_Module.pyMLSInterface(MLS_confFile, FSIInterface.globalFluidCoordinates, 
                                               FSIInterface.globalSolidCoordinates)
        except TypeError as exception:
            print('ERROR building the MLS Interpolation: ', exception)

    else:
        MLS = None

    if have_MPI:
        comm.barrier()

    # Run the solver
    if myid == 0:
        print("\n------------------------------ Begin Solver -----------------------------\n")
    sys.stdout.flush()
    if have_MPI:
        comm.Barrier()

    FSIInterface.SteadyFSI(FSI_config, FluidSolver, SolidSolver, MLS)

    # Postprocess the solver and exit cleanly
    FluidSolver.Postprocessing()

    if FluidSolver is not None:
        del FluidSolver


    return


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
