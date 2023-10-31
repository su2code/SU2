#!/usr/bin/env python

## \file flatPlate_rigidMotion.py
#  \brief Python script to launch SU2_CFD with customized unsteady boundary conditions using the Python wrapper.
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

import sys
from optparse import OptionParser	# use a parser for configuration
import pysu2			            # imports the SU2 wrapped module
from math import *
import numpy as np

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

  # Command line options
  parser=OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
  parser.add_option("--parallel", action="store_true",
                    help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  (options, args) = parser.parse_args()
  options.nDim  = int(2)
  options.nZone = int(1)

  # Import mpi4py for parallel run
  if options.with_MPI == True:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
  else:
    comm = 0
    rank = 0

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  try:
      SU2Driver = pysu2.CSinglezoneDriver(options.filename, options.nZone, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return


  MovingMarkerID = None
  MovingMarker = 'plate'       #specified by the user

  # Get all the tags with the moving option
  MovingMarkerList =  SU2Driver.GetMarkerTags()

  # Get all the markers defined on this rank and their associated indices.
  allMarkerIDs = SU2Driver.GetMarkerIndices()

  # Check if the specified marker has a moving option and if it exists on this rank.
  if MovingMarker in MovingMarkerList and MovingMarker in allMarkerIDs.keys():
    MovingMarkerID = allMarkerIDs[MovingMarker]

  # Number of vertices on the specified marker (per rank)
  nVertex_MovingMarker = 0         #total number of vertices (physical + halo)

  if MovingMarkerID != None:
    nVertex_MovingMarker = SU2Driver.GetNumberMarkerNodes(MovingMarkerID)\

  # Retrieve some control parameters from the driver
  deltaT = SU2Driver.GetUnsteadyTimeStep()
  TimeIter = SU2Driver.GetTimeIter()
  nTimeIter = SU2Driver.GetNumberTimeIter()
  time = TimeIter*deltaT

  # Extract the initial position of each node on the moving marker
  CoordX = np.zeros(nVertex_MovingMarker)
  CoordY = np.zeros(nVertex_MovingMarker)
  for iVertex in range(nVertex_MovingMarker):
    CoordX[iVertex], CoordY[iVertex] = SU2Driver.MarkerInitialCoordinates(MovingMarkerID).Get(iVertex)

  # Time loop is defined in Python so that we have acces to SU2 functionalities at each time step
  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  while (TimeIter < nTimeIter):
    # Define the rigid body displacement and set the new coords of each node on the marker
    value = 0.0, 0.0175*sin(2*pi*time)

    for iVertex in range(nVertex_MovingMarker):
      SU2Driver.SetMarkerCustomDisplacement(MovingMarkerID, int(iVertex), value)

    # Time iteration preprocessing
    SU2Driver.Preprocess(TimeIter)
    # Run one time iteration (e.g. dual-time)
    SU2Driver.Run()
    # Postprocess the solver
    SU2Driver.Postprocess()
    # Update the solver for the next time iteration
    SU2Driver.Update()
    # Monitor the solver and output solution to file if required
    stopCalc = SU2Driver.Monitor(TimeIter)
    SU2Driver.Output(TimeIter)
    if (stopCalc == True):
      break
    # Update control parameters
    TimeIter += 1
    time += deltaT

  # Finalize the solver and exit cleanly
  SU2Driver.Finalize()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
