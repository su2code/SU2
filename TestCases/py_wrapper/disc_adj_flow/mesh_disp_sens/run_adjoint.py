#!/usr/bin/env python

## \file run_adjoint.py
#  \brief Python script to launch SU2_CFD_AD
#  \author Ruben Sanchez
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
import pysu2ad as pysu2          # imports the SU2 adjoint-wrapped module

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
  options.nDim  = 2
  options.nZone = 1

  print(args)

  # Import mpi4py for parallel run
  if options.with_MPI == True:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
  else:
    comm = 0
    rank = 0

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  SU2Driver = pysu2.CDiscAdjSinglezoneDriver(options.filename, options.nZone, comm);

  MarkerID = None
  MarkerName = 'wallF'       # Specified by the user

  # Get all the boundary tags
  MarkerList = SU2Driver.GetMarkerTags()

  # Get all the markers defined on this rank and their associated indices.
  allMarkerIDs = SU2Driver.GetMarkerIndices()

  # Check if the specified marker exists and if it belongs to this rank.
  if MarkerName in MarkerList and MarkerName in allMarkerIDs.keys():
    MarkerID = allMarkerIDs[MarkerName]

  # Number of vertices on the specified marker (per rank)
  nVertex_Marker = 0         #total number of vertices (physical + halo)

  if MarkerID != None:
    nVertex_Marker = SU2Driver.GetNumberMarkerNodes(MarkerID)

  # Time loop is defined in Python so that we have acces to SU2 functionalities at each time step
  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  # Time iteration preprocessing
  SU2Driver.Preprocess(0)

  # Run one time-step (static: one simulation)
  SU2Driver.Run()

  # Postprocess
  SU2Driver.Postprocess()

  # Update the solver for the next time iteration
  SU2Driver.Update()

  # Monitor the solver and output solution to file if required
  SU2Driver.Monitor(0)

  # Output the solution to file
  SU2Driver.Output(0)

  # Sensitivities of the marker
  print("\n------------------------------ Sensitivities -----------------------------\n")
  for iVertex in range(nVertex_Marker):
    sensX, sensY = SU2Driver.GetMarkerDisplacementSensitivity(MarkerID, iVertex)

    if (iVertex == 30) and rank == 0:
      print(1000, 1000, iVertex, sensX, sensY, 0.0)

  # Finalize the solver and exit cleanly
  SU2Driver.Finalize()



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
