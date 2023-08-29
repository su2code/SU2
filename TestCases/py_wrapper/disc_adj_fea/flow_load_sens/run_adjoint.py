#!/usr/bin/env python

## \file run_adjoint.py
#  \brief Python script to launch SU2_CFD_AD and compute the sensitivity of the FEA problem respect to flow loads.
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
  else:
    comm = 0

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  SU2Driver = pysu2.CDiscAdjSinglezoneDriver(options.filename, options.nZone, comm);

  MarkerID = None
  MarkerName = 'RightBeamS'       # Specified by the user

  # Get all the boundary tags
  MarkerList = SU2Driver.GetMarkerTags()

  # Get all the markers defined on this rank and their associated indices.
  allMarkerIDs = SU2Driver.GetMarkerIndices()

  #Check if the specified marker exists and if it belongs to this rank.
  if MarkerName in MarkerList and MarkerName in allMarkerIDs.keys():
    MarkerID = allMarkerIDs[MarkerName]

  # Only print on the rank to which the marker belongs.
  # WARNING: We only do this for the regression test, there is no guarantee that a marker will only belong to one rank.

  # Time loop is defined in Python so that we have acces to SU2 functionalities at each time step
  if MarkerID != None:
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  # Define the load at the target vertex
  if MarkerID != None:
    SU2Driver.SetMarkerCustomFEALoad(MarkerID, 5, (0, -0.005, 0))

  # Time iteration preprocessing
  SU2Driver.Preprocess(0)

  # Run one time-step (static: one simulation)
  SU2Driver.Run()

  # Update the solver for the next time iteration
  SU2Driver.Update()

  # Monitor the solver and output solution to file if required
  SU2Driver.Monitor(0)

  # Output the solution to file
  SU2Driver.Output(0)

  if MarkerID != None:
    sens=[]
    disp=[]

    # Recover the sensitivity
    sens.append(SU2Driver.GetMarkerFEALoadSensitivity(MarkerID, 5))

    fea_sol = SU2Driver.GetSolverIndices()["FEA"]
    marker_disp = SU2Driver.MarkerSolution(fea_sol, MarkerID)
    disp.append(marker_disp.Get(5))

    print("Sens[0]\tSens[1]\tDisp[0]\tDisp[1]\t")
    print(100, 100, sens[0][0], sens[0][1], disp[0][0], disp[0][1])

  # Finalize the solver and exit cleanly
  SU2Driver.Finalize()



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
