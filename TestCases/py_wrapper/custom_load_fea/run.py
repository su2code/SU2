#!/usr/bin/env python

## \file run.py
#  \brief FEA case with custom load.
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

import pysu2
from mpi4py import MPI

def main():
  comm = MPI.COMM_WORLD

  # Initialize the corresponding driver of SU2, this includes solver preprocessing.
  try:
    SU2Driver = pysu2.CSinglezoneDriver('config.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ', exception)
    raise

  # Get the ID of the marker we want to deform.
  AllMarkerIDs = SU2Driver.GetMarkerIndices()
  MarkerName = 'y_minus'
  MarkerID = AllMarkerIDs[MarkerName] if MarkerName in AllMarkerIDs else -1

  # Number of vertices on the specified marker (per rank).
  nVertex = SU2Driver.GetNumberMarkerNodes(MarkerID) if MarkerID >= 0 else 0

  # Apply a load based on the coordinates.
  if nVertex > 0:
    MarkerCoords = SU2Driver.MarkerCoordinates(MarkerID)
    L = 0.5
    dx = L / 16 # known from mesh settings in this case.
    for iVertex in range(nVertex):
      x = MarkerCoords(iVertex, 0)
      nodalForce = (2 * x / L) * dx
      # Half load due to half dx on first and last node.
      if abs(x) < 1e-6 or abs(x - L) < 1e-6:
        nodalForce = nodalForce / 2
      SU2Driver.SetMarkerCustomFEALoad(MarkerID, iVertex, (0, nodalForce))

  # Solve.
  SU2Driver.StartSolver()

  # Find the tip displacement.
  MarkerName = 'x_plus'
  MarkerID = AllMarkerIDs[MarkerName] if MarkerName in AllMarkerIDs else -1
  nVertex = SU2Driver.GetNumberMarkerNodes(MarkerID) if MarkerID >= 0 else 0
  Disp = 0
  NodeFound = False

  if nVertex > 0:
    MarkerCoords = SU2Driver.MarkerCoordinates(MarkerID)
    SolverID = SU2Driver.GetSolverIndices()["FEA"]
    Solution = SU2Driver.MarkerSolution(SolverID, MarkerID)
    DispID = SU2Driver.GetFEASolutionIndices()["DISPLACEMENT_Y"]

    for iVertex in range(nVertex):
      y = MarkerCoords(iVertex, 1)
      if abs(y - 0.025) < 1e-6:
        Disp = Solution(iVertex, DispID)
        NodeFound = True

  if NodeFound:
    print(f"Vertical displacement of tip: {Disp}")
    # Test the value against expected.
    assert abs(Disp / 0.095439 - 1) < 1e-5, "Test FAILED"

  # Finalize the solver and exit cleanly.
  SU2Driver.Finalize()


if __name__ == '__main__':
  main()
