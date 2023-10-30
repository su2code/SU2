#!/usr/bin/env python

## \file run.py
#  \brief Deforming bump in channel.
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
import numpy as np

def main():
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  # Initialize the corresponding driver of SU2, this includes solver preprocessing.
  try:
    SU2Driver = pysu2.CSinglezoneDriver('config.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ', exception)
    raise

  # Get the ID of the marker we want to deform.
  AllMarkerIDs = SU2Driver.GetMarkerIndices()
  MarkerName = 'interface'
  MarkerID = AllMarkerIDs[MarkerName] if MarkerName in AllMarkerIDs else -1

  # Number of vertices on the specified marker (per rank).
  nVertex = SU2Driver.GetNumberMarkerNodes(MarkerID) if MarkerID >= 0 else 0

  # Retrieve some control parameters from the driver.
  deltaT = SU2Driver.GetUnsteadyTimeStep()
  TimeIter = SU2Driver.GetTimeIter()
  nTimeIter = SU2Driver.GetNumberTimeIter()
  time = TimeIter * deltaT

  # Extract the initial position of each node on the moving marker.
  CoordX = np.zeros(nVertex)
  CoordY = np.zeros(nVertex)
  for iVertex in range(nVertex):
    CoordX[iVertex], CoordY[iVertex] = SU2Driver.MarkerInitialCoordinates(MarkerID).Get(iVertex)

  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------\n")

  # The time loop is defined in Python so that we have acces to SU2 functionalities at each time step.
  while (TimeIter < nTimeIter):
    # Apply the surface deformation.
    for iVertex in range(nVertex):
      dy = np.real(DeformFunction(CoordX[iVertex] - 0.9, time))
      SU2Driver.SetMarkerCustomDisplacement(MarkerID, iVertex, (0.0, dy))

    # Time iteration preprocessing.
    SU2Driver.Preprocess(TimeIter)

    # Run one time iteration (e.g. dual-time).
    SU2Driver.Run()
    SU2Driver.Postprocess()
    SU2Driver.Update()

    # Monitor the solver and output solution to file if required.
    stopCalc = SU2Driver.Monitor(TimeIter)
    SU2Driver.Output(TimeIter)

    if (stopCalc == True):
      break

    # Update control parameters
    TimeIter += 1
    time += deltaT

  # Check the value of an output to cover the functionality in a regression test.
  assert 'DRAG' in SU2Driver.GetOutputNames()
  assert abs(SU2Driver.GetOutputValue('DRAG') -
             SU2Driver.GetMarkerMonitoringOutputValue('DRAG_ON_SURFACE', MarkerName)) < np.finfo(float).eps

  # Finalize the solver and exit cleanly
  SU2Driver.Finalize()


# Imposed deformation
def DeformFunction(x,t):
  A1 = 0.016
  L = 1.0
  k1 = 4.730/L
  return A1*(np.cosh(k1*x) - np.cos(k1*x) + ((np.cosh(k1*L)-np.cos(k1*L))*(np.sin(k1*x)-np.sinh(k1*x)))/(np.sinh(k1*L) - np.sin(k1*L)))*(np.exp(1j*t*2) + np.exp(-1j*t*2))


if __name__ == '__main__':
  main()

