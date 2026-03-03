#!/usr/bin/env python

## \file run.py
#  \brief Unsteady inlet boundary conditions.
#  \version 8.2.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

import sys
import pysu2
import math
# from mpi4py import MPI

def main():
  """
  Run the flow solver with a custom inlet (function of time and space).
  """
  # comm = MPI.COMM_WORLD
  comm = 0

  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2.CSinglezoneDriver('lam_flatplate.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CSinglezoneDriver : ', exception)
    raise

  # Get the ID of the inlet marker.
  all_marker_ids = driver.GetMarkerIndices()
  marker_name = 'x_minus'
  marker_id = all_marker_ids[marker_name] if marker_name in all_marker_ids else -1

  # Run the time loop in python to vary the inlet conditions.
  dt = driver.GetUnsteadyTimeStep()

  print("\n------------------------------ Begin Solver -----------------------------")
  sys.stdout.flush()

  for time_iter in range(driver.GetNumberTimeIter()):
    # Change the total pressure.
    if marker_id >= 0:
      for i_vertex in range(driver.GetNumberMarkerNodes(marker_id)):
        y = driver.MarkerCoordinates(marker_id)(i_vertex, 1)
        t = time_iter * dt
        pt = 1e5 + 2e4 * y / 0.01 * (1 - math.cos(2 * math.pi * t / 0.1))
        driver.SetMarkerCustomInletFlowVar1(marker_id, i_vertex, pt)
    driver.BoundaryConditionsUpdate()

    driver.Preprocess(time_iter)

    # Run one time iteration.
    driver.Run()
    driver.Postprocess()
    driver.Update()

    # Monitor the solver and output solution to file if required.
    driver.Monitor(time_iter)
    driver.Output(time_iter)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

if __name__ == '__main__':
  main()
