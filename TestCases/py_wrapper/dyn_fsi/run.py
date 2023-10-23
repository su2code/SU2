#!/usr/bin/env python

## \file run.py
#  \brief Unsteady FSI case with custom load.
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

import sys
import pysu2
import math
from mpi4py import MPI

def main():
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  # Initialize the corresponding driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2.CMultizoneDriver('config.cfg', 2, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ', exception)
    raise

  # By default the first zone is selected (flow in the case).
  # Select the structural zone (second zone) to work with the FEA solver.
  driver.SelectZone(1)

  # Get the ID of the marker where we will apply a force.
  # This marker cannot be used by the fluid-solid interface otherwise the imposed
  # load will be cleared when interpolating forces.
  all_marker_ids = driver.GetMarkerIndices()
  ctrl_id = all_marker_ids['internal'] if 'internal' in all_marker_ids else -1

  # Number of vertices on the specified markers (per rank).
  n_vertex_ctrl = driver.GetNumberMarkerNodes(ctrl_id) if ctrl_id >= 0 else 0

  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------")
    sys.stdout.flush()

  for time_iter in range(driver.GetNumberTimeIter()):
    # Apply a custom load and then solve the time step.
    time = time_iter * driver.GetUnsteadyTimeStep()
    for i_vertex in range(n_vertex_ctrl):
      i_point = driver.GetMarkerNode(ctrl_id, i_vertex)
      if driver.GetNodeDomain(i_point):
        driver.SetMarkerCustomFEALoad(ctrl_id, i_vertex, (-0.002 + 0.002 * math.cos(2 * math.pi * time / 0.02), 0))

    driver.Preprocess(time_iter)

    driver.Run()
    driver.Postprocess()
    driver.Update()

    stop = driver.Monitor(time_iter)
    driver.Output(time_iter)
    if stop: break

  # Finalize the solver and exit cleanly.
  driver.Finalize()


if __name__ == '__main__':
  main()
