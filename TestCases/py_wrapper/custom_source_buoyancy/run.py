#!/usr/bin/env python

## \file run.py
#  \brief Buoyancy force using user defines source term
#  \version 8.1.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

# with mpi:
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
# without mpi:
#  comm = 0

def main():


  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2.CSinglezoneDriver('lam_buoyancy_cavity.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CSinglezoneDriver : ', exception)
    raise

  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------")
    sys.stdout.flush()

  # we need to add a source term to the energy equation. For this, we need to get the solver and the variable first.
  # we then loop over all points and for these points, we add the source term
  nDim = driver.GetNumberDimensions()

  # index to the flow solver
  iSOLVER = driver.GetSolverIndices()['INC.FLOW']

  # all the indices and the map to the names of the primitives
  primindex = driver.GetPrimitiveIndices()

  iDENSITY = primindex.get("DENSITY")

  Body_Force_Vector = [0.0, -9.81, 0.0]
  DensityInc_0 = driver.GetDensityFreeStreamND()
  Force_Ref = driver.GetForceRef()

  # super important to actually push the commands.
  sys.stdout.flush()

  # run N iterations
  for inner_iter in range(2):
    if (rank==0):
      print("python iteration ", inner_iter)

    # set the source term, per point
    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
      DensityInc_i = driver.Primitives()(i_node,iDENSITY)

      for iDim in range(nDim):
        custom_source_vector = (DensityInc_i - DensityInc_0) * Body_Force_Vector[iDim] / Force_Ref
        driver.UserDefinedSource(iSOLVER).Set(i_node,iDim+1,custom_source_vector)

    driver.Preprocess(inner_iter)
    driver.Run()

    driver.Postprocess()
    driver.Update()

    # Monitor the solver and output solution to file if required.
    driver.Output(inner_iter)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

if __name__ == '__main__':
  main()
