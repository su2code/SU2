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
from mpi4py import MPI


def main():
  """
  custom source to add buoyancy term.
  """
  # parallel 
  comm = MPI.COMM_WORLD
  # serial
  #comm = 0

  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2.CSinglezoneDriver('lam_buoyancy_cavity.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CSinglezoneDriver : ', exception)
    raise

  print("\n------------------------------ Begin Solver -----------------------------")
  sys.stdout.flush()

  # we need to add a source term to the energy equation. For this, we need to get the solver and the variable first.
  # we then loop over all points and for these points, we add the source term
  nDim = driver.GetNumberDimensions()

  # index to the flow solver
  iSOLVER = driver.GetSolverIndices()['INC.FLOW']
  #print("index of flow solver = ",iSOLVER)

  # all the indices and the map to the names of the primitives
  primindex = driver.GetPrimitiveIndices()
  #print("indices of primitives=",primindex)
  #print("number of primitives:",len(primindex))

  #print("number of elements:",driver.GetNumberElements())

  nVars = driver.GetNumberSolverVars(iSOLVER)
  #print("number of solver variables:",nVars)
  varindex = primindex.copy()
  for prim in varindex.copy():
    if varindex[prim] >=nVars:
      del varindex[prim]
  varindex = dict(sorted(varindex.items(), key=lambda item: item[1]))



  print("solver variable names:",varindex)
  iDENSITY = primindex.get("DENSITY")
  print("index of density = ",iDENSITY)

  index_Vel = varindex.get("VELOCITY_X")
  print("index of velocity = ",index_Vel)
  custom_source_vector = [0.0 for i in range(nVars)]
  print("custom source vector = ", custom_source_vector)

  #print("max. number of inner iterations: ",driver.GetNumberInnerIter());
  #print("max nr of outer iterations: ",driver.GetNumberOuterIter());

  # is in domain: isdomain = driver.GetNodeDomain(iPoint)  
  #for i_vertex in range(n_vertex)
  #AllSolutions = driver.GetAllSolutions(iSOLVER)
  Body_Force_Vector = [0.0, -9.81, 0.0]
  DensityInc_0 = driver.GetDensity_FreeStreamND()
  #print("rho freestream = ",DensityInc_0)
  Force_Ref = driver.GetForce_Ref()
  #print("reference force = ",Force_Ref)

  Iter = driver.GetNumberInnerIter()
  print("1. inner iterations = ",Iter)
  # set the inner iterations to 1
  driver.SetNumberInnerIter(1)

  for inner_iter in range(Iter):
    # set the source term, per point
    print(driver.GetNumberNodes() - driver.GetNumberHaloNodes())
    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
      #SolutionVector =  driver.GetSolutionVector(iSOLVER,i_node)
      PrimitiveVector =  driver.GetPrimitiveVector(iSOLVER,i_node)
      DensityInc_i = PrimitiveVector[iDENSITY] 

      for iDim in range(nDim):
        custom_source_vector[iDim+1] = -(DensityInc_i - DensityInc_0) * Body_Force_Vector[iDim] / Force_Ref

      driver.SetPointCustomSource(iSOLVER, i_node, custom_source_vector)

    print("   ***   inner iteration:",inner_iter)
    driver.Preprocess(inner_iter)

    # Run one time iteration.
    driver.Run()
    driver.Postprocess()
    driver.Update()

    # Monitor the solver and output solution to file if required.
    #driver.Monitor(inner_iter)
    driver.Output(inner_iter)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

if __name__ == '__main__':
  main()
