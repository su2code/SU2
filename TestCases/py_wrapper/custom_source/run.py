#!/usr/bin/env python

## \file run.py
#  \brief Unsteady inlet boundary conditions.
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
import math
import itertools
# from mpi4py import MPI

def main():
  """
  Run the flow solver with a custom inlet (function of time and space).
  """
  # comm = MPI.COMM_WORLD
  comm = 0

  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2.CSinglezoneDriver('lam_buoyancy_cavity.cfg', 1, comm)
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

  # we need to add a source term to the energy equation. For this, we need to get the solver, and the variable first.
  # we then loop over all points and for these points, we add the source term
  nDim = driver.GetNumberDimensions()
  nNodes = driver.GetNumberNodes()
  #solverindex = driver.GetSolverIndices()
  #primindex = driver.GetPrimitiveIndices()

  # index to the flow solver
  iSOLVER = driver.GetSolverIndices()['INC.FLOW']
  print("index of flow solver = ",iSOLVER)
  coords = driver.Coordinates()
  # all the indices and the map to the names of the primitives
  primindex = driver.GetPrimitiveIndices()
  print("indices of primitives=",primindex)
  print("number of primitives:",len(primindex))

  nElem = driver.GetNumberElements()
  print("number of elements:",nElem)

  nVars = driver.GetNumberSolverVars(iSOLVER)
  print("number of solver variables:",nVars)
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

  inner_iter = driver.GetNumberInnerIter();
  print("max. number of inner iterations: ",driver.GetNumberInnerIter());
  print("max nr of outer iterations: ",driver.GetNumberOuterIter());

  # is in domain: isdomain = driver.GetNodeDomain(iPoint)  
  #for i_vertex in range(n_vertex)
  #AllSolutions = driver.GetAllSolutions(iSOLVER)
  Body_Force_Vector = [0.0, -9.81, 0.0]
  DensityInc_0 = driver.GetDensity_FreeStreamND()
  print("rho freestream = ",DensityInc_0)
  Force_Ref = driver.GetForce_Ref()
  print("reference force = ",Force_Ref)
  # get the density from the solution
  #residual[iDim+1] = -(DensityInc_i - DensityInc_0) * Body_Force_Vector[iDim] / Force_Ref;

  #for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
    #print("(x,y) = ( ",coords.Get(i_node,0)," , ",coords.Get(i_node,1)," )")
    #custom_source = -9.81
    # we need to set the custom source to the correct primitive equation, in this case let us add it to the y-momentum equation
    #i_momy = 2
    #driver.SetCustomSource(iSOLVER, iVar, i_node,custom_source)
    #DensityInc_i =  driver.GetSolution(iSOLVER,iPoint,2)
    #driver.SetPointCustomSource(iSOLVER, i_node,custom_source_vector)



  # we get the numer of iterations:

  DensityInc_i = 1
  # we set the actual inner iterations to 1
  #for time_iter in range(inner_iter):
  for inner_iter in range(500):

    # set the source term, per point
    print("loop over nodes and set the source term field")
    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
      #print("node = ",i_node)
      SolutionVector =  driver.GetSolutionVector(iSOLVER,i_node)
      #print("solutionvector=",SolutionVector)
      PrimitiveVector =  driver.GetPrimitiveVector(iSOLVER,i_node)
      #print("primitivevector=",PrimitiveVector)
      DensityInc_i = PrimitiveVector[iDENSITY] 
      #residual[iDim+1] = -Volume * (DensityInc_i - DensityInc_0) * Body_Force_Vector[iDim] / Force_Ref;
      for iDim in range(nDim):
        custom_source_vector[iDim+1] = -(DensityInc_i - DensityInc_0) * Body_Force_Vector[iDim] / Force_Ref
      #print("density=",DensityInc_i)
      driver.SetPointCustomSource(iSOLVER, i_node,custom_source_vector)
    print("end setting custom source term")

    # Change the total pressure.
    #if marker_id >= 0:
    #  for i_vertex in range(driver.GetNumberMarkerNodes(marker_id)):
    #    y = driver.MarkerCoordinates(marker_id)(i_vertex, 1)
    #    t = time_iter * dt
    #    pt = 1e5 + 2e4 * y / 0.01 * (1 - math.cos(2 * math.pi * t / 0.1))
    #    driver.SetMarkerCustomInletFlowVar1(marker_id, i_vertex, pt)
    #driver.BoundaryConditionsUpdate()
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
