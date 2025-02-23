#!/usr/bin/env python

## \file run.py
#  \brief turbulent premixed dump combustor simulation (PSI flame)
# phi=0.5, methane-air, U=40 m/s
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
import numpy as np
#from mpi4py import MPI

# unburnt temperature of the propane-air mixture
# flame temperature of the methane-air mixture (phi=0.5, P=5)
Tf = 1777

Tu = 673.0
Pu = 5.0
phi = 0.5
# unburnt density at P=5
rho_u = 2.52
# unburnt thermal conductivity of methane-air at phi=0.5 (P=5)
k_u = 0.0523
# unburnt heat capacity of methane-air at phi=0.5 (P=1)
cp_u = 1311.0

# P = rho*R*T
# 5 = 2.55 * R * 673
# R = 0.0029

# ################################################################## #
# create a function for the initial progress variable                #
# ################################################################## #
def initC(coord):
    x = coord[0]
    #y = coord[1]
    #z = coord[2]
    #print("x,y = ",x," ",y)
    # location where the flame should be
    flame_x = 0.012
    if (x < flame_x):
      C = 0.0
    else:
      C = 1.0

    return C

# ################################################################## #
# loop over all vertices and set the species progress variable       #
# ################################################################## #
def SetInitialSpecies(SU2Driver):
    allCoords = SU2Driver.Coordinates()
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    print("index of species solver = ",iSPECIESSOLVER)
    nVarsSpecies = SU2Driver.GetNumberSolverVars(iSPECIESSOLVER)
    print("number of species solver variables:",nVarsSpecies)
    for iPoint in range(SU2Driver.GetNumberNodes() - SU2Driver.GetNumberHaloNodes()):
      coord = allCoords.Get(iPoint)
      C = initC(coord)
      # now update the initial condition
      SU2Driver.SetSolutionVector(iSPECIESSOLVER, iPoint, [C])

def SetInitialVelocity(SU2Driver):
    allCoords = SU2Driver.Coordinates()
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    print("index of FLOW solver = ",iFLOWSOLVER)
    nVarsFlow = SU2Driver.GetNumberSolverVars(iFLOWSOLVER)
    print("number of flow solver variables:",nVarsFlow)
    for iPoint in range(SU2Driver.GetNumberNodes() - SU2Driver.GetNumberHaloNodes()):
      coord = allCoords.Get(iPoint)
      C = initC(coord)
      # now update the initial condition
      SU2Driver.SetSolutionVector(iFLOWSOLVER, iPoint, [C])

def update_temperature(SU2Driver, iPoint):
    # first, get the progress variable
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    # returns a list
    C = SU2Driver.GetSolutionVector(iSPECIESSOLVER, iPoint)
    T = Tu*(1-C[0]) + Tf*C[0]
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    solvar = list(SU2Driver.GetSolutionVector(iFLOWSOLVER, iPoint))
    # the list with names
    solindex = getsolvar(SU2Driver)
    primindex = SU2Driver.GetPrimitiveIndices()
    iTEMP = solindex.get("TEMPERATURE")
    solvar[iTEMP] = T
    SU2Driver.SetSolutionVector(iFLOWSOLVER, iPoint, solvar)


def zimont(SU2Driver, iPoint):
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    iSSTSOLVER = SU2Driver.GetSolverIndices()['SST']
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    primindex = SU2Driver.GetPrimitiveIndices()
    primvar = list(SU2Driver.GetPrimitiveVector(iFLOWSOLVER, iPoint))

    iDENSITY = primindex.get("DENSITY")
    iMU = primindex.get("LAMINAR_VISCOSITY")

    # laminar burning velocity of methane-air at phi=0.5, P=5
    Slu = 0.232

    rho = primvar[iDENSITY]
    mu = primvar[iMU]
    nu=mu/rho
    tke, dissipation = SU2Driver.GetSolutionVector(iSSTSOLVER,iPoint)
    gradc = SU2Driver.GetGradient(iSPECIESSOLVER,iPoint,0)
    # Turbulent Flamespeed Closure with Dinkelacker correction
    up = np.sqrt((2.0/3.0) * tke )
    lt = (0.09**0.75) * (tke**1.5) / dissipation
    Re = up*lt/nu
    Le = 1.0
    Ut = Slu * (1.0 + (0.46/Le)*np.power(Re,0.25)*np.power(up/Slu,0.3)*np.power(Pu,0.2))

    norm_gradc = np.sqrt(gradc[0]*gradc[0] + gradc[1]*gradc[1])

    Sc = rho_u * Ut * norm_gradc

    return Sc

def getsolvar(SU2Driver):
    primindex = SU2Driver.GetPrimitiveIndices()
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    nVars = SU2Driver.GetNumberSolverVars(iFLOWSOLVER)
    varindex = primindex.copy()
    for prim in varindex.copy():
      if varindex[prim] >=nVars:
        del varindex[prim]
    varindex = dict(sorted(varindex.items(), key=lambda item: item[1]))
    return varindex


def main():
  """
  Run the flow solver with a custom inlet (function of time and space).
  """
  # comm = MPI.COMM_WORLD
  comm = 0

  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2.CSinglezoneDriver('psi.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CSinglezoneDriver : ', exception)
    raise

  print("\n------------------------------ Begin Solver -----------------------------")
  sys.stdout.flush()

  nDim = driver.GetNumberDimensions()
  print("Dimensions of the problem = ",nDim)

  # index to the flow solver
  # C.FLOW
  # INC.FLOW
  # HEAT
  # FLAMELET
  # SPECIES
  # SA
  # SST
  iFLOWSOLVER = driver.GetSolverIndices()['INC.FLOW']
  print("index of flow solver = ",iFLOWSOLVER)
  iSPECIESSOLVER = driver.GetSolverIndices()['SPECIES']
  print("index of species solver = ",iSPECIESSOLVER)
  iSSTSOLVER = driver.GetSolverIndices()['SST']
  print("index of turbulence solver = ",iSSTSOLVER)


  # all the indices and the map to the names of the primitives
  primindex = driver.GetPrimitiveIndices()
  print("indices of primitives=",primindex)
  print("number of primitives:",len(primindex))

  nElem = driver.GetNumberElements()
  print("number of elements:",nElem)

  nVars = driver.GetNumberSolverVars(iFLOWSOLVER)
  print("number of flow solver variables:",nVars)

  nVarsSpecies = driver.GetNumberSolverVars(iSPECIESSOLVER)
  print("number of species solver variables:",nVarsSpecies)
  nVarsTurb = driver.GetNumberSolverVars(iSSTSOLVER)
  print("number of turbulence solver variables:",nVarsTurb)

  varindex = primindex.copy()
  for prim in varindex.copy():
    if varindex[prim] >=nVars:
      del varindex[prim]
  varindex = dict(sorted(varindex.items(), key=lambda item: item[1]))


# it is possible to get the solver type by doing
#  iFLOWSOLVER = driver.GetSolverIndices()['INC.FLOW']
# with the solver type we then get the solver variables using:

#  nVars = driver.GetNumberSolverVars(iFLOWSOLVER)
# and the solver variable names:
#  print("solver variable names:",varindex)

# we can overwrite the solution using:
# driver.SetSolutionVector(iSolver, iPoint, solutionVector)
#

  #print("solver variable names:",varindex)
  iDENSITY = primindex.get("DENSITY")
  #print("index of density = ",iDENSITY)

  index_Vel = varindex.get("VELOCITY_X")
  #print("index of velocity = ",index_Vel)

  custom_source_vector = [0.0 for i in range(nVars)]
  #print("custom source vector = ", custom_source_vector)

  #print("max. number of inner iterations: ",driver.GetNumberInnerIter());
  #print("max nr of outer iterations: ",driver.GetNumberOuterIter());

  # We can set an initial condition by calling this function:
  #print("Start calling SetInitialSpecies")
  #SetInitialSpecies(driver)
  #print("End calling SetInitialSpecies")

  # super important to actually push the commands.
  sys.stdout.flush()

  # run 5 iterations
  for inner_iter in range(5):

    driver.Preprocess(inner_iter)
    driver.Run()

    # set the source term, per point
    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
      # add source term:
      # default TFC of Zimont: rho*Sc = rho_u * U_t * grad(c)
      S = [zimont(driver,i_node)]
      driver.SetPointCustomSource(iSPECIESSOLVER, i_node,S)
      # at this point we also need to update the temperature based on the progress variable:
      # set the temperature to T = c*Tf + (1-c)*Tu
      update_temperature(driver, i_node)

    driver.Postprocess()
    driver.Update()
    driver.Monitor(inner_iter)
    driver.Output(inner_iter)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

if __name__ == '__main__':
  main()
