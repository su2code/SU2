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

# with mpi:
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
# without mpi:
#  comm = 0

# flame temperature of the methane-air mixture (phi=0.5, P=5)
Tf = 1777

# unburnt temperature of the methane-air mixture (phi=0.5, P=5)
Tu = 673.0
Pu = 5.0
phi = 0.5
# unburnt density at P=5
rho_u = 2.52
# unburnt thermal conductivity of methane-air (phi=0.5, P=5)
k_u = 0.0523
# unburnt heat capacity of methane-air (phi=0.5, P=5)
cp_u = 1311.0

# P = rho*R*T
# 5 = 2.55 * R * 673
# R = 0.0029


# ################################################################## #
# create a function for the initial progress variable c              #
# ################################################################## #
def initC(coord):
    x = coord[0]
    #y = coord[1]
    #z = coord[2]
    # location where the flame should be
    flame_x = 0.012
    if (x < flame_x):
      C = 0.0
    else:
      C = 1.0

    return C

# ################################################################## #
# loop over all vertices and set the species progress variable c     #
# ################################################################## #
def SetInitialSpecies(SU2Driver):
    allCoords = SU2Driver.Coordinates()
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    for iPoint in range(SU2Driver.GetNumberNodes() - SU2Driver.GetNumberHaloNodes()):
      coord = allCoords.Get(iPoint)
      C = initC(coord)
      # now update the initial condition for the species
      SU2Driver.Solution(iSPECIESSOLVER).Set(iPoint,0,C)

# ################################################################## #
# Temperature is an algebraic function of c
# ################################################################## #
def update_temperature(SU2Driver, iPoint):
    # first, get the progress variable
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    # Note: returns a list
    C = SU2Driver.Solution(iSPECIESSOLVER)(iPoint,0)

    T = Tu*(1-C) + Tf*C

    prim_indices = SU2Driver.GetPrimitiveIndices()
    iTemp = prim_indices['TEMPERATURE']
    SU2Driver.Primitives().Set(iPoint,iTemp, T)

# ################################################################## #
# Source term according to Zimont
# ################################################################## #
def zimont(SU2Driver, iPoint):

    iSSTSOLVER = SU2Driver.GetSolverIndices()['SST']
    tke, dissipation = SU2Driver.Solution(iSSTSOLVER)(iPoint)

    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    # get the gradient of species_0
    gradc = SU2Driver.Gradient(iSPECIESSOLVER)(iPoint,0)
    primindex = SU2Driver.GetPrimitiveIndices()
    iDENSITY = primindex.get("DENSITY")
    iMU = primindex.get("LAMINAR_VISCOSITY")

    # laminar burning velocity of methane-air at phi=0.5, P=5
    Slu = 0.232

    rho = SU2Driver.Primitives()(iPoint,iDENSITY)
    mu = SU2Driver.Primitives()(iPoint,iMU)
    nu=mu/rho
    # Turbulent Flamespeed Closure with Dinkelacker correction
    up = np.sqrt((2.0/3.0) * tke )
    lt = (0.09**0.75) * (tke**1.5) / dissipation
    Re = up*lt/nu
    Le = 1.0
    Ut = Slu * (1.0 + (0.46/Le) * np.power(Re,0.25) * np.power(up/Slu,0.3) * np.power(Pu,0.2) )
    norm_gradc = np.sqrt(gradc[0]*gradc[0] + gradc[1]*gradc[1])
    Sc = rho_u * Ut * norm_gradc

    return Sc

# ################################################################## #
# Get the list of solver variable names
# ################################################################## #
def getsolvar(SU2Driver):
    primindex = SU2Driver.GetPrimitiveIndices()
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    nVars = SU2Driver.Solution(iFLOWSOLVER).Shape()[1]
    varindex = primindex.copy()
    for prim in varindex.copy():
      if varindex[prim] >=nVars:
        del varindex[prim]
    varindex = dict(sorted(varindex.items(), key=lambda item: item[1]))
    return varindex

# ################################################################## #
# Main routine
# ################################################################## #
def main():

  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2.CSinglezoneDriver('psi.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CSinglezoneDriver : ', exception)
    raise

  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------")
    sys.stdout.flush()

  nDim = driver.GetNumberDimensions()

  # index to the flow solver
  # C.FLOW
  # INC.FLOW
  # HEAT
  # FLAMELET
  # SPECIES
  # SA
  # SST
  iFLOWSOLVER = driver.GetSolverIndices()['INC.FLOW']
  iSPECIESSOLVER = driver.GetSolverIndices()['SPECIES']
  iSSTSOLVER = driver.GetSolverIndices()['SST']
  # all the indices and the map to the names of the primitives
  primindex = driver.GetPrimitiveIndices()
  nElem = driver.GetNumberElements()
  nVars = driver.Solution(iFLOWSOLVER).Shape()[1]
  nVarsSpecies = driver.Solution(iSPECIESSOLVER).Shape()[1]
  nVarsTurb = driver.Solution(iSSTSOLVER).Shape()[1]

  if rank == 0:
    print("Dimensions of the problem = ",nDim)
    print("index of flow solver = ",iFLOWSOLVER)
    print("index of turbulence solver = ",iSSTSOLVER)
    print("indices of primitives=",primindex)
    print("number of primitives:",len(primindex))
    print("number of elements:",nElem)
    print("number of flow solver variables:",nVars)
    print("number of species solver variables:",nVarsSpecies)
    print("number of turbulence solver variables:",nVarsTurb)
    sys.stdout.flush()



  # ### Check if we do a restart or not. ###
  with open('psi.cfg') as f:
    if 'RESTART_SOL= YES' in f.read():
      if rank == 0:
        print("restarting from file")
    else:
        # We can set an initial condition by calling this function:
        if rank == 0:
          print("Using user defined initial condition.")
        SetInitialSpecies(driver)

  # super important to actually push the commands.
  sys.stdout.flush()

  # run N iterations
  for inner_iter in range(2):
    if (rank==0):
      print("python iteration ", inner_iter)
    driver.Preprocess(inner_iter)
    driver.Run()

    Source = driver.UserDefinedSource(iSPECIESSOLVER)

    # set the source term, per point
    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
      # add source term:
      # default TFC of Zimont: rho*Sc = rho_u * U_t * grad(c)
      S = zimont(driver,i_node)
      Source.Set(i_node,0,S)

    # for the update of temperature, we need to update also the halo nodes
    for i_node in range(driver.GetNumberNodes()):
      # set the temperature to T = c*Tf + (1-c)*Tu
      update_temperature(driver, i_node)

    driver.Postprocess()
    driver.Update()
    # Monitor the solver and output solution to file if required.
    #driver.Monitor(inner_iter)
    # Output the solution to file
    driver.Output(inner_iter)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

if __name__ == '__main__':
  main()
