#!/usr/bin/env python

## \file launch_species_init.py
#  \brief Python script to launch SU2_CFD with customized initial conditions using the Python wrapper.
#  \author Nijso Beishuizen
#  \version 7.5.1 "Blackbird"
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

# ### INSTRUCTIONS
#
# run this script using:
# $ python launch_species_init.py --parallel -f species_init.cfg
#
# ###

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import sys
from optparse import OptionParser	# use a parser for configuration
import pysu2			                # imports the SU2 wrapped module
from math import *

# create a function for initial species distribution
# we want the vertical inlet to be filled with species_0=1
# and the horizontal part with species_0=0
def initSpecies(coord):
    #x = coord[0] # not used
    y = coord[1]
    #z = coord[2] # only for 3D

    species_0 = 0.0
    line = 0.016

    if (y>line):
      species_0 = 1.0

    return species_0

# create initial velocity field
def initVelocity(coord):
    #x = coord[0] # not used
    y = coord[1]
    #z = coord[2] # only for 3D

    velx = 1.0
    vely = 0.0
    velz = 0.0

    line = 0.016

    if (y>line):
        velx = 0.0
        vely = -1.0
        velz = 0.0

    return [velx, vely, velz]


# now loop over all vertices in the domain and set the species progress variable
def SetInitialSpecies(driver):
    # ### get the mesh coordinates ###
    coords = driver.Coordinates()
    print("Dimension of the problem: N = ",driver.GetNumberDimensions())
    print("nodes = ",driver.GetNumberNodes())
    print("elements = ",driver.GetNumberElements())
    # ### get the list of solver indices ###
    solverIndices = driver.GetSolverIndices()
    print("solver indices = ",solverIndices)

    speciesIndex = solverIndices["SPECIES"]
    print("species index = ",speciesIndex)
    # ### get the species solution ###
    solution = driver.Solution(speciesIndex)

    flowIndex = solverIndices["INC.FLOW"]
    print("flow index = ",flowIndex)
    # ### get the flow solution ###
    flowsolution = driver.Solution(flowIndex)
    primitiveIndices = driver.GetPrimitiveIndices()
    print("flow solver indices = ",primitiveIndices)
    velxIndex = primitiveIndices["VELOCITY_X"]
    print("index of X-velocity = ",velxIndex)
    velyIndex = primitiveIndices["VELOCITY_Y"]
    print("index of Y-velocity = ",velyIndex)

    print("loop over vertices")
    # number of points on mesh 0
    nPoints = driver.GetNumberNodes()
    print("number of points = ",nPoints)
    # now loop from 0..nPoints and get the coordinates on the (point,imesh)
    # loop over all points and get the progress variable
    for iPoint in range(nPoints):
      species_0 = initSpecies(coords.Get(iPoint))
      solution.Set(iPoint, (species_0, 1.0-species_0))
      vel = initVelocity(coords.Get(iPoint))

      flowsolution.Set(iPoint, velxIndex, vel[0] )
      flowsolution.Set(iPoint, velyIndex, vel[1] )

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

  #  ### Command line options ###
  parser=OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
  parser.add_option("--parallel", action="store_true",
                    help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  (options, args) = parser.parse_args()

  # number of spatial dimensions of the problem
  options.nDim = int(2)
  # number of (multi-)zones
  options.nZone = int(1)

  # Import mpi4py for parallel run
  if options.with_MPI == True:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
  else:
    comm = 0
    rank = 0

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  try:
      SU2Driver = pysu2.CSinglezoneDriver(options.filename, options.nZone, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return

  print("Start calling SetInitialSpecies")
  SetInitialSpecies(SU2Driver)
  print("End calling SetInitialSpecies")

  # Time loop is defined in Python so that we have acces to SU2 functionalities at each time step
  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  # number of inner terations per Run (config option ITER= )
  nIter = 10

  # starting iteration
  Iter = 0

  # total number of outer iterations: N=1 for steady cases
  maxIter = 1

  while (Iter < maxIter):
    print("Iter = ",Iter)

    # this routine sets the number of inner iterations ITER= for steady cases
    SU2Driver.Preprocess(nIter)

    # This runs for the number of iterations nIter
    SU2Driver.Run()

    SU2Driver.Postprocess()
    SU2Driver.Update()

    # stopping criterion, based on convergence and nr of iterations
    # we stop after 1 outer iteration
    stopCalc = SU2Driver.Monitor(1)
    if (stopCalc == True):
      break

    SU2Driver.Output(1)

  # ### We are finished, let's exit ###
  SU2Driver.Output(1)
  if SU2Driver != None:
    del SU2Driver

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
