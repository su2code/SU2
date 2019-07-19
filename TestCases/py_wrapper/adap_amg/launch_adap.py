#!/usr/bin/env python

## \file run_adaptation.py
#  \brief Python script to launch SU2-AMG interface
#  \author Brian Munguía
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import sys
import shutil
from optparse import OptionParser	# use a parser for configuration
import pysu2            # imports the SU2 wrapped module
import pysu2ad          # imports the SU2 adjoint-wrapped module
import numpy as np
from math import *

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

  # Command line options
  parser=OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
  parser.add_option("--nDim", dest="nDim", default=2, help="Define the number of DIMENSIONS",
                    metavar="DIMENSIONS")
  parser.add_option("--nZone", dest="nZone", default=1, help="Define the number of ZONES",
                    metavar="ZONES")
  parser.add_option("--parallel", action="store_true",
                    help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  parser.add_option("--fsi", dest="fsi", default="False", help="Launch the FSI driver", metavar="FSI")

  parser.add_option("--fem", dest="fem", default="False", help="Launch the FEM driver (General driver)", metavar="FEM")

  parser.add_option("--harmonic_balance", dest="harmonic_balance", default="False",
                    help="Launch the Harmonic Balance (HB) driver", metavar="HB")


  (options, args) = parser.parse_args()
  options.nDim  = int( options.nDim )
  options.nZone = int( options.nZone )

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
    SU2Driver = pysu2ad.CDiscAdjSinglezoneDriver(options.filename, options.nZone, options.nDim, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return

  MarkerID = None
  MarkerName = 'airfoil'       # Specified by the user

  # Get all the boundary tags
  MarkerList = SU2Driver.GetAllBoundaryMarkersTag()

  # Get all the markers defined on this rank and their associated indices.
  allMarkerIDs = SU2Driver.GetAllBoundaryMarkers()

  #Check if the specified marker exists and if it belongs to this rank.
  if MarkerName in MarkerList and MarkerName in allMarkerIDs.keys():
    MarkerID = allMarkerIDs[MarkerName]

  # Number of vertices on the specified marker (per rank)
  nVertex_Marker = 0         #total number of vertices (physical + halo)
  nVertex_Marker_HALO = 0    #number of halo vertices
  nVertex_Marker_PHYS = 0    #number of physical vertices

  if MarkerID != None:
    nVertex_Marker = SU2Driver.GetNumberVertices(MarkerID)
    nVertex_Marker_HALO = SU2Driver.GetNumberHaloVertices(MarkerID)
    nVertex_Marker_PHYS = nVertex_Marker - nVertex_Marker_HALO

  # Retrieve some control parameters from the driver
  TimeIter = SU2Driver.GetExtIter()
  nTimeIter = SU2Driver.GetnExtIter()

  # Time loop is defined in Python so that we have access to SU2 functionalities at each time step
  if rank == 0:
    print("\n--------------------------- Begin Flow Solver ---------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  while (TimeIter < nTimeIter):
    # Time iteration preprocessing
    stopCalc = SU2Driver.DirectIteration(TimeIter)
    if (stopCalc == True):
      break
    # Update control parameters
    TimeIter += 1

  # Retrieve some control parameters from the driver
  TimeIter = SU2Driver.GetExtIter()
  nTimeIter = SU2Driver.GetnExtIter()

  if rank == 0:
    print("\n-------------------------- Begin Adjoint Solver -------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  while (TimeIter < nTimeIter):
    # Time iteration preprocessing
    SU2Driver.Preprocess(TimeIter)
    # Run one time-step (static: one simulation)
    SU2Driver.Run()
    # Postprocess
    SU2Driver.Postprocess()  
    # Update the solver for the next time iteration
    SU2Driver.Update()
    # Monitor the solver and output solution to file if required
    stopCalc = SU2Driver.Monitor(TimeIter)
    
    # Output the solution to file
    SU2Driver.Output(TimeIter)
    if (stopCalc == True):
      break
    # Update control parameters
    TimeIter += 1

  # Initialize the error estimation driver
  try:
    SU2Error = pysu2ad.CErrorEstimationDriver(SU2Driver, options.nZone, options.nDim, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return

  if rank == 0:
    print("\n------------------------- Begin Error Estimation ------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  # Compute the metric
  SU2Error.ComputeMetric()

  # Sort the solution data for AMG
  SU2Error.SetAdaptationData()

  # Retrieve the solution data
  nPoint_Local = SU2Error.GetnPoinPar()
  nVar_Par     = SU2Error.GetnVarPar()

  Sol = np.zeros((nPoint_Local, nVar_Par))
  for iPoint in range(0, nPoint_Local):
    for iVar in range(0, nVar_Par):
      Sol[iPoint,iVar] = SU2Error.GetAdaptationData(iVar, iPoint)

  # Sort the connectivity data for AMG
  SU2Error.SetConnectivityData()

  # Retrieve the connectivity data
  nEdg = SU2Error.GetnEdgPar()
  nTri = SU2Error.GetnTriPar(options.nDim)
  nTet = SU2Error.GetnTetPar()

  Edg = np.zeros((nEdg, 3), int)
  Tri = np.zeros((nTri, 4), int)
  Tet = np.zeros((nTet, 5), int)

  for iEdg in range(0, nEdg):
    e1, e2 = SU2Error.GetConnectivityEdg(iEdg)
    Edg[iEdg, 0:2] = np.array([e1, e2])

  for iTri in range(0, nTri):
    t1, t2, t3 = SU2Error.GetConnectivityTri(iTri, options.nDim)
    Tri[iTri, 0:3] = np.array([t1, t2, t3])

  for iTet in range(0, nTet):
    t1, t2, t3, t4 = SU2Error.GetConnectivityTet(iTet)
    Tet[iTet, 0:4] = np.array([t1, t2, t3, t4])

  # Gather data to rank 0
  if options.with_MPI == True:
    sendSolCounts = np.array(comm.gather(len(Sol)*nVar_Par, root=0))
    sendEdgCounts = np.array(comm.gather(len(Edg)*3,        root=0))
    sendTriCounts = np.array(comm.gather(len(Tri)*4,        root=0))
    sendTetCounts = np.array(comm.gather(len(Tet)*5,        root=0))

    if rank == 0:
      recvSolBuf = np.empty(sum(sendSolCounts), Sol.dtype)
      recvEdgBuf = np.empty(sum(sendEdgCounts), Edg.dtype)
      recvTriBuf = np.empty(sum(sendTriCounts), Tri.dtype)
      recvTetBuf = np.empty(sum(sendTetCounts), Tet.dtype)

    else:
      recvSolBuf = None
      recvEdgBuf = None
      recvTriBuf = None
      recvTetBuf = None
    
    comm.Gatherv(sendbuf=Sol, recvbuf=(recvSolBuf, sendSolCounts), root=0)
    comm.Gatherv(sendbuf=Edg, recvbuf=(recvEdgBuf, sendEdgCounts), root=0)
    comm.Gatherv(sendbuf=Tri, recvbuf=(recvTriBuf, sendTriCounts), root=0)
    comm.Gatherv(sendbuf=Tet, recvbuf=(recvTetBuf, sendTetCounts), root=0)

    if rank == 0:
      recvSolBuf = np.array(recvSolBuf).reshape(recvSolBuf.size/nVar_Par, nVar_Par)
      recvEdgBuf = np.array(recvEdgBuf).reshape(recvEdgBuf.size/3,        3       )
      recvTriBuf = np.array(recvTriBuf).reshape(recvTriBuf.size/4,        4       )
      recvTetBuf = np.array(recvTetBuf).reshape(recvTetBuf.size/5,        5       )

  # Postprocess the solver and exit cleanly
  SU2Driver.Postprocessing()

  if SU2Driver != None:
    del SU2Driver

  

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()  
