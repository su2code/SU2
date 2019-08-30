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
import pysu2ad          # imports the SU2 adjoint-wrapped module
import numpy as np
import pyamg
from SU2 import amginria as su2amg
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
    SU2Driver = pysu2ad.CDiscAdjSinglezoneDriver(options.filename, options.nZone, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return

  # Retrieve some control parameters from the driver
  TimeIter = SU2Driver.GetIter()
  nTimeIter = SU2Driver.GetnIter()

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
  TimeAdjIter = 0
  nTimeAdjIter = 1

  if rank == 0:
    print("\n-------------------------- Begin Adjoint Solver -------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  while (TimeAdjIter < nTimeAdjIter):
    # Time iteration preprocessing
    SU2Driver.Preprocess(TimeAdjIter)
    # Run one time-step (static: one simulation)
    SU2Driver.Run()
    # Postprocess
    SU2Driver.Postprocess()  
    # Update the solver for the next time iteration
    SU2Driver.Update()
    # Monitor the solver and output solution to file if required
    stopCalc = SU2Driver.Monitor(TimeAdjIter)
    
    # Output the solution to file
    SU2Driver.Output(TimeAdjIter)
    if (stopCalc == True):
      break
    # Update control parameters
    TimeAdjIter += 1

  # Initialize the error estimation driver
  try:
    SU2Error = pysu2ad.CErrorEstimationDriver(SU2Driver, options.nZone, comm);
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
  Sol = np.array(SU2Error.GetAdaptationData(), float)
  nVar = Sol.shape[1]

  # Free memory
  SU2Error.CleanAdaptationData()

  # Sort the connectivity data for AMG
  SU2Error.SetConnectivityData()

  # Retrieve the connectivity data
  iZone = 0
  iInst = 0

  if options.nDim == 2:
    Edg = np.array(SU2Error.GetConnectivityEdg(iZone, iInst), int)
    Tri = np.array(SU2Error.GetConnectivityTri(iZone, iInst), int)
    Tet = np.empty(0, int)
  else:
    Edg = np.empty(0, int)
    Tri = np.array(SU2Error.GetConnectivityTri(iZone, iInst), int)
    Tet = np.array(SU2Error.GetConnectivityTet(iZone, iInst), int)

  # Free memory
  SU2Error.CleanConnectivityData()

  # Gather data to rank 0
  if options.with_MPI == True:
    sendSolCounts = np.array(comm.gather(len(Sol)*nVar, root=0))
    sendEdgCounts = np.array(comm.gather(len(Edg)*3,    root=0))
    sendTriCounts = np.array(comm.gather(len(Tri)*4,    root=0))
    sendTetCounts = np.array(comm.gather(len(Tet)*5,    root=0))

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
      recvSolBuf = np.array(recvSolBuf).reshape(recvSolBuf.size/nVar, nVar)
      recvEdgBuf = np.array(recvEdgBuf).reshape(recvEdgBuf.size/3,    3   )
      recvTriBuf = np.array(recvTriBuf).reshape(recvTriBuf.size/4,    4   )
      recvTetBuf = np.array(recvTetBuf).reshape(recvTetBuf.size/5,    5   )

      # Transfer all data structures to mesh dict for pyAMG
      print("Transferring mesh and solution data to pyAMG.")
      mesh = dict()

      mesh['dimension']  = options.nDim
      mesh['Edges']      = recvEdgBuf
      mesh['Triangles']  = recvTriBuf
      mesh['Tetrahedra'] = recvTetBuf
      mesh['Corners']    = []

      mesh['solution']   = recvSolBuf[:,options.nDim:]

      if options.nDim == 2:
        mesh['xy']       = recvSolBuf[:,0:2]
        mesh['metric']   = recvSolBuf[:,-3:]
      else:
        mesh['xyz']      = recvSolBuf[:,0:3]
        mesh['metric']   = recvSolBuf[:,-6:]

      # Get markers
      nMarker_All = SU2Error.GetnMarker_All()
      mesh['markers'] = np.empty(nMarker_All, 'string')
      for iMarker in range(nMarker_All):
        mesh['markers'][iMarker] = SU2Error.GetMarker_All_TagBound(iMarker)

      # Remesh options
      remesh_options                = {}
      remesh_options['Lp']          = 1
      remesh_options['gradation']   = 1.8
      remesh_options['hmax']        = 10.0
      remesh_options['hmin']        = 1.0e-8
      remesh_options['target']      = 2500
      remesh_options['logfile']     = "amg.log"

      # Run pyAMG
      print("Running pyAMG.")
      try:
        mesh_new = pyamg.adapt_mesh(mesh, remesh_options)        
      except:
          sys.stderr("## ERROR : pyamg failed.\n")
          raise

      current_mesh = "mesh_new.meshb"
      current_solution = "mesh_new.solb" 

      su2amg.write_mesh(current_mesh, current_solution, mesh_new)

      if options.nDim == 2:
        SolAdap = mesh_new['solution'].tolist()
        PoiAdap = mesh_new['xy'].tolist()
        EdgAdap = mesh_new['Edges'].tolist()
        TriAdap = mesh_new['Triangles'].tolist()
        TetAdap = [[]]
      else:
        SolAdap = mesh_new['solution'].tolist()
        EdgAdap = [[]]
        PoiAdap = mesh_new['xyz'].tolist()
        TriAdap = mesh_new['Triangles'].tolist()
        TetAdap = mesh_new['Tetrahedra'].tolist()

      del [mesh, mesh_new]

      SU2Driver.Adapted_Input_Preprocessing(comm, options.filename,
                                            SolAdap, PoiAdap, EdgAdap, TriAdap, TetAdap,
                                            options.nDim, iZone, options.nZone)

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
