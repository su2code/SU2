#!/usr/bin/env python

## \file FSIInterface.py
#  \brief FSI interface class that handles fluid/solid solvers synchronisation and communication.
#  \author David Thomas
#  \version 7.0.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

import os, sys, shutil, copy
import numpy as np
import scipy as sp
import scipy.spatial.distance as spdist
from math import *
from rtree import index
from petsc4py import PETSc

# ----------------------------------------------------------------------
#  FSI Interface Class
# ----------------------------------------------------------------------

class Interface:
    """ 
    FSI interface class that handles fluid/solid solvers synchronisation and communication
    """
   
    def __init__(self, FSI_config, FluidSolver, SolidSolver, have_MPI):
	""" 
	Class constructor. Declare some variables and do some screen outputs.
	"""
        
        if have_MPI == True:
          from mpi4py import MPI
          self.MPI = MPI
          self.comm = MPI.COMM_WORLD			#MPI World communicator
          self.have_MPI = True
          myid = self.comm.Get_rank()
        else:
          self.comm = 0
          self.have_MPI = False
          myid = 0

	self.rootProcess = 0				#the root process is chosen to be MPI rank = 0

	self.nDim = FSI_config['NDIM']			#problem dimension

        self.haveFluidSolver = False			#True if the fluid solver is initialized on the current rank
        self.haveSolidSolver = False			#True if the solid solver is initialized on the current rank
        self.haveFluidInterface = False			#True if the current rank owns at least one fluid interface node
        self.haveSolidInterface = False			#True if the current rank owns at least one solid interface node

        self.fluidSolverProcessors = list()		#list of partitions where the fluid solver is initialized
        self.solidSolverProcessors = list()		#list of partitions where the solid solver is initialized
        self.fluidInterfaceProcessors = list()          #list of partitions where there are fluid interface nodes
        self.solidInterfaceProcessors = list() 	        #list of partitions where there are solid interface nodes

        self.fluidInterfaceIdentifier = None		#object that can identify the f/s interface within the fluid solver
        self.solidInterfaceIdentifier = None		#object that can identify the f/s interface within the solid solver

        self.fluidGlobalIndexRange = {}			#contains the global FSI indexing of each fluid interface node for all partitions
        self.solidGlobalIndexRange = {}			#contains the global FSI indexing of each solid interface node for all partitions

        self.FluidHaloNodeList = {}			#contains the the indices (fluid solver indexing) of the halo nodes for each partition
        self.fluidIndexing = {}				#links between the fluid solver indexing and the FSI indexing for the interface nodes
        self.SolidHaloNodeList = {}			#contains the the indices (solid solver indexing) of the halo nodes for each partition
        self.solidIndexing = {}				#links between the solid solver indexing and the FSI indexing for the interface nodes

	self.nLocalFluidInterfaceNodes = 0		#number of nodes (halo nodes included) on the fluid interface, on each partition
        self.nLocalFluidInterfaceHaloNode = 0		#number of halo nodes on the fluid intrface, on each partition
        self.nLocalFluidInterfacePhysicalNodes = 0	#number of physical (= non halo) nodes on the fluid interface, on each partition
	self.nFluidInterfaceNodes = 0			#number of nodes on the fluid interface, sum over all the partitions
        self.nFluidInterfacePhysicalNodes = 0		#number of physical nodes on the fluid interface, sum over all partitions

        self.nLocalSolidInterfaceNodes = 0     		#number of physical nodes on the solid interface, on each partition
        self.nLocalSolidInterfaceHaloNode = 0		#number of halo nodes on the solid intrface, on each partition
        self.nLocalSolidInterfacePhysicalNodes = 0	#number of physical (= non halo) nodes on the solid interface, on each partition
	self.nSolidInterfaceNodes = 0			#number of nodes on the solid interface, sum over all partitions
        self.nSolidInterfacePhysicalNodes = 0		#number of physical nodes on the solid interface, sum over all partitions

        if FSI_config['MATCHING_MESH'] == 'NO' and (FSI_config['MESH_INTERP_METHOD'] == 'RBF' or FSI_config['MESH_INTERP_METHOD'] == 'TPS'):
          self.MappingMatrixA = None
          self.MappingMatrixA_T = None
          self.MappingMatrixB = None
          self.MappingMatrixB_T = None
          self.d_RBF = self.nDim+1
        else:
          self.MappingMatrix = None			#interpolation/mapping matrix for meshes interpolation/mapping
          self.MappingMatrix_T = None			#transposed interpolation/mapping matrix for meshes interpolation/mapping
          self.d_RBF = 0

        self.localFluidInterface_array_X_init = None	#initial fluid interface position on each partition (used for the meshes mapping)
        self.localFluidInterface_array_Y_init = None
        self.localFluidInterface_array_Z_init = None

        self.haloNodesPositionsInit = {}		#initial position of the halo nodes (fluid side only)

        self.solidInterface_array_DispX = None		#solid interface displacement
        self.solidInterface_array_DispY = None
        self.solidInterface_array_DispZ = None

        self.solidInterfaceResidual_array_X = None	#solid interface position residual
        self.solidInterfaceResidual_array_Y = None
        self.solidInterfaceResidual_array_Z = None

        self.solidInterfaceResidualnM1_array_X = None	#solid interface position residual at the previous BGS iteration
        self.solidInterfaceResidualnM1_array_Y = None
        self.solidInterfaceResidualnM1_array_Z = None
       
        self.fluidInterface_array_DispX = None		#fluid interface displacement
        self.fluidInterface_array_DispY = None
        self.fluidInterface_array_DispZ = None

        self.fluidLoads_array_X = None			#loads on the fluid side of the f/s interface
        self.fluidLoads_array_Y = None
        self.fluidLoads_array_Z = None

        self.solidLoads_array_X = None			#loads on the solid side of the f/s interface
        self.solidLoads_array_Y = None
        self.solidLoads_array_Z = None

	self.aitkenParam = FSI_config['AITKEN_PARAM']			#relaxation parameter for the BGS method
	self.FSIIter = 0				#current FSI iteration
        self.unsteady = False				#flag for steady or unsteady simulation (default is steady)

	# ---Some screen output ---
	self.MPIPrint('Fluid solver : SU2_CFD')
	self.MPIPrint('Solid solver : {}'.format(FSI_config['CSD_SOLVER']))

	if FSI_config['TIME_MARCHING'] == 'YES':
          self.MPIPrint('Unsteady coupled simulation with physical time step : {} s'.format(FSI_config['UNST_TIMESTEP']))
          self.unsteady = True
	else:
	  self.MPIPrint('Steady coupled simulation')

	if FSI_config['MATCHING_MESH'] == 'YES':
	  self.MPIPrint('Matching fluid-solid interface')
	else:
          if FSI_config['MESH_INTERP_METHOD'] == 'TPS':
	    self.MPIPrint('Non matching fluid-solid interface with Thin Plate Spline interpolation')
          elif FSI_config['MESH_INTERP_METHOD'] == 'RBF':
            self.MPIPrint('Non matching fluid-solid interface with Radial Basis Function interpolation')
            self.RBF_rad = FSI_config['RBF_RADIUS']
            self.MPIPrint('Radius value : {}'.format(self.RBF_rad))           
          else:
	    self.MPIPrint('Non matching fluid-solid interface with Nearest Neighboor interpolation')

	self.MPIPrint('Solid predictor : {}'.format(FSI_config['DISP_PRED']))

	self.MPIPrint('Maximum number of FSI iterations : {}'.format(FSI_config['NB_FSI_ITER']))

	self.MPIPrint('FSI tolerance : {}'.format(FSI_config['FSI_TOLERANCE']))

	if FSI_config['AITKEN_RELAX'] == 'STATIC':
	  self.MPIPrint('Static Aitken under-relaxation with constant parameter {}'.format(FSI_config['AITKEN_PARAM']))
	elif FSI_config['AITKEN_RELAX'] == 'DYNAMIC':
	  self.MPIPrint('Dynamic Aitken under-relaxation with initial parameter {}'.format(FSI_config['AITKEN_PARAM']))
	else:
	  self.MPIPrint('No Aitken under-relaxation')

        self.MPIPrint('FSI interface is set')

    def MPIPrint(self, message):
      """ 
      Print a message on screen only from the master process.
      """

      if self.have_MPI == True:
        myid = self.comm.Get_rank()
      else:
        myid = 0

      if myid == self.rootProcess:
        print(message)

    def MPIBarrier(self):
      """
      Perform a synchronization barrier in case of parallel run with MPI.
      """
      
      if self.have_MPI == True:
        self.comm.barrier()

    def connect(self, FSI_config, FluidSolver, SolidSolver):
	"""
	Connection between solvers. 
	Creates the communication support between the two solvers.
	Gets information about f/s interfaces from the two solvers.
	"""
        if self.have_MPI == True:
          myid = self.comm.Get_rank()
	  MPIsize = self.comm.Get_size()
        else:
          myid = 0
          MPIsize = 1
	
	# --- Identify the fluid and solid interfaces and store the number of nodes on both sides (and for each partition) ---
        self.fluidInterfaceIdentifier = None
        self.nLocalFluidInterfaceNodes = 0
        if FluidSolver != None:
	    print('Fluid solver is initialized on process {}'.format(myid))
            self.haveFluidSolver = True
            allMovingMarkersTags = FluidSolver.GetAllMovingMarkersTag()
            allMarkersID = FluidSolver.GetAllBoundaryMarkers()
            if not allMovingMarkersTags:
                raise Exception('No interface for FSI was defined.')
            else:
                if allMovingMarkersTags[0] in allMarkersID.keys():
                    self.fluidInterfaceIdentifier = allMarkersID[allMovingMarkersTags[0]]
            if self.fluidInterfaceIdentifier != None:
	        self.nLocalFluidInterfaceNodes = FluidSolver.GetNumberVertices(self.fluidInterfaceIdentifier)
	    if self.nLocalFluidInterfaceNodes != 0:
              self.haveFluidInterface = True
	      print('Number of interface fluid nodes (halo nodes included) on proccess {} : {}'.format(myid,self.nLocalFluidInterfaceNodes))
	else:
	    pass

	if SolidSolver != None:
	    print('Solid solver is initialized on process {}'.format(myid))
            self.haveSolidSolver = True
	    self.solidInterfaceIdentifier = SolidSolver.getFSIMarkerID()
	    self.nLocalSolidInterfaceNodes = SolidSolver.getNumberOfSolidInterfaceNodes(self.solidInterfaceIdentifier)
	    if self.nLocalSolidInterfaceNodes != 0:
              self.haveSolidInterface = True
              print('Number of interface solid nodes (halo nodes included) on proccess {} : {}'.format(myid,self.nLocalSolidInterfaceNodes))
	else:
	    pass

        # --- Exchange information about processors on which the solvers are defined and where the interface nodes are lying ---
        if self.have_MPI == True:
          if self.haveFluidSolver == True:
            sendBufFluid = np.array(int(1))
          else:
            sendBufFluid = np.array(int(0))
          if self.haveSolidSolver == True:
            sendBufSolid = np.array(int(1))
          else:
            sendBufSolid = np.array(int(0))
          if self.haveFluidInterface == True:
            sendBufFluidInterface = np.array(int(1))
          else:
            sendBufFluidInterface = np.array(int(0))
          if self.haveSolidInterface == True:
            sendBufSolidInterface = np.array(int(1))
          else:
            sendBufSolidInterface = np.array(int(0))
          rcvBufFluid = np.zeros(MPIsize, dtype = int)
	  rcvBufSolid = np.zeros(MPIsize, dtype = int)
          rcvBufFluidInterface = np.zeros(MPIsize, dtype = int)
	  rcvBufSolidInterface = np.zeros(MPIsize, dtype = int)
          self.comm.Allgather(sendBufFluid, rcvBufFluid)
          self.comm.Allgather(sendBufSolid, rcvBufSolid)
          self.comm.Allgather(sendBufFluidInterface, rcvBufFluidInterface)
          self.comm.Allgather(sendBufSolidInterface, rcvBufSolidInterface)
          for iProc in range(MPIsize):
	    if rcvBufFluid[iProc] == 1:
              self.fluidSolverProcessors.append(iProc)
            if rcvBufSolid[iProc] == 1:
	      self.solidSolverProcessors.append(iProc)
            if rcvBufFluidInterface[iProc] == 1:
              self.fluidInterfaceProcessors.append(iProc)
            if rcvBufSolidInterface[iProc] == 1:
              self.solidInterfaceProcessors.append(iProc)
          del sendBufFluid, sendBufSolid, rcvBufFluid, rcvBufSolid, sendBufFluidInterface, sendBufSolidInterface, rcvBufFluidInterface, rcvBufSolidInterface
        else:
          self.fluidSolverProcessors.append(0)
	  self.solidSolverProcessors.append(0)
          self.fluidInterfaceProcessors.append(0)
          self.solidInterfaceProcessors.append(0)

	self.MPIBarrier()
	
	# --- Calculate the total number of nodes at the fluid interface (sum over all the partitions) ---
        # Calculate the number of halo nodes on each partition
        self.nLocalFluidInterfaceHaloNode = 0
	for iVertex in range(self.nLocalFluidInterfaceNodes):
            if FluidSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == True:
              GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
	      self.FluidHaloNodeList[GlobalIndex] = iVertex
              self.nLocalFluidInterfaceHaloNode += 1
        # Calculate the number of physical (= not halo) nodes on each partition
        self.nLocalFluidInterfacePhysicalNodes = self.nLocalFluidInterfaceNodes - self.nLocalFluidInterfaceHaloNode
        if self.have_MPI == True:
          self.FluidHaloNodeList = self.comm.allgather(self.FluidHaloNodeList)
        else:
          self.FluidHaloNodeList = [{}]

        # Same thing for the solid part
        self.nLocalSolidInterfaceHaloNode = 0
	#for iVertex in range(self.nLocalSolidInterfaceNodes):
            #if SoliddSolver.IsAHaloNode(self.fluidInterfaceIdentifier, iVertex) == True:
              #GlobalIndex = SolidSolver.GetVertexGlobalIndex(self.solidInterfaceIdentifier, iVertex)
	      #self.SolidHaloNodeList[GlobalIndex] = iVertex
              #self.nLocalSolidInterfaceHaloNode += 1
        self.nLocalSolidInterfacePhysicalNodes = self.nLocalSolidInterfaceNodes - self.nLocalSolidInterfaceHaloNode
        if self.have_MPI == True:
          self.SolidHaloNodeList = self.comm.allgather(self.SolidHaloNodeList)
        else:
          self.SolidHaloNodeList = [{}]


        # --- Calculate the total number of nodes (with and without halo) at the fluid interface (sum over all the partitions) and broadcast the number accross all processors ---
        sendBuffHalo = np.array(int(self.nLocalFluidInterfaceNodes))
        sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
	rcvBuffHalo = np.zeros(1, dtype=int)
        rcvBuffPhysical = np.zeros(1, dtype=int)
        if self.have_MPI == True:			
          self.comm.barrier()
	  self.comm.Allreduce(sendBuffHalo,rcvBuffHalo,op=self.MPI.SUM)
          self.comm.Allreduce(sendBuffPhysical,rcvBuffPhysical,op=self.MPI.SUM)
          self.nFluidInterfaceNodes = rcvBuffHalo[0]
          self.nFluidInterfacePhysicalNodes = rcvBuffPhysical[0]
        else:
          self.nFluidInterfaceNodes = np.copy(sendBuffHalo)
          self.nFluidInterfacePhysicalNodes = np.copy(sendBuffPhysical)
        del sendBuffHalo, rcvBuffHalo, sendBuffPhysical, rcvBuffPhysical

        # Same thing for the solid part
        sendBuffHalo = np.array(int(self.nLocalSolidInterfaceNodes))
        sendBuffPhysical = np.array(int(self.nLocalSolidInterfacePhysicalNodes))
	rcvBuffHalo = np.zeros(1, dtype=int)
        rcvBuffPhysical = np.zeros(1, dtype=int)
        if self.have_MPI == True:
	  self.comm.barrier()
	  self.comm.Allreduce(sendBuffHalo,rcvBuffHalo,op=self.MPI.SUM)
          self.comm.Allreduce(sendBuffPhysical,rcvBuffPhysical,op=self.MPI.SUM)
          self.nSolidInterfaceNodes = rcvBuffHalo[0]
          self.nSolidInterfacePhysicalNodes = rcvBuffPhysical[0]
        else:
          self.nSolidInterfaceNodes = np.copy(sendBuffHalo)
          self.nSolidInterfacePhysicalNodes = np.copy(sendBuffPhysical)
        del sendBuffHalo, rcvBuffHalo, sendBuffPhysical, rcvBuffPhysical

        # --- Store the number of physical interface nodes on each processor and allgather the information ---
        self.fluidPhysicalInterfaceNodesDistribution = np.zeros(MPIsize, dtype=int)
        if self.have_MPI == True:
          sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
          self.comm.Allgather(sendBuffPhysical,self.fluidPhysicalInterfaceNodesDistribution)
          del sendBuffPhysical
        else:
          self.fluidPhysicalInterfaceNodesDistribution[0] = self.nFluidInterfacePhysicalNodes

        # Same thing for the solid part
        self.solidPhysicalInterfaceNodesDistribution = np.zeros(MPIsize, dtype=int)
        if self.have_MPI == True:
          sendBuffPhysical = np.array(int(self.nLocalSolidInterfaceNodes))
          self.comm.Allgather(sendBuffPhysical,self.solidPhysicalInterfaceNodesDistribution)
          del sendBuffPhysical
        else:
          self.solidPhysicalInterfaceNodesDistribution[0] = self.nSolidInterfacePhysicalNodes

        # --- Calculate and store the global indexing of interface physical nodes on each processor and allgather the information ---
        if self.have_MPI == True:
          if myid in self.fluidInterfaceProcessors:
            globalIndexStart = 0
            for iProc in range(myid):
	        globalIndexStart += self.fluidPhysicalInterfaceNodesDistribution[iProc]
            globalIndexStop = globalIndexStart + self.nLocalFluidInterfacePhysicalNodes-1
          else:
            globalIndexStart = 0
            globalIndexStop = 0
          self.fluidGlobalIndexRange[myid] = [globalIndexStart,globalIndexStop]
          self.fluidGlobalIndexRange = self.comm.allgather(self.fluidGlobalIndexRange)
        else:
          temp = {}
          temp[0] = [0,self.nLocalFluidInterfacePhysicalNodes-1]
          self.fluidGlobalIndexRange = list()
          self.fluidGlobalIndexRange.append(temp)
          
	# Same thing for the solid part
        if self.have_MPI == True:
          if myid in self.solidInterfaceProcessors:
            globalIndexStart = 0
            for iProc in range(myid):
              globalIndexStart += self.solidPhysicalInterfaceNodesDistribution[iProc]
            globalIndexStop = globalIndexStart + self.nLocalSolidInterfaceNodes-1
          else:
            globalIndexStart = 0
            globalIndexStop = 0
          self.solidGlobalIndexRange[myid] = [globalIndexStart,globalIndexStop]
          self.solidGlobalIndexRange = self.comm.allgather(self.solidGlobalIndexRange)
        else:
          temp = {}
          temp[0] = [0,self.nSolidInterfacePhysicalNodes-1]
          self.solidGlobalIndexRange = list()
          self.solidGlobalIndexRange.append(temp)        

	self.MPIPrint('Total number of fluid interface nodes (halo nodes included) : {}'.format(self.nFluidInterfaceNodes))
	self.MPIPrint('Total number of solid interface nodes (halo nodes included) : {}'.format(self.nSolidInterfaceNodes))
        self.MPIPrint('Total number of fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))
        self.MPIPrint('Total number of solid interface nodes : {}'.format(self.nSolidInterfacePhysicalNodes))

	self.MPIBarrier()

        # --- Create all the PETSc vectors required for parallel communication and parallel mesh mapping/interpolation (working for serial too) ---
        if self.have_MPI == True:
          self.solidInterface_array_DispX = PETSc.Vec().create(self.comm)
          self.solidInterface_array_DispY = PETSc.Vec().create(self.comm)
          self.solidInterface_array_DispZ = PETSc.Vec().create(self.comm)
          self.solidInterface_array_DispX.setType('mpi')
          self.solidInterface_array_DispY.setType('mpi')
          self.solidInterface_array_DispZ.setType('mpi')
        else:
          self.solidInterface_array_DispX = PETSc.Vec().create()
          self.solidInterface_array_DispY = PETSc.Vec().create()
          self.solidInterface_array_DispZ = PETSc.Vec().create()
          self.solidInterface_array_DispX.setType('seq')
          self.solidInterface_array_DispY.setType('seq')
          self.solidInterface_array_DispZ.setType('seq')
        self.solidInterface_array_DispX.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidInterface_array_DispY.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidInterface_array_DispZ.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidInterface_array_DispX.set(0.0)
	self.solidInterface_array_DispY.set(0.0)
	self.solidInterface_array_DispZ.set(0.0)

        if self.have_MPI == True:
          self.fluidInterface_array_DispX = PETSc.Vec().create(self.comm)
          self.fluidInterface_array_DispY = PETSc.Vec().create(self.comm)
          self.fluidInterface_array_DispZ = PETSc.Vec().create(self.comm)
          self.fluidInterface_array_DispX.setType('mpi')
          self.fluidInterface_array_DispY.setType('mpi')
          self.fluidInterface_array_DispZ.setType('mpi')
        else:
          self.fluidInterface_array_DispX = PETSc.Vec().create()
          self.fluidInterface_array_DispY = PETSc.Vec().create()
          self.fluidInterface_array_DispZ = PETSc.Vec().create()
          self.fluidInterface_array_DispX.setType('seq')
          self.fluidInterface_array_DispY.setType('seq')
          self.fluidInterface_array_DispZ.setType('seq')
        self.fluidInterface_array_DispX.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidInterface_array_DispY.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidInterface_array_DispZ.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidInterface_array_DispX.set(0.0)
        self.fluidInterface_array_DispY.set(0.0)
        self.fluidInterface_array_DispZ.set(0.0)

        if self.have_MPI == True:
          self.fluidLoads_array_X = PETSc.Vec().create(self.comm)
          self.fluidLoads_array_Y = PETSc.Vec().create(self.comm)
          self.fluidLoads_array_Z = PETSc.Vec().create(self.comm)
          self.fluidLoads_array_X.setType('mpi')
          self.fluidLoads_array_Y.setType('mpi')
          self.fluidLoads_array_Z.setType('mpi')
        else:
          self.fluidLoads_array_X = PETSc.Vec().create()
          self.fluidLoads_array_Y = PETSc.Vec().create()
          self.fluidLoads_array_Z = PETSc.Vec().create()
          self.fluidLoads_array_X.setType('seq')
          self.fluidLoads_array_Y.setType('seq')
          self.fluidLoads_array_Z.setType('seq')
        self.fluidLoads_array_X.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidLoads_array_Y.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidLoads_array_Z.setSizes(self.nFluidInterfacePhysicalNodes)

        if self.have_MPI == True:
          self.solidLoads_array_X = PETSc.Vec().create(self.comm)
          self.solidLoads_array_Y = PETSc.Vec().create(self.comm)
          self.solidLoads_array_Z = PETSc.Vec().create(self.comm)
          self.solidLoads_array_X.setType('mpi')
          self.solidLoads_array_Y.setType('mpi')
          self.solidLoads_array_Z.setType('mpi')
        else:
          self.solidLoads_array_X = PETSc.Vec().create()
          self.solidLoads_array_Y = PETSc.Vec().create()
          self.solidLoads_array_Z = PETSc.Vec().create()
          self.solidLoads_array_X.setType('seq')
          self.solidLoads_array_Y.setType('seq')
          self.solidLoads_array_Z.setType('seq')
        self.solidLoads_array_X.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidLoads_array_Y.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidLoads_array_Z.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidLoads_array_X.set(0.0)
        self.solidLoads_array_Y.set(0.0)
        self.solidLoads_array_Z.set(0.0)

        # --- Create the PETSc vectors required for parallel relaxed BGS algo (working for serial too) ---
        if self.have_MPI == True:
          self.solidInterfaceResidual_array_X = PETSc.Vec().create(self.comm)
          self.solidInterfaceResidual_array_Y = PETSc.Vec().create(self.comm)
          self.solidInterfaceResidual_array_Z = PETSc.Vec().create(self.comm)
          self.solidInterfaceResidual_array_X.setType('mpi')
          self.solidInterfaceResidual_array_Y.setType('mpi')
          self.solidInterfaceResidual_array_Z.setType('mpi')
        else:
          self.solidInterfaceResidual_array_X = PETSc.Vec().create()
          self.solidInterfaceResidual_array_Y = PETSc.Vec().create()
          self.solidInterfaceResidual_array_Z = PETSc.Vec().create()
          self.solidInterfaceResidual_array_X.setType('seq')
          self.solidInterfaceResidual_array_Y.setType('seq')
          self.solidInterfaceResidual_array_Z.setType('seq')
        self.solidInterfaceResidual_array_X.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidInterfaceResidual_array_Y.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidInterfaceResidual_array_Z.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)

        if self.have_MPI == True:
          self.solidInterfaceResidualnM1_array_X = PETSc.Vec().create(self.comm)
          self.solidInterfaceResidualnM1_array_Y = PETSc.Vec().create(self.comm)
          self.solidInterfaceResidualnM1_array_Z = PETSc.Vec().create(self.comm)
          self.solidInterfaceResidualnM1_array_X.setType('mpi')
          self.solidInterfaceResidualnM1_array_Y.setType('mpi')
          self.solidInterfaceResidualnM1_array_Z.setType('mpi')
        else:
          self.solidInterfaceResidualnM1_array_X = PETSc.Vec().create()
          self.solidInterfaceResidualnM1_array_Y = PETSc.Vec().create()
          self.solidInterfaceResidualnM1_array_Z = PETSc.Vec().create()
          self.solidInterfaceResidualnM1_array_X.setType('seq')
          self.solidInterfaceResidualnM1_array_Y.setType('seq')
          self.solidInterfaceResidualnM1_array_Z.setType('seq')
        self.solidInterfaceResidualnM1_array_X.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidInterfaceResidualnM1_array_Y.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidInterfaceResidualnM1_array_Z.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        self.solidInterfaceResidualnM1_array_X.set(0.0)
        self.solidInterfaceResidualnM1_array_Y.set(0.0)
        self.solidInterfaceResidualnM1_array_Z.set(0.0)

    def interfaceMapping(self,FluidSolver, SolidSolver, FSI_config):
	""" 
	Creates the one-to-one mapping between interfaces in case of matching meshes.
	Creates the interpolation rules between interfaces in case of non-matching meshes.
	"""
	if self.have_MPI == True:
          myid = self.comm.Get_rank()
	  MPIsize = self.comm.Get_size()
        else:
          myid = 0
          MPIsize = 1

	# --- Get the fluid interface from fluid solver on each partition ---
	GlobalIndex = int()
        localIndex = 0
        fluidIndexing_temp = {}
        self.localFluidInterface_array_X_init = np.zeros((self.nLocalFluidInterfacePhysicalNodes))
        self.localFluidInterface_array_Y_init = np.zeros((self.nLocalFluidInterfacePhysicalNodes))
        self.localFluidInterface_array_Z_init = np.zeros((self.nLocalFluidInterfacePhysicalNodes))
        for iVertex in range(self.nLocalFluidInterfaceNodes):
	    GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
	    posx = FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex)
	    posy = FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex)
	    posz = FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex)
	    if GlobalIndex in self.FluidHaloNodeList[myid].keys():
              self.haloNodesPositionsInit[GlobalIndex] = (posx, posy, posz)
            else:
              fluidIndexing_temp[GlobalIndex] = self.__getGlobalIndex('fluid', myid, localIndex)
              self.localFluidInterface_array_X_init[localIndex] = posx
              self.localFluidInterface_array_Y_init[localIndex] = posy
              self.localFluidInterface_array_Z_init[localIndex] = posz
              localIndex += 1
        if self.have_MPI == True:
          fluidIndexing_temp = self.comm.allgather(fluidIndexing_temp)
          for ii in range(len(fluidIndexing_temp)):
            for key, value in fluidIndexing_temp[ii].items():
              self.fluidIndexing[key] = value
        else:
          self.fluidIndexing = fluidIndexing_temp.copy()
        del fluidIndexing_temp

	# --- Get the solid interface from solid solver on each partition ---
        localIndex = 0
        solidIndexing_temp = {}
	self.localSolidInterface_array_X = np.zeros(self.nLocalSolidInterfaceNodes)
        self.localSolidInterface_array_Y = np.zeros(self.nLocalSolidInterfaceNodes)
        self.localSolidInterface_array_Z = np.zeros(self.nLocalSolidInterfaceNodes)
        for iVertex in range(self.nLocalSolidInterfaceNodes):
          GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
	  posx = SolidSolver.getInterfaceNodePosX(self.solidInterfaceIdentifier, iVertex)
	  posy = SolidSolver.getInterfaceNodePosY(self.solidInterfaceIdentifier, iVertex)
	  posz = SolidSolver.getInterfaceNodePosZ(self.solidInterfaceIdentifier, iVertex)
          if GlobalIndex in self.SolidHaloNodeList[myid].keys():
            pass
          else:
            solidIndexing_temp[GlobalIndex] = self.__getGlobalIndex('solid', myid, localIndex)
            self.localSolidInterface_array_X[localIndex] = posx
            self.localSolidInterface_array_Y[localIndex] = posy
            self.localSolidInterface_array_Z[localIndex] = posz
            localIndex += 1
        if self.have_MPI == True:
          solidIndexing_temp = self.comm.allgather(solidIndexing_temp)
          for ii in range(len(solidIndexing_temp)):
            for key, value in solidIndexing_temp[ii].items():
              self.solidIndexing[key] = value
        else:
          self.solidIndexing = solidIndexing_temp.copy()
        del solidIndexing_temp


	# --- Create the PETSc parallel interpolation matrix ---
        if FSI_config['MATCHING_MESH'] == 'NO' and (FSI_config['MESH_INTERP_METHOD'] == 'RBF' or FSI_config['MESH_INTERP_METHOD'] == 'TPS'):
          if self.have_MPI == True:
            self.MappingMatrixA = PETSc.Mat().create(self.comm)
            self.MappingMatrixB = PETSc.Mat().create(self.comm)
            self.MappingMatrixA_T = PETSc.Mat().create(self.comm)
            self.MappingMatrixB_T = PETSc.Mat().create(self.comm)
            if FSI_config['MESH_INTERP_METHOD'] == 'RBF' : 
              self.MappingMatrixA.setType('mpiaij')
              self.MappingMatrixB.setType('mpiaij')
              self.MappingMatrixA_T.setType('mpiaij')
              self.MappingMatrixB_T.setType('mpiaij')
            else:
              self.MappingMatrixA.setType('mpiaij')
              self.MappingMatrixB.setType('mpiaij')
              self.MappingMatrixA_T.setType('mpiaij')
              self.MappingMatrixB_T.setType('mpiaij')
          else:
            self.MappingMatrixA = PETSc.Mat().create()
            self.MappingMatrixB = PETSc.Mat().create()
            self.MappingMatrixA_T = PETSc.Mat().create()
            self.MappingMatrixB_T = PETSc.Mat().create()
            if FSI_config['MESH_INTERP_METHOD'] == 'RBF' : 
              self.MappingMatrixA.setType('aij')
              self.MappingMatrixB.setType('aij')
              self.MappingMatrixA_T.setType('aij')
              self.MappingMatrixB_T.setType('aij')
            else:
              self.MappingMatrixA.setType('aij')
              self.MappingMatrixB.setType('aij')
              self.MappingMatrixA_T.setType('aij')
              self.MappingMatrixB_T.setType('aij')
	  self.MappingMatrixA.setSizes((self.nSolidInterfacePhysicalNodes+self.d_RBF, self.nSolidInterfacePhysicalNodes+self.d_RBF))
          self.MappingMatrixA.setUp()
          self.MappingMatrixA.setOption(PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False)
	  self.MappingMatrixB.setSizes((self.nFluidInterfacePhysicalNodes, self.nSolidInterfacePhysicalNodes+self.d_RBF))
          self.MappingMatrixB.setUp()
          self.MappingMatrixB.setOption(PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False)
	  self.MappingMatrixA_T.setSizes((self.nSolidInterfacePhysicalNodes+self.d_RBF, self.nSolidInterfacePhysicalNodes+self.d_RBF))
          self.MappingMatrixA_T.setUp()
          self.MappingMatrixA_T.setOption(PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False)
	  self.MappingMatrixB_T.setSizes((self.nSolidInterfacePhysicalNodes+self.d_RBF, self.nFluidInterfacePhysicalNodes))
          self.MappingMatrixB_T.setUp()
          self.MappingMatrixB_T.setOption(PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False)
        else:
          if self.have_MPI == True:
            self.MappingMatrix = PETSc.Mat().create(self.comm)
            self.MappingMatrix_T = PETSc.Mat().create(self.comm)
            self.MappingMatrix.setType('mpiaij')
            self.MappingMatrix_T.setType('mpiaij')
          else:
            self.MappingMatrix = PETSc.Mat().create()
            self.MappingMatrix_T = PETSc.Mat().create()
            self.MappingMatrix.setType('aij')
            self.MappingMatrix_T.setType('aij')
	  self.MappingMatrix.setSizes((self.nFluidInterfacePhysicalNodes, self.nSolidInterfacePhysicalNodes))
          self.MappingMatrix.setUp()
          self.MappingMatrix.setOption(PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False)
	  self.MappingMatrix_T.setSizes((self.nSolidInterfacePhysicalNodes, self.nFluidInterfacePhysicalNodes))
          self.MappingMatrix_T.setUp()
          self.MappingMatrix_T.setOption(PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False)
                  
	
        # --- Fill the interpolation matrix in parallel (working in serial too) ---
        if FSI_config['MATCHING_MESH'] == 'NO' and (FSI_config['MESH_INTERP_METHOD'] == 'RBF' or FSI_config['MESH_INTERP_METHOD'] == 'TPS'):
          self.MPIPrint('Building interpolation matrices...')
          if self.have_MPI == True:
            for iProc in self.solidInterfaceProcessors:
              if myid == iProc:
                for jProc in self.solidInterfaceProcessors:  
                  self.comm.Send(self.localSolidInterface_array_X, dest=jProc, tag=1)
                  self.comm.Send(self.localSolidInterface_array_Y, dest=jProc, tag=2)
                  self.comm.Send(self.localSolidInterface_array_Z, dest=jProc, tag=3)
              if myid in self.solidInterfaceProcessors:
                sizeOfBuff = self.solidPhysicalInterfaceNodesDistribution[iProc]
                solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
                solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
                solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
                self.comm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
                self.comm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
                self.comm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
                if FSI_config['MESH_INTERP_METHOD'] == 'RBF':
                  self.RBFMeshMapping_A(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc, self.RBF_rad)
                else:
                  self.TPSMeshMapping_A(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
          else:
            if FSI_config['MESH_INTERP_METHOD'] == 'RBF':
              self.RBFMeshMapping_A(self.localSolidInterface_array_X, self.localSolidInterface_array_Y, self.localSolidInterface_array_Z, 0, self.RBF_rad)
            else:
              self.TPSMeshMapping_A(self.localSolidInterface_array_X, self.localSolidInterface_array_Y, self.localSolidInterface_array_Z, 0)
          self.MappingMatrixA.assemblyBegin()
          self.MappingMatrixA.assemblyEnd()
          self.MappingMatrixA_T.assemblyBegin()
          self.MappingMatrixA_T.assemblyEnd()
          self.MPIPrint('Matrix A is built.')
        else:
          self.MPIPrint("Building interpolation matrix...")
        self.MPIBarrier()

        if self.have_MPI == True:
          for iProc in self.solidInterfaceProcessors:
            if myid == iProc:
              for jProc in self.fluidInterfaceProcessors:
                self.comm.Send(self.localSolidInterface_array_X, dest=jProc, tag=1)
                self.comm.Send(self.localSolidInterface_array_Y, dest=jProc, tag=2)
                self.comm.Send(self.localSolidInterface_array_Z, dest=jProc, tag=3)
            if myid in self.fluidInterfaceProcessors:
              sizeOfBuff = self.solidPhysicalInterfaceNodesDistribution[iProc]
              solidInterfaceBuffRcv_X = np.zeros(sizeOfBuff)
              solidInterfaceBuffRcv_Y = np.zeros(sizeOfBuff)
              solidInterfaceBuffRcv_Z = np.zeros(sizeOfBuff)
              self.comm.Recv(solidInterfaceBuffRcv_X, iProc, tag=1)
              self.comm.Recv(solidInterfaceBuffRcv_Y, iProc, tag=2)
              self.comm.Recv(solidInterfaceBuffRcv_Z, iProc, tag=3)
              if FSI_config['MATCHING_MESH'] == 'NO':
                if FSI_config['MESH_INTERP_METHOD'] == 'RBF':
                  self.RBFMeshMapping_B(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc, self.RBF_rad)
                elif FSI_config['MESH_INTERP_METHOD'] == 'TPS':
                  self.TPSMeshMapping_B(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
                else:
                  self.NearestNeighboorMeshMapping(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
	      else:
                self.matchingMeshMapping(solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc)
        else:
          if FSI_config['MATCHING_MESH'] == 'NO':
            if FSI_config['MESH_INTERP_METHOD'] == 'RBF':
              self.RBFMeshMapping_B(self.localSolidInterface_array_X, self.localSolidInterface_array_Y, self.localSolidInterface_array_Z, 0, self.RBF_rad)
            elif FSI_config['MESH_INTERP_METHOD'] == 'TPS' :
              self.TPSMeshMapping_B(self.localSolidInterface_array_X, self.localSolidInterface_array_Y, self.localSolidInterface_array_Z, 0)
            else:
              self.NearestNeighboorMeshMapping(self.localSolidInterface_array_X, self.localSolidInterface_array_Y, self.localSolidInterface_array_Z, 0) 
	  else:	    
            self.matchingMeshMapping(self.localSolidInterface_array_X, self.localSolidInterface_array_Y, self.localSolidInterface_array_Z, 0)
        
        if FSI_config['MATCHING_MESH'] == 'NO' and (FSI_config['MESH_INTERP_METHOD'] == 'RBF' or FSI_config['MESH_INTERP_METHOD'] == 'TPS'):
          self.MappingMatrixB.assemblyBegin()
          self.MappingMatrixB.assemblyEnd()
          self.MappingMatrixB_T.assemblyBegin()
          self.MappingMatrixB_T.assemblyEnd()
          self.MPIPrint('Matrix B is built.')
        else:
          self.MappingMatrix.assemblyBegin()
          self.MappingMatrix.assemblyEnd()
          self.MappingMatrix_T.assemblyBegin()
          self.MappingMatrix_T.assemblyEnd()
          self.MPIPrint("Interpolation matrix is built.")
  
        self.MPIBarrier()
  
        del self.localSolidInterface_array_X
        del self.localSolidInterface_array_Y
        del self.localSolidInterface_array_Z

    def matchingMeshMapping(self,solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Fill the mapping matrix in case of matching meshes at the f/s interface.
        """
        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0

        # --- Instantiate the spatial indexing ---
	prop_index = index.Property()
	prop_index.dimension = self.nDim
	SolidSpatialTree = index.Index(properties=prop_index)
        
        nSolidNodes = solidInterfaceBuffRcv_X.shape[0]

        for jVertex in range(nSolidNodes):
          posX = solidInterfaceBuffRcv_X[jVertex]
          posY = solidInterfaceBuffRcv_Y[jVertex]
          posZ = solidInterfaceBuffRcv_Z[jVertex]
	  if self.nDim == 2 :
	    SolidSpatialTree.add(jVertex, (posX, posY))
   	  else :
	    SolidSpatialTree.add(jVertex, (posX, posY, posZ))

        if self.nFluidInterfacePhysicalNodes != self.nSolidInterfacePhysicalNodes:
          raise Exception("Fluid and solid interface must have the same number of nodes for matching meshes ! ")

        # --- For each fluid interface node, find the nearest solid interface node and fill the boolean mapping matrix ---
        for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
          posX = self.localFluidInterface_array_X_init[iVertexFluid]
          posY = self.localFluidInterface_array_Y_init[iVertexFluid]
          posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
          if self.nDim == 2:
            neighboors = list(SolidSpatialTree.nearest((posX, posY),1))
          elif self.nDim == 3:
            neighboors = list(SolidSpatialTree.nearest((posX, posY, posZ),1))
          jVertexSolid = neighboors[0]
          # Check if the distance is small enough to ensure coincidence
          NodeA = np.array([posX, posY, posZ])
          NodeB = np.array([solidInterfaceBuffRcv_X[jVertexSolid], solidInterfaceBuffRcv_Y[jVertexSolid], solidInterfaceBuffRcv_Z[jVertexSolid]])
          distance = spdist.euclidean(NodeA, NodeB)
          iGlobalVertexFluid = self.__getGlobalIndex('fluid', myid, iVertexFluid)
          jGlobalVertexSolid = self.__getGlobalIndex('solid', iProc, jVertexSolid)
          if distance > 1e-6:
            print("WARNING : Tolerance for matching meshes is not matched between node F{} and S{} : ({}, {}, {})<-->({}, {}, {}) , DISTANCE : {} !".format(iGlobalVertexFluid,jGlobalVertexSolid,posX, posY, posZ,solidInterfaceBuffRcv_X[jVertexSolid], solidInterfaceBuffRcv_Y[jVertexSolid], solidInterfaceBuffRcv_Z[jVertexSolid], distance))
          self.MappingMatrix.setValue(iGlobalVertexFluid,jGlobalVertexSolid,1.0)
          self.MappingMatrix_T.setValue(jGlobalVertexSolid, iGlobalVertexFluid,1.0)

        del solidInterfaceBuffRcv_X
        del solidInterfaceBuffRcv_Y
        del solidInterfaceBuffRcv_Z

    def NearestNeighboorMeshMapping(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description
        """

        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0

        # --- Instantiate the spatial indexing ---
	prop_index = index.Property()
	prop_index.dimension = self.nDim
	SolidSpatialTree = index.Index(properties=prop_index)
        
        nSolidNodes = solidInterfaceBuffRcv_X.shape[0]

        for jVertex in range(nSolidNodes):
          posX = solidInterfaceBuffRcv_X[jVertex]
          posY = solidInterfaceBuffRcv_Y[jVertex]
          posZ = solidInterfaceBuffRcv_Z[jVertex]
	  if self.nDim == 2 :
	    SolidSpatialTree.add(jVertex, (posX, posY))
   	  else :
	    SolidSpatialTree.add(jVertex, (posX, posY, posZ))

        # --- For each fluid interface node, find the nearest solid interface node and fill the boolean mapping matrix ---
        for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
          posX = self.localFluidInterface_array_X_init[iVertexFluid]
          posY = self.localFluidInterface_array_Y_init[iVertexFluid]
          posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
          if self.nDim == 2:
            neighboors = list(SolidSpatialTree.nearest((posX, posY),1))
          elif self.nDim == 3:
            neighboors = list(SolidSpatialTree.nearest((posX, posY, posZ),1))
          jVertexSolid = neighboors[0]
          iGlobalVertexFluid = self.__getGlobalIndex('fluid', myid, iVertexFluid)
          jGlobalVertexSolid = self.__getGlobalIndex('solid', iProc, jVertexSolid)
          self.MappingMatrix.setValue(iGlobalVertexFluid,jGlobalVertexSolid,1.0)
          self.MappingMatrix_T.setValue(jGlobalVertexSolid, iGlobalVertexFluid,1.0)

    def RBFMeshMapping_A(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc, rad):
        """
        Description
        """

        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0

        # --- Instantiate the spatial indexing ---
	prop_index = index.Property()
	prop_index.dimension = self.nDim
	SolidSpatialTree = index.Index(properties=prop_index)
        
        nSolidNodes = solidInterfaceBuffRcv_X.shape[0]

        for jVertex in range(nSolidNodes):
          posX = solidInterfaceBuffRcv_X[jVertex]
          posY = solidInterfaceBuffRcv_Y[jVertex]
          posZ = solidInterfaceBuffRcv_Z[jVertex]
	  if self.nDim == 2 :
	    SolidSpatialTree.add(jVertex, (posX, posY))
   	  else :
	    SolidSpatialTree.add(jVertex, (posX, posY, posZ))

        for iVertexSolid in range(self.nLocalSolidInterfaceNodes):
          posX = self.localSolidInterface_array_X[iVertexSolid]
          posY = self.localSolidInterface_array_Y[iVertexSolid]
          posZ = self.localSolidInterface_array_Z[iVertexSolid]
          NodeA = np.array([posX, posY, posZ])
          iGlobalVertexSolid = self.__getGlobalIndex('solid', myid, iVertexSolid)
          if self.nDim == 2:
            neighboors = list(SolidSpatialTree.intersection((posX-rad, posY-rad, posX+rad, posY+rad)))
          elif self.nDim == 3:
            neighboors = list(SolidSpatialTree.intersection((posX-rad, posY-rad, posZ-rad, posX+rad, posY+rad, posZ+rad)))
          for jVertexSolid in neighboors:
            NodeB = np.array([solidInterfaceBuffRcv_X[jVertexSolid], solidInterfaceBuffRcv_Y[jVertexSolid], solidInterfaceBuffRcv_Z[jVertexSolid]])
            distance = spdist.euclidean(NodeA, NodeB)
            phi = self.__CPC2(distance, rad)
            jGlobalVertexSolid = self.__getGlobalIndex('solid', iProc, jVertexSolid)
            self.MappingMatrixA.setValue(iGlobalVertexSolid, jGlobalVertexSolid, phi)
            self.MappingMatrixA_T.setValue(jGlobalVertexSolid, iGlobalVertexSolid, phi)
          self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes, 1.0)
          self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes+1, posX)
          self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes+2, posY)
          self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes+3, posZ)
          self.MappingMatrixA_T.setValue(nSolidNodes, iGlobalVertexSolid, 1.0)
          self.MappingMatrixA_T.setValue(nSolidNodes+1, iGlobalVertexSolid, posX)
          self.MappingMatrixA_T.setValue(nSolidNodes+2, iGlobalVertexSolid, posY)
          self.MappingMatrixA_T.setValue(nSolidNodes+3, iGlobalVertexSolid, posZ)

    def RBFMeshMapping_B(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc, rad):
        """
        Description
        """

        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0

        # --- Instantiate the spatial indexing ---
	prop_index = index.Property()
	prop_index.dimension = self.nDim
	SolidSpatialTree = index.Index(properties=prop_index)
        
        nSolidNodes = solidInterfaceBuffRcv_X.shape[0]

        for jVertex in range(nSolidNodes):
          posX = solidInterfaceBuffRcv_X[jVertex]
          posY = solidInterfaceBuffRcv_Y[jVertex]
          posZ = solidInterfaceBuffRcv_Z[jVertex]
	  if self.nDim == 2 :
	    SolidSpatialTree.add(jVertex, (posX, posY))
   	  else :
	    SolidSpatialTree.add(jVertex, (posX, posY, posZ))

        for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
          posX = self.localFluidInterface_array_X_init[iVertexFluid]
          posY = self.localFluidInterface_array_Y_init[iVertexFluid]
          posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
          NodeA = np.array([posX, posY, posZ])
          iGlobalVertexFluid = self.__getGlobalIndex('fluid', myid, iVertexFluid)
          if self.nDim == 2:
            neighboors = list(SolidSpatialTree.intersection((posX-rad, posY-rad, posX+rad, posY+rad)))
          elif self.nDim == 3:
            neighboors = list(SolidSpatialTree.intersection((posX-rad, posY-rad, posZ-rad, posX+rad, posY+rad, posZ+rad)))
          for jVertexSolid in neighboors:
            NodeB = np.array([solidInterfaceBuffRcv_X[jVertexSolid], solidInterfaceBuffRcv_Y[jVertexSolid], solidInterfaceBuffRcv_Z[jVertexSolid]])
            distance = spdist.euclidean(NodeA, NodeB)
            phi = self.__CPC2(distance, rad)
            jGlobalVertexSolid = self.__getGlobalIndex('solid', iProc, jVertexSolid)
            self.MappingMatrixB.setValue(iGlobalVertexFluid, jGlobalVertexSolid, phi)
            self.MappingMatrixB_T.setValue(jGlobalVertexSolid, iGlobalVertexFluid, phi)
          self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes, 1.0)
          self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes+1, posX)
          self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes+2, posY)
          self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes+3, posZ)
          self.MappingMatrixB_T.setValue(nSolidNodes, iGlobalVertexFluid, 1.0)
          self.MappingMatrixB_T.setValue(nSolidNodes+1, iGlobalVertexFluid, posX)
          self.MappingMatrixB_T.setValue(nSolidNodes+2, iGlobalVertexFluid, posY)
          self.MappingMatrixB_T.setValue(nSolidNodes+3, iGlobalVertexFluid, posZ)

    def TPSMeshMapping_A(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description
        """

        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0
        
        nSolidNodes = solidInterfaceBuffRcv_X.shape[0]

        for iVertexSolid in range(self.nLocalSolidInterfaceNodes):
          posX = self.localSolidInterface_array_X[iVertexSolid]
          posY = self.localSolidInterface_array_Y[iVertexSolid]
          posZ = self.localSolidInterface_array_Z[iVertexSolid]
          NodeA = np.array([posX, posY, posZ])
          iGlobalVertexSolid = self.__getGlobalIndex('solid', myid, iVertexSolid)
          for jVertexSolid in range(nSolidNodes):
            NodeB = np.array([solidInterfaceBuffRcv_X[jVertexSolid], solidInterfaceBuffRcv_Y[jVertexSolid], solidInterfaceBuffRcv_Z[jVertexSolid]])
            distance = spdist.euclidean(NodeA, NodeB)
            phi = self.__TPS(distance)
            jGlobalVertexSolid = self.__getGlobalIndex('solid', iProc, jVertexSolid)
            self.MappingMatrixA.setValue(iGlobalVertexSolid, jGlobalVertexSolid, phi)
            self.MappingMatrixA_T.setValue(jGlobalVertexSolid, iGlobalVertexSolid, phi)
          self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes, 1.0)
          self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes+1, posX)
          self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes+2, posY)
          self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes+3, posZ)
          self.MappingMatrixA_T.setValue(nSolidNodes, iGlobalVertexSolid, 1.0)
          self.MappingMatrixA_T.setValue(nSolidNodes+1, iGlobalVertexSolid, posX)
          self.MappingMatrixA_T.setValue(nSolidNodes+2, iGlobalVertexSolid, posY)
          self.MappingMatrixA_T.setValue(nSolidNodes+3, iGlobalVertexSolid, posZ)

    def TPSMeshMapping_B(self, solidInterfaceBuffRcv_X, solidInterfaceBuffRcv_Y, solidInterfaceBuffRcv_Z, iProc):
        """
        Description
        """

        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0
        
        nSolidNodes = solidInterfaceBuffRcv_X.shape[0]

        for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
          posX = self.localFluidInterface_array_X_init[iVertexFluid]
          posY = self.localFluidInterface_array_Y_init[iVertexFluid]
          posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
          NodeA = np.array([posX, posY, posZ])
          iGlobalVertexFluid = self.__getGlobalIndex('fluid', myid, iVertexFluid)
          for jVertexSolid in range(nSolidNodes):
            NodeB = np.array([solidInterfaceBuffRcv_X[jVertexSolid], solidInterfaceBuffRcv_Y[jVertexSolid], solidInterfaceBuffRcv_Z[jVertexSolid]])
            distance = spdist.euclidean(NodeA, NodeB)
            phi = self.__TPS(distance)
            jGlobalVertexSolid = self.__getGlobalIndex('solid', iProc, jVertexSolid)
            self.MappingMatrixB.setValue(iGlobalVertexFluid, jGlobalVertexSolid, phi)
            self.MappingMatrixB_T.setValue(jGlobalVertexSolid, iGlobalVertexFluid, phi)
          self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes, 1.0)
          self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes+1, posX)
          self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes+2, posY)
          self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes+3, posZ)
          self.MappingMatrixB_T.setValue(nSolidNodes, iGlobalVertexFluid, 1.0)
          self.MappingMatrixB_T.setValue(nSolidNodes+1, iGlobalVertexFluid, posX)
          self.MappingMatrixB_T.setValue(nSolidNodes+2, iGlobalVertexFluid, posY)
          self.MappingMatrixB_T.setValue(nSolidNodes+3, iGlobalVertexFluid, posZ)


    def __CPC2(self, distance, rad):
        """
        Description.
        """
        phi = 0.0
        eps = distance/rad
 
        if eps < 1:
          phi = ((1.0-eps)**4)*(4.0*eps+1.0)
        else:
          phi = 0.0

        return phi

    def __TPS(self, distance):
        """
        Description
        """
        phi = 0.0
     
        if distance > 0.0:
          phi = (distance**2)*np.log10(distance)
        else:
          phi = 0.0

        return phi            


    def interpolateSolidPositionOnFluidMesh(self, FSI_config):
	"""
	Applies the one-to-one mapping or the interpolaiton rules from solid to fluid mesh.
	"""
	if self.have_MPI == True:
          myid = self.comm.Get_rank()
          MPIsize = self.comm.Get_size()
        else:
          myid = 0
          MPIsize = 1


        # --- Interpolate (or map) in parallel the solid interface displacement on the fluid interface ---
        if FSI_config['MATCHING_MESH'] == 'NO' and (FSI_config['MESH_INTERP_METHOD'] == 'RBF' or FSI_config['MESH_INTERP_METHOD'] == 'TPS'):
          if self.have_MPI == True:
            gamma_array_DispX = PETSc.Vec().create(self.comm)
            gamma_array_DispY = PETSc.Vec().create(self.comm)
            gamma_array_DispZ = PETSc.Vec().create(self.comm)
            gamma_array_DispX.setType('mpi')
            gamma_array_DispY.setType('mpi')
            gamma_array_DispZ.setType('mpi')
            KSP_solver = PETSc.KSP().create(self.comm)
          else:
            gamma_array_DispX = PETSc.Vec().create()
            gamma_array_DispY = PETSc.Vec().create()
            gamma_array_DispZ = PETSc.Vec().create()
            gamma_array_DispX.setType('seq')
            gamma_array_DispY.setType('seq')
            gamma_array_DispZ.setType('seq')
            KSP_solver = PETSc.KSP().create()
          gamma_array_DispX.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
          gamma_array_DispY.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
          gamma_array_DispZ.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
          gamma_array_DispX.set(0.0)
          gamma_array_DispY.set(0.0)
          gamma_array_DispZ.set(0.0)
          KSP_solver.setType('fgmres')
          KSP_solver.getPC().setType('jacobi')
          KSP_solver.setOperators(self.MappingMatrixA)
          KSP_solver.setFromOptions()
          #print(KSP_solver.getInitialGuessNonzero())
          KSP_solver.setInitialGuessNonzero(True)
          #print(KSP_solver.getInitialGuessNonzero())
          KSP_solver.solve(self.solidInterface_array_DispX, gamma_array_DispX)
          KSP_solver.solve(self.solidInterface_array_DispY, gamma_array_DispY)
          KSP_solver.solve(self.solidInterface_array_DispZ, gamma_array_DispZ)
          self.MappingMatrixB.mult(gamma_array_DispX, self.fluidInterface_array_DispX)
          self.MappingMatrixB.mult(gamma_array_DispY, self.fluidInterface_array_DispY)
          self.MappingMatrixB.mult(gamma_array_DispZ, self.fluidInterface_array_DispZ)
          gamma_array_DispX.destroy()
          gamma_array_DispY.destroy()
          gamma_array_DispZ.destroy()
          KSP_solver.destroy()
          del gamma_array_DispX
          del gamma_array_DispY
          del gamma_array_DispZ
          del KSP_solver
        else:  
          self.MappingMatrix.mult(self.solidInterface_array_DispX, self.fluidInterface_array_DispX)
          self.MappingMatrix.mult(self.solidInterface_array_DispY, self.fluidInterface_array_DispY)
          self.MappingMatrix.mult(self.solidInterface_array_DispZ, self.fluidInterface_array_DispZ)

        # --- Checking conservation --- 
        WSX = self.solidLoads_array_X.dot(self.solidInterface_array_DispX)
        WSY = self.solidLoads_array_Y.dot(self.solidInterface_array_DispY)
        WSZ = self.solidLoads_array_Z.dot(self.solidInterface_array_DispZ)

        WFX = self.fluidLoads_array_X.dot(self.fluidInterface_array_DispX)
        WFY = self.fluidLoads_array_Y.dot(self.fluidInterface_array_DispY)
        WFZ = self.fluidLoads_array_Z.dot(self.fluidInterface_array_DispZ)

	self.MPIPrint("Checking f/s interface conservation...")
	self.MPIPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ))        
	self.MPIPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ))

   
        # --- Redistribute the interpolated fluid interface according to the partitions that own the fluid interface ---
        # Gather the fluid interface on the master process
        if self.have_MPI == True:
          sendBuff_X = None
          sendBuff_Y = None
          sendBuff_Z = None
          self.fluidInterface_array_DispX_recon = None
          self.fluidInterface_array_DispY_recon = None
          self.fluidInterface_array_DispZ_recon = None

          if myid == self.rootProcess:
            self.fluidInterface_array_DispX_recon = np.zeros(self.nFluidInterfacePhysicalNodes)
            self.fluidInterface_array_DispY_recon = np.zeros(self.nFluidInterfacePhysicalNodes)
            self.fluidInterface_array_DispZ_recon = np.zeros(self.nFluidInterfacePhysicalNodes)

          myNumberOfNodes = self.fluidInterface_array_DispX.getArray().shape[0]
          sendBuffNumber = np.array([myNumberOfNodes], dtype=int)
          rcvBuffNumber = np.zeros(MPIsize, dtype=int)
          self.comm.Allgather(sendBuffNumber, rcvBuffNumber)

          counts = tuple(rcvBuffNumber)
          displ = np.zeros(MPIsize, dtype=int)
          for ii in range(rcvBuffNumber.shape[0]):
            displ[ii] = rcvBuffNumber[0:ii].sum()
          displ = tuple(displ)

          del sendBuffNumber, rcvBuffNumber
              
          #print("DEBUG MESSAGE From proc {}, counts = {}".format(myid, counts))
          #print("DEBUG MESSAGE From proc {}, displ = {}".format(myid, displ))

          self.comm.Gatherv(self.fluidInterface_array_DispX.getArray(), [self.fluidInterface_array_DispX_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)
          self.comm.Gatherv(self.fluidInterface_array_DispY.getArray(), [self.fluidInterface_array_DispY_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)
          self.comm.Gatherv(self.fluidInterface_array_DispZ.getArray(), [self.fluidInterface_array_DispZ_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)

          # Send the partitioned interface to the right fluid partitions
          if myid == self.rootProcess:
            for iProc in self.fluidInterfaceProcessors:
              sendBuff_X = np.zeros(self.fluidPhysicalInterfaceNodesDistribution[iProc])
              sendBuff_Y = np.zeros(self.fluidPhysicalInterfaceNodesDistribution[iProc])
              sendBuff_Z = np.zeros(self.fluidPhysicalInterfaceNodesDistribution[iProc])
              globalIndex = self.fluidGlobalIndexRange[iProc][iProc][0]
              for iVertex in range(self.fluidPhysicalInterfaceNodesDistribution[iProc]):
                sendBuff_X[iVertex] = self.fluidInterface_array_DispX_recon[globalIndex]
                sendBuff_Y[iVertex] = self.fluidInterface_array_DispY_recon[globalIndex]
                sendBuff_Z[iVertex] = self.fluidInterface_array_DispZ_recon[globalIndex]
                globalIndex += 1
              self.comm.Send(sendBuff_X, dest=iProc, tag = 1)
              self.comm.Send(sendBuff_Y, dest=iProc, tag = 2)
              self.comm.Send(sendBuff_Z, dest=iProc, tag = 3)
          if myid in self.fluidInterfaceProcessors:
            self.localFluidInterface_array_DispX = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
            self.localFluidInterface_array_DispY = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
            self.localFluidInterface_array_DispZ = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
            self.comm.Recv(self.localFluidInterface_array_DispX, source=self.rootProcess, tag = 1)
            self.comm.Recv(self.localFluidInterface_array_DispY, source=self.rootProcess, tag = 2)
            self.comm.Recv(self.localFluidInterface_array_DispZ, source=self.rootProcess, tag = 3)
          del sendBuff_X
          del sendBuff_Y
          del sendBuff_Z
        else:
          self.localFluidInterface_array_DispX = self.fluidInterface_array_DispX.getArray().copy()
          self.localFluidInterface_array_DispY = self.fluidInterface_array_DispY.getArray().copy()
          self.localFluidInterface_array_DispZ = self.fluidInterface_array_DispZ.getArray().copy()

        # Special treatment for the halo nodes on the fluid interface
        self.haloNodesDisplacements = {}
        sendBuff = {}
        if self.have_MPI == True:
          if myid == self.rootProcess:
            for iProc in self.fluidInterfaceProcessors:
              sendBuff = {}
              for key in self.FluidHaloNodeList[iProc].keys():
                globalIndex = self.fluidIndexing[key]
                DispX = self.fluidInterface_array_DispX_recon[globalIndex]
                DispY = self.fluidInterface_array_DispY_recon[globalIndex]
                DispZ = self.fluidInterface_array_DispZ_recon[globalIndex]
                sendBuff[key] = (DispX, DispY, DispZ)
              self.comm.send(sendBuff, dest = iProc, tag=4)
          if myid in self.fluidInterfaceProcessors:
            self.haloNodesDisplacements = self.comm.recv(source = self.rootProcess, tag = 4)
          del sendBuff

    def interpolateFluidLoadsOnSolidMesh(self, FSI_config):
	"""
	Applies the one-to-one mapping or the interpolaiton rules from fluid to solid mesh.
	"""
	if self.have_MPI == True:
          myid = self.comm.Get_rank()
          MPIsize = self.comm.Get_size()
        else:
          myid = 0
          MPIsize = 1
	
        # --- Interpolate (or map) in parallel the fluid interface loads on the solid interface ---
	#self.MappingMatrix.transpose()
        if FSI_config['MATCHING_MESH'] == 'NO' and (FSI_config['MESH_INTERP_METHOD'] == 'RBF' or FSI_config['MESH_INTERP_METHOD'] == 'TPS'):
          if self.have_MPI == True:
            gamma_array_LoadX = PETSc.Vec().create(self.comm)
            gamma_array_LoadY = PETSc.Vec().create(self.comm)
            gamma_array_LoadZ = PETSc.Vec().create(self.comm)
            gamma_array_LoadX.setType('mpi')
            gamma_array_LoadY.setType('mpi')
            gamma_array_LoadZ.setType('mpi')
            KSP_solver = PETSc.KSP().create(self.comm)
          else:
            gamma_array_LoadX = PETSc.Vec().create()
            gamma_array_LoadY = PETSc.Vec().create()
            gamma_array_LoadZ = PETSc.Vec().create()
            gamma_array_LoadX.setType('seq')
            gamma_array_LoadY.setType('seq')
            gamma_array_LoadZ.setType('seq')
            KSP_solver = PETSc.KSP().create()
          gamma_array_LoadX.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
          gamma_array_LoadY.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
          gamma_array_LoadZ.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
          gamma_array_LoadX.set(0.0)
          gamma_array_LoadY.set(0.0)
          gamma_array_LoadZ.set(0.0)
          KSP_solver.setType('fgmres')
          KSP_solver.getPC().setType('jacobi')
          KSP_solver.setOperators(self.MappingMatrixA_T)
          KSP_solver.setFromOptions()
          self.MappingMatrixB_T.mult(self.fluidLoads_array_X, gamma_array_LoadX)
          self.MappingMatrixB_T.mult(self.fluidLoads_array_Y, gamma_array_LoadY)
          self.MappingMatrixB_T.mult(self.fluidLoads_array_Z, gamma_array_LoadZ)
          KSP_solver.solve(gamma_array_LoadX, self.solidLoads_array_X)
          KSP_solver.solve(gamma_array_LoadY, self.solidLoads_array_Y)
          KSP_solver.solve(gamma_array_LoadZ, self.solidLoads_array_Z)
          gamma_array_LoadX.destroy()
          gamma_array_LoadY.destroy()
          gamma_array_LoadZ.destroy()
          KSP_solver.destroy()
          del gamma_array_LoadX
          del gamma_array_LoadY
          del gamma_array_LoadZ
          del KSP_solver
        else:
          self.MappingMatrix_T.mult(self.fluidLoads_array_X, self.solidLoads_array_X)
          self.MappingMatrix_T.mult(self.fluidLoads_array_Y, self.solidLoads_array_Y)
          self.MappingMatrix_T.mult(self.fluidLoads_array_Z, self.solidLoads_array_Z)

        # --- Redistribute the interpolated solid loads according to the partitions that own the solid interface ---
        # Gather the solid loads on the master process
        if self.have_MPI:
          sendBuff_X = None
          sendBuff_Y = None
          sendBuff_Z = None
          self.solidLoads_array_X_recon = None
          self.solidLoads_array_Y_recon = None
          self.solidLoads_array_Z_recon = None
	  if myid == self.rootProcess:
	    self.solidLoads_array_X_recon = np.zeros(self.nSolidInterfacePhysicalNodes+self.d_RBF)
	    self.solidLoads_array_Y_recon = np.zeros(self.nSolidInterfacePhysicalNodes+self.d_RBF)
	    self.solidLoads_array_Z_recon = np.zeros(self.nSolidInterfacePhysicalNodes+self.d_RBF)
          myNumberOfNodes =  self.solidLoads_array_X.getArray().shape[0]
          sendBuffNumber = np.array([myNumberOfNodes], dtype=int)
          rcvBuffNumber = np.zeros(MPIsize, dtype=int)
          self.comm.Allgather(sendBuffNumber, rcvBuffNumber)

          counts = tuple(rcvBuffNumber)
          displ = np.zeros(MPIsize, dtype=int)
          for ii in range(rcvBuffNumber.shape[0]):
            displ[ii] = rcvBuffNumber[0:ii].sum()
          displ = tuple(displ)      

          del sendBuffNumber, rcvBuffNumber   

          self.comm.Gatherv(self.solidLoads_array_X.getArray(), [self.solidLoads_array_X_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)
          self.comm.Gatherv(self.solidLoads_array_Y.getArray(), [self.solidLoads_array_Y_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)
          self.comm.Gatherv(self.solidLoads_array_Z.getArray(), [self.solidLoads_array_Z_recon, counts, displ, self.MPI.DOUBLE], root=self.rootProcess)

          # Send the partitioned loads to the right solid partitions
          if myid == self.rootProcess:
            for iProc in self.solidInterfaceProcessors:
              sendBuff_X = np.zeros(self.solidPhysicalInterfaceNodesDistribution[iProc])
              sendBuff_Y = np.zeros(self.solidPhysicalInterfaceNodesDistribution[iProc])
              sendBuff_Z = np.zeros(self.solidPhysicalInterfaceNodesDistribution[iProc])
              globalIndex = self.solidGlobalIndexRange[iProc][iProc][0]
              for iVertex in range(self.solidPhysicalInterfaceNodesDistribution[iProc]):
                sendBuff_X[iVertex] = self.solidLoads_array_X_recon[globalIndex]
                sendBuff_Y[iVertex] = self.solidLoads_array_Y_recon[globalIndex]
                sendBuff_Z[iVertex] = self.solidLoads_array_Z_recon[globalIndex]
                globalIndex += 1
              self.comm.Send(sendBuff_X, dest=iProc, tag = 1)
              self.comm.Send(sendBuff_Y, dest=iProc, tag = 2)
              self.comm.Send(sendBuff_Z, dest=iProc, tag = 3)
          if myid in self.solidInterfaceProcessors:
            self.localSolidLoads_array_X = np.zeros(self.nLocalSolidInterfaceNodes)
            self.localSolidLoads_array_Y = np.zeros(self.nLocalSolidInterfaceNodes)
            self.localSolidLoads_array_Z = np.zeros(self.nLocalSolidInterfaceNodes)
            self.comm.Recv(self.localSolidLoads_array_X, source=self.rootProcess, tag = 1)
            self.comm.Recv(self.localSolidLoads_array_Y, source=self.rootProcess, tag = 2)
            self.comm.Recv(self.localSolidLoads_array_Z, source=self.rootProcess, tag = 3)
          del sendBuff_X
          del sendBuff_Y
          del sendBuff_Z
        else:
          self.localSolidLoads_array_X = self.solidLoads_array_X.getArray().copy()
          self.localSolidLoads_array_Y = self.solidLoads_array_Y.getArray().copy()
          self.localSolidLoads_array_Z = self.solidLoads_array_Z.getArray().copy()

        # Special treatment for the halo nodes on the fluid interface
        # TODO when we will use parallel solid solver !!


    '''def getSolidInterfacePosition(self, SolidSolver):
	"""
	Gets the current solid interface position from the solid solver.
	"""
        if self.have_MPI == True:
	  myid = self.comm.Get_rank()
        else:
          myid = 0
	
        # --- Get the solid interface position from the solid solver and directly fill the corresponding PETSc vector ---
        GlobalIndex = int()
        localIndex = 0
	for iVertex in range(self.nLocalSolidInterfaceNodes):
          GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
          if GlobalIndex in self.SolidHaloNodeList[myid].keys():
            pass
          else:
	    newPosx = SolidSolver.getInterfaceNodePosX(self.solidInterfaceIdentifier, iVertex)
	    newPosy = SolidSolver.getInterfaceNodePosY(self.solidInterfaceIdentifier, iVertex)
	    newPosz = SolidSolver.getInterfaceNodePosZ(self.solidInterfaceIdentifier, iVertex)
            iGlobalVertex = self.__getGlobalIndex('solid', myid, localIndex)
            self.solidInterface_array_X.setValues([iGlobalVertex],newPosx)
            self.solidInterface_array_Y.setValues([iGlobalVertex],newPosy)
            self.solidInterface_array_Z.setValues([iGlobalVertex],newPosz)
            localIndex += 1
          #print("DEBUG MESSAGE From proc {} : In loop !".format(myid))

        #print("DEBUG MESSAGE From proc {} : Prepare for assembly !".format(myid))

        self.solidInterface_array_X.assemblyBegin()
        self.solidInterface_array_X.assemblyEnd()
        self.solidInterface_array_Y.assemblyBegin()
        self.solidInterface_array_Y.assemblyEnd()
        self.solidInterface_array_Z.assemblyBegin()
        self.solidInterface_array_Z.assemblyEnd()

        #print("DEBUG MESSAGE From PROC {} : Assembly is done !".format(myid))
        #print("DEBUG MESSAGE From PROC {} : array_X = {}".format(myid, self.solidInterface_array_X.getArray()))'''

    def getSolidInterfaceDisplacement(self, SolidSolver):
	"""
	Gets the current solid interface position from the solid solver.
	"""
        if self.have_MPI == True:
	  myid = self.comm.Get_rank()
        else:
          myid = 0
	
        # --- Get the solid interface position from the solid solver and directly fill the corresponding PETSc vector ---
        GlobalIndex = int()
        localIndex = 0
	for iVertex in range(self.nLocalSolidInterfaceNodes):
          GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
          if GlobalIndex in self.SolidHaloNodeList[myid].keys():
            pass
          else:
	    newDispx = SolidSolver.getInterfaceNodeDispX(self.solidInterfaceIdentifier, iVertex)
	    newDispy = SolidSolver.getInterfaceNodeDispY(self.solidInterfaceIdentifier, iVertex)
	    newDispz = SolidSolver.getInterfaceNodeDispZ(self.solidInterfaceIdentifier, iVertex)
            iGlobalVertex = self.__getGlobalIndex('solid', myid, localIndex)
            self.solidInterface_array_DispX.setValues([iGlobalVertex],newDispx)
            self.solidInterface_array_DispY.setValues([iGlobalVertex],newDispy)
            self.solidInterface_array_DispZ.setValues([iGlobalVertex],newDispz)
            localIndex += 1

        self.solidInterface_array_DispX.assemblyBegin()
        self.solidInterface_array_DispX.assemblyEnd()
        self.solidInterface_array_DispY.assemblyBegin()
        self.solidInterface_array_DispY.assemblyEnd()
        self.solidInterface_array_DispZ.assemblyBegin()
        self.solidInterface_array_DispZ.assemblyEnd()

    def getFluidInterfaceNodalForce(self, FSI_config, FluidSolver):
	"""
	Gets the fluid interface loads from the fluid solver.
	"""
        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0

        localIndex = 0
        FX = 0.0
        FY = 0.0
        FZ = 0.0

        # --- Get the fluid interface loads from the fluid solver and directly fill the corresponding PETSc vector ---
	for iVertex in range(self.nLocalFluidInterfaceNodes):
	    halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex) # !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
	    if halo==False:
		if FSI_config['CSD_SOLVER'] == 'GETDP':
		    newFx = FluidSolver.GetVertexForceDensityX(self.fluidInterfaceIdentifier, iVertex)
	            newFy = FluidSolver.GetVertexForceDensityY(self.fluidInterfaceIdentifier, iVertex)
	            newFz = FluidSolver.GetVertexForceDensityZ(self.fluidInterfaceIdentifier, iVertex)
		else:
	            newFx = FluidSolver.GetVertexForceX(self.fluidInterfaceIdentifier, iVertex)
	            newFy = FluidSolver.GetVertexForceY(self.fluidInterfaceIdentifier, iVertex)
	            newFz = FluidSolver.GetVertexForceZ(self.fluidInterfaceIdentifier, iVertex)
                iGlobalVertex = self.__getGlobalIndex('fluid', myid, localIndex)
                self.fluidLoads_array_X.setValues([iGlobalVertex], newFx)
                self.fluidLoads_array_Y.setValues([iGlobalVertex], newFy)
                self.fluidLoads_array_Z.setValues([iGlobalVertex], newFz)
                FX += newFx
                FY += newFy
                FZ += newFz
                localIndex += 1

        if self.have_MPI == True:
          FX = self.comm.allreduce(FX)
          FY = self.comm.allreduce(FY)
          FZ = self.comm.allreduce(FZ)

        self.fluidLoads_array_X.assemblyBegin()
        self.fluidLoads_array_X.assemblyEnd()
        self.fluidLoads_array_Y.assemblyBegin()
        self.fluidLoads_array_Y.assemblyEnd()
        self.fluidLoads_array_Z.assemblyBegin()
        self.fluidLoads_array_Z.assemblyEnd()

        FX_b = self.fluidLoads_array_X.sum()
        FY_b = self.fluidLoads_array_Y.sum()
        FZ_b = self.fluidLoads_array_Z.sum()
        

    def setFluidInterfaceVarCoord(self, FluidSolver):
	"""
	Communicate the change of coordinates of the fluid interface to the fluid solver.
	Prepare the fluid solver for mesh deformation.
	"""
        if self.have_MPI == True:
	  myid = self.comm.Get_rank()
        else:
          myid = 0
	
        # --- Send the new fluid interface position to the fluid solver (on each partition, halo nodes included) ---
        localIndex = 0
	for iVertex in range(self.nLocalFluidInterfaceNodes):
	    GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
            if GlobalIndex in self.FluidHaloNodeList[myid].keys():
              posX0, posY0, posZ0 = self.haloNodesPositionsInit[GlobalIndex]
              DispX, DispY, DispZ = self.haloNodesDisplacements[GlobalIndex]
              posX = posX0 + DispX
              posY = posY0 + DispY
              posZ = posZ0 + DispZ
              FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, posX)
              FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, posY)
              FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, posZ)
            else:
              posX = self.localFluidInterface_array_DispX[localIndex] + self.localFluidInterface_array_X_init[localIndex]
              posY = self.localFluidInterface_array_DispY[localIndex] + self.localFluidInterface_array_Y_init[localIndex]
              posZ = self.localFluidInterface_array_DispZ[localIndex] + self.localFluidInterface_array_Z_init[localIndex]
              FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, posX)
              FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, posY)
              FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, posZ)
              localIndex += 1
            # Prepares the mesh deformation in the fluid solver
	    nodalVarCoordNorm = FluidSolver.SetVertexVarCoord(self.fluidInterfaceIdentifier, iVertex)

	    
    def setSolidInterfaceLoads(self, SolidSolver, FSI_config, time):
	"""
	Communicates the new solid interface loads to the solid solver.
	In case of rigid body motion, calculates the new resultant forces (lift, drag, ...).
	"""
        if self.have_MPI == True:
	  myid = self.comm.Get_rank()
        else:
          myid = 0

	FY = 0.0 # solid-side resultant forces
        FX = 0.0
        FZ = 0.0
	FFX = 0.0 # fluid-side resultant forces
	FFY = 0.0
	FFZ = 0.0

        # --- Check for total force conservation after interpolation
        FFX = self.fluidLoads_array_X.sum()
        FFY = self.fluidLoads_array_Y.sum()
        FFZ = self.fluidLoads_array_Z.sum()

 	
        for iVertex in range(self.nLocalSolidInterfaceNodes):
          FX += self.localSolidLoads_array_X[iVertex]
          FY += self.localSolidLoads_array_Y[iVertex]
          FZ += self.localSolidLoads_array_Z[iVertex]

        if self.have_MPI == True:
          FX = self.comm.allreduce(FX)
          FY = self.comm.allreduce(FY)
          FZ = self.comm.allreduce(FZ)

	self.MPIPrint("Checking f/s interface total force...")
	self.MPIPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FX, FY, FZ))        
	self.MPIPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ))

        # --- Send the new solid interface loads to the solid solver (on each partition, halo nodes included) ---
        GlobalIndex = int()
        localIndex = 0
        if myid in self.solidInterfaceProcessors:
          for iVertex in range(self.nLocalSolidInterfaceNodes):
            GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
            if GlobalIndex in self.SolidHaloNodeList[myid].keys():
              pass
            else:
              Fx = self.localSolidLoads_array_X[localIndex]
	      Fy = self.localSolidLoads_array_Y[localIndex]
	      Fz = self.localSolidLoads_array_Z[localIndex]
              SolidSolver.applyload(iVertex, Fx, Fy, Fz, time)
              localIndex += 1
          if FSI_config['CSD_SOLVER'] == 'NATIVE':  
            SolidSolver.setGeneralisedForce()
	    SolidSolver.setGeneralisedMoment()

    def computeSolidInterfaceResidual(self, SolidSolver):
	"""
	Computes the solid interface FSI displacement residual.
	"""

        if self.have_MPI == True:
	  myid = self.comm.Get_rank()
        else:
          myid = 0

	normInterfaceResidualSquare = 0.0

        # --- Create and fill the PETSc vector for the predicted solid interface position (predicted by the solid computation) ---
        if self.have_MPI == True:
          predDisp_array_X = PETSc.Vec().create(self.comm)
          predDisp_array_X.setType('mpi')
          predDisp_array_Y = PETSc.Vec().create(self.comm)
          predDisp_array_Y.setType('mpi')
          predDisp_array_Z = PETSc.Vec().create(self.comm)
          predDisp_array_Z.setType('mpi')
        else:
          predDisp_array_X = PETSc.Vec().create()
          predDisp_array_X.setType('seq')
          predDisp_array_Y = PETSc.Vec().create()
          predDisp_array_Y.setType('seq')
          predDisp_array_Z = PETSc.Vec().create()
          predDisp_array_Z.setType('seq') 
        predDisp_array_X.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        predDisp_array_Y.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        predDisp_array_Z.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
	
	if myid in self.solidSolverProcessors:
	  for iVertex in range(self.nLocalSolidInterfaceNodes):
	    predDispx = SolidSolver.getInterfaceNodeDispX(self.solidInterfaceIdentifier, iVertex)
	    predDispy = SolidSolver.getInterfaceNodeDispY(self.solidInterfaceIdentifier, iVertex)
	    predDispz = SolidSolver.getInterfaceNodeDispZ(self.solidInterfaceIdentifier, iVertex)
            iGlobalVertex = self.__getGlobalIndex('solid', myid, iVertex)
            predDisp_array_X.setValues([iGlobalVertex], predDispx)
            predDisp_array_Y.setValues([iGlobalVertex], predDispy)
            predDisp_array_Z.setValues([iGlobalVertex], predDispz)
	
	predDisp_array_X.assemblyBegin()
	predDisp_array_X.assemblyEnd()
	predDisp_array_Y.assemblyBegin()
	predDisp_array_Y.assemblyEnd()
	predDisp_array_Z.assemblyBegin()
	predDisp_array_Z.assemblyEnd()

        # --- Calculate the residual (vector and norm) ---
        self.solidInterfaceResidual_array_X = predDisp_array_X - self.solidInterface_array_DispX
        self.solidInterfaceResidual_array_Y = predDisp_array_Y - self.solidInterface_array_DispY
        self.solidInterfaceResidual_array_Z = predDisp_array_Z - self.solidInterface_array_DispZ

        normInterfaceResidual_X = self.solidInterfaceResidual_array_X.norm()
        normInterfaceResidual_Y = self.solidInterfaceResidual_array_Y.norm()
        normInterfaceResidual_Z = self.solidInterfaceResidual_array_Z.norm()

        normInterfaceResidualSquare = normInterfaceResidual_X**2 + normInterfaceResidual_Y**2 + normInterfaceResidual_Z**2

        predDisp_array_X.destroy()
        predDisp_array_Y.destroy()
        predDisp_array_Z.destroy()
        del predDisp_array_X
        del predDisp_array_Y
        del predDisp_array_Z

	return sqrt(normInterfaceResidualSquare)

    def relaxSolidPosition(self,FSI_config):
	"""
	Apply solid displacement under-relaxation.
	"""
        if self.have_MPI == True:
	  myid = self.comm.Get_rank()
        else:
          myid = 0

        # --- Set the Aitken coefficient for the relaxation ---
	if FSI_config['AITKEN_RELAX'] == 'STATIC':
	    self.aitkenParam = FSI_config['AITKEN_PARAM']
	elif FSI_config['AITKEN_RELAX'] == 'DYNAMIC':	
	    self.setAitkenCoefficient(FSI_config)
	else:
	    self.aitkenParam = 1.0

	self.MPIPrint('Aitken under-relaxation step with parameter {}'.format(self.aitkenParam))

        # --- Relax the solid interface position ---
        self.solidInterface_array_DispX += self.aitkenParam*self.solidInterfaceResidual_array_X
        self.solidInterface_array_DispY += self.aitkenParam*self.solidInterfaceResidual_array_Y
        self.solidInterface_array_DispZ += self.aitkenParam*self.solidInterfaceResidual_array_Z
	

    def setAitkenCoefficient(self, FSI_config):
	"""
	Computes the Aitken coefficients for solid displacement under-relaxation.
	"""

	deltaResNormSquare = 0.0
	prodScalRes = 0.0
	
        # --- Create the PETSc vector for the difference between the residuals (current and previous FSI iter) ---
	if self.FSIIter == 0:
	    self.aitkenParam = max(FSI_config['AITKEN_PARAM'], self.aitkenParam)
	else:
            if self.have_MPI:
              deltaResx_array_X = PETSc.Vec().create(self.comm)
              deltaResx_array_X.setType('mpi')
              deltaResx_array_Y = PETSc.Vec().create(self.comm)
              deltaResx_array_Y.setType('mpi')
              deltaResx_array_Z = PETSc.Vec().create(self.comm)
              deltaResx_array_Z.setType('mpi')
            else:
              deltaResx_array_X = PETSc.Vec().create()
              deltaResx_array_X.setType('seq')
              deltaResx_array_Y = PETSc.Vec().create()
              deltaResx_array_Y.setType('seq')
              deltaResx_array_Z = PETSc.Vec().create()
              deltaResx_array_Z.setType('seq')
            deltaResx_array_X.setSizes(self.nSolidInterfacePhysicalNodes)
            deltaResx_array_X.set(0.0)
            deltaResx_array_Y.setSizes(self.nSolidInterfacePhysicalNodes)
            deltaResx_array_Y.set(0.0)
            deltaResx_array_Z.setSizes(self.nSolidInterfacePhysicalNodes)
            deltaResx_array_Z.set(0.0)

            # --- Compute the dynamic Aitken coefficient ---
            deltaResx_array_X = self.solidInterfaceResidual_array_X - self.solidInterfaceResidualnM1_array_X
            deltaResx_array_Y = self.solidInterfaceResidual_array_Y - self.solidInterfaceResidualnM1_array_Y
            deltaResx_array_Z = self.solidInterfaceResidual_array_Z - self.solidInterfaceResidualnM1_array_Z

            prodScalRes_X = deltaResx_array_X.dot(self.solidInterfaceResidualnM1_array_X)
            prodScalRes_Y = deltaResx_array_Y.dot(self.solidInterfaceResidualnM1_array_Y)
            prodScalRes_Z = deltaResx_array_Z.dot(self.solidInterfaceResidualnM1_array_Z)
            prodScalRes = prodScalRes_X + prodScalRes_Y + prodScalRes_Z

            deltaResNormSquare_X = (deltaResx_array_X.norm())**2
            deltaResNormSquare_Y = (deltaResx_array_Y.norm())**2
            deltaResNormSquare_Z = (deltaResx_array_Z.norm())**2
	    deltaResNormSquare = deltaResNormSquare_X + deltaResNormSquare_Y + deltaResNormSquare_Z

	    self.aitkenParam *= -prodScalRes/deltaResNormSquare

            deltaResx_array_X.destroy()
            deltaResx_array_Y.destroy()
            deltaResx_array_Z.destroy()
            del deltaResx_array_X
            del deltaResx_array_Y
            del deltaResx_array_Z

        self.aitkenParam = min(self.aitkenParam, 1.0)
        self.aitkenParam = max(self.aitkenParam, 0.0)

        # --- Update the value of the residual for the next FSI iteration ---
        self.solidInterfaceResidual_array_X.copy(self.solidInterfaceResidualnM1_array_X)
        self.solidInterfaceResidual_array_Y.copy(self.solidInterfaceResidualnM1_array_Y)
        self.solidInterfaceResidual_array_Z.copy(self.solidInterfaceResidualnM1_array_Z)

    def displacementPredictor(self, FSI_config , SolidSolver, deltaT):
	"""
	Calculates a prediciton for the solid interface position for the next time step.
	"""

        if self.have_MPI == True:
	  myid = self.comm.Get_rank()
        else:
          myid = 0

	if FSI_config['DISP_PRED'] == 'FIRST_ORDER':
	    self.MPIPrint("First order predictor")	
	    alpha_0 = 1.0
	    alpha_1 = 0.0
	elif FSI_config['DISP_PRED'] == 'SECOND_ORDER':
	    self.MPIPrint("Second order predictor")
	    alpha_0 = 1.0
	    alpha_1 = 0.5
	else:
	    self.MPIPrint("No predictor")
	    alpha_0 = 0.0
	    alpha_1 = 0.0

        # --- Create the PETSc vectors to store the solid interface velocity ---
        if self.have_MPI == True:
          Vel_array_X = PETSc.Vec().create(self.comm)
          Vel_array_X.setType('mpi')
          Vel_array_Y = PETSc.Vec().create(self.comm)
          Vel_array_Y.setType('mpi')
          Vel_array_Z = PETSc.Vec().create(self.comm)
          Vel_array_Z.setType('mpi')
          VelnM1_array_X = PETSc.Vec().create(self.comm)
          VelnM1_array_X.setType('mpi')
          VelnM1_array_Y = PETSc.Vec().create(self.comm)
          VelnM1_array_Y.setType('mpi')
          VelnM1_array_Z = PETSc.Vec().create(self.comm)
          VelnM1_array_Z.setType('mpi')
        else:
          Vel_array_X = PETSc.Vec().create()
          Vel_array_X.setType('seq')
          Vel_array_Y = PETSc.Vec().create()
          Vel_array_Y.setType('seq')
          Vel_array_Z = PETSc.Vec().create()
          Vel_array_Z.setType('seq')
          VelnM1_array_X = PETSc.Vec().create()
          VelnM1_array_X.setType('seq')
          VelnM1_array_Y = PETSc.Vec().create()
          VelnM1_array_Y.setType('seq')
          VelnM1_array_Z = PETSc.Vec().create()
          VelnM1_array_Z.setType('seq')
        Vel_array_X.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        Vel_array_Y.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        Vel_array_Z.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        Vel_array_X.set(0.0)
        Vel_array_Y.set(0.0)
        Vel_array_Z.set(0.0)
        VelnM1_array_X.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        VelnM1_array_Y.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        VelnM1_array_Z.setSizes(self.nSolidInterfacePhysicalNodes+self.d_RBF)
        VelnM1_array_X.set(0.0)
        VelnM1_array_Y.set(0.0)
        VelnM1_array_Z.set(0.0)


        # --- Fill the PETSc vectors ---
        GlobalIndex = int()
        localIndex = 0
	for iVertex in range(self.nLocalSolidInterfaceNodes):
	    GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
            if GlobalIndex in self.SolidHaloNodeList[myid].keys():
              pass
            else:
              iGlobalVertex = self.__getGlobalIndex('solid', myid, localIndex)
	      velx = SolidSolver.getInterfaceNodeVelX(self.solidInterfaceIdentifier, iVertex)
	      vely = SolidSolver.getInterfaceNodeVelY(self.solidInterfaceIdentifier, iVertex)
	      velz = SolidSolver.getInterfaceNodeVelZ(self.solidInterfaceIdentifier, iVertex)
	      velxNm1 = SolidSolver.getInterfaceNodeVelXNm1(self.solidInterfaceIdentifier, iVertex)
	      velyNm1 = SolidSolver.getInterfaceNodeVelYNm1(self.solidInterfaceIdentifier, iVertex)
	      velzNm1 = SolidSolver.getInterfaceNodeVelZNm1(self.solidInterfaceIdentifier, iVertex)
              Vel_array_X.setValues([iGlobalVertex],velx)
              Vel_array_Y.setValues([iGlobalVertex],vely)
              Vel_array_Z.setValues([iGlobalVertex],velz)
              VelnM1_array_X.setValues([iGlobalVertex],velxNm1)
              VelnM1_array_Y.setValues([iGlobalVertex],velyNm1)
              VelnM1_array_Z.setValues([iGlobalVertex],velzNm1)
              localIndex += 1

        Vel_array_X.assemblyBegin()
        Vel_array_X.assemblyEnd()
        Vel_array_Y.assemblyBegin()
        Vel_array_Y.assemblyEnd()
        Vel_array_Z.assemblyBegin()
        Vel_array_Z.assemblyEnd()
        VelnM1_array_X.assemblyBegin()
        VelnM1_array_X.assemblyEnd()
        VelnM1_array_Y.assemblyBegin()
        VelnM1_array_Y.assemblyEnd()
        VelnM1_array_Z.assemblyBegin()
        VelnM1_array_Z.assemblyEnd()

        # --- Predict the solid position for the next time step ---
        self.solidInterface_array_DispX += alpha_0*deltaT*Vel_array_X + alpha_1*deltaT*(Vel_array_X - VelnM1_array_X)
        self.solidInterface_array_DispY += alpha_0*deltaT*Vel_array_Y + alpha_1*deltaT*(Vel_array_Y - VelnM1_array_Y)
        self.solidInterface_array_DispZ += alpha_0*deltaT*Vel_array_Z + alpha_1*deltaT*(Vel_array_Z - VelnM1_array_Z)

        Vel_array_X.destroy()
        Vel_array_Y.destroy()
        Vel_array_Z.destroy()
        VelnM1_array_X.destroy()
        VelnM1_array_Y.destroy()
        VelnM1_array_Z.destroy()
        del Vel_array_X, Vel_array_Y, Vel_array_Z
        del VelnM1_array_X, VelnM1_array_Y, VelnM1_array_Z

    def writeFSIHistory(self, TimeIter, time, varCoordNorm, FSIConv):
	"""
	Write the FSI history file of the computaion.
	"""

        if self.have_MPI == True:
          myid = self.comm.Get_rank()
        else:
          myid = 0
        
        if myid == self.rootProcess:
          if self.unsteady:	
	    if TimeIter == 0:
	      histFile = open('FSIhistory.dat', "w")
              histFile.write("TimeIter\tTime\tFSIRes\tFSINbIter\n")
	    else:
	      histFile = open('FSIhistory.dat', "a")
	    if FSIConv:
	      histFile.write(str(TimeIter) + '\t' + str(time) + '\t' + str(varCoordNorm) + '\t' + str(self.FSIIter+1) + '\n')
	    else:
	      histFile.write(str(TimeIter) + '\t' + str(time) + '\t' + str(varCoordNorm) + '\t' + str(self.FSIIter) + '\n')
	    histFile.close()
          else:
            if self.FSIIter == 0:
              histFile = open('FSIhistory.dat', "w")
              histFile.write("FSI Iter\tFSIRes\n")
            else :
              histFile = open('FSIhistory.dat', "a")
            histFile.write(str(self.FSIIter) + '\t' + str(varCoordNorm) + '\n')
            histFile.close()
          

        self.MPIBarrier()

    def __getGlobalIndex(self, physics, iProc, iLocalVertex):
        """
        Calculate the global indexing of interface nodes accross all the partitions. This does not include halo nodes.
        """

        if physics == 'fluid':
          globalStartIndex = self.fluidGlobalIndexRange[iProc][iProc][0]
        elif physics == 'solid':
          globalStartIndex = self.solidGlobalIndexRange[iProc][iProc][0]

        globalIndex = globalStartIndex + iLocalVertex

        return globalIndex
	    

    def UnsteadyFSI(self,FSI_config, FluidSolver, SolidSolver):
	  """ 
	  Run the unsteady FSI computation by synchronizing the fluid and solid solvers.
	  F/s interface data are exchanged through interface mapping and interpolation (if non mathcing meshes).
	  """

          if self.have_MPI == True:
	    myid = self.comm.Get_rank()
	    numberPart = self.comm.Get_size()
          else:
            myid = 0
            numberPart = 1

	  # --- Set some general variables for the unsteady computation --- #
  	  deltaT = FSI_config['UNST_TIMESTEP']		# physical time step
	  totTime = FSI_config['UNST_TIME']		# physical simulation time
	  NbFSIIterMax = FSI_config['NB_FSI_ITER']	# maximum number of FSI iteration (for each time step)
	  FSITolerance = FSI_config['FSI_TOLERANCE']	# f/s interface tolerance
	  TimeIterTreshold = 0				# time iteration from which we allow the solid to deform

	  if FSI_config['RESTART_SOL'] == 'YES':
	    startTime = FSI_config['START_TIME']
	    NbTimeIter = ((totTime)/deltaT)-1		
	    time = startTime
	    TimeIter = FSI_config['RESTART_ITER']
	  else:
	    NbTimeIter = (totTime/deltaT)-1		# number of time iterations
	    time = 0.0					# initial time
	    TimeIter = 0				# initial time iteration

	  NbTimeIter = int(NbTimeIter)			# be sure that NbTimeIter is an integer

	  varCoordNorm = 0.0				# FSI residual
	  FSIConv = False				# FSI convergence flag

	  self.MPIPrint('\n**********************************')
	  self.MPIPrint('* Begin unsteady FSI computation *')
	  self.MPIPrint('**********************************\n')
	  
	  # --- Initialize the coupled solution --- #
	  #If restart (DOES NOT WORK YET)
	  if FSI_config['RESTART_SOL'] == 'YES':
	    TimeIterTreshold = -1
	    FluidSolver.setTemporalIteration(TimeIter)
	    if myid == self.rootProcess:
	      SolidSolver.outputDisplacements(FluidSolver.getInterRigidDispArray(), True)
            if self.have_MPI == True:
	      self.comm.barrier()
	    FluidSolver.setInitialMesh(True)
	    if myid == self.rootProcess:
	      SolidSolver.displacementPredictor(FluidSolver.getInterRigidDispArray())
	    if self.have_MPI == True:
              self.comm.barrier()								
	    if myid == self.rootProcess:
	      SolidSolver.updateSolution()
	  #If no restart
	  else:
	    self.MPIPrint('Setting FSI initial conditions')
            if myid in self.solidSolverProcessors:
	      SolidSolver.setInitialDisplacements()
            self.getSolidInterfaceDisplacement(SolidSolver)
	    self.interpolateSolidPositionOnFluidMesh(FSI_config)
	    self.setFluidInterfaceVarCoord(FluidSolver)
	    FluidSolver.SetInitialMesh()	# if there is an initial deformation in the solid, it has to be communicated to the fluid solver
	    self.MPIPrint('\nFSI initial conditions are set')
	    self.MPIPrint('Beginning time integration\n')

	  # --- External temporal loop --- #
	  while TimeIter <= NbTimeIter:

		if TimeIter > TimeIterTreshold:
		  NbFSIIter = NbFSIIterMax
		  self.MPIPrint('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI on time iteration {} ***************'.format(TimeIter))  
		else:
		  NbFSIIter = 1

		self.FSIIter = 0
		FSIConv = False
		FluidSolver.PreprocessExtIter(TimeIter)	# set some parameters before temporal fluid iteration
	
		# --- Internal FSI loop --- #
		while self.FSIIter <= (NbFSIIter-1):

			self.MPIPrint("\n>>>> Time iteration {} / FSI iteration {} <<<<".format(TimeIter,self.FSIIter))

			# --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
			self.interpolateSolidPositionOnFluidMesh(FSI_config)
                        self.MPIPrint('\nPerforming dynamic mesh deformation (ALE)...\n')
                        self.setFluidInterfaceVarCoord(FluidSolver)
                        FluidSolver.DynamicMeshUpdate(TimeIter)
			
			# --- Fluid solver call for FSI subiteration --- #
		        self.MPIPrint('\nLaunching fluid solver for one single dual-time iteration...')
                        self.MPIBarrier()
			FluidSolver.ResetConvergence()
			FluidSolver.Run()
                        self.MPIBarrier()

			# --- Surface fluid loads interpolation and communication --- #
		        self.MPIPrint('\nProcessing interface fluid loads...\n')
                        self.MPIBarrier()
		        self.getFluidInterfaceNodalForce(FSI_config, FluidSolver)
                        self.MPIBarrier()
			if TimeIter > TimeIterTreshold:
			  self.interpolateFluidLoadsOnSolidMesh(FSI_config)
			  self.setSolidInterfaceLoads(SolidSolver, FSI_config, time)

			  # --- Solid solver call for FSI subiteration --- #
			  self.MPIPrint('\nLaunching solid solver for a single time iteration...\n')
                          if myid in self.solidSolverProcessors:
			    if FSI_config['CSD_SOLVER'] == 'NATIVE':
			        SolidSolver.timeIteration(time)
			    elif FSI_config['CSD_SOLVER'] == 'METAFOR' or FSI_config['CSD_SOLVER'] == 'GETDP' or FSI_config['CSD_SOLVER'] == 'TESTER':
				SolidSolver.run(time-deltaT, time)

			  # --- Compute and monitor the FSI residual --- #
			  varCoordNorm = self.computeSolidInterfaceResidual(SolidSolver)
			  self.MPIPrint('\nFSI displacement norm : {}\n'.format(varCoordNorm))
			  if varCoordNorm < FSITolerance:		
			    FSIConv = True
			    break

			  # --- Relaxe the solid position --- #
                          self.MPIPrint('\nProcessing interface displacements...\n')
			  self.relaxSolidPosition(FSI_config)
			
			self.FSIIter += 1
		# --- End OF FSI loop --- #

                self.MPIBarrier()

		# --- Update the FSI history file --- # 
		if TimeIter > TimeIterTreshold:
		  self.MPIPrint('\nBGS is converged (strong coupling)')
		self.writeFSIHistory(TimeIter, time, varCoordNorm, FSIConv)
		
		# --- Update, monitor and output the fluid solution before the next time step  ---#
		FluidSolver.Update()
		FluidSolver.Monitor(TimeIter)
		FluidSolver.Output(TimeIter)

	   	if TimeIter >= TimeIterTreshold:
		  if myid in self.solidSolverProcessors:
		    # --- Output the solid solution before thr next time step --- #
		    SolidSolver.writeSolution(time, self.FSIIter, TimeIter, NbTimeIter)
		
		  # --- Displacement predictor for the next time step and update of the solid solution --- #
		  self.MPIPrint('\nSolid displacement prediction for next time step')
		  self.displacementPredictor(FSI_config, SolidSolver, deltaT)
                  if myid in self.solidSolverProcessors:
		    SolidSolver.updateSolution()
		
		TimeIter += 1
		time += deltaT
	  #--- End of the temporal loop --- #

          self.MPIBarrier()

	  self.MPIPrint('\n*************************')
	  self.MPIPrint('*  End FSI computation  *')
	  self.MPIPrint('*************************\n')

    def SteadyFSI(self, FSI_config,FluidSolver, SolidSolver):
	  """
	  Runs the steady FSI computation by synchronizing the fluid and solid solver with data exchange at the f/s interface.
	  """

          if self.have_MPI == True:
	    myid = self.comm.Get_rank()
	    numberPart = self.comm.Get_size()
          else:
            myid = 0
            numberPart = 1

	  # --- Set some general variables for the steady computation --- #
	  NbIter = FSI_config['NB_EXT_ITER']		# number of fluid iteration at each FSI step
	  NbFSIIterMax = FSI_config['NB_FSI_ITER']	# maximum number of FSI iteration (for each time step)
	  FSITolerance = FSI_config['FSI_TOLERANCE']	# f/s interface tolerance
	  varCoordNorm = 0.0

	  self.MPIPrint('\n********************************')
	  self.MPIPrint('* Begin steady FSI computation *')
	  self.MPIPrint('********************************\n')
	  self.MPIPrint('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI ***************')

          self.getSolidInterfaceDisplacement(SolidSolver)

	  # --- External FSI loop --- #
	  self.FSIIter = 0
	  while self.FSIIter < NbFSIIterMax:
	    self.MPIPrint("\n>>>> FSI iteration {} <<<<".format(self.FSIIter))
	    self.MPIPrint('\nLaunching fluid solver for a steady computation...')
	    # --- Fluid solver call for FSI subiteration ---#
	    Iter = 0
	    FluidSolver.ResetConvergence()          
	    while Iter < NbIter:
	      FluidSolver.PreprocessExtIter(Iter)
	      FluidSolver.Run()				
	      StopIntegration = FluidSolver.Monitor(Iter)
	      FluidSolver.Output(Iter)
	      if StopIntegration:
		break;
	      Iter += 1
	    
	    # --- Surface fluid loads interpolation and communication ---#
	    self.MPIPrint('\nProcessing interface fluid loads...\n')
            self.MPIBarrier()
	    self.getFluidInterfaceNodalForce(FSI_config, FluidSolver)
            self.MPIBarrier()
	    self.interpolateFluidLoadsOnSolidMesh(FSI_config)
	    self.setSolidInterfaceLoads(SolidSolver, FSI_config, 0.05)
	     
	    # --- Solid solver call for FSI subiteration --- #
	    self.MPIPrint('\nLaunching solid solver for a static computation...\n')
            if myid in self.solidSolverProcessors:
	      if FSI_config['CSD_SOLVER'] == 'NATIVE':
	          SolidSolver.staticComputation()
              else:
                  SolidSolver.run(0.0, 0.05)
	      SolidSolver.writeSolution(0.0, self.FSIIter, Iter, NbIter)		

	    # --- Compute and monitor the FSI residual --- #
	    varCoordNorm = self.computeSolidInterfaceResidual(SolidSolver)
	    self.MPIPrint('\nFSI displacement norm : {}\n'.format(varCoordNorm))
            self.writeFSIHistory(0, 0.0, varCoordNorm, False)
	    if varCoordNorm < FSITolerance:			
	      break

            # --- Relaxe the solid displacement and update the solid solution --- #
            self.MPIPrint('\nProcessing interface displacements...\n')
	    self.relaxSolidPosition(FSI_config)
            if myid in self.solidSolverProcessors:
              SolidSolver.updateSolution()
	
	    # --- Mesh morphing step (displacement interpolation, displacements communication, and mesh morpher call) --- #
	    self.interpolateSolidPositionOnFluidMesh(FSI_config)						
	    self.MPIPrint('\nPerforming static mesh deformation...\n')
	    self.setFluidInterfaceVarCoord(FluidSolver)									
	    FluidSolver.StaticMeshUpdate()
	    self.FSIIter += 1

          self.MPIBarrier()

	  self.MPIPrint('\nBGS is converged (strong coupling)')
	  self.MPIPrint(' ')
	  self.MPIPrint('*************************')
	  self.MPIPrint('*  End FSI computation  *')
	  self.MPIPrint('*************************')
	  self.MPIPrint(' ')
