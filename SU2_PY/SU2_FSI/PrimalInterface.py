#!/usr/bin/env python

## \file Interface.py
#  \brief Interface class that handles fluid/solid solvers synchronisation and communication.
#  \author Rocco Bombardieri and Ruben Sanchez based on previous work by David Thomas.
#  \version 7.0.0
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
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


import numpy as np
import shutil
import os
# ----------------------------------------------------------------------
#  FSI Interface Class
# ----------------------------------------------------------------------


class Interface:
    """ 
    FSI interface class that handles fluid/solid solvers synchronisation and communication
    """

    def __init__(self, FSI_config, FluidSolver, SolidSolver, MLS_Solver, have_MPI):
        """
        Class constructor. Declare some variables and do some screen outputs.
        """

        if have_MPI:
            from mpi4py import MPI
            self.MPI = MPI
            self.comm = MPI.COMM_WORLD  # MPI World communicator
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

        self.localFluidInterface_vertex_indices = None

        self.globalFluidInterfaceXcoor = None
        self.globalFluidInterfaceYcoor = None
        self.globalFluidInterfaceZcoor = None

        self.globalFluidCoordinates = None
        self.globalSolidCoordinates = None

        self.sendCounts = None
        self.globalFluidDispX = None
        self.globalFluidDispY = None
        self.globalFluidDispZ = None

        self.globalSolidDispX = None
        self.globalSolidDispY = None
        self.globalSolidDispZ = None

        self.globalSolidDispXOld = None
        self.globalSolidDispYOld = None
        self.globalSolidDispZOld = None

        self.globalFluidLoadX = None
        self.globalFluidLoadY = None
        self.globalFluidLoadZ = None

        self.globalSolidLoadX = None
        self.globalSolidLoadY = None
        self.globalSolidLoadZ = None

        self.haloNodesPositionsInit = {}  # initial position of the halo nodes (fluid side only)

        self.solidInterface_array_DispX = None  # solid interface displacement
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

        self.FSIIter = 0  # current FSI iteration
        self.unsteady = False  # flag for steady or unsteady simulation (default is steady)

        # ---Some screen output ---
        self.MPIPrint('Fluid solver : SU2_CFD')
        self.MPIPrint('Solid solver : pyBeam')
        self.MPIPrint('Steady coupled simulation')
        self.MPIPrint('Matching fluid-solid interface using Moving Least Squares method')
        self.MPIPrint('Maximum number of FSI iterations : {}'.format(FSI_config['NB_FSI_ITER']))
        self.MPIPrint('FSI tolerance : {}'.format(FSI_config['FSI_TOLERANCE']))
        self.MPIPrint('Static under-relaxation with constant parameter {}'.format(FSI_config['RELAX_PARAM']))
        self.MPIPrint('FSI interface is set')

    def MPIPrint(self, message):
        """
        Print a message on screen only from the master process.
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        if myid == self.rootProcess:
            print(message)

    def MPIBarrier(self):
        """
        Perform a synchronization barrier in case of parallel run with MPI.
        """

        if self.have_MPI:
            self.comm.barrier()

    def checkMPI(self):
        """
        Return the MPI characteristics of the problem
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

        return myid, MPIsize

    def connect(self, FSI_config, FluidSolver, SolidSolver):
        """
        Connection between solvers.
        Creates the communication support between the two solvers.
        Gets information about f/s interfaces from the two solvers.
        """

        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        # --- Identify the fluid interface and store the number of nodes for each partition ---#
        self.fluidInterfaceIdentifier = None
        self.nLocalFluidInterfaceNodes = 0
        if FluidSolver != None:
            print('Fluid solver is initialized on process {}'.format(myid))
            self.haveFluidSolver = True
            allInterfaceMarkersTags = FluidSolver.GetAllDeformMeshMarkersTag()
            allMarkersID = FluidSolver.GetAllBoundaryMarkers()
            if not allInterfaceMarkersTags:
                raise Exception('No moving marker was defined in SU2.')
            else:
                if allInterfaceMarkersTags[0] in allMarkersID.keys():
                    self.fluidInterfaceIdentifier = allMarkersID[allInterfaceMarkersTags[0]]
            if self.fluidInterfaceIdentifier != None:
                self.nLocalFluidInterfaceNodes = FluidSolver.GetNumberVertices(self.fluidInterfaceIdentifier)
            if self.nLocalFluidInterfaceNodes != 0:
                self.haveFluidInterface = True
                print('Number of interface fluid nodes (halo nodes included) on proccess {} and marker {}: {}'\
                      .format(myid,allInterfaceMarkersTags[0],self.nLocalFluidInterfaceNodes))
        else:
            pass

        # --- Identify the solid interface and store the number of nodes (single core) ---#
        if SolidSolver != None:
            print('Solid solver is initialized on process {}'.format(myid))
            self.haveSolidSolver = True
            self.nSolidInterfaceNodes = SolidSolver.nPoint
            self.nSolidInterfacePhysicalNodes = SolidSolver.nPoint
            self.nLocalSolidInterfaceNodes = SolidSolver.nPoint
            self.globalSolidCoordinates = np.zeros((SolidSolver.nPoint, 3))
            for iPoint in range(0, SolidSolver.nPoint):
                coordX, coordY, coordZ = SolidSolver.GetInitialCoordinates(iPoint)
                self.globalSolidCoordinates[iPoint, 0] = coordX
                self.globalSolidCoordinates[iPoint, 1] = coordY
                self.globalSolidCoordinates[iPoint, 2] = coordZ
        else:
            pass

        # --- Exchange information about processors --- #
        if self.have_MPI:
            if self.haveFluidSolver:
                sendBufFluid = np.array(int(1))
            else:
                sendBufFluid = np.array(int(0))
            if self.haveFluidInterface:
                sendBufFluidInterface = np.array(int(1))
            else:
                sendBufFluidInterface = np.array(int(0))
            rcvBufFluid = np.zeros(MPIsize, dtype=int)
            rcvBufFluidInterface = np.zeros(MPIsize, dtype=int)
            self.comm.Allgather(sendBufFluid, rcvBufFluid)
            self.comm.Allgather(sendBufFluidInterface, rcvBufFluidInterface)

            for iProc in range(MPIsize):
                if rcvBufFluid[iProc] == 1:
                    self.fluidSolverProcessors.append(iProc)
                if rcvBufFluidInterface[iProc] == 1:
                    self.fluidInterfaceProcessors.append(iProc)

            del sendBufFluid, rcvBufFluid, sendBufFluidInterface, rcvBufFluidInterface
        else:
            self.fluidSolverProcessors.append(0)
            self.fluidInterfaceProcessors.append(0)

        self.MPIBarrier()

        # --- Calculate the total number of nodes at the fluid interface (sum over all the partitions)
        # --- Calculate the number of halo nodes on each partition
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

        # --- Calculate the total number of nodes (with and without halo) at the fluid interface (sum over all the partitions) and broadcast the number accross all processors ---
        sendBuffHalo = np.array(int(self.nLocalFluidInterfaceNodes))
        sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
        rcvBuffHalo = np.zeros(1, dtype=int)
        rcvBuffPhysical = np.zeros(1, dtype=int)
        if self.have_MPI:
            self.comm.barrier()
            self.comm.Allreduce(sendBuffHalo, rcvBuffHalo, op=self.MPI.SUM)
            self.comm.Allreduce(sendBuffPhysical, rcvBuffPhysical, op=self.MPI.SUM)
            self.nFluidInterfaceNodes = rcvBuffHalo[0]
            self.nFluidInterfacePhysicalNodes = rcvBuffPhysical[0]
        else:
            self.nFluidInterfaceNodes = np.copy(sendBuffHalo)
            self.nFluidInterfacePhysicalNodes = np.copy(sendBuffPhysical)
        del sendBuffHalo, rcvBuffHalo, sendBuffPhysical, rcvBuffPhysical

        # --- Store the number of physical interface nodes on each processor and allgather the information ---
        self.fluidPhysicalInterfaceNodesDistribution = np.zeros(MPIsize, dtype=int)
        if self.have_MPI:
            sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
            self.comm.Allgather(sendBuffPhysical, self.fluidPhysicalInterfaceNodesDistribution)
            del sendBuffPhysical
        else:
            self.fluidPhysicalInterfaceNodesDistribution[0] = self.nFluidInterfacePhysicalNodes

        # --- Calculate and store the global indexing of interface physical nodes on each processor and allgather the information ---
        if self.have_MPI:
            if myid in self.fluidInterfaceProcessors:
                globalIndexStart = 0
                for iProc in range(myid):
                    globalIndexStart += self.fluidPhysicalInterfaceNodesDistribution[iProc]
                globalIndexStop = globalIndexStart + self.nLocalFluidInterfacePhysicalNodes - 1
            else:
                globalIndexStart = 0
                globalIndexStop = 0
            self.fluidGlobalIndexRange[myid] = [globalIndexStart, globalIndexStop]
            self.fluidGlobalIndexRange = self.comm.allgather(self.fluidGlobalIndexRange)
        else:
            temp = {}
            temp[0] = [0, self.nLocalFluidInterfacePhysicalNodes - 1]
            self.fluidGlobalIndexRange = list()
            self.fluidGlobalIndexRange.append(temp)


        self.MPIPrint(
            'Total number of fluid interface nodes (halo nodes included) : {}'.format(self.nFluidInterfaceNodes))
        self.MPIPrint(
            'Total number of physical fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))
        self.MPIPrint(
            'Total number of beam interface nodes : {}'.format(self.nSolidInterfaceNodes))
        self.MPIPrint('Total number of fluid interface nodes : {}'.format(self.nFluidInterfacePhysicalNodes))
        self.MPIPrint('Total number of solid interface nodes : {}'.format(self.nSolidInterfacePhysicalNodes))

        self.MPIBarrier()

        # --- Get, for the fluid interface on each partition:
        # --- The vertex indices, which are stored locally on each processor
        # --- The coordinates X, Y and Z, which are stored in a global list in processor 0
        GlobalIndex = int()
        localIndex = 0
        fluidIndexing_temp = {}
        self.localFluidInterface_vertex_indices = np.zeros(self.nLocalFluidInterfacePhysicalNodes, dtype=int)
        localFluidInterface_array_X_init = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidInterface_array_Y_init = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidInterface_array_Z_init = np.zeros(self.nLocalFluidInterfacePhysicalNodes)

        for iVertex in range(self.nLocalFluidInterfaceNodes):
            GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
            posx, posy, posz = FluidSolver.GetVertex_UndeformedCoord(self.fluidInterfaceIdentifier, iVertex)

            if GlobalIndex in self.FluidHaloNodeList[myid].keys():
                self.haloNodesPositionsInit[GlobalIndex] = (posx, posy, posz)
            else:
                fluidIndexing_temp[GlobalIndex] = self.__getGlobalIndex('fluid', myid, localIndex)
                localFluidInterface_array_X_init[localIndex] = posx
                localFluidInterface_array_Y_init[localIndex] = posy
                localFluidInterface_array_Z_init[localIndex] = posz
                self.localFluidInterface_vertex_indices[localIndex] = int(iVertex)
                localIndex += 1

        #print("rank: {}, local_vertex_indices: {}".format(myid, self.localFluidInterface_vertex_indices))

        if self.have_MPI:
            fluidIndexing_temp = self.comm.allgather(fluidIndexing_temp)
            for ii in range(len(fluidIndexing_temp)):
                for key, value in fluidIndexing_temp[ii].items():
                    self.fluidIndexing[key] = value

             # --- Collect local array sizes using gather in root process
            bufXCoor = np.array(localFluidInterface_array_X_init)
            bufYCoor = np.array(localFluidInterface_array_Y_init)
            bufZCoor = np.array(localFluidInterface_array_Z_init)
            self.sendCounts = np.array(self.comm.gather(self.nLocalFluidInterfacePhysicalNodes, 0))

            if myid == self.rootProcess:
                print("sendCounts: {}, total: {}".format(self.sendCounts, sum(self.sendCounts)))
                self.globalFluidInterfaceXcoor = np.empty(sum(self.sendCounts))
                self.globalFluidInterfaceYcoor = np.empty(sum(self.sendCounts))
                self.globalFluidInterfaceZcoor = np.empty(sum(self.sendCounts))

            self.comm.Gatherv(sendbuf=bufXCoor, recvbuf=(self.globalFluidInterfaceXcoor, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufYCoor, recvbuf=(self.globalFluidInterfaceYcoor, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufZCoor, recvbuf=(self.globalFluidInterfaceZcoor, self.sendCounts), root=0)
            #if myid == 0:
                #print("Gathered array X: {}".format(self.globalFluidInterfaceXcoor))
                #print("Gathered array Y: {}".format(self.globalFluidInterfaceYcoor))
                #print("Gathered array Z: {}".format(self.globalFluidInterfaceZcoor))

        else:
            self.fluidIndexing = fluidIndexing_temp.copy()
            self.globalFluidInterfaceXcoor = localFluidInterface_array_X_init.copy()
            self.globalFluidInterfaceYcoor = localFluidInterface_array_Y_init.copy()
            self.globalFluidInterfaceZcoor = localFluidInterface_array_Z_init.copy()

        # Store the global fluid coordinates
        if myid == self.rootProcess:
            self.globalFluidCoordinates = np.zeros((self.nFluidInterfacePhysicalNodes, 3))
            for i in range(0, self.nFluidInterfacePhysicalNodes):
                self.globalFluidCoordinates[i][0] = self.globalFluidInterfaceXcoor[i]
                self.globalFluidCoordinates[i][1] = self.globalFluidInterfaceYcoor[i]
                self.globalFluidCoordinates[i][2] = self.globalFluidInterfaceZcoor[i]

            print(self.globalFluidCoordinates.shape)
            print(self.globalSolidCoordinates.shape)

        del fluidIndexing_temp, localFluidInterface_array_X_init, \
            localFluidInterface_array_Y_init, localFluidInterface_array_Z_init

        ###################################################################################################
        # Initialize the local load array
        # The initial displacements are set to 0 (this might be changed when we have restart capabilities)
        ###################################################################################################

        self.globalFluidLoadX = np.zeros(self.nFluidInterfacePhysicalNodes)
        self.globalFluidLoadY = np.zeros(self.nFluidInterfacePhysicalNodes)
        self.globalFluidLoadZ = np.zeros(self.nFluidInterfacePhysicalNodes)

        if self.haveSolidSolver:

            self.globalSolidLoadX = np.zeros(self.nSolidInterfaceNodes)
            self.globalSolidLoadY = np.zeros(self.nSolidInterfaceNodes)
            self.globalSolidLoadZ = np.zeros(self.nSolidInterfaceNodes)

            self.globalSolidDispX = np.zeros(self.nSolidInterfaceNodes)
            self.globalSolidDispY = np.zeros(self.nSolidInterfaceNodes)
            self.globalSolidDispZ = np.zeros(self.nSolidInterfaceNodes)

    def transferFluidTractions(self, FluidSolver, SolidSolver, MLSSolver):
        """
        Transfer fluid tractions.
        Gathers the fluid tractions from the interface into the root process.
        Interpolates the tractions using the transposed matrix.
        Applies the tractions into the solid solver.
        """

        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        ################################################################################################################
        # --- STEP 1: Retrieve the fluid loads
        # --- Get, for the fluid interface on each partition:
        # --- The vertex indices, which are stored locally on each processor
        # --- The coordinates X, Y and Z, which are stored in a global list in processor 0
        ################################################################################################################

        # Initialize the local load array
        localFluidLoadX = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidLoadY = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidLoadZ = np.zeros(self.nLocalFluidInterfacePhysicalNodes)

        localIndex = 0
        # For the vertices that belong to the interface
        for iVertex in self.localFluidInterface_vertex_indices:
            # Compute the vertex forces on the fluid solver
            loadX, loadY, loadZ = FluidSolver.GetFlowLoad(self.fluidInterfaceIdentifier, int(iVertex))
            # Store them in the local load array
            localFluidLoadX[localIndex] = loadX
            localFluidLoadY[localIndex] = loadY
            localFluidLoadZ[localIndex] = loadZ
            localIndex += 1

        if self.have_MPI:

            # Store the local loads in buffers in the form of numpy arrays
            bufXLoad = np.array(localFluidLoadX)
            bufYLoad = np.array(localFluidLoadY)
            bufZLoad = np.array(localFluidLoadZ)

            # Initialize the global load array
            if myid == self.rootProcess:
                print("sendCounts: {}, total: {}".format(self.sendCounts, sum(self.sendCounts)))
                self.globalFluidLoadX = np.empty(sum(self.sendCounts))
                self.globalFluidLoadY = np.empty(sum(self.sendCounts))
                self.globalFluidLoadZ = np.empty(sum(self.sendCounts))

            # Gatherv using self.sendCounts maintains the ordering of the coordinates
            self.comm.Gatherv(sendbuf=bufXLoad, recvbuf=(self.globalFluidLoadX, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufYLoad, recvbuf=(self.globalFluidLoadY, self.sendCounts), root=0)
            self.comm.Gatherv(sendbuf=bufZLoad, recvbuf=(self.globalFluidLoadZ, self.sendCounts), root=0)

        else:
            self.globalFluidLoadX = localFluidLoadX.copy()
            self.globalFluidLoadY = localFluidLoadY.copy()
            self.globalFluidLoadZ = localFluidLoadZ.copy()

        # Delete local variables
        del localFluidLoadX, localFluidLoadY, localFluidLoadZ

        ################################################################################################################
        # --- STEP 2: Interpolate
        ################################################################################################################

        # ---> Input: self.globalFluidLoadX, self.globalFluidLoadY, self.globalFluidLoadZ

        if myid == self.rootProcess:
            self.globalSolidLoadX = MLSSolver.interpolation_matrix.transpose().dot(self.globalFluidLoadX)
            self.globalSolidLoadY = MLSSolver.interpolation_matrix.transpose().dot(self.globalFluidLoadY)
            self.globalSolidLoadZ = MLSSolver.interpolation_matrix.transpose().dot(self.globalFluidLoadZ)

        # ---> Output: self.globalSolidLoadX, self.globalSolidLoadY, self.globalSolidLoadZ

        ################################################################################################################
        # --- STEP 3: Check conservation
        ################################################################################################################

        # --- Check for total force conservation after interpolation
        if myid == self.rootProcess:

            # Total loads before interpolation, fluid side
            FFX = self.globalFluidLoadX.sum()
            FFY = self.globalFluidLoadY.sum()
            FFZ = self.globalFluidLoadZ.sum()

            # Total loads after interpolation, solid side
            FX = self.globalSolidLoadX.sum()
            FY = self.globalSolidLoadY.sum()
            FZ = self.globalSolidLoadZ.sum()

            self.MPIPrint("Checking f/s interface total force...")
            self.MPIPrint('Solid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FX, FY, FZ))
            self.MPIPrint('Fluid side (Fx, Fy, Fz) = ({}, {}, {})'.format(FFX, FFY, FFZ))

            hist_file = open("historyFSI.dat", "a")
            hist_file.write(str(FFX) + "\t" + str(FFY) + "\t" + str(FFZ) + "\t" + str(FX) + "\t" + str(FY) + "\t" + str(FZ) + "\t")
            hist_file.close()

        ################################################################################################################
        # --- STEP 4: Transfer to the structural solver
        # --- pyBeam runs in single core, so there is no need to deal with the parallelization
        ################################################################################################################

        if self.haveSolidSolver:
            # For the vertices that belong to the interface
            for iVertex in range(0, self.nSolidInterfaceNodes):
                # Store them in the solid solver directly
                SolidSolver.SetLoads(iVertex, self.globalSolidLoadX[iVertex],
                                              self.globalSolidLoadY[iVertex],
                                              self.globalSolidLoadZ[iVertex])

    def transferStructuralDisplacements(self, FSIConfig, FluidSolver, SolidSolver, MLSSolver):
        """
        Transfer structural displacements.
        Gathers the structural displacements from the interface.
        Interpolates the displacements using the transposed matrix.
        Applies the fluid displacements by scattering them into the correct position for the fluid solver.
        """

        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        ################################################################################################################
        # --- STEP 1: Retrieve the structural displacements from pyBeam and apply the relaxation
        # --- pyBeam runs in single core, so there is no need to deal with the parallelization
        ################################################################################################################

        if self.haveSolidSolver:

            # Recover the relaxation parameter
            relaxParam = float(FSIConfig["RELAX_PARAM"])

            # Store the old displacements
            self.globalSolidDispXOld = self.globalSolidDispX.copy()
            self.globalSolidDispYOld = self.globalSolidDispY.copy()
            self.globalSolidDispZOld = self.globalSolidDispZ.copy()

            # Initialize local vectors, to store the relaxed displacements
            relaxedSolidDispX = np.zeros(self.nSolidInterfaceNodes)
            relaxedSolidDispY = np.zeros(self.nSolidInterfaceNodes)
            relaxedSolidDispZ = np.zeros(self.nSolidInterfaceNodes)

            # For the vertices that belong to the interface
            for iVertex in range(0, self.nSolidInterfaceNodes):

                # Store the new displacements in the global load array directly
                dispX, dispY, dispZ = SolidSolver.ExtractDisplacements(iVertex)
                self.globalSolidDispX[iVertex] = dispX
                self.globalSolidDispY[iVertex] = dispY
                self.globalSolidDispZ[iVertex] = dispZ

                # Compute the relaxed, interface displacements
                relaxedSolidDispX[iVertex] = relaxParam * dispX + (1.0 - relaxParam) * self.globalSolidDispXOld[iVertex]
                relaxedSolidDispY[iVertex] = relaxParam * dispY + (1.0 - relaxParam) * self.globalSolidDispYOld[iVertex]
                relaxedSolidDispZ[iVertex] = relaxParam * dispZ + (1.0 - relaxParam) * self.globalSolidDispZOld[iVertex]

        ################################################################################################################
        # --- STEP 2: Interpolate
        ################################################################################################################

        # ---> Input: relaxedSolidDispX, relaxedSolidDispY, relaxedSolidDispZ

        if myid == self.rootProcess:

            self.globalFluidDispX = MLSSolver.interpolation_matrix.dot(relaxedSolidDispX)
            self.globalFluidDispY = MLSSolver.interpolation_matrix.dot(relaxedSolidDispY)
            self.globalFluidDispZ = MLSSolver.interpolation_matrix.dot(relaxedSolidDispZ)

        # ---> Output: self.globalFluidDispX, self.globalFluidDispY, self.globalFluidDispZ

        ################################################################################################################
        # --- STEP 3: Check conservation
        ################################################################################################################

        # --- Checking conservation ---

        if myid == self.rootProcess:

            WSX = self.globalSolidLoadX.dot(relaxedSolidDispX)
            WSY = self.globalSolidLoadY.dot(relaxedSolidDispY)
            WSZ = self.globalSolidLoadZ.dot(relaxedSolidDispZ)

            WFX = self.globalFluidLoadX.dot(self.globalFluidDispX)
            WFY = self.globalFluidLoadY.dot(self.globalFluidDispY)
            WFZ = self.globalFluidLoadZ.dot(self.globalFluidDispZ)

            self.MPIPrint("Checking f/s interface conservation...")
            self.MPIPrint('Solid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WSX, WSY, WSZ))
            self.MPIPrint('Fluid side (Wx, Wy, Wz) = ({}, {}, {})'.format(WFX, WFY, WFZ))

        ################################################################################################################
        # --- STEP 4: Transfer to the fluid solver
        ################################################################################################################

        # --- Recover them from the interpolated vectors
        localFluidDispX = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidDispY = np.zeros(self.nLocalFluidInterfacePhysicalNodes)
        localFluidDispZ = np.zeros(self.nLocalFluidInterfacePhysicalNodes)

        if self.have_MPI:

            self.comm.Scatterv(sendbuf=(self.globalFluidDispX, self.sendCounts), recvbuf=localFluidDispX, root=0)
            self.comm.Scatterv(sendbuf=(self.globalFluidDispY, self.sendCounts), recvbuf=localFluidDispY, root=0)
            self.comm.Scatterv(sendbuf=(self.globalFluidDispZ, self.sendCounts), recvbuf=localFluidDispZ, root=0)

            # print("rank: {}, local_array X: {}".format(myid, localFluidDispX))
            # print("rank: {}, local_array Y: {}".format(myid, localFluidDispY))
            # print("rank: {}, local_array Z: {}".format(myid, localFluidDispZ))

        else:
            localFluidDispX = self.globalFluidDispX.copy()
            localFluidDispY = self.globalFluidDispY.copy()
            localFluidDispZ = self.globalFluidDispZ.copy()

        # For the vertices that belong to the interface
        localIndex = 0
        for iVertex in self.localFluidInterface_vertex_indices:
            # Store them in the mesh displacement routine
            FluidSolver.SetMeshDisplacement(self.fluidInterfaceIdentifier, int(iVertex), localFluidDispX[localIndex],
                                            localFluidDispY[localIndex], localFluidDispZ[localIndex])
            # Increment the local index
            localIndex += 1

        # Delete local variables
        del localFluidDispX, localFluidDispY, localFluidDispZ
        if myid == self.rootProcess:
            del relaxedSolidDispX, relaxedSolidDispY, relaxedSolidDispZ

    def SteadyFSI(self, FSIconfig, FluidSolver, SolidSolver, MLSSolver):
        """
        Runs the steady FSI computation
        Synchronizes the fluid and solid solver with data exchange at the f/s interface.
        """

        # Recover the process and the size of the parallelization
        myid, MPIsize = self.checkMPI()

        # --- Set some general variables for the steady computation --- #
        nFSIIter = FSIconfig['NB_FSI_ITER']  # maximum number of FSI iteration (for each time step)

        if myid is 0:
            hist_file = open("historyFSI.dat", "w")
            hist_file.write("FF(X) \t FF(Y) \t FF(Z) \t FS(X) \t FS(Y) \t FS(Z) \t CD \t CL \n")
            hist_file.close()


        self.MPIPrint('\n********************************')
        self.MPIPrint('* Begin steady FSI computation *')
        self.MPIPrint('********************************\n')

        self.MPIPrint('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI ***************')

        # --- External FSI loop --- #
        self.FSIIter = 0

        # For the number of iterations allowed
        while self.FSIIter < nFSIIter:

            self.MPIPrint("\n>>>> FSI iteration {} <<<<".format(self.FSIIter))

            if self.FSIIter > 0:
                # --- Surface displacements interpolation and communication ---#
                self.MPIPrint('\n##### Transferring displacements\n')
                self.MPIBarrier()
                self.transferStructuralDisplacements(FSIconfig, FluidSolver, SolidSolver, MLSSolver)

            # --- Fluid solver call for FSI subiteration --- #
            self.MPIPrint('\n##### Launching fluid solver for a steady computation\n')
            self.MPIBarrier()
            FluidSolver.ResetConvergence()     # Make sure the solver starts convergence from 0
            FluidSolver.Preprocess(0)          # Time iteration pre-processing
            FluidSolver.Run()                  # Run one time-step (static: one simulation)
            FluidSolver.Postprocess()          # Run one time-step (static: one simulation)
            FluidSolver.Update()               # Update the solver for the next time iteration
            FluidSolver.Monitor(0)             # Monitor the solver and output solution to file if required
            FluidSolver.Output(0)              # Output the solution to file

            # --- Surface fluid loads interpolation and communication ---#
            self.MPIPrint('\n##### Transferring fluid tractions to the beam solver\n')
            self.MPIBarrier()
            self.transferFluidTractions(FluidSolver, SolidSolver, MLSSolver)

            # --- Solid solver call for FSI subiteration --- #
            self.MPIPrint('\n##### Launching solid solver for a static computation\n')
            self.MPIBarrier()
            if self.haveSolidSolver:
                SolidSolver.run()

            self.FSIIter += 1

            # Store the surface flow history
            if myid is 0:
                new_name_surf = "./surface_flow_" + str("{:04d}".format(self.FSIIter)) + ".vtk"
                shutil.move("surface_flow.vtk", new_name_surf)

                hist_file = open("historyFSI.dat", "a")
                hist_file.write(str(FluidSolver.Get_DragCoeff()) + "\t")
                hist_file.write(str(FluidSolver.Get_LiftCoeff()) + "\n")
                hist_file.close()

        self.MPIBarrier()

        self.MPIPrint('\nBGS is converged (strong coupling)')
        self.MPIPrint(' ')
        self.MPIPrint('*************************')
        self.MPIPrint('*  End FSI computation  *')
        self.MPIPrint('*************************')
        self.MPIPrint(' ')

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
