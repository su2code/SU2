#!/usr/bin/env python

## \file FSIInterface.py
#  \brief FSI interface class that handles fluid/solid solvers synchronisation and communication.
#  \authors Nicola Fonzi, Vittorio Cavalieri based on the work of David Thomas
#  \version 8.0.0 "Harrier"
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os
import csv
import numpy as np
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

        self.rootProcess = 0  # the root process is chosen to be MPI rank = 0

        self.nDim = FSI_config["NDIM"]  # problem dimension

        self.haveFluidSolver = (
            False  # True if the fluid solver is initialized on the current rank
        )
        self.haveSolidSolver = (
            False  # True if the solid solver is initialized on the current rank
        )
        self.haveFluidInterface = (
            False  # True if the current rank owns at least one fluid interface node
        )
        self.haveSolidInterface = (
            False  # True if the current rank owns at least one solid interface node
        )

        self.fluidSolverProcessors = (
            list()
        )  # list of partitions where the fluid solver is initialized
        self.solidSolverProcessors = (
            list()
        )  # list of partitions where the solid solver is initialized
        self.fluidInterfaceProcessors = (
            list()
        )  # list of partitions where there are fluid interface nodes
        self.solidInterfaceProcessors = (
            list()
        )  # list of partitions where there are solid interface nodes

        self.fluidInterfaceIdentifier = (
            None  # object that can identify the f/s interface within the fluid solver
        )
        self.solidInterfaceIdentifier = (
            None  # object that can identify the f/s interface within the solid solver
        )

        self.fluidGlobalIndexRange = (
            {}
        )  # contains the global FSI indexing of each fluid interface node for all partitions
        self.solidGlobalIndexRange = (
            {}
        )  # contains the global FSI indexing of each solid interface node for all partitions

        self.FluidHaloNodeList = (
            {}
        )  # contains the the indices (fluid solver indexing) of the halo nodes for each partition
        self.fluidIndexing = (
            {}
        )  # links between the fluid solver indexing and the FSI indexing for the interface nodes
        self.SolidHaloNodeList = (
            {}
        )  # contains the the indices (solid solver indexing) of the halo nodes for each partition
        self.solidIndexing = (
            {}
        )  # links between the solid solver indexing and the FSI indexing for the interface nodes

        self.nLocalFluidInterfaceNodes = 0  # number of nodes (halo nodes included) on the fluid interface, on each partition
        self.nLocalFluidInterfaceHaloNode = (
            0  # number of halo nodes on the fluid intrface, on each partition
        )
        self.nLocalFluidInterfacePhysicalNodes = 0  # number of physical (= non halo) nodes on the fluid interface, on each partition
        self.nFluidInterfaceNodes = np.array(
            int(0)
        )  # number of nodes on the fluid interface, sum over all the partitions
        self.nFluidInterfacePhysicalNodes = np.array(
            int(0)
        )  # number of physical nodes on the fluid interface, sum over all partitions

        self.nLocalSolidInterfaceNodes = (
            0  # number of physical nodes on the solid interface, on each partition
        )
        self.nLocalSolidInterfaceHaloNode = (
            0  # number of halo nodes on the solid intrface, on each partition
        )
        self.nLocalSolidInterfacePhysicalNodes = 0  # number of physical (= non halo) nodes on the solid interface, on each partition
        self.nSolidInterfaceNodes = np.array(
            int(0)
        )  # number of nodes on the solid interface, sum over all partitions
        self.nSolidInterfacePhysicalNodes = np.array(
            int(0)
        )  # number of physical nodes on the solid interface, sum over all partitions

        if FSI_config["MATCHING_MESH"] == "NO" and (
            FSI_config["MESH_INTERP_METHOD"] == "RBF"
            or FSI_config["MESH_INTERP_METHOD"] == "TPS"
        ):
            self.MappingMatrixA = None
            self.MappingMatrixA_T = None
            self.MappingMatrixB = None
            self.MappingMatrixB_T = None
            self.d_RBF = self.nDim + 1
        else:
            self.MappingMatrix = (
                None  # interpolation/mapping matrix for meshes interpolation/mapping
            )
            self.MappingMatrix_T = None  # transposed interpolation/mapping matrix for meshes interpolation/mapping
            self.d_RBF = 0

        self.localFluidInterface_array_X_init = None  # initial fluid interface position on each partition (used for the meshes mapping)
        self.localFluidInterface_array_Y_init = None
        self.localFluidInterface_array_Z_init = None

        self.localSolidInterface_array_X_init = None  # initial solid interface position on each partition (used for mesh mapping)
        self.localSolidInterface_array_Y_init = None
        self.localSolidInterface_array_Z_init = None

        self.solidInterface_array_DispX = None  # solid interface displacement
        self.solidInterface_array_DispY = None
        self.solidInterface_array_DispZ = None

        self.solidInterfaceResidual_array_X = None  # solid interface position residual
        self.solidInterfaceResidual_array_Y = None
        self.solidInterfaceResidual_array_Z = None

        self.solidInterfaceResidualnM1_array_X = (
            None  # solid interface position residual at the previous BGS iteration
        )
        self.solidInterfaceResidualnM1_array_Y = None
        self.solidInterfaceResidualnM1_array_Z = None

        self.fluidInterface_array_DispX = None  # fluid interface displacement
        self.fluidInterface_array_DispY = None
        self.fluidInterface_array_DispZ = None

        self.fluidLoads_array_X = None  # loads on the fluid side of the f/s interface
        self.fluidLoads_array_Y = None
        self.fluidLoads_array_Z = None

        self.solidLoads_array_X = None  # loads on the solid side of the f/s interface
        self.solidLoads_array_Y = None
        self.solidLoads_array_Z = None

        self.aitkenParam = FSI_config[
            "AITKEN_PARAM"
        ]  # relaxation parameter for the BGS method
        self.FSIIter = 0  # current FSI iteration
        self.unsteady = (
            False  # flag for steady or unsteady simulation (default is steady)
        )
        if FSI_config["IMPOSED_MOTION"] == "YES":
            self.ImposedMotion = True
        else:
            self.ImposedMotion = False

        # ---Some screen output ---
        self.MPIPrint("Fluid solver : SU2_CFD")
        self.MPIPrint("Solid solver : {}".format(FSI_config["CSD_SOLVER"]))

        if FSI_config["TIME_MARCHING"] == "YES":
            self.MPIPrint(
                "Unsteady coupled simulation with physical time step : {} s".format(
                    FSI_config["UNST_TIMESTEP"]
                )
            )
            self.unsteady = True
        else:
            self.MPIPrint("Steady coupled simulation")

        if FSI_config["MATCHING_MESH"] == "YES":
            self.MPIPrint("Matching fluid-solid interface")
        else:
            if FSI_config["MESH_INTERP_METHOD"] == "TPS":
                self.MPIPrint(
                    "Non matching fluid-solid interface with Thin Plate Spline interpolation"
                )
            elif FSI_config["MESH_INTERP_METHOD"] == "RBF":
                self.MPIPrint(
                    "Non matching fluid-solid interface with Radial Basis Function interpolation"
                )
                self.RBF_rad = FSI_config["RBF_RADIUS"]
                self.MPIPrint("Radius value : {}".format(self.RBF_rad))
            else:
                self.MPIPrint(
                    "Non matching fluid-solid interface with Nearest Neighboor interpolation"
                )

        self.MPIPrint("Solid predictor : {}".format(FSI_config["DISP_PRED"]))

        self.MPIPrint(
            "Maximum number of FSI iterations : {}".format(FSI_config["NB_FSI_ITER"])
        )

        self.MPIPrint("FSI tolerance : {}".format(FSI_config["FSI_TOLERANCE"]))

        if FSI_config["AITKEN_RELAX"] == "STATIC":
            self.MPIPrint(
                "Static Aitken under-relaxation with constant parameter {}".format(
                    FSI_config["AITKEN_PARAM"]
                )
            )
        elif FSI_config["AITKEN_RELAX"] == "DYNAMIC":
            self.MPIPrint(
                "Dynamic Aitken under-relaxation with initial parameter {}".format(
                    FSI_config["AITKEN_PARAM"]
                )
            )
        else:
            self.MPIPrint("No Aitken under-relaxation")

        self.MPIPrint("FSI interface is set")

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

    def connect(self, FSI_config, FluidSolver, SolidSolver):
        """
        Connection between solvers.
        Creates the communication support between the two solvers.
        Gets information about f/s interfaces from the two solvers.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

        # --- Identify the fluid and solid interfaces and store the number of nodes on both sides (and for each partition) ---
        self.fluidInterfaceIdentifier = None
        self.nLocalFluidInterfaceNodes = 0
        if FluidSolver is not None:
            print("Fluid solver is initialized on process {}".format(myid))
            self.haveFluidSolver = True
            allMovingMarkersTags = FluidSolver.GetDeformableMarkerTags()
            allMarkersID = FluidSolver.GetMarkerTags()
            if not allMovingMarkersTags:
                raise Exception("No interface for FSI was defined.")
            else:
                if allMovingMarkersTags[0] in allMarkersID.keys():
                    self.fluidInterfaceIdentifier = allMarkersID[
                        allMovingMarkersTags[0]
                    ]
            if self.fluidInterfaceIdentifier is not None:
                self.nLocalFluidInterfaceNodes = FluidSolver.GetNumberMarkerNodes(
                    self.fluidInterfaceIdentifier
                )
            if self.nLocalFluidInterfaceNodes != 0:
                self.haveFluidInterface = True
                print(
                    "Number of interface fluid nodes (halo nodes included) on proccess {} : {}".format(
                        myid, self.nLocalFluidInterfaceNodes
                    )
                )
        else:
            pass

        if SolidSolver is not None:
            print("Solid solver is initialized on process {}".format(myid))
            self.haveSolidSolver = True
            self.solidInterfaceIdentifier = SolidSolver.getFSIMarkerID()
            self.nLocalSolidInterfaceNodes = SolidSolver.getNumberOfSolidInterfaceNodes(
                self.solidInterfaceIdentifier
            )
            if self.nLocalSolidInterfaceNodes != 0:
                self.haveSolidInterface = True
                print(
                    "Number of interface solid nodes (halo nodes included) on proccess {} : {}".format(
                        myid, self.nLocalSolidInterfaceNodes
                    )
                )
        else:
            pass

        # --- Exchange information about processors on which the solvers are defined and where the interface nodes are lying ---
        if self.have_MPI:
            if self.haveFluidSolver:
                sendBufFluid = np.array(int(1))
            else:
                sendBufFluid = np.array(int(0))
            if self.haveSolidSolver:
                sendBufSolid = np.array(int(1))
            else:
                sendBufSolid = np.array(int(0))
            if self.haveFluidInterface:
                sendBufFluidInterface = np.array(int(1))
            else:
                sendBufFluidInterface = np.array(int(0))
            if self.haveSolidInterface:
                sendBufSolidInterface = np.array(int(1))
            else:
                sendBufSolidInterface = np.array(int(0))
            rcvBufFluid = np.zeros(MPIsize, dtype=int)
            rcvBufSolid = np.zeros(MPIsize, dtype=int)
            rcvBufFluidInterface = np.zeros(MPIsize, dtype=int)
            rcvBufSolidInterface = np.zeros(MPIsize, dtype=int)
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
            del (
                sendBufFluid,
                sendBufSolid,
                rcvBufFluid,
                rcvBufSolid,
                sendBufFluidInterface,
                sendBufSolidInterface,
                rcvBufFluidInterface,
                rcvBufSolidInterface,
            )
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
            iPoint = FluidSolver.GetMarkerNode(self.fluidInterfaceIdentifier, iVertex)
            if not FluidSolver.GetNodeDomain(iPoint):
                GlobalIndex = FluidSolver.GetNodeGlobalIndex(iPoint)
                self.FluidHaloNodeList[GlobalIndex] = iVertex
                self.nLocalFluidInterfaceHaloNode += 1
        # Calculate the number of physical (= not halo) nodes on each partition
        self.nLocalFluidInterfacePhysicalNodes = (
            self.nLocalFluidInterfaceNodes - self.nLocalFluidInterfaceHaloNode
        )
        if self.have_MPI:
            self.FluidHaloNodeList = self.comm.allgather(self.FluidHaloNodeList)
        else:
            self.FluidHaloNodeList = [{}]

        # Same thing for the solid part
        self.nLocalSolidInterfaceHaloNode = 0
        for iVertex in range(self.nLocalSolidInterfaceNodes):
            if SolidSolver.IsAHaloNode(self.solidInterfaceIdentifier, iVertex):
                GlobalIndex = SolidSolver.getVertexGlobalIndex(
                    self.solidInterfaceIdentifier, iVertex
                )
                self.SolidHaloNodeList[GlobalIndex] = iVertex
                self.nLocalSolidInterfaceHaloNode += 1
        self.nLocalSolidInterfacePhysicalNodes = (
            self.nLocalSolidInterfaceNodes - self.nLocalSolidInterfaceHaloNode
        )
        if self.have_MPI:
            self.SolidHaloNodeList = self.comm.allgather(self.SolidHaloNodeList)
        else:
            self.SolidHaloNodeList = [{}]

        # --- Calculate the total number of nodes (with and without halo) at the fluid interface (sum over all the partitions) and broadcast the number accross all processors ---
        sendBuffTotal = np.array(int(self.nLocalFluidInterfaceNodes))
        sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
        rcvBuffTotal = np.zeros(1, dtype=int)
        rcvBuffPhysical = np.zeros(1, dtype=int)
        if self.have_MPI:
            self.comm.barrier()
            self.comm.Allreduce(
                sendBuffTotal, self.nFluidInterfaceNodes, op=self.MPI.SUM
            )
            self.comm.Allreduce(
                sendBuffPhysical, self.nFluidInterfacePhysicalNodes, op=self.MPI.SUM
            )
        else:
            self.nFluidInterfaceNodes = np.copy(sendBuffTotal)
            self.nFluidInterfacePhysicalNodes = np.copy(sendBuffPhysical)
        del sendBuffTotal, rcvBuffTotal, sendBuffPhysical, rcvBuffPhysical

        # Same thing for the solid part
        sendBuffTotal = np.array(int(self.nLocalSolidInterfaceNodes))
        sendBuffPhysical = np.array(int(self.nLocalSolidInterfacePhysicalNodes))
        rcvBuffTotal = np.zeros(1, dtype=int)
        rcvBuffPhysical = np.zeros(1, dtype=int)
        if self.have_MPI:
            self.comm.barrier()
            self.comm.Allreduce(
                sendBuffTotal, self.nSolidInterfaceNodes, op=self.MPI.SUM
            )
            self.comm.Allreduce(
                sendBuffPhysical, self.nSolidInterfacePhysicalNodes, op=self.MPI.SUM
            )
        else:
            self.nSolidInterfaceNodes = np.copy(sendBuffTotal)
            self.nSolidInterfacePhysicalNodes = np.copy(sendBuffPhysical)
        del sendBuffTotal, rcvBuffTotal, sendBuffPhysical, rcvBuffPhysical

        # --- Store the number of physical interface nodes on each processor and allgather the information ---
        self.fluidPhysicalInterfaceNodesDistribution = np.zeros(MPIsize, dtype=int)
        if self.have_MPI:
            sendBuffPhysical = np.array(int(self.nLocalFluidInterfacePhysicalNodes))
            self.comm.Allgather(
                sendBuffPhysical, self.fluidPhysicalInterfaceNodesDistribution
            )
            del sendBuffPhysical
        else:
            self.fluidPhysicalInterfaceNodesDistribution[
                0
            ] = self.nFluidInterfacePhysicalNodes

        # Same thing for the solid part
        self.solidPhysicalInterfaceNodesDistribution = np.zeros(MPIsize, dtype=int)
        if self.have_MPI:
            sendBuffPhysical = np.array(int(self.nLocalSolidInterfacePhysicalNodes))
            self.comm.Allgather(
                sendBuffPhysical, self.solidPhysicalInterfaceNodesDistribution
            )
            del sendBuffPhysical
        else:
            self.solidPhysicalInterfaceNodesDistribution[
                0
            ] = self.nSolidInterfacePhysicalNodes

        # --- Calculate and store the global indexing of interface physical nodes on each processor and allgather the information ---
        if self.have_MPI:
            if myid in self.fluidInterfaceProcessors:
                globalIndexStart = 0
                for iProc in range(myid):
                    globalIndexStart += self.fluidPhysicalInterfaceNodesDistribution[
                        iProc
                    ]
                globalIndexStop = (
                    globalIndexStart + self.nLocalFluidInterfacePhysicalNodes - 1
                )
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

        # Same thing for the solid part
        if self.have_MPI:
            if myid in self.solidInterfaceProcessors:
                globalIndexStart = 0
                for iProc in range(myid):
                    globalIndexStart += self.solidPhysicalInterfaceNodesDistribution[
                        iProc
                    ]
                globalIndexStop = (
                    globalIndexStart + self.nLocalSolidInterfacePhysicalNodes - 1
                )
            else:
                globalIndexStart = 0
                globalIndexStop = 0
            self.solidGlobalIndexRange[myid] = [globalIndexStart, globalIndexStop]
            self.solidGlobalIndexRange = self.comm.allgather(self.solidGlobalIndexRange)
        else:
            temp = {}
            temp[0] = [0, self.nSolidInterfacePhysicalNodes - 1]
            self.solidGlobalIndexRange = list()
            self.solidGlobalIndexRange.append(temp)

        self.MPIPrint(
            "Total number of fluid interface nodes (halo nodes included) : {}".format(
                self.nFluidInterfaceNodes
            )
        )
        self.MPIPrint(
            "Total number of solid interface nodes (halo nodes included) : {}".format(
                self.nSolidInterfaceNodes
            )
        )
        self.MPIPrint(
            "Total number of fluid interface nodes : {}".format(
                self.nFluidInterfacePhysicalNodes
            )
        )
        self.MPIPrint(
            "Total number of solid interface nodes : {}".format(
                self.nSolidInterfacePhysicalNodes
            )
        )

        self.MPIBarrier()

        # --- Create all the PETSc vectors required for parallel communication and parallel mesh mapping/interpolation (working for serial too) ---
        if self.have_MPI:
            self.solidInterface_array_DispX = PETSc.Vec().create(self.comm)
            self.solidInterface_array_DispY = PETSc.Vec().create(self.comm)
            self.solidInterface_array_DispZ = PETSc.Vec().create(self.comm)
            self.solidInterface_array_DispX.setType("mpi")
            self.solidInterface_array_DispY.setType("mpi")
            self.solidInterface_array_DispZ.setType("mpi")
        else:
            self.solidInterface_array_DispX = PETSc.Vec().create()
            self.solidInterface_array_DispY = PETSc.Vec().create()
            self.solidInterface_array_DispZ = PETSc.Vec().create()
            self.solidInterface_array_DispX.setType("seq")
            self.solidInterface_array_DispY.setType("seq")
            self.solidInterface_array_DispZ.setType("seq")
        self.solidInterface_array_DispX.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterface_array_DispY.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterface_array_DispZ.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterface_array_DispX.set(0.0)
        self.solidInterface_array_DispY.set(0.0)
        self.solidInterface_array_DispZ.set(0.0)

        if self.have_MPI:
            self.fluidInterface_array_DispX = PETSc.Vec().create(self.comm)
            self.fluidInterface_array_DispY = PETSc.Vec().create(self.comm)
            self.fluidInterface_array_DispZ = PETSc.Vec().create(self.comm)
            self.fluidInterface_array_DispX.setType("mpi")
            self.fluidInterface_array_DispY.setType("mpi")
            self.fluidInterface_array_DispZ.setType("mpi")
        else:
            self.fluidInterface_array_DispX = PETSc.Vec().create()
            self.fluidInterface_array_DispY = PETSc.Vec().create()
            self.fluidInterface_array_DispZ = PETSc.Vec().create()
            self.fluidInterface_array_DispX.setType("seq")
            self.fluidInterface_array_DispY.setType("seq")
            self.fluidInterface_array_DispZ.setType("seq")
        self.fluidInterface_array_DispX.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidInterface_array_DispY.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidInterface_array_DispZ.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidInterface_array_DispX.set(0.0)
        self.fluidInterface_array_DispY.set(0.0)
        self.fluidInterface_array_DispZ.set(0.0)

        if self.have_MPI:
            self.fluidLoads_array_X = PETSc.Vec().create(self.comm)
            self.fluidLoads_array_Y = PETSc.Vec().create(self.comm)
            self.fluidLoads_array_Z = PETSc.Vec().create(self.comm)
            self.fluidLoads_array_X.setType("mpi")
            self.fluidLoads_array_Y.setType("mpi")
            self.fluidLoads_array_Z.setType("mpi")
        else:
            self.fluidLoads_array_X = PETSc.Vec().create()
            self.fluidLoads_array_Y = PETSc.Vec().create()
            self.fluidLoads_array_Z = PETSc.Vec().create()
            self.fluidLoads_array_X.setType("seq")
            self.fluidLoads_array_Y.setType("seq")
            self.fluidLoads_array_Z.setType("seq")
        self.fluidLoads_array_X.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidLoads_array_Y.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidLoads_array_Z.setSizes(self.nFluidInterfacePhysicalNodes)
        self.fluidLoads_array_X.set(0.0)
        self.fluidLoads_array_Y.set(0.0)
        self.fluidLoads_array_Z.set(0.0)

        if self.have_MPI:
            self.solidLoads_array_X = PETSc.Vec().create(self.comm)
            self.solidLoads_array_Y = PETSc.Vec().create(self.comm)
            self.solidLoads_array_Z = PETSc.Vec().create(self.comm)
            self.solidLoads_array_X.setType("mpi")
            self.solidLoads_array_Y.setType("mpi")
            self.solidLoads_array_Z.setType("mpi")
        else:
            self.solidLoads_array_X = PETSc.Vec().create()
            self.solidLoads_array_Y = PETSc.Vec().create()
            self.solidLoads_array_Z = PETSc.Vec().create()
            self.solidLoads_array_X.setType("seq")
            self.solidLoads_array_Y.setType("seq")
            self.solidLoads_array_Z.setType("seq")
        self.solidLoads_array_X.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        self.solidLoads_array_Y.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        self.solidLoads_array_Z.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        self.solidLoads_array_X.set(0.0)
        self.solidLoads_array_Y.set(0.0)
        self.solidLoads_array_Z.set(0.0)

        # --- Create the PETSc vectors required for parallel relaxed BGS algo (working for serial too) ---
        if self.have_MPI:
            self.solidInterfaceResidual_array_X = PETSc.Vec().create(self.comm)
            self.solidInterfaceResidual_array_Y = PETSc.Vec().create(self.comm)
            self.solidInterfaceResidual_array_Z = PETSc.Vec().create(self.comm)
            self.solidInterfaceResidual_array_X.setType("mpi")
            self.solidInterfaceResidual_array_Y.setType("mpi")
            self.solidInterfaceResidual_array_Z.setType("mpi")
        else:
            self.solidInterfaceResidual_array_X = PETSc.Vec().create()
            self.solidInterfaceResidual_array_Y = PETSc.Vec().create()
            self.solidInterfaceResidual_array_Z = PETSc.Vec().create()
            self.solidInterfaceResidual_array_X.setType("seq")
            self.solidInterfaceResidual_array_Y.setType("seq")
            self.solidInterfaceResidual_array_Z.setType("seq")
        self.solidInterfaceResidual_array_X.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterfaceResidual_array_Y.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterfaceResidual_array_Z.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterfaceResidual_array_X.set(0.0)
        self.solidInterfaceResidual_array_Y.set(0.0)
        self.solidInterfaceResidual_array_Z.set(0.0)

        if self.have_MPI:
            self.solidInterfaceResidualnM1_array_X = PETSc.Vec().create(self.comm)
            self.solidInterfaceResidualnM1_array_Y = PETSc.Vec().create(self.comm)
            self.solidInterfaceResidualnM1_array_Z = PETSc.Vec().create(self.comm)
            self.solidInterfaceResidualnM1_array_X.setType("mpi")
            self.solidInterfaceResidualnM1_array_Y.setType("mpi")
            self.solidInterfaceResidualnM1_array_Z.setType("mpi")
        else:
            self.solidInterfaceResidualnM1_array_X = PETSc.Vec().create()
            self.solidInterfaceResidualnM1_array_Y = PETSc.Vec().create()
            self.solidInterfaceResidualnM1_array_Z = PETSc.Vec().create()
            self.solidInterfaceResidualnM1_array_X.setType("seq")
            self.solidInterfaceResidualnM1_array_Y.setType("seq")
            self.solidInterfaceResidualnM1_array_Z.setType("seq")
        self.solidInterfaceResidualnM1_array_X.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterfaceResidualnM1_array_Y.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterfaceResidualnM1_array_Z.setSizes(
            self.nSolidInterfacePhysicalNodes + self.d_RBF
        )
        self.solidInterfaceResidualnM1_array_X.set(0.0)
        self.solidInterfaceResidualnM1_array_Y.set(0.0)
        self.solidInterfaceResidualnM1_array_Z.set(0.0)

    def interfaceMapping(self, FluidSolver, SolidSolver, FSI_config):
        """
        Creates the one-to-one mapping between interfaces in case of matching meshes.
        Creates the interpolation rules between interfaces in case of non-matching meshes.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

        # --- Get the fluid interface from fluid solver on each partition ---
        GlobalIndex = int()
        localIndex = 0
        fluidIndexing_temp = {}
        self.localFluidInterface_array_X_init = np.zeros(
            (self.nLocalFluidInterfacePhysicalNodes)
        )
        self.localFluidInterface_array_Y_init = np.zeros(
            (self.nLocalFluidInterfacePhysicalNodes)
        )
        self.localFluidInterface_array_Z_init = np.zeros(
            (self.nLocalFluidInterfacePhysicalNodes)
        )
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            # Note that the fluid solver is separated in more processors outside the python script
            # thus when, from a core, we request for the vertices on the interface, we only obtain
            # those in that core
            iPoint = FluidSolver.GetMarkerNode(self.fluidInterfaceIdentifier, iVertex)
            GlobalIndex = FluidSolver.GetNodeGlobalIndex(iPoint)
            if self.nDim == 2:
                posx, posy = FluidSolver.InitialCoordinates().Get(iPoint)
                posz = 0
            else:
                posx, posy, posz = FluidSolver.InitialCoordinates().Get(iPoint)
            if GlobalIndex not in self.FluidHaloNodeList[myid].keys():
                fluidIndexing_temp[GlobalIndex] = self.__getGlobalIndex(
                    "fluid", myid, localIndex
                )
                self.localFluidInterface_array_X_init[localIndex] = posx
                self.localFluidInterface_array_Y_init[localIndex] = posy
                self.localFluidInterface_array_Z_init[localIndex] = posz
                localIndex += 1
        if self.have_MPI:
            fluidIndexing_temp = self.comm.allgather(fluidIndexing_temp)
            for ii in range(len(fluidIndexing_temp)):
                for key, value in fluidIndexing_temp[ii].items():
                    # This contains the link between the global index in python and that in SU2
                    self.fluidIndexing[key] = value
        else:
            self.fluidIndexing = fluidIndexing_temp.copy()
        del fluidIndexing_temp

        # --- Get the solid interface from solid solver on each partition ---
        localIndex = 0
        solidIndexing_temp = {}
        self.localSolidInterface_array_X_init = np.zeros(self.nLocalSolidInterfaceNodes)
        self.localSolidInterface_array_Y_init = np.zeros(self.nLocalSolidInterfaceNodes)
        self.localSolidInterface_array_Z_init = np.zeros(self.nLocalSolidInterfaceNodes)
        for iVertex in range(self.nLocalSolidInterfaceNodes):
            GlobalIndex = SolidSolver.getVertexGlobalIndex(
                self.solidInterfaceIdentifier, iVertex
            )
            posx, posy, posz = SolidSolver.getInterfaceNodePosInit(
                self.solidInterfaceIdentifier, iVertex
            )
            if GlobalIndex not in self.SolidHaloNodeList[myid].keys():
                solidIndexing_temp[GlobalIndex] = self.__getGlobalIndex(
                    "solid", myid, localIndex
                )
                self.localSolidInterface_array_X_init[localIndex] = posx
                self.localSolidInterface_array_Y_init[localIndex] = posy
                self.localSolidInterface_array_Z_init[localIndex] = posz
                localIndex += 1
        if self.have_MPI:
            solidIndexing_temp = self.comm.allgather(solidIndexing_temp)
            for ii in range(len(solidIndexing_temp)):
                for key, value in solidIndexing_temp[ii].items():
                    self.solidIndexing[key] = value
        else:
            self.solidIndexing = solidIndexing_temp.copy()
        del solidIndexing_temp

        # --- Create the PETSc parallel interpolation matrix ---
        if FSI_config["MATCHING_MESH"] == "NO" and (
            FSI_config["MESH_INTERP_METHOD"] == "RBF"
            or FSI_config["MESH_INTERP_METHOD"] == "TPS"
        ):
            if self.have_MPI:
                self.MappingMatrixA = PETSc.Mat().create(self.comm)
                self.MappingMatrixB = PETSc.Mat().create(self.comm)
                self.MappingMatrixA_T = PETSc.Mat().create(self.comm)
                self.MappingMatrixB_T = PETSc.Mat().create(self.comm)
                if FSI_config["MESH_INTERP_METHOD"] == "RBF":
                    self.MappingMatrixA.setType("mpiaij")
                    self.MappingMatrixB.setType("mpiaij")
                    self.MappingMatrixA_T.setType("mpiaij")
                    self.MappingMatrixB_T.setType("mpiaij")
                else:
                    self.MappingMatrixA.setType("mpiaij")
                    self.MappingMatrixB.setType("mpiaij")
                    self.MappingMatrixA_T.setType("mpiaij")
                    self.MappingMatrixB_T.setType("mpiaij")
            else:
                self.MappingMatrixA = PETSc.Mat().create()
                self.MappingMatrixB = PETSc.Mat().create()
                self.MappingMatrixA_T = PETSc.Mat().create()
                self.MappingMatrixB_T = PETSc.Mat().create()
                if FSI_config["MESH_INTERP_METHOD"] == "RBF":
                    self.MappingMatrixA.setType("aij")
                    self.MappingMatrixB.setType("aij")
                    self.MappingMatrixA_T.setType("aij")
                    self.MappingMatrixB_T.setType("aij")
                else:
                    self.MappingMatrixA.setType("aij")
                    self.MappingMatrixB.setType("aij")
                    self.MappingMatrixA_T.setType("aij")
                    self.MappingMatrixB_T.setType("aij")
            self.MappingMatrixA.setSizes(
                (
                    self.nSolidInterfacePhysicalNodes + self.d_RBF,
                    self.nSolidInterfacePhysicalNodes + self.d_RBF,
                )
            )
            self.MappingMatrixA.setUp()
            self.MappingMatrixA.setOption(
                PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False
            )
            self.MappingMatrixB.setSizes(
                (
                    self.nFluidInterfacePhysicalNodes,
                    self.nSolidInterfacePhysicalNodes + self.d_RBF,
                )
            )
            self.MappingMatrixB.setUp()
            self.MappingMatrixB.setOption(
                PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False
            )
            self.MappingMatrixA_T.setSizes(
                (
                    self.nSolidInterfacePhysicalNodes + self.d_RBF,
                    self.nSolidInterfacePhysicalNodes + self.d_RBF,
                )
            )
            self.MappingMatrixA_T.setUp()
            self.MappingMatrixA_T.setOption(
                PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False
            )
            self.MappingMatrixB_T.setSizes(
                (
                    self.nSolidInterfacePhysicalNodes + self.d_RBF,
                    self.nFluidInterfacePhysicalNodes,
                )
            )
            self.MappingMatrixB_T.setUp()
            self.MappingMatrixB_T.setOption(
                PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False
            )
        else:
            if self.have_MPI:
                self.MappingMatrix = PETSc.Mat().create(self.comm)
                self.MappingMatrix_T = PETSc.Mat().create(self.comm)
                self.MappingMatrix.setType("mpiaij")
                self.MappingMatrix_T.setType("mpiaij")
            else:
                self.MappingMatrix = PETSc.Mat().create()
                self.MappingMatrix_T = PETSc.Mat().create()
                self.MappingMatrix.setType("aij")
                self.MappingMatrix_T.setType("aij")
            self.MappingMatrix.setSizes(
                (self.nFluidInterfacePhysicalNodes, self.nSolidInterfacePhysicalNodes)
            )
            self.MappingMatrix.setUp()
            self.MappingMatrix.setOption(
                PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False
            )
            self.MappingMatrix_T.setSizes(
                (self.nSolidInterfacePhysicalNodes, self.nFluidInterfacePhysicalNodes)
            )
            self.MappingMatrix_T.setUp()
            self.MappingMatrix_T.setOption(
                PETSc.Mat().Option.NEW_NONZERO_ALLOCATION_ERR, False
            )

        # --- Fill the interpolation matrix in parallel (working in serial too) ---
        if FSI_config["MATCHING_MESH"] == "NO" and (
            FSI_config["MESH_INTERP_METHOD"] == "RBF"
            or FSI_config["MESH_INTERP_METHOD"] == "TPS"
        ):
            self.MPIPrint("Building interpolation matrices...")
            if self.have_MPI:
                for iProc in self.solidInterfaceProcessors:
                    if myid == iProc:
                        for jProc in self.solidInterfaceProcessors:
                            if jProc != iProc:
                                self.comm.Send(
                                    self.localSolidInterface_array_X_init,
                                    dest=jProc,
                                    tag=1,
                                )
                                self.comm.Send(
                                    self.localSolidInterface_array_Y_init,
                                    dest=jProc,
                                    tag=2,
                                )
                                self.comm.Send(
                                    self.localSolidInterface_array_Z_init,
                                    dest=jProc,
                                    tag=3,
                                )
                            else:
                                solidInterfaceBuffRcv_X = np.copy(
                                    self.localSolidInterface_array_X_init
                                )
                                solidInterfaceBuffRcv_Y = np.copy(
                                    self.localSolidInterface_array_Y_init
                                )
                                solidInterfaceBuffRcv_Z = np.copy(
                                    self.localSolidInterface_array_Z_init
                                )
                    if myid in self.solidInterfaceProcessors:
                        if myid != iProc:
                            sizeOfBuff = self.solidPhysicalInterfaceNodesDistribution[
                                iProc
                            ]
                            solidInterfaceBuffRcv_X = np.empty(
                                sizeOfBuff, dtype=np.float64
                            )
                            solidInterfaceBuffRcv_Y = np.empty(
                                sizeOfBuff, dtype=np.float64
                            )
                            solidInterfaceBuffRcv_Z = np.empty(
                                sizeOfBuff, dtype=np.float64
                            )
                            self.comm.Recv(solidInterfaceBuffRcv_X, source=iProc, tag=1)
                            self.comm.Recv(solidInterfaceBuffRcv_Y, source=iProc, tag=2)
                            self.comm.Recv(solidInterfaceBuffRcv_Z, source=iProc, tag=3)
                        if FSI_config["MESH_INTERP_METHOD"] == "RBF":
                            self.RBFMeshMapping_A(
                                solidInterfaceBuffRcv_X,
                                solidInterfaceBuffRcv_Y,
                                solidInterfaceBuffRcv_Z,
                                iProc,
                                self.RBF_rad,
                            )
                        else:
                            self.TPSMeshMapping_A(
                                solidInterfaceBuffRcv_X,
                                solidInterfaceBuffRcv_Y,
                                solidInterfaceBuffRcv_Z,
                                iProc,
                            )
            else:
                if FSI_config["MESH_INTERP_METHOD"] == "RBF":
                    self.RBFMeshMapping_A(
                        self.localSolidInterface_array_X_init,
                        self.localSolidInterface_array_Y_init,
                        self.localSolidInterface_array_Z_init,
                        0,
                        self.RBF_rad,
                    )
                else:
                    self.TPSMeshMapping_A(
                        self.localSolidInterface_array_X_init,
                        self.localSolidInterface_array_Y_init,
                        self.localSolidInterface_array_Z_init,
                        0,
                    )
            self.MappingMatrixA.assemblyBegin()
            self.MappingMatrixA.assemblyEnd()
            self.MappingMatrixA_T.assemblyBegin()
            self.MappingMatrixA_T.assemblyEnd()
            self.MPIPrint("Matrix A is built.")
        else:
            self.MPIPrint("Building interpolation matrix...")
        self.MPIBarrier()
        if self.have_MPI:
            for iProc in self.solidInterfaceProcessors:
                if myid == iProc:
                    for jProc in self.fluidInterfaceProcessors:
                        if jProc != iProc:
                            self.comm.Send(
                                self.localSolidInterface_array_X_init, dest=jProc, tag=1
                            )
                            self.comm.Send(
                                self.localSolidInterface_array_Y_init, dest=jProc, tag=2
                            )
                            self.comm.Send(
                                self.localSolidInterface_array_Z_init, dest=jProc, tag=3
                            )
                        else:
                            solidInterfaceBuffRcv_X = np.copy(
                                self.localSolidInterface_array_X_init
                            )
                            solidInterfaceBuffRcv_Y = np.copy(
                                self.localSolidInterface_array_Y_init
                            )
                            solidInterfaceBuffRcv_Z = np.copy(
                                self.localSolidInterface_array_Z_init
                            )
                if myid in self.fluidInterfaceProcessors:
                    if myid != iProc:
                        sizeOfBuff = self.solidPhysicalInterfaceNodesDistribution[iProc]
                        solidInterfaceBuffRcv_X = np.empty(sizeOfBuff, dtype=np.float64)
                        solidInterfaceBuffRcv_Y = np.empty(sizeOfBuff, dtype=np.float64)
                        solidInterfaceBuffRcv_Z = np.empty(sizeOfBuff, dtype=np.float64)
                        self.comm.Recv(solidInterfaceBuffRcv_X, source=iProc, tag=1)
                        self.comm.Recv(solidInterfaceBuffRcv_Y, source=iProc, tag=2)
                        self.comm.Recv(solidInterfaceBuffRcv_Z, source=iProc, tag=3)
                    if FSI_config["MATCHING_MESH"] == "NO":
                        if FSI_config["MESH_INTERP_METHOD"] == "RBF":
                            self.RBFMeshMapping_B(
                                solidInterfaceBuffRcv_X,
                                solidInterfaceBuffRcv_Y,
                                solidInterfaceBuffRcv_Z,
                                iProc,
                                self.RBF_rad,
                            )
                        elif FSI_config["MESH_INTERP_METHOD"] == "TPS":
                            self.TPSMeshMapping_B(
                                solidInterfaceBuffRcv_X,
                                solidInterfaceBuffRcv_Y,
                                solidInterfaceBuffRcv_Z,
                                iProc,
                            )
                        else:
                            self.NearestNeighboorMeshMapping(
                                solidInterfaceBuffRcv_X,
                                solidInterfaceBuffRcv_Y,
                                solidInterfaceBuffRcv_Z,
                                iProc,
                            )
                    else:
                        self.matchingMeshMapping(
                            solidInterfaceBuffRcv_X,
                            solidInterfaceBuffRcv_Y,
                            solidInterfaceBuffRcv_Z,
                            iProc,
                        )
        else:
            if FSI_config["MATCHING_MESH"] == "NO":
                if FSI_config["MESH_INTERP_METHOD"] == "RBF":
                    self.RBFMeshMapping_B(
                        self.localSolidInterface_array_X_init,
                        self.localSolidInterface_array_Y_init,
                        self.localSolidInterface_array_Z_init,
                        0,
                        self.RBF_rad,
                    )
                elif FSI_config["MESH_INTERP_METHOD"] == "TPS":
                    self.TPSMeshMapping_B(
                        self.localSolidInterface_array_X_init,
                        self.localSolidInterface_array_Y_init,
                        self.localSolidInterface_array_Z_init,
                        0,
                    )
                else:
                    self.NearestNeighboorMeshMapping(
                        self.localSolidInterface_array_X_init,
                        self.localSolidInterface_array_Y_init,
                        self.localSolidInterface_array_Z_init,
                        0,
                    )
            else:
                self.matchingMeshMapping(
                    self.localSolidInterface_array_X_init,
                    self.localSolidInterface_array_Y_init,
                    self.localSolidInterface_array_Z_init,
                    0,
                )

        if FSI_config["MATCHING_MESH"] == "NO" and (
            FSI_config["MESH_INTERP_METHOD"] == "RBF"
            or FSI_config["MESH_INTERP_METHOD"] == "TPS"
        ):
            self.MappingMatrixB.assemblyBegin()
            self.MappingMatrixB.assemblyEnd()
            self.MappingMatrixB_T.assemblyBegin()
            self.MappingMatrixB_T.assemblyEnd()
            self.MPIPrint("Matrix B is built.")
        else:
            self.MappingMatrix.assemblyBegin()
            self.MappingMatrix.assemblyEnd()
            self.MappingMatrix_T.assemblyBegin()
            self.MappingMatrix_T.assemblyEnd()
            self.MPIPrint("Interpolation matrix is built.")

        self.MPIBarrier()

        del self.localSolidInterface_array_X_init
        del self.localSolidInterface_array_Y_init
        del self.localSolidInterface_array_Z_init
        del self.localFluidInterface_array_X_init
        del self.localFluidInterface_array_Y_init
        del self.localFluidInterface_array_Z_init

    def matchingMeshMapping(
        self,
        solidInterfaceBuffRcv_X,
        solidInterfaceBuffRcv_Y,
        solidInterfaceBuffRcv_Z,
        iProc,
    ):
        """
        Fill the mapping matrix in case of matching meshes at the f/s interface.
        """
        if self.have_MPI:
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
            if self.nDim == 2:
                SolidSpatialTree.add(jVertex, (posX, posY))
            else:
                SolidSpatialTree.add(jVertex, (posX, posY, posZ))

        if self.nFluidInterfacePhysicalNodes != self.nSolidInterfacePhysicalNodes:
            raise Exception(
                "Fluid and solid interface must have the same number of nodes for matching meshes ! "
            )

        # --- For each fluid interface node, find the nearest solid interface node and fill the boolean mapping matrix ---
        for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
            posX = self.localFluidInterface_array_X_init[iVertexFluid]
            posY = self.localFluidInterface_array_Y_init[iVertexFluid]
            posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
            if self.nDim == 2:
                neighboors = list(SolidSpatialTree.nearest((posX, posY), 1))
            elif self.nDim == 3:
                neighboors = list(SolidSpatialTree.nearest((posX, posY, posZ), 1))
            jVertexSolid = neighboors[0]
            # Check if the distance is small enough to ensure coincidence
            NodeA = np.array([posX, posY, posZ])
            NodeB = np.array(
                [
                    solidInterfaceBuffRcv_X[jVertexSolid],
                    solidInterfaceBuffRcv_Y[jVertexSolid],
                    solidInterfaceBuffRcv_Z[jVertexSolid],
                ]
            )
            distance = spdist.euclidean(NodeA, NodeB)
            iGlobalVertexFluid = self.__getGlobalIndex("fluid", myid, iVertexFluid)
            jGlobalVertexSolid = self.__getGlobalIndex("solid", iProc, jVertexSolid)
            if distance > 1e-6:
                print(
                    "WARNING : Tolerance for matching meshes is not matched between node F{} and S{} : ({}, {}, {})<-->({}, {}, {}) , DISTANCE : {} !".format(
                        iGlobalVertexFluid,
                        jGlobalVertexSolid,
                        posX,
                        posY,
                        posZ,
                        solidInterfaceBuffRcv_X[jVertexSolid],
                        solidInterfaceBuffRcv_Y[jVertexSolid],
                        solidInterfaceBuffRcv_Z[jVertexSolid],
                        distance,
                    )
                )
            self.MappingMatrix.setValue(iGlobalVertexFluid, jGlobalVertexSolid, 1.0)
            self.MappingMatrix_T.setValue(jGlobalVertexSolid, iGlobalVertexFluid, 1.0)

    def NearestNeighboorMeshMapping(
        self,
        solidInterfaceBuffRcv_X,
        solidInterfaceBuffRcv_Y,
        solidInterfaceBuffRcv_Z,
        iProc,
    ):
        """
        Interpolation based on the nearest neighboor.
        For each node, the mesh is scanned to find the closed node to the first
        one.
        """

        if self.have_MPI:
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
            if self.nDim == 2:
                SolidSpatialTree.add(jVertex, (posX, posY))
            else:
                SolidSpatialTree.add(jVertex, (posX, posY, posZ))

        # --- For each fluid interface node, find the nearest solid interface node and fill the boolean mapping matrix ---
        for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
            posX = self.localFluidInterface_array_X_init[iVertexFluid]
            posY = self.localFluidInterface_array_Y_init[iVertexFluid]
            posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
            if self.nDim == 2:
                neighboors = list(SolidSpatialTree.nearest((posX, posY), 1))
            elif self.nDim == 3:
                neighboors = list(SolidSpatialTree.nearest((posX, posY, posZ), 1))
            jVertexSolid = neighboors[0]
            iGlobalVertexFluid = self.__getGlobalIndex("fluid", myid, iVertexFluid)
            jGlobalVertexSolid = self.__getGlobalIndex("solid", iProc, jVertexSolid)
            self.MappingMatrix.setValue(iGlobalVertexFluid, jGlobalVertexSolid, 1.0)
            self.MappingMatrix_T.setValue(jGlobalVertexSolid, iGlobalVertexFluid, 1.0)

    def RBFMeshMapping_A(
        self,
        solidInterfaceBuffRcv_X,
        solidInterfaceBuffRcv_Y,
        solidInterfaceBuffRcv_Z,
        iProc,
        rad,
    ):
        """
        First part of the RBF mapping. This method provides the matrix required to
        obtain, from the structural displacements, the loadings of the kernel
        functions.
        """

        if self.have_MPI:
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
            if self.nDim == 2:
                SolidSpatialTree.add(jVertex, (posX, posY))
            else:
                SolidSpatialTree.add(jVertex, (posX, posY, posZ))

        for iVertexSolid in range(self.nLocalSolidInterfacePhysicalNodes):
            posX = self.localSolidInterface_array_X_init[iVertexSolid]
            posY = self.localSolidInterface_array_Y_init[iVertexSolid]
            posZ = self.localSolidInterface_array_Z_init[iVertexSolid]
            NodeA = np.array([posX, posY, posZ])
            iGlobalVertexSolid = self.__getGlobalIndex("solid", myid, iVertexSolid)
            if self.nDim == 2:
                neighboors = list(
                    SolidSpatialTree.intersection(
                        (posX - rad, posY - rad, posX + rad, posY + rad)
                    )
                )
            elif self.nDim == 3:
                neighboors = list(
                    SolidSpatialTree.intersection(
                        (
                            posX - rad,
                            posY - rad,
                            posZ - rad,
                            posX + rad,
                            posY + rad,
                            posZ + rad,
                        )
                    )
                )
            for jVertexSolid in neighboors:
                NodeB = np.array(
                    [
                        solidInterfaceBuffRcv_X[jVertexSolid],
                        solidInterfaceBuffRcv_Y[jVertexSolid],
                        solidInterfaceBuffRcv_Z[jVertexSolid],
                    ]
                )
                distance = spdist.euclidean(NodeA, NodeB)
                phi = self.__CPC2(distance, rad)
                jGlobalVertexSolid = self.__getGlobalIndex("solid", iProc, jVertexSolid)
                self.MappingMatrixA.setValue(
                    iGlobalVertexSolid, jGlobalVertexSolid, phi
                )
                self.MappingMatrixA_T.setValue(
                    jGlobalVertexSolid, iGlobalVertexSolid, phi
                )
            self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes, 1.0)
            self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes + 1, posX)
            self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes + 2, posY)
            if self.nDim == 3:
                self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes + 3, posZ)
            self.MappingMatrixA_T.setValue(nSolidNodes, iGlobalVertexSolid, 1.0)
            self.MappingMatrixA_T.setValue(nSolidNodes + 1, iGlobalVertexSolid, posX)
            self.MappingMatrixA_T.setValue(nSolidNodes + 2, iGlobalVertexSolid, posY)
            if self.nDim == 3:
                self.MappingMatrixA_T.setValue(
                    nSolidNodes + 3, iGlobalVertexSolid, posZ
                )

    def RBFMeshMapping_B(
        self,
        solidInterfaceBuffRcv_X,
        solidInterfaceBuffRcv_Y,
        solidInterfaceBuffRcv_Z,
        iProc,
        rad,
    ):
        """
        Second part of the RBF mapping. This method provides the matrix required to
        obtain, from the kernel function loadings, the fluid nodes displacements.
        """

        if self.have_MPI:
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
            if self.nDim == 2:
                SolidSpatialTree.add(jVertex, (posX, posY))
            else:
                SolidSpatialTree.add(jVertex, (posX, posY, posZ))

        for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
            posX = self.localFluidInterface_array_X_init[iVertexFluid]
            posY = self.localFluidInterface_array_Y_init[iVertexFluid]
            posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
            NodeA = np.array([posX, posY, posZ])
            iGlobalVertexFluid = self.__getGlobalIndex("fluid", myid, iVertexFluid)
            if self.nDim == 2:
                neighboors = list(
                    SolidSpatialTree.intersection(
                        (posX - rad, posY - rad, posX + rad, posY + rad)
                    )
                )
            elif self.nDim == 3:
                neighboors = list(
                    SolidSpatialTree.intersection(
                        (
                            posX - rad,
                            posY - rad,
                            posZ - rad,
                            posX + rad,
                            posY + rad,
                            posZ + rad,
                        )
                    )
                )
            for jVertexSolid in neighboors:
                NodeB = np.array(
                    [
                        solidInterfaceBuffRcv_X[jVertexSolid],
                        solidInterfaceBuffRcv_Y[jVertexSolid],
                        solidInterfaceBuffRcv_Z[jVertexSolid],
                    ]
                )
                distance = spdist.euclidean(NodeA, NodeB)
                phi = self.__CPC2(distance, rad)
                jGlobalVertexSolid = self.__getGlobalIndex("solid", iProc, jVertexSolid)
                self.MappingMatrixB.setValue(
                    iGlobalVertexFluid, jGlobalVertexSolid, phi
                )
                self.MappingMatrixB_T.setValue(
                    jGlobalVertexSolid, iGlobalVertexFluid, phi
                )
            self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes, 1.0)
            self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes + 1, posX)
            self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes + 2, posY)
            if self.nDim == 3:
                self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes + 3, posZ)
            self.MappingMatrixB_T.setValue(nSolidNodes, iGlobalVertexFluid, 1.0)
            self.MappingMatrixB_T.setValue(nSolidNodes + 1, iGlobalVertexFluid, posX)
            self.MappingMatrixB_T.setValue(nSolidNodes + 2, iGlobalVertexFluid, posY)
            if self.nDim == 3:
                self.MappingMatrixB_T.setValue(
                    nSolidNodes + 3, iGlobalVertexFluid, posZ
                )

    def TPSMeshMapping_A(
        self,
        solidInterfaceBuffRcv_X,
        solidInterfaceBuffRcv_Y,
        solidInterfaceBuffRcv_Z,
        iProc,
    ):
        """
        First part of the RBF mapping. This method provides the matrix required to
        obtain, from the structural displacements, the loadings of the kernel
        functions.
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        nSolidNodes = solidInterfaceBuffRcv_X.shape[0]

        for iVertexSolid in range(self.nLocalSolidInterfacePhysicalNodes):
            posX = self.localSolidInterface_array_X_init[iVertexSolid]
            posY = self.localSolidInterface_array_Y_init[iVertexSolid]
            posZ = self.localSolidInterface_array_Z_init[iVertexSolid]
            NodeA = np.array([posX, posY, posZ])
            iGlobalVertexSolid = self.__getGlobalIndex("solid", myid, iVertexSolid)
            for jVertexSolid in range(nSolidNodes):
                NodeB = np.array(
                    [
                        solidInterfaceBuffRcv_X[jVertexSolid],
                        solidInterfaceBuffRcv_Y[jVertexSolid],
                        solidInterfaceBuffRcv_Z[jVertexSolid],
                    ]
                )
                distance = spdist.euclidean(NodeA, NodeB)
                phi = self.__TPS(distance)
                jGlobalVertexSolid = self.__getGlobalIndex("solid", iProc, jVertexSolid)
                self.MappingMatrixA.setValue(
                    iGlobalVertexSolid, jGlobalVertexSolid, phi
                )
                self.MappingMatrixA_T.setValue(
                    jGlobalVertexSolid, iGlobalVertexSolid, phi
                )
            self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes, 1.0)
            self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes + 1, posX)
            self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes + 2, posY)
            if self.nDim == 3:
                self.MappingMatrixA.setValue(iGlobalVertexSolid, nSolidNodes + 3, posZ)
            self.MappingMatrixA_T.setValue(nSolidNodes, iGlobalVertexSolid, 1.0)
            self.MappingMatrixA_T.setValue(nSolidNodes + 1, iGlobalVertexSolid, posX)
            self.MappingMatrixA_T.setValue(nSolidNodes + 2, iGlobalVertexSolid, posY)
            if self.nDim == 3:
                self.MappingMatrixA_T.setValue(
                    nSolidNodes + 3, iGlobalVertexSolid, posZ
                )

    def TPSMeshMapping_B(
        self,
        solidInterfaceBuffRcv_X,
        solidInterfaceBuffRcv_Y,
        solidInterfaceBuffRcv_Z,
        iProc,
    ):
        """
        Second part of the TPS mapping. This method provides the matrix required to
        obtain, from the kernel function loadings, the fluid nodes displacements.
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        nSolidNodes = solidInterfaceBuffRcv_X.shape[0]

        for iVertexFluid in range(self.nLocalFluidInterfacePhysicalNodes):
            posX = self.localFluidInterface_array_X_init[iVertexFluid]
            posY = self.localFluidInterface_array_Y_init[iVertexFluid]
            posZ = self.localFluidInterface_array_Z_init[iVertexFluid]
            NodeA = np.array([posX, posY, posZ])
            iGlobalVertexFluid = self.__getGlobalIndex("fluid", myid, iVertexFluid)
            for jVertexSolid in range(nSolidNodes):
                NodeB = np.array(
                    [
                        solidInterfaceBuffRcv_X[jVertexSolid],
                        solidInterfaceBuffRcv_Y[jVertexSolid],
                        solidInterfaceBuffRcv_Z[jVertexSolid],
                    ]
                )
                distance = spdist.euclidean(NodeA, NodeB)
                phi = self.__TPS(distance)
                jGlobalVertexSolid = self.__getGlobalIndex("solid", iProc, jVertexSolid)
                self.MappingMatrixB.setValue(
                    iGlobalVertexFluid, jGlobalVertexSolid, phi
                )
                self.MappingMatrixB_T.setValue(
                    jGlobalVertexSolid, iGlobalVertexFluid, phi
                )
            self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes, 1.0)
            self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes + 1, posX)
            self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes + 2, posY)
            if self.nDim == 3:
                self.MappingMatrixB.setValue(iGlobalVertexFluid, nSolidNodes + 3, posZ)
            self.MappingMatrixB_T.setValue(nSolidNodes, iGlobalVertexFluid, 1.0)
            self.MappingMatrixB_T.setValue(nSolidNodes + 1, iGlobalVertexFluid, posX)
            self.MappingMatrixB_T.setValue(nSolidNodes + 2, iGlobalVertexFluid, posY)
            if self.nDim == 3:
                self.MappingMatrixB_T.setValue(
                    nSolidNodes + 3, iGlobalVertexFluid, posZ
                )

    def __CPC2(self, distance, rad):
        """
        This method provides the value of the kernel function given the euclidean
        distance. The kernel function is the one used for RBF.
        """
        phi = 0.0
        eps = distance / rad

        if eps < 1:
            phi = ((1.0 - eps) ** 4) * (4.0 * eps + 1.0)
        else:
            phi = 0.0

        return phi

    def __TPS(self, distance):
        """
        This method provides the value of the kernel function given the euclidean
        distance. The kernel function is the one used for TPS.
        """
        phi = 0.0

        if distance > 0.0:
            phi = (distance**2) * np.log10(distance)
        else:
            phi = 0.0

        return phi

    def interpolateSolidPositionOnFluidMesh(self, FSI_config):
        """
        Applies the one-to-one mapping or the interpolaiton rules from solid to fluid mesh.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

        # --- Interpolate (or map) in parallel the solid interface displacement on the fluid interface ---
        if FSI_config["MATCHING_MESH"] == "NO" and (
            FSI_config["MESH_INTERP_METHOD"] == "RBF"
            or FSI_config["MESH_INTERP_METHOD"] == "TPS"
        ):
            if self.have_MPI:
                gamma_array_DispX = PETSc.Vec().create(self.comm)
                gamma_array_DispY = PETSc.Vec().create(self.comm)
                gamma_array_DispZ = PETSc.Vec().create(self.comm)
                gamma_array_DispX.setType("mpi")
                gamma_array_DispY.setType("mpi")
                gamma_array_DispZ.setType("mpi")
                KSP_solver = PETSc.KSP().create(self.comm)
            else:
                gamma_array_DispX = PETSc.Vec().create()
                gamma_array_DispY = PETSc.Vec().create()
                gamma_array_DispZ = PETSc.Vec().create()
                gamma_array_DispX.setType("seq")
                gamma_array_DispY.setType("seq")
                gamma_array_DispZ.setType("seq")
                KSP_solver = PETSc.KSP().create()
            gamma_array_DispX.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            gamma_array_DispY.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            gamma_array_DispZ.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            gamma_array_DispX.set(0.0)
            gamma_array_DispY.set(0.0)
            gamma_array_DispZ.set(0.0)
            KSP_solver.setType("fgmres")
            KSP_solver.getPC().setType("jacobi")
            KSP_solver.setOperators(self.MappingMatrixA)
            KSP_solver.setFromOptions()
            KSP_solver.setInitialGuessNonzero(True)
            KSP_solver.solve(self.solidInterface_array_DispX, gamma_array_DispX)
            KSP_solver.solve(self.solidInterface_array_DispY, gamma_array_DispY)
            if self.nDim == 3:
                KSP_solver.solve(self.solidInterface_array_DispZ, gamma_array_DispZ)
            self.MappingMatrixB.mult(gamma_array_DispX, self.fluidInterface_array_DispX)
            self.MappingMatrixB.mult(gamma_array_DispY, self.fluidInterface_array_DispY)
            if self.nDim == 3:
                self.MappingMatrixB.mult(
                    gamma_array_DispZ, self.fluidInterface_array_DispZ
                )
            gamma_array_DispX.destroy()
            gamma_array_DispY.destroy()
            gamma_array_DispZ.destroy()
            KSP_solver.destroy()
            del gamma_array_DispX
            del gamma_array_DispY
            del gamma_array_DispZ
            del KSP_solver
        else:
            self.MappingMatrix.mult(
                self.solidInterface_array_DispX, self.fluidInterface_array_DispX
            )
            self.MappingMatrix.mult(
                self.solidInterface_array_DispY, self.fluidInterface_array_DispY
            )
            if self.nDim == 3:
                self.MappingMatrix.mult(
                    self.solidInterface_array_DispZ, self.fluidInterface_array_DispZ
                )

        # --- Checking conservation ---
        WSX = self.solidLoads_array_X.dot(self.solidInterface_array_DispX)
        WSY = self.solidLoads_array_Y.dot(self.solidInterface_array_DispY)
        WSZ = self.solidLoads_array_Z.dot(self.solidInterface_array_DispZ)

        WFX = self.fluidLoads_array_X.dot(self.fluidInterface_array_DispX)
        WFY = self.fluidLoads_array_Y.dot(self.fluidInterface_array_DispY)
        WFZ = self.fluidLoads_array_Z.dot(self.fluidInterface_array_DispZ)

        self.MPIPrint("Checking f/s interface conservation...")
        self.MPIPrint("Solid side (Wx, Wy, Wz) = ({}, {}, {})".format(WSX, WSY, WSZ))
        self.MPIPrint("Fluid side (Wx, Wy, Wz) = ({}, {}, {})".format(WFX, WFY, WFZ))

        # --- Redistribute the interpolated fluid interface according to the partitions that own the fluid interface ---
        # Gather the fluid interface on the master process
        # This is required because PETSc redistributes evenly in the cores, and does not use the same division
        # of SU2, thus we need to redistribute
        if self.have_MPI:
            sendBuff_X = None
            sendBuff_Y = None
            sendBuff_Z = None
            self.fluidInterface_array_DispX_recon = None
            self.fluidInterface_array_DispY_recon = None
            self.fluidInterface_array_DispZ_recon = None

            if myid == self.rootProcess:
                self.fluidInterface_array_DispX_recon = np.zeros(
                    self.nFluidInterfacePhysicalNodes
                )
                self.fluidInterface_array_DispY_recon = np.zeros(
                    self.nFluidInterfacePhysicalNodes
                )
                self.fluidInterface_array_DispZ_recon = np.zeros(
                    self.nFluidInterfacePhysicalNodes
                )

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

            self.comm.Gatherv(
                self.fluidInterface_array_DispX.getArray(),
                [self.fluidInterface_array_DispX_recon, counts, displ, self.MPI.DOUBLE],
                root=self.rootProcess,
            )
            self.comm.Gatherv(
                self.fluidInterface_array_DispY.getArray(),
                [self.fluidInterface_array_DispY_recon, counts, displ, self.MPI.DOUBLE],
                root=self.rootProcess,
            )
            self.comm.Gatherv(
                self.fluidInterface_array_DispZ.getArray(),
                [self.fluidInterface_array_DispZ_recon, counts, displ, self.MPI.DOUBLE],
                root=self.rootProcess,
            )

            # Send the partitioned interface to the right fluid partitions
            if myid == self.rootProcess:
                for iProc in self.fluidInterfaceProcessors:
                    sendBuff_X = np.empty(
                        self.fluidPhysicalInterfaceNodesDistribution[iProc],
                        dtype=np.float64,
                    )
                    sendBuff_Y = np.empty(
                        self.fluidPhysicalInterfaceNodesDistribution[iProc],
                        dtype=np.float64,
                    )
                    sendBuff_Z = np.empty(
                        self.fluidPhysicalInterfaceNodesDistribution[iProc],
                        dtype=np.float64,
                    )
                    globalIndex = self.__getGlobalIndex("fluid", iProc, 0)
                    for iVertex in range(
                        self.fluidPhysicalInterfaceNodesDistribution[iProc]
                    ):
                        sendBuff_X[iVertex] = self.fluidInterface_array_DispX_recon[
                            globalIndex
                        ]
                        sendBuff_Y[iVertex] = self.fluidInterface_array_DispY_recon[
                            globalIndex
                        ]
                        sendBuff_Z[iVertex] = self.fluidInterface_array_DispZ_recon[
                            globalIndex
                        ]
                        globalIndex += 1
                    if iProc == self.rootProcess:
                        self.localFluidInterface_array_DispX = np.copy(sendBuff_X)
                        self.localFluidInterface_array_DispY = np.copy(sendBuff_Y)
                        self.localFluidInterface_array_DispZ = np.copy(sendBuff_Z)
                    else:
                        self.comm.Send(sendBuff_X, dest=iProc, tag=1)
                        self.comm.Send(sendBuff_Y, dest=iProc, tag=2)
                        self.comm.Send(sendBuff_Z, dest=iProc, tag=3)
            if myid in self.fluidInterfaceProcessors:
                if myid != self.rootProcess:
                    self.localFluidInterface_array_DispX = np.empty(
                        self.nLocalFluidInterfacePhysicalNodes, dtype=np.float64
                    )
                    self.localFluidInterface_array_DispY = np.empty(
                        self.nLocalFluidInterfacePhysicalNodes, dtype=np.float64
                    )
                    self.localFluidInterface_array_DispZ = np.empty(
                        self.nLocalFluidInterfacePhysicalNodes, dtype=np.float64
                    )
                    self.comm.Recv(
                        self.localFluidInterface_array_DispX,
                        source=self.rootProcess,
                        tag=1,
                    )
                    self.comm.Recv(
                        self.localFluidInterface_array_DispY,
                        source=self.rootProcess,
                        tag=2,
                    )
                    self.comm.Recv(
                        self.localFluidInterface_array_DispZ,
                        source=self.rootProcess,
                        tag=3,
                    )
            del sendBuff_X
            del sendBuff_Y
            del sendBuff_Z
            self.comm.barrier()
        else:
            self.localFluidInterface_array_DispX = (
                self.fluidInterface_array_DispX.getArray().copy()
            )
            self.localFluidInterface_array_DispY = (
                self.fluidInterface_array_DispY.getArray().copy()
            )
            self.localFluidInterface_array_DispZ = (
                self.fluidInterface_array_DispZ.getArray().copy()
            )

        # Special treatment for the halo nodes on the fluid interface
        self.haloNodesDisplacements = {}
        sendBuff = {}
        if self.have_MPI:
            if myid == self.rootProcess:
                for iProc in self.fluidInterfaceProcessors:
                    sendBuff = {}
                    for key in self.FluidHaloNodeList[
                        iProc
                    ].keys():  # The keys are the SU2 global IDs of the interface nodes
                        globalIndex = self.fluidIndexing[
                            key
                        ]  # These are the interface global IDs, not the SU2 global IDs
                        DispX = self.fluidInterface_array_DispX_recon[globalIndex]
                        DispY = self.fluidInterface_array_DispY_recon[globalIndex]
                        DispZ = self.fluidInterface_array_DispZ_recon[globalIndex]
                        sendBuff[key] = (DispX, DispY, DispZ)
                    if iProc == self.rootProcess:
                        self.haloNodesDisplacements = sendBuff
                    else:
                        self.comm.send(sendBuff, dest=iProc, tag=4)
            if myid in self.fluidInterfaceProcessors:
                if myid != self.rootProcess:
                    self.haloNodesDisplacements = self.comm.recv(
                        source=self.rootProcess, tag=4
                    )
            self.comm.barrier()
            del self.fluidInterface_array_DispX_recon
            del self.fluidInterface_array_DispY_recon
            del self.fluidInterface_array_DispZ_recon
        del sendBuff

    def interpolateFluidLoadsOnSolidMesh(self, FSI_config):
        """
        Applies the one-to-one mapping or the interpolaiton rules from fluid to solid mesh.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
            MPIsize = self.comm.Get_size()
        else:
            myid = 0
            MPIsize = 1

        # --- Interpolate (or map) in parallel the fluid interface loads on the solid interface ---
        # self.MappingMatrix.transpose()
        if FSI_config["MATCHING_MESH"] == "NO" and (
            FSI_config["MESH_INTERP_METHOD"] == "RBF"
            or FSI_config["MESH_INTERP_METHOD"] == "TPS"
        ):
            if self.have_MPI:
                gamma_array_LoadX = PETSc.Vec().create(self.comm)
                gamma_array_LoadY = PETSc.Vec().create(self.comm)
                gamma_array_LoadZ = PETSc.Vec().create(self.comm)
                gamma_array_LoadX.setType("mpi")
                gamma_array_LoadY.setType("mpi")
                gamma_array_LoadZ.setType("mpi")
                KSP_solver = PETSc.KSP().create(self.comm)
            else:
                gamma_array_LoadX = PETSc.Vec().create()
                gamma_array_LoadY = PETSc.Vec().create()
                gamma_array_LoadZ = PETSc.Vec().create()
                gamma_array_LoadX.setType("seq")
                gamma_array_LoadY.setType("seq")
                gamma_array_LoadZ.setType("seq")
                KSP_solver = PETSc.KSP().create()
            gamma_array_LoadX.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            gamma_array_LoadY.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            gamma_array_LoadZ.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            gamma_array_LoadX.set(0.0)
            gamma_array_LoadY.set(0.0)
            gamma_array_LoadZ.set(0.0)
            KSP_solver.setType("fgmres")
            KSP_solver.getPC().setType("jacobi")
            KSP_solver.setOperators(self.MappingMatrixA_T)
            KSP_solver.setFromOptions()
            self.MappingMatrixB_T.mult(self.fluidLoads_array_X, gamma_array_LoadX)
            self.MappingMatrixB_T.mult(self.fluidLoads_array_Y, gamma_array_LoadY)
            if self.nDim == 3:
                self.MappingMatrixB_T.mult(self.fluidLoads_array_Z, gamma_array_LoadZ)
            KSP_solver.solve(gamma_array_LoadX, self.solidLoads_array_X)
            KSP_solver.solve(gamma_array_LoadY, self.solidLoads_array_Y)
            if self.nDim == 3:
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
            if self.nDim == 3:
                self.MappingMatrix_T.mult(
                    self.fluidLoads_array_Z, self.solidLoads_array_Z
                )

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
                self.solidLoads_array_X_recon = np.zeros(
                    self.nSolidInterfacePhysicalNodes + self.d_RBF
                )
                self.solidLoads_array_Y_recon = np.zeros(
                    self.nSolidInterfacePhysicalNodes + self.d_RBF
                )
                self.solidLoads_array_Z_recon = np.zeros(
                    self.nSolidInterfacePhysicalNodes + self.d_RBF
                )
            myNumberOfNodes = self.solidLoads_array_X.getArray().shape[0]
            sendBuffNumber = np.array([myNumberOfNodes], dtype=int)
            rcvBuffNumber = np.zeros(MPIsize, dtype=int)
            self.comm.Allgather(sendBuffNumber, rcvBuffNumber)

            counts = tuple(rcvBuffNumber)
            displ = np.zeros(MPIsize, dtype=int)
            for ii in range(rcvBuffNumber.shape[0]):
                displ[ii] = rcvBuffNumber[0:ii].sum()
            displ = tuple(displ)

            del sendBuffNumber, rcvBuffNumber

            self.comm.Gatherv(
                self.solidLoads_array_X.getArray(),
                [self.solidLoads_array_X_recon, counts, displ, self.MPI.DOUBLE],
                root=self.rootProcess,
            )
            self.comm.Gatherv(
                self.solidLoads_array_Y.getArray(),
                [self.solidLoads_array_Y_recon, counts, displ, self.MPI.DOUBLE],
                root=self.rootProcess,
            )
            self.comm.Gatherv(
                self.solidLoads_array_Z.getArray(),
                [self.solidLoads_array_Z_recon, counts, displ, self.MPI.DOUBLE],
                root=self.rootProcess,
            )

            # Send the partitioned loads to the right solid partitions
            if myid == self.rootProcess:
                for iProc in self.solidInterfaceProcessors:
                    sendBuff_X = np.empty(
                        self.solidPhysicalInterfaceNodesDistribution[iProc],
                        dtype=np.float64,
                    )
                    sendBuff_Y = np.empty(
                        self.solidPhysicalInterfaceNodesDistribution[iProc],
                        dtype=np.float64,
                    )
                    sendBuff_Z = np.empty(
                        self.solidPhysicalInterfaceNodesDistribution[iProc],
                        dtype=np.float64,
                    )
                    globalIndex = self.__getGlobalIndex("solid", iProc, 0)
                    for iVertex in range(
                        self.solidPhysicalInterfaceNodesDistribution[iProc]
                    ):
                        sendBuff_X[iVertex] = self.solidLoads_array_X_recon[globalIndex]
                        sendBuff_Y[iVertex] = self.solidLoads_array_Y_recon[globalIndex]
                        sendBuff_Z[iVertex] = self.solidLoads_array_Z_recon[globalIndex]
                        globalIndex += 1
                    if iProc != myid:
                        self.comm.Send(sendBuff_X, dest=iProc, tag=1)
                        self.comm.Send(sendBuff_Y, dest=iProc, tag=2)
                        self.comm.Send(sendBuff_Z, dest=iProc, tag=3)
                    else:
                        self.localSolidLoads_array_X = np.copy(sendBuff_X)
                        self.localSolidLoads_array_Y = np.copy(sendBuff_Y)
                        self.localSolidLoads_array_Z = np.copy(sendBuff_Z)
            if myid in self.solidInterfaceProcessors:
                if myid != self.rootProcess:
                    self.localSolidLoads_array_X = np.empty(
                        self.nLocalSolidInterfacePhysicalNodes, dtype=np.float64
                    )
                    self.localSolidLoads_array_Y = np.empty(
                        self.nLocalSolidInterfacePhysicalNodes, dtype=np.float64
                    )
                    self.localSolidLoads_array_Z = np.empty(
                        self.nLocalSolidInterfacePhysicalNodes, dtype=np.float64
                    )
                    self.comm.Recv(
                        self.localSolidLoads_array_X, source=self.rootProcess, tag=1
                    )
                    self.comm.Recv(
                        self.localSolidLoads_array_Y, source=self.rootProcess, tag=2
                    )
                    self.comm.Recv(
                        self.localSolidLoads_array_Z, source=self.rootProcess, tag=3
                    )
            del sendBuff_X
            del sendBuff_Y
            del sendBuff_Z
            self.comm.barrier()
        else:
            self.localSolidLoads_array_X = self.solidLoads_array_X.getArray().copy()
            self.localSolidLoads_array_Y = self.solidLoads_array_Y.getArray().copy()
            self.localSolidLoads_array_Z = self.solidLoads_array_Z.getArray().copy()

        # Special treatment for the halo nodes on the solid interface
        self.haloNodesLoads = {}
        sendBuff = {}
        if self.have_MPI:
            if myid == self.rootProcess:
                for iProc in self.solidInterfaceProcessors:
                    sendBuff = {}
                    for key in self.SolidHaloNodeList[iProc].keys():
                        globalIndex = self.solidIndexing[key]
                        DispX = self.solidLoads_array_X_recon[globalIndex]
                        DispY = self.solidLoads_array_Y_recon[globalIndex]
                        DispZ = self.solidLoads_array_Z_recon[globalIndex]
                        sendBuff[key] = (DispX, DispY, DispZ)
                    if iProc == self.rootProcess:
                        self.haloNodesLoads = sendBuff
                    else:
                        self.comm.send(sendBuff, dest=iProc, tag=4)
            if myid in self.solidInterfaceProcessors:
                if myid != self.rootProcess:
                    self.haloNodesLoads = self.comm.recv(source=self.rootProcess, tag=4)
            self.comm.barrier()
            del self.solidLoads_array_X_recon
            del self.solidLoads_array_Y_recon
            del self.solidLoads_array_Z_recon
        del sendBuff

    def getSolidInterfaceDisplacement(self, SolidSolver):
        """
        Gets the current solid interface position from the solid solver.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        # --- Get the solid interface position from the solid solver and directly fill the corresponding PETSc vector ---
        GlobalIndex = int()
        localIndex = 0
        for iVertex in range(self.nLocalSolidInterfaceNodes):
            GlobalIndex = SolidSolver.getVertexGlobalIndex(
                self.solidInterfaceIdentifier, iVertex
            )
            if GlobalIndex not in self.SolidHaloNodeList[myid].keys():
                newDispx, newDispy, newDispz = SolidSolver.getInterfaceNodeDisp(
                    self.solidInterfaceIdentifier, iVertex
                )
                iGlobalVertex = self.__getGlobalIndex("solid", myid, localIndex)
                self.solidInterface_array_DispX.setValues([iGlobalVertex], newDispx)
                self.solidInterface_array_DispY.setValues([iGlobalVertex], newDispy)
                self.solidInterface_array_DispZ.setValues([iGlobalVertex], newDispz)
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
        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        GlobalIndex = int()
        localIndex = 0

        # --- Get the fluid interface loads from the fluid solver and directly fill the corresponding PETSc vector ---
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            GlobalIndex = FluidSolver.GetNodeGlobalIndex(
                FluidSolver.GetMarkerNode(self.fluidInterfaceIdentifier, iVertex)
            )
            if GlobalIndex not in self.FluidHaloNodeList[myid].keys():
                load = FluidSolver.GetMarkerFlowLoad(
                    self.fluidInterfaceIdentifier, iVertex
                )
                iGlobalVertex = self.__getGlobalIndex("fluid", myid, localIndex)
                self.fluidLoads_array_X.setValues([iGlobalVertex], load[0])
                self.fluidLoads_array_Y.setValues([iGlobalVertex], load[1])
                self.fluidLoads_array_Z.setValues(
                    [iGlobalVertex], load[2] if len(load) == 3 else 0.0
                )
                localIndex += 1

        self.fluidLoads_array_X.assemblyBegin()
        self.fluidLoads_array_X.assemblyEnd()
        self.fluidLoads_array_Y.assemblyBegin()
        self.fluidLoads_array_Y.assemblyEnd()
        self.fluidLoads_array_Z.assemblyBegin()
        self.fluidLoads_array_Z.assemblyEnd()

    def setFluidInterfaceVarCoord(self, FluidSolver):
        """
        Communicate the change of coordinates of the fluid interface to the fluid solver.
        Prepare the fluid solver for mesh deformation.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        # --- Send the new fluid interface position to the fluid solver (on each partition, halo nodes included) ---
        localIndex = 0
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            GlobalIndex = FluidSolver.GetNodeGlobalIndex(
                FluidSolver.GetMarkerNode(self.fluidInterfaceIdentifier, iVertex)
            )
            if GlobalIndex in self.FluidHaloNodeList[myid].keys():
                DispX, DispY, DispZ = self.haloNodesDisplacements[GlobalIndex]
                FluidSolver.SetMarkerCustomDisplacement(
                    self.fluidInterfaceIdentifier,
                    int(iVertex),
                    np.array([DispX, DispY, DispZ]),
                )
            else:
                DispX = self.localFluidInterface_array_DispX[localIndex]
                DispY = self.localFluidInterface_array_DispY[localIndex]
                DispZ = self.localFluidInterface_array_DispZ[localIndex]
                FluidSolver.SetMarkerCustomDisplacement(
                    self.fluidInterfaceIdentifier,
                    int(iVertex),
                    np.array([DispX, DispY, DispZ]),
                )
                localIndex += 1

    def setSolidInterfaceLoads(self, SolidSolver, FSI_config):
        """
        Communicates the new solid interface loads to the solid solver.
        Calculates the new resultant forces (lift, drag, ...).
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        FX = np.array(0.0, dtype=np.float64)
        FY = np.array(0.0, dtype=np.float64)  # solid-side resultant forces
        FZ = np.array(0.0, dtype=np.float64)
        FXSendBuff = np.array(0.0, dtype=np.float64)
        FYSendBuff = np.array(0.0, dtype=np.float64)
        FZSendBuff = np.array(0.0, dtype=np.float64)
        FFX = 0.0  # fluid-side resultant forces
        FFY = 0.0
        FFZ = 0.0

        # --- Check for total force conservation after interpolation
        FFX = self.fluidLoads_array_X.sum()
        FFY = self.fluidLoads_array_Y.sum()
        FFZ = self.fluidLoads_array_Z.sum()

        for iVertex in range(self.nLocalSolidInterfacePhysicalNodes):
            FXSendBuff += self.localSolidLoads_array_X[iVertex]
            FYSendBuff += self.localSolidLoads_array_Y[iVertex]
            FZSendBuff += self.localSolidLoads_array_Z[iVertex]

        if self.have_MPI:
            self.comm.Allreduce(FXSendBuff, FX, op=self.MPI.SUM)
            self.comm.Allreduce(FYSendBuff, FY, op=self.MPI.SUM)
            self.comm.Allreduce(FZSendBuff, FZ, op=self.MPI.SUM)
        else:
            FX = np.copy(FXSendBuff)
            FY = np.copy(FYSendBuff)
            FZ = np.copy(FZSendBuff)

        del FXSendBuff
        del FYSendBuff
        del FZSendBuff

        self.MPIPrint("Checking f/s interface total force...")
        self.MPIPrint("Solid side (Fx, Fy, Fz) = ({}, {}, {})".format(FX, FY, FZ))
        self.MPIPrint("Fluid side (Fx, Fy, Fz) = ({}, {}, {})".format(FFX, FFY, FFZ))

        # --- Send the new solid interface loads to the solid solver (on each partition, halo nodes included) ---
        GlobalIndex = int()
        localIndex = 0
        for iVertex in range(self.nLocalSolidInterfaceNodes):
            GlobalIndex = SolidSolver.getVertexGlobalIndex(
                self.solidInterfaceIdentifier, iVertex
            )
            if GlobalIndex in self.SolidHaloNodeList[myid].keys():
                pass  # TODO here, when the solid solver will run in parallel, we will need to pass the halo loads
            else:
                Fx = self.localSolidLoads_array_X[localIndex]
                Fy = self.localSolidLoads_array_Y[localIndex]
                Fz = self.localSolidLoads_array_Z[localIndex]
                SolidSolver.applyload(iVertex, Fx, Fy, Fz)
                localIndex += 1

    def computeSolidInterfaceResidual(self, SolidSolver):
        """
        Computes the solid interface FSI displacement residual.
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        normInterfaceResidualSquare = 0.0

        # --- Create and fill the PETSc vector for the predicted solid interface position (predicted by the solid computation) ---
        if self.have_MPI:
            predDisp_array_X = PETSc.Vec().create(self.comm)
            predDisp_array_X.setType("mpi")
            predDisp_array_Y = PETSc.Vec().create(self.comm)
            predDisp_array_Y.setType("mpi")
            predDisp_array_Z = PETSc.Vec().create(self.comm)
            predDisp_array_Z.setType("mpi")
        else:
            predDisp_array_X = PETSc.Vec().create()
            predDisp_array_X.setType("seq")
            predDisp_array_Y = PETSc.Vec().create()
            predDisp_array_Y.setType("seq")
            predDisp_array_Z = PETSc.Vec().create()
            predDisp_array_Z.setType("seq")
        predDisp_array_X.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        predDisp_array_Y.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        predDisp_array_Z.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        predDisp_array_X.set(0.0)
        predDisp_array_Y.set(0.0)
        predDisp_array_Z.set(0.0)

        for iVertex in range(self.nLocalSolidInterfaceNodes):
            predDispx, predDispy, predDispz = SolidSolver.getInterfaceNodeDisp(
                self.solidInterfaceIdentifier, iVertex
            )
            iGlobalVertex = self.__getGlobalIndex("solid", myid, iVertex)
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
        self.solidInterfaceResidual_array_X = (
            predDisp_array_X - self.solidInterface_array_DispX
        )
        self.solidInterfaceResidual_array_Y = (
            predDisp_array_Y - self.solidInterface_array_DispY
        )
        self.solidInterfaceResidual_array_Z = (
            predDisp_array_Z - self.solidInterface_array_DispZ
        )

        normInterfaceResidual_X = self.solidInterfaceResidual_array_X.norm()
        normInterfaceResidual_Y = self.solidInterfaceResidual_array_Y.norm()
        normInterfaceResidual_Z = self.solidInterfaceResidual_array_Z.norm()

        normInterfaceResidualSquare = (
            normInterfaceResidual_X**2
            + normInterfaceResidual_Y**2
            + normInterfaceResidual_Z**2
        )

        predDisp_array_X.destroy()
        predDisp_array_Y.destroy()
        predDisp_array_Z.destroy()
        del predDisp_array_X
        del predDisp_array_Y
        del predDisp_array_Z

        return sqrt(normInterfaceResidualSquare)

    def relaxSolidPosition(self, FSI_config):
        """
        Apply solid displacement under-relaxation.
        """
        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        # --- Set the Aitken coefficient for the relaxation ---
        if FSI_config["AITKEN_RELAX"] == "STATIC":
            self.aitkenParam = FSI_config["AITKEN_PARAM"]
        elif FSI_config["AITKEN_RELAX"] == "DYNAMIC":
            self.setAitkenCoefficient(FSI_config)
        else:
            self.aitkenParam = 1.0

        self.MPIPrint(
            "Aitken under-relaxation step with parameter {}".format(self.aitkenParam)
        )

        # --- Relax the solid interface position ---
        self.solidInterface_array_DispX += (
            self.aitkenParam * self.solidInterfaceResidual_array_X
        )
        self.solidInterface_array_DispY += (
            self.aitkenParam * self.solidInterfaceResidual_array_Y
        )
        self.solidInterface_array_DispZ += (
            self.aitkenParam * self.solidInterfaceResidual_array_Z
        )

    def setAitkenCoefficient(self, FSI_config):
        """
        Computes the Aitken coefficients for solid displacement under-relaxation.
        """

        deltaResNormSquare = 0.0
        prodScalRes = 0.0

        # --- Create the PETSc vector for the difference between the residuals (current and previous FSI iter) ---
        if self.FSIIter == 0:
            self.aitkenParam = max(FSI_config["AITKEN_PARAM"], self.aitkenParam)
        else:
            if self.have_MPI:
                deltaResx_array_X = PETSc.Vec().create(self.comm)
                deltaResx_array_X.setType("mpi")
                deltaResx_array_Y = PETSc.Vec().create(self.comm)
                deltaResx_array_Y.setType("mpi")
                deltaResx_array_Z = PETSc.Vec().create(self.comm)
                deltaResx_array_Z.setType("mpi")
            else:
                deltaResx_array_X = PETSc.Vec().create()
                deltaResx_array_X.setType("seq")
                deltaResx_array_Y = PETSc.Vec().create()
                deltaResx_array_Y.setType("seq")
                deltaResx_array_Z = PETSc.Vec().create()
                deltaResx_array_Z.setType("seq")
            deltaResx_array_X.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            deltaResx_array_X.set(0.0)
            deltaResx_array_Y.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            deltaResx_array_Y.set(0.0)
            deltaResx_array_Z.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
            deltaResx_array_Z.set(0.0)

            # --- Compute the dynamic Aitken coefficient ---
            deltaResx_array_X = (
                self.solidInterfaceResidual_array_X
                - self.solidInterfaceResidualnM1_array_X
            )
            deltaResx_array_Y = (
                self.solidInterfaceResidual_array_Y
                - self.solidInterfaceResidualnM1_array_Y
            )
            deltaResx_array_Z = (
                self.solidInterfaceResidual_array_Z
                - self.solidInterfaceResidualnM1_array_Z
            )

            prodScalRes_X = deltaResx_array_X.dot(
                self.solidInterfaceResidualnM1_array_X
            )
            prodScalRes_Y = deltaResx_array_Y.dot(
                self.solidInterfaceResidualnM1_array_Y
            )
            prodScalRes_Z = deltaResx_array_Z.dot(
                self.solidInterfaceResidualnM1_array_Z
            )
            prodScalRes = prodScalRes_X + prodScalRes_Y + prodScalRes_Z

            deltaResNormSquare_X = (deltaResx_array_X.norm()) ** 2
            deltaResNormSquare_Y = (deltaResx_array_Y.norm()) ** 2
            deltaResNormSquare_Z = (deltaResx_array_Z.norm()) ** 2
            deltaResNormSquare = (
                deltaResNormSquare_X + deltaResNormSquare_Y + deltaResNormSquare_Z
            )

            self.aitkenParam *= -prodScalRes / deltaResNormSquare

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

    def displacementPredictor(self, FSI_config, SolidSolver, deltaT):
        """
        Calculates a prediciton for the solid interface position for the next time step.
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        if FSI_config["DISP_PRED"] == "FIRST_ORDER":
            self.MPIPrint("First order predictor")
            alpha_0 = 1.0
            alpha_1 = 0.0
        elif FSI_config["DISP_PRED"] == "SECOND_ORDER":
            self.MPIPrint("Second order predictor")
            alpha_0 = 1.0
            alpha_1 = 0.5
        else:
            self.MPIPrint("No predictor")
            alpha_0 = 0.0
            alpha_1 = 0.0

        # --- Create the PETSc vectors to store the solid interface velocity ---
        if self.have_MPI:
            Vel_array_X = PETSc.Vec().create(self.comm)
            Vel_array_X.setType("mpi")
            Vel_array_Y = PETSc.Vec().create(self.comm)
            Vel_array_Y.setType("mpi")
            Vel_array_Z = PETSc.Vec().create(self.comm)
            Vel_array_Z.setType("mpi")
            VelnM1_array_X = PETSc.Vec().create(self.comm)
            VelnM1_array_X.setType("mpi")
            VelnM1_array_Y = PETSc.Vec().create(self.comm)
            VelnM1_array_Y.setType("mpi")
            VelnM1_array_Z = PETSc.Vec().create(self.comm)
            VelnM1_array_Z.setType("mpi")
        else:
            Vel_array_X = PETSc.Vec().create()
            Vel_array_X.setType("seq")
            Vel_array_Y = PETSc.Vec().create()
            Vel_array_Y.setType("seq")
            Vel_array_Z = PETSc.Vec().create()
            Vel_array_Z.setType("seq")
            VelnM1_array_X = PETSc.Vec().create()
            VelnM1_array_X.setType("seq")
            VelnM1_array_Y = PETSc.Vec().create()
            VelnM1_array_Y.setType("seq")
            VelnM1_array_Z = PETSc.Vec().create()
            VelnM1_array_Z.setType("seq")
        Vel_array_X.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        Vel_array_Y.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        Vel_array_Z.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        Vel_array_X.set(0.0)
        Vel_array_Y.set(0.0)
        Vel_array_Z.set(0.0)
        VelnM1_array_X.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        VelnM1_array_Y.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        VelnM1_array_Z.setSizes(self.nSolidInterfacePhysicalNodes + self.d_RBF)
        VelnM1_array_X.set(0.0)
        VelnM1_array_Y.set(0.0)
        VelnM1_array_Z.set(0.0)

        # --- Fill the PETSc vectors ---
        GlobalIndex = int()
        localIndex = 0
        for iVertex in range(self.nLocalSolidInterfaceNodes):
            GlobalIndex = SolidSolver.getVertexGlobalIndex(
                self.solidInterfaceIdentifier, iVertex
            )
            if GlobalIndex not in self.SolidHaloNodeList[myid].keys():
                iGlobalVertex = self.__getGlobalIndex("solid", myid, localIndex)
                velx, vely, velz = SolidSolver.getInterfaceNodeVel(
                    self.solidInterfaceIdentifier, iVertex
                )
                velxNm1, velyNm1, velzNm1 = SolidSolver.getInterfaceNodeVelNm1(
                    self.solidInterfaceIdentifier, iVertex
                )
                Vel_array_X.setValues([iGlobalVertex], velx)
                Vel_array_Y.setValues([iGlobalVertex], vely)
                Vel_array_Z.setValues([iGlobalVertex], velz)
                VelnM1_array_X.setValues([iGlobalVertex], velxNm1)
                VelnM1_array_Y.setValues([iGlobalVertex], velyNm1)
                VelnM1_array_Z.setValues([iGlobalVertex], velzNm1)
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
        self.solidInterface_array_DispX += (
            alpha_0 * deltaT * Vel_array_X
            + alpha_1 * deltaT * (Vel_array_X - VelnM1_array_X)
        )
        self.solidInterface_array_DispY += (
            alpha_0 * deltaT * Vel_array_Y
            + alpha_1 * deltaT * (Vel_array_Y - VelnM1_array_Y)
        )
        self.solidInterface_array_DispZ += (
            alpha_0 * deltaT * Vel_array_Z
            + alpha_1 * deltaT * (Vel_array_Z - VelnM1_array_Z)
        )

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

        if self.have_MPI:
            myid = self.comm.Get_rank()
        else:
            myid = 0

        if myid == self.rootProcess:
            if self.unsteady:
                if TimeIter == 0:
                    histFile = open("FSIhistory.dat", "w")
                    histFile.write("TimeIter\tTime\tFSIRes\tFSINbIter\n")
                else:
                    histFile = open("FSIhistory.dat", "a")
                if FSIConv:
                    histFile.write(
                        str(TimeIter)
                        + "\t"
                        + str(time)
                        + "\t"
                        + str(varCoordNorm)
                        + "\t"
                        + str(self.FSIIter + 1)
                        + "\n"
                    )
                else:
                    histFile.write(
                        str(TimeIter)
                        + "\t"
                        + str(time)
                        + "\t"
                        + str(varCoordNorm)
                        + "\t"
                        + str(self.FSIIter)
                        + "\n"
                    )
                histFile.close()
            else:
                if self.FSIIter == 0:
                    histFile = open("FSIhistory.dat", "w")
                    histFile.write("FSI Iter\tFSIRes\n")
                else:
                    histFile = open("FSIhistory.dat", "a")
                histFile.write(str(self.FSIIter) + "\t" + str(varCoordNorm) + "\n")
                histFile.close()

        self.MPIBarrier()

    def __getGlobalIndex(self, physics, iProc, iLocalVertex):
        """
        Calculate the global indexing of interface nodes accross all the partitions. This does not include halo nodes.
        This is needed because the global index of the fluid solver takes into account all the nodes, thus also those
        in the volume mesh, not only on the interface. Here, we compute the global index of all the nodes on the
        interface.
        """

        if physics == "fluid":
            globalStartIndex = self.fluidGlobalIndexRange[iProc][iProc][0]
        elif physics == "solid":
            globalStartIndex = self.solidGlobalIndexRange[iProc][iProc][0]

        globalIndex = globalStartIndex + iLocalVertex

        return globalIndex

    def UnsteadyFSI(self, FSI_config, FluidSolver, SolidSolver):
        """
        Run the unsteady FSI computation by synchronizing the fluid and solid solvers.
        F/s interface data are exchanged through interface mapping and interpolation (if non mathcing meshes).
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
            numberPart = self.comm.Get_size()
        else:
            myid = 0
            numberPart = 1

        # --- Set some general variables for the unsteady computation --- #
        deltaT = FSI_config["UNST_TIMESTEP"]  # physical time step
        totTime = FSI_config["UNST_TIME"]  # physical simulation time
        NbFSIIterMax = FSI_config[
            "NB_FSI_ITER"
        ]  # maximum number of FSI iteration (for each time step)
        FSITolerance = FSI_config["FSI_TOLERANCE"]  # f/s interface tolerance
        TimeIterTreshold = FSI_config[
            "TIME_TRESHOLD"
        ]  # time iteration from which we allow the solid to deform
        self.MPIPrint(
            "The FSI coupling will start after {} iterations".format(TimeIterTreshold)
        )

        if FSI_config["RESTART_SOL"] == "YES":
            NbTimeIter = ((totTime) / deltaT) - 1
            time = (FSI_config["RESTART_ITER"]) * deltaT
            TimeIter = FSI_config["RESTART_ITER"]
        else:
            NbTimeIter = (totTime / deltaT) - 1  # number of time iterations
            time = 0.0  # initial time
            TimeIter = 0  # initial time iteration

        NbTimeIter = int(NbTimeIter)  # be sure that NbTimeIter is an integer

        varCoordNorm = 0.0  # FSI residual
        FSIConv = False  # FSI convergence flag

        self.MPIPrint("\n**********************************")
        self.MPIPrint("* Begin unsteady FSI computation *")
        self.MPIPrint("**********************************\n")

        # --- Initialize the coupled solution --- #
        # If restart
        if FSI_config["RESTART_SOL"] == "YES":
            self.getSolidInterfaceDisplacement(SolidSolver)
            self.displacementPredictor(FSI_config, SolidSolver, deltaT)
            if myid in self.solidSolverProcessors:
                SolidSolver.updateSolution()
        # If no restart
        else:
            self.MPIPrint("Setting FSI initial conditions")
            if myid in self.solidSolverProcessors:
                SolidSolver.setInitialDisplacements()
            self.getSolidInterfaceDisplacement(SolidSolver)
            self.interpolateSolidPositionOnFluidMesh(FSI_config)
            self.setFluidInterfaceVarCoord(FluidSolver)
            self.MPIPrint(
                "\nPerforming static mesh deformation (ALE) of initial mesh...\n"
            )
            if myid in self.fluidSolverProcessors:
                FluidSolver.SetInitialMesh()  # if there is an initial deformation in the solid, it has to be communicated to the fluid solver
            self.MPIPrint("\nFSI initial conditions are set")
            self.MPIPrint("Beginning time integration\n")

        # --- External temporal loop --- #
        while TimeIter <= NbTimeIter:

            if TimeIter > TimeIterTreshold:
                NbFSIIter = NbFSIIterMax
                self.MPIPrint("\n")
                self.MPIPrint(
                    " Enter Block Gauss Seidel (BGS) method for strong coupling FSI on time iteration {} ".format(
                        TimeIter
                    ).center(
                        80, "*"
                    )
                )
            else:
                NbFSIIter = 1

            self.FSIIter = 0
            FSIConv = False

            # --- Internal FSI loop --- #
            while self.FSIIter <= (NbFSIIter - 1):

                self.MPIPrint(
                    "\n>>>> Time iteration {} / FSI iteration {} <<<<".format(
                        TimeIter, self.FSIIter
                    )
                )

                # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
                self.interpolateSolidPositionOnFluidMesh(FSI_config)
                self.MPIPrint("\nPerforming dynamic mesh deformation (ALE)...\n")
                self.setFluidInterfaceVarCoord(FluidSolver)
                if myid in self.fluidSolverProcessors:
                    if self.FSIIter == 0:
                        FluidSolver.Preprocess(
                            TimeIter
                        )  # set some parameters before temporal fluid iteration and dynamic mesh update
                    else:
                        FluidSolver.DynamicMeshUpdate(TimeIter)
                # --- Fluid solver call for FSI subiteration --- #
                self.MPIPrint(
                    "\nLaunching fluid solver for one single dual-time iteration..."
                )
                self.MPIBarrier()
                if myid in self.fluidSolverProcessors:
                    FluidSolver.Run()
                    self.MPIBarrier()
                    FluidSolver.Postprocess()
                    self.MPIBarrier()

                # --- Surface fluid loads interpolation and communication --- #
                if TimeIter > TimeIterTreshold:
                    if not self.ImposedMotion:
                        self.MPIPrint("\nProcessing interface fluid loads...\n")
                        self.MPIBarrier()
                        self.getFluidInterfaceNodalForce(FSI_config, FluidSolver)
                        self.MPIBarrier()
                        self.interpolateFluidLoadsOnSolidMesh(FSI_config)
                        self.setSolidInterfaceLoads(SolidSolver, FSI_config)

                    # --- Solid solver call for FSI subiteration --- #
                    self.MPIPrint(
                        "\nLaunching solid solver for a single time iteration...\n"
                    )
                    if myid in self.solidSolverProcessors:
                        SolidSolver.run(time)

                    # --- Compute and monitor the FSI residual --- #
                    varCoordNorm = self.computeSolidInterfaceResidual(SolidSolver)
                    self.MPIPrint("\nFSI displacement norm : {}\n".format(varCoordNorm))
                    if varCoordNorm < FSITolerance:
                        FSIConv = True
                        break

                    # --- Relaxe the solid position --- #
                    self.MPIPrint("\nProcessing interface displacements...\n")
                    self.relaxSolidPosition(FSI_config)

                self.FSIIter += 1
            # --- End OF FSI loop --- #

            self.MPIBarrier()

            # --- Update the FSI history file --- #
            if TimeIter > TimeIterTreshold:
                self.MPIPrint("\nBGS is converged (strong coupling)")
            self.writeFSIHistory(TimeIter, time, varCoordNorm, FSIConv)

            # --- Update, monitor and output the fluid solution before the next time step  ---#
            if myid in self.fluidSolverProcessors:
                FluidSolver.Update()
                FluidSolver.Monitor(TimeIter)
                FluidSolver.Output(TimeIter)

            if TimeIter >= TimeIterTreshold:
                if myid in self.solidSolverProcessors:
                    # --- Output the solid solution before thr next time step --- #
                    SolidSolver.writeSolution(time, TimeIter, self.FSIIter)

            if TimeIter > TimeIterTreshold:
                # --- Displacement predictor for the next time step and update of the solid solution --- #
                self.MPIPrint("\nSolid displacement prediction for next time step")
                self.displacementPredictor(FSI_config, SolidSolver, deltaT)
                if myid in self.solidSolverProcessors:
                    SolidSolver.updateSolution()

            TimeIter += 1
            time += deltaT
        # --- End of the temporal loop --- #

        self.MPIBarrier()

        self.MPIPrint("\n*************************")
        self.MPIPrint("*  End FSI computation  *")
        self.MPIPrint("*************************\n")

    def SteadyFSI(self, FSI_config, FluidSolver, SolidSolver):
        """
        Runs the steady FSI computation by synchronizing the fluid and solid solver with data exchange at the f/s interface.
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
            numberPart = self.comm.Get_size()
        else:
            myid = 0
            numberPart = 1

        # --- Set some general variables for the steady computation --- #
        NbFSIIterMax = FSI_config[
            "NB_FSI_ITER"
        ]  # maximum number of FSI iteration (for each time step)
        FSITolerance = FSI_config["FSI_TOLERANCE"]  # f/s interface tolerance
        varCoordNorm = 0.0

        self.MPIPrint("\n********************************")
        self.MPIPrint("* Begin steady FSI computation *")
        self.MPIPrint("********************************\n")
        self.MPIPrint("\n")
        self.MPIPrint(
            " Enter Block Gauss Seidel (BGS) method for strong coupling FSI ".center(
                80, "*"
            )
        )

        self.MPIPrint("Setting initial deformed mesh")
        if myid in self.solidSolverProcessors:
            SolidSolver.setInitialDisplacements()
        self.getSolidInterfaceDisplacement(SolidSolver)
        self.interpolateSolidPositionOnFluidMesh(FSI_config)
        self.setFluidInterfaceVarCoord(FluidSolver)
        self.MPIPrint("\nFSI initial conditions are set")

        # --- External FSI loop --- #
        self.FSIIter = 0
        while self.FSIIter < NbFSIIterMax:
            self.MPIPrint("\n>>>> FSI iteration {} <<<<".format(self.FSIIter))
            self.MPIPrint("\nLaunching fluid solver for a steady computation...")
            # --- Fluid solver call for FSI subiteration ---#

            if myid in self.fluidSolverProcessors:
                # The mesh will be deformed in the context of the preprocessor, there is no need to set the initial
                # mesh pushing back the solution to avoid spurious velocities, as the velocity is not computed at all
                self.MPIPrint("\nPerforming static mesh deformation...\n")
                FluidSolver.Preprocess(
                    0
                )  # This will attempt to always set the initial condition, but there is a flag on the unsteady computation that will avoid it
                FluidSolver.Run()
                FluidSolver.Postprocess()
                FluidSolver.Monitor(
                    0
                )  # This is actually not needed, it only saves the fact that the fluid solver converged innerly or reached max iterations
                FluidSolver.Output(0)

            # --- Surface fluid loads interpolation and communication ---#
            if not self.ImposedMotion:
                self.MPIPrint("\nProcessing interface fluid loads...\n")
                self.MPIBarrier()
                self.getFluidInterfaceNodalForce(FSI_config, FluidSolver)
                self.MPIBarrier()
                self.interpolateFluidLoadsOnSolidMesh(FSI_config)
                self.setSolidInterfaceLoads(SolidSolver, FSI_config)
                # --- Solid solver call for FSI subiteration --- #
                self.MPIPrint("\nLaunching solid solver for a static computation...\n")
                if myid in self.solidSolverProcessors:
                    SolidSolver.run(0.0)
                    SolidSolver.writeSolution(0.0, 0, self.FSIIter)

            # --- Compute and monitor the FSI residual --- #
            varCoordNorm = self.computeSolidInterfaceResidual(SolidSolver)
            self.MPIPrint("\nFSI displacement norm : {}\n".format(varCoordNorm))
            self.writeFSIHistory(0, 0.0, varCoordNorm, False)
            if varCoordNorm < FSITolerance:
                break

            # --- Relaxe the solid displacement and update the solid solution --- #
            self.MPIPrint("\nProcessing interface displacements...\n")
            self.relaxSolidPosition(FSI_config)
            if myid in self.solidSolverProcessors:
                SolidSolver.updateSolution()

            # --- Mesh morphing step (displacement interpolation, displacements communication, and mesh morpher call) --- #
            self.interpolateSolidPositionOnFluidMesh(FSI_config)
            self.setFluidInterfaceVarCoord(FluidSolver)
            self.FSIIter += 1

        self.MPIBarrier()

        self.MPIPrint("\nBGS is converged (strong coupling)")
        self.MPIPrint(" ")
        self.MPIPrint("*************************")
        self.MPIPrint("*  End FSI computation  *")
        self.MPIPrint("*************************")
        self.MPIPrint(" ")

    def MapModes(self, FSI_config, FluidSolver, SolidSolver):
        """
        Runs nothing, just extract the structural modes mapped on the fluid mesh
        """

        if self.have_MPI:
            myid = self.comm.Get_rank()
            numberPart = self.comm.Get_size()
        else:
            myid = 0
            numberPart = 1

        nodeNormals = {}
        for iVertex in range(self.nLocalFluidInterfaceNodes):
            if self.nDim == 2:
                nx, ny = FluidSolver.GetMarkerVertexNormals(
                    self.fluidInterfaceIdentifier, iVertex, False
                )
                nz = 0
            else:
                nx, ny, nz = FluidSolver.GetMarkerVertexNormals(
                    self.fluidInterfaceIdentifier, iVertex, False
                )
            GlobalIndex = FluidSolver.GetNodeGlobalIndex(
                FluidSolver.GetMarkerNode(self.fluidInterfaceIdentifier, iVertex)
            )
            nodeNormals[GlobalIndex] = [nx, ny, nz]

        nodeNormals = self.comm.gather(nodeNormals, root=self.rootProcess)
        if myid == self.rootProcess:
            normalsToPrint = {}
            for iDictionary in range(numberPart):
                for key, value in nodeNormals[iDictionary].items():
                    normalsToPrint[key] = value
            normalsToPrint = dict(sorted(normalsToPrint.items()))
            with open("Normals.csv", "w") as f:
                writer = csv.writer(f)
                for key, value in normalsToPrint.items():
                    writer.writerow([key, value])

        SurfaceFileName = FluidSolver.GetSurfaceFileName()

        self.MPIPrint("\n********************************")
        self.MPIPrint("* Begin mapping the modes *")
        self.MPIPrint("********************************\n")
        self.MPIPrint("\n")

        if (
            myid == self.rootProcess
        ):  # The root process contains the solid solver for sure
            modesNumber = np.array(int(SolidSolver.getNumberOfModes()))
        else:
            modesNumber = np.empty(1, dtype=np.int)

        self.comm.Bcast(modesNumber, root=self.rootProcess)

        for mode in range(np.asscalar(modesNumber)):
            self.MPIPrint("Setting mode {} active".format(mode))
            if myid in self.solidSolverProcessors:
                SolidSolver.activateMode(mode)
            self.MPIBarrier()
            self.getSolidInterfaceDisplacement(SolidSolver)
            self.interpolateSolidPositionOnFluidMesh(FSI_config)
            self.setFluidInterfaceVarCoord(FluidSolver)

            self.MPIPrint("\nPerforming mesh deformation...\n")
            FluidSolver.DynamicMeshUpdate(0)
            FluidSolver.Output(0)
            self.MPIBarrier()

            if myid == self.rootProcess:
                AllFiles = os.listdir()
                for FileNumber, FileName in enumerate(AllFiles):
                    if SurfaceFileName in FileName:
                        file = FileName.split(".")[0]
                        extension = FileName.split(".")[1]
                        os.rename(
                            file + "." + extension, "Mode{}.".format(mode) + extension
                        )

        self.MPIPrint("\n*************************")
        self.MPIPrint("*  Mapping completed  *")
        self.MPIPrint("*************************\n")
