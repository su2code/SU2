#!/usr/bin/env python
# -*-coding:utf-8 -* 

# \file FSIInterface.py
#  \brief FSI interface class that handles fluid/solid solvers synchronisation and communication.
#  \author D. THOMAS, University of Liege, Belgium. Department of Mechanical and Aerospace Engineering.
#  \version BETA

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
from mpi4py import MPI
import numpy as np
from math import *

# ----------------------------------------------------------------------
#  FSI Interface Class
# ----------------------------------------------------------------------

class Interface:
    """ 
    FSI interface class that handles fluid/solid solvers synchronisation and communication
    """
   
    def __init__(self, FSI_config, FluidSolver, SolidSolver):
	""" 
	Class constructor. Define some variables and do some screen outputs.
	"""
	self.rootProcess = 0			#the root process is imposed
        self.comm = MPI.COMM_WORLD		#MPI communicator
        self.fluidInterfaceIdentifier = None	#object that can identify the f/s interface within the fluid solver
        self.solidInterfaceIdentifier = None	#object that can identify the f/s interface within the solid solver
	self.nLocalFluidInterfaceNodes = 0	#number of interface fluid nodes for the current partition
	self.nLocalSolidInterfaceNodes = 0      #number of interface solid nodes for the current partition
	self.nFluidInterfaceNodes = 0		#number of interface fluid nodes (sum over all the partitions)
	self.nSolidInterfaceNodes = 0		#number of interface solid nodes
        self.localFluidInterface = {}		#fluid interface for the current partition, localFluidInterface[GlobalIndex] = (posx, posy, posz, Fx, Fy, Fz)
	self.globalFluidInterface= {}		#fluid interface after MPI reconstruction on the root process
	self.globalSolidInterface= {}		#solid interface on the root process (assume the solid solver is not parallelized), globalSolidInterface[GlobalIndex] = (posx, posy, posz, Fx, Fy, Fz)
	self.globalSolidInterfaceResidual = {}
	self.globalSolidInterfaceResidualnM1 = {}
	self.aitkenParam = float()
	self.FSIIter = 0
	self.mappingInterface = {}		#interface mapping/interpolation "matrix"
	self.Center = [0,0,0]			#in case of rigid body structural model, the center of rotation has to be tracked

        myid = self.comm.Get_rank()
        if myid == self.rootProcess:
	    print('Fluid solver : SU2_CFD')
	    print('Solid solver : {}'.format(FSI_config['CSD_SOLVER']))

	    if FSI_config['UNSTEADY_SIMULATION'] == 'YES':
		print('Unsteady coupled simulation with physical time step : {} s'.format(FSI_config['UNST_TIMESTEP']))
	    else:
		print('Steady coupled simulation')

	    if FSI_config['MATCHING_MESH'] == 'YES':
		print('Matching fluid-solid interface')
	    else:
		print('Non matching fluid-solid interface')

	    print('Solid predictor : {}'.format(FSI_config['DISP_PRED']))

	    print('Maximum number of FSI iterations : {}'.format(FSI_config['NB_FSI_ITER']))

	    print('FSI tolerance : {}'.format(FSI_config['FSI_TOLERANCE']))

	    if FSI_config['AITKEN_RELAX'] == 'STATIC':
		print('Static Aitken under-relaxation with constant parameter {}'.format(FSI_config['AITKEN_PARAM']))
	    elif FSI_config['AITKEN_RELAX'] == 'DYNAMIC':
		print('Dynamic Aitken under-relaxation with initial parameter {}'.format(FSI_config['AITKEN_PARAM']))
	    else:
		print('No Aitken under-relaxation')

            print('FSI interface is set')

    def connect(self, FluidSolver, SolidSolver):
	"""
	Connection between solvers. 
	Creates the communication support between the two solvers.
	Gets information about f/s interfaces from the two solvers.
	"""
        myid = self.comm.Get_rank()
	MPIsize = self.comm.Get_size()
	
	# identifies the fluid and solid interfaces and store the number of nodes on both sides (and for each partition)
        if FluidSolver != None:
	    print('Fluid solver is initialized on process {}'.format(myid))
	    self.fluidInterfaceIdentifier = FluidSolver.GetMovingMarker()
	    self.nLocalFluidInterfaceNodes = FluidSolver.GetNumberVertices(self.fluidInterfaceIdentifier)
	    if self.nLocalFluidInterfaceNodes != 0:
	      print('Number of interface fluid nodes (halo nodes included) on proccess {} : {}'.format(myid,self.nLocalFluidInterfaceNodes))
	else:
	    pass

	if SolidSolver != None:
	    print('Solid solver is initialized on process {}'.format(myid))
	    self.solidInterfaceIdentifier = SolidSolver.getFSIMarkerID()
	    self.nLocalSolidInterfaceNodes = SolidSolver.getNumberOfSolidInterfaceNodes(self.solidInterfaceIdentifier)
	else:
	    pass

	self.comm.barrier()
	
	# calculates the total number of nodes at the interface (sum over all the partitions)
	# fluid side
        sendBuff = np.array(int(self.nLocalFluidInterfaceNodes))
	rcvbuff = np.zeros(1, dtype=int)			
        self.comm.barrier()
	self.comm.Reduce(sendBuff,rcvbuff,op=MPI.SUM, root=self.rootProcess)
	self.nFluidInterfaceNodes = rcvbuff[0]
	# solid side
        sendBuff = np.array(int(self.nLocalSolidInterfaceNodes))
	rcvbuff = np.zeros(1, dtype=int)
	self.comm.barrier()
	self.comm.Reduce(sendBuff,rcvbuff,op=MPI.SUM, root=self.rootProcess)
	self.nSolidInterfaceNodes = rcvbuff[0]
        del sendBuff
        del rcvbuff
	if myid==self.rootProcess:
	    print('Total number of fluid nodes (halo nodes included) : {}'.format(self.nFluidInterfaceNodes))
	    print('Total number of solid nodes (halo nodes included) : {}'.format(self.nSolidInterfaceNodes))

	self.comm.barrier()

	# gets the fluid interface from fluid solver (!on each process!)
        self.getFluidInterfacePosition(FluidSolver)

	# gets the solid interface from solid solver
	self.getSolidInterfacePosition(SolidSolver)

	# gathers the fluid interface on the root process
        buff = self.localFluidInterface
	self.__ReconstructFluidInterface(buff)
          
	if myid==self.rootProcess:
	  print('Total number of physical fluid nodes : {}'.format(len(self.globalFluidInterface)))
	  # we assume the solid is only running on one process (we skip the gather)
	  print('Total number of physical solid nodes : {}'.format(self.nSolidInterfaceNodes))

	self.comm.barrier()

    def interfaceMapping(self,FluidSolver, SolidSolver, FSI_config):
	""" 
	Creates the one-to-one mapping between interfaces in case of matching meshes.
	Creates the interpolation rules between interfaces in case of non-matching meshes.
	"""
	myid = self.comm.Get_rank()
	MPIsize = self.comm.Get_size()

	# the f/s interfaces mapping/interpolation is performed on the root process only
	if myid==self.rootProcess:
	    if FSI_config['MATCHING_MESH'] == 'YES':			# in case of matching meshes, the mapping is one-to-one relation between node identifiers
	        print('Mapping for matching meshes')
		if len(self.globalFluidInterface) != len(self.globalSolidInterface):
		    raise Exception("Fluid and solid interface must have the same number of nodes for matching meshes ! ")
		for fluidNode in self.globalFluidInterface.iterkeys():
		    minDist = 1e10
		    self.mappingInterface[fluidNode] = None
		    for solidNode in self.globalSolidInterface.iterkeys():
  		        tempDist = self.__computeDistance(self.globalFluidInterface[fluidNode], self.globalSolidInterface[solidNode])	# mapping based on the distance between fluid and solid nodes
			if tempDist < minDist:
			    minDist = tempDist
			    self.mappingInterface[fluidNode] = solidNode
		    if minDist > 1e-10:
			print('Tolerance for matching meshes is not matched !')				# checking that the distance is small enough to consider the two nodes as matching
		print('FSI interface is mapped for {} matching nodes.'.format(len(self.globalFluidInterface)))
	    else:
		print('Interpolating non matching meshes')

    def interpolateSolidPositionOnFluidMesh(self, FSI_config):
	"""
	Applies the one-to-one mapping or the interpolaiton rules from solid to fluid mesh.
	"""
	myid = self.comm.Get_rank()

	# the f/s interfaces interpolation is performed on the root process only
	if myid==self.rootProcess:	
	    if FSI_config['MATCHING_MESH'] == 'YES':			# in case of matching meshes, simply apply the one-to-one mapping
	        for fluidNode in self.globalFluidInterface.iterkeys():
	            self.globalFluidInterface[fluidNode] = self.globalSolidInterface[self.mappingInterface[fluidNode]]

    def interpolateFluidLoadsOnSolidMesh(self, FSI_config):
	"""
	Applies the one-to-one mapping or the interpolaiton rules from fluid to solid mesh.
	"""
	myid = self.comm.Get_rank()

	#The f/s interfaces interpolation is performed on the root process only
	if myid==self.rootProcess:	
	    if FSI_config['MATCHING_MESH'] == 'YES':			#in case of matching meshes, simply apply the one-to-one mapping
	        for fluidNode in self.globalFluidInterface.iterkeys():
		     self.globalSolidInterface[self.mappingInterface[fluidNode]] = self.globalFluidInterface[fluidNode]

    def __computeDistance(self, nodeA, nodeB):
	"""
	Computes the distance between two nodes (x,y,z)
	"""
	myid = self.comm.Get_rank()
	MPIsize = self.comm.Get_size()

	distance = sqrt((nodeB[0]-nodeA[0])**2 + (nodeB[1]-nodeA[1])**2 + (nodeB[2]-nodeA[2])**2)
	
	return distance

    def getSolidInterfacePosition(self, SolidSolver):
	"""
	Gets the current solid interface position from the solid solver.
	"""
	myid = self.comm.Get_rank()
	
	GlobalIndex = int()
	normVarPosSquare = 0.0
	if myid==self.rootProcess:
		for iVertex in range(self.nSolidInterfaceNodes):
		    GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
		    newPosx = SolidSolver.getInterfaceNodePosX(self.solidInterfaceIdentifier, iVertex)
		    newPosy = SolidSolver.getInterfaceNodePosY(self.solidInterfaceIdentifier, iVertex)
		    newPosz = SolidSolver.getInterfaceNodePosZ(self.solidInterfaceIdentifier, iVertex)
		    self.globalSolidInterface[GlobalIndex] = (newPosx,newPosy,newPosz, 0.0, 0.0, 0.0)
 	    	self.Center[0] = SolidSolver.getRotationCenterPosX()
            	self.Center[1] = SolidSolver.getRotationCenterPosY()
            	self.Center[2] = SolidSolver.getRotationCenterPosZ()
	    	print("Center of rotation : {};{};{}".format(self.Center[0],self.Center[1],self.Center[2]))
	return sqrt(normVarPosSquare)

    def computeSolidInterfaceResidual(self, SolidSolver):
	"""
	Computes the solid interface displacement residual.
	"""

	myid = self.comm.Get_rank()

	normInterfaceResidualSquare = 0.0
	
	if myid == self.rootProcess:
	    for iVertex in range(self.nSolidInterfaceNodes):
		GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
		predPosx = SolidSolver.getInterfaceNodePosX(self.solidInterfaceIdentifier, iVertex)
		predPosy = SolidSolver.getInterfaceNodePosY(self.solidInterfaceIdentifier, iVertex)
		predPosz = SolidSolver.getInterfaceNodePosZ(self.solidInterfaceIdentifier, iVertex)
		varPosx = predPosx - self.globalSolidInterface[GlobalIndex][0]
		varPosy = predPosy - self.globalSolidInterface[GlobalIndex][1]
		varPosz = predPosz - self.globalSolidInterface[GlobalIndex][2]
		normInterfaceResidualSquare += (varPosx**2+varPosy**2+varPosz**2)
		self.globalSolidInterfaceResidual[GlobalIndex] = (varPosx,varPosy,varPosz,0.0,0.0,0.0)

	return sqrt(normInterfaceResidualSquare)

    def relaxSolidPosition(self,FSI_config, SolidSolver):
	"""
	Apply solid displacement under-relaxation.
	"""
	myid = self.comm.Get_rank()

	if FSI_config['AITKEN_RELAX'] == 'STATIC':
	    self.aitkenParam = FSI_config['AITKEN_PARAM']
	elif FSI_config['AITKEN_RELAX'] == 'DYNAMIC':	
	    self.setAitkenCoefficient(FSI_config, SolidSolver)
	else:
	    self.aitkenParam = 1.0

	if myid == self.rootProcess:
	    print('Aitken under-relaxation step with parameter {}'.format(self.aitkenParam))
	    for iVertex in range(self.nSolidInterfaceNodes):
	        GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
	        resx, resy, resz, fx, fy, fz = self.globalSolidInterfaceResidual[GlobalIndex]
	        newPosx = self.globalSolidInterface[GlobalIndex][0] + self.aitkenParam*resx
	        newPosy = self.globalSolidInterface[GlobalIndex][1] + self.aitkenParam*resy
	        newPosz = self.globalSolidInterface[GlobalIndex][2] + self.aitkenParam*resz
	        self.globalSolidInterface[GlobalIndex] = (newPosx,newPosy,newPosz, 0.0, 0.0, 0.0)
	    if FSI_config['CSD_SOLVER'] == 'GETDP':
	        print('Computed position of the node 7 : ({} , {} , {})'.format(self.globalSolidInterface[7][0],self.globalSolidInterface[7][1],self.globalSolidInterface[7][2]))
	    if FSI_config['CSD_SOLVER'] == 'METAFOR':
	        print('Computed position of the node 7 : ({} , {} , {})'.format(self.globalSolidInterface[7][0],self.globalSolidInterface[7][1],self.globalSolidInterface[7][2]))
	

    def setAitkenCoefficient(self, FSI_config, SolidSolver):
	"""
	Computes the Aitken coefficients for solid displacement under-relaxation.
	"""

	deltaResNormSquare = 0.0
	prodScalRes = 0.0
	
	if self.FSIIter == 0:
	    self.aitkenParam = FSI_config['AITKEN_PARAM']
	else:
	    for iVertex in range(self.nSolidInterfaceNodes):
		GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
	        resx, resy, resz, fx, fy, fz = self.globalSolidInterfaceResidual[GlobalIndex]
		resxnM1, resynM1, resznM1, fx, fy, fz = self.globalSolidInterfaceResidualnM1[GlobalIndex]
		deltaResx = resx - resxnM1
		deltaResy = resy - resynM1
		deltaResz = resz - resznM1
		#prodScalRes += (deltaResx*resx+deltaResy*resy+deltaResz*resz)
		prodScalRes += (deltaResx*resxnM1+deltaResy*resynM1+deltaResz*resznM1)
		deltaResNormSquare += deltaResx**2 + deltaResy**2 + deltaResz**2
	    self.aitkenParam *= -prodScalRes/deltaResNormSquare
	self.globalSolidInterfaceResidualnM1 = self.globalSolidInterfaceResidual.copy()

    def getFluidInterfacePosition(self, FluidSolver):
	"""
	Gets the current fluid interface position from the fluid solver.
	"""

	GlobalIndex = int()
        for iVertex in range(self.nLocalFluidInterfaceNodes):
	    GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
	    posx = FluidSolver.GetVertexCoordX(self.fluidInterfaceIdentifier, iVertex)
	    posy = FluidSolver.GetVertexCoordY(self.fluidInterfaceIdentifier, iVertex)
	    posz = FluidSolver.GetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex)
            self.localFluidInterface[GlobalIndex] = (posx,posy,posz, 0.0 , 0.0, 0.0)

    def getFluidInterfaceNodalForce(self, FSI_config, FluidSolver):
	"""
	Gets the fluid interface loads from the fluid solver.
	"""
	buff = {}	# this will be gathered

	for iVertex in range(self.nLocalFluidInterfaceNodes):
	    GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
	    halo = FluidSolver.ComputeVertexForces(self.fluidInterfaceIdentifier, iVertex)		# !!we have to ignore halo node coming from mesh partitioning because they introduice non-physical forces
	    if halo==False:
		if FSI_config['CSD_SOLVER'] == 'GETDP':
		    newFx = FluidSolver.GetVertexForceDensityX(self.fluidInterfaceIdentifier, iVertex)
	            newFy = FluidSolver.GetVertexForceDensityY(self.fluidInterfaceIdentifier, iVertex)
	            newFz = FluidSolver.GetVertexForceDensityZ(self.fluidInterfaceIdentifier, iVertex)
		else:
	            newFx = FluidSolver.GetVertexForceX(self.fluidInterfaceIdentifier, iVertex)
	            newFy = FluidSolver.GetVertexForceY(self.fluidInterfaceIdentifier, iVertex)
	            newFz = FluidSolver.GetVertexForceZ(self.fluidInterfaceIdentifier, iVertex)
		posx,posy,posz,fx,fy,fz = self.localFluidInterface[GlobalIndex]
		self.localFluidInterface[GlobalIndex] = (posx,posy,posz,newFx,newFy,newFz)
		buff[GlobalIndex] = self.localFluidInterface[GlobalIndex]
	return buff

    def setFluidInterfaceVarCoord(self, FluidSolver):
	"""
	Communicates the change of coordinates of the fluid interface to the fluid solver.
	Prepares the fluid solver for mesh deformation.
	"""
	myid = self.comm.Get_rank()
	
	varCoordNormSquare = 0.0	# this can be used for FSI convergence criterion
	for iVertex in range(self.nLocalFluidInterfaceNodes):
	    GlobalIndex = FluidSolver.GetVertexGlobalIndex(self.fluidInterfaceIdentifier, iVertex)
	    FluidSolver.SetVertexCoordX(self.fluidInterfaceIdentifier, iVertex, self.globalFluidInterface[GlobalIndex][0])
	    FluidSolver.SetVertexCoordY(self.fluidInterfaceIdentifier, iVertex, self.globalFluidInterface[GlobalIndex][1])
	    FluidSolver.SetVertexCoordZ(self.fluidInterfaceIdentifier, iVertex, self.globalFluidInterface[GlobalIndex][2])
	    nodalVarCoordNorm = FluidSolver.SetVertexVarCoord(self.fluidInterfaceIdentifier, iVertex)		# prepares the mesh deformation in the fluid solver
	    varCoordNormSquare += nodalVarCoordNorm**2
	    posx,posy,posz,fx,fy,fz = self.localFluidInterface[GlobalIndex]
	    newposx,newposy,newposz,FX,FY,FZ = self.globalFluidInterface[GlobalIndex]
	    self.localFluidInterface[GlobalIndex] = (newposx,newposy,newposz,fx,fy,fz)
	return sqrt(varCoordNormSquare)

    def setSolidInterfaceLoads(self, SolidSolver, FSI_config):
	"""
	Communicates the new solid interface loads to the solid solver.
	In case of rigid body motion, calculates the new resultant forces (lift, drag, ...).
	"""
	myid = self.comm.Get_rank()

	FY = 0.0 # solid-side resultant forces
        FX = 0.0
        FZ = 0.0
        Mz = 0.0
	FFX = 0.0 # fluid-side resultant forces
	FFY = 0.0
	FFZ = 0.0
 	# checking force conservation after interpolation    
        for node in self.globalSolidInterface.iterkeys():
  	    FX += self.globalSolidInterface[node][3]
	    FY += self.globalSolidInterface[node][4]
	    FZ += self.globalSolidInterface[node][5]				
	    Mz += (self.globalSolidInterface[node][4]*(self.globalSolidInterface[node][0]-self.Center[0])-self.globalSolidInterface[node][3]*(self.globalSolidInterface[node][1]-self.Center[1]));
	for node in self.globalFluidInterface.iterkeys():
	    FFX += self.globalFluidInterface[node][3]
	    FFY += self.globalFluidInterface[node][4]
	    FFZ += self.globalFluidInterface[node][5]
	print("Checking f/s interface conservation...")
	print('Solid-side FX = {}'.format(FX))        
	print('Solid-side FY = {}'.format(FY))        
        print('Solid-side FZ = {}'.format(FZ))
	print('Fluid-side FX = {}'.format(FFX))        
	print('Fluid-side FY = {}'.format(FFY))        
        print('Fluid-side FZ = {}'.format(FFZ))
	# in case of NativeSolid, we communicate only global fluid loads
	if FSI_config['CSD_SOLVER'] == 'NATIVE':
	    SolidSolver.setGeneralisedForce(FX, FY)
	    SolidSolver.setGeneralisedMoment(Mz)
	# in other cases, we communicate nodal fluid loads
	if FSI_config['CSD_SOLVER'] == 'METAFOR' or FSI_config['CSD_SOLVER'] == 'GETDP':
	    for solidNode in self.globalSolidInterface.iterkeys():
	        Fx = self.globalSolidInterface[solidNode][3]
	        Fy = self.globalSolidInterface[solidNode][4]
	        SolidSolver.applyload(solidNode, Fx, Fy)

    def displacementPredictor(self, FSI_config , SolidSolver, deltaT):
	"""
	Calculates a prediciton for the solid interface deformation for the next time step.
	"""

	if FSI_config['DISP_PRED'] == 'FIRST_ORDER':
	    print("First order predictor")	
	    alpha_0 = 1.0
	    alpha_1 = 0.0
	elif FSI_config['DISP_PRED'] == 'SECOND_ORDER':
	    print("Second order predictor")
	    alpha_0 = 1.0
	    alpha_1 = 0.5
	else:
	    print("No predictor")
	    alpha_0 = 0.0
	    alpha_1 = 0.0

	for iVertex in range(self.nSolidInterfaceNodes):
	    GlobalIndex = SolidSolver.getInterfaceNodeGlobalIndex(self.solidInterfaceIdentifier, iVertex)
	    posx,posy,posz,fx,fy,fz = self.globalSolidInterface[GlobalIndex]
	    velx = SolidSolver.getInterfaceNodeVelX(self.solidInterfaceIdentifier, iVertex)
	    vely = SolidSolver.getInterfaceNodeVelY(self.solidInterfaceIdentifier, iVertex)
	    velz = SolidSolver.getInterfaceNodeVelZ(self.solidInterfaceIdentifier, iVertex)
	    velxNm1 = SolidSolver.getInterfaceNodeVelXNm1(self.solidInterfaceIdentifier, iVertex)
	    velyNm1 = SolidSolver.getInterfaceNodeVelYNm1(self.solidInterfaceIdentifier, iVertex)
	    velzNm1 = SolidSolver.getInterfaceNodeVelZNm1(self.solidInterfaceIdentifier, iVertex)
	    posxPred = posx + alpha_0*deltaT*velx + alpha_1*deltaT*(velx-velxNm1)
	    posyPred = posy + alpha_0*deltaT*vely + alpha_1*deltaT*(vely-velyNm1)
	    poszPred = posz + alpha_0*deltaT*velz + alpha_1*deltaT*(velz-velzNm1)
            self.globalSolidInterface[GlobalIndex] = (posxPred,posyPred,poszPred,fx,fy,fz)
	
   	if FSI_config['CSD_SOLVER'] == 'GETDP':
	    print('Predicted position of the node 7 : ({} , {} , {})'.format(self.globalSolidInterface[7][0],self.globalSolidInterface[7][1],self.globalSolidInterface[7][2]))
	if FSI_config['CSD_SOLVER'] == 'METAFOR':
	    print('Predicted position of the node 7 : ({} , {} , {})'.format(self.globalSolidInterface[7][0],self.globalSolidInterface[7][1],self.globalSolidInterface[7][2]))
	    

    def __BroadcastFluidInterface(self):
	"""
	Distributes the fluid interface over all the partitions.
	"""
	self.globalFluidInterface = self.comm.bcast(self.globalFluidInterface, root=self.rootProcess)

    def __ReconstructFluidInterface(self, buff):
	"""
	Reconstructs the fluid interface on the root process from all the partitions.
	"""
	
	myid = self.comm.Get_rank()	
	MPIsize = self.comm.Get_size()

	buff  = self.comm.gather(buff, root=self.rootProcess)	# use of a buffer because it is cleared after the gather
	if myid ==self.rootProcess:
          for i in range(MPIsize):
            for key, value in buff[i].items():
              self.globalFluidInterface[key] = value
	del buff

    def writeFSIHistory(self, TimeIter, time, varCoordNorm, FSIConv):
	"""
	Write the fsi history file of the computaion.
	"""
	
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
	    

    def UnsteadyFSI(self,FSI_config, FluidSolver, SolidSolver):
	  """ 
	  Run the unsteady FSI computation by synchronizing the fluid and solid solvers.
	  F/s interface data are exchanged through interface mapping and interpolation (if non mathcing meshes).
	  """

	  myid = self.comm.Get_rank()
	  numberPart = self.comm.Get_size()

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

	  if myid == self.rootProcess:
	    print('\n**********************************')
	    print('* Begin unsteady FSI computation *')
	    print('**********************************\n')
	  
	  # --- Initialize the coupled solution --- #
	  #If restart (DOES NOT WORK YET)
	  if FSI_config['RESTART_SOL'] == 'YES':
	    TimeIterTreshold = -1
	    FluidSolver.setTemporalIteration(TimeIter)
	    if myid == self.rootProcess:
	      SolidSolver.outputDisplacements(FluidSolver.getInterRigidDispArray(), True)
	    self.comm.barrier()
	    FluidSolver.setInitialMesh(True)
	    if myid == self.rootProcess:
	      SolidSolver.displacementPredictor(FluidSolver.getInterRigidDispArray())
	    self.comm.barrier()								
	    if myid == self.rootProcess:
	      SolidSolver.updateSolution()
	  #If no restart
	  else:
	    if myid ==self.rootProcess:
	      print('Setting FSI initial conditions')
	      SolidSolver.setInitialDisplacements()
	      self.getSolidInterfacePosition(SolidSolver)
	      self.interpolateSolidPositionOnFluidMesh(FSI_config)
	    self.comm.barrier()
	    self.__BroadcastFluidInterface()
	    self.setFluidInterfaceVarCoord(FluidSolver)
	    FluidSolver.SetInitialMesh()	# if there is an initial deformation in the solid, it has to be communicated to the fluid solver
	    if myid == self.rootProcess:
	      print('\nFSI initial conditions are set')
	      print('Beginning time integration\n')

	  # --- External temporal loop --- #
	  while TimeIter <= NbTimeIter:

		if TimeIter > TimeIterTreshold:
		  NbFSIIter = NbFSIIterMax
		  if myid == self.rootProcess:
		    print('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI on time iteration {} ***************'.format(TimeIter))
		  
		else:
		  NbFSIIter = 1

		self.FSIIter = 0
		FSIConv = False
		FluidSolver.PreprocessExtIter(TimeIter)	# set some parameters before temporal fluid iteration
	
		# --- Internal FSI loop --- #
		while self.FSIIter <= (NbFSIIter-1):

			if myid == self.rootProcess:
			  print("\n>>>> Time iteration {} / FSI iteration {} <<<<".format(TimeIter,self.FSIIter))
			  # --- Mesh morphing step (displacements interpolation, displacements communication, and mesh morpher call) --- #
			  self.interpolateSolidPositionOnFluidMesh(FSI_config)
                          print('\nPerforming dynamic mesh deformation (ALE)...\n')
                        self.comm.barrier()
                        self.__BroadcastFluidInterface()
                        self.setFluidInterfaceVarCoord(FluidSolver)
                        FluidSolver.DynamicMeshUpdate(TimeIter)
			
			# --- Fluid solver call for FSI subiteration --- #
			if myid == self.rootProcess:
		          print('\nLaunching fluid solver for one single dual-time iteration...')
			self.comm.barrier()
			FluidSolver.ResetConvergence()
			FluidSolver.Run()

			# --- Surface fluid loads interpolation and communication --- #
		        if myid == self.rootProcess:
		          print('\nProcessing surface fluid loads...')
		 	self.comm.barrier()
		        buff = self.getFluidInterfaceNodalForce(FSI_config, FluidSolver)
		        self.__ReconstructFluidInterface(buff)
		   	self.comm.barrier()
			if TimeIter > TimeIterTreshold:
			  if myid == self.rootProcess:
			    self.interpolateFluidLoadsOnSolidMesh(FSI_config)
			    self.setSolidInterfaceLoads(SolidSolver, FSI_config)

			    # --- Solid solver call for FSI subiteration --- #
			    print('\nLaunching solid solver for a single time iteration...\n')
			    if FSI_config['CSD_SOLVER'] == 'NATIVE':
			        SolidSolver.timeIteration(time)
			    elif FSI_config['CSD_SOLVER'] == 'METAFOR' or FSI_config['CSD_SOLVER'] == 'GETDP':
				SolidSolver.run(time-deltaT, time)

			    # --- Compute and monitor the FSI residual --- #
			    varCoordNorm = self.computeSolidInterfaceResidual(SolidSolver)
			    print('\nFSI displacement norm : {}\n'.format(varCoordNorm))
			  self.comm.barrier()
			  varCoordNorm = self.comm.bcast(varCoordNorm, root=0)
			  if varCoordNorm < FSITolerance:		
			    FSIConv = True
			    break

			  # --- Relaxe the solid position --- #
			  if myid == self.rootProcess:
			    self.relaxSolidPosition(FSI_config, SolidSolver)			
			self.FSIIter += 1
		# --- End OF FSI loop --- #

		# --- Update the FSI history file --- # 
		if myid == self.rootProcess:
		  if TimeIter > TimeIterTreshold:
		    print('\nBGS is converged (strong coupling)')
		  self.writeFSIHistory(TimeIter, time, varCoordNorm, FSIConv)
		
		# --- Update, monitor and output the fluid solution before the next time step  ---#
		FluidSolver.Update()
		FluidSolver.Monitor()
		FluidSolver.Output()

	   	if TimeIter >= TimeIterTreshold:
		  if myid == self.rootProcess:
		    # --- Output the solid solution before thr next time step --- #
		    SolidSolver.writeSolution(time, self.FSIIter, TimeIter, NbTimeIter)
		
		    # --- Displacement predictor for the next time step and update of the solid solution --- #
		    print('\nSolid displacement prediction for next time step')
		    if FSI_config['CSD_SOLVER'] == 'NATIVE':
			SolidSolver.updateGeometry()
		        SolidSolver.displacementPredictor()
			self.getSolidInterfacePosition(SolidSolver)
		    elif FSI_config['CSD_SOLVER'] == 'METAFOR' or FSI_config['CSD_SOLVER'] == 'GETDP':
			self.displacementPredictor(FSI_config, SolidSolver, deltaT)
		    SolidSolver.updateSolution()
		
		TimeIter += 1
		time += deltaT
	  #--- End of the temporal loop --- #

	  self.comm.barrier()

	  if myid == self.rootProcess:
	    print('\n*************************')
	    print('*  End FSI computation  *')
	    print('*************************\n')

    def SteadyFSI(self, FSI_config,FluidSolver, SolidSolver):
	  """
	  Runs the steady FSI computation by synchronizing the fluid and solid solver with data exchange at the f/s interface.
	  """

	  myid = self.comm.Get_rank()
	  numberPart = self.comm.Get_size()

	  # --- Set some general variables for the steady computation --- #
	  NbIter = FSI_config['NB_EXT_ITER']		# number of fluid iteration at each FSI step
	  NbFSIIterMax = FSI_config['NB_FSI_ITER']	# maximum number of FSI iteration (for each time step)
	  FSITolerance = FSI_config['FSI_TOLERANCE']	# f/s interface tolerance
	  varCoordNorm = 0.0

	  if myid == self.rootProcess:
	    print('\n********************************')
	    print('* Begin steady FSI computation *')
	    print('********************************\n')
	    print('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI ***************')

	  # --- External FSI loop --- #
	  self.FSIIter = 0
	  while self.FSIIter < NbFSIIterMax:
	    if myid == self.rootProcess:
	      print("\n>>>> FSI iteration {} <<<<".format(self.FSIIter))
	      print('\nLaunching fluid solver for a steady computation...')
	    # --- Fluid solver call for FSI subiteration ---#
	    Iter = 0
	    FluidSolver.ResetConvergence()          
	    while Iter < NbIter:
	      FluidSolver.PreprocessExtIter(Iter)
	      FluidSolver.Run()				
	      StopIntegration = FluidSolver.Monitor()
	      FluidSolver.Output()
	      if StopIntegration:
		break;
	      Iter += 1
	    
	    # --- Surface fluid loads interpolation and communication ---#
            if myid == self.rootProcess:
		print('\nProcessing surface fluid loads...')
	    self.comm.barrier()
	    buff = self.getFluidInterfaceNodalForce(FSI_config, FluidSolver)
	    self.__ReconstructFluidInterface(buff)
	    self.comm.barrier()
	    if myid == self.rootProcess:
	      self.interpolateFluidLoadsOnSolidMesh(FSI_config)
	      self.setSolidInterfaceLoads(SolidSolver, FSI_config)
	     
	      # --- Solid solver call for FSO subiteration --- #
	      print('\nLaunching solid solver for a static computation...\n')
	      if FSI_config['CSD_SOLVER'] == 'NATIVE':
	          SolidSolver.staticComputation()	
	      SolidSolver.writeSolution(0.0, self.FSIIter, Iter, NbIter)		

	      # --- Compute and monitor the FSI residual --- #
	      varCoordNorm = self.computeSolidInterfaceResidual(SolidSolver)
	      print('\n>>>> FSI Residual : {} <<<<'.format(varCoordNorm))
	    self.comm.barrier()
	    varCoordNorm = self.comm.bcast(varCoordNorm, root=0)
	    if varCoordNorm < FSITolerance:			
	      break

            # --- Relaxe the solid displacement and update the solid solution --- #
	    if myid == self.rootProcess:
	      self.relaxSolidPosition(FSI_config, SolidSolver)
              SolidSolver.updateSolution()
	
	      #--- Mesh morphing step (displacement interpolation, displacements communication, and mesh morpher call) --- #
	      self.interpolateSolidPositionOnFluidMesh(FSI_config)						
	      print('\nPerforming static mesh deformation...\n')
	    self.comm.barrier()
	    self.__BroadcastFluidInterface()
	    self.setFluidInterfaceVarCoord(FluidSolver)									
	    FluidSolver.StaticMeshUpdate()
	    self.FSIIter += 1

	  self.comm.barrier()

	  if myid == self.rootProcess:
	    print('\nBGS is converged (strong coupling)')
	    print(' ')
	    print('*************************')
	    print('*  End FSI computation  *')
	    print('*************************')
	    print(' ')
