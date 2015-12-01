# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pypar
from ..io import MPIComm
from math import *

# -------------------------------------------------------------------
#  Unsteady FSI computation 
# -------------------------------------------------------------------

def UnsteadyFSI(FSI_config, FluidSolver, SolidSolver, time, TimeIterTreshold, NbTimeIter, NbFSIIterMax, FSITolerance, TimeIter, deltaT):

  MASTER_NODE = 0
  myid = pypar.rank()
  DispNorm = 0.0

  if myid == MASTER_NODE:
    print(' ')
    print('*************************')
    print('* Begin FSI computation *')
    print('*************************')
    print(' ')
  
  #If restart
  if FSI_config['RESTART_SOL'] == 'YES':
    TimeIterTreshold = -1
    FluidSolver.setTemporalIteration(TimeIter)
    if myid == MASTER_NODE:
      SolidSolver.outputDisplacements(FluidSolver.getInterRigidDispArray(), True)
    pypar.barrier()
    FluidSolver.setInitialMesh(True)
    if myid == MASTER_NODE:
      SolidSolver.displacementPredictor(FluidSolver.getInterRigidDispArray())		#Compute the solid displacement predictor for the next time step
    pypar.barrier()								
    if myid == MASTER_NODE:
      SolidSolver.updateSolution()
  #If no restart
  else:
    if myid == MASTER_NODE:
      print('Setting FSI initial conditions')
      FluidSolver.getSolidDisplacement(SolidSolver.getSolidInterface(), SolidSolver.getCenterCoordinate())
    pypar.barrier()
    FluidSolver.setInitialMesh(False)
    if myid == MASTER_NODE:
      print('\nFSI initial conditions are set')
      print('Beginning time integration\n')
  
  #Initialise the FSIResidual object
  StopComp = FSIResidual()

  #TIME LOOP
  while TimeIter <= NbTimeIter:
        if TimeIter > TimeIterTreshold:
          NbFSIIter = NbFSIIterMax
          if myid == MASTER_NODE:
	    print('\n*************** Enter Block Gauss Seidel (BGS) method for strong coupling FSI on time iteration {} ***************'.format(TimeIter))
          
        else:
          NbFSIIter = 1
	FSIIter = 0
	FluidSolver.setTemporalIteration(TimeIter)	#Set some parameters before temporal fluid iteration (may disappear in the future)
	
	#FSI LOOP
	while FSIIter <= (NbFSIIter-1):
		if myid == MASTER_NODE:
		  print("\n>>>> Time iteration {} / FSI iteration {} <<<<".format(TimeIter,FSIIter))
                  print('\nLaunching fluid solver for one single dual-time iteration...')
		FluidSolver.dualTimeInnerLoop(TimeIter)		#Fluid sub-iteration (dual-time integration)
                Drag = FluidSolver.outputFluidLoads_Drag()	#Output global aero load : drag (in [N])
		Lift = FluidSolver.outputFluidLoads_Lift()	#Output global aero load : lift (in [N]) only to compute FSI residual (based on loads)
                if myid == MASTER_NODE:
                  print('Processing surface fluid loads...')
                FluidSolver.mergeSurfaceLoads()
                pypar.barrier()
		if TimeIter > TimeIterTreshold:
		  if myid == MASTER_NODE:
                    FluidSolver.outputSurfaceLoads(SolidSolver.getSolidSurfaceLoads())          #Transfer the global fluid loads (lift or drag) to the solid solver
		    SolidSolver.applyGlobalFluidLoads() 						#Apply the global fluid loads to the solid structure
		    #??.mapFluidLoadsOnSolidMesh()
		    print('\nLaunching solid solver for a single time iteration...\n')
		    SolidSolver.timeIteration(time)	#Solid sub-iteration, the current solid solution is automatically reset
		    DispNorm = SolidSolver.getVarCoordNorm()
                    FluidSolver.getSolidDisplacement(SolidSolver.getSolidInterface(), SolidSolver.getCenterCoordinate())	#Communicate the new position of the FSI interface (solid side)
                    print('\nPerforming dynamic mesh deformation (ALE)...\n')
                  pypar.barrier()
                  DispNorm = MPIComm.BroadcastOneDouble(DispNorm)
                  FluidSolver.mapRigidDisplacementOnFluidMesh()		#Interpolate the FSI interface displacement from solid to fluid mesh
		  FluidSolver.dynamicMeshUpdate(TimeIter)		#Deform the volume mesh, calculate grid velocity, update multi-grid structure
		  #print("Relaxation step")
		  #SolidSolver.setAitkenCoefficient(FSIIter)		#Implemented but not functional yet
                  #print(DispNorm)
  		  res = StopComp.setValue(FSIIter, DispNorm) 		#Calculate the residual based on fluid global loads (lift or drag)		
		  if myid == MASTER_NODE: 
		    print(' ')
		    #print('>>>> FSI Residual : {} <<<<'.format(res))
		    print(DispNorm)
		    print(' ')
 		  #if res < FSITolerance:					#Exit the FSI loop if the specified tolerance is reached
                  if DispNorm < FSITolerance:
		    break	
		FSIIter += 1
	#END OF FSI LOOP
        if myid == MASTER_NODE:
          if TimeIter > TimeIterTreshold:
            print('\nBGS is converged (strong coupling)')
        
	FluidSolver.updateDualTime()				#Update the dual-time solution of the fluid solver
	StopIntegration = FluidSolver.writeSolution(TimeIter)	#Write the solution in the history file once the FSI loop is converged

   	if TimeIter >= TimeIterTreshold:
	  if myid == MASTER_NODE:
            SolidSolver.writeSolution(time, FSIIter, TimeIter, NbTimeIter)				#Write the solid solution in the history file once the FSI loop is converged
            SolidSolver.updateGeometry()
	    print('\nHigh order solid displacement prediction for next time step')
            SolidSolver.displacementPredictor()
            FluidSolver.getSolidDisplacement(SolidSolver.getSolidInterface(), SolidSolver.getCenterCoordinate())
            SolidSolver.updateSolution()				#Update the solid solution for the next time-step
            print('\nPerforming dynamic mesh deformation (ALE)...\n')
          pypar.barrier()
          FluidSolver.mapRigidDisplacementOnFluidMesh()
	  FluidSolver.dynamicMeshUpdate(TimeIter)		#Deform the volume mesh, calculate grid velocity, update multi-grid structure

	TimeIter += 1
	time += deltaT
  #END OF TIME LOOP

  pypar.barrier()

  if myid == MASTER_NODE:
    print(' ')
    print('*************************')
    print('*  End FSI computation  *')
    print('*************************')
    print(' ')

#--------------------------------------------------------------

def UnsteadyFSI_Old(FSI_config, FluidSolver, SolidSolver, time, TimeIterTreshold, NbTimeIter, NbFSIIterMax, FSITolerance, TimeIter, deltaT):

  MASTER_NODE = 0
  myid = pypar.rank()

  if myid == MASTER_NODE:
    print(' ')
    print('*************************')
    print('* Begin FSI computation *')
    print('*************************')
    print(' ')

  #TIME LOOP
  #TimeIter = 0
  #time = 0.0
  

  #If restart
  if FSI_config['RESTART_SOL'] == 'YES':
    TimeIterTreshold = -1
    FluidSolver.setTemporalIteration(TimeIter)
    if myid == MASTER_NODE:
      SolidSolver.outputDisplacements(FluidSolver.getInterRigidDispArray(), True)
    pypar.barrier()
    FluidSolver.broadcastInterfDisp()
    pypar.barrier()
    FluidSolver.setInitialMesh(True)
    if myid == MASTER_NODE:
      SolidSolver.displacementPredictor_Old(FluidSolver.getInterRigidDispArray())		#Compute the solid displacement predictor for the next time step
    pypar.barrier()								
    FluidSolver.broadcastInterfDisp()			#Broadcast the predicted displacement from MASTER_NODE to other processes (use of a buffer)
    pypar.barrier()
    if myid == MASTER_NODE:
      SolidSolver.updateSolution()
  #If no restart
  else:
    if myid == MASTER_NODE:
      print('Setting FSI initial conditions')
      SolidSolver.outputDisplacements(FluidSolver.getInterRigidDispArray(), True) 
    pypar.barrier()							
    FluidSolver.broadcastInterfDisp()		#Broadcast the displacement from MASTER_NODE to other processes (use of a buffer)
    pypar.barrier()
    FluidSolver.setInitialMesh_Old(False)
    if myid == MASTER_NODE:
      print('\nFSI initial conditions are set')
      print('Beginning time integration\n')
      #SolidSolver.updateSolution()
  

  StopComp = FSIResidual()
  while TimeIter <= NbTimeIter:
        if TimeIter > TimeIterTreshold:
          NbFSIIter = NbFSIIterMax
        else:
          NbFSIIter = 1
	FSIIter = 0
	FluidSolver.setTemporalIteration(TimeIter)
	
	#FSI LOOP
	while FSIIter <= (NbFSIIter-1):
		if myid == MASTER_NODE:
		  print(' ')
		  print(">>>> Time iteration {} / FSI iteration {} <<<<".format(TimeIter,FSIIter))
		  print(' ')
		FluidSolver.resetConvergence()			#Reset the convergence flag of the fluid solution to FALSE
		#if TimeIter > TimeIterTreshold:
		FluidSolver.mapRigidDisplacementOnFluidMesh_Old()  #The displacement of the FSI interface (fluid side) is applied
		if myid == MASTER_NODE:
  		  print('\nPerforming dynamic mesh deformation (ALE)...\n')
		FluidSolver.dynamicMeshUpdate(TimeIter)		#Deform the volume mesh, calculate grid velocity, update multi-grid structure
		if myid == MASTER_NODE:
  		  print('\nLaunching fluid solver for one single dual-time iteration...')
		FluidSolver.dualTimeInnerLoop(TimeIter)		#Fluid sub-iteration (dual-time integration)
		Drag = FluidSolver.outputFluidLoads_Drag()	#Output global aero load : drag (in [N])
		Lift = FluidSolver.outputFluidLoads_Lift()	#Output global aero load : lift (in [N])
		if myid == 0 and TimeIter > TimeIterTreshold:
		  FluidSolver.outputGlobalFluidLoads(SolidSolver.getGlobalFluidLoadsArray())	#Transfer the global fluid loads (lift or drag) to the solid solver
		  SolidSolver.applyGlobalFluidLoads() 						#Apply the global fluid loads to the solid structure
		  #??.mapFluidLoadsOnSolidMesh()
		  print('\nLaunching solid solver for a single time iteration...\n')
		  SolidSolver.timeIteration(time)	#Solid sub-iteration, the current solid solution is automatically reset
		  #print("Relaxation step")
		  #SolidSolver.setAitkenCoefficient(FSIIter)		#Implemented but not functional yet
		  SolidSolver.outputDisplacements(FluidSolver.getInterRigidDispArray(), False)	#Output the displacement of the FSI interface (solid side)
		pypar.barrier()							
		FluidSolver.broadcastInterfDisp()		#Broadcast the displacement from MASTER_NODE to other processes (use of a buffer)
		pypar.barrier()
  		res = StopComp.setValue(FSIIter, Lift) 		#Calculate the residual based on fluid global loads (lift or drag)		
		if myid == MASTER_NODE:
		  print(' ')
		  print('>>>> FSI Residual : {} <<<<'.format(res))
		  print(' ')
 		if res < FSITolerance:					#Exit the FSI loop if the specified tolerance is reached
		  break	
		FSIIter += 1
	#END OF FSI LOOP
        
	FluidSolver.updateDualTime()				#Update the dual-time solution of the fluid solver
	StopIntegration = FluidSolver.writeSolution(TimeIter)	#Write the solution in the history file once the FSI loop is converged

   	if myid == MASTER_NODE and TimeIter >= TimeIterTreshold:
          SolidSolver.writeSolution(time, FSIIter, TimeIter, NbTimeIter)				#Write the solid solution in the history file once the FSI loop is converged
	  SolidSolver.displacementPredictor_Old(FluidSolver.getInterRigidDispArray())		#Compute the solid displacement predictor for the next time step
	pypar.barrier()								
        FluidSolver.broadcastInterfDisp()			#Broadcast the predicted displacement from MASTER_NODE to other processes (use of a buffer)
	pypar.barrier()
	if myid == 0 and TimeIter >= TimeIterTreshold:
	  SolidSolver.updateSolution()				#Update the solid solution for the next time-step

	TimeIter += 1
	time += deltaT
  #END OF TIME LOOP

  pypar.barrier()

  if myid == MASTER_NODE:
    print(' ')
    print('*************************')
    print('*  End FSI computation  *')
    print('*************************')
    print(' ')

#---------------------------------------------------------------
class FSIResidual:

  def __init__(self):
    self.critFSI0 = 0.0
    self.critFSIPrev = 0.0
    self.critFSI = 0.0
    self.res = 1.0
    self.logres = log10(1.0)

  def setValue(self,FSIIter,currentCritFSI):
    self.critFSI = currentCritFSI
    if FSIIter == 0:
      self.critFSI0 = currentCritFSI
      self.critFSIPrev = 0.0
    self.res = abs(self.critFSI-self.critFSIPrev)/abs(self.critFSI0)
    self.logres = log10(self.res)
    self.critFSIPrev = self.critFSI
    return self.res
