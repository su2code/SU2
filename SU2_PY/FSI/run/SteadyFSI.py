# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pypar
from math import *

# -------------------------------------------------------------------
#  Unsteady FSI computation 
# -------------------------------------------------------------------

def SteadyFSI(FluidSolver, SolidSolver, NbIter, NbFSIIterMax, FSITolerance):

  MASTER_NODE = 0
  myid = pypar.rank()

  if myid == MASTER_NODE:
    print(' ')
    print('*************************')
    print('* Begin FSI computation *')
    print('*************************')
    print(' ')

  Iter = 0
  StopComp = SteadyFSIResidual()
  NbFSIIter = NbFSIIterMax
  FSIIter = 0
  
  #FSI LOOP
  while FSIIter < NbFSIIterMax:
    if myid == MASTER_NODE:
      print(' ')
      print(">>>> FSI iteration {} <<<<".format(FSIIter))
      print(' ')
    FluidSolver.resetConvergence() 		#Reset the convergence flag of the fluid solution to FALSE

    #STEADY LOOP
    if myid == MASTER_NODE:
        print('\nLaunching fluid solver for a steady computation...')
    Iter = 0
    while Iter < NbIter:
      FluidSolver.setTemporalIteration(Iter)
      FluidSolver.steadyFluidIteration(Iter)		#Fluid iteration (steady iteration)
      StopIntegration = FluidSolver.writeSolution(Iter)	#Write the solution in the history file once the FSI loop is converged
      if StopIntegration:
        break;
      Iter += 1
    #END OF STEADY LOOP

    Drag = FluidSolver.outputFluidLoads_Drag()	#Output global aero load : drag (in [N])
    Lift = FluidSolver.outputFluidLoads_Lift()	#Output global aero load : lift (in [N])
    if myid == 0:
      FluidSolver.outputGlobalFluidLoads(SolidSolver.getGlobalFluidLoadsArray())	#Transfer the fluid loads (lift or drag) to the solid solver
      SolidSolver.applyGlobalFluidLoads()						#Apply the global fluid loads to the solid structure
      #??.mapFluidLoadsOnSolidMesh()
      print('\nLaunching solid solver for a static computation...\n')
      SolidSolver.staticComputation()	#Solid static computation, the current solid solution is automatically reset
      #print("Relaxation step")
      #SolidSolver.setAitkenCoefficient(FSIIter)				#Implemented but not functional yet
      FluidSolver.getSolidDisplacement(SolidSolver.getSolidInterface(), SolidSolver.getCenterCoordinate())	#Communicate the new position of the FSI interface (solid side)
      SolidSolver.writeSolution(0.0, FSIIter, Iter, NbIter)				#Write the solid solution in the history file
      SolidSolver.updateSolution()
      print('\nPerforming static mesh deformation...\n')							
    FluidSolver.mapRigidDisplacementOnFluidMesh()		#Interpolate the FSI interface displacement from solid to fluid mesh
    FluidSolver.staticMeshUpdate()		#Deform the volume mesh, update multi-grid structure			
    pypar.barrier()

    res = StopComp.setValue(FSIIter, Drag) 		#Calculate the residual based on fluid global loads		
    if myid == MASTER_NODE:
      print(' ')
      print('>>>> FSI Residual : {} <<<<'.format(res))
      print(' ')
    if res < FSITolerance:					#Exit the FSI loop if the specified tolerance os reached
      break	
    FSIIter += 1
  #END OF FSI LOOP	

  pypar.barrier()

  if myid == MASTER_NODE:
    print(' ')
    print('*************************')
    print('*  End FSI computation  *')
    print('*************************')
    print(' ')

#--------------------------------------------------------------

def SteadyFSI_Old(FluidSolver, SolidSolver, NbIter, NbFSIIterMax, FSITolerance):

  MASTER_NODE = 0
  myid = pypar.rank()

  if myid == MASTER_NODE:
    print(' ')
    print('*************************')
    print('* Begin FSI computation *')
    print('*************************')
    print(' ')

  Iter = 0
  StopComp = SteadyFSIResidual()
  NbFSIIter = NbFSIIterMax
  FSIIter = 0
  
  #FSI LOOP
  while FSIIter < NbFSIIterMax:
    if myid == MASTER_NODE:
      print(' ')
      print(">>>> FSI iteration {} <<<<".format(FSIIter))
      print(' ')
    FluidSolver.resetConvergence() 		#Reset the convergence flag of the fluid solution to FALSE
    FluidSolver.mapRigidDisplacementOnFluidMesh()  #The displacement of the FSI interface (fluid side) is applied
    if myid == MASTER_NODE:
      print('\nPerforming static mesh deformation...\n')
    FluidSolver.staticMeshUpdate()		#Deform the volume mesh, update multi-grid structure

    #STEADY LOOP
    if myid == MASTER_NODE:
        print('\nLaunching fluid solver for a steady computation...')
    Iter = 0
    while Iter < NbIter:
      FluidSolver.setTemporalIteration(Iter)
      FluidSolver.steadyFluidIteration(Iter)		#Fluid iteration (steady iteration)
      StopIntegration = FluidSolver.writeSolution(Iter)	#Write the solution in the history file once the FSI loop is converged
      if StopIntegration:
        break;
      Iter += 1
    #END OF STEADY LOOP

    Drag = FluidSolver.outputFluidLoads_Drag()	#Output global aero load : drag (in [N])
    Lift = FluidSolver.outputFluidLoads_Lift()	#Output global aero load : lift (in [N])
    if myid == 0:
      FluidSolver.outputGlobalFluidLoads(SolidSolver.getGlobalFluidLoadsArray())	#Transfer the fluid loads (lift or drag) to the solid solver
      SolidSolver.applyGlobalFluidLoads()						#Apply the global fluid loads to the solid structure
      #??.mapFluidLoadsOnSolidMesh()
      print('\nLaunching solid solver for a static computation...\n')
      SolidSolver.staticComputation()	#Solid static computation, the current solid solution is automatically reset
      #print("Relaxation step")
      #SolidSolver.setAitkenCoefficient(FSIIter)				#Implemented but not functional yet
      SolidSolver.outputDisplacements(FluidSolver.getInterRigidDispArray())	#Output the displacement of the FSI interface (solid side)
      SolidSolver.writeSolution(0.0, FSIIter)				#Write the solid solution in the history file
      SolidSolver.updateSolution()
    pypar.barrier()							
    FluidSolver.broadcastInterfDisp()				#Broadcast the displacement from MASTER_NODE to other processes (use of a buffer)
    pypar.barrier()

    res = StopComp.setValue(FSIIter, Drag) 		#Calculate the residual based on fluid global loads		
    if myid == MASTER_NODE:
      print(' ')
      print('>>>> FSI Residual : {} <<<<'.format(res))
      print(' ')
    if res < FSITolerance:					#Exit the FSI loop if the specified tolerance os reached
      break	
    FSIIter += 1
  #END OF FSI LOOP	

  pypar.barrier()

  if myid == MASTER_NODE:
    print(' ')
    print('*************************')
    print('*  End FSI computation  *')
    print('*************************')
    print(' ')

#--------------------------------------------------------------

class SteadyFSIResidual:

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
