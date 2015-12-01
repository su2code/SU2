#!/usr/bin/env python
# -*-coding:utf-8 -* 

# \file LaunchFSIComputation.py
#  \brief Python wrapper code for FSI computation using C++ fluid and solid solvers
#  \author THOMAS David, University of Liege, Belgium. Aerospace and Mechanical Engineering Department
#  \version BETA

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pypar	#MPI is initialized from now by python and can be continued in C++ !

import os, sys, shutil, copy
from math import *	#Use mathematical expressions
from optparse import OptionParser	#Use a parser for configuration

import SU2	#For SU2 python parser
import FSI	#For FSI python parser


#Import the CFD and CSD modules for FSI computation
import NativeSolid
import SU2Solver

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

  myid = pypar.rank()
  numberPart = pypar.size()
  MASTER_NODE = 0

  # Command Line Options
  parser=OptionParser()
  parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
  parser.add_option("-d", "--divide_grid",dest="divide_grid",default="True",
                      help="DIVIDE_GRID the numerical grid", metavar="DIVIDE_GRID")

  (options, args)=parser.parse_args()

  confFile = str(options.filename)
  options.divide_grid = options.divide_grid.upper() == 'TRUE'

  # FSI configuration file
  FSI_config = FSI.io.FSIConfig(confFile)
  CFD_ConFile = FSI_config['CFD_CONFIG_FILE_NAME']
  CSD_ConFile = FSI_config['CSD_CONFIG_FILE_NAME']

  # Fluid solver configuration file (case of SU2)
  config = SU2.io.Config(CFD_ConFile)
  config.NUMBER_PART = numberPart
  config.DECOMPOSED = False

  # Redefine the fluid solver configuration file if we are running in parallel
  if numberPart > 1:
    config.DECOMPOSED  = True
    config.NUMBER_PART = numberPart
  else:
    config.DECOMPOSED  = False

  
  CFD_ConFile = 'config_CFD.cfg'

  if myid == MASTER_NODE:
    konfig = copy.deepcopy(config)
    konfig.dump(CFD_ConFile)

  pypar.barrier()
  

  #Compute general variables
  deltaT = FSI_config['UNST_TIMESTEP']
  totTime = FSI_config['UNST_TIME']
  NbFSIIterMax = FSI_config['NB_FSI_ITER']
  FSITolerance = FSI_config['FSI_TOLERANCE']
  TimeIterTreshold = 0
  StopIntegration = False

  if FSI_config['RESTART_SOL'] == 'YES':
    startTime = FSI_config['START_TIME']
    NbTimeIter = ((totTime)/deltaT)-1
    time = startTime
    TimeIter = FSI_config['RESTART_ITER']
  else:
    NbTimeIter = (totTime/deltaT)-1
    time = 0.0
    TimeIter = 0

  NbTimeIter = int(NbTimeIter)

#--------------------------------

  #Defition of the components of the computation
  #Each method that is available is coded in C++ into the API's
  FluidSolver = SU2Solver.SU2Solver(CFD_ConFile)
  SolidSolver = NativeSolid.NativeSolidSolver(CSD_ConFile)

  #Initialization of the solvers
  FluidSolver.initialize(True)
  if myid == MASTER_NODE:
    SolidSolver.initialize(True)
    print(' ')
    print('***************************** Connect fluid and solid solver *****************************')
    FluidSolver.connectToSolidSolver(SolidSolver.getnSolidInterfaceVertex())
    print(' ')
    print('***************************** Setting common FSI interface *****************************')
  
  pypar.barrier()
  FluidSolver.setFSIInterface()

  #Launch steady or unsteady FSI simulation (with some python errors interception)
  if FSI_config['UNSTEADY_SIMULATION'] == "YES":
    try:
      FSI.run.UnsteadyFSI(FSI_config, FluidSolver, SolidSolver, time,TimeIterTreshold, NbTimeIter, NbFSIIterMax, FSITolerance, TimeIter,deltaT)
    except NameError as exception:
      if myid == MASTER_NODE:
        print('An NameError occured in FSI.run.UnsteadyFSI : ',exception)
    except TypeError as exception:
      if myid == MASTER_NODE:
        print('A TypeError occured in FSI.run.UnsteadyFSI : ',exception)
    except KeyboardInterrupt as exception :
      if myid == MASTER_NODE:
        print('A KeyboardInterrupt occured in FSI.run.UnsteadyFSI : ',exception)
  else:
    try:
      NbExtIter = FSI_config['NB_EXT_ITER']
      FSI.run.SteadyFSI(FluidSolver, SolidSolver, NbExtIter, NbFSIIterMax, FSITolerance)
    except NameError as exception:
      if myid == MASTER_NODE:
        print('An NameError occured in FSI.run.UnsteadyFSI : ',exception)
    except TypeError as exception:
      if myid == MASTER_NODE:
        print('A TypeError occured in FSI.run.UnsteadyFSI : ',exception)
    except KeyboardInterrupt as exception :
      if myid == MASTER_NODE:
        print('A KeyboardInterrupt occured in FSI.run.UnsteadyFSI : ',exception)

  pypar.barrier()

  #Exit cleanly the fluid and solid solvers
  FluidSolver.exit()
  if myid == 0:
    SolidSolver.exit()

  pypar.finalize()	#Exit MPI

  return

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
