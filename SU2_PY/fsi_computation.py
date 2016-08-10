#!/usr/bin/env python
# -*-coding:utf-8 -* 

# \file LaunchFSIComputation.py
#  \brief Python wrapper code for FSI computation using C++ fluid and solid solvers.
#  \author D. THOMAS, University of Liege, Belgium. Department of Mechanical and Aerospace Engineering.
#  \version BETA

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from mpi4py import MPI  # MPI is initialized from now by python and can be continued in C++ !

import os, sys, shutil, copy
import time as timer
from math import *	# use mathematical expressions
from optparse import OptionParser	# use a parser for configuration

import SU2	# imports SU2 python tools
import FSI	# imports FSI python tools


# imports the CFD (SU2) module for FSI computation
import SU2Solver

# -------------------------------------------------------------------
#  Main 
# -------------------------------------------------------------------

def main():

  comm = MPI.COMM_WORLD
  myid = comm.Get_rank()
  numberPart = comm.Get_size()
 
  rootProcess = 0

  # --- Set the working directory --- #
  if myid == rootProcess:
      if os.getcwd() not in sys.path:
          sys.path.append(os.getcwd())
	  print("Setting working directory : {}".format(os.getcwd()))
      else: 
	  print ("Working directory is set to {}".format(os.getcwd()))

  # starts timer
  start = timer.time()

  # --- Get the FSI conig file name form the command line options --- #
  parser=OptionParser()
  parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
  parser.add_option("-d", "--divide_grid",dest="divide_grid",default="True",
                      help="DIVIDE_GRID the numerical grid", metavar="DIVIDE_GRID")

  (options, args)=parser.parse_args()

  confFile = str(options.filename)
  options.divide_grid = options.divide_grid.upper() == 'TRUE'

  FSI_config = FSI.io.FSIConfig(confFile) 		# FSI configuration file
  CFD_ConFile = FSI_config['CFD_CONFIG_FILE_NAME']	# CFD configuration file
  CSD_ConFile = FSI_config['CSD_CONFIG_FILE_NAME']	# CSD configuration file

  CSD_Solver = FSI_config['CSD_SOLVER']			# CSD solver

  # Redefine the fluid solver configuration file if we are running in parallel
  config = SU2.io.Config(CFD_ConFile)
  config.NUMBER_PART = numberPart
  config.DECOMPOSED = False
  if numberPart > 1:
    config.DECOMPOSED  = True
    config.NUMBER_PART = numberPart
  else:
    config.DECOMPOSED  = False
  CFD_ConFile = 'config_CFD.cfg'
  if myid == rootProcess:
    konfig = copy.deepcopy(config)
    konfig.dump(CFD_ConFile)

  comm.barrier()

  # --- Initialize the fluid solver --- #
  if myid == rootProcess:
    print('\n***************************** Initializing fluid solver *****************************')
  FluidSolver = SU2Solver.CSingleZoneDriver(CFD_ConFile, 1, 2)

  comm.barrier()
  
  # --- Initialize the solid solver --- # (!! for now we are using only serial solid solvers)
  if myid == rootProcess:
    print('\n***************************** Initializing solid solver *****************************')
    if CSD_Solver == 'METAFOR':
      from MetaforSolver import MtfSolver
      SolidSolver = MtfSolver(CSD_ConFile)
    elif CSD_Solver == 'NATIVE':
      import NativeSolid
      SolidSolver = NativeSolid.NativeSolidSolver(CSD_ConFile, True)
    elif CSD_Solver == 'GETDP':
      import GetDPSolver
      SolidSolver = GetDPSolver.GetDPSolver(CSD_ConFile, True)
  else:
    SolidSolver = None

  comm.barrier()

  # --- Initialize and set the FSI interface (coupling environement) --- #
  if myid == rootProcess:
    print('\n***************************** Initializing FSI interface *****************************')
  comm.barrier()
  FSIInterface = FSI.Interface(FSI_config, FluidSolver, SolidSolver)
  
  if myid == rootProcess:
    print('\n***************************** Connect fluid and solid solvers *****************************')
  comm.barrier()
  FSIInterface.connect(FluidSolver, SolidSolver)

  if myid == rootProcess:
    print('\n***************************** Mapping fluid-solid interfaces *****************************')
  comm.barrier()
  FSIInterface.interfaceMapping(FluidSolver, SolidSolver, FSI_config)
  
  comm.barrier()

  # --- Launch a steady or unsteady FSI computation --- #
  if FSI_config['UNSTEADY_SIMULATION'] == "YES":
    try:
      FSIInterface.UnsteadyFSI(FSI_config, FluidSolver, SolidSolver)
    except NameError as exception:
      if myid == rootProcess:
        print('An NameError occured in FSIInterface.UnsteadyFSI : ',exception)
    except TypeError as exception:
      if myid == rootProcess:
        print('A TypeError occured in FSIInterface.UnsteadyFSI : ',exception)
    except KeyboardInterrupt as exception :
      if myid == rootProcess:
        print('A KeyboardInterrupt occured in FSIInterface.UnsteadyFSI : ',exception)
  else:
    try:
      NbExtIter = FSI_config['NB_EXT_ITER']
      FSIInterface.SteadyFSI(FSI_config, FluidSolver, SolidSolver)
    except NameError as exception:
      if myid == rootProcess:
        print('An NameError occured in FSIInterface.SteadyFSI : ',exception)
    except TypeError as exception:
      if myid == rootProcess:
        print('A TypeError occured in FSIInterface.SteadyFSI : ',exception)
    except KeyboardInterrupt as exception :
      if myid == rootProcess:
        print('A KeyboardInterrupt occured in FSIInterface.SteadyFSI : ',exception)
  
  comm.barrier()

  # --- Exit cleanly the fluid and solid solvers --- #
  FluidSolver.Postprocessing()
  if myid == rootProcess:
      SolidSolver.exit()

  comm.barrier()

  # stops timer
  stop = timer.time()
  elapsedTime = stop-start
  
  if myid == rootProcess:
    print("\n Computation successfully performed in {} seconds.".format(elapsedTime))

  return

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# --- This is only accessed if running from command prompt --- #
if __name__ == '__main__':
    main()
