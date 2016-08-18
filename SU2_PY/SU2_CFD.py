#!/usr/bin/env python

## \file test_pyWrapper.py
#  \brief Python script to launch SU2_CFD through the Python Wrapper
#  \author David Thomas
#  \version 4.2.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

from optparse import OptionParser	# use a parser for configuration
import SU2				# imports SU2 python tools

import SU2Solver			# imports the SU2 wrapped module

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
  parser.add_option("--fsi", dest="fsi", action="store_true", default="False", 
                    help="Launch the FSI driver", metavar="FSI")
  parser.add_option("--spectral", action="store_true", dest="time_spectral", default="False",
                    help="Launch the time SPECTRAL driver", metavar="SPECTRAL")
  parser.add_option("--parallel", action="store_true",
                    help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  (options, args) = parser.parse_args()
  options.nDim  = int( options.nDim )
  options.nZone = int( options.nZone )

  if options.filename == None:
    raise Exception("No config file provided. Use -f flag")

  if options.with_MPI == True:
    from mpi4py import MPI			# use mpi4py for parallel run (also valid for serial)
    comm = MPI.COMM_WORLD
  else:
    comm = 0 

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  if options.nZone == 1:
    SU2Driver = SU2Solver.CSingleZoneDriver(options.filename, options.nZone, options.nDim, comm);
  elif options.time_spectral == True:
    print options.time_spectral
    SU2Driver = SU2Solver.CSpectralDriver(options.filename, options.nZone, options.nDim, comm);
  elif (options.nZone == 2) and (options.fsi == True):
    SU2Driver = SU2Solver.CFSIDriver(options.filename, options.nZone, options.nDim, comm);
  else:
    SU2Driver = SU2Solver.CMultiZoneDriver(options.filename, options.nZone, options.nDim, comm);

  # Launch the solver for the entire computation
  SU2Driver.StartSolver()

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
