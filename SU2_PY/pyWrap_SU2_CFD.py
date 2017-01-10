#!/usr/bin/env python

## \file test_pyWrapper.py
#  \brief Python script to launch SU2_CFD through the Python Wrapper
#  \author David Thomas
#  \version 4.3.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
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

from mpi4py import MPI			# use mpi4py for parallel run (also valid for serial) 
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

  parser.add_option("--fsi", dest="fsi", default="False", help="Launch the FSI driver", metavar="FSI")

  parser.add_option("--fem", dest="fem", default="False", help="Launch the FEM driver (General driver)", metavar="FEM")

  parser.add_option("--harmonic_balance", dest="harmonic_balance", default="False",
                    help="Launch the Harmonic Balance (HB) driver", metavar="HB")

  parser.add_option("--poisson_equation", dest="poisson_equation", default="False",
                    help="Launch the poisson equation driver (General driver)", metavar="POIS_EQ")

  parser.add_option("--wave_equation", dest="wave_equation", default="False",
                    help="Launch the wave equation driver (General driver)", metavar="WAVE_EQ")

  parser.add_option("--heat_equation", dest="heat_equation", default="False",
                    help="Launch the heat equation driver (General driver)", metavar="HEAT_EQ")

  (options, args) = parser.parse_args()
  options.nDim  = int( options.nDim )
  options.nZone = int( options.nZone )
  options.fsi = options.fsi.upper() == 'TRUE'
  options.fem = options.fem.upper() == 'TRUE'
  options.harmonic_balance = options.harmonic_balance.upper() == 'TRUE'
  options.poisson_equation = options.poisson_equation.upper() == 'TRUE'
  options.wave_equation    = options.wave_equation.upper()    == 'TRUE'
  options.heat_equation    = options.heat_equation.upper()    == 'TRUE'

  if options.filename == None:
    raise Exception("No config file provided. Use -f flag")

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  if (options.nZone == 1) and ( options.fem or options.poisson_equation or options.wave_equation or options.heat_equation ):
    SU2Driver = SU2Solver.CGeneralDriver(options.filename, options.nZone, options.nDim);
  elif options.harmonic_balance:
    SU2Driver = SU2Solver.CHBDriver(options.filename, options.nZone, options.nDim);
  elif (options.nZone == 2) and (options.fsi):
    SU2Driver = SU2Solver.CFSIDriver(options.filename, options.nZone, options.nDim);
  else:
    SU2Driver = SU2Solver.CFluidDriver(options.filename, options.nZone, options.nDim);

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
