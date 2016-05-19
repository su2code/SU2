#!/usr/bin/env python

## \file set_ffd_design_var.py.py
#  \brief Python script for automatically generating a list of FFD variables.
#  \author T. Economon, F. Palacios
#  \version 4.1.2 "Cardinal"
#
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

import os, time
from optparse import OptionParser
from numpy import *

parser = OptionParser()
parser.add_option("-i", "--iDegree", dest="iDegree", default=4,
                  help="i degree of the FFD box", metavar="IDEGREE")
parser.add_option("-j", "--jDegree", dest="jDegree", default=4,
                  help="j degree of the FFD box", metavar="JDEGREE")
parser.add_option("-k", "--kDegree", dest="kDegree", default=1,
                  help="k degree of the FFD box", metavar="KDEGREE")
parser.add_option("-b", "--ffdid", dest="ffd_id", default=0,
                  help="ID of the FFD box", metavar="FFD_ID")
parser.add_option("-m", "--marker", dest="marker",
                  help="marker name of the design surface", metavar="MARKER")
parser.add_option("-s", "--scale", dest="scale", default=1.0,
                  help="scale factor for the bump functions", metavar="SCALE")
parser.add_option("-d", "--dimension", dest="dimension", default=3.0,
                  help="dimension of the problem", metavar="DIMENSION")
parser.add_option("-x", "--xMovement", dest="xMove", default=0.0,
                  help="movement in x direction", metavar="XMOVEMENT")
parser.add_option("-y", "--yMovement", dest="yMove", default=0.0,
                  help="movement in y direction", metavar="YMOVEMENT")
parser.add_option("-z", "--zMovement", dest="zMove", default=0.0,
                  help="movement in z direction", metavar="ZMOVEMENT")


(options, args)=parser.parse_args()

# Process options
options.iOrder  = int(options.iDegree) + 1
options.jOrder  = int(options.jDegree) + 1
options.kOrder  = int(options.kDegree) + 1
options.ffd_id  = str(options.ffd_id)
options.marker = str(options.marker)
options.scale  = float(options.scale)
options.dim  = int(options.dimension)
options.xMove = float(options.xMove)
options.yMove = float(options.yMove)
options.zMove = float(options.zMove)

if options.dim == 3:
  
  print " "
  print "FFD_CONTROL_POINT"

  iVariable = 0
  dvList = "DEFINITION_DV= "
  for kIndex in range(options.kOrder):
    for jIndex in range(options.jOrder):
      for iIndex in range(options.iOrder):
        iVariable = iVariable + 1
        dvList = dvList + "( 7, " + str(options.scale) + " | " + options.marker + " | "
        dvList = dvList + options.ffd_id + ", " + str(iIndex) + ", " + str(jIndex) + ", " + str(kIndex) + ", " + str(options.xMove) +", " + str(options.yMove) + ", " + str(options.zMove) + " )"
        if iVariable < (options.iOrder*(options.jOrder)*options.kOrder):
          dvList = dvList + "; "


  print dvList

  print " "
  print "FFD_CAMBER & FFD_THICKNESS"

  iVariable = 0
  dvList = "DEFINITION_DV= "
  for jIndex in range(options.jOrder):
    for iIndex in range(options.iOrder):
      iVariable = iVariable + 1
      dvList = dvList + "( 11, " + str(options.scale) + " | " + options.marker + " | "
      dvList = dvList + options.ffd_id + ", " + str(iIndex) + ", " + str(jIndex) + " )"
      dvList = dvList + "; "
  iVariable = 0
  for jIndex in range(options.jOrder):
    for iIndex in range(options.iOrder):
      iVariable = iVariable + 1
      dvList = dvList + "( 12, " + str(options.scale) + " | " + options.marker + " | "
      dvList = dvList + options.ffd_id + ", " + str(iIndex) + ", " + str(jIndex) + " )"
      if iVariable < (options.iOrder*(options.jOrder)):
        dvList = dvList + "; "

  print dvList

if options.dim == 2:

  iVariable = 0
  dvList = "DEFINITION_DV= "
  for jIndex in range(options.jOrder):
    for iIndex in range(options.iOrder):
      iVariable = iVariable + 1
      dvList = dvList + "( 15, " + str(options.scale) + " | " + options.marker + " | "
      dvList = dvList + options.ffd_id + ", " + str(iIndex) + ", " + str(jIndex) + ", " + str(options.xMove) +", " + str(options.yMove) + " )"
      if iVariable < (options.iOrder*options.jOrder):
        dvList = dvList + "; "

  print dvList

  print " "
  print "FFD_CAMBER & FFD_THICKNESS"

  iVariable = 0
  dvList = "DEFINITION_DV= "
  for iIndex in range(options.iOrder):
    iVariable = iVariable + 1
    dvList = dvList + "( 16, " + str(options.scale) + " | " + options.marker + " | "
    dvList = dvList + options.ffd_id + ", " + str(iIndex) + " )"
    dvList = dvList + "; "
  iVariable = 0
  for iIndex in range(options.iOrder):
    iVariable = iVariable + 1
    dvList = dvList + "( 17, " + str(options.scale) + " | " + options.marker + " | "
    dvList = dvList + options.ffd_id + ", " + str(iIndex) + " )"
    if iVariable < (options.iOrder):
      dvList = dvList + "; "

  print dvList





