#!/usr/bin/env python

## \file wmles.py
#  \brief Python script for WMLES regression
#  \author E. Molina
#  \version 7.0.3 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

import os, sys
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

def channel_mesh():

    PI = 3.141592653589793
    KindElem = 12
    KindBound = 9
    nNode = 48 + 1
    mNode = 24 + 1
    lNode = 24 + 1
    xLength = 2. * PI
    yLength = 1. * PI
    zLength = 2.0

    Mesh_File = open("channel.su2","w")

    Mesh_File.write( "%\n" )
    Mesh_File.write( "% Problem dimension\n" )
    Mesh_File.write( "%\n" )
    Mesh_File.write( "NDIME=3\n" )
    Mesh_File.write( "%\n" )
    Mesh_File.write( "% Inner elements\n" )
    Mesh_File.write( "%\n" )
    Mesh_File.write( "NELEM=%s\n" % ((lNode-1)*(nNode-1)*(mNode-1)))


    iElem = 0
    for kNode in range(lNode-1):
        for jNode in range(mNode-1):
            for iNode in range(nNode-1):
                Point0 = kNode*mNode*nNode + jNode*nNode + iNode
                Point1 = kNode*mNode*nNode + jNode*nNode + iNode + 1
                Point2 = kNode*mNode*nNode + (jNode+1)*nNode + (iNode+1)
                Point3 = kNode*mNode*nNode + (jNode+1)*nNode + iNode
                Point4 = (kNode+1)*mNode*nNode + jNode*nNode + iNode
                Point5 = (kNode+1)*mNode*nNode + jNode*nNode + iNode + 1
                Point6 = (kNode+1)*mNode*nNode + (jNode + 1)*nNode + (iNode + 1)
                Point7 = (kNode+1)*mNode*nNode + (jNode + 1)*nNode + iNode
                Mesh_File.write( "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n" % (KindElem, Point0, Point1, Point2, Point3, Point4, Point5, Point6, Point7, iElem) )
                iElem = iElem + 1

    nPoint = (nNode)*(mNode)*(lNode)
    Mesh_File.write( "%\n" )
    Mesh_File.write( "NPOIN=%s\n" % ((nNode)*(mNode)*(lNode)))
    iPoint = 0
    for kNode in range(lNode):
        for jNode in range(mNode):
            for iNode in range(nNode):
                Mesh_File.write( "%15.14f \t %15.14f \t %15.14f \t %s\n" % (xLength*float(iNode)/float(nNode-1), zLength*float(kNode)/float(lNode-1), yLength*float(jNode)/float(mNode-1), iPoint) )
                iPoint = iPoint + 1

    Mesh_File.write( "%\n" )
    Mesh_File.write( "% Boundary elements\n" )
    Mesh_File.write( "%\n" )
    Mesh_File.write( "NMARK=6\n" )

    Mesh_File.write( "MARKER_TAG= lower\n" )
    elem = (nNode-1)*(mNode-1);
    Mesh_File.write( "MARKER_ELEMS=%s\n" % elem)
    for jNode in range(mNode-1):
        for iNode in range(nNode-1):
            Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, jNode*nNode + iNode, jNode*nNode + (iNode+1), (jNode + 1)*nNode + (iNode + 1), (jNode + 1)*nNode + iNode) )

    Mesh_File.write( "MARKER_TAG= upper\n" )
    elem = (nNode-1)*(mNode-1)
    Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
    for jNode in range(mNode-1):
        for iNode in range(nNode-1):
            Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, nNode*mNode*(lNode - 1) + jNode*nNode + iNode, nNode*mNode*(lNode - 1) + jNode*nNode + iNode + 1, nNode*mNode*(lNode - 1) + (jNode + 1)*nNode + (iNode + 1), nNode*mNode*(lNode - 1) + (jNode + 1)*nNode + iNode) )

    Mesh_File.write( "MARKER_TAG= left\n" )
    elem = (nNode-1)*(lNode-1)
    Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
    for iNode in range(nNode-1):
        for kNode in range(lNode-1):
            Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, iNode + kNode*nNode*mNode, iNode + (kNode+1)*nNode*mNode, iNode + 1 + (kNode+1)*nNode*mNode, iNode + 1 + kNode*nNode*mNode) )

    Mesh_File.write( "MARKER_TAG= right\n" )
    elem = (nNode-1)*(lNode-1)
    Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
    for iNode in range(nNode-1):
        for kNode in range(lNode-1):
            Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, (nNode*mNode - 1) - iNode + kNode*nNode*mNode,  (nNode*mNode - 1) - iNode + (kNode+1)*nNode*mNode, (nNode*mNode - 1) - (iNode + 1) + (kNode+1)*nNode*mNode, (nNode*mNode - 1) - (iNode + 1) + kNode*nNode*mNode) )

    Mesh_File.write( "MARKER_TAG= outlet\n" )
    elem = (mNode-1)*(lNode-1)
    Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
    for jNode in range(mNode-1):
        for kNode in range(lNode-1):
            Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, jNode*nNode + (nNode - 1) + kNode*nNode*mNode,  (jNode + 1)*nNode + (nNode - 1) + kNode*nNode*mNode,  (jNode + 1)*nNode + (nNode - 1)+ (kNode+1)*nNode*mNode, jNode*nNode + (nNode - 1)+ (kNode+1)*nNode*mNode ) )

    Mesh_File.write( "MARKER_TAG= inlet\n" )
    elem = (mNode-1)*(lNode-1)
    Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
    for jNode in range(mNode-2, -1, -1):
        for kNode in range(lNode-1):
            Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, (jNode + 1)*nNode + kNode*nNode*mNode, jNode*nNode + kNode*nNode*mNode, jNode*nNode+ (kNode+1)*nNode*mNode, (jNode + 1)*nNode+ (kNode+1)*nNode*mNode ) )

    Mesh_File.close()

    return None
# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

    # Command Line Options
    parser=OptionParser()
    parser.add_option("-f", "--file",       dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=2,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-c", "--compute",    dest="compute",    default="True",
                      help="COMPUTE direct and adjoint problem", metavar="COMPUTE")

    (options, args)=parser.parse_args()
    options.partitions  = int( options.partitions )
    options.compute     = options.compute.upper() == 'TRUE'

    if options.filename == None:
        raise Exception("No config file provided. Use -f flag")

    channel_mesh()
    parallel_computation( options.filename    ,
                          options.partitions  ,
                          options.compute      )

#: def main()


# -------------------------------------------------------------------
#  CFD Solution
# -------------------------------------------------------------------

def parallel_computation( filename           ,
                          partitions  = 0    ,
                          compute     = True  ):

    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions

    if config.SOLVER == "MULTIPHYSICS":
        print("Parallel computation script not compatible with MULTIPHYSICS solver.")
        exit(1)

    # State
    state = SU2.io.State()

    # check for existing files
    if not compute:
        state.find_files(config)
    else:
        state.FILES.MESH = config.MESH_FILENAME

    # CFD Solution (direct or adjoint)
    info = SU2.run.CFD(config)
    state.update(info)

    # Solution merging
    # if config.MATH_PROBLEM == 'DIRECT':
    #     config.SOLUTION_FILENAME = config.RESTART_FILENAME
    # elif config.MATH_PROBLEM in ['CONTINUOUS_ADJOINT', 'DISCRETE_ADJOINT']:
    #     config.SOLUTION_ADJ_FILENAME = config.RESTART_ADJ_FILENAME
    # info = SU2.run.merge(config)
    # state.update(info)
  
    return state

#: parallel_computation()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
