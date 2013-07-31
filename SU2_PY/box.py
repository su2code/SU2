#!/usr/bin/env python 

## \file box.py
#  \brief Python script for box meshing
#  \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.6
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from optparse import OptionParser

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename", default="channel.su2",
                  help="write mesh to FILE", metavar="FILE")
parser.add_option("-n", "--nNode", dest="nNode", default=5,
                  help="use this NNODE in x direction", metavar="NNODE")
parser.add_option("-m", "--mNode", dest="mNode", default=5,
                  help="use this MNODE in y direction", metavar="MNODE")
parser.add_option("-l", "--lNode", dest="lNode", default=5,
                  help="use this LNODE in z direction", metavar="LNODE")
parser.add_option("-x", "--xLength", dest="xLength", default=1.0,
                  help="use this XLENGTH", metavar="XLENGTH")
parser.add_option("-y", "--yLength", dest="yLength", default=1.0,
                  help="use this YLENGTH", metavar="YLENGTH")
parser.add_option("-z", "--zLength", dest="zLength", default=1.0,
                  help="use this ZLENGTH", metavar="ZLENGTH")
(options, args)=parser.parse_args()

KindElem = 12
KindBound = 9
nNode = int(options.nNode)
mNode = int(options.mNode)
lNode = int(options.lNode)
xLength = float(options.xLength)
yLength = float(options.yLength)
zLength = float(options.zLength)

Mesh_File = open(options.filename,"w")

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

Mesh_File.write( "MARKER_TAG= left\n" )
elem = (nNode-1)*(mNode-1);
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem)
for jNode in range(mNode-1):
    for iNode in range(nNode-1):
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, jNode*nNode + iNode, jNode*nNode + (iNode+1), (jNode + 1)*nNode + (iNode + 1), (jNode + 1)*nNode + iNode) )

Mesh_File.write( "MARKER_TAG= right\n" )
elem = (nNode-1)*(mNode-1)
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
for jNode in range(mNode-1):
    for iNode in range(nNode-1):
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, nNode*mNode*(lNode - 1) + jNode*nNode + iNode, nNode*mNode*(lNode - 1) + jNode*nNode + iNode + 1, nNode*mNode*(lNode - 1) + (jNode + 1)*nNode + (iNode + 1), nNode*mNode*(lNode - 1) + (jNode + 1)*nNode + iNode) )

Mesh_File.write( "MARKER_TAG= lower\n" )
elem = (nNode-1)*(lNode-1)
Mesh_File.write( "MARKER_ELEMS=%s\n" % elem )
for iNode in range(nNode-1):
    for kNode in range(lNode-1):
        Mesh_File.write( "%s \t %s \t %s \t %s \t %s\n" % (KindBound, iNode + kNode*nNode*mNode, iNode + (kNode+1)*nNode*mNode, iNode + 1 + (kNode+1)*nNode*mNode, iNode + 1 + kNode*nNode*mNode) )

Mesh_File.write( "MARKER_TAG= upper\n" )
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

Mesh_File.write( "NCHUNK=1\n")
Mesh_File.write( "CHUNK_TAG=0\n")
Mesh_File.write( "CHUNK_DEGREE_I=6\n")
Mesh_File.write( "CHUNK_DEGREE_J=6\n")
Mesh_File.write( "CHUNK_DEGREE_K=1\n")
Mesh_File.write( "CHUNK_CORNER_POINTS=8\n")
Mesh_File.write( "4.0	0	-0.1\n")
Mesh_File.write( "6.0	0	-0.1\n")
Mesh_File.write( "6.0	2.0	-0.1\n")
Mesh_File.write( "4.0	2.0	-0.1\n")
Mesh_File.write( "4.0	0	0.1\n")
Mesh_File.write( "6.0	0	0.1\n")
Mesh_File.write( "6.0	2.0	0.1\n")
Mesh_File.write( "4.0	2.0	0.1\n")
Mesh_File.write( "CHUNK_CONTROL_POINTS=0\n")
Mesh_File.write( "CHUNK_SURFACE_POINTS=0\n")


    
Mesh_File.close()
