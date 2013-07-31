#!/usr/bin/env python 

## \file tools.py
#  \brief mesh functions
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
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

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import numpy as np

# ---------------------------------------------------------------------- 
#  Read SU2 Mesh File
# ---------------------------------------------------------------------- 
def read(filename,scale=1.0):
    ''' imports mesh and builds python dictionary structure 
        input: filename
               scale: apply scaling factor (optional)
        output:
           meshdata            mesh data dictionary
           meshdata['NDIME']   number of dimensions
           meshdata['NELEM']   number of elements
           meshdata['ELEM']    element array [ type, nodes, index ]
           meshdata['NPOIN']   number of points
           meshdata['NMARK']   number of markers
           meshdata['MARKS']   marker data dictionary
           meshdata['MARKS']['tag_name']           marker data for 'tag_name'
           meshdata['MARKS']['tag_name']['NELEM']  number of elements
           meshdata['MARKS']['tag_name']['ELEM']   element array [type,nodes]
    '''

    # initialize variables
    data  = {} 
    marks = {}

    # open meshfile
    meshfile = open(filename,'r')

    # scan file until end of file
    keepon = True
    while keepon:

        # read line
        line = meshfile.readline()
        line = line.replace('\t',' ')
        line = line.replace('\n',' ')
        
        # stop if line is empty
        if not line: 
            keepon = False

        # skip comments
        elif line[0] == "%":
            pass

        # number of dimensions
        elif "NDIME=" in line:
            # save to SU2_MESH data
            data['NDIME'] = int( line.split("=")[1].strip() )
        #:if NDIME

        # elements
        elif "NELEM=" in line:
            # number of elements
            nelem = long( line.split("=")[1].strip() )
            # save to SU2_MESH data
            data['NELEM'] = nelem
            
            # element data list
            elem = []

            # scan next lines for element data
            for ielem in range(nelem):
                # read line
                line = meshfile.readline()
                # split line, convert to long ints
                thiselem = map(long, line.split() )
                # add to element list
                elem = elem + [thiselem]

            # save to SU2_MESH data
            data['ELEM'] = elem
        #: if NELEM

        # points
        elif "NPOIN=" in line:
            npoin = long( line.split("=")[1].strip().split(' ')[0] )
            # save to SU2_MESH data
            data['NPOIN'] = npoin

            # point data list
            poin = []
            
            #scan next lines for point data
            for ipoin in range(npoin):
                # read line
                line = meshfile.readline()
                # split line, convert to long ints
                thispoin = [float(x)*scale for x in line.split() ]
                # add to point list
                poin = poin + [thispoin]

            # save to SU2_MESH data
            data['POIN'] = poin
        #:if NPOIN

        # number of markers
        elif "NMARK=" in line:
            nmark = long( line.split("=")[1].strip() )
            # save to SU2_MESH data
            data['NMARK'] = nmark
        #:if NMARK

        # a marker
        elif "MARKER_TAG=" in line:
            # marker tag
            thistag = line.split("=")[1].strip()
            # start SU2_MARK dictionary
            thismark = {} 
            # save to SU2_MARK data
            thismark['TAG'] = thistag

            # read number of marker elements
            line = meshfile.readline()
            if not "MARKER_ELEMS=" in line:
                raise Exception("Marker Specification Error")
            # convert string to long int
            thisnelem = long( line.split("=")[1].strip() )
            # save to SU2_MARK data
            thismark['NELEM'] = thisnelem

            # marker element data list
            markelem = []
            # scan next lines for marker elements
            for ielem in range(thisnelem):
                # read line
                line = meshfile.readline()
                # split and convert to long ints
                thiselem = map(long, line.split() )
                # add to marker element list
                markelem = markelem + [thiselem]
            # save to SU2_MARK data
            thismark['ELEM'] = markelem
            # add to marker list
            marks[thismark['TAG']] = thismark
        #:if MARKER_TAG

    #:while not end of file

    # save to SU2_MESH data
    data['MARKS'] = marks
    
    return data
#: def read


# ---------------------------------------------------------------------- 
#  Write SU2 Mesh File
# ---------------------------------------------------------------------- 
def write(filename,meshdata,scale=1.0):
    ''' writes meshdata to file
        inputs: filename, meshdata 
    '''

    # open file for writing
    outputfile = open(filename,'w')

    # numbers
    ndime = meshdata['NDIME']

    # write dimension
    outputfile.write("% \n% Problem Dimension \n% \n")
    outputfile.write("NDIME= %i\n" % meshdata['NDIME'])

    # write elements
    outputfile.write("% \n% Inner element connectivity \n% \n")
    outputfile.write("NELEM= %i\n" % meshdata['NELEM'])
    for elem in meshdata['ELEM']:
        for num in elem:
            outputfile.write("%i " % num)
        outputfile.write("\n")

    # write nodes
    outputfile.write("% \n% Node coordinates \n% \n")
    outputfile.write("NPOIN= %i\n" % meshdata['NPOIN'])
    for poin in meshdata['POIN']:
        for inum in range(ndime):
            outputfile.write("%#18.10e " % (poin[inum]*scale))
        outputfile.write( "%i\n" % (long(poin[inum+1])) )

    # write markers 
    outputfile.write("% \n% Boundary elements \n% \n")
    outputfile.write( "NMARK= %i\n" % meshdata['NMARK'] )
    for mark_tag in meshdata['MARKS'].keys():
        this_mark = meshdata['MARKS'][mark_tag]
        outputfile.write( "MARKER_TAG= %s\n" % this_mark['TAG'] )
        outputfile.write( "MARKER_ELEMS= %i\n" % this_mark['NELEM'] )
        for elem in this_mark['ELEM']:
            for num in elem:
                outputfile.write("%i " % num)
            outputfile.write("\n")

    # close file
    outputfile.close()

    return
#: def write


# ---------------------------------------------------------------------- 
#  Get Marker Mesh Points
# ---------------------------------------------------------------------- 
def get_markerPoints(meshdata,mark_tags):
    ''' pulls all mesh nodes on markers 
        checks for duplicates (from edges) '''

    # marker tags should be a list
    if not isinstance(mark_tags,list):
        mark_tags = [mark_tags]

    # some numbers
    nmark = meshdata['NMARK']
    ndim  = meshdata['NDIME']
    
    # list for marker node numbers
    markernodes  = []

    # scan each marker
    for this_tag in mark_tags:
        # current mark
        this_mark = meshdata['MARKS'][this_tag]
        # marker elements
        markelems = this_mark['ELEM']
        # list for marker nodes
        marknodes = []
        # pull all marker nodes, there will be duplicates
        for row in markelems:
            # ignore first marker element entry (element type)
            marknodes = marknodes + row[1:]
        # find unique node points
        marknodes = dict(map(lambda i:(i,1),marknodes)).keys()
        # add to mesh node list
        markernodes  = markernodes  + marknodes
    #: for each marker

    # one more unique check
    markernodes = dict(map(lambda i:(i,1),markernodes)).keys()

    # list for marker points
    markerpoints = []
    # pull all nodes on markers
    for inode in markernodes:
        markerpoints = markerpoints + [ tuple( meshdata['POIN'][inode][0:ndim] ) ]

    return markerpoints, markernodes

#: def get_markerPoints()


# ---------------------------------------------------------------------- 
#  Set Mesh Points
# ---------------------------------------------------------------------- 
def set_meshPoints(meshdata,meshnodes,meshpoints):
    ''' stores array of meshpoints in the meshdata structure
        note: will operate on the input meshdata by pointer
              if a new mesh is needed make a deep copy 
              before calling this function
    '''

    n_nodes = len(meshnodes)
    n_dim   = meshdata['NDIME']

    # for each given node, update meshdata['POIN']
    for ipoint in range(n_nodes):
        inode = meshnodes[ipoint]
        for iDim in range(n_dim):
            meshdata['POIN'][inode][iDim] = meshpoints[ipoint][iDim]

    return meshdata

#def: set_meshPoints

# ---------------------------------------------------------------------- 
#  Sort Airfoil
# ---------------------------------------------------------------------- 
def sort_airfoil(mesh_data,marker_name):
    ''' sorts xy airfoil points in clockwise loop from trailing edge 
        returns list of mesh point indeces
        assumes:
          - airfoil oriented nearly parallel with x-axis
          - oriented from leading to trailing edge in the +x-direction
          - one airfoil element with name 'marker_name' 
    '''
    
    # find airfoil elements and points
    airfoil_elems  = mesh_data['MARKS'][marker_name]['ELEM']
    airfoil_elems  = np.array(airfoil_elems)    
    airfoil_points = mesh_data['POIN']
    airfoil_points = np.array(airfoil_points)
    airfoil_points = airfoil_points[airfoil_elems[:,1],:]
    n_P,_ = airfoil_elems.shape
    
    # build transfer arrays
    EP = airfoil_elems[:,1:3]    # edge to point
    PX = airfoil_points[:,0:2]   # point to coord
    IP = np.arange(0,n_P)        # loop index to point
    
    # sorted airfoil point indeces tobe
    Psort = np.zeros(n_P,long)
    Isort = np.arange(0,n_P)
    
    # find trailing edge
    iP0 = np.argmax(PX[:,0])
    P0  = EP[iP0,0]
    I0  = IP[iP0]
    Psort[0] = P0
    
    # build loop
    for this_iP in range(1,n_P):
        P0 = EP[EP[:,0]==P0,1]
        I0 = IP[EP[:,0]==P0]
        Psort[this_iP] = P0
        Isort[this_iP] = I0
    
    
    # check for clockwise
    D1 = PX[Isort[1],1]  - PX[Isort[0],1]
    D2 = PX[Isort[-1],1] - PX[Isort[0],1]
    if D1>D2:   
        Psort = Psort[-1::-1]
    
    # done
    points_sorted = Psort
    loop_sorted   = Isort
    
    return points_sorted,loop_sorted


#def: sort_airfoil()
