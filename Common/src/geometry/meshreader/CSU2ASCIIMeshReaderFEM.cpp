/*!
 * \file CSU2ASCIIMeshReaderFVM.cpp
 * \brief Reads a native SU2 ASCII grid into linear partitions for the
 *        finite element solver (FEM).
 * \author T. Economon, E. van der Weide
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../../include/geometry/meshreader/CSU2ASCIIMeshReaderFEM.hpp"
#include "../../../include/fem/CFEMStandardElementBase.hpp"

CSU2ASCIIMeshReaderFEM::CSU2ASCIIMeshReaderFEM(CConfig        *val_config,
                                               unsigned short val_iZone,
                                               unsigned short val_nZone)
: CSU2ASCIIMeshReaderBase(val_config, val_iZone, val_nZone) {
  
  /* Read the basic metadata and perform some basic error checks. */
  ReadMetadata(true, val_config);
  
  /*--- Read the volume connectivity and distribute it
        linearly over the MPI ranks. ---*/
  ReadVolumeElementConnectivity();

  /*--- Read the coordinates of the points that are needed
        on this MPI rank. ---*/
  ReadPointCoordinates();

  /*--- Read the surface connectivity and store the surface elements whose
        corresponding volume element is stored on this MPI rank. ---*/
  ReadSurfaceElementConnectivity();
}

CSU2ASCIIMeshReaderFEM::~CSU2ASCIIMeshReaderFEM(void) { }

void CSU2ASCIIMeshReaderFEM::ReadPointCoordinates() {

  /*--- Loop over the local elements to determine the global
        point IDs to be stored on this rank. --*/
  unsigned long ind = 0;
  for(unsigned long i=0; i<numberOfLocalElements; ++i) {

    /*--- Store the number of grid DOFs for this element and
          skip the meta data for this element (5 entries). ---*/
    const unsigned long nDOFsGrid = localVolumeElementConnectivity[ind+3];
    ind += 5;

    /*--- Copy the connectivity to globalPointIDs. ---*/
    unsigned long *conn = localVolumeElementConnectivity.data() + ind;
    ind += nDOFsGrid;
    globalPointIDs.insert(globalPointIDs.end(), conn, conn+nDOFsGrid);
  }

  /*--- Sort the globalPointIDs and remove the duplicate entries. ---*/
  sort(globalPointIDs.begin(), globalPointIDs.end());
  vector<unsigned long>::iterator lastNode;
  lastNode = unique(globalPointIDs.begin(), globalPointIDs.end());
  globalPointIDs.erase(lastNode, globalPointIDs.end());

  /*--- Determine the number of locally stored points. ---*/
  numberOfLocalPoints = globalPointIDs.size();

  /*--- Prepare our data structure for the point coordinates. ---*/
  localPointCoordinates.resize(dimension);
  for(int k=0; k<dimension; ++k)
    localPointCoordinates[k].reserve(numberOfLocalPoints);

  /*--- Open the mesh file and jump to our zone. ---*/
  mesh_file.open(meshFilename, ios::in);
  FastForwardToMyZone();

  /*--- Find the section containing the coordinates. ---*/
  string text_line;
  while (getline (mesh_file, text_line)) {

    string::size_type position = text_line.find ("NPOIN=",0);
    if (position != string::npos) break;
  }

  /*--- Loop over the global number of points in the grid. ---*/
  for(unsigned long i=0; i<numberOfGlobalPoints; ++i) {

    /*--- Read the line in the grid file. This line must always
          be read, even if the point is not stored on this rank. ---*/
    getline(mesh_file, text_line);

    /*--- Determine whether this point must be stored on this rank. ---*/
    if( binary_search(globalPointIDs.begin(), globalPointIDs.end(), i) ) {

      /*--- Read the coordinates from the string and store them
            in localPointCoordinates. ---*/
      istringstream point_line(text_line);
      for (unsigned short iDim = 0; iDim < dimension; ++iDim) {
        passivedouble Coord;
        point_line >> Coord;
        localPointCoordinates[iDim].push_back(Coord);
      }
    }
  }

  /*--- Close the mesh file again. ---*/
  mesh_file.close();
}

void CSU2ASCIIMeshReaderFEM::ReadVolumeElementConnectivity() {

  /* Get a partitioner to help with linear partitioning. */
  CLinearPartitioner elemPartitioner(numberOfGlobalElements,0);

  /*--- Open the mesh file and jump to our zone. ---*/
  mesh_file.open(meshFilename, ios::in);
  FastForwardToMyZone();

  /*--- Find the section containing the interior elements. ---*/
  string text_line;
  while (getline (mesh_file, text_line)) {

    string::size_type position = text_line.find ("NELEM=",0);
    if (position != string::npos) break;
  }

  /*--- Skip the elements, which are read by ranks lower than my rank. ---*/
  const unsigned long firstIndex = elemPartitioner.GetFirstIndexOnRank(rank);
  for(unsigned long i=0; i<firstIndex; ++i) getline(mesh_file, text_line);

  /*--- Determine the number of local elements. ---*/
  numberOfLocalElements = elemPartitioner.GetSizeOnRank(rank);

  /*--- Loop over the elements that must be stored on this rank. ---*/
  for(unsigned long i=0; i<numberOfLocalElements; ++i) {

    /*--- Read the line for this element and couple it to an istringstream
          to enable the actual reading of the data. ---*/
    getline(mesh_file, text_line);
    istringstream elem_line(text_line);

    /*--- Read the value that defines the element type and the polynomial degree
          of the geometry and the solution. Extract this info as well. ---*/
    unsigned long typeRead; elem_line >> typeRead;
    unsigned long typeReadErrorMessage = typeRead;

    unsigned short nPolySol, nPolyGrid;
    if(typeRead > 10000) {
      nPolySol  = typeRead/10000 -1;
      typeRead  = typeRead%10000;
      nPolyGrid = typeRead/100 + 1;
    }
    else {
      nPolyGrid = typeRead/100 + 1;
      nPolySol  = nPolyGrid;
    }

    unsigned short VTK_Type = typeRead%100;

    /*--- Determine the number of grid DOFs for this element. ---*/
    unsigned short nDOFsGrid = CFEMStandardElementBase::GetNDOFsStatic(VTK_Type, nPolyGrid);
    if(nDOFsGrid == 0) {
      ostringstream message;
      message << "Unknown FEM element type, " << typeReadErrorMessage
              << ", encountered.";
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
    }

    /*--- Store the data just created in localVolumeElementConnectivity. ---*/
    localVolumeElementConnectivity.push_back(VTK_Type);
    localVolumeElementConnectivity.push_back(nPolyGrid);
    localVolumeElementConnectivity.push_back(nPolySol);
    localVolumeElementConnectivity.push_back(nDOFsGrid);
    localVolumeElementConnectivity.push_back(firstIndex+i);  // Global elem ID.

    /*--- Store the current index, as this is needed when the
          connectivity of linear elements is swapped. ---*/
    const unsigned long ind = localVolumeElementConnectivity.size();

    /*--- Read the connectivity from elem_line and store it in
          localVolumeElementConnectivity. ---*/
    for(unsigned short j=0; j<nDOFsGrid; ++j) {
      unsigned long nodeID;
      elem_line >> nodeID;
      localVolumeElementConnectivity.push_back(nodeID);
    }

    /*--- If a linear element is used, the node numbering for non-simplices
          must be adapted. The reason is that compatability with the original
          SU2 format is maintained for linear elements, but for the FEM solver
          the nodes of the elements are stored row-wise. ---*/
    if(nPolyGrid == 1) {
      switch( VTK_Type ) {

      case QUADRILATERAL:
        swap(localVolumeElementConnectivity[ind+2],
             localVolumeElementConnectivity[ind+3]);
        break;

      case HEXAHEDRON:
        swap(localVolumeElementConnectivity[ind+2],
             localVolumeElementConnectivity[ind+3]);
        swap(localVolumeElementConnectivity[ind+6],
             localVolumeElementConnectivity[ind+7]);
        break;

      case PYRAMID:
        swap(localVolumeElementConnectivity[ind+2],
             localVolumeElementConnectivity[ind+3]);
        break;
      }
    }
  }

  /*--- Close the mesh file again. ---*/
  mesh_file.close();
}

void CSU2ASCIIMeshReaderFEM::ReadSurfaceElementConnectivity() {

  /*--- Determine the vector to hold the faces of the local elements. ---*/
  vector<CFaceOfElement> localFaces;
  DetermineFacesVolumeElements(localFaces);

  /*--- We already read in the number of markers with the metadata.
        Allocate the memory for the marker names, number of local surface
        elements and the first index of the surface element connectivity. ---*/
  surfaceElementConnectivity.resize(numberOfMarkers);
  markerNames.resize(numberOfMarkers);
  numberOfLocalSurfaceElements.resize(numberOfMarkers, 0);

  /*--- Open the mesh file and jump to our zone. ---*/
  mesh_file.open(meshFilename, ios::in);
  FastForwardToMyZone();

  /*--- Find the section containing the markers. ---*/
  string text_line;
  while (getline (mesh_file, text_line)) {
    string::size_type position = text_line.find ("NMARK=",0);
    if (position != string::npos) break;
  }

  /*--- Loop over the number of boundary markers. ---*/
  for(unsigned long iMarker=0; iMarker<numberOfMarkers; ++iMarker) {

    /*--- Find the section containing the marker name. ---*/
    while (getline (mesh_file, text_line)) {
      string::size_type position = text_line.find ("MARKER_TAG=",0);
      if (position != string::npos) break;
    }

    /*--- Extract the marker name. Remove spaces returns and tabs
          and store the name in markerNames. ---*/
    text_line.erase(0,11);
    string::size_type position;
        
    for (unsigned short iChar = 0; iChar < 20; iChar++) {
      position = text_line.find( " ", 0 );
      if (position != string::npos) text_line.erase (position,1);
      position = text_line.find( "\r", 0 );
      if (position != string::npos) text_line.erase (position,1);
      position = text_line.find( "\n", 0 );
      if (position != string::npos) text_line.erase (position,1);
    }
    markerNames[iMarker] = text_line.c_str();

    /*--- Find the section containing the number of surface elements
          for this marker and determine this number. ---*/
    while (getline (mesh_file, text_line)) {
      string::size_type position = text_line.find ("MARKER_ELEMS=",0);
      if (position != string::npos) break;
    }

    text_line.erase (0,13);
    unsigned long nElem_Bound = atoi(text_line.c_str());

    /*--- Loop over the surface elements for this marker. ---*/
    for(unsigned long i=0; i<nElem_Bound; ++i) {

      /*--- Read the line for this element and couple it to an istringstream
            to enable the actual reading of the data. ---*/
      getline(mesh_file, text_line);
      istringstream bound_line(text_line);

      /*--- Determine the element type and polynomial degree. ---*/
      unsigned long  typeRead; bound_line >> typeRead;

      const unsigned short nPolyGrid = typeRead/100 + 1;
      const unsigned short VTK_Type  = typeRead%100;

      /*--- Make a distinction between the possible element surface types and
            determine the corner points in local numbering of the element. ---*/
      const unsigned short nDOFEdgeGrid = nPolyGrid + 1;

      unsigned short nDOFsGrid = 0;
      CFaceOfElement thisFace;
      thisFace.cornerPoints[0] = 0; thisFace.cornerPoints[1] = nPolyGrid;     

      switch( VTK_Type ) {
        case LINE:
          nDOFsGrid = nDOFEdgeGrid;
          thisFace.nCornerPoints = 2;
          break;

        case TRIANGLE:
          nDOFsGrid = nDOFEdgeGrid*(nDOFEdgeGrid+1)/2;
          thisFace.nCornerPoints = 3;
          thisFace.cornerPoints[2] = nDOFsGrid -1;
          break;

        case QUADRILATERAL:
          nDOFsGrid = nDOFEdgeGrid*nDOFEdgeGrid;
          thisFace.nCornerPoints = 4;
          thisFace.cornerPoints[2] = nPolyGrid*nDOFEdgeGrid;
          thisFace.cornerPoints[3] = nDOFsGrid -1;
          break;

        default:
          ostringstream message;
          message << "Unknown FEM boundary element value, " << typeRead
                  << ", in " << meshFilename;
          SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /*--- Read the connectivity information. ---*/
      vector<unsigned long> connFace(nDOFsGrid);
      for(unsigned short j=0; j<nDOFsGrid; ++j) bound_line >> connFace[j];

      /*--- If a linear quadrilateral is used, the node numbering must be adapted.
            The reason is that compatability with the original SU2 format is
            maintained for linear elements, but for the FEM solver the nodes
            of the elements are stored row-wise. ---*/
      if((nPolyGrid == 1) && (VTK_Type == QUADRILATERAL))
        swap(connFace[2], connFace[3]);

      /*--- Convert the local numbering of thisFace to global numbering
            and create a unique numbering of corner points. ---*/
      for(unsigned short j=0; j<thisFace.nCornerPoints; ++j)
        thisFace.cornerPoints[j] = connFace[thisFace.cornerPoints[j]];
      thisFace.CreateUniqueNumbering();

      /*--- Check if this boundary face must be stored on this rank. ---*/ 
      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);
      if(low != localFaces.end()) {
        if( !(thisFace < *low) ) {

          /*--- Update the counter for this boundary marker and store
                the meta data in surfaceElementConnectivity. ---*/
          ++numberOfLocalSurfaceElements[iMarker];

          surfaceElementConnectivity[iMarker].push_back(VTK_Type);
          surfaceElementConnectivity[iMarker].push_back(nPolyGrid);
          surfaceElementConnectivity[iMarker].push_back(nDOFsGrid);
          surfaceElementConnectivity[iMarker].push_back(i);            // Global surface elem ID.
          surfaceElementConnectivity[iMarker].push_back(low->elemID0); // Global volume elem ID.

          /*--- Copy the connectivity to surfaceElementConnectivity. ---*/
          surfaceElementConnectivity[iMarker].insert(surfaceElementConnectivity[iMarker].end(),
                                                     connFace.begin(), connFace.end());
        }
      }
    }
  }

  /*--- Close the mesh file again. ---*/
  mesh_file.close();
}
