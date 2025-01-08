/*!
 * \file geometry_structure_fem_part.cpp
 * \brief Main subroutines for distributin the grid for the Fluid FEM solver.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/geometry/CPhysicalGeometry.hpp"
#include "../../include/fem/fem_standard_element.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridFEM.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"

#ifdef HAVE_CGNS
#include "../../include/fem/fem_cgns_elements.hpp"
#endif
#include "../../include/adt/CADTElemClass.hpp"

#include "../../include/linear_algebra/blas_structure.hpp"
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
/*--- Epsilon definition ---*/

CFaceOfElement::CFaceOfElement() {
  nCornerPoints = 0;
  cornerPoints[0] = cornerPoints[1] = cornerPoints[2] = cornerPoints[3] = ULONG_MAX;
  elemID0 = elemID1 = ULONG_MAX;
  nPolyGrid0 = nPolyGrid1 = 0;
  nPolySol0 = nPolySol1 = 0;
  nDOFsElem0 = nDOFsElem1 = 0;
  elemType0 = elemType1 = 0;
  faceID0 = faceID1 = 0;
  periodicIndex = periodicIndexDonor = 0;
  faceIndicator = 0;

  JacFaceIsConsideredConstant = false;
  elem0IsOwner = false;
}

CFaceOfElement::CFaceOfElement(const unsigned short VTK_Type, const unsigned short nPoly, const unsigned long* Nodes) {
  /* Set the default values of the member variables. */
  nCornerPoints = 0;
  cornerPoints[0] = cornerPoints[1] = cornerPoints[2] = cornerPoints[3] = ULONG_MAX;
  elemID0 = elemID1 = ULONG_MAX;
  nPolyGrid0 = nPolyGrid1 = 0;
  nPolySol0 = nPolySol1 = 0;
  nDOFsElem0 = nDOFsElem1 = 0;
  elemType0 = elemType1 = 0;
  faceID0 = faceID1 = 0;
  periodicIndex = periodicIndexDonor = 0;
  faceIndicator = 0;

  JacFaceIsConsideredConstant = false;
  elem0IsOwner = false;

  /* Determine the face element type and set the corner points accordingly. */
  switch (VTK_Type) {
    case LINE: {
      nCornerPoints = 2;
      cornerPoints[0] = Nodes[0];
      cornerPoints[1] = Nodes[nPoly];
      break;
    }

    case TRIANGLE: {
      const unsigned short ind2 = (nPoly + 1) * (nPoly + 2) / 2 - 1;
      nCornerPoints = 3;
      cornerPoints[0] = Nodes[0];
      cornerPoints[1] = Nodes[nPoly];
      cornerPoints[2] = Nodes[ind2];
      break;
    }

    case QUADRILATERAL: {
      const unsigned short ind2 = nPoly * (nPoly + 1);
      const unsigned short ind3 = (nPoly + 1) * (nPoly + 1) - 1;
      nCornerPoints = 4;
      cornerPoints[0] = Nodes[0];
      cornerPoints[1] = Nodes[nPoly];
      cornerPoints[2] = Nodes[ind2];
      cornerPoints[3] = Nodes[ind3];
      break;
    }

    default: {
      ostringstream message;
      message << "Unknown VTK surface element type, " << VTK_Type;
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
    }
  }
}

bool CFaceOfElement::operator<(const CFaceOfElement& other) const {
  if (nCornerPoints != other.nCornerPoints) return nCornerPoints < other.nCornerPoints;

  for (unsigned short i = 0; i < nCornerPoints; ++i)
    if (cornerPoints[i] != other.cornerPoints[i]) return cornerPoints[i] < other.cornerPoints[i];

  return false;
}

bool CFaceOfElement::operator==(const CFaceOfElement& other) const {
  if (nCornerPoints != other.nCornerPoints) return false;
  for (unsigned short i = 0; i < nCornerPoints; ++i)
    if (cornerPoints[i] != other.cornerPoints[i]) return false;

  return true;
}

void CFaceOfElement::Copy(const CFaceOfElement& other) {
  nCornerPoints = other.nCornerPoints;
  for (unsigned short i = 0; i < nCornerPoints; ++i) cornerPoints[i] = other.cornerPoints[i];
  for (unsigned short i = nCornerPoints; i < 4; ++i) cornerPoints[i] = ULONG_MAX;

  elemID0 = other.elemID0;
  elemID1 = other.elemID1;

  nPolyGrid0 = other.nPolyGrid0;
  nPolyGrid1 = other.nPolyGrid1;
  nPolySol0 = other.nPolySol0;
  nPolySol1 = other.nPolySol1;

  nDOFsElem0 = other.nDOFsElem0;
  nDOFsElem1 = other.nDOFsElem1;
  elemType0 = other.elemType0;
  elemType1 = other.elemType1;
  faceID0 = other.faceID0;
  faceID1 = other.faceID1;

  periodicIndex = other.periodicIndex;
  periodicIndexDonor = other.periodicIndexDonor;
  faceIndicator = other.faceIndicator;

  JacFaceIsConsideredConstant = other.JacFaceIsConsideredConstant;
  elem0IsOwner = other.elem0IsOwner;
}

void CFaceOfElement::CreateUniqueNumberingWithOrientation() {
  /*--- Determine the element type and create the unique numbering accordingly. ---*/
  bool swapElements = false;
  switch (nCornerPoints) {
    case 2: {
      /* Element is a line. Check if the node numbering must be swapped. If so
         also the element information must be swapped, because element 0 is to
         the left of the face and element 1 to the right. */
      if (cornerPoints[1] < cornerPoints[0]) {
        swap(cornerPoints[0], cornerPoints[1]);
        swapElements = true;
      }
      break;
    }

    case 3: {
      /* Element is a triangle. The vertices are sorted in increasing order.
         If the sequence of the new numbering is opposite to the current
         numbering, the element information must be exchanged, because
         element 0 is to the left of the face and element 1 to the right. */
      unsigned long nn[] = {cornerPoints[0], cornerPoints[1], cornerPoints[2]};
      unsigned short ind = 0;
      if (nn[1] < nn[ind]) ind = 1;
      if (nn[2] < nn[ind]) ind = 2;

      unsigned short indm1 = ind == 0 ? 2 : ind - 1;  // Next lower index.
      unsigned short indp1 = ind == 2 ? 0 : ind + 1;  // Next upper index.

      if (nn[indp1] < nn[indm1]) {
        /* The orientation of the triangle remains the same.
           Store the new sorted node numbering. */
        cornerPoints[0] = nn[ind];
        cornerPoints[1] = nn[indp1];
        cornerPoints[2] = nn[indm1];
      } else {
        /* The orientation of the triangle changes. Store the new
           sorted node numbering and set swapElements to true. */
        cornerPoints[0] = nn[ind];
        cornerPoints[1] = nn[indm1];
        cornerPoints[2] = nn[indp1];
        swapElements = true;
      }

      break;
    }

    case 4: {
      /* Element is a quadrilateral. The vertices are sorted in increasing order
         under the condition neighboring vertices remain neighbors. If the
         sequence of the new numbering is opposite to the current
         numbering, the element information must be exchanged, because
         element 0 is to the left of the face and element 1 to the right. */
      unsigned long nn[] = {cornerPoints[0], cornerPoints[1], cornerPoints[2], cornerPoints[3]};
      unsigned short ind = 0;
      if (nn[1] < nn[ind]) ind = 1;
      if (nn[2] < nn[ind]) ind = 2;
      if (nn[3] < nn[ind]) ind = 3;

      unsigned short indm1 = ind == 0 ? 3 : ind - 1;        // Next lower index.
      unsigned short indp1 = ind == 3 ? 0 : ind + 1;        // Next upper index.
      unsigned short indp2 = ind >= 2 ? ind - 2 : ind + 2;  // Opposite index.

      if (nn[indp1] < nn[indm1]) {
        /* The orientation of the quadrilateral remains the same.
           Store the new sorted node numbering. */
        cornerPoints[0] = nn[ind];
        cornerPoints[1] = nn[indp1];
        cornerPoints[2] = nn[indp2];
        cornerPoints[3] = nn[indm1];
      } else {
        /* The orientation of the quadrilateral changes. Store the new
           sorted node numbering and set swapElements to true. */
        cornerPoints[0] = nn[ind];
        cornerPoints[1] = nn[indm1];
        cornerPoints[2] = nn[indp2];
        cornerPoints[3] = nn[indp1];
        swapElements = true;
      }

      break;
    }

    default: {
      ostringstream message;
      message << "Unknown surface element type with " << nCornerPoints << " corners." << endl;
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
    }
  }

  /* Swap the element information, if needed. */
  if (swapElements) {
    swap(elemID0, elemID1);
    swap(nPolyGrid0, nPolyGrid1);
    swap(nPolySol0, nPolySol1);
    swap(nDOFsElem0, nDOFsElem1);
    swap(elemType0, elemType1);
    swap(faceID0, faceID1);
  }
}

void CBoundaryFace::Copy(const CBoundaryFace& other) {
  VTK_Type = other.VTK_Type;
  nPolyGrid = other.nPolyGrid;
  nDOFsGrid = other.nDOFsGrid;
  globalBoundElemID = other.globalBoundElemID;
  domainElementID = other.domainElementID;
  Nodes = other.Nodes;
}

CMatchingFace::CMatchingFace() {
  nCornerPoints = 0;
  nDim = 0;
  nPoly = 0;
  nDOFsElem = 0;
  elemType = 0;
  elemID = 0;

  cornerCoor[0][0] = cornerCoor[0][1] = cornerCoor[0][2] = 0.0;
  cornerCoor[1][0] = cornerCoor[1][1] = cornerCoor[1][2] = 0.0;
  cornerCoor[2][0] = cornerCoor[2][1] = cornerCoor[2][2] = 0.0;
  cornerCoor[3][0] = cornerCoor[3][1] = cornerCoor[3][2] = 0.0;

  tolForMatching = 0.0;
}

bool CMatchingFace::operator<(const CMatchingFace& other) const {
  /* First compare the number of corner points. ---*/
  if (nCornerPoints != other.nCornerPoints) return nCornerPoints < other.nCornerPoints;

  /*--- Determine the tolerance for comparing both objects. ---*/
  const su2double tol = min(tolForMatching, other.tolForMatching);

  /*--- Loop over the number of corner points and dimensions and compare the
        coordinates. If considered different, return true if the current face
        is considered smaller and false otherwise.            ---*/
  for (unsigned short k = 0; k < nCornerPoints; ++k) {
    for (unsigned short l = 0; l < nDim; ++l) {
      if (fabs(cornerCoor[k][l] - other.cornerCoor[k][l]) > tol) return cornerCoor[k][l] < other.cornerCoor[k][l];
    }
  }

  /*--- Both objects are considered the same. Return false. ---*/
  return false;
}

void CMatchingFace::SortFaceCoordinates() {
  /*--- Determine the tolerance for a matching point for this face. This is
        accomplished by computing the minimum distance between the points of
        the face, multiplied by a relative tolerance.        ---*/
  for (unsigned short k = 0; k < nCornerPoints; ++k) {
    for (unsigned short j = (k + 1); j < nCornerPoints; ++j) {
      su2double dist = 0.0;
      for (unsigned short l = 0; l < nDim; ++l) {
        su2double ds = cornerCoor[k][l] - cornerCoor[j][l];
        dist += ds * ds;
      }
      dist = sqrt(dist);

      if (k == 0 && j == 1)
        tolForMatching = dist;
      else
        tolForMatching = min(tolForMatching, dist);
    }
  }

  tolForMatching *= 1.e-6;

  /*--- Sort the points in increasing order based on the coordinates.
        An insertion sort algorithm is used, which is quite efficient
        for at most four corner points.                        ---*/
  for (unsigned short k = 1; k < nCornerPoints; ++k) {
    for (unsigned short j = k; j > 0; --j) {
      /* Check if cornerCoor[j] is considered less than cornerCoor[j-1]. */
      bool lessThan = false;
      for (unsigned short l = 0; l < nDim; ++l) {
        if (fabs(cornerCoor[j][l] - cornerCoor[j - 1][l]) > tolForMatching) {
          lessThan = cornerCoor[j][l] < cornerCoor[j - 1][l];
          break;
        }
      }

      /* If cornerCoor[j] is less than cornerCoor[j-1] they must be swapped.
         Otherwise an exit can be made from the j-loop.         */
      if (lessThan) {
        for (unsigned short l = 0; l < nDim; ++l) swap(cornerCoor[j][l], cornerCoor[j - 1][l]);
      } else
        break;
    }
  }
}

void CMatchingFace::Copy(const CMatchingFace& other) {
  nCornerPoints = other.nCornerPoints;
  nDim = other.nDim;
  nPoly = other.nPoly;
  nDOFsElem = other.nDOFsElem;
  elemType = other.elemType;
  elemID = other.elemID;

  for (unsigned short k = 0; k < nCornerPoints; ++k) {
    for (unsigned l = 0; l < nDim; ++l) {
      cornerCoor[k][l] = other.cornerCoor[k][l];
    }
  }

  tolForMatching = other.tolForMatching;
}

void CPhysicalGeometry::Read_SU2_Format_Parallel_FEM(CConfig* config, const string& val_mesh_filename,
                                                     unsigned short val_iZone, unsigned short val_nZone) {
  string text_line, Marker_Tag;
  ifstream mesh_file;
  string::size_type position;
  unsigned long nDOFsGrid_Local = 0, loc_element_count = 0;
  bool domain_flag = false;
  bool time_spectral = config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE;
  unsigned short nMarker_Max = config->GetnMarker_Max();
  nZone = val_nZone;

  /*--- Initialize counters for local/global points & elements ---*/
  Global_nPoint = 0;
  Global_nPointDomain = 0;
  Global_nElem = 0;
  nelem_edge = 0;
  Global_nelem_edge = 0;
  nelem_triangle = 0;
  Global_nelem_triangle = 0;
  nelem_quad = 0;
  Global_nelem_quad = 0;
  nelem_tetra = 0;
  Global_nelem_tetra = 0;
  nelem_hexa = 0;
  Global_nelem_hexa = 0;
  nelem_prism = 0;
  Global_nelem_prism = 0;
  nelem_pyramid = 0;
  Global_nelem_pyramid = 0;

  /*--- Allocate memory for the linear partition of the elements of the mesh.
        These arrays are the size of the number of ranks. ---*/

  beg_node = new unsigned long[size];
  end_node = new unsigned long[size];
  nPointLinear = new unsigned long[size];

  /*--- Open grid file ---*/

  mesh_file.open(val_mesh_filename.c_str(), ios::in);

  /*--- Check the grid ---*/

  if (mesh_file.fail())
    SU2_MPI::Error(string("There is no mesh file (CPhysicalGeometry)!! ") + val_mesh_filename, CURRENT_FUNCTION);

  /*--- If more than one, find the zone in the mesh file ---*/

  if (val_nZone > 1 || time_spectral) {
    if (time_spectral) {
      if (rank == MASTER_NODE) cout << "Reading time spectral instance " << val_iZone + 1 << ":" << endl;
    } else {
      while (getline(mesh_file, text_line)) {
        /*--- Search for the current domain ---*/
        position = text_line.find("IZONE=", 0);
        if (position != string::npos) {
          text_line.erase(0, 6);
          unsigned short jDomain = atoi(text_line.c_str());
          if (jDomain == val_iZone + 1) {
            if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone + 1 << " points:" << endl;
            break;
          }
        }
      }
    }
  }

  /*--- Read grid file with format SU2 ---*/

  while (getline(mesh_file, text_line)) {
    /*--- Read the dimension of the problem ---*/

    position = text_line.find("NDIME=", 0);
    if (position != string::npos) {
      if (!domain_flag) {
        text_line.erase(0, 6);
        nDim = atoi(text_line.c_str());
        if (rank == MASTER_NODE) {
          if (nDim == 2) cout << "Two dimensional problem." << endl;
          if (nDim == 3) cout << "Three dimensional problem." << endl;
        }
        domain_flag = true;
      } else {
        break;
      }
    }

    /*--- Read the information about inner elements ---*/

    position = text_line.find("NELEM=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      stringstream stream_line(text_line);
      stream_line >> Global_nElem;

      if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
        cout << Global_nElem << " interior elements before parallel partitioning." << endl;

      /*--- Check if the number of cores used is larger than the number of
            elements. Terminate if this is the case, because it does not make
            sense to do this. ---*/
      unsigned long nCores = size;  // Correct for the number of cores per rank.
      if (nCores > Global_nElem) {
        ostringstream message;
        message << "The number of cores, " << nCores;
        message << ", is larger than the number of elements, " << Global_nElem << "." << endl;
        message << "This is not an efficient use of the resources and therefore "
                << "SU2 will terminate.";

        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /*--- Compute the number of elements that will be on each processor.
            This is a linear partitioning with the addition of a simple load
            balancing for any remainder elements. ---*/
      unsigned long total_elem_accounted = 0;
      for (unsigned long i = 0; i < (unsigned long)size; i++) {
        nPointLinear[i] = Global_nElem / size;
        total_elem_accounted = total_elem_accounted + nPointLinear[i];
      }

      /*--- Get the number of remainder elements after the even division ---*/
      unsigned long rem_elem = Global_nElem - total_elem_accounted;
      for (unsigned long i = 0; i < rem_elem; i++) {
        ++nPointLinear[i];
      }

      /*--- Store the local number of elements and the beginning/end index ---*/
      nElem = nPointLinear[rank];
      beg_node[0] = 0;
      end_node[0] = beg_node[0] + nPointLinear[0];
      for (unsigned long i = 1; i < (unsigned long)size; i++) {
        beg_node[i] = end_node[i - 1];
        end_node[i] = beg_node[i] + nPointLinear[i];
      }

      /*--- Allocate space for elements ---*/
      elem = new CPrimalGrid*[nElem]();

      /*--- Loop over all the elements and store the elements to be stored on
            this rank. Furthermore, determine the total amount of DOFs for
            the solution (which may differ from the number of DOFS for the
            grid).                                                          ---*/
      unsigned long nDOFs_tot = 0;
      for (unsigned long i = 0; i < Global_nElem; i++) {
        getline(mesh_file, text_line);
        istringstream elem_line(text_line);

        /*--- Read the value that defines the element type, the polynomial
              degree of the geometry and the polynomial degree of the
              solution. Extract this info as well.                     ---*/
        unsigned long typeRead;
        elem_line >> typeRead;
        unsigned long typeReadErrorMessage = typeRead;
        unsigned short nPolySol, nPolyGrid;
        if (typeRead > 10000) {
          nPolySol = typeRead / 10000 - 1;
          typeRead = typeRead % 10000;
          nPolyGrid = typeRead / 100 + 1;
        } else {
          nPolyGrid = typeRead / 100 + 1;
          nPolySol = nPolyGrid;
        }

        unsigned short VTK_Type = typeRead % 100;

        unsigned short nDOFsGrid = CFEMStandardElementBase::GetNDOFsStatic(VTK_Type, nPolyGrid, typeReadErrorMessage);
        unsigned short nDOFsSol = CFEMStandardElementBase::GetNDOFsStatic(VTK_Type, nPolySol, typeReadErrorMessage);

        /*--- Allocate the memory for a new primary grid FEM element if
              this element must be stored on this rank.                 ---*/
        if ((i >= beg_node[rank]) && (i < end_node[rank])) {
          elem[loc_element_count] =
              new CPrimalGridFEM(i, VTK_Type, nPolyGrid, nPolySol, nDOFsGrid, nDOFsSol, nDOFs_tot, elem_line);
          nDOFsGrid_Local += nDOFsGrid;
          loc_element_count++;
        }

        /*--- Update the total value of the number of DOFs. ---*/
        nDOFs_tot += nDOFsSol;
      }

      /*--- Break from the loop to find the element information. ---*/
      break;
    }
  }

  mesh_file.close();

  /*--- Create a vector, which contains the global node IDs of the local elements. ---*/
  vector<unsigned long> nodeIDsElemLoc;
  nodeIDsElemLoc.reserve(nDOFsGrid_Local);
  for (unsigned long i = 0; i < loc_element_count; i++) {
    unsigned short nDOFsElem = elem[i]->GetnNodes();
    for (unsigned short j = 0; j < nDOFsElem; ++j) nodeIDsElemLoc.push_back(elem[i]->GetNode(j));
  }

  sort(nodeIDsElemLoc.begin(), nodeIDsElemLoc.end());
  vector<unsigned long>::iterator lastNode;
  lastNode = unique(nodeIDsElemLoc.begin(), nodeIDsElemLoc.end());
  nodeIDsElemLoc.erase(lastNode, nodeIDsElemLoc.end());

  /*--- Allocate the memory for the coordinates to be stored on this rank. ---*/
  nPoint = nodeIDsElemLoc.size();
  nodes = new CPoint(nPoint, nDim);

  /*--- Open the grid file again and go to the position where
        the correct zone is stored.                           ---*/
  mesh_file.open(val_mesh_filename.c_str(), ios::in);

  if (val_nZone > 1 && !time_spectral) {
    while (getline(mesh_file, text_line)) {
      position = text_line.find("IZONE=", 0);
      if (position != string::npos) {
        text_line.erase(0, 6);
        unsigned short jDomain = atoi(text_line.c_str());
        if (jDomain == val_iZone + 1) break;
      }
    }
  }

  /*--- While loop to read the point information. ---*/
  while (getline(mesh_file, text_line)) {
    position = text_line.find("NPOIN=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      stringstream stream_line(text_line);
      stream_line >> Global_nPoint;

      /*--- Loop over the global number of points and store the
            ones that are needed on this processor.             ---*/
      unsigned long ii = 0;
      for (unsigned long i = 0; i < Global_nPoint; ++i) {
        getline(mesh_file, text_line);

        if (binary_search(nodeIDsElemLoc.begin(), nodeIDsElemLoc.end(), i)) {
          istringstream point_line(text_line);
          su2double Coord[3] = {0.0};
          point_line >> Coord[0];
          point_line >> Coord[1];
          if (nDim == 3) point_line >> Coord[2];
          nodes->SetCoord(ii, Coord);
          nodes->SetGlobalIndex(ii, i);
          ++ii;
        }
      }

      break;
    }
  }

  mesh_file.close();

  /*--- Determine the faces of the local elements. --- */
  vector<CFaceOfElement> localFaces;
  for (unsigned long k = 0; k < loc_element_count; k++) {
    /*--- Get the global IDs of the corner points of all the faces of this elements. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long faceConn[6][4];

    elem[k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Loop over the faces and add them to localFaces. ---*/
    for (unsigned short i = 0; i < nFaces; ++i) {
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for (unsigned short j = 0; j < nPointsPerFace[i]; ++j) thisFace.cornerPoints[j] = faceConn[i][j];
      thisFace.elemID0 = k + beg_node[rank];

      thisFace.CreateUniqueNumbering();
      localFaces.push_back(thisFace);
    }
  }

  /*--- Sort localFaces in increasing order and remove the double entities,
        such that the binary search later on is a bit more efficient.       ---*/
  sort(localFaces.begin(), localFaces.end());
  vector<CFaceOfElement>::iterator lastFace;
  lastFace = unique(localFaces.begin(), localFaces.end());
  localFaces.erase(lastFace, localFaces.end());

  /*--- Open the grid file again and go to the position where
        the correct zone is stored.                           ---*/
  mesh_file.open(val_mesh_filename.c_str(), ios::in);

  if (val_nZone > 1 && !time_spectral) {
    while (getline(mesh_file, text_line)) {
      position = text_line.find("IZONE=", 0);
      if (position != string::npos) {
        text_line.erase(0, 6);
        unsigned short jDomain = atoi(text_line.c_str());
        if (jDomain == val_iZone + 1) break;
      }
    }
  }

  /*--- While loop to read the boundary information. ---*/
  while (getline(mesh_file, text_line)) {
    /*--- Read number of markers ---*/
    position = text_line.find("NMARK=", 0);
    if (position != string::npos) {
      text_line.erase(0, 6);
      istringstream stream_line(text_line);
      stream_line >> nMarker;

      if (rank == MASTER_NODE) cout << nMarker << " surface markers." << endl;
      config->SetnMarker_All(nMarker);
      bound = new CPrimalGrid**[nMarker];
      nElem_Bound = new unsigned long[nMarker];
      Tag_to_Marker = new string[nMarker_Max];

      for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++) {
        getline(mesh_file, text_line);
        text_line.erase(0, 11);
        text_line.erase(remove(text_line.begin(), text_line.end(), ' '), text_line.end());
        text_line.erase(remove(text_line.begin(), text_line.end(), '\r'), text_line.end());
        text_line.erase(remove(text_line.begin(), text_line.end(), '\n'), text_line.end());

        Marker_Tag = text_line;

        /*--- Read the number of elements for this marker. ---*/
        getline(mesh_file, text_line);
        text_line.erase(0, 13);
        istringstream nmark_line(text_line);
        unsigned long nElem_Bound_Global;
        nmark_line >> nElem_Bound_Global;
        if (rank == MASTER_NODE)
          cout << nElem_Bound_Global << " boundary elements in index " << iMarker << " (Marker = " << Marker_Tag << ")."
               << endl;

        /*--- Define a vector of FEM boundary elements to store the local boundary faces. ---*/
        vector<CBoundaryFace> boundElems;

        /*--- Loop over the global boundary faces. ---*/
        for (unsigned long i = 0; i < nElem_Bound_Global; ++i) {
          getline(mesh_file, text_line);
          istringstream bound_line(text_line);

          /*--- Determine the element type, its number of DOFs and read
                its node IDs.                                           ---*/
          unsigned long typeRead;
          bound_line >> typeRead;
          unsigned short nPolyGrid = typeRead / 100 + 1;
          unsigned short VTK_Type = typeRead % 100;
          unsigned short nDOFEdgeGrid = nPolyGrid + 1;

          unsigned short nDOFsGrid = 0;
          CFaceOfElement thisFace;
          thisFace.cornerPoints[0] = 0;
          thisFace.cornerPoints[1] = nPolyGrid;
          switch (VTK_Type) {
            case LINE:
              nDOFsGrid = nDOFEdgeGrid;
              thisFace.nCornerPoints = 2;
              break;

            case TRIANGLE:
              nDOFsGrid = nDOFEdgeGrid * (nDOFEdgeGrid + 1) / 2;
              thisFace.nCornerPoints = 3;
              thisFace.cornerPoints[2] = nDOFsGrid - 1;
              break;

            case QUADRILATERAL:
              nDOFsGrid = nDOFEdgeGrid * nDOFEdgeGrid;
              thisFace.nCornerPoints = 4;
              thisFace.cornerPoints[2] = nPolyGrid * nDOFEdgeGrid;
              thisFace.cornerPoints[3] = nDOFsGrid - 1;
              break;

            default:
              ostringstream message;
              message << "Unknown FEM boundary element value, " << typeRead << ", in " << val_mesh_filename;
              SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
          }

          vector<unsigned long> nodeIDs(nDOFsGrid);
          for (unsigned short j = 0; j < nDOFsGrid; ++j) bound_line >> nodeIDs[j];

          /*--- Convert the local numbering of thisFace to global numbering
                and create a unique numbering of these nodes.               ---*/
          for (unsigned short j = 0; j < thisFace.nCornerPoints; ++j)
            thisFace.cornerPoints[j] = nodeIDs[thisFace.cornerPoints[j]];
          thisFace.CreateUniqueNumbering();

          /*--- Check if this boundary face must be stored on this rank.
                If so, create an object of CBoundaryFace and add it
                to boundElems.                                          ---*/
          vector<CFaceOfElement>::iterator low;
          low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);
          if (low != localFaces.end()) {
            if (!(thisFace < *low)) {
              CBoundaryFace thisBoundFace;
              thisBoundFace.VTK_Type = VTK_Type;
              thisBoundFace.nPolyGrid = nPolyGrid;
              thisBoundFace.nDOFsGrid = nDOFsGrid;
              thisBoundFace.globalBoundElemID = i;
              thisBoundFace.domainElementID = low->elemID0;
              thisBoundFace.Nodes = nodeIDs;

              boundElems.push_back(thisBoundFace);
            }
          }
        }

        /*--- Allocate space for the boundary elements and store the ones
              whose parent element is stored on this rank. ---*/
        nElem_Bound[iMarker] = boundElems.size();
        bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];

        for (unsigned long i = 0; i < nElem_Bound[iMarker]; ++i)
          bound[iMarker][i] = new CPrimalGridBoundFEM(boundElems[i].globalBoundElemID, boundElems[i].domainElementID,
                                                      boundElems[i].VTK_Type, boundElems[i].nPolyGrid,
                                                      boundElems[i].nDOFsGrid, boundElems[i].Nodes);

        /*--- Update config information storing the boundary information in the right place ---*/
        Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
        config->SetMarker_All_TagBound(iMarker, Marker_Tag);
        config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
        config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
        config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
        config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
        config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
        config->SetMarker_All_Analyze(iMarker, config->GetMarker_CfgFile_Analyze(Marker_Tag));
        config->SetMarker_All_ZoneInterface(iMarker, config->GetMarker_CfgFile_ZoneInterface(Marker_Tag));
        config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
        config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
        config->SetMarker_All_SobolevBC(iMarker, config->GetMarker_CfgFile_SobolevBC(Marker_Tag));
        config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
        config->SetMarker_All_SendRecv(iMarker, NONE);
      }

      break;
    }
  }

  mesh_file.close();
}

void CPhysicalGeometry::Read_CGNS_Format_Parallel_FEM(CConfig* config, const string& val_mesh_filename,
                                                      unsigned short val_iZone, unsigned short val_nZone) {
#ifdef HAVE_CGNS

  /*--- For proper support of the high order elements, at least version 3.3
        of CGNS must be used. Check this. ---*/
#if CGNS_VERSION >= 3300

  /*--- Check whether the supplied file is truly a CGNS file. ---*/
  int file_type;
  if (cg_is_cgns(val_mesh_filename.c_str(), &file_type) != CG_OK)
    SU2_MPI::Error(val_mesh_filename + string(" is not a CGNS file that can be read."), CURRENT_FUNCTION);

  /*--- Initialize counters for local/global points & elements ---*/
  Global_nPoint = 0;
  Global_nPointDomain = 0;
  Global_nElem = 0;
  nelem_edge = 0;
  Global_nelem_edge = 0;
  nelem_triangle = 0;
  Global_nelem_triangle = 0;
  nelem_quad = 0;
  Global_nelem_quad = 0;
  nelem_tetra = 0;
  Global_nelem_tetra = 0;
  nelem_hexa = 0;
  Global_nelem_hexa = 0;
  nelem_prism = 0;
  Global_nelem_prism = 0;
  nelem_pyramid = 0;
  Global_nelem_pyramid = 0;

  /*--------------------------------------------------------------------------*/
  /*---       Checking of the file, determine the dimensions, etc.         ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Allocate memory for the linear partition of the elements of the mesh.
        These arrays are the size of the number of ranks. ---*/
  beg_node = new unsigned long[size];
  end_node = new unsigned long[size];
  nPointLinear = new unsigned long[size];

  /* Open the CGNS file for reading and check if it went OK. */
  int fn;
  if (cg_open(val_mesh_filename.c_str(), CG_MODE_READ, &fn) != CG_OK) cg_error_exit();
  if (rank == MASTER_NODE) {
    cout << "Reading the CGNS file: " << val_mesh_filename << "." << endl;
  }

  /* Get the number of databases. This is the highest node in the CGNS heirarchy.
     The current implementation assumes that there is only one database. */
  int nbases;
  if (cg_nbases(fn, &nbases) != CG_OK) cg_error_exit();
  if (nbases > 1) {
    ostringstream message;
    message << "CGNS file contains " << nbases << " databases." << endl;
    message << "CGNS reader can handle only 1 at the moment." << endl;

    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /* Read the information of the base, especially the number of dimensions. */
  char cgnsname[CGNS_STRING_SIZE];
  int cellDim, physDim;
  const int iBase = 1;
  if (cg_base_read(fn, iBase, cgnsname, &cellDim, &physDim)) cg_error_exit();
  nDim = physDim;

  if (cellDim != physDim) {
    ostringstream message;
    message << "The element dimension, " << cellDim << ", differs from the physical dimension, " << physDim << "."
            << endl;
    message << "These should be the same for the DG-FEM solver." << endl;

    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /* Write the info about the number of dimensions. */
  if (rank == MASTER_NODE) {
    if (nDim == 2) cout << "Two dimensional problem." << endl;
    if (nDim == 3) cout << "Three dimensional problem." << endl;
  }

  /* Set the zone ID for CGNS, CGNS indices starts at 1, and check
     if the requested zone ID is a valid one. */
  const int iZone = val_iZone + 1;

  int nzones;
  if (cg_nzones(fn, iBase, &nzones) != CG_OK) cg_error_exit();
  if (iZone < 1 || iZone > nzones) {
    ostringstream message;
    message << "Zone " << iZone << " requested for reading, but there are only " << nzones
            << " zones present in the CGNS file." << endl;

    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /* Determine the zone type for the requested zone and check if it is
     unstructured. */
  ZoneType_t zoneType;
  if (cg_zone_type(fn, iBase, iZone, &zoneType) != CG_OK) cg_error_exit();
  if (zoneType != Unstructured)
    SU2_MPI::Error("Structured CGNS zone found while unstructured expected.", CURRENT_FUNCTION);

  /* Determine the number of sections for the connectivities in this zone. */
  int nsections;
  if (cg_nsections(fn, iBase, iZone, &nsections) != CG_OK) cg_error_exit();

  /*--------------------------------------------------------------------------*/
  /*---            Reading and distributing the volume elements.           ---*/
  /*--------------------------------------------------------------------------*/

  /* Loop over the sections to store the meta data in CGNSElemTypes
     and to determine the number of volume elements.
     Note that the indices start at 1 in CGNS. */
  vector<CCGNSElementType> CGNSElemTypes(nsections);

  for (int iConn = 1; iConn <= nsections; ++iConn) {
    CGNSElemTypes[iConn - 1].DetermineMetaData(nDim, fn, iBase, iZone, iConn);
    if (CGNSElemTypes[iConn - 1].volumeConn) Global_nElem += CGNSElemTypes[iConn - 1].nElem;
  }

  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << Global_nElem << " interior elements before parallel partitioning." << endl;

  /*--- Check if the number of cores used is larger than the number of
        elements. Terminate if this is the case, because it does not make
        sense to do this. ---*/
  unsigned long nCores = size;  // Correct for the number of cores per rank.
  if (nCores > Global_nElem) {
    ostringstream message;
    message << "The number of cores, " << nCores;
    message << ", is larger than the number of elements, " << Global_nElem << "." << endl;
    message << "This is not an efficient use of the resources and therefore "
            << "SU2 will terminate.";

    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /*--- Compute the number of elements that will be stored on each rank.
        This is a linear partitioning with the addition of a simple load
        balancing for any remainder elements. ---*/
  unsigned long total_elem_accounted = 0;
  for (unsigned long i = 0; i < (unsigned long)size; i++) {
    nPointLinear[i] = Global_nElem / size;
    total_elem_accounted = total_elem_accounted + nPointLinear[i];
  }

  /*--- Get the number of remainder elements after the even division ---*/
  const unsigned long rem_elem = Global_nElem - total_elem_accounted;
  for (unsigned long i = 0; i < rem_elem; i++) {
    ++nPointLinear[i];
  }

  /*--- Store the local number of elements and the beginning/end index ---*/
  nElem = nPointLinear[rank];
  beg_node[0] = 0;
  end_node[0] = beg_node[0] + nPointLinear[0];
  for (int i = 1; i < size; i++) {
    beg_node[i] = end_node[i - 1];
    end_node[i] = beg_node[i] + nPointLinear[i];
  }

  /*--- Allocate space for elements ---*/
  elem = new CPrimalGrid*[nElem]();

  /*--- Loop over over the connectivities and read the elements to be stored on
        this rank. Furthermore, determine the local amount of DOFs for the
        solution (which may differ from the number of DOFS for the grid). ---*/
  unsigned long nDOFsLoc = 0, elemCount = 0, locElemCount = 0;

  for (int iConn = 0; iConn < nsections; ++iConn) {
    if (CGNSElemTypes[iConn].volumeConn) {
      /* Determine the global volume element range for this connectivity. */
      const unsigned long elemCountOld = elemCount;
      elemCount += CGNSElemTypes[iConn].nElem;

      /* Check for overlap with the element range this rank is responsible for. */
      const unsigned long indBegOverlap = max(elemCountOld, beg_node[rank]);
      const unsigned long indEndOverlap = min(elemCount, end_node[rank]);

      if (indEndOverlap > indBegOverlap) {
        /* This rank must read element data from this connectivity section.
           Determine the offset relative to the start of this section and the
           number of elements to be read by this rank. */
        const unsigned long offsetRank = indBegOverlap - elemCountOld;
        const unsigned long nElemRank = indEndOverlap - indBegOverlap;

        /* Read the connectivity range determined above. */
        CGNSElemTypes[iConn].ReadConnectivityRange(fn, iBase, iZone, offsetRank, nElemRank, beg_node[rank], elem,
                                                   locElemCount, nDOFsLoc);
      }
    }
  }

#ifdef HAVE_MPI
  /* The global offset of the DOFs must be corrected when running in
     parallel. Therefore gather the number of DOFs of all the ranks. */
  vector<unsigned long> nDOFsPerRank(size);
  SU2_MPI::Allgather(&nDOFsLoc, 1, MPI_UNSIGNED_LONG, nDOFsPerRank.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  /* Determine the offset for the DOFs on this rank. */
  unsigned long offsetRank = 0;
  for (int i = 0; i < rank; ++i) offsetRank += nDOFsPerRank[i];

  /* Loop over the local elements to correct the global offset of the DOFs. */
  for (unsigned long i = 0; i < nElem; ++i) elem[i]->AddOffsetGlobalDOFs(offsetRank);
#endif

  /*--------------------------------------------------------------------------*/
  /*---            Reading and distributing the coordinates.               ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the global number of vertices in the requested zone.
     The other size information is not used. */
  cgsize_t sizes[3];
  if (cg_zone_read(fn, iBase, iZone, cgnsname, sizes) != CG_OK) cg_error_exit();
  Global_nPoint = sizes[0];

  /*--- Determine the number of points per rank in cumulative storage format.
        This is done to avoid that every rank reads all the coordinates.
        The required coordinates for each rank are later obtained via
        communication.  ---*/
  unsigned long totalPointsAccounted = 0;
  vector<cgsize_t> nPointsPerRank(size + 1);
  for (int i = 1; i <= size; ++i) {
    nPointsPerRank[i] = Global_nPoint / size;
    totalPointsAccounted += nPointsPerRank[i];
  }

  const unsigned long nPointsRem = Global_nPoint - totalPointsAccounted;
  for (unsigned long i = 1; i <= nPointsRem; ++i) ++nPointsPerRank[i];

  nPointsPerRank[0] = 0;
  for (int i = 0; i < size; ++i) nPointsPerRank[i + 1] += nPointsPerRank[i];

  /* Determine the number of points that must be read by this rank and
     allocate the memory for the coordinate buffers. */
  const cgsize_t nPointsRead = nPointsPerRank[rank + 1] - nPointsPerRank[rank];
  vector<vector<su2double> > coorBuf(nDim, vector<su2double>(nPointsRead));

  /* Loop over the number of dimensions to read the coordinates. Note that
     the loop starts at 1 and ends at nDim because CGNS requires this. */
  for (unsigned short iDim = 1; iDim <= nDim; ++iDim) {
    /* Determine the data type and name of the coordinate. Copy the name
       of the coordinate in a string for easier comparison. */
    DataType_t datatype;
    if (cg_coord_info(fn, iBase, iZone, iDim, &datatype, cgnsname) != CG_OK) cg_error_exit();
    string coorname = cgnsname;

    /* Check the name of the coordinate and determine the index in coorBuf
       where to store this coordinate. Normally this should be iDim-1. */
    unsigned short indC = 0;
    if (coorname == "CoordinateX")
      indC = 0;
    else if (coorname == "CoordinateY")
      indC = 1;
    else if (coorname == "CoordinateZ")
      indC = 2;
    else
      SU2_MPI::Error(string("Unknown coordinate name, ") + coorname + string(", encountered in the CGNS file."),
                     CURRENT_FUNCTION);

    /* Easier storage of the range in the CGNS file. */
    cgsize_t range_min = nPointsPerRank[rank] + 1;
    cgsize_t range_max = nPointsPerRank[rank + 1];

    /*--- Read the coordinate with the required precision and copy
          this data to the correct index in coorBuf. ---*/
    switch (datatype) {
      case RealSingle: {
        /* Single precision used. */
        vector<float> buf(nPointsRead);
        if (cg_coord_read(fn, iBase, iZone, cgnsname, datatype, &range_min, &range_max, buf.data()) != CG_OK)
          cg_error_exit();

        for (cgsize_t i = 0; i < nPointsRead; ++i) coorBuf[indC][i] = buf[i];
        break;
      }

      case RealDouble: {
        /* Double precision used. */
        vector<double> buf(nPointsRead);
        if (cg_coord_read(fn, iBase, iZone, cgnsname, datatype, &range_min, &range_max, buf.data()) != CG_OK)
          cg_error_exit();

        for (cgsize_t i = 0; i < nPointsRead; ++i) coorBuf[indC][i] = buf[i];
        break;
      }

      default: {
        ostringstream message;
        message << "Datatype for coordinates must be RealSingle or RealDouble, "
                << "not " << datatype << endl;
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }
    }
  }

  /* Make a distinction between sequential and parallel mode. */
#ifdef HAVE_MPI
  /*--- Parallel mode. Create a vector, which contains the global
        node IDs of the local elements. ---*/
  vector<unsigned long> nodeIDsElemLoc;
  nodeIDsElemLoc.reserve(nDOFsLoc);
  for (unsigned long i = 0; i < locElemCount; ++i) {
    unsigned short nDOFsElem = elem[i]->GetnNodes();
    for (unsigned short j = 0; j < nDOFsElem; ++j) nodeIDsElemLoc.push_back(elem[i]->GetNode(j));
  }

  sort(nodeIDsElemLoc.begin(), nodeIDsElemLoc.end());
  vector<unsigned long>::iterator lastNode;
  lastNode = unique(nodeIDsElemLoc.begin(), nodeIDsElemLoc.end());
  nodeIDsElemLoc.erase(lastNode, nodeIDsElemLoc.end());

  /*--- Allocate the memory for the coordinates to be stored on this rank. ---*/
  nPoint = nodeIDsElemLoc.size();
  nodes = new CPoint(nPoint, nDim);

  /*--- Store the global ID's of the nodes in such a way that they can
        be sent to the rank that actually stores the coordinates.. ---*/
  vector<vector<unsigned long> > nodeBuf(size, vector<unsigned long>(0));
  for (unsigned long i = 0; i < nodeIDsElemLoc.size(); ++i) {
    const cgsize_t nodeID = nodeIDsElemLoc[i];
    vector<cgsize_t>::iterator low;
    low = lower_bound(nPointsPerRank.begin(), nPointsPerRank.end(), nodeID);
    cgsize_t rankNode = low - nPointsPerRank.begin();
    if (*low > nodeID) --rankNode;

    nodeBuf[rankNode].push_back(nodeIDsElemLoc[i]);
  }

  /*--- Determine the total number of ranks to which this rank will send
        a message and also determine the number of ranks from which this
        rank will receive a message. Furthermore, determine the starting
        indices where data from the different ranks should be stored in
        node. ---*/
  int nRankSend = 0;
  vector<int> sendToRank(size, 0);
  vector<unsigned long> startingIndRanksInNode(size + 1);
  startingIndRanksInNode[0] = 0;

  for (int i = 0; i < size; ++i) {
    startingIndRanksInNode[i + 1] = startingIndRanksInNode[i] + nodeBuf[i].size();

    if (!nodeBuf[i].empty()) {
      ++nRankSend;
      sendToRank[i] = 1;
    }
  }

  int nRankRecv;
  vector<int> sizeRecv(size, 1);
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /*--- Send out the messages with the global node numbers. Use nonblocking
        sends to avoid deadlock. ---*/
  vector<SU2_MPI::Request> sendReqs(nRankSend);
  nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (!nodeBuf[i].empty()) {
      SU2_MPI::Isend(nodeBuf[i].data(), nodeBuf[i].size(), MPI_UNSIGNED_LONG, i, i, SU2_MPI::GetComm(),
                     &sendReqs[nRankSend]);
      ++nRankSend;
    }
  }

  /* Define the communication buffer for the coordinates and the vector
     for the return communication requests. */
  vector<SU2_MPI::Request> returnReqs(nRankRecv);
  vector<vector<su2double> > coorReturnBuf(nRankRecv, vector<su2double>(0));

  /*--- Loop over the number of ranks from which this rank receives global
        point numbers that should be stored on this rank. ---*/
  for (int i = 0; i < nRankRecv; ++i) {
    /* Block until a message arrives. Determine the source and size
       of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for a buffer to receive this message and also
       for the buffer to return to coordinates. */
    vector<unsigned long> nodeRecvBuf(sizeMess);
    coorReturnBuf[i].resize(nDim * sizeMess);

    /* Receive the message using a blocking receive. */
    SU2_MPI::Recv(nodeRecvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, source, rank, SU2_MPI::GetComm(), &status);

    /*--- Loop over the nodes just received and fill the return communication
          buffer with the coordinates of the requested nodes. ---*/
    for (int j = 0; j < sizeMess; ++j) {
      const int jj = nDim * j;
      const long kk = nodeRecvBuf[j] - nPointsPerRank[rank];
      if (kk < 0 || kk >= nPointsRead)
        SU2_MPI::Error("Invalid point requested. This should not happen.", CURRENT_FUNCTION);

      for (unsigned short k = 0; k < nDim; ++k) coorReturnBuf[i][jj + k] = coorBuf[k][kk];
    }

    /* Send the buffer just filled back to the requesting rank.
       Use a non-blocking send to avoid deadlock. */
    SU2_MPI::Isend(coorReturnBuf[i].data(), coorReturnBuf[i].size(), MPI_DOUBLE, source, source + 1, SU2_MPI::GetComm(),
                   &returnReqs[i]);
  }

  /* Loop over the ranks from which this rank has requested coordinates. */
  for (int i = 0; i < nRankSend; ++i) {
    /* Block until a message arrives. Determine the source of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank + 1, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    /* Allocate the memory for the coordinate receive buffer. */
    vector<su2double> coorRecvBuf(nDim * nodeBuf[source].size());

    /* Receive the message using a blocking receive. */
    SU2_MPI::Recv(coorRecvBuf.data(), coorRecvBuf.size(), MPI_DOUBLE, source, rank + 1, SU2_MPI::GetComm(), &status);

    /*--- Make a distinction between 2D and 3D to store the data of the nodes.
          This data is created by taking the offset of the source rank into
          account. In this way the nodes are numbered with increading
          global node ID. ---*/
    for (unsigned long j = 0; j < nodeBuf[source].size(); ++j) {
      const unsigned long jj = nDim * j;
      const unsigned long kk = startingIndRanksInNode[source] + j;

      nodes->SetCoord(kk, &coorRecvBuf[jj]);
      nodes->SetGlobalIndex(kk, nodeBuf[source][j]);
    }
  }

  /* Complete the non-blocking sends of both rounds. */
  SU2_MPI::Waitall(sendReqs.size(), sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Waitall(returnReqs.size(), returnReqs.data(), MPI_STATUSES_IGNORE);

  /* Wild cards have been used in the communication,
     so synchronize the ranks to avoid problems.    */
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#else
  /*--- Sequential mode. Create the data for the points. The global
        number of points equals the local number of points. ---*/
  nPoint = Global_nPoint;
  nodes = new CPoint(nPoint, nDim);

  for (unsigned long i = 0; i < nPoint; ++i) {
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) nodes->SetCoord(i, iDim, coorBuf[iDim][i]);
    nodes->SetGlobalIndex(i, i);
  }

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Determine and distribute the single faces of the elements. These   ---*/
  /*--- faces are distributed over the ranks such that later the boundary  ---*/
  /*--- can retrieve the info on which rank they must be stored without    ---*/
  /*--- each rank having to read the entire connectivity data.             ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Determine the faces of the local elements. --- */
  vector<CFaceOfElement> localFaces;
  for (unsigned long k = 0; k < nElem; ++k) {
    /*--- Get the global IDs of the corner points of all the faces of this elements. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long faceConn[6][4];

    elem[k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Loop over the faces and add them to localFaces. For consistency
          between sequential and parallel mode the rank is stored at the
          position for the second element ID. ---*/
    for (unsigned short i = 0; i < nFaces; ++i) {
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for (unsigned short j = 0; j < nPointsPerFace[i]; ++j) thisFace.cornerPoints[j] = faceConn[i][j];
      thisFace.elemID0 = k + beg_node[rank];
      thisFace.elemID1 = rank;

      thisFace.CreateUniqueNumbering();
      localFaces.push_back(thisFace);
    }
  }

  /*--- Sort localFaces in increasing order and remove the double entities,
        such that unnecessary data is not communicated later on. ---*/
  sort(localFaces.begin(), localFaces.end());
  vector<CFaceOfElement>::iterator lastFace;
  lastFace = unique(localFaces.begin(), localFaces.end());
  localFaces.erase(lastFace, localFaces.end());

#ifdef HAVE_MPI

  /*--- In parallel mode these faces must be distributed over the ranks.
        A face is stored on the rank where its first node ID is located
        based on nPointsPerRank. Define the communication buffers and
        determine their contents. ---*/
  vector<vector<unsigned long> > faceBuf(size, vector<unsigned long>(0));
  for (unsigned long i = 0; i < localFaces.size(); ++i) {
    const cgsize_t nodeID = localFaces[i].cornerPoints[0];
    vector<cgsize_t>::iterator low;
    low = lower_bound(nPointsPerRank.begin(), nPointsPerRank.end(), nodeID);
    cgsize_t rankNode = low - nPointsPerRank.begin();
    if (*low > nodeID) --rankNode;

    faceBuf[rankNode].push_back(localFaces[i].nCornerPoints);
    for (unsigned short j = 0; j < localFaces[i].nCornerPoints; ++j)
      faceBuf[rankNode].push_back(localFaces[i].cornerPoints[j]);
    faceBuf[rankNode].push_back(localFaces[i].elemID0);
  }

  /* Delete the memory of localFaces again, because its contents will be
     build from the messages that this rank will receive. */
  localFaces.clear();

  /*--- Determine the number of messages this rank will send and receive. ---*/
  nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (!faceBuf[i].empty()) {
      ++nRankSend;
      sendToRank[i] = 1;
    } else {
      sendToRank[i] = 0;
    }
  }

  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /*--- Send the messages using non-blocking sends to avoid deadlock. ---*/
  sendReqs.resize(nRankSend);
  nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (!faceBuf[i].empty()) {
      SU2_MPI::Isend(faceBuf[i].data(), faceBuf[i].size(), MPI_UNSIGNED_LONG, i, i + 4, SU2_MPI::GetComm(),
                     &sendReqs[nRankSend]);
      ++nRankSend;
    }
  }

  /* Loop over the number of ranks from which this rank will receive data. */
  for (int i = 0; i < nRankRecv; ++i) {
    /* Block until a message arrives and determine the source and size
       of the message. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank + 4, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for the receive buffer and receive the
       message using a non-blocking receive. */
    vector<unsigned long> faceRecvBuf(sizeMess);
    SU2_MPI::Recv(faceRecvBuf.data(), faceRecvBuf.size(), MPI_UNSIGNED_LONG, source, rank + 4, SU2_MPI::GetComm(),
                  &status);

    /* Loop to extract the data from the receive buffer. */
    int ii = 0;
    while (ii < sizeMess) {
      /* Store the data for this face in localFaces. The rank where the
         corresponding element is physically present is stored in the
         second element ID. Note that it is not necessary to create a unique
         numbering anymore, because this has already been done before the
         communication buffer was created. */
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = (unsigned short)faceRecvBuf[ii++];
      for (unsigned short j = 0; j < thisFace.nCornerPoints; ++j, ++ii) thisFace.cornerPoints[j] = faceRecvBuf[ii];
      thisFace.elemID0 = faceRecvBuf[ii++];

      thisFace.elemID1 = source;
      localFaces.push_back(thisFace);
    }
  }

  /*--- Sort localFaces in increasing order and remove the double entities,
       such that searching is a bit more efficient later on. ---*/
  sort(localFaces.begin(), localFaces.end());
  lastFace = unique(localFaces.begin(), localFaces.end());
  localFaces.erase(lastFace, localFaces.end());

  /* Complete the non-blocking sends. Afterwards, synchronize the ranks,
     because wild cards have been used. */
  SU2_MPI::Waitall(sendReqs.size(), sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif

  /*--------------------------------------------------------------------------*/
  /*---           Reading and distributing the surface elements.           ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine the number of families in this base and read their names.
     Note that when multiple zones are present, this step is repeated for
     every zone. */
  int nFamilies;
  if (cg_nfamilies(fn, iBase, &nFamilies) != CG_OK) cg_error_exit();

  vector<string> familyNames(nFamilies);
  for (int i = 1; i <= nFamilies; ++i) {
    int nFamBC, nGeo;
    if (cg_family_read(fn, iBase, i, cgnsname, &nFamBC, &nGeo) != CG_OK) cg_error_exit();
    familyNames[i - 1] = cgnsname;
  }

  /* Determine the number of boundary conditions for this zone. */
  int nBCs;
  if (cg_nbocos(fn, iBase, iZone, &nBCs) != CG_OK) cg_error_exit();

  /* Read the names of the boundary conditions and determine their family names.
     If not family name is specified for a boundary condition, the family name
     is set to the name of the boundary condition. */
  vector<string> BCNames(nBCs), BCFamilyNames(nBCs);
  for (int i = 1; i <= nBCs; ++i) {
    /* Read the info for this boundary condition. */
    BCType_t BCType;
    PointSetType_t ptsetType;
    cgsize_t npnts, NormalListSize;
    int NormalIndex, nDataSet;
    DataType_t NormalDataType;
    if (cg_boco_info(fn, iBase, iZone, i, cgnsname, &BCType, &ptsetType, &npnts, &NormalIndex, &NormalListSize,
                     &NormalDataType, &nDataSet) != CG_OK)
      cg_error_exit();
    BCNames[i - 1] = cgnsname;

    /* Read the possibly family name and set it. If not present, it is
       equal to BCName. */
    if (cg_goto(fn, iBase, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t", i, "end") != CG_OK) cg_error_exit();

    int ierr = cg_famname_read(cgnsname);
    if (ierr == CG_ERROR)
      cg_error_exit();
    else if (ierr == CG_OK)
      BCFamilyNames[i - 1] = cgnsname;
    else
      BCFamilyNames[i - 1] = BCNames[i - 1];
  }

  /*--- Determine the number of different surface connectivities. It is
        possible to specify a family name for a surface connectivity. If that
        family name is the same for multiple surface connectivities, these
        connectivities are merged together for the boundary condition
        treatment in the DG-FEM solver. ---*/
  vector<string> surfaceNames;
  vector<vector<int> > surfaceConnIDs;

  for (int i = 0; i < nsections; ++i) {
    if (CGNSElemTypes[i].surfaceConn) {
      /*--- Determine the surface name to use for this connectivity.
            This name is determined as follows (in terms of importance).
            1) Family name specified for this connectivity.
            2) Family name of the corresponding boundary condition,
               if present.
            3) Name of the connectivity. ---*/
      string thisSurfaceName;
      const int connID = CGNSElemTypes[i].connID;

      /* First try to read the family name for the connectivity. */
      if (cg_goto(fn, iBase, "Zone_t", iZone, "Elements_t", connID, "end") != CG_OK) cg_error_exit();
      int ierr = cg_famname_read(cgnsname);
      if (ierr == CG_ERROR)
        cg_error_exit();
      else if (ierr == CG_OK)
        thisSurfaceName = cgnsname;
      else {
        /* No family name. Check the boundary conditions. It is assumed that
           the boundary conditions have the same name as the connectivities. */
        int j;
        for (j = 0; j < nBCs; ++j) {
          if (BCNames[j] == CGNSElemTypes[i].connName) {
            thisSurfaceName = BCFamilyNames[j];
            break;
          }
        }

        /* If the name is not found in the boundary conditions, set the name
           of this surface connectivity to the name of this connectivity. */
        if (j == nBCs) thisSurfaceName = CGNSElemTypes[i].connName;
      }

      /* Loop over the previously stored surface names and check if this
         surface name is already present. */
      unsigned long j;
      for (j = 0; j < surfaceNames.size(); ++j) {
        if (thisSurfaceName == surfaceNames[j]) {
          surfaceConnIDs[j].push_back(i);
        }
      }

      /* If the surface name is not stored yet, create new entries in
          surfaceNames and surfaceConnIDs. */
      if (j == surfaceNames.size()) {
        surfaceNames.push_back(thisSurfaceName);
        vector<int> thisSurfaceConn(1, i);
        surfaceConnIDs.push_back(thisSurfaceConn);
      }
    }
  }

  /* Write a message about the number of surface markers and allocate the
     memory for the data structures to store the required information. */
  nMarker = surfaceNames.size();
  if (rank == MASTER_NODE) cout << nMarker << " surface markers." << endl;
  config->SetnMarker_All(nMarker);

  unsigned short nMarker_Max = config->GetnMarker_Max();

  bound = new CPrimalGrid**[nMarker];
  nElem_Bound = new unsigned long[nMarker];
  Tag_to_Marker = new string[nMarker_Max];

  /* Loop over the number of markers to read and distribute the connectivities. */
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    /* Easier storage of the entries in CGNSElemTypes that contribute
       to this boundary marker. */
    const int nEntries = surfaceConnIDs[iMarker].size();
    const int* entries = surfaceConnIDs[iMarker].data();

    /* Determine the global number of elements for this boundary marker. */
    cgsize_t nElem_Bound_Global = 0;
    for (int iConn = 0; iConn < nEntries; ++iConn) nElem_Bound_Global += CGNSElemTypes[entries[iConn]].nElem;

    /* Write a message about the global number of surface elements
       present in this marker. */
    string Marker_Tag = surfaceNames[iMarker];
    if (rank == MASTER_NODE)
      cout << nElem_Bound_Global << " boundary elements in index " << iMarker << " (Marker = " << Marker_Tag << ")."
           << endl;

    /* Determine the number of surface elements per rank in cumulative storage
       format. This is done to avoid that every rank reads all the elements.
       This is to avoid that every rank reads all the elements. The correct
       rank for storage is determined later via communication. */
    unsigned long totalBoundElemAccounted = 0;
    vector<unsigned long> nBoundElemPerRank(size + 1);
    for (int i = 1; i <= size; ++i) {
      nBoundElemPerRank[i] = nElem_Bound_Global / size;
      totalBoundElemAccounted += nBoundElemPerRank[i];
    }

    const unsigned long nBoundElemRem = nElem_Bound_Global - totalBoundElemAccounted;
    for (unsigned long i = 1; i <= nBoundElemRem; ++i) ++nBoundElemPerRank[i];

    nBoundElemPerRank[0] = 0;
    for (int i = 0; i < size; ++i) nBoundElemPerRank[i + 1] += nBoundElemPerRank[i];

    /* Define a vector of FEM boundary elements to store the local
       boundary faces to be read. */
    vector<CBoundaryFace> boundElems;

    /* Loop over the connectivity sections that contribute. */
    elemCount = locElemCount = 0;
    for (int iConn = 0; iConn < nEntries; ++iConn) {
      /* Determine the global range for this connectivity. */
      const unsigned long elemCountOld = elemCount;
      elemCount += CGNSElemTypes[entries[iConn]].nElem;

      /* Check for overlap with the element range this rank is responsible for. */
      const unsigned long indBegOverlap = max(elemCountOld, nBoundElemPerRank[rank]);
      const unsigned long indEndOverlap = min(elemCount, nBoundElemPerRank[rank + 1]);

      if (indEndOverlap > indBegOverlap) {
        /* This rank must read boundary element data from this connectivity
           section. Determine the offset relative to the start of this section
           and the number of elements to be read by this rank. */
        const unsigned long offsetRank = indBegOverlap - elemCountOld;
        const unsigned long nElemRank = indEndOverlap - indBegOverlap;

        /* Read the connectivity range determined above. */
        CGNSElemTypes[entries[iConn]].ReadBoundaryConnectivityRange(fn, iBase, iZone, offsetRank, nElemRank,
                                                                    nBoundElemPerRank[rank], locElemCount, boundElems);
      }
    }

    /* Make a distinction between sequential and parallel mode. */
#ifdef HAVE_MPI
    /*--- Parallel mode. The information stored in boundElems must be
          communicated to find out where it must be stored. First clear the
          contents of the faceBuf, such that it can be used again to send
          data to the appropiate rank. ---*/
    for (int i = 0; i < size; ++i) faceBuf[i].clear();

    /*-- Loop over the locally read boundary elements and store its contents
         in the appropiate location in faceBuf. ---*/
    for (unsigned long i = 0; i < boundElems.size(); ++i) {
      /* Create an object of the class CFaceOfElement to determine the
         node number that determines on which rank the corresponding face
         in localFaces is stored. */
      CFaceOfElement thisFace(boundElems[i].VTK_Type, boundElems[i].nPolyGrid, boundElems[i].Nodes.data());
      thisFace.CreateUniqueNumbering();

      /* Determine the rank to which this face must be sent to. */
      const cgsize_t nodeID = thisFace.cornerPoints[0];
      vector<cgsize_t>::iterator low;
      low = lower_bound(nPointsPerRank.begin(), nPointsPerRank.end(), nodeID);
      cgsize_t rankNode = low - nPointsPerRank.begin();
      if (*low > nodeID) --rankNode;

      /*--- Copy the relevant data of this boundary element into faceBuf. ---*/
      faceBuf[rankNode].push_back(boundElems[i].VTK_Type);
      faceBuf[rankNode].push_back(boundElems[i].nPolyGrid);
      faceBuf[rankNode].push_back(boundElems[i].nDOFsGrid);
      faceBuf[rankNode].push_back(boundElems[i].globalBoundElemID);
      faceBuf[rankNode].insert(faceBuf[rankNode].end(), boundElems[i].Nodes.begin(), boundElems[i].Nodes.end());
    }

    /* The contents of boundElems is copied into faceBuf, so it can
       be released. */
    boundElems.clear();

    /*--- Determine the number of messages this rank will send and receive. ---*/
    nRankSend = 0;
    for (int i = 0; i < size; ++i) {
      if (!faceBuf[i].empty()) {
        ++nRankSend;
        sendToRank[i] = 1;
      } else {
        sendToRank[i] = 0;
      }
    }

    SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

    /*--- Send the messages using non-blocking sends to avoid deadlock. ---*/
    sendReqs.resize(nRankSend);
    nRankSend = 0;
    for (int i = 0; i < size; ++i) {
      if (!faceBuf[i].empty()) {
        SU2_MPI::Isend(faceBuf[i].data(), faceBuf[i].size(), MPI_UNSIGNED_LONG, i, i + 5, SU2_MPI::GetComm(),
                       &sendReqs[nRankSend]);
        ++nRankSend;
      }
    }

    /* Use nodeBuf as storage for the sending of the surface element data
       to the correct rank. First clear its contents. */
    for (int i = 0; i < size; ++i) nodeBuf[i].clear();

    /*--- Loop over the number of ranks from which this rank receives
          surface elements to be processed. ---*/
    for (int i = 0; i < nRankRecv; ++i) {
      /* Block until a message arrives. Determine the source and size
         of the message. */
      SU2_MPI::Status status;
      SU2_MPI::Probe(MPI_ANY_SOURCE, rank + 5, SU2_MPI::GetComm(), &status);
      int source = status.MPI_SOURCE;

      int sizeMess;
      SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

      /* Allocate the memory for the receive buffer and receive the message
         using a blocking send. */
      vector<unsigned long> boundElemRecvBuf(sizeMess);
      SU2_MPI::Recv(boundElemRecvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, source, rank + 5, SU2_MPI::GetComm(),
                    &status);

      /* Loop to extract the data from the receive buffer. */
      int ii = 0;
      while (ii < sizeMess) {
        /* Store the data for this boundary element. */
        const auto VTK_Type = (unsigned short)boundElemRecvBuf[ii++];
        const auto nPolyGrid = (unsigned short)boundElemRecvBuf[ii++];
        const auto nDOFsGrid = (unsigned short)boundElemRecvBuf[ii++];

        const unsigned long globalBoundElemID = boundElemRecvBuf[ii++];
        const unsigned long* Nodes = boundElemRecvBuf.data() + ii;
        ii += nDOFsGrid;

        /* Determine the corner nodes and store them in an object of
           the class CFaceOfElement to carry out the search in localFaces. */
        CFaceOfElement thisFace(VTK_Type, nPolyGrid, Nodes);

        /* Check if the face is actually present. If not, print an error message
           and exit. */
        thisFace.CreateUniqueNumbering();
        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        bool thisFaceFound = false;
        if (low != localFaces.end()) {
          if (!(thisFace < *low)) thisFaceFound = true;
        }

        if (!thisFaceFound)
          SU2_MPI::Error("Boundary element not found in list of faces. This is a bug.", CURRENT_FUNCTION);

        /* Determine the domain element and the rank where
           this boundary element should be sent to.. */
        const unsigned long domainElementID = low->elemID0;
        const int rankBoundElem = (int)low->elemID1;

        /*--- Store the data for this element in the communication buffer
              for rankBoundElem. ---*/
        nodeBuf[rankBoundElem].push_back(VTK_Type);
        nodeBuf[rankBoundElem].push_back(nPolyGrid);
        nodeBuf[rankBoundElem].push_back(nDOFsGrid);
        nodeBuf[rankBoundElem].push_back(globalBoundElemID);
        nodeBuf[rankBoundElem].push_back(domainElementID);

        for (unsigned short j = 0; j < nDOFsGrid; ++j) nodeBuf[rankBoundElem].push_back(Nodes[j]);
      }
    }

    /* Complete the non-blocking sends. */
    SU2_MPI::Waitall(sendReqs.size(), sendReqs.data(), MPI_STATUSES_IGNORE);

    /*--- Determine the number of messages this rank will send and receive. ---*/
    nRankSend = 0;
    for (int i = 0; i < size; ++i) {
      if (!nodeBuf[i].empty()) {
        ++nRankSend;
        sendToRank[i] = 1;
      } else {
        sendToRank[i] = 0;
      }
    }

    SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeRecv.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

    /*--- Send the messages using non-blocking sends to avoid deadlock. ---*/
    sendReqs.resize(nRankSend);
    nRankSend = 0;
    for (int i = 0; i < size; ++i) {
      if (!nodeBuf[i].empty()) {
        SU2_MPI::Isend(nodeBuf[i].data(), nodeBuf[i].size(), MPI_UNSIGNED_LONG, i, i + 6, SU2_MPI::GetComm(),
                       &sendReqs[nRankSend]);
        ++nRankSend;
      }
    }

    /*--- Loop over the number of ranks from which this rank receives
          surface elements to be stored on this rank. ---*/
    for (int i = 0; i < nRankRecv; ++i) {
      /* Block until a message arrives. Determine the source and size
         of the message. */
      SU2_MPI::Status status;
      SU2_MPI::Probe(MPI_ANY_SOURCE, rank + 6, SU2_MPI::GetComm(), &status);
      int source = status.MPI_SOURCE;

      int sizeMess;
      SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

      /* Allocate the memory for the receive buffer and receive the message
         using a blocking send. */
      vector<unsigned long> boundElemRecvBuf(sizeMess);
      SU2_MPI::Recv(boundElemRecvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, source, rank + 6, SU2_MPI::GetComm(),
                    &status);

      /* Loop to extract the data from the receive buffer. */
      int ii = 0;
      while (ii < sizeMess) {
        /* Store the data for this boundary element. */
        const auto VTK_Type = (unsigned short)boundElemRecvBuf[ii++];
        const auto nPolyGrid = (unsigned short)boundElemRecvBuf[ii++];
        const auto nDOFsGrid = (unsigned short)boundElemRecvBuf[ii++];

        const unsigned long globalBoundElemID = boundElemRecvBuf[ii++];
        const unsigned long domainElementID = boundElemRecvBuf[ii++];
        const unsigned long* Nodes = boundElemRecvBuf.data() + ii;
        ii += nDOFsGrid;

        /* Create an object of CBoundaryFace and store it in boundElems. */
        CBoundaryFace thisBoundFace;
        thisBoundFace.VTK_Type = VTK_Type;
        thisBoundFace.nPolyGrid = nPolyGrid;
        thisBoundFace.nDOFsGrid = nDOFsGrid;
        thisBoundFace.globalBoundElemID = globalBoundElemID;
        thisBoundFace.domainElementID = domainElementID;

        thisBoundFace.Nodes.resize(nDOFsGrid);
        for (unsigned short j = 0; j < nDOFsGrid; ++j) thisBoundFace.Nodes[j] = Nodes[j];

        boundElems.push_back(thisBoundFace);
      }
    }

    /* Sort boundElems in increasing order, where the less than operator
       is a comparison between the globalBoundElemID's. */
    sort(boundElems.begin(), boundElems.end());

    /* Complete the non-blocking sends and synchronize the processors,
       because wild cards have been used. */
    SU2_MPI::Waitall(sendReqs.size(), sendReqs.data(), MPI_STATUSES_IGNORE);

    SU2_MPI::Barrier(SU2_MPI::GetComm());

#else
    /*--- Sequential mode. All boundary elements read must be stored on this
          rank. The only information missing is the global ID of the domain
          element of the boundary elements. This is created below. ---*/
    for (unsigned long i = 0; i < boundElems.size(); ++i) {
      /* Determine the corner nodes and store them in an object of
         the class CFaceOfElement to carry out the search in localFaces. */
      CFaceOfElement thisFace(boundElems[i].VTK_Type, boundElems[i].nPolyGrid, boundElems[i].Nodes.data());

      /* Check if the face is actually present. If not, print an error message
         and exit. */
      thisFace.CreateUniqueNumbering();
      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

      bool thisFaceFound = false;
      if (low != localFaces.end()) {
        if (!(thisFace < *low)) thisFaceFound = true;
      }

      if (!thisFaceFound)
        SU2_MPI::Error("Boundary element not found in list of faces. This is a bug.", CURRENT_FUNCTION);

      /* Set the domain element. */
      boundElems[i].domainElementID = low->elemID0;
    }
#endif

    /*--- Allocate space for the local boundary elements and copy the data
          from boundElems into bound. ---*/
    nElem_Bound[iMarker] = boundElems.size();
    bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];

    for (unsigned long i = 0; i < nElem_Bound[iMarker]; ++i)
      bound[iMarker][i] = new CPrimalGridBoundFEM(boundElems[i].globalBoundElemID, boundElems[i].domainElementID,
                                                  boundElems[i].VTK_Type, boundElems[i].nPolyGrid,
                                                  boundElems[i].nDOFsGrid, boundElems[i].Nodes);

    /*--- Update config information storing the boundary information in the right place ---*/
    Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
    config->SetMarker_All_TagBound(iMarker, Marker_Tag);
    config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
    config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
    config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
    config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
    config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
    config->SetMarker_All_Analyze(iMarker, config->GetMarker_CfgFile_Analyze(Marker_Tag));
    config->SetMarker_All_ZoneInterface(iMarker, config->GetMarker_CfgFile_ZoneInterface(Marker_Tag));
    config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
    config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
    config->SetMarker_All_SobolevBC(iMarker, config->GetMarker_CfgFile_SobolevBC(Marker_Tag));
    config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
    config->SetMarker_All_SendRecv(iMarker, NONE);
  }

  /* Close the CGNS file again. */
  if (cg_close(fn) != CG_OK) cg_error_exit();
  if (rank == MASTER_NODE) cout << "Successfully closed the CGNS file." << endl;

#else /* CGNS_VERSION >= 3300 */

  SU2_MPI::Error("CGNS version 3.3 or higher is necessary for the DG FEM solver", CURRENT_FUNCTION);

#endif /* CGNS_VERSION >= 3300 */

#else /* HAVE_CGNS. */

  SU2_MPI::Error("SU2 built without CGNS support!!\nTo use CGNS, build SU2 accordingly.", CURRENT_FUNCTION);

#endif /* HAVE_CGNS. */
}

void CPhysicalGeometry::SetColorFEMGrid_Parallel(CConfig* config) {
  /*--- Initialize the color vector of the elements. ---*/
  for (unsigned long i = 0; i < nElem; ++i) elem[i]->SetColor(0);

  /*--- Determine the faces of the elements. ---*/
  vector<CFaceOfElement> localFaces;
  for (unsigned long k = 0; k < nElem; k++) {
    /*--- Get the global IDs of the corner points of all the faces of this element. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long faceConn[6][4];

    elem[k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Loop over the faces and add them to localFaces. ---*/
    for (unsigned short i = 0; i < nFaces; ++i) {
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for (unsigned short j = 0; j < nPointsPerFace[i]; ++j) thisFace.cornerPoints[j] = faceConn[i][j];
      thisFace.elemID0 = beg_node[rank] + k;
      thisFace.nPolySol0 = elem[k]->GetNPolySol();
      thisFace.nDOFsElem0 = elem[k]->GetNDOFsSol();
      thisFace.elemType0 = elem[k]->GetVTK_Type();

      thisFace.CreateUniqueNumbering();
      localFaces.push_back(thisFace);
    }
  }

  /*--- Sort localFaces in increasing order. */
  sort(localFaces.begin(), localFaces.end());

  /*--- Loop over the boundary markers and its faces in order to flag the
        physical boundary faces in localFaces. Note that use is made of the
        overloaded function GetCornerPointsAllFaces, which explains the
        dimensions of the variables used in the function call. Also note that
        the periodic boundaries are excluded, because they are not physical.  ---*/
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    if (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {
      for (unsigned long k = 0; k < nElem_Bound[iMarker]; ++k) {
        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long faceConn[6][4];
        bound[iMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for (unsigned short j = 0; j < nPointsPerFace[0]; ++j) thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        /*--- Store the corresponding element ID for the boundary element and
              invalidate the face by setting the element ID to an invalid value. ---*/
        bound[iMarker][k]->SetDomainElement(low->elemID0);
        low->elemID0 = Global_nElem + 10;
      }
    }
  }

  /*--- Loop over the faces and check for double entries. If a double entry is
        found, the elemID from the second entry is copied to the first entry,
        the polynomial degree is updated, and the second entry is invalidated. ---*/
  unsigned long nFacesLoc = localFaces.size();
  for (unsigned long i = 1; i < nFacesLoc; ++i) {
    if (localFaces[i] == localFaces[i - 1]) {
      localFaces[i - 1].elemID1 = localFaces[i].elemID0;
      localFaces[i - 1].nPolySol1 = localFaces[i].nPolySol0;
      localFaces[i - 1].nDOFsElem1 = localFaces[i].nDOFsElem0;
      localFaces[i - 1].elemType1 = localFaces[i].elemType0;
      localFaces[i].elemID0 = Global_nElem + 10;
    }
  }

  /*--- Remove the invalidated faces. This is accomplished by giving the
        face four points a global node ID that is larger than the largest
        point ID in the grid. In this way the sorting operator puts
        these faces at the end of the vector, see also the < operator
        of CFaceOfElement.                                         ---*/
  unsigned long nFacesLocOr = nFacesLoc;
  for (unsigned long i = 0; i < nFacesLocOr; ++i) {
    if (localFaces[i].elemID0 > Global_nElem) {
      localFaces[i].nCornerPoints = 4;
      localFaces[i].cornerPoints[0] = Global_nPoint;
      localFaces[i].cornerPoints[1] = Global_nPoint;
      localFaces[i].cornerPoints[2] = Global_nPoint;
      localFaces[i].cornerPoints[3] = Global_nPoint;
      --nFacesLoc;
    }
  }

  sort(localFaces.begin(), localFaces.end());
  localFaces.resize(nFacesLoc);

  /*--- The following part of this function should only be present if we have
        parallel support with MPI. ---*/

#ifdef HAVE_MPI

  /*--- Determine and store the faces, for which only one neighboring element
        was found. For these faces the other neighbor might be stored on
        a different rank (unless there are non-matching interfaces).     ---*/
  vector<CFaceOfElement> localFacesComm;
  for (unsigned long i = 0; i < nFacesLoc; ++i)
    if (localFaces[i].elemID1 > Global_nElem) localFacesComm.push_back(localFaces[i]);

  /*--- Determine the maximum global point ID that occurs in localFacesComm
        of all ranks. Note that only the first point is taken into account,
        because this point determines the rank where the face is sent to. ---*/
  unsigned long nFacesLocComm = localFacesComm.size();
  unsigned long maxPointIDLoc = 0;
  for (unsigned long i = 0; i < nFacesLocComm; ++i)
    maxPointIDLoc = max(maxPointIDLoc, localFacesComm[i].cornerPoints[0]);

  unsigned long maxPointID;
  SU2_MPI::Allreduce(&maxPointIDLoc, &maxPointID, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
  ++maxPointID;

  /*--- Create a vector with a linear distribution over the ranks for
        the points that occur in the faces of localFacesComm.        ---*/
  vector<unsigned long> facePointsProc(size + 1, 0);
  unsigned long total_point_accounted = 0;
  for (unsigned long i = 1; i <= (unsigned long)size; ++i) {
    facePointsProc[i] = maxPointID / size;
    total_point_accounted += facePointsProc[i];
  }

  unsigned long rem_point = maxPointID - total_point_accounted;
  for (unsigned long i = 1; i <= rem_point; ++i) ++facePointsProc[i];

  for (unsigned long i = 0; i < (unsigned long)size; ++i) facePointsProc[i + 1] += facePointsProc[i];

  /*--- Determine the number of faces that has to be sent to each rank.
        Note that the rank is stored in elemID1, such that the search
        does not have to be repeated below.               ---*/
  vector<unsigned long> nFacesComm(size, 0);
  for (unsigned long i = 0; i < nFacesLocComm; ++i) {
    vector<unsigned long>::iterator low;
    low = lower_bound(facePointsProc.begin(), facePointsProc.end(), localFacesComm[i].cornerPoints[0]);
    unsigned long rankFace = low - facePointsProc.begin();
    if (*low > localFacesComm[i].cornerPoints[0]) --rankFace;

    ++nFacesComm[rankFace];
    localFacesComm[i].elemID1 = rankFace;
  }

  /*--- Create the send buffer for the faces to be communicated. ---*/
  vector<unsigned long> sendBufFace(9 * nFacesLocComm);
  vector<unsigned long> counter(size);
  counter[0] = 0;
  for (unsigned long i = 1; i < (unsigned long)size; ++i) counter[i] = counter[i - 1] + 9 * nFacesComm[i - 1];

  for (unsigned long i = 0; i < nFacesLocComm; ++i) {
    unsigned long rankFace = localFacesComm[i].elemID1;
    unsigned long ii = counter[rankFace];
    counter[rankFace] += 9;

    sendBufFace[ii] = localFacesComm[i].nCornerPoints;
    sendBufFace[ii + 1] = localFacesComm[i].cornerPoints[0];
    sendBufFace[ii + 2] = localFacesComm[i].cornerPoints[1];
    sendBufFace[ii + 3] = localFacesComm[i].cornerPoints[2];
    sendBufFace[ii + 4] = localFacesComm[i].cornerPoints[3];
    sendBufFace[ii + 5] = localFacesComm[i].elemID0;
    sendBufFace[ii + 6] = localFacesComm[i].nPolySol0;
    sendBufFace[ii + 7] = localFacesComm[i].nDOFsElem0;
    sendBufFace[ii + 8] = localFacesComm[i].elemType0;
  }

  /*--- Determine the number of ranks from which I receive a message. */
  unsigned long nMessSend = 0;
  for (unsigned long i = 0; i < (unsigned long)size; ++i) {
    if (nFacesComm[i]) {
      counter[i] = 1;
      ++nMessSend;
    } else
      counter[i] = 0;
  }
  vector<int> sizeRecv(size, 1);

  unsigned long nMessRecv;
  SU2_MPI::Reduce_scatter(counter.data(), &nMessRecv, sizeRecv.data(), MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  /*--- Send the data using nonblocking sends. ---*/
  vector<SU2_MPI::Request> commReqs(max(nMessSend, nMessRecv));

  nMessSend = 0;
  unsigned long indSend = 0;
  for (unsigned long i = 0; i < (unsigned long)size; ++i) {
    if (nFacesComm[i]) {
      unsigned long count = 9 * nFacesComm[i];
      SU2_MPI::Isend(&sendBufFace[indSend], count, MPI_UNSIGNED_LONG, i, i, SU2_MPI::GetComm(), &commReqs[nMessSend]);
      ++nMessSend;
      indSend += count;
    }
  }

  /*--- Loop over the number of ranks from which faces are received.
        Receive the messages and store them in facesRecv.         ---*/
  vector<CFaceOfElement> facesRecv;
  vector<unsigned long> nFacesRecv(nMessRecv + 1);
  vector<int> rankRecv(nMessRecv);
  nFacesRecv[0] = 0;
  for (unsigned long i = 0; i < nMessRecv; ++i) {
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    rankRecv[i] = status.MPI_SOURCE;
    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, rankRecv[i], rank, SU2_MPI::GetComm(), &status);

    nFacesRecv[i + 1] = nFacesRecv[i] + sizeMess / 9;
    facesRecv.resize(nFacesRecv[i + 1]);
    unsigned long ii = 0;
    for (unsigned long j = nFacesRecv[i]; j < nFacesRecv[i + 1]; ++j, ii += 9) {
      facesRecv[j].nCornerPoints = recvBuf[ii];
      facesRecv[j].cornerPoints[0] = recvBuf[ii + 1];
      facesRecv[j].cornerPoints[1] = recvBuf[ii + 2];
      facesRecv[j].cornerPoints[2] = recvBuf[ii + 3];
      facesRecv[j].cornerPoints[3] = recvBuf[ii + 4];
      facesRecv[j].elemID0 = recvBuf[ii + 5];
      facesRecv[j].nPolySol0 = recvBuf[ii + 6];
      facesRecv[j].nDOFsElem0 = recvBuf[ii + 7];
      facesRecv[j].elemType0 = recvBuf[ii + 8];
    }
  }

  /*--- The exact contents of facesRecv is needed when communicating back
        the data to the original ranks, so a copy is made of facesRecv.
        As localFacesComm is not needed anymore, this vector is used
        for the sorting and searching of the data of facesRecv. ---*/
  localFacesComm = facesRecv;
  sort(localFacesComm.begin(), localFacesComm.end());

  nFacesLocComm = localFacesComm.size();
  for (unsigned long i = 1; i < localFacesComm.size(); ++i) {
    if (localFacesComm[i] == localFacesComm[i - 1]) {
      localFacesComm[i - 1].elemID1 = localFacesComm[i].elemID0;
      localFacesComm[i - 1].nPolySol1 = localFacesComm[i].nPolySol0;
      localFacesComm[i - 1].nDOFsElem1 = localFacesComm[i].nDOFsElem0;
      localFacesComm[i - 1].elemType1 = localFacesComm[i].elemType0;

      localFacesComm[i].nCornerPoints = 4;
      localFacesComm[i].cornerPoints[0] = Global_nPoint;
      localFacesComm[i].cornerPoints[1] = Global_nPoint;
      localFacesComm[i].cornerPoints[2] = Global_nPoint;
      localFacesComm[i].cornerPoints[3] = Global_nPoint;
      --nFacesLocComm;
    }
  }

  sort(localFacesComm.begin(), localFacesComm.end());
  localFacesComm.resize(nFacesLocComm);

  /*--- Complete the first round of non-blocking sends. ---*/
  SU2_MPI::Waitall(nMessSend, commReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Send the data back to the requesting ranks. ---*/
  sendBufFace.resize(9 * nFacesRecv[nMessRecv]);
  indSend = 0;
  for (unsigned long i = 0; i < nMessRecv; ++i) {
    unsigned long ii = indSend;
    for (unsigned long j = nFacesRecv[i]; j < nFacesRecv[i + 1]; ++j, ii += 9) {
      sendBufFace[ii] = facesRecv[j].nCornerPoints;
      sendBufFace[ii + 1] = facesRecv[j].cornerPoints[0];
      sendBufFace[ii + 2] = facesRecv[j].cornerPoints[1];
      sendBufFace[ii + 3] = facesRecv[j].cornerPoints[2];
      sendBufFace[ii + 4] = facesRecv[j].cornerPoints[3];

      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFacesComm.begin(), localFacesComm.end(), facesRecv[j]);
      if (facesRecv[j].elemID0 == low->elemID0) {
        sendBufFace[ii + 5] = low->elemID1;
        sendBufFace[ii + 6] = low->nPolySol1;
        sendBufFace[ii + 7] = low->nDOFsElem1;
        sendBufFace[ii + 8] = low->elemType1;
      } else {
        sendBufFace[ii + 5] = low->elemID0;
        sendBufFace[ii + 6] = low->nPolySol0;
        sendBufFace[ii + 7] = low->nDOFsElem0;
        sendBufFace[ii + 8] = low->elemType0;
      }
    }

    unsigned long count = ii - indSend;
    SU2_MPI::Isend(&sendBufFace[indSend], count, MPI_UNSIGNED_LONG, rankRecv[i], rankRecv[i] + 1, SU2_MPI::GetComm(),
                   &commReqs[i]);
    indSend = ii;
  }

  /*--- Loop over the ranks to which I originally sent my face data.
        The return data contains information about the neighboring element. ---*/
  for (unsigned long i = 0; i < nMessSend; ++i) {
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank + 1, SU2_MPI::GetComm(), &status);
    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, status.MPI_SOURCE, rank + 1, SU2_MPI::GetComm(),
                  &status);

    sizeMess /= 9;
    unsigned long jj = 0;
    for (unsigned long j = 0; j < (unsigned long)sizeMess; ++j, jj += 9) {
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = recvBuf[jj];
      thisFace.cornerPoints[0] = recvBuf[jj + 1];
      thisFace.cornerPoints[1] = recvBuf[jj + 2];
      thisFace.cornerPoints[2] = recvBuf[jj + 3];
      thisFace.cornerPoints[3] = recvBuf[jj + 4];

      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);
      low->elemID1 = recvBuf[jj + 5];
      low->nPolySol1 = recvBuf[jj + 6];
      low->nDOFsElem1 = recvBuf[jj + 7];
      low->elemType1 = recvBuf[jj + 8];
    }
  }

  /*--- Complete the second round of non-blocking sends. ---*/
  SU2_MPI::Waitall(nMessRecv, commReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Wild cards have been used in the communication, so
        synchronize the ranks to avoid problems.          ---*/
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif

  /*--- In the procedure above the periodic boundaries are not found.
        A different treatment must be used in order to find these. ---*/
  DeterminePeriodicFacesFEMGrid(config, localFaces);

  /*--- Determine the total number of non-matching faces in the grid. ---*/
  nFacesLoc = localFaces.size();
  nFacesLocOr = nFacesLoc;
  for (unsigned long i = 0; i < nFacesLocOr; ++i) {
    if (localFaces[i].elemID1 > Global_nElem) {
      localFaces[i].nCornerPoints = 4;
      localFaces[i].cornerPoints[0] = Global_nPoint;
      localFaces[i].cornerPoints[1] = Global_nPoint;
      localFaces[i].cornerPoints[2] = Global_nPoint;
      localFaces[i].cornerPoints[3] = Global_nPoint;
      --nFacesLoc;
    }
  }

  nFacesLocOr -= nFacesLoc;
  if (nFacesLocOr) {
    sort(localFaces.begin(), localFaces.end());
    localFaces.resize(nFacesLoc);
  }

  unsigned long nNonMatchingFaces = nFacesLocOr;

#ifdef HAVE_MPI
  SU2_MPI::Reduce(&nFacesLocOr, &nNonMatchingFaces, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
#endif
  if (rank == MASTER_NODE && nNonMatchingFaces) {
    cout << "There are " << nNonMatchingFaces << " non-matching faces in the grid. "
         << "These are ignored in the partitioning." << endl;
  }

  /*--- Determine whether or not the Jacobians of the elements and faces
        can be considered constant. ---*/
  DetermineFEMConstantJacobiansAndLenScale(config);

  /*--- Determine the donor elements for the wall function treatment. ---*/
  DetermineDonorElementsWallFunctions(config);

  /*--- Determine the time levels of the elements. This is only relevant
        when time accurate local time stepping is employed. ---*/
  map<unsigned long, CUnsignedShort2T> mapExternalElemIDToTimeLevel;
  DetermineTimeLevelElements(config, localFaces, mapExternalElemIDToTimeLevel);

  /*--- Determine the ownership of the internal faces, i.e. which adjacent
        element is responsible for computing the fluxes through the face. ---*/
  for (unsigned long i = 0; i < nFacesLoc; ++i) {
    /* Check for a matching face. */
    if (localFaces[i].elemID1 < Global_nElem) {
      /* Determine the time level for both elements. Elem0 is always owned, while
         elem1 is either owned or external. The data for external elements is
         stored in mapExternalElemIDToTimeLevel. */
      unsigned long elemID0 = localFaces[i].elemID0 - beg_node[rank];
      unsigned short timeLevel0 = elem[elemID0]->GetTimeLevel();

      unsigned short timeLevel1;
      if (localFaces[i].elemID1 >= beg_node[rank] && localFaces[i].elemID1 < beg_node[rank] + nElem) {
        unsigned long elemID1 = localFaces[i].elemID1 - beg_node[rank];
        timeLevel1 = elem[elemID1]->GetTimeLevel();
      } else {
        map<unsigned long, CUnsignedShort2T>::const_iterator MI;
        MI = mapExternalElemIDToTimeLevel.find(localFaces[i].elemID1);
        if (MI == mapExternalElemIDToTimeLevel.end())
          SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
        timeLevel1 = MI->second.short0;
      }

      /* Check if both elements have the same time level. */
      if (timeLevel0 == timeLevel1) {
        /* Same time level, hence both elements can own the face. First check whether
           elemID0 == elemID1 (which happens for periodic problems with only one
           element in the periodic direction), because this is a special case. */
        if (localFaces[i].elemID0 == localFaces[i].elemID1) {
          /* This face occurs twice, but should be owned only once. Base this
             decision on the periodic index. */
          localFaces[i].elem0IsOwner = localFaces[i].periodicIndex < localFaces[i].periodicIndexDonor;
        } else {
          /* Different elements on both sides of the face. The ownership decision
             below makes an attempt to spread the workload evenly. */
          const unsigned long sumElemID = localFaces[i].elemID0 + localFaces[i].elemID1;
          if (sumElemID % 2)
            localFaces[i].elem0IsOwner = localFaces[i].elemID0 < localFaces[i].elemID1;
          else
            localFaces[i].elem0IsOwner = localFaces[i].elemID0 > localFaces[i].elemID1;
        }
      } else {
        /* The time level of both elements differ. The element with the smallest
           time level must be the owner of the element. */
        localFaces[i].elem0IsOwner = timeLevel0 < timeLevel1;
      }
    } else {
      /* Non-matching face. Give the ownership to element 0. */
      localFaces[i].elem0IsOwner = true;
    }
  }

  /*--- All the matching face information is known now, including periodic
        faces. Store the information of the neighbors in the data structure
        for the local elements.     ---*/
  for (unsigned long k = 0; k < nElem; ++k) {
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long faceConn[6][4];

    elem[k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);
    elem[k]->InitializeNeighbors(nFaces);

    for (unsigned short i = 0; i < nFaces; ++i) {
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for (unsigned short j = 0; j < nPointsPerFace[i]; ++j) thisFace.cornerPoints[j] = faceConn[i][j];
      thisFace.elemID0 = beg_node[rank] + k;
      thisFace.nPolySol0 = elem[k]->GetNPolySol();
      thisFace.nDOFsElem0 = elem[k]->GetNDOFsSol();

      thisFace.CreateUniqueNumbering();

      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

      if (low != localFaces.end()) {
        if (!(thisFace < *low)) {
          if (low->elemID0 == thisFace.elemID0) {
            elem[k]->SetNeighbor_Elements(low->elemID1, i);
            elem[k]->SetOwnerFace(low->elem0IsOwner, i);
          } else {
            elem[k]->SetNeighbor_Elements(low->elemID0, i);
            elem[k]->SetOwnerFace(!(low->elem0IsOwner), i);
          }

          if (low->periodicIndex > 0) elem[k]->SetPeriodicIndex(low->periodicIndex - 1, i);
        }
      }
    }
  }

  /*--- Create the vector of vectors that describe the connectivity
        of the graph. First the faces. ---*/
  vector<vector<unsigned long> > adjacency(nElem, vector<unsigned long>(0));
  for (unsigned long i = 0; i < nFacesLoc; ++i) {
    long ii = localFaces[i].elemID0 - beg_node[rank];
    adjacency[ii].push_back(localFaces[i].elemID1);

    if (localFaces[i].periodicIndex == 0) {
      ii = localFaces[i].elemID1 - beg_node[rank];
      if (ii >= 0 && ii < (long)nElem) adjacency[ii].push_back(localFaces[i].elemID0);
    }
  }

  /* It is possible that some neighbors appear multiple times due to e.g.
     periodic boundary conditions. ParMETIS is not able to deal with this
     situation, hence these multiple entries must be removed. */
  for (unsigned long i = 0; i < nElem; ++i) {
    sort(adjacency[i].begin(), adjacency[i].end());
    vector<unsigned long>::iterator lastEntry;
    lastEntry = unique(adjacency[i].begin(), adjacency[i].end());
    adjacency[i].erase(lastEntry, adjacency[i].end());
  }

  /* Due to periodic boundary conditions it is also possible that self entries
     are present. ParMETIS is not able to deal with self entries, hence
     they must be removed as well. */
  for (unsigned long i = 0; i < nElem; ++i) {
    const unsigned long globalElemID = i + beg_node[rank];
    unsigned long nEntriesNew = adjacency[i].size();

    for (unsigned long j = 0; j < adjacency[i].size(); ++j) {
      if (adjacency[i][j] == globalElemID) {
        adjacency[i][j] = ULONG_MAX;
        --nEntriesNew;
      }
    }

    sort(adjacency[i].begin(), adjacency[i].end());
    adjacency[i].resize(nEntriesNew);
  }

  /*--- Possibly add the connectivities in the graph from the wall function
        treatment. As these connectivities are one-sided, they must be stored
        and communicated later. ---*/
  vector<unsigned long> additionalExternalEntriesGraph;

  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    for (unsigned long l = 0; l < nElem_Bound[iMarker]; ++l) {
      /* Get the global and local element ID adjacent to this boundary face. */
      const unsigned long globalElemID = bound[iMarker][l]->GetDomainElement();
      const unsigned long elemID = globalElemID - beg_node[rank];

      /* Get the number of donor elements for the wall function treatment
         and the pointer to the array which stores this info. */
      const unsigned short nDonors = bound[iMarker][l]->GetNDonorsWallFunctions();
      const unsigned long* donors = bound[iMarker][l]->GetDonorsWallFunctions();

      /* Loop over the number of donors and add the entry in the graph,
         if not already present. */
      for (unsigned short i = 0; i < nDonors; ++i) {
        if (donors[i] != globalElemID) {
          if (!binary_search(adjacency[elemID].begin(), adjacency[elemID].end(), donors[i])) {
            /* Donor not present in the graph for elemID. Add it and sort it
               afterwards, such that a binary search can be applied later. */
            adjacency[elemID].push_back(donors[i]);
            sort(adjacency[elemID].begin(), adjacency[elemID].end());

            /* Check if the donor element is stored locally. */
            if (donors[i] >= beg_node[rank] && donors[i] < beg_node[rank] + nElem) {
              /* Donor is stored locally. Add the entry to the graph
                 and sort it afterwards. */
              const unsigned long localDonorID = donors[i] - beg_node[rank];
              adjacency[localDonorID].push_back(globalElemID);
              sort(adjacency[localDonorID].begin(), adjacency[localDonorID].end());
            } else {
              /* Donor is stored externally. Store the graph entry in
                 additionalExternalEntriesGraph. */
              additionalExternalEntriesGraph.push_back(donors[i]);
              additionalExternalEntriesGraph.push_back(globalElemID);
            }
          }
        }
      }
    }
  }

#ifdef HAVE_MPI
  /*--- Create the send buffers with the additional graph data and determine
        to which ranks data is sent. ---*/
  vector<vector<unsigned long> > sendBufsGraphData(size, vector<unsigned long>(0));
  vector<int> sendToRank(size, 0);

  for (unsigned long i = 0; i < additionalExternalEntriesGraph.size(); i += 2) {
    /* Determine the rank where this external is stored. */
    const unsigned long elemID = additionalExternalEntriesGraph[i];
    int rankElem;
    if (elemID >= beg_node[size - 1])
      rankElem = size - 1;
    else {
      const unsigned long* low;
      low = lower_bound(beg_node, beg_node + size, elemID);

      rankElem = low - beg_node;
      if (*low > elemID) --rankElem;
    }

    sendBufsGraphData[rankElem].push_back(additionalExternalEntriesGraph[i]);
    sendBufsGraphData[rankElem].push_back(additionalExternalEntriesGraph[i + 1]);
    sendToRank[rankElem] = 1;
  }

  /*-- Determine to how many ranks this rank will send data and from how
       many ranks it will receive data. ---*/
  int nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (sendToRank[i]) ++nRankSend;
  }

  int nRankRecv;
  vector<int> sizeSend(size, 1);
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeSend.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /* Send the data using non-blocking sends. */
  vector<SU2_MPI::Request> sendReqs(nRankSend);
  nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (sendToRank[i])
      SU2_MPI::Isend(sendBufsGraphData[i].data(), sendBufsGraphData[i].size(), MPI_UNSIGNED_LONG, i, i,
                     SU2_MPI::GetComm(), &sendReqs[nRankSend++]);
  }

  /* Loop over the number of ranks from which this rank receives data. */
  for (int i = 0; i < nRankRecv; ++i) {
    /* Block until a message with unsigned longs arrives from any processor.
       Determine the source and the size of the message.   */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /* Allocate the memory for the receive buffer and receive the message
       using a blocking receive. */
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, source, rank, SU2_MPI::GetComm(), &status);

    /* Loop over the contents of the receive buffer and update the
       graph accordingly. */
    for (int j = 0; j < sizeMess; j += 2) {
      const unsigned long elemID = recvBuf[j] - beg_node[rank];
      adjacency[elemID].push_back(recvBuf[j + 1]);
      sort(adjacency[elemID].begin(), adjacency[elemID].end());
    }
  }

  /* Complete the non-blocking sends amd synchronize the ranks, because
     wild cards have been used in the above communication. */
  SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif

  /* Allocate the memory to store the weights of the graph. */
  vector<su2double> vwgt(2 * nElem);
  vector<vector<su2double> > adjwgt(nElem, vector<su2double>(0));

  for (unsigned long i = 0; i < nElem; ++i) adjwgt[i].resize(adjacency[i].size());

  /* Compute the weigts of the graph. */
  ComputeFEMGraphWeights(config, localFaces, adjacency, mapExternalElemIDToTimeLevel, vwgt, adjwgt);

  /*--- The remainder of this function should only be called if we have parallel
        support with MPI and have the ParMETIS library compiled and linked. ---*/

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS

  /*--- Only call ParMETIS if we have more than one rank to avoid errors ---*/
  if (size > SINGLE_NODE) {
    /*--- The scalar variables and the options array for the call to ParMETIS. ---*/
    idx_t wgtflag = 3;              // Weights on both the vertices and edges.
    idx_t numflag = 0;              // C-numbering.
    idx_t ncon = 2;                 // Number of constraints.
    real_t ubvec[] = {1.05, 1.05};  // Tolerances for the vertex weights, recommended value is 1.05.
    auto nparts = (idx_t)size;      // Number of subdomains. Must be number of MPI ranks.
    idx_t options[METIS_NOPTIONS];  // Just use the default options.
    METIS_SetDefaultOptions(options);
    options[1] = 0;

    /*--- Determine the array, which stores the distribution of the graph nodes
          over the ranks.     ---*/
    vector<idx_t> vtxdist(size + 1);
    vtxdist[0] = 0;
    for (int i = 0; i < size; ++i) vtxdist[i + 1] = (idx_t)end_node[i];

    /* Create the array xadjPar, which contains the number of edges for each
       vertex of the graph in ParMETIS format. */
    vector<idx_t> xadjPar(nElem + 1);
    xadjPar[0] = 0;
    for (unsigned long i = 0; i < nElem; ++i) xadjPar[i + 1] = xadjPar[i] + (idx_t)adjacency[i].size();

    /* Create the adjacency information in ParMETIS format. */
    vector<idx_t> adjacencyPar(xadjPar[nElem]);
    unsigned long ii = 0;
    for (unsigned long i = 0; i < nElem; ++i) {
      for (unsigned long j = 0; j < adjacency[i].size(); ++j, ++ii) adjacencyPar[ii] = adjacency[i][j];
    }

    /* Create the vertex weights in ParMETIS format. */
    vector<idx_t> vwgtPar(nElem * ncon);
    for (unsigned long i = 0; i < nElem * ncon; ++i) vwgtPar[i] = (idx_t)ceil(vwgt[i]);

    /* Create the adjacency weight in ParMETIS format. */
    vector<idx_t> adjwgtPar(xadjPar[nElem]);
    ii = 0;
    for (unsigned long i = 0; i < nElem; ++i) {
      for (unsigned long j = 0; j < adjwgt[i].size(); ++j, ++ii) adjwgtPar[ii] = (idx_t)ceil(adjwgt[i][j]);
    }

    /* Make sure that an equal distribution is obtained. */
    vector<real_t> tpwgts(size * ncon, 1.0 / ((real_t)size));

    /*--- Calling ParMETIS ---*/
    vector<idx_t> part(nElem);
    if (rank == MASTER_NODE) cout << "Calling ParMETIS...";

    idx_t edgecut;
    MPI_Comm comm = SU2_MPI::GetComm();
    ParMETIS_V3_PartKway(vtxdist.data(), xadjPar.data(), adjacencyPar.data(), vwgtPar.data(), adjwgtPar.data(),
                         &wgtflag, &numflag, &ncon, &nparts, tpwgts.data(), ubvec, options, &edgecut, part.data(),
                         &comm);
    if (rank == MASTER_NODE) {
      cout << " graph partitioning complete (";
      cout << edgecut << " edge cuts)." << endl;
    }

    /*--- Set the color of the elements to the outcome of ParMETIS. ---*/
    for (unsigned long i = 0; i < nElem; ++i) elem[i]->SetColor(part[i]);
  }

#else /* HAVE_PARMETIS */

  if (size > SINGLE_NODE) {
    if (rank == MASTER_NODE) {
      cout << endl;
      cout << "--------------------- WARNING -------------------------------" << endl;
      cout << "SU2 compiled without PARMETIS. A linear distribution is used." << endl;
      cout << "This is very inefficient" << endl;
      cout << "-------------------------------------------------------------" << endl;
      cout << endl;
    }

    /* Set the color to the current rank. */
    for (unsigned long i = 0; i < nElem; ++i) elem[i]->SetColor(rank);
  }

#endif /* HAVE_PARMETIS */

#endif /* HAVE_MPI */
}

void CPhysicalGeometry::DeterminePeriodicFacesFEMGrid(CConfig* config, vector<CFaceOfElement>& localFaces) {
  /*--- Determine a mapping from the global point ID to the local index
        of the points.            ---*/
  map<unsigned long, unsigned long> globalPointIDToLocalInd;
  for (unsigned i = 0; i < nPoint; ++i) {
    globalPointIDToLocalInd[nodes->GetGlobalIndex(i)] = i;
  }

  /*--- Loop over the number of markers present in the grid and check for a periodic one. ---*/
  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); ++iMarker) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      /*--- Determine the donor marker and the transformation from the
            current marker to the donor marker.    ---*/
      unsigned short jMarker = config->GetMarker_Periodic_Donor(config->GetMarker_All_TagBound(iMarker));

      auto center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
      auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
      auto trans = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));

      /*--- Store (center+trans) as it is constant and will be added on. ---*/
      su2double translation[] = {center[0] + trans[0], center[1] + trans[1], center[2] + trans[2]};

      /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
      su2double theta = angles[0];
      su2double phi = angles[1];
      su2double psi = angles[2];

      su2double cosTheta = cos(theta), cosPhi = cos(phi), cosPsi = cos(psi);
      su2double sinTheta = sin(theta), sinPhi = sin(phi), sinPsi = sin(psi);

      /*--- Compute the rotation matrix. Note that the implicit
            ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
      su2double rotMatrix[3][3];
      rotMatrix[0][0] = cosPhi * cosPsi;
      rotMatrix[1][0] = cosPhi * sinPsi;
      rotMatrix[2][0] = -sinPhi;

      rotMatrix[0][1] = sinTheta * sinPhi * cosPsi - cosTheta * sinPsi;
      rotMatrix[1][1] = sinTheta * sinPhi * sinPsi + cosTheta * cosPsi;
      rotMatrix[2][1] = sinTheta * cosPhi;

      rotMatrix[0][2] = cosTheta * sinPhi * cosPsi + sinTheta * sinPsi;
      rotMatrix[1][2] = cosTheta * sinPhi * sinPsi - sinTheta * cosPsi;
      rotMatrix[2][2] = cosTheta * cosPhi;

      /*--- Define the vector to store the faces of the donor. Initialize its
            size to the number of local donor faces.    ---*/
      vector<CMatchingFace> facesDonor(nElem_Bound[jMarker]);

      /*------------------------------------------------------------------*/
      /*--- Step 1: Store the information of the local faces of the donor
                    marker in the variables defined above.             ---*/
      /*------------------------------------------------------------------*/

      /*--- Loop over the local elements of the donor marker. ---*/
      for (unsigned long k = 0; k < nElem_Bound[jMarker]; ++k) {
        /*--- Get the connectivity of this face. The reason for the used
              function arguments for GetCornerPointsAllFaces, is that
              GetCornerPointsAllFaces is an overloaded function.   ---*/
        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long faceConn[6][4];
        bound[jMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        /*--- Search for this face in localFaces. It must be present. ---*/
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for (unsigned short j = 0; j < nPointsPerFace[0]; ++j) thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        /*--- Store the relevant data in facesDonor. ---*/
        facesDonor[k].nDim = nDim;
        facesDonor[k].nCornerPoints = nPointsPerFace[0];
        facesDonor[k].elemID = low->elemID0;
        facesDonor[k].nPoly = low->nPolySol0;
        facesDonor[k].nDOFsElem = low->nDOFsElem0;
        facesDonor[k].elemType = low->elemType0;

        for (unsigned short j = 0; j < nPointsPerFace[0]; ++j) {
          map<unsigned long, unsigned long>::const_iterator MI;
          MI = globalPointIDToLocalInd.find(faceConn[0][j]);
          unsigned long ind = MI->second;

          for (unsigned l = 0; l < nDim; ++l) facesDonor[k].cornerCoor[j][l] = nodes->GetCoord(ind, l);
        }

        /*--- Create the tolerance for this face and sort the coordinates. ---*/
        facesDonor[k].SortFaceCoordinates();
      }

      /*------------------------------------------------------------------*/
      /*--- Step 2: In parallel mode the data of the donor marker is
                    gathered on all ranks.                             ---*/
      /*------------------------------------------------------------------*/

#ifdef HAVE_MPI

      /*--- Check if this is indeed a parallel simulation. ---*/
      if (size > 1) {
        /*--- Allocate the memory for the size arrays in Allgatherv. ---*/
        vector<int> recvCounts(size), displs(size);

        /*--- Create the values of recvCounts for the gather of the facesDonor. ---*/
        int sizeLocal = facesDonor.size();

        SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, SU2_MPI::GetComm());

        /*--- Create the data for the vector displs from the known values of
              recvCounts. Also determine the total size of the data.   ---*/
        displs[0] = 0;
        for (int i = 1; i < size; ++i) displs[i] = displs[i - 1] + recvCounts[i - 1];

        int sizeGlobal = displs.back() + recvCounts.back();

        /*--- SU2_MPI does not support the communication of derived data types,
              at least not at the moment when this was coded. Therefore the data
              from facesDonor is copied into buffers of primitive types, which
              are communicated. ---*/
        vector<unsigned short> shortLocBuf(5 * sizeLocal);
        vector<unsigned long> longLocBuf(sizeLocal);
        vector<su2double> doubleLocBuf(13 * sizeLocal);

        unsigned long cS = 0, cL = 0, cD = 0;
        for (vector<CMatchingFace>::const_iterator fIt = facesDonor.begin(); fIt != facesDonor.end(); ++fIt) {
          shortLocBuf[cS++] = fIt->nCornerPoints;
          shortLocBuf[cS++] = fIt->nDim;
          shortLocBuf[cS++] = fIt->nPoly;
          shortLocBuf[cS++] = fIt->nDOFsElem;
          shortLocBuf[cS++] = fIt->elemType;

          longLocBuf[cL++] = fIt->elemID;

          doubleLocBuf[cD++] = fIt->cornerCoor[0][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[0][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[0][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[1][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[1][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[1][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[2][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[2][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[2][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[3][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[3][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[3][2];

          doubleLocBuf[cD++] = fIt->tolForMatching;
        }

        /*--- Gather the faces from all ranks to all ranks. Use Allgatherv
              for this purpose.                  ---*/
        vector<unsigned short> shortGlobBuf(5 * sizeGlobal);
        vector<unsigned long> longGlobBuf(sizeGlobal);
        vector<su2double> doubleGlobBuf(13 * sizeGlobal);

        SU2_MPI::Allgatherv(longLocBuf.data(), longLocBuf.size(), MPI_UNSIGNED_LONG, longGlobBuf.data(),
                            recvCounts.data(), displs.data(), MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

        for (int i = 0; i < size; ++i) {
          recvCounts[i] *= 5;
          displs[i] *= 5;
        }

        SU2_MPI::Allgatherv(shortLocBuf.data(), shortLocBuf.size(), MPI_UNSIGNED_SHORT, shortGlobBuf.data(),
                            recvCounts.data(), displs.data(), MPI_UNSIGNED_SHORT, SU2_MPI::GetComm());

        for (int i = 0; i < size; ++i) {
          recvCounts[i] /= 5;
          displs[i] /= 5;
          recvCounts[i] *= 13;
          displs[i] *= 13;
        }

        SU2_MPI::Allgatherv(doubleLocBuf.data(), doubleLocBuf.size(), MPI_DOUBLE, doubleGlobBuf.data(),
                            recvCounts.data(), displs.data(), MPI_DOUBLE, SU2_MPI::GetComm());

        /*--- Copy the data back into facesDonor, which will contain the
              global information after the copies. ---*/
        facesDonor.resize(sizeGlobal);
        cS = cL = cD = 0;

        for (auto fIt = facesDonor.begin(); fIt != facesDonor.end(); ++fIt) {
          fIt->nCornerPoints = shortGlobBuf[cS++];
          fIt->nDim = shortGlobBuf[cS++];
          fIt->nPoly = shortGlobBuf[cS++];
          fIt->nDOFsElem = shortGlobBuf[cS++];
          fIt->elemType = shortGlobBuf[cS++];

          fIt->elemID = longGlobBuf[cL++];

          fIt->cornerCoor[0][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[0][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[0][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[1][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[1][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[1][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[2][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[2][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[2][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[3][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[3][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[3][2] = doubleGlobBuf[cD++];

          fIt->tolForMatching = doubleGlobBuf[cD++];
        }
      }
#endif

      /*--- Sort facesDonor in increasing order. ---*/
      sort(facesDonor.begin(), facesDonor.end());

      /*------------------------------------------------------------------*/
      /*--- Step 3: Find for the current marker the required data in the
                    variables for the donor marker.                    ---*/
      /*------------------------------------------------------------------*/

      /*--- Loop over the local faces of this boundary marker. ---*/
      for (unsigned long k = 0; k < nElem_Bound[iMarker]; ++k) {
        /*--- Get the connectivity of this face. The reason for the used
              function arguments for GetCornerPointsAllFaces, is that
              GetCornerPointsAllFaces is an overloaded function.   ---*/
        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long faceConn[6][4];
        bound[iMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        /*--- Search for this face in localFaces. It must be present. ---*/
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for (unsigned short j = 0; j < nPointsPerFace[0]; ++j) thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        /*--- Indicate that this face is a periodic face. This is accomplished
              by setting periodicIndex to iMarker + 1.       ---*/
        low->periodicIndex = iMarker + 1;

        /*--- Store the data for this boundary element also in thisMatchingFace,
              such that a search can be carried out in donorFaces. Note that the
              periodic transformation must be applied to the coordinates. ---*/
        CMatchingFace thisMatchingFace;
        thisMatchingFace.nDim = nDim;
        thisMatchingFace.nCornerPoints = nPointsPerFace[0];
        thisMatchingFace.elemID = low->elemID0;
        thisMatchingFace.nPoly = low->nPolySol0;
        thisMatchingFace.nDOFsElem = low->nDOFsElem0;
        thisMatchingFace.elemType = low->elemType0;

        for (unsigned short j = 0; j < nPointsPerFace[0]; ++j) {
          map<unsigned long, unsigned long>::const_iterator MI;
          MI = globalPointIDToLocalInd.find(faceConn[0][j]);
          unsigned long ind = MI->second;
          const su2double* coor = nodes->GetCoord(ind);

          const su2double dx = coor[0] - center[0];
          const su2double dy = coor[1] - center[1];
          const su2double dz = nDim == 3 ? coor[2] - center[2] : su2double(0.0);

          thisMatchingFace.cornerCoor[j][0] =
              rotMatrix[0][0] * dx + rotMatrix[0][1] * dy + rotMatrix[0][2] * dz + translation[0];
          thisMatchingFace.cornerCoor[j][1] =
              rotMatrix[1][0] * dx + rotMatrix[1][1] * dy + rotMatrix[1][2] * dz + translation[1];
          thisMatchingFace.cornerCoor[j][2] =
              rotMatrix[2][0] * dx + rotMatrix[2][1] * dy + rotMatrix[2][2] * dz + translation[2];
        }

        /*--- Create the tolerance for this face and sort the coordinates. ---*/
        thisMatchingFace.SortFaceCoordinates();

        /*--- Check if thisMatchingFace is present in facesDonor. If so, set the
              missing information in the face corresponding to the iterator low. ---*/
        vector<CMatchingFace>::const_iterator donorLow;
        donorLow = lower_bound(facesDonor.begin(), facesDonor.end(), thisMatchingFace);

        if (donorLow != facesDonor.end()) {
          if (!(thisMatchingFace < *donorLow)) {
            low->elemID1 = donorLow->elemID;
            low->nPolySol1 = donorLow->nPoly;
            low->nDOFsElem1 = donorLow->nDOFsElem;
            low->elemType1 = donorLow->elemType;
            low->periodicIndexDonor = jMarker + 1;
          }
        }
      }
    }
  }
}

void CPhysicalGeometry::DetermineFEMConstantJacobiansAndLenScale(CConfig* config) {
  /* Definition of the object that is used to carry out the BLAS calls. */
  CBlasStructure blasFunctions;

  /*--- Determine a mapping from the global point ID to the local index
        of the points.    ---*/
  map<unsigned long, unsigned long> globalPointIDToLocalInd;
  for (unsigned long i = 0; i < nPoint; ++i) {
    globalPointIDToLocalInd[nodes->GetGlobalIndex(i)] = i;
  }

  /*--- Define the vectors to store the standard elements for the volume elements
        and surface faces. These standard elements will be created based on the
        polynomial degree of the grid.      ---*/
  vector<CFEMStandardElement> standardVolumeElements, standardFaceElements;

  /*--- Loop over the local volume elements. ---*/
  for (unsigned long i = 0; i < nElem; ++i) {
    /*------------------------------------------------------------------------*/
    /*--- Compute the Jacobians of the volume element and determine        ---*/
    /*--- whether the Jacobians can be considered constant.                ---*/
    /*------------------------------------------------------------------------*/

    /*--- Determine the standard element of the volume element. If it does not
          exist yet, it will be created. Note that it suffices to indicate that
          the Jacobians are constant, because the only task here is to determine
          whether or not the Jacobians are constant. ---*/
    unsigned short VTK_Type = elem[i]->GetVTK_Type();
    unsigned short nPolyGrid = elem[i]->GetNPolyGrid();

    unsigned long ii;
    for (ii = 0; ii < standardVolumeElements.size(); ++ii) {
      if (standardVolumeElements[ii].SameStandardElement(VTK_Type, nPolyGrid, true)) break;
    }

    if (ii == standardVolumeElements.size()) standardVolumeElements.emplace_back(VTK_Type, nPolyGrid, true, config);

    /*--- Allocate the memory for some help vectors to carry out the matrix
          product to determine the derivatives of the coordinates w.r.t.
          the parametric coordinates. ---*/
    unsigned short nDOFs = standardVolumeElements[ii].GetNDOFs();
    unsigned short nIntegration = standardVolumeElements[ii].GetNIntegration();

    unsigned long sizeResult = nDim * nDim * nIntegration;
    unsigned long sizeRHS = nDim * nDOFs;

    vector<su2double> vecResult(sizeResult), vecRHS(sizeRHS);

    /* Store the coordinates in vecRHS. */
    unsigned long jj = 0;
    for (unsigned short j = 0; j < nDOFs; ++j) {
      unsigned long nodeID = elem[i]->GetNode(j);

      map<unsigned long, unsigned long>::const_iterator MI = globalPointIDToLocalInd.find(nodeID);
      unsigned long ind = MI->second;
      for (unsigned short k = 0; k < nDim; ++k, ++jj) vecRHS[jj] = nodes->GetCoord(ind, k);
    }

    /*--- Get the pointer to the matrix storage of the basis functions and its
          derivatives. The first nDOFs*nIntegration entries of this matrix
          correspond to the interpolation data to the integration points
          and are not needed. Hence this part is skipped. ---*/
    const su2double* matBasisInt = standardVolumeElements[ii].GetMatBasisFunctionsIntegration();
    const su2double* matDerBasisInt = &matBasisInt[nDOFs * nIntegration];

    /* Carry out the matrix matrix product. The last argument is NULL, such
       that this gemm call is ignored in the profiling. Replace by config if
       if should be included. */
    blasFunctions.gemm(nDim * nIntegration, nDim, nDOFs, matDerBasisInt, vecRHS.data(), vecResult.data(), nullptr);

    /*--- Compute the Jacobians in the integration points and determine
          the minimum and maximum values. Make a distinction between
          a 2D and 3D element. ---*/
    su2double jacMin = 1.e+25, jacMax = -1.e+25;

    switch (nDim) {
      case 2: {
        /* 2D computation. Store the offset between the r and s derivatives. */
        const unsigned int off = 2 * nIntegration;

        /*--- Loop over the integration points to compute the Jacobians. ---*/
        for (unsigned short j = 0; j < nIntegration; ++j) {
          const unsigned int jx = 2 * j;
          const unsigned int jy = jx + 1;
          const su2double dxdr = vecResult[jx], dydr = vecResult[jy];
          const su2double dxds = vecResult[jx + off], dyds = vecResult[jy + off];

          const su2double Jac = dxdr * dyds - dxds * dydr;

          jacMin = min(jacMin, Jac);
          jacMax = max(jacMax, Jac);
        }

        break;
      }

      case 3: {
        /* 3D computation. Store the offset between the r and s and r and t derivatives. */
        const unsigned int offS = 3 * nIntegration, offT = 6 * nIntegration;

        /*--- Loop over the integration points to compute the Jacobians. ---*/
        for (unsigned short j = 0; j < nIntegration; ++j) {
          const unsigned int jx = 3 * j;
          const unsigned int jy = jx + 1, jz = jx + 2;
          const su2double dxdr = vecResult[jx], dydr = vecResult[jy], dzdr = vecResult[jz];
          const su2double dxds = vecResult[jx + offS], dyds = vecResult[jy + offS], dzds = vecResult[jz + offS];
          const su2double dxdt = vecResult[jx + offT], dydt = vecResult[jy + offT], dzdt = vecResult[jz + offT];

          const su2double Jac = dxdr * (dyds * dzdt - dzds * dydt) - dxds * (dydr * dzdt - dzdr * dydt) +
                                dxdt * (dydr * dzds - dzdr * dyds);

          jacMin = min(jacMin, Jac);
          jacMax = max(jacMax, Jac);
        }

        break;
      }
    }

    /*--- Check for negative Jacobians. ---*/
    if (jacMin <= 0.0) SU2_MPI::Error("Negative Jacobian found", CURRENT_FUNCTION);

    /*--- Determine the ratio of the maximum and minimum value of the Jacobian.
          From this ratio, determine whether or not the element is considered to
          have a constant Jacobian. ---*/
    bool constJacobian = (jacMax / jacMin) <= 1.000001;
    elem[i]->SetJacobianConsideredConstant(constJacobian);

    /*------------------------------------------------------------------------*/
    /*--- Compute the normal vectors of the faces of the element and check ---*/
    /*--- for constant Jacobians of the faces.                             ---*/
    /*------------------------------------------------------------------------*/

    /*--- Get the global IDs of the corner points of all the faces
          of this element. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long faceConn[6][4];

    elem[i]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Reset the array, which stores whether or not the faces are
          considered to have a constant Jacobian, to false. ---*/
    elem[i]->ResetJacobianConstantFaces();

    /*--- Loop over the number of faces of this element. ---*/
    su2double jacFaceMax = 0.0;
    for (unsigned short j = 0; j < nFaces; ++j) {
      /*--- Determine the index jj in the standard face elements for this face.
            If not present, create the standard element. ---*/
      switch (nPointsPerFace[j]) {
        case 2:
          VTK_Type = LINE;
          break;
        case 3:
          VTK_Type = TRIANGLE;
          break;
        case 4:
          VTK_Type = QUADRILATERAL;
          break;
      }

      for (jj = 0; jj < standardFaceElements.size(); ++jj) {
        if (standardFaceElements[jj].SameStandardElement(VTK_Type, nPolyGrid, true)) break;
      }

      if (jj == standardFaceElements.size()) standardFaceElements.emplace_back(VTK_Type, nPolyGrid, true, config);

      /*--- Set the pointer to store the face connectivity of this face. ---*/
      unsigned short* connFace = nullptr;
      switch (j) {
        case 0:
          connFace = standardVolumeElements[ii].GetConnFace0();
          break;
        case 1:
          connFace = standardVolumeElements[ii].GetConnFace1();
          break;
        case 2:
          connFace = standardVolumeElements[ii].GetConnFace2();
          break;
        case 3:
          connFace = standardVolumeElements[ii].GetConnFace3();
          break;
        case 4:
          connFace = standardVolumeElements[ii].GetConnFace4();
          break;
        case 5:
          connFace = standardVolumeElements[ii].GetConnFace5();
          break;
      }

      /*--- Store the relevant derivative vectors of the standard element of
            the face as well as the number of DOFs and integration points. ---*/
      nDOFs = standardFaceElements[jj].GetNDOFs();
      nIntegration = standardFaceElements[jj].GetNIntegration();

      const su2double* dr = standardFaceElements[jj].GetDrBasisFunctionsIntegration();
      const su2double* ds = standardFaceElements[jj].GetDsBasisFunctionsIntegration();

      /*--- Allocate the memory for the vector to store the normals in the
            integration points of the face. */
      vector<su2double> normalsFace((nDim + 1) * nIntegration);

      /*--- Compute the unit normals in the integration points. Make a
            distinction between two and three dimensions.  ---*/
      switch (nDim) {
        case 2: {
          /*--- Two dimensional case, for which the faces are edges.
                The normal is the vector normal to the tangent vector
                of the edge. Loop over the integration points. ---*/
          for (unsigned short k = 0; k < nIntegration; ++k) {
            const su2double* drr = &dr[k * nDOFs];
            su2double dxdr = 0.0, dydr = 0.0;
            for (unsigned short l = 0; l < nDOFs; ++l) {
              dxdr += vecRHS[2 * connFace[l]] * drr[l];
              dydr += vecRHS[2 * connFace[l] + 1] * drr[l];
            }

            normalsFace[3 * k] = dydr;
            normalsFace[3 * k + 1] = -dxdr;
            normalsFace[3 * k + 2] = sqrt(dxdr * dxdr + dydr * dydr);
          }

          break;
        }

        case 3: {
          /*--- Three dimensional case, for which the faces are triangles or
                quadrilaterals. The normal is the vector obtained via the
                cross product of two tangent vectors.
                Loop over the integration points.               ---*/
          for (unsigned short k = 0; k < nIntegration; ++k) {
            const su2double *drr = &dr[k * nDOFs], *dss = &ds[k * nDOFs];
            su2double dxdr = 0.0, dydr = 0.0, dzdr = 0.0, dxds = 0.0, dyds = 0.0, dzds = 0.0;

            for (unsigned short l = 0; l < nDOFs; ++l) {
              dxdr += vecRHS[3 * connFace[l]] * drr[l];
              dxds += vecRHS[3 * connFace[l]] * dss[l];
              dydr += vecRHS[3 * connFace[l] + 1] * drr[l];
              dyds += vecRHS[3 * connFace[l] + 1] * dss[l];
              dzdr += vecRHS[3 * connFace[l] + 2] * drr[l];
              dzds += vecRHS[3 * connFace[l] + 2] * dss[l];
            }

            su2double nx = dydr * dzds - dyds * dzdr;
            su2double ny = dxds * dzdr - dxdr * dzds;
            su2double nz = dxdr * dyds - dxds * dydr;

            normalsFace[4 * k] = nx;
            normalsFace[4 * k + 1] = ny;
            normalsFace[4 * k + 2] = nz;
            normalsFace[4 * k + 3] = sqrt(nx * nx + ny * ny + nz * nz);
          }

          break;
        }
      }

      /*--- Double loop over the integration points to determine the minimum
            value of the cosine between the normals. Also the minimum and
            maximum length of the normal is determined.     ---*/
      su2double normLenMin = 1.0e-15, normLenMax = 0.0;
      su2double minCosAngleFaceNormals = 1.0;
      unsigned short kk = (nDim + 1);
      for (unsigned short k = 0; k < nIntegration; ++k) {
        if (k == 0)
          normLenMin = normLenMax = normalsFace[nDim];
        else {
          normLenMin = min(normLenMin, normalsFace[k * kk + nDim]);
          normLenMax = max(normLenMax, normalsFace[k * kk + nDim]);
        }

        for (unsigned short l = k + 1; l < nIntegration; ++l) {
          su2double dot = 0.0;
          for (unsigned short m = 0; m < nDim; ++m) dot += normalsFace[k * kk + m] * normalsFace[l * kk + m];
          dot /= normalsFace[k * kk + nDim] * normalsFace[l * kk + nDim];

          minCosAngleFaceNormals = min(minCosAngleFaceNormals, dot);
        }
      }

      /*--- Compute the ratio of normLenMin and normLenMax and determine
            whether the face is considered to have a constant Jacobian. ---*/
      su2double maxRatioLenFaceNormals = normLenMax / normLenMin;

      constJacobian = minCosAngleFaceNormals >= 0.999999 && maxRatioLenFaceNormals <= 1.000001;

      elem[i]->SetJacobianConstantFace(constJacobian, j);

      /* Update the maximum value of the Jacobian of all faces surrounding
         the current element. */
      jacFaceMax = max(jacFaceMax, normLenMax);
    }

    /*------------------------------------------------------------------------*/
    /*--- Determine the length scale of the element. This is needed to     ---*/
    /*--- compute the computational weights of an element when time        ---*/
    /*--- accurate local time stepping is employed and to determine a      ---*/
    /*--- tolerance when periodic transformations are present. Note that a ---*/
    /*--- factor 2 must be taken into account, which is the length scale   ---*/
    /*--- of all the reference elements used in this code.                 ---*/
    /*------------------------------------------------------------------------*/

    const su2double lenScale = 2.0 * jacMin / jacFaceMax;
    elem[i]->SetLengthScale(lenScale);
  }
}

void CPhysicalGeometry::DetermineDonorElementsWallFunctions(CConfig* config) {
  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Check whether wall functions are used at all.              ---*/
  /*--------------------------------------------------------------------------*/

  bool wallFunctions = false;
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX: {
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) wallFunctions = true;
        break;
      }
      default: /* Just to avoid a compiler warning. */
        break;
    }
  }

  /* If no wall functions are used, nothing needs to be done and a
     return can be made. */
  if (!wallFunctions) return;

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Build the ADT of the linear sub-elements of the locally    ---*/
  /*---         stored volume elements.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /* Determine a mapping from the global point ID to the local index
     of the points. */
  map<unsigned long, unsigned long> globalPointIDToLocalInd;
  for (unsigned long i = 0; i < nPoint; ++i) {
    globalPointIDToLocalInd[nodes->GetGlobalIndex(i)] = i;
  }

  /* Define the vector to store the standard element for the volume elements. */
  vector<CFEMStandardElement> standardVolumeElements;

  /* Define the vectors, which store the mapping from the subelement to the
     parent element, subelement ID within the parent element, the element
     type and the connectivity of the subelements. */
  vector<unsigned long> parentElement;
  vector<unsigned short> subElementIDInParent;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned long> elemConn;

  /* Loop over the local volume elements to create the connectivity of
     the linear sub-elements. */
  for (unsigned long l = 0; l < nElem; ++l) {
    /* Determine the standard element of the volume element. If it does not
       exist yet, it will be created. Note that it suffices to create a
       standard element for a constant Jacobian, because the element is
       split into its linear sub-elements for the ADT. */
    unsigned short VTK_Parent = elem[l]->GetVTK_Type();
    unsigned short nPolyGrid = elem[l]->GetNPolyGrid();

    unsigned long ii;
    for (ii = 0; ii < standardVolumeElements.size(); ++ii) {
      if (standardVolumeElements[ii].SameStandardElement(VTK_Parent, nPolyGrid, true)) break;
    }

    if (ii == standardVolumeElements.size()) standardVolumeElements.emplace_back(VTK_Parent, nPolyGrid, true, config);

    /* Determine the necessary data for splitting the element in its linear
       sub-elements. */
    unsigned short VTK_Type[] = {standardVolumeElements[ii].GetVTK_Type1(), standardVolumeElements[ii].GetVTK_Type2()};
    unsigned short nSubElems[] = {0, 0};
    unsigned short nDOFsPerSubElem[] = {0, 0};
    const unsigned short* connSubElems[] = {nullptr, nullptr};

    if (VTK_Type[0] != NONE) {
      nSubElems[0] = standardVolumeElements[ii].GetNSubElemsType1();
      nDOFsPerSubElem[0] = standardVolumeElements[ii].GetNDOFsPerSubElem(VTK_Type[0]);
      connSubElems[0] = standardVolumeElements[ii].GetSubConnType1();
    }

    if (VTK_Type[1] != NONE) {
      nSubElems[1] = standardVolumeElements[ii].GetNSubElemsType2();
      nDOFsPerSubElem[1] = standardVolumeElements[ii].GetNDOFsPerSubElem(VTK_Type[1]);
      connSubElems[1] = standardVolumeElements[ii].GetSubConnType2();
    }

    /* Store the connectivity of the sub-elements. Note that local node
       numbering must be used for these sub-elements. */
    unsigned short jj = 0;
    for (unsigned short i = 0; i < 2; ++i) {
      unsigned short kk = 0;
      for (unsigned short j = 0; j < nSubElems[i]; ++j, ++jj) {
        parentElement.push_back(elem[l]->GetGlobalElemID());
        subElementIDInParent.push_back(jj);
        VTK_TypeElem.push_back(VTK_Type[i]);

        for (unsigned short k = 0; k < nDOFsPerSubElem[i]; ++k, ++kk) {
          unsigned long nodeID = elem[l]->GetNode(connSubElems[i][kk]);
          map<unsigned long, unsigned long>::const_iterator MI;
          MI = globalPointIDToLocalInd.find(nodeID);

          elemConn.push_back(MI->second);
        }
      }
    }
  }

  /* Store the coordinates of the locally stored nodes in the format
     expected by the ADT. */
  vector<su2double> volCoor(nDim * nPoint);

  unsigned long jj = 0;
  for (unsigned long l = 0; l < nPoint; ++l) {
    for (unsigned short k = 0; k < nDim; ++k, ++jj) volCoor[jj] = nodes->GetCoord(l, k);
  }

  /* Build the local ADT. */
  CADTElemClass localVolumeADT(nDim, volCoor, elemConn, VTK_TypeElem, subElementIDInParent, parentElement, false);

  /* Release the memory of the vectors used to build the ADT. To make sure
     that all the memory is deleted, the swap function is used. */
  vector<unsigned short>().swap(subElementIDInParent);
  vector<unsigned short>().swap(VTK_TypeElem);
  vector<unsigned long>().swap(parentElement);
  vector<unsigned long>().swap(elemConn);
  vector<su2double>().swap(volCoor);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3. Search for donor elements at the exchange locations in     ---*/
  /*---         the local elements.                                        ---*/
  /*--------------------------------------------------------------------------*/

  /* Define the vectors, which store the boundary marker, boundary element ID
     and exchange coordinates for the integration points for which no donor
     element was found in the locally stored volume elements. */
  vector<unsigned short> markerIDGlobalSearch;
  vector<unsigned long> boundaryElemIDGlobalSearch;
  vector<su2double> coorExGlobalSearch;

  /* Define the standard boundary faces for the solution and the grid. */
  vector<CFEMStandardBoundaryFace> standardBoundaryFacesSol, standardBoundaryFacesGrid;

  /* Loop over the markers and select the ones for which a wall function
     treatment must be carried out. */
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX: {
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) {
          /* Retrieve the floating point information for this boundary marker.
             The exchange location is the first element of this array. */
          const su2double* doubleInfo = config->GetWallFunction_DoubleInfo(Marker_Tag);

          /* Loop over the local boundary elements for this marker. */
          for (unsigned long l = 0; l < nElem_Bound[iMarker]; ++l) {
            /* Easier storage of the element type, the corresponding volume
               element and the polynomial degree for the solution and grid. */
            const unsigned short VTK_Type = bound[iMarker][l]->GetVTK_Type();
            const unsigned long elemID = bound[iMarker][l]->GetDomainElement() - beg_node[rank];
            const unsigned short nPolyGrid = bound[iMarker][l]->GetNPolyGrid();
            const unsigned short nPolySol = elem[elemID]->GetNPolySol();
            const unsigned short VTK_Elem = elem[elemID]->GetVTK_Type();

            /* Get the corner points of the boundary element. Note that this
               is an overloaded function, hence the arguments that allow for
               multiple faces. */
            unsigned short nFaces;
            unsigned short nPointsPerFace[6];
            unsigned long faceConn[6][4];
            bound[iMarker][l]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

            /* Create an object of CFaceOfElement to store the information. */
            CFaceOfElement boundaryFace;
            boundaryFace.nCornerPoints = nPointsPerFace[0];
            for (unsigned short i = 0; i < nPointsPerFace[0]; ++i) boundaryFace.cornerPoints[i] = faceConn[0][i];
            boundaryFace.elemID0 = elemID;

            /* Renumber the corner points of the face, but keep the orientation. */
            boundaryFace.CreateUniqueNumberingWithOrientation();
            const bool boundaryFaceSwapped = boundaryFace.elemID1 == elemID;

            /* Get the corner points of all the faces from the corresponding
               volume element. */
            elem[elemID]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

            /* Loop over the faces of the element to determine which one is
               the boundary face. */
            unsigned long ii;
            bool outwardPointing = true; /* Initialized to avoid compiler warning. */
            for (ii = 0; ii < nFaces; ++ii) {
              /* Create an object of CFaceOfElement to store the information.
                 Renumber the corner points, but keep the orientation. */
              CFaceOfElement thisFace;
              thisFace.nCornerPoints = nPointsPerFace[ii];
              for (unsigned short i = 0; i < nPointsPerFace[ii]; ++i) thisFace.cornerPoints[i] = faceConn[ii][i];
              thisFace.elemID0 = elemID;

              thisFace.CreateUniqueNumberingWithOrientation();
              const bool thisFaceSwapped = thisFace.elemID1 == elemID;

              /* Check if this face is the boundary face. */
              if (boundaryFace == thisFace) {
                /* Determine whether the orientation of the boundary face is
                   pointing out of the adjacent element. */
                outwardPointing = boundaryFaceSwapped == thisFaceSwapped;

                /* Break the loop over the faces of the element. */
                break;
              }
            }

            /* Additional check, just to be sure. */
            if (ii == nFaces) SU2_MPI::Error("Boundary face not found in faces of element", CURRENT_FUNCTION);

            /* Determine whether or not the Jacobian of the boundary face
               is constant. */
            const bool constJac = elem[elemID]->GetJacobianConstantFace(ii);
            bound[iMarker][l]->SetJacobianConsideredConstant(constJac);

            /* Determine the corresponding standard element. If it does not
               exist yet, it will be created. Note that the polynomial degree
               of both the solution and grid must match. If a new standard
               element must be created, make sure that the integration points
               of the solution element are used for the grid element. */
            for (ii = 0; ii < standardBoundaryFacesSol.size(); ++ii) {
              if (standardBoundaryFacesSol[ii].SameStandardBoundaryFace(VTK_Type, constJac, VTK_Elem, nPolySol,
                                                                        false) &&
                  standardBoundaryFacesGrid[ii].SameStandardBoundaryFace(VTK_Type, constJac, VTK_Elem, nPolyGrid,
                                                                         false))
                break;
            }

            if (ii == standardBoundaryFacesSol.size()) {
              standardBoundaryFacesSol.emplace_back(VTK_Type, VTK_Elem, nPolySol, constJac, false, config);
              standardBoundaryFacesGrid.emplace_back(VTK_Type, VTK_Elem, nPolyGrid, constJac, false, config,
                                                     standardBoundaryFacesSol[ii].GetOrderExact());
            }

            /* Get the required information from the standard element. */
            const unsigned short nDOFs = standardBoundaryFacesGrid[ii].GetNDOFsFace();
            const unsigned short nInt = standardBoundaryFacesGrid[ii].GetNIntegration();
            const su2double* lag = standardBoundaryFacesGrid[ii].GetBasisFaceIntegration();
            const su2double* dr = standardBoundaryFacesGrid[ii].GetDrBasisFaceIntegration();
            const su2double* ds = standardBoundaryFacesGrid[ii].GetDsBasisFaceIntegration();

            /* Create the vector of coordinates for this boundary face. */
            vector<su2double> coorBoundFace(nDim * nDOFs);
            ii = 0;
            for (unsigned short j = 0; j < nDOFs; ++j) {
              unsigned long nodeID = bound[iMarker][l]->GetNode(j);
              map<unsigned long, unsigned long>::const_iterator MI;
              MI = globalPointIDToLocalInd.find(nodeID);
              nodeID = MI->second;
              for (unsigned short k = 0; k < nDim; ++k, ++ii) coorBoundFace[ii] = nodes->GetCoord(nodeID, k);
            }

            /* Set the multiplication factor for the normal, such that
               an inward pointing normal is obtained. It is scaled with
               the exchange distance. */
            const su2double factNorm = outwardPointing ? -doubleInfo[0] : doubleInfo[0];

            /* Allocate the memory for the exchange locations corresponding to
               the integration points of this face. */
            vector<su2double> coorExchange(nDim * nInt);

            /* Make a distinction between 2D and 3D to compute the actual
               coordinates of the exchange locations. */
            switch (nDim) {
              case 2: {
                /* Two dimensional computation. Loop over the integration
                   points to compute the corresponding exchange coordinates. */
                for (unsigned short i = 0; i < nInt; ++i) {
                  /* Compute the coordinates and the tangential vector
                     in the integration point. */
                  su2double xInt = 0.0, yInt = 0.0, dxdr = 0.0, dydr = 0.0;
                  const su2double* lagInt = lag + i * nDOFs;
                  const su2double* drr = dr + i * nDOFs;

                  for (unsigned short k = 0; k < nDOFs; ++k) {
                    const su2double* coorDOF = coorBoundFace.data() + 2 * k;
                    xInt += lagInt[k] * coorDOF[0];
                    yInt += lagInt[k] * coorDOF[1];
                    dxdr += drr[k] * coorDOF[0];
                    dydr += drr[k] * coorDOF[1];
                  }

                  /* Compute the vector from the integration point to
                     the coordinates of the exchange location. This is the
                     unit inward normal multiplied by the exchange distance. */
                  const su2double lenNorm = sqrt(dxdr * dxdr + dydr * dydr);
                  const su2double invLenNorm = lenNorm < su2double(1.e-35) ? su2double(1.e+35) : 1.0 / lenNorm;

                  const su2double dx = factNorm * dydr * invLenNorm;
                  const su2double dy = -factNorm * dxdr * invLenNorm;

                  /* Compute the coordinates of the exchange location. */
                  su2double* coorEx = coorExchange.data() + nDim * i;
                  coorEx[0] = xInt + dx;
                  coorEx[1] = yInt + dy;
                }

                break;
              }

              case 3: {
                /* Three dimensional computation. Loop over the integration
                   points to compute the corresponding exchange coordinates. */
                for (unsigned short i = 0; i < nInt; ++i) {
                  /* Compute the coordinates and the tangential vectors
                     in the integration point. */
                  const su2double* lagInt = lag + i * nDOFs;
                  const su2double *drr = dr + i * nDOFs, *dss = ds + i * nDOFs;
                  su2double xInt = 0.0, yInt = 0.0, zInt = 0.0, dxdr = 0.0, dydr = 0.0, dzdr = 0.0, dxds = 0.0,
                            dyds = 0.0, dzds = 0.0;

                  for (unsigned short k = 0; k < nDOFs; ++k) {
                    const su2double* coorDOF = coorBoundFace.data() + 3 * k;
                    xInt += lagInt[k] * coorDOF[0];
                    yInt += lagInt[k] * coorDOF[1];
                    zInt += lagInt[k] * coorDOF[2];
                    dxdr += drr[k] * coorDOF[0];
                    dydr += drr[k] * coorDOF[1];
                    dzdr += drr[k] * coorDOF[2];
                    dxds += dss[k] * coorDOF[0];
                    dyds += dss[k] * coorDOF[1];
                    dzds += dss[k] * coorDOF[2];
                  }

                  /* Compute the outward pointing normal. */
                  const su2double nx = -(dydr * dzds - dyds * dzdr);
                  const su2double ny = -(dxds * dzdr - dxdr * dzds);
                  const su2double nz = -(dxdr * dyds - dxds * dydr);

                  /* Compute the vector from the integration point to
                     the coordinates of the exchange location. This is the
                     unit inward normal multiplied by the exchange distance. */
                  const su2double lenNorm = sqrt(nx * nx + ny * ny + nz * nz);
                  const su2double invLenNorm = lenNorm < su2double(1.e-35) ? su2double(1.e+35) : 1.0 / lenNorm;

                  const su2double dx = factNorm * nx * invLenNorm;
                  const su2double dy = factNorm * ny * invLenNorm;
                  const su2double dz = factNorm * nz * invLenNorm;

                  /* Compute the coordinates of the exchange location. */
                  su2double* coorEx = coorExchange.data() + nDim * i;
                  coorEx[0] = xInt + dx;
                  coorEx[1] = yInt + dy;
                  coorEx[2] = zInt + dz;
                }

                break;
              }
            }

            /* Loop over the integration points and determine the donor elements
               for the corresponding exchange locations. */
            vector<unsigned long> donorElementsFace;
            for (unsigned short i = 0; i < nInt; ++i) {
              /* Search for the element, which contains the exchange location. */
              unsigned short subElem;
              unsigned long parElem;
              int mpirank;
              su2double parCoor[3], weightsInterpol[8];
              su2double* coorEx = coorExchange.data() + nDim * i;

              if (localVolumeADT.DetermineContainingElement(coorEx, subElem, parElem, mpirank, parCoor,
                                                            weightsInterpol)) {
                /* Donor element found. Store it in donorElementsFace. */
                donorElementsFace.push_back(parElem);
              } else {
                /* Donor element not found in the local ADT. Store the exchange
                   coordinates, boundary marker and boundary element ID. */
                markerIDGlobalSearch.push_back(iMarker);
                boundaryElemIDGlobalSearch.push_back(l);
                for (unsigned short iDim = 0; iDim < nDim; ++iDim) coorExGlobalSearch.push_back(coorEx[iDim]);
              }
            }

            /* Sort donorElementsFace in increasing order and remove the
               the double entities. */
            sort(donorElementsFace.begin(), donorElementsFace.end());
            vector<unsigned long>::iterator lastEntry;
            lastEntry = unique(donorElementsFace.begin(), donorElementsFace.end());
            donorElementsFace.erase(lastEntry, donorElementsFace.end());

            /* Store the donor elements in the data structure for
               this boundary element. */
            bound[iMarker][l]->SetDonorsWallFunctions(donorElementsFace);
          }
        }

        break;
      }
      default: /* Just to avoid a compiler warning. */
        break;
    }
  }

#ifndef HAVE_MPI
  /* In sequential mode it should not be possible that there are points for which
     no donor elements are found. Check this. */
  if (markerIDGlobalSearch.size()) {
    ostringstream message;
    message << markerIDGlobalSearch.size() << " integration points for which "
            << "no donor elements for the wall functions were found." << endl
            << "Something is seriously wrong";
    SU2_MPI::Error("message.str()", CURRENT_FUNCTION);
  }
#endif

  /* The remaining part of this function only needs to be carried out in parallel mode. */
#ifdef HAVE_MPI

  /* Determine the number of search points for which a global search must be
     carried out for each rank and store them in such a way that the info can
     be used directly in Allgatherv. */
  vector<int> recvCounts(size), displs(size);
  int nLocalSearchPoints = (int)markerIDGlobalSearch.size();

  SU2_MPI::Allgather(&nLocalSearchPoints, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, SU2_MPI::GetComm());
  displs[0] = 0;
  for (int i = 1; i < size; ++i) displs[i] = displs[i - 1] + recvCounts[i - 1];

  int nGlobalSearchPoints = displs.back() + recvCounts.back();

  /* Check if there actually are global searches to be carried out. */
  if (nGlobalSearchPoints > 0) {
    /* Create a cumulative storage version of recvCounts. */
    vector<int> nSearchPerRank(size + 1);
    nSearchPerRank[0] = 0;

    for (int i = 0; i < size; ++i) nSearchPerRank[i + 1] = nSearchPerRank[i] + recvCounts[i];

    /* Gather the data of the search points for which a global search must
       be carried out on all ranks. */
    vector<unsigned short> bufMarkerIDGlobalSearch(nGlobalSearchPoints);
    SU2_MPI::Allgatherv(markerIDGlobalSearch.data(), nLocalSearchPoints, MPI_UNSIGNED_SHORT,
                        bufMarkerIDGlobalSearch.data(), recvCounts.data(), displs.data(), MPI_UNSIGNED_SHORT,
                        SU2_MPI::GetComm());

    vector<unsigned long> bufBoundaryElemIDGlobalSearch(nGlobalSearchPoints);
    SU2_MPI::Allgatherv(boundaryElemIDGlobalSearch.data(), nLocalSearchPoints, MPI_UNSIGNED_LONG,
                        bufBoundaryElemIDGlobalSearch.data(), recvCounts.data(), displs.data(), MPI_UNSIGNED_LONG,
                        SU2_MPI::GetComm());

    for (int i = 0; i < size; ++i) {
      recvCounts[i] *= nDim;
      displs[i] *= nDim;
    }
    vector<su2double> bufCoorExGlobalSearch(nDim * nGlobalSearchPoints);
    SU2_MPI::Allgatherv(coorExGlobalSearch.data(), nDim * nLocalSearchPoints, MPI_DOUBLE, bufCoorExGlobalSearch.data(),
                        recvCounts.data(), displs.data(), MPI_DOUBLE, SU2_MPI::GetComm());

    /* Buffers to store the return information. */
    vector<unsigned short> markerIDReturn;
    vector<unsigned long> boundaryElemIDReturn;
    vector<unsigned long> volElemIDDonorReturn;

    /* Loop over the number of global search points to check if these points
       are contained in the volume elements of this rank. The loop is carried
       out as a double loop, such that the rank where the point resides is
       known as well. Furthermore, it is not necessary to search the points
       that were not found earlier on this rank. The vector recvCounts is used
       as storage for the number of search items that must be returned to the
       other ranks. */
    for (int rankID = 0; rankID < size; ++rankID) {
      recvCounts[rankID] = 0;
      if (rankID != rank) {
        for (int i = nSearchPerRank[rankID]; i < nSearchPerRank[rankID + 1]; ++i) {
          /* Search the local ADT for the coordinate of the exchange point
             and check if it is found. */
          unsigned short subElem;
          unsigned long parElem;
          int rankDonor;
          su2double parCoor[3], weightsInterpol[8];
          if (localVolumeADT.DetermineContainingElement(bufCoorExGlobalSearch.data() + i * nDim, subElem, parElem,
                                                        rankDonor, parCoor, weightsInterpol)) {
            /* Store the required data in the return buffers. */
            ++recvCounts[rankID];
            markerIDReturn.push_back(bufMarkerIDGlobalSearch[i]);
            boundaryElemIDReturn.push_back(bufBoundaryElemIDGlobalSearch[i]);
            volElemIDDonorReturn.push_back(parElem);
          }
        }
      }
    }

    /* Create a cumulative version of recvCounts. */
    for (int i = 0; i < size; ++i) nSearchPerRank[i + 1] = nSearchPerRank[i] + recvCounts[i];

    /* Determine the number of return messages this rank has to receive.
       Use displs and recvCounts as temporary storage. */
    int nRankSend = 0;
    for (int i = 0; i < size; ++i) {
      if (recvCounts[i]) {
        recvCounts[i] = 1;
        ++nRankSend;
      }
      displs[i] = 1;
    }

    int nRankRecv;
    SU2_MPI::Reduce_scatter(recvCounts.data(), &nRankRecv, displs.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

    /* Send the data using nonblocking sends to avoid deadlock. */
    vector<SU2_MPI::Request> commReqs(3 * nRankSend);
    nRankSend = 0;
    for (int i = 0; i < size; ++i) {
      if (recvCounts[i]) {
        const int sizeMessage = nSearchPerRank[i + 1] - nSearchPerRank[i];
        SU2_MPI::Isend(markerIDReturn.data() + nSearchPerRank[i], sizeMessage, MPI_UNSIGNED_SHORT, i, i,
                       SU2_MPI::GetComm(), &commReqs[nRankSend++]);
        SU2_MPI::Isend(boundaryElemIDReturn.data() + nSearchPerRank[i], sizeMessage, MPI_UNSIGNED_LONG, i, i + 1,
                       SU2_MPI::GetComm(), &commReqs[nRankSend++]);
        SU2_MPI::Isend(volElemIDDonorReturn.data() + nSearchPerRank[i], sizeMessage, MPI_UNSIGNED_LONG, i, i + 2,
                       SU2_MPI::GetComm(), &commReqs[nRankSend++]);
      }
    }

    /* Loop over the number of ranks from which I receive return data. */
    for (int i = 0; i < nRankRecv; ++i) {
      /* Block until a message with unsigned shorts arrives from any processor.
         Determine the source and the size of the message.   */
      SU2_MPI::Status status;
      SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
      int source = status.MPI_SOURCE;

      int sizeMess;
      SU2_MPI::Get_count(&status, MPI_UNSIGNED_SHORT, &sizeMess);

      /* Allocate the memory for the receive buffers. */
      vector<unsigned short> bufMarkerIDReturn(sizeMess);
      vector<unsigned long> bufBoundaryElemIDReturn(sizeMess);
      vector<unsigned long> bufVolElemIDDonorReturn(sizeMess);

      /* Receive the three messages using blocking receives. */
      SU2_MPI::Recv(bufMarkerIDReturn.data(), sizeMess, MPI_UNSIGNED_SHORT, source, rank, SU2_MPI::GetComm(), &status);

      SU2_MPI::Recv(bufBoundaryElemIDReturn.data(), sizeMess, MPI_UNSIGNED_LONG, source, rank + 1, SU2_MPI::GetComm(),
                    &status);

      SU2_MPI::Recv(bufVolElemIDDonorReturn.data(), sizeMess, MPI_UNSIGNED_LONG, source, rank + 2, SU2_MPI::GetComm(),
                    &status);

      /* Loop over the data just received and add it to the wall function
         donor information of the corresponding boundary element. */
      for (int j = 0; j < sizeMess; ++j) {
        const unsigned short iMarker = bufMarkerIDReturn[j];
        const unsigned long l = bufBoundaryElemIDReturn[j];
        const unsigned long volID = bufVolElemIDDonorReturn[j];

        bound[iMarker][l]->AddDonorWallFunctions(volID);
      }
    }

    /* Complete the non-blocking sends. */
    SU2_MPI::Waitall(nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

    /* Wild cards have been used in the communication,
       so synchronize the ranks to avoid problems. */
    SU2_MPI::Barrier(SU2_MPI::GetComm());

    /* Loop again over the boundary elements of the marker for which a wall
       function treatment must be used and make remove the multiple entries
       of the donor information. */
    for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case ISOTHERMAL:
        case HEAT_FLUX: {
          const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) {
            for (unsigned long l = 0; l < nElem_Bound[iMarker]; ++l)
              bound[iMarker][l]->RemoveMultipleDonorsWallFunctions();
          }

          break;
        }

        default: /* Just to avoid a compiler warning. */
          break;
      }
    }
  }

#endif
}

void CPhysicalGeometry::DetermineTimeLevelElements(CConfig* config, const vector<CFaceOfElement>& localFaces,
                                                   map<unsigned long, CUnsignedShort2T>& mapExternalElemIDToTimeLevel) {
  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Initialize the map mapExternalElemIDToTimeLevel.           ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Initialize the time level of external elements to zero.
        First the externals from the faces. ---*/
  for (auto FI = localFaces.begin(); FI != localFaces.end(); ++FI) {
    if (FI->elemID1 < Global_nElem) {  // Safeguard against non-matching faces.

      if (FI->elemID1 < beg_node[rank] || FI->elemID1 >= beg_node[rank] + nElem) {
        /* This element is an external element. Store it in the map
           mapExternalElemIDToTimeLevel if not already done so. */
        map<unsigned long, CUnsignedShort2T>::iterator MI;
        MI = mapExternalElemIDToTimeLevel.find(FI->elemID1);
        if (MI == mapExternalElemIDToTimeLevel.end())
          mapExternalElemIDToTimeLevel[FI->elemID1] = CUnsignedShort2T(0, 0);
      }
    }
  }

  /* Define the communication buffers to send additional externals
     to other ranks. Also define the buffer recvFromRank which is used
     in Reduce_scatter later on. It indicates whether or not a message
     is sent to a certain rank. */
  vector<vector<unsigned long> > sendBufAddExternals(size, vector<unsigned long>(0));
  vector<int> recvFromRank(size, 0);

  /*--- Add the externals from the wall function donors. Loop over the
        boundary elements of all markers. ---*/
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    for (unsigned long l = 0; l < nElem_Bound[iMarker]; ++l) {
      /* Get the number of donor elements for the wall function treatment
         and the pointer to the array which stores this info. */
      const unsigned short nDonors = bound[iMarker][l]->GetNDonorsWallFunctions();
      const unsigned long* donors = bound[iMarker][l]->GetDonorsWallFunctions();

      /* Loop over the number of donors and add the externals to
         mapExternalElemIDToTimeLevel, if not already present. */
      for (unsigned short i = 0; i < nDonors; ++i) {
        if (donors[i] < beg_node[rank] || donors[i] >= beg_node[rank] + nElem) {
          map<unsigned long, CUnsignedShort2T>::iterator MI;
          MI = mapExternalElemIDToTimeLevel.find(donors[i]);
          if (MI == mapExternalElemIDToTimeLevel.end()) {
            /* Element not present in external. Add it. */
            mapExternalElemIDToTimeLevel[donors[i]] = CUnsignedShort2T(0, 0);
          }

          /* The reverse connection may not be present either. Store the global
             ID of this element in the send buffers for the additional
             externals. */
          int rankDonor;
          if (donors[i] >= beg_node[size - 1])
            rankDonor = size - 1;
          else {
            const unsigned long* low;
            low = lower_bound(beg_node, beg_node + size, donors[i]);

            rankDonor = (int)(low - beg_node);
            if (*low > donors[i]) --rankDonor;
          }

          sendBufAddExternals[rankDonor].push_back(bound[iMarker][l]->GetDomainElement());
          recvFromRank[rankDonor] = 1;
        }
      }
    }
  }

#ifdef HAVE_MPI

  /* Determine the number of messages this rank will receive with additional
     externals to be stored. */
  int nRankRecv;
  vector<int> sizeSend(size, 1);
  SU2_MPI::Reduce_scatter(recvFromRank.data(), &nRankRecv, sizeSend.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /* Determine the number of messages this rank will send. */
  int nRankSend = 0;
  for (int i = 0; i < size; ++i) nRankSend += recvFromRank[i];

  /* Send the data using non-blocking sends to avoid deadlock. */
  vector<SU2_MPI::Request> sendReqs(nRankSend);
  nRankSend = 0;
  for (int i = 0; i < size; ++i) {
    if (recvFromRank[i]) {
      sort(sendBufAddExternals[i].begin(), sendBufAddExternals[i].end());
      auto lastElem = unique(sendBufAddExternals[i].begin(), sendBufAddExternals[i].end());
      sendBufAddExternals[i].erase(lastElem, sendBufAddExternals[i].end());

      SU2_MPI::Isend(sendBufAddExternals[i].data(), sendBufAddExternals[i].size(), MPI_UNSIGNED_LONG, i, i,
                     SU2_MPI::GetComm(), &sendReqs[nRankSend++]);
    }
  }

  /* Loop over the number of ranks from which this rank will receive data
     to be stored in mapExternalElemIDToTimeLevel. */
  for (int i = 0; i < nRankRecv; ++i) {
    /* Block until a message arrives and determine the source and size
       of the message. Allocate the memory for a receive buffer. */
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);
    vector<unsigned long> recvBuf(sizeMess);

    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG, source, rank, SU2_MPI::GetComm(), &status);

    /* Loop over the entries of recvBuf and add them to
       mapExternalElemIDToTimeLevel, if not present already. */
    for (int j = 0; j < sizeMess; ++j) {
      map<unsigned long, CUnsignedShort2T>::iterator MI;
      MI = mapExternalElemIDToTimeLevel.find(recvBuf[j]);
      if (MI == mapExternalElemIDToTimeLevel.end()) mapExternalElemIDToTimeLevel[recvBuf[j]] = CUnsignedShort2T(0, 0);
    }
  }

  /* Complete the non-blocking sends. Synchronize the processors afterwards,
     because wild cards have been used in the communication. */
  SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Initialize the time level of the owned elements.           ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Get the number of time levels used and check whether or not
        time accurate local time stepping is used. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

  if (nTimeLevels == 1) {
    /*--- No time accurate local time stepping. Set the time level of
          all elements to zero. ---*/
    for (unsigned long i = 0; i < nElem; ++i) elem[i]->SetTimeLevel(0);
  } else {
    /*--- Time accurate local time stepping is used. The estimate of the time
          step is based on free stream values at the moment, but this is easy
          to change, if needed. ---*/
    const su2double Gamma = config->GetGamma();
    const su2double Prandtl = config->GetPrandtl_Lam();

    const su2double Density = config->GetDensity_FreeStreamND();
    const su2double* Vel = config->GetVelocity_FreeStreamND();
    const su2double Viscosity = config->GetViscosity_FreeStreamND();

    su2double VelMag = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) VelMag += Vel[iDim] * Vel[iDim];
    VelMag = sqrt(VelMag);

    const su2double SoundSpeed = VelMag / config->GetMach();

    /* In the estimate of the time step the spectral radius of the inviscid terms is
       needed. As the current estimate is based on the free stream, the value of this
       spectral radius can be computed beforehand. Note that this is a rather
       conservative estimate. */
    su2double charVel2 = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      const su2double rad = fabs(Vel[iDim]) + SoundSpeed;
      charVel2 += rad * rad;
    }

    const su2double charVel = sqrt(charVel2);

    /* Also the viscous contribution to the time step is constant. Compute it. */
    const su2double factHeatFlux = Gamma / Prandtl;
    const su2double lambdaOverMu = -TWO3;
    const su2double radVisc = max(max(1.0, 2.0 + lambdaOverMu), factHeatFlux) * Viscosity / Density;

    /*--- Allocate the memory for time step estimate of the local elements and
          determine the values. Also keep track of the minimum value. ---*/
    su2double minDeltaT = 1.e25;
    vector<su2double> timeStepElements(nElem);

    for (unsigned long i = 0; i < nElem; ++i) {
      unsigned short nPoly = elem[i]->GetNPolySol();
      if (nPoly == 0) nPoly = 1;
      const su2double lenScaleInv = nPoly / elem[i]->GetLengthScale();
      const su2double lenScale = 1.0 / lenScaleInv;

      timeStepElements[i] = lenScale / (charVel + lenScaleInv * radVisc);
      minDeltaT = min(minDeltaT, timeStepElements[i]);
    }

    /* Determine the minimum value of all elements in the grid.
       Only needed for a parallel implementation. */
#ifdef HAVE_MPI
    su2double locVal = minDeltaT;
    SU2_MPI::Allreduce(&locVal, &minDeltaT, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
#endif

    /* Initial estimate of the time level of the owned elements. */
    for (unsigned long i = 0; i < nElem; ++i) {
      unsigned short timeLevel;
      su2double deltaT = minDeltaT;
      for (timeLevel = 0; timeLevel < (nTimeLevels - 1); ++timeLevel) {
        deltaT *= 2;
        if (timeStepElements[i] < deltaT) break;
      }

      elem[i]->SetTimeLevel(timeLevel);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Set up the variables to carry out the MPI communication    ---*/
  /*---         of the external element data.                              ---*/
  /*--------------------------------------------------------------------------*/

  map<unsigned long, CUnsignedShort2T>::iterator MI;
#ifdef HAVE_MPI

  /*--- Determine the ranks from which I receive element data during
        the actual exchange. ---*/
  recvFromRank.assign(size, 0);

  for (MI = mapExternalElemIDToTimeLevel.begin(); MI != mapExternalElemIDToTimeLevel.end(); ++MI) {
    /* Determine the rank where this external is stored. */
    const unsigned long elemID = MI->first;
    int rankElem;
    if (elemID >= beg_node[size - 1])
      rankElem = size - 1;
    else {
      const unsigned long* low;
      low = lower_bound(beg_node, beg_node + size, elemID);

      rankElem = low - beg_node;
      if (*low > elemID) --rankElem;
    }

    /* Set the corresponding index of recvFromRank to 1. */
    recvFromRank[rankElem] = 1;
  }

  map<int, int> mapRankToIndRecv;
  for (int i = 0; i < size; ++i) {
    if (recvFromRank[i]) {
      int ind = mapRankToIndRecv.size();
      mapRankToIndRecv[i] = ind;
    }
  }

  /* Determine the number of ranks from which I will receive data and to
     which I will send data. */
  nRankRecv = mapRankToIndRecv.size();
  SU2_MPI::Reduce_scatter(recvFromRank.data(), &nRankSend, sizeSend.data(), MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /*--- Create the vector of vectors of the global element ID's that
        will be received from other ranks. ---*/
  vector<vector<unsigned long> > recvElem(nRankRecv, vector<unsigned long>(0));

  for (MI = mapExternalElemIDToTimeLevel.begin(); MI != mapExternalElemIDToTimeLevel.end(); ++MI) {
    const unsigned long elemID = MI->first;
    int rankElem;
    if (elemID >= beg_node[size - 1])
      rankElem = size - 1;
    else {
      const unsigned long* low;
      low = lower_bound(beg_node, beg_node + size, elemID);

      rankElem = low - beg_node;
      if (*low > elemID) --rankElem;
    }

    map<int, int>::const_iterator MRI = mapRankToIndRecv.find(rankElem);
    recvElem[MRI->second].push_back(elemID);
  }

  /*--- Loop over the ranks from which I receive data during the actual
        exchange and send over the global element ID's. In order to avoid
        unnecessary communication, multiple entries are filtered out. ---*/
  sendReqs.resize(nRankRecv);
  map<int, int>::const_iterator MRI = mapRankToIndRecv.begin();

  for (int i = 0; i < nRankRecv; ++i, ++MRI) {
    sort(recvElem[i].begin(), recvElem[i].end());
    auto lastElem = unique(recvElem[i].begin(), recvElem[i].end());
    recvElem[i].erase(lastElem, recvElem[i].end());

    SU2_MPI::Isend(recvElem[i].data(), recvElem[i].size(), MPI_UNSIGNED_LONG, MRI->first, MRI->first,
                   SU2_MPI::GetComm(), &sendReqs[i]);
  }

  /*--- Receive the messages in arbitrary sequence and store the requested
        element ID's, which are converted to local ID's. Furthermore, store
        the processors from which the requested element ID's came from. ---*/
  vector<vector<unsigned long> > sendElem(nRankSend, vector<unsigned long>(0));
  vector<int> sendRank(nRankSend);

  for (int i = 0; i < nRankSend; ++i) {
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    sendRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);
    sendElem[i].resize(sizeMess);

    SU2_MPI::Recv(sendElem[i].data(), sizeMess, MPI_UNSIGNED_LONG, sendRank[i], rank, SU2_MPI::GetComm(), &status);

    for (int j = 0; j < sizeMess; ++j) sendElem[i][j] -= beg_node[rank];
  }

  /* Complete the non-blocking sends. Synchronize the processors afterwards,
     because wild cards have been used in the communication. */
  SU2_MPI::Waitall(nRankRecv, sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Communicate the data of the externals. This data is the    ---*/
  /*---         time level and the number of DOFs of the element.          ---*/
  /*--------------------------------------------------------------------------*/

#ifdef HAVE_MPI

  /* Define the send buffers and resize the vector of send requests. */
  vector<vector<unsigned short> > sendBuf(nRankSend, vector<unsigned short>(0));
  sendReqs.resize(nRankSend);

  /* Define the return buffers and the vector of return requests. */
  vector<vector<unsigned short> > returnBuf(nRankRecv, vector<unsigned short>(0));
  vector<SU2_MPI::Request> returnReqs(nRankRecv);

  /*--- Copy the information of the time level and number of DOFs into the
        send buffers and send the data using non-blocking sends. ---*/
  for (int i = 0; i < nRankSend; ++i) {
    sendBuf[i].resize(2 * sendElem[i].size());
    for (unsigned long j = 0; j < sendElem[i].size(); ++j) {
      sendBuf[i][2 * j] = elem[sendElem[i][j]]->GetTimeLevel();
      sendBuf[i][2 * j + 1] = elem[sendElem[i][j]]->GetNDOFsSol();
    }

    SU2_MPI::Isend(sendBuf[i].data(), sendBuf[i].size(), MPI_UNSIGNED_SHORT, sendRank[i], sendRank[i],
                   SU2_MPI::GetComm(), &sendReqs[i]);
  }

  /*--- Receive the data for the externals. As this data is needed immediately,
        blocking communication is used. The time level and the number of DOFs
        of the externals is stored in the second entry of the
        mapExternalElemIDToTimeLevel, which is set accordingly. ---*/
  MRI = mapRankToIndRecv.begin();
  for (int i = 0; i < nRankRecv; ++i, ++MRI) {
    returnBuf[i].resize(2 * recvElem[i].size());
    SU2_MPI::Status status;
    SU2_MPI::Recv(returnBuf[i].data(), returnBuf[i].size(), MPI_UNSIGNED_SHORT, MRI->first, rank, SU2_MPI::GetComm(),
                  &status);

    for (unsigned long j = 0; j < recvElem[i].size(); ++j) {
      MI = mapExternalElemIDToTimeLevel.find(recvElem[i][j]);
      if (MI == mapExternalElemIDToTimeLevel.end())
        SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
      MI->second.short0 = returnBuf[i][2 * j];
      MI->second.short1 = returnBuf[i][2 * j + 1];
    }
  }

  /* Complete the nonblocking sends. */
  SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Only when time accurate local time stepping is employed.   ---*/
  /*---         Iterative algorithm to make sure that neighboring elements ---*/
  /*---         differ at most one time level. This is no fundamental      ---*/
  /*---         issue, but it makes the parallel implementation easier,    ---*/
  /*---         while the penalty in efficiency should be small for most   ---*/
  /*---         cases. Furthermore, also for practical reasons, make sure  ---*/
  /*---         that the donor elements for the wall function treatment    ---*/
  /*---         have the same time level as the element to which the       ---*/
  /*---         boundary face belongs.                                     ---*/
  /*--------------------------------------------------------------------------*/

  /* Test for time accurate local time stepping. */
  if (nTimeLevels > 1) {
    /* Infinite loop for the iterative algorithm. */
    for (;;) {
      /* Variable to indicate whether the local situation has changed. */
      unsigned short localSituationChanged = 0;

      /*--- Loop over the boundary markers and make sure that the time level
            of the element adjacent to the boundary and the donor element for
            the wall function treatment is the same. This is not much of a
            restriction, because in the vast majority of cases, the donor
            element is this adjacent element itself. Requiring it to be the
            same makes the implementation of the wall functions easier. ---*/
      for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
        for (unsigned long l = 0; l < nElem_Bound[iMarker]; ++l) {
          /* Determine the ID of the adjacent element. */
          const unsigned long elemID = bound[iMarker][l]->GetDomainElement() - beg_node[rank];

          /* Get the number of donor elements for the wall function treatment
             and the pointer to the array which stores this info. */
          const unsigned short nDonors = bound[iMarker][l]->GetNDonorsWallFunctions();
          const unsigned long* donors = bound[iMarker][l]->GetDonorsWallFunctions();

          /* Loop over the number of donors and check the time levels. */
          for (unsigned short i = 0; i < nDonors; ++i) {
            /* Determine the status of the donor element. */
            if (donors[i] >= beg_node[rank] && donors[i] < beg_node[rank] + nElem) {
              /* Donor is stored locally. Determine its local ID and
                 get the time levels of both elements. */
              const unsigned long donorID = donors[i] - beg_node[rank];
              const unsigned short timeLevelB = elem[elemID]->GetTimeLevel();
              const unsigned short timeLevelD = elem[donorID]->GetTimeLevel();
              const unsigned short timeLevel = min(timeLevelB, timeLevelD);

              /* If the time level of either element is larger than timeLevel,
                 adapt the time levels and indicate that the local situation
                 has changed. */
              if (timeLevelB > timeLevel) {
                elem[elemID]->SetTimeLevel(timeLevel);
                localSituationChanged = 1;
              }

              if (timeLevelD > timeLevel) {
                elem[donorID]->SetTimeLevel(timeLevel);
                localSituationChanged = 1;
              }
            } else {
              /* The donor element is stored on a different processor.
                 Retrieve its time level from mapExternalElemIDToTimeLevel
                 and determine the minimum time level. */
              const unsigned short timeLevelB = elem[elemID]->GetTimeLevel();
              MI = mapExternalElemIDToTimeLevel.find(donors[i]);
              if (MI == mapExternalElemIDToTimeLevel.end())
                SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
              const unsigned short timeLevel = min(timeLevelB, MI->second.short0);

              /* If the time level of either element is larger than timeLevel,
                 adapt the time levels and indicate that the local situation
                 has changed. */
              if (timeLevelB > timeLevel) {
                elem[elemID]->SetTimeLevel(timeLevel);
                localSituationChanged = 1;
              }

              if (MI->second.short0 > timeLevel) {
                MI->second.short0 = timeLevel;
                localSituationChanged = 1;
              }
            }
          }
        }
      }

      /*--- Loop over the matching faces and update the time levels of the
            adjacent elements, if needed. ---*/
      for (auto FI = localFaces.begin(); FI != localFaces.end(); ++FI) {
        /* Safeguard against non-matching faces. */
        if (FI->elemID1 < Global_nElem) {
          /* Local element ID of the first element. Per definition this is
             always a locally stored element. Also store its time level. */
          const unsigned long elemID0 = FI->elemID0 - beg_node[rank];
          const unsigned short timeLevel0 = elem[elemID0]->GetTimeLevel();

          /* Determine the status of the second element. */
          if (FI->elemID1 >= beg_node[rank] && FI->elemID1 < beg_node[rank] + nElem) {
            /* Both elements are stored locally. Determine the local
               element of the second element and determine the minimum
               time level. */
            const unsigned long elemID1 = FI->elemID1 - beg_node[rank];
            const unsigned short timeLevel1 = elem[elemID1]->GetTimeLevel();
            const unsigned short timeLevel = min(timeLevel0, timeLevel1);

            /* If the time level of either element is larger than timeLevel+1,
               adapt the time levels and indicate that the local situation
               has changed. */
            if (timeLevel0 > timeLevel + 1) {
              elem[elemID0]->SetTimeLevel(timeLevel + 1);
              localSituationChanged = 1;
            }

            if (timeLevel1 > timeLevel + 1) {
              elem[elemID1]->SetTimeLevel(timeLevel + 1);
              localSituationChanged = 1;
            }
          } else {
            /* The second element is stored on a different processor.
               Retrieve its time level from mapExternalElemIDToTimeLevel
               and determine the minimum time level. */
            MI = mapExternalElemIDToTimeLevel.find(FI->elemID1);
            if (MI == mapExternalElemIDToTimeLevel.end())
              SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
            const unsigned short timeLevel = min(timeLevel0, MI->second.short0);

            /* If the time level of either element is larger than timeLevel+1,
               adapt the time levels and indicate that the local situation
               has changed. */
            if (timeLevel0 > timeLevel + 1) {
              elem[elemID0]->SetTimeLevel(timeLevel + 1);
              localSituationChanged = 1;
            }

            if (MI->second.short0 > timeLevel + 1) {
              MI->second.short0 = timeLevel + 1;
              localSituationChanged = 1;
            }
          }
        }
      }

      /* Determine whether or not the global situation changed. If not
         the infinite loop can be terminated. */
      unsigned short globalSituationChanged = localSituationChanged;

#ifdef HAVE_MPI
      SU2_MPI::Allreduce(&localSituationChanged, &globalSituationChanged, 1, MPI_UNSIGNED_SHORT, MPI_MAX,
                         SU2_MPI::GetComm());
#endif
      if (!globalSituationChanged) break;

        /*--- Communicate the information of the externals, if needed. ---*/
#ifdef HAVE_MPI

      /*--- Copy the information of the time level into the send buffers
            and send the data using non-blocking sends. Note that the size
            of sendElem is used and not sendBuf, because the size of sendBuf
            is twice the size, see step 4.  ---*/
      for (int i = 0; i < nRankSend; ++i) {
        for (unsigned long j = 0; j < sendElem[i].size(); ++j) sendBuf[i][j] = elem[sendElem[i][j]]->GetTimeLevel();

        SU2_MPI::Isend(sendBuf[i].data(), sendElem[i].size(), MPI_UNSIGNED_SHORT, sendRank[i], sendRank[i],
                       SU2_MPI::GetComm(), &sendReqs[i]);
      }

      /*--- Receive the data for the externals. As this data is needed
            immediately, blocking communication is used. The time level of the
            externals is stored in mapExternalElemIDToTimeLevel, whose second
            element, which contains the time level, is updated accordingly.
            Note that the minimum of the two levels is taken, such that changes
            for the external are also incorporated. As such this information
            must be returned to the original rank. Also note that the size
            of recvElem is used and not returnBuf, because the latter is twice
            as large, see step 4. ---*/
      MRI = mapRankToIndRecv.begin();
      for (int i = 0; i < nRankRecv; ++i, ++MRI) {
        SU2_MPI::Status status;
        SU2_MPI::Recv(returnBuf[i].data(), recvElem[i].size(), MPI_UNSIGNED_SHORT, MRI->first, rank, SU2_MPI::GetComm(),
                      &status);

        for (unsigned long j = 0; j < recvElem[i].size(); ++j) {
          MI = mapExternalElemIDToTimeLevel.find(recvElem[i][j]);
          if (MI == mapExternalElemIDToTimeLevel.end())
            SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
          MI->second.short0 = min(returnBuf[i][j], MI->second.short0);
          returnBuf[i][j] = MI->second.short0;
        }

        SU2_MPI::Isend(returnBuf[i].data(), recvElem[i].size(), MPI_UNSIGNED_SHORT, MRI->first, MRI->first + 1,
                       SU2_MPI::GetComm(), &returnReqs[i]);
      }

      /* Complete the first round of nonblocking sends, such that the
         original send buffers can be used as receive buffers for the
         second round of messages. */
      SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);

      /* Loop again over the original sending processors to receive
         the updated time level of elements that may have been updated
         on other ranks. */
      for (int i = 0; i < nRankSend; ++i) {
        SU2_MPI::Status status;
        SU2_MPI::Recv(sendBuf[i].data(), sendElem[i].size(), MPI_UNSIGNED_SHORT, sendRank[i], rank + 1,
                      SU2_MPI::GetComm(), &status);

        for (unsigned long j = 0; j < sendElem[i].size(); ++j) elem[sendElem[i][j]]->SetTimeLevel(sendBuf[i][j]);
      }

      /* Complete the second round of nonblocking sends. */
      SU2_MPI::Waitall(nRankRecv, returnReqs.data(), MPI_STATUSES_IGNORE);
#endif
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 6: Print a nice output about the number of elements per     ---*/
    /*---         time level.                                              ---*/
    /*------------------------------------------------------------------------*/

    if (rank == MASTER_NODE)
      cout << endl << "------- Element distribution for time accurate local time stepping ------" << endl;

    /* Determine the local number of elements per time level. */
    vector<unsigned long> nLocalElemPerLevel(nTimeLevels, 0);
    for (unsigned long i = 0; i < nElem; ++i) ++nLocalElemPerLevel[elem[i]->GetTimeLevel()];

    /* Determine the global version of nLocalElemPerLevel. This only needs to
       be known on the master node. */
    vector<unsigned long> nGlobalElemPerLevel = nLocalElemPerLevel;

#ifdef HAVE_MPI
    SU2_MPI::Reduce(nLocalElemPerLevel.data(), nGlobalElemPerLevel.data(), nTimeLevels, MPI_UNSIGNED_LONG, MPI_SUM,
                    MASTER_NODE, SU2_MPI::GetComm());
#endif

    /* Write the output. */
    if (rank == MASTER_NODE) {
      for (unsigned short i = 0; i < nTimeLevels; ++i) {
        if (nGlobalElemPerLevel[i])
          cout << "Number of elements time level " << i << ": " << nGlobalElemPerLevel[i] << endl;
      }
    }
  }
}

void CPhysicalGeometry::ComputeFEMGraphWeights(CConfig* config, const vector<CFaceOfElement>& localFaces,
                                               const vector<vector<unsigned long> >& adjacency,
                                               const map<unsigned long, CUnsignedShort2T>& mapExternalElemIDToTimeLevel,
                                               vector<su2double>& vwgt, vector<vector<su2double> >& adjwgt) {
  /*--- Determine the maximum time level that occurs in the grid. ---*/
  unsigned short maxTimeLevel = 0;
  for (unsigned long i = 0; i < nElem; ++i) maxTimeLevel = max(maxTimeLevel, elem[i]->GetTimeLevel());

#ifdef HAVE_MPI
  unsigned short maxTimeLevelLocal = maxTimeLevel;
  SU2_MPI::Allreduce(&maxTimeLevelLocal, &maxTimeLevel, 1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

  /*--------------------------------------------------------------------------*/
  /* Step 1: Determine the vertex weights of the graph. Per element two       */
  /*         weights are determined. The first weight is proportional to the  */
  /*         amount of work for the volume element. The second weight is the  */
  /*         number of DOFs of the element, such that the number of DOFs per  */
  /*         rank will also be the same.                                      */
  /*--------------------------------------------------------------------------*/

  /*--- Define the standard elements for the volume elements, the boundary faces
        and the matching internal faces. ---*/
  vector<CFEMStandardElement> standardElements;
  vector<CFEMStandardBoundaryFace> standardBoundaryFaces;
  vector<CFEMStandardInternalFace> standardMatchingFaces;

  /*--- Loop over the elements to determine the amount of computational work.
        This amount has a contribution from both the volume integral and
        surface integral to allow for Discontinuous and Continuous Galerkin
        schemes. For the latter the contribution of the surface integral will
        be negligible to the total amount of work.         ---*/
  for (unsigned long i = 0; i < nElem; ++i) {
    /*------------------------------------------------------------------------*/
    /*--- Determine the computational weight of the volume integral in the ---*/
    /*--- (DG-)FEM formulation.                                            ---*/
    /*------------------------------------------------------------------------*/

    /* Easier storage of the two indices for the vertex weights of this element. */
    const unsigned long ind0 = 2 * i;
    const unsigned long ind1 = ind0 + 1;

    /*--- Determine the corresponding standard element for this volume element.
          Create it, if it does not exist. ---*/
    unsigned short VTK_Type = elem[i]->GetVTK_Type();
    unsigned short nPolySol = elem[i]->GetNPolySol();
    bool JacIsConstant = elem[i]->GetJacobianConsideredConstant();

    unsigned long ii;
    for (ii = 0; ii < standardElements.size(); ++ii) {
      if (standardElements[ii].SameStandardElement(VTK_Type, nPolySol, JacIsConstant)) break;
    }

    if (ii == standardElements.size()) standardElements.emplace_back(VTK_Type, nPolySol, JacIsConstant, config);

    /* Initialize the computational work for this element, which is stored
       in the 1st vertex weight. */
    vwgt[ind0] = standardElements[ii].WorkEstimateMetis(config);

    /*------------------------------------------------------------------------*/
    /*--- Determine the computational weight of the surface integral in    ---*/
    /*--- the (DG-)FEM formulation.                                        ---*/
    /*------------------------------------------------------------------------*/

    /* Determine the global element ID of this element. */
    unsigned long elemID = beg_node[rank] + i;

    /*--- Get the global IDs of the corner points of all the faces of this element. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long faceConn[6][4];

    elem[i]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Loop over the number of faces of this element. ---*/
    for (unsigned short j = 0; j < nFaces; ++j) {
      /* Determine the VTK type of the face. */
      unsigned short VTK_Type_Face = 0;
      switch (nPointsPerFace[j]) {
        case 2:
          VTK_Type_Face = LINE;
          break;
        case 3:
          VTK_Type_Face = TRIANGLE;
          break;
        case 4:
          VTK_Type_Face = QUADRILATERAL;
          break;
      }

      /* Easier storage whether or not the Jacobian is constant. */
      JacIsConstant = elem[i]->GetJacobianConstantFace(j);

      /*--- Check if this is a matching internal face. ---*/
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[j];
      for (unsigned short k = 0; k < nPointsPerFace[j]; ++k) thisFace.cornerPoints[k] = faceConn[j][k];
      thisFace.CreateUniqueNumbering();

      vector<CFaceOfElement>::const_iterator low;
      low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

      bool thisFaceFound = false;
      if (low != localFaces.end()) {
        if (!(thisFace < *low)) thisFaceFound = true;
      }

      if (thisFaceFound) {
        /* Check if this internal matching face is owned by this element. */
        bool faceIsOwned;
        if (elemID == low->elemID0)
          faceIsOwned = low->elem0IsOwner;
        else
          faceIsOwned = !low->elem0IsOwner;

        if (faceIsOwned) {
          /*--- Determine the index of the corresponding standard matching
                face. If it does not exist, create it. ---*/
          for (ii = 0; ii < standardMatchingFaces.size(); ++ii) {
            if (standardMatchingFaces[ii].SameStandardMatchingFace(VTK_Type_Face, JacIsConstant, low->elemType0,
                                                                   low->nPolySol0, low->elemType1, low->nPolySol1,
                                                                   false, false))
              break;
          }

          if (ii == standardMatchingFaces.size())
            standardMatchingFaces.emplace_back(VTK_Type_Face, low->elemType0, low->nPolySol0, low->elemType1,
                                               low->nPolySol1, JacIsConstant, false, false, config);

          /* Update the computational work for this element, i.e. the 1st
             vertex weight. */
          vwgt[ind0] += standardMatchingFaces[ii].WorkEstimateMetis(config);
        }
      } else {
        /*--- This is a boundary face, which is owned by definition. Determine
              the index of the corresponding standard boundary face. If it
              does not exist, create it. ---*/
        for (ii = 0; ii < standardBoundaryFaces.size(); ++ii) {
          if (standardBoundaryFaces[ii].SameStandardBoundaryFace(VTK_Type_Face, JacIsConstant, VTK_Type, nPolySol,
                                                                 false))
            break;
        }

        if (ii == standardBoundaryFaces.size())
          standardBoundaryFaces.emplace_back(VTK_Type_Face, VTK_Type, nPolySol, JacIsConstant, false, config);

        /* Update the computational work for this element, i.e. the 1st
           vertex weight. */
        vwgt[ind0] += standardBoundaryFaces[ii].WorkEstimateMetis(config);
      }
    }

    /* Set the value of the second vertex weight to the number of DOFs. */
    vwgt[ind1] = elem[i]->GetNDOFsSol();
  }

  /*--------------------------------------------------------------------------*/
  /*--- Check for boundary faces for which a wall function treatment       ---*/
  /*--- must be used. This type of boundary condition treatment is         ---*/
  /*--- computationally intensive and must be added to the vertex weight   ---*/
  /*--- of the corresponding element.                                      ---*/
  /*--------------------------------------------------------------------------*/

  /* Loop over the markers and select the ones for which a wall function
     treatment must be carried out. */
  for (unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX: {
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) {
          /* Retrieve the integer information for this boundary marker.
             The number of points in normal direction for the wall function
             treatment is the first element of this array. */
          const unsigned short* shortInfo = config->GetWallFunction_IntInfo(Marker_Tag);

          /* Loop over the local boundary elements for this marker. */
          for (unsigned long l = 0; l < nElem_Bound[iMarker]; ++l) {
            /* Get the required data to define the standard boundary element. */

            /* Easier storage of the element type, the corresponding volume
               element, the polynomial degree for the solution and whether or
               not the Jacobian can be considered constant. */
            const unsigned short VTK_Type_Face = bound[iMarker][l]->GetVTK_Type();
            const unsigned long elemID = bound[iMarker][l]->GetDomainElement() - beg_node[rank];
            const unsigned short nPolySol = elem[elemID]->GetNPolySol();
            const unsigned short VTK_Type_Elem = elem[elemID]->GetVTK_Type();
            const bool JacIsConstant = bound[iMarker][l]->GetJacobianConsideredConstant();

            /* Determine the corresponding entry in the standard elements for
               boundary faces. This entry must be found. */
            unsigned long ii;
            for (ii = 0; ii < standardBoundaryFaces.size(); ++ii) {
              if (standardBoundaryFaces[ii].SameStandardBoundaryFace(VTK_Type_Face, JacIsConstant, VTK_Type_Elem,
                                                                     nPolySol, false))
                break;
            }

            if (ii == standardBoundaryFaces.size())
              SU2_MPI::Error("No matching standard boundary element found", CURRENT_FUNCTION);

            /* Update the computational work for the corresponding volume element,
               i.e. the 1st vertex weight. */
            vwgt[2 * elemID] += standardBoundaryFaces[ii].WorkEstimateMetisWallFunctions(config, shortInfo[0]);
          }
        }

        break;
      }

      default: /* Just to avoid a compiler warning. */
        break;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- The final weight is obtained by taking the amount of work in time  ---*/
  /*--- into account. Note that this correction is only relevant when time ---*/
  /*--- accurate local time stepping is employed.                          ---*/
  /*--------------------------------------------------------------------------*/

  for (unsigned long i = 0; i < nElem; ++i) {
    const unsigned short diffLevel = maxTimeLevel - elem[i]->GetTimeLevel();
    vwgt[2 * i] *= pow(2, diffLevel);
  }

  /*--- Determine the minimum of the workload of the elements, i.e. 1st vertex
        weight, over the entire domain. ---*/
  su2double minvwgt = vwgt[0];
  for (unsigned long i = 0; i < nElem; ++i) minvwgt = min(minvwgt, vwgt[2 * i]);

#ifdef HAVE_MPI
  su2double locminvwgt = minvwgt;
  SU2_MPI::Allreduce(&locminvwgt, &minvwgt, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
#endif

  /*--- Scale the workload of the elements, the 1st vertex weight, with the
        minimum value and multiply by 100, such that the conversion to integer
        weights (ParMETIS is using integers for the weights) does not lead to
        a significant increase in the load imbalance. ---*/
  minvwgt = 100.0 / minvwgt;
  for (unsigned long i = 0; i < nElem; ++i) vwgt[2 * i] *= minvwgt;

  /*--------------------------------------------------------------------------*/
  /* Step 2: Determine the adjacency weights, which are proportional to the   */
  /*         amount of communication needed when the two neighboring          */
  /*         elements are stored on different ranks.                          */
  /*--------------------------------------------------------------------------*/

  /*--- Determine the edge weights by a double loop over the local elements
        and its edges. The map mapElemIDsToFaceInd is used to determine
        the corresponding index in localFaces.   ---*/
  for (unsigned long i = 0; i < nElem; ++i) {
    /* Easier storage of the time level of the current element and the
       number of solution DOFs. */
    const unsigned short timeLevel0 = elem[i]->GetTimeLevel();
    const unsigned short nDOFs0 = elem[i]->GetNDOFsSol();

    /* Loop over the number of entries in the graph for this element. */
    for (unsigned long j = 0; j < adjacency[i].size(); ++j) {
      /* Define the variables for the time level and number of solution
         DOFs for the neighboring element. */
      unsigned short timeLevel1, nDOFs1;

      /* Check if the neighor is stored locally. */
      if (adjacency[i][j] >= beg_node[rank] && adjacency[i][j] < beg_node[rank] + nElem) {
        /* Locally stored element. Determine its local ID and set the
           time level and number of solution DOFs. */
        unsigned long elemID1 = adjacency[i][j] - beg_node[rank];
        timeLevel1 = elem[elemID1]->GetTimeLevel();
        nDOFs1 = elem[elemID1]->GetNDOFsSol();
      } else {
        /* The neighbor is an external element. Find it in mapExternalElemIDToTimeLevel
           and set the time level and number of solution DOFs accordingly. */
        map<unsigned long, CUnsignedShort2T>::const_iterator MI;
        MI = mapExternalElemIDToTimeLevel.find(adjacency[i][j]);
        if (MI == mapExternalElemIDToTimeLevel.end())
          SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
        timeLevel1 = MI->second.short0;
        nDOFs1 = MI->second.short1;
      }

      /* Determine the difference of the maximum time level that occurs and
         the minimum of the time level of the current and the adjacent element.
         This value can only be nonzero when time accurate local time stepping
         is employed. */
      const unsigned short diffLevel = maxTimeLevel - min(timeLevel0, timeLevel1);

      /* Set the edge weight. As ParMetis expects an undirected graph, set the edge weight
         to the sum of the number of DOFs on both sides, multiplied by the weight factor
         to account for different time levels. */
      adjwgt[i][j] = pow(2, diffLevel) * (nDOFs0 + nDOFs1);
    }
  }
}
