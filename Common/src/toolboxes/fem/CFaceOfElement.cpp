/*!
 * \file CFaceOfElement.cpp
 * \brief Helper class used in distributing the surface elements and
 *        creating the surface elements for the DG solver.
 * \author E. van der Weide
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

#include "../../../include/toolboxes/fem/CFaceOfElement.hpp" 


CFaceOfElement::CFaceOfElement() {

  /*--- Set the default values. ---*/
  nCornerPoints   = 0;
  cornerPoints[0] = cornerPoints[1] = cornerPoints[2] = cornerPoints[3] = ULONG_MAX;
  elemID0         = elemID1 = ULONG_MAX;
  nPolyGrid0      = nPolyGrid1 = 0;
  nPolySol0       = nPolySol1  = 0;
  nDOFsElem0      = nDOFsElem1 = 0;
  elemType0       = elemType1  = 0;
  faceID0         = faceID1    = 0;
  periodicIndex   = periodicIndexDonor = 0;
  faceIndicator   = 0;

  JacFaceIsConsideredConstant = false;
  elem0IsOwner                = false;
}

CFaceOfElement::CFaceOfElement(const unsigned short VTK_Type,
                               const unsigned short nPoly,
                               const unsigned long  *Nodes) {

  /*--- Set the default values of the member variables. ---*/
  nCornerPoints   = 0;
  cornerPoints[0] = cornerPoints[1] = cornerPoints[2] = cornerPoints[3] = ULONG_MAX;
  elemID0         = elemID1 = ULONG_MAX;
  nPolyGrid0      = nPolyGrid1 = 0;
  nPolySol0       = nPolySol1  = 0;
  nDOFsElem0      = nDOFsElem1 = 0;
  elemType0       = elemType1  = 0;
  faceID0         = faceID1    = 0;
  periodicIndex   = periodicIndexDonor = 0;
  faceIndicator   = 0;

  JacFaceIsConsideredConstant = false;
  elem0IsOwner                = false;

  /*--- Determine the face element type and set the corner points accordingly. ---*/
  switch( VTK_Type ) {
    case LINE: {
      nCornerPoints   = 2;
      cornerPoints[0] = Nodes[0];
      cornerPoints[1] = Nodes[nPoly];
      break;
    }

    case TRIANGLE: {
      const unsigned short ind2 = (nPoly+1)*(nPoly+2)/2 -1;
      nCornerPoints   = 3;
      cornerPoints[0] = Nodes[0];
      cornerPoints[1] = Nodes[nPoly];
      cornerPoints[2] = Nodes[ind2];
      break;
    }

    case QUADRILATERAL: {
      const unsigned short ind2 = nPoly*(nPoly+1);
      const unsigned short ind3 = (nPoly+1)*(nPoly+1) -1;
      nCornerPoints   = 4;
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

bool CFaceOfElement::operator<(const CFaceOfElement &other) const {

  /*--- First compare the number of corner points. ---*/
  if(nCornerPoints != other.nCornerPoints) return nCornerPoints < other.nCornerPoints;

  /*--- Compare the corner points. ---*/
  for(unsigned short i=0; i<nCornerPoints; ++i)
    if(cornerPoints[i] != other.cornerPoints[i]) return cornerPoints[i] < other.cornerPoints[i];

  /*--- Faces are identical, hence return false. ---*/
  return false;
}

bool CFaceOfElement::operator ==(const CFaceOfElement &other) const {

  /*--- First compare the number of corner points. ---*/
  if(nCornerPoints != other.nCornerPoints) return false;

  /*--- Compare the corner points. ---*/
  for(unsigned short i=0; i<nCornerPoints; ++i)
    if(cornerPoints[i] != other.cornerPoints[i]) return false;

  /*--- Faces are identical, hence return true. ---*/
  return true;
}

void CFaceOfElement::Copy(const CFaceOfElement &other) {

  /*--- Simply copy the member variables. ---*/
  nCornerPoints = other.nCornerPoints;
  for(unsigned short i=0; i<nCornerPoints; ++i) cornerPoints[i] = other.cornerPoints[i];
  for(unsigned short i=nCornerPoints; i<4; ++i) cornerPoints[i] = ULONG_MAX;

  elemID0 = other.elemID0;
  elemID1 = other.elemID1;

  nPolyGrid0 = other.nPolyGrid0;
  nPolyGrid1 = other.nPolyGrid1;
  nPolySol0  = other.nPolySol0;
  nPolySol1  = other.nPolySol1;

  nDOFsElem0 = other.nDOFsElem0;
  nDOFsElem1 = other.nDOFsElem1;
  elemType0  = other.elemType0;
  elemType1  = other.elemType1;
  faceID0    = other.faceID0;
  faceID1    = other.faceID1;

  periodicIndex      = other.periodicIndex;
  periodicIndexDonor = other.periodicIndexDonor;
  faceIndicator      = other.faceIndicator;

  JacFaceIsConsideredConstant = other.JacFaceIsConsideredConstant;
  elem0IsOwner                = other.elem0IsOwner;
}

void CFaceOfElement::CreateUniqueNumberingWithOrientation(void) {

  /*--- Determine the element type and create the unique numbering accordingly. ---*/
  bool swapElements = false;
  switch( nCornerPoints ) {

    case 2: {
      /*--- Element is a line. Check if the node numbering must be swapped. If so
            also the element information must be swapped, because element 0 is to
            the left of the face and element 1 to the right. ---*/
      if(cornerPoints[1] < cornerPoints[0]) {
        swap(cornerPoints[0], cornerPoints[1]);
        swapElements = true;
      }
      break;
    }

    case 3: {
      /*--- Element is a triangle. The vertices are sorted in increasing order.
            If the sequence of the new numbering is opposite to the current
            numbering, the element information must be exchanged, because
            element 0 is to the left of the face and element 1 to the right. ---*/
      unsigned long nn[] = {cornerPoints[0], cornerPoints[1], cornerPoints[2]};
      unsigned short ind = 0;
      if(nn[1] < nn[ind]) ind = 1;
      if(nn[2] < nn[ind]) ind = 2;

      unsigned short indm1 = ind==0 ? 2 : ind - 1;   // Next lower index.
      unsigned short indp1 = ind==2 ? 0 : ind + 1;   // Next upper index.

      if(nn[indp1] < nn[indm1]) {

        /*--- The orientation of the triangle remains the same.
              Store the new sorted node numbering. ---*/
        cornerPoints[0] = nn[ind];
        cornerPoints[1] = nn[indp1];
        cornerPoints[2] = nn[indm1];
      }
      else {

        /*--- The orientation of the triangle changes. Store the new
              sorted node numbering and set swapElements to true. ---*/
        cornerPoints[0] = nn[ind];
        cornerPoints[1] = nn[indm1];
        cornerPoints[2] = nn[indp1];
        swapElements    = true;
      }

      break;
    }

    case 4: {
      /*--- Element is a quadrilateral. The vertices are sorted in increasing order
            under the condition neighboring vertices remain neighbors. If the
            sequence of the new numbering is opposite to the current
            numbering, the element information must be exchanged, because
            element 0 is to the left of the face and element 1 to the right. ---*/
      unsigned long nn[] = {cornerPoints[0], cornerPoints[1],
                            cornerPoints[2], cornerPoints[3]};
      unsigned short ind = 0;
      if(nn[1] < nn[ind]) ind = 1;
      if(nn[2] < nn[ind]) ind = 2;
      if(nn[3] < nn[ind]) ind = 3;

      unsigned short indm1 = ind==0 ?       3 : ind - 1; // Next lower index.
      unsigned short indp1 = ind==3 ?       0 : ind + 1; // Next upper index.
      unsigned short indp2 = ind>=2 ? ind - 2 : ind + 2; // Opposite index.

      if(nn[indp1] < nn[indm1]) {

        /*--- The orientation of the quadrilateral remains the same.
              Store the new sorted node numbering. ---*/
        cornerPoints[0] = nn[ind];
        cornerPoints[1] = nn[indp1];
        cornerPoints[2] = nn[indp2];
        cornerPoints[3] = nn[indm1];
      }
      else {

        /*--- The orientation of the quadrilateral changes. Store the new
              sorted node numbering and set swapElements to true. ---*/
        cornerPoints[0] = nn[ind];
        cornerPoints[1] = nn[indm1];
        cornerPoints[2] = nn[indp2];
        cornerPoints[3] = nn[indp1];
        swapElements    = true;
      }

      break;
    }

    default: {
      ostringstream message;
      message << "Unknown surface element type with " << nCornerPoints
              << " corners." << endl;
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
    }
  }

  /*--- Swap the element information, if needed. ---*/
  if( swapElements ) {
    swap(elemID0,    elemID1);
    swap(nPolyGrid0, nPolyGrid1);
    swap(nPolySol0,  nPolySol1);
    swap(nDOFsElem0, nDOFsElem1);
    swap(elemType0,  elemType1);
    swap(faceID0,    faceID1);
  }
}

unsigned short CFaceOfElement::DetermineOrientationElemSide1(
                               const vector<CVolumeElementFEM_DG> &volElem) const {

  /*--- Determine the corner points of all the faces of the element on side 1. ---*/
  unsigned short nFaces;
  unsigned short nPointsPerFace[6];
  unsigned long  faceConn[6][4];

  volElem[elemID1].GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

  /*--- Define the return variable orientation and initialize it to avoid a
        compiler warning. Also initialize the bool wrongOrientation to false. ---*/
  unsigned short orientation      = 0;
  bool           wrongOrientation = false;

  /*--- Make a distinction between the possible surface elements. ---*/
  switch( nCornerPoints ) {

    case 2: {

      /*---- Element is a line. Abbreviate the points of the face. ---*/
      const unsigned long n0 = faceConn[faceID1][0], n1 = faceConn[faceID1][1];

      /*--- Determine the situation and set orientation accordingly. ---*/
      if((n0 == cornerPoints[1]) && (n1 == cornerPoints[0])) orientation = 1;
      else wrongOrientation = true;

      break;
    }

    case 3: {

      /*---- Element is a triangle. Abbreviate the points of the face. ---*/
      const unsigned long n0 = faceConn[faceID1][0], n1 = faceConn[faceID1][1],
                          n2 = faceConn[faceID1][2];

      /*--- Determine the situation and set orientation accordingly. ---*/
      if(     (n0 == cornerPoints[0]) && (n1 == cornerPoints[2]) &&
              (n2 == cornerPoints[1])) orientation = 1;
      else if((n0 == cornerPoints[1]) && (n1 == cornerPoints[0]) &&
              (n2 == cornerPoints[2])) orientation = 2;
      else if((n0 == cornerPoints[2]) && (n1 == cornerPoints[1]) &&
              (n2 == cornerPoints[0])) orientation = 3;
      else wrongOrientation = true;

      break;
    }

    case 4: {

      /*---- Element is a quadrilateral. Abbreviate the points of the face. ---*/
      const unsigned long n0 = faceConn[faceID1][0], n1 = faceConn[faceID1][1],
                          n2 = faceConn[faceID1][2], n3 = faceConn[faceID1][3];

      /*--- Determine the situation and set orientation accordingly. ---*/
      if(     (n0 == cornerPoints[0]) && (n1 == cornerPoints[3]) &&
              (n2 == cornerPoints[2]) && (n3 == cornerPoints[1])) orientation = 1;
      else if((n0 == cornerPoints[1]) && (n1 == cornerPoints[0]) &&
              (n2 == cornerPoints[3]) && (n3 == cornerPoints[2])) orientation = 2;
      else if((n0 == cornerPoints[2]) && (n1 == cornerPoints[1]) &&
              (n2 == cornerPoints[0]) && (n3 == cornerPoints[3])) orientation = 3;
      else if((n0 == cornerPoints[3]) && (n1 == cornerPoints[2]) &&
              (n2 == cornerPoints[1]) && (n3 == cornerPoints[0])) orientation = 4;
      else wrongOrientation = true;

      break;
    }

    default: {
      ostringstream message;
      message << "Unknown surface element type with " << nCornerPoints
              << " corners." << endl;
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
    }
  }

  /*--- Check if something went wrong. ---*/
  if( wrongOrientation )
    SU2_MPI::Error(string("No matching orientation found"), CURRENT_FUNCTION);

  /*--- Return the value of orientation. ---*/
  return orientation;
}

void CFaceOfElement::MatchOrientationElemSide0(const vector<CVolumeElementFEM_DG> &volElem) {

  /*--- Determine the corner points of all the faces of the element on side 0. ---*/
  unsigned short nFaces;
  unsigned short nPointsPerFace[6];
  unsigned long  faceConn[6][4];

  volElem[elemID0].GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

  /*--- Reset the node number of this face to the node numbering of the
        corresponding face of the element. ---*/
  for(unsigned short j=0; j<nCornerPoints; ++j)
    cornerPoints[j] = faceConn[faceID0][j];
}

void CFaceOfElement::SwapSidesIfNeeded(const unsigned long nVolElemOwned,
                                       const unsigned long nVolElemTot) {
 
  /*--- Definition of the boolean, which defines whether or not the adjacent
        elements must be swapped. ---*/
  bool swapElements;

  /*--- Check if this is an internal matching face. ---*/
  if(elemID0 < nVolElemTot && elemID1 < nVolElemTot) {

    /*--- Internal matching face. Make sure that the element with
          the highest polynomial degree of the grid is on side 0.
          In this way the coordinates and the normals can be
          computed from the info from the element on side 0. ---*/
    if(nPolyGrid0 != nPolyGrid1) {
      swapElements = nPolyGrid1 > nPolyGrid0;
    }
    else {

      /*--- The polynomial degree of the grid is the same. The following
            criteria for swapping are such that the least amount of
            standard elements for matching faces need to be created. ---*/
      if(elemType0 == elemType1) {

        /*--- Elements are of the same type. The elements are swapped if the
              face ID of element 1 is lower. In this way the element with the
              lowest face ID is always on side 0 of the face. If the face ID's
              are the same, look at the polynomial degree. ---*/
        if(faceID0 != faceID1) swapElements = faceID1 < faceID0;
        else                   swapElements = nPolySol1 > nPolySol0;
      }
      else {

        /*--- Adjacent elements are of a different type. Below follow some criteria
              to swap the elements. ---*/
        if(nCornerPoints == 2) {

          /*--- 2D simulation. The quadrilateral is always on side 0. ---*/
          swapElements = elemType1 == QUADRILATERAL;
        }
        else if(nCornerPoints == 3) {

          /*--- Triangular face. The prism is always on side 0 and
                the tetrahedron never. ---*/
          swapElements = (elemType0 == TETRAHEDRON) || (elemType1 == PRISM);
        }
        else {

          /*--- Quadrilateral face. The pyramid is always on side 0 and
                the prism never. ---*/
          swapElements = (elemType0 == PRISM) || (elemType1 == PYRAMID);
        }
      }
    }
  }
  else {

    /*--- Either a boundary face or a non-matching face. It must be swapped
          if the element is currently on side 1 of the face. ---*/
    swapElements = elemID1 < nVolElemTot;
  }

  /*--- Swap the adjacent elements of the face, if needed. Note that
        also the sequence of the corner points must be altered in order
        to obey the right hand rule. ---*/
  if( swapElements ) {
    swap(elemID0,    elemID1);
    swap(nPolyGrid0, nPolyGrid1);
    swap(nPolySol0,  nPolySol1);
    swap(nDOFsElem0, nDOFsElem1);
    swap(elemType0,  elemType1);
    swap(faceID0,    faceID1);

    if(nCornerPoints == 2) swap(cornerPoints[0], cornerPoints[1]);
    else                   swap(cornerPoints[0], cornerPoints[2]);

    elem0IsOwner = !elem0IsOwner;
  }
}
