/*!
 * \file CFaceOfElement.cpp
 * \brief Helper class used in distributing the surface elements and
 *        creating the surface elements for the DG solver.
 * \author E. van der Weide
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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
      std::ostringstream message;
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
        std::swap(cornerPoints[0], cornerPoints[1]);
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
      std::ostringstream message;
      message << "Unknown surface element type with " << nCornerPoints << " corners." << std::endl;
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
    }
  }

  /* Swap the element information, if needed. */
  if (swapElements) {
    std::swap(elemID0, elemID1);
    std::swap(nPolyGrid0, nPolyGrid1);
    std::swap(nPolySol0, nPolySol1);
    std::swap(nDOFsElem0, nDOFsElem1);
    std::swap(elemType0, elemType1);
    std::swap(faceID0, faceID1);
  }
}
