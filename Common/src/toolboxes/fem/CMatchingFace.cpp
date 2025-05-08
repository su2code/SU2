/*!
 * \file CMatchingFace.cpp
 * \brief Helper class used to determine whether or not faces of
 *        surface elements for the DG solver are matching faces.
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

#include "../../../include/toolboxes/fem/CMatchingFace.hpp"

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
