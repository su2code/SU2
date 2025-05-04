/*!
 * \file CFaceOfElement.hpp
 * \brief Header file for the class CFaceOfElement.
 *        The implementations are in the <i>CFaceOfElement.cpp</i> file.
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

#pragma once

#include "../../code_config.hpp"
#include "../../option_structure.hpp"
#include <iostream>
#include <algorithm>
#include <climits>

/*!
 * \class CFaceOfElement
 * \brief Class used in the partitioning of the FEM grid as well as the building of
          the faces of DG. It stores a face of an element.
 */
class CFaceOfElement {
 public:
  unsigned short nCornerPoints;          /*!< \brief Number of corner points of the face. */
  unsigned long cornerPoints[4];         /*!< \brief Global ID's of ther corner points. */
  unsigned long elemID0, elemID1;        /*!< \brief Element ID's to the left and right. */
  unsigned short nPolyGrid0, nPolyGrid1; /*!< \brief Polynomial degrees of the grid of the elements
                                                     to the left and right. */
  unsigned short nPolySol0, nPolySol1;   /*!< \brief Polynomial degrees of the solution of the elements
                                                     to the left and right. */
  unsigned short nDOFsElem0, nDOFsElem1; /*!< \brief Number of DOFs of the elements to the left and right. */
  unsigned short elemType0, elemType1;   /*!< \brief Type of the elements to the left and right. */
  unsigned short faceID0, faceID1;       /*!< \brief The local face ID in the corresponding elements
                                                     to the left and right of the face. */
  unsigned short periodicIndex;          /*!< \brief Periodic indicator of the face. A value of 0 means no
                                                     periodic face. A value larger than 0 gives the index of
                                                     the periodic boundary + 1. */
  unsigned short periodicIndexDonor;     /*!< \brief Periodic indicator of the donor face. A value of 0 means no
                                                     periodic donor face. A value larger than 0 gives the index of
                                                     the periodic donor boundary + 1. */
  short faceIndicator;                   /*!< \brief The corresponding boundary marker if this face is on a
                                                     boundary. For an internal face the value is -1,
                                                     while an invalidated face has the value -2. */
  bool JacFaceIsConsideredConstant;      /*!< \brief Whether or not the Jacobian of the transformation
                                                     to the standard element is considered constant. */
  bool elem0IsOwner;                     /*!< \brief Whether or not the neighboring element 0 is the owner
                                                     of the face. If false, element 1 is the owner. */

  /* Constructor. Initialize the member variables. */
  CFaceOfElement();

  /* Alternative constructor to set the corner points. */
  CFaceOfElement(const unsigned short VTK_Type, const unsigned short nPoly, const unsigned long* Nodes);

  /* Copy constructor and assignment operator. */
  inline CFaceOfElement(const CFaceOfElement& other) { Copy(other); }

  inline CFaceOfElement& operator=(const CFaceOfElement& other) {
    Copy(other);
    return (*this);
  }

  /* Less than operator. Needed for the sorting and searching. */
  bool operator<(const CFaceOfElement& other) const;

  /* Equal operator. Needed for removing double entities. */
  bool operator==(const CFaceOfElement& other) const;

  /*--- Member function, which creates a unique numbering for the corner points.
        A sort in increasing order is OK for this purpose.                       ---*/
  inline void CreateUniqueNumbering() { std::sort(cornerPoints, cornerPoints + nCornerPoints); }

  /*--- Member function, which creates a unique numbering for the corner points
        while the orientation is taken into account. ---*/
  void CreateUniqueNumberingWithOrientation();

 private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CFaceOfElement& other);
};
