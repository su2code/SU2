/*!
 * \file geometry_structure_fem_part.hpp
 * \brief Helper classes for the Fluid FEM solver.
 * \author E. van der Weide
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

#pragma once

#include "../basic_types/datatype_structure.hpp"

#include <vector>
#include <algorithm>

/*!
 * \struct CUnsignedLong2T
 * \brief Helper struct used to store two integral types as one entity.
 */
struct CUnsignedLong2T {
  unsigned long long0; /*!< \brief First integer to store in this class. */
  unsigned long long1; /*!< \brief Second integer to store in this class. */

  CUnsignedLong2T(unsigned long a = 0, unsigned long b = 0) : long0(a), long1(b) {}

  inline bool operator<(const CUnsignedLong2T& other) const {
    if (long0 != other.long0) return (long0 < other.long0);
    return (long1 < other.long1);
  }

  inline bool operator==(const CUnsignedLong2T& other) const {
    return (long0 == other.long0) && (long1 == other.long1);
  }
};

/*!
 * \struct CUnsignedShort2T
 * \brief Help struct used to store two integral types as one entity.
 */
struct CUnsignedShort2T {
  unsigned short short0; /*!< \brief First integer to store in this class. */
  unsigned short short1; /*!< \brief Second integer to store in this class. */

  CUnsignedShort2T(unsigned short a = 0, unsigned short b = 0) : short0(a), short1(b) {}

  inline bool operator<(const CUnsignedShort2T& other) const {
    if (short0 != other.short0) return (short0 < other.short0);
    return (short1 < other.short1);
  }

  inline bool operator==(const CUnsignedShort2T& other) const {
    return (short0 == other.short0) && (short1 == other.short1);
  }
};

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

  /* Standard constructor and destructor. */
  CFaceOfElement();
  ~CFaceOfElement() {}

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
  inline void CreateUniqueNumbering(void) { std::sort(cornerPoints, cornerPoints + nCornerPoints); }

  /*--- Member function, which creates a unique numbering for the corner points
        while the orientation is taken into account. ---*/
  void CreateUniqueNumberingWithOrientation(void);

 private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CFaceOfElement& other);
};

/*!
 * \class CBoundaryFace
 * \brief Help class used in the partitioning of the FEM grid.
          It stores a boundary element.
 */
class CBoundaryFace {
 public:
  unsigned short VTK_Type, nPolyGrid, nDOFsGrid;
  unsigned long globalBoundElemID, domainElementID;
  std::vector<unsigned long> Nodes;

  /* Standard constructor and destructor. Nothing to be done. */
  CBoundaryFace() {}
  ~CBoundaryFace() {}

  /* Copy constructor and assignment operator. */
  inline CBoundaryFace(const CBoundaryFace& other) { Copy(other); }

  inline CBoundaryFace& operator=(const CBoundaryFace& other) {
    Copy(other);
    return (*this);
  }

  /* Less than operator. Needed for the sorting. */
  inline bool operator<(const CBoundaryFace& other) const { return (globalBoundElemID < other.globalBoundElemID); }

 private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CBoundaryFace& other);
};

/*!
 * \class CMatchingFace
 * \brief Help class used to determine whether or not (periodic) faces match.
 */
class CMatchingFace {
 public:
  unsigned short nCornerPoints; /*!< \brief Number of corner points of the face. */
  unsigned short nDim;          /*!< \brief Number of spatial dimensions. */
  unsigned short nPoly;         /*!< \brief Polynomial degree of the face. */
  unsigned short nDOFsElem;     /*!< \brief Number of DOFs of the relevant adjacent element. */
  unsigned short elemType;      /*!< \brief Type of the adjacent element. */
  unsigned long elemID;         /*!< \brief The relevant adjacent element ID. */
  su2double cornerCoor[4][3];   /*!< \brief Coordinates of the corner points of the face. */
  su2double tolForMatching;     /*!< \brief Tolerance for this face for matching points. */

  /* Standard constructor. */
  CMatchingFace();

  /* Destructor, nothing to be done. */
  ~CMatchingFace() {}

  /* Copy constructor and assignment operator. */
  inline CMatchingFace(const CMatchingFace& other) { Copy(other); }

  inline CMatchingFace& operator=(const CMatchingFace& other) {
    Copy(other);
    return (*this);
  }

  /* Less than operator. Needed for the sorting and searching. */
  bool operator<(const CMatchingFace& other) const;

  /*--- Member function, which sorts the coordinates of the face. ---*/
  void SortFaceCoordinates(void);

 private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CMatchingFace& other);
};
