/*!
 * \file geometry_structure_fem_part.hpp
 * \brief Helper classes for the Fluid FEM solver.
 * \author E. van der Weide
 * \version 7.0.5 "Blackbird"
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

#pragma once

#include "../basic_types/datatype_structure.hpp"

#include <vector>
#include <algorithm>

/*!
 * \struct CUnsignedLong2T
 * \brief Helper struct used to store two integral types as one entity.
 */
struct CUnsignedLong2T {

  unsigned long long0;  /*!< \brief First integer to store in this class. */
  unsigned long long1;  /*!< \brief Second integer to store in this class. */

  CUnsignedLong2T(unsigned long a = 0, unsigned long b = 0) : long0(a), long1(b) {}

  inline bool operator<(const CUnsignedLong2T &other) const {
    if(long0 != other.long0)
      return (long0 < other.long0);
    return (long1 < other.long1);
  }

  inline bool operator==(const CUnsignedLong2T &other) const {
    return (long0 == other.long0) && (long1 == other.long1);
  }
};

/*!
 * \struct CUnsignedShort2T
 * \brief Help struct used to store two integral types as one entity.
 */
struct CUnsignedShort2T {

  unsigned short short0;  /*!< \brief First integer to store in this class. */
  unsigned short short1;  /*!< \brief Second integer to store in this class. */

  CUnsignedShort2T(unsigned short a = 0, unsigned short b = 0) : short0(a), short1(b) {}

  inline bool operator<(const CUnsignedShort2T &other) const {
    if(short0 != other.short0)
      return (short0 < other.short0);
    return (short1 < other.short1);
  }

  inline bool operator==(const CUnsignedShort2T &other) const {
    return (short0 == other.short0) && (short1 == other.short1);
  }
};

/*!
 * \class CBoundaryFace
 * \brief Help class used in the partitioning of the FEM grid.
          It stores a boundary element.
 */
class CBoundaryFace {
public:
  unsigned short VTK_Type, nPolyGrid, nDOFsGrid;
  unsigned long  globalBoundElemID, domainElementID;
  std::vector<unsigned long>  Nodes;

  /* Standard constructor and destructor. Nothing to be done. */
  CBoundaryFace(){}
  ~CBoundaryFace(){}

  /* Copy constructor and assignment operator. */
  inline CBoundaryFace(const CBoundaryFace &other) { Copy(other); }

  inline CBoundaryFace& operator=(const CBoundaryFace &other) { Copy(other); return (*this); }

  /* Less than operator. Needed for the sorting. */
  inline bool operator<(const CBoundaryFace &other) const {
    return (globalBoundElemID < other.globalBoundElemID);
  }

private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CBoundaryFace &other);
};

/*!
 * \class CMatchingFace
 * \brief Help class used to determine whether or not (periodic) faces match.
 */
class CMatchingFace {
public:
  unsigned short nCornerPoints;          /*!< \brief Number of corner points of the face. */
  unsigned short nDim;                   /*!< \brief Number of spatial dimensions. */
  unsigned short nPoly;                  /*!< \brief Polynomial degree of the face. */
  unsigned short nDOFsElem;              /*!< \brief Number of DOFs of the relevant adjacent element. */
  unsigned short elemType;               /*!< \brief Type of the adjacent element. */
  unsigned long  elemID;                 /*!< \brief The relevant adjacent element ID. */
  su2double cornerCoor[4][3];            /*!< \brief Coordinates of the corner points of the face. */
  su2double tolForMatching;              /*!< \brief Tolerance for this face for matching points. */

  /* Standard constructor. */
  CMatchingFace();

  /* Destructor, nothing to be done. */
  ~CMatchingFace(){}

  /* Copy constructor and assignment operator. */
  inline CMatchingFace(const CMatchingFace &other) { Copy(other); }

  inline CMatchingFace& operator=(const CMatchingFace &other) { Copy(other); return (*this); }

  /* Less than operator. Needed for the sorting and searching. */
  bool operator<(const CMatchingFace &other) const;

  /*--- Member function, which sorts the coordinates of the face. ---*/
  void SortFaceCoordinates(void);

private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CMatchingFace &other);
};

