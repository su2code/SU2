/*!
 * \file CADTBaseClass.hpp
 * \brief Base class for storing an ADT in an arbitrary number of dimensions.
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

#include <vector>
#include <array>

#include "../basic_types/datatype_structure.hpp"
#include "./CADTNodeClass.hpp"
#include "../parallelization/omp_structure.hpp"

using namespace std;

/*!
 * \class CADTBaseClass
 * \ingroup ADT
 * \brief  Base class for storing an ADT in an arbitrary number of dimensions.
 * \author E. van der Weide
 */
class CADTBaseClass {
 protected:
  unsigned long nLeaves;  /*!< \brief Number of leaves in the ADT. */
  unsigned short nDimADT; /*!< \brief Number of dimensions of the ADT. */
  bool isEmpty;           /*!< \brief Whether or not the ADT is empty. */

  vector<CADTNodeClass> leaves; /*!< \brief Vector, which contains all the leaves of the ADT. */

#ifdef HAVE_OMP
  vector<vector<unsigned long> > FrontLeaves;    /*!< \brief Vector used in the tree traversal. */
  vector<vector<unsigned long> > FrontLeavesNew; /*!< \brief Vector used in the tree traversal. */
#else
  array<vector<unsigned long>, 1> FrontLeaves;
  array<vector<unsigned long>, 1> FrontLeavesNew;
#endif
 private:
  vector<su2double> coorMinLeaves; /*!< \brief Vector, which contains all the minimum coordinates
                                               of the leaves. */
  vector<su2double> coorMaxLeaves; /*!< \brief Vector, which contains all the maximum coordinates
                                               of the leaves. */
 protected:
  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  CADTBaseClass() = default;

  /*--- Disable copy operations ---*/
  CADTBaseClass(const CADTBaseClass&) = delete;
  CADTBaseClass& operator=(const CADTBaseClass&) = delete;

  /*!
   * \brief Function, which builds the ADT of the given coordinates.
   * \param[in] nDim    Number of dimensions of the ADT.
   * \param[in] nPoints Number of points present in the ADT.
   * \param[in] coor    Coordinates of the points.
   */
  void BuildADT(unsigned short nDim, unsigned long nPoints, const su2double* coor);

 public:
  /*!
   * \brief Function, which returns whether or not the ADT is empty.
   * \return  Whether or not the ADT is empty.
   */
  inline bool IsEmpty(void) const { return isEmpty; }
};
