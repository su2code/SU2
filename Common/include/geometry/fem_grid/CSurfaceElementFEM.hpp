/*!
 * \file CSurfaceElementFEM.hpp
 * \brief Class definition for a surface element for the FEM solver.
 *        The implementations are in the <i>CSurfaceElementFEM.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "../../fem/CFEMStandardElementBase.hpp"

using namespace std;

/*!
 * \class CSurfaceeSurfaceFEM
 * \brief Class to store a surface element for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
struct CSurfaceElementFEM {
  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */

  unsigned long volElemID;         /*!< \brief ID of the corresponding volume element. */
  unsigned long boundElemIDGlobal; /*!< \brief Global ID of this surface element inside
                                               the boundary to which it belongs. */

  vector<unsigned long> nodeIDsGrid; /*!< \brief Vector with the node IDs of the grid for this element. */

  CFEMStandardElementBase *standardElemGrid = nullptr; /*!< \brief Pointer to the standard element for the grid. */
  CFEMStandardElementBase *standardElemFlow = nullptr; /*!< \brief Pointer to the standard element for the
                                                                   standard flow solution variables. */
  CFEMStandardElementBase *standardElemP    = nullptr; /*!< \brief Pointer to the standard element for the
                                                                   pressure for an incompressible flow. */

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
            The criterion for comparison is the corresponding (local) volume ID.
   */
  bool operator<(const CSurfaceElementFEM &other) const { return volElemID < other.volElemID; }

  /*!
   *  \brief Function, which determines the corner points of this surface element.
   *  \param[out] nPointsPerFace - Number of corner points of the face.
   *  \param[out] faceConn       - The corner points of the face.
   */
  void GetCornerPointsFace(unsigned short &nPointsPerFace,
                           unsigned long  faceConn[]);
};
