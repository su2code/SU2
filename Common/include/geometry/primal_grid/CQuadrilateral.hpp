/*!
 * \file CQuadrilateral.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>CQuadrilateral.cpp</i> file.
 * \author F. Palacios
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

#include "CPrimalGrid.hpp"

/*! \class CQuadrilateralConnectivity
 * Defines the connectivity of a quadrilateral element.
 * See CPrimalGridWithConnectivity.
 */
struct CQuadrilateralConnectivity {
  enum { nNodes = N_POINTS_QUADRILATERAL };
  enum { nFaces = N_FACES_QUADRILATERAL };
  enum { maxNodesFace = N_POINTS_LINE };
  enum { VTK_Type = QUADRILATERAL };
  static constexpr unsigned short nNodesFace[4] = {2, 2, 2, 2};
  static constexpr unsigned short Faces[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  static constexpr unsigned short nNeighbor_Nodes[4] = {2, 2, 2, 2};
  static constexpr unsigned short Neighbor_Nodes[4][2] = {{1, 3}, {2, 0}, {3, 1}, {0, 2}};
};

/*!
 * \class CQuadrilateral
 * \brief Class for quadrilateral element definition.
 * \author F. Palacios
 */
class CQuadrilateral final : public CPrimalGridWithConnectivity<CQuadrilateralConnectivity> {
 public:
  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th point read from the grid file.
   * \param[in] val_point_3 - Index of the 4th point read from the grid file.
   */
  CQuadrilateral(unsigned long val_point_0, unsigned long val_point_1, unsigned long val_point_2,
                 unsigned long val_point_3);

  /*!
   * \brief Change the orientation of an element.
   */
  void Change_Orientation() override;
};
