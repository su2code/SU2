/*!
 * \file CLine.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>CLine.cpp</i> file.
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

/*! \class CLineConnectivity
 * Defines the connectivity of a line element.
 * See CPrimalGridWithConnectivity.
 */
struct CLineConnectivity {
  enum { nNodes = N_POINTS_LINE };
  enum { nFaces = N_FACES_LINE };
  enum { maxNodesFace = N_POINTS_LINE };
  enum { VTK_Type = LINE };
  static constexpr unsigned short nNodesFace[1] = {2};
  static constexpr unsigned short Faces[1][2] = {{0, 1}};
  static constexpr unsigned short nNeighbor_Nodes[2] = {1, 1};
  static constexpr unsigned short Neighbor_Nodes[2][1] = {{1}, {0}};
};

/*!
 * \class CLine
 * \brief Class for line element definition.
 * \author F. Palacios
 */
class CLine final : public CPrimalGridWithConnectivity<CLineConnectivity> {
 public:
  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st triangle point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd triangle point read from the grid file.
   */
  CLine(unsigned long val_point_0, unsigned long val_point_1);

  /*!
   * \brief Change the orientation of an element.
   */
  inline void Change_Orientation() override { std::swap(Nodes[0], Nodes[1]); }
};
