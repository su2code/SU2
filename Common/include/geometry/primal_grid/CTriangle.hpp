/*!
 * \file CTriangle.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>CTriangle.cpp</i> file.
 * \author F. Palacios
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

struct CTriangleConnectivity {
  static constexpr unsigned short Faces[3][2] = {{0,1},{1,2},{2,0}};
  static constexpr unsigned short Neighbor_Nodes[3][2] = {{1,2},{2,0},{0,1}};
  static constexpr unsigned short nNodesFace[3] = {2,2,2};
  static constexpr unsigned short nNeighbor_Nodes[3] = {2,2,2};
  static constexpr unsigned short nFaces = N_POINTS_TRIANGLE;
  static constexpr unsigned short nNodes = N_FACES_TRIANGLE;
  static constexpr unsigned short nNeighbor_Elements = 3;
  static constexpr unsigned short VTK_Type = 5;
  static constexpr unsigned short maxNodesFace = 2;
};

/*!
 * \class CTriangle
 * \brief Class for triangle element definition.
 * \author F. Palacios
 */
class CTriangle final: public CPrimalGridWithConnectivity<CTriangleConnectivity> {
public:

  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st triangle point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd triangle point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th triangle point read from the grid file.
   */
  CTriangle(unsigned long val_point_0, unsigned long val_point_1,
            unsigned long val_point_2);

  /*!
   * \brief Destructor of the class.
   */
  ~CTriangle(void) override;

  /*!
   * \brief Change the orientation of an element.
   */
  inline void Change_Orientation(void) override {
    swap(Nodes[0], Nodes[2]);
  }
};
