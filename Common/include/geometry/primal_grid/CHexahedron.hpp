/*!
 * \file CHexahedron.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>CHexahedron.cpp</i> file.
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

/*! \class CHexahedronConnectivity
 * Defines the connectivity of a hexahedron element.
 * See CPrimalGridWithConnectivity.
 */
struct CHexahedronConnectivity {
  enum { nNodes = N_POINTS_HEXAHEDRON };
  enum { nFaces = N_FACES_HEXAHEDRON };
  enum { maxNodesFace = N_POINTS_QUADRILATERAL };
  enum { VTK_Type = HEXAHEDRON };
  static constexpr unsigned short nNodesFace[6] = {4, 4, 4, 4, 4, 4};
  static constexpr unsigned short Faces[6][4] = {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6},
                                                 {3, 0, 4, 7}, {0, 3, 2, 1}, {4, 5, 6, 7}};
  static constexpr unsigned short nNeighbor_Nodes[8] = {3, 3, 3, 3, 3, 3, 3, 3};
  static constexpr unsigned short Neighbor_Nodes[8][3] = {{1, 3, 4}, {0, 2, 5}, {1, 3, 6}, {0, 2, 7},
                                                          {0, 5, 7}, {4, 6, 1}, {2, 5, 7}, {4, 3, 6}};
};

/*!
 * \class CHexahedron
 * \brief Class for hexahedron element definition.
 * \author F. Palacios
 */
class CHexahedron : public CPrimalGridWithConnectivity<CHexahedronConnectivity> {
 public:
  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th point read from the grid file.
   * \param[in] val_point_3 - Index of the 4th point read from the grid file.
   * \param[in] val_point_4 - Index of the 5td point read from the grid file.
   * \param[in] val_point_5 - Index of the 6th point read from the grid file.
   * \param[in] val_point_6 - Index of the 7th point read from the grid file.
   * \param[in] val_point_7 - Index of the 8th point read from the grid file.
   */
  CHexahedron(unsigned long val_point_0, unsigned long val_point_1, unsigned long val_point_2,
              unsigned long val_point_3, unsigned long val_point_4, unsigned long val_point_5,
              unsigned long val_point_6, unsigned long val_point_7);

  /*!
   * \brief Change the orientation of an element.
   */
  void Change_Orientation() override;
};
