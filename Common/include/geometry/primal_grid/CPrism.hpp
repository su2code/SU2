/*!
 * \file CPrism.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>CPrism.cpp</i> file.
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

/*! \class CPrismConnectivity
 * Defines the connectivity of a prism element.
 * See CPrimalGridWithConnectivity.
 */
struct CPrismConnectivity {
  enum { nNodes = N_POINTS_PRISM };
  enum { nFaces = N_FACES_PRISM };
  enum { maxNodesFace = N_POINTS_QUADRILATERAL };
  enum { VTK_Type = PRISM };
  static constexpr unsigned short nNodesFace[5] = {4, 4, 4, 3, 3};
  static constexpr unsigned short Faces[5][4] = {{3, 4, 1, 0}, {5, 2, 1, 4}, {2, 5, 3, 0}, {0, 1, 2, 2}, {5, 4, 3, 3}};
  static constexpr unsigned short nNeighbor_Nodes[6] = {3, 3, 3, 3, 3, 3};
  static constexpr unsigned short Neighbor_Nodes[6][3] = {{1, 2, 3}, {0, 2, 4}, {1, 0, 5},
                                                          {0, 4, 5}, {3, 5, 1}, {4, 3, 2}};
};

/*!
 * \class CPrism
 * \brief Class for prism element definition.
 * \author F. Palacios
 */
class CPrism final : public CPrimalGridWithConnectivity<CPrismConnectivity> {
 public:
  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point_0 - Index of the 1st point read from the grid file.
   * \param[in] val_point_1 - Index of the 2nd point read from the grid file.
   * \param[in] val_point_2 - Index of the 3th point read from the grid file.
   * \param[in] val_point_3 - Index of the 4th point read from the grid file.
   * \param[in] val_point_4 - Index of the 5th point read from the grid file.
   * \param[in] val_point_5 - Index of the 6th point read from the grid file.
   */
  CPrism(unsigned long val_point_0, unsigned long val_point_1, unsigned long val_point_2, unsigned long val_point_3,
         unsigned long val_point_4, unsigned long val_point_5);

  /*!
   * \brief Change the orientation of an element.
   */
  void Change_Orientation() override;
};
