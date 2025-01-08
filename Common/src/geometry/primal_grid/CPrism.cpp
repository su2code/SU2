/*!
 * \file CPrism.cpp
 * \brief Main classes for defining the primal grid elements
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

#include "../../../include/geometry/primal_grid/CPrism.hpp"
#include "../../../include/option_structure.hpp"

constexpr unsigned short CPrismConnectivity::nNodesFace[5];
constexpr unsigned short CPrismConnectivity::Faces[5][4];
constexpr unsigned short CPrismConnectivity::nNeighbor_Nodes[6];
constexpr unsigned short CPrismConnectivity::Neighbor_Nodes[6][3];

CPrism::CPrism(unsigned long val_point_0, unsigned long val_point_1, unsigned long val_point_2,
               unsigned long val_point_3, unsigned long val_point_4, unsigned long val_point_5)
    : CPrimalGridWithConnectivity<CPrismConnectivity>(false) {
  /*--- Define face structure of the element ---*/
  Nodes[0] = val_point_0;
  Nodes[1] = val_point_1;
  Nodes[2] = val_point_2;
  Nodes[3] = val_point_3;
  Nodes[4] = val_point_4;
  Nodes[5] = val_point_5;
}

void CPrism::Change_Orientation() {
  std::swap(Nodes[0], Nodes[1]);
  std::swap(Nodes[3], Nodes[4]);
}
