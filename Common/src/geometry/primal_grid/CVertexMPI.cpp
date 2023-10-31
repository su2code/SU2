/*!
 * \file CVertexMPI.cpp
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

#include "../../../include/geometry/primal_grid/CVertexMPI.hpp"

constexpr unsigned short CVertexMPIConnectivity::nNodesFace[1];
constexpr unsigned short CVertexMPIConnectivity::Faces[1][1];
constexpr unsigned short CVertexMPIConnectivity::nNeighbor_Nodes[1];
constexpr unsigned short CVertexMPIConnectivity::Neighbor_Nodes[1][1];

CVertexMPI::CVertexMPI(unsigned long val_point) : CPrimalGridWithConnectivity<CVertexMPIConnectivity>(false) {
  /*--- Define face structure of the element ---*/
  Nodes[0] = val_point;

  /*--- By default, no rotation in the solution ---*/
  Rotation_Type = 0;
}
