/*!
 * \file CPrimalGrid.cpp
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

#include "../../../include/geometry/primal_grid/CPrimalGrid.hpp"

CPrimalGrid::CPrimalGrid(bool FEM, unsigned short nNodes, unsigned short nNeighbor_Elements)
    : Nodes(new unsigned long[nNodes]), Neighbor_Elements(new long[nNeighbor_Elements]), FEM(FEM) {
  GlobalIndex_DomainElement = 0;
  for (unsigned short i = 0; i < nNeighbor_Elements; i++) Neighbor_Elements[i] = -1;
}

void CPrimalGrid::InitializeNeighbors(unsigned short val_nFaces) {
  /*--- Initialize arrays to -1/false to indicate that no neighbor is present and
        that no periodic transformation is needed to the neighbor. ---*/
  for (size_t i = 0; i < val_nFaces; i++) {
    Neighbor_Elements[i] = -1;
    PeriodIndexNeighbors[i] = -1;
  }

  for (auto i = 0; i < N_FACES_MAXIMUM; ++i) ElementOwnsFace[i] = false;
}
