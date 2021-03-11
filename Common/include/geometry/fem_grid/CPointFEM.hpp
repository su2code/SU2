/*!
 * \file CPointFEM.hpp
 * \brief Class definition for a grid point for the FEM solver.
 *        The implementations are in the <i>CPointFEM.cpp</i> file.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

#include "../../parallelization/mpi_structure.hpp"

using namespace std;

/*!
 * \class CPointFEM
 * \brief Class to define a point for the FEM solver.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
struct CPointFEM {
  unsigned long globalID;    /*!< \brief The global ID of this point in the grid. */
  short periodIndexToDonor;  /*!< \brief The index of the periodic transformation to the donor
                                  element. Only for halo elements. A -1 indicates no periodic
                                  transformation. */
  su2double coor[3] = {0.0}; /*!< \brief Array with the coordinates of the node. */

  /*!
   * \brief Less than operator of the class. Needed for the sorting.
   */
  bool operator<(const CPointFEM &other) const;

  /*!
   * \brief Equal operator of the class. Needed for the removal of double entities.
   */
  bool operator==(const CPointFEM &other) const;
};
