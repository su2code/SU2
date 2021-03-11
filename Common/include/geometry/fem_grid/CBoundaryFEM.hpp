/*!
 * \file CBoundaryFEM.hpp
 * \brief Class definition for a boundary for the FEM solver.
 *        The implementations are in the <i>CBoundaryFEM.cpp</i> file.
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

#include "CSurfaceElementFEM.hpp"
#include "../../wall_model.hpp"

using namespace std;

/*!
 * \class CBoundaryFEM
 * \brief Class to store a boundary for the FEM solver.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
 */
struct CBoundaryFEM {
  string markerTag;  /*!< \brief Marker tag of this boundary. */

  bool periodicBoundary = false;     /*!< \brief Whether or not this boundary is a periodic boundary. */
  bool haloInfoNeededForBC = false;  /*!< \brief Whether or not information of halo elements
                                                 is needed to impose the boundary conditions. */

  vector<unsigned long> nSurfElem; /*!< \brief Number of surface elements per time level,
                                               cumulative storage format. */

  vector<CSurfaceElementFEM> surfElem; /*!< \brief Vector of the local surface elements. */

  CWallModel *wallModel = nullptr;     /*!< \brief Wall model for LES. */

  ~CBoundaryFEM(void) { delete wallModel; }
};
