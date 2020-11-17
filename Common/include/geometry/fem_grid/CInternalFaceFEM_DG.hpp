/*!
 * \file CInternalFaceFEM_DG.hpp
 * \brief Class for an internal face element for the DG-FEM solver.
 *        The implementations are in the <i>CInternalFaceFEM_DG.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "../../fem/CFEMStandardInternalFaceGrid.hpp"
#include "../../fem/CFEMStandardInternalFaceSol.hpp"

using namespace std;

/*!
 * \class CInternalFaceElementFEM
 * \brief Class to store an internal face for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
class CInternalFaceFEM_DG final {
public:
  unsigned long  elemID0;         /*!< \brief Element ID adjacent to side 0 of the face. */
  unsigned long  elemID1;         /*!< \brief Element ID adjacent to side 1 of the face. */

  ColMajorMatrix<su2double> metricNormalsFace;     /*!< \brief The normals in the integration points of the face.
                                                               The normals point from side 0 to side 1. */
  ColMajorMatrix<su2double> metricCoorDerivFace0;  /*!< \brief The terms drdx, dsdx, etc. of side 0 in the
                                                               integration points of the face. */
  ColMajorMatrix<su2double> metricCoorDerivFace1;  /*!< \brief The terms drdx, dsdx, etc. of side 1 in the
                                                               integration points of the face. */

  ColMajorMatrix<su2double> coorIntegrationPoints;  /*!< \brief Coordinates of the integration points of this face. */
  ColMajorMatrix<su2double> gridVelocities;         /*!< \brief Grid velocities of the integration points of this face. */
  ColMajorMatrix<su2double> wallDistance;           /*!< \brief The wall distance to the viscous walls for
                                                                the integration points of this face. */

  CFEMStandardInternalFaceGrid *standardElemGrid = nullptr; /*!< \brief Pointer to the standard element for the grid. */
  CFEMStandardInternalFaceSol  *standardElemFlow = nullptr; /*!< \brief Pointer to the standard element for the
                                                                        standard flow solution variables. */
  CFEMStandardInternalFaceSol  *standardElemP    = nullptr; /*!< \brief Pointer to the standard element for the
                                                                        pressure for an incompressible flow. */
};
