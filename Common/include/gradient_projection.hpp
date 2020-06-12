/*!
 * \file gradient_projection.hpp
 * \brief Include files and headers of the functions to carry out
 *        projections between mesh coordinate sensitivities and design values
 *        similar to the functions in SU2_DOT
 * \author T.Dick
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "option_structure.hpp"
#include "CConfig.hpp"
#include "grid_movement_structure.hpp"
#include "eigen_structure.hpp"

/*!
 * \brief GetParameterizationJacobianForward
 * \param geometry
 * \param config
 * \param surface_movement
 * \param Jacobian of the parameterization
 */
void GetParameterizationJacobianForward(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double *Jacobian);

/*!
 * \brief GetParameterizationJacobianReverse
 * \param geometry
 * \param config
 * \param surface_movement
 * \param Jacobian of the parameterization
 */
void GetParameterizationJacobianReverse(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double *Jacobian);

/*!
 * \brief GetParameterizationJacobianPreaccumulation
 * \param geometry
 * \param config
 * \param surface_movement
 * \param Jacobian of the parameterization
 */
void GetParameterizationJacobianPreaccumulation(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double *Jacobian);

/*!
 * \brief Cast2Eigenmatrix
 * \param Jacobian of the parameterization
 */
MatrixType Cast2Eigenmatrix(CGeometry *geometry, CConfig *config, su2double *Jacobian);
