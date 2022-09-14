/*!
 * \file SU2_DOT.hpp
 * \brief Headers of the main subroutines of the code SU2_DOT.
 *        The subroutines and functions are in the <i>SU2_DOT.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#define ENABLE_MAPS
#include "../../Common/include/CConfig.hpp"
#undef ENABLE_MAPS

#include "../../Common/include/parallelization/mpi_structure.hpp"
#include "../../Common/include/parallelization/omp_structure.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/fem/fem_geometry_structure.hpp"
#include "../../SU2_CFD/include/output/CBaselineOutput.hpp"
#include "../../SU2_CFD/include/solvers/CBaselineSolver.hpp"

#include "../../SU2_CFD/include/solvers/CGradientSmoothingSolver.hpp"
#include "../../SU2_CFD/include/numerics/CGradSmoothing.hpp"

using namespace std;


/*!
 * \brief Projection of the surface sensitivity using finite differences (FD).
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement class of the problem.
 * \param[in] Gradient_file - Output file to store the gradient data.
 */

void SetProjection_FD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double **Gradient);

/*!
 * \brief Projection of the surface sensitivity using algorithmic differentiation (AD).
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] surface_movement - Surface movement class of the problem.
 * \param[in] Gradient_file - Output file to store the gradient data.
 */

void SetProjection_AD(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double **Gradient);

/*!
 * \brief Prints the gradient information to a file.
 * \param[in] Gradient - The gradient data.
 * \param[in] config - Definition of the particular problem.
 * \param[in] Gradient_file - Output file to store the gradient data.
 */

void OutputGradient(su2double** Gradient, CConfig* config, ofstream& Gradient_file);

/*!
 * \brief Write the sensitivity (including mesh sensitivity) computed with the discrete adjoint method
 *  on the surface and in the volume to a file.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] val_nZone - Number of Zones.
 */

void SetSensitivity_Files(CGeometry ***geometry, CConfig **config, unsigned short val_nZone);

/*!
 * \brief Treatment of derivatives with the Sobolev smoothing solver.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] grid_movement - Volumetric movement class of the problem.
 */

void DerivativeTreatment_MeshSensitivity(CGeometry *geometry, CConfig *config, CVolumetricMovement *grid_movement);

/*!
 * \brief Treatment of derivatives with the Sobolev smoothing solver.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] grid_movement - Volumetric movement class of the problem.
 * \param[in] surface_movement - Surface movement class of the problem.
 * \param[in] Gradient - Output array to store the gradient data.
 */

void DerivativeTreatment_Gradient(CGeometry *geometry, CConfig *config, CVolumetricMovement *grid_movement, CSurfaceMovement *surface_movement, su2double **Gradient);
