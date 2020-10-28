/*!
 * \file CMeshFEM_DG.cpp
 * \brief Implementations of the member functions of CMeshFEM_DG.
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

#include "../../../include/geometry/fem_grid/CMeshFEM_DG.hpp"

CMeshFEM_DG::CMeshFEM_DG(CGeometry *geometry, CConfig *config)
  : CMeshFEM(geometry, config) {
}

void CMeshFEM_DG::CoordinatesIntegrationPoints(void) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::CoordinatesSolDOFs(void) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::CreateFaces(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::CreateStandardVolumeElements(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::InitStaticMeshMovement(CConfig              *config,
                                         const unsigned short Kind_Grid_Movement,
                                         const unsigned short iZone) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::LengthScaleVolumeElements(void) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::MetricTermsSurfaceElements(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::MetricTermsVolumeElements(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::SetGlobal_to_Local_Point(void) {
  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CMeshFEM_DG::WallFunctionPreprocessing(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
