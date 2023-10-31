/*!
 * \file CFVMOutput.cpp
 * \brief Main subroutines for Finite Volume Method output
 * \author T. Kattmann
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

#include "../../include/output/CFVMOutput.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"


CFVMOutput::CFVMOutput(const CConfig *config, unsigned short nDim, bool fem_output) : COutput (config, nDim, fem_output){ }

void CFVMOutput::AddCoordinates() {

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");
}


void CFVMOutput::AddCommonFVMOutputs(const CConfig *config) {

  // Mesh quality metrics
  AddVolumeOutput("ORTHOGONALITY", "Orthogonality", "MESH_QUALITY", "Orthogonality Angle (deg.)");
  AddVolumeOutput("ASPECT_RATIO",  "Aspect_Ratio",  "MESH_QUALITY", "CV Face Area Aspect Ratio");
  AddVolumeOutput("VOLUME_RATIO",  "Volume_Ratio",  "MESH_QUALITY", "CV Sub-Volume Ratio");

  AddVolumeOutput("RANK", "rank", "MPI", "Rank of the MPI-partition");

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); ++iMesh) {
    stringstream key, name;
    key << "MG_" << iMesh;
    name << "Coarse_Grid_" << iMesh;
    AddVolumeOutput(key.str(), name.str(), "MULTIGRID", "Coarse mesh");
  }

  if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    AddVolumeOutput("LINELET", "Linelet", "LINELET", "Mesh lines built for the line implicit preconditioner");
  }
}

void CFVMOutput::LoadCommonFVMOutputs(const CConfig* config, const CGeometry* geometry, unsigned long iPoint) {

  // Mesh quality metrics, computed in CPhysicalGeometry::ComputeMeshQualityStatistics.
  if (config->GetWrt_MeshQuality()) {
    SetVolumeOutputValue("ORTHOGONALITY", iPoint, geometry->Orthogonality[iPoint]);
    SetVolumeOutputValue("ASPECT_RATIO",  iPoint, geometry->Aspect_Ratio[iPoint]);
    SetVolumeOutputValue("VOLUME_RATIO",  iPoint, geometry->Volume_Ratio[iPoint]);
  }

  SetVolumeOutputValue("RANK", iPoint, rank);

  if (config->GetWrt_MultiGrid()) {
    for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); ++iMesh) {
      stringstream key;
      key << "MG_" << iMesh;
      SetVolumeOutputValue(key.str(), iPoint, geometry->CoarseGridColor(iPoint,iMesh-1));
    }
  }

  if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    SetVolumeOutputValue("LINELET", iPoint, geometry->GetLineletInfo(config).lineletColor[iPoint]);
  }
}
