/*!
 * \file CFVMOutput.cpp
 * \brief Main subroutines for Finite Volume Method output
 * \author T. Kattmann
 * \version 7.5.0 "Blackbird"
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

#include "../../include/output/CFVMOutput.hpp"
#include <cstdint>
#include <limits>
#include <unordered_set>
#include "../../include/solvers/CSolver.hpp"
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

  AddVolumeOutput("LINELET", "Linelet", "LINELET", "Mesh lines built for the line implicit preconditioner");
}

void CFVMOutput::LoadCommonFVMOutputs(const CConfig* config, const CGeometry* geometry, const CSolver* solver,
                                      unsigned long iPoint) {

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
    /*--- Lazy build. ---*/
    if (lineletPointColor.empty()) {
      /*--- First color the lines based on connections between then. ---*/
      const auto& linelets = solver->Jacobian.GetLinelets();
      const auto nLine = linelets.size();

      /*--- Adjacency between lines, computed from point neighbors. ---*/
      std::vector<std::vector<unsigned long>> adjacency(nLine);
      for (auto iLine = 0ul; iLine < nLine; ++iLine) {
        std::unordered_set<unsigned long> neighbors;
        for (const auto iPoint : linelets[iLine]) {
          neighbors.insert(iPoint);
          for (const auto jPoint : geometry->nodes->GetPoints(iPoint)) {
            neighbors.insert(jPoint);
          }
        }
        adjacency[iLine].reserve(neighbors.size());
        for (const auto iPoint : neighbors) {
          adjacency[iLine].push_back(iPoint);
        }
      }

      std::vector<uint8_t> lineletColor;
      const unsigned long nColors = colorSparsePattern<uint8_t, std::numeric_limits<uint8_t>::max()>(
          CCompressedSparsePatternUL(adjacency), 1, true, &lineletColor).getOuterSize();

      /*--- Offset colors to avoid coloring across ranks. ---*/
      std::vector<unsigned long> allNColors(size);
      SU2_MPI::Allgather(&nColors, 1, MPI_UNSIGNED_LONG, allNColors.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());
      unsigned long offset = 0;
      for (int i = 0; i < rank; ++i) offset += allNColors[i];

      /*--- Finally, transfer colors to points. ---*/
      lineletPointColor.resize(geometry->GetnPoint(), 0);
      for (auto iLine = 0ul; iLine < nLine; ++iLine) {
        for (const auto iPoint : linelets[iLine]) {
          lineletPointColor[iPoint] = 1 + offset + lineletColor[iLine];
        }
      }
    }

    SetVolumeOutputValue("LINELET", iPoint, lineletPointColor[iPoint]);
  }
}
