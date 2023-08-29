/*!
 * \file CNearestNeighbor.cpp
 * \brief Implementation of nearest neighbor interpolation.
 * \author H. Kline
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

#include "../../include/interface_interpolation/CNearestNeighbor.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

CNearestNeighbor::CNearestNeighbor(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone,
                                   unsigned int jZone)
    : CInterpolator(geometry_container, config, iZone, jZone) {
  SetTransferCoeff(config);
}

void CNearestNeighbor::PrintStatistics() const {
  if (rank != MASTER_NODE) return;
  cout << "  Avg/max distance to closest donor point: " << AvgDistance << "/" << MaxDistance << endl;
}

void CNearestNeighbor::SetTransferCoeff(const CConfig* const* config) {
  /*--- Desired number of donor points. ---*/
  const auto nDonor = max<unsigned long>(config[donorZone]->GetNumNearestNeighbors(), 1);

  /*--- Epsilon used to avoid division by zero. ---*/
  const su2double eps = numeric_limits<passivedouble>::epsilon();

  const int nProcessor = size;
  const auto nMarkerInt = config[donorZone]->GetMarker_n_ZoneInterface() / 2;
  const auto nDim = donor_geometry->GetnDim();

  Buffer_Receive_nVertex_Donor = new unsigned long[nProcessor];

  targetVertices.resize(config[targetZone]->GetnMarker_All());

  vector<vector<DonorInfo> > DonorInfoVec(omp_get_max_threads());

  /*--- Cycle over nMarkersInt interface to determine communication pattern. ---*/

  AvgDistance = MaxDistance = 0.0;
  unsigned long totalTargetPoints = 0;

  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; iMarkerInt++) {
    /*--- On the donor side: find the tag of the boundary sharing the interface. ---*/
    const auto markDonor = config[donorZone]->FindInterfaceMarker(iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface. ---*/
    const auto markTarget = config[targetZone]->FindInterfaceMarker(iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step. ---*/
    if (!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    unsigned long nVertexDonor = 0, nVertexTarget = 0;
    if (markDonor != -1) nVertexDonor = donor_geometry->GetnVertex(markDonor);
    if (markTarget != -1) nVertexTarget = target_geometry->GetnVertex(markTarget);

    /*--- Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor. ---*/
    Determine_ArraySize(markDonor, markTarget, nVertexDonor, nDim);
    if (nVertexTarget) targetVertices[markTarget].resize(nVertexTarget);

    const auto nPossibleDonor =
        accumulate(Buffer_Receive_nVertex_Donor, Buffer_Receive_nVertex_Donor + nProcessor, 0ul);

    Buffer_Send_Coord.resize(MaxLocalVertex_Donor, nDim);
    Buffer_Send_GlobalPoint.resize(MaxLocalVertex_Donor);
    Buffer_Receive_Coord.resize(nProcessor * MaxLocalVertex_Donor, nDim);
    Buffer_Receive_GlobalPoint.resize(nProcessor * MaxLocalVertex_Donor);

    /*--- Collect coordinates and global point indices. ---*/
    Collect_VertexInfo(markDonor, markTarget, nVertexDonor, nDim);

    /*--- Find the closest donor points to each target. ---*/
    SU2_OMP_PARALLEL {
      /*--- Working array for this thread. ---*/
      auto& donorInfo = DonorInfoVec[omp_get_thread_num()];
      donorInfo.resize(nPossibleDonor);

      su2double avgDist = 0.0, maxDist = 0.0;
      unsigned long numTarget = 0;

      SU2_OMP_FOR_DYN(roundUpDiv(nVertexTarget, 2 * omp_get_max_threads()))
      for (auto iVertexTarget = 0ul; iVertexTarget < nVertexTarget; iVertexTarget++) {
        auto& target_vertex = targetVertices[markTarget][iVertexTarget];
        const auto Point_Target = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();

        if (!target_geometry->nodes->GetDomain(Point_Target)) continue;

        /*--- Coordinates of the target point. ---*/
        const su2double* Coord_i = target_geometry->nodes->GetCoord(Point_Target);

        /*--- Compute all distances. ---*/
        for (int iProcessor = 0, iDonor = 0; iProcessor < nProcessor; ++iProcessor) {
          for (auto jVertex = 0ul; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; ++jVertex) {
            const auto idx = iProcessor * MaxLocalVertex_Donor + jVertex;
            const auto pGlobalPoint = Buffer_Receive_GlobalPoint[idx];
            const su2double* Coord_j = Buffer_Receive_Coord[idx];
            const auto dist2 = GeometryToolbox::SquaredDistance(nDim, Coord_i, Coord_j);

            donorInfo[iDonor++] = DonorInfo(dist2, pGlobalPoint, iProcessor);
          }
        }

        /*--- Find k closest points. ---*/
        partial_sort(donorInfo.begin(), donorInfo.begin() + nDonor, donorInfo.end(),
                     [](const DonorInfo& a, const DonorInfo& b) {
                       /*--- Global index is used as tie-breaker to make sorted order independent of initial. ---*/
                       return (a.dist != b.dist) ? (a.dist < b.dist) : (a.pidx < b.pidx);
                     });

        /*--- Update stats. ---*/
        numTarget += 1;
        su2double d = sqrt(donorInfo[0].dist);
        avgDist += d;
        maxDist = max(maxDist, d);

        /*--- Compute interpolation numerators and denominator. ---*/
        su2double denom = 0.0;
        for (auto iDonor = 0ul; iDonor < nDonor; ++iDonor) {
          donorInfo[iDonor].dist = 1.0 / (donorInfo[iDonor].dist + eps);
          denom += donorInfo[iDonor].dist;
        }

        /*--- Set interpolation coefficients. ---*/
        target_vertex.resize(nDonor);

        for (auto iDonor = 0ul; iDonor < nDonor; ++iDonor) {
          target_vertex.globalPoint[iDonor] = donorInfo[iDonor].pidx;
          target_vertex.processor[iDonor] = donorInfo[iDonor].proc;
          target_vertex.coefficient[iDonor] = donorInfo[iDonor].dist / denom;
        }
      }
      END_SU2_OMP_FOR
      SU2_OMP_CRITICAL {
        totalTargetPoints += numTarget;
        AvgDistance += avgDist;
        MaxDistance = max(MaxDistance, maxDist);
      }
      END_SU2_OMP_CRITICAL
    }
    END_SU2_OMP_PARALLEL
  }

  delete[] Buffer_Receive_nVertex_Donor;

  unsigned long tmp = totalTargetPoints;
  SU2_MPI::Allreduce(&tmp, &totalTargetPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  su2double tmp1 = AvgDistance, tmp2 = MaxDistance;
  SU2_MPI::Allreduce(&tmp1, &AvgDistance, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&tmp2, &MaxDistance, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
  AvgDistance /= totalTargetPoints;
}
