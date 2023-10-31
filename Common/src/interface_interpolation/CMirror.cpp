/*!
 * \file CMirror.cpp
 * \brief Implementation of mirror interpolation (conservative approach in FSI problems).
 * \author P. Gomes
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

#include "../../include/interface_interpolation/CMirror.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/toolboxes/printing_toolbox.hpp"

CMirror::CMirror(CGeometry**** geometry_container, const CConfig* const* config, const CInterpolator* interpolator,
                 unsigned int iZone, unsigned int jZone)
    : CInterpolator(geometry_container, config, iZone, jZone), transpInterpolator(interpolator) {
  using PrintingToolbox::to_string;
  if (jZone < iZone) {
    SU2_MPI::Error(string("The order of the zones does not allow conservative interpolation to be setup.\n"
                          "Swap zones ") +
                       to_string(iZone) + string(" and ") + to_string(jZone) + string("."),
                   CURRENT_FUNCTION);
  }
  SetTransferCoeff(config);
}

void CMirror::SetTransferCoeff(const CConfig* const* config) {
  const int nProcessor = size;

  vector<unsigned long> allNumVertexTarget(nProcessor);
  vector<unsigned long> allNumVertexDonor(nProcessor);
  vector<unsigned long> allNumNodeDonor(nProcessor);

  /*--- The target vertex information of the transpose interpolator. ---*/
  const auto& donorVertices = transpInterpolator->targetVertices;

  targetVertices.resize(config[targetZone]->GetnMarker_All());

  /*--- Number of markers on the interface ---*/
  const auto nMarkerInt = (config[targetZone]->GetMarker_n_ZoneInterface()) / 2;

  /*--- For the number of markers on the interface... ---*/
  for (unsigned short iMarkerInt = 0; iMarkerInt < nMarkerInt; iMarkerInt++) {
    /* High level procedure:
     * - Gather the interpolation matrix of the donor geometry;
     * - Set the interpolation matrix of the target as the transpose.
     */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    const auto markDonor = config[donorZone]->FindInterfaceMarker(iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    const auto markTarget = config[targetZone]->FindInterfaceMarker(iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if (!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    unsigned long nVertexDonor = 0, nVertexTarget = 0;
    if (markDonor != -1) nVertexDonor = donor_geometry->GetnVertex(markDonor);
    if (markTarget != -1) nVertexTarget = target_geometry->GetnVertex(markTarget);

    /*--- Count the number of donor nodes on the donor geometry. ---*/
    unsigned long nVertexDonorLocal = 0;
    unsigned long nNodeDonorLocal = 0;
    for (auto iVertex = 0ul; iVertex < nVertexDonor; iVertex++) {
      const auto iPoint = donor_geometry->vertex[markDonor][iVertex]->GetNode();
      if (donor_geometry->nodes->GetDomain(iPoint)) {
        nNodeDonorLocal += donorVertices[markDonor][iVertex].nDonor();
        nVertexDonorLocal++;
      }
    }

    /*--- Communicate vertex and donor node counts. ---*/
    SU2_MPI::Allgather(&nVertexTarget, 1, MPI_UNSIGNED_LONG, allNumVertexTarget.data(), 1, MPI_UNSIGNED_LONG,
                       SU2_MPI::GetComm());
    SU2_MPI::Allgather(&nVertexDonorLocal, 1, MPI_UNSIGNED_LONG, allNumVertexDonor.data(), 1, MPI_UNSIGNED_LONG,
                       SU2_MPI::GetComm());
    SU2_MPI::Allgather(&nNodeDonorLocal, 1, MPI_UNSIGNED_LONG, allNumNodeDonor.data(), 1, MPI_UNSIGNED_LONG,
                       SU2_MPI::GetComm());

    /*--- Copy donor interpolation matrix (triplet format). ---*/
    vector<long> sendGlobalIndex(nNodeDonorLocal);
    vector<long> sendDonorIndex(nNodeDonorLocal);
    vector<su2double> sendDonorCoeff(nNodeDonorLocal);

    for (auto iVertex = 0ul, iDonor = 0ul; iVertex < nVertexDonor; ++iVertex) {
      auto& donor_vertex = donorVertices[markDonor][iVertex];
      const auto iPoint = donor_geometry->vertex[markDonor][iVertex]->GetNode();

      if (!donor_geometry->nodes->GetDomain(iPoint)) continue;

      const auto nDonor = donor_vertex.nDonor();
      const auto donorGlobalIndex = donor_geometry->nodes->GetGlobalIndex(iPoint);

      for (auto i = 0u; i < nDonor; ++i) {
        sendGlobalIndex[iDonor] = donorGlobalIndex;
        sendDonorIndex[iDonor] = donor_vertex.globalPoint[i];
        sendDonorCoeff[iDonor] = donor_vertex.coefficient[i];
        ++iDonor;
      }
    }

    /*--- Sort the matrix by donor index, effectively transposing the triplets. ---*/
    vector<int> order(nNodeDonorLocal);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&sendDonorIndex](int i, int j) { return sendDonorIndex[i] < sendDonorIndex[j]; });
    for (int i = 0; i < int(nNodeDonorLocal); ++i) {
      int j = order[i];
      while (j < i) j = order[j];
      if (i == j) continue;
      swap(sendGlobalIndex[i], sendGlobalIndex[j]);
      swap(sendDonorIndex[i], sendDonorIndex[j]);
      swap(sendDonorCoeff[i], sendDonorCoeff[j]);
    }
    vector<int>().swap(order);  // no longer needed

    /*--- Communicate donor interpolation matrix and info. We only gather the
     *    matrix in ranks that need it, i.e. have target vertices, to avoid
     *    potentially huge memory usage when running in parallel. This is an
     *    in-place Gatherv (done manually due to AD issues in the past). ---*/

    vector<int> iSendProcessor;
    for (int iProcessor = 0; iProcessor < nProcessor; ++iProcessor)
      if (allNumVertexDonor[iProcessor] != 0) iSendProcessor.push_back(iProcessor);

    const int nSend = iSendProcessor.size();

    vector<long*> GlobalIndex(nSend, nullptr), DonorIndex(nSend, nullptr);
    vector<su2double*> DonorCoeff(nSend, nullptr);

    /*--- For each "target processor" that needs the interpolation matrix. ---*/
    for (int iProcessor = 0; iProcessor < nProcessor; ++iProcessor) {
      if (allNumVertexTarget[iProcessor] == 0) continue;

      /*--- For each "donor processor that contain a part of the matrix. ---*/
      for (int iSend = 0; iSend < nSend; ++iSend) {
        const auto jProcessor = iSendProcessor[iSend];
        const auto numCoeff = allNumNodeDonor[jProcessor];

        if ((rank == iProcessor) && (rank == jProcessor)) {
          /*--- "Self" communication. ---*/
          GlobalIndex[iSend] = sendGlobalIndex.data();
          DonorIndex[iSend] = sendDonorIndex.data();
          DonorCoeff[iSend] = sendDonorCoeff.data();
        } else if (rank == iProcessor) {
          /*--- "I'm" the target, allocate and receive. ---*/
          GlobalIndex[iSend] = new long[numCoeff];
          DonorIndex[iSend] = new long[numCoeff];
          DonorCoeff[iSend] = new su2double[numCoeff];
          SU2_MPI::Recv(GlobalIndex[iSend], numCoeff, MPI_LONG, jProcessor, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
          SU2_MPI::Recv(DonorIndex[iSend], numCoeff, MPI_LONG, jProcessor, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
          SU2_MPI::Recv(DonorCoeff[iSend], numCoeff, MPI_DOUBLE, jProcessor, 0, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
        } else if (rank == jProcessor) {
          /*--- "I'm" the donor, send. ---*/
          SU2_MPI::Send(sendGlobalIndex.data(), numCoeff, MPI_LONG, iProcessor, 0, SU2_MPI::GetComm());
          SU2_MPI::Send(sendDonorIndex.data(), numCoeff, MPI_LONG, iProcessor, 0, SU2_MPI::GetComm());
          SU2_MPI::Send(sendDonorCoeff.data(), numCoeff, MPI_DOUBLE, iProcessor, 0, SU2_MPI::GetComm());
        }
      }
    }

    if (nVertexTarget) targetVertices[markTarget].resize(nVertexTarget);

    /*--- Loop over the vertices on the target marker, define one row of the transpose matrix. ---*/

    SU2_OMP_PARALLEL {
      SU2_OMP_FOR_DYN(roundUpDiv(nVertexTarget, 2 * omp_get_max_threads()))
      for (auto iVertex = 0ul; iVertex < nVertexTarget; ++iVertex) {
        auto& target_vertex = targetVertices[markTarget][iVertex];
        const auto iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();

        if (!target_geometry->nodes->GetDomain(iPoint)) continue;

        /*--- Any point of the donor geometry, that has this target point as a donor, becomes a donor. ---*/
        const long targetGlobalIndex = target_geometry->nodes->GetGlobalIndex(iPoint);

        /*--- Count donors and safe the binary search results (this is why we sorted the matrix). ---*/
        auto nDonor = 0ul;
        vector<pair<long*, long*> > ranges(nSend);
        for (int iSend = 0; iSend < nSend; ++iSend) {
          const auto iProcessor = iSendProcessor[iSend];
          const auto numCoeff = allNumNodeDonor[iProcessor];
          auto p = equal_range(DonorIndex[iSend], DonorIndex[iSend] + numCoeff, targetGlobalIndex);
          nDonor += (p.second - p.first);
          ranges[iSend] = p;
        }

        target_vertex.resize(nDonor);

        /*--- Use the search results to set the interpolation coefficients. ---*/
        for (int iSend = 0, iDonor = 0; iSend < nSend; ++iSend) {
          const auto iProcessor = iSendProcessor[iSend];

          const auto first = ranges[iSend].first - DonorIndex[iSend];
          const auto last = ranges[iSend].second - DonorIndex[iSend];

          for (auto iCoeff = first; iCoeff < last; ++iCoeff) {
            target_vertex.processor[iDonor] = iProcessor;
            target_vertex.coefficient[iDonor] = DonorCoeff[iSend][iCoeff];
            target_vertex.globalPoint[iDonor] = GlobalIndex[iSend][iCoeff];
            ++iDonor;
          }
        }
      }
      END_SU2_OMP_FOR
    }
    END_SU2_OMP_PARALLEL

    /*--- Free the heap allocations. ---*/
    for (auto ptr : GlobalIndex)
      if (ptr != sendGlobalIndex.data()) delete[] ptr;
    for (auto ptr : DonorIndex)
      if (ptr != sendDonorIndex.data()) delete[] ptr;
    for (auto ptr : DonorCoeff)
      if (ptr != sendDonorCoeff.data()) delete[] ptr;

  }  // end marker loop
}
