/*!
 * \file CMirror.cpp
 * \brief Implementation of mirror interpolation (conservative approach in FSI problems).
 * \author H. Kline, P. Gomes
 * \version 7.0.3 "Blackbird"
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

#include "../../include/interface_interpolation/CMirror.hpp"
#include "../../include/CConfig.hpp"
#include "../../include/geometry/CGeometry.hpp"
#include "../../include/toolboxes/printing_toolbox.hpp"


CMirror::CMirror(CGeometry ****geometry_container, const CConfig* const* config,  unsigned int iZone,
                 unsigned int jZone) : CInterpolator(geometry_container, config, iZone, jZone) {
  using PrintingToolbox::to_string;
  if (jZone < iZone) {
    SU2_MPI::Error(string("The order of the zones does not allow conservative interpolation to be setup.\n"
      "Swap zones ") + to_string(iZone) + string(" and ") + to_string(jZone) + string("."),CURRENT_FUNCTION);
  }
  SetTransferCoeff(config);
}

void CMirror::SetTransferCoeff(const CConfig* const* config) {

  const int nProcessor = size;

  vector<unsigned long> allNumVertexTarget(nProcessor);
  vector<unsigned long> allNumVertexDonor(nProcessor);
  vector<unsigned long> allNumNodeDonor(nProcessor);

  /*--- Number of markers on the interface ---*/
  const auto nMarkerInt = (config[targetZone]->GetMarker_n_ZoneInterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (unsigned short iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

   /* High level procedure:
    * - Gather the interpolation matrix of the donor geometry;
    * - Set the interpolation matrix of the target as the transpose.
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    const auto markDonor = FindInterfaceMarker(config[donorZone], iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    const auto markTarget = FindInterfaceMarker(config[targetZone], iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if (!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    unsigned long nVertexDonor = 0, nVertexTarget = 0;
    if (markDonor != -1) nVertexDonor = donor_geometry->GetnVertex( markDonor );
    if (markTarget != -1) nVertexTarget = target_geometry->GetnVertex( markTarget );

    /*--- Count the number of donor nodes on the donor geometry. ---*/
    unsigned long nVertexDonorLocal = 0;
    unsigned long nNodeDonorLocal = 0;
    for (auto iVertex = 0ul; iVertex < nVertexDonor; iVertex++) {

      auto donor_vertex = donor_geometry->vertex[markDonor][iVertex];

      if (donor_geometry->node[donor_vertex->GetNode()]->GetDomain()) {
        nNodeDonorLocal += donor_vertex->GetnDonorPoints();
        nVertexDonorLocal++;
      }
    }

    /*--- Communicate vertex and donor node counts. ---*/
    SU2_MPI::Allgather(&nVertexTarget, 1, MPI_UNSIGNED_LONG,
                       allNumVertexTarget.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(&nVertexDonorLocal, 1, MPI_UNSIGNED_LONG,
                       allNumVertexDonor.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(&nNodeDonorLocal, 1, MPI_UNSIGNED_LONG,
                       allNumNodeDonor.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    /*--- Copy donor interpolation matrix and info. ---*/
    const auto sizeRowPtr = nVertexDonorLocal + 1;
    vector<long> sendGlobalIndex(nVertexDonorLocal);
    vector<long> sendDonorRowPtr(sizeRowPtr);
    vector<long> sendDonorIndex(nNodeDonorLocal);
    vector<su2double> sendDonorCoeff(nNodeDonorLocal);

    sendDonorRowPtr[0] = 0;
    for (auto iVertex = 0ul, iVertexLocal = 0ul, iDonor = 0ul; iVertex < nVertexDonor; ++iVertex) {

      auto donor_vertex = donor_geometry->vertex[markDonor][iVertex];
      const auto iPoint = donor_vertex->GetNode();

      if (!donor_geometry->node[iPoint]->GetDomain()) continue;

      const auto nDonor = donor_vertex->GetnDonorPoints();
      sendGlobalIndex[iVertexLocal] = donor_geometry->node[iPoint]->GetGlobalIndex();
      sendDonorRowPtr[iVertexLocal+1] = sendDonorRowPtr[iVertexLocal] + nDonor;

      for (auto i = 0u; i < nDonor; ++i) {
        sendDonorIndex[iDonor] = donor_vertex->GetInterpDonorPoint(i);
        sendDonorCoeff[iDonor] = donor_vertex->GetDonorCoeff(i);
        ++iDonor;
      }
      ++iVertexLocal;
    }

    /*--- Communicate donor interpolation matrix and info. We only gather the
     *    matrix in ranks that need it, i.e. have target vertices, to avoid
     *    potentially huge memory usage when running in parallel. This is an
     *    in-place Gatherv (done manually due to AD issues in the past). ---*/

    vector<int> iSendProcessor;
    for (int iProcessor = 0; iProcessor < nProcessor; ++iProcessor)
      if (allNumVertexDonor[iProcessor] != 0)
        iSendProcessor.push_back(iProcessor);

    const int nSend = iSendProcessor.size();

    vector<long*> GlobalIndex(nSend,nullptr), DonorRowPtr(nSend,nullptr), DonorIndex(nSend,nullptr);
    vector<su2double*> DonorCoeff(nSend,nullptr);

    /*--- For each "target processor" that needs the interpolation matrix. ---*/
    for (int iProcessor = 0; iProcessor < nProcessor; ++iProcessor) {
      if (allNumVertexTarget[iProcessor] == 0) continue;

      /*--- For each "donor processor that contain a part of the matrix. ---*/
      for (int iSend = 0; iSend < nSend; ++iSend) {
        const auto jProcessor = iSendProcessor[iSend];

        const auto numRows = allNumVertexDonor[jProcessor];
        const auto numCoeff = allNumNodeDonor[jProcessor];

        if ((rank == iProcessor) && (rank == jProcessor)) {
          /*--- "Self" communication. ---*/
          GlobalIndex[iSend] = sendGlobalIndex.data();
          DonorRowPtr[iSend] = sendDonorRowPtr.data();
          DonorIndex[iSend] = sendDonorIndex.data();
          DonorCoeff[iSend] = sendDonorCoeff.data();
        }
        else if (rank == iProcessor) {
          /*--- "I'm" the target, allocate and receive. ---*/
          GlobalIndex[iSend] = new long [numRows];
          DonorRowPtr[iSend] = new long [numRows+1];
          DonorIndex[iSend] = new long [numCoeff];
          DonorCoeff[iSend] = new su2double [numCoeff];
          SU2_MPI::Recv(GlobalIndex[iSend], numRows,   MPI_LONG, jProcessor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          SU2_MPI::Recv(DonorRowPtr[iSend], numRows+1, MPI_LONG, jProcessor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          SU2_MPI::Recv(DonorIndex[iSend],  numCoeff,  MPI_LONG, jProcessor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          SU2_MPI::Recv(DonorCoeff[iSend], numCoeff, MPI_DOUBLE, jProcessor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (rank == jProcessor) {
          /*--- "I'm" the donor, send. ---*/
          SU2_MPI::Send(sendGlobalIndex.data(), numRows,   MPI_LONG, iProcessor, 0, MPI_COMM_WORLD);
          SU2_MPI::Send(sendDonorRowPtr.data(), numRows+1, MPI_LONG, iProcessor, 0, MPI_COMM_WORLD);
          SU2_MPI::Send(sendDonorIndex.data(),  numCoeff,  MPI_LONG, iProcessor, 0, MPI_COMM_WORLD);
          SU2_MPI::Send(sendDonorCoeff.data(), numCoeff, MPI_DOUBLE, iProcessor, 0, MPI_COMM_WORLD);
        }
      }
    }

    /*--- Loop over the vertices on the target marker, define one row of the transpose matrix. ---*/

    SU2_OMP_PARALLEL_(for schedule(dynamic,roundUpDiv(nVertexTarget, 2*omp_get_max_threads())))
    for (auto iVertex = 0ul; iVertex < nVertexTarget; ++iVertex) {

      auto target_vertex = target_geometry->vertex[markTarget][iVertex];
      const auto iPoint = target_vertex->GetNode();

      if (!target_geometry->node[iPoint]->GetDomain()) continue;

      /*--- Any point of the donor geometry, that has this target point as a donor, becomes a donor. ---*/
      const long targetGlobalIndex = target_geometry->node[iPoint]->GetGlobalIndex();

      /*--- Traverse entire matrix to count donors. ---*/
      auto nDonor = 0ul;
      for (int iSend = 0; iSend < nSend; ++iSend) {
        const auto iProcessor = iSendProcessor[iSend];
        for (auto iCoeff = 0ul; iCoeff < allNumNodeDonor[iProcessor]; ++iCoeff)
          nDonor += (targetGlobalIndex == DonorIndex[iSend][iCoeff]);
      }

      target_vertex->Allocate_DonorInfo(nDonor);

      /*--- Traverse matrix again to set interpolation coefficients. ---*/
      auto iDonor = 0ul;
      for (int iSend = 0; iSend < nSend; ++iSend) {
        const auto iProcessor = iSendProcessor[iSend];

        for (auto iVertexDonor = 0ul; iVertexDonor < allNumVertexDonor[iProcessor]; ++iVertexDonor) {
          for (auto iCoeff = DonorRowPtr[iSend][iVertexDonor];
               iCoeff < DonorRowPtr[iSend][iVertexDonor+1]; ++iCoeff)
          {
            if (targetGlobalIndex != DonorIndex[iSend][iCoeff]) continue;

            target_vertex->SetInterpDonorProcessor(iDonor, iProcessor);
            target_vertex->SetDonorCoeff(iDonor, DonorCoeff[iSend][iCoeff]);
            target_vertex->SetInterpDonorPoint(iDonor, GlobalIndex[iSend][iVertexDonor]);
            ++iDonor;
            /*--- we could break here on the assumption that donors do not repeat, but... ---*/
          }
          if (iDonor == nDonor) break;
        }
      }

    } // end target loop

    /*--- Free the heap allocations. ---*/
    for (auto ptr : GlobalIndex) if (ptr != sendGlobalIndex.data()) delete [] ptr;
    for (auto ptr : DonorRowPtr) if (ptr != sendDonorRowPtr.data()) delete [] ptr;
    for (auto ptr : DonorIndex) if (ptr != sendDonorIndex.data()) delete [] ptr;
    for (auto ptr : DonorCoeff) if (ptr != sendDonorCoeff.data()) delete [] ptr;

  } // end marker loop

}
