/*!
 * \file CMirror.cpp
 * \brief Implementation of mirror interpolation (conservative approach in FSI problems).
 * \author H. Kline
 * \version 7.0.2 "Blackbird"
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


CMirror::CMirror(CGeometry ****geometry_container, const CConfig* const* config,  unsigned int iZone,
                 unsigned int jZone) : CInterpolator(geometry_container, config, iZone, jZone) {
  Set_TransferCoeff(config);
}

void CMirror::Set_TransferCoeff(const CConfig* const* config) {

  const int nProcessor = size;

  Buffer_Receive_nFace_Donor = new unsigned long [nProcessor];
  Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

  /*--- Number of markers on the interface ---*/
  const auto nMarkerInt = (config[targetZone]->GetMarker_n_ZoneInterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (unsigned short iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
   /*--- Procedure:
    * - Loop through vertices of the aero grid
    * - Find nearest element and allocate enough space in the aero grid donor point info
    * - Set the transfer coefficient values
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    const auto markDonor = Find_InterfaceMarker(config[donorZone], iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    const auto markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    unsigned long nVertexDonor = 0, nVertexTarget = 0;
    if (markDonor != -1) nVertexDonor = donor_geometry->GetnVertex( markDonor );
    if (markTarget != -1) nVertexTarget = target_geometry->GetnVertex( markTarget );

    /*-- Collect the number of donor nodes: re-use 'Face' containers --*/
    auto nLocalFace_Donor = 0ul;
    auto nLocalFaceNodes_Donor = 0ul;
    for (auto jVertex = 0ul; jVertex<nVertexDonor; jVertex++) {
      auto Point_Donor = donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex

      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        auto nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        nLocalFaceNodes_Donor += nNodes;
        nLocalFace_Donor++;
      }
    }

    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;

    /*--- Send Interface vertex information --*/
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG,
                       Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG,
                       Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;

    /*-- Send donor info and init to 0 --*/
    Buffer_Send_FaceIndex = new unsigned long[MaxFace_Donor] ();
    Buffer_Send_FaceNodes = new unsigned long[MaxFaceNodes_Donor] ();
    Buffer_Send_GlobalPoint = new long[MaxFaceNodes_Donor] ();
    auto Buffer_Send_Coeff = new su2double[MaxFaceNodes_Donor] ();

    Buffer_Receive_FaceIndex = new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_GlobalPoint = new long[MaxFaceNodes_Donor*nProcessor];
    auto Buffer_Receive_Coeff = new su2double[MaxFaceNodes_Donor*nProcessor];

    Buffer_Send_FaceIndex[0]=rank*MaxFaceNodes_Donor;
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    for (auto jVertex = 0ul; jVertex < nVertexDonor; jVertex++) {

      auto Point_Donor = donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex

      if (!donor_geometry->node[Point_Donor]->GetDomain()) continue;

      auto nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
      for (auto iDonor = 0; iDonor < nNodes; iDonor++) {
        Buffer_Send_FaceNodes[nLocalFaceNodes_Donor] = donor_geometry->node[Point_Donor]->GetGlobalIndex();
        Buffer_Send_GlobalPoint[nLocalFaceNodes_Donor] =
            donor_geometry->vertex[markDonor][jVertex]->GetInterpDonorPoint(iDonor);
        Buffer_Send_Coeff[nLocalFaceNodes_Donor] =
            donor_geometry->vertex[markDonor][jVertex]->GetDonorCoeff(iDonor);
        nLocalFaceNodes_Donor++;
      }
      Buffer_Send_FaceIndex[nLocalFace_Donor+1] =Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
      nLocalFace_Donor++;
    }

    SU2_MPI::Allgather(Buffer_Send_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG,
                       Buffer_Receive_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_GlobalPoint, MaxFaceNodes_Donor, MPI_LONG,
                       Buffer_Receive_GlobalPoint, MaxFaceNodes_Donor, MPI_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Coeff, MaxFaceNodes_Donor, MPI_DOUBLE,
                       Buffer_Receive_Coeff, MaxFaceNodes_Donor, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG,
                       Buffer_Receive_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    /*--- Loop over the vertices on the target Marker ---*/
    SU2_OMP_PARALLEL_(for schedule(dynamic,roundUpDiv(nVertexTarget,2*omp_get_max_threads())))
    for (auto iVertex = 0ul; iVertex < nVertexTarget; iVertex++) {

      auto target_vertex = target_geometry->vertex[markTarget][iVertex];
      auto iPoint = target_vertex->GetNode();

      if (!target_geometry->node[iPoint]->GetDomain()) continue;

      const long Global_Point = target_geometry->node[iPoint]->GetGlobalIndex();

      auto nNodes = 0;
      for (int iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (auto iFace = 0ul; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {
          auto faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
          auto iNodes = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1] - faceindex;
          for (auto iTarget = 0ul; iTarget<iNodes; iTarget++)
            nNodes += (Global_Point == long(Buffer_Receive_GlobalPoint[faceindex+iTarget]));
        }
      }

      target_vertex->Allocate_DonorInfo(nNodes);

      for (int iProcessor = 0, iDonor = 0; iProcessor < nProcessor; iProcessor++) {
        for (auto iFace = 0ul; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {

          auto faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
          auto iNodes = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1] - faceindex;

          for (auto iTarget = 0ul; iTarget < iNodes; iTarget++) {
            if (Global_Point == long(Buffer_Receive_GlobalPoint[faceindex+iTarget])) {
              auto coeff = Buffer_Receive_Coeff[faceindex+iTarget];
              auto pGlobalPoint = Buffer_Receive_FaceNodes[faceindex+iTarget];
              target_vertex->SetInterpDonorPoint(iDonor,pGlobalPoint);
              target_vertex->SetDonorCoeff(iDonor,coeff);
              target_vertex->SetInterpDonorProcessor(iDonor, iProcessor);
              iDonor++;
            }
          }
        }
      }
    }

    delete[] Buffer_Send_FaceIndex;
    delete[] Buffer_Send_FaceNodes;
    delete[] Buffer_Send_GlobalPoint;
    delete[] Buffer_Send_Coeff;

    delete[] Buffer_Receive_FaceIndex;
    delete[] Buffer_Receive_FaceNodes;
    delete[] Buffer_Receive_GlobalPoint;
    delete[] Buffer_Receive_Coeff;
  }

  delete[] Buffer_Receive_nFace_Donor;
  delete[] Buffer_Receive_nFaceNodes_Donor;

}
