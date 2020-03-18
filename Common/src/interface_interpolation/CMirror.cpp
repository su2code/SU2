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


CMirror::CMirror(CGeometry ****geometry_container, CConfig **config,  unsigned int iZone,
                 unsigned int jZone) : CInterpolator(geometry_container, config, iZone, jZone) {
  Set_TransferCoeff(config);
}

void CMirror::Set_TransferCoeff(CConfig **config) {

  unsigned long iVertex, jVertex;
  unsigned long iPoint;
  unsigned short iDonor=0, iFace=0, iTarget=0;

  unsigned short nMarkerInt;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;

  unsigned int nNodes=0, iNodes=0;
  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Donor = 0;
  unsigned long pGlobalPoint = 0;
  int iProcessor;

  unsigned long nLocalFace_Donor = 0, nLocalFaceNodes_Donor=0;

  unsigned long faceindex;

  int nProcessor = size;

  su2double *Buffer_Send_Coeff, *Buffer_Receive_Coeff;
  su2double coeff;

  Buffer_Send_nFace_Donor= new unsigned long [1];
  Buffer_Send_nFaceNodes_Donor= new unsigned long [1];

  Buffer_Receive_nFace_Donor = new unsigned long [nProcessor];
  Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

  /*--- Number of markers on the interface ---*/
  nMarkerInt = (config[targetZone]->GetMarker_n_ZoneInterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
   /*--- Procedure:
    * - Loop through vertices of the aero grid
    * - Find nearest element and allocate enough space in the aero grid donor point info
    * - Set the transfer coefficient values
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor = Find_InterfaceMarker(config[donorZone], iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(!CheckInterfaceBoundary(markDonor, markTarget)) continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;

    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;

    /*-- Collect the number of donor nodes: re-use 'Face' containers --*/
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;
    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {
      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex

      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        nLocalFaceNodes_Donor+=nNodes;
        nLocalFace_Donor++;
      }
    }

    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;

    /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG,
                       Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG,
                       Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
#else
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];
#endif

    /*-- Send donor info --*/
    Buffer_Send_FaceIndex = new unsigned long[MaxFace_Donor];
    Buffer_Send_FaceNodes = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_GlobalPoint = new long[MaxFaceNodes_Donor];
    Buffer_Send_Coeff = new su2double[MaxFaceNodes_Donor];

    Buffer_Receive_FaceIndex = new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_GlobalPoint = new long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_Coeff = new su2double[MaxFaceNodes_Donor*nProcessor];

    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex]=0;
    }
    for (iVertex=0; iVertex<MaxFaceNodes_Donor; iVertex++) {
      Buffer_Send_FaceNodes[iVertex]=0;
      Buffer_Send_GlobalPoint[iVertex]=0;
      Buffer_Send_Coeff[iVertex]=0.0;
    }
    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex]=0;
    }

    Buffer_Send_FaceIndex[0]=rank*MaxFaceNodes_Donor;
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {

      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex
      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        for (iDonor=0; iDonor<nNodes; iDonor++) {
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
    for (iVertex = 0; iVertex<nVertexTarget; iVertex++) {

      iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();
      if (!target_geometry->node[iPoint]->GetDomain()) continue;

      long Global_Point = target_geometry->node[iPoint]->GetGlobalIndex();
      nNodes = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {
          faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
          iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
          for (iTarget=0; iTarget<iNodes; iTarget++) {
            if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget])
              nNodes++;
          }
        }
      }

      target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
      target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

      iDonor = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {

          faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
          iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
          for (iTarget=0; iTarget<iNodes; iTarget++) {
            if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget]) {
              coeff = Buffer_Receive_Coeff[faceindex+iTarget];
              pGlobalPoint = Buffer_Receive_FaceNodes[faceindex+iTarget];
              target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(iDonor,pGlobalPoint);
              target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor,coeff);
              target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, iProcessor);
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

  delete[] Buffer_Send_nFace_Donor;
  delete[] Buffer_Send_nFaceNodes_Donor;

  delete[] Buffer_Receive_nFace_Donor;
  delete[] Buffer_Receive_nFaceNodes_Donor;

}
