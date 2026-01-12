/*!
 * \file CInterface.cpp
 * \brief Main subroutines for MPI transfer of information between zones
 * \author R. Sanchez
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/interfaces/CInterface.hpp"
#include "../../../Common/include/interface_interpolation/CInterpolator.hpp"
#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CInterface::CInterface() :
  rank(SU2_MPI::GetRank()),
  size(SU2_MPI::GetSize()) {
}

CInterface::CInterface(unsigned short val_nVar, unsigned short val_nConst) :
  rank(SU2_MPI::GetRank()),
  size(SU2_MPI::GetSize()),
  nVar(val_nVar) {

  Physical_Constants = new su2double[val_nConst] ();
  Donor_Variable     = new su2double[val_nVar] ();
  Target_Variable    = new su2double[val_nVar] ();

  /*--- By default, the value is aggregated in the transfer routine ---*/
  valAggregated      = true;
}

CInterface::~CInterface() {

  delete [] Physical_Constants;
  delete [] Donor_Variable;
  delete [] Target_Variable;

  delete[] SpanValueCoeffTarget;
  delete[] SpanLevelDonor;
}

void CInterface::BroadcastData(const CInterpolator& interpolator,
                               CSolver *donor_solution, CSolver *target_solution,
                               CGeometry *donor_geometry, CGeometry *target_geometry,
                               const CConfig *donor_config, const CConfig *target_config) {
  static_assert(su2activematrix::Storage == StorageType::RowMajor,"");

  GetPhysical_Constants(donor_solution, target_solution, donor_geometry, target_geometry,
                        donor_config, target_config);

  /*--- Loop over interface markers. ---*/

  for (auto iMarkerInt = 0u; iMarkerInt < donor_config->GetMarker_n_ZoneInterface()/2; iMarkerInt++) {

    /*--- Check if this interface connects the two zones, if not continue. ---*/

    const auto markDonor = donor_config->FindInterfaceMarker(iMarkerInt);
    const auto markTarget = target_config->FindInterfaceMarker(iMarkerInt);

    if(!CInterpolator::CheckInterfaceBoundary(markDonor, markTarget)) continue;

    /*--- Count donor vertices on this rank. ---*/

    int nLocalVertexDonor = 0;
    if (markDonor >= 0) {
      for (auto iVertex = 0ul; iVertex < donor_geometry->GetnVertex(markDonor); iVertex++) {
        auto Point_Donor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
        /*--- Only domain points are donors. ---*/
        nLocalVertexDonor += donor_geometry->nodes->GetDomain(Point_Donor);
      }
    }

    /*--- Gather donor counts and compute total sizes, and displacements (cumulative
     * sums) to perform an Allgatherv of donor indices and variables. ---*/

    vector<int> nAllVertexDonor(size), nAllVarCounts(size), displIdx(size,0), displVar(size);
    SU2_MPI::Allgather(&nLocalVertexDonor, 1, MPI_INT, nAllVertexDonor.data(), 1, MPI_INT, SU2_MPI::GetComm());

    for (int i = 0; i < size; ++i) {
      nAllVarCounts[i] = nAllVertexDonor[i] * nVar;
      if(i) displIdx[i] = displIdx[i-1] + nAllVertexDonor[i-1];
      displVar[i] = displIdx[i] * nVar;
    }

    /*--- Fill send buffers. ---*/

    vector<unsigned long> sendDonorIdx(nLocalVertexDonor);
    su2activematrix sendDonorVar(nLocalVertexDonor, nVar);

    if (markDonor >= 0) {

      /*--- Apply contact resistance if specified. ---*/
      
      SetContactResistance(donor_config->GetContactResistance(iMarkerInt));

      for (auto iVertex = 0ul, iSend = 0ul; iVertex < donor_geometry->GetnVertex(markDonor); iVertex++) {
        const auto iPoint = donor_geometry->vertex[markDonor][iVertex]->GetNode();

        /*--- If this processor owns the node. ---*/
        if (donor_geometry->nodes->GetDomain(iPoint)) {

          GetDonor_Variable(donor_solution, donor_geometry, donor_config, markDonor, iVertex, iPoint);
          for (auto iVar = 0u; iVar < nVar; iVar++) sendDonorVar(iSend, iVar) = Donor_Variable[iVar];

          sendDonorIdx[iSend] = donor_geometry->nodes->GetGlobalIndex(iPoint);
          ++iSend;
        }
      }
    }

    /*--- Gather data. ---*/

    const auto nGlobalVertexDonor = displIdx.back() + nAllVertexDonor.back();

    vector<unsigned long> donorIdx(nGlobalVertexDonor);
    su2activematrix donorVar(nGlobalVertexDonor, nVar);

    SU2_MPI::Allgatherv(sendDonorIdx.data(), sendDonorIdx.size(), MPI_UNSIGNED_LONG, donorIdx.data(),
                        nAllVertexDonor.data(), displIdx.data(), MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

    SU2_MPI::Allgatherv(sendDonorVar.data(), sendDonorVar.size(), MPI_DOUBLE, donorVar.data(),
                        nAllVarCounts.data(), displVar.data(), MPI_DOUBLE, SU2_MPI::GetComm());

    /*--- This rank does not need to do more work. ---*/
    if (markTarget < 0) continue;

    /*--- Sort the donor information by index to then use binary searches. ---*/

    vector<size_t> order(donorIdx.size());
    iota(order.begin(), order.end(), 0ul);
    sort(order.begin(), order.end(), [&donorIdx](size_t i, size_t j) {return donorIdx[i] < donorIdx[j];} );

    /*--- inplace permutation. ---*/
    for (size_t i = 0; i < order.size(); ++i) {
      auto j = order[i];
      while (j < i) j = order[j];
      if (i == j) continue;
      swap(donorIdx[i], donorIdx[j]);
      for (auto iVar = 0u; iVar < nVar; ++iVar)
        swap(donorVar(i,iVar), donorVar(j,iVar));
    }

    /*--- Loop over target vertices. ---*/

    for (auto iVertex = 0ul; iVertex < target_geometry->GetnVertex(markTarget); iVertex++) {
      const auto iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();

      if (!target_geometry->nodes->GetDomain(iPoint)) continue;

      auto& targetVertex = interpolator.targetVertices[markTarget][iVertex];
      const auto nDonorPoints = targetVertex.nDonor();

      InitializeTarget_Variable(target_solution, markTarget, iVertex, nDonorPoints);

      /*--- For the number of donor points. ---*/
      for (auto iDonorPoint = 0ul; iDonorPoint < nDonorPoints; iDonorPoint++) {

        /*--- Get the global index of the donor and the interpolation coefficient. ---*/

        const auto donorGlobalIndex = targetVertex.globalPoint[iDonorPoint];
        const auto donorCoeff = targetVertex.coefficient[iDonorPoint];

        /*--- Find the index of the global donor point in the donor data. ---*/

        const auto idx = lower_bound(donorIdx.begin(), donorIdx.end(), donorGlobalIndex) - donorIdx.begin();
        assert(idx < static_cast<long>(donorIdx.size()));

        /*--- Recover the Target_Variable from the buffer of variables. ---*/
        RecoverTarget_Variable(donorVar[idx], donorCoeff);

        /*--- If the value is not directly aggregated in the previous function. ---*/
        if (!valAggregated)
          SetTarget_Variable(target_solution, target_geometry, target_config, markTarget, iVertex, iPoint);
      }

      /*--- If we have aggregated the values in the function RecoverTarget_Variable, the set is outside the loop. ---*/
      if (valAggregated)
        SetTarget_Variable(target_solution, target_geometry, target_config, markTarget, iVertex, iPoint);
    }
  }
}
