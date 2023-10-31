/*!
 * \file CInterface.cpp
 * \brief Main subroutines for MPI transfer of information between zones
 * \author R. Sanchez
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

void CInterface::PreprocessAverage(CGeometry *donor_geometry, CGeometry *target_geometry,
                                   const CConfig *donor_config, const CConfig *target_config,
                                   unsigned short iMarkerInt){

  unsigned short  nMarkerDonor, nMarkerTarget;    // Number of markers on the interface, donor and target side
  unsigned short  iMarkerDonor, iMarkerTarget;    // Variables for iteration over markers
  unsigned short iSpan,jSpan, tSpan = 0, kSpan = 0, nSpanDonor, nSpanTarget, Donor_Flag = 0, Target_Flag = 0;
  int Marker_Donor = -1, Marker_Target = -1;

  const su2double *SpanValuesDonor, *SpanValuesTarget;
  su2double dist, test, dist2, test2;

  nMarkerDonor   = donor_geometry->GetnMarker();
  nMarkerTarget  = target_geometry->GetnMarker();
  //TODO turbo this approach only works if all the turboamchinery marker
  //    of all zones have the same amount of span wise sections.
  //TODO turbo initialization needed for the MPI routine should be place somewhere else.
  nSpanDonor     = donor_config->GetnSpanWiseSections();
  nSpanTarget    = target_config->GetnSpanWiseSections();

  /*--- On the donor side ---*/
  for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
    /*--- If the tag GetMarker_All_MixingPlaneInterface equals the index we are looping at ---*/
    if ( donor_config->GetMarker_All_MixingPlaneInterface(iMarkerDonor) == iMarkerInt ){
      /*--- We have identified the local index of the Donor marker ---*/
      /*--- Now we are going to store the average values that belong to Marker_Donor on each processor ---*/
      /*--- Store the identifier for the structural marker ---*/
      Marker_Donor = iMarkerDonor;
      Donor_Flag = donor_config->GetMarker_All_TurbomachineryFlag(iMarkerDonor);
      /*--- Exit the for loop: we have found the local index for Mixing-Plane interface ---*/
      break;
    }
          /*--- If the tag hasn't matched any tag within the donor markers ---*/
      Marker_Donor = -1;
      Donor_Flag   = -1;
   
  }

#ifdef HAVE_MPI
  auto BuffMarkerDonor = new int[size];
  auto BuffDonorFlag = new int[size];
  for (int iSize=0; iSize<size; iSize++){
    BuffMarkerDonor[iSize] = -1;
    BuffDonorFlag[iSize] = -1;
  }

  SU2_MPI::Allgather(&Marker_Donor, 1 , MPI_INT, BuffMarkerDonor, 1, MPI_INT, SU2_MPI::GetComm());
  SU2_MPI::Allgather(&Donor_Flag, 1 , MPI_INT, BuffDonorFlag, 1, MPI_INT, SU2_MPI::GetComm());

  Marker_Donor= -1;
  Donor_Flag= -1;

  for (int iSize=0; iSize<size; iSize++){
    if(BuffMarkerDonor[iSize] > 0.0){
      Marker_Donor = BuffMarkerDonor[iSize];
      Donor_Flag = BuffDonorFlag[iSize];
      break;
    }
  }
  delete [] BuffMarkerDonor;
  delete [] BuffDonorFlag;
#endif

  /*--- On the target side we have to identify the marker as well ---*/

  for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
    /*--- If the tag GetMarker_All_MixingPlaneInterface(iMarkerTarget) equals the index we are looping at ---*/
    if ( target_config->GetMarker_All_MixingPlaneInterface(iMarkerTarget) == iMarkerInt ){
      /*--- Store the identifier for the fluid marker ---*/

      // here i should then store it in the target zone

      Marker_Target = iMarkerTarget;
      Target_Flag = target_config->GetMarker_All_TurbomachineryFlag(iMarkerTarget);
      /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
      break;
    }
          /*--- If the tag hasn't matched any tag within the Flow markers ---*/
      Marker_Target = -1;
   
  }

  if (Marker_Target != -1 && Marker_Donor != -1){

    SpanValuesDonor  = donor_geometry->GetSpanWiseValue(Donor_Flag);
    SpanValuesTarget = target_geometry->GetSpanWiseValue(Target_Flag);


    for(iSpan = 1; iSpan <nSpanTarget-1; iSpan++){
      dist  = 10E+06;
      dist2 = 10E+06;
      for(jSpan = 0; jSpan < nSpanDonor;jSpan++){
        test = abs(SpanValuesTarget[iSpan] - SpanValuesDonor[jSpan]);
        test2 = abs(SpanValuesTarget[iSpan] - SpanValuesDonor[jSpan]);
        if(test < dist && SpanValuesTarget[iSpan] > SpanValuesDonor[jSpan]){
          dist = test;
          kSpan = jSpan;
        }
        if(test2 < dist2){
          dist2 = test2;
          tSpan =jSpan;
        }

      }
      switch(donor_config->GetKind_MixingPlaneInterface()){
        case MATCHING:
          SpanLevelDonor[iSpan]        = iSpan;
          SpanValueCoeffTarget[iSpan]  = 0.0;
          break;
        case NEAREST_SPAN:
          SpanLevelDonor[iSpan]        = tSpan;
          SpanValueCoeffTarget[iSpan]  = 0.0;
          break;
        case LINEAR_INTERPOLATION:
          SpanLevelDonor[iSpan]        = kSpan;
          SpanValueCoeffTarget[iSpan]  = (SpanValuesTarget[iSpan] - SpanValuesDonor[kSpan])
                                         /(SpanValuesDonor[kSpan + 1] - SpanValuesDonor[kSpan]);
          break;
        default:
          SU2_MPI::Error("MixingPlane interface option not implemented yet", CURRENT_FUNCTION);
          break;

      }
    }
  }

}


void CInterface::AllgatherAverage(CSolver *donor_solution, CSolver *target_solution,
                                  CGeometry *donor_geometry, CGeometry *target_geometry,
                                  const CConfig *donor_config, const CConfig *target_config, unsigned short iMarkerInt){

  unsigned short  nMarkerDonor, nMarkerTarget;    // Number of markers on the interface, donor and target side
  unsigned short  iMarkerDonor, iMarkerTarget;    // Variables for iteration over markers
  unsigned short iSpan, nSpanDonor, nSpanTarget;
  int Marker_Donor = -1, Marker_Target = -1;
  su2double *avgPressureDonor = nullptr, *avgDensityDonor = nullptr, *avgNormalVelDonor = nullptr,
      *avgTangVelDonor = nullptr, *avg3DVelDonor = nullptr, *avgNuDonor = nullptr,
      *avgOmegaDonor = nullptr, *avgKineDonor = nullptr;
  su2double *avgPressureTarget = nullptr, *avgDensityTarget = nullptr, *avgNormalVelTarget = nullptr,
      *avg3DVelTarget = nullptr, *avgTangVelTarget = nullptr, *avgNuTarget = nullptr,
      *avgOmegaTarget = nullptr, *avgKineTarget = nullptr;

#ifdef HAVE_MPI
  int iSize;
  su2double *BuffAvgPressureDonor = nullptr, *BuffAvgDensityDonor = nullptr, *BuffAvgNormalVelDonor = nullptr,
      *BuffAvg3DVelDonor = nullptr, *BuffAvgTangVelDonor = nullptr, *BuffAvgNuDonor = nullptr,
      *BuffAvgKineDonor = nullptr, *BuffAvgOmegaDonor = nullptr;
  int nSpanSize, *BuffMarkerDonor;
#endif


  nMarkerTarget  = target_geometry->GetnMarker();
  nMarkerDonor   = donor_geometry->GetnMarker();
  nSpanDonor     = donor_config->GetnSpanWiseSections() +1;
  nSpanTarget    = target_config->GetnSpanWiseSections() +1;


  avgDensityDonor                  = new su2double[nSpanDonor];
  avgPressureDonor                 = new su2double[nSpanDonor];
  avgNormalVelDonor                = new su2double[nSpanDonor];
  avgTangVelDonor                  = new su2double[nSpanDonor];
  avg3DVelDonor                    = new su2double[nSpanDonor];
  avgNuDonor                       = new su2double[nSpanDonor];
  avgKineDonor                     = new su2double[nSpanDonor];
  avgOmegaDonor                    = new su2double[nSpanDonor];

  for (iSpan = 0; iSpan < nSpanDonor; iSpan++){
    avgDensityDonor[iSpan]         = -1.0;
    avgPressureDonor[iSpan]        = -1.0;
    avgNormalVelDonor[iSpan]       = -1.0;
    avgTangVelDonor[iSpan]         = -1.0;
    avg3DVelDonor[iSpan]           = -1.0;
    avgNuDonor[iSpan]              = -1.0;
    avgKineDonor[iSpan]            = -1.0;
    avgOmegaDonor[iSpan]           = -1.0;
  }

  avgDensityTarget                 = new su2double[nSpanTarget];
  avgPressureTarget                = new su2double[nSpanTarget];
  avgNormalVelTarget               = new su2double[nSpanTarget];
  avgTangVelTarget                 = new su2double[nSpanTarget];
  avg3DVelTarget                   = new su2double[nSpanTarget];
  avgNuTarget                      = new su2double[nSpanTarget];
  avgKineTarget                    = new su2double[nSpanTarget];
  avgOmegaTarget                   = new su2double[nSpanTarget];


  for (iSpan = 0; iSpan < nSpanTarget; iSpan++){
    avgDensityTarget[iSpan]        = -1.0;
    avgPressureTarget[iSpan]       = -1.0;
    avgNormalVelTarget[iSpan]      = -1.0;
    avgTangVelTarget[iSpan]        = -1.0;
    avg3DVelTarget[iSpan]          = -1.0;
    avgNuTarget[iSpan]             = -1.0;
    avgKineTarget[iSpan]           = -1.0;
    avgOmegaTarget[iSpan]          = -1.0;
  }

  /*--- Outer loop over the markers on the Mixing-Plane interface: compute one by one ---*/
  /*--- The tags are always an integer greater than 1: loop from 1 to nMarkerMixingPlane ---*/
  Marker_Donor = -1;
  Marker_Target = -1;

  /*--- The donor and target markers are tagged with the same index.
   *--- This is independent of the MPI domain decomposition.
   *--- We need to loop over all markers on both sides  ---*/

  /*--- On the donor side ---*/

  for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++){
    /*--- If the tag GetMarker_All_MixingPlaneInterface equals the index we are looping at ---*/
    if ( donor_config->GetMarker_All_MixingPlaneInterface(iMarkerDonor) == iMarkerInt ){
      /*--- We have identified the local index of the Donor marker ---*/
      /*--- Now we are going to store the average values that belong to Marker_Donor on each processor ---*/
      /*--- Store the identifier for the structural marker ---*/
      Marker_Donor = iMarkerDonor;
      /*--- Exit the for loop: we have found the local index for Mixing-Plane interface ---*/
      break;
    }
          /*--- If the tag hasn't matched any tag within the donor markers ---*/
      Marker_Donor = -1;
   
  }
  /*--- Here we want to make available the quantities for all the processors and collect them in a buffer
   * for each span of the donor the span-wise height vector also so
   * that then we can interpolate on the target side  ---*/
  if (Marker_Donor != -1){
    for(iSpan = 0; iSpan < nSpanDonor; iSpan++){
      GetDonor_Variable(donor_solution, donor_geometry, donor_config, Marker_Donor, iSpan, rank);
      avgDensityDonor[iSpan]          = Donor_Variable[0];
      avgPressureDonor[iSpan]         = Donor_Variable[1];
      avgNormalVelDonor[iSpan]        = Donor_Variable[2];
      avgTangVelDonor[iSpan]          = Donor_Variable[3];
      avg3DVelDonor[iSpan]            = Donor_Variable[4];
      avgNuDonor[iSpan]               = Donor_Variable[5];
      avgKineDonor[iSpan]             = Donor_Variable[6];
      avgOmegaDonor[iSpan]            = Donor_Variable[7];
    }
  }

#ifdef HAVE_MPI
  nSpanSize = size*nSpanDonor;
  BuffAvgDensityDonor                 = new su2double[nSpanSize];
  BuffAvgPressureDonor                = new su2double[nSpanSize];
  BuffAvgNormalVelDonor               = new su2double[nSpanSize];
  BuffAvgTangVelDonor                 = new su2double[nSpanSize];
  BuffAvg3DVelDonor                   = new su2double[nSpanSize];
  BuffAvgNuDonor                      = new su2double[nSpanSize];
  BuffAvgKineDonor                    = new su2double[nSpanSize];
  BuffAvgOmegaDonor                   = new su2double[nSpanSize];
  BuffMarkerDonor                     = new int[size];

  for (iSpan=0;iSpan<nSpanSize;iSpan++){
    BuffAvgDensityDonor[iSpan]        = -1.0;
    BuffAvgPressureDonor[iSpan]       = -1.0;
    BuffAvgNormalVelDonor[iSpan]      = -1.0;
    BuffAvgTangVelDonor[iSpan]        = -1.0;
    BuffAvg3DVelDonor[iSpan]          = -1.0;
    BuffAvgNuDonor[iSpan]             = -1.0;
    BuffAvgKineDonor[iSpan]           = -1.0;
    BuffAvgOmegaDonor[iSpan]          = -1.0;
  }

  for (iSize=0; iSize<size;iSize++){
    BuffMarkerDonor[iSize]            = -1;
  }

  SU2_MPI::Allgather(avgDensityDonor, nSpanDonor , MPI_DOUBLE, BuffAvgDensityDonor,
                     nSpanDonor, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(avgPressureDonor, nSpanDonor , MPI_DOUBLE, BuffAvgPressureDonor,
                     nSpanDonor, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(avgNormalVelDonor, nSpanDonor , MPI_DOUBLE, BuffAvgNormalVelDonor,
                     nSpanDonor, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(avgTangVelDonor, nSpanDonor , MPI_DOUBLE, BuffAvgTangVelDonor,
                     nSpanDonor, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(avg3DVelDonor, nSpanDonor , MPI_DOUBLE, BuffAvg3DVelDonor,
                     nSpanDonor, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(avgNuDonor, nSpanDonor , MPI_DOUBLE, BuffAvgNuDonor,
                     nSpanDonor, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(avgKineDonor, nSpanDonor , MPI_DOUBLE, BuffAvgKineDonor,
                     nSpanDonor, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(avgOmegaDonor, nSpanDonor , MPI_DOUBLE, BuffAvgOmegaDonor,
                     nSpanDonor, MPI_DOUBLE, SU2_MPI::GetComm());
  SU2_MPI::Allgather(&Marker_Donor, 1 , MPI_INT, BuffMarkerDonor, 1, MPI_INT, SU2_MPI::GetComm());

  for (iSpan = 0; iSpan < nSpanDonor; iSpan++){
    avgDensityDonor[iSpan]            = -1.0;
    avgPressureDonor[iSpan]           = -1.0;
    avgNormalVelDonor[iSpan]          = -1.0;
    avgTangVelDonor[iSpan]            = -1.0;
    avg3DVelDonor[iSpan]              = -1.0;
    avgNuDonor[iSpan]                 = -1.0;
    avgKineDonor[iSpan]               = -1.0;
    avgOmegaDonor[iSpan]              = -1.0;
  }

  Marker_Donor= -1;

  for (iSize=0; iSize<size;iSize++){
    if(BuffAvgDensityDonor[nSpanDonor*iSize] > 0.0){
      for (iSpan = 0; iSpan < nSpanDonor; iSpan++){
        avgDensityDonor[iSpan]        = BuffAvgDensityDonor[nSpanDonor*iSize + iSpan];
        avgPressureDonor[iSpan]       = BuffAvgPressureDonor[nSpanDonor*iSize + iSpan];
        avgNormalVelDonor[iSpan]      = BuffAvgNormalVelDonor[nSpanDonor*iSize + iSpan];
        avgTangVelDonor[iSpan]        = BuffAvgTangVelDonor[nSpanDonor*iSize + iSpan];
        avg3DVelDonor[iSpan]          = BuffAvg3DVelDonor[nSpanDonor*iSize + iSpan];
        avgNuDonor[iSpan]             = BuffAvgNuDonor[nSpanDonor*iSize + iSpan];
        avgKineDonor[iSpan]           = BuffAvgKineDonor[nSpanDonor*iSize + iSpan];
        avgOmegaDonor[iSpan]          = BuffAvgOmegaDonor[nSpanDonor*iSize + iSpan];
      }
      Marker_Donor                    = BuffMarkerDonor[iSize];
      break;
    }
  }
  delete [] BuffAvgDensityDonor;
  delete [] BuffAvgPressureDonor;
  delete [] BuffAvgNormalVelDonor;
  delete [] BuffAvgTangVelDonor;
  delete [] BuffAvg3DVelDonor;
  delete [] BuffAvgNuDonor;
  delete [] BuffAvgKineDonor;
  delete [] BuffAvgOmegaDonor;
  delete [] BuffMarkerDonor;
#endif

  /*--- On the target side we have to identify the marker as well ---*/
  for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++){
    /*--- If the tag GetMarker_All_MixingPlaneInterface(iMarkerTarget) equals the index we are looping at ---*/
    if ( target_config->GetMarker_All_MixingPlaneInterface(iMarkerTarget) == iMarkerInt ){
      /*--- Store the identifier for the fluid marker ---*/
      Marker_Target = iMarkerTarget;
      /*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
      break;
    }
          /*--- If the tag hasn't matched any tag within the Flow markers ---*/
      Marker_Target = -1;
   
  }


  if (Marker_Target != -1 && Marker_Donor != -1){

    /*--- linear interpolation of the average value of for the internal span-wise levels ---*/
    for(iSpan = 1; iSpan < nSpanTarget -2 ; iSpan++){
      avgDensityTarget[iSpan]             = SpanValueCoeffTarget[iSpan]*(avgDensityDonor[SpanLevelDonor[iSpan] + 1]
                                            - avgDensityDonor[SpanLevelDonor[iSpan]]);
      avgDensityTarget[iSpan]            += avgDensityDonor[SpanLevelDonor[iSpan]];
      avgPressureTarget[iSpan]            = SpanValueCoeffTarget[iSpan]*(avgPressureDonor[SpanLevelDonor[iSpan] + 1]
                                            - avgPressureDonor[SpanLevelDonor[iSpan]]);
      avgPressureTarget[iSpan]           += avgPressureDonor[SpanLevelDonor[iSpan]];
      avgNormalVelTarget[iSpan]           = SpanValueCoeffTarget[iSpan]*(avgNormalVelDonor[SpanLevelDonor[iSpan] + 1]
                                            - avgNormalVelDonor[SpanLevelDonor[iSpan]]);
      avgNormalVelTarget[iSpan]          += avgNormalVelDonor[SpanLevelDonor[iSpan]];
      avgTangVelTarget[iSpan]             = SpanValueCoeffTarget[iSpan]*(avgTangVelDonor[SpanLevelDonor[iSpan] + 1]
                                            - avgTangVelDonor[SpanLevelDonor[iSpan]]);
      avgTangVelTarget[iSpan]            += avgTangVelDonor[SpanLevelDonor[iSpan]];
      avg3DVelTarget[iSpan]               = SpanValueCoeffTarget[iSpan]*(avg3DVelDonor[SpanLevelDonor[iSpan] + 1]
                                            - avg3DVelDonor[SpanLevelDonor[iSpan]]);
      avg3DVelTarget[iSpan]              += avg3DVelDonor[SpanLevelDonor[iSpan]];
      avgNuTarget[iSpan]                  = SpanValueCoeffTarget[iSpan]*(avgNuDonor[SpanLevelDonor[iSpan] + 1]
                                            - avgNuDonor[SpanLevelDonor[iSpan]]);
      avgNuTarget[iSpan]                 += avgNuDonor[SpanLevelDonor[iSpan]];
      avgKineTarget[iSpan]                = SpanValueCoeffTarget[iSpan]*(avgKineDonor[SpanLevelDonor[iSpan] + 1]
                                            - avgKineDonor[SpanLevelDonor[iSpan]]);
      avgKineTarget[iSpan]               += avgKineDonor[SpanLevelDonor[iSpan]];
      avgOmegaTarget[iSpan]               = SpanValueCoeffTarget[iSpan]*(avgOmegaDonor[SpanLevelDonor[iSpan] + 1]
                                            - avgOmegaDonor[SpanLevelDonor[iSpan] ]);
      avgOmegaTarget[iSpan]              += avgOmegaDonor[SpanLevelDonor[iSpan]];
    }


    /*--- transfer values at the hub ---*/
    avgDensityTarget[0]                      = avgDensityDonor[0];
    avgPressureTarget[0]                     = avgPressureDonor[0];
    avgNormalVelTarget[0]                    = avgNormalVelDonor[0];
    avgTangVelTarget[0]                      = avgTangVelDonor[0];
    avg3DVelTarget[0]                        = avg3DVelDonor[0];
    avgNuTarget[0]                           = avgNuDonor[0];
    avgKineTarget[0]                         = avgKineDonor[0];
    avgOmegaTarget[0]                        = avgOmegaDonor[0];

    /*--- transfer values at the shroud ---*/
    avgDensityTarget[nSpanTarget - 2]        = avgDensityDonor[nSpanDonor - 2];
    avgPressureTarget[nSpanTarget - 2]       = avgPressureDonor[nSpanDonor - 2];
    avgNormalVelTarget[nSpanTarget - 2]      = avgNormalVelDonor[nSpanDonor - 2];
    avgTangVelTarget[nSpanTarget - 2]        = avgTangVelDonor[nSpanDonor - 2];
    avg3DVelTarget[nSpanTarget - 2]          = avg3DVelDonor[nSpanDonor - 2];
    avgNuTarget[nSpanTarget - 2]             = avgNuDonor[nSpanDonor - 2];
    avgKineTarget[nSpanTarget - 2]           = avgKineDonor[nSpanDonor - 2];
    avgOmegaTarget[nSpanTarget - 2]          = avgOmegaDonor[nSpanDonor - 2];

    /*--- transfer 1D values ---*/
    avgDensityTarget[nSpanTarget - 1]        = avgDensityDonor[nSpanDonor - 1];
    avgPressureTarget[nSpanTarget - 1]       = avgPressureDonor[nSpanDonor - 1];
    avgNormalVelTarget[nSpanTarget - 1]      = avgNormalVelDonor[nSpanDonor - 1];
    avgTangVelTarget[nSpanTarget - 1]        = avgTangVelDonor[nSpanDonor - 1];
    avg3DVelTarget[nSpanTarget - 1]          = avg3DVelDonor[nSpanDonor - 1];
    avgNuTarget[nSpanTarget - 1]             = avgNuDonor[nSpanDonor - 1];
    avgKineTarget[nSpanTarget - 1]           = avgKineDonor[nSpanDonor - 1];
    avgOmegaTarget[nSpanTarget - 1]          = avgOmegaDonor[nSpanDonor - 1];


    /*---finally, the interpolated value is sent  to the target zone ---*/
    for(iSpan = 0; iSpan < nSpanTarget ; iSpan++){
      Target_Variable[0]                     = avgDensityTarget[iSpan];
      Target_Variable[1]                     = avgPressureTarget[iSpan];
      Target_Variable[2]                     = avgNormalVelTarget[iSpan];
      Target_Variable[3]                     = avgTangVelTarget[iSpan];
      Target_Variable[4]                     = avg3DVelTarget[iSpan];
      Target_Variable[5]                     = avgNuTarget[iSpan];
      Target_Variable[6]                     = avgKineTarget[iSpan];
      Target_Variable[7]                     = avgOmegaTarget[iSpan];


      SetTarget_Variable(target_solution, target_geometry, target_config, Marker_Target, iSpan, rank);
    }
  }

  delete [] avgDensityDonor;
  delete [] avgPressureDonor;
  delete [] avgNormalVelDonor;
  delete [] avgTangVelDonor;
  delete [] avg3DVelDonor;
  delete [] avgNuDonor;
  delete [] avgKineDonor;
  delete [] avgOmegaDonor;


  delete [] avgDensityTarget;
  delete [] avgPressureTarget;
  delete [] avgNormalVelTarget;
  delete [] avgTangVelTarget;
  delete [] avg3DVelTarget;
  delete [] avgNuTarget;
  delete [] avgKineTarget;
  delete [] avgOmegaTarget;


}

void CInterface::GatherAverageValues(CSolver *donor_solution, CSolver *target_solution, unsigned short donorZone){


  /*--- here we made the strong assumption that the mesh zone order
   * follows the same order of the turbomachinery markers ---*/
  SetAverageValues(donor_solution, target_solution, donorZone);

}

void CInterface::GatherAverageTurboGeoValues(CGeometry *donor_geometry, CGeometry *target_geometry,
                                             unsigned short donorZone){


  /*--- here we made the strong assumption that the mesh zone order
   * follows the same order of the turbomachinery markers ---*/
  SetAverageTurboGeoValues(donor_geometry, target_geometry, donorZone);

}
