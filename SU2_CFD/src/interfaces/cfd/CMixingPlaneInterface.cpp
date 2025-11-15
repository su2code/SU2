/*!
 * \file CMixingPlaneInterface.cpp
 * \brief Declaration and inlines of the class to transfer average variables
 *        needed for MixingPlane computation from a generic zone into another one.
 * \author S. Vitale
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

#include "../../../include/interfaces/cfd/CMixingPlaneInterface.hpp"
#include "../../../Common/include/interface_interpolation/CInterpolator.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

CMixingPlaneInterface::CMixingPlaneInterface(unsigned short val_nVar, unsigned short val_nConst){
  nVar = val_nVar;
  Donor_Variable     = new su2double[nVar + 5]();
  Target_Variable    = new su2double[nVar + 5]();
  InterfaceType      = ENUM_TRANSFER::MIXING_PLANE;
}

void CMixingPlaneInterface::SetSpanWiseLevels(const CConfig *donor_config, const CConfig *target_config){

  unsigned short iSpan;
  nSpanMaxAllZones = donor_config->GetnSpanMaxAllZones();


  SpanValueCoeffTarget = new su2double[target_config->GetnSpanWiseSections() + 1];
  SpanLevelDonor       = new unsigned short[target_config->GetnSpanWiseSections() + 1];


  for (iSpan = 0; iSpan < target_config->GetnSpanWiseSections() + 1;iSpan++){
    SpanValueCoeffTarget[iSpan] = 0.0;
    SpanLevelDonor[iSpan]       = 1;
  }
}

void CMixingPlaneInterface::BroadcastData_MixingPlane(const CInterpolator& interpolator,
                               CSolver *donor_solution, CSolver *target_solution,
                               CGeometry *donor_geometry, CGeometry *target_geometry,
                               const CConfig *donor_config, const CConfig *target_config) {
  static_assert(su2activematrix::Storage == StorageType::RowMajor,"");

  /*--- Loop over interface markers. ---*/

  for (auto iMarkerInt = 0u; iMarkerInt < donor_config->GetMarker_n_ZoneInterface()/2; iMarkerInt++) {

    /*--- Check if this interface connects the two zones, if not continue. ---*/

    const auto markDonor = donor_config->FindInterfaceMarker(iMarkerInt);
    const auto markTarget = target_config->FindInterfaceMarker(iMarkerInt);

    if(!CInterpolator::CheckInterfaceBoundary(markDonor, markTarget)) continue;

    // The number of spans is available on every rank
    const auto nSpanDonor = donor_config->GetnSpanWiseSections();

    /*--- Fill send buffers. ---*/

    vector<unsigned long> sendDonorIdx(nSpanDonor);
    su2activematrix sendDonorVar(nSpanDonor, 8);

    if (markDonor >= 0) {
      for (auto iSpan = 0ul, iSend = 0ul; iSpan < nSpanDonor; iSpan++) {
        GetDonor_Variable(donor_solution, donor_geometry, donor_config, markDonor, iSpan, 0);
        for (auto iVar = 0u; iVar < nVar; iVar++) sendDonorVar(iSend, iVar) = Donor_Variable[iVar];
        sendDonorIdx[iSend] = iSpan;
        ++iSend;
      }
    }
    /*--- Gather data. ---*/
    const int nTotalDonors = nSpanDonor * size;
    const int nSpanDonorVars = nSpanDonor * 8;
    vector<unsigned long> donorIdx(nTotalDonors);
    su2activematrix donorVar(nTotalDonors, 8);

    SU2_MPI::Allgather(sendDonorIdx.data(), nSpanDonor, MPI_UNSIGNED_LONG,
                      donorIdx.data(), nSpanDonor, MPI_UNSIGNED_LONG,
                      SU2_MPI::GetComm());

    SU2_MPI::Allgather(sendDonorVar.data(), nSpanDonorVars, MPI_DOUBLE,
                      donorVar.data(), nSpanDonorVars, MPI_DOUBLE,
                      SU2_MPI::GetComm());

    /*--- This rank does not need to do more work. ---*/
    if (markTarget < 0) continue;

    /*--- Loop over target spans. ---*/
    auto nTargetSpan = target_config->GetnSpanWiseSections() + 1;

    for (auto iTargetSpan = 0ul; iTargetSpan < nTargetSpan; iTargetSpan++) {

      auto& targetSpan = interpolator.targetSpans[iMarkerInt][markTarget];
      const auto nDonorSpan = targetSpan.nDonor();

      InitializeTarget_Variable(target_solution, markTarget, iTargetSpan, nDonorSpan);

      if ((iTargetSpan == 0) || (iTargetSpan < nTargetSpan - 3)) {
        /*--- Transfer values at hub, shroud and 1D values ---*/
        unsigned long iDonorSpan;
        if (iTargetSpan == 0) iDonorSpan = 0;
        else if (iTargetSpan == nTargetSpan - 2) iDonorSpan = nDonorSpan - 2;
        else if (iTargetSpan == nTargetSpan - 1) iDonorSpan = nDonorSpan - 1;

        const auto idx = lower_bound(donorIdx.begin(), donorIdx.end(), iDonorSpan) - donorIdx.begin();
        assert(idx < static_cast<long>(donorIdx.size()));

        RecoverTarget_Variable(donorVar[idx]);

        SetTarget_Variable(target_solution, target_geometry, target_config, markTarget, iTargetSpan, 0);
      }
      else {
        /*--- Get the global index of the donor and the interpolation coefficient. ---*/

        const auto donorSpan = targetSpan.globalSpan[iTargetSpan];
        const auto donorCoeff = targetSpan.coefficient[iTargetSpan];

        /*--- Find the index of the global donor point in the donor data. ---*/

        const auto idx = lower_bound(donorIdx.begin(), donorIdx.end(), donorSpan) - donorIdx.begin();
        assert(idx < static_cast<long>(donorIdx.size()));

        /*--- Recover the Target_Variable from the buffer of variables. ---*/
        RecoverTarget_Variable(donorVar[idx], donorVar[idx+1], donorCoeff);

        SetTarget_Variable(target_solution, target_geometry, target_config, markTarget, iTargetSpan, 0);
      }
    }
  }
}

void CMixingPlaneInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                                          const CConfig *donor_config, unsigned long Marker_Donor,
                                          unsigned long Span_Donor, unsigned long Point_Donor) {

  Donor_Variable[0] = donor_solution->GetAverageDensity(Marker_Donor, Span_Donor);
  Donor_Variable[1] = donor_solution->GetAveragePressure(Marker_Donor, Span_Donor);
  Donor_Variable[2] = donor_solution->GetAverageTurboVelocity(Marker_Donor, Span_Donor)[0];
  Donor_Variable[3] = donor_solution->GetAverageTurboVelocity(Marker_Donor, Span_Donor)[1];

  if(donor_geometry->GetnDim() == 3){
    Donor_Variable[4] = donor_solution->GetAverageTurboVelocity(Marker_Donor, Span_Donor)[2];
  }
  else{
    Donor_Variable[4] = -1.0;
  }

  if(donor_config->GetKind_Turb_Model() != TURB_MODEL::NONE){
    Donor_Variable[5] = donor_solution->GetAverageNu(Marker_Donor, Span_Donor);
    Donor_Variable[6] = donor_solution->GetAverageKine(Marker_Donor, Span_Donor);
    Donor_Variable[7] = donor_solution->GetAverageOmega(Marker_Donor, Span_Donor);
  }
  else{
    Donor_Variable[5] = -1.0;
    Donor_Variable[6] = -1.0;
    Donor_Variable[7] = -1.0;
  }
}

void CMixingPlaneInterface::InitializeTarget_Variable(CSolver *target_solution, unsigned long Marker_Target,
                                                  unsigned long Span_Target, unsigned short nSpanDonor) {

  target_solution->SetnMixingStates(Marker_Target, Span_Target, nSpanDonor); // This is to allocate
  target_solution->SetMixingStateStructure(Marker_Target, Span_Target);
  target_solution->SetnMixingStates(Marker_Target, Span_Target, 0); // Reset counter to 0

}

void CMixingPlaneInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                                           const CConfig *target_config, unsigned long Marker_Target,
                                           unsigned long Span_Target, unsigned long Point_Target) {

  unsigned short iVar, iDonorSpan, nTargetVar;
  nTargetVar = 8;
  /*--- Set the mixing plane solution with the value of the Target Variable ---*/

  iDonorSpan = target_solution->GetnMixingStates(Marker_Target, Span_Target);

  for (iVar = 0; iVar < nTargetVar; iVar++)
    target_solution->SetMixingState(Marker_Target, Span_Target, iVar, Target_Variable[iVar]);

  target_solution->SetnMixingStates( Marker_Target, Span_Target, iDonorSpan + 1 );
}

void CMixingPlaneInterface::SetAverageValues(CSolver *donor_solution, CSolver *target_solution,
                                             unsigned short donorZone){
  
  unsigned short iSpan;

  for(iSpan = 0; iSpan<nSpanMaxAllZones +1; iSpan++){
    /*--- transfer inviscid quantities ---*/
    target_solution->SetDensityIn(donor_solution->GetDensityIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetPressureIn(donor_solution->GetPressureIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetTurboVelocityIn(donor_solution->GetTurboVelocityIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetDensityOut(donor_solution->GetDensityOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetPressureOut(donor_solution->GetPressureOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetTurboVelocityOut(donor_solution->GetTurboVelocityOut(donorZone, iSpan), donorZone, iSpan);

    /*--- transfer turbulent quantities ---*/
    target_solution->SetKineIn(donor_solution->GetKineIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetOmegaIn(donor_solution->GetOmegaIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetNuIn(donor_solution->GetNuIn(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetKineOut(donor_solution->GetKineOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetOmegaOut(donor_solution->GetOmegaOut(donorZone, iSpan), donorZone, iSpan);
    target_solution->SetNuOut(donor_solution->GetNuOut(donorZone, iSpan), donorZone, iSpan);

  }
}
